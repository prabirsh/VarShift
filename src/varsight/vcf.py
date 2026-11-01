"""Streaming VCF reader. Nothing is loaded whole, so file size does not change memory use."""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from .builds import detect_build_from_lengths, parse_build
from .errors import InputError

MAX_LINE_BYTES = 64 * 1024 * 1024
_CONTIG_LINE_RE = re.compile(r"^##contig=<ID=([^,>]+)(?:,.*?length=(\d+))?", re.IGNORECASE)
_REFERENCE_RE = re.compile(r"^##reference=(.*)$", re.IGNORECASE)
_BUILD_WORD_RE = re.compile(r"(hg19|hg38|grch37|grch38|b37|b38)", re.IGNORECASE)
_PLAIN_ALLELE_RE = re.compile(r"^[ACGTNacgtn]+$")


@dataclass
class VcfHeader:
    lines: list[str]
    columns: list[str]
    contig_lengths: dict[str, int]
    build: str | None
    wrong_contigs: list[str]
    named_build: str | None


@dataclass
class VcfRecord:
    line_number: int
    fields: list[str]
    alt: str
    allele_index: int
    allele_count: int

    @property
    def chrom(self) -> str:
        return self.fields[0]

    @property
    def pos(self) -> int:
        return int(self.fields[1])

    @property
    def ref(self) -> str:
        return self.fields[3]

    @property
    def is_plain(self) -> bool:
        return bool(_PLAIN_ALLELE_RE.match(self.alt)) and bool(_PLAIN_ALLELE_RE.match(self.ref))


def _open_text(path: Path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="\n")
    return open(path, encoding="utf-8", errors="replace", newline="\n")


class VcfReader:
    def __init__(self, path: Path):
        self.path = Path(path)
        self._handle = _open_text(self.path)
        self._line_number = 0
        self.header = self._read_header()

    def _readline(self) -> str:
        line = self._handle.readline(MAX_LINE_BYTES)
        if line and not line.endswith("\n") and len(line) >= MAX_LINE_BYTES:
            raise InputError(f"{self.path}: line {self._line_number + 1} is longer than {MAX_LINE_BYTES} bytes")
        self._line_number += 1
        return line

    def _read_header(self) -> VcfHeader:
        lines: list[str] = []
        lengths: dict[str, int] = {}
        named = None
        while True:
            line = self._readline()
            if not line:
                raise InputError(f"{self.path} ended before a #CHROM line; is it a VCF?")
            line = line.rstrip("\r\n")
            if line.startswith("##"):
                lines.append(line)
                match = _CONTIG_LINE_RE.match(line)
                if match and match.group(2):
                    lengths[match.group(1)] = int(match.group(2))
                match = _REFERENCE_RE.match(line)
                if match:
                    word = _BUILD_WORD_RE.search(match.group(1))
                    if word:
                        named = parse_build(word.group(1))
                continue
            if line.startswith("#CHROM"):
                columns = line[1:].split("\t")
                if columns[:5] != ["CHROM", "POS", "ID", "REF", "ALT"]:
                    raise InputError(f"{self.path}: the #CHROM line does not have the standard VCF columns")
                build, wrong = detect_build_from_lengths(lengths)
                return VcfHeader(lines, columns, lengths, build, wrong, named)
            raise InputError(f"{self.path}: expected a header line at line {self._line_number}")

    def records(self) -> Iterator[VcfRecord]:
        """Yield one record per ALT allele."""
        while True:
            line = self._readline()
            if not line:
                return
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise InputError(f"{self.path}: line {self._line_number} has fewer than 5 columns")
            if not fields[1].isdigit():
                raise InputError(f"{self.path}: line {self._line_number} has a non-numeric POS")
            alts = fields[4].split(",")
            for index, alt in enumerate(alts, 1):
                yield VcfRecord(self._line_number, fields, alt, index, len(alts))

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "VcfReader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def looks_like_vcf(text: str) -> bool:
    lowered = text.lower()
    return lowered.endswith((".vcf", ".vcf.gz", ".vcf.bgz")) and Path(text).is_file()
