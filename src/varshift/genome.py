"""Random access to reference FASTA files through a memory map.

Only the pages that are read come into memory, which is what lets a whole genome be used on a
machine with less RAM than the genome is long.

A FASTA record named like ``chr7:140400001-140630000`` is treated as a window of that contig. The
bundled demo data uses this, so it can hold real sequence at real coordinates without 3 GB of genome.
"""

from __future__ import annotations

import mmap
import os
import re
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

from .errors import DataError

_WINDOW_RE = re.compile(r"^(.+):(\d+)-(\d+)$")


@dataclass(frozen=True)
class _Segment:
    start: int  # 0-based position of the first base on the contig
    length: int
    offset: int  # byte offset of the first base in the file
    line_bases: int
    line_bytes: int


class FastaFile:
    """One FASTA file plus its .fai index."""

    def __init__(self, fasta: Path):
        self.path = Path(fasta)
        index = Path(str(fasta) + ".fai")
        if not self.path.is_file() or not index.is_file():
            raise DataError(f"reference {self.path} or its index is missing; run: varshift setup")
        self._segments: dict[str, list[_Segment]] = {}
        self._starts: dict[str, list[int]] = {}
        self._load_index(index)
        self._file = open(self.path, "rb")
        self._map = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

    def _load_index(self, index: Path) -> None:
        with open(index, encoding="ascii") as handle:
            for number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise DataError(f"{index}: line {number} is not a valid .fai record")
                name, length, offset, line_bases, line_bytes = fields[:5]
                start = 0
                window = _WINDOW_RE.match(name)
                if window:
                    name = window.group(1)
                    start = int(window.group(2)) - 1
                segment = _Segment(start, int(length), int(offset), int(line_bases), int(line_bytes))
                self._segments.setdefault(name, []).append(segment)
        for name, segments in self._segments.items():
            segments.sort(key=lambda s: s.start)
            self._starts[name] = [s.start for s in segments]

    def has(self, contig: str) -> bool:
        return contig in self._segments

    def contigs(self) -> list[str]:
        return list(self._segments)

    def covers(self, contig: str, start: int, end: int) -> bool:
        segments = self._segments.get(contig)
        if not segments or start < 0 or end < start:
            return False
        index = bisect_right(self._starts[contig], start) - 1
        if index < 0:
            return False
        segment = segments[index]
        return end <= segment.start + segment.length

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Uppercase bases for the 0-based half-open range, clipped to the sequence ends."""
        segments = self._segments.get(contig)
        if not segments:
            raise DataError(f"contig {contig} is not in {self.path.name}")
        index = bisect_right(self._starts[contig], max(start, 0)) - 1
        if index < 0:
            return ""
        segment = segments[index]
        begin = max(start - segment.start, 0)
        stop = min(end - segment.start, segment.length)
        if begin >= stop:
            return ""
        first = segment.offset + begin // segment.line_bases * segment.line_bytes + begin % segment.line_bases
        last = (
            segment.offset
            + (stop - 1) // segment.line_bases * segment.line_bytes
            + (stop - 1) % segment.line_bases
        )
        raw = self._map[first : last + 1]
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()

    def close(self) -> None:
        self._map.close()
        self._file.close()


class Genome:
    """Reference bases for one build, including the rCRS special case on hg19."""

    def __init__(self, build: str, fasta: Path, rcrs: FastaFile | None = None, owns_rcrs: bool = False):
        self.build = build
        self._fasta = FastaFile(fasta)
        self._rcrs = rcrs
        self._owns_rcrs = owns_rcrs

    def has(self, contig: str) -> bool:
        if contig == "chrMT":
            return self._rcrs is not None and self._rcrs.has("chrM")
        return self._fasta.has(contig)

    def covers(self, contig: str, start: int, end: int) -> bool:
        if contig == "chrMT":
            return self._rcrs is not None and self._rcrs.covers("chrM", start, end)
        return self._fasta.covers(contig, start, end)

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig == "chrMT":
            if self._rcrs is None:
                raise DataError("the rCRS sequence needs the hg38 reference; run: varshift setup build=hg38")
            return self._rcrs.fetch("chrM", start, end)
        return self._fasta.fetch(contig, start, end)

    def close(self) -> None:
        self._fasta.close()
        if self._owns_rcrs and self._rcrs is not None:
            self._rcrs.close()


def open_build_genome(genome_dir: Path, build: str) -> Genome:
    """Open one build for a builder. On hg19 the hg38 file, if present, supplies the rCRS sequence."""
    rcrs = None
    hg38 = Path(genome_dir) / "hg38.fa"
    if build == "hg19" and hg38.is_file():
        rcrs = FastaFile(hg38)
    return Genome(build, Path(genome_dir) / f"{build}.fa", rcrs, owns_rcrs=True)


def write_fai(fasta: Path) -> int:
    """Write the .fai index for a FASTA file and return the number of records.

    The file is scanned in binary with a small buffer, so memory stays flat however large it is.
    """
    index_path = Path(str(fasta) + ".fai")
    tmp_path = index_path.with_suffix(".fai.tmp")
    records = 0
    with open(fasta, "rb") as handle, open(tmp_path, "w", encoding="ascii", newline="\n") as out:
        name = None
        length = offset = line_bases = line_bytes = 0
        position = 0
        short_line_seen = False
        for line in handle:
            size = len(line)
            if line.startswith(b">"):
                if name is not None:
                    out.write(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_bytes}\n")
                    records += 1
                name = line[1:].split(None, 1)[0].decode("ascii")
                length = line_bases = line_bytes = 0
                offset = position + size
                short_line_seen = False
            elif name is None:
                raise DataError(f"{fasta}: sequence found before the first header")
            else:
                bases = len(line.rstrip(b"\r\n"))
                if bases:
                    if short_line_seen:
                        raise DataError(f"{fasta}: uneven line length in {name}, cannot index")
                    if line_bases == 0:
                        line_bases, line_bytes = bases, size
                    elif bases > line_bases:
                        raise DataError(f"{fasta}: uneven line length in {name}, cannot index")
                    elif bases < line_bases:
                        short_line_seen = True
                    length += bases
            position += size
        if name is not None:
            out.write(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_bytes}\n")
            records += 1
    os.replace(tmp_path, index_path)
    return records
