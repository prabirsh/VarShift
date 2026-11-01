"""Annotate a VCF or a plain list of variants, one record at a time."""

from __future__ import annotations

import logging
import time
from collections import Counter
from pathlib import Path

from .annotate import Annotator
from .builds import canonical_contig, parse_build
from .errors import InputError, ReferenceMismatch, VarSightError
from .output import ResultWriter
from .resolver import Candidate
from .variant import normalize
from .vcf import VcfReader

log = logging.getLogger("varsight")


def _error_record(raw: str, kind: str, build: str | None, message: str, flag: str) -> dict:
    return {
        "input": {"raw": raw, "type": kind, "build": build},
        "candidate": {"rank": None, "of": 0},
        "flags": [flag],
        "error": message,
    }


def annotate_vcf(annotator: Annotator, path: Path, writer: ResultWriter, build: str | None) -> Counter:
    counts: Counter = Counter()
    started = time.time()
    with VcfReader(path) as reader:
        chosen = _vcf_build(reader, build)
        columns = reader.header.columns
        for record in reader.records():
            counts["alleles"] += 1
            prefix = dict(zip(columns, record.fields))
            prefix["ALT"] = record.alt
            prefix["ALT_ORIGINAL"] = record.fields[4]
            raw = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
            if not record.is_plain:
                counts["symbolic_skipped"] += 1
                writer.write(
                    _error_record(raw, "vcf", chosen, f"allele {record.alt} is not a plain sequence", "symbolic_allele_skipped"),
                    prefix,
                )
                continue
            try:
                contig = canonical_contig(record.chrom, chosen)
                variant = normalize(annotator.store.genome(chosen), contig, record.pos, record.ref, record.alt)
                result = annotator.annotate_candidate(Candidate(variant, via="vcf"))
                result["input"] = {"raw": raw, "type": "vcf", "build": chosen}
                if record.allele_count > 1:
                    result["flags"].append("split_from_multiallelic")
                writer.write(result, prefix)
                counts["annotated"] += 1
                if result["evidence"]["clinvar"].get("status") == "found":
                    counts["in_clinvar"] += 1
            except ReferenceMismatch as exc:
                counts["ref_mismatch"] += 1
                writer.write(_error_record(raw, "vcf", chosen, str(exc), "ref_mismatch"), prefix)
            except VarSightError as exc:
                counts["not_annotated"] += 1
                writer.write(_error_record(raw, "vcf", chosen, str(exc), "not_annotated"), prefix)
    counts["seconds"] = int(time.time() - started)
    return counts


def _vcf_build(reader: VcfReader, build: str | None) -> str:
    header = reader.header
    chosen = parse_build(build) if build else None
    if chosen and header.build and chosen != header.build:
        raise InputError(
            f"--build says {chosen}, but the contig lengths in {reader.path.name} are {header.build}. "
            "Annotating against the wrong build gives wrong results, so VarSight stopped."
        )
    chosen = chosen or header.build or header.named_build
    if chosen is None:
        raise InputError(f"cannot tell the build of {reader.path.name} from its header; add --build hg19 or --build hg38")
    if header.build and header.named_build and header.build != header.named_build:
        log.warning("the header names %s but the contig lengths are %s; using %s", header.named_build, header.build, chosen)
    if header.wrong_contigs:
        names = ", ".join(header.wrong_contigs[:5])
        log.warning("these contigs have unexpected lengths for %s: %s", chosen, names)
    return chosen


def vcf_columns(path: Path) -> list[str]:
    with VcfReader(path) as reader:
        return reader.header.columns + ["ALT_ORIGINAL"]


def annotate_list(annotator: Annotator, path: Path, writer: ResultWriter, build: str | None, gene: str = "") -> Counter:
    """One description per line. Blank lines and lines starting with # are skipped."""
    counts: Counter = Counter()
    with open(path, encoding="utf-8", errors="replace") as handle:
        for number, line in enumerate(handle, 1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            counts["lines"] += 1
            try:
                for record in annotator.annotate_text(text, build, gene):
                    writer.write(record)
                    counts["results"] += 1
            except VarSightError as exc:
                counts["errors"] += 1
                writer.write(_error_record(text, "list", build, str(exc), "not_annotated"))
                log.warning("line %d: %s", number, exc)
    return counts
