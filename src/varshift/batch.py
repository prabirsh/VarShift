"""Annotate a VCF or a plain list of variants, one record at a time."""

from __future__ import annotations

import logging
import multiprocessing
import threading
import time
from collections import Counter
from pathlib import Path

from .annotate import Annotator
from .builds import canonical_contig, parse_build
from .errors import DataError, InputError, ReferenceMismatch, VarShiftError
from .resolver import Candidate
from .store import DataStore
from .tidy import tidy
from .variant import normalize
from .vcf import VcfReader

log = logging.getLogger("varshift")
_CHUNK = 200


def _error_record(raw: str, kind: str, build: str | None, message: str, flag: str, extra: dict | None = None) -> dict:
    return {
        "input": {"raw": raw, "type": kind, "build": build, **(extra or {})},
        "candidate": {"rank": None, "of": 0},
        "flags": [flag],
        "error": message,
    }


def _annotate_record(annotator: Annotator, record, columns: list[str], build: str, mode: str | None = None) -> tuple[dict, dict, str]:
    """One VCF allele. Returns the result, the original columns, and how it went."""
    prefix = dict(zip(columns, record.fields))
    prefix["ALT"] = record.alt
    prefix["ALT_ORIGINAL"] = record.fields[4]
    raw = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
    extra = {"mode": mode} if mode else None
    if not record.is_plain:
        message = f"allele {record.alt} is not a plain sequence"
        return _error_record(raw, "vcf", build, message, "symbolic_allele_skipped", extra), prefix, "symbolic_skipped"
    try:
        contig = canonical_contig(record.chrom, build)
        variant = normalize(annotator.store.genome(build), contig, record.pos, record.ref, record.alt)
        result = annotator.annotate_candidate(Candidate(variant, via="vcf"))
    except ReferenceMismatch as exc:
        return _error_record(raw, "vcf", build, str(exc), "ref_mismatch", extra), prefix, "ref_mismatch"
    except VarShiftError as exc:
        return _error_record(raw, "vcf", build, str(exc), "not_annotated", extra), prefix, "not_annotated"
    result["input"] = {"raw": raw, "type": "vcf", "build": build, **(extra or {})}
    if record.allele_count > 1:
        result["flags"].append("split_from_multiallelic")
    return result, prefix, "annotated"


_worker: dict = {}


def _start_worker(data_dir: str, all_transcripts: bool, columns: list[str], build: str, mode: str | None = None) -> None:
    try:
        _worker.update(annotator=Annotator(DataStore(Path(data_dir)), all_transcripts), columns=columns, build=build, mode=mode)
    except VarShiftError as exc:
        # An error raised here would make the pool restart this worker for ever, so report it with the first chunk.
        _worker["failed"] = str(exc)


def _annotate_chunk(records: list) -> list:
    if "failed" in _worker:
        raise DataError(_worker["failed"])
    return [_annotate_record(_worker["annotator"], record, _worker["columns"], _worker["build"], _worker["mode"]) for record in records]


def _in_parallel(annotator: Annotator, records, columns: list[str], build: str, workers: int, mode: str | None = None):
    """Annotate in worker processes, each with its own read-only connections. Results keep the input order."""
    waiting = workers * 4
    # The reader runs ahead on its own thread, so this stops a huge VCF piling up in memory.
    gate = threading.Semaphore(waiting)
    stopping = threading.Event()

    def chunks():
        chunk = []
        for record in records:
            chunk.append(record)
            if len(chunk) == _CHUNK:
                gate.acquire()
                if stopping.is_set():
                    return
                yield chunk
                chunk = []
        if chunk:
            gate.acquire()
            if not stopping.is_set():
                yield chunk

    settings = (str(annotator.store.data_dir), annotator.all_transcripts, columns, build, mode)
    # spawn on every system: a forked child would inherit open SQLite connections, which isn't safe.
    with multiprocessing.get_context("spawn").Pool(workers, _start_worker, settings) as pool:
        try:
            for done in pool.imap(_annotate_chunk, chunks()):
                gate.release()
                yield from done
        finally:
            stopping.set()
            gate.release(waiting + 1)


def annotate_vcf(annotator: Annotator, path: Path, writer, build: str | None, workers: int = 1, mode: str | None = None) -> Counter:
    counts: Counter = Counter()
    started = time.time()
    with VcfReader(path) as reader:
        chosen = _vcf_build(reader, build)
        columns = reader.header.columns
        if workers > 1:
            results = _in_parallel(annotator, reader.records(), columns, chosen, workers, mode)
        else:
            results = (_annotate_record(annotator, record, columns, chosen, mode) for record in reader.records())
        try:
            for result, prefix, outcome in results:
                writer.write(result, prefix)
                counts["alleles"] += 1
                counts[outcome] += 1
                if outcome == "annotated" and result["evidence"]["clinvar"].get("status") == "found":
                    counts["in_clinvar"] += 1
        finally:
            results.close()
    counts["seconds"] = int(time.time() - started)
    return counts


def _vcf_build(reader: VcfReader, build: str | None) -> str:
    header = reader.header
    chosen = parse_build(build) if build else None
    if chosen and header.build and chosen != header.build:
        raise InputError(
            f"build= says {chosen}, but the contig lengths in {reader.path.name} are {header.build}. "
            "Annotating against the wrong build gives wrong results, so VarShift stopped."
        )
    if chosen and not header.build and header.named_build and chosen != header.named_build:
        raise InputError(
            f"build= says {chosen}, but the ##reference line of {reader.path.name} names {header.named_build}. "
            "Annotating against the wrong build gives wrong results, so VarShift stopped."
        )
    chosen = chosen or header.build or header.named_build
    if chosen is None:
        raise InputError(f"cannot tell the build of {reader.path.name} from its header; add build=hg19 or build=hg38")
    if header.build and header.named_build and header.build != header.named_build:
        log.warning("the header names %s but the contig lengths are %s; using %s", header.named_build, header.build, chosen)
    if header.wrong_contigs:
        names = ", ".join(header.wrong_contigs[:5])
        log.warning("these contigs have unexpected lengths for %s: %s", chosen, names)
    return chosen


def vcf_columns(path: Path) -> list[str]:
    with VcfReader(path) as reader:
        return reader.header.columns + ["ALT_ORIGINAL"]


def annotate_list(annotator: Annotator, path: Path, writer, build: str | None, gene: str = "", mode: str | None = None) -> Counter:
    """One description per line. Blank lines and lines starting with # are skipped."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        stripped = ((number, line.strip()) for number, line in enumerate(handle, 1))
        wanted = ((number, text) for number, text in stripped if text and not text.startswith("#"))
        return annotate_lines(annotator, wanted, writer, build, gene, mode)


def annotate_lines(
    annotator: Annotator, lines, writer, build: str | None, gene: str = "", mode: str | None = None, label: str = "line"
) -> Counter:
    """Annotate (number, description) pairs, from a list file or a ; list. A failing one becomes an error row."""
    counts: Counter = Counter()
    for number, text in lines:
        counts["lines"] += 1
        try:
            for record in annotator.annotate_text(text, build, gene, mode):
                writer.write(record)
                counts["results"] += 1
        except VarShiftError as exc:
            counts["errors"] += 1
            reading = tidy(text, gene)
            extra = {"read_as": reading.text if reading.changed else "", "read_as_notes": list(reading.notes), "mode": mode}
            writer.write(_error_record(text, "list", build, str(exc), "not_annotated", extra))
            log.warning("%s %d: %s", label, number, exc)
    return counts


def split_queries(text: str) -> list[str]:
    """Split a description on ; into queries. A ; inside square brackets belongs to an HGVS change."""
    parts, current, depth = [], [], 0
    for character in text:
        if character == "[":
            depth += 1
        elif character == "]":
            depth = max(depth - 1, 0)
        if character == ";" and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(character)
    parts.append("".join(current))
    return [part.strip() for part in parts if part.strip()]
