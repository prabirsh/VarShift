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
from .errors import DataError, InputError, ReferenceMismatch, VarSightError
from .resolver import Candidate
from .store import DataStore
from .variant import normalize
from .vcf import VcfReader

log = logging.getLogger("varsight")
_CHUNK = 200


def _error_record(raw: str, kind: str, build: str | None, message: str, flag: str) -> dict:
    return {
        "input": {"raw": raw, "type": kind, "build": build},
        "candidate": {"rank": None, "of": 0},
        "flags": [flag],
        "error": message,
    }


def _annotate_record(annotator: Annotator, record, columns: list[str], build: str) -> tuple[dict, dict, str]:
    """One VCF allele. Returns the result, the original columns, and how it went."""
    prefix = dict(zip(columns, record.fields))
    prefix["ALT"] = record.alt
    prefix["ALT_ORIGINAL"] = record.fields[4]
    raw = f"{record.chrom}:{record.pos}:{record.ref}:{record.alt}"
    if not record.is_plain:
        return _error_record(raw, "vcf", build, f"allele {record.alt} is not a plain sequence", "symbolic_allele_skipped"), prefix, "symbolic_skipped"
    try:
        contig = canonical_contig(record.chrom, build)
        variant = normalize(annotator.store.genome(build), contig, record.pos, record.ref, record.alt)
        result = annotator.annotate_candidate(Candidate(variant, via="vcf"))
    except ReferenceMismatch as exc:
        return _error_record(raw, "vcf", build, str(exc), "ref_mismatch"), prefix, "ref_mismatch"
    except VarSightError as exc:
        return _error_record(raw, "vcf", build, str(exc), "not_annotated"), prefix, "not_annotated"
    result["input"] = {"raw": raw, "type": "vcf", "build": build}
    if record.allele_count > 1:
        result["flags"].append("split_from_multiallelic")
    return result, prefix, "annotated"


_worker: dict = {}


def _start_worker(data_dir: str, all_transcripts: bool, columns: list[str], build: str) -> None:
    try:
        _worker.update(annotator=Annotator(DataStore(Path(data_dir)), all_transcripts), columns=columns, build=build)
    except VarSightError as exc:
        # An error raised here would make the pool restart this worker for ever, so report it with the first chunk.
        _worker["failed"] = str(exc)


def _annotate_chunk(records: list) -> list:
    if "failed" in _worker:
        raise DataError(_worker["failed"])
    return [_annotate_record(_worker["annotator"], record, _worker["columns"], _worker["build"]) for record in records]


def _in_parallel(annotator: Annotator, records, columns: list[str], build: str, workers: int):
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

    settings = (str(annotator.store.data_dir), annotator.all_transcripts, columns, build)
    # spawn on every system: a forked child would inherit open SQLite connections, which isn't safe.
    with multiprocessing.get_context("spawn").Pool(workers, _start_worker, settings) as pool:
        try:
            for done in pool.imap(_annotate_chunk, chunks()):
                gate.release()
                yield from done
        finally:
            stopping.set()
            gate.release(waiting + 1)


def annotate_vcf(annotator: Annotator, path: Path, writer, build: str | None, workers: int = 1) -> Counter:
    counts: Counter = Counter()
    started = time.time()
    with VcfReader(path) as reader:
        chosen = _vcf_build(reader, build)
        columns = reader.header.columns
        if workers > 1:
            results = _in_parallel(annotator, reader.records(), columns, chosen, workers)
        else:
            results = (_annotate_record(annotator, record, columns, chosen) for record in reader.records())
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
            "Annotating against the wrong build gives wrong results, so VarSight stopped."
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


def annotate_list(annotator: Annotator, path: Path, writer, build: str | None, gene: str = "") -> Counter:
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
