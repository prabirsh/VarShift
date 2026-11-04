"""Measure speed and memory on this computer, using public ClinVar records as the test VCF."""

from __future__ import annotations

import os
import platform
import shutil
import statistics
import sys
import tempfile
import time
from pathlib import Path

from .. import __version__
from ..annotate import Annotator
from ..batch import annotate_vcf, vcf_columns
from ..builds import contig_length, key_contig
from ..errors import DataError, VarShiftError
from ..output import ExcelSink, HtmlSink, JsonLinesSink, ReportContext, ReportFolder, TsvSink
from ..store import DataStore


def peak_memory_mb() -> float | None:
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t)] + [
                (f"unused{n}", ctypes.c_size_t) for n in range(7)
            ]

        counters = Counters()
        counters.cb = ctypes.sizeof(Counters)
        process = ctypes.windll.kernel32.GetCurrentProcess
        process.restype = wintypes.HANDLE
        query = ctypes.windll.psapi.GetProcessMemoryInfo
        query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if query(process(), ctypes.byref(counters), counters.cb):
            return counters.PeakWorkingSetSize / 1024**2
        return None
    import resource

    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def _Reports(folder: Path, context: ReportContext) -> ReportFolder:
    return ReportFolder(folder, context)


def _size_of(path: Path) -> int:
    return sum(item.stat().st_size for item in path.iterdir()) if path.is_dir() else path.stat().st_size


def _write_test_vcf(store: DataStore, build: str, path: Path, records: int) -> int:
    db = store._db(f"clinvar.{build}")
    if db is None:
        raise DataError(f"ClinVar for {build} is not installed; run: varshift setup")
    total = db.execute("SELECT count(*) FROM clinvar").fetchone()[0]
    step = max(total // records, 1)
    rows = db.execute("SELECT contig, pos, ref, alt FROM clinvar WHERE variation_id % ? = 0 ORDER BY contig, pos", (step,)).fetchmany(records)
    contigs = sorted({row["contig"] for row in rows})
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        out.write(f"##fileformat=VCFv4.2\n##reference={build}\n")
        for contig in contigs:
            length = contig_length(contig, build)
            if length:
                out.write(f"##contig=<ID={contig},length={length}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST\n")
        for row in rows:
            contig = f"chr{key_contig(row['contig'])}" if row["contig"] != "chrMT" else "MT"
            out.write(f"{contig}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t50\tPASS\t.\tGT\t0/1\n")
    return len(rows)


def run_benchmark(data_dir: Path, records: int = 3000) -> str:
    store = DataStore(data_dir)
    build = "hg38" if "hg38" in store.builds else store.builds[0]
    work = Path(tempfile.mkdtemp(prefix="varshift-benchmark-", dir=data_dir))
    try:
        vcf = work / "test.vcf"
        written = _write_test_vcf(store, build, vcf, records)
        annotator = Annotator(store)
        rows = []
        cores = min(os.cpu_count() or 1, 4)
        sinks = [("TSV", TsvSink, "results.tsv", 1), ("Excel", ExcelSink, "results.xlsx", 1), ("HTML report", HtmlSink, "report.html", 1), ("JSON Lines", JsonLinesSink, "results.jsonl", 1)]
        if cores > 1:
            sinks += [(f"TSV, {cores} workers", TsvSink, "workers.tsv", cores), (f"All three reports, {cores} workers", _Reports, "reports", cores)]
        for label, sink_class, file_name, workers in sinks:
            context = ReportContext(input_name="benchmark", prefix=vcf_columns(vcf), installed=store.installed(), versions=store.versions())
            sink = sink_class(work / file_name, context)
            started = time.perf_counter()
            counts = annotate_vcf(annotator, vcf, sink, build, workers)
            sink.close()
            took = time.perf_counter() - started
            rows.append((label, took, counts["alleles"] / took if took else 0, _size_of(work / file_name) / 1024**2))
        keys = [line.split("\t")[:5] for line in vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")][:300]
        timings = []
        for chrom, pos, _, ref, alt in keys:
            started = time.perf_counter()
            try:
                annotator.annotate_text(f"{chrom}-{pos}-{ref}-{alt}", build)
            except VarShiftError:  # a few ClinVar records have N in an allele, which text input refuses
                continue
            timings.append((time.perf_counter() - started) * 1000)
        memory = peak_memory_mb()
    finally:
        store.close()
        shutil.rmtree(work, ignore_errors=True)
    lines = [
        "# VarShift benchmark",
        "",
        f"VarShift {__version__} on {platform.system()} {platform.release()}, Python {platform.python_version()}, {os.cpu_count()} logical CPUs.",
        f"Data: {data_dir}. Test VCF: {written:,} public ClinVar records on {build}.",
        "",
        "| Output | Time | Records per second | File size |",
        "| - | -: | -: | -: |",
    ]
    for label, took, rate, size in rows:
        lines.append(f"| {label} | {took:.1f} s | {rate:,.0f} | {size:.1f} MB |")
    lines += ["", "| One variant | Time |", "| - | -: |", f"| Typical (median of {len(timings)}) | {statistics.median(timings):.1f} ms |", f"| Slowest | {max(timings):.1f} ms |"]
    if memory:
        lines += ["", f"Peak memory: {memory:,.0f} MB."]
    lines.append("")
    sys.stdout.flush()
    return "\n".join(lines)
