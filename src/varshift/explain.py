"""Show exactly what happens to one variant: each step, each file read, each query, each timing."""

from __future__ import annotations

import time

from .annotate import Annotator
from .inputs import parse_input
from .store import DataStore

_WHAT = {
    "genes": "gene names and IDs (HGNC, ClinGen, gnomAD constraint, HPO)",
    "proteins": "protein domains (UniProt, Pfam)",
    "civic": "cancer evidence (CIViC)",
}
_FOR = {
    "transcripts": "gene structure (RefSeq, GENCODE)",
    "clinvar": "ClinVar records",
    "liftover": "hg19 to hg38 mapping",
}


def _describe(database: str) -> str:
    if database in _WHAT:
        return _WHAT[database]
    return _FOR.get(database.split(".")[0], database)


def _print_queries(entries: list[tuple[str, str]], started: float) -> None:
    for database, statement in entries:
        text = " ".join(statement.split())
        if len(text) > 150:
            text = text[:147] + "..."
        print(f"    {database}.sqlite  ({_describe(database)})")
        print(f"      {text}")
    print(f"    Time: {(time.perf_counter() - started) * 1000:.2f} ms")


def run_explain(store: DataStore, text: str, build: str | None, gene: str, mode: str | None = None) -> int:
    store.sql_log = []
    annotator = Annotator(store)
    total = time.perf_counter()

    print(f"Input: {text}\n")
    started = time.perf_counter()
    parsed = parse_input(text, gene=gene)
    print("Step 1. Read the text. No database is needed: patterns decide what kind of input this is.")
    details = {
        "kind": parsed.kind,
        "gene": parsed.gene,
        "transcript": parsed.transcript,
        "position": f"{parsed.contig}:{parsed.pos}" if parsed.pos else "",
        "change": (parsed.protein.clinvar_style if parsed.protein else ""),
    }
    print("    Understood as: " + ", ".join(f"{key} {value}" for key, value in details.items() if value))
    print(f"    Time: {(time.perf_counter() - started) * 1000:.2f} ms\n")

    chosen = annotator._pick_build(parsed, build)
    print(f"Step 2. Place it on the {chosen} genome.")
    print("    SQLite finds the transcript, then the reference file gives the bases to check against.")
    started = time.perf_counter()
    candidates = annotator.resolver.resolve(parsed, chosen)
    _print_queries(store.sql_log, started)
    store.sql_log.clear()
    print(f"    Result: {len(candidates)} candidate(s): {', '.join(c.variant.key for c in candidates[:6])}")
    print(f"    The REF base was read from genome/{chosen}.fa through its .fai index and matched.\n")

    best = candidates[0]
    print("Step 3. Annotate the best candidate. Each line below is one indexed lookup, not a scan of the file.")
    started = time.perf_counter()
    record = annotator.annotate_candidate(best)
    _print_queries(store.sql_log, started)
    store.sql_log.clear()

    clinvar = record["evidence"]["clinvar"]
    print("\nResult")
    print(f"    hg19      {record['variant']['hg19']['key']}")
    print(f"    hg38      {record['variant']['hg38']['key']}")
    print(f"    Names     {record['hgvs']['c']}  {record['hgvs']['p']}")
    print(f"    Effect    {record['effect']['type']}")
    print(f"    ClinVar   {clinvar.get('id')}: {clinvar.get('classification') or clinvar.get('status')}")
    print(f"    Priority  {record['priority']['label']}: {'; '.join(record['priority']['reasons'])}")
    print(f"    Mode      {mode or 'none'}")
    print(f"\nTotal: {(time.perf_counter() - total) * 1000:.1f} ms. Nothing left this computer.")
    return 0
