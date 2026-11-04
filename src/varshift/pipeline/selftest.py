"""Fixed known variants, checked after every setup and before any update goes live."""

from __future__ import annotations

import json
from importlib import resources

from ..annotate import Annotator
from ..errors import VarShiftError
from ..store import DataStore


def known_variants() -> list[dict]:
    text = resources.files("varshift").joinpath("data/known_variants.json").read_text(encoding="utf-8")
    return json.loads(text)["variants"]


def run_selftest(store: DataStore, cases: list[dict] | None = None) -> list[dict]:
    annotator = Annotator(store)
    results = []
    for case in cases if cases is not None else known_variants():
        if case["build"] not in store.builds:
            continue
        results.append(_run_case(annotator, case))
    return results


def _run_case(annotator: Annotator, case: dict) -> dict:
    name = f"{case['input']} [{case['build']}]"
    try:
        records = annotator.annotate_text(case["input"], case["build"])
    except VarShiftError as exc:
        return {"name": name, "passed": False, "detail": f"could not annotate: {exc}"}
    expected = case["expect"]
    record = records[0]
    problems = []
    for build in ("hg19", "hg38"):
        want = expected.get(f"key_{build}")
        if want is None:
            continue
        got = record["variant"][build]["key"]
        # The other build's key can be absent when only one build is installed.
        if got != want and (got is not None or build == case["build"]):
            problems.append(f"{build} key {got}, expected {want}")
    if "key_any_candidate" in expected:
        found = {r["variant"][case["build"]]["key"] for r in records}
        if expected["key_any_candidate"] not in found:
            problems.append(f"{expected['key_any_candidate']} is not among the candidates {sorted(found)}")
    _compare(problems, "gene", record["gene"]["symbol"], expected.get("gene"))
    _compare(problems, "effect", record["effect"]["type"], expected.get("effect"))
    _compare_suffix(problems, "c.", record["hgvs"]["c"], expected.get("c"))
    _compare_suffix(problems, "p.", record["hgvs"]["p"], expected.get("p"))
    if "clinvar_id" in expected and record["evidence"]["clinvar"].get("id") != expected["clinvar_id"]:
        problems.append(f"ClinVar id {record['evidence']['clinvar'].get('id')}, expected {expected['clinvar_id']}")
    return {"name": name, "passed": not problems, "detail": "; ".join(problems) or "ok"}


def _compare(problems: list[str], label: str, got, want) -> None:
    if want is not None and got != want:
        problems.append(f"{label} {got}, expected {want}")


def _compare_suffix(problems: list[str], label: str, got: str | None, want: str | None) -> None:
    """Transcript versions change between releases, so only the change itself is compared."""
    if want is None:
        return
    body = got.split(":", 1)[-1] if got else None
    if body != want:
        problems.append(f"{label} {body}, expected {want}")
