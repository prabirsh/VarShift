"""Measure accuracy against ClinVar's own names, in both directions and across builds.

ClinVar and VarSight can each be wrong, so a mismatch is something to review rather than proof of a
bug. The report lists examples of every mismatch group for that reason.
"""

from __future__ import annotations

import re
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from ..consequence import SequenceCache, annotate
from ..errors import VarSightError
from ..hgvs import parse_dna_change, parse_protein_change, protein_style
from ..inputs import ParsedInput
from ..liftover import lift
from ..resolver import Resolver, _apply_on_transcript
from ..store import DataStore
from ..variant import Variant

_INS_TER_RE = re.compile(r"^p\.[A-Z][a-z]{2}\d+_([A-Z][a-z]{2}\d+)insTer$")
_DELINS_TER_RE = re.compile(r"^p\.([A-Z][a-z]{2}\d+)_[A-Z][a-z]{2}\d+delinsTer$")
_EXT_RE = re.compile(r"^p\.(Ter\d+[A-Z][a-z]{2})(?:extTer(?:\d+|\?))?$")
_OLD_STYLE_RE = re.compile(r"(del|dup)[ACGT]+$")  # older ClinVar names repeat the deleted bases
_EXAMPLES = 25
_MC_EQUIVALENT = {
    # ClinVar calls a frameshift that stops at once "nonsense"; VarSight keeps frameshift_variant.
    "nonsense": {"stop_gained", "frameshift_variant"},
    "missense_variant": {"missense_variant", "protein_altering_variant"},
    "synonymous_variant": {"synonymous_variant", "stop_retained_variant", "start_retained_variant"},
    "frameshift_variant": {"frameshift_variant"},
    "splice_donor_variant": {"splice_donor_variant"},
    "splice_acceptor_variant": {"splice_acceptor_variant"},
    "intron_variant": {"intron_variant", "splice_region_variant"},
    "5_prime_UTR_variant": {"5_prime_UTR_variant"},
    "3_prime_UTR_variant": {"3_prime_UTR_variant"},
    "inframe_deletion": {"inframe_deletion"},
    "inframe_insertion": {"inframe_insertion"},
    "inframe_indel": {"protein_altering_variant", "inframe_deletion", "inframe_insertion"},
    "initiator_codon_variant": {"start_lost", "start_retained_variant"},
    "stop_lost": {"stop_lost"},
    "non-coding_transcript_variant": {"non_coding_transcript_exon_variant", "non_coding_transcript_variant", "intron_variant"},
}


@dataclass
class Score:
    passed: int = 0
    total: int = 0
    examples: list[str] = field(default_factory=list)

    def add(self, ok: bool, example: str = "") -> None:
        self.total += 1
        if ok:
            self.passed += 1
        elif example and len(self.examples) < _EXAMPLES:
            self.examples.append(example)

    @property
    def rate(self) -> str:
        return f"{100 * self.passed / self.total:.2f}%" if self.total else "n/a"


@dataclass
class BuildResult:
    build: str
    rows: int = 0
    seconds: float = 0.0
    scores: dict[str, Score] = field(default_factory=lambda: defaultdict(Score))
    counts: Counter = field(default_factory=Counter)


def comparable(protein: str | None) -> str | None:
    """Put a p. name in one form, so that only real differences count against VarSight.

    ClinVar and HGVS allow several spellings of the same protein change. These are reduced to one form
    on both sides: p.Glu23fs for any frameshift, p.Met1? for any change of the start codon, p.Ser5Ter for
    a stop written as an insertion before residue 5 (p.Glu4_Ser5insTer), and p.Ter9Gln for an extension
    with or without its length.
    """
    if not protein:
        return None
    body = protein.split(":", 1)[-1]
    if body.startswith("p.Met1") and body != "p.Met1=":
        return "p.Met1?"
    match = _INS_TER_RE.match(body)
    if match:
        return f"p.{match.group(1)}Ter"
    match = _DELINS_TER_RE.match(body)
    if match:
        return f"p.{match.group(1)}Ter"
    match = _EXT_RE.match(body)
    if match:
        return f"p.{match.group(1)}"
    return protein_style(body)


def _class_of(ref: str, alt: str) -> str:
    return "snv" if len(ref) == 1 and len(alt) == 1 else "indel"


def run_validation(store: DataStore, builds: list[str], sample: int = 20000) -> list[BuildResult]:
    return [_validate_build(store, build, sample) for build in builds]


def _validate_build(store: DataStore, build: str, sample: int) -> BuildResult:
    result = BuildResult(build)
    db = store._db(f"clinvar.{build}")
    if db is None:
        return result
    total = db.execute("SELECT count(*) FROM clinvar_names").fetchone()[0]
    step = max(total // sample, 1) if sample else 1
    genome = store.genome(build)
    cache = SequenceCache(genome)
    resolver = Resolver(store)
    other = "hg38" if build == "hg19" else "hg19"
    started = time.time()
    rows = db.execute(
        "SELECT n.variation_id, n.transcript, n.gene, n.hgvs_c, n.hgvs_p, c.contig, c.pos, c.ref, c.alt, "
        "c.molecular_consequence FROM clinvar_names n JOIN clinvar c ON c.variation_id=n.variation_id "
        "WHERE n.rowid % ? = 0",
        (step,),
    )
    for row in rows:
        result.rows += 1
        variant = Variant(build, row["contig"], row["pos"], row["ref"], row["alt"])
        kind = _class_of(row["ref"], row["alt"])
        label = f"{row['transcript']}:{row['hgvs_c']} {variant.key} (ClinVar {row['variation_id']})"
        tx = _find_transcript(store, build, variant, row["transcript"], result)
        if tx is not None:
            _forward(genome, cache, variant, tx, row, kind, label, result)
        _reverse(resolver, build, variant, row, kind, label, result)
        if other in store.builds:
            _cross_build(store, variant, other, row, kind, label, result)
    result.seconds = time.time() - started
    return result


def _find_transcript(store, build, variant, named, result):
    candidates = store.transcript_versions(build, named)
    if not candidates:
        result.counts["transcript_not_in_annotation"] += 1
        return None
    exact = [t for t in candidates if t.accession == named and t.contig == variant.contig]
    if exact:
        result.counts["transcript_exact_version"] += 1
        return exact[0]
    same_contig = [t for t in candidates if t.contig == variant.contig]
    if not same_contig:
        result.counts["transcript_on_other_contig"] += 1
        return None
    result.counts["transcript_other_version"] += 1
    return same_contig[0]


def _forward(genome, cache, variant, tx, row, kind, label, result) -> None:
    try:
        effect = annotate(genome, variant, tx, cache)
    except VarSightError as exc:
        result.scores[f"forward_c_{kind}"].add(False, f"{label}: {exc}")
        return
    ours_c = effect.hgvs_c.split(":", 1)[1] if effect.hgvs_c else None
    same_version = tx.accession == row["transcript"]
    bucket = f"forward_c_{kind}" if same_version else f"forward_c_{kind}_other_version"
    theirs_c = _OLD_STYLE_RE.sub(r"\1", row["hgvs_c"])
    agrees = ours_c == theirs_c
    if not agrees and "[" in theirs_c:
        # ClinVar's repeat form (c.3901ACT[1]) names the same change; check it reaches the same variant.
        agrees = _same_variant(genome, tx, theirs_c, variant)
        if agrees:
            result.counts["repeat_notation_matched_by_position"] += 1
    result.scores[bucket].add(agrees, f"{label}: VarSight {ours_c}")
    theirs_p = row["hgvs_p"]
    if theirs_p and "[" in theirs_p:
        result.counts["protein_repeat_notation_not_compared"] += 1
    elif theirs_p:
        ours_p = comparable(effect.hgvs_p)
        bucket = f"forward_p_{kind}" if same_version else f"forward_p_{kind}_other_version"
        result.scores[bucket].add(ours_p == comparable(theirs_p), f"{label} {theirs_p}: VarSight {effect.hgvs_p}")
    consequence = row["molecular_consequence"]
    if consequence:
        wanted = set()
        known = True
        for term in consequence.split(","):
            if term in _MC_EQUIVALENT:
                wanted |= _MC_EQUIVALENT[term]
            else:
                known = False
        if wanted and known:
            result.scores[f"consequence_{kind}"].add(bool(wanted & set(effect.terms)), f"{label}: ClinVar {consequence}, VarSight {effect.terms}")


def _same_variant(genome, tx, name, variant) -> bool:
    try:
        return _apply_on_transcript(genome, tx, parse_dna_change(name.partition(".")[2])).key == variant.key
    except VarSightError:
        return False


def _reverse(resolver, build, variant, row, kind, label, result) -> None:
    store = resolver.store
    c_body = row["hgvs_c"]
    level, _, change = c_body.partition(".")
    versions = store.transcript_versions(build, row["transcript"])
    if not versions:
        result.counts["reverse_transcript_not_installed"] += 1
    elif level in ("c", "n") and not change.endswith("="):
        exact = any(t.accession == row["transcript"] for t in versions)
        bucket = f"reverse_c_{kind}" if exact else f"reverse_c_{kind}_other_version"
        try:
            parsed = ParsedInput(c_body, "coding", transcript=row["transcript"], dna=parse_dna_change(change), level=level)
            candidates = resolver.resolve(parsed, build)
            result.scores[bucket].add(
                any(c.variant.key == variant.key for c in candidates),
                f"{label}: got {[c.variant.key[:60] for c in candidates]}",
            )
        except VarSightError as exc:
            result.scores[bucket].add(False, f"{label}: {exc}")
    protein = row["hgvs_p"]
    if protein and "[" in protein:
        return
    if protein and protein.startswith("p.") and store.gene_symbol(row["gene"]) is None:
        result.counts["reverse_gene_not_installed"] += 1
    elif protein and protein.startswith("p."):
        try:
            change = parse_protein_change(protein[2:])
            parsed = ParsedInput(protein, "protein", gene=row["gene"], protein=change)
            candidates = resolver.resolve(parsed, build)
            result.scores[f"reverse_p_{kind}"].add(
                any(c.variant.key == variant.key for c in candidates),
                f"{row['gene']} {protein} {variant.key} (ClinVar {row['variation_id']}): got {[c.variant.key for c in candidates][:4]}",
            )
        except VarSightError as exc:
            result.scores[f"reverse_p_{kind}"].add(False, f"{row['gene']} {protein} {variant.key}: {exc}")


def _cross_build(store, variant, other, row, kind, label, result) -> None:
    db = store._db(f"clinvar.{other}")
    if db is None:
        return
    match = db.execute("SELECT contig, pos, ref, alt FROM clinvar WHERE variation_id=?", (row["variation_id"],)).fetchone()
    if match is None:
        result.counts["cross_build_no_clinvar_record"] += 1
        return
    expected = f"{match['contig']}-{match['pos']}-{match['ref']}-{match['alt']}"
    lifted = lift(store, variant, other)
    if lifted.variant is None and "other_build_region_not_installed" in lifted.flags:
        result.counts["cross_build_region_not_installed"] += 1
        return
    got = lifted.variant.key if lifted.variant else f"none ({','.join(lifted.flags)})"
    result.scores[f"cross_build_{kind}"].add(got == expected, f"{label}: expected {expected}, got {got}")


def render_report(results: list[BuildResult]) -> str:
    lines = ["# VarSight validation against ClinVar", ""]
    for result in results:
        lines += [f"## {result.build}", ""]
        rate = result.rows / result.seconds if result.seconds else 0
        lines += [f"{result.rows:,} ClinVar names tested in {result.seconds:.1f} s ({rate:,.0f} per second).", ""]
        lines += ["| Test | Passed | Total | Rate |", "| - | -: | -: | -: |"]
        for name in sorted(result.scores):
            score = result.scores[name]
            lines.append(f"| {name} | {score.passed:,} | {score.total:,} | {score.rate} |")
        lines.append("")
        if result.counts:
            lines += ["| Not scored | Count |", "| - | -: |"]
            for name, count in sorted(result.counts.items()):
                lines.append(f"| {name} | {count:,} |")
            lines.append("")
        for name in sorted(result.scores):
            score = result.scores[name]
            if score.examples:
                lines += [f"### Examples of mismatches: {name}", "", "```text"]
                lines += score.examples
                lines += ["```", ""]
    return "\n".join(lines)
