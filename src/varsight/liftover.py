"""Map a variant between hg19 and hg38 using a chain file stored in SQLite."""

from __future__ import annotations

from dataclasses import dataclass, field

from .builds import HG19, is_mitochondrial
from .errors import InputError, ReferenceMismatch
from .sequence import reverse_complement
from .store import DataStore
from .variant import Variant, normalize


@dataclass
class LiftResult:
    variant: Variant | None = None
    ref_check: str = "not_run"
    flags: list[str] = field(default_factory=list)


def _map_base(db, contig: str, position: int) -> tuple[str, int, int] | None:
    """Target contig, 0-based target base and strand for a 0-based source base."""
    row = db.execute("SELECT value FROM meta WHERE contig=?", (contig,)).fetchone()
    if row is None:
        return None
    longest = row["value"]
    rows = db.execute(
        "SELECT * FROM blocks WHERE contig=? AND start<=? AND start>=? AND end>? ORDER BY score DESC, start DESC",
        (contig, position, position - longest, position),
    ).fetchall()
    if not rows:
        return None
    block = rows[0]
    within = position - block["start"]
    if block["target_strand"] > 0:
        return block["target_contig"], block["target_start"] + within, 1
    return block["target_contig"], block["target_end"] - 1 - within, -1


def lift(store: DataStore, variant: Variant, target_build: str) -> LiftResult:
    """Map to the other build, then repeat the reference check there."""
    if variant.build == target_build:
        return LiftResult(variant, "pass")
    result = LiftResult()
    if variant.contig == "chrM" and variant.build == HG19:
        result.flags.append("chrM_not_rCRS")
        return result
    try:
        target_genome = store.genome(target_build)
    except Exception:
        result.flags.append("other_build_not_installed")
        return result

    if is_mitochondrial(variant.contig):
        contig = "chrM" if target_build != HG19 else "chrMT"
        try:
            result.variant = normalize(target_genome, contig, variant.pos, variant.ref, variant.alt)
            result.ref_check = "pass"
        except (ReferenceMismatch, InputError):
            result.ref_check = "fail"
            result.flags.append("ref_mismatch_on_other_build")
        return result

    db = store.liftover_db(variant.build, target_build)
    if db is None:
        result.flags.append("liftover_data_not_installed")
        return result
    first = _map_base(db, variant.contig, variant.pos - 1)
    last = _map_base(db, variant.contig, variant.pos + len(variant.ref) - 2)
    if first is None or last is None:
        result.flags.append("no_liftover")
        return result
    if first[0] != last[0] or first[2] != last[2]:
        result.flags.append("no_liftover")
        return result
    low, high = sorted((first[1], last[1]))
    if high - low != len(variant.ref) - 1:
        result.flags.append("liftover_span_changed")
        return result
    ref, alt = variant.ref, variant.alt
    if first[2] < 0:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
        result.flags.append("liftover_strand_flipped")
    if not target_genome.covers(first[0], low, high + 1):
        result.flags.append("other_build_region_not_installed")
        return result
    found = target_genome.fetch(first[0], low, high + 1)
    if found == ref:
        try:
            result.variant = normalize(target_genome, first[0], low + 1, ref, alt)
            result.ref_check = "pass"
        except (ReferenceMismatch, InputError):
            result.ref_check = "fail"
            result.flags.append("ref_mismatch_on_other_build")
        return result
    result.ref_check = "fail"
    if found == alt:
        result.flags.append("ref_alt_swapped_on_other_build")
    else:
        result.flags.append("ref_mismatch_on_other_build")
    return result
