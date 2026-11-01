"""A sorting aid: which results deserve a look first.

This is not a classification. It never says a variant is pathogenic or benign. It only puts the
results most likely to matter at the top of a long list, and always states why, so the reader can
disagree. Score thresholds are the ClinGen-calibrated ones where those exist: REVEL 0.644 (Pejaver et
al., 2022) and SpliceAI 0.2 (Walker et al., 2023).
"""

from __future__ import annotations

import re

HIGH, MEDIUM, REVIEW, LOW = 1, 2, 3, 4
LABELS = {HIGH: "High", MEDIUM: "Medium", REVIEW: "Review", LOW: "Low"}

_TRUNCATING = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost", "stop_lost"}
_PROTEIN_CHANGING = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "splice_region_variant"}
_REVEL_SUPPORTING = 0.644
_SPLICEAI_SUPPORTING = 0.2
_COMMON_FREQUENCY = 0.05


def _population_frequency(advanced: dict, clinvar: dict) -> tuple[float | None, str]:
    """The highest population frequency available, with where it came from."""
    best, source = None, ""
    candidates = [
        ("gnomAD popmax", (advanced.get("gnomad") or {}).get("af_popmax")),
        ("gnomAD", (advanced.get("gnomad") or {}).get("af")),
        ("TOPMed", (advanced.get("dbsnp") or {}).get("af_topmed")),
        ("gnomAD genomes", (advanced.get("dbsnp") or {}).get("af_gnomad_genomes")),
        ("gnomAD exomes", (advanced.get("dbsnp") or {}).get("af_gnomad_exomes")),
        ("1000 Genomes", (advanced.get("thousand_genomes") or {}).get("af")),
    ]
    candidates += [(f"{name.upper()} via ClinVar", value) for name, value in (clinvar.get("frequencies") or {}).items()]
    for name, value in candidates:
        if value is not None and (best is None or value > best):
            best, source = value, name
    return best, source


def assign(record: dict) -> dict:
    """Work out the priority for one result. Returns its level, label and reasons."""
    evidence = record.get("evidence", {})
    clinvar = evidence.get("clinvar") or {}
    advanced = evidence.get("advanced") or {}
    terms = set((record.get("effect") or {}).get("all") or [])
    classification = (clinvar.get("classification") or "").lower()
    oncogenicity = (clinvar.get("oncogenicity") or "").lower()
    somatic = (clinvar.get("somatic_impact") or "").lower()
    high: list[str] = []
    medium: list[str] = []
    review: list[str] = []
    low: list[str] = []

    if "pathogenic" in classification and "conflicting" not in classification and "benign" not in classification:
        high.append(f"ClinVar: {clinvar['classification']}")
    elif "conflicting" in classification:
        medium.append("ClinVar: conflicting classifications")
    elif "uncertain" in classification:
        review.append("ClinVar: uncertain significance")
    elif "benign" in classification:
        low.append(f"ClinVar: {clinvar['classification']}")
    elif classification and classification != "not provided":
        review.append(f"ClinVar: {clinvar['classification']}")
    if "oncogenic" in oncogenicity:
        high.append(f"ClinVar oncogenicity: {clinvar['oncogenicity']}")
    if re.search(r"tier ii? -", somatic):  # AMP/ASCO/CAP tiers I and II: strong or potential significance
        high.append(f"ClinVar somatic impact: {clinvar['somatic_impact']}")

    clingen = advanced.get("clingen_variants") or {}
    assertion = (clingen.get("assertion") or "").lower()
    if "pathogenic" in assertion:
        high.append(f"ClinGen expert panel: {clingen['assertion']}")
    elif "benign" in assertion:
        low.append(f"ClinGen expert panel: {clingen['assertion']}")

    civic = [item for item in evidence.get("civic") or [] if item.get("evidence_items")]
    if civic:
        high.append(f"CIViC: {sum(item['evidence_items'] for item in civic)} evidence items")
    hotspot = advanced.get("hotspots") or {}
    if hotspot.get("label"):
        high.append(f"Cancer hotspot: {hotspot['label']}")
    cosmic = advanced.get("cosmic") or {}
    if (cosmic.get("samples") or 0) >= 10:
        medium.append(f"COSMIC: {cosmic['samples']} tumour samples")

    truncating = sorted(terms & _TRUNCATING)
    if truncating:
        medium.append(f"Predicted to disrupt the protein: {truncating[0]}")
    elif terms & _PROTEIN_CHANGING:
        review.append(f"Changes the protein or nearby splicing: {sorted(terms & _PROTEIN_CHANGING)[0]}")

    revel = (advanced.get("revel") or {}).get("score")
    if revel is not None and revel >= _REVEL_SUPPORTING:
        medium.append(f"REVEL {revel:.3f} (at least {_REVEL_SUPPORTING})")
    alpha = (advanced.get("alphamissense") or {}).get("am_class")
    if alpha and "pathogenic" in alpha:
        medium.append(f"AlphaMissense: {alpha.replace('_', ' ')}")
    splice = (advanced.get("spliceai") or {}).get("max_delta")
    if splice is not None and splice >= _SPLICEAI_SUPPORTING:
        medium.append(f"SpliceAI {splice:.2f} (at least {_SPLICEAI_SUPPORTING})")

    frequency, frequency_source = _population_frequency(advanced, clinvar)
    common = frequency is not None and frequency >= _COMMON_FREQUENCY
    if common:
        low.append(f"Common in the population: {frequency:.3f} in {frequency_source}")

    if high:
        level, reasons = HIGH, high + medium
    elif low and "ClinVar" in low[0] or common:
        level, reasons = LOW, low + medium + review
    elif medium:
        level, reasons = MEDIUM, medium + review
    elif review:
        level, reasons = REVIEW, review
    else:
        level, reasons = LOW, low or ["No protein change and no database record of concern"]
    return {"level": level, "label": LABELS[level], "reasons": reasons}
