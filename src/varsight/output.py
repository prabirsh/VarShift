"""Write results as JSON, JSON Lines, TSV or Excel.

A file being written carries an .incomplete suffix until the run finishes, so a stopped run can't be
mistaken for a finished one. Output files are readable only by their owner because they can hold
patient variants.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from .errors import InputError
from .links import link_names
from .xlsx import XlsxWriter

BASE_COLUMNS = [
    "input",
    "input_type",
    "build",
    "candidate_rank",
    "candidates_total",
    "hg19_key",
    "hg38_key",
    "gene",
    "transcript",
    "mane",
    "hgvs_g",
    "hgvs_c",
    "hgvs_p",
    "protein_short",
    "exon",
    "intron",
    "effect",
    "clinvar_id",
    "clinvar_classification",
    "clinvar_review_status",
    "clinvar_condition",
    "clinvar_oncogenicity",
    "clinvar_somatic_impact",
    "rsid",
    "uniprot",
    "pfam_at_residue",
    "pfam_gene",
    "civic_variant",
    "civic_evidence_items",
    "clingen_validity",
    "clingen_dosage",
    "gnomad_pli",
    "gnomad_loeuf",
    "gnomad_af",
    "alphamissense",
    "flags",
    "error",
]
_FORMULA_START = ("=", "+", "@", "\t", "\r")
_FORMATS = {".json": "json", ".jsonl": "jsonl", ".tsv": "tsv", ".txt": "tsv", ".xlsx": "xlsx"}


def format_for(path: str | None, requested: str | None) -> str:
    if requested:
        return requested
    if path is None:
        return "json"
    fmt = _FORMATS.get(Path(path).suffix.lower())
    if fmt is None:
        raise InputError(f"cannot tell the output format from {path}; use .json, .jsonl, .tsv or .xlsx, or add --format")
    return fmt


def flatten(record: dict) -> dict:
    """One flat row from one result, in the shape used by TSV and Excel."""
    variant = record.get("variant", {})
    hgvs = record.get("hgvs", {})
    effect = record.get("effect", {})
    protein = record.get("protein", {})
    evidence = record.get("evidence", {})
    clinvar = evidence.get("clinvar") or {}
    transcript = record.get("transcript") or {}
    gene = record.get("gene") or {}
    civic = evidence.get("civic") or []
    constraint = evidence.get("constraint") or {}
    dosage = evidence.get("clingen_dosage") or {}
    gnomad = evidence.get("gnomad") or {}
    am = evidence.get("alphamissense") or {}
    row = {
        "input": record.get("input", {}).get("raw"),
        "input_type": record.get("input", {}).get("type"),
        "build": record.get("input", {}).get("build"),
        "candidate_rank": record.get("candidate", {}).get("rank"),
        "candidates_total": record.get("candidate", {}).get("of"),
        "hg19_key": (variant.get("hg19") or {}).get("key"),
        "hg38_key": (variant.get("hg38") or {}).get("key"),
        "gene": gene.get("symbol"),
        "transcript": transcript.get("accession"),
        "mane": transcript.get("mane"),
        "hgvs_g": hgvs.get("g"),
        "hgvs_c": hgvs.get("c"),
        "hgvs_p": hgvs.get("p"),
        "protein_short": hgvs.get("protein_short"),
        "exon": effect.get("exon"),
        "intron": effect.get("intron"),
        "effect": effect.get("type"),
        "clinvar_id": clinvar.get("id"),
        "clinvar_classification": clinvar.get("classification") or _clinvar_status(clinvar),
        "clinvar_review_status": clinvar.get("review_status"),
        "clinvar_condition": clinvar.get("condition"),
        "clinvar_oncogenicity": clinvar.get("oncogenicity"),
        "clinvar_somatic_impact": clinvar.get("somatic_impact"),
        "rsid": evidence.get("dbsnp"),
        "uniprot": protein.get("uniprot"),
        "pfam_at_residue": ";".join(d["accession"] for d in protein.get("domains_at_residue", [])),
        "pfam_gene": ";".join(p["accession"] for p in protein.get("pfam", [])),
        "civic_variant": ";".join(c["variant"] for c in civic),
        "civic_evidence_items": sum(c["evidence_items"] or 0 for c in civic) or None,
        "clingen_validity": ";".join(
            f"{v['classification']} ({v['disease']})" for v in evidence.get("clingen_validity", [])
        ),
        "clingen_dosage": dosage.get("haploinsufficiency"),
        "gnomad_pli": constraint.get("pli"),
        "gnomad_loeuf": constraint.get("loeuf"),
        "gnomad_af": gnomad.get("af"),
        "alphamissense": am.get("am_class"),
        "flags": ";".join(record.get("flags", [])),
        "error": record.get("error"),
    }
    for link in record.get("links", []):
        row[f"link_{link['name']}"] = link.get("url")
    row["data_versions"] = ";".join(f"{k}={v}" for k, v in sorted(record.get("versions", {}).items()))
    return row


def _clinvar_status(clinvar: dict) -> str | None:
    if clinvar.get("status") == "not_found":
        return f"not in ClinVar (version {clinvar.get('version')})"
    if clinvar.get("status") == "not_installed":
        return "ClinVar not installed"
    return None


def columns(prefix: list[str] | None = None) -> list[str]:
    return (prefix or []) + BASE_COLUMNS + [f"link_{name}" for name in link_names()] + ["data_versions"]


def _guard_cell(value) -> str:
    """Stop a spreadsheet opening a TSV from running text as a formula, and keep rows intact."""
    if value is None:
        return ""
    text = str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    if text.startswith(_FORMULA_START) or (text.startswith("-") and text[1:2].isalpha()):
        return "'" + text
    return text


class ResultWriter:
    """Write rows to a file, or to standard output when path is None."""

    def __init__(self, path: str | None, fmt: str, prefix_columns: list[str] | None = None, force: bool = False):
        self.fmt = fmt
        self.prefix = prefix_columns or []
        self.columns = columns(self.prefix)
        self.final_path = Path(path) if path else None
        self.rows = 0
        self._records: list[dict] = []
        self._xlsx: XlsxWriter | None = None
        self._text = None
        if fmt == "xlsx" and path is None:
            raise InputError("Excel output needs a file; add --out results.xlsx")
        self._path = None
        if self.final_path is not None:
            if self.final_path.exists() and not force:
                raise InputError(f"{self.final_path} already exists; add --force to replace it")
            if self.final_path.is_symlink():
                raise InputError(f"{self.final_path} is a symbolic link; refusing to write through it")
            self.final_path.parent.mkdir(parents=True, exist_ok=True)
            self._path = self.final_path.with_name(self.final_path.name + ".incomplete")
            _touch_private(self._path)
        if fmt == "xlsx":
            self._xlsx = XlsxWriter(self._path)
            self._xlsx.add_sheet("VarSight", [24] * len(self.columns))
            self._xlsx.write_row(self.columns, header=True)
        elif self._path is not None:
            self._text = open(self._path, "w", encoding="utf-8", newline="\n")
        else:
            self._text = sys.stdout
        if fmt == "tsv":
            self._text.write("\t".join(self.columns) + "\n")

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        if self.fmt == "json":
            self._records.append(record if not prefix else {"vcf": prefix, **record})
            return
        if self.fmt == "jsonl":
            payload = record if not prefix else {"vcf": prefix, **record}
            self._text.write(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n")
            return
        row = flatten(record)
        if prefix:
            row.update(prefix)
        if self.fmt == "tsv":
            self._text.write("\t".join(_guard_cell(row.get(name)) for name in self.columns) + "\n")
        else:
            self._xlsx.write_row([row.get(name) for name in self.columns])

    def close(self) -> None:
        if self.fmt == "json":
            json.dump(self._records, self._text, ensure_ascii=False, indent=2)
            self._text.write("\n")
        if self._xlsx is not None:
            self._xlsx.close()
        elif self._path is not None:
            self._text.close()
        else:
            self._text.flush()
        if self._path is not None:
            os.replace(self._path, self.final_path)


def _touch_private(path: Path) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    descriptor = os.open(path, flags, 0o600)
    os.close(descriptor)
