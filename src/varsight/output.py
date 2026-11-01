"""Write results: on the screen, or as Excel, an interactive HTML report, TSV, JSON or JSON Lines.

A file being written carries an .incomplete suffix until the run finishes, so a stopped run can't be
mistaken for a finished one. Output files are readable only by their owner because they can hold
patient variants.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from . import xlsx
from .advanced import flatten_advanced, output_columns
from .errors import InputError
from .flags import FLAG_HELP
from .links import link_names
from .report import HtmlCollector, write_html

BASE_COLUMNS = [
    "priority",
    "priority_reasons",
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
]
_TAIL_COLUMNS = ["flags", "error"]
# The sheet people read first: what matters, then where to check it.
_PRIORITY_SHEET = [
    "priority",
    "priority_reasons",
    "gene",
    "hgvs_c",
    "hgvs_p",
    "effect",
    "exon",
    "clinvar_id",
    "clinvar_classification",
    "clinvar_review_status",
    "rsid",
    "hg19_key",
    "hg38_key",
    "flags",
    "link_ClinVar",
    "link_Franklin",
    "link_VarSome",
    "link_OncoKB",
    "link_COSMIC",
]
_LEVEL_STYLE = {"High": xlsx.HIGH, "Medium": xlsx.MEDIUM, "Review": xlsx.REVIEW, "Low": xlsx.LOW}
_MAX_PRIORITY_ROWS = 50_000
_FORMULA_START = ("=", "+", "@", "\t", "\r")
_FORMATS = {".json": "json", ".jsonl": "jsonl", ".tsv": "tsv", ".txt": "tsv", ".xlsx": "xlsx", ".html": "html"}


@dataclass
class ReportContext:
    """What every report needs to know besides the results."""

    input_name: str = ""
    prefix: list[str] = field(default_factory=list)
    installed: list[str] = field(default_factory=list)
    versions: dict[str, str] = field(default_factory=dict)
    notices: list[str] = field(default_factory=list)
    generated: str = field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M"))


def columns(context: ReportContext) -> list[str]:
    links = [f"link_{name}" for name in link_names()]
    return context.prefix + BASE_COLUMNS + output_columns(context.installed) + _TAIL_COLUMNS + links + ["data_versions"]


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
    priority = record.get("priority") or {}
    row = {
        "priority": priority.get("label"),
        "priority_reasons": "; ".join(priority.get("reasons", [])),
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
        "clingen_validity": ";".join(f"{v['classification']} ({v['disease']})" for v in evidence.get("clingen_validity", [])),
        "clingen_dosage": dosage.get("haploinsufficiency"),
        "gnomad_pli": constraint.get("pli"),
        "gnomad_loeuf": constraint.get("loeuf"),
        "flags": ";".join(record.get("flags", [])),
        "error": record.get("error"),
    }
    row.update(flatten_advanced(evidence))
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


def _guard_cell(value) -> str:
    """Stop a spreadsheet opening a TSV from running text as a formula, and keep rows intact."""
    if value is None:
        return ""
    text = str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
    if text.startswith(_FORMULA_START) or (text.startswith("-") and text[1:2].isalpha()):
        return "'" + text
    return text


class _FileSink:
    """A result file that stays .incomplete, and private, until it is closed."""

    def __init__(self, path: str | Path, replace: bool):
        self.final_path = Path(path)
        if self.final_path.is_symlink():
            raise InputError(f"{self.final_path} is a symbolic link; refusing to write through it")
        if self.final_path.exists() and not replace:
            raise InputError(f"{self.final_path} already exists; add the word replace to overwrite it")
        self.final_path.parent.mkdir(parents=True, exist_ok=True)
        self.path = self.final_path.with_name(self.final_path.name + ".incomplete")
        descriptor = os.open(self.path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        os.close(descriptor)
        self.rows = 0

    def _finish(self) -> None:
        os.replace(self.path, self.final_path)


class TsvSink(_FileSink):
    def __init__(self, path, context: ReportContext, replace: bool = False):
        super().__init__(path, replace)
        self.columns = columns(context)
        self._handle = open(self.path, "w", encoding="utf-8", newline="\n")
        self._handle.write("\t".join(self.columns) + "\n")

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        row = flatten(record)
        if prefix:
            row.update(prefix)
        self._handle.write("\t".join(_guard_cell(row.get(name)) for name in self.columns) + "\n")

    def close(self) -> None:
        self._handle.close()
        self._finish()


class JsonLinesSink(_FileSink):
    def __init__(self, path, context: ReportContext, replace: bool = False):
        super().__init__(path, replace)
        self._handle = open(self.path, "w", encoding="utf-8", newline="\n")

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        payload = record if not prefix else {"vcf": prefix, **record}
        self._handle.write(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n")

    def close(self) -> None:
        self._handle.close()
        self._finish()


class JsonSink(_FileSink):
    def __init__(self, path, context: ReportContext, replace: bool = False):
        super().__init__(path, replace)
        self._records: list[dict] = []

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        self._records.append(record if not prefix else {"vcf": prefix, **record})

    def close(self) -> None:
        with open(self.path, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(self._records, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        self._finish()


class ExcelSink(_FileSink):
    """A workbook with a summary, a sorted priority sheet, every result, and what the warnings mean."""

    def __init__(self, path, context: ReportContext, replace: bool = False):
        super().__init__(path, replace)
        self.context = context
        self.columns = columns(context)
        self._book = xlsx.XlsxWriter(self.path)
        self._book.add_sheet("All results", [18] * len(self.columns), table=True)
        self._book.write_row(self.columns, style=xlsx.HEADER)
        self._priority_at = self.columns.index("priority")
        self._top: list[dict] = []
        self._levels: Counter = Counter()
        self._effects: Counter = Counter()
        self._genes: Counter = Counter()
        self._flags: Counter = Counter()

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        row = flatten(record)
        if prefix:
            row.update(prefix)
        label = row.get("priority")
        styles = {self._priority_at: _LEVEL_STYLE[label]} if label in _LEVEL_STYLE else None
        self._book.write_row([row.get(name) for name in self.columns], cell_styles=styles)
        self._levels[label or "Not annotated"] += 1
        self._effects[row.get("effect") or "not annotated"] += 1
        if row.get("gene") and label in ("High", "Medium"):
            self._genes[row["gene"]] += 1
        for flag in record.get("flags", []):
            self._flags[flag] += 1
        if label in ("High", "Medium", "Review") and len(self._top) < _MAX_PRIORITY_ROWS:
            self._top.append({name: row.get(name) for name in _PRIORITY_SHEET})

    def _write_priority(self) -> None:
        order = {"High": 0, "Medium": 1, "Review": 2}
        widths = [10, 46, 12, 30, 28, 24, 6, 11, 32, 26, 14, 22, 22, 22, 18, 18, 18, 18, 18]
        self._book.add_sheet("Priority", widths, table=True)
        self._book.write_row(_PRIORITY_SHEET, style=xlsx.HEADER)
        for row in sorted(self._top, key=lambda item: order[item["priority"]]):
            self._book.write_row([row.get(name) for name in _PRIORITY_SHEET], cell_styles={0: _LEVEL_STYLE[row["priority"]]})

    def _write_summary(self) -> None:
        book = self._book
        book.add_sheet("Summary", [30, 16, 60])
        book.write_row(["VarSight report"], style=xlsx.TITLE)
        book.write_row(["Input", self.context.input_name])
        book.write_row(["Created", self.context.generated])
        book.write_row(["Results", self.rows])
        for line in self.context.notices:
            book.write_row([line], style=xlsx.BOLD)
        book.write_row([])
        book.write_row(["Review priority", "Results", "What it means"], style=xlsx.HEADER)
        meanings = {
            "High": "A database reports this as pathogenic, oncogenic, or with cancer evidence",
            "Medium": "Predicted to disrupt the protein, conflicting reports, or a high prediction score",
            "Review": "Changes the protein, or uncertain significance",
            "Low": "No protein change, reported benign, or common in the population",
        }
        for label, meaning in meanings.items():
            book.write_row([label, self._levels.get(label, 0), meaning], cell_styles={0: _LEVEL_STYLE[label]})
        if self._levels.get("Not annotated"):
            book.write_row(["Not annotated", self._levels["Not annotated"], "See the error column in All results"])
        book.write_row(["Priority only sorts the list. It is not a classification, and every row says why it was placed there."])
        book.write_row([])
        book.write_row(["Effect", "Results"], style=xlsx.HEADER)
        for effect, count in self._effects.most_common():
            book.write_row([effect, count])
        if self._genes:
            book.write_row([])
            book.write_row(["Genes with High or Medium results", "Results"], style=xlsx.HEADER)
            for gene, count in self._genes.most_common(25):
                book.write_row([gene, count])

    def _write_sources(self) -> None:
        book = self._book
        book.add_sheet("Data versions", [32, 24], table=True)
        book.write_row(["Database", "Version"], style=xlsx.HEADER)
        for name, version in sorted(self.context.versions.items()):
            book.write_row([name, version])

    def _write_flags(self) -> None:
        book = self._book
        book.add_sheet("Warnings explained", [34, 10, 60, 60], table=True)
        book.write_row(["Warning", "Results", "What it means", "What to do"], style=xlsx.HEADER)
        for flag, count in self._flags.most_common():
            meaning, action = FLAG_HELP.get(flag, ("", ""))
            book.write_row([flag, count, meaning, action])

    def close(self) -> None:
        self._write_priority()
        self._write_summary()
        self._write_sources()
        self._write_flags()
        self._book.close(order=["Summary", "Priority", "All results", "Warnings explained", "Data versions"])
        self._finish()


class HtmlSink(_FileSink):
    def __init__(self, path, context: ReportContext, replace: bool = False):
        super().__init__(path, replace)
        self.context = context
        self._collected = HtmlCollector()

    def write(self, record: dict, prefix: dict | None = None) -> None:
        self.rows += 1
        self._collected.add(record, prefix)

    def close(self) -> None:
        write_html(self.path, self._collected, self.context)
        self._finish()


_SINKS = {"tsv": TsvSink, "json": JsonSink, "jsonl": JsonLinesSink, "xlsx": ExcelSink, "html": HtmlSink}


def format_of(path: str) -> str | None:
    return _FORMATS.get(Path(path).suffix.lower())


def open_file_sink(path: str | Path, context: ReportContext, replace: bool = False, fmt: str | None = None):
    chosen = fmt or format_of(str(path))
    if chosen not in _SINKS:
        raise InputError(f"cannot tell the output format from {path}; use .xlsx, .html, .tsv, .json or .jsonl")
    return _SINKS[chosen](path, context, replace)


class ReportFolder:
    """The default for a file input: the report, the workbook and a TSV, side by side."""

    def __init__(self, folder: str | Path, context: ReportContext, replace: bool = False):
        self.folder = Path(folder)
        self.sinks = [
            HtmlSink(self.folder / "report.html", context, replace),
            ExcelSink(self.folder / "results.xlsx", context, replace),
            TsvSink(self.folder / "results.tsv", context, replace),
        ]

    def write(self, record: dict, prefix: dict | None = None) -> None:
        for sink in self.sinks:
            sink.write(record, prefix)

    def close(self) -> None:
        for sink in self.sinks:
            sink.close()

    @property
    def paths(self) -> list[Path]:
        return [sink.final_path for sink in self.sinks]


class ScreenSink:
    """Print to the terminal: a readable card, or JSON or TSV for piping into something else."""

    def __init__(self, fmt: str, context: ReportContext, stream=None):
        self.fmt = fmt
        self.context = context
        self.stream = stream or sys.stdout
        self.columns = columns(context)
        self._records: list[dict] = []
        if fmt == "tsv":
            self.stream.write("\t".join(self.columns) + "\n")

    def write(self, record: dict, prefix: dict | None = None) -> None:
        if self.fmt == "json":
            self._records.append(record)
        elif self.fmt == "jsonl":
            self.stream.write(json.dumps(record, separators=(",", ":")) + "\n")
        elif self.fmt == "tsv":
            row = flatten(record)
            self.stream.write("\t".join(_guard_cell(row.get(name)) for name in self.columns) + "\n")
        else:
            self.stream.write(render_card(record) + "\n")

    def close(self) -> None:
        if self.fmt == "json":
            json.dump(self._records, self.stream, indent=2)
            self.stream.write("\n")
        self.stream.flush()


def render_card(record: dict) -> str:
    """One result, laid out to be read in a terminal."""
    if record.get("error"):
        return f"{record['input']['raw']}\n  Could not annotate: {record['error']}"
    candidate = record["candidate"]
    variant = record["variant"]
    transcript = record.get("transcript") or {}
    effect = record["effect"]
    hgvs = record["hgvs"]
    evidence = record["evidence"]
    clinvar = evidence["clinvar"]
    priority = record.get("priority") or {}
    gene = record["gene"]["symbol"] or "no gene here"
    heading = f"{record['input']['raw']}"
    if candidate.get("of", 1) > 1:
        heading += f"   (candidate {candidate['rank']} of {candidate['of']})"
    where = f"exon {effect['exon']}" if effect.get("exon") else (f"intron {effect['intron']}" if effect.get("intron") else "")
    transcript_text = " ".join(filter(None, [transcript.get("accession"), f"({transcript['mane']})" if transcript.get("mane") else None]))
    lines = [
        heading,
        f"  Priority    {priority.get('label')}: {'; '.join(priority.get('reasons', []))}",
        f"  Position    hg19 {variant['hg19']['key'] or 'no mapping'}    hg38 {variant['hg38']['key'] or 'no mapping'}",
        f"  Gene        {gene}   {transcript_text}",
        f"  HGVS        {_short_name(hgvs['c'])}   {_short_name(hgvs['p'])}",
        f"  Effect      {' '.join(filter(None, [effect.get('type'), where])) or 'not determined'}",
    ]
    if clinvar.get("status") == "found":
        lines.append(f"  ClinVar     {clinvar['id']}: {clinvar['classification']} ({clinvar.get('review_status')})")
    else:
        lines.append(f"  ClinVar     {_clinvar_status(clinvar)}")
    domains = record["protein"].get("domains_at_residue") or []
    if domains:
        lines.append(f"  Domain      {domains[0]['name']} ({domains[0]['accession']})")
    if evidence.get("civic"):
        items = sum(item["evidence_items"] or 0 for item in evidence["civic"])
        lines.append(f"  CIViC       {evidence['civic'][0]['variant']}: {items} evidence items")
    for key, value in flatten_advanced(evidence).items():
        if value not in (None, ""):
            lines.append(f"  {key:<20} {value}")
    if record.get("flags"):
        lines.append(f"  Warnings    {', '.join(record['flags'])}")
    lines.append("  Check it")
    for link in record["links"]:
        if link.get("url"):
            lines.append(f"    {link['name']:<20} {link['url']}")
    return "\n".join(lines) + "\n"


def _short_name(name: str | None) -> str:
    return name or "none"
