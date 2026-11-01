"""The interactive HTML report: one self-contained file that opens with a double click.

It needs no server and no internet: the page and its data travel together, so it can be emailed or
kept with a case. All text reaches the page through textContent, and only https links are made
clickable, so nothing in an input file can run as code.
"""

from __future__ import annotations

import json
from collections import Counter
from importlib import resources
from pathlib import Path

from .advanced import DATABASES, flatten_advanced
from .flags import FLAG_HELP
from .links import link_names

# A browser stays comfortable with this many rows. Beyond it, Low results are left out, and the page says so.
MAX_ROWS = 20_000
_LABELS = {"confirmed_id": "c", "built": "b"}


class HtmlCollector:
    """Keeps a compact form of each row, so a very large VCF can't use up the memory.

    Gene details and link addresses repeat across thousands of rows, so each is stored once.
    """

    def __init__(self):
        self.names = {name: index for index, name in enumerate(link_names())}
        self.rows: list[dict] = []
        self.genes: dict[str, dict] = {}
        self.texts: dict[str, int] = {}
        self.levels: Counter = Counter()
        self.total = 0
        self.left_out = 0

    def add(self, record: dict, prefix: dict | None) -> None:
        self.total += 1
        level = (record.get("priority") or {}).get("level", 0)
        self.levels[level] += 1
        if len(self.rows) >= MAX_ROWS:
            if level == 4:
                self.left_out += 1
                return
            for index in range(len(self.rows) - 1, -1, -1):
                if self.rows[index].get("p") == 4:
                    del self.rows[index]
                    self.left_out += 1
                    break
            else:
                self.left_out += 1
                return
        self.rows.append(_compact(record, prefix, self))

    def text_id(self, text: str) -> int:
        return self.texts.setdefault(text, len(self.texts))


def _compact(record: dict, prefix: dict | None, shared: HtmlCollector) -> dict:
    if record.get("error"):
        return {"i": record["input"]["raw"], "err": record["error"], "fl": record.get("flags", []), "p": 0, **_vcf_part(prefix)}
    variant = record["variant"]
    hgvs = record["hgvs"]
    effect = record["effect"]
    evidence = record["evidence"]
    clinvar = evidence.get("clinvar") or {}
    transcript = record.get("transcript") or {}
    priority = record.get("priority") or {}
    links = [
        [shared.names[link["name"]], _LABELS.get(link["label"], "b"), shared.text_id(link.get("url") or ""), shared.text_id(link.get("note") or "")]
        for link in record.get("links", [])
    ]
    symbol = record["gene"].get("symbol")
    if symbol and symbol not in shared.genes:
        details = {"gn": record["gene"].get("name"), "cgv": evidence.get("clingen_validity"), "cgd": evidence.get("clingen_dosage"), "con": evidence.get("constraint")}
        shared.genes[symbol] = {key: value for key, value in details.items() if value not in (None, "", [], {})}
    row = {
        "i": record["input"]["raw"],
        "r": record["candidate"].get("rank"),
        "n": record["candidate"].get("of"),
        "p": priority.get("level", 4),
        "pr": priority.get("reasons", []),
        "g": symbol,
        "tx": transcript.get("accession"),
        "mn": transcript.get("mane"),
        "c": hgvs.get("c"),
        "pn": hgvs.get("p"),
        "ps": hgvs.get("protein_short"),
        "gh": hgvs.get("g"),
        "e": effect.get("type"),
        "ea": effect.get("all"),
        "ex": effect.get("exon"),
        "in": effect.get("intron"),
        "k19": (variant.get("hg19") or {}).get("key"),
        "k38": (variant.get("hg38") or {}).get("key"),
        "cs": clinvar.get("status"),
        "cid": clinvar.get("id"),
        "cc": clinvar.get("classification"),
        "cr": clinvar.get("review_status"),
        "cd": clinvar.get("condition"),
        "co": clinvar.get("oncogenicity"),
        "csi": clinvar.get("somatic_impact"),
        "cv": clinvar.get("version"),
        "cf": clinvar.get("frequencies"),
        "cn": (clinvar.get("names") or [])[:4],
        "rs": evidence.get("dbsnp"),
        "dom": record["protein"].get("domains_at_residue"),
        "ft": record["protein"].get("features_at_residue"),
        "up": record["protein"].get("uniprot"),
        "civ": evidence.get("civic"),
        "adv": _advanced_part(evidence),
        "fl": record.get("flags", []),
        "l": links,
        "ot": record.get("other_transcripts", [])[:8],
    }
    row.update(_vcf_part(prefix))
    return {key: value for key, value in row.items() if value not in (None, "", [], {})}


def _advanced_part(evidence: dict) -> dict:
    row = flatten_advanced(evidence)
    return {key: value for key, value in row.items() if value not in (None, "")}


def _vcf_part(prefix: dict | None) -> dict:
    if not prefix:
        return {}
    part = {"vf": prefix.get("FILTER"), "vq": prefix.get("QUAL"), "vp": f"{prefix.get('CHROM')}:{prefix.get('POS')} {prefix.get('REF')}>{prefix.get('ALT')}"}
    return {key: value for key, value in part.items() if value}


def write_html(path: Path, collected: HtmlCollector, context) -> None:
    installed = [
        database.title
        for database in DATABASES.values()
        if any(a == database.name or a.startswith(f"{database.name}.") for a in context.installed)
    ]
    payload = {
        "input": context.input_name,
        "generated": context.generated,
        "notices": context.notices,
        "versions": context.versions,
        "advanced": installed,
        "links": list(collected.names),
        "genes": collected.genes,
        "texts": list(collected.texts),
        "flagHelp": FLAG_HELP,
        "total": collected.total,
        "leftOut": collected.left_out,
        "levels": {str(level): count for level, count in collected.levels.items()},
        "rows": collected.rows,
    }
    # Inside JSON a < only appears within strings, so this stops the data ever closing its script tag.
    data = json.dumps(payload, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c")
    template = resources.files("varsight").joinpath("data/report_template.html").read_text(encoding="utf-8")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(template.replace("/*VARSIGHT_DATA*/null", data))
