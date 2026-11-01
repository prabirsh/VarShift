"""Work out what kind of variant description a piece of text is.

Colleagues send these in many shapes, so the patterns are deliberately forgiving. Every pattern is
anchored and uses only simple repeats, and input length is capped before matching.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .builds import build_of_accession
from .errors import InputError
from .hgvs import DnaChange, ProteinChange, parse_dna_change, parse_protein_change

MAX_INPUT_LENGTH = 2000

_TX = r"(?:NM|NR|XM|XR|ENST)_?\d+(?:\.\d+)?"
_GENE = r"[A-Za-z][A-Za-z0-9_.\-]{0,30}"
_CONTIG = r"(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)"

_RSID_RE = re.compile(r"^rs(\d{1,12})$", re.IGNORECASE)
_GENOMIC_HGVS_RE = re.compile(rf"^(?P<contig>NC_\d{{6}}\.\d+|{_CONTIG})\s*:\s*(?P<level>[gm])\.\s*(?P<change>.+)$", re.IGNORECASE)
_KEY_RE = re.compile(
    rf"^(?P<contig>{_CONTIG})\s*[-:_\s]\s*(?P<pos>\d{{1,10}})\s*[-:_\s]\s*(?P<ref>[ACGTacgt]+)\s*(?:>|[-:_/\s])\s*(?P<alt>[ACGTacgt]+)$",
    re.IGNORECASE,
)
_TX_HGVS_RE = re.compile(
    rf"^(?P<tx>{_TX})(?:\s*\(\s*(?P<gene>{_GENE})\s*\))?\s*:\s*(?P<level>[cn])\.\s*(?P<change>[^()]+?)(?:\s*\(\s*p\.[^)]*\)+)?\s*$",
    re.IGNORECASE,
)
_GENE_CODING_RE = re.compile(
    rf"^(?P<gene>{_GENE})\s*[:\s]\s*\(?(?P<level>[cn])\.\s*(?P<change>[^()]+?)\)?\s*$"
)
_BARE_CODING_RE = re.compile(r"^(?P<level>[cn])\.\s*.+$")
_NP_HGVS_RE = re.compile(r"^(?P<np>(?:NP|XP|ENSP)_?\d+(?:\.\d+)?)\s*:\s*p\.\s*(?P<change>.+)$", re.IGNORECASE)
_GENE_PROTEIN_RE = re.compile(rf"^(?P<gene>{_GENE})[\s:,]+\(?\s*(?:p\.)?\s*(?P<change>\(?[A-Za-z*][^\s]*\)?)\s*\)?$")
_BARE_PROTEIN_RE = re.compile(r"^(?:p\.\s*)?\(?(?:[A-Z][a-z]{2}|[A-Z*])\d+[^\s]*$")


@dataclass(frozen=True)
class ParsedInput:
    raw: str
    kind: str  # rsid, genomic, position, coding or protein
    gene: str = ""
    transcript: str = ""
    protein_accession: str = ""
    contig: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    rsid: str = ""
    build_hint: str = ""
    dna: DnaChange | None = None
    protein: ProteinChange | None = None
    level: str = ""


def parse_input(text: str, gene: str = "") -> ParsedInput:
    """Classify one description. Raises InputError naming what is missing."""
    raw = text
    text = text.strip().strip('"').strip("'").strip()
    if not text:
        raise InputError("the variant description is empty")
    if len(text) > MAX_INPUT_LENGTH:
        raise InputError(f"the variant description is longer than {MAX_INPUT_LENGTH} characters")

    match = _RSID_RE.match(text)
    if match:
        return ParsedInput(raw, "rsid", rsid=match.group(1))

    match = _GENOMIC_HGVS_RE.match(text)
    if match:
        contig = match.group("contig")
        return ParsedInput(
            raw,
            "genomic",
            contig=contig,
            build_hint=(build_of_accession(contig) or "") if contig.upper().startswith("NC_") else "",
            dna=parse_dna_change(match.group("change")),
            level=match.group("level").lower(),
        )

    match = _KEY_RE.match(text)
    if match:
        return ParsedInput(
            raw,
            "position",
            contig=match.group("contig"),
            pos=int(match.group("pos")),
            ref=match.group("ref").upper(),
            alt=match.group("alt").upper(),
        )

    match = _TX_HGVS_RE.match(text)
    if match:
        return ParsedInput(
            raw,
            "coding",
            transcript=_clean_accession(match.group("tx")),
            gene=match.group("gene") or gene,
            dna=parse_dna_change(match.group("change")),
            level=match.group("level").lower(),
        )

    match = _NP_HGVS_RE.match(text)
    if match:
        return ParsedInput(
            raw,
            "protein",
            protein_accession=_clean_accession(match.group("np")),
            gene=gene,
            protein=parse_protein_change(match.group("change")),
        )

    match = _GENE_CODING_RE.match(text)
    if match and _plausible_gene(match.group("gene")):
        return ParsedInput(
            raw,
            "coding",
            gene=match.group("gene"),
            dna=parse_dna_change(match.group("change")),
            level=match.group("level").lower(),
        )

    match = _BARE_CODING_RE.match(text)
    if match:
        if not gene:
            raise InputError(f"{text!r} gives no gene or transcript; add one, as in BRAF {text}")
        return ParsedInput(
            raw, "coding", gene=gene, dna=parse_dna_change(text[2:]), level=match.group("level").lower()
        )

    match = _BARE_PROTEIN_RE.match(text)
    if match:
        if not gene:
            raise InputError(f"{text!r} gives no gene; add one, as in BRAF {text}")
        return ParsedInput(raw, "protein", gene=gene, protein=parse_protein_change(_strip_p(text)))

    match = _GENE_PROTEIN_RE.match(text)
    if match and _plausible_gene(match.group("gene")):
        return ParsedInput(
            raw, "protein", gene=match.group("gene"), protein=parse_protein_change(_strip_p(match.group("change")))
        )

    raise InputError(
        f"cannot read {text!r}. Expected a form such as: BRAF V600E, NM_004333.6:c.1799T>A, "
        "chr7:g.140453136A>T, chr7-140453136-A-T or rs113488022"
    )


def _strip_p(text: str) -> str:
    text = text.strip()
    if text.endswith(")") and "(" not in text:
        text = text[:-1]
    if text.lower().startswith("p."):
        return text[2:]
    return text


def _plausible_gene(token: str) -> bool:
    return token.lower() not in {"c", "n", "p", "g", "m", "chr"}


def _clean_accession(accession: str) -> str:
    accession = accession.upper()
    for prefix in ("NM", "NR", "XM", "XR", "NP", "XP"):
        if accession.startswith(prefix) and not accession.startswith(prefix + "_"):
            return f"{prefix}_{accession[len(prefix):]}"
    return accession
