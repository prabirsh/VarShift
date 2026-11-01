"""Parsers for HGVS change descriptions at DNA and protein level."""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import InputError
from .sequence import AA1, AA3
from .transcripts import TxPosition

_POS = r"[-*]?\d+(?:[+-]\d+)?"
_CHANGE_RE = re.compile(
    rf"^(?P<first>{_POS})(?:_(?P<last>{_POS}))?"
    r"(?:(?P<ref>[ACGTN]+)?>(?P<sub>[ACGTN]+)"
    r"|(?P<syn>(?P<same>[ACGTN]*)=)"
    r"|del(?P<del>[ACGTN]*)ins(?P<delins>[ACGTN]+)"
    r"|del(?:\d+|(?P<delonly>[ACGTN]*))"
    r"|dup(?:\d+|(?P<dup>[ACGTN]*))"
    r"|ins(?P<ins>[ACGTN]+)"
    r"|inv(?P<inv>[ACGTN]*)"
    r"|(?P<unit>[ACGTN]+)\[(?P<copies>\d+)\])$"
)
_POS_RE = re.compile(r"^(?P<region>[-*]?)(?P<base>\d+)(?P<offset>[+-]\d+)?$")
_FS_FORM_RE = re.compile(r"^p\.([A-Z][a-z]{2}\d+)(?:[A-Z][a-z]{2})?fs(?:Ter(?:\d+|\?))?$")

_AA_3 = "|".join(sorted((three for three in AA1 if three != "Xaa"), key=len, reverse=True))
_AA = rf"(?:{_AA_3}|[A-Z*])"
_AA_SEQ = rf"(?:{_AA})+"
_PROTEIN_RE = re.compile(
    rf"^(?P<ref>{_AA})(?P<pos>\d+)"
    rf"(?:_(?P<ref2>{_AA})(?P<pos2>\d+))?"
    r"(?:"
    r"(?P<syn>=)"
    r"|(?P<unknown>\?)"
    r"|(?P<del_>del)(?!ins)"
    r"|(?P<dup_>dup)"
    rf"|delins(?P<delins>{_AA_SEQ})"
    rf"|ins(?P<ins>{_AA_SEQ})"
    rf"|(?P<fs_alt>{_AA})?fs(?:(?:Ter|\*|X)(?P<fs_len>\d+|\?))?"
    rf"|(?P<ext_alt>{_AA})?ext(?:(?:Ter|\*|X)?(?P<ext_len>\d+|\?))?"
    rf"|(?P<alt>{_AA})"
    r")$"
)


@dataclass(frozen=True)
class DnaChange:
    """A change on a transcript (c., n.) or chromosome (g., m.)."""

    kind: str  # sub, del, dup, ins, delins, inv, syn or repeat
    first: TxPosition
    last: TxPosition
    ref: str = ""  # the repeated unit, for a repeat
    alt: str = ""  # the number of copies as text, for a repeat


@dataclass(frozen=True)
class ProteinChange:
    kind: str  # sub, syn, del, dup, ins, delins, fs, ext or start
    ref: str
    pos: int
    ref2: str = ""
    pos2: int = 0
    alt: str = ""
    fs_length: str = ""

    @property
    def clinvar_style(self) -> str:
        """The form ClinVar puts in variant names, used to match against its index."""
        first = f"{AA3[self.ref]}{self.pos}"
        span = f"{first}_{AA3[self.ref2]}{self.pos2}" if self.pos2 else first
        if self.kind == "sub":
            return f"p.{first}{AA3[self.alt]}"
        if self.kind == "syn":
            return f"p.{first}="
        if self.kind in ("del", "dup"):
            return f"p.{span}{self.kind}"
        if self.kind == "ins":
            return f"p.{span}ins{three_letter(self.alt)}"
        if self.kind == "delins":
            return f"p.{span}delins{three_letter(self.alt)}"
        if self.kind == "fs":
            return f"p.{first}fs"
        if self.kind == "start":
            return f"p.{first}?"
        return f"p.{first}{AA3[self.alt] if self.alt else ''}ext"


def three_letter(protein: str) -> str:
    return "".join(AA3[aa] for aa in protein)


def protein_style(name: str | None) -> str | None:
    """Reduce a p. name to the form used for matching: every frameshift becomes p.Glu23fs."""
    if not name:
        return None
    match = _FS_FORM_RE.match(name)
    return f"p.{match.group(1)}fs" if match else name


def parse_position(text: str) -> TxPosition:
    match = _POS_RE.match(text)
    if not match:
        raise InputError(f"cannot read the position {text!r}")
    return TxPosition(match.group("region"), int(match.group("base")), int(match.group("offset") or 0))


def parse_dna_change(text: str) -> DnaChange:
    """Parse the part after 'c.' or 'g.', for example 1799T>A or 68_69del."""
    compact = re.sub(r"\s+", "", text).strip("()")
    match = _CHANGE_RE.match(compact)
    if not match:
        raise InputError(f"cannot read the DNA change {text!r}")
    first = parse_position(match.group("first"))
    last = parse_position(match.group("last")) if match.group("last") else first
    if match.group("unit") is not None:
        return DnaChange("repeat", first, last, match.group("unit"), match.group("copies"))
    if match.group("sub") is not None:
        ref = match.group("ref") or ""
        alt = match.group("sub")
        if len(ref) == 1 and len(alt) == 1:
            return DnaChange("sub", first, last, ref, alt)
        return DnaChange("delins", first, last, ref, alt)
    if match.group("syn") is not None:
        return DnaChange("syn", first, last, match.group("same") or "")
    if match.group("delins") is not None:
        return DnaChange("delins", first, last, match.group("del") or "", match.group("delins"))
    if match.group("ins") is not None:
        if first == last:
            raise InputError(f"an insertion needs two flanking positions, as in 100_101insA: {text!r}")
        return DnaChange("ins", first, last, "", match.group("ins"))
    if match.group("inv") is not None:
        return DnaChange("inv", first, last, match.group("inv"))
    # Older names can give a count instead of the bases (del26, dup22), which leaves both groups empty.
    if "dup" in compact:
        return DnaChange("dup", first, last, match.group("dup") or "")
    return DnaChange("del", first, last, match.group("delonly") or "")


def _aa_one(token: str) -> str:
    if token == "X":
        return "*"
    if len(token) == 1:
        if token not in AA3:
            raise InputError(f"{token!r} is not an amino acid code")
        return token
    return AA1[token]


def _aa_seq(text: str) -> str:
    return "".join(_aa_one(token) for token in re.findall(_AA, text))


def parse_protein_change(text: str) -> ProteinChange:
    """Parse the part after 'p.', in three-letter or one-letter form."""
    compact = re.sub(r"\s+", "", text)
    if compact.startswith("(") and compact.endswith(")"):
        compact = compact[1:-1]
    match = _PROTEIN_RE.match(compact)
    if not match:
        raise InputError(f"cannot read the protein change {text!r}")
    ref = _aa_one(match.group("ref"))
    pos = int(match.group("pos"))
    if pos < 1:
        raise InputError(f"protein position must be 1 or more: {text!r}")
    ref2 = _aa_one(match.group("ref2")) if match.group("ref2") else ""
    pos2 = int(match.group("pos2") or 0)
    if pos2 and pos2 <= pos:
        raise InputError(f"protein range runs backwards: {text!r}")
    if match.group("syn"):
        return ProteinChange("syn", ref, pos)
    if match.group("unknown"):
        return ProteinChange("start", ref, pos)
    if match.group("del_"):
        return ProteinChange("del", ref, pos, ref2, pos2)
    if match.group("dup_"):
        return ProteinChange("dup", ref, pos, ref2, pos2)
    if match.group("delins"):
        return ProteinChange("delins", ref, pos, ref2, pos2, _aa_seq(match.group("delins")))
    if match.group("ins"):
        if not pos2:
            raise InputError(f"a protein insertion needs two flanking residues: {text!r}")
        return ProteinChange("ins", ref, pos, ref2, pos2, _aa_seq(match.group("ins")))
    if "fs" in compact:
        alt = _aa_one(match.group("fs_alt")) if match.group("fs_alt") else ""
        return ProteinChange("fs", ref, pos, alt=alt, fs_length=match.group("fs_len") or "")
    if "ext" in compact:
        alt = _aa_one(match.group("ext_alt")) if match.group("ext_alt") else ""
        return ProteinChange("ext", ref, pos, alt=alt, fs_length=match.group("ext_len") or "")
    return ProteinChange("sub", ref, pos, alt=_aa_one(match.group("alt")))
