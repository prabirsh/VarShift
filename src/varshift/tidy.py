"""Step 1 of the input assistant: fix the slips that have only one sensible reading, and say so.

Text that already reads is never changed. Each rule below matches one anchored pattern, and the rules
stop as soon as the text reads. Anything that needs a guess, such as a missing amino acid, is left for
step 2 (suggest.py), which asks rather than assumes.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import VarShiftError
from .inputs import MAX_INPUT_LENGTH, parse_input
from .sequence import AA1


@dataclass(frozen=True)
class Reading:
    text: str
    notes: tuple[str, ...] = ()

    @property
    def changed(self) -> bool:
        return bool(self.notes)


# Characters that arrive when a variant is copied from a document, an email or a web page.
_PASTED = str.maketrans(
    {
        "\u201c": '"', "\u201d": '"', "\u2018": "'", "\u2019": "'",
        "\u2013": "-", "\u2014": "-", "\u2212": "-", "\uff1e": ">",
        "\u00a0": " ", "\u2002": " ", "\u2003": " ", "\u2009": " ", "\u202f": " ", "\u3000": " ",
        "\u200b": None, "\ufeff": None,
    }
)
_THOUSANDS_RE = re.compile(r"(?<![\w.,])(\d{1,3}(?:,\d{3})+)(?![\d,])")
_GENE_BRACKETS_RE = re.compile(r"^(?P<gene>[A-Za-z][A-Za-z0-9_.\-]{0,30})\((?P<change>(?:[pPcC]\.)?[A-Za-z0-9*_=>]+)\)$")
_DNA_CONTEXT_RE = re.compile(r"(?:\b[cngm]\.|:)\s*[-*]?\d", re.IGNORECASE)
_ARROW_RE = re.compile(r"(\d)([ACGTacgt])\s*(?:->|=>|>>|-|/|\u2192)\s*([ACGTacgt])\s*$")
_TRANSCRIPT_NO_LEVEL_RE = re.compile(
    r"^(?P<acc>(?:NM|XM|NR|XR)_?\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?)(?P<gene>\s*\([^)]*\))?\s*:\s*(?=[-*]?\d)", re.IGNORECASE
)
_CONTIG_NO_LEVEL_RE = re.compile(r"^(?P<contig>NC_\d{6}\.\d+|(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT))\s*:\s*(?=\d+\s*[ACGTacgt])", re.IGNORECASE)
_PROTEIN_PART_RE = re.compile(r"^(?P<head>.*?[\s:,(]|)(?P<prefix>[pP]\.\s*\(?)?(?P<change>[A-Za-z*][A-Za-z0-9_*=?]*?)(?P<tail>\)?)$")
_LOOKS_LIKE_PROTEIN_RE = re.compile(r"^[A-Za-z*]{1,3}\d")
_THREE_LETTER = {three.lower(): three for three in AA1 if three != "Xaa"}
_KEYWORDS_RE = re.compile(r"delins|del|dup|ins|ext|fs", re.IGNORECASE)


def tidy(text: str, gene: str = "") -> Reading:
    """Return the text to annotate, and a note for each change made to it."""
    original = text.strip()
    if len(original) > MAX_INPUT_LENGTH or _reads(original, gene):
        return Reading(original)
    current = original
    notes: list[str] = []
    for rule in _RULES:
        changed, note = rule(current)
        if changed != current:
            current = changed
            notes.append(note)
            if _reads(current, gene):
                break
    return Reading(current, tuple(notes))


def _reads(text: str, gene: str) -> bool:
    try:
        parse_input(text, gene=gene)
    except (VarShiftError, ValueError):
        return False
    return True


# The rules, in the order they are tried. Each returns the new text and what it did.
def _pasted_characters(text: str) -> tuple[str, str]:
    cleaned = text.translate(_PASTED).strip().strip("\"'").strip().rstrip(".,").rstrip()
    return cleaned, "replaced characters pasted from a document"


def _thousands(text: str) -> tuple[str, str]:
    return _THOUSANDS_RE.sub(lambda match: match.group(1).replace(",", ""), text), "removed thousands separators"


def _gene_brackets(text: str) -> tuple[str, str]:
    match = _GENE_BRACKETS_RE.match(text)
    if not match:
        return text, ""
    return f"{match.group('gene')} {match.group('change')}", "gene and change written as GENE CHANGE"


def _substitution_arrow(text: str) -> tuple[str, str]:
    if not _DNA_CONTEXT_RE.search(text):
        return text, ""
    return _ARROW_RE.sub(r"\1\2>\3", text), "substitution written with >"


def _missing_level(text: str) -> tuple[str, str]:
    match = _TRANSCRIPT_NO_LEVEL_RE.match(text)
    if match:
        accession = match.group("acc")
        level = "n" if accession.upper().startswith(("NR", "XR")) else "c"
        rest = text[match.end():]
        return f"{accession}{match.group('gene') or ''}:{level}.{rest}", f"added the missing {level}."
    match = _CONTIG_NO_LEVEL_RE.match(text)
    if match:
        contig = match.group("contig")
        mitochondrial = contig.upper().startswith("NC_012920") or contig.upper().removeprefix("CHR") in ("M", "MT")
        level = "m" if mitochondrial else "g"
        return f"{contig}:{level}.{text[match.end():]}", f"added the missing {level}."
    return text, ""


def _protein_case(text: str) -> tuple[str, str]:
    match = _PROTEIN_PART_RE.match(text)
    if not match or not _LOOKS_LIKE_PROTEIN_RE.match(match.group("change")):
        return text, ""
    change = re.sub(r"[A-Za-z]{3}", _standard_three_letter, match.group("change"))
    change = _KEYWORDS_RE.sub(lambda word: word.group(0).lower(), change)
    # One-letter codes: a lone letter before a position, and a lone letter at the very end.
    change = re.sub(r"(?<![A-Za-z])([a-z])(?=\d)", lambda letter: letter.group(1).upper(), change)
    change = re.sub(r"(?<=\d)([a-z])$", lambda letter: letter.group(1).upper(), change)
    prefix = (match.group("prefix") or "").replace("P.", "p.")
    return f"{match.group('head')}{prefix}{change}{match.group('tail')}", "amino acids written in standard case"


def _standard_three_letter(word: re.Match) -> str:
    return _THREE_LETTER.get(word.group(0).lower(), word.group(0))


_RULES = (_pasted_characters, _thousands, _gene_brackets, _substitution_arrow, _missing_level, _protein_case)
