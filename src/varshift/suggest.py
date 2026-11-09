"""Step 2 of the input assistant: when a description fails, offer what it probably meant.

Step 1 (tidy.py) fixes slips that have only one reading. This step is for guesses, so it never applies
one: it turns the error into a question with numbered options, each built from the installed data, and
the person chooses. Where nobody can choose, as in a file, the options are only reported.

Each rule below looks at the error and the text, and either returns a question or leaves it to the
next rule. A suggestion must never stop a run, so a rule that fails gives no suggestion.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import NamedTuple

from .builds import HG19, HG38, parse_build
from .errors import InputError, ReferenceMismatch, ResolveError, VarShiftError
from .hgvs import parse_protein_change
from .inputs import parse_input
from .resolver import Resolver
from .sequence import AA1, AA3
from .tidy import reads, replace_word, tidy

MAX_OPTIONS = 10
MAX_SPELLINGS = 4
# Two edits are tried for longer symbols only, and not past this length, where there are too many to try
# in the time a suggestion may take: about 150 ms for 10 letters, 250 ms for 11.
TWO_EDITS_FROM, TWO_EDITS_UP_TO = 6, 10
_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-"
_DOT = " \u00b7 "

_RESIDUE_RE = re.compile(r"^(?P<gene>[A-Za-z][A-Za-z0-9_.\-]{0,30})[\s:]+(?:[pP]\.\s*)?\(?(?P<ref>[A-Z][a-z]{2}|[A-Z*])(?P<pos>\d+)\)?$")
_THREE_LETTER_RE = re.compile(r"[A-Z][a-z]{2}")
_SUBSTITUTION_END_RE = re.compile(r"([ACGTacgt])\s*>\s*([ACGTacgt])\s*$")


@dataclass(frozen=True)
class Option:
    text: str  # a whole description to annotate
    label: str  # short, for the numbered list
    detail: str  # the evidence shown beside it
    build: str | None = None  # set when the option is the same text on another build


@dataclass(frozen=True)
class Question:
    text: str
    options: tuple[Option, ...] = ()


def suggest(store, raw: str, error: Exception, build: str | None = None, fast_only: bool = False, gene: str = "") -> Question | None:
    """A question for the error that annotating raw gave, or None when there is nothing to offer.

    fast_only leaves out the one slow search (which genes have a protein change), for files.
    """
    text = tidy(raw, gene).text
    try:
        chosen = parse_build(build) if build else (HG38 if HG38 in store.builds else (store.builds or [None])[0])
    except InputError:
        return None
    if chosen is None:
        return None
    for rule in _RULES:
        try:
            question = rule(store, text, error, chosen, fast_only, gene)
        except (VarShiftError, ValueError):
            question = None
        if question is not None:
            return question
    return None


def apply_choice(records: list[dict], typed: str, option: Option) -> list[dict]:
    """Records annotated from a chosen option: the input stays as typed, and Read as says what was chosen."""
    source = records[0]["input"]  # shared by every record of one query
    source["read_as"] = source["read_as"] or option.text
    chosen = ["chosen from the suggestions"] + ([f"on {option.build}"] if option.build else [])
    source["read_as_notes"] = list(source.get("read_as_notes") or []) + chosen
    source["raw"] = typed.strip()
    return records


def option_texts(question: Question | None) -> list[str]:
    """What to type for each option, with build= when it is on the other build."""
    if question is None:
        return []
    return [option.text + (f" build={option.build}" if option.build else "") for option in question.options]


# 1. A protein change with no ALT: BRAF V600
def _missing_alt(store, text, error, build, fast_only, gene) -> Question | None:
    match = _RESIDUE_RE.match(text)
    if not isinstance(error, InputError) or not match:
        return None
    ref = match.group("ref")
    ref = AA1.get(ref) if len(ref) == 3 else ref
    symbol = store.gene_symbol(match.group("gene"))
    if ref not in AA3 or not symbol:
        return None
    pos = int(match.group("pos"))
    found: dict[str, dict] = {}
    prefix = f"p.{AA3[ref]}{pos}"
    for style, classification in store.clinvar_changes_at(build, symbol, prefix):
        if style[len(prefix) : len(prefix) + 1].isdigit():
            continue  # p.Val6001, not p.Val600
        entry = found.setdefault(one_letter(style), {"style": style, "classes": [], "civic": None})
        entry["classes"].append(_pretty(classification))
    for variant, items in store.civic_variants_at(symbol, f"{ref}{pos}"):
        change = _protein_change(variant)
        if change is None or change.pos != pos:
            continue
        entry = found.setdefault(variant.upper(), {"style": None, "classes": [], "civic": None})
        entry["civic"] = (entry["civic"] or 0) + items
    options = []
    for short, entry in found.items():
        description = f"{symbol} {short}"
        if not reads(description) and entry["style"]:
            description = f"{symbol} {entry['style']}"
        if not reads(description):
            continue
        classes = sorted(set(c for c in entry["classes"] if c))
        clinvar = "; ".join(classes[:2]) if classes else ("no classification" if entry["classes"] else "not in ClinVar")
        civic = f"CIViC: {entry['civic']} evidence items" if entry["civic"] is not None else "CIViC: none"
        rank = (-(entry["civic"] or 0), not any(_pathogenic(c) for c in classes), -len(entry["classes"]), short)
        options.append((rank, Option(description, short, f"ClinVar: {clinvar}{_DOT}{civic}")))
    if not options:
        return None
    options.sort(key=lambda pair: pair[0])
    return Question(f"{AA3[ref]}{pos} changes to what? Changes known at this residue:", tuple(o for _, o in options[:MAX_OPTIONS]))


# 2. A protein change with no gene: V600E
def _missing_gene(store, text, error, build, fast_only, gene) -> Question | None:
    if fast_only or not isinstance(error, InputError) or "gives no gene" not in str(error):
        return None
    change = _protein_change(re.sub(r"^[pP]\.\s*", "", text))
    if change is None:
        return None
    short = one_letter(change.clinvar_style)
    found: dict[str, list[int]] = {}
    for symbol, names in store.genes_with_change(build, change.clinvar_style):
        found.setdefault(symbol, [0, 0])[1] += names
    for symbol, items in store.civic_genes_with(short):
        found.setdefault(symbol, [0, 0])[0] += items
    options = []
    for symbol, (items, names) in sorted(found.items(), key=lambda pair: (-pair[1][0], -pair[1][1], pair[0])):
        description = f"{symbol} {text}"
        if reads(description):
            options.append(Option(description, symbol, f"ClinVar: {_records(names)}{_DOT}CIViC: {items} evidence items"))
    if not options:
        return None
    shown = f" ({MAX_OPTIONS} of {len(options)}, most evidence first)" if len(options) > MAX_OPTIONS else ""
    return Question(f"Which gene? Genes where ClinVar or CIViC records {short}{shown}:", tuple(options[:MAX_OPTIONS]))


# 3. A gene symbol that isn't one: BRFA V600E
def _unknown_gene(store, text, error, build, fast_only, gene) -> Question | None:
    if not isinstance(error, ResolveError) or "is not a known gene symbol or alias" not in str(error):
        return None
    typed = parse_input(text, gene=gene).gene
    counts: dict[str, int] = {}

    def count(symbol: str) -> int:
        if symbol not in counts:
            counts[symbol] = store.clinvar_count(build, symbol)
        return counts[symbol]

    word = typed.upper()
    # Most slips are one edit away: those few hundred names, aliases included, are looked up through the
    # indexes. Only when no approved symbol is one edit away are approved symbols two edits away tried.
    spellings = close_symbols(word, store.gene_symbols_for(edits(word)), count, two_edits=False)
    if not any(found.direct for found in spellings) and TWO_EDITS_FROM <= len(word) <= TWO_EDITS_UP_TO:
        found_symbols = {found.symbol for found in spellings}
        farther = [found for found in close_symbols(word, store.approved_symbols(), count) if found.symbol not in found_symbols]
        spellings = (spellings + farther)[:MAX_SPELLINGS]
    options = []
    for found in spellings:
        description = replace_word(text, typed, found.symbol)
        if description != text:  # a gene given with gene= isn't in the text, so it can't be rewritten
            through = "" if found.direct else f"{_DOT}through its alias {found.matched}"
            options.append(Option(description, found.symbol, f"ClinVar: {_records(count(found.symbol))}{through}"))
    if not options:
        return None
    return Question(f"{typed} isn't a gene symbol. Did you mean:", tuple(options))


# 4. A residue that isn't what the protein has: BRAF A600E
def _wrong_residue(store, text, error, build, fast_only, gene) -> Question | None:
    cause = _first_mismatch(error)
    if cause is None or not isinstance(cause.position, int) or cause.found not in AA3 or cause.expected not in AA3:
        return None
    parsed = parse_input(text, gene=gene)
    if parsed.kind != "protein":
        return None
    start = text.find(parsed.gene) + len(parsed.gene) if parsed.gene and parsed.gene in text else 0
    pattern = re.compile(rf"(?<![A-Za-z])({re.escape(AA3[cause.expected])}|{re.escape(cause.expected)})(?={cause.position}(?!\d))")
    match = pattern.search(text, start)
    if match is None:
        return None
    letter = AA3[cause.found] if len(match.group(1)) == 3 else cause.found
    description = text[: match.start()] + letter + text[match.end() :]
    question = f"Residue {cause.position} is {AA3[cause.found]}, not {AA3[cause.expected]}. Did you mean {description}?"
    return Question(question, (Option(description, description, f"residue {cause.position} of {cause.transcript}"),))


# 5. A c. or n. substitution whose REF isn't the transcript's base: BRAF c.1799A>T
def _wrong_base(store, text, error, build, fast_only, gene) -> Question | None:
    cause = _first_mismatch(error)
    if cause is None or not isinstance(cause.position, str) or len(cause.found) != 1:
        return None
    parsed = parse_input(text, gene=gene)
    if parsed.kind != "coding" or parsed.dna.kind != "sub":
        return None
    ref, alt, found = parsed.dna.ref.upper(), parsed.dna.alt.upper(), cause.found.upper()
    if alt == found:
        new_ref, new_alt, how = found, ref, "REF and ALT swapped"
    else:
        new_ref, new_alt, how = found, alt, "REF corrected"
    if new_ref == new_alt or not _SUBSTITUTION_END_RE.search(text):
        return None
    description = _SUBSTITUTION_END_RE.sub(f"{new_ref}>{new_alt}", text)
    name = f"{parsed.level}.{cause.position}"
    question = f"{name} is {found}. Did you mean {name}{new_ref}>{new_alt} ({how})?"
    return Question(question, (Option(description, f"{name}{new_ref}>{new_alt}", f"{how}, on {cause.transcript}"),))


# 6. A position that fits the other build: chr7-140753336-A-T build=hg19
def _other_build(store, text, error, build, fast_only, gene) -> Question | None:
    message = str(error)
    off_this_build = isinstance(error, ReferenceMismatch) or "is not covered by the installed" in message or "is outside the contig" in message
    if not off_this_build:
        return None
    parsed = parse_input(text, gene=gene)
    if parsed.kind not in ("position", "genomic") or parsed.build_hint:
        return None
    other = HG38 if build == HG19 else HG19
    why = f"The REF doesn't match {build}" if isinstance(error, ReferenceMismatch) else f"That position isn't in the installed {build} reference"
    if other not in store.builds:
        return Question(f"{why}. It may be on {other}, which isn't installed here. To add it: varshift setup build={other}")
    try:
        Resolver(store).resolve(parsed, other)
    except VarShiftError:
        return None
    option = Option(text, f"{text} build={other}", f"the REF matches {other}", build=other)
    return Question(f"{why}, but it does match {other}. Use build={other}?", (option,))


_RULES: tuple[Callable, ...] = (_missing_alt, _missing_gene, _unknown_gene, _wrong_residue, _wrong_base, _other_build)


# spelling
class Spelling(NamedTuple):
    symbol: str  # the approved symbol
    matched: str  # the symbol, or the alias, that was near the word
    edits: int  # 1 or 2

    @property
    def direct(self) -> bool:
        return self.matched == self.symbol.upper()


def edits(word: str) -> set[str]:
    """Every string one edit from the word: a letter removed, two neighbours swapped, replaced or added."""
    return set(_edit_list(word))


def _edit_list(word: str) -> list[str]:
    splits = [(word[:i], word[i:]) for i in range(len(word) + 1)]
    removed = [a + b[1:] for a, b in splits if b]
    swapped = [a + b[1] + b[0] + b[2:] for a, b in splits if len(b) > 1]
    replaced = [a + c + b[1:] for a, b in splits if b for c in _LETTERS]
    added = [a + c + b for a, b in splits for c in _LETTERS]
    return removed + swapped + replaced + added


def close_symbols(
    word: str, symbols: Mapping[str, str], count: Callable[[str], int], limit: int = MAX_SPELLINGS, two_edits: bool = True
) -> list[Spelling]:
    """Approved symbols near the word: one edit away, or two for longer words when two_edits allows.

    symbols maps upper-case symbols and aliases to the approved symbol. Fewer edits come first, then a
    symbol before one reached through an alias, then the gene with more ClinVar records. The standard
    library's close match search isn't used: for BRFA it offers RBFA, BRF2 and BRF1, and misses BRAF.
    """
    word = word.upper()
    best: dict[str, Spelling] = {}

    def keep(candidate: str, distance: int) -> None:
        found = Spelling(symbols[candidate], candidate, distance)
        held = best.get(found.symbol)
        if held is None or (found.edits, not found.direct) < (held.edits, not held.direct):
            best[found.symbol] = found

    near = _edit_list(word)
    for candidate in near:
        if candidate in symbols:
            keep(candidate, 1)
    if two_edits and TWO_EDITS_FROM <= len(word) <= TWO_EDITS_UP_TO:
        for first in set(near):
            for candidate in _edit_list(first):
                if candidate in symbols and candidate != word:
                    keep(candidate, 2)
    best.pop(symbols.get(word, ""), None)
    ranked = sorted(best.values(), key=lambda found: (found.edits, not found.direct, -(count(found.symbol) or 0), found.symbol))
    return ranked[:limit]


# helpers
def one_letter(style: str) -> str:
    """V600E from p.Val600Glu: three-letter amino acids as one letter, without the p."""
    body = re.sub(r"^p\.", "", style)
    return _THREE_LETTER_RE.sub(lambda match: AA1.get(match.group(0), match.group(0)), body)


def _protein_change(text: str):
    try:
        return parse_protein_change(text)
    except (VarShiftError, ValueError, KeyError):
        return None


def _first_mismatch(error: Exception) -> ReferenceMismatch | None:
    """The mismatch on the first transcript tried, which is the default one."""
    if isinstance(error, ResolveError) and error.causes and isinstance(error.causes[0], ReferenceMismatch):
        return error.causes[0]
    return None


def _records(count: int) -> str:
    # The counts are of ClinVar names, and on the hg19 data each record has one name per gene.
    return f"{count:,} record{'' if count == 1 else 's'}"


def _pretty(value: str | None) -> str:
    return (value or "").replace("_", " ").replace(" ,", ",")


def _pathogenic(classification: str) -> bool:
    return classification.lower().startswith(("pathogenic", "likely pathogenic"))
