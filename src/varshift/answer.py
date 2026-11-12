"""Sort the candidates for one query into one answer.

One protein change can come from several DNA changes, and an rs number can cover several alleles.
Rather than show every candidate with the same weight, each query gets one of three kinds of answer:

- single: there is one candidate.
- clear: one candidate stands out, because it is the only one in ClinVar, or, with none in ClinVar,
  the only one that changes the fewest bases. It is shown first, with why, and the others are named.
- choice: nothing decides between them, so they are listed for the user to choose.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass, field

from .advanced import flatten_advanced
from .flags import FLAG_HELP
from .output import LINK_LABELS, _clinvar_status
from .suggest import typo_note

ALSO_SHOWN = 3
LINKS_SHOWN = 5
MAX_CHOICES = 20
# The links a card shows first, in order of preference, for each mode.
KEY_LINKS = {
    None: ("ClinVar", "gnomAD", "UCSC Genome Browser", "dbSNP", "VarSome", "Franklin", "NCBI Gene"),
    "tumour": ("ClinVar", "OncoKB", "CIViC", "COSMIC", "cBioPortal", "gnomAD", "UCSC Genome Browser"),
    "inherited": ("ClinVar", "gnomAD", "ClinGen", "OMIM", "Franklin", "UCSC Genome Browser", "dbSNP"),
}
_DOT = " \u00b7 "
_LABEL = 12  # the width of a field name; values start after it and two spaces of indent
_INDENT = " " * (2 + _LABEL)
_NAME = 21  # the width of a link name on the Check it lines
# Also says it better than the warning would.
_NOT_NOTED = {"multiple_candidates"}


@dataclass
class Answer:
    raw: str
    read_as: str
    kind: str  # single, clear or choice
    records: list[dict]  # the best first, then rank order
    best: dict | None = None  # set for single and clear
    why: list[str] = field(default_factory=list)
    also: list[dict] = field(default_factory=list)  # the other records, for clear


def decide(records: list[dict]) -> Answer:
    """Decide the kind of answer for the results of one query, in the resolver's rank order."""
    if not records:
        raise ValueError("there are no results to decide between")
    source = records[0].get("input") or {}
    raw, read_as = source.get("raw", ""), source.get("read_as", "")
    if len(records) == 1:
        return Answer(raw, read_as, "single", list(records), best=records[0])
    # Each allele of an rs number is a different variant, so there is nothing to prefer.
    if source.get("type") == "rsid":
        return Answer(raw, read_as, "choice", list(records))

    in_clinvar = [r for r in records if r["candidate"].get("in_clinvar")]
    if len(in_clinvar) == 1:
        best, why = in_clinvar[0], ["in ClinVar; the others aren't"]
    elif in_clinvar:
        return Answer(raw, read_as, "choice", list(records))
    else:
        sizes = sorted(r["candidate"].get("bases_changed") or 0 for r in records)
        if sizes[0] == sizes[1]:
            return Answer(raw, read_as, "choice", list(records))
        best = next(r for r in records if (r["candidate"].get("bases_changed") or 0) == sizes[0])
        why = [f"fewest bases changed ({sizes[0]}; the others need {sizes[1]} or more)"]
    if ((best.get("transcript") or {}).get("mane") or "").startswith("MANE Select"):
        why.append("MANE Select transcript")
    others = [r for r in records if r is not best]
    return Answer(raw, read_as, "clear", [best] + others, best=best, why=why, also=others)


def also_text(answer: Answer) -> str:
    """The other candidates of a clear answer, as in "c.1799_1800delinsAA (2 bases, not in ClinVar)"."""
    parts = [_describe(r) for r in answer.also[:ALSO_SHOWN]]
    if len(answer.also) > ALSO_SHOWN:
        parts.append(f"and {len(answer.also) - ALSO_SHOWN} more")
    return "; ".join(parts)


def _describe(record: dict) -> str:
    c_name = (record.get("hgvs") or {}).get("c")
    if c_name:
        name = c_name.split(":", 1)[-1]
    else:
        build = (record.get("input") or {}).get("build")
        name = ((record.get("variant") or {}).get(build) or {}).get("key") or "another change"
    candidate = record["candidate"]
    bases = candidate.get("bases_changed") or 0
    status = "in ClinVar" if candidate.get("in_clinvar") else "not in ClinVar"
    return f"{name} ({bases} base{'' if bases == 1 else 's'}, {status})"


def apply(answer: Answer) -> list[dict]:
    """Write the answer into its records, best first, and return them ranked in that order."""
    for rank, record in enumerate(answer.records, 1):
        record["input"]["answer"] = answer.kind
        record["candidate"]["rank"] = rank
        record["candidate"]["best"] = record is answer.best
        record["candidate"]["why"] = list(answer.why) if record is answer.best else []
    if answer.kind == "clear":
        answer.best["also"] = also_text(answer)
    return answer.records


# The card: one answer, short enough to read at a glance, with the links to check it.
def render_answer(answer: Answer, term, mode: str | None = None, links_all: bool = False, also_hint: str = "show=all") -> str:
    """The best candidate's card. A choice has no best candidate, so it gets the choice table."""
    if answer.best is None:
        return render_choices(answer, term)
    record = answer.best
    source = record.get("input") or {}
    mode = _mode(mode if mode is not None else source.get("mode"))
    lines = _heading(answer, record, term)
    if answer.read_as:
        notes = source.get("read_as_notes") or []
        lines += _field(term, "Read as", answer.read_as + (f" ({'; '.join(notes)})" if notes else ""), "strong")
    priority = record.get("priority") or {}
    if priority.get("label"):
        lines += _field(term, priority["label"].upper(), _DOT.join(priority.get("reasons") or []), priority["label"].lower())
    lines += _field(term, "Gene", _gene_text(record))
    lines += _field(term, "Names", _names_text(record))
    lines += _field(term, "Effect", (record.get("effect") or {}).get("type") or "not determined")
    for build_line in _EVIDENCE[mode]:
        found = build_line(record, mode)
        if found:
            lines += _field(term, *found)
    advanced = [value for value in flatten_advanced(record.get("evidence") or {}).values() if value not in (None, "")]
    if advanced:
        lines += _field(term, "Advanced", f"{len(advanced)} values from the advanced databases. To see them: show=all")
    if answer.kind == "clear":
        lines += _field(term, "Why this", _DOT.join(answer.why))
        lines += _field(term, "Also", _also_line(answer, record, also_hint))
    lines += _link_lines(term, record.get("links") or [], mode, links_all)
    flags = record.get("flags") or []
    for flag in flags:
        # The typo note says which transcript was used, and why, so the general one would only repeat it.
        if flag in _NOT_NOTED or (flag == "other_transcript_used" and "default_transcript_differs" in flags):
            continue
        typo = typo_note([record]) if flag == "default_transcript_differs" else None
        lines += _field(term, "Note", typo or _note(flag))
    return "\n".join(lines) + "\n"


def _mode(mode: str | None) -> str | None:
    if mode in (None, "", "none"):
        return None
    if mode not in KEY_LINKS:
        raise ValueError(f"mode must be tumour, inherited or none, not {mode!r}")
    return mode


def _heading(answer: Answer, record: dict, term) -> list[str]:
    build = (record.get("input") or {}).get("build")
    variant = record.get("variant") or {}
    places = [f"{(variant.get(build) or {}).get('key')} ({build})"]
    places += [f"{position['key']} ({other})" for other, position in variant.items() if other != build and position.get("key")]
    return _wrap(f"{answer.raw}  ->  " + ", ".join(places), term.width, "    ")


def _field(term, label: str, text: str, role: str = "label") -> list[str]:
    """A field name, then its value wrapped at the terminal width under itself."""
    parts = _wrap(text, max(term.width - len(_INDENT), 20)) or [""]
    return ["  " + term.style(label.ljust(_LABEL), role) + parts[0]] + [_INDENT + part for part in parts[1:]]


def _wrap(text: str, width: int, indent: str = "") -> list[str]:
    # Names and links are never split, even when one is longer than the line.
    return textwrap.wrap(text, width, subsequent_indent=indent, break_long_words=False, break_on_hyphens=False)


def _short(name: str | None) -> str:
    """c.1799T>A from NM_004333.6:c.1799T>A."""
    return name.split(":", 1)[-1] if name else ""


def _gene_text(record: dict) -> str:
    transcript = record.get("transcript") or {}
    effect = record.get("effect") or {}
    where = f"exon {effect['exon']}" if effect.get("exon") else (f"intron {effect['intron']}" if effect.get("intron") else "")
    accession = transcript.get("accession") or ""
    if accession and transcript.get("mane"):
        accession += f" ({transcript['mane']})"
    parts = [(record.get("gene") or {}).get("symbol") or "no gene here", accession, where]
    return "   ".join(part for part in parts if part)


def _names_text(record: dict) -> str:
    hgvs = record.get("hgvs") or {}
    return "   ".join(name for name in (_short(hgvs.get("c")), _short(hgvs.get("p")), hgvs.get("g")) if name) or "none"


def _clinvar_line(record: dict, mode: str | None):
    clinvar = (record.get("evidence") or {}).get("clinvar") or {}
    if clinvar.get("status") != "found":
        return "ClinVar", _clinvar_status(clinvar) or "not looked up"
    called = clinvar.get("classification") or "no classification"
    if clinvar.get("review_status"):
        called += f" ({clinvar['review_status']})"
    parts = [str(clinvar["id"])]
    if mode == "tumour":
        if clinvar.get("oncogenicity"):
            parts.append(f"oncogenicity {clinvar['oncogenicity']}")
        if clinvar.get("somatic_impact"):
            parts.append(f"somatic impact {clinvar['somatic_impact']}")
    return "ClinVar", _DOT.join(parts + [called])


def _civic_line(record: dict, mode: str | None):
    civic = (record.get("evidence") or {}).get("civic") or []
    if not civic:
        return None
    items = sum(item.get("evidence_items") or 0 for item in civic)
    return "CIViC", f"{civic[0]['variant']}{_DOT}{items} evidence items"


def _domain_line(record: dict, mode: str | None):
    domains = (record.get("protein") or {}).get("domains_at_residue") or []
    return ("Domain", f"{domains[0]['name']} ({domains[0]['accession']})") if domains else None


def _clingen_line(record: dict, mode: str | None):
    validity = (record.get("evidence") or {}).get("clingen_validity") or []
    if not validity:
        return None
    parts = [f"{v['disease']} ({v['moi']}): {v['classification']}" for v in validity[:2]]
    if len(validity) > 2:
        parts.append(f"and {len(validity) - 2} more")
    return "ClinGen", "; ".join(parts)


def _constraint_line(record: dict, mode: str | None):
    constraint = (record.get("evidence") or {}).get("constraint")
    if not constraint:
        return None
    values = [f"{name} {constraint[key]:.2f}" for name, key in (("pLI", "pli"), ("LOEUF", "loeuf")) if constraint.get(key) is not None]
    return ("Constraint", _DOT.join(values) + f" ({constraint.get('source')})") if values else None


# The evidence lines on the card, in each mode's order.
_EVIDENCE = {
    None: (_clinvar_line, _civic_line, _domain_line),
    "tumour": (_civic_line, _clinvar_line, _domain_line),
    "inherited": (_clinvar_line, _clingen_line, _constraint_line, _domain_line),
}


def _also_line(answer: Answer, record: dict, hint: str) -> str:
    count = len(answer.also)
    gives = _short((record.get("hgvs") or {}).get("p")) or answer.raw
    which = "it" if count == 1 else "them"
    return f"{count} other DNA change{'s give' if count > 1 else ' gives'} {gives}: {also_text(answer)}. To see {which}: {hint}"


def _link_lines(term, links: list[dict], mode: str | None, links_all: bool) -> list[str]:
    if links_all:
        rows = []
        for link in links:
            rows.append(f"{link['name']:<{_NAME}}{link.get('url') or 'no link'}")
            about = "; ".join(part for part in (LINK_LABELS.get(link.get("label"), link.get("label")), link.get("note")) if part)
            rows += [" " * _NAME + part for part in _wrap(about, max(term.width - len(_INDENT) - _NAME, 20))]
    else:
        by_name = {link["name"]: link for link in links}
        chosen = [by_name[name] for name in KEY_LINKS[mode] if (by_name.get(name) or {}).get("url")][:LINKS_SHOWN]
        rows = [f"{link['name']:<{_NAME}}{link['url']}" for link in chosen]
        if len(links) > len(chosen):
            rows.append(f"{len(links) - len(chosen)} more: links=all")
    if not rows:
        return []
    return ["  " + term.style("Check it".ljust(_LABEL), "label") + rows[0]] + [_INDENT + row for row in rows[1:]]


def _note(flag: str) -> str:
    sentence, advice = FLAG_HELP.get(flag, (flag.replace("_", " "), ""))
    # Advice that is a command is worth the space; the rest is in the reports and the docs.
    return f"{sentence} {advice}" if "varshift " in advice else sentence


# The choice table: candidates that nothing decides between, one row each.
def choice_table(answer: Answer, term) -> tuple[str, list[str]]:
    """The heading and one row per candidate, without numbers, so the rows can also be offered as options."""
    total = len(answer.records)
    first = answer.records[0]
    if (first.get("input") or {}).get("type") == "rsid":
        where = "ClinVar" if (first.get("candidate") or {}).get("via") == "clinvar_rsid" else "dbSNP"
        heading = f"{answer.raw} covers {total} different variants in {where}:"
    else:
        heading = f"{answer.raw} fits {total} DNA changes:"
    shown = answer.records[:MAX_CHOICES]
    words = [_priority(record) for record in shown]
    # Room for the number in front ("  20  ") and the priority word at the end, which is never cut.
    room = max(term.width - 6 - 2 - max(len(word) for word, _ in words), 20)
    cells = [_choice_cells(record) for record in shown]
    # Columns: gene and protein, c. name, key, ClinVar. On a narrow terminal ClinVar goes first, then the
    # c. name, so the key that tells the variants apart is kept whole whenever it can be.
    for keep in ((0, 1, 2, 3), (0, 1, 2), (0, 2)):
        bodies = _aligned([[row[column] for column in keep] for row in cells])
        if max(len(body) for body in bodies) <= room:
            break
    column = min(max(len(body) for body in bodies), room)
    rows = [_fit(body, column).ljust(column) + "  " + term.style(word, role) for body, (word, role) in zip(bodies, words)]
    return heading, rows


def _aligned(cells: list[list[str]]) -> list[str]:
    widths = [max(len(row[column]) for row in cells) for column in range(len(cells[0]))]
    return ["  ".join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in cells]


def render_choices(answer: Answer, term) -> str:
    heading, rows = choice_table(answer, term)
    lines = [heading] + [f"{number:>3}  {row}" for number, row in enumerate(rows, 1)]
    hidden = len(answer.records) - len(rows)
    if hidden > 0:
        lines.append(f"     {hidden} more: be more specific")
    return "\n".join(lines) + "\n"


def _choice_cells(record: dict) -> list[str]:
    gene = (record.get("gene") or {}).get("symbol") or ""
    protein = (record.get("hgvs") or {}).get("protein_short") or ""
    return [f"{gene} {protein}".strip(), _short((record.get("hgvs") or {}).get("c")), _key(record), _classification(record)]


# The summary line: one query of several, on one line.
def summary_line(number: int, answer: Answer, term, input_width: int = 18) -> str:
    """Number, input, the best match (or how many to choose from), priority, and the main evidence.

    Lines for several queries line up when input_width is the longest of their inputs.
    """
    if answer.best is not None:
        key, (word, role), evidence = _key(answer.best), _priority(answer.best), _evidence_short(answer.best)
    else:
        known = sum(1 for record in answer.records if (record.get("candidate") or {}).get("in_clinvar"))
        key = f"{len(answer.records)} choices"
        word, role = min((_priority(record) for record in answer.records), key=lambda pair: _LEVELS.index(pair[0]))
        evidence = f"{known} of {len(answer.records)} in ClinVar"
    raw_room = max(term.width - (3 + 2 + 2 + max(len(key), 22) + 2 + len(word)), 8)
    raw = _fit(answer.raw, raw_room)
    head = f"{number:>3}  {raw:<{min(input_width, raw_room)}}  {key:<22}  "
    room = term.width - len(head) - len(word) - 2
    tail = "  " + _fit(evidence, room) if room >= 8 else ""
    return head + term.style(word, role) + tail


def error_line(number: int, text: str, message: str, term, input_width: int = 18) -> str:
    """A query of several that could not be annotated, lined up with the summary lines."""
    raw_room = max(term.width - 30, 8)
    raw = _fit(text, raw_room)
    head = f"{number:>3}  {raw:<{min(input_width, raw_room)}}  "
    word = "could not annotate"
    return head + term.style(word, "error") + ": " + _fit(message, max(term.width - len(head) - len(word) - 2, 8))


_LEVELS = ["HIGH", "MEDIUM", "REVIEW", "LOW", ""]


def _priority(record: dict) -> tuple[str, str]:
    label = (record.get("priority") or {}).get("label") or ""
    return label.upper(), (label.lower() or "label")


def _key(record: dict) -> str:
    build = (record.get("input") or {}).get("build")
    return ((record.get("variant") or {}).get(build) or {}).get("key") or ""


def _classification(record: dict) -> str:
    clinvar = (record.get("evidence") or {}).get("clinvar") or {}
    if clinvar.get("status") != "found":
        return "not in ClinVar"
    called = clinvar.get("classification") or "no classification"
    return "Conflicting" if called.lower().startswith("conflicting") else called


def _evidence_short(record: dict) -> str:
    called = _classification(record)
    parts = [called if called == "not in ClinVar" else f"ClinVar {called}"]
    civic = (record.get("evidence") or {}).get("civic") or []
    if civic:
        parts.append(f"CIViC {sum(item.get('evidence_items') or 0 for item in civic)} items")
    return _DOT.join(parts)


def _fit(text: str, room: int) -> str:
    """The text, cut to the room with ... when it is longer."""
    if len(text) <= room:
        return text
    return text[: max(room - 3, 1)].rstrip() + "..."
