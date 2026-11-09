"""How to write a variant: every form VarShift reads, each with an example that works.

The same text serves varshift examples, the guided start and its help, so they can't disagree. A test
annotates every example on the demo data, so none of them can stop working unnoticed.
"""

from __future__ import annotations

import textwrap

KINDS = ("protein", "cdna", "genomic", "rsid", "file", "several")
TITLES = {
    "protein": "Protein change",
    "cdna": "cDNA change",
    "genomic": "Genomic position",
    "rsid": "rs number",
    "file": "A file",
    "several": "Several at once",
}

# (kind, example, what it means)
EXAMPLES: list[tuple[str, str, str]] = [
    ("protein", "BRAF V600E", "Gene and change, with one-letter amino acids"),
    ("protein", "BRAF p.Val600Glu", "Three-letter amino acids; the p. can be left out"),
    ("protein", "NP_004324.2:p.Val600Glu", "A RefSeq protein and its change"),
    ("protein", "EGFR E746_A750del", "Deletions, duplications and insertions of residues"),
    ("protein", "BRCA1 E23fs", "A frameshift: ClinVar says which DNA changes give it"),
    ("cdna", "NM_004333.6:c.1799T>A", "A transcript and its c. change; the version can be left out"),
    ("cdna", "NM_000546.6(TP53):c.524G>A", "With the gene in brackets, as ClinVar writes it"),
    ("cdna", "BRAF c.1799T>A", "A gene and c. change, tried on its MANE Select transcript first"),
    ("cdna", "BRCA1 c.68_69del", "Deletions, duplications, insertions and delins too"),
    ("cdna", "NR_027676.2:n.200G>A", "A non-coding transcript, with n. numbering"),
    ("genomic", "chr7:g.140453136A>T", "A chromosome and g. change; add build=hg19 or build=hg38"),
    ("genomic", "NC_000007.13:g.140453136A>T", "A RefSeq chromosome, whose version also says the build"),
    ("genomic", "chr7-140453136-A-T", "Chromosome, position, REF and ALT, as gnomAD writes them"),
    ("genomic", "7:140453136:A>T", "The same with colons, with or without chr"),
    ("genomic", "chr7 140453136 A T", "The same with spaces, as copied from a table"),
    ("genomic", "NC_012920.1:m.3243A>G", "The mitochondrial genome (rCRS), with m. numbering; needs the hg38 data"),
    ("rsid", "rs113488022", "A dbSNP number, found through ClinVar (or dbSNP, once added)"),
    ("file", "demo/variants.txt", "A list: one description per line; lines starting with # are skipped"),
    ("file", "demo/demo.hg19.vcf", "A VCF, .vcf or .vcf.gz; its header says the build"),
    ("several", "BRAF V600E; KRAS G12D; TP53 R175H", "Several descriptions at once, separated by ;"),
]

# The five-line menu at the top of the guided start.
_MENU = (
    ("protein", ("BRAF V600E",)),
    ("cdna", ("NM_004333.6:c.1799T>A",)),
    ("genomic", ("chr7:g.140453136A>T", "chr7-140453136-A-T")),
    ("rsid", ("rs113488022",)),
    ("file", ("variants.txt", "sample.vcf.gz")),
)
MENU_KINDS = tuple(kind for kind, _ in _MENU)

_INTRO = (
    "VarShift reads a variant written in any of these ways. Give it to varshift annotate, "
    "with build= for a position, as in:"
)
_TRY = '  varshift annotate "BRAF V600E" build=hg19'


def menu_lines() -> list[str]:
    """One numbered line per kind of input, with its main example."""
    return [f"  {number}  {TITLES[kind]:<17}  {'   or   '.join(examples)}" for number, (kind, examples) in enumerate(_MENU, 1)]


def render_examples(term, kind: str | None = None) -> str:
    """Every example, grouped by kind, or those of one kind, wrapped at the terminal's width."""
    chosen = [entry for entry in EXAMPLES if kind is None or entry[0] == kind]
    if not chosen:
        raise ValueError(f"no examples of kind {kind!r}")
    width = term.width
    column = min(max(len(example) for _, example, _ in chosen), width // 2)
    room = max(width - 4 - column, 20)
    lines = _wrap(_INTRO, width) + [_TRY] if kind is None else []
    for group in KINDS:
        entries = [entry for entry in chosen if entry[0] == group]
        if not entries:
            continue
        lines += ["", term.style(TITLES[group], "strong")]
        for _, example, meaning in entries:
            parts = _wrap(meaning, room)
            if len(example) > column:
                lines.append(f"  {example}")
            else:
                lines.append(f"  {example:<{column}}  {parts.pop(0)}")
            lines += [" " * (4 + column) + part for part in parts]
    return "\n".join(lines).lstrip("\n") + "\n"


def _wrap(text: str, width: int) -> list[str]:
    # Never split "one-letter" or "build=" across lines.
    return textwrap.wrap(text, width, break_long_words=False, break_on_hyphens=False)
