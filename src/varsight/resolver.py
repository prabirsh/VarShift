"""Turn g., c. and p. descriptions into normalised genome variants, keeping every candidate.

One protein change can come from several DNA changes, so a p. input usually gives several
candidates. They are ranked by what can be defended: known to ClinVar first, then the fewest
changed bases.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .builds import canonical_contig, is_mitochondrial
from .errors import InputError, ReferenceMismatch, ResolveError
from .genome import Genome
from .hgvs import DnaChange, ProteinChange
from .inputs import ParsedInput
from .sequence import AA3, codons_for, reverse_complement, translate
from .store import DataStore
from .transcripts import Transcript, TxPosition
from .variant import Variant, normalize

MAX_CANDIDATES = 32


@dataclass
class Candidate:
    variant: Variant
    transcript: Transcript | None = None
    via: str = ""
    bases_changed: int = 0
    flags: list[str] = field(default_factory=list)
    in_clinvar: bool = False


class Resolver:
    def __init__(self, store: DataStore):
        self.store = store

    def resolve(self, parsed: ParsedInput, build: str) -> list[Candidate]:
        if parsed.kind == "position":
            contig = canonical_contig(parsed.contig, build)
            return [Candidate(normalize(self.store.genome(build), contig, parsed.pos, parsed.ref, parsed.alt), via="genomic")]
        if parsed.kind == "genomic":
            return [self._from_genomic(parsed, build)]
        if parsed.kind == "rsid":
            return self._from_rsid(parsed, build)
        if parsed.kind == "coding":
            return self._from_coding(parsed, build)
        return self._from_protein(parsed, build)

    # g.
    def _from_genomic(self, parsed: ParsedInput, build: str) -> Candidate:
        contig = canonical_contig(parsed.contig, build)
        if parsed.level == "m" and not is_mitochondrial(contig):
            raise InputError(f"m. numbering only applies to the mitochondrial genome, not {parsed.contig}")
        genome = self.store.genome(build)
        dna = parsed.dna
        for position in (dna.first, dna.last):
            if position.region or position.offset:
                raise InputError(f"a genomic position cannot use {position}; use a plain number")
        if dna.kind == "repeat":
            return Candidate(_apply_repeat(genome, contig, dna.first.base - 1, 1, dna.ref, int(dna.alt)), via="genomic")
        return Candidate(_apply_forward(genome, contig, dna.first.base - 1, dna.last.base, dna), via="genomic")

    # rs
    def _from_rsid(self, parsed: ParsedInput, build: str) -> list[Candidate]:
        rows = self.store.clinvar_by_rsid(build, parsed.rsid)
        known_to_clinvar = bool(rows)
        if not rows:
            rows = self.store.dbsnp_by_rsid(build, parsed.rsid)
        if not rows:
            raise ResolveError(
                f"rs{parsed.rsid} is not in the databases installed for {build}. An rs number needs ClinVar to "
                "have it, or dbSNP installed (varsight add dbsnp). You can also give the position or the HGVS name."
            )
        genome = self.store.genome(build)
        found = []
        for row in rows:
            variant = normalize(genome, row["contig"], row["pos"], row["ref"], row["alt"])
            source = "clinvar_rsid" if known_to_clinvar else "dbsnp_rsid"
            found.append(Candidate(variant, via=source, in_clinvar=known_to_clinvar, bases_changed=_size(variant)))
        return _rank(found)

    # c.
    def _from_coding(self, parsed: ParsedInput, build: str) -> list[Candidate]:
        genome = self.store.genome(build)
        transcripts, flags = self._transcripts_for(parsed, build)
        tried = []
        for number, tx in enumerate(transcripts):
            if parsed.level == "c" and not tx.coding:
                tried.append(f"{tx.accession}: non-coding, so c. numbering does not apply")
                continue
            if parsed.level == "n" and tx.coding:
                tried.append(f"{tx.accession}: protein-coding, so n. numbering does not apply")
                continue
            try:
                variant = _apply_on_transcript(genome, tx, parsed.dna)
            except (ReferenceMismatch, ResolveError, InputError) as exc:
                tried.append(f"{tx.accession}: {exc}")
                continue
            candidate = Candidate(variant, tx, via="coding", bases_changed=_size(variant), flags=list(flags))
            if number and not parsed.transcript:
                candidate.flags.append("other_transcript_used")
            return [candidate]
        raise ResolveError(f"no transcript fits {parsed.raw.strip()}", tried)

    def _transcripts_for(self, parsed: ParsedInput, build: str) -> tuple[list[Transcript], list[str]]:
        if parsed.transcript:
            return self._versions(parsed.transcript, build)
        if parsed.protein_accession:
            matches = self.store.transcripts_for_protein(build, parsed.protein_accession)
            if not matches:
                raise ResolveError(f"{parsed.protein_accession} is not among the {build} transcripts")
            exact = [t for t in matches if t.protein == parsed.protein_accession]
            if exact:
                return exact, []
            return matches, ["protein_version_substituted"]
        symbol = self._approved_symbol(parsed.gene)
        found = self.store.transcripts_for_gene(build, symbol)
        if not found:
            raise ResolveError(f"no {build} transcripts are installed for {symbol}")
        return found, []

    def _versions(self, accession: str, build: str) -> tuple[list[Transcript], list[str]]:
        matches = self.store.transcript_versions(build, accession)
        if not matches:
            raise ResolveError(f"{accession} is not among the installed {build} transcripts")
        if "." not in accession:
            newest = matches[0].accession
            return [t for t in matches if t.accession == newest], []
        exact = [t for t in matches if t.accession == accession]
        if exact:
            return exact, []
        # The genome-specific annotation may only hold another version of the same transcript.
        return matches, ["transcript_version_substituted"]

    def _approved_symbol(self, gene: str) -> str:
        if not gene:
            raise InputError("a gene or transcript is needed for c. and p. input")
        symbol = self.store.gene_symbol(gene)
        if symbol is None:
            raise ResolveError(f"{gene} is not a known gene symbol or alias")
        return symbol

    # p.
    def _from_protein(self, parsed: ParsedInput, build: str) -> list[Candidate]:
        genome = self.store.genome(build)
        try:
            transcripts, flags = self._transcripts_for(parsed, build)
        except ResolveError:
            # Mitochondrial genes have no transcript model here, but ClinVar can still name the variant.
            if not parsed.gene or self.store.gene_symbol(parsed.gene) is None:
                raise
            transcripts, flags = [], []
        coding = [t for t in transcripts if t.coding]
        change = parsed.protein
        tried: list[str] = []
        found: list[Candidate] = []
        gene = coding[0].gene if coding else self.store.gene_symbol(parsed.gene)
        for number, tx in enumerate(coding):
            try:
                found = _protein_candidates(genome, tx, change)
            except (ResolveError, ReferenceMismatch, InputError) as exc:
                tried.append(f"{tx.accession}: {exc}")
                continue
            for candidate in found:
                candidate.flags.extend(flags)
                if number and not parsed.protein_accession:
                    candidate.flags.append("other_transcript_used")
            break

        known = {c.variant.key: c for c in found}
        for row in self.store.clinvar_by_protein(build, gene, change.clinvar_style):
            try:
                variant = normalize(genome, row["contig"], row["pos"], row["ref"], row["alt"])
            except (ReferenceMismatch, InputError):
                continue
            if variant.key not in known:
                known[variant.key] = Candidate(variant, via="clinvar_protein", bases_changed=_size(variant))
            known[variant.key].in_clinvar = True
        for candidate in known.values():
            if not candidate.in_clinvar:
                v = candidate.variant
                candidate.in_clinvar = self.store.clinvar_row(build, v.contig, v.pos, v.ref, v.alt) is not None
        if not known:
            if change.kind in ("fs", "ext", "ins", "delins"):
                tried.insert(
                    0,
                    f"{change.clinvar_style} does not define one DNA change, and ClinVar has no record of it "
                    "for this gene. Give the c. name instead.",
                )
            raise ResolveError(f"no transcript fits {parsed.raw.strip()}", tried)
        return _rank(list(known.values()))


def _protein_candidates(genome: Genome, tx: Transcript, change: ProteinChange) -> list[Candidate]:
    seq = tx.spliced(genome)
    first, stop = tx.cds_tx_range
    cds = seq[first:stop]
    protein = translate(cds, to_stop=False)
    _check_residue(tx, protein, change.ref, change.pos)
    if change.pos2:
        _check_residue(tx, protein, change.ref2, change.pos2)

    if change.kind in ("sub", "syn"):
        target = change.alt if change.kind == "sub" else change.ref
        return _codon_candidates(genome, tx, cds, change.pos, target)
    if change.kind in ("del", "dup"):
        last = change.pos2 or change.pos
        first_c = TxPosition("", (change.pos - 1) * 3 + 1)
        last_c = TxPosition("", last * 3)
        variant = _apply_on_transcript(genome, tx, DnaChange(change.kind, first_c, last_c))
        return [Candidate(variant, tx, via="protein", bases_changed=(last - change.pos + 1) * 3)]
    return []


def _check_residue(tx: Transcript, protein: str, residue: str, position: int) -> None:
    if position > len(protein):
        raise ResolveError(f"residue {position} is past the end of the {len(protein)} residue protein")
    found = protein[position - 1]
    if found != residue:
        raise ReferenceMismatch(
            f"residue {position} is {AA3[found]}, not {AA3[residue]}", expected=residue, found=found
        )


def _codon_candidates(genome: Genome, tx: Transcript, cds: str, residue: int, target: str) -> list[Candidate]:
    start = (residue - 1) * 3
    codon = cds[start : start + 3]
    found = []
    for new_codon in codons_for(target):
        if new_codon == codon:
            continue
        differs = [i for i in range(3) if codon[i] != new_codon[i]]
        first, last = differs[0], differs[-1]
        first_pos = TxPosition("", start + first + 1)
        last_pos = TxPosition("", start + last + 1)
        if first == last:
            dna = DnaChange("sub", first_pos, last_pos, codon[first], new_codon[first])
        else:
            dna = DnaChange("delins", first_pos, last_pos, codon[first : last + 1], new_codon[first : last + 1])
        low, high = sorted((tx.genomic_of(first_pos), tx.genomic_of(last_pos)))
        if high - low != last - first:
            continue  # the codon is split by an intron; no single contiguous DNA change makes this
        variant = _apply_on_transcript(genome, tx, dna)
        found.append(Candidate(variant, tx, via="protein", bases_changed=len(differs)))
    return found


def _apply_repeat(genome: Genome, contig: str, first: int, strand: int, unit: str, copies: int) -> Variant:
    """Read a tandem repeat such as c.3901ACT[1]: count the copies present, then remove or add."""
    size = len(unit)
    forward_unit = unit if strand > 0 else reverse_complement(unit)
    present = 0
    while present < 1000:
        if strand > 0:
            piece = genome.fetch(contig, first + present * size, first + (present + 1) * size)
        else:
            piece = genome.fetch(contig, first + 1 - (present + 1) * size, first + 1 - present * size)
        if piece != forward_unit:
            break
        present += 1
    if present == 0:
        raise ReferenceMismatch(f"the repeat unit {unit} is not at {contig}:{first + 1}", expected=unit)
    if copies == present:
        raise InputError(f"the reference already has {present} copies of {unit} here, so nothing changes")
    low = first if strand > 0 else first + 1 - present * size
    if copies < present:
        removed = (present - copies) * size
        return _apply_forward(genome, contig, low, low + removed, DnaChange("del", TxPosition("", 0), TxPosition("", 0)))
    anchor = genome.fetch(contig, low - 1, low)
    return normalize(genome, contig, low, anchor, anchor + forward_unit * (copies - present))


def _apply_on_transcript(genome: Genome, tx: Transcript, dna: DnaChange) -> Variant:
    """Convert a transcript-level change to a normalised genome variant."""
    g_first = tx.genomic_of(dna.first)
    if dna.kind == "repeat":
        return _apply_repeat(genome, tx.contig, g_first, tx.strand, dna.ref, int(dna.alt))
    g_last = tx.genomic_of(dna.last)
    low, high = min(g_first, g_last), max(g_first, g_last)
    if tx.strand > 0 and g_first > g_last or tx.strand < 0 and g_first < g_last:
        raise InputError(f"the positions in {dna.first}_{dna.last} run backwards along {tx.accession}")
    if dna.kind in ("del", "dup") and _crosses_junction(tx, dna, g_first, g_last):
        moved = _move_off_junction(genome, tx, g_first, g_last)
        if moved is None:
            raise ResolveError(
                f"{dna.first}_{dna.last} joins two neighbouring exons of {tx.accession}, which can't be written as "
                "one change on the genome without also removing the intron between them"
            )
        low, high = moved
    forward = DnaChange(
        dna.kind,
        dna.first,
        dna.last,
        dna.ref if tx.strand > 0 else reverse_complement(dna.ref),
        dna.alt if tx.strand > 0 else reverse_complement(dna.alt),
    )
    return _apply_forward(genome, tx.contig, low, high + 1, forward)


def _crosses_junction(tx: Transcript, dna: DnaChange, g_first: int, g_last: int) -> bool:
    """True when two exonic positions are next-door in the transcript but an intron apart on the genome.

    ClinVar's 3' rule works on the spliced sequence, so a small deletion at the end of an exon can be
    named across the junction (c.8331_8332del). Read literally on the genome it would delete a whole intron.
    """
    if dna.first.offset or dna.last.offset:
        return False
    first_index, last_index = tx.exon_index_of(g_first), tx.exon_index_of(g_last)
    if first_index is None or last_index is None:
        return False
    return abs(last_index - first_index) != abs(g_last - g_first) and abs(last_index - first_index) < 100


def _move_off_junction(genome: Genome, tx: Transcript, g_first: int, g_last: int) -> tuple[int, int] | None:
    """Slide the change 5' on the spliced sequence until it fits inside one exon."""
    seq = tx.spliced(genome)
    first = tx.exon_index_of(g_first)
    last = tx.exon_index_of(g_last)
    while first > 0 and seq[first - 1] == seq[last]:
        first, last = first - 1, last - 1
        low, high = sorted((tx.genomic_of_index(first), tx.genomic_of_index(last)))
        if high - low == last - first:
            return low, high
    return None


def _apply_forward(genome: Genome, contig: str, start: int, end: int, dna: DnaChange) -> Variant:
    """Apply a forward-strand change to the 0-based half-open span [start, end)."""
    if not genome.covers(contig, max(start - 1, 0), end):
        raise InputError(f"{contig}:{start + 1}-{end} is not covered by the installed {genome.build} reference")
    span = genome.fetch(contig, start, end)
    if dna.ref and dna.kind != "ins" and dna.ref != span:
        raise ReferenceMismatch(
            f"the reference at {contig}:{start + 1} is {span[:24]}, not {dna.ref[:24]}", expected=dna.ref, found=span
        )
    if dna.kind == "sub":
        return normalize(genome, contig, start + 1, span, dna.alt)
    if dna.kind == "syn":
        raise InputError("this description says nothing changes, so there is no variant to annotate")
    if dna.kind == "delins":
        return normalize(genome, contig, start + 1, span, dna.alt)
    if dna.kind == "inv":
        if len(span) < 2:
            raise InputError("an inversion needs at least two bases")
        return normalize(genome, contig, start + 1, span, reverse_complement(span))
    if dna.kind == "del":
        base = genome.fetch(contig, start - 1, start)
        return normalize(genome, contig, start, base + span, base)
    if dna.kind == "dup":
        base = genome.fetch(contig, end - 1, end)
        return normalize(genome, contig, end, base, base + span)
    if dna.kind == "ins":
        if end - start != 2:
            raise InputError("an insertion must lie between two neighbouring positions")
        anchor = span[0]
        return normalize(genome, contig, start + 1, anchor, anchor + dna.alt)
    raise InputError(f"unsupported change type {dna.kind}")


def _size(variant: Variant) -> int:
    return max(len(variant.ref), len(variant.alt))


def _rank(candidates: list[Candidate]) -> list[Candidate]:
    ordered = sorted(candidates, key=lambda c: (not c.in_clinvar, c.bases_changed, c.variant.pos, c.variant.alt))
    return ordered[:MAX_CANDIDATES]
