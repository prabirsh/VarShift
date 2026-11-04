"""Transcript model and the coordinate arithmetic between genome and c./n. positions.

All genomic coordinates here are 0-based half-open unless a name says otherwise. Exons are stored
in transcript order, so exon 0 is the 5' exon on either strand.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

from .errors import ResolveError
from .genome import Genome
from .sequence import reverse_complement


class TxPosition(NamedTuple):
    """A position in c. or n. numbering. region is '' (CDS or n.), '-' (5' of CDS) or '*' (3' of it)."""

    region: str
    base: int
    offset: int = 0

    def __str__(self) -> str:
        offset = f"{self.offset:+d}" if self.offset else ""
        return f"{self.region}{self.base}{offset}"


class Location(NamedTuple):
    kind: str  # exon, intron, upstream or downstream
    index: int  # exon or intron number, 0-based in transcript order


@dataclass(frozen=True)
class Transcript:
    accession: str
    gene: str
    contig: str
    strand: int
    exons: tuple[tuple[int, int], ...]
    cds_start: int | None = None
    cds_end: int | None = None
    biotype: str = ""
    protein: str = ""
    source: str = "RefSeq"
    gene_id: str = ""
    hgnc_id: str = ""
    mane: str = ""
    mane_accession: str = ""
    refseq_select: bool = False
    genome_note: str = ""
    translation_exception: bool = False
    # RefSeq's alignment, only where the transcript doesn't match the genome base for base.
    blocks: tuple[tuple[int, int, int], ...] = ()

    @property
    def coding(self) -> bool:
        return self.cds_start is not None

    @property
    def base_accession(self) -> str:
        return self.accession.split(".")[0]

    @cached_property
    def tx_low(self) -> int:
        return min(start for start, _ in self.exons)

    @cached_property
    def tx_high(self) -> int:
        return max(end for _, end in self.exons)

    @cached_property
    def aligned(self) -> tuple[tuple[int, int, int], ...]:
        """Stretches of genome matching the transcript base for base: (start, end, index of the 5' base)."""
        if self.blocks:
            return self.blocks
        stretches, total = [], 0
        for start, end in self.exons:
            stretches.append((start, end, total))
            total += end - start
        return tuple(stretches)

    @cached_property
    def tx_length(self) -> int:
        start, end, index = self.aligned[-1]
        return index + end - start

    @cached_property
    def introns(self) -> tuple[tuple[int, int], ...]:
        """Genomic spans of the introns; intron i lies between exon i and exon i + 1."""
        spans = []
        for first, second in zip(self.exons, self.exons[1:]):
            if self.strand > 0:
                spans.append((first[1], second[0]))
            else:
                spans.append((second[1], first[0]))
        return tuple(spans)

    @cached_property
    def cds_tx_range(self) -> tuple[int, int]:
        """Transcript indices of the first CDS base and one past the last, as (start, end)."""
        if not self.coding:
            raise ValueError(f"{self.accession} has no CDS")
        low = self.exon_index_of(self.cds_start)
        high = self.exon_index_of(self.cds_end - 1)
        if low is None or high is None:
            raise ResolveError(f"{self.accession}: CDS ends do not fall on exons")
        return (low, high + 1) if self.strand > 0 else (high, low + 1)

    def exon_index_of(self, genomic: int) -> int | None:
        """Transcript index of a genomic base, or None when the transcript has no such base."""
        for start, end, index in self.aligned:
            if start <= genomic < end:
                return index + (genomic - start if self.strand > 0 else end - 1 - genomic)
        return None

    def genomic_of_index(self, index: int) -> int:
        """Genomic base for a transcript index. Indices past either end extend straight out."""
        if index < 0:
            five_prime = self.exons[0][0] if self.strand > 0 else self.exons[0][1] - 1
            return five_prime + (index - self.aligned[0][2]) * self.strand
        for start, end, first in self.aligned:
            if first <= index < first + end - start:
                within = index - first
                return start + within if self.strand > 0 else end - 1 - within
        if index < self.tx_length:
            raise ResolveError(f"base {index + 1} of {self.accession} isn't in this genome build; RefSeq's alignment has a gap there")
        three_prime = self.exons[-1][1] - 1 if self.strand > 0 else self.exons[-1][0]
        return three_prime + (index - self.tx_length + 1) * self.strand

    def _exon_indices(self, number: int) -> tuple[int | None, int | None]:
        """Transcript indices of an exon's 5' and 3' bases."""
        start, end = self.exons[number]
        five, three = (start, end - 1) if self.strand > 0 else (end - 1, start)
        return self.exon_index_of(five), self.exon_index_of(three)

    def locate(self, genomic: int) -> Location:
        for number, (start, end) in enumerate(self.exons):
            if start <= genomic < end:
                return Location("exon", number)
        for number, (start, end) in enumerate(self.introns):
            if start <= genomic < end:
                return Location("intron", number)
        before = genomic < self.tx_low
        if (before and self.strand > 0) or (not before and self.strand < 0):
            return Location("upstream", 0)
        return Location("downstream", 0)

    def _label_index(self, index: int) -> TxPosition | None:
        if not self.coding:
            if 0 <= index < self.tx_length:
                return TxPosition("", index + 1)
            return None
        cds_first, cds_stop = self.cds_tx_range
        if index < cds_first:
            return TxPosition("-", cds_first - index)
        if index < cds_stop:
            return TxPosition("", index - cds_first + 1)
        return TxPosition("*", index - cds_stop + 1)

    def position_of(self, genomic: int) -> TxPosition | None:
        """c. or n. position of a genomic base, using the HGVS rule for intronic numbering."""
        where = self.locate(genomic)
        if where.kind == "exon":
            index = self.exon_index_of(genomic)
            return None if index is None else self._label_index(index)
        if where.kind == "intron":
            intron_start, intron_end = self.introns[where.index]
            if self.strand > 0:
                from_donor = genomic - intron_start + 1
                from_acceptor = intron_end - genomic
            else:
                from_donor = intron_end - genomic
                from_acceptor = genomic - intron_start + 1
            if from_donor <= from_acceptor:
                label = self._label_index(self._exon_indices(where.index)[1])
                return TxPosition(label.region, label.base, from_donor)
            label = self._label_index(self._exon_indices(where.index + 1)[0])
            return TxPosition(label.region, label.base, -from_acceptor)
        distance = (self.tx_low - genomic) if genomic < self.tx_low else (genomic - self.tx_high + 1)
        if where.kind == "upstream":
            return self._label_index(self.aligned[0][2] - distance)
        return self._label_index(self.tx_length - 1 + distance)

    def genomic_of(self, position: TxPosition) -> int:
        """Genomic base of a c. or n. position. Raises ResolveError if the position cannot exist."""
        if self.coding:
            cds_first, cds_stop = self.cds_tx_range
            if position.region == "-":
                index = cds_first - position.base
            elif position.region == "*":
                index = cds_stop + position.base - 1
            else:
                if position.base > cds_stop - cds_first:
                    raise ResolveError(
                        f"c.{position} is past the end of the {cds_stop - cds_first} bp CDS of {self.accession}"
                    )
                index = cds_first + position.base - 1
        else:
            if position.region:
                raise ResolveError(f"{self.accession} is non-coding; n. numbering has no {position.region} region")
            if position.base > self.tx_length:
                raise ResolveError(f"n.{position} is past the end of {self.accession}")
            index = position.base - 1
        genomic = self.genomic_of_index(index)
        if not position.offset:
            return genomic
        if not 0 <= index < self.tx_length:
            raise ResolveError(f"{position} has an intron offset but lies outside {self.accession}")
        exon_number = next(
            number
            for number, (start, end) in enumerate(self.exons)
            if start <= genomic < end
        )
        first_of_exon, last_of_exon = self._exon_indices(exon_number)
        if position.offset > 0 and index != last_of_exon:
            raise ResolveError(f"{position}: an offset of +N is only valid on the last base of an exon")
        if position.offset < 0 and index != first_of_exon:
            raise ResolveError(f"{position}: an offset of -N is only valid on the first base of an exon")
        target = genomic + position.offset * self.strand
        if self.locate(target).kind != "intron":
            raise ResolveError(f"{position} runs past the intron on {self.accession}")
        return target

    def spliced(self, genome: Genome) -> str:
        """The transcript's sequence. Bases the genome doesn't have are N, which keeps later codons in frame."""
        parts, made = [], 0
        for start, end, index in self.aligned:
            part = genome.fetch(self.contig, start, end)
            if len(part) != end - start:
                raise ResolveError(f"{self.accession}: exon {start}-{end} is not in the installed reference")
            parts.append("N" * (index - made))
            parts.append(reverse_complement(part) if self.strand < 0 else part)
            made = index + end - start
        return "".join(parts)

    def default_rank(self) -> tuple:
        """Sort key for choosing one transcript to show first; lower is better."""
        cds_length = 0
        if self.coding:
            first, stop = self.cds_tx_range
            cds_length = stop - first
        return (
            {"MANE Select": 0, "MANE Plus Clinical": 1}.get(self.mane, 2),
            0 if self.refseq_select else 1,
            0 if self.source == "RefSeq" else 1,
            0 if self.coding else 1,
            -cds_length,
            -self.tx_length,
            self.accession,
        )
