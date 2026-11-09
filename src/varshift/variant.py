"""The genomic variant record and the normalisation that produces it."""

from __future__ import annotations

from dataclasses import dataclass

from .builds import contig_length
from .errors import InputError, ReferenceMismatch
from .genome import Genome
from .sequence import is_dna

MAX_ALLELE_LENGTH = 100_000
_LEFT_ALIGN_BLOCK = 64


@dataclass(frozen=True)
class Variant:
    """A normalised variant: VCF-style 1-based position, left-aligned and minimal."""

    build: str
    contig: str
    pos: int
    ref: str
    alt: str

    @property
    def key(self) -> str:
        return f"{self.contig}-{self.pos}-{self.ref}-{self.alt}"

    @property
    def end(self) -> int:
        return self.pos + len(self.ref) - 1

    def edit(self) -> tuple[int, int, str]:
        """The change without padding: 0-based half-open reference span and the inserted bases."""
        ref, alt = self.ref, self.alt
        suffix = 0
        while suffix < min(len(ref), len(alt)) and ref[-1 - suffix] == alt[-1 - suffix]:
            suffix += 1
        if suffix:
            ref, alt = ref[:-suffix], alt[:-suffix]
        prefix = 0
        while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
            prefix += 1
        start = self.pos - 1 + prefix
        return start, start + len(ref) - prefix, alt[prefix:]


def normalize(genome: Genome, contig: str, pos: int, ref: str, alt: str) -> Variant:
    """Check REF against the genome, then trim and left-align.

    Follows Tan, Abecasis and Kang (2015), the algorithm behind vt normalize and bcftools norm.
    """
    ref = ref.upper()
    alt = alt.upper()
    _check_allele(ref, "reference")
    _check_allele(alt, "alternate")
    length = contig_length(contig, genome.build)
    if pos < 1 or (length is not None and pos + len(ref) - 1 > length):
        raise InputError(f"position {contig}:{pos} is outside the contig on {genome.build}")
    genome.require(contig)
    if not genome.covers(contig, pos - 1, pos - 1 + len(ref)):
        raise InputError(f"{contig}:{pos} is not covered by the installed {genome.build} reference")
    found = genome.fetch(contig, pos - 1, pos - 1 + len(ref))
    if found != ref:
        raise ReferenceMismatch(
            f"REF {_shorten(ref)} does not match {genome.build} {contig}:{pos}, which has {_shorten(found)}",
            expected=ref,
            found=found,
        )
    if ref == alt:
        raise InputError(f"REF and ALT are the same at {contig}:{pos}, so nothing changes")

    pos, ref, alt = _left_align(genome, contig, pos, ref, alt)
    return Variant(genome.build, contig, pos, ref, alt)


def _left_align(genome: Genome, contig: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    block = ""
    block_start = pos - 1
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif (not ref or not alt) and pos > 1:
            if pos - 2 < block_start or not block:
                block_start = max(pos - 1 - _LEFT_ALIGN_BLOCK, 0)
                block = genome.fetch(contig, block_start, pos - 1)
                if not block:
                    break  # at the start of a demo window; left-aligned as far as data allows
            base = block[pos - 2 - block_start]
            ref, alt, pos = base + ref, base + alt, pos - 1
        else:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _check_allele(allele: str, label: str) -> None:
    if len(allele) > MAX_ALLELE_LENGTH:
        raise InputError(f"{label} allele is longer than {MAX_ALLELE_LENGTH} bases")
    if not is_dna(allele):
        raise InputError(f"{label} allele {_shorten(allele)!r} is not plain DNA (A, C, G, T, N)")


def _shorten(seq: str, limit: int = 24) -> str:
    return seq if len(seq) <= limit else f"{seq[:limit]}...({len(seq)} bases)"
