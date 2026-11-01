"""A small gene written by hand, so every expected HGVS name in the tests can be checked by eye.

    M   A   E   K   L  | intron |  G   W   P   *
    ATG GCT GAA AAA CTG GTAAGT.. GGC TGG CCC TAA

The same gene is available on either strand. A reverse-strand copy is the exact reverse complement,
so a correct engine must give identical c. and p. names for both.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from varsight.genome import Genome, write_fai
from varsight.sequence import reverse_complement
from varsight.transcripts import Transcript

UPSTREAM = "ACGTTGCAACGTAGCATCGA"
UTR5 = "TTCCA"
CDS1 = "ATGGCTGAAAAACTG"
INTRON = "GTAAGTCCCCCTTTTAAACAG"
CDS2 = "GGCTGGCCCTAA"
UTR3 = "GGATCCTAAAT"
DOWNSTREAM = "TCGATGCTACGTTGCAACGT"
FORWARD = UPSTREAM + UTR5 + CDS1 + INTRON + CDS2 + UTR3 + DOWNSTREAM
CONTIG = "chrT"


def write_fasta(folder: Path, name: str, records: dict[str, str], width: int = 10) -> Path:
    path = folder / name
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        for contig, seq in records.items():
            handle.write(f">{contig}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")
    write_fai(path)
    return path


class ToyGene:
    """The toy gene and a genome holding it, on the chosen strand."""

    def __init__(self, strand: int):
        self.strand = strand
        self._folder = tempfile.TemporaryDirectory()
        sequence = FORWARD if strand > 0 else reverse_complement(FORWARD)
        fasta = write_fasta(Path(self._folder.name), "toy.fa", {CONTIG: sequence})
        self.genome = Genome("hg38", fasta)
        length = len(FORWARD)

        def place(start: int, end: int) -> tuple[int, int]:
            return (start, end) if strand > 0 else (length - end, length - start)

        exon1_start = len(UPSTREAM)
        exon1_end = exon1_start + len(UTR5) + len(CDS1)
        exon2_start = exon1_end + len(INTRON)
        exon2_end = exon2_start + len(CDS2) + len(UTR3)
        cds_start = exon1_start + len(UTR5)
        cds_end = exon2_start + len(CDS2)
        low, high = place(cds_start, cds_end)
        self.transcript = Transcript(
            accession="NM_TOY.1",
            gene="TOY",
            contig=CONTIG,
            strand=strand,
            exons=(place(exon1_start, exon1_end), place(exon2_start, exon2_end)),
            cds_start=low,
            cds_end=high,
            biotype="mRNA",
            protein="NP_TOY.1",
        )
        self.length = length

    def close(self) -> None:
        self.genome.close()
        self._folder.cleanup()
