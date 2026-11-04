"""A small gene written by hand, so every expected HGVS name in the tests can be checked by eye.

    M   A   E   K   L  | intron |  G   W   P   *
    ATG GCT GAA AAA CTG GTAAGT.. GGC TGG CCC TAA

The same gene is available on either strand. A reverse-strand copy is the exact reverse complement,
so a correct engine must give identical c. and p. names for both.
"""

from __future__ import annotations

import struct
import tempfile
import zlib
from pathlib import Path

from varshift.genome import Genome, write_fai
from varshift.sequence import reverse_complement
from varshift.transcripts import Transcript

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


def tiny_bigwig() -> bytes:
    """A real little-endian bigWig holding one block of each kind, with a two-level index."""
    header_size = 64
    tree_at = header_size
    key = b"chr7\0\0"
    tree = struct.pack("<IIIIQQ", 0x78CA8C91, 256, len(key), 8, 1, 0) + struct.pack("<BBH", 1, 0, 1) + key + struct.pack("<II", 0, 1000)
    data_at = tree_at + len(tree)
    body = bytearray(struct.pack("<Q", 3))

    def block(start, end, step, span, kind, items, packed):
        section = struct.pack("<IIIIIBBH", 0, start, end, step, span, kind, 0, items) + packed
        compressed = zlib.compress(section)
        offset = data_at + len(body)
        body.extend(compressed)
        return (0, start, 0, end, offset, len(compressed))

    leaves = [
        block(100, 103, 1, 1, 3, 3, struct.pack("<fff", 1.5, 2.5, 3.5)),
        block(200, 206, 0, 1, 2, 2, struct.pack("<IfIf", 200, -1.0, 205, 7.25)),
        block(300, 310, 0, 0, 1, 1, struct.pack("<IIf", 300, 310, 9.0)),
    ]
    index_at = data_at + len(body)
    index = struct.pack("<IIQIIIIQII", 0x2468ACE0, 256, 3, 0, 100, 0, 310, index_at, 1024, 0)
    leaf_at = index_at + len(index) + 4 + 24
    root = struct.pack("<BBH", 0, 0, 1) + struct.pack("<IIIIQ", 0, 100, 0, 310, leaf_at)
    leaf = struct.pack("<BBH", 1, 0, len(leaves)) + b"".join(struct.pack("<IIIIQQ", *item) for item in leaves)
    header = struct.pack("<IHHQQQHHQQIQ", 0x888FFC26, 4, 0, tree_at, data_at, index_at, 0, 0, 0, 0, 4096, 0)
    return header + tree + bytes(body) + index + root + leaf
