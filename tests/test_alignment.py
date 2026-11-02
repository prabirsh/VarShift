"""Transcripts that don't match the genome base for base: numbering, translation, and RefSeq's real rows.

The toy transcript is ATG GCT CAT GAA TGG TAA, so Met Ala His Glu Trp stop. Each toy genome differs from it in
one way, and every expected name was worked out by hand. Each case is checked on both strands.
"""

import tempfile
import unittest
from pathlib import Path

from tests.helpers import write_fasta
from varsight.consequence import annotate
from varsight.errors import ResolveError, SetupError
from varsight.genome import Genome
from varsight.hgvs import parse_dna_change
from varsight.pipeline import builders
from varsight.resolver import _apply_on_transcript
from varsight.sequence import reverse_complement
from varsight.transcripts import Transcript, TxPosition
from varsight.variant import Variant

FIXTURES = Path(__file__).resolve().parent / "fixtures"
FLANK = "CCCCCCCCCC"
TRANSCRIPT = "ATGGCTCATGAATGGTAA"


class Toy:
    """A one-exon coding transcript on a small genome, described as it reads on the transcript's own strand."""

    def __init__(self, test, strand, pieces):
        """pieces are (bases on the genome, transcript index of the first base), or (bases, None) if it has none."""
        folder = tempfile.TemporaryDirectory()
        test.addCleanup(folder.cleanup)
        along_transcript = "".join(bases for bases, _ in pieces)
        self.strand = strand
        self.span = len(along_transcript)
        forward = FLANK + along_transcript + FLANK
        if strand < 0:
            forward = reverse_complement(forward)
        self.genome = Genome("hg19", write_fasta(Path(folder.name), "g.fa", {"chr1": forward}))
        test.addCleanup(self.genome.close)
        blocks = []
        made = 0
        for bases, index in pieces:
            if index is not None:
                blocks.append((*self.forward_span(made, made + len(bases)), index))
            made += len(bases)
        low, high = self.forward_span(0, self.span)
        self.tx = Transcript(
            "NM_1.1", "TOY", "chr1", strand, ((low, high),), cds_start=low, cds_end=high, protein="NP_1.1", blocks=tuple(blocks)
        )

    def forward_span(self, first: int, last: int) -> tuple[int, int]:
        """Genomic span of bases first to last, counted along the transcript's strand on the genome."""
        if self.strand > 0:
            return len(FLANK) + first, len(FLANK) + last
        return len(FLANK) + self.span - last, len(FLANK) + self.span - first

    def genomic(self, along: int) -> int:
        low, high = self.forward_span(along, along + 1)
        return low

    def change(self, along: int, ref: str, alt: str):
        """A base change at a place counted along the transcript's strand on the genome."""
        if self.strand < 0:
            ref, alt = reverse_complement(ref), reverse_complement(alt)
        return annotate(self.genome, Variant("hg19", "chr1", self.genomic(along) + 1, ref, alt), self.tx)


def on_both_strands(test_method):
    def run(self):
        for strand in (1, -1):
            with self.subTest(strand=strand):
                test_method(self, strand)

    run.__name__ = test_method.__name__
    return run


class ThreeExtraBasesInTheTranscript(unittest.TestCase):
    """The genome lacks CAT, the third codon, so it reads ATG GCT, then GAA TGG TAA."""

    def toy(self, strand):
        return Toy(self, strand, [("ATGGCT", 0), ("GAATGGTAA", 9)])

    @on_both_strands
    def test_numbering_after_the_gap_counts_the_missing_bases(self, strand):
        toy = self.toy(strand)
        effect = toy.change(9, "T", "A")  # T of TGG is transcript base 13 on a genome that has it as base 10
        self.assertEqual(effect.hgvs_c, "NM_1.1:c.13T>A")
        self.assertEqual(effect.hgvs_p, "NP_1.1:p.Trp5Arg")
        self.assertEqual(effect.effect, "missense_variant")

    @on_both_strands
    def test_numbering_before_the_gap_is_unchanged(self, strand):
        effect = self.toy(strand).change(4, "C", "A")  # GCT to GAT
        self.assertEqual(effect.hgvs_c, "NM_1.1:c.5C>A")
        self.assertEqual(effect.hgvs_p, "NP_1.1:p.Ala2Asp")

    @on_both_strands
    def test_a_name_on_the_far_side_of_the_gap_finds_its_position(self, strand):
        toy = self.toy(strand)
        variant = _apply_on_transcript(toy.genome, toy.tx, parse_dna_change("13T>A"))
        self.assertEqual(variant.pos, toy.genomic(9) + 1)

    @on_both_strands
    def test_a_name_inside_the_missing_bases_is_refused(self, strand):
        toy = self.toy(strand)
        with self.assertRaises(ResolveError) as caught:
            _apply_on_transcript(toy.genome, toy.tx, parse_dna_change("8A>G"))
        self.assertIn("isn't in this genome build", str(caught.exception))

    @on_both_strands
    def test_every_base_goes_there_and_back(self, strand):
        toy = self.toy(strand)
        for along in range(toy.span):
            genomic = toy.genomic(along)
            self.assertEqual(toy.tx.genomic_of_index(toy.tx.exon_index_of(genomic)), genomic)

    @on_both_strands
    def test_the_sequence_keeps_the_missing_bases_as_unknown(self, strand):
        toy = self.toy(strand)
        self.assertEqual(toy.tx.spliced(toy.genome), "ATGGCTNNNGAATGGTAA")
        self.assertEqual(toy.tx.tx_length, len(TRANSCRIPT))


class OneExtraBaseInTheTranscript(unittest.TestCase):
    """The genome lacks the A of CAT, which on its own would shift the reading frame."""

    def toy(self, strand):
        return Toy(self, strand, [("ATGGCTC", 0), ("TGAATGGTAA", 8)])

    @on_both_strands
    def test_the_reading_frame_is_kept_after_the_gap(self, strand):
        effect = self.toy(strand).change(11, "T", "A")
        self.assertEqual(effect.hgvs_c, "NM_1.1:c.13T>A")
        self.assertEqual(effect.hgvs_p, "NP_1.1:p.Trp5Arg")

    @on_both_strands
    def test_a_change_in_the_broken_codon_gets_no_protein_name(self, strand):
        effect = self.toy(strand).change(6, "C", "G")
        self.assertEqual(effect.hgvs_c, "NM_1.1:c.7C>G")
        self.assertIsNone(effect.hgvs_p)
        self.assertIn("transcript_gap_here", effect.flags)


class TwoExtraBasesInTheGenome(unittest.TestCase):
    """The genome has GG that the transcript lacks: ATG GCT CA[GG]T GAA TGG TAA."""

    def toy(self, strand):
        return Toy(self, strand, [("ATGGCTCA", 0), ("GG", None), ("TGAATGGTAA", 8)])

    @on_both_strands
    def test_numbering_after_the_gap_skips_those_bases(self, strand):
        effect = self.toy(strand).change(14, "T", "A")
        self.assertEqual(effect.hgvs_c, "NM_1.1:c.13T>A")
        self.assertEqual(effect.hgvs_p, "NP_1.1:p.Trp5Arg")

    @on_both_strands
    def test_a_change_in_those_bases_has_no_transcript_name(self, strand):
        effect = self.toy(strand).change(8, "G", "C")
        self.assertIsNone(effect.hgvs_c)
        self.assertIsNone(effect.hgvs_p)
        self.assertIn("transcript_gap_here", effect.flags)

    @on_both_strands
    def test_they_are_not_mistaken_for_an_intron(self, strand):
        toy = self.toy(strand)
        self.assertEqual(toy.tx.introns, ())
        self.assertNotIn("splice_donor_variant", toy.change(7, "A", "G").terms)


class TranscriptsWithoutGaps(unittest.TestCase):
    def test_a_transcript_without_an_alignment_counts_along_its_exons(self):
        tx = Transcript("NM_2.1", "PLAIN", "chr1", 1, ((10, 20), (30, 40)), cds_start=12, cds_end=38)
        self.assertEqual(tx.tx_length, 20)
        self.assertEqual(tx.exon_index_of(31), 11)
        self.assertEqual(tx.genomic_of_index(11), 31)
        self.assertEqual(tx.position_of(21), TxPosition("", 8, 2))
        self.assertEqual(tx.position_of(29), TxPosition("", 9, -1))
        self.assertEqual(tx.genomic_of(TxPosition("", 9, -1)), 29)


class RealRefSeqRows(unittest.TestCase):
    """RefSeq's own rows for four transcripts, checked against ClinVar's names, which the old numbering missed."""

    @classmethod
    def setUpClass(cls):
        sample = FIXTURES / "refseq_gaps.sample.gff"
        with open(sample, encoding="utf-8") as handle:
            alignments = builders.read_refseq_alignments(handle)
        cls.transcripts = {}
        with open(sample, encoding="utf-8") as handle:
            for record in builders.iter_refseq_transcripts(handle, "hg19"):
                row = builders._transcript_row(record, {}, alignments)
                exons, blocks = row[11], row[-1]
                strand = row[7]
                cls.transcripts[record["accession"]] = Transcript(
                    record["accession"],
                    record["gene"],
                    record["contig"],
                    strand,
                    tuple(tuple(int(n) for n in item.split("-")) for item in exons.split(",")),
                    cds_start=row[12],
                    cds_end=row[13],
                    blocks=tuple(tuple(int(n) for n in item.split("-")) for item in blocks.split(",")),
                )

    def test_positions_match_clinvar_in_both_directions(self):
        cases = [
            ("NM_001378454.1", 73675519, TxPosition("", 1865)),  # ALMS1, three extra bases, plus strand
            ("NM_001378454.1", 73800485, TxPosition("", 11481)),
            ("NM_000348.4", 31756486, TxPosition("", 502)),  # SRD5A2, one extra base, minus strand
            ("NM_000348.4", 31754536, TxPosition("", 548, -9)),
            ("NM_001372044.2", 51112850, TxPosition("", 8)),  # SHANK3, two bases missing from the transcript
            ("NM_001372044.2", 51150043, TxPosition("", 2336)),  # and after 39 bases that hg19 lacks
            ("NM_001145026.2", 80888844, TxPosition("", 1673)),  # a stretch of the transcript that hg19 lacks
        ]
        for accession, position, expected in cases:
            with self.subTest(accession=accession, position=position):
                tx = self.transcripts[accession]
                self.assertEqual(tx.position_of(position - 1), expected)
                self.assertEqual(tx.genomic_of(expected), position - 1)

    def test_skipped_bases_in_the_exon_rows_are_not_an_intron(self):
        tx = self.transcripts["NM_000348.4"]
        self.assertEqual(len(tx.exons), 5)
        self.assertEqual(tx.locate(31805881).kind, "exon")

    def test_the_real_rows_follow_the_transcript_on_the_minus_strand(self):
        tx = self.transcripts["NM_000348.4"]
        self.assertEqual(tx.aligned[0], (31805880, 31806007, 0))  # M127 comes first, from the 5' end
        self.assertEqual(tx.aligned[1], (31805689, 31805880, 128))  # then I1, then M191


class AlignmentRows(unittest.TestCase):
    ROW = "NC_000001.10\tRefSeq\tcDNA_match\t{start}\t{end}\t.\t{strand}\t.\tID=a;Target=NM_1.1 {first} {last} +{gap}\n"

    def read(self, **values):
        fields = {"start": 101, "end": 110, "strand": "+", "first": 1, "last": 10, "gap": ""} | values
        return builders.read_refseq_alignments([self.ROW.format(**fields)])

    def test_an_exact_alignment_is_not_kept(self):
        self.assertEqual(self.read(), {})

    def test_gaps_split_a_row_on_the_plus_strand(self):
        found = self.read(end=112, last=11, gap=";Gap=M4 I1 M3 D2 M3")
        exons, blocks = found[("NC_000001.10", "NM_1.1")][0]
        self.assertEqual(exons, [(100, 112)])
        self.assertEqual(blocks, [(100, 104, 0), (104, 107, 5), (109, 112, 8)])

    def test_gaps_split_a_row_on_the_minus_strand(self):
        found = self.read(end=112, last=11, strand="-", gap=";Gap=M4 I1 M3 D2 M3")
        _, blocks = found[("NC_000001.10", "NM_1.1")][0]
        self.assertEqual(blocks, [(108, 112, 0), (105, 108, 5), (100, 103, 8)])

    def test_a_stretch_of_the_transcript_with_no_alignment_is_kept(self):
        lines = [
            self.ROW.format(start=101, end=110, strand="+", first=1, last=10, gap=""),
            self.ROW.format(start=201, end=210, strand="+", first=41, last=50, gap=""),
        ]
        _, blocks = builders.read_refseq_alignments(lines)[("NC_000001.10", "NM_1.1")][0]
        self.assertEqual(blocks, [(100, 110, 0), (200, 210, 40)])

    def test_rows_on_other_contigs_are_left_out(self):
        line = self.ROW.format(start=101, end=110, strand="+", first=5, last=14, gap="").replace("NC_000001.10", "NT_113878.1")
        self.assertEqual(builders.read_refseq_alignments([line]), {})

    def test_a_changed_format_stops_the_build(self):
        with self.assertRaises(SetupError):
            self.read(gap=";Gap=M4 X1 M5")
        with self.assertRaises(SetupError):
            self.read(gap=";Gap=M4 I1 M3")  # doesn't add up to the row it describes
        with self.assertRaises(SetupError):
            builders.read_refseq_alignments([self.ROW.format(start=101, end=110, strand="+", first=1, last=10, gap="").replace(" +", " -")])

    def test_the_alignment_nearest_the_transcript_is_used(self):
        record = {"seqid": "chr", "accession": "NM_1.1", "exons": [(1000, 1100)]}
        far = ([(10, 110)], [(10, 110, 0)])
        near = ([(1000, 1100)], [(1000, 1050, 0), (1052, 1100, 50)])
        self.assertEqual(builders._alignment_for(record, {("chr", "NM_1.1"): [far, near]}), near)
        self.assertIsNone(builders._alignment_for(record, {("chr", "NM_1.1"): [far]}))
        self.assertIsNone(builders._alignment_for(record, {}))


if __name__ == "__main__":
    unittest.main()
