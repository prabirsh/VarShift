"""HGVS naming and effect, checked on both strands against names worked out by hand."""

import unittest
from dataclasses import replace

from tests.helpers import ToyGene
from varsight.consequence import annotate
from varsight.hgvs import parse_dna_change
from varsight.resolver import _apply_on_transcript
from varsight.variant import Variant

# Input c. change, expected c. name after normalising, expected p. name, expected top effect.
CASES = [
    ("5C>A", "5C>A", "Ala2Asp", "missense_variant"),
    ("6T>C", "6T>C", "Ala2=", "synonymous_variant"),
    ("21G>A", "21G>A", "Trp7Ter", "stop_gained"),
    ("1A>G", "1A>G", "Met1?", "start_lost"),
    ("25T>C", "25T>C", "Ter9GlnextTer3", "stop_lost"),
    ("27A>G", "27A>G", "Ter9=", "stop_retained_variant"),
    ("4del", "4del", "Ala2LeufsTer?", "frameshift_variant"),
    ("7_9del", "7_9del", "Glu3del", "inframe_deletion"),
    # A run of five As (c.8 to c.12): a deletion anywhere in it is named at the 3' end.
    ("8_10del", "10_12del", "Lys4del", "inframe_deletion"),
    ("12dup", "12dup", "Leu5ThrfsTer?", "frameshift_variant"),
    ("8dup", "12dup", "Leu5ThrfsTer?", "frameshift_variant"),
    ("10_12dup", "10_12dup", "Lys4dup", "inframe_insertion"),
    ("9_10insGGG", "9_10insGGG", "Glu3_Lys4insGly", "inframe_insertion"),
    ("5_6delinsAA", "5_6delinsAA", "Ala2Glu", "missense_variant"),
    ("7_9inv", "7_9inv", "Glu3Phe", "missense_variant"),
    ("15+1G>A", "15+1G>A", None, "splice_donor_variant"),
    ("15+2T>C", "15+2T>C", None, "splice_donor_variant"),
    ("15+4A>T", "15+4A>T", None, "splice_region_variant"),
    ("15+11C>T", "15+11C>T", None, "intron_variant"),
    ("16-10T>G", "16-10T>G", None, "intron_variant"),
    ("16-1G>C", "16-1G>C", None, "splice_acceptor_variant"),
    ("16-2A>T", "16-2A>T", None, "splice_acceptor_variant"),
    # The 3' rule crosses an exon-intron junction when the bases allow it, as ClinVar's names do.
    ("15del", "15+1del", None, "splice_donor_variant"),
    ("15dup", "15+1dup", None, "splice_donor_variant"),
    ("16-1del", "17del", "Gly6AlafsTer?", "splice_acceptor_variant"),
    # Repeat notation: the run of five As reduced to four, and increased to six.
    ("8A[4]", "12del", "Lys4AsnfsTer?", "frameshift_variant"),
    ("8A[6]", "12dup", "Leu5ThrfsTer?", "frameshift_variant"),
    # A stop codon inside an insertion is a truncation, not an in-frame insertion.
    ("9_10insTAA", "9_10insTAA", "Glu3_Lys4insTer", "stop_gained"),
    ("-1A>G", "-1A>G", None, "5_prime_UTR_variant"),
    ("-5T>C", "-5T>C", None, "5_prime_UTR_variant"),
    ("*1G>A", "*1G>A", None, "3_prime_UTR_variant"),
    ("*11T>A", "*11T>A", None, "3_prime_UTR_variant"),
]


class NamingBothStrands(unittest.TestCase):
    def _check(self, strand):
        gene = ToyGene(strand)
        self.addCleanup(gene.close)
        for change, want_c, want_p, want_effect in CASES:
            with self.subTest(change=change, strand=strand):
                variant = _apply_on_transcript(gene.genome, gene.transcript, parse_dna_change(change))
                effect = annotate(gene.genome, variant, gene.transcript)
                self.assertEqual(effect.hgvs_c, f"NM_TOY.1:c.{want_c}")
                self.assertEqual(effect.hgvs_p, f"NP_TOY.1:p.{want_p}" if want_p else None)
                self.assertEqual(effect.effect, want_effect)

    def test_forward_strand(self):
        self._check(1)

    def test_reverse_strand(self):
        self._check(-1)


class Regions(unittest.TestCase):
    def setUp(self):
        self.gene = ToyGene(1)
        self.addCleanup(self.gene.close)

    def _effect(self, change):
        variant = _apply_on_transcript(self.gene.genome, self.gene.transcript, parse_dna_change(change))
        return annotate(self.gene.genome, variant, self.gene.transcript)

    def test_exon_and_intron_numbers(self):
        self.assertEqual(self._effect("5C>A").exon, 1)
        self.assertEqual(self._effect("21G>A").exon, 2)
        self.assertEqual(self._effect("15+11C>T").intron, 1)

    def test_last_exon_base_is_also_splice_region(self):
        effect = self._effect("15G>A")
        self.assertIn("splice_region_variant", effect.terms)
        self.assertEqual(effect.hgvs_p, "NP_TOY.1:p.Leu5=")

    def test_amino_acid_position(self):
        self.assertEqual(self._effect("5C>A").aa_position, 2)
        self.assertEqual(self._effect("21G>A").aa_position, 7)

    def test_short_protein_form(self):
        self.assertEqual(self._effect("5C>A").protein_short, "A2D")
        self.assertEqual(self._effect("21G>A").protein_short, "W7*")

    def test_upstream_position_names(self):
        # One base before the transcript starts: c.-6 on this gene.
        variant = Variant("hg38", "chrT", 20, "A", "C")
        result = annotate(self.gene.genome, variant, self.gene.transcript)
        self.assertEqual(result.effect, "upstream_gene_variant")
        self.assertEqual(result.hgvs_c, "NM_TOY.1:c.-6A>C")

    def test_downstream_position_names(self):
        last = self.gene.transcript.tx_high
        base = self.gene.genome.fetch("chrT", last, last + 1)
        alt = "A" if base != "A" else "C"
        result = annotate(self.gene.genome, Variant("hg38", "chrT", last + 1, base, alt), self.gene.transcript)
        self.assertEqual(result.effect, "downstream_gene_variant")
        self.assertEqual(result.hgvs_c, f"NM_TOY.1:c.*12{base}>{alt}")

    def test_note_on_transcript_becomes_flag(self):
        tx = replace(self.gene.transcript, genome_note="The RefSeq transcript has 1 substitution")
        result = annotate(self.gene.genome, Variant("hg38", "chrT", 30, *self._alleles(30)), tx)
        self.assertIn("transcript_differs_from_genome", result.flags)

    def _alleles(self, pos):
        base = self.gene.genome.fetch("chrT", pos - 1, pos)
        return base, ("A" if base != "A" else "C")


if __name__ == "__main__":
    unittest.main()
