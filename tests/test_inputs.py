"""Reading variant descriptions in the many shapes people actually send them."""

import unittest

from varsight.errors import InputError
from varsight.hgvs import parse_dna_change, parse_protein_change
from varsight.inputs import MAX_INPUT_LENGTH, parse_input
from varsight.transcripts import TxPosition


class Protein(unittest.TestCase):
    def test_the_same_change_in_many_forms(self):
        for text in (
            "BRAF V600E",
            "BRAF:V600E",
            "BRAF p.V600E",
            "BRAF:p.V600E",
            "BRAF p.Val600Glu",
            "BRAF p.(Val600Glu)",
            "BRAF (V600E)",
            "  BRAF   V600E  ",
            "braf V600E",
        ):
            with self.subTest(text=text):
                parsed = parse_input(text)
                self.assertEqual(parsed.kind, "protein")
                self.assertEqual(parsed.gene.upper(), "BRAF")
                self.assertEqual((parsed.protein.kind, parsed.protein.ref, parsed.protein.pos, parsed.protein.alt), ("sub", "V", 600, "E"))

    def test_protein_accession(self):
        parsed = parse_input("NP_004324.2:p.Val600Glu")
        self.assertEqual((parsed.kind, parsed.protein_accession), ("protein", "NP_004324.2"))

    def test_gene_can_come_from_the_option(self):
        self.assertEqual(parse_input("V600E", gene="BRAF").gene, "BRAF")

    def test_no_gene_says_what_is_missing(self):
        with self.assertRaises(InputError) as caught:
            parse_input("V600E")
        self.assertIn("gene", str(caught.exception))

    def test_change_types(self):
        cases = {
            "R132*": ("sub", "*"),
            "R132X": ("sub", "*"),
            "Arg132Ter": ("sub", "*"),
            "V600=": ("syn", ""),
            "E23fs": ("fs", ""),
            "Glu23ValfsTer17": ("fs", "V"),
            "E23Vfs*17": ("fs", "V"),
            "E746_A750del": ("del", ""),
            "Glu746del": ("del", ""),
            "L5dup": ("dup", ""),
            "Lys10_Ala11insGly": ("ins", "G"),
            "K10delinsNG": ("delins", "NG"),
            "Ter110GlnextTer17": ("ext", "Q"),
            "M1?": ("start", ""),
        }
        for text, (kind, alt) in cases.items():
            with self.subTest(text=text):
                change = parse_protein_change(text)
                self.assertEqual((change.kind, change.alt), (kind, alt))

    def test_clinvar_style(self):
        self.assertEqual(parse_protein_change("V600E").clinvar_style, "p.Val600Glu")
        self.assertEqual(parse_protein_change("R132*").clinvar_style, "p.Arg132Ter")
        self.assertEqual(parse_protein_change("E23Vfs*17").clinvar_style, "p.Glu23fs")
        self.assertEqual(parse_protein_change("E746_A750del").clinvar_style, "p.Glu746_Ala750del")

    def test_nonsense_is_rejected(self):
        for text in ("Z", "V0E", "V600", "600E", "A10_B5del"):
            with self.subTest(text=text), self.assertRaises(InputError):
                parse_protein_change(text)


class Coding(unittest.TestCase):
    def test_forms(self):
        cases = [
            ("NM_004333.6:c.1799T>A", "NM_004333.6", ""),
            ("NM_004333:c.1799T>A", "NM_004333", ""),
            ("NM_004333.6(BRAF):c.1799T>A", "NM_004333.6", "BRAF"),
            ("NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)", "NM_004333.6", "BRAF"),
            ("NM_004333.6: c.1799 T>A", "NM_004333.6", ""),
            ("BRAF c.1799T>A", "", "BRAF"),
            ("BRAF:c.1799T>A", "", "BRAF"),
            ("ENST00000288602.11:c.1799T>A", "ENST00000288602.11", ""),
        ]
        for text, tx, gene in cases:
            with self.subTest(text=text):
                parsed = parse_input(text)
                self.assertEqual((parsed.kind, parsed.transcript, parsed.gene), ("coding", tx, gene))
                self.assertEqual((parsed.dna.kind, parsed.dna.ref, parsed.dna.alt), ("sub", "T", "A"))

    def test_bare_c_needs_a_gene(self):
        with self.assertRaises(InputError):
            parse_input("c.1799T>A")
        self.assertEqual(parse_input("c.1799T>A", gene="BRAF").gene, "BRAF")

    def test_change_types(self):
        cases = {
            "1799T>A": ("sub", TxPosition("", 1799), TxPosition("", 1799)),
            "68_69del": ("del", TxPosition("", 68), TxPosition("", 69)),
            "68_69delAG": ("del", TxPosition("", 68), TxPosition("", 69)),
            "100dup": ("dup", TxPosition("", 100), TxPosition("", 100)),
            "100_101insAT": ("ins", TxPosition("", 100), TxPosition("", 101)),
            "100_102delinsAT": ("delins", TxPosition("", 100), TxPosition("", 102)),
            "100delTinsAT": ("delins", TxPosition("", 100), TxPosition("", 100)),
            "100_103inv": ("inv", TxPosition("", 100), TxPosition("", 103)),
            "-21G>A": ("sub", TxPosition("-", 21), TxPosition("-", 21)),
            "*5C>T": ("sub", TxPosition("*", 5), TxPosition("*", 5)),
            "1234+1G>A": ("sub", TxPosition("", 1234, 1), TxPosition("", 1234, 1)),
            "1235-2A>G": ("sub", TxPosition("", 1235, -2), TxPosition("", 1235, -2)),
            "-21-5G>A": ("sub", TxPosition("-", 21, -5), TxPosition("-", 21, -5)),
            "*5+3G>A": ("sub", TxPosition("*", 5, 3), TxPosition("*", 5, 3)),
        }
        for text, (kind, first, last) in cases.items():
            with self.subTest(text=text):
                change = parse_dna_change(text)
                self.assertEqual((change.kind, change.first, change.last), (kind, first, last))

    def test_an_insertion_needs_two_positions(self):
        with self.assertRaises(InputError):
            parse_dna_change("100insA")

    def test_rubbish_is_rejected(self):
        for text in ("", "abc", "100", "100T>", "T>A", "100T>Z", "1_2_3del"):
            with self.subTest(text=text), self.assertRaises(InputError):
                parse_dna_change(text)


class Genomic(unittest.TestCase):
    def test_hgvs_forms(self):
        for text, contig in (
            ("chr7:g.140453136A>T", "chr7"),
            ("7:g.140453136A>T", "7"),
            ("NC_000007.13:g.140453136A>T", "NC_000007.13"),
        ):
            with self.subTest(text=text):
                parsed = parse_input(text)
                self.assertEqual((parsed.kind, parsed.contig, parsed.dna.first.base), ("genomic", contig, 140453136))

    def test_refseq_accession_decides_the_build(self):
        self.assertEqual(parse_input("NC_000007.13:g.1A>T").build_hint, "hg19")
        self.assertEqual(parse_input("NC_000007.14:g.1A>T").build_hint, "hg38")
        self.assertEqual(parse_input("chr7:g.1A>T").build_hint, "")

    def test_position_forms(self):
        for text in (
            "chr7-140453136-A-T",
            "7-140453136-A-T",
            "chr7:140453136:A:T",
            "chr7:140453136 A>T",
            "chr7:140453136:A>T",
            "7 140453136 A T",
            "chr7_140453136_A_T",
            "chr7-140453136-a-t",
        ):
            with self.subTest(text=text):
                parsed = parse_input(text)
                self.assertEqual((parsed.kind, parsed.pos, parsed.ref, parsed.alt), ("position", 140453136, "A", "T"))

    def test_rsid(self):
        parsed = parse_input("rs113488022")
        self.assertEqual((parsed.kind, parsed.rsid), ("rsid", "113488022"))
        self.assertEqual(parse_input("RS113488022").rsid, "113488022")

    def test_mitochondrial(self):
        parsed = parse_input("NC_012920.1:m.3243A>G")
        self.assertEqual((parsed.kind, parsed.level), ("genomic", "m"))


class Limits(unittest.TestCase):
    def test_empty(self):
        for text in ("", "   ", '""'):
            with self.assertRaises(InputError):
                parse_input(text)

    def test_too_long(self):
        with self.assertRaises(InputError):
            parse_input("A" * (MAX_INPUT_LENGTH + 1))

    def test_hostile_text_is_only_ever_an_input_error(self):
        for text in (
            "'; DROP TABLE clinvar; " + "-" * 2,
            "../../etc/passwd",
            "$(rm -rf /)",
            "BRAF V600E; rm -rf /",
            "<script>alert(1)</script>",
            "A" * 1500 + ">" + "T" * 400,
            "\x00\x01\x02",
        ):
            with self.subTest(text=text[:30]), self.assertRaises(InputError):
                parse_input(text)


if __name__ == "__main__":
    unittest.main()
