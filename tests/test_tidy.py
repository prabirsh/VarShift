"""Step 1 of the input assistant: obvious slips are fixed, and nothing that already reads is changed."""

import unittest
from pathlib import Path

from varshift.inputs import parse_input
from varshift.tidy import Reading, tidy

DEMO_LIST = Path(__file__).resolve().parents[1] / "demo" / "variants.txt"


class Fixes(unittest.TestCase):
    CASES = [
        # typed, read as, words the note must contain
        ("braf val600glu", "braf Val600Glu", "standard case"),
        ("BRAF VAL600GLU", "BRAF Val600Glu", "standard case"),
        ("BRAF p.v600e", "BRAF p.V600E", "standard case"),
        ("BRAF v600e", "BRAF V600E", "standard case"),
        ("BRAF val600=", "BRAF Val600=", "standard case"),
        ("TP53 arg196*", "TP53 Arg196*", "standard case"),
        ("EGFR e746_a750del", "EGFR E746_A750del", "standard case"),
        ("EGFR GLU746_ALA750DEL", "EGFR Glu746_Ala750del", "standard case"),
        ("BRCA1 e23fs", "BRCA1 E23fs", "standard case"),
        ("NP_004324.2:p.val600glu", "NP_004324.2:p.Val600Glu", "standard case"),
        ("BRAF(V600E)", "BRAF V600E", "gene and change"),
        ("BRAF c.1799T-A", "BRAF c.1799T>A", "substitution"),
        ("BRAF c.1799T->A", "BRAF c.1799T>A", "substitution"),
        ("BRAF c.1799T\u2192A", "BRAF c.1799T>A", "substitution"),
        ("BRAF c.1799T/A", "BRAF c.1799T>A", "substitution"),
        ("BRAF c.1799T\uff1eA", "BRAF c.1799T>A", "pasted"),
        ("NM_004333.6:1799T>A", "NM_004333.6:c.1799T>A", "missing c."),
        ("NR_046018.2:100A>G", "NR_046018.2:n.100A>G", "missing n."),
        ("NC_000007.13:140453136A>T", "NC_000007.13:g.140453136A>T", "missing g."),
        ("chr7:140453136A>T", "chr7:g.140453136A>T", "missing g."),
        ("chrM:3243A>G", "chrM:m.3243A>G", "missing m."),
        ("NC_012920.1:3243A>G", "NC_012920.1:m.3243A>G", "missing m."),
        ("chr7:140,453,136 A>T", "chr7:140453136 A>T", "thousands"),
        ("\u201cBRAF V600E\u201d", "BRAF V600E", "pasted"),
        ("BRAF V600E.", "BRAF V600E", "pasted"),
        ("chr7\u2013140453136\u2013A\u2013T", "chr7-140453136-A-T", "pasted"),
    ]

    def test_each_slip_is_fixed_and_said(self):
        for typed, expected, note in self.CASES:
            with self.subTest(typed=ascii(typed)):
                reading = tidy(typed)
                self.assertEqual(reading.text, expected)
                self.assertTrue(reading.changed)
                self.assertTrue(any(note in line for line in reading.notes), reading.notes)

    def test_every_fix_now_reads(self):
        for typed, _, _ in self.CASES:
            with self.subTest(typed=ascii(typed)):
                parse_input(tidy(typed).text)

    def test_two_slips_at_once_give_two_notes(self):
        reading = tidy("NM_004333.6:1799T-A")
        self.assertEqual(reading.text, "NM_004333.6:c.1799T>A")
        self.assertEqual(len(reading.notes), 2)

    def test_a_gene_given_separately_counts(self):
        self.assertEqual(tidy("v600e", gene="BRAF").text, "V600E")
        self.assertEqual(tidy("v600e").text, "V600E")


class LeftAlone(unittest.TestCase):
    def assert_unchanged(self, text):
        reading = tidy(text)
        self.assertEqual(reading, Reading(text.strip()))
        self.assertFalse(reading.changed)

    def test_every_demo_input_is_left_alone(self):
        lines = [line.strip() for line in DEMO_LIST.read_text(encoding="utf-8").splitlines()]
        inputs = [line for line in lines if line and not line.startswith("#")]
        self.assertGreater(len(inputs), 15)
        for text in inputs:
            with self.subTest(text=text):
                self.assert_unchanged(text)

    def test_forms_that_already_read_are_left_alone(self):
        for text in (
            "BRAF V600E",
            "BRAF:V600E",
            "BRAF (p.V600E)",
            "BRAF P.V600E",
            "BRAF:p.V600E",
            "chr7-140453136-A-T",
            "7:140453136:A>T",
            "chr7:140453136 A>T",
            "rs113488022",
            "HLA-A V600E",
            "C4orf3 R10H",
            "TP53 R196*",
            "BRAF Val600=",
            "  BRAF V600E  ",
        ):
            with self.subTest(text=text):
                self.assert_unchanged(text)

    def test_what_cannot_be_fixed_is_not_pretended(self):
        """Missing information is for step 2 to ask about. Step 1 changes nothing it can't be sure of."""
        for text in ("BRAF V600", "nonsense words", "V600E", "chr7-140453136A-T", "BRAF c.1799"):
            with self.subTest(text=text):
                self.assertFalse(tidy(text).changed)

    def test_case_fixes_still_help_when_something_else_is_missing(self):
        reading = tidy("braf val600")
        self.assertEqual(reading.text, "braf Val600")
        self.assertTrue(reading.changed)


if __name__ == "__main__":
    unittest.main()
