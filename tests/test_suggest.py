"""Step 2 of the input assistant: an error becomes a question with options built from the installed data."""

import time
import unittest
from unittest import mock

from tests.test_integration import WithData
from varshift.errors import ReferenceMismatch, ResolveError, VarShiftError
from varshift.store import DataStore
from varshift.suggest import Question, close_symbols, edits, suggest, typo_question


class Suggestions(WithData):
    def question(self, text, build="hg19", fast_only=False, gene=""):
        try:
            self.annotator.annotate_text(text, build, gene)
        except VarShiftError as exc:
            return suggest(self.store, text, exc, build, fast_only=fast_only, gene=gene)
        self.fail(f"{text} annotated without an error")

    def texts(self, question):
        self.assertIsInstance(question, Question)
        return [option.text for option in question.options]

    def test_a_change_with_no_alt_lists_the_changes_known_at_that_residue(self):
        question = self.question("BRAF V600")
        self.assertTrue(question.text.startswith("Val600 changes to what?"), question.text)
        self.assertIn("BRAF V600E", self.texts(question))
        option = next(o for o in question.options if o.text == "BRAF V600E")
        self.assertEqual(option.label, "V600E")
        self.assertIn("ClinVar:", option.detail)
        self.assertIn("CIViC:", option.detail)

    def test_every_option_for_a_residue_annotates(self):
        for option in self.question("BRAF V600").options:
            with self.subTest(option=option.text):
                self.assertTrue(self.annotator.annotate_text(option.text, "hg19"))

    def test_the_tidied_text_is_what_is_suggested_from(self):
        self.assertIn("BRAF V600E", self.texts(self.question("braf v600")))
        self.assertIn("BRAF V600E", self.texts(self.question("BRAF p.Val600")))

    def test_a_change_with_no_gene_lists_the_genes_that_have_it(self):
        question = self.question("V600E")
        self.assertTrue(question.text.startswith("Which gene?"), question.text)
        self.assertEqual(self.texts(question)[0], "BRAF V600E")

    def test_the_slow_search_is_skipped_in_files(self):
        self.assertIsNone(self.question("V600E", fast_only=True))

    def test_an_unknown_gene_is_corrected_by_spelling(self):
        question = self.question("BRFA V600E")
        self.assertEqual(question.text, "BRFA isn't a gene symbol. Did you mean:")
        self.assertEqual(self.texts(question)[0], "BRAF V600E")
        self.assertLessEqual(len(question.options), 4)
        self.assertEqual(self.texts(self.question("brfa V600E"))[0], "BRAF V600E")

    def test_names_are_looked_up_with_the_same_rules_as_gene_symbol(self):
        names = ["BRAF", "braf1", "KRAS2", "K-RAS4B", "P53", "ERBB", "NOTAGENE"]
        expected = {name.upper(): self.store.gene_symbol(name) for name in names if self.store.gene_symbol(name)}
        self.assertEqual(self.store.gene_symbols_for(names), expected)
        self.assertEqual(expected["KRAS2"], "KRAS")
        approved = self.store.approved_symbols()
        self.assertEqual(approved["BRAF"], "BRAF")
        self.assertNotIn("KRAS2", approved)

    def test_the_whole_symbol_list_is_read_only_when_one_edit_finds_nothing(self):
        for text, reads_all in (("BRFA V600E", False), ("NOTAGENE V600E", True)):
            with self.subTest(text=text), mock.patch.object(self.store, "approved_symbols", wraps=self.store.approved_symbols) as full:
                self.question(text)
                self.assertEqual(full.called, reads_all)

    def test_a_gene_given_with_gene_is_not_rewritten(self):
        self.assertIsNone(self.question("V600E", gene="BRFA"))

    def test_a_wrong_residue_is_corrected(self):
        # Not A600E: another BRAF transcript has Ala at 600, so that one annotates, with other_transcript_used.
        question = self.question("BRAF W600E")
        self.assertTrue(question.text.startswith("Residue 600 is Val, not Trp."), question.text)
        self.assertEqual(self.texts(question), ["BRAF V600E"])
        self.assertEqual(self.texts(self.question("BRAF p.Trp600Glu")), ["BRAF p.Val600Glu"])

    def test_swapped_ref_and_alt_in_a_c_name(self):
        question = self.question("BRAF c.1799A>T")
        self.assertEqual(self.texts(question), ["BRAF c.1799T>A"])
        self.assertIn("swapped", question.text)

    def test_a_wrong_ref_in_a_c_name(self):
        # Not c.1799C>G: another BRAF transcript has C at c.1799, so that one annotates.
        question = self.question("BRAF c.1799A>C")
        self.assertEqual(self.texts(question), ["BRAF c.1799T>C"])
        self.assertIn("corrected", question.text)

    def test_a_position_on_the_other_build(self):
        question = self.question("chr7-140753336-A-T", "hg19")
        self.assertEqual(len(question.options), 1)
        self.assertEqual(question.options[0].text, "chr7-140753336-A-T")
        self.assertEqual(question.options[0].build, "hg38")
        self.assertIn("build=hg38", question.text)
        self.assertTrue(self.annotator.annotate_text("chr7-140753336-A-T", question.options[0].build))

    def test_the_other_build_not_installed_gives_a_hint_only(self):
        with mock.patch.object(DataStore, "builds", new_callable=mock.PropertyMock, return_value=["hg19"]):
            question = suggest(self.store, "chr7-140753336-A-T", ReferenceMismatch("REF A does not match"), "hg19")
        self.assertEqual(question.options, ())
        self.assertIn("varshift setup build=hg38", question.text)

    def test_nothing_to_suggest(self):
        for text in ("xyz123", "BRAF c.1799_1800delinsAA garbage"):
            with self.subTest(text=text):
                self.assertIsNone(self.question(text))

    def test_the_resolver_keeps_what_went_wrong(self):
        with self.assertRaises(ResolveError) as caught:
            self.annotator.annotate_text("BRAF W600E", "hg19")
        cause = caught.exception.causes[0]
        self.assertIsInstance(cause, ReferenceMismatch)
        self.assertEqual((cause.expected, cause.found, cause.position), ("W", "V", 600))
        self.assertEqual(len(caught.exception.causes), len(caught.exception.tried))
        self.assertTrue(cause.transcript.startswith("NM_004333"))
        with self.assertRaises(ResolveError) as caught:
            self.annotator.annotate_text("BRAF c.1799A>T", "hg19")
        cause = caught.exception.causes[0]
        self.assertEqual((cause.expected, cause.found, cause.position), ("A", "T", "1799"))

    def test_each_suggestion_is_quick(self):
        for text in ("BRAF V600", "V600E", "BRFA V600E", "NOTAGENE V600E", "BRAF W600E", "BRAF c.1799A>T"):
            with self.subTest(text=text):
                try:
                    self.annotator.annotate_text(text, "hg19")
                except VarShiftError as exc:
                    started = time.perf_counter()
                    suggest(self.store, text, exc, "hg19")
                    self.assertLess(time.perf_counter() - started, 0.3)


class Spelling(unittest.TestCase):
    SYMBOLS = {"BRAF": "BRAF", "RBFA": "RBFA", "BRF1": "BRF1", "BRF2": "BRF2", "BRAP": "BRAP"}
    COUNTS = {"BRAF": 900, "RBFA": 3, "BRF1": 40, "BRF2": 20, "BRAP": 10}

    def test_a_swapped_pair_counts_as_one_edit(self):
        self.assertIn("BRAF", edits("BRFA"))

    @staticmethod
    def names(found):
        return [spelling.symbol for spelling in found]

    def test_the_gene_with_most_clinvar_records_comes_first(self):
        found = close_symbols("BRFA", self.SYMBOLS, self.COUNTS.get)
        self.assertEqual(found[0].symbol, "BRAF")
        self.assertLessEqual(len(found), 4)

    def test_aliases_lead_to_the_approved_symbol_once(self):
        symbols = {"KRAS": "KRAS", "KRAS2": "KRAS", "KRAS1": "KRAS"}
        self.assertEqual(self.names(close_symbols("KRAS3", symbols, lambda symbol: 1)), ["KRAS"])

    def test_fewer_edits_then_a_real_symbol_come_before_more_records(self):
        symbols = {"EGFR": "EGFR", "EGR1": "EGR1", "PTK3A": "DDR1", "PIK3CA": "PIK3CA", "PIK3CD": "PIK3CD"}
        counts = {"EGFR": 10, "EGR1": 5000, "DDR1": 9000, "PIK3CA": 1, "PIK3CD": 2}
        self.assertEqual(self.names(close_symbols("PIK3AC", symbols, counts.get)), ["PIK3CA", "PIK3CD", "DDR1"])
        found = close_symbols("PIK3AC", symbols, counts.get)
        self.assertEqual([(f.edits, f.direct) for f in found], [(1, True), (2, True), (2, False)])
        self.assertEqual(found[2].matched, "PTK3A")

    def test_short_symbols_get_one_edit_and_longer_ones_two(self):
        self.assertEqual(close_symbols("BXXF", {"BRAF": "BRAF"}, lambda s: 1), [])
        self.assertEqual(self.names(close_symbols("PIK3XX", {"PIK3CA": "PIK3CA"}, lambda s: 1)), ["PIK3CA"])
        self.assertEqual(close_symbols("ABCDEFGHIXX", {"ABCDEFGHIJK": "ABCDEFGHIJK"}, lambda s: 1), [])

    def test_the_word_itself_is_not_offered(self):
        self.assertEqual(close_symbols("BRAF", {"BRAF": "BRAF"}, lambda s: 1), [])


if __name__ == "__main__":
    unittest.main()


class TypoOnTheDefaultTranscript(WithData):
    """A change that fails on the default transcript but fits another is often a typo, so it says so."""

    def test_a_residue_typo_that_fits_another_transcript_is_noted(self):
        records = self.annotator.annotate_text("BRAF A600E", "hg19")
        source = records[0]["input"]
        default = source["default_transcript"]
        self.assertTrue(default["accession"].startswith("NM_004333"))
        self.assertIn("has Val at residue 600, not Ala", default["says"])
        self.assertEqual(default["suggestion"], "BRAF V600E")
        for record in records:
            self.assertIn("default_transcript_differs", record["flags"])
            self.assertEqual(record["suggestions"], ["BRAF V600E"])
            self.assertNotEqual(record["transcript"]["accession"], default["accession"])

    def test_a_base_typo_that_fits_another_transcript_is_noted(self):
        source = self.annotator.annotate_text("BRAF c.1799C>G", "hg19")[0]["input"]
        self.assertIn("has T at c.1799, not C", source["default_transcript"]["says"])
        self.assertEqual(source["default_transcript"]["suggestion"], "BRAF c.1799T>G")

    def test_nothing_is_noted_when_the_default_fits_or_a_transcript_is_named(self):
        # The c. name of what BRAF A600E found, with its transcript named: that transcript was meant.
        named_c = self.annotator.annotate_text("BRAF A600E", "hg19")[0]["hgvs"]["c"]
        for text in ("BRAF V600E", "NM_004333.6:c.1799T>A", "KRAS G12D", named_c):
            with self.subTest(text=text):
                record = self.annotator.annotate_text(text, "hg19")[0]
                self.assertNotIn("default_transcript", record["input"])
                self.assertNotIn("default_transcript_differs", record["flags"])

    def test_the_question_offers_the_default_transcript_first_and_keeping_it(self):
        records = self.annotator.annotate_text("BRAF A600E", "hg19")
        question = typo_question(records)
        self.assertIn("has Val at residue 600, not Ala", question.text)
        self.assertEqual([o.text for o in question.options], ["BRAF V600E", "BRAF A600E"])
        self.assertEqual([o.keep for o in question.options], [False, True])
        self.assertIsNone(typo_question(self.annotator.annotate_text("BRAF V600E", "hg19")))

    def test_the_card_says_it_in_words(self):
        from varshift.answer import decide, render_answer
        from varshift.flags import FLAG_HELP
        from varshift.term import QuietTerminal

        # On the demo data A600E fits two DNA changes, a choice, so one candidate's card is drawn.
        text = render_answer(decide(self.annotator.annotate_text("BRAF A600E", "hg19")[:1]), QuietTerminal())
        self.assertIn("If you meant the default transcript: BRAF V600E", text)
        self.assertNotIn(FLAG_HELP["other_transcript_used"][0], text)
        self.assertNotIn("default_transcript_differs", text)
