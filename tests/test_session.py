"""The guided start: type variants, see answers, keep them, save them."""

import re
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from unittest import mock

from tests.helpers import ScriptedTerminal
from tests.test_integration import DEMO, WithData
from varshift import session as session_module
from varshift.errors import DataError
from varshift.session import DEMO_NOTICE, Session
from varshift.store import DataStore


class SessionTest(WithData):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.cwd = Path(folder.name)

    def run_session(self, *answers, pasted=(), opener=None, **options):
        self.term = ScriptedTerminal(answers=list(answers), pasted=pasted)
        self.session = Session(opener or (lambda: DataStore(self.data)), self.term, self.cwd, build="hg19", **options)
        self.code = self.session.run()
        return self.term.output

    def headings(self, text):
        return [line for line in text.splitlines() if "  ->  " in line]


class Queries(SessionTest):
    def test_the_header_names_the_data_and_the_kinds_of_input(self):
        out = self.run_session()
        self.assertIn("VarShift: offline variant annotation. Data: hg19, hg38", out)
        self.assertIn("ClinVar", out.splitlines()[0])
        self.assertIn("What do you have? Type one or more, separated by ;", out)
        self.assertIn("1  Protein change", out)
        self.assertEqual(self.code, 0)

    def test_one_query_shows_its_card_and_a_number_shows_it_again(self):
        out = self.run_session("BRAF V600E", "1", "list", "", "n")
        self.assertEqual(len(self.headings(out)), 2)
        self.assertEqual(out.count("Why this"), 2)
        self.assertRegex(out, r"\n\s+1  BRAF V600E\s+chr7-140453136-A-T")
        self.assertEqual(len(self.session.items), 1)

    def test_several_queries_get_one_line_each(self):
        out = self.run_session("BRAF V600E; KRAS G12D; TP53 R175H", "", "n")
        for number, text in enumerate(("BRAF V600E", "KRAS G12D", "TP53 R175H"), 1):
            self.assertRegex(out, rf"\n\s+{number}  {text}\s")
        self.assertIn("Type a number to see its full answer. 3 results in this session; save writes them to Excel.", out)
        self.assertEqual(self.headings(out), [])

    def test_the_same_query_twice_is_two_results(self):
        out = self.run_session("BRAF V600E; BRAF V600E", "list", "", "n")
        self.assertEqual(len(self.session.items), 2)
        self.assertEqual(len(re.findall(r"\n\s+2  BRAF V600E\s", out)), 2, "numbered 2 in the summary and in list")

    def test_pasted_lines_are_answered_together(self):
        out = self.run_session("BRAF V600E", "KRAS G12D", "", "n", pasted={1})
        self.assertRegex(out, r"\n\s+2  KRAS G12D\s")
        self.assertEqual(self.headings(out), [])
        self.assertEqual(len(self.session.items), 2)

    def test_a_slip_is_asked_about_and_the_choice_is_kept(self):
        out = self.run_session("BRAF V600", "1", "", "n")
        self.assertIn("Val600 changes to what?", out)
        read_as = next(line for line in out.splitlines() if line.lstrip().startswith("Read as"))
        self.assertIn("BRAF V600E", read_as)
        self.assertIn("chosen from the suggestions", read_as)
        self.assertEqual(self.session.items[0].answer.raw, "BRAF V600")

    def test_choosing_every_allele_adds_each_one(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        self.run_session("rs113488022", "a", "", "n")
        self.assertEqual(len(self.session.items), len(records))
        self.assertEqual([item.chosen[0]["variant"]["hg19"]["key"] for item in self.session.items], [r["variant"]["hg19"]["key"] for r in records])

    def test_skipping_a_question_adds_nothing(self):
        out = self.run_session("BRAF V600", "", "")
        self.assertIn("Skipped BRAF V600", out)
        self.assertEqual(self.session.items, [])

    def test_nothing_to_suggest_says_where_the_examples_are(self):
        out = self.run_session("xyz123", "")
        self.assertIn("could not annotate", out)
        self.assertIn("examples", out)
        self.assertEqual(self.session.items, [])

    def test_a_file_is_annotated_into_its_own_folder(self):
        out = self.run_session(str(DEMO / "variants.txt"), "")
        self.assertTrue((self.cwd / "variants_varshift" / "report.html").is_file())
        self.assertIn("results", out)
        self.assertEqual(self.session.items, [])


class Words(SessionTest):
    def test_a_number_before_any_result_shows_examples_of_that_kind(self):
        out = self.run_session("1", "")
        self.assertIn("BRAF p.Val600Glu", out)
        self.assertNotIn("NM_000546.6(TP53):c.524G>A", out)

    def test_help_shows_the_words_and_the_examples(self):
        out = self.run_session("help", "")
        for word in ("save", "list", "remove N", "clear", "mode", "examples"):
            self.assertIn(word, out)
        self.assertIn("NR_027676.2:n.200G>A", out)

    def test_remove_and_clear(self):
        out = self.run_session("BRAF V600E; KRAS G12D", "remove 1", "list", "clear", "list", "remove 4", "")
        self.assertIn("Removed 1: BRAF V600E", out)
        self.assertRegex(out, r"\n\s+1  KRAS G12D\s")
        self.assertIn("Cleared 1 result.", out)
        self.assertIn("No results yet.", out)
        self.assertIn("There is no result 4.", out)
        self.assertEqual(self.session.items, [])

    def test_the_mode_can_change_during_the_session(self):
        out = self.run_session("mode tumour", "KRAS G12D", "mode off", "", "n")
        self.assertIn("Mode: tumour", out)
        self.assertIn("OncoKB", out.split("Check it")[1])
        self.assertIn("Mode: none", out)
        self.assertEqual(self.session.items[0].chosen[0]["input"]["mode"], "tumour")
        self.assertIn("mode must be", self.run_session("mode sideways", ""))


class Saving(SessionTest):
    def test_save_writes_a_dated_workbook_and_never_overwrites(self):
        moment = datetime(2026, 11, 8, 14, 30)
        with mock.patch.object(session_module, "_now", return_value=moment):
            out = self.run_session("BRAF V600E", "save", "save", "")
        first = self.cwd / "varshift-session-2026-11-08-1430.xlsx"
        second = self.cwd / "varshift-session-2026-11-08-1430-2.xlsx"
        self.assertTrue(first.is_file())
        self.assertTrue(second.is_file())
        self.assertIn(str(first), out)
        self.assertIn("Links", out)
        with zipfile.ZipFile(first) as archive:
            self.assertIn("Links", archive.read("xl/workbook.xml").decode())
        self.assertTrue(all(item.saved for item in self.session.items))

    def test_save_to_a_named_file_or_a_folder(self):
        out = self.run_session("BRAF V600E", "save x.tsv", "save out/", "")
        header = (self.cwd / "x.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
        self.assertEqual(header[-7:], ["read_as", "answer", "best", "why", "also", "suggestion", "mode"])
        rows = (self.cwd / "x.tsv").read_text(encoding="utf-8").splitlines()[1:]
        self.assertEqual(len(rows), 1, "only the answer shown, not every candidate")
        for name in ("report.html", "results.xlsx", "results.tsv"):
            self.assertTrue((self.cwd / "out" / name).is_file(), name)
        self.assertIn("Saved 1 result", out)

    def test_saving_over_a_file_asks_first(self):
        (self.cwd / "keep.tsv").write_text("keep me", encoding="utf-8")
        out = self.run_session("BRAF V600E", "save keep.tsv", "n", "", "n")
        self.assertEqual((self.cwd / "keep.tsv").read_text(encoding="utf-8"), "keep me")
        self.assertIn("Replace it?", out)
        self.run_session("BRAF V600E", "save keep.tsv", "y", "")
        self.assertNotEqual((self.cwd / "keep.tsv").read_text(encoding="utf-8"), "keep me")

    def test_a_format_it_cant_write_keeps_the_results(self):
        out = self.run_session("BRAF V600E", "save x.pdf", "", "n")
        self.assertIn(".xlsx", out)
        self.assertFalse(self.session.items[0].saved)

    def test_a_failed_save_keeps_the_results(self):
        with mock.patch.object(session_module, "open_file_sink", side_effect=OSError("disk full")):
            out = self.run_session("BRAF V600E", "save x.tsv", "", "n")
        self.assertIn("disk full", out)
        self.assertFalse(self.session.items[0].saved)

    def test_nothing_to_save(self):
        self.assertIn("Nothing to save yet.", self.run_session("save", ""))


class Quitting(SessionTest):
    def test_quitting_with_unsaved_results_offers_to_save(self):
        with mock.patch.object(session_module, "_now", return_value=datetime(2026, 11, 8, 9, 5)):
            out = self.run_session("BRAF V600E", "", "y")
        self.assertIn("Save the 1 result to Excel first?", out)
        self.assertTrue((self.cwd / "varshift-session-2026-11-08-0905.xlsx").is_file())
        self.assertEqual(self.code, 0)

    def test_answering_no_saves_nothing(self):
        self.run_session("BRAF V600E", "", "n")
        self.assertEqual(list(self.cwd.glob("*.xlsx")), [])
        self.assertEqual(self.code, 0)

    def test_end_of_input_quits_cleanly(self):
        self.run_session("BRAF V600E")
        self.assertEqual(self.code, 0)

    def test_ctrl_c_at_the_prompt_stops_at_once(self):
        out = self.run_session("BRAF V600E", KeyboardInterrupt())
        self.assertEqual(self.code, 130)
        self.assertIn("Stopped. 1 result was not saved.", out)

    def test_ctrl_c_during_a_query_stops_only_that_query(self):
        out = self.run_session("BRAF V600", KeyboardInterrupt(), "KRAS G12D", "", "n")
        self.assertIn("Stopped that query.", out)
        self.assertEqual([item.answer.raw for item in self.session.items], ["KRAS G12D"])
        self.assertEqual(self.code, 0)


class NoData(SessionTest):
    def test_the_demo_data_can_be_tried_and_is_removed_after(self):
        def no_data():
            raise DataError("no VarShift data here")

        with mock.patch.object(session_module, "_now", return_value=datetime(2026, 11, 8, 10, 0)):
            out = self.run_session("1", "BRAF V600E", "save", "", opener=no_data)
        self.assertIn("Try the 12-gene demo data now", out)
        self.assertIn("Data: 12-gene demo data", out)
        self.assertEqual(len(self.headings(out)), 1)
        self.assertEqual(list(self.cwd.glob("varshift-demo-session-*")), [])
        with zipfile.ZipFile(self.cwd / "varshift-session-2026-11-08-1000.xlsx") as archive:
            texts = "".join(archive.read(name).decode() for name in archive.namelist() if name.endswith(".xml"))
        self.assertIn(DEMO_NOTICE, texts)

    def test_setting_up_runs_setup_then_opens_the_data(self):
        opened = []

        def opener():
            if not opened:
                opened.append("tried")
                raise DataError("no VarShift data here")
            return DataStore(self.data)

        setup = mock.Mock(return_value=0)
        out = self.run_session("2", "BRAF V600E", "", "n", opener=opener, setup=setup)
        setup.assert_called_once()
        self.assertEqual(len(self.headings(out)), 1)

    def test_choosing_nothing_ends_the_session(self):
        def no_data():
            raise DataError("no VarShift data here")

        out = self.run_session("", opener=no_data)
        self.assertIn("varshift setup", out)
        self.assertEqual(self.code, 1)
