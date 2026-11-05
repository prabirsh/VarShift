"""The whole pipeline on demo/sources: real sequence, real transcripts, real ClinVar records."""

import contextlib
import io
import json
import os
import stat
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

from tests.helpers import ScriptedTerminal
from varshift import cli
from varshift.annotate import Annotator
from varshift.batch import annotate_lines
from varshift.cli import main
from varshift.config import is_windows_mount
from varshift.errors import DataError, InputError, ReferenceMismatch, ResolveError
from varshift.pipeline.demo import build_demo_store
from varshift.pipeline.selftest import run_selftest
from varshift.store import DataStore

DEMO = Path(__file__).resolve().parents[1] / "demo"


class WithData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        cls.data = Path(cls.folder.name) / "data"
        build_demo_store(DEMO, cls.data)
        cls.store = DataStore(cls.data)
        cls.annotator = Annotator(cls.store)

    @classmethod
    def tearDownClass(cls):
        cls.store.close()
        cls.folder.cleanup()

    def one(self, text, build=None):
        return self.annotator.annotate_text(text, build)[0]


class KnownVariant(WithData):
    def test_the_self_test_passes(self):
        results = run_selftest(self.store)
        self.assertGreater(len(results), 0)
        self.assertEqual([r for r in results if not r["passed"]], [])

    def test_braf_v600e_from_protein(self):
        record = self.one("BRAF V600E", "hg19")
        self.assertEqual(record["variant"]["hg19"]["key"], "chr7-140453136-A-T")
        self.assertEqual(record["variant"]["hg38"]["key"], "chr7-140753336-A-T")
        self.assertEqual(record["gene"]["symbol"], "BRAF")
        self.assertEqual(record["gene"]["hgnc"], "HGNC:1097")
        self.assertEqual(record["gene"]["ncbi_gene"], "673")
        self.assertEqual(record["gene"]["ensembl"], "ENSG00000157764")
        self.assertTrue(record["hgvs"]["c"].endswith(":c.1799T>A"))
        self.assertTrue(record["hgvs"]["p"].endswith(":p.Val600Glu"))
        self.assertEqual(record["effect"]["type"], "missense_variant")
        self.assertEqual(record["effect"]["aa_position"], 600)
        self.assertEqual(record["protein"]["uniprot"], "P15056")
        self.assertIn("PF07714", [p["accession"] for p in record["protein"]["pfam"]])
        self.assertIn("PF07714", [d["accession"] for d in record["protein"]["domains_at_residue"]])

    def test_clinvar_is_quoted_as_it_stands(self):
        clinvar = self.one("BRAF V600E", "hg19")["evidence"]["clinvar"]
        self.assertEqual(clinvar["id"], 13961)
        self.assertEqual(clinvar["classification"], "Conflicting classifications of pathogenicity")
        self.assertIn("names", clinvar)

    def test_one_protein_change_gives_every_dna_candidate(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19")
        keys = [r["variant"]["hg19"]["key"] for r in records]
        self.assertEqual(keys[0], "chr7-140453136-A-T")
        self.assertIn("chr7-140453135-CA-TT", keys)
        self.assertEqual(records[0]["candidate"]["of"], len(records))
        self.assertTrue(records[0]["candidate"]["in_clinvar"])
        self.assertIn("multiple_candidates", records[0]["flags"])

    def test_in_clinvar_is_true_for_a_position_in_clinvar(self):
        for text in ("chr7-140453136-A-T", "BRAF c.1799T>A"):
            with self.subTest(text=text):
                self.assertTrue(self.one(text, "hg19")["candidate"]["in_clinvar"])

    def test_results_carry_the_answer_fields(self):
        record = self.one("NM_004333.6:c.1799T>A", "hg19")
        self.assertEqual(record["input"]["answer"], "single")
        self.assertEqual(record["input"]["read_as"], "")
        self.assertIn("mode", record["input"])
        self.assertEqual(record["candidate"]["bases_changed"], 1)
        self.assertIn("best", record["candidate"])
        self.assertIn("why", record["candidate"])
        changed = {r["variant"]["hg19"]["key"]: r["candidate"]["bases_changed"] for r in self.annotator.annotate_text("BRAF V600E", "hg19")}
        self.assertEqual(changed["chr7-140453136-A-T"], 1)
        self.assertEqual(changed["chr7-140453135-CA-TT"], 2)

    def test_every_input_form_reaches_the_same_variant(self):
        for text in (
            "BRAF V600E",
            "BRAF p.Val600Glu",
            "BRAF c.1799T>A",
            "NC_000007.13:g.140453136A>T",
            "chr7:g.140453136A>T",
            "chr7-140453136-A-T",
            "7:140453136:A>T",
        ):
            with self.subTest(text=text):
                self.assertEqual(self.one(text, "hg19")["variant"]["hg19"]["key"], "chr7-140453136-A-T")

    def test_an_rs_number_lists_each_allele_clinvar_has(self):
        keys = [r["variant"]["hg19"]["key"] for r in self.annotator.annotate_text("rs113488022", "hg19")]
        self.assertIn("chr7-140453136-A-T", keys)
        self.assertGreater(len(keys), 1)

    def test_same_answer_from_either_build(self):
        from_hg19 = self.one("chr7-140453136-A-T", "hg19")
        from_hg38 = self.one("chr7-140753336-A-T", "hg38")
        self.assertEqual(from_hg19["variant"], {**from_hg19["variant"], "hg38": from_hg38["variant"]["hg38"]})
        self.assertEqual(from_hg19["evidence"]["clinvar"]["id"], from_hg38["evidence"]["clinvar"]["id"])

    def test_refseq_accession_sets_the_build(self):
        self.assertEqual(self.one("NC_000007.14:g.140753336A>T")["input"]["build"], "hg38")

    def test_all_requested_links_are_present(self):
        links = {link["name"]: link for link in self.one("BRAF V600E", "hg19")["links"]}
        for name in ("ClinVar", "NCBI Gene", "NCBI RefSeq", "HGNC", "Ensembl", "UniProt", "Pfam", "ProteinPaint", "Franklin", "COSMIC", "VarSome", "OncoKB", "dbSNP"):
            self.assertTrue(links[name]["url"], name)
        self.assertEqual(links["ClinVar"]["url"], "https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/")
        self.assertEqual(links["ClinVar"]["label"], "confirmed_id")
        self.assertEqual(links["Franklin"]["label"], "built")

    def test_every_result_records_data_versions(self):
        self.assertIn("clinvar.hg19", self.one("BRAF V600E", "hg19")["versions"])

    def test_oncology_evidence(self):
        evidence = self.one("BRAF V600E", "hg19")["evidence"]
        self.assertTrue(evidence["civic"], "CIViC has BRAF V600E")
        self.assertGreater(evidence["civic"][0]["evidence_items"], 0)


class Collect:
    """A writer that keeps what it is given."""

    def __init__(self):
        self.records = []

    def write(self, record, prefix=None):
        self.records.append(record)


class TidiedInput(WithData):
    def test_a_slip_is_fixed_and_the_reading_is_recorded(self):
        for typed, read_as in (
            ("braf val600glu", "BRAF Val600Glu"),
            ("BRAF c.1799T-A", "BRAF c.1799T>A"),
            ("NM_004333.6:1799T>A", "NM_004333.6:c.1799T>A"),
            ("chr7:140453136A>T", "chr7:g.140453136A>T"),
        ):
            with self.subTest(typed=typed):
                record = self.one(typed, "hg19")
                self.assertEqual(record["variant"]["hg19"]["key"], "chr7-140453136-A-T")
                self.assertEqual(record["input"]["raw"], typed)
                self.assertEqual(record["input"]["read_as"], read_as)
                self.assertTrue(record["input"]["read_as_notes"])

    def test_text_that_reads_is_left_alone(self):
        for text in ("BRAF V600E", "NM_004333.6:c.1799T>A", "chr7-140453136-A-T", "rs113488022"):
            with self.subTest(text=text):
                source = self.one(text, "hg19")["input"]
                self.assertEqual(source["read_as"], "")
                self.assertEqual(source["read_as_notes"], [])

    def test_the_gene_symbol_is_written_as_hgnc_writes_it(self):
        source = self.one("braf V600E", "hg19")["input"]
        self.assertEqual(source["raw"], "braf V600E")
        self.assertEqual(source["read_as"], "BRAF V600E")
        self.assertEqual(source["read_as_notes"], ["gene symbol written as HGNC writes it"])

    def test_an_alias_is_named(self):
        for typed, read_as in (("KRAS2 G12D", "KRAS G12D (KRAS2 is an alias)"), ("kras2 g12d", "KRAS G12D (kras2 is an alias)")):
            with self.subTest(typed=typed):
                record = self.one(typed, "hg19")
                self.assertEqual(record["gene"]["symbol"], "KRAS")
                self.assertEqual(record["input"]["read_as"], read_as)

    def test_an_alias_given_with_gene_leaves_the_text_alone(self):
        source = self.annotator.annotate_text("G12D", "hg19", gene="KRAS2")[0]["input"]
        self.assertEqual(source["read_as"], "")

    def test_one_query_does_not_carry_the_last_gene_into_the_next(self):
        self.one("KRAS2 G12D", "hg19")
        self.assertEqual(self.one("chr7-140453136-A-T", "hg19")["input"]["read_as"], "")

    def test_the_mode_is_recorded(self):
        self.assertEqual(self.annotator.annotate_text("BRAF V600E", "hg19", mode="tumour")[0]["input"]["mode"], "tumour")
        self.assertIsNone(self.one("BRAF V600E", "hg19")["input"]["mode"])

    def test_a_list_keeps_the_reading_of_a_line_that_still_fails(self):
        writer = Collect()
        counts = annotate_lines(self.annotator, [(1, "BRAF(V600)"), (2, "braf v600e")], writer, "hg19", mode="inherited")
        self.assertEqual((counts["lines"], counts["errors"]), (2, 1))
        failed, first = writer.records[0], writer.records[1]
        self.assertTrue(failed["error"])
        self.assertEqual(failed["input"]["raw"], "BRAF(V600)")
        self.assertEqual(failed["input"]["read_as"], "BRAF V600")
        self.assertEqual(failed["input"]["mode"], "inherited")
        self.assertEqual(first["input"]["read_as"], "BRAF V600E")
        self.assertEqual(first["input"]["mode"], "inherited")


class NoTyping:
    """Standard input that says it is a terminal but must never be read."""

    def isatty(self):
        return True

    def fileno(self):
        raise OSError("not a real stream")

    def read(self, *args):
        raise AssertionError("standard input was read")

    readline = read


class ScreenAnswers(WithData):
    """What varshift annotate prints for one description, or a ; list, with no out= or format=."""

    def _cli(self, *args, terminal=None):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.ExitStack() as stack:
            stack.enter_context(contextlib.redirect_stdout(stdout))
            stack.enter_context(contextlib.redirect_stderr(stderr))
            stack.enter_context(mock.patch("sys.stdin", NoTyping()))
            if terminal is not None:
                stack.enter_context(mock.patch.object(cli, "_terminal", return_value=terminal))
            code = main(list(args) + [f"data={self.data}"])
        return code, stdout.getvalue(), stderr.getvalue()

    def test_one_clear_answer_with_why(self):
        code, out, err = self._cli("annotate", "BRAF V600E", "build=hg19")
        self.assertEqual(code, 0, err)
        self.assertIn("Why this", out)
        self.assertIn("chr7-140453136-A-T", out.splitlines()[0])
        self.assertLess(len(out.splitlines()), 25)
        self.assertNotIn("candidate 2 of 2", out)
        self.assertNotIn("\x1b", out)

    def test_show_all_prints_every_candidate_in_full(self):
        code, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19", "show=all")
        self.assertEqual(code, 0)
        self.assertIn("candidate 2 of 2", out)

    def test_a_choice_is_listed_without_waiting_when_no_one_can_answer(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        code, out, err = self._cli("annotate", "rs113488022", "build=hg19")
        self.assertEqual(code, 0, err)
        self.assertIn(f"covers {len(records)} different variants", out)
        for record in records:
            self.assertIn(record["variant"]["hg19"]["key"], out)
        self.assertIn("Each is a different variant. To see all in full: show=all", out)

    def test_a_choice_made_at_the_terminal_shows_that_variant(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        terminal = ScriptedTerminal(answers=["2"])
        code, _, err = self._cli("annotate", "rs113488022", "build=hg19", terminal=terminal)
        self.assertEqual(code, 0, err)
        cards = [line for line in terminal.output.splitlines() if line.startswith("rs113488022  ->")]
        self.assertEqual(len(cards), 1, terminal.output)
        self.assertIn(records[1]["variant"]["hg19"]["key"], cards[0])
        self.assertIn("Choose 1 to 2, a for all, or Enter for the first", terminal.output)

    def test_all_can_be_chosen(self):
        terminal = ScriptedTerminal(answers=["a"])
        self.assertEqual(self._cli("annotate", "rs113488022", "build=hg19", terminal=terminal)[0], 0)
        self.assertEqual(len([line for line in terminal.output.splitlines() if line.startswith("rs113488022  ->")]), 2)

    def test_no_answer_at_the_end_of_input_chooses_nothing(self):
        terminal = ScriptedTerminal(answers=[])
        code, _, err = self._cli("annotate", "rs113488022", "build=hg19", terminal=terminal)
        self.assertEqual(code, 1)
        self.assertIn("Nothing chosen", err)

    def test_the_mode_is_checked_and_recorded(self):
        code, _, err = self._cli("annotate", "BRAF V600E", "build=hg19", "mode=bogus")
        self.assertEqual(code, 1)
        for word in ("tumour", "inherited", "none"):
            self.assertIn(word, err)
        for given, recorded in (("tumour", "tumour"), ("tumor", "tumour"), ("inherited", "inherited"), ("none", None)):
            with self.subTest(given=given):
                code, out, err = self._cli("annotate", "BRAF V600E", "build=hg19", f"mode={given}", "format=json")
                self.assertEqual(code, 0, err)
                self.assertEqual(json.loads(out)[0]["input"]["mode"], recorded)

    def test_the_mode_chooses_the_links_on_the_card(self):
        _, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19", "mode=tumour")
        self.assertIn("OncoKB", out.split("Check it")[1])

    def test_colour_is_checked_and_only_used_when_asked(self):
        code, _, err = self._cli("annotate", "BRAF V600E", "build=hg19", "color=bogus")
        self.assertEqual(code, 1)
        for word in ("auto", "yes", "no"):
            self.assertIn(word, err)
        with mock.patch.dict(os.environ):
            os.environ.pop("NO_COLOR", None)
            _, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19", "color=yes")
        self.assertIn("\x1b[1;31mHIGH", out)
        with mock.patch.dict(os.environ, {"NO_COLOR": "1"}):
            self.assertNotIn("\x1b", self._cli("annotate", "BRAF V600E", "build=hg19", "color=yes")[1])

    def test_links_all_and_show_take_only_all(self):
        _, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19", "links=all")
        self.assertIn("NCBI RefSeq", out)
        self.assertNotIn("links=all", out)
        for word in ("show=some", "links=few"):
            with self.subTest(word=word):
                code, _, err = self._cli("annotate", "BRAF V600E", "build=hg19", word)
                self.assertEqual(code, 1)
                self.assertIn("all", err)

    def test_several_queries_get_one_line_each(self):
        code, out, err = self._cli("annotate", "BRAF V600E; KRAS G12D; NOTAGENE V600E", "build=hg19")
        self.assertEqual(code, 0, err)
        lines = out.splitlines()
        self.assertTrue(lines[0].lstrip().startswith("1  BRAF V600E"), lines[0])
        self.assertIn("chr7-140453136-A-T", lines[0])
        self.assertTrue(lines[1].lstrip().startswith("2  KRAS G12D"), lines[1])
        self.assertTrue(lines[2].lstrip().startswith("3  NOTAGENE V600E"), lines[2])
        self.assertIn("could not annotate", lines[2])
        self.assertEqual(lines[0].index("chr7"), lines[1].index("chr12"))
        self.assertIn('varshift annotate "BRAF V600E" build=hg19', out)
        self.assertIn("out=results.xlsx", out)

    def test_a_list_file_and_a_vcf_record_the_mode(self):
        listing = Path(self.folder.name) / "modes.txt"
        listing.write_text("BRAF V600E\nNOTAGENE V600E\n", encoding="utf-8")
        target = Path(self.folder.name) / "modes.jsonl"
        code, _, err = self._cli("annotate", str(listing), "build=hg19", "mode=inherited", f"out={target}", "replace")
        self.assertEqual(code, 0, err)
        rows = [json.loads(line) for line in target.read_text(encoding="utf-8").splitlines()]
        self.assertEqual({row["input"]["mode"] for row in rows}, {"inherited"})
        tsv = Path(self.folder.name) / "modes.tsv"
        code, _, err = self._cli("annotate", str(DEMO / "demo.hg19.vcf"), "mode=tumour", f"out={tsv}", "replace")
        self.assertEqual(code, 0, err)
        lines = tsv.read_text(encoding="utf-8").splitlines()
        column = lines[0].split("\t").index("mode")
        self.assertEqual({line.split("\t")[column] for line in lines[1:]}, {"tumour"})


class OtherClasses(WithData):
    def test_frameshift_through_clinvar_names(self):
        record = self.one("BRCA1 c.68_69del", "hg19")
        self.assertEqual(record["variant"]["hg19"]["key"], "chr17-41276044-ACT-A")
        self.assertRegex(record["hgvs"]["p"], r":p\.Glu23ValfsTer17$")
        self.assertEqual(record["effect"]["type"], "frameshift_variant")

    def test_minus_strand_gene(self):
        record = self.one("TP53 R175H", "hg38")
        self.assertEqual(record["transcript"]["strand"], "-")
        self.assertEqual(record["variant"]["hg38"]["key"], "chr17-7675088-C-T")
        self.assertEqual(record["variant"]["hg19"]["key"], "chr17-7578406-C-T")

    def test_inframe_deletion(self):
        records = self.annotator.annotate_text("EGFR E746_A750del", "hg19")
        self.assertEqual(records[0]["effect"]["type"], "inframe_deletion")
        self.assertRegex(records[0]["hgvs"]["p"], r":p\.Glu746_Ala750del$")

    def test_alias_finds_the_gene(self):
        self.assertEqual(self.one("ERBB1 L858R", "hg19")["gene"]["symbol"], "EGFR")

    def test_nonsense_candidates(self):
        records = self.annotator.annotate_text("TP53 R196*", "hg19")
        self.assertEqual(records[0]["effect"]["type"], "stop_gained")

    def test_variant_not_in_clinvar_is_not_called_benign(self):
        ref = self.store.genome("hg19").fetch("chr7", 140499999, 140500000)
        alt = "A" if ref != "A" else "C"
        record = self.one(f"chr7-140500000-{ref}-{alt}", "hg19")
        clinvar = record["evidence"]["clinvar"]
        self.assertEqual(clinvar["status"], "not_found")
        self.assertTrue(clinvar["version"])
        self.assertNotIn("classification", clinvar)


class Rejections(WithData):
    def test_wrong_reference_base_at_a_position(self):
        with self.assertRaises(ReferenceMismatch) as caught:
            self.one("chr7-140453136-G-T", "hg19")
        self.assertIn("which has A", str(caught.exception))

    def test_wrong_reference_base_in_c_name(self):
        with self.assertRaises(ResolveError) as caught:
            self.one("NM_004333.6:c.1799G>A", "hg19")
        self.assertTrue(caught.exception.tried)

    def test_wrong_residue_in_p_name(self):
        with self.assertRaises(ResolveError) as caught:
            self.one("BRAF W600E", "hg19")
        self.assertIn("Val", " ".join(caught.exception.tried))

    def test_another_isoform_is_tried_and_flagged(self):
        record = self.one("BRAF A600E", "hg19")
        self.assertIn("other_transcript_used", record["flags"])

    def test_build_is_required_for_a_bare_position(self):
        with self.assertRaises(InputError):
            self.one("chr7-140453136-A-T")

    def test_build_that_contradicts_the_accession(self):
        with self.assertRaises(InputError):
            self.one("NC_000007.14:g.140753336A>T", "hg19")

    def test_unknown_gene(self):
        with self.assertRaises(ResolveError):
            self.one("NOTAGENE V600E", "hg19")

    def test_unknown_transcript(self):
        with self.assertRaises(ResolveError):
            self.one("NM_999999999.1:c.1A>T", "hg19")

    def test_frameshift_not_in_clinvar_explains_why(self):
        with self.assertRaises(ResolveError) as caught:
            self.one("BRAF V600fs", "hg19")
        self.assertIn("does not define one DNA change", " ".join(caught.exception.tried))

    def test_sql_text_in_a_gene_name_is_harmless(self):
        with self.assertRaises((InputError, ResolveError)):
            self.one("x' OR '1'='1 V600E", "hg19")
        self.assertGreater(self.store.clinvar_row("hg19", "chr7", 140453136, "A", "T")["variation_id"], 0)


class Mitochondria(WithData):
    def test_hg19_chrm_is_not_matched_to_clinvar(self):
        genome = self.store.genome("hg19")
        ref = genome.fetch("chrM", 100, 101)
        alt = "A" if ref != "A" else "C"
        record = self.one(f"chrM-101-{ref}-{alt}", "hg19")
        self.assertIn("chrM_not_rCRS", record["flags"])
        self.assertEqual(record["evidence"]["clinvar"]["status"], "not_found")
        self.assertIsNone(record["variant"]["hg38"]["key"])

    def test_rcrs_is_shared_between_builds(self):
        record = self.one("NC_012920.1:m.3243A>G", "hg19")
        self.assertEqual(record["variant"]["hg19"]["key"], "chrMT-3243-A-G")
        self.assertEqual(record["variant"]["hg38"]["key"], "chrM-3243-A-G")


class FileRuns(WithData):
    def setUp(self):
        self.out = tempfile.TemporaryDirectory()
        self.addCleanup(self.out.cleanup)
        self.root = Path(self.out.name)
        self.vcf = DEMO / "demo.hg19.vcf"

    def _cli(self, *args):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = main(list(args) + [f"data={self.data}"]) if args[0] != "demo" else main(list(args))
        return code, stdout.getvalue(), stderr.getvalue()

    def test_single_variant_prints_json(self):
        code, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19", "format=json")
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(out)[0]["evidence"]["clinvar"]["id"], 13961)

    def test_single_variant_prints_a_card_by_default(self):
        code, out, _ = self._cli("annotate", "BRAF V600E", "build=hg19")
        self.assertEqual(code, 0)
        for expected in ("HIGH", "Why this", "NM_004333", "13961", "https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/"):
            self.assertIn(expected, out)

    def test_dashes_typed_out_of_habit_still_work(self):
        code, out, _ = self._cli("annotate", "BRAF V600E", "-build", "hg19", "-format=json")
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(out)[0]["input"]["build"], "hg19")

    def test_file_input_gets_a_folder_with_every_report(self):
        folder = self.root / "reports"
        code, _, err = self._cli("annotate", str(self.vcf), f"out={folder}")
        self.assertEqual(code, 0, err)
        for name in ("report.html", "results.xlsx", "results.tsv"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertIn(str(folder / "report.html"), err)
        self.assertEqual(list(folder.glob("*.incomplete")), [])

    def test_existing_output_needs_the_word_replace(self):
        target = self.root / "out.tsv"
        target.write_text("keep me")
        code, _, err = self._cli("annotate", str(self.vcf), f"out={target}")
        self.assertEqual(code, 1)
        self.assertIn("replace", err)
        self.assertEqual(target.read_text(), "keep me")
        self.assertEqual(self._cli("annotate", str(self.vcf), f"out={target}", "replace")[0], 0)

    def test_workers_give_exactly_the_same_file_in_the_same_order(self):
        one, many = self.root / "one.tsv", self.root / "many.tsv"
        self.assertEqual(self._cli("annotate", str(self.vcf), f"out={one}", "workers=1")[0], 0)
        code, _, err = self._cli("annotate", str(self.vcf), f"out={many}", "workers=2")
        self.assertEqual(code, 0, err)
        self.assertEqual(one.read_text(encoding="utf-8"), many.read_text(encoding="utf-8"))

    def test_a_large_file_is_split_across_chunks_and_stays_in_order(self):
        lines = self.vcf.read_text(encoding="utf-8").splitlines(keepends=True)
        header = [line for line in lines if line.startswith("#")]
        records = [line for line in lines if not line.startswith("#")]
        large = self.root / "large.vcf"
        large.write_text("".join(header + records * 40), encoding="utf-8", newline="\n")
        one, many = self.root / "one.tsv", self.root / "many.tsv"
        self.assertEqual(self._cli("annotate", str(large), f"out={one}", "workers=1")[0], 0)
        code, _, err = self._cli("annotate", str(large), f"out={many}", "workers=2")
        self.assertEqual(code, 0, err)
        self.assertEqual(one.read_text(encoding="utf-8"), many.read_text(encoding="utf-8"))

    def test_error_gives_exit_code_and_message(self):
        code, _, err = self._cli("annotate", "BRAF W600E", "build=hg19")
        self.assertEqual(code, 1)
        self.assertIn("no transcript fits", err)
        self.assertIn("tried", err)

    def test_vcf_keeps_columns_and_reports_bad_lines(self):
        target = self.root / "out.tsv"
        code, _, err = self._cli("annotate", str(self.vcf), f"out={target}")
        self.assertEqual(code, 0, err)
        lines = target.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        self.assertEqual(header[:5], ["CHROM", "POS", "ID", "REF", "ALT"])
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
        braf = next(r for r in rows if r["POS"] == "140453136")
        self.assertEqual(braf["clinvar_id"], "13961")
        self.assertTrue(any(r["flags"] == "ref_mismatch" for r in rows))
        self.assertTrue(any(r["flags"] == "symbolic_allele_skipped" for r in rows))
        self.assertTrue(any("split_from_multiallelic" in r["flags"] for r in rows))

    def test_vcf_build_comes_from_the_header_and_a_contradiction_stops(self):
        code, _, err = self._cli("annotate", str(self.vcf), "build=hg38", f"out={self.root / 'x.tsv'}")
        self.assertEqual(code, 1)
        self.assertIn("wrong build", err)

    def test_excel_output_opens_as_a_valid_workbook(self):
        target = self.root / "out.xlsx"
        code, _, err = self._cli("annotate", str(self.vcf), f"out={target}")
        self.assertEqual(code, 0, err)
        with zipfile.ZipFile(target) as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn("xl/worksheets/sheet1.xml", archive.namelist())

    def test_list_file_reports_each_failure_and_carries_on(self):
        listing = self.root / "list.txt"
        listing.write_text("BRAF V600E\n# comment\n\nNOTAGENE V600E\nchr7-140453136-A-T\nTP53 R175H\n", encoding="utf-8")
        target = self.root / "out.jsonl"
        code, _, err = self._cli("annotate", str(listing), "build=hg19", f"out={target}")
        self.assertEqual(code, 0, err)
        rows = [json.loads(line) for line in target.read_text(encoding="utf-8").splitlines()]
        self.assertTrue(any(r.get("error") for r in rows))
        self.assertTrue(any(r.get("gene", {}).get("symbol") == "TP53" for r in rows if "gene" in r))

    def test_several_queries_separated_by_semicolons_run_as_a_list(self):
        target = self.root / "x.tsv"
        code, _, err = self._cli("annotate", "BRAF V600E; KRAS G12D", "build=hg19", f"out={target}")
        self.assertEqual(code, 0, err)
        lines = target.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
        self.assertEqual({r["input"] for r in rows}, {"BRAF V600E", "KRAS G12D"})
        self.assertEqual({r["gene"] for r in rows}, {"BRAF", "KRAS"})
        self.assertIn("lines: 2", err)

    def test_a_failing_query_in_a_semicolon_list_does_not_stop_the_others(self):
        with self.assertLogs("varshift", "WARNING") as logged:
            code, out, err = self._cli("annotate", "BRAF(V600); TP53 R175H", "build=hg19", "format=jsonl")
        self.assertEqual(code, 0, err)
        rows = [json.loads(line) for line in out.splitlines()]
        self.assertEqual(rows[0]["input"]["read_as"], "BRAF V600")
        self.assertTrue(rows[0]["error"])
        self.assertEqual(rows[1]["gene"]["symbol"], "TP53")
        self.assertIn("lines: 2, errors: 1", err)
        self.assertIn("query 1: ", "\n".join(logged.output))

    def test_one_query_with_a_trailing_semicolon_is_one_query(self):
        code, out, err = self._cli("annotate", "BRAF V600E;", "build=hg19", "format=json")
        self.assertEqual(code, 0, err)
        self.assertEqual(json.loads(out)[0]["variant"]["hg19"]["key"], "chr7-140453136-A-T")

    def test_status_and_doctor_and_selftest(self):
        self.assertEqual(self._cli("status")[0], 0)
        self.assertEqual(self._cli("selftest")[0], 0)

    def test_missing_data_says_what_to_do(self):
        with self.assertRaises(DataError) as caught:
            DataStore(self.root / "nothing")
        self.assertIn("varshift setup", str(caught.exception))

    @unittest.skipIf(os.name == "nt" or is_windows_mount(Path(tempfile.gettempdir())), "POSIX file modes")
    def test_vcf_results_are_private(self):
        target = self.root / "private.tsv"
        self._cli("annotate", str(self.vcf), f"out={target}")
        self.assertEqual(stat.S_IMODE(target.stat().st_mode), 0o600)

    def test_demo_command_runs(self):
        code, out, err = self._cli("demo", f"demodir={DEMO}", f"out={self.root / 'demo'}", "replace")
        self.assertEqual(code, 0, out + err)
        for folder in ("variants", "demo.hg19"):
            for name in ("report.html", "results.xlsx", "results.tsv"):
                self.assertTrue((self.root / "demo" / folder / name).exists(), f"{folder}/{name}")


if __name__ == "__main__":
    unittest.main()
