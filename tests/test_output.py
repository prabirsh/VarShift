"""Reports: safe against spreadsheet formulas, private, never mistaken for complete, and readable."""

import csv
import io
import json
import os
import re
import stat
import tempfile
import unittest
import xml.etree.ElementTree as ElementTree
import zipfile
from pathlib import Path
from unittest import mock

from varshift import report
from varshift.config import is_windows_mount
from varshift.errors import InputError
from varshift.links import build_links, link_names
from varshift.output import ExcelSink, ReportContext, ReportFolder, ScreenSink, _guard_cell, columns, flatten, format_of, open_file_sink, render_card
from varshift.xlsx import XlsxWriter

MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RECORD = {
    "input": {"raw": "=HYPERLINK(\"https://evil.example\")", "type": "protein", "build": "hg19"},
    "candidate": {"rank": 1, "of": 2},
    "variant": {"hg19": {"key": "chr7-140453136-A-T"}, "hg38": {"key": "chr7-140753336-A-T"}},
    "gene": {"symbol": "BRAF"},
    "transcript": {"accession": "NM_004333.6"},
    "hgvs": {"c": "NM_004333.6:c.1799T>A", "p": "NP_004324.2:p.Val600Glu"},
    "effect": {"type": "missense_variant", "exon": 15},
    "protein": {"uniprot": "P15056", "pfam": [{"accession": "PF07714"}], "domains_at_residue": []},
    "evidence": {"clinvar": {"status": "found", "id": 13961, "classification": "Conflicting"}, "dbsnp": "rs113488022"},
    "links": [{"name": "ClinVar", "label": "confirmed_id", "url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/"}],
    "versions": {"clinvar.hg19": "2026-10-25"},
    "flags": ["multiple_candidates"],
    "priority": {"level": 2, "label": "Medium", "reasons": ["ClinVar: conflicting classifications"]},
}


def with_priority(level, label, gene="BRAF"):
    return {**RECORD, "gene": {"symbol": gene}, "priority": {"level": level, "label": label, "reasons": [f"because {label}"]}}


class Folder(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.context = ReportContext(input_name="test.vcf", versions={"clinvar.hg19": "2026-10-25"})


class FileHandling(Folder):
    def test_format_comes_from_the_extension(self):
        self.assertEqual(format_of("a.xlsx"), "xlsx")
        self.assertEqual(format_of("a.HTML"), "html")
        self.assertEqual(format_of("a.tsv"), "tsv")
        self.assertIsNone(format_of("a_folder"))
        with self.assertRaises(InputError):
            open_file_sink(self.root / "a.pdf", self.context)

    def test_incomplete_until_closed(self):
        target = self.root / "out.tsv"
        sink = open_file_sink(target, self.context)
        sink.write(RECORD)
        self.assertTrue((self.root / "out.tsv.incomplete").exists())
        self.assertFalse(target.exists())
        sink.close()
        self.assertTrue(target.exists())
        self.assertFalse((self.root / "out.tsv.incomplete").exists())

    def test_an_existing_file_is_not_replaced_by_accident(self):
        target = self.root / "out.tsv"
        target.write_text("precious")
        with self.assertRaises(InputError):
            open_file_sink(target, self.context)
        self.assertEqual(target.read_text(), "precious")
        open_file_sink(target, self.context, replace=True).close()
        self.assertNotEqual(target.read_text(), "precious")

    @unittest.skipIf(os.name == "nt" or is_windows_mount(Path(tempfile.gettempdir())), "POSIX file modes")
    def test_every_report_is_private(self):
        folder = ReportFolder(self.root / "reports", self.context)
        folder.write(RECORD)
        folder.close()
        for path in folder.paths:
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600, path.name)

    @unittest.skipIf(os.name == "nt", "symlinks need privileges on Windows")
    def test_a_symbolic_link_is_not_followed(self):
        real = self.root / "real.tsv"
        real.write_text("x")
        link = self.root / "link.tsv"
        link.symlink_to(real)
        with self.assertRaises(InputError):
            open_file_sink(link, self.context, replace=True)


class Tsv(Folder):
    def test_formula_text_is_neutralised(self):
        self.assertEqual(_guard_cell("=1+1"), "'=1+1")
        self.assertEqual(_guard_cell("+cmd"), "'+cmd")
        self.assertEqual(_guard_cell("@SUM(A1)"), "'@SUM(A1)")
        self.assertEqual(_guard_cell("-abc"), "'-abc")
        self.assertEqual(_guard_cell("-0.5"), "-0.5")
        self.assertEqual(_guard_cell("BRAF"), "BRAF")
        self.assertEqual(_guard_cell(None), "")

    def test_tabs_and_newlines_cannot_break_a_row(self):
        self.assertEqual(_guard_cell("a\tb\nc"), "a b c")

    def test_a_row_has_one_value_per_column(self):
        target = self.root / "out.tsv"
        sink = open_file_sink(target, self.context)
        sink.write(RECORD)
        sink.close()
        with open(target, encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(rows[0]), len(rows[1]))
        row = dict(zip(rows[0], rows[1]))
        self.assertEqual(row["hg19_key"], "chr7-140453136-A-T")
        self.assertEqual(row["clinvar_id"], "13961")
        self.assertTrue(row["input"].startswith("'="))

    def test_every_database_gets_a_link_column(self):
        for name in ("ClinVar", "NCBI Gene", "HGNC", "Ensembl", "UniProt", "Pfam", "ProteinPaint", "Franklin", "COSMIC", "VarSome", "OncoKB"):
            self.assertIn(f"link_{name}", columns(self.context))

    def test_advanced_columns_appear_only_when_installed(self):
        self.assertNotIn("revel", columns(ReportContext()))
        self.assertNotIn("topmed_af", columns(ReportContext()))
        with_both = columns(ReportContext(installed=["revel.hg19", "dbsnp.hg38"]))
        for name in ("revel", "topmed_af", "dbsnp_rsid"):
            self.assertIn(name, with_both)
        self.assertNotIn("cosmic_id", with_both)


class ExcelCells(Folder):
    def _write(self, rows):
        path = self.root / "t.xlsx"
        writer = XlsxWriter(path)
        writer.add_sheet("Sheet")
        for row in rows:
            writer.write_row(row)
        writer.close()
        return path

    def _cells(self, path):
        with zipfile.ZipFile(path) as archive:
            root = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        return root.findall(f".//{MAIN}c")

    def test_file_is_well_formed(self):
        path = self._write([["a", 1, 2.5, None, True]])
        with zipfile.ZipFile(path) as archive:
            self.assertIsNone(archive.testzip())
            for name in archive.namelist():
                ElementTree.fromstring(archive.read(name))
            self.assertIn("xl/workbook.xml", archive.namelist())

    def test_a_formula_stays_as_text(self):
        cells = self._cells(self._write([['=HYPERLINK("https://evil.example")', "@SUM(A1)"]]))
        for cell in cells:
            self.assertEqual(cell.get("t"), "inlineStr")
            self.assertIsNone(cell.find(f"{MAIN}f"))

    def test_markup_is_escaped(self):
        cells = self._cells(self._write([["<b>&\"x\"</b>"]]))
        self.assertEqual(cells[0].find(f"{MAIN}is/{MAIN}t").text, '<b>&"x"</b>')

    def test_control_characters_do_not_break_the_file(self):
        cells = self._cells(self._write([["a\x00b\x1fc"]]))
        self.assertEqual(cells[0].find(f"{MAIN}is/{MAIN}t").text, "abc")

    def test_only_https_becomes_a_link(self):
        path = self._write([["https://www.ncbi.nlm.nih.gov/gene/673", "http://plain.example", "javascript:alert(1)", "file:///etc/passwd"]])
        with zipfile.ZipFile(path) as archive:
            rels = archive.read("xl/worksheets/_rels/sheet1.xml.rels").decode()
        self.assertIn("ncbi.nlm.nih.gov", rels)
        for bad in ("plain.example", "javascript", "passwd"):
            self.assertNotIn(bad, rels)

    def test_many_rows_stay_cheap(self):
        path = self._write([[f"row{i}", i] for i in range(5000)])
        self.assertEqual(len(self._cells(path)), 10000)


class Workbook(Folder):
    def setUp(self):
        super().setUp()
        self.path = self.root / "results.xlsx"
        sink = open_file_sink(self.path, self.context)
        for level, label, gene in ((4, "Low", "AAA"), (1, "High", "BBB"), (3, "Review", "CCC"), (2, "Medium", "DDD"), (1, "High", "BBB")):
            sink.write(with_priority(level, label, gene))
        sink.close()

    def _sheets(self):
        with zipfile.ZipFile(self.path) as archive:
            book = ElementTree.fromstring(archive.read("xl/workbook.xml"))
            rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
            targets = {rel.get("Id"): rel.get("Target") for rel in rels}
            sheets = {}
            key = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            for sheet in book.find(f"{MAIN}sheets"):
                root = ElementTree.fromstring(archive.read("xl/" + targets[sheet.get(key)].lstrip("/").removeprefix("xl/")))
                rows = []
                for row in root.findall(f".//{MAIN}row"):
                    values = []
                    for cell in row:
                        text = cell.find(f"{MAIN}is/{MAIN}t")
                        number = cell.find(f"{MAIN}v")
                        values.append(text.text if text is not None else number.text if number is not None else None)
                    rows.append(values)
                sheets[sheet.get("name")] = rows
        return sheets

    def test_sheets_come_in_the_order_people_read_them(self):
        self.assertEqual(list(self._sheets()), ["Summary", "Priority", "All results", "Warnings explained", "Data versions"])

    def test_priority_sheet_is_sorted_and_leaves_out_low(self):
        rows = self._sheets()["Priority"]
        self.assertEqual([row[0] for row in rows[1:]], ["High", "High", "Medium", "Review"])
        self.assertEqual(rows[1][1], "because High")

    def test_summary_counts_each_level(self):
        rows = self._sheets()["Summary"]
        counts = {row[0]: row[1] for row in rows if len(row) > 1}
        self.assertEqual((counts["High"], counts["Medium"], counts["Review"], counts["Low"]), ("2", "1", "1", "1"))

    def test_all_results_keeps_every_row(self):
        self.assertEqual(len(self._sheets()["All results"]), 6)

    def test_warnings_are_explained(self):
        row = self._sheets()["Warnings explained"][1]
        self.assertEqual(row[0], "multiple_candidates")
        self.assertTrue(row[2])
        self.assertTrue(row[3])

    def test_past_excel_row_limit_the_run_finishes_and_the_summary_says_so(self):
        self.path = self.root / "big.xlsx"
        with mock.patch.object(ExcelSink, "row_limit", 3):
            sink = open_file_sink(self.path, self.context)
            for _ in range(5):
                sink.write(with_priority(3, "Review"))
            sink.close()
        sheets = self._sheets()
        self.assertEqual(len(sheets["All results"]), 3, "the header and two results fit")
        summary = [row[0] for row in sheets["Summary"] if row]
        self.assertTrue(any(text and text.startswith("3 results are not in All results") for text in summary), summary)
        counts = {row[0]: row[1] for row in sheets["Summary"] if len(row) > 1}
        self.assertEqual(counts["Review"], "5", "the counts still include every result")


class HtmlReport(Folder):
    def _payload(self, path):
        text = path.read_text(encoding="utf-8")
        match = re.search(r'<script type="application/json" id="data">(.*?)</script>', text, re.S)
        return text, json.loads(match.group(1)) if match else None

    def _write(self, records):
        path = self.root / "report.html"
        sink = open_file_sink(path, self.context)
        for record in records:
            sink.write(record)
        sink.close()
        return path

    def test_the_data_is_inside_the_page(self):
        text, payload = self._payload(self._write([with_priority(1, "High")]))
        self.assertNotIn("/*VARSHIFT_DATA*/", text)
        self.assertIsNotNone(payload, "report data block not found")
        self.assertEqual(payload["total"], 1)
        self.assertEqual(payload["rows"][0]["g"], "BRAF")
        self.assertEqual(payload["input"], "test.vcf")

    def test_the_page_needs_nothing_from_the_internet(self):
        text, _ = self._payload(self._write([RECORD]))
        page = text[: text.index('id="data"')]
        self.assertIsNone(re.search(r'<(script|link|img)[^>]+(src|href)="https?:', page))
        self.assertNotIn("@import", page)

    def test_input_text_cannot_close_the_script_or_add_markup(self):
        hostile = {**RECORD, "input": {"raw": "</script><script>alert(1)</script>", "type": "x", "build": "hg19"}}
        text, payload = self._payload(self._write([hostile]))
        self.assertNotIn("<script>alert(1)", text)
        self.assertEqual(payload["rows"][0]["i"], "</script><script>alert(1)</script>")

    def test_the_page_writes_text_safely(self):
        template = (Path(report.__file__).parent / "data" / "report_template.html").read_text(encoding="utf-8")
        self.assertNotIn("innerHTML", template)
        self.assertNotIn("document.write", template)
        self.assertNotIn("eval(", template)

    def test_when_too_large_low_results_leave_first_and_the_page_says_so(self):
        with mock.patch.object(report, "MAX_ROWS", 3):
            path = self._write([with_priority(4, "Low"), with_priority(1, "High"), with_priority(4, "Low"), with_priority(2, "Medium"), with_priority(3, "Review")])
        _, payload = self._payload(path)
        self.assertEqual(sorted(row["p"] for row in payload["rows"]), [1, 2, 3])
        self.assertEqual((payload["total"], payload["leftOut"]), (5, 2))

    def test_when_full_a_high_result_still_replaces_a_lower_one(self):
        with mock.patch.object(report, "MAX_ROWS", 3):
            path = self._write([with_priority(3, "Review"), with_priority(3, "Review"), with_priority(2, "Medium"), with_priority(1, "High"), with_priority(4, "Low")])
        _, payload = self._payload(path)
        self.assertEqual(sorted(row["p"] for row in payload["rows"]), [1, 2, 3])
        self.assertEqual((payload["total"], payload["leftOut"]), (5, 2))

    def test_an_error_row_is_kept(self):
        failed = {"input": {"raw": "NOTAGENE V600E", "type": "protein", "build": "hg19"}, "error": "no such gene", "flags": []}
        _, payload = self._payload(self._write([failed]))
        self.assertEqual(payload["rows"][0]["err"], "no such gene")


class Screen(Folder):
    def test_json_is_valid(self):
        stream = io.StringIO()
        sink = ScreenSink("json", self.context, stream)
        sink.write(RECORD)
        sink.close()
        self.assertEqual(json.loads(stream.getvalue())[0]["gene"]["symbol"], "BRAF")

    def test_tsv_has_a_header(self):
        stream = io.StringIO()
        sink = ScreenSink("tsv", self.context, stream)
        sink.write(RECORD)
        sink.close()
        lines = stream.getvalue().splitlines()
        self.assertTrue(lines[0].startswith("priority\t"))
        self.assertEqual(len(lines), 2)

    def test_a_build_that_is_not_installed_is_not_called_unmappable(self):
        record = json.loads(json.dumps(RECORD))
        record["variant"]["hg38"]["key"] = None
        record["flags"] = ["other_build_not_installed"]
        self.assertIn("hg38 not installed", render_card(record))
        record["flags"] = []
        self.assertIn("hg38 no mapping", render_card(record))

    def test_a_failed_input_says_why(self):
        card = render_card({"input": {"raw": "x"}, "error": "no such gene"})
        self.assertIn("no such gene", card)


class Links(unittest.TestCase):
    VALUES = {
        "build": "hg19",
        "contig": "chr7",
        "contig_chr": "chr7",
        "contig_bare": "7",
        "pos": 140453136,
        "ref": "A",
        "alt": "T",
        "gene": "BRAF",
        "protein_short": "V600E",
        "clinvar_id": 13961,
        "rsid": "113488022",
        "ncbi_gene": "673",
        "hgnc_id": "HGNC:1097",
        "ensembl_gene": "ENSG00000157764",
        "omim": "164757",
        "uniprot": "P15056",
        "pfam": "PF07714",
        "refseq": "NM_004333.6",
        "civic_url": "https://civicdb.org/links/variants/12",
        "cosmic_query": "BRAF V600E",
        "pubmed_query": "BRAF V600E",
    }

    def _urls(self, values):
        return {link["name"]: link for link in build_links(values)}

    def test_every_requested_database_has_a_link(self):
        links = self._urls(self.VALUES)
        for name in ("ClinVar", "NCBI Gene", "NCBI RefSeq", "HGNC", "Ensembl", "UniProt", "Pfam", "ProteinPaint", "Franklin", "COSMIC", "VarSome", "OncoKB"):
            self.assertTrue(links[name]["url"], name)
            self.assertTrue(links[name]["url"].startswith("https://"), name)

    def test_exact_urls_on_hg19(self):
        links = self._urls(self.VALUES)
        self.assertEqual(links["ClinVar"]["url"], "https://www.ncbi.nlm.nih.gov/clinvar/variation/13961/")
        self.assertEqual(links["Franklin"]["url"], "https://franklin.genoox.com/clinical-db/variant/snp/chr7-140453136-A-T")
        self.assertEqual(links["VarSome"]["url"], "https://varsome.com/variant/hg19/7-140453136-A-T")
        self.assertEqual(links["OncoKB"]["url"], "https://www.oncokb.org/gene/BRAF/somatic/V600E")
        self.assertEqual(links["Ensembl"]["url"], "https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000157764")

    def test_hgnc_id_keeps_its_colon(self):
        links = self._urls(self.VALUES)
        self.assertEqual(links["HGNC"]["url"], "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1097")
        self.assertEqual(links["ClinGen"]["url"], "https://search.clinicalgenome.org/kb/genes/HGNC:1097")

    def test_hg38_forms(self):
        links = self._urls({**self.VALUES, "build": "hg38", "pos": 140753336})
        self.assertEqual(links["Franklin"]["url"], "https://franklin.genoox.com/clinical-db/variant/snp/chr7-140753336-A-T-hg38")
        self.assertEqual(links["Ensembl"]["url"], "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000157764")

    def test_missing_id_gives_a_reason_and_no_invented_link(self):
        links = self._urls({**self.VALUES, "rsid": None, "pfam": None})
        self.assertIsNone(links["dbSNP"]["url"])
        self.assertIn("rs number", links["dbSNP"]["note"])
        self.assertIsNone(links["Pfam"]["url"])

    def test_not_in_clinvar_is_labelled_as_a_search(self):
        link = self._urls({**self.VALUES, "clinvar_id": None})["ClinVar"]
        self.assertEqual(link["label"], "built")
        self.assertIn("Search only", link["note"])

    def test_input_cannot_change_the_address(self):
        hostile = {**self.VALUES, "gene": "../../evil?x=1#", "protein_short": "a/b c&d=e"}
        link = self._urls(hostile)["OncoKB"]["url"]
        self.assertTrue(link.startswith("https://www.oncokb.org/gene/"))
        self.assertNotIn("../", link)
        self.assertNotIn(" ", link)
        self.assertIn("%2F", link)

    def test_a_non_https_url_from_data_is_dropped(self):
        link = self._urls({**self.VALUES, "civic_url": "javascript:alert(1)"})["CIViC"]
        self.assertIsNone(link["url"])

    def test_names_are_unique(self):
        names = link_names()
        self.assertEqual(len(names), len(set(names)))


class Flatten(unittest.TestCase):
    def test_not_in_clinvar_is_stated_with_its_version(self):
        row = flatten({"evidence": {"clinvar": {"status": "not_found", "version": "2026-10-25"}}})
        self.assertEqual(row["clinvar_classification"], "not in ClinVar (version 2026-10-25)")

    def test_an_error_record_still_flattens(self):
        row = flatten({"input": {"raw": "x", "type": "vcf", "build": "hg19"}, "flags": ["ref_mismatch"], "error": "bad"})
        self.assertEqual((row["flags"], row["error"]), ("ref_mismatch", "bad"))

    def test_advanced_values_reach_their_columns(self):
        row = flatten({"evidence": {"clinvar": {}, "advanced": {"revel": {"status": "found", "score": 0.93}, "dbsnp": {"rsid": "1", "af_topmed": 0.01}}}})
        self.assertEqual((row["revel"], row["topmed_af"]), (0.93, 0.01))


if __name__ == "__main__":
    unittest.main()
