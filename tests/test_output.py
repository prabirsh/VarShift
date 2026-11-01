"""Output files: safe against spreadsheet formulas, private, and never mistaken for complete."""

import csv
import os
import stat
import tempfile
import unittest
import xml.etree.ElementTree as ElementTree
import zipfile
from pathlib import Path

from varsight.errors import InputError
from varsight.links import build_links, link_names
from varsight.output import ResultWriter, _guard_cell, columns, flatten, format_for
from varsight.xlsx import XlsxWriter

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
}


class Folder(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)


class FileHandling(Folder):
    def test_format_comes_from_the_extension(self):
        self.assertEqual(format_for("a.xlsx", None), "xlsx")
        self.assertEqual(format_for("a.tsv", None), "tsv")
        self.assertEqual(format_for(None, None), "json")
        self.assertEqual(format_for("a.xlsx", "tsv"), "tsv")
        with self.assertRaises(InputError):
            format_for("a.pdf", None)

    def test_incomplete_until_closed(self):
        target = self.root / "out.tsv"
        writer = ResultWriter(str(target), "tsv")
        writer.write(RECORD)
        self.assertTrue((self.root / "out.tsv.incomplete").exists())
        self.assertFalse(target.exists())
        writer.close()
        self.assertTrue(target.exists())
        self.assertFalse((self.root / "out.tsv.incomplete").exists())

    def test_an_existing_file_is_not_replaced_by_accident(self):
        target = self.root / "out.tsv"
        target.write_text("precious")
        with self.assertRaises(InputError):
            ResultWriter(str(target), "tsv")
        self.assertEqual(target.read_text(), "precious")
        writer = ResultWriter(str(target), "tsv", force=True)
        writer.close()
        self.assertNotEqual(target.read_text(), "precious")

    @unittest.skipIf(os.name == "nt", "POSIX file modes")
    def test_results_are_private(self):
        target = self.root / "out.tsv"
        writer = ResultWriter(str(target), "tsv")
        writer.close()
        self.assertEqual(stat.S_IMODE(target.stat().st_mode), 0o600)

    @unittest.skipIf(os.name == "nt", "symlinks need privileges on Windows")
    def test_a_symbolic_link_is_not_followed(self):
        real = self.root / "real.tsv"
        real.write_text("x")
        link = self.root / "link.tsv"
        link.symlink_to(real)
        with self.assertRaises(InputError):
            ResultWriter(str(link), "tsv", force=True)


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
        writer = ResultWriter(str(target), "tsv")
        writer.write(RECORD)
        writer.close()
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
            self.assertIn(f"link_{name}", columns())


class Excel(Folder):
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


if __name__ == "__main__":
    unittest.main()
