"""The advanced databases on samples cut from each real public file, or built to the provider's documented
columns where the file needs an account (COSMIC, OMIM, SpliceAI, PrimateAI). Long INFO fields were shortened.
"""

import contextlib
import gzip
import io
import sqlite3
import tempfile
import unittest
import zipfile
from pathlib import Path

from tests.helpers import tiny_bigwig
from varsight.advanced import ALIAS, DATABASES, OWN_FILE, artifact_names, flatten_advanced, get_database, output_columns
from varsight.annotate import Annotator
from varsight.bigwig import BigWig
from varsight.errors import DataError, SetupError
from varsight.explain import run_explain
from varsight.genome import open_build_genome
from varsight.pipeline import advanced_builders as builders
from varsight.pipeline.benchmark import run_benchmark
from varsight.pipeline.demo import build_demo_store
from varsight.store import DB_DIR, GENOME_DIR, DataStore

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def lines_of(name):
    return (FIXTURES / name).read_text(encoding="utf-8").splitlines(keepends=True)


def rows_of(path, sql, *parameters):
    connection = sqlite3.connect(path)
    try:
        return connection.execute(sql, parameters).fetchall()
    finally:
        connection.close()


class WithDemo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        cls.root = Path(cls.folder.name)
        cls.data = cls.root / "data"
        build_demo_store(ROOT / "demo", cls.data)
        cls.genome = open_build_genome(cls.data / GENOME_DIR, "hg19")
        cls.regions = builders.CodingRegions(cls.data / DB_DIR / "transcripts.hg19.sqlite")

    @classmethod
    def tearDownClass(cls):
        cls.genome.close()
        cls.folder.cleanup()

    def out(self, name):
        path = self.root / name
        path.unlink(missing_ok=True)
        return path


class Regions(WithDemo):
    def test_keeps_exons_and_their_edges_and_drops_the_rest(self):
        self.assertTrue(self.regions.holds("chr7", 140453136))
        self.assertFalse(self.regions.holds("chr7", 1000))
        self.assertFalse(self.regions.holds("chrZ", 1000))
        self.assertLess(self.regions.total_bases(), 10_000_000)

    def test_padding_widens_each_exon_by_exactly_that_much(self):
        path = self.out("tiny.sqlite")
        connection = sqlite3.connect(path)
        connection.execute("CREATE TABLE transcripts(contig TEXT, exons TEXT)")
        connection.execute("INSERT INTO transcripts VALUES('chr1', '999-1100,2000-2100')")
        connection.commit()
        connection.close()
        regions = builders.CodingRegions(path, pad=50)
        for position, expected in ((949, False), (950, True), (1150, True), (1151, False), (1951, True), (2150, True), (2151, False)):
            self.assertEqual(regions.holds("chr1", position), expected, position)


class Dbsnp(WithDemo):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.path = cls.root / "dbsnp.sqlite"
        cls.stats = builders.build_dbsnp(lines_of("dbsnp.hg19.sample.txt"), cls.path, "hg19", cls.genome, cls.regions)

    def _row(self, pos, alt):
        return rows_of(self.path, "SELECT rsid, af_topmed, af_alfa, af_gnomad_exomes, af_gnomad_genomes FROM variants WHERE contig='chr7' AND pos=? AND alt=?", pos, alt)

    def test_the_accession_becomes_a_chromosome(self):
        self.assertGreater(self.stats["rows"], 20)
        self.assertEqual(rows_of(self.path, "SELECT DISTINCT contig FROM variants"), [("chr7",)])

    def test_rs_number_of_braf_v600e(self):
        self.assertEqual(self._row(140453136, "T")[0][0], "113488022")

    def test_each_alt_takes_its_own_frequency(self):
        # FREQ lists the reference first: GnomAD_exomes:1,7.141e-07,.,. and dbGaP_PopFreq:1,0,3.724e-05,.
        first = self._row(140453132, "A")[0]
        second = self._row(140453132, "C")[0]
        self.assertEqual((first[0], first[2], first[3]), ("121913365", 0.0, 7.141e-07))
        self.assertEqual((second[2], second[3]), (3.724e-05, None))

    def test_an_indel_is_stored_left_aligned(self):
        self.assertTrue(rows_of(self.path, "SELECT 1 FROM variants WHERE rsid='121913226' AND length(ref) > 1"))

    def test_rows_outside_the_regions_are_dropped(self):
        stats = builders.build_dbsnp(["NC_000007.13\t1000\trs1\tA\tT\t.\t.\tRS=1\n"] + lines_of("dbsnp.hg19.sample.txt"), self.out("d2.sqlite"), "hg19", self.genome, self.regions)
        self.assertEqual(stats["outside_scope"], 1)

    def test_the_wrong_build_gives_no_rows_and_a_clear_error(self):
        with self.assertRaises(SetupError) as caught:
            builders.build_dbsnp(lines_of("dbsnp.hg19.sample.txt"), self.out("d3.sqlite"), "hg38", self.genome, None)
        self.assertIn("no rows", str(caught.exception))


class Frequencies(WithDemo):
    def test_thousand_genomes_real_format(self):
        path = self.out("g1k.sqlite")
        stats = builders.build_thousand_genomes(lines_of("thousand_genomes.sample.txt"), path, "hg19", self.genome, None)
        self.assertGreaterEqual(stats["rows"], 1)
        row = rows_of(path, "SELECT af, af_afr, af_eas FROM variants WHERE pos=10505")[0]
        self.assertEqual(row, (0.000199681, 0.0008, 0.0))

    def test_gnomad_real_format(self):
        path = self.out("gnomad.sqlite")
        stats = builders.build_gnomad([("7", lines_of("gnomad.sample.txt"))], path, "hg19", self.genome, self.regions)
        self.assertEqual(stats["rows"], 4)
        row = rows_of(path, "SELECT filter, ac, an, af, nhomalt, popmax, af_popmax FROM variants WHERE pos=140453136")[0]
        self.assertEqual(row, ("PASS", 1, 251260, 3.97994e-06, 0, "sas", 3.26669e-05))
        self.assertEqual(rows_of(path, "SELECT filter, popmax FROM variants WHERE pos=140453132")[0], ("AC0;RF", None))

    def test_gnomad_sv_real_format(self):
        path = self.out("sv.sqlite")
        stats = builders.build_gnomad_sv(lines_of("gnomad_sv.sample.txt"), path, "hg19")
        self.assertEqual(stats["rows"], 5)
        row = rows_of(path, "SELECT start, end, sv_type, sv_length, af, filter FROM intervals WHERE sv_id='gnomAD-SV_v2.1_DEL_1_1'")[0]
        self.assertEqual(row, (21000, 26000, "DEL", 5000, 0.016073, "UNSTABLE_AF_PCRMINUS;LOW_CALL_RATE"))


class Scores(WithDemo):
    def test_revel_reads_the_column_for_the_build(self):
        archive_path = self.root / "revel.zip"
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.write(FIXTURES / "revel.sample.csv", "revel.csv")
        path = self.out("revel.sqlite")
        builders.build_revel(archive_path, path, "hg19")
        self.assertEqual(rows_of(path, "SELECT contig, pos, score FROM variants WHERE ref='G' AND alt='A' LIMIT 1")[0], ("chr1", 35142, 0.027))
        path = self.out("revel38.sqlite")
        self.assertGreater(builders.build_revel(archive_path, path, "hg38")["rows"], 5)

    def test_alphamissense_real_format(self):
        path = self.out("am.sqlite")
        stats = builders.build_alphamissense(FIXTURES / "alphamissense.sample.tsv", path, "hg19")
        self.assertGreater(stats["rows"], 10)
        row = rows_of(path, "SELECT contig, am_class FROM variants LIMIT 1")[0]
        self.assertEqual(row[0], "chr1")
        self.assertIn(row[1], ("benign", "pathogenic", "ambiguous"))

    def test_spliceai_keeps_the_highest_delta(self):
        lines = [
            "##fileformat=VCFv4.0\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "7\t140453136\t.\tA\tT\t.\t.\tSpliceAI=T|BRAF|0.00|0.03|0.31|0.00|-2|34|-1|8\n",
            "7\t140453136\t.\tA\tC\t.\t.\tSpliceAI=C|BRAF|0.00|0.00|0.01|0.00|-2|34|-1|8\n",
        ]
        path = self.out("splice.sqlite")
        stats = builders.build_spliceai(lines, path, "hg19", self.genome, self.regions)
        self.assertEqual(stats["rows"], 2)
        self.assertEqual(rows_of(path, "SELECT max_delta, gene FROM variants WHERE alt='T'")[0], (0.31, "BRAF"))

    def test_primateai_documented_columns(self):
        file = self.root / "primate.tsv"
        file.write_text("# PrimateAI\nchr\tpos\tref\talt\tprimateDL_score\nchr1\t100\tA\tG\t0.7\n", encoding="utf-8")
        path = self.out("primate.sqlite")
        self.assertEqual(builders.build_primateai(file, path, "hg19")["rows"], 1)
        self.assertEqual(rows_of(path, "SELECT score FROM variants")[0], (0.7,))

    def test_primateai_wrong_file_is_refused(self):
        file = self.root / "wrong.tsv"
        file.write_text("a\tb\n1\t2\n", encoding="utf-8")
        with self.assertRaises(SetupError):
            builders.build_primateai(file, self.out("w.sqlite"), "hg19")


class Cancer(WithDemo):
    def test_hotspots_real_format(self):
        path = self.out("hot.sqlite")
        stats = builders.build_hotspots(FIXTURES / "hotspots.sample.json", path)
        self.assertEqual(stats["rows"], 5)
        gene, residue, start, tumours = rows_of(path, "SELECT gene, residue, start, tumor_count FROM hotspots WHERE gene='NRAS'")[0]
        self.assertEqual((gene, residue, start, tumours), ("NRAS", "Q61", 61, 422))
        self.assertIn("R:204", rows_of(path, "SELECT variants FROM hotspots WHERE gene='NRAS'")[0][0])

    def test_hotspots_format_change_is_loud(self):
        file = self.root / "hot.json"
        file.write_text('[{"unexpected": 1}]', encoding="utf-8")
        with self.assertRaises(SetupError):
            builders.build_hotspots(file, self.out("h2.sqlite"))

    def test_cosmic_counts_samples_for_one_change(self):
        file = self.root / "cosmic.tsv"
        file.write_text(
            "GENE_SYMBOL\tSAMPLE_NAME\tGENOMIC_MUTATION_ID\tMUTATION_CDS\tMUTATION_AA\tHGVSG\n"
            "BRAF\tS1\tCOSV56056643\tc.1799T>A\tp.V600E\t7:g.140453136A>T\n"
            "BRAF\tS2\tCOSV56056643\tc.1799T>A\tp.V600E\t7:g.140453136A>T\n"
            "BRAF\tS3\tCOSV1\tc.1798G>A\tp.V600M\t7:g.140453137C>T\n"
            "BRAF\tS4\tCOSV2\t\t\t\n"
            "BRAF\tS5\tCOSV3\tc.1T>A\tp.M1K\t7:g.140453136G>T\n",
            encoding="utf-8",
        )
        path = self.out("cosmic.sqlite")
        stats = builders.build_cosmic(file, path, "hg19", self.genome)
        self.assertEqual(rows_of(path, "SELECT cosmic_id, samples FROM variants WHERE pos=140453136 AND alt='T'")[0], ("COSV56056643", 2))
        self.assertEqual((stats["no_position"], stats["skipped"]), (1, 1))

    def test_cosmic_wrong_file_is_refused(self):
        file = self.root / "wrong.tsv"
        file.write_text("a\tb\n1\t2\n", encoding="utf-8")
        with self.assertRaises(SetupError):
            builders.build_cosmic(file, self.out("c2.sqlite"), "hg19", self.genome)

    def test_census_documented_columns(self):
        file = self.root / "census.tsv"
        file.write_text("GENE_SYMBOL\tTIER\tROLE_IN_CANCER\tTUMOUR_TYPES_SOMATIC\nBRAF\t1\toncogene, fusion\tmelanoma, colorectal\n", encoding="utf-8")
        path = self.out("census.sqlite")
        builders.build_cosmic_census(file, path)
        self.assertEqual(rows_of(path, "SELECT * FROM genes")[0], ("BRAF", "1", "oncogene, fusion", "melanoma, colorectal"))

    def test_fusions_count_both_directions(self):
        file = self.root / "fusions.csv"
        file.write_text("FIVE_PRIME_GENE_SYMBOL,THREE_PRIME_GENE_SYMBOL\nKIAA1549_ENST1,BRAF_ENST2\nKIAA1549,BRAF\n", encoding="utf-8")
        path = self.out("fusions.sqlite")
        builders.build_cosmic_fusions(file, path)
        self.assertEqual(rows_of(path, "SELECT partner, samples FROM fusions WHERE symbol='BRAF'")[0], ("KIAA1549", 2))
        self.assertEqual(rows_of(path, "SELECT partner FROM fusions WHERE symbol='KIAA1549'")[0], ("BRAF",))


class Genes(WithDemo):
    def test_decipher_and_ddg2p_real_format(self):
        path = self.out("decipher.sqlite")
        stats = builders.build_decipher(FIXTURES / "decipher_hi.sample.txt", FIXTURES / "ddg2p.sample.csv", path)
        self.assertEqual(stats["hi_genes"], 6)
        self.assertEqual(rows_of(path, "SELECT hi_score, hi_percent FROM genes WHERE symbol='HYPM'")[0], (0.000043663, "100%"))
        self.assertEqual(rows_of(path, "SELECT inheritance, confidence FROM ddg2p WHERE symbol='HMX1'")[0], ("biallelic_autosomal", "definitive"))

    def test_clingen_real_format(self):
        path = self.out("clingen.sqlite")
        stats = builders.build_clingen_variants(FIXTURES / "clingen_erepo.sample.txt", path)
        self.assertEqual(stats["rows"], 3)
        gene, assertion = rows_of(path, "SELECT gene, assertion FROM classifications WHERE variation_id=586")[0]
        self.assertEqual(gene, "PAH")
        self.assertIn("athogenic", assertion)

    def test_omim_documented_header(self):
        file = self.root / "genemap2.txt"
        file.write_text(
            "# Copyright (c) OMIM\n# Chromosome\tGenomic Position Start\tMIM Number\tApproved Gene Symbol\tPhenotypes\n"
            "chr7\t140719327\t164757\tBRAF\tCardiofaciocutaneous syndrome, 115150 (3)\nchr1\t1\t100000\t\tnone\n",
            encoding="utf-8",
        )
        path = self.out("omim.sqlite")
        self.assertEqual(builders.build_omim(file, path)["rows"], 1)
        self.assertEqual(rows_of(path, "SELECT * FROM genes")[0], ("BRAF", "164757", "Cardiofaciocutaneous syndrome, 115150 (3)"))

    def test_omim_wrong_file_is_refused(self):
        file = self.root / "wrong.txt"
        file.write_text("not a genemap\n", encoding="utf-8")
        with self.assertRaises(SetupError):
            builders.build_omim(file, self.out("o2.sqlite"))


class BigWigFile(unittest.TestCase):
    def setUp(self):
        self.reader = BigWig(io.BytesIO(tiny_bigwig()))

    def test_fixed_step_values(self):
        self.assertEqual([self.reader.value("chr7", p) for p in (100, 101, 102)], [1.5, 2.5, 3.5])

    def test_variable_step_values(self):
        self.assertEqual(self.reader.value("chr7", 200), -1.0)
        self.assertEqual(self.reader.value("chr7", 205), 7.25)
        self.assertIsNone(self.reader.value("chr7", 202))

    def test_bedgraph_values(self):
        self.assertEqual(self.reader.value("chr7", 300), 9.0)
        self.assertEqual(self.reader.value("chr7", 309), 9.0)

    def test_outside_the_data_gives_nothing(self):
        for contig, position in (("chr7", 50), ("chr7", 150), ("chr7", 999), ("chr1", 100)):
            self.assertIsNone(self.reader.value(contig, position))

    def test_contig_names(self):
        self.assertTrue(self.reader.has("chr7"))
        self.assertFalse(self.reader.has("chr8"))

    def test_another_file_type_is_refused(self):
        with self.assertRaises(DataError):
            BigWig(io.BytesIO(b"\0" * 64))


class Catalog(unittest.TestCase):
    def test_every_database_the_user_named_is_listed(self):
        for name in ("dbsnp", "gnomad", "thousand_genomes", "revel", "topmed", "cosmic", "primateai", "clingen_variants", "spliceai", "decipher", "gnomad_sv", "cosmic_fusions", "omim", "cosmic_census", "phylop", "hotspots"):
            self.assertIsNotNone(get_database(name), name)

    def test_hyphens_and_case_in_a_name_are_fine(self):
        self.assertEqual(get_database("Gnomad-SV").name, "gnomad_sv")

    def test_licensed_databases_are_never_fetched_automatically(self):
        for name in ("cosmic", "cosmic_census", "cosmic_fusions", "omim", "spliceai", "primateai"):
            self.assertEqual(DATABASES[name].access, OWN_FILE, name)

    def test_topmed_points_at_dbsnp(self):
        self.assertEqual(DATABASES["topmed"].access, ALIAS)
        self.assertEqual(DATABASES["topmed"].alias_of, "dbsnp")
        self.assertEqual(artifact_names(DATABASES["topmed"], ["hg19"]), [])

    def test_artifact_names_follow_the_builds(self):
        self.assertEqual(artifact_names(DATABASES["revel"], ["hg19", "hg38"]), ["revel.hg19", "revel.hg38"])
        self.assertEqual(artifact_names(DATABASES["omim"], ["hg19", "hg38"]), ["omim"])

    def test_every_artifact_is_in_the_download_catalog(self):
        from varsight.pipeline.manager import load_catalog

        catalog = load_catalog()
        for database in DATABASES.values():
            for name in artifact_names(database, ["hg19", "hg38"]):
                self.assertIn(name, catalog["artifacts"], name)

    def test_every_column_can_be_filled(self):
        evidence = {"advanced": {"revel": {"score": 0.5}}}
        self.assertEqual(flatten_advanced(evidence), {"revel": 0.5})
        self.assertEqual(output_columns(["revel.hg19"]), ["revel"])
        self.assertEqual(output_columns([]), [])


class InTheResults(unittest.TestCase):
    """Install advanced databases into a data folder and read them back through the annotator."""

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        root = Path(cls.folder.name)
        cls.data = root / "data"
        manifest = build_demo_store(ROOT / "demo", cls.data)
        db = cls.data / DB_DIR
        genome = open_build_genome(cls.data / GENOME_DIR, "hg19")
        regions = builders.CodingRegions(db / "transcripts.hg19.sqlite")
        builders.build_dbsnp(lines_of("dbsnp.hg19.sample.txt"), db / "dbsnp.hg19.sqlite", "hg19", genome, regions)
        builders.build_gnomad([("7", lines_of("gnomad.sample.txt"))], db / "gnomad.hg19.sqlite", "hg19", genome, regions)
        builders.build_hotspots(FIXTURES / "hotspots.sample.json", db / "hotspots.sqlite")
        genome.close()
        census = root / "census.tsv"
        census.write_text("GENE_SYMBOL\tTIER\tROLE_IN_CANCER\tTUMOUR_TYPES_SOMATIC\nBRAF\t1\toncogene\tmelanoma\n", encoding="utf-8")
        builders.build_cosmic_census(census, db / "cosmic_census.sqlite")
        (cls.data / "scores").mkdir()
        (cls.data / "scores" / "phylop.hg19.bw").write_bytes(tiny_bigwig())
        for name in ("dbsnp.hg19", "gnomad.hg19", "hotspots", "cosmic_census", "phylop.hg19"):
            manifest["files"][name] = {"policy": "demo", "kind": name.split(".")[0], "source": {"version": "test"}, "files": {}}
        import json

        (cls.data / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        cls.store = DataStore(cls.data)
        cls.annotator = Annotator(cls.store)
        cls.record = cls.annotator.annotate_text("chr7-140453136-A-T", "hg19")[0]
        cls.advanced = cls.record["evidence"]["advanced"]

    @classmethod
    def tearDownClass(cls):
        cls.store.close()
        cls.folder.cleanup()

    def test_each_installed_database_answers(self):
        self.assertEqual(self.advanced["dbsnp"]["rsid"], "113488022")
        self.assertEqual(self.advanced["gnomad"]["af_popmax"], 3.26669e-05)
        self.assertEqual(self.advanced["cosmic_census"]["tier"], "1")

    def test_a_database_that_is_not_installed_adds_nothing(self):
        self.assertNotIn("revel", self.advanced)
        self.assertNotIn("omim", self.advanced)

    def test_no_match_says_not_found(self):
        record = self.annotator.annotate_text("TP53 R175H", "hg19")[0]
        self.assertEqual(record["evidence"]["advanced"]["hotspots"]["status"], "not_found")
        self.assertEqual(self.advanced["phylop"]["status"], "not_found")

    def test_a_hotspot_on_the_right_gene_and_residue(self):
        record = self.annotator.annotate_text("NRAS Q61K", "hg19")[0]
        hotspot = record["evidence"]["advanced"]["hotspots"]
        self.assertEqual(hotspot["status"], "found")
        self.assertIn("Q61", hotspot["label"])
        self.assertEqual(record["priority"]["label"], "High")

    def test_a_wide_indel_region_does_not_hide_the_single_residue_hotspot(self):
        """The real list has BRAF V600 and an indel region 592-604, and both cover residue 600."""
        record = self.annotator.annotate_text("BRAF V600E", "hg19")[0]
        self.assertEqual(record["evidence"]["advanced"]["hotspots"]["label"], "BRAF V600: 897 tumours")

    def test_a_missense_change_inside_an_indel_region_is_not_a_hotspot(self):
        record = self.annotator.annotate_text("BRAF F595L", "hg19")[0]
        self.assertEqual(record["evidence"]["advanced"]["hotspots"]["status"], "not_found")

    def test_an_in_frame_change_inside_an_indel_region_is_reported_as_such(self):
        record = self.annotator.annotate_text("BRAF T599dup", "hg19")[0]
        self.assertIn("inframe_insertion", record["effect"]["all"])
        self.assertIn("592-604 in-frame indels: 8 tumours", record["evidence"]["advanced"]["hotspots"]["label"])

    def test_a_different_residue_letter_is_not_trusted(self):
        record = self.annotator.annotate_text("NRAS Q61K", "hg19")[0]
        self.assertIn("Q61", record["evidence"]["advanced"]["hotspots"]["label"])
        store_rows = self.store.advanced_hotspots("NRAS", 61)
        self.assertEqual(store_rows[0]["residue"], "Q61")
        self.assertEqual(self.store.advanced_hotspots("NRAS", 62), [])

    def test_the_card_shows_advanced_values(self):
        from varsight.output import render_card

        self.assertIn("gnomad_popmax", render_card(self.record))

    def test_an_rs_number_works_through_dbsnp(self):
        keys = [r["variant"]["hg19"]["key"] for r in self.annotator.annotate_text("rs121913365", "hg19")]
        self.assertIn("chr7-140453132-T-A", keys)

    def test_explain_shows_the_real_queries(self):
        output = io.StringIO()
        store = DataStore(self.data)
        self.addCleanup(store.close)
        with contextlib.redirect_stdout(output):
            code = run_explain(store, "BRAF V600E", "hg19", "")
        text = output.getvalue()
        self.assertEqual(code, 0)
        for expected in ("Step 1", "Step 2", "Step 3", "clinvar.hg19.sqlite", "transcripts.hg19.sqlite", "SELECT", "genome/hg19.fa", "13961", "Nothing left this computer"):
            self.assertIn(expected, text)

    def test_benchmark_reports_every_format(self):
        report = run_benchmark(self.data, 40)
        for expected in ("TSV", "Excel", "HTML report", "JSON Lines", "One variant", "Peak memory"):
            self.assertIn(expected, report)
        self.assertNotIn("-" * 2, report)


class Compressed(unittest.TestCase):
    def test_builders_read_gzip_too(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            file = root / "census.tsv.gz"
            with gzip.open(file, "wt", encoding="utf-8") as handle:
                handle.write("GENE_SYMBOL\tTIER\tROLE_IN_CANCER\nBRAF\t1\toncogene\n")
            self.assertEqual(builders.build_cosmic_census(file, root / "c.sqlite")["rows"], 1)


if __name__ == "__main__":
    unittest.main()
