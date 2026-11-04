"""Setup, weekly update and rollback, end to end against a real HTTP server on 127.0.0.1."""

import gzip
import hashlib
import http.server
import json
import logging
import os
import shutil
import tempfile
import threading
import time
import unittest
from datetime import datetime, timezone
from email.utils import format_datetime
from pathlib import Path
from unittest import mock

from tests.helpers import tiny_bigwig
from varshift.annotate import Annotator
from varshift.errors import DataError, SetupError
from varshift.pipeline import manager as manager_module
from varshift.pipeline import updatecheck
from varshift.pipeline.fetch import Fetcher
from varshift.pipeline.manager import DataManager, load_catalog
from varshift.store import DataStore
from varshift.updates import read_status

DEMO_SOURCES = Path(__file__).resolve().parents[1] / "demo" / "sources"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
MAPPING = {
    "hg19_fasta": "hg19.windows.fa.gz",
    "hg38_fasta": "hg38.windows.fa.gz",
    "chain_hg19_hg38": "chain.hg19_hg38.gz",
    "chain_hg38_hg19": "chain.hg38_hg19.gz",
    "refseq_hg19": "refseq.hg19.gff.gz",
    "refseq_hg38": "refseq.hg38.gff.gz",
    "mane": "mane.summary.txt.gz",
    "clinvar_hg19": "clinvar.hg19.vcf.gz",
    "clinvar_hg38": "clinvar.hg38.vcf.gz",
    "clinvar_summary": "variant_summary.txt.gz",
    "hgnc": "hgnc.txt.gz",
    "clingen_validity": "clingen_validity.csv.gz",
    "clingen_dosage": "clingen_dosage.tsv.gz",
    "gnomad_constraint": "gnomad_constraint.tsv.gz",
    "hpo_diseases": "hpo_diseases.txt.gz",
    "uniprot": "uniprot.dat.gz",
    "pfam_regions": "pfam_regions.tsv.gz",
    "civic_variants": "civic_variants.tsv.gz",
    "civic_evidence": "civic_evidence.tsv.gz",
}
ADVANCED = {
    "dbsnp_hg19": "dbsnp.hg19.gz",
    "dbsnp_hg38": "dbsnp.hg38.gz",
    "hotspots": "hotspots.sample.json",
    "clingen_variants": "clingen_erepo.sample.txt",
    "phylop_hg19": "phylop.bw",
    "phylop_hg38": "phylop.bw",
}


class CoreSources(unittest.TestCase):
    def test_setup_reads_no_huge_file_off_the_network(self):
        """Reading Pfam's 4.7 GB file once made setup take most of a day, on a server giving 100 KB/s."""
        catalog = load_catalog()
        plan = DataManager(Path("unused"), catalog, Fetcher([]))._plan(["hg19", "hg38"])
        core = {source for name in plan for source in catalog["artifacts"][name]["sources"]}
        self.assertIn("pfam_regions", core)
        for source in core:
            self.assertFalse(catalog["sources"][source].get("stream"), source)
        self.assertLessEqual(catalog["sources"]["pfam_regions"]["max_mb"], 20)

    def test_revel_comes_from_zenodo_with_its_checksum(self):
        """REVEL's own server refuses scripts, and its download page names Zenodo as the other source."""
        revel = load_catalog()["sources"]["revel"]
        self.assertTrue(revel["url"].startswith("https://zenodo.org/api/records/7072866/"))
        self.assertRegex(revel["md5"], r"^[0-9a-f]{32}$")


class Reporting(unittest.TestCase):
    def test_a_long_build_keeps_saying_it_is_going(self):
        with mock.patch.object(manager_module, "STILL_GOING_EVERY", 0.01), self.assertLogs("varshift.data", "INFO") as logs:
            with manager_module._still_going("building something"):
                time.sleep(0.1)
        self.assertIn("still building something", logs.output[0])


class Server(http.server.BaseHTTPRequestHandler):
    root: Path
    get_requests: list = []

    def log_message(self, *args):
        pass

    def do_HEAD(self):
        self._serve(False)

    def do_GET(self):
        type(self).get_requests.append(self.path)
        self._serve(True)

    def _serve(self, body):
        path = self.root / self.path.lstrip("/")
        if not path.is_file():
            self.send_error(404)
            return
        data = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Length", str(len(data)))
        # The content hash stands in for the ETag a real server sends when a file changes.
        self.send_header("ETag", f'"{hashlib.sha1(data, usedforsecurity=False).hexdigest()}"')
        self.send_header("Last-Modified", "Mon, 26 Oct 2026 13:51:56 GMT")
        self.end_headers()
        if body:
            self.wfile.write(data)


class WithServer(unittest.TestCase):
    """Serves the demo sources, and a few samples for the advanced databases, from 127.0.0.1."""

    @classmethod
    def setUpClass(cls):
        logging.disable(logging.CRITICAL)
        cls.folder = tempfile.TemporaryDirectory()
        root = Path(cls.folder.name)
        cls.www = root / "www"
        cls.www.mkdir()
        for name in set(MAPPING.values()):
            shutil.copy(DEMO_SOURCES / name, cls.www / name)
        shutil.copy(FIXTURES / "gencode.sample.gff3", cls.www / "gencode.gff3")
        for name in ("hotspots.sample.json", "clingen_erepo.sample.txt"):
            shutil.copy(FIXTURES / name, cls.www / name)
        for build in ("hg19", "hg38"):
            with open(FIXTURES / f"dbsnp.{build}.sample.txt", "rb") as source, gzip.open(cls.www / f"dbsnp.{build}.gz", "wb") as target:
                shutil.copyfileobj(source, target)
        (cls.www / "phylop.bw").write_bytes(tiny_bigwig())
        handler = type("Handler", (Server,), {"root": cls.www, "get_requests": []})
        cls.handler = handler
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

        catalog = load_catalog()
        catalog["allowed_hosts"] = ["127.0.0.1"]
        for source_id, source in catalog["sources"].items():
            for key in ("list_url", "pattern", "md5_url", "md5_list", "md5", "url_template"):
                source.pop(key, None)
            name = {**MAPPING, **ADVANCED}.get(source_id, "gencode.gff3" if source_id.startswith("gencode") else "absent")
            source["url"] = f"{cls.base}/{name}"
        cls.catalog = catalog
        cls.data = root / "data"
        cls.manager = cls._make_manager()
        cls.installed = cls.manager.setup(["hg19", "hg38"])

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.folder.cleanup()
        logging.disable(logging.NOTSET)

    @classmethod
    def _make_manager(cls):
        fetcher = Fetcher(["127.0.0.1"], retry_delay=0.01, allow_http_hosts=("127.0.0.1",))
        return DataManager(cls.data, json.loads(json.dumps(cls.catalog)), fetcher)


class Lifecycle(WithServer):
    def _clinvar_rows(self, build="hg19"):
        store = DataStore(self.data)
        try:
            return store._db(f"clinvar.{build}").execute("SELECT count(*) FROM clinvar").fetchone()[0]
        finally:
            store.close()

    def _rewrite_clinvar(self, build, keep=None, append=""):
        """Publish a changed ClinVar file on the test server."""
        path = self.www / f"clinvar.{build}.vcf.gz"
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            lines = handle.readlines()
        header = [line for line in lines if line.startswith("#")]
        records = [line for line in lines if not line.startswith("#")]
        if keep is not None:
            records = records[: int(len(records) * keep)]
        with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
            handle.writelines(header + records)
            handle.write(append)

    def test_1_setup_installs_everything_and_answers_correctly(self):
        names = {item["artifact"] for item in self.installed}
        for expected in ("genome.hg19", "genome.hg38", "transcripts.hg19", "transcripts.hg38", "clinvar.hg19", "clinvar.hg38", "genes", "proteins", "civic", "liftover.hg19_hg38", "liftover.hg38_hg19"):
            self.assertIn(expected, names)
        manifest = self.manager.load_manifest()
        self.assertEqual(sorted(manifest["builds"]), ["hg19", "hg38"])
        store = DataStore(self.data)
        try:
            record = Annotator(store).annotate_text("BRAF V600E", "hg19")[0]
        finally:
            store.close()
        self.assertEqual(record["variant"]["hg38"]["key"], "chr7-140753336-A-T")
        self.assertEqual(record["evidence"]["clinvar"]["id"], 13961)

    def test_2_downloads_are_removed_after_the_build(self):
        leftovers = [p.name for p in (self.data / "downloads").iterdir()]
        self.assertEqual(leftovers, [])

    def test_3_setup_again_rebuilds_nothing(self):
        self.assertEqual(self._make_manager().setup(["hg19", "hg38"]), [])

    def test_4_checksums_recorded_and_verified(self):
        results = self.manager.verify()
        self.assertTrue(results)
        self.assertTrue(all(ok for _, ok, _ in results), results)

    def test_5_nothing_new_means_no_change(self):
        before = self.manager.load_manifest()["files"]["clinvar.hg19"]["built_at"]
        items = self._make_manager().update()
        self.assertTrue(all(i.action == "none" for i in items), [(i.artifact, i.action, i.changed) for i in items])
        self.assertEqual(self.manager.load_manifest()["files"]["clinvar.hg19"]["built_at"], before)

    def test_5b_a_file_made_on_request_is_not_new_every_time(self):
        # ClinGen sends no ETag and dates its gene-validity file with the time of the request.
        manager = self._make_manager()
        real_head = manager.fetcher.head

        def head(url):
            remote = real_head(url)
            if url.endswith("/clingen_validity.csv.gz"):
                remote.etag, remote.last_modified = "", format_datetime(datetime.now(timezone.utc), usegmt=True)
            return remote

        with mock.patch.object(manager.fetcher, "head", side_effect=head):
            items = {i.artifact: i for i in manager.check_updates()}
        self.assertEqual(items["genes"].changed, [])
        self.assertFalse(manager_module._just_now("Mon, 26 Oct 2026 13:51:56 GMT"))
        self.assertFalse(manager_module._just_now(""))

    def test_6_new_clinvar_is_applied_automatically_then_can_be_rolled_back(self):
        rows_before = self._clinvar_rows()
        store = DataStore(self.data)
        ref = store.genome("hg19").fetch("chr7", 140499999, 140500000)
        store.close()
        alt = "A" if ref != "A" else "C"
        extra = f"7\t140500000\t999999999\t{ref}\t{alt}\t.\t.\tALLELEID=1;CLNSIG=Benign;GENEINFO=BRAF:673\n"
        self._rewrite_clinvar("hg19", append=extra)
        manager = self._make_manager()

        items = {i.artifact: i for i in manager.update(check_only=True)}
        self.assertEqual(items["clinvar.hg19"].action, "available")
        self.assertEqual(self._clinvar_rows(), rows_before, "a check must change nothing")

        items = {i.artifact: i for i in manager.update()}
        self.assertEqual(items["clinvar.hg19"].action, "updated", items["clinvar.hg19"].detail)
        self.assertEqual(self._clinvar_rows(), rows_before + 1)
        self.assertTrue(manager.list_snapshots())

        restored = manager.rollback()
        self.assertIn("clinvar.hg19", restored)
        self.assertEqual(self._clinvar_rows(), rows_before)

        # Serve the file again without the extra row, so later tests start from a clean state.
        path = self.www / "clinvar.hg19.vcf.gz"
        shutil.copy(DEMO_SOURCES / "clinvar.hg19.vcf.gz", path)

    def test_7_a_shrunken_file_is_refused_and_current_data_stays(self):
        rows_before = self._clinvar_rows("hg38")
        self._rewrite_clinvar("hg38", keep=0.5)
        manager = self._make_manager()
        items = {i.artifact: i for i in manager.update()}
        self.assertEqual(items["clinvar.hg38"].action, "failed")
        self.assertIn("shrink", items["clinvar.hg38"].detail)
        self.assertEqual(self._clinvar_rows("hg38"), rows_before)
        shutil.copy(DEMO_SOURCES / "clinvar.hg38.vcf.gz", self.www / "clinvar.hg38.vcf.gz")

    def test_8_a_file_that_breaks_the_self_test_is_refused(self):
        path = self.www / "clinvar.hg19.vcf.gz"
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            text = handle.read()
        self.assertIn("\t140453136\t13961\t", text)
        with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(text.replace("\t140453136\t13961\t", "\t140453136\t99999\t"))
        manager = self._make_manager()
        items = {i.artifact: i for i in manager.update()}
        self.assertEqual(items["clinvar.hg19"].action, "failed")
        self.assertIn("self-test", items["clinvar.hg19"].detail)
        store = DataStore(self.data)
        try:
            self.assertIsNotNone(store.clinvar_row("hg19", "chr7", 140453136, "A", "T"))
        finally:
            store.close()
        shutil.copy(DEMO_SOURCES / "clinvar.hg19.vcf.gz", path)

    def test_9_transcript_changes_wait_for_approval(self):
        path = self.www / "refseq.hg19.gff.gz"
        original = path.read_bytes()
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            text = handle.read()
        with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(text + "# a new release\n")
        manager = self._make_manager()
        built_at = manager.load_manifest()["files"]["transcripts.hg19"]["built_at"]

        items = {i.artifact: i for i in manager.update()}
        self.assertEqual(items["transcripts.hg19"].action, "needs_approval")
        self.assertIn("approve=", items["transcripts.hg19"].detail)
        self.assertEqual(manager.load_manifest()["files"]["transcripts.hg19"]["built_at"], built_at)

        items = {i.artifact: i for i in manager.update(approve=["transcripts.hg19"])}
        self.assertEqual(items["transcripts.hg19"].action, "updated", items["transcripts.hg19"].detail)
        manager.rollback()
        path.write_bytes(original)

    def test_a_changed_source_that_cannot_be_read_is_reported(self):
        manager = self._make_manager()
        manager.catalog["sources"]["civic_variants"]["url"] = f"{self.base}/gone.tsv"
        items = {i.artifact: i for i in manager.update()}
        self.assertEqual(items["civic"].action, "error")

    def test_a_tampered_file_is_detected(self):
        target = self.data / "db" / "civic.sqlite"
        original = target.read_bytes()
        target.chmod(0o644)
        try:
            target.write_bytes(original + b"tampered")
            bad = [name for name, ok, _ in self.manager.verify() if not ok]
            self.assertEqual(bad, ["civic.sqlite"])
        finally:
            target.write_bytes(original)

    def test_rollback_that_fails_part_way_changes_nothing(self):
        store = DataStore(self.data)
        ref = store.genome("hg19").fetch("chr7", 140500099, 140500100)
        store.close()
        alt = "A" if ref != "A" else "C"
        self._rewrite_clinvar("hg19", append=f"7\t140500100\t999999998\t{ref}\t{alt}\t.\t.\tALLELEID=2;CLNSIG=Benign;GENEINFO=BRAF:673\n")
        manager = self._make_manager()
        try:
            items = {i.artifact: i for i in manager.update()}
            self.assertEqual(items["clinvar.hg19"].action, "updated", items["clinvar.hg19"].detail)
            rows_now = self._clinvar_rows()
            real_replace, calls = os.replace, []

            def in_use_on_the_second_move(source, target):
                calls.append(source)
                if len(calls) == 2:
                    raise PermissionError("the file is in use by another process")
                return real_replace(source, target)

            with mock.patch("varshift.pipeline.manager.os.replace", side_effect=in_use_on_the_second_move):
                with self.assertRaises(DataError):
                    manager.rollback()
            self.assertEqual(self._clinvar_rows(), rows_now, "the newer data must still be in place")
            self.assertTrue(all(ok for _, ok, _ in manager.verify()), "files and manifest must still agree")
            self.assertIn("clinvar.hg19", manager.rollback(), "the snapshot must still be usable")
            self.assertEqual(self._clinvar_rows(), rows_now - 1)
        finally:
            shutil.copy(DEMO_SOURCES / "clinvar.hg19.vcf.gz", self.www / "clinvar.hg19.vcf.gz")

    def test_rollback_with_nothing_saved_says_so(self):
        empty = tempfile.TemporaryDirectory()
        self.addCleanup(empty.cleanup)
        with self.assertRaises(DataError):
            DataManager(Path(empty.name), self.catalog, Fetcher(["127.0.0.1"])).rollback()

    def test_a_second_run_at_once_is_refused(self):
        lock = self.data / ".lock"
        lock.write_text(str(__import__("os").getpid()))
        try:
            with self.assertRaises(SetupError):
                self._make_manager().update()
        finally:
            lock.unlink()

    def test_a_stale_lock_from_a_dead_process_is_cleared(self):
        lock = self.data / ".lock"
        lock.write_text("999999999")
        self._make_manager().update(check_only=True)
        self.assertFalse(lock.exists())


class AdvancedDatabases(WithServer):
    def _record(self, text="BRAF V600E"):
        store = DataStore(self.data)
        try:
            return store, Annotator(store).annotate_text(text, "hg19")[0]
        finally:
            store.close()

    def _installed(self):
        return set(self.manager.load_manifest()["files"])

    def test_a_small_database_is_added_used_and_removed(self):
        done = self.manager.add("hotspots")
        self.assertEqual([item["artifact"] for item in done], ["hotspots"])
        self.assertTrue((self.data / "db" / "hotspots.sqlite").is_file())
        self.assertEqual(list((self.data / "downloads").iterdir()), [])
        _, record = self._record("NRAS Q61K")
        self.assertEqual(record["evidence"]["advanced"]["hotspots"]["status"], "found")

        self.assertEqual(self.manager.remove("hotspots"), ["hotspots"])
        self.assertFalse((self.data / "db" / "hotspots.sqlite").exists())
        _, record = self._record("NRAS Q61K")
        self.assertNotIn("hotspots", record["evidence"]["advanced"])

    def test_a_large_file_is_read_as_it_arrives_for_each_build(self):
        done = self.manager.add("dbsnp")
        self.assertEqual(sorted(item["artifact"] for item in done), ["dbsnp.hg19", "dbsnp.hg38"])
        self.addCleanup(self.manager.remove, "dbsnp")
        _, record = self._record()
        self.assertEqual(record["evidence"]["advanced"]["dbsnp"]["rsid"], "113488022")
        self.assertFalse([p for p in (self.data / "downloads").iterdir()])

    def test_conservation_scores_are_kept_as_the_original_file(self):
        self.manager.add("phylop")
        self.addCleanup(self.manager.remove, "phylop")
        self.assertTrue((self.data / "scores" / "phylop.hg19.bw").is_file())
        store = DataStore(self.data)
        try:
            self.assertEqual(store.phylop("hg19", "chr7", 101), 2.5)
        finally:
            store.close()

    def test_terms_must_be_accepted_first(self):
        with self.assertRaises(SetupError) as caught:
            self.manager.add("revel")
        self.assertIn("accept", str(caught.exception))

    def test_a_checksum_kept_in_the_catalog_is_checked(self):
        manager = self._make_manager()
        manager.catalog["sources"]["hotspots"]["md5"] = "0" * 32
        with self.assertRaises(SetupError) as caught:
            manager.add("hotspots")
        self.assertIn("MD5", str(caught.exception))
        self.assertFalse(any((self.data / "downloads").glob("hotspots__*")))

    def test_a_licensed_database_needs_the_users_own_file(self):
        with self.assertRaises(SetupError) as caught:
            self.manager.add("omim", {"accepted": True})
        self.assertIn("file=PATH", str(caught.exception))
        with self.assertRaises(SetupError) as caught:
            self.manager.add("omim", {"accepted": True, "file": str(self.data / "absent.txt")})
        self.assertIn("cannot find", str(caught.exception))

    def test_an_own_file_is_installed_and_left_alone_by_updates(self):
        file = Path(self.folder.name) / "census.tsv"
        file.write_text("GENE_SYMBOL\tTIER\tROLE_IN_CANCER\nBRAF\t1\toncogene\n", encoding="utf-8")
        self.manager.add("cosmic_census", {"accepted": True, "file": str(file)})
        self.addCleanup(self.manager.remove, "cosmic_census")
        self.assertIn("cosmic_census", self._installed())
        self.assertNotIn("cosmic_census", [item.artifact for item in self._make_manager().update()])
        _, record = self._record()
        self.assertEqual(record["evidence"]["advanced"]["cosmic_census"]["tier"], "1")

    def test_a_per_build_licensed_file_must_say_its_build(self):
        file = Path(self.folder.name) / "cosmic.tsv"
        file.write_text("GENOMIC_MUTATION_ID\tHGVSG\nCOSV1\t7:g.140453136A>T\n", encoding="utf-8")
        with self.assertRaises(SetupError) as caught:
            self.manager.add("cosmic", {"accepted": True, "file": str(file)})
        self.assertIn("build=hg19", str(caught.exception))
        self.manager.add("cosmic", {"accepted": True, "file": str(file), "build": "hg19"})
        self.addCleanup(self.manager.remove, "cosmic")
        self.assertIn("cosmic.hg19", self._installed())
        self.assertNotIn("cosmic.hg38", self._installed())

    def test_new_public_database_version_is_noticed(self):
        self.manager.add("clingen_variants")
        self.addCleanup(self.manager.remove, "clingen_variants")
        path = self.www / "clingen_erepo.sample.txt"
        original = path.read_text(encoding="utf-8")
        self.addCleanup(path.write_text, original, encoding="utf-8")
        self.assertEqual({i.artifact: i.changed for i in self.manager.check_updates()}["clingen_variants"], [])
        path.write_text(original + "NM_0:c.1A>G\t999\tCA1\tx\tPAH\tphenylketonuria\t\t\tPathogenic\n", encoding="utf-8")
        self.assertEqual({i.artifact: i.changed for i in self.manager.check_updates()}["clingen_variants"], ["clingen_variants"])

    def test_the_background_check_saves_its_answer_for_next_time(self):
        self.manager.add("clingen_variants")
        self.addCleanup(self.manager.remove, "clingen_variants")
        path = self.www / "clingen_erepo.sample.txt"
        original = path.read_text(encoding="utf-8")
        self.addCleanup(path.write_text, original, encoding="utf-8")
        path.write_text(original + "NM_0:c.1A>G\t999\tCA1\tx\tPAH\tphenylketonuria\t\t\tPathogenic\n", encoding="utf-8")
        with mock.patch.object(updatecheck, "load_catalog", return_value=json.loads(json.dumps(self.catalog))), mock.patch.object(updatecheck, "Fetcher", lambda hosts, **options: self._make_manager().fetcher):
            self.assertEqual(updatecheck.run(self.data), 0)
        waiting = [item["artifact"] for item in read_status(self.data)["items"]]
        self.assertIn("clingen_variants", waiting)
        self.assertEqual(self.manager.load_manifest()["files"]["clingen_variants"]["stats"]["rows"], 3, "a check changes nothing")

    def test_removing_what_is_not_installed_says_so(self):
        with self.assertRaises(SetupError) as caught:
            self.manager.remove("gnomad")
        self.assertIn("not installed", str(caught.exception))

    def test_core_databases_cannot_be_removed_by_this_command(self):
        with self.assertRaises(SetupError):
            self.manager.remove("clinvar")
        self.assertIn("clinvar.hg19", self._installed())

    def test_topmed_points_to_dbsnp(self):
        with self.assertRaises(SetupError) as caught:
            self.manager.add("topmed")
        self.assertIn("dbsnp", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
