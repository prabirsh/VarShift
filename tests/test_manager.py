"""Setup, weekly update and rollback, end to end against a real HTTP server on 127.0.0.1."""

import gzip
import hashlib
import http.server
import json
import logging
import shutil
import tempfile
import threading
import unittest
from pathlib import Path

from varsight.annotate import Annotator
from varsight.errors import DataError, SetupError
from varsight.pipeline.fetch import Fetcher
from varsight.pipeline.manager import DataManager, load_catalog
from varsight.store import DataStore

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


class Lifecycle(unittest.TestCase):
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
        handler = type("Handler", (Server,), {"root": cls.www, "get_requests": []})
        cls.handler = handler
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

        catalog = load_catalog()
        catalog["allowed_hosts"] = ["127.0.0.1"]
        for source_id, source in catalog["sources"].items():
            for key in ("list_url", "pattern", "md5_url", "md5_list", "url_template"):
                source.pop(key, None)
            name = MAPPING.get(source_id, "gencode.gff3" if source_id.startswith("gencode") else "absent")
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
        self.assertEqual(self._clinvar_rows(), rows_before, "--check must change nothing")

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
        self.assertIn("--approve", items["transcripts.hg19"].detail)
        self.assertEqual(manager.load_manifest()["files"]["transcripts.hg19"]["built_at"], built_at)

        items = {i.artifact: i for i in manager.update(approve=["transcripts.hg19"])}
        self.assertEqual(items["transcripts.hg19"].action, "updated", items["transcripts.hg19"].detail)
        manager.rollback()
        path.write_bytes(original)

    def test_skipping_pfam_at_setup_is_remembered_by_updates(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        fetcher = Fetcher(["127.0.0.1"], retry_delay=0.01, allow_http_hosts=("127.0.0.1",))
        manager = DataManager(Path(folder.name), json.loads(json.dumps(self.catalog)), fetcher)
        manager._prepare_folders()
        manifest = manager.load_manifest()
        built = manager._build("proteins", {"skip_pfam_regions": True})
        manager._publish(manifest, built, None)
        self.handler.get_requests.clear()

        items = {i.artifact: i for i in manager.update()}

        self.assertEqual(items["proteins"].action, "none", items["proteins"].changed)
        self.assertNotIn("/pfam_regions.tsv.gz", self.handler.get_requests)

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


if __name__ == "__main__":
    unittest.main()
