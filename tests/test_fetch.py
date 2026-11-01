"""The downloader against a real HTTP server on 127.0.0.1, including a connection that drops."""

import gzip
import hashlib
import http.server
import tempfile
import threading
import unittest
from pathlib import Path

from varsight.errors import SetupError
from varsight.pipeline.fetch import Fetcher, md5_from_list

BODY = b"".join(b"line %06d of the source file\n" % i for i in range(60000))
GZIPPED = gzip.compress(BODY, mtime=0)


class Handler(http.server.BaseHTTPRequestHandler):
    drop_after: int | None = None
    requests_seen: list = []

    def log_message(self, *args):
        pass

    def _files(self):
        return {"/plain.txt": BODY, "/data.gz": GZIPPED, "/broken.gz": GZIPPED[: len(GZIPPED) // 2] + b"x" * 200}

    def do_HEAD(self):
        self._respond(head=True)

    def do_GET(self):
        self._respond(head=False)

    def _respond(self, head):
        type(self).requests_seen.append((self.command, self.path, self.headers.get("Range")))
        if self.path == "/redirect-bad":
            self.send_response(302)
            self.send_header("Location", "https://evil.example/file")
            self.end_headers()
            return
        if self.path == "/redirect-plain":
            self.send_response(302)
            self.send_header("Location", "http://127.0.0.1:1/file")
            self.end_headers()
            return
        body = self._files().get(self.path)
        if body is None:
            self.send_error(404)
            return
        start = 0
        requested = self.headers.get("Range")
        if requested:
            start = int(requested.split("=")[1].split("-")[0])
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(body) - start))
        self.send_header("ETag", '"v1"')
        self.send_header("Last-Modified", "Mon, 26 Oct 2026 13:51:56 GMT")
        self.end_headers()
        if head:
            return
        payload = body[start:]
        limit = type(self).drop_after
        if limit is not None:
            type(self).drop_after = None
            self.wfile.write(payload[:limit])
            self.wfile.flush()
            self.connection.close()
            return
        self.wfile.write(payload)


class Downloads(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def setUp(self):
        Handler.drop_after = None
        Handler.requests_seen = []
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.fetcher = Fetcher(["127.0.0.1"], retry_delay=0.01, allow_http_hosts=("127.0.0.1",))

    def test_a_simple_download_records_checksums(self):
        meta = self.fetcher.download(f"{self.base}/plain.txt", self.root / "a", 10**8)
        self.assertEqual((self.root / "a").read_bytes(), BODY)
        self.assertEqual(meta["sha256"], hashlib.sha256(BODY).hexdigest())
        self.assertEqual(meta["etag"], '"v1"')

    def test_a_dropped_connection_resumes_where_it_stopped(self):
        Handler.drop_after = 500_000
        self.fetcher.download(f"{self.base}/plain.txt", self.root / "b", 10**8)
        self.assertEqual((self.root / "b").read_bytes(), BODY)
        ranges = [r[2] for r in Handler.requests_seen if r[0] == "GET"]
        self.assertTrue(any(r and r != "bytes=0-" for r in ranges), ranges)

    def test_an_interrupted_file_resumes_on_the_next_run(self):
        partial = self.root / "c.part"
        partial.write_bytes(BODY[:1000])
        (self.root / "c.part.json").write_text('{"fingerprint": "\\"v1\\""}')
        self.fetcher.download(f"{self.base}/plain.txt", self.root / "c", 10**8)
        self.assertEqual((self.root / "c").read_bytes(), BODY)
        self.assertIn("bytes=1000-", [r[2] for r in Handler.requests_seen])

    def test_a_partial_file_of_an_older_version_is_discarded(self):
        partial = self.root / "d.part"
        partial.write_bytes(b"old version bytes")
        (self.root / "d.part.json").write_text('{"fingerprint": "\\"v0\\""}')
        self.fetcher.download(f"{self.base}/plain.txt", self.root / "d", 10**8)
        self.assertEqual((self.root / "d").read_bytes(), BODY)

    def test_a_wrong_md5_removes_the_file(self):
        with self.assertRaises(SetupError):
            self.fetcher.download(f"{self.base}/plain.txt", self.root / "e", 10**8, "0" * 32)
        self.assertFalse((self.root / "e").exists())
        self.assertFalse((self.root / "e.part").exists())

    def test_the_right_md5_passes(self):
        digest = hashlib.md5(BODY, usedforsecurity=False).hexdigest()
        meta = self.fetcher.download(f"{self.base}/plain.txt", self.root / "f", 10**8, digest)
        self.assertTrue(meta["md5_verified"])

    def test_a_damaged_gzip_is_caught(self):
        with self.assertRaises(SetupError):
            self.fetcher.download(f"{self.base}/broken.gz", self.root / "g", 10**8)
        self.assertFalse((self.root / "g").exists())

    def test_a_good_gzip_passes(self):
        self.fetcher.download(f"{self.base}/data.gz", self.root / "h", 10**8)
        self.assertEqual(gzip.decompress((self.root / "h").read_bytes()), BODY)

    def test_a_file_over_its_size_limit_is_refused(self):
        with self.assertRaises(SetupError):
            self.fetcher.download(f"{self.base}/plain.txt", self.root / "i", 1000)
        self.assertFalse((self.root / "i").exists())

    def test_plain_http_is_refused(self):
        fetcher = Fetcher(["127.0.0.1"], retry_delay=0.01)
        with self.assertRaises(SetupError) as caught:
            fetcher.head(f"{self.base}/plain.txt")
        self.assertIn("non-HTTPS", str(caught.exception))

    def test_a_host_off_the_list_is_refused(self):
        with self.assertRaises(SetupError) as caught:
            self.fetcher.head("https://evil.example/file")
        self.assertIn("allow-list", str(caught.exception))

    def test_a_redirect_to_a_host_off_the_list_is_refused(self):
        with self.assertRaises(SetupError):
            self.fetcher.head(f"{self.base}/redirect-bad")

    def test_a_missing_file_fails_clearly(self):
        with self.assertRaises(SetupError):
            self.fetcher.download(f"{self.base}/absent", self.root / "j", 10**8)

    def test_streaming_reads_the_whole_file_through_a_dropped_connection(self):
        Handler.drop_after = 300_000
        with self.fetcher.stream(f"{self.base}/plain.txt", 10**8) as stream:
            collected = b""
            while True:
                block = stream.read(65536)
                if not block:
                    break
                collected += block
        self.assertEqual(collected, BODY)
        self.assertEqual(stream.metadata()["sha256"], hashlib.sha256(BODY).hexdigest())

    def test_a_directory_listing_picks_the_newest_version(self):
        class Lister(Fetcher):
            def read_text(self, url, max_bytes=0):
                return '<a href="gencode.v9.basic.annotation.gff3.gz">x</a><a href="gencode.v50.basic.annotation.gff3.gz">y</a>'

        lister = Lister(["127.0.0.1"], allow_http_hosts=("127.0.0.1",))
        url = lister.resolve_listing(f"{self.base}/", r"gencode\.v\d+\.basic\.annotation\.gff3\.gz")
        self.assertTrue(url.endswith("gencode.v50.basic.annotation.gff3.gz"))

    def test_a_listing_without_a_match_fails_clearly(self):
        class Lister(Fetcher):
            def read_text(self, url, max_bytes=0):
                return "<html>moved</html>"

        with self.assertRaises(SetupError):
            Lister(["127.0.0.1"], allow_http_hosts=("127.0.0.1",)).resolve_listing(f"{self.base}/", "x")


class Md5List(unittest.TestCase):
    TEXT = "2c85c9e94fc9b52a7c66bd954b855da4  chromAgp.tar.gz\nbcdbfbe9da62f19bee88b74dabef8cd3  hg19.2bit\n"

    def test_finds_the_named_file(self):
        self.assertEqual(md5_from_list(self.TEXT, "hg19.2bit"), "bcdbfbe9da62f19bee88b74dabef8cd3")

    def test_absent_file_gives_none(self):
        self.assertIsNone(md5_from_list(self.TEXT, "hg19.fa.gz"))


if __name__ == "__main__":
    unittest.main()
