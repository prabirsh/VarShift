"""Careful file transfer: HTTPS only, a host allow-list, resume, retries and integrity checks."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import socket
import ssl
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urljoin, urlsplit

from .. import __version__
from ..errors import SetupError

log = logging.getLogger("varsight.fetch")

CHUNK = 1024 * 1024
ATTEMPTS = 6
TIMEOUT = 60
_USER_AGENT = f"VarSight/{__version__} (reference data setup)"
_MD5_RE = re.compile(r"^[0-9a-f]{32}$")


@dataclass
class RemoteInfo:
    url: str
    size: int | None
    etag: str
    last_modified: str

    def fingerprint(self) -> str:
        return self.etag or self.last_modified or (str(self.size) if self.size else "")


class Fetcher:
    def __init__(
        self,
        allowed_hosts: list[str],
        retry_delay: float = 2.0,
        allow_http_hosts: tuple[str, ...] = (),
        attempts: int = ATTEMPTS,
        timeout: int = TIMEOUT,
    ):
        self.allowed_hosts = set(allowed_hosts)
        self.retry_delay = retry_delay
        self.attempts = attempts
        self.timeout = timeout
        # Only the test suite passes this, for a server it runs on 127.0.0.1.
        self.allow_http_hosts = set(allow_http_hosts)
        self._context = ssl.create_default_context()
        self._opener = urllib.request.build_opener(
            _CheckedRedirect(self._check_url),
            urllib.request.HTTPSHandler(context=self._context),
        )

    def _check_url(self, url: str) -> None:
        parts = urlsplit(url)
        host = parts.hostname or ""
        if parts.scheme == "http" and host in self.allow_http_hosts:
            return
        if parts.scheme != "https":
            raise SetupError(f"refusing a non-HTTPS address: {url}")
        if host not in self.allowed_hosts:
            raise SetupError(f"refusing a host that is not on the allow-list: {host}")

    def _request(self, url: str, headers: dict | None = None, method: str = "GET"):
        self._check_url(url)
        request = urllib.request.Request(url, method=method, headers={"User-Agent": _USER_AGENT, **(headers or {})})
        return self._opener.open(request, timeout=self.timeout)

    def head(self, url: str) -> RemoteInfo:
        for attempt in range(1, self.attempts + 1):
            try:
                with self._request(url, method="HEAD") as response:
                    return _remote_info(url, response)
            except urllib.error.HTTPError as exc:
                if exc.code in (403, 405, 501):
                    return self._head_by_range(url)
                if attempt == self.attempts:
                    raise SetupError(f"cannot reach {url}: {exc}") from exc
                self._wait(attempt)
            except (urllib.error.URLError, socket.timeout, ConnectionError, ssl.SSLError) as exc:
                if attempt == self.attempts:
                    raise SetupError(f"cannot reach {url}: {exc}") from exc
                self._wait(attempt)
        raise SetupError(f"cannot reach {url}")

    def _head_by_range(self, url: str) -> RemoteInfo:
        """Some servers refuse HEAD. Ask for one byte instead, and read the same headers."""
        try:
            with self._request(url, headers={"Range": "bytes=0-0"}) as response:
                info = _remote_info(url, response)
                total = response.headers.get("Content-Range", "").rpartition("/")[2]
                info.size = int(total) if total.isdigit() else None
                return info
        except (urllib.error.URLError, socket.timeout, ConnectionError, ssl.SSLError) as exc:
            raise SetupError(f"cannot reach {url}: {exc}") from exc

    def read_text(self, url: str, max_bytes: int = 4 * 1024 * 1024) -> str:
        for attempt in range(1, self.attempts + 1):
            try:
                with self._request(url) as response:
                    body = response.read(max_bytes + 1)
                if len(body) > max_bytes:
                    raise SetupError(f"{url} is larger than expected")
                return body.decode("utf-8", "replace")
            except (urllib.error.URLError, socket.timeout, ConnectionError, ssl.SSLError) as exc:
                if attempt == self.attempts:
                    raise SetupError(f"cannot read {url}: {exc}") from exc
                self._wait(attempt)
        raise SetupError(f"cannot read {url}")

    def resolve_listing(self, list_url: str, pattern: str) -> str:
        """Find the file in a directory listing that matches pattern; the newest name wins."""
        listing = self.read_text(list_url, 8 * 1024 * 1024)
        names = sorted(set(re.findall(rf'href="({pattern})"', listing)), key=_version_key)
        if not names:
            raise SetupError(f"no file matching {pattern} at {list_url}; the source layout may have changed")
        return urljoin(list_url, names[-1])

    def download(self, url: str, dest: Path, max_bytes: int, expected_md5: str | None = None) -> dict:
        """Fetch url to dest, resuming a partial file. Returns what was fetched and its checksums."""
        partial = dest.with_name(dest.name + ".part")
        meta_path = dest.with_name(dest.name + ".part.json")
        remote = self.head(url)
        if remote.size is not None and remote.size > max_bytes:
            raise SetupError(f"{url} is {remote.size:,} bytes, over the {max_bytes:,} byte limit set for it")
        _check_free_space(dest.parent, remote.size or 0)

        if partial.exists() and _load_json(meta_path).get("fingerprint") != remote.fingerprint():
            partial.unlink()  # the file changed on the server; a resume would mix two versions
        meta_path.write_text(json.dumps({"fingerprint": remote.fingerprint(), "url": url}), encoding="utf-8")

        for attempt in range(1, self.attempts + 1):
            have = partial.stat().st_size if partial.exists() else 0
            if remote.size is not None and have >= remote.size:
                break
            headers = {"Range": f"bytes={have}-"} if have else {}
            try:
                with self._request(url, headers=headers) as response:
                    if have and response.status != 206:
                        have = 0  # the server ignored the range request, so start over
                    mode = "ab" if have else "wb"
                    written = have
                    with open(partial, mode) as out:
                        while True:
                            block = response.read(CHUNK)
                            if not block:
                                break
                            written += len(block)
                            if written > max_bytes:
                                raise SetupError(f"{url} grew past the {max_bytes:,} byte limit set for it")
                            out.write(block)
                if remote.size is None or written >= remote.size:
                    break
            except (urllib.error.URLError, socket.timeout, ConnectionError, ssl.SSLError, OSError) as exc:
                if isinstance(exc, urllib.error.HTTPError) and exc.code == 416:
                    partial.unlink(missing_ok=True)
                if attempt == self.attempts:
                    raise SetupError(f"download of {url} failed after {self.attempts} tries: {exc}") from exc
                log.warning("download interrupted (%s); retrying", exc)
                self._wait(attempt)

        size = partial.stat().st_size
        if remote.size is not None and size != remote.size:
            raise SetupError(f"{url}: expected {remote.size:,} bytes, got {size:,}")
        md5, sha256 = _checksums(partial)
        if expected_md5 and md5 != expected_md5:
            partial.unlink()
            meta_path.unlink(missing_ok=True)
            raise SetupError(f"{url} failed its MD5 check (expected {expected_md5}, got {md5}); the file was removed")
        if _is_gzip(partial):
            _test_gzip(partial)
        os.replace(partial, dest)
        meta_path.unlink(missing_ok=True)
        return {
            "url": url,
            "size": size,
            "etag": remote.etag,
            "last_modified": remote.last_modified,
            "md5": md5,
            "sha256": sha256,
            "md5_verified": bool(expected_md5),
        }

    def stream(self, url: str, max_bytes: int) -> "ResumingStream":
        """A readable stream for a file too large to keep; it reconnects where it left off."""
        remote = self.head(url)
        if remote.size is not None and remote.size > max_bytes:
            raise SetupError(f"{url} is {remote.size:,} bytes, over the {max_bytes:,} byte limit set for it")
        return ResumingStream(self, remote)

    def _wait(self, attempt: int) -> None:
        time.sleep(self.retry_delay * 2 ** (attempt - 1))


class ResumingStream(io.RawIOBase):
    def __init__(self, fetcher: Fetcher, remote: RemoteInfo):
        super().__init__()
        self._fetcher = fetcher
        self.remote = remote
        self._offset = 0
        self._response = None
        self._failures = 0
        self._digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def _connect(self) -> None:
        headers = {"Range": f"bytes={self._offset}-"} if self._offset else {}
        if self._offset and self.remote.etag:
            headers["If-Range"] = self.remote.etag
        self._response = self._fetcher._request(self.remote.url, headers=headers)
        if self._offset and self._response.status != 206:
            raise SetupError(f"{self.remote.url} changed during the transfer; run the step again")

    def readinto(self, buffer) -> int:
        while True:
            try:
                if self._response is None:
                    self._connect()
                count = self._response.readinto(buffer)
                self._failures = 0
            except (urllib.error.URLError, socket.timeout, ConnectionError, ssl.SSLError, OSError) as exc:
                self._failures += 1
                self._close_response()
                if self._failures >= self._fetcher.attempts:
                    raise SetupError(f"transfer of {self.remote.url} failed: {exc}") from exc
                log.warning("transfer interrupted at %s bytes (%s); reconnecting", f"{self._offset:,}", exc)
                self._fetcher._wait(self._failures)
                continue
            if count == 0 and self.remote.size is not None and self._offset < self.remote.size:
                self._failures += 1
                self._close_response()
                if self._failures >= self._fetcher.attempts:
                    raise SetupError(f"{self.remote.url} ended early at {self._offset:,} bytes")
                self._fetcher._wait(self._failures)
                continue
            self._offset += count
            self._digest.update(memoryview(buffer)[:count])
            return count

    def _close_response(self) -> None:
        if self._response is not None:
            try:
                self._response.close()
            except OSError:
                pass
            self._response = None

    def close(self) -> None:
        self._close_response()
        super().close()

    def metadata(self) -> dict:
        return {
            "url": self.remote.url,
            "size": self._offset,
            "etag": self.remote.etag,
            "last_modified": self.remote.last_modified,
            "sha256": self._digest.hexdigest(),
            "streamed": True,
        }


class _CheckedRedirect(urllib.request.HTTPRedirectHandler):
    """Hold a redirect to the same rules as the first address."""

    def __init__(self, check):
        self._check = check

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        self._check(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _remote_info(url: str, response) -> RemoteInfo:
    length = response.headers.get("Content-Length")
    return RemoteInfo(
        url,
        int(length) if length else None,
        response.headers.get("ETag", "") or "",
        response.headers.get("Last-Modified", "") or "",
    )


def md5_from_list(text: str, file_name: str) -> str | None:
    """Pull one file's MD5 from the 'hash  name' lists that UCSC and NCBI publish."""
    for line in text.splitlines():
        fields = line.split()
        if len(fields) >= 2 and _MD5_RE.match(fields[0].lower()):
            if Path(fields[-1].lstrip("*")).name == file_name:
                return fields[0].lower()
    return None


def _version_key(name: str) -> list:
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def _check_free_space(folder: Path, needed: int) -> None:
    free = shutil.disk_usage(folder).free
    margin = max(needed // 2, 16 * 1024**2)
    if free < needed + margin:
        raise SetupError(
            f"not enough disk space in {folder}: {free / 1024**3:.1f} GB free, about "
            f"{(needed + margin) / 1024**3:.1f} GB needed for the next file"
        )


def _checksums(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(CHUNK):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _test_gzip(path: Path) -> None:
    """Read the whole archive once, which makes gzip verify its stored CRC."""
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(CHUNK):
                pass
    except (OSError, EOFError) as exc:
        path.unlink(missing_ok=True)
        raise SetupError(f"{path.name} is a damaged gzip file ({exc}); it was removed, run the step again") from exc


def _load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
