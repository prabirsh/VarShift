"""First setup, weekly update and rollback of the local data.

The rule for every change: build beside the current data, test the new file, and only then switch.
A failed download, build or test leaves what is already installed untouched.
"""

from __future__ import annotations

import io
import json
import logging
import os
import shutil
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from importlib import resources
from pathlib import Path
from urllib.parse import urlsplit

from .. import __version__
from ..builds import BUILDS, primary_contigs
from ..errors import DataError, SetupError
from ..genome import Genome, open_build_genome
from ..store import DB_DIR, GENOME_DIR, MANIFEST, SCHEMA_VERSION, DataStore
from . import builders
from .fetch import CHUNK, Fetcher, md5_from_list

log = logging.getLogger("varsight.data")

SNAPSHOT_DIR = "snapshots"
DOWNLOAD_DIR = "downloads"
STAGING_DIR = "staging"
REPORT_DIR = "reports"
CACHE_DIR = "cache"
KEEP_SNAPSHOTS = 2
# A new database may not shrink past this share of the one it replaces.
SHRINK_LIMIT = 0.9
_GATE_METRIC = {"clinvar": "rows", "transcripts": "refseq", "genes": "genes", "proteins": "reviewed_proteins", "civic": "variants"}


def load_catalog() -> dict:
    text = resources.files("varsight").joinpath("data/sources.json").read_text(encoding="utf-8")
    return json.loads(text)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class UpdateItem:
    artifact: str
    policy: str
    changed: list[str] = field(default_factory=list)
    error: str = ""
    action: str = "none"
    detail: str = ""


class _Lock:
    """One setup or update at a time per data folder."""

    def __init__(self, path: Path):
        self.path = path

    def __enter__(self):
        try:
            descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError:
            if self._stale():
                self.path.unlink(missing_ok=True)
                descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
            else:
                raise SetupError(f"another setup or update is running in this data folder ({self.path})") from None
        with os.fdopen(descriptor, "w") as handle:
            handle.write(str(os.getpid()))
        return self

    def _stale(self) -> bool:
        try:
            pid = int(self.path.read_text().strip())
        except (OSError, ValueError):
            return True
        return not _process_alive(pid) or time.time() - self.path.stat().st_mtime > 3 * 86400

    def __exit__(self, *exc):
        self.path.unlink(missing_ok=True)


def _process_alive(pid: int) -> bool:
    if os.name == "nt":
        # On Windows os.kill(pid, 0) would send a console Ctrl+C rather than test the process.
        import ctypes

        query_limited_information = 0x1000
        handle = ctypes.windll.kernel32.OpenProcess(query_limited_information, False, pid)
        if not handle:
            return False
        ctypes.windll.kernel32.CloseHandle(handle)
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True
    return True


class DataManager:
    def __init__(self, data_dir: Path, catalog: dict | None = None, fetcher: Fetcher | None = None):
        self.data_dir = Path(data_dir).expanduser().resolve()
        self.catalog = catalog or load_catalog()
        self.fetcher = fetcher or Fetcher(self.catalog["allowed_hosts"])
        self.manifest_path = self.data_dir / MANIFEST

    # manifest
    def load_manifest(self) -> dict:
        if self.manifest_path.is_file():
            with open(self.manifest_path, encoding="utf-8") as handle:
                return json.load(handle)
        return {"schema": SCHEMA_VERSION, "builds": [], "created": utc_now(), "files": {}}

    def _save_manifest(self, manifest: dict) -> None:
        manifest["updated"] = utc_now()
        manifest["varsight_version"] = __version__
        temp = self.manifest_path.with_suffix(".json.tmp")
        with open(temp, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temp, self.manifest_path)

    def _prepare_folders(self) -> None:
        self.data_dir.mkdir(parents=True, exist_ok=True)
        try:
            self.data_dir.chmod(0o700)
        except OSError:
            pass
        for name in (DB_DIR, GENOME_DIR, SNAPSHOT_DIR, DOWNLOAD_DIR, STAGING_DIR, REPORT_DIR, CACHE_DIR):
            (self.data_dir / name).mkdir(exist_ok=True)

    # planning
    def _plan(self, builds: list[str], optional: list[str]) -> list[str]:
        order = ["genes", "civic", "proteins"]
        for build in builds:
            order.insert(0, f"genome.{build}")
            order.append(f"transcripts.{build}")
            order.append(f"clinvar.{build}")
        if set(builds) >= set(BUILDS):
            order.extend(["liftover.hg19_hg38", "liftover.hg38_hg19"])
        for item in optional:
            matched = [name for name in self.catalog["artifacts"] if name == item or name.startswith(f"{item}.")]
            if not matched:
                raise SetupError(f"unknown optional source {item}")
            for name in matched:
                artifact = self.catalog["artifacts"][name]
                if artifact.get("build") in builds:
                    order.append(name)
        return order

    def _artifact_paths(self, name: str) -> list[Path]:
        kind = self.catalog["artifacts"][name]["kind"]
        if kind == "genome":
            build = self.catalog["artifacts"][name]["build"]
            fasta = self.data_dir / GENOME_DIR / f"{build}.fa"
            return [fasta, Path(str(fasta) + ".fai")]
        return [self.data_dir / DB_DIR / f"{name}.sqlite"]

    def _source_url(self, source_id: str, chrom: str = "") -> str:
        source = self.catalog["sources"][source_id]
        if "url" in source:
            return source["url"]
        if "list_url" in source:
            return self.fetcher.resolve_listing(source["list_url"], source["pattern"])
        return source["url_template"].format(chrom=chrom)

    # downloads
    def _download(self, source_id: str) -> tuple[Path, dict]:
        source = self.catalog["sources"][source_id]
        url = self._source_url(source_id)
        file_name = Path(urlsplit(url).path).name or source_id
        dest = self.data_dir / DOWNLOAD_DIR / f"{source_id}__{file_name}"
        meta_path = dest.with_name(dest.name + ".meta.json")
        remote = self.fetcher.head(url)
        if dest.is_file() and meta_path.is_file():
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if (meta.get("etag"), meta.get("last_modified")) == (remote.etag, remote.last_modified):
                log.info("  already downloaded: %s", source["name"])
                return dest, meta
        expected = None
        if "md5_url" in source:
            fields = self.fetcher.read_text(source["md5_url"]).split()
            expected = fields[0].lower() if fields else None
        elif "md5_list" in source:
            expected = md5_from_list(self.fetcher.read_text(source["md5_list"]), file_name)
        log.info("  downloading %s (%s)", source["name"], _size_label(remote.size))
        meta = self.fetcher.download(url, dest, source["max_mb"] * 1024 * 1024, expected)
        meta["license"] = source.get("license", "")
        meta_path.write_text(json.dumps(meta), encoding="utf-8")
        return dest, meta

    def _cleanup_downloads(self, shared_with: list[str]) -> None:
        folder = self.data_dir / DOWNLOAD_DIR
        needed = {source for name in shared_with for source in self.catalog["artifacts"][name]["sources"]}
        for path in folder.iterdir():
            source_id = path.name.split("__", 1)[0]
            if source_id not in needed:
                path.unlink(missing_ok=True)

    # build one artifact
    def _build(self, name: str, options: dict) -> dict:
        """Download inputs and build into staging. Returns the record for the manifest."""
        artifact = self.catalog["artifacts"][name]
        kind = artifact["kind"]
        staging = self.data_dir / STAGING_DIR
        started = time.time()
        metas: dict[str, dict] = {}
        paths: dict[str, Path] = {}
        for source_id in artifact["sources"]:
            if self.catalog["sources"][source_id].get("stream"):
                continue
            paths[source_id], metas[source_id] = self._download(source_id)

        log.info("  building %s", name)
        build = artifact.get("build")
        if kind == "genome":
            out = staging / f"{build}.fa"
            stats = builders.build_genome(paths[artifact["sources"][0]], out)
            files = [out, Path(str(out) + ".fai")]
        else:
            out = staging / f"{name}.sqlite"
            files = [out]
            if kind == "liftover":
                source_build, target_build = name.split(".")[1].split("_")
                stats = builders.build_liftover(paths[artifact["sources"][0]], out, source_build, target_build)
            elif kind == "transcripts":
                refseq, gencode, mane = artifact["sources"]
                inputs = {"refseq": paths[refseq], "gencode": paths[gencode], "mane": paths[mane]}
                stats = builders.build_transcripts(inputs, out, build)
            elif kind == "clinvar":
                vcf, summary = artifact["sources"]
                genome = self._genome_for(build)
                try:
                    stats = builders.build_clinvar(paths[vcf], paths[summary], out, build, genome)
                finally:
                    genome.close()
            elif kind == "genes":
                stats = builders.build_genes(
                    {
                        "hgnc": paths["hgnc"],
                        "clingen_validity": paths["clingen_validity"],
                        "clingen_dosage": paths["clingen_dosage"],
                        "gnomad_constraint": paths["gnomad_constraint"],
                        "hpo_diseases": paths["hpo_diseases"],
                    },
                    out,
                )
            elif kind == "proteins":
                stats, metas["pfam_regions"] = self._build_proteins(paths["uniprot"], out, options)
            elif kind == "civic":
                stats = builders.build_civic(paths["civic_variants"], paths["civic_evidence"], out)
            elif kind == "alphamissense":
                stats = builders.build_alphamissense(paths[artifact["sources"][0]], out, build)
            elif kind == "gnomad":
                stats, metas["gnomad"] = self._build_gnomad(out, build, options)
            else:
                raise SetupError(f"no builder for {kind}")
        return {
            "name": name,
            "files": files,
            "stats": stats,
            "metas": metas,
            "seconds": round(time.time() - started, 1),
        }

    def _genome_for(self, build: str) -> Genome:
        fasta = self.data_dir / GENOME_DIR / f"{build}.fa"
        if not fasta.is_file():
            raise SetupError(f"the {build} reference must be installed first")
        return open_build_genome(self.data_dir / GENOME_DIR, build)

    def _build_proteins(self, uniprot: Path, out: Path, options: dict):
        source = self.catalog["sources"]["pfam_regions"]
        cache = self.data_dir / CACHE_DIR / "pfam_regions.kept.tsv.gz"
        cache_meta = self.data_dir / CACHE_DIR / "pfam_regions.kept.json"
        if options.get("skip_pfam_regions"):
            log.warning("  skipping Pfam regions as requested; domain positions will be missing")
            return builders.build_proteins(uniprot, None, out), {"skipped": True}
        remote = self.fetcher.head(source["url"])
        if cache.is_file() and cache_meta.is_file():
            meta = json.loads(cache_meta.read_text(encoding="utf-8"))
            if (meta.get("etag"), meta.get("last_modified")) == (remote.etag, remote.last_modified):
                log.info("  Pfam unchanged; using the rows kept from the last read")
                with builders.open_text(cache) as handle:
                    return builders.build_proteins(uniprot, handle, out), meta
        log.info("  reading Pfam regions (%s) without saving the whole file", _size_label(remote.size))
        stream = self.fetcher.stream(source["url"], source["max_mb"] * 1024 * 1024)
        temp_cache = cache.with_name(cache.name + ".tmp")
        with stream as raw, _gzip_text(raw) as lines:
            stats = builders.build_proteins(uniprot, lines, out, temp_cache)
        meta = stream.metadata()
        meta["license"] = source.get("license", "")
        os.replace(temp_cache, cache)
        cache_meta.write_text(json.dumps(meta), encoding="utf-8")
        return stats, meta

    def _build_gnomad(self, out: Path, build: str, options: dict):
        source = self.catalog["sources"][f"gnomad_{build}"]
        chroms = options.get("gnomad_chroms") or [c[3:] for c in primary_contigs()]
        genome = self._genome_for(build)
        metas = []

        def chromosomes():
            for chrom in chroms:
                log.info("  reading gnomAD chromosome %s", chrom)
                url = source["url_template"].format(chrom=chrom)
                stream = self.fetcher.stream(url, source["max_mb"] * 1024 * 1024)
                with stream as raw, _gzip_text(raw) as lines:
                    yield chrom, lines
                metas.append(stream.metadata())

        try:
            stats = builders.build_gnomad(chromosomes(), out, build, genome)
        finally:
            genome.close()
        return stats, {"chromosomes": chroms, "files": metas, "license": source.get("license", "")}

    # publish
    def _entry(self, built: dict) -> dict:
        name = built["name"]
        artifact = self.catalog["artifacts"][name]
        primary_id = artifact["sources"][0]
        primary = built["metas"].get(primary_id, {})
        version = built["stats"].get("file_date") or _date_of(primary.get("last_modified", "")) or utc_now()[:10]
        return {
            "policy": artifact["policy"],
            "kind": artifact["kind"],
            "built_at": utc_now(),
            "source": {"version": version, "name": self.catalog["sources"][primary_id]["name"]},
            "sources": built["metas"],
            "stats": built["stats"],
            "build_seconds": built["seconds"],
        }

    def _publish(self, manifest: dict, built: dict, snapshot: Path | None) -> None:
        name = built["name"]
        targets = self._artifact_paths(name)
        if snapshot is not None:
            old_entry = manifest["files"].get(name)
            if old_entry is not None:
                snapshot.mkdir(parents=True, exist_ok=True)
                for target in targets:
                    if target.exists():
                        shutil.move(str(target), str(snapshot / target.name))
                record_path = snapshot / "entries.json"
                record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.exists() else {}
                record[name] = old_entry
                record_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        for source, target in zip(built["files"], targets):
            os.replace(source, target)
            if os.name != "nt":
                target.chmod(0o444 if target.suffix == ".sqlite" else 0o644)
        entry = self._entry(built)
        entry["files"] = {t.name: {"size": t.stat().st_size, "sha256": _sha256(t)} for t in targets}
        manifest["files"][name] = entry
        self._save_manifest(manifest)

    # setup
    def setup(self, builds: list[str], optional: list[str] | None = None, options: dict | None = None, force: bool = False) -> list[dict]:
        options = options or {}
        self._prepare_folders()
        with _Lock(self.data_dir / ".lock"):
            manifest = self.load_manifest()
            manifest["builds"] = sorted(set(manifest.get("builds", [])) | set(builds))
            self._save_manifest(manifest)
            plan = self._plan(sorted(set(manifest["builds"])), optional or [])
            done = []
            for number, name in enumerate(plan, 1):
                if name in manifest["files"] and not force and all(p.exists() for p in self._artifact_paths(name)):
                    log.info("[%d/%d] %s is already installed", number, len(plan), name)
                    continue
                log.info("[%d/%d] %s", number, len(plan), name)
                built = self._build(name, options)
                self._publish(manifest, built, None)
                done.append({"artifact": name, "stats": built["stats"], "seconds": built["seconds"]})
                self._cleanup_downloads(plan[number:])
            self._cleanup_downloads([])
            self._write_report("setup", {"builds": manifest["builds"], "installed": done})
            return done

    # update
    def check_updates(self) -> list[UpdateItem]:
        manifest = self.load_manifest()
        if not manifest["files"]:
            raise DataError(f"nothing is installed in {self.data_dir}; run: varsight setup")
        items = []
        for name, entry in sorted(manifest["files"].items()):
            artifact = self.catalog["artifacts"].get(name)
            if artifact is None:
                continue
            item = UpdateItem(name, artifact["policy"])
            for source_id in artifact["sources"]:
                stored = entry.get("sources", {}).get(source_id)
                if stored is None or stored.get("skipped") or self.catalog["sources"][source_id].get("url_template"):
                    continue
                try:
                    url = self._source_url(source_id)
                    remote = self.fetcher.head(url)
                except SetupError as exc:
                    item.error = str(exc)
                    break
                same_file = url == stored.get("url")
                same_content = (remote.etag, remote.last_modified) == (stored.get("etag"), stored.get("last_modified"))
                if not (same_file and same_content):
                    item.changed.append(source_id)
            items.append(item)
        return items

    def update(self, approve: list[str] | None = None, check_only: bool = False, options: dict | None = None) -> list[UpdateItem]:
        """Apply new source versions. Automatic sources go ahead; the rest wait for approval."""
        approve = approve or []
        options = options or {}
        self._prepare_folders()
        with _Lock(self.data_dir / ".lock"):
            items = self.check_updates()
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            manifest = self.load_manifest()
            snapshot = self.data_dir / SNAPSHOT_DIR / stamp
            proteins = manifest["files"].get("proteins", {}).get("sources", {})
            if proteins.get("pfam_regions", {}).get("skipped"):
                # Someone who chose to skip the 4.7 GB Pfam read at setup keeps that choice.
                options = {**options, "skip_pfam_regions": True}
            order = {"genome": 0, "liftover": 1, "transcripts": 2, "clinvar": 3}
            for item in sorted(items, key=lambda i: order.get(self.catalog["artifacts"][i.artifact]["kind"], 9)):
                if item.error:
                    item.action = "error"
                    item.detail = item.error
                    continue
                if not item.changed:
                    item.detail = "up to date"
                    continue
                approved = "all" in approve or item.artifact in approve or item.artifact.split(".")[0] in approve
                if item.policy == "once" and not approved:
                    item.action = "skipped"
                    item.detail = "this source is not expected to change; approve it to rebuild"
                    continue
                if item.policy in ("approval", "optional") and not approved:
                    item.action = "needs_approval"
                    item.detail = f"can change results; run: varsight data update --approve {item.artifact}"
                    continue
                if check_only:
                    item.action = "available"
                    continue
                self._apply_update(manifest, item, snapshot, options)
            self._cleanup_downloads([])
            self._prune_snapshots()
            self._write_report("update", {"items": [vars(i) for i in items]})
            return items

    def _apply_update(self, manifest: dict, item: UpdateItem, snapshot: Path, options: dict) -> None:
        name = item.artifact
        log.info("updating %s (new: %s)", name, ", ".join(item.changed))
        try:
            built = self._build(name, options)
            self._quality_gate(name, manifest["files"].get(name, {}), built)
            self._self_test_new(name, built)
            self._publish(manifest, built, snapshot)
            item.action = "updated"
            item.detail = f"now version {manifest['files'][name]['source']['version']}"
            log.info("  %s %s", name, item.detail)
        except Exception as exc:  # noqa: BLE001 - any failure must leave the current data in place
            item.action = "failed"
            item.detail = f"kept the current data: {exc}"
            log.error("  %s %s", name, item.detail)
            for path in (self.data_dir / STAGING_DIR).iterdir():
                path.unlink(missing_ok=True)

    def _quality_gate(self, name: str, old_entry: dict, built: dict) -> None:
        kind = self.catalog["artifacts"][name]["kind"]
        metric = _GATE_METRIC.get(kind)
        if metric is None or not old_entry:
            return
        before = old_entry.get("stats", {}).get(metric, 0)
        after = built["stats"].get(metric, 0)
        if before and after < before * SHRINK_LIMIT:
            raise SetupError(f"{name} would shrink from {before:,} to {after:,} {metric}, so it looks wrong")

    def _self_test_new(self, name: str, built: dict) -> None:
        if not self.manifest_path.is_file() or built["files"][0].suffix != ".sqlite":
            return
        from .selftest import run_selftest

        store = DataStore(self.data_dir, overrides={name: built["files"][0]})
        try:
            failures = [r for r in run_selftest(store) if not r["passed"]]
        finally:
            store.close()
        if failures:
            first = failures[0]
            raise SetupError(f"self-test failed on {first['name']}: {first['detail']}")

    # rollback
    def list_snapshots(self) -> list[str]:
        folder = self.data_dir / SNAPSHOT_DIR
        return sorted(p.name for p in folder.iterdir() if p.is_dir()) if folder.is_dir() else []

    def rollback(self, stamp: str | None = None) -> list[str]:
        """Bring back the databases saved in a snapshot, newest by default."""
        snapshots = self.list_snapshots()
        if not snapshots:
            raise DataError("there is no older data to go back to yet")
        chosen = stamp or snapshots[-1]
        if chosen not in snapshots:
            raise DataError(f"no snapshot named {chosen}; available: {', '.join(snapshots)}")
        folder = self.data_dir / SNAPSHOT_DIR / chosen
        with _Lock(self.data_dir / ".lock"):
            manifest = self.load_manifest()
            record = json.loads((folder / "entries.json").read_text(encoding="utf-8"))
            restored = []
            for name, entry in record.items():
                for target in self._artifact_paths(name):
                    saved = folder / target.name
                    if not saved.exists():
                        raise DataError(f"snapshot {chosen} is missing {target.name}")
                    target.unlink(missing_ok=True)
                    shutil.move(str(saved), str(target))
                manifest["files"][name] = entry
                restored.append(name)
            self._save_manifest(manifest)
            shutil.rmtree(folder)
            self._write_report("rollback", {"snapshot": chosen, "restored": restored})
            return restored

    def _prune_snapshots(self) -> None:
        snapshots = self.list_snapshots()
        for stamp in snapshots[:-KEEP_SNAPSHOTS]:
            shutil.rmtree(self.data_dir / SNAPSHOT_DIR / stamp, ignore_errors=True)

    # verify
    def verify(self) -> list[tuple[str, bool, str]]:
        """Compare every installed file with the checksum recorded when it was built."""
        results = []
        manifest = self.load_manifest()
        for name, entry in sorted(manifest["files"].items()):
            for target in self._artifact_paths(name):
                recorded = entry.get("files", {}).get(target.name, {})
                if not target.is_file():
                    results.append((target.name, False, "missing"))
                elif _sha256(target) != recorded.get("sha256"):
                    results.append((target.name, False, "checksum differs from the one recorded at build"))
                else:
                    results.append((target.name, True, "ok"))
        return results

    def _write_report(self, action: str, body: dict) -> Path:
        folder = self.data_dir / REPORT_DIR
        folder.mkdir(exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        path = folder / f"{action}-{stamp}.json"
        path.write_text(json.dumps({"action": action, "time": utc_now(), **body}, indent=2, default=str), encoding="utf-8")
        return path


def _gzip_text(raw: io.RawIOBase) -> io.TextIOBase:
    import gzip

    buffered = io.BufferedReader(raw, buffer_size=CHUNK)
    return gzip.open(buffered, "rt", encoding="utf-8", errors="replace", newline="\n")


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(CHUNK):
            digest.update(block)
    return digest.hexdigest()


def _date_of(http_date: str) -> str:
    if not http_date:
        return ""
    try:
        return parsedate_to_datetime(http_date).strftime("%Y-%m-%d")
    except (TypeError, ValueError):
        return ""


def _size_label(size: int | None) -> str:
    return f"{size / 1024**2:,.1f} MB" if size else "size unknown"

