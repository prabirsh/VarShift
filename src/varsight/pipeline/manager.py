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
import threading
import time
from contextlib import contextmanager
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
from ..advanced import ACCEPT, ALIAS, OWN_FILE, artifact_names, get_database
from . import advanced_builders, builders
from .fetch import CHUNK, Fetcher, duration_label, md5_from_list

log = logging.getLogger("varsight.data")

SCORE_DIR = "scores"
SNAPSHOT_DIR = "snapshots"
DOWNLOAD_DIR = "downloads"
STAGING_DIR = "staging"
REPORT_DIR = "reports"
STILL_GOING_EVERY = 60
KEEP_SNAPSHOTS = 2
# A new database may not shrink past this share of the one it replaces.
SHRINK_LIMIT = 0.9
_GATE_METRIC = {
    "clinvar": "rows",
    "transcripts": "refseq",
    "genes": "genes",
    "proteins": "reviewed_proteins",
    "civic": "variants",
    "dbsnp": "rows",
    "thousand_genomes": "rows",
    "gnomad": "rows",
    "gnomad_sv": "rows",
    "revel": "rows",
    "alphamissense": "rows",
    "hotspots": "rows",
    "clingen_variants": "rows",
}
# A source that sends no version headers is simply read again after this long.
RECHECK_DAYS = 7


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
        for name in (DB_DIR, GENOME_DIR, SCORE_DIR, SNAPSHOT_DIR, DOWNLOAD_DIR, STAGING_DIR, REPORT_DIR):
            (self.data_dir / name).mkdir(exist_ok=True)

    # planning
    def _plan(self, builds: list[str]) -> list[str]:
        order = ["genes", "civic", "proteins"]
        for build in builds:
            order.insert(0, f"genome.{build}")
            order.append(f"transcripts.{build}")
            order.append(f"clinvar.{build}")
        if set(builds) >= set(BUILDS):
            order.extend(["liftover.hg19_hg38", "liftover.hg38_hg19"])
        return order

    def _artifact_paths(self, name: str) -> list[Path]:
        kind = self.catalog["artifacts"][name]["kind"]
        if kind == "genome":
            build = self.catalog["artifacts"][name]["build"]
            fasta = self.data_dir / GENOME_DIR / f"{build}.fa"
            return [fasta, Path(str(fasta) + ".fai")]
        if kind == "phylop":
            return [self.data_dir / SCORE_DIR / f"{name}.bw"]
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
        started = time.time()
        metas: dict[str, dict] = {}
        paths: dict[str, Path] = {}
        for source_id in artifact["sources"]:
            if self.catalog["sources"][source_id].get("stream"):
                continue
            paths[source_id], metas[source_id] = self._download(source_id)

        log.info("  building %s", name)
        with _still_going(f"building {name}"):
            files, stats = self._build_files(name, paths, metas, options)
        return {
            "name": name,
            "files": files,
            "stats": stats,
            "metas": metas,
            "seconds": round(time.time() - started, 1),
        }

    def _build_files(self, name: str, paths: dict, metas: dict, options: dict) -> tuple[list[Path], dict]:
        artifact = self.catalog["artifacts"][name]
        kind = artifact["kind"]
        staging = self.data_dir / STAGING_DIR
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
                with builders.open_text(paths["pfam_regions"]) as pfam:
                    stats = builders.build_proteins(paths["uniprot"], pfam, out)
            elif kind == "civic":
                stats = builders.build_civic(paths["civic_variants"], paths["civic_evidence"], out)
            elif kind == "phylop":
                out = staging / f"{name}.bw"
                files = [out]
                os.replace(paths[artifact["sources"][0]], out)
                stats = {"bytes": out.stat().st_size}
            else:
                stats = self._build_advanced(name, kind, build, out, paths, metas, options)
        return files, stats

    def _genome_for(self, build: str) -> Genome:
        fasta = self.data_dir / GENOME_DIR / f"{build}.fa"
        if not fasta.is_file():
            raise SetupError(f"the {build} reference must be installed first")
        return open_build_genome(self.data_dir / GENOME_DIR, build)

    def _build_advanced(self, name: str, kind: str, build: str | None, out: Path, paths: dict, metas: dict, options: dict) -> dict:
        own_file = Path(options["file"]) if options.get("file") else None
        if kind == "revel":
            return advanced_builders.build_revel(paths["revel"], out, build)
        if kind == "alphamissense":
            return advanced_builders.build_alphamissense(paths[f"alphamissense_{build}"], out, build)
        if kind == "hotspots":
            return advanced_builders.build_hotspots(paths["hotspots"], out)
        if kind == "decipher":
            return advanced_builders.build_decipher(paths["decipher_hi"], paths["ddg2p"], out)
        if kind == "clingen_variants":
            return advanced_builders.build_clingen_variants(paths["clingen_variants"], out)
        if kind == "cosmic_census":
            return advanced_builders.build_cosmic_census(own_file, out)
        if kind == "cosmic_fusions":
            return advanced_builders.build_cosmic_fusions(own_file, out)
        if kind == "omim":
            return advanced_builders.build_omim(own_file, out)
        if kind == "primateai":
            return advanced_builders.build_primateai(own_file, out, build)
        if kind == "gnomad_sv":
            source_id = f"gnomad_sv_{build}"
            return self._read_stream(source_id, metas, lambda lines: advanced_builders.build_gnomad_sv(lines, out, build))
        genome = self._genome_for(build)
        try:
            if kind == "cosmic":
                return advanced_builders.build_cosmic(own_file, out, build, genome)
            regions = self._regions(build, options, pad=500 if kind == "spliceai" else 50)
            if kind == "spliceai":
                with builders.open_text(own_file) as lines:
                    return advanced_builders.build_spliceai(lines, out, build, genome, regions)
            if kind == "dbsnp":
                build_step = advanced_builders.build_dbsnp
                source_id = f"dbsnp_{build}"
            elif kind == "thousand_genomes":
                build_step = advanced_builders.build_thousand_genomes
                source_id = f"g1k_{build}"
            elif kind == "gnomad":
                stats, metas[f"gnomad_{build}"] = self._build_gnomad(out, build, genome, regions, options)
                return stats
            else:
                raise SetupError(f"no builder for {kind}")
            return self._read_stream(source_id, metas, lambda lines: build_step(lines, out, build, genome, regions))
        finally:
            genome.close()

    def _regions(self, build: str, options: dict, pad: int):
        if options.get("scope") == "all":
            return None
        return advanced_builders.CodingRegions(self.data_dir / DB_DIR / f"transcripts.{build}.sqlite", pad)

    def _read_stream(self, source_id: str, metas: dict, step) -> dict:
        """Read a very large file straight off the network, without saving it."""
        source = self.catalog["sources"][source_id]
        url = self._source_url(source_id)
        remote = self.fetcher.head(url)
        log.info("  reading %s (%s) without saving the whole file", source["name"], _size_label(remote.size))
        stream = self.fetcher.stream(url, source["max_mb"] * 1024 * 1024)
        with stream as raw, _gzip_text(raw) as lines:
            stats = step(lines)
        metas[source_id] = {**stream.metadata(), "license": source.get("license", "")}
        return stats

    def _build_gnomad(self, out: Path, build: str, genome, regions, options: dict):
        source = self.catalog["sources"][f"gnomad_{build}"]
        chroms = options.get("chrom") or [c[3:] for c in primary_contigs()]
        metas = []

        def chromosomes():
            for chrom in chroms:
                log.info("  reading gnomAD chromosome %s", chrom)
                url = source["url_template"].format(chrom=chrom)
                stream = self.fetcher.stream(url, source["max_mb"] * 1024 * 1024)
                with stream as raw, _gzip_text(raw) as lines:
                    yield chrom, lines
                metas.append(stream.metadata())

        stats = advanced_builders.build_gnomad(chromosomes(), out, build, genome, regions)
        return stats, {"chromosomes": chroms, "files": metas, "license": source.get("license", "")}

    # publish
    def _entry(self, built: dict) -> dict:
        name = built["name"]
        artifact = self.catalog["artifacts"][name]
        primary_id = artifact["sources"][0] if artifact["sources"] else None
        primary = built["metas"].get(primary_id, {})
        version = built["stats"].get("file_date") or _date_of(primary.get("last_modified", "")) or utc_now()[:10]
        source_name = self.catalog["sources"][primary_id]["name"] if primary_id else artifact.get("title", name)
        return {
            "policy": artifact["policy"],
            "kind": artifact["kind"],
            "built_at": utc_now(),
            "source": {"version": version, "name": source_name},
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
    def setup(self, builds: list[str], options: dict | None = None, force: bool = False) -> list[dict]:
        options = options or {}
        self._prepare_folders()
        with _Lock(self.data_dir / ".lock"):
            manifest = self.load_manifest()
            manifest["builds"] = sorted(set(manifest.get("builds", [])) | set(builds))
            self._save_manifest(manifest)
            plan = self._plan(sorted(set(manifest["builds"])))
            began = time.monotonic()
            done = []
            for number, name in enumerate(plan, 1):
                if name in manifest["files"] and not force and all(p.exists() for p in self._artifact_paths(name)):
                    log.info("[%d/%d] %s is already installed", number, len(plan), name)
                    continue
                log.info("[%d/%d] %s", number, len(plan), name)
                built = self._build(name, options)
                self._publish(manifest, built, None)
                log.info("  finished %s in %s", name, duration_label(built["seconds"]))
                done.append({"artifact": name, "stats": built["stats"], "seconds": built["seconds"]})
                self._cleanup_downloads(plan[number:])
            self._cleanup_downloads([])
            self._write_report("setup", {"builds": manifest["builds"], "installed": done})
            log.info("Setup finished in %s", duration_label(time.monotonic() - began))
            return done

    # advanced databases

    def add(self, name: str, options: dict | None = None) -> list[dict]:
        """Install one advanced database for the builds set up here."""
        options = dict(options or {})
        database = get_database(name)
        if database is None:
            raise SetupError(f"no database named {name}; see: varsight databases")
        if database.access == ALIAS:
            raise SetupError(f"{database.title} comes with {database.alias_of}. To get it, run: varsight add {database.alias_of}")
        if database.access in (ACCEPT, OWN_FILE) and not options.get("accepted"):
            raise SetupError(f"{database.title}: {database.license}. If that fits your use, add the word accept to the command")
        if database.access == OWN_FILE:
            file = options.get("file")
            if not file:
                raise SetupError(f"{database.title} needs your own licensed file. {database.get_it}. Then run: varsight add {database.name} file=PATH accept")
            if not Path(file).is_file():
                raise SetupError(f"cannot find the file {file}")
        self._prepare_folders()
        with _Lock(self.data_dir / ".lock"):
            manifest = self.load_manifest()
            builds = manifest.get("builds", [])
            if not builds:
                raise DataError(f"nothing is installed in {self.data_dir} yet; run: varsight setup")
            if database.access == OWN_FILE and database.per_build:
                builds = [_build_of_file(options, builds, database)]
            done = []
            for artifact in artifact_names(database, builds):
                log.info("adding %s", artifact)
                built = self._build(artifact, options)
                if database.access == OWN_FILE:
                    built["metas"] = {"own_file": {"name": Path(options["file"]).name, "license": database.license}}
                self._publish(manifest, built, None)
                done.append({"artifact": artifact, "stats": built["stats"], "seconds": built["seconds"]})
            self._cleanup_downloads([])
            self._write_report("add", {"database": database.name, "installed": done})
            return done

    def remove(self, name: str) -> list[str]:
        database = get_database(name)
        if database is None or database.access == ALIAS:
            raise SetupError(f"{name} is not an advanced database you can remove; see: varsight databases")
        with _Lock(self.data_dir / ".lock"):
            manifest = self.load_manifest()
            gone = []
            for artifact in list(manifest.get("files", {})):
                if artifact == database.name or artifact.startswith(f"{database.name}."):
                    for target in self._artifact_paths(artifact):
                        target.unlink(missing_ok=True)
                    del manifest["files"][artifact]
                    gone.append(artifact)
            if not gone:
                raise SetupError(f"{database.title} is not installed")
            self._save_manifest(manifest)
            return gone

    # update

    def check_updates(self) -> list[UpdateItem]:
        manifest = self.load_manifest()
        if not manifest["files"]:
            raise DataError(f"nothing is installed in {self.data_dir}; run: varsight setup")
        items = []
        for name, entry in sorted(manifest["files"].items()):
            artifact = self.catalog["artifacts"].get(name)
            if artifact is None or artifact["policy"] == "own_file":
                continue
            item = UpdateItem(name, artifact["policy"])
            for source_id in artifact["sources"]:
                stored = entry.get("sources", {}).get(source_id)
                if stored is None or self.catalog["sources"][source_id].get("url_template"):
                    continue
                try:
                    url = self._source_url(source_id)
                    remote = self.fetcher.head(url)
                except SetupError as exc:
                    item.error = str(exc)
                    break
                same_file = url == stored.get("url")
                same_content = (remote.etag, remote.last_modified) == (stored.get("etag"), stored.get("last_modified"))
                if not (remote.etag or remote.last_modified):
                    # No version headers, so a small file like this is just read again after a while.
                    same_content = not _older_than(entry.get("built_at", ""), RECHECK_DAYS)
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
                    item.detail = f"can change results; to apply: varsight update approve={item.artifact}"
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


def _build_of_file(options: dict, builds: list[str], database) -> str:
    """Which build a licensed file is for: from build=, its name, or the only build installed."""
    chosen = options.get("build")
    file_name = Path(options["file"]).name.lower()
    if not chosen:
        if "grch37" in file_name or "hg19" in file_name or "b37" in file_name:
            chosen = "hg19"
        elif "grch38" in file_name or "hg38" in file_name or "b38" in file_name:
            chosen = "hg38"
        elif len(builds) == 1:
            chosen = builds[0]
    if chosen is None:
        raise SetupError(f"cannot tell which genome build {Path(options['file']).name} is for; add build=hg19 or build=hg38")
    if chosen not in builds:
        raise SetupError(f"{database.title} is for {chosen}, but only {', '.join(builds)} is set up here")
    return chosen


@contextmanager
def _still_going(what: str):
    """Some builds run for minutes with nothing to report. A line each minute shows they haven't stopped."""
    began = time.monotonic()
    finished = threading.Event()

    def speak():
        while not finished.wait(STILL_GOING_EVERY):
            log.info("    still %s, %s so far", what, duration_label(time.monotonic() - began))

    speaker = threading.Thread(target=speak, daemon=True)
    speaker.start()
    try:
        yield
    finally:
        finished.set()
        speaker.join()


def _older_than(stamp: str, days: int) -> bool:
    try:
        built = datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except ValueError:
        return True
    return (datetime.now(timezone.utc) - built).days >= days


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

