"""The varsight command."""

from __future__ import annotations

import logging
import os
import shutil
import sqlite3
import sys
from pathlib import Path

from . import __version__
from .advanced import ACCEPT, ALIAS, DATABASES, DOWNLOAD, OWN_FILE, get_database
from .annotate import Annotator
from .batch import annotate_list, annotate_vcf, vcf_columns
from .commandline import COMMANDS, Values, parse, render_help
from .config import default_data_dir, is_windows_mount, load_settings, resolve_data_dir, save_setting
from .errors import InputError, ResolveError, SetupError, VarSightError
from .output import ReportContext, ReportFolder, ScreenSink, format_of, open_file_sink
from .store import DataStore
from .updates import check_in_background, clear_applied, notice_lines, read_status
from .vcf import looks_like_vcf

log = logging.getLogger("varsight")
_PARALLEL_FROM = 1024 * 1024


def main(argv: list[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        # A Windows console can't show every character. Replace those rather than stop.
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)
    try:
        values = parse(sys.argv[1:] if argv is None else argv)
        return HANDLERS[values.command](values)
    except ResolveError as exc:
        print(f"error: {exc}", file=sys.stderr)
        for line in exc.tried:
            print(f"  tried {line}", file=sys.stderr)
        return 1
    except VarSightError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("stopped", file=sys.stderr)
        return 130
    except BrokenPipeError:
        return 0
    except OSError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _help(values: Values) -> int:
    print(render_help(values.positional))
    return 0


def _version(values: Values) -> int:
    print(f"VarSight {__version__}")
    return 0


# annotate
def _annotate(values: Values) -> int:
    source = values.positional
    data_dir = resolve_data_dir(values.get("data"))
    store = DataStore(data_dir)
    try:
        notices = notice_lines(data_dir, store.manifest)
        for line in notices:
            print(line, file=sys.stderr)
        check_in_background(data_dir, store.manifest)
        context = ReportContext(
            input_name=Path(source).name if Path(source).is_file() else source,
            installed=store.installed(),
            versions=store.versions(),
            notices=[line.strip() for line in notices],
        )
        annotator = Annotator(store, all_transcripts=values.get("transcripts") == "all")
        if looks_like_vcf(source):
            context.prefix = vcf_columns(Path(source))
            sink = _sink_for_file(source, values, context)
            counts = annotate_vcf(annotator, Path(source), sink, values.get("build"), _worker_count(values, source))
            return _finish_file(sink, counts)
        if Path(source).is_file():
            sink = _sink_for_file(source, values, context)
            counts = annotate_list(annotator, Path(source), sink, values.get("build"), values.get("gene") or "")
            return _finish_file(sink, counts)
        records = annotator.annotate_text(source, values.get("build"), values.get("gene") or "")
        sink = _sink_for_single(values, context)
        for record in records:
            sink.write(record)
        sink.close()
        if hasattr(sink, "final_path"):
            print(f"Wrote {len(records)} result(s) to {sink.final_path}", file=sys.stderr)
        return 0
    finally:
        store.close()


def _worker_count(values: Values, source: str) -> int:
    if values.get("workers"):
        return max(_whole_number(values, "workers", 1), 1)
    # Starting processes takes a moment, so a small file is quicker in just one.
    return 1 if Path(source).stat().st_size < _PARALLEL_FROM else min(os.cpu_count() or 1, 4)


def _sink_for_single(values: Values, context: ReportContext):
    out = values.get("out")
    if out:
        if format_of(out) is None:
            return ReportFolder(out, context, values.has("replace"))
        return open_file_sink(out, context, values.has("replace"))
    return ScreenSink(_screen_format(values), context)


def _screen_format(values: Values) -> str:
    fmt = values.get("format") or "card"
    if fmt not in ("card", "json", "jsonl", "tsv"):
        raise InputError(f"format must be card, json, jsonl or tsv, not {fmt!r}")
    return fmt


def _sink_for_file(source: str, values: Values, context: ReportContext):
    out = values.get("out")
    if values.get("format"):
        if out:
            raise InputError("use out= or format=, not both")
        return ScreenSink(_screen_format(values), context)
    if out and format_of(out) is not None:
        return open_file_sink(out, context, values.has("replace"))
    return ReportFolder(out or _default_folder(source), context, values.has("replace"))


def _default_folder(source: str) -> str:
    name = Path(source).name
    for ending in (".vcf.gz", ".vcf.bgz", ".vcf", ".txt", ".tsv"):
        if name.lower().endswith(ending):
            name = name[: -len(ending)]
            break
    return f"{name}_varsight"


def _finish_file(sink, counts) -> int:
    sink.close()
    summary = ", ".join(f"{name.replace('_', ' ')}: {value:,}" for name, value in counts.items())
    print(summary, file=sys.stderr)
    if isinstance(sink, ReportFolder):
        print(f"Open this to explore the results: {sink.paths[0]}", file=sys.stderr)
        for path in sink.paths[1:]:
            print(f"Also written: {path}", file=sys.stderr)
    elif hasattr(sink, "final_path"):
        print(f"Wrote {sink.final_path}", file=sys.stderr)
    return 0


# setup
_BUILD_CHOICES = {"1": "hg19", "2": "hg38", "3": "both", "hg19": "hg19", "hg38": "hg38", "both": "both"}
_CORE_NEEDS = {
    "hg19": ("about 2.3 GB, plus 4.7 GB of Pfam read without saving", "about 5 GB", "about 15 GB"),
    "hg38": ("about 2.4 GB, plus 4.7 GB of Pfam read without saving", "about 5 GB", "about 15 GB"),
    "both": ("about 3.3 GB, plus 4.7 GB of Pfam read without saving", "about 10 GB", "about 25 GB"),
}


def _ask(question: str, default: str = "") -> str:
    prompt = f"{question} [{default}] " if default else f"{question} "
    try:
        answer = input(prompt).strip()
    except EOFError:
        return default
    return answer or default


def _ask_build() -> str:
    print("\nWhich genome build do you work with?")
    print("  1  hg19 (GRCh37)")
    print("  2  hg38 (GRCh38)")
    print("  3  Both. Results then always show the position on both builds")
    while True:
        answer = _ask("Your choice, 1, 2 or 3:", "3").lower()
        if answer in _BUILD_CHOICES:
            return _BUILD_CHOICES[answer]
        print("Please type 1, 2 or 3.")


def _ask_advanced() -> list[str]:
    choices = [db for db in DATABASES.values() if db.access in (DOWNLOAD, ACCEPT)]
    print("\nAdvanced databases are optional, and can be added any time later.")
    for number, database in enumerate(choices, 1):
        print(f"  {number:>2}  {database.title:<34} {database.transfer}")
    print("Databases that need your own account (COSMIC, OMIM, SpliceAI, PrimateAI) are added later: varsight databases")
    answer = _ask("Numbers to add now, separated by commas, or press Enter for none:")
    picked = []
    for item in answer.split(","):
        item = item.strip()
        if not item:
            continue
        if item.isdigit() and 1 <= int(item) <= len(choices):
            picked.append(choices[int(item) - 1].name)
        else:
            print(f"Skipping {item!r}: not one of the numbers shown.")
    return picked


def _setup(values: Values) -> int:
    from .pipeline.manager import DataManager

    interactive = sys.stdin.isatty() and not values.has("yes")
    print("VarSight setup")
    build = values.get("build")
    if build is None:
        build = _ask_build() if interactive else "both"
    build = _BUILD_CHOICES.get(build.lower())
    if build is None:
        raise InputError("build must be hg19, hg38 or both")
    data_text = values.get("data")
    if data_text is None and interactive:
        data_text = _ask("\nWhere should the data go?", str(resolve_data_dir()))
    data_dir = resolve_data_dir(data_text)
    with_names = values.items("with")
    if interactive and not with_names:
        with_names = _ask_advanced()
    accepted = set(values.items("accept"))
    for name in with_names:
        database = get_database(name)
        if database is None:
            raise InputError(f"no advanced database named {name}; see: varsight databases")
        if database.access == OWN_FILE:
            raise InputError(f"{database.title} needs your own file. After setup, run: varsight add {database.name} file=PATH accept")
        if database.access == ACCEPT and database.name not in accepted:
            print(f"\n{database.title}: {database.license}")
            if interactive and _ask("Does that fit your use? yes or no:", "no").lower() in ("y", "yes"):
                accepted.add(database.name)
            else:
                raise InputError(f"{database.title} needs you to accept its terms; add accept={database.name} if they fit your use")
    if is_windows_mount(data_dir):
        log.warning(
            "\nNote: %s is a Windows drive seen from WSL. Lookups measured about 50 to 65 times slower there than on the Linux disk.",
            data_dir,
        )
    transfer, kept, free = _CORE_NEEDS[build]
    print(f"\nThis will install: {'hg19 and hg38' if build == 'both' else build} into {data_dir}")
    print(f"  Download   {transfer}")
    print(f"  Disk       {kept} kept, and keep {free} free during setup (estimates)")
    if with_names:
        print(f"  Advanced   {', '.join(with_names)}")
    if interactive and _ask("Start? yes or no:", "yes").lower() not in ("y", "yes"):
        print("Nothing was changed.")
        return 1
    builds = ["hg19", "hg38"] if build == "both" else [build]
    manager = DataManager(data_dir)
    manager.setup(builds, {"skip_pfam_regions": values.has("nopfam")}, values.has("rebuild"))
    for name in with_names:
        manager.add(name, {"accepted": True})
    save_setting("data_dir", str(data_dir))
    print(f"\nInstalled into {data_dir}. Checking it with known variants.", file=sys.stderr)
    values.options["data"] = str(data_dir)
    return _selftest(values)


def _add(values: Values) -> int:
    from .pipeline.manager import DataManager

    data_dir = resolve_data_dir(values.get("data"))
    database = get_database(values.positional)
    if database is None:
        raise InputError(f"no advanced database named {values.positional}; see: varsight databases")
    if database.access == ALIAS:
        print(f"{database.title}: {database.provides}.")
        print(f"To get it, run: varsight add {database.alias_of}")
        return 0
    accepted = values.has("accept")
    if database.access in (ACCEPT, OWN_FILE) and not accepted and sys.stdin.isatty():
        print(f"{database.title}: {database.license}")
        accepted = _ask("Does that fit your use? yes or no:", "no").lower() in ("y", "yes")
    if database.access == DOWNLOAD and sys.stdin.isatty() and not values.has("yes"):
        print(f"{database.title}: {database.provides}")
        print(f"  Transfer  {database.transfer}\n  Kept      {database.kept}")
        if _ask("Start? yes or no:", "yes").lower() not in ("y", "yes"):
            print("Nothing was changed.")
            return 1
    options = {"accepted": accepted, "file": values.get("file"), "scope": values.get("scope"), "chrom": values.items("chrom"), "build": values.get("build")}
    done = DataManager(data_dir).add(database.name, options)
    for item in done:
        print(f"{item['artifact']}: {item['stats'].get('rows', 'installed')} rows, {item['seconds']} s")
    return 0


def _remove(values: Values) -> int:
    from .pipeline.manager import DataManager

    gone = DataManager(resolve_data_dir(values.get("data"))).remove(values.positional)
    print("Removed: " + ", ".join(gone))
    return 0


def _databases(values: Values) -> int:
    data_dir = resolve_data_dir(values.get("data"))
    installed: list[str] = []
    try:
        store = DataStore(data_dir)
        installed = store.installed()
        store.close()
    except VarSightError:
        print(f"Nothing is installed in {data_dir} yet. Start with: varsight setup\n")
    print("Core databases, installed by setup: UCSC reference, RefSeq, GENCODE, MANE, ClinVar, HGNC, UniProt, Pfam, CIViC,")
    print("ClinGen gene curation, gnomAD gene constraint, HPO, and the hg19 to hg38 mapping.\n")
    print("Advanced databases, optional:\n")
    for database in DATABASES.values():
        have = any(a == database.name or a.startswith(f"{database.name}.") for a in installed)
        if database.access == ALIAS:
            status = f"comes with {database.alias_of}"
        else:
            status = "INSTALLED" if have else "not installed"
        print(f"{database.title} ({database.name}): {status}")
        print(f"  {database.provides}")
        print(f"  Transfer: {database.transfer}.  Kept: {database.kept}")
        print(f"  Licence: {database.license}")
        if database.access == OWN_FILE:
            print(f"  How to get it: {database.get_it}")
            print(f"  Then run: varsight add {database.name} file=PATH accept")
        elif database.access == ACCEPT:
            print(f"  To add: varsight add {database.name} accept")
        elif database.access == DOWNLOAD:
            print(f"  To add: varsight add {database.name}")
        print()
    return 0


# update and rollback
def _update(values: Values) -> int:
    from .pipeline.manager import DataManager

    data_dir = resolve_data_dir(values.get("data"))
    manager = DataManager(data_dir)
    items = manager.update(values.items("approve"), check_only=values.has("check"))
    failed = False
    updated = 0
    waiting = 0
    for item in items:
        label = {"none": "up to date", "updated": "UPDATED", "failed": "FAILED", "available": "new version"}.get(item.action, item.action.replace("_", " "))
        detail = item.detail if item.action != "none" else ""
        print(f"{item.artifact:24} {label:16} {detail}")
        failed = failed or item.action in ("failed", "error")
        updated += item.action == "updated"
        waiting += item.action in ("needs_approval", "available")
    clear_applied(data_dir, [item.artifact for item in items if item.action in ("updated", "none")])
    print(f"\n{updated} updated, {waiting} waiting for you.")
    return 1 if failed else 0


def _rollback(values: Values) -> int:
    from .pipeline.manager import DataManager

    restored = DataManager(resolve_data_dir(values.get("data"))).rollback(values.get("to"))
    print("Restored: " + ", ".join(restored))
    return 0


# information
def _status(values: Values) -> int:
    data_dir = resolve_data_dir(values.get("data"))
    store = DataStore(data_dir)
    try:
        manifest = store.manifest
        print(f"Data folder   {data_dir}")
        print(f"Builds        {', '.join(manifest.get('builds', []))}")
        print(f"Last updated  {manifest.get('updated')}")
        print(f"Update check  {read_status(data_dir).get('checked_at', 'not run yet')}")
        for line in notice_lines(data_dir, manifest):
            print(line)
        print(f"\n{'Database':24} {'Version':12} {'Updates':10} {'Size':>10}")
        for name, entry in sorted(manifest.get("files", {}).items()):
            size = sum(item.get("size", 0) for item in entry.get("files", {}).values())
            print(f"{name:24} {entry['source']['version']:12} {entry['policy']:10} {size / 1024**2:>7,.0f} MB")
    finally:
        store.close()
    return 0


def _doctor(values: Values) -> int:
    data_dir = resolve_data_dir(values.get("data"))
    problems = 0
    version = sys.version_info
    python_ok = version >= (3, 10)
    problems += not python_ok
    print(f"{'ok' if python_ok else 'PROBLEM':8} Python {version.major}.{version.minor}.{version.micro} (3.10 or newer needed)")
    print(f"{'ok':8} SQLite {sqlite3.sqlite_version}")
    print(f"{'ok':8} Data folder {data_dir}")
    if is_windows_mount(data_dir):
        print("WARNING  This is a Windows drive seen from WSL, so lookups will be much slower than on the Linux disk")
    probe = data_dir if data_dir.exists() else data_dir.parent
    if probe.exists():
        free = shutil.disk_usage(probe).free / 1024**3
        print(f"{'ok' if free > 30 else 'WARNING':8} {free:.1f} GB free (a first setup of both builds needs about 25 GB)")
    if not (data_dir / "manifest.json").is_file():
        print("Not yet  No data installed. To install: varsight setup")
        return 1 if problems else 0
    from .pipeline.manager import DataManager

    manager = DataManager(data_dir)
    manifest = manager.load_manifest()
    print(f"{'ok':8} Installed builds: {', '.join(manifest.get('builds', []))}")
    for name in sorted(manifest.get("files", {})):
        for target in manager._artifact_paths(name):
            if not target.is_file():
                print(f"PROBLEM  {target} is missing")
                problems += 1
    if values.has("verify"):
        for file_name, passed, detail in manager.verify():
            print(f"{'ok' if passed else 'PROBLEM':8} {file_name}: {detail}")
            problems += not passed
    return 1 if problems else 0


def _selftest(values: Values) -> int:
    from .pipeline.selftest import run_selftest

    store = DataStore(resolve_data_dir(values.get("data")))
    try:
        results = run_selftest(store)
    finally:
        store.close()
    failed = 0
    for result in results:
        print(f"{'ok' if result['passed'] else 'FAIL':5} {result['name']}: {result['detail']}")
        failed += not result["passed"]
    print(f"{len(results) - failed} of {len(results)} known variants passed")
    return 1 if failed or not results else 0


def _validate(values: Values) -> int:
    from .pipeline.validate import render_report, run_validation

    store = DataStore(resolve_data_dir(values.get("data")))
    try:
        builds = [values.get("build")] if values.get("build") else store.builds
        results = run_validation(store, builds, _whole_number(values, "sample", 20000))
    finally:
        store.close()
    report = render_report(results)
    if values.get("report"):
        Path(values.get("report")).write_text(report, encoding="utf-8")
        print(f"Wrote {values.get('report')}", file=sys.stderr)
    else:
        print(report)
    return 0


def _whole_number(values: Values, name: str, default: int) -> int:
    text = values.get(name)
    if text is None:
        return default
    if not text.isdigit():
        raise InputError(f"{name} must be a whole number, not {text!r}")
    return int(text)


def _explain(values: Values) -> int:
    from .explain import run_explain

    data_dir = resolve_data_dir(values.get("data"))
    store = DataStore(data_dir)
    try:
        return run_explain(store, values.positional, values.get("build"), values.get("gene") or "")
    finally:
        store.close()


def _benchmark(values: Values) -> int:
    from .pipeline.benchmark import run_benchmark

    report = run_benchmark(resolve_data_dir(values.get("data")), _whole_number(values, "records", 3000))
    if values.get("report"):
        Path(values.get("report")).write_text(report, encoding="utf-8")
        print(f"Wrote {values.get('report')}", file=sys.stderr)
    print(report)
    return 0


def _demo(values: Values) -> int:
    from .pipeline.demo import run_demo

    return run_demo(values.get("demodir"), Path(values.get("out") or "demo_output"), values.has("replace"))


def _config(values: Values) -> int:
    updates = values.get("updates")
    if updates is not None:
        if updates not in ("on", "off"):
            raise InputError("updates must be on or off")
        save_setting("updates", updates)
    if values.get("data"):
        save_setting("data_dir", str(Path(values.get("data")).expanduser().resolve()))
    settings = load_settings()
    print(f"Data folder            {settings.get('data_dir', default_data_dir())}")
    print(f"Quick update check     {settings.get('updates', 'on')}")
    return 0


HANDLERS = {
    "annotate": _annotate,
    "setup": _setup,
    "databases": _databases,
    "add": _add,
    "remove": _remove,
    "update": _update,
    "rollback": _rollback,
    "status": _status,
    "doctor": _doctor,
    "explain": _explain,
    "selftest": _selftest,
    "validate": _validate,
    "benchmark": _benchmark,
    "demo": _demo,
    "config": _config,
    "help": _help,
    "version": _version,
}
if set(HANDLERS) != set(COMMANDS):
    raise SetupError("every command needs a handler")
