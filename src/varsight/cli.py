"""Command line for VarSight."""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .annotate import Annotator
from .batch import annotate_list, annotate_vcf, vcf_columns
from .config import is_windows_mount, resolve_data_dir, save_data_dir
from .errors import ReferenceMismatch, ResolveError, VarSightError
from .output import ResultWriter, format_for
from .store import DataStore
from .vcf import looks_like_vcf

log = logging.getLogger("varsight")

_STALE_DAYS = 30


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(message)s",
        stream=sys.stderr,
    )
    if not getattr(args, "run", None):
        parser.print_help()
        return 2
    try:
        return args.run(args)
    except ResolveError as exc:
        print(f"error: {exc}", file=sys.stderr)
        for line in exc.tried:
            print(f"  tried {line}", file=sys.stderr)
        return 1
    except ReferenceMismatch as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except VarSightError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("stopped", file=sys.stderr)
        return 130
    except BrokenPipeError:
        return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="varsight",
        description="Offline variant annotation for hg19 and hg38. Give it p., c., g. or a VCF.",
    )
    parser.add_argument("--version", action="version", version=f"VarSight {__version__}")
    commands = parser.add_subparsers(title="commands", metavar="command")

    annotate = commands.add_parser(
        "annotate",
        help="annotate one variant, a list file or a VCF",
        description="INPUT is a description such as 'BRAF V600E', or the path of a .vcf, .vcf.gz or .txt list.",
    )
    annotate.add_argument("input", help="variant description, VCF path or list file")
    annotate.add_argument("--build", help="hg19 or hg38. Needed for bare genomic positions; read from the header for a VCF")
    annotate.add_argument("--gene", default="", help="gene to use when the description has none")
    annotate.add_argument("--out", help="output file; the extension sets the format (.json .jsonl .tsv .xlsx)")
    annotate.add_argument("--format", choices=["json", "jsonl", "tsv", "xlsx"], help="override the output format")
    annotate.add_argument("--force", action="store_true", help="replace the output file if it exists")
    annotate.add_argument("--all-transcripts", action="store_true", help="report every transcript, not the best five per gene")
    _add_common(annotate)
    annotate.set_defaults(run=_cmd_annotate)

    setup = commands.add_parser("setup", help="download and build the local data (first run)")
    setup.add_argument("--build", choices=["hg19", "hg38", "both"], default="both")
    setup.add_argument(
        "--with",
        dest="optional",
        action="append",
        default=[],
        choices=["alphamissense", "gnomad"],
        help="add an optional source; see docs/data-sources.md for size and licence",
    )
    setup.add_argument("--accept-license", action="append", default=[], metavar="SOURCE", help="confirm you accept the licence of an optional source")
    setup.add_argument("--gnomad-chrom", action="append", default=[], metavar="CHROM", help="limit gnomAD to these chromosomes")
    setup.add_argument("--skip-pfam-regions", action="store_true", help="skip the 4.7 GB Pfam read; domain positions will be missing")
    setup.add_argument("--force", action="store_true", help="rebuild even what is already installed")
    setup.add_argument("--yes", action="store_true", help="do not ask for confirmation")
    _add_common(setup)
    setup.set_defaults(run=_cmd_setup)

    status = commands.add_parser("status", help="show what is installed and how old it is")
    _add_common(status)
    status.set_defaults(run=_cmd_status)

    doctor = commands.add_parser("doctor", help="check this machine and the data folder")
    doctor.add_argument("--verify", action="store_true", help="also compare every data file with its recorded checksum")
    _add_common(doctor)
    doctor.set_defaults(run=_cmd_doctor)

    selftest = commands.add_parser("selftest", help="run the fixed list of known variants")
    _add_common(selftest)
    selftest.set_defaults(run=_cmd_selftest)

    validate = commands.add_parser("validate", help="measure accuracy against ClinVar's own names")
    validate.add_argument("--build", choices=["hg19", "hg38"], help="limit to one build")
    validate.add_argument("--sample", type=int, default=20000, help="ClinVar names to test per build (0 for all)")
    validate.add_argument("--report", help="write the report to this Markdown file")
    _add_common(validate)
    validate.set_defaults(run=_cmd_validate)

    demo = commands.add_parser("demo", help="run the bundled demo with no download")
    demo.add_argument("--out", default="demo_output", help="folder for the demo results")
    demo.add_argument("--demo-dir", help="location of the demo folder if it is not beside the code")
    demo.add_argument("--force", action="store_true", help="replace existing results")
    demo.set_defaults(run=_cmd_demo, verbose=False)

    data = commands.add_parser("data", help="update or roll back the local data")
    data_commands = data.add_subparsers(title="actions", metavar="action")
    update = data_commands.add_parser("update", help="apply new source versions")
    update.add_argument("--check", action="store_true", help="only report what is new")
    update.add_argument("--approve", action="append", default=[], metavar="NAME", help="approve an update that can change results, or 'all'")
    _add_common(update)
    update.set_defaults(run=_cmd_update)
    rollback = data_commands.add_parser("rollback", help="go back to the previous data")
    rollback.add_argument("--to", help="snapshot name; the newest by default")
    _add_common(rollback)
    rollback.set_defaults(run=_cmd_rollback)
    return parser


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-dir", help="data folder (default: $VARSIGHT_DATA, the saved setting or ~/varsight-data)")
    parser.add_argument("--verbose", action="store_true")


# annotate
def _cmd_annotate(args) -> int:
    store = DataStore(resolve_data_dir(args.data_dir))
    try:
        _warn_if_stale(store)
        annotator = Annotator(store, all_transcripts=args.all_transcripts)
        return _run_annotation(annotator, args)
    finally:
        store.close()


def _run_annotation(annotator: Annotator, args) -> int:
    fmt = format_for(args.out, args.format)
    source = args.input
    if looks_like_vcf(source):
        writer = ResultWriter(args.out, fmt, vcf_columns(Path(source)), args.force)
        counts = annotate_vcf(annotator, Path(source), writer, args.build)
        writer.close()
        _report_batch(counts, args.out)
        return 0
    if Path(source).is_file():
        writer = ResultWriter(args.out, fmt, force=args.force)
        counts = annotate_list(annotator, Path(source), writer, args.build, args.gene)
        writer.close()
        _report_batch(counts, args.out)
        return 0
    records = annotator.annotate_text(source, args.build, args.gene)
    writer = ResultWriter(args.out, fmt, force=args.force)
    for record in records:
        writer.write(record)
    writer.close()
    if args.out:
        print(f"wrote {len(records)} result(s) to {args.out}", file=sys.stderr)
    return 0


def _report_batch(counts, out) -> None:
    summary = ", ".join(f"{name.replace('_', ' ')}: {value:,}" for name, value in counts.items())
    print(summary, file=sys.stderr)
    if out:
        print(f"wrote {out}", file=sys.stderr)


def _warn_if_stale(store: DataStore) -> None:
    updated = store.manifest.get("updated")
    if not updated:
        return
    try:
        age = datetime.now(timezone.utc) - datetime.strptime(updated, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except ValueError:
        return
    if age.days > _STALE_DAYS:
        log.warning("data was last updated %d days ago; run: varsight data update", age.days)


# setup and data
def _cmd_setup(args) -> int:
    from .pipeline.manager import DataManager

    data_dir = resolve_data_dir(args.data_dir)
    builds = ["hg19", "hg38"] if args.build == "both" else [args.build]
    if "alphamissense" in args.optional and "alphamissense" not in args.accept_license:
        print(
            "error: the AlphaMissense file header states CC BY-NC-SA 4.0 (non-commercial). Check it fits your use, "
            "then add: --accept-license alphamissense",
            file=sys.stderr,
        )
        return 2
    if is_windows_mount(data_dir):
        log.warning(
            "warning: %s is a Windows drive seen from WSL. Lookups measured about 50 times slower there "
            "than on the Linux disk. Consider: --data-dir ~/varsight-data",
            data_dir,
        )
    manager = DataManager(data_dir)
    if not args.yes and sys.stdin.isatty():
        print(f"VarSight will download public reference data into {data_dir}.")
        print("About 3.3 GB is saved and processed, and 4.7 GB of Pfam is read through without saving it.")
        print("It needs about 25 GB free during setup and about 12 GB afterwards (estimate).")
        if input("Continue? [y/N] ").strip().lower() not in ("y", "yes"):
            print("Nothing was changed.")
            return 1
    options = {"skip_pfam_regions": args.skip_pfam_regions, "gnomad_chroms": args.gnomad_chrom}
    manager.setup(builds, args.optional, options, args.force)
    save_data_dir(data_dir)
    print(f"Installed into {data_dir}. Running the self-test.", file=sys.stderr)
    args.data_dir = str(data_dir)
    return _cmd_selftest(args)


def _cmd_update(args) -> int:
    from .pipeline.manager import DataManager

    manager = DataManager(resolve_data_dir(args.data_dir))
    items = manager.update(args.approve, check_only=args.check)
    failed = False
    for item in items:
        label = item.action if item.action != "none" else "ok"
        detail = item.detail or ", ".join(item.changed)
        print(f"{item.artifact:24} {label:15} {detail}")
        failed = failed or item.action in ("failed", "error")
    return 1 if failed else 0


def _cmd_rollback(args) -> int:
    from .pipeline.manager import DataManager

    restored = DataManager(resolve_data_dir(args.data_dir)).rollback(args.to)
    print("restored: " + ", ".join(restored))
    return 0


def _cmd_status(args) -> int:
    data_dir = resolve_data_dir(args.data_dir)
    store = DataStore(data_dir)
    try:
        manifest = store.manifest
        print(f"data folder   {data_dir}")
        print(f"builds        {', '.join(manifest.get('builds', []))}")
        print(f"last updated  {manifest.get('updated')}")
        print(f"\n{'database':24} {'version':12} {'policy':10} {'size':>10}")
        for name, entry in sorted(manifest.get("files", {}).items()):
            size = sum(item.get("size", 0) for item in entry.get("files", {}).values())
            print(f"{name:24} {entry['source']['version']:12} {entry['policy']:10} {size / 1024**2:>7,.0f} MB")
    finally:
        store.close()
    return 0


def _cmd_doctor(args) -> int:
    import shutil
    import sqlite3

    data_dir = resolve_data_dir(args.data_dir)
    problems = 0
    version = sys.version_info
    ok = version >= (3, 10)
    problems += not ok
    print(f"{'ok' if ok else 'PROBLEM':8} Python {version.major}.{version.minor}.{version.micro} (3.10 or newer needed)")
    print(f"{'ok':8} SQLite {sqlite3.sqlite_version}")
    print(f"{'ok':8} data folder {data_dir}")
    if is_windows_mount(data_dir):
        print("WARNING  this is a Windows drive seen from WSL, so lookups will be much slower than on the Linux disk")
    probe = data_dir if data_dir.exists() else data_dir.parent
    if probe.exists():
        free = shutil.disk_usage(probe).free / 1024**3
        print(f"{'ok' if free > 30 else 'WARNING':8} {free:.1f} GB free (setup needs about 25 GB)")
    if not (data_dir / "manifest.json").is_file():
        print("not yet  no data installed; run: varsight setup")
        return 1 if problems else 0
    from .pipeline.manager import DataManager

    manager = DataManager(data_dir)
    manifest = manager.load_manifest()
    print(f"{'ok':8} installed builds: {', '.join(manifest.get('builds', []))}")
    for name in sorted(manifest.get("files", {})):
        for target in manager._artifact_paths(name):
            if not target.is_file():
                print(f"PROBLEM  {target} is missing")
                problems += 1
    if args.verify:
        for file_name, passed, detail in manager.verify():
            print(f"{'ok' if passed else 'PROBLEM':8} {file_name}: {detail}")
            problems += not passed
    return 1 if problems else 0


def _cmd_selftest(args) -> int:
    from .pipeline.selftest import run_selftest

    store = DataStore(resolve_data_dir(args.data_dir))
    try:
        results = run_selftest(store)
    finally:
        store.close()
    failed = 0
    for result in results:
        mark = "ok" if result["passed"] else "FAIL"
        print(f"{mark:5} {result['name']}: {result['detail']}")
        failed += not result["passed"]
    print(f"{len(results) - failed} of {len(results)} known variants passed")
    return 1 if failed or not results else 0


def _cmd_validate(args) -> int:
    from .pipeline.validate import render_report, run_validation

    store = DataStore(resolve_data_dir(args.data_dir))
    try:
        builds = [args.build] if args.build else store.builds
        results = run_validation(store, builds, args.sample)
    finally:
        store.close()
    report = render_report(results)
    if args.report:
        Path(args.report).write_text(report, encoding="utf-8")
        print(f"wrote {args.report}", file=sys.stderr)
    else:
        print(report)
    return 0


def _cmd_demo(args) -> int:
    from .pipeline.demo import run_demo

    return run_demo(args.demo_dir, Path(args.out), args.force)

