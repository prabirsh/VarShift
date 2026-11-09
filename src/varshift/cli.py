"""The varshift command."""

from __future__ import annotations

import errno
import logging
import os
import shutil
import sqlite3
import sys
from pathlib import Path

from . import __version__
from .advanced import ACCEPT, ALIAS, DATABASES, DOWNLOAD, OWN_FILE, get_database
from .annotate import Annotator
from .answer import Answer, choice_table, decide, error_line, render_answer, render_choices, summary_line
from .batch import annotate_lines, annotate_list, annotate_vcf, split_queries, vcf_columns
from .builds import parse_build
from .commandline import COMMANDS, Values, parse, render_help
from .config import default_data_dir, forget_setting, is_windows_mount, load_settings, resolve_data_dir, save_setting
from .errors import DataError, InputError, ReferenceMismatch, ResolveError, SetupError, VarShiftError
from .examples import render_examples
from .output import ReportContext, ReportFolder, ScreenSink, default_folder, format_of, open_file_sink, render_card
from .session import Session
from .store import DataStore
from .suggest import apply_choice, option_texts, suggest
from .term import SETTINGS, Terminal
from .updates import check_in_background, clear_applied, notice_lines, read_status
from .vcf import looks_like_vcf

log = logging.getLogger("varshift")
_OUT_OF_MEMORY = (
    "The system ran out of free memory. Close other large programs and run the same command again: finished steps are kept. "
    "In WSL, this happens when reading a Windows drive under /mnt, so keep the data on the Linux disk, or use PowerShell."
)
_PARALLEL_FROM = 1024 * 1024


def main(argv: list[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        # A Windows console can't show every character. Replace those rather than stop.
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)
    try:
        arguments = sys.argv[1:] if argv is None else argv
        if not arguments and _at_a_terminal():
            return _start(Values("start"))
        values = parse(arguments)
        return HANDLERS[values.command](values)
    except VarShiftError as exc:
        _report(exc)
        return 1
    except KeyboardInterrupt:
        print("stopped", file=sys.stderr)
        return 130
    except BrokenPipeError:
        return 0
    except OSError as exc:
        print(f"error: {exc}", file=sys.stderr)
        if exc.errno == errno.ENOMEM:
            print(_OUT_OF_MEMORY, file=sys.stderr)
        return 1


def _report(exc: VarShiftError) -> None:
    print(f"error: {exc}", file=sys.stderr)
    for line in getattr(exc, "tried", []):
        print(f"  tried {line}", file=sys.stderr)


def _help(values: Values) -> int:
    print(render_help(values.positional))
    return 0


def _version(values: Values) -> int:
    print(f"VarShift {__version__}")
    return 0


def _examples(values: Values) -> int:
    term = _terminal(values)
    term.write(render_examples(term))
    return 0


def _at_a_terminal() -> bool:
    """For varshift with nothing after it: open the guided start only when a person is there."""
    try:
        return _terminal(Values("start")).interactive
    except InputError:
        return False  # a bad saved colour; the command list still shows, and start would say what to fix


def _start(values: Values) -> int:
    term = _terminal(values)
    if not term.interactive:
        print("varshift start needs a terminal. For scripts and pipes, use: varshift annotate", file=sys.stderr)
        return 1
    mode = _mode(values)
    build = parse_build(values.get("build")) if values.get("build") else None
    data_dir = resolve_data_dir(values.get("data"))
    try:
        store = DataStore(data_dir)
    except DataError:
        pass  # the session offers the demo data or setup
    else:
        try:
            for line in notice_lines(data_dir, store.manifest):
                term.write(line + "\n")
            check_in_background(data_dir, store.manifest)
        finally:
            store.close()
    session = Session(
        # Opened again after setup, which may have chosen another folder.
        lambda: DataStore(resolve_data_dir(values.get("data"))),
        term,
        Path.cwd(),
        build=build,
        mode=mode,
        setup=lambda: _setup(Values("setup")),
    )
    return session.run()


# annotate
def _annotate(values: Values) -> int:
    source = values.positional
    mode = _mode(values)
    for name in ("show", "links"):
        if values.get(name) not in (None, "all"):
            raise InputError(f"{name}= takes all, not {values.get(name)!r}")
    term = _terminal(values)
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
        build, gene = values.get("build"), values.get("gene") or ""
        if looks_like_vcf(source):
            context.prefix = vcf_columns(Path(source))
            sink = _sink_for_file(source, values, context)
            counts = annotate_vcf(annotator, Path(source), sink, build, _worker_count(values, source), mode)
            return _finish_file(sink, counts)
        if Path(source).is_file():
            sink = _sink_for_file(source, values, context)
            counts = annotate_list(annotator, Path(source), sink, build, gene, mode)
            return _finish_file(sink, counts)
        on_screen = not values.get("out") and (values.get("format") or "card") == "card"
        queries = split_queries(source)
        if len(queries) > 1 and on_screen:
            return _show_queries(annotator, queries, values, mode, term)
        if len(queries) > 1:
            # Several queries are handled as a list file with these lines, so one failing doesn't stop the rest.
            sink = _sink_for_single(values, context)
            counts = annotate_lines(annotator, enumerate(queries, 1), sink, build, gene, mode, label="query")
            return _finish_file(sink, counts)
        text = queries[0] if queries else source
        try:
            records = annotator.annotate_text(text, build, gene, mode)
        except (InputError, ResolveError, ReferenceMismatch) as exc:
            records = _suggest_instead(annotator, text, exc, values, mode, term, on_screen)
            if records is None:
                return 1
        if on_screen:
            return _show_answer(records, values, mode, term)
        sink = _sink_for_single(values, context)
        for record in records:
            sink.write(record)
        sink.close()
        if hasattr(sink, "final_path"):
            print(f"Wrote {len(records)} result(s) to {sink.final_path}", file=sys.stderr)
        return 0
    finally:
        store.close()


_MODES = {"tumour": "tumour", "tumor": "tumour", "inherited": "inherited", "none": None}


def _mode(values: Values) -> str | None:
    """mode=, else the saved setting. None when there is no mode."""
    if values.get("mode"):
        return _mode_named(values.get("mode"))
    saved = load_settings().get("mode") or "none"
    try:
        return _mode_named(saved)
    except InputError:
        raise InputError(f"the saved mode {saved!r} isn't tumour, inherited or none; to clear it: varshift config mode=none") from None


def _mode_named(word: str) -> str | None:
    if word.lower() not in _MODES:
        raise InputError(f"mode must be tumour, inherited or none, not {word!r}")
    return _MODES[word.lower()]


def _terminal(values: Values) -> Terminal:
    """The screen VarShift writes to, and asks on only when a person is there. Tests patch this."""
    given = values.get("color")
    setting = given or load_settings().get("color") or "auto"
    try:
        return Terminal.detect(setting)
    except ValueError:
        if given:
            raise InputError(f"color must be auto, yes or no, not {setting!r}") from None
        raise InputError(f"the saved colour {setting!r} isn't auto, yes or no; to reset it: varshift config color=auto") from None


FORMAT_HINT = "Not sure of the format? See: varshift examples"


def _suggest_instead(annotator: Annotator, text: str, error: VarShiftError, values: Values, mode, term: Terminal, on_screen: bool):
    """A description failed. Offer what it probably meant: ask when someone can answer, else say it.

    Returns the records of the chosen suggestion, or None when nothing was annotated (the reason is printed).
    """
    build, gene = values.get("build"), values.get("gene") or ""
    question = suggest(annotator.store, text, error, build, gene=gene)
    if not (on_screen and term.interactive and question is not None and question.options):
        _report(error)
        if question is not None and question.options:
            print(f"Did you mean: {', '.join(option_texts(question)[:5])}", file=sys.stderr)
        else:
            print(question.text if question is not None else FORMAT_HINT, file=sys.stderr)
        return None
    width = max(len(option.label) for option in question.options)
    picked = term.choose(term.style(question.text, "strong"), [f"{o.label:<{width}}  {o.detail}" for o in question.options])
    if picked is None:
        print("Nothing annotated.", file=sys.stderr)
        return None
    option = question.options[picked[0]]
    # Annotated once, as if it had been typed. If it fails too, that error is shown as it is.
    return apply_choice(annotator.annotate_text(option.text, option.build or build, gene, mode), text, option)


def _show_answer(records: list[dict], values: Values, mode: str | None, term: Terminal) -> int:
    """One description on the screen: the best match's card, or a choice."""
    if values.get("show") == "all":
        for record in records:
            term.write(render_card(record) + "\n")
        return 0
    answer = decide(records)
    links_all = values.get("links") == "all"
    if answer.kind != "choice":
        term.write(render_answer(answer, term, mode, links_all))
        return 0
    if not term.interactive:
        term.write(render_choices(answer, term))
        term.write("Each is a different variant. To see all in full: show=all\n")
        return 0
    heading, rows = choice_table(answer, term)
    if len(rows) < len(answer.records):
        heading += f"\n(the first {len(rows)} of {len(answer.records)}; a more specific description narrows them)"
    picked = term.choose(heading, rows, allow_all=True)
    if picked is None:
        print("Nothing chosen.", file=sys.stderr)
        return 1
    for index in picked:
        record = answer.records[index]
        chosen = Answer(answer.raw, answer.read_as, "single", [record], best=record)
        term.write("\n" + render_answer(chosen, term, mode, links_all))
    return 0


def _show_queries(annotator: Annotator, queries: list[str], values: Values, mode: str | None, term: Terminal) -> int:
    """Several descriptions on the screen: one numbered line each, then how to see one in full."""
    build, gene = values.get("build"), values.get("gene") or ""
    width = min(max(len(query) for query in queries), 30)
    example = None
    for number, text in enumerate(queries, 1):
        try:
            answer = decide(annotator.annotate_text(text, build, gene, mode))
        except VarShiftError as exc:
            texts = option_texts(suggest(annotator.store, text, exc, build, fast_only=True, gene=gene))
            message = f"did you mean {texts[0]}? ({exc})" if texts else str(exc)
            term.write(error_line(number, text, message, term, width) + "\n")
            continue
        example = example or text
        term.write(summary_line(number, answer, term, width) + "\n")
    words = "".join(f" {name}={values.get(name)}" for name in ("build", "gene", "mode") if values.get(name))
    term.write(f'To see one in full: varshift annotate "{example or queries[0]}"{words}\n')
    term.write("To save them all, with every link: add out=results.xlsx\n")
    return 0


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
    return ReportFolder(out or default_folder(source), context, values.has("replace"))


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
    "hg19": ("about 1.7 GB", "about 5 GB", "about 15 GB", "about 45 minutes"),
    "hg38": ("about 1.9 GB", "about 5 GB", "about 15 GB", "about 45 minutes"),
    "both": ("about 3.0 GB", "about 10 GB", "about 25 GB", "about 1.5 hours"),
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
    print("Databases that need your own account (COSMIC, OMIM, SpliceAI, PrimateAI) are added later: varshift databases")
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


def _refuse_windows_drive(values: Values, data_dir: Path, instead: str) -> None:
    """Setup, update and add build databases, which from WSL on a Windows drive can run out of memory."""
    if is_windows_mount(data_dir) and not values.has("windowsdrive"):
        raise InputError(
            f"{data_dir} is a Windows drive seen from WSL. Building databases there can stop with 'Cannot allocate memory', "
            f"and measured about 30 times slower. {instead} To build here anyway, add the word windowsdrive"
        )


def _setup(values: Values) -> int:
    from .pipeline.manager import DataManager

    interactive = sys.stdin.isatty() and not values.has("yes")
    print("VarShift setup")
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
    _refuse_windows_drive(
        values,
        data_dir,
        f"Use a folder on the Linux disk, such as data={default_data_dir()}, or run setup from PowerShell, after which WSL can read the same folder.",
    )
    with_names = values.items("with")
    if interactive and not with_names:
        with_names = _ask_advanced()
    accepted = set(values.items("accept"))
    for name in with_names:
        database = get_database(name)
        if database is None:
            raise InputError(f"no advanced database named {name}; see: varshift databases")
        if database.access == OWN_FILE:
            raise InputError(f"{database.title} needs your own file. After setup, run: varshift add {database.name} file=PATH accept")
        if database.access == ACCEPT and database.name not in accepted:
            print(f"\n{database.title}: {database.license}")
            if interactive and _ask("Does that fit your use? yes or no:", "no").lower() in ("y", "yes"):
                accepted.add(database.name)
            else:
                raise InputError(f"{database.title} needs you to accept its terms; add accept={database.name} if they fit your use")
    transfer, kept, free, duration = _CORE_NEEDS[build]
    print(f"\nThis will install: {'hg19 and hg38' if build == 'both' else build} into {data_dir}")
    print(f"  Download   {transfer}")
    print(f"  Disk       {kept} kept, and keep {free} free during setup")
    print(f"  Time       {duration} on a 2 MB/s connection. If it stops, run it again and it carries on")
    if with_names:
        print(f"  Advanced   {', '.join(with_names)}")
    if interactive and _ask("Start? yes or no:", "yes").lower() not in ("y", "yes"):
        print("Nothing was changed.")
        return 1
    builds = ["hg19", "hg38"] if build == "both" else [build]
    manager = DataManager(data_dir)
    # Saved first, so that running setup again after an interruption offers this same folder.
    save_setting("data_dir", str(data_dir))
    manager.setup(builds, {}, values.has("rebuild"))
    for name in with_names:
        manager.add(name, {"accepted": True})
    print(f"\nInstalled into {data_dir}. Checking it with known variants.", file=sys.stderr)
    values.options["data"] = str(data_dir)
    return _selftest(values)


def _add(values: Values) -> int:
    from .pipeline.manager import DataManager

    data_dir = resolve_data_dir(values.get("data"))
    database = get_database(values.positional)
    if database is None:
        raise InputError(f"no advanced database named {values.positional}; see: varshift databases")
    if database.access == ALIAS:
        print(f"{database.title}: {database.provides}.")
        print(f"To get it, run: varshift add {database.alias_of}")
        return 0
    _refuse_windows_drive(values, data_dir, f"Run varshift add {database.name} from PowerShell instead, which uses the same folder.")
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
    except VarShiftError:
        print(f"Nothing is installed in {data_dir} yet. Start with: varshift setup\n")
    print("Core databases, installed by setup: UCSC reference, RefSeq, GENCODE, MANE, ClinVar, HGNC, UniProt, Pfam, CIViC,")
    print("ClinGen gene curation, gnomAD gene constraint, HPO, and, when both builds are set up, the hg19 to hg38 mapping.\n")
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
            print(f"  Then run: varshift add {database.name} file=PATH accept")
        elif database.access == ACCEPT:
            print(f"  To add: varshift add {database.name} accept")
        elif database.access == DOWNLOAD:
            print(f"  To add: varshift add {database.name}")
        print()
    return 0


# update and rollback
def _update(values: Values) -> int:
    from .pipeline.manager import DataManager

    data_dir = resolve_data_dir(values.get("data"))
    if not values.has("check"):
        _refuse_windows_drive(
            values, data_dir, "Run varshift update from PowerShell instead, which uses the same folder. varshift update check works here, as it changes nothing."
        )
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
    if values.has("check"):
        print(f"\n{waiting} new, nothing changed.")
    else:
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
        print("Not yet  No data installed. To install: varshift setup")
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

    mode = _mode(values)
    data_dir = resolve_data_dir(values.get("data"))
    store = DataStore(data_dir)
    try:
        return run_explain(store, values.positional, values.get("build"), values.get("gene") or "", mode)
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
    # Every value is checked before any is saved, so a mistake changes nothing.
    updates, mode, color = values.get("updates"), values.get("mode"), values.get("color")
    if updates is not None and updates not in ("on", "off"):
        raise InputError("updates must be on or off")
    mode_word = _mode_named(mode) if mode is not None else None
    if color is not None and color.lower() not in SETTINGS:
        raise InputError(f"color must be auto, yes or no, not {color!r}")
    if updates is not None:
        save_setting("updates", updates)
    if values.get("data"):
        save_setting("data_dir", str(resolve_data_dir(values.get("data"))))
    if mode is not None:
        if mode_word is None:
            forget_setting("mode")
        else:
            save_setting("mode", mode_word)
    if color is not None:
        save_setting("color", color.lower())
    settings = load_settings()
    print(f"Data folder            {settings.get('data_dir', default_data_dir())}")
    print(f"Quick update check     {settings.get('updates', 'on')}")
    print(f"Mode                   {settings.get('mode') or 'none'}")
    print(f"Colour                 {settings.get('color') or 'auto'}")
    return 0


HANDLERS = {
    "annotate": _annotate,
    "start": _start,
    "examples": _examples,
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
