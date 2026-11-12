"""The guided start: type variants the way a colleague wrote them, see each answer, and save them.

The data stays open between queries, so after the first answer each one is immediate. Results stay in
memory until they are saved, and nothing about a query is written to disk unless the person saves.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .annotate import Annotator
from .answer import Answer, choice_table, decide, render_answer, summary_line
from .batch import annotate_list, annotate_vcf, split_queries, vcf_columns
from .errors import DataError, InputError, ReferenceMismatch, ResolveError, VarShiftError
from .examples import MENU_KINDS, menu_lines, render_examples
from .output import ReportContext, ReportFolder, default_folder, format_of, open_file_sink, render_card
from .pipeline.demo import build_demo_store, find_demo_dir
from .store import DataStore
from .suggest import apply_choice, suggest, typo_question
from .vcf import looks_like_vcf

DEMO_NOTICE = "Made with the 12-gene demo data: not for real results"
SHEETS = "Summary, Priority, All results, Links, Warnings explained, Data versions"
_MODES = {"tumour": "tumour", "tumor": "tumour", "inherited": "inherited", "off": None, "none": None}
_QUIT = {"quit", "exit", "q"}
_WORDS = """Words you can type at the prompt:
  save               Save the results to Excel: varshift-session-DATE-TIME.xlsx, here
  save NAME.xlsx     That file. .html, .tsv, .json and .jsonl write those formats
  save FOLDER/       report.html, results.xlsx and results.tsv in that folder
  list               The results so far, numbered
  a number           That result in full
  all N              Every candidate of result N in full, not only the best
  remove N, clear    Drop result N, or all of them
  mode tumour        Evidence and links for tumour work first; also inherited, or off
  examples           How to write a variant
  help               This list
  Enter, or quit     Quit; it asks first if there are results not saved
"""


def _now() -> datetime:
    return datetime.now()


@dataclass(eq=False)  # the same query typed twice gives two results, so equal content isn't the same result
class Item:
    answer: Answer  # what is shown for this result
    chosen: list[dict]  # what save writes for it
    saved: bool = False


class Session:
    def __init__(
        self,
        store_opener: Callable[[], DataStore],
        terminal,
        cwd: Path,
        build: str | None = None,
        mode: str | None = None,
        gene: str = "",
        setup: Callable[[], int] | None = None,
    ):
        """store_opener raises DataError when there is no data; setup, when given, runs varshift setup."""
        self.store_opener = store_opener
        self.term = terminal
        self.cwd = Path(cwd)
        self.build = build
        self.mode = mode
        self.gene = gene
        self.setup = setup
        self.items: list[Item] = []
        self.demo = False
        self.store: DataStore | None = None
        self.annotator: Annotator | None = None

    # the whole session
    def run(self) -> int:
        work = None
        try:
            try:
                store = self.store_opener()
            except DataError:
                store, work = self._no_data()
                if store is None:
                    return 1
            self.store = store
            try:
                self.annotator = Annotator(store)
                self._header()
                return self._loop()
            finally:
                store.close()
        finally:
            if work is not None:
                shutil.rmtree(work, ignore_errors=True)

    def _loop(self) -> int:
        while True:
            try:
                lines = self._read()
            except KeyboardInterrupt:
                count = len(self._unsaved())
                self._say(f"\nStopped. {count} result{' was' if count == 1 else 's were'} not saved.")
                return 130
            if lines is None or (len(lines) == 1 and lines[0].strip().lower() in {""} | _QUIT):
                return self._quit()
            self._dispatch(lines)

    def _read(self) -> list[str] | None:
        """One line, or several when they were pasted together."""
        first = self.term.ask("> ")
        if first is None:
            return None
        lines = [first]
        while self.term.input_waiting():
            more = self.term.ask("")
            if more is None:
                break
            lines.append(more)
        if len(lines) > 1:
            lines = [line for line in lines if line.strip()] or [""]
        return lines

    def _dispatch(self, lines: list[str]) -> None:
        if len(lines) > 1:
            self._queries([query for line in lines for query in split_queries(line)])
            return
        text = lines[0].strip()
        word, _, rest = text.partition(" ")
        word, rest = word.lower(), rest.strip()
        if text.isdigit():
            self._number(int(text))
        elif word == "save":
            self._save(rest or None)
        elif word == "list" and not rest:
            self._list()
        elif word == "remove":
            self._remove(rest)
        elif word == "clear" and not rest:
            self._clear()
        elif word == "help" and not rest:
            self._say(render_examples(self.term) + "\n" + _WORDS)
        elif word == "all":
            self._all(rest)
        elif word == "examples":
            self._examples(rest)
        elif word == "mode":
            self._set_mode(rest)
        elif self._is_file(text):
            self._file(self._path(text))
        else:
            self._queries(split_queries(text))

    # queries
    def _queries(self, queries: list[str]) -> None:
        added: list[Item] = []
        for text in queries:
            try:
                added += self._one(text)
            except KeyboardInterrupt:
                self._say("\nStopped that query.")
        if not added:
            return
        if len(queries) == 1:
            for item in added:
                self._card(item)
            count = len(self.items)
            self._say(f"{count} result{'' if count == 1 else 's'} in this session; save writes {'it' if count == 1 else 'them'} to Excel.")
            return
        self._summary(added)
        self._say(f"Type a number to see its full answer. {len(self.items)} results in this session; save writes them to Excel.")

    def _one(self, text: str) -> list[Item]:
        """Annotate one query, asking about a slip or a choice. Returns the results it added."""
        try:
            records = self.annotator.annotate_text(text, self.build, self.gene, self.mode)
        except (InputError, ResolveError, ReferenceMismatch) as exc:
            records = self._suggested(text, exc)
            if records is None:
                return []
        except VarShiftError as exc:
            self._say(f"{text}: could not annotate: {exc}")
            return []
        records = self._typo_checked(text, records)
        if records is None:
            return []
        answer = decide(records)
        if answer.kind != "choice":
            items = [Item(answer, [answer.best])]
        else:
            heading, rows = choice_table(answer, self.term)
            picked = self.term.choose(self.term.style(heading, "strong"), rows, allow_all=True)
            if picked is None:
                self._say(f"Skipped {text}.")
                return []
            items = [self._chosen(answer, answer.records[index]) for index in picked]
        self.items += items
        return items

    def _suggested(self, text: str, error: VarShiftError) -> list[dict] | None:
        question = suggest(self.store, text, error, self.build, gene=self.gene)
        if question is None or not question.options:
            self._say(f"{text}: could not annotate: {error}")
            self._say(question.text if question is not None else "Not sure of the format? Type examples.")
            return None
        width = max(len(option.label) for option in question.options)
        rows = [f"{option.label:<{width}}  {option.detail}" for option in question.options]
        picked = self.term.choose(self.term.style(question.text, "strong"), rows)
        if picked is None:
            self._say(f"Skipped {text}.")
            return None
        option = question.options[picked[0]]
        try:
            records = self.annotator.annotate_text(option.text, option.build or self.build, self.gene, self.mode)
        except VarShiftError as exc:
            self._say(f"{option.text}: could not annotate: {exc}")
            return None
        return apply_choice(records, text, option)

    def _typo_checked(self, text: str, records: list[dict]) -> list[dict] | None:
        """It fitted only a transcript other than the default, often a typo: ask which was meant."""
        question = typo_question(records)
        if question is None:
            return records
        width = max(len(option.label) for option in question.options)
        rows = [f"{option.label:<{width}}  {option.detail}" for option in question.options]
        picked = self.term.choose(self.term.style(question.text, "strong"), rows)
        if picked is None:
            self._say(f"Skipped {text}.")
            return None
        option = question.options[picked[0]]
        if option.keep:
            return records
        try:
            corrected = self.annotator.annotate_text(option.text, self.build, self.gene, self.mode)
        except VarShiftError as exc:
            self._say(f"{option.text}: could not annotate: {exc}")
            return None
        return apply_choice(corrected, text, option)

    @staticmethod
    def _chosen(answer: Answer, record: dict) -> Item:
        """One of several candidates, picked by the person, kept as a result of its own."""
        record["candidate"]["best"] = True
        record["candidate"]["why"] = ["chosen in the session"]
        return Item(Answer(answer.raw, answer.read_as, "single", [record], best=record), [record])

    # words
    def _number(self, number: int) -> None:
        if not self.items:
            if 1 <= number <= len(MENU_KINDS):
                self.term.write(render_examples(self.term, MENU_KINDS[number - 1]))
            else:
                self._say("Type a variant, or help for what else you can type.")
        elif 1 <= number <= len(self.items):
            self._card(self.items[number - 1])
        else:
            self._say(f"There is no result {number}. Type list to see them.")

    def _card(self, item: Item) -> None:
        number = next(n for n, held in enumerate(self.items, 1) if held is item)
        self.term.write(render_answer(item.answer, self.term, also_hint=f"all {number}"))

    def _all(self, rest: str) -> None:
        """Every candidate of one result in full, as show=all does for varshift annotate."""
        if not rest.isdigit():
            self._say("Type all and the number of a result, as in all 1.")
            return
        if not 1 <= int(rest) <= len(self.items):
            self._say(f"There is no result {rest}. Type list to see them.")
            return
        for record in self.items[int(rest) - 1].answer.records:
            self.term.write(render_card(record) + "\n")

    def _list(self) -> None:
        if not self.items:
            self._say("No results yet.")
            return
        self._summary(self.items)
        unsaved = len(self._unsaved())
        self._say(f"{len(self.items)} result{'' if len(self.items) == 1 else 's'}; {unsaved} not saved yet.")

    def _summary(self, items: list[Item]) -> None:
        width = min(max(len(item.answer.raw) for item in self.items), 30)
        numbers = {id(item): number for number, item in enumerate(self.items, 1)}
        for item in items:
            self._say(summary_line(numbers[id(item)], item.answer, self.term, width))

    def _remove(self, rest: str) -> None:
        if not rest.isdigit():
            self._say("Type remove and a number, as in remove 2.")
            return
        number = int(rest)
        if not 1 <= number <= len(self.items):
            self._say(f"There is no result {number}. Type list to see them.")
            return
        item = self.items.pop(number - 1)
        left = len(self.items)
        self._say(f"Removed {number}: {item.answer.raw}. {left} result{'' if left == 1 else 's'} left.")

    def _clear(self) -> None:
        count = len(self.items)
        self.items.clear()
        self._say(f"Cleared {count} result{'' if count == 1 else 's'}.")

    def _examples(self, rest: str) -> None:
        if rest.isdigit() and 1 <= int(rest) <= len(MENU_KINDS):
            self.term.write(render_examples(self.term, MENU_KINDS[int(rest) - 1]))
        else:
            self.term.write(render_examples(self.term))

    def _set_mode(self, rest: str) -> None:
        if rest.lower() not in _MODES:
            self._say(f"mode must be tumour, inherited or off, not {rest!r}")
            return
        self.mode = _MODES[rest.lower()]
        self._say(f"Mode: {self.mode or 'none'}. It applies to the queries you type from now on.")

    # files
    def _is_file(self, text: str) -> bool:
        try:
            return self._path(text).is_file()
        except (OSError, ValueError):
            return False

    def _path(self, text: str) -> Path:
        return self.cwd / text.strip().strip('"').strip("'")

    def _file(self, path: Path) -> None:
        """A file typed at the prompt is annotated as varshift annotate FILE would, into its own folder."""
        folder = self.cwd / default_folder(path.name)
        replace = False
        if folder.is_dir() and any(folder.iterdir()):
            if not self.term.confirm(f"{folder} already has results. Replace them?"):
                self._say("Left as it is.")
                return
            replace = True
        context = self._context(path.name)
        try:
            if looks_like_vcf(str(path)):
                context.prefix = vcf_columns(path)
                sink = ReportFolder(folder, context, replace)
                counts = annotate_vcf(self.annotator, path, sink, self.build, 1, self.mode)
            else:
                sink = ReportFolder(folder, context, replace)
                counts = annotate_list(self.annotator, path, sink, self.build, self.gene, self.mode)
            sink.close()
        except (VarShiftError, OSError) as exc:
            self._say(f"{path.name}: {exc}")
            return
        summary = ", ".join(f"{name.replace('_', ' ')}: {value:,}" for name, value in counts.items())
        self._say(f"{path.name}: {summary}")
        self._say(f"Open this to explore the results: {sink.paths[0]}")

    # saving and quitting
    def _save(self, target: str | None) -> None:
        if not self.items:
            self._say("Nothing to save yet.")
            return
        context = self._context("VarShift session")
        try:
            if target is None:
                shown = self._free_name()
                sink = open_file_sink(shown, context)
            elif target.endswith(("/", "\\")) or (self.cwd / target).is_dir():
                shown = self.cwd / target
                replace = shown.is_dir() and any(shown.iterdir())
                if replace and not self.term.confirm(f"{shown} already has files. Replace them?"):
                    self._say("Not saved.")
                    return
                sink = ReportFolder(shown, context, replace)
            else:
                shown = self.cwd / target
                if format_of(str(shown)) is None:
                    self._say("save writes .xlsx, .html, .tsv, .json or .jsonl files, or a folder ending in /")
                    return
                replace = shown.exists()
                if replace and not self.term.confirm(f"{shown} exists. Replace it?"):
                    self._say("Not saved.")
                    return
                sink = open_file_sink(shown, context, replace)
            for item in self.items:
                for record in item.chosen:
                    sink.write(record)
            sink.close()
        except (OSError, VarShiftError) as exc:
            self._say(f"Could not save: {exc}. The results are kept.")
            return
        for item in self.items:
            item.saved = True
        count = sum(len(item.chosen) for item in self.items)
        self._say(f"Saved {count} result{'' if count == 1 else 's'} to {shown}")
        if str(shown).lower().endswith(".xlsx"):
            links = sum(len(record.get("links") or []) for item in self.items for record in item.chosen)
            self._say(f"Sheets: {SHEETS}. Links: {links:,}.")

    def _free_name(self) -> Path:
        """varshift-session-DATE-TIME.xlsx here, with -2, -3 and so on: a save never overwrites."""
        base = f"varshift-session-{_now():%Y-%m-%d-%H%M}"
        path, number = self.cwd / f"{base}.xlsx", 2
        while path.exists():
            path, number = self.cwd / f"{base}-{number}.xlsx", number + 1
        return path

    def _quit(self) -> int:
        count = len(self._unsaved())
        if count and self.term.confirm(f"Save the {count} result{'' if count == 1 else 's'} to Excel first?"):
            self._save(None)
        return 0

    def _unsaved(self) -> list[Item]:
        return [item for item in self.items if not item.saved]

    def _context(self, name: str) -> ReportContext:
        return ReportContext(
            input_name=name,
            installed=self.store.installed(),
            versions=self.store.versions(),
            notices=[DEMO_NOTICE] if self.demo else [],
        )

    # the start
    def _header(self) -> None:
        if self.demo:
            data = "12-gene demo data"
        else:
            versions = self.store.versions()
            clinvar = next((versions[f"clinvar.{build}"] for build in self.store.builds if versions.get(f"clinvar.{build}")), None)
            data = ", ".join(self.store.builds) + (f", ClinVar {clinvar}" if clinvar else "")
        lines = [
            f"VarShift: offline variant annotation. Data: {data}",
            "",
            "What do you have? Type one or more, separated by ;  (or paste several lines)",
            *menu_lines(),
            "Type a number for examples of that kind, help for what else you can type, or Enter to quit.",
            "",
        ]
        self._say("\n".join(lines))

    def _no_data(self) -> tuple[DataStore | None, Path | None]:
        self._say("VarShift: offline variant annotation. No data is set up here yet.")
        options = ["Try the 12-gene demo data now (no download)", "Set up hg19, hg38 or both"]
        picked = self.term.choose(self.term.style("What would you like to do?", "strong"), options)
        if picked is None:
            self._say("Nothing set up. When you are ready: varshift setup")
            return None, None
        if picked[0] == 0:
            work = Path(tempfile.mkdtemp(prefix="varshift-demo-session-", dir=self.cwd))
            try:
                self._say("Building the demo data set: a few seconds, and nothing is downloaded. It covers 12 genes only.")
                build_demo_store(find_demo_dir(None), work)
                store = DataStore(work)
            except BaseException:
                shutil.rmtree(work, ignore_errors=True)
                raise
            self.demo = True
            return store, work
        if self.setup is None or self.setup() != 0:
            self._say("Setup didn't finish. Run it again any time: varshift setup")
            return None, None
        try:
            return self.store_opener(), None
        except DataError as exc:
            self._say(f"error: {exc}")
            return None, None

    def _say(self, text: str) -> None:
        self.term.write(text + "\n")
