"""Reads the command line.

Options are written as words, the way make and dd do it: build=hg19 out=results.xlsx. That reads
well for people who don't live in a terminal. Anyone used to dashed options can still type them, since
leading dashes are ignored.
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass, field

from .errors import InputError


@dataclass(frozen=True)
class Command:
    summary: str
    positional: str = ""
    positional_help: str = ""
    required: bool = False
    options: dict[str, str] = field(default_factory=dict)
    flags: dict[str, str] = field(default_factory=dict)
    examples: tuple[str, ...] = ()


DATA_OPTION = {"data": "data folder (default: the folder chosen at setup, or ~/varsight-data)"}

COMMANDS: dict[str, Command] = {
    "annotate": Command(
        "Annotate one variant, a list file or a VCF",
        "INPUT",
        "a description such as 'BRAF V600E', or the path of a .vcf, .vcf.gz or .txt file",
        required=True,
        options={
            "build": "hg19 or hg38; needed for a bare genomic position, read from the header for a VCF",
            "gene": "gene to use when the description has none",
            "out": "output file (.xlsx .html .tsv .json .jsonl) or folder. A file input gets a folder with all reports",
            "format": "card, json, jsonl or tsv, for output on the screen",
            "transcripts": "best (the top five per gene, the default) or all",
            **DATA_OPTION,
        },
        flags={"replace": "replace output that already exists"},
        examples=(
            'varsight annotate "BRAF V600E" build=hg19',
            'varsight annotate "NM_004333.6:c.1799T>A" build=hg38',
            "varsight annotate sample.vcf.gz",
            "varsight annotate variants.txt build=hg19 out=results.xlsx",
        ),
    ),
    "setup": Command(
        "Choose hg19, hg38 or both, then download and build the databases",
        options={
            "build": "hg19, hg38 or both. If left out, you are asked",
            "with": "advanced databases to add, separated by commas. See: varsight databases",
            "accept": "advanced databases whose licence you accept, separated by commas",
            **DATA_OPTION,
        },
        flags={
            "yes": "don't ask anything; use the options given",
            "rebuild": "rebuild even what is already installed",
            "nopfam": "skip the 4.7 GB Pfam read; domain positions will be missing",
        },
        examples=("varsight setup", "varsight setup build=hg38 data=/data/varsight yes"),
    ),
    "databases": Command("List every database: what is installed, its size, licence and how to add it", options=DATA_OPTION),
    "add": Command(
        "Add an advanced database",
        "NAME",
        "a name from: varsight databases",
        required=True,
        options={
            "file": "your own licensed download, for databases that need an account (COSMIC, OMIM, SpliceAI, PrimateAI)",
            "scope": "coding (the default: coding exons and 50 bases around them, far less disk) or all",
            "chrom": "limit to these chromosomes, separated by commas (gnomAD only)",
            **DATA_OPTION,
        },
        flags={"accept": "confirm the licence shown fits your use", "yes": "don't ask for confirmation"},
        examples=("varsight add dbsnp", "varsight add cosmic file=Cosmic_GenomeScreensMutant_v101_GRCh37.tsv.gz accept"),
    ),
    "remove": Command("Remove an advanced database and free its disk space", "NAME", "a name from: varsight databases", required=True, options=DATA_OPTION),
    "update": Command(
        "Bring the databases up to date",
        options={"approve": "also apply updates that can change names, such as transcripts (or: all)", **DATA_OPTION},
        flags={"check": "only report what is new; change nothing"},
        examples=("varsight update check", "varsight update", "varsight update approve=transcripts"),
    ),
    "rollback": Command("Go back to the previous data", options={"to": "snapshot name; the newest by default", **DATA_OPTION}),
    "status": Command("Show what is installed, how old it is, and whether updates are waiting", options=DATA_OPTION),
    "doctor": Command("Check this computer and the data folder", options=DATA_OPTION, flags={"verify": "also compare every data file with its recorded checksum"}),
    "explain": Command(
        "Show each step for one variant: which file was read, which query ran, how long it took",
        "INPUT",
        "a description such as 'BRAF V600E'",
        required=True,
        options={"build": "hg19 or hg38", "gene": "gene to use when the description has none", **DATA_OPTION},
        examples=('varsight explain "BRAF V600E" build=hg19',),
    ),
    "selftest": Command("Run the fixed list of known variants", options=DATA_OPTION),
    "validate": Command(
        "Measure accuracy against ClinVar's own names",
        options={"build": "hg19 or hg38; both by default", "sample": "names to test per build, 0 for all (default 20000)", "report": "write the report to this Markdown file", **DATA_OPTION},
    ),
    "benchmark": Command(
        "Measure speed and memory on this computer",
        options={"records": "VCF records to test with (default 3000)", "report": "write the result to this Markdown file", **DATA_OPTION},
    ),
    "demo": Command("Run the bundled demo with no download", options={"out": "folder for the results (default demo_output)", "demodir": "the demo folder, if it is not beside the code"}, flags={"replace": "replace existing results"}),
    "config": Command("Show or change settings", options={"updates": "on or off: the quick check for new database versions", "data": "the data folder to use from now on"}),
    "help": Command("Show help for a command", "COMMAND", "the command to explain"),
    "version": Command("Show the version"),
}


class Values:
    """The parsed options of one command."""

    def __init__(self, command: str):
        self.command = command
        self.positional: str | None = None
        self.options: dict[str, str] = {}
        self.flags: set[str] = set()

    def get(self, name: str, default: str | None = None) -> str | None:
        return self.options.get(name, default)

    def has(self, flag: str) -> bool:
        return flag in self.flags

    def items(self, name: str) -> list[str]:
        """A comma-separated option as a list."""
        return [item.strip() for item in (self.get(name) or "").split(",") if item.strip()]


def parse(argv: list[str]) -> Values:
    if not argv:
        return Values("help")
    first = argv[0].lstrip("-").lower()
    if first in ("h", "help"):
        values = Values("help")
        values.positional = argv[1].lstrip("-").lower() if len(argv) > 1 else None
        return values
    if first in ("v", "version"):
        return Values("version")
    if first not in COMMANDS:
        raise InputError(f"unknown command {argv[0]!r}{_suggest(first, COMMANDS)}. See: varsight help")
    command = COMMANDS[first]
    values = Values(first)
    tokens = argv[1:]
    index = 0
    while index < len(tokens):
        token = tokens[index]
        index += 1
        word = token.lstrip("-")
        if word.lower() in ("h", "help") and word != token:
            values.command = "help"
            values.positional = first
            return values
        key, equals, value = word.partition("=")
        key = key.lower()
        if equals and key in command.options:
            values.options[key] = value
        elif not equals and key in command.flags:
            values.flags.add(key)
        elif word != token and not equals and key in command.options:
            # Dashed habit: -build hg19 as well as build=hg19.
            if index >= len(tokens):
                raise InputError(f"{key} needs a value, as in {key}=VALUE")
            values.options[key] = tokens[index]
            index += 1
        elif command.positional and values.positional is None:
            values.positional = token
        else:
            known = list(command.options) + list(command.flags)
            raise InputError(f"{first} does not understand {token!r}{_suggest(key, known)}. See: varsight help {first}")
    if command.required and values.positional is None:
        raise InputError(f"{first} needs {command.positional}: {command.positional_help}. See: varsight help {first}")
    return values


def _suggest(word: str, choices) -> str:
    close = difflib.get_close_matches(word, list(choices), n=1, cutoff=0.6)
    return f"; did you mean {close[0]}?" if close else ""


def render_help(name: str | None = None) -> str:
    if name is None:
        lines = ["VarSight: offline variant annotation for hg19 and hg38.", "", "Usage: varsight COMMAND [INPUT] [option=value ...]", "", "Commands:"]
        lines += [f"  {command:11} {info.summary}" for command, info in COMMANDS.items()]
        lines += ["", 'Try: varsight annotate "BRAF V600E" build=hg19', "More: varsight help annotate"]
        return "\n".join(lines)
    if name not in COMMANDS:
        raise InputError(f"unknown command {name!r}{_suggest(name, COMMANDS)}")
    info = COMMANDS[name]
    usage = f"varsight {name}" + (f" {info.positional}" if info.positional else "") + (" [option=value ...]" if info.options else "")
    lines = [info.summary, "", f"Usage: {usage}"]
    if info.positional:
        lines += ["", f"  {info.positional:12} {info.positional_help}"]
    if info.options:
        lines += ["", "Options:"] + [f"  {key + '=':12} {text}" for key, text in info.options.items()]
    if info.flags:
        lines += ["", "Words you can add:"] + [f"  {key:12} {text}" for key, text in info.flags.items()]
    if info.examples:
        lines += ["", "Examples:"] + [f"  {example}" for example in info.examples]
    return "\n".join(lines)
