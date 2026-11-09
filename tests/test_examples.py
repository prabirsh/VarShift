"""varshift examples: every accepted form, each with an example that really works."""

import contextlib
import io
import os
import re
import unittest
from pathlib import Path
from unittest import mock

from tests.test_integration import DEMO, WithData
from varshift import cli
from varshift.batch import split_queries
from varshift.commandline import render_help
from varshift.examples import EXAMPLES, KINDS, menu_lines, render_examples
from varshift.inputs import parse_input
from varshift.term import QuietTerminal, Terminal

ROOT = DEMO.parent


def descriptions():
    """The examples that are descriptions, not files."""
    return [example for kind, example, _ in EXAMPLES if kind != "file"]


class EveryExampleWorks(WithData):
    def test_every_description_annotates_on_the_demo_data(self):
        for example in descriptions():
            for query in split_queries(example):
                with self.subTest(example=query):
                    self.assertTrue(self.annotator.annotate_text(query, "hg19"))

    def test_the_file_examples_are_the_demo_files(self):
        files = [example for kind, example, _ in EXAMPLES if kind == "file"]
        self.assertEqual(files, ["demo/variants.txt", "demo/demo.hg19.vcf"])
        for name in files:
            self.assertTrue((ROOT / name).is_file(), name)


class EveryFormIsShown(unittest.TestCase):
    def test_each_form_in_the_spec_has_an_example(self):
        parsed = {example: parse_input(example) for example in descriptions() if ";" not in example}
        kinds = [result.kind for result in parsed.values()]
        self.assertGreaterEqual(kinds.count("position"), 3, "positions in three spellings")
        self.assertIn("rsid", kinds)
        forms = {
            "one-letter protein": r"^[A-Z0-9]+ [A-Z]\d+[A-Z]$",
            "three-letter protein": r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$",
            "NP_ name": r"^NP_",
            "c. with a transcript": r"^NM_.*:c\.",
            "c. with a gene": r"^[A-Z0-9]+ c\.",
            "n.": r":n\.",
            "g. with NC_": r"^NC_0000.*:g\.",
            "g. with chr": r"^chr.*:g\.",
            "m.": r":m\.",
        }
        for form, pattern in forms.items():
            with self.subTest(form=form):
                self.assertTrue(any(re.search(pattern, example) for example in parsed), form)
        self.assertTrue(any(len(split_queries(example)) > 1 for example in descriptions()), "the ; list")

    def test_every_example_says_what_it_means(self):
        for kind, example, meaning in EXAMPLES:
            with self.subTest(example=example):
                self.assertIn(kind, KINDS)
                self.assertTrue(meaning and meaning[0].isupper() and not meaning.endswith("."), meaning)


class Layout(unittest.TestCase):
    def test_the_page_lists_every_example_and_how_to_use_them(self):
        text = render_examples(QuietTerminal())
        for _, example, meaning in EXAMPLES:
            self.assertIn(example, text)
        self.assertIn('varshift annotate "BRAF V600E" build=hg19', text)
        self.assertNotIn("\x1b", text)

    def test_lines_stay_inside_the_width(self):
        for width in (60, 100):
            for line in render_examples(QuietTerminal(width=width)).splitlines():
                self.assertLessEqual(len(line), width, line)

    def test_headings_are_bold_with_colour(self):
        self.assertIn("\x1b[1m", render_examples(Terminal(interactive=False, color=True, width=100)))

    def test_the_menu_has_five_kinds_numbered_from_one(self):
        lines = menu_lines()
        self.assertEqual(len(lines), 5)
        for number, line in enumerate(lines, 1):
            self.assertTrue(line.lstrip().startswith(str(number)), line)
        self.assertIn("BRAF V600E", lines[0])
        self.assertIn("rs113488022", lines[3])


class Command(unittest.TestCase):
    def _main(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = cli.main(list(args))
        return code, out.getvalue(), err.getvalue()

    def test_varshift_examples_needs_no_data(self):
        with mock.patch.dict(os.environ, {"VARSHIFT_DATA": str(Path(__file__).parent / "no such folder")}):
            code, out, _ = self._main("examples")
        self.assertEqual(code, 0)
        self.assertIn("NM_004333.6:c.1799T>A", out)

    def test_one_kind_can_be_shown_alone(self):
        text = render_examples(QuietTerminal(), "cdna")
        self.assertIn("NR_027676.2:n.200G>A", text)
        self.assertNotIn("BRAF V600E", text)
        self.assertNotIn("varshift annotate", text)

    def test_help_examples_works(self):
        self.assertIn("examples", render_help("examples"))
        self.assertEqual(self._main("help", "examples")[0], 0)

    def test_an_input_nothing_fits_points_at_the_examples(self):
        self.assertIn("varshift examples", cli.FORMAT_HINT)


if __name__ == "__main__":
    unittest.main()
