"""The documentation stays honest: every warning is explained, and links and commands really exist."""

import re
import shlex
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "varsight"
DOCS = sorted(list(ROOT.glob("*.md")) + list((ROOT / "docs").glob("*.md")) + list((ROOT / "demo").glob("*.md")))
_LINK_RE = re.compile(r"\]\(([^)#\s]+)(#[^)\s]*)?\)")
_HEADING_RE = re.compile(r"^#{1,6} (.+)$", re.MULTILINE)


def anchors_of(path: Path) -> set[str]:
    found = set()
    for heading in _HEADING_RE.findall(path.read_text(encoding="utf-8")):
        slug = re.sub(r"[^\w\- ]", "", heading.lower()).replace(" ", "-")
        found.add(slug)
    return found


class Flags(unittest.TestCase):
    def test_every_flag_in_the_code_is_documented(self):
        documented = (ROOT / "docs" / "input-output.md").read_text(encoding="utf-8")
        found = set()
        for path in SRC.rglob("*.py"):
            text = path.read_text(encoding="utf-8")
            found.update(re.findall(r'flags\.append\("([a-zA-Z0-9_]+)"\)', text))
            found.update(re.findall(r', \["([a-z_]+)"\]', text))
        found.update(re.findall(r'"(ref_mismatch|not_annotated|symbolic_allele_skipped)"', (SRC / "batch.py").read_text(encoding="utf-8")))
        self.assertGreater(len(found), 20)
        for flag in sorted(found):
            self.assertIn(f"`{flag}`", documented, f"{flag} is raised by the code but not explained")


class Links(unittest.TestCase):
    def test_relative_links_point_at_real_files_and_headings(self):
        for doc in DOCS:
            for target, anchor in _LINK_RE.findall(doc.read_text(encoding="utf-8")):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                path = (doc.parent / target).resolve()
                self.assertTrue(path.exists(), f"{doc.name} links to a missing file: {target}")
                if anchor and path.suffix == ".md":
                    self.assertIn(anchor[1:], anchors_of(path), f"{doc.name} links to a missing heading: {target}{anchor}")

    def test_links_within_a_page_point_at_real_headings(self):
        for doc in DOCS:
            text = doc.read_text(encoding="utf-8")
            for anchor in re.findall(r"\]\(#([^)\s]+)\)", text):
                self.assertIn(anchor, anchors_of(doc), f"{doc.name} links to a missing heading: #{anchor}")


class Commands(unittest.TestCase):
    def test_commands_in_the_docs_are_really_accepted(self):
        """Every line starting with varsight must parse, so the docs can't drift from the code."""
        from varsight.commandline import parse
        from varsight.errors import InputError

        checked = 0
        for doc in DOCS:
            for line in doc.read_text(encoding="utf-8").splitlines():
                command = line.strip()
                if not re.match(r"^varsight [a-z]+", command):
                    continue
                command = command.split("  #")[0].strip()
                command = re.sub(r"\s{2,}\S.*$", "", command)
                with self.subTest(doc=doc.name, line=line):
                    try:
                        parse(shlex.split(command)[1:])
                    except InputError as exc:
                        self.fail(f"{doc.name}: {command}: {exc}")
                    checked += 1
        self.assertGreater(checked, 20)

    def test_every_warning_has_plain_help(self):
        from varsight.flags import FLAG_HELP

        documented = (ROOT / "docs" / "input-output.md").read_text(encoding="utf-8")
        for flag, (meaning, action) in FLAG_HELP.items():
            self.assertTrue(meaning and action, flag)
            self.assertIn(f"`{flag}`", documented, flag)

    def test_the_readme_names_every_requested_database(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        for name in ("ClinVar", "NCBI", "HGNC", "Ensembl", "UniProt", "Pfam", "ProteinPaint", "Franklin", "COSMIC", "VarSome", "OncoKB"):
            self.assertIn(name, readme)

    def test_demo_inputs_are_all_understood(self):
        from varsight.inputs import parse_input

        for line in (ROOT / "demo" / "variants.txt").read_text(encoding="utf-8").splitlines():
            if line.strip() and not line.startswith("#"):
                parse_input(line)


class Version(unittest.TestCase):
    def test_version_is_the_same_everywhere(self):
        from varsight import __version__

        self.assertIn(f"version: {__version__}", (ROOT / "CITATION.cff").read_text(encoding="utf-8"))
        self.assertIn(f"[{__version__}]", (ROOT / "CHANGELOG.md").read_text(encoding="utf-8"))


class NoPatientData(unittest.TestCase):
    def test_no_variant_files_outside_the_demo(self):
        header_only = ROOT / "tests" / "fixtures" / "gnomad.header.vcf"
        self.assertTrue(all(line.startswith("#") for line in header_only.read_text(encoding="utf-8").splitlines()))
        for pattern in ("*.vcf", "*.vcf.gz", "*.vcf.bgz"):
            for path in ROOT.rglob(pattern):
                if ".git" in path.parts or path == header_only:
                    continue
                self.assertEqual(path.relative_to(ROOT).parts[0], "demo", f"unexpected variant file: {path}")

    def test_the_demo_vcf_says_it_is_invented(self):
        text = (ROOT / "demo" / "demo.hg19.vcf").read_text(encoding="utf-8")
        self.assertIn("Not from a patient", text)


if __name__ == "__main__":
    unittest.main()
