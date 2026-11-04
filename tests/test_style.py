"""House rules for the text of the package itself. Third-party samples are left as they came."""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP_FOLDERS = {".git", "__pycache__", "fixtures", "sources", "demo_output", "build", "dist"}
TEXT_SUFFIXES = {".py", ".md", ".sh", ".ps1", ".bat", ".yml", ".yaml", ".toml", ".cff", ".json", ".html", ".svg", ".txt", ".vcf", ""}
DOUBLE_HYPHEN = "-" * 2
LONG_DASHES = (chr(0x2013), chr(0x2014))
# Google's public bucket for gnomAD really is named this way, so a download address can't avoid it.
ALLOWED = "gcp-public-data" + DOUBLE_HYPHEN + "gnomad"


def package_files():
    for path in ROOT.rglob("*"):
        parts = path.relative_to(ROOT).parts
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES or SKIP_FOLDERS & set(parts):
            continue
        # A VCF or a results folder in the working copy is usually someone's sample, not part of the package.
        if any(part.endswith("_varshift") for part in parts) or (path.suffix.lower() == ".vcf" and parts != ("demo", "demo.hg19.vcf")):
            continue
        yield path


class Text(unittest.TestCase):
    def test_no_double_hyphens_anywhere(self):
        found = []
        for path in package_files():
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if DOUBLE_HYPHEN in line.replace(ALLOWED, ""):
                    found.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()[:80]}")
        self.assertEqual(found, [])

    def test_no_long_dashes_anywhere(self):
        found = []
        for path in package_files():
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if any(dash in line for dash in LONG_DASHES):
                    found.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()[:80]}")
        self.assertEqual(found, [])

    def test_this_check_actually_looks_at_the_package(self):
        names = {path.name for path in package_files()}
        for expected in ("README.md", "cli.py", "install.sh", "install.ps1", "install.bat", "ci.yml", "report_template.html"):
            self.assertIn(expected, names)
        self.assertNotIn("gnomad.header.vcf", names)
        self.assertIn("demo.hg19.vcf", names)


if __name__ == "__main__":
    unittest.main()
