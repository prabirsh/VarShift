"""Security properties enforced by reading the source, so a later change cannot quietly break them."""

import ast
import unittest
from pathlib import Path


SRC = Path(__file__).resolve().parents[1] / "src" / "varshift"
NETWORK_MODULES = {"socket", "ssl", "http", "urllib.request", "urllib.error", "ftplib", "smtplib", "requests", "httpx"}
FORBIDDEN_MODULES = {"pickle", "marshal", "shelve", "telnetlib", "dill"}
FORBIDDEN_CALLS = {"eval", "exec", "compile", "__import__"}


def imports_of(tree: ast.AST) -> set[str]:
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            found.add(node.module)
            found.update(f"{node.module}.{alias.name}" for alias in node.names)
    return found


def sources():
    for path in SRC.rglob("*.py"):
        yield path, ast.parse(path.read_text(encoding="utf-8"))


class SourceRules(unittest.TestCase):
    def test_annotation_code_cannot_reach_the_network(self):
        for path, tree in sources():
            if "pipeline" in path.parts:
                continue
            bad = imports_of(tree) & NETWORK_MODULES
            self.assertFalse(bad, f"{path.name} imports {bad}")

    def test_only_fetch_module_talks_to_the_network(self):
        for path, tree in sources():
            if path.name == "fetch.py":
                continue
            bad = imports_of(tree) & NETWORK_MODULES
            self.assertFalse(bad, f"{path.name} imports {bad}")

    def test_no_dangerous_modules(self):
        for path, tree in sources():
            bad = {name for name in imports_of(tree) if name.split(".")[0] in FORBIDDEN_MODULES}
            self.assertFalse(bad, f"{path.name} imports {bad}")

    def test_only_the_update_check_and_the_windows_root_request_start_a_process(self):
        for path, tree in sources():
            if path.name not in ("updates.py", "fetch.py"):
                self.assertNotIn("subprocess", {name.split(".")[0] for name in imports_of(tree)}, path.name)

    def test_the_windows_root_request_has_a_fixed_command_line(self):
        """curl.exe from System32 is given only the host of a download that has already passed the allow-list."""
        text = (SRC / "pipeline" / "fetch.py").read_text(encoding="utf-8")
        self.assertEqual(text.count("subprocess.run("), 1)
        self.assertIn('subprocess.run([str(curl), "-sI", "-m", "30", f"https://{host}/"]', text)
        self.assertIn('"System32" / "curl.exe"', text)

    def test_no_process_ever_goes_through_a_shell(self):
        for path, tree in sources():
            for node in ast.walk(tree):
                if isinstance(node, ast.keyword) and node.arg == "shell":
                    self.fail(f"{path.name}:{node.value.lineno} passes shell=")

    def test_the_update_check_cannot_carry_variants(self):
        """It is started with fixed arguments: the interpreter, the module and the data folder."""
        text = (SRC / "updates.py").read_text(encoding="utf-8")
        self.assertIn('[sys.executable, "-m", "varshift.pipeline.updatecheck", str(data_dir)]', text)
        check = (SRC / "pipeline" / "updatecheck.py").read_text(encoding="utf-8")
        for word in ("variant", "vcf", "annotate"):
            self.assertNotIn(word, check.lower())

    def test_no_code_execution_calls(self):
        for path, tree in sources():
            for node in ast.walk(tree):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                    self.assertNotIn(node.func.id, FORBIDDEN_CALLS, f"{path.name}:{node.lineno}")

    def test_no_shell_helpers(self):
        for path, tree in sources():
            for node in ast.walk(tree):
                if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "os":
                    self.assertNotIn(node.attr, {"system", "popen", "spawnl", "startfile"}, f"{path.name}:{node.lineno}")

    def test_no_third_party_runtime_imports(self):
        allowed = set(__import__("sys").stdlib_module_names) | {"varshift"}
        for path, tree in sources():
            for name in imports_of(tree):
                root = name.split(".")[0]
                if root:
                    self.assertIn(root, allowed, f"{path.name} imports {name}")

    def test_sql_is_never_built_from_formatted_user_text(self):
        """f-strings near execute are limited to fixed table names and placeholder lists."""
        allowed = {
            "store.py": {"marks", "table"},
            "builders.py": {"folder"},
        }
        for path, tree in sources():
            for node in ast.walk(tree):
                if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)):
                    continue
                if node.func.attr not in ("execute", "executemany", "executescript") or not node.args:
                    continue
                argument = node.args[0]
                if isinstance(argument, ast.JoinedStr):
                    names = {n.id for v in argument.values if isinstance(v, ast.FormattedValue) for n in ast.walk(v.value) if isinstance(n, ast.Name)}
                    self.assertLessEqual(names, allowed.get(path.name, set()), f"{path.name}:{node.lineno}")

    def test_connections_in_the_reader_are_read_only(self):
        text = (SRC / "store.py").read_text(encoding="utf-8")
        self.assertIn("mode=ro", text)
        self.assertIn("PRAGMA query_only=ON", text)
        self.assertIn("PRAGMA trusted_schema=OFF", text)

    def test_only_ucsc_own_amazon_bucket_is_allowed(self):
        """UCSC redirects its downloads there. Allowing all of Amazon would let anyone's bucket through."""
        import json

        hosts = json.loads((SRC / "data" / "sources.json").read_text(encoding="utf-8"))["allowed_hosts"]
        self.assertIn("genome-browser.s3.amazonaws.com", hosts)
        self.assertEqual([host for host in hosts if "amazonaws" in host], ["genome-browser.s3.amazonaws.com"])

    def test_every_download_host_is_on_the_allow_list_and_https(self):
        import json
        from urllib.parse import urlsplit

        catalog = json.loads((SRC / "data" / "sources.json").read_text(encoding="utf-8"))
        hosts = set(catalog["allowed_hosts"])
        for name, source in catalog["sources"].items():
            for key in ("url", "list_url", "md5_url", "md5_list", "url_template"):
                if key in source:
                    parts = urlsplit(source[key])
                    self.assertEqual(parts.scheme, "https", name)
                    self.assertIn(parts.hostname, hosts, name)
            self.assertGreater(source["max_mb"], 0, name)
            self.assertTrue(source.get("license"), f"{name} has no licence note")

    def test_every_artifact_source_exists(self):
        import json

        catalog = json.loads((SRC / "data" / "sources.json").read_text(encoding="utf-8"))
        for name, artifact in catalog["artifacts"].items():
            for source in artifact["sources"]:
                self.assertIn(source, catalog["sources"], name)
            self.assertIn(artifact["policy"], ("once", "auto", "approval", "own_file"), name)


if __name__ == "__main__":
    unittest.main()
