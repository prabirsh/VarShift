"""varshift alone, and varshift start: the guided start opens only when a person is at a terminal."""

import contextlib
import io
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest import mock

from tests.helpers import ScriptedTerminal
from tests.test_integration import NoTyping, WithData, isolate_settings
from varshift import cli

ROOT = Path(__file__).resolve().parents[1]


class Start(WithData):
    def setUp(self):
        isolate_settings(self)
        patcher = mock.patch.dict(os.environ, {"VARSHIFT_DATA": str(self.data)})
        patcher.start()
        self.addCleanup(patcher.stop)
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.cwd = Path(folder.name)

    def _main(self, argv, terminal=None):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.ExitStack() as stack:
            stack.enter_context(contextlib.redirect_stdout(out))
            stack.enter_context(contextlib.redirect_stderr(err))
            stack.enter_context(mock.patch("sys.stdin", NoTyping()))
            stack.enter_context(mock.patch.object(cli.Path, "cwd", return_value=self.cwd))
            if terminal is not None:
                stack.enter_context(mock.patch.object(cli, "_terminal", return_value=terminal))
            code = cli.main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_varshift_alone_without_a_terminal_prints_the_commands(self):
        code, out, _ = self._main([])
        self.assertEqual(code, 0)
        self.assertIn("Commands:", out)
        self.assertIn("start", out)

    def test_start_without_a_terminal_says_it_needs_one(self):
        code, _, err = self._main(["start"])
        self.assertEqual(code, 1)
        self.assertIn("varshift start needs a terminal", err)
        self.assertIn("varshift annotate", err)

    def test_varshift_alone_at_a_terminal_opens_the_session(self):
        terminal = ScriptedTerminal(answers=[""])
        code, out, _ = self._main([], terminal)
        self.assertEqual(code, 0)
        self.assertIn("What do you have?", terminal.output)
        self.assertNotIn("Commands:", out)

    def test_start_takes_the_mode_and_the_data_folder(self):
        terminal = ScriptedTerminal(answers=["KRAS G12D", "", "n"])
        code, _, err = self._main(["start", f"data={self.data}", "mode=tumour"], terminal)
        self.assertEqual(code, 0, err)
        self.assertIn("OncoKB", terminal.output.split("Check it")[1])

    def test_start_checks_its_options_before_opening(self):
        for word in ("build=hg18", "mode=bogus"):
            with self.subTest(word=word):
                terminal = ScriptedTerminal(answers=["BRAF V600E"])
                code, _, _ = self._main(["start", word], terminal)
                self.assertEqual(code, 1)
                self.assertEqual(terminal.position, 0, "nothing was read")

    def test_a_card_in_the_session_says_how_to_see_the_other_candidates(self):
        terminal = ScriptedTerminal(answers=["BRAF V600E", "all 1", "", "n"])
        self.assertEqual(self._main(["start"], terminal)[0], 0)
        self.assertIn("To see it: all 1", terminal.output)
        self.assertNotIn("show=all", terminal.output)
        self.assertIn("candidate 2 of 2", terminal.output)


class NoPromptThroughAPipe(WithData):
    """Run the real program with standard input a closed pipe: it must never wait, and never ask."""

    def _run(self, *args):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        environment = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1", XDG_CONFIG_HOME=folder.name)
        environment.pop("VARSHIFT_DATA", None)
        process = subprocess.Popen(
            [sys.executable, "-m", "varshift", *args],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=folder.name,
            env=environment,
        )
        process.stdin.close()
        try:
            out, err = process.communicate(timeout=60)
        except subprocess.TimeoutExpired:
            process.kill()
            self.fail(f"varshift {' '.join(args)} waited for input through a pipe")
        out, err = out.decode("utf-8", "replace"), err.decode("utf-8", "replace")
        self.assertNotIn("Choose 1 to", out + err)
        return process.returncode, out, err

    def test_varshift_alone(self):
        code, out, _ = self._run()
        self.assertEqual(code, 0)
        self.assertIn("Commands:", out)

    def test_start(self):
        code, _, err = self._run("start", f"data={self.data}")
        self.assertEqual(code, 1)
        self.assertIn("needs a terminal", err)

    def test_a_slip(self):
        code, _, err = self._run("annotate", "BRAF V600", f"data={self.data}", "build=hg19")
        self.assertEqual(code, 1)
        self.assertIn("Did you mean: BRAF V600E", err)

    def test_a_semicolon_list_with_a_slip(self):
        code, out, err = self._run("annotate", "BRAF V600E; BRAF V600", f"data={self.data}", "build=hg19")
        self.assertEqual(code, 0, err)
        self.assertIn("did you mean BRAF V600E", out)

    def test_a_choice(self):
        code, out, err = self._run("annotate", "rs113488022", f"data={self.data}", "build=hg19")
        self.assertEqual(code, 0, err)
        self.assertIn("To see all in full: show=all", out)

    def test_a_list_file(self):
        listing = Path(tempfile.mkdtemp(dir=self.data.parent)) / "slips.txt"
        listing.write_text("BRAF V600\nKRAS G12D\n", encoding="utf-8")
        code, _, err = self._run("annotate", str(listing), f"data={self.data}", "build=hg19", f"out={listing.parent / 'out.tsv'}")
        self.assertEqual(code, 0, err)
        self.assertIn("with a suggestion: 1", err)
