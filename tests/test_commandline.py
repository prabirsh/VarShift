"""The command line, the setup questions, the update notice, and the priority sorting."""

import contextlib
import errno
import io
import json
import os
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest import mock

from varshift import cli, priority, updates
from varshift.commandline import COMMANDS, parse, render_help
from varshift.config import is_windows_mount, load_settings, resolve_data_dir, save_setting
from varshift.errors import InputError


class Parsing(unittest.TestCase):
    def test_options_are_words(self):
        values = parse(["annotate", "BRAF V600E", "build=hg19", "out=x.xlsx", "replace"])
        self.assertEqual(values.command, "annotate")
        self.assertEqual(values.positional, "BRAF V600E")
        self.assertEqual(values.get("build"), "hg19")
        self.assertEqual(values.get("out"), "x.xlsx")
        self.assertTrue(values.has("replace"))

    def test_leading_dashes_are_forgiven(self):
        double = "-" * 2
        values = parse(["annotate", "x", f"{double}build", "hg19", f"{double}out=a.tsv", f"{double}replace"])
        self.assertEqual((values.get("build"), values.get("out"), values.has("replace")), ("hg19", "a.tsv", True))
        self.assertEqual(parse([f"{double}version"]).command, "version")
        self.assertEqual(parse([f"{double}help"]).command, "help")

    def test_a_variant_with_equals_or_a_leading_dash_is_still_a_variant(self):
        self.assertEqual(parse(["annotate", "BRAF p.Val600="]).positional, "BRAF p.Val600=")
        self.assertEqual(parse(["annotate", "NM_004333.6:c.-5A>T"]).positional, "NM_004333.6:c.-5A>T")

    def test_a_typo_gets_a_suggestion(self):
        with self.assertRaises(InputError) as caught:
            parse(["annotate", "x", "buid=hg19"])
        self.assertIn("build", str(caught.exception))
        with self.assertRaises(InputError) as caught:
            parse(["anotate", "x"])
        self.assertIn("annotate", str(caught.exception))

    def test_a_missing_input_says_what_is_needed(self):
        with self.assertRaises(InputError) as caught:
            parse(["annotate"])
        self.assertIn("INPUT", str(caught.exception))

    def test_a_list_option_splits_on_commas(self):
        self.assertEqual(parse(["setup", "with=dbsnp, revel"]).items("with"), ["dbsnp", "revel"])
        self.assertEqual(parse(["setup"]).items("with"), [])

    def test_no_arguments_shows_help(self):
        self.assertEqual(parse([]).command, "help")

    def test_help_exists_for_every_command(self):
        for name in COMMANDS:
            self.assertIn(COMMANDS[name].summary, render_help(name))
        self.assertIn("annotate", render_help())


class Settings(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        patcher = mock.patch.dict(os.environ, {"XDG_CONFIG_HOME": str(self.root)}, clear=False)
        patcher.start()
        self.addCleanup(patcher.stop)
        os.environ.pop("VARSHIFT_DATA", None)


class FakeTty(io.StringIO):
    def isatty(self):
        return True


class Config(Settings):
    def _config(self, *words):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = cli.main(["config", *words])
        return code, out.getvalue(), err.getvalue()

    def test_a_setting_is_remembered(self):
        save_setting("updates", "off")
        self.assertEqual(load_settings()["updates"], "off")

    def test_the_mode_is_saved_checked_and_cleared(self):
        code, out, _ = self._config("mode=tumour")
        self.assertEqual(code, 0)
        self.assertEqual(load_settings()["mode"], "tumour")
        self.assertRegex(out, r"Mode\s+tumour")
        self.assertEqual(self._config("mode=tumor")[0], 0)
        self.assertEqual(load_settings()["mode"], "tumour")
        code, _, err = self._config("mode=bogus")
        self.assertEqual(code, 1)
        for word in ("tumour", "inherited", "none"):
            self.assertIn(word, err)
        self.assertEqual(load_settings()["mode"], "tumour")
        code, out, _ = self._config("mode=none")
        self.assertEqual(code, 0)
        self.assertNotIn("mode", load_settings())
        self.assertRegex(out, r"Mode\s+none")

    def test_the_colour_is_saved_and_checked(self):
        self.assertEqual(self._config("color=no")[0], 0)
        self.assertEqual(load_settings()["color"], "no")
        code, _, err = self._config("color=bogus")
        self.assertEqual(code, 1)
        self.assertIn("auto, yes or no", err)
        self.assertEqual(load_settings()["color"], "no")
        code, out, _ = self._config()
        self.assertRegex(out, r"Colour\s+no")

    def test_the_saved_colour_is_used_and_color_on_a_command_wins(self):
        environ = {key: value for key, value in os.environ.items() if key != "NO_COLOR"}
        with mock.patch.dict(os.environ, environ, clear=True), mock.patch("sys.stdin", FakeTty()), mock.patch("sys.stdout", FakeTty()):
            save_setting("color", "no")
            self.assertFalse(cli._terminal(parse(["annotate", "BRAF V600E"])).color)
            save_setting("color", "yes")
            self.assertTrue(cli._terminal(parse(["annotate", "BRAF V600E"])).color)
            self.assertFalse(cli._terminal(parse(["annotate", "BRAF V600E", "color=no"])).color)

    def test_a_bad_saved_value_says_how_to_put_it_right(self):
        save_setting("mode", "weird")
        with self.assertRaises(InputError) as caught:
            cli._mode(parse(["annotate", "BRAF V600E"]))
        self.assertIn("varshift config mode=none", str(caught.exception))
        self.assertEqual(cli._mode(parse(["annotate", "BRAF V600E", "mode=inherited"])), "inherited")
        save_setting("color", "purple")
        with self.assertRaises(InputError) as caught:
            cli._terminal(parse(["annotate", "BRAF V600E"]))
        self.assertIn("varshift config color=auto", str(caught.exception))

    def test_data_folder_order(self):
        save_setting("data_dir", str(self.root / "saved"))
        self.assertEqual(resolve_data_dir(), self.root / "saved")
        with mock.patch.dict(os.environ, {"VARSHIFT_DATA": str(self.root / "env")}):
            self.assertEqual(resolve_data_dir(), (self.root / "env").resolve())
            self.assertEqual(resolve_data_dir(str(self.root / "given")), (self.root / "given").resolve())

    def test_setup_refuses_a_windows_drive_from_wsl_unless_told(self):
        windows_drive = self.root / "drive"
        with mock.patch.object(cli, "is_windows_mount", return_value=True), mock.patch("varshift.pipeline.manager.DataManager") as manager:
            self.assertEqual(cli.main(["setup", "build=hg19", f"data={windows_drive}", "yes"]), 1)
            manager.assert_not_called()
            with mock.patch.object(cli, "_selftest", return_value=0):
                self.assertEqual(cli.main(["setup", "build=hg19", f"data={windows_drive}", "yes", "windowsdrive"]), 0)
            manager.assert_called_once()

    def test_update_refuses_a_windows_drive_from_wsl_but_still_checks(self):
        windows_drive = str(self.root / "drive")
        with mock.patch.object(cli, "is_windows_mount", return_value=True), mock.patch("varshift.pipeline.manager.DataManager") as manager:
            manager.return_value.update.return_value = []
            with mock.patch("sys.stderr", new=io.StringIO()) as stderr, contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["update", f"data={windows_drive}"]), 1)
            self.assertIn("update check works here", stderr.getvalue())
            manager.assert_not_called()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["update", "check", f"data={windows_drive}"]), 0)
                self.assertEqual(cli.main(["update", f"data={windows_drive}", "windowsdrive"]), 0)
            self.assertEqual(manager.call_count, 2)

    def test_add_refuses_a_windows_drive_from_wsl_unless_told(self):
        windows_drive = str(self.root / "drive")
        with mock.patch.object(cli, "is_windows_mount", return_value=True), mock.patch("varshift.pipeline.manager.DataManager") as manager:
            manager.return_value.add.return_value = []
            with mock.patch("sys.stderr", new=io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["add", "hotspots", f"data={windows_drive}", "yes"]), 1)
            manager.assert_not_called()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["add", "hotspots", f"data={windows_drive}", "yes", "windowsdrive"]), 0)
            manager.assert_called_once()

    def test_running_out_of_memory_says_what_to_do(self):
        with mock.patch.object(cli, "HANDLERS", {"version": mock.Mock(side_effect=OSError(errno.ENOMEM, "Cannot allocate memory"))}):
            with mock.patch("sys.stderr", new=io.StringIO()) as stderr:
                self.assertEqual(cli.main(["version"]), 1)
        self.assertIn("Cannot allocate memory", stderr.getvalue())
        self.assertIn("run the same command again", stderr.getvalue())

    @unittest.skipIf(os.name == "nt", "a Windows folder is correct on Windows")
    def test_a_windows_folder_is_refused_on_linux(self):
        """Bash removes the backslash, so I:\\data used to become a new folder called I:data right here."""
        for text, hint in (("I:VIP-Annotations", "/mnt/i/VIP-Annotations"), ("C:\\data\\run", "/mnt/c/data/run"), ("d:/x", "/mnt/d/x")):
            with self.assertRaises(InputError) as caught:
                resolve_data_dir(text)
            self.assertIn(hint, str(caught.exception))
        with mock.patch.dict(os.environ, {"VARSHIFT_DATA": "I:\\data"}), self.assertRaises(InputError):
            resolve_data_dir()
        self.assertEqual(resolve_data_dir("/mnt/i/data"), Path("/mnt/i/data").resolve())

    @unittest.skipUnless(os.name == "nt", "only Windows has drive letters")
    def test_a_windows_folder_is_accepted_on_windows(self):
        self.assertEqual(resolve_data_dir("I:\\data"), Path("I:\\data").resolve())

    def test_a_damaged_settings_file_is_ignored(self):
        save_setting("updates", "off")
        (self.root / "varshift" / "config.json").write_text("{not json", encoding="utf-8")
        self.assertEqual(load_settings(), {})

    @unittest.skipIf(os.name == "nt", "WSL paths")
    def test_windows_drive_seen_from_wsl(self):
        self.assertTrue(is_windows_mount(Path("/mnt/c/data")))
        self.assertFalse(is_windows_mount(Path("/home/me/data")))
        self.assertFalse(is_windows_mount(Path("/mnt/shared/data")))


class FakeManager:
    """Stands in for the real downloader, and records what setup asked it to do."""

    calls: list = []

    def __init__(self, data_dir):
        type(self).calls.append(("data", data_dir))

    def setup(self, builds, options, force):
        type(self).calls.append(("setup", builds, options, force))
        return []

    def add(self, name, options):
        type(self).calls.append(("add", name, options))
        return []


class Setup(Settings):
    def setUp(self):
        super().setUp()
        # These test the questions, wherever the temporary folder happens to be. The refusal has its own test.
        patcher = mock.patch.object(cli, "is_windows_mount", return_value=False)
        patcher.start()
        self.addCleanup(patcher.stop)

    def _run(self, *arguments, answers=None, tty=True):
        FakeManager.calls = []
        output = io.StringIO()
        stdin = mock.Mock()
        stdin.isatty.return_value = tty
        with (
            mock.patch("varshift.pipeline.manager.DataManager", FakeManager),
            mock.patch.object(cli, "_selftest", return_value=0),
            mock.patch.object(cli.sys, "stdin", stdin),
            mock.patch("builtins.input", side_effect=list(answers or [])),
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(output),
        ):
            code = cli.main(["setup", *arguments])
        return code, output.getvalue(), FakeManager.calls

    def _builds(self, calls):
        return next(call[1] for call in calls if call[0] == "setup")

    def test_asks_which_build_and_installs_only_that(self):
        for answer, expected in (("1", ["hg19"]), ("2", ["hg38"]), ("3", ["hg19", "hg38"]), ("hg38", ["hg38"]), ("", ["hg19", "hg38"])):
            with self.subTest(answer=answer):
                code, text, calls = self._run(answers=[answer, str(self.root / "data"), "", "yes"])
                self.assertEqual(code, 0)
                self.assertIn("Which genome build", text)
                self.assertEqual(self._builds(calls), expected)

    def test_a_wrong_answer_asks_again(self):
        code, text, calls = self._run(answers=["7", "1", str(self.root / "data"), "", "yes"])
        self.assertEqual(code, 0)
        self.assertIn("Please type 1, 2 or 3", text)
        self.assertEqual(self._builds(calls), ["hg19"])

    def test_saying_no_changes_nothing(self):
        code, text, calls = self._run(answers=["1", str(self.root / "data"), "", "no"])
        self.assertEqual(code, 1)
        self.assertIn("Nothing was changed", text)
        self.assertFalse([call for call in calls if call[0] == "setup"])
        self.assertNotIn("data_dir", load_settings())

    def test_shows_the_download_size_before_starting(self):
        _, text, _ = self._run(answers=["2", str(self.root / "data"), "", "yes"])
        self.assertIn("hg38", text)
        self.assertIn("Download", text)
        self.assertIn("GB", text)

    def test_the_chosen_folder_is_remembered(self):
        self._run(answers=["1", str(self.root / "mine"), "", "yes"])
        self.assertEqual(Path(load_settings()["data_dir"]), (self.root / "mine").resolve())

    def test_advanced_choices_are_installed_after_the_core(self):
        # The first menu number is dbSNP. Thereafter the core runs first, then each choice.
        _, _, calls = self._run(answers=["1", str(self.root / "data"), "1", "yes"])
        kinds = [call[0] for call in calls]
        self.assertLess(kinds.index("setup"), kinds.index("add"))
        self.assertEqual(next(call[1] for call in calls if call[0] == "add"), "dbsnp")

    def test_no_questions_when_told_not_to_ask(self):
        code, text, calls = self._run("build=hg38", f"data={self.root / 'data'}", "yes")
        self.assertEqual(code, 0)
        self.assertNotIn("Which genome build", text)
        self.assertEqual(self._builds(calls), ["hg38"])

    def test_no_terminal_means_no_questions(self):
        code, _, calls = self._run(f"data={self.root / 'data'}", tty=False)
        self.assertEqual(code, 0)
        self.assertEqual(self._builds(calls), ["hg19", "hg38"])

    def test_a_bad_build_is_refused(self):
        code, text, _ = self._run("build=hg18", f"data={self.root / 'data'}", "yes")
        self.assertEqual(code, 1)
        self.assertIn("hg19, hg38 or both", text)

    def test_a_licensed_database_cannot_be_selected_without_its_file(self):
        code, text, calls = self._run("build=hg19", f"data={self.root / 'data'}", "with=cosmic", "yes")
        self.assertEqual(code, 1)
        self.assertIn("file=PATH", text)
        self.assertFalse([call for call in calls if call[0] == "setup"])

    def test_terms_must_be_accepted(self):
        code, text, _ = self._run("build=hg19", f"data={self.root / 'data'}", "with=revel", "yes")
        self.assertEqual(code, 1)
        self.assertIn("accept=revel", text)
        code, _, calls = self._run("build=hg19", f"data={self.root / 'data'}", "with=revel", "accept=revel", "yes")
        self.assertEqual(code, 0)
        self.assertEqual(next(call[1] for call in calls if call[0] == "add"), "revel")


class UpdateNotice(Settings):
    def setUp(self):
        super().setUp()
        self.data = self.root / "data"
        self.data.mkdir()
        self.manifest = {"updated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}

    def _status(self, items, hours=1):
        stamp = (datetime.now(timezone.utc) - timedelta(hours=hours)).strftime("%Y-%m-%dT%H:%M:%SZ")
        (self.data / updates.STATUS_FILE).write_text(json.dumps({"checked_at": stamp, "items": items}), encoding="utf-8")

    def test_nothing_new_means_no_message(self):
        self._status([])
        self.assertEqual(updates.notice_lines(self.data, self.manifest), [])

    def test_a_new_clinvar_is_announced_with_the_command(self):
        self._status([{"artifact": "clinvar.hg19", "title": "ClinVar", "policy": "auto"}, {"artifact": "clinvar.hg38", "title": "ClinVar", "policy": "auto"}])
        lines = updates.notice_lines(self.data, self.manifest)
        self.assertEqual(lines[0].count("ClinVar"), 1)
        self.assertIn("varshift update", " ".join(lines))
        self.assertIn("installed versions", lines[0])

    def test_a_change_that_can_rename_things_asks_for_approval(self):
        self._status([{"artifact": "transcripts.hg19", "title": "RefSeq transcripts", "policy": "approval"}])
        self.assertIn("approve=transcripts", " ".join(updates.notice_lines(self.data, self.manifest)))

    def test_very_old_data_is_mentioned_when_no_check_has_worked(self):
        old = {"updated": (datetime.now(timezone.utc) - timedelta(days=45)).strftime("%Y-%m-%dT%H:%M:%SZ")}
        self.assertIn("45 days", " ".join(updates.notice_lines(self.data, old)))

    def test_it_can_be_switched_off(self):
        self._status([{"artifact": "clinvar.hg19", "title": "ClinVar", "policy": "auto"}])
        save_setting("updates", "off")
        self.assertEqual(updates.notice_lines(self.data, self.manifest), [])
        with mock.patch.object(updates.subprocess, "Popen") as started:
            self.assertFalse(updates.check_in_background(self.data, self.manifest))
        started.assert_not_called()

    def test_demo_data_never_checks(self):
        with mock.patch.object(updates.subprocess, "Popen") as started:
            self.assertFalse(updates.check_in_background(self.data, {"demo": True}))
        started.assert_not_called()

    def test_a_check_runs_only_when_the_last_is_over_a_day_old(self):
        self._status([], hours=2)
        with mock.patch.object(updates.subprocess, "Popen") as started:
            self.assertFalse(updates.check_in_background(self.data, self.manifest))
        started.assert_not_called()
        self._status([], hours=25)
        with mock.patch.object(updates.subprocess, "Popen") as started:
            self.assertTrue(updates.check_in_background(self.data, self.manifest))
        arguments = started.call_args.args[0]
        self.assertEqual(arguments[1:3], ["-m", "varshift.pipeline.updatecheck"])
        self.assertEqual(arguments[3], str(self.data))

    def test_a_loop_of_commands_starts_only_one_check(self):
        with mock.patch.object(updates.subprocess, "Popen") as started:
            results = [updates.check_in_background(self.data, self.manifest) for _ in range(5)]
        self.assertEqual(results, [True, False, False, False, False])
        self.assertEqual(started.call_count, 1)

    def test_a_check_is_tried_again_after_an_hour(self):
        marker = self.data / updates.STARTED_FILE
        with mock.patch.object(updates.subprocess, "Popen") as started:
            updates.check_in_background(self.data, self.manifest)
            past = marker.stat().st_mtime - 2 * 3600
            os.utime(marker, (past, past))
            self.assertTrue(updates.check_in_background(self.data, self.manifest))
        self.assertEqual(started.call_count, 2)

    def test_a_folder_that_cannot_be_written_means_no_check(self):
        with mock.patch.object(updates.subprocess, "Popen") as started:
            self.assertFalse(updates.check_in_background(self.root / "absent", self.manifest))
        started.assert_not_called()

    def test_a_failure_to_start_never_stops_annotation(self):
        with mock.patch.object(updates.subprocess, "Popen", side_effect=OSError):
            self.assertFalse(updates.check_in_background(self.data, self.manifest))

    def test_applied_updates_leave_the_notice(self):
        self._status([{"artifact": "clinvar.hg19", "title": "ClinVar", "policy": "auto"}, {"artifact": "civic", "title": "CIViC", "policy": "auto"}])
        updates.clear_applied(self.data, ["clinvar.hg19"])
        self.assertEqual([item["artifact"] for item in updates.read_status(self.data)["items"]], ["civic"])

    def test_a_damaged_status_file_is_ignored(self):
        (self.data / updates.STATUS_FILE).write_text("{not json", encoding="utf-8")
        self.assertEqual(updates.notice_lines(self.data, self.manifest), [])


class Priority(unittest.TestCase):
    def _level(self, clinvar=None, terms=(), advanced=None, civic=None):
        record = {"evidence": {"clinvar": clinvar or {}, "advanced": advanced or {}, "civic": civic or []}, "effect": {"all": list(terms)}}
        return priority.assign(record)

    def test_pathogenic_comes_first(self):
        result = self._level({"classification": "Pathogenic"})
        self.assertEqual(result["label"], "High")
        self.assertIn("Pathogenic", result["reasons"][0])

    def test_conflicting_is_not_treated_as_pathogenic(self):
        self.assertEqual(self._level({"classification": "Conflicting classifications of pathogenicity"})["label"], "Medium")

    def test_benign_is_low_even_when_it_is_missense(self):
        self.assertEqual(self._level({"classification": "Benign"}, ["missense_variant"])["label"], "Low")
        self.assertEqual(self._level({"classification": "Likely benign"}, ["missense_variant"])["label"], "Low")

    def test_an_expert_panel_benign_call_lowers_it(self):
        result = self._level(terms=["missense_variant"], advanced={"clingen_variants": {"assertion": "Benign"}})
        self.assertEqual(result["label"], "Low")
        self.assertIn("ClinGen expert panel: Benign", result["reasons"])

    def test_truncating_with_no_record_is_medium(self):
        self.assertEqual(self._level(terms=["stop_gained"])["label"], "Medium")

    def test_uncertain_and_plain_missense_need_review(self):
        self.assertEqual(self._level({"classification": "Uncertain significance"})["label"], "Review")
        self.assertEqual(self._level(terms=["missense_variant"])["label"], "Review")

    def test_synonymous_with_nothing_known_is_low(self):
        self.assertEqual(self._level(terms=["synonymous_variant"])["label"], "Low")

    def test_cancer_evidence_raises_it(self):
        self.assertEqual(self._level(civic=[{"evidence_items": 12}])["label"], "High")
        self.assertEqual(self._level(advanced={"hotspots": {"label": "V600 in 800 tumours"}})["label"], "High")

    def test_only_somatic_tiers_one_and_two_raise_it(self):
        self.assertEqual(self._level({"somatic_impact": "Tier I - Strong"})["label"], "High")
        self.assertEqual(self._level({"somatic_impact": "Tier II - Potential"})["label"], "High")
        self.assertNotEqual(self._level({"somatic_impact": "Tier III - Unknown"})["label"], "High")
        self.assertNotEqual(self._level({"somatic_impact": "Tier IV - Benign/Likely benign"})["label"], "High")

    def test_score_thresholds_are_the_calibrated_ones(self):
        self.assertEqual(self._level(terms=["missense_variant"], advanced={"revel": {"score": 0.644}})["label"], "Medium")
        self.assertEqual(self._level(terms=["missense_variant"], advanced={"revel": {"score": 0.643}})["label"], "Review")
        self.assertEqual(self._level(advanced={"spliceai": {"max_delta": 0.2}})["label"], "Medium")

    def test_common_variants_are_low(self):
        result = self._level(terms=["missense_variant"], advanced={"gnomad": {"af": 0.2}})
        self.assertEqual(result["label"], "Low")
        self.assertIn("gnomAD", result["reasons"][0])

    def test_a_common_variant_reported_pathogenic_stays_high(self):
        result = self._level({"classification": "Pathogenic"}, advanced={"gnomad": {"af": 0.2}})
        self.assertEqual(result["label"], "High")

    def test_a_reason_is_always_given(self):
        self.assertTrue(self._level()["reasons"])


if __name__ == "__main__":
    unittest.main()
