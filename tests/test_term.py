"""The terminal: whether it is a person at a keyboard, colour, width, and asking."""

import io
import os
import unittest
from unittest import mock

from tests.helpers import ScriptedTerminal
from varshift import term
from varshift.term import QuietTerminal, Terminal


class FakeTty(io.StringIO):
    def isatty(self):
        return True


def typed(*lines: str) -> Terminal:
    """An interactive terminal whose keyboard holds these lines."""
    return Terminal(interactive=True, color=False, width=100, stdin=io.StringIO("".join(f"{line}\n" for line in lines)), stdout=io.StringIO())


class Detect(unittest.TestCase):
    def test_pipes_are_not_a_person(self):
        found = Terminal.detect(stdin=io.StringIO(), stdout=io.StringIO(), environ={})
        self.assertFalse(found.interactive)
        self.assertFalse(found.color)

    def test_both_ends_must_be_a_terminal(self):
        self.assertFalse(Terminal.detect(stdin=FakeTty(), stdout=io.StringIO(), environ={}).interactive)
        self.assertFalse(Terminal.detect(stdin=io.StringIO(), stdout=FakeTty(), environ={}).interactive)
        self.assertTrue(Terminal.detect(stdin=FakeTty(), stdout=FakeTty(), environ={}).interactive)

    def test_colour_follows_what_the_terminal_supports(self):
        for supported in (True, False):
            with self.subTest(supported=supported), mock.patch.object(term, "_supports_color", return_value=supported):
                self.assertEqual(Terminal.detect(stdin=FakeTty(), stdout=FakeTty(), environ={}).color, supported)

    def test_no_color_wins_over_everything(self):
        with mock.patch.object(term, "_supports_color", return_value=True):
            found = Terminal.detect("yes", stdin=FakeTty(), stdout=FakeTty(), environ={"NO_COLOR": "1"})
        self.assertFalse(found.color)

    def test_an_empty_no_color_is_ignored(self):
        """no-color.org: only a value that is present and not empty turns colour off."""
        with mock.patch.object(term, "_supports_color", return_value=True):
            self.assertTrue(Terminal.detect(stdin=FakeTty(), stdout=FakeTty(), environ={"NO_COLOR": ""}).color)

    def test_the_setting_can_force_colour_off_or_on(self):
        with mock.patch.object(term, "_supports_color", return_value=True):
            self.assertFalse(Terminal.detect("no", stdin=FakeTty(), stdout=FakeTty(), environ={}).color)
        self.assertTrue(Terminal.detect("yes", stdin=io.StringIO(), stdout=io.StringIO(), environ={}).color)

    def test_an_unknown_setting_is_refused(self):
        with self.assertRaises(ValueError):
            Terminal.detect("sometimes", stdin=io.StringIO(), stdout=io.StringIO(), environ={})

    def test_width_stays_readable(self):
        for columns, expected in ((200, 120), (30, 60), (100, 100)):
            with self.subTest(columns=columns), mock.patch.object(term.shutil, "get_terminal_size", return_value=os.terminal_size((columns, 40))):
                self.assertEqual(Terminal.detect(stdin=io.StringIO(), stdout=io.StringIO(), environ={}).width, expected)

    def test_a_dumb_terminal_gets_no_colour(self):
        if os.name == "nt":
            self.skipTest("TERM is a Linux convention")
        self.assertFalse(term._supports_color(FakeTty(), {"TERM": "dumb"}))
        self.assertTrue(term._supports_color(FakeTty(), {"TERM": "xterm-256color"}))

    def test_a_windows_console_that_refuses_colour_means_none(self):
        self.assertFalse(term._enable_windows_vt(io.StringIO()))


class Style(unittest.TestCase):
    def test_plain_without_colour(self):
        self.assertEqual(Terminal(interactive=True, color=False, width=100).style("HIGH", "high"), "HIGH")

    def test_codes_with_colour(self):
        styled = Terminal(interactive=True, color=True, width=100).style("HIGH", "high")
        self.assertEqual(styled, "\x1b[1;31mHIGH\x1b[0m")

    def test_every_role_has_a_style_and_an_unknown_role_is_an_error(self):
        painted = Terminal(interactive=True, color=True, width=100)
        for role in ("high", "medium", "review", "low", "label", "strong", "error"):
            self.assertIn("\x1b[", painted.style("x", role))
        with self.assertRaises(KeyError):
            painted.style("x", "sparkly")


class Choose(unittest.TestCase):
    def test_a_number_picks_that_option(self):
        self.assertEqual(typed("2").choose("Which?", ["V600E", "V600K"]), [1])

    def test_a_for_all_when_allowed(self):
        self.assertEqual(typed("a").choose("Which?", ["x", "y", "z"], allow_all=True), [0, 1, 2])

    def test_enter_takes_the_first_when_all_are_allowed_and_cancels_otherwise(self):
        self.assertEqual(typed("").choose("Which?", ["x", "y"], allow_all=True), [0])
        self.assertIsNone(typed("").choose("Which?", ["x", "y"]))

    def test_it_asks_again_then_gives_up(self):
        self.assertEqual(typed("9", "x", "1").choose("Which?", ["only"]), [0])
        self.assertIsNone(typed("9", "x", "nope").choose("Which?", ["only"]))

    def test_end_of_input_cancels(self):
        self.assertIsNone(typed().choose("Which?", ["x"], allow_all=True))

    def test_the_options_are_shown_numbered(self):
        asking = typed("1")
        asking.choose("Val600 changes to what?", ["V600E", "V600K"])
        shown = asking.stdout.getvalue()
        self.assertIn("Val600 changes to what?", shown)
        self.assertIn("  1  V600E", shown)
        self.assertIn("  2  V600K", shown)
        self.assertIn("Choose 1 to 2, or Enter to cancel", shown)


class Confirm(unittest.TestCase):
    def test_yes_no_and_the_default(self):
        self.assertTrue(typed("y").confirm("Save first?"))
        self.assertTrue(typed("YES").confirm("Save first?"))
        self.assertFalse(typed("n").confirm("Save first?", default=True))
        self.assertFalse(typed("").confirm("Save first?"))
        self.assertTrue(typed("").confirm("Save first?", default=True))
        self.assertFalse(typed().confirm("Save first?"))

    def test_the_default_is_shown_in_capitals(self):
        asking = typed("")
        asking.confirm("Save the 3 results first?")
        self.assertIn("Save the 3 results first? [y/N]", asking.stdout.getvalue())


class Quiet(unittest.TestCase):
    def test_it_never_asks(self):
        quiet = QuietTerminal()
        self.assertFalse(quiet.interactive)
        self.assertIsNone(quiet.ask("> "))
        self.assertIsNone(quiet.choose("Which?", ["x"], allow_all=True))
        self.assertTrue(quiet.confirm("Sure?", default=True))
        self.assertFalse(quiet.input_waiting())
        self.assertEqual(quiet.style("HIGH", "high"), "HIGH")

    def test_nothing_is_waiting_when_not_interactive(self):
        self.assertFalse(Terminal(interactive=False, color=False, width=100, stdin=io.StringIO("more\n")).input_waiting())


class Scripted(unittest.TestCase):
    def test_answers_in_order_and_records_what_was_shown(self):
        script = ScriptedTerminal(["BRAF V600E", "save"])
        self.assertEqual(script.ask("> "), "BRAF V600E")
        script.write("card\n")
        self.assertEqual(script.ask("> "), "save")
        self.assertIsNone(script.ask("> "))
        self.assertEqual(script.output, "> BRAF V600E\ncard\n> save\n> ")

    def test_pasted_answers_say_more_is_waiting(self):
        script = ScriptedTerminal(["a", "b", "c"], pasted={1, 2})
        script.ask("> ")
        self.assertTrue(script.input_waiting())
        script.ask("")
        self.assertTrue(script.input_waiting())
        script.ask("")
        self.assertFalse(script.input_waiting())

    def test_an_exception_in_the_script_is_raised_at_that_point(self):
        script = ScriptedTerminal(["x", KeyboardInterrupt()])
        script.ask("> ")
        with self.assertRaises(KeyboardInterrupt):
            script.ask("> ")


if __name__ == "__main__":
    unittest.main()
