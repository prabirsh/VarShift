"""What kind of terminal VarShift is talking to, and how to ask the person at it.

VarShift only waits for an answer when a person is there: standard input and standard output must
both be terminals. Anything piped, scheduled or redirected gets a QuietTerminal, which never reads.
"""

from __future__ import annotations

import os
import shutil
import sys

SETTINGS = ("auto", "yes", "no")
MIN_WIDTH = 60
MAX_WIDTH = 120
_ATTEMPTS = 3
_RESET = "\x1b[0m"
_STYLES = {
    "high": "1;31",
    "medium": "1;33",
    "review": "1;34",
    "low": "2",
    "label": "2",
    "strong": "1",
    "error": "31",
}


class Terminal:
    def __init__(self, interactive: bool, color: bool, width: int, stdin=None, stdout=None):
        self.interactive = interactive
        self.color = color
        self.width = width
        self.stdin = stdin if stdin is not None else sys.stdin
        self.stdout = stdout if stdout is not None else sys.stdout

    @classmethod
    def detect(cls, color_setting: str | None = "auto", stdin=None, stdout=None, environ=None) -> Terminal:
        """Look at the real streams. color_setting is auto, yes or no; NO_COLOR overrides all three."""
        setting = (color_setting or "auto").lower()
        if setting not in SETTINGS:
            raise ValueError(f"colour must be auto, yes or no, not {color_setting!r}")
        stdin = stdin if stdin is not None else sys.stdin
        stdout = stdout if stdout is not None else sys.stdout
        environ = os.environ if environ is None else environ
        interactive = _isatty(stdin) and _isatty(stdout)
        if environ.get("NO_COLOR"):
            color = False
        elif setting == "yes":
            color = True
            _supports_color(stdout, environ)  # still switch a Windows console into colour mode if it can be
        elif setting == "no":
            color = False
        else:
            color = interactive and _supports_color(stdout, environ)
        columns = shutil.get_terminal_size((100, 24)).columns
        return cls(interactive, color, min(max(columns, MIN_WIDTH), MAX_WIDTH), stdin, stdout)

    # showing
    def style(self, text: str, role: str) -> str:
        code = _STYLES[role]
        return f"\x1b[{code}m{text}{_RESET}" if self.color else text

    def write(self, text: str) -> None:
        self.stdout.write(text)
        self.stdout.flush()

    # asking
    def ask(self, prompt: str) -> str | None:
        """One line from the person, without its line ending. None at the end of input."""
        if self.stdin is sys.stdin and self.stdout is sys.stdout:
            try:
                import readline  # noqa: F401  # gives arrow keys and history on Linux; Windows has its own
            except ImportError:
                pass
            try:
                return input(prompt)
            except EOFError:
                return None
        self.write(prompt)
        line = self.stdin.readline()
        if not line:
            return None
        return line.rstrip("\r\n")

    def choose(self, question: str, options: list[str], allow_all: bool = False) -> list[int] | None:
        """Numbered options. Returns the chosen positions from 0, or None when cancelled."""
        count = len(options)
        lines = [question] + [f"{number:>3}  {option}" for number, option in enumerate(options, 1)]
        self.write("\n".join(lines) + "\n")
        if allow_all:
            prompt = f"Choose 1 to {count}, a for all, or Enter for the first: "
        else:
            prompt = f"Choose 1 to {count}, or Enter to cancel: "
        for _ in range(_ATTEMPTS):
            answer = self.ask(self.style(prompt, "strong"))
            if answer is None:
                return None
            answer = answer.strip().lower()
            if not answer:
                return [0] if allow_all else None
            if allow_all and answer in ("a", "all"):
                return list(range(count))
            if answer.isdigit() and 1 <= int(answer) <= count:
                return [int(answer) - 1]
            self.write(f"Please type a number from 1 to {count}.\n")
        return None

    def confirm(self, question: str, default: bool = False) -> bool:
        hint = "[Y/n]" if default else "[y/N]"
        for _ in range(_ATTEMPTS):
            answer = self.ask(self.style(f"{question} {hint} ", "strong"))
            if answer is None:
                return default
            answer = answer.strip().lower()
            if not answer:
                return default
            if answer in ("y", "yes"):
                return True
            if answer in ("n", "no"):
                return False
            self.write("Please answer y or n.\n")
        return default

    def input_waiting(self) -> bool:
        """True when more pasted text is already waiting, so pasted lines can be answered together."""
        if not self.interactive:
            return False
        try:
            if os.name == "nt":
                import msvcrt

                return bool(msvcrt.kbhit())
            import select

            ready, _, _ = select.select([self.stdin], [], [], 0.05)
            return bool(ready)
        except (OSError, ValueError, ImportError):
            return False


class QuietTerminal(Terminal):
    """For scripts and files: never reads, never colours."""

    def __init__(self, stdout=None, width: int = 100):
        super().__init__(interactive=False, color=False, width=width, stdin=None, stdout=stdout)

    def ask(self, prompt: str) -> str | None:
        return None

    def choose(self, question: str, options: list[str], allow_all: bool = False) -> list[int] | None:
        return None

    def confirm(self, question: str, default: bool = False) -> bool:
        return default

    def input_waiting(self) -> bool:
        return False


def _isatty(stream) -> bool:
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError, OSError):
        return False


def _supports_color(stream, environ) -> bool:
    if os.name == "nt":
        return _enable_windows_vt(stream)
    return environ.get("TERM", "") != "dumb"


def _enable_windows_vt(stream) -> bool:
    """Switch a Windows console into the mode that understands colour codes. False if it won't."""
    try:
        import ctypes
        import msvcrt

        kernel32 = ctypes.windll.kernel32
        handle = msvcrt.get_osfhandle(stream.fileno())
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        if mode.value & 0x0004:
            return True
        return bool(kernel32.SetConsoleMode(handle, mode.value | 0x0004))
    except (AttributeError, ImportError, OSError, ValueError):
        return False
