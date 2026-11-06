"""Errors that are shown to the user. Anything else is a bug and keeps its traceback."""

from __future__ import annotations


class VarShiftError(Exception):
    """Base class for expected failures."""


class InputError(VarShiftError):
    """The input could not be understood, or something needed is missing from it."""


class ReferenceMismatch(VarShiftError):
    """The stated reference allele is not what the genome or transcript has.

    position and transcript say where, when the mismatch is on a transcript: a residue number for a
    protein change, or the c. or n. position as text for a DNA change. Suggestions read them.
    """

    def __init__(self, message: str, expected: str = "", found: str = "", position: int | str | None = None, transcript: str = ""):
        super().__init__(message)
        self.expected = expected
        self.found = found
        self.position = position
        self.transcript = transcript


class ResolveError(VarShiftError):
    """A c. or p. description could not be placed on the genome.

    tried says, in words, what failed on each transcript; causes keeps those errors, in the same order.
    """

    def __init__(self, message: str, tried: list[str] | None = None, causes: list[VarShiftError] | None = None):
        super().__init__(message)
        self.tried = tried or []
        self.causes = causes or []


class DataError(VarShiftError):
    """Local data is missing, damaged or the wrong version."""


class SetupError(VarShiftError):
    """Download or build of the local data failed."""
