"""Errors that are shown to the user. Anything else is a bug and keeps its traceback."""

from __future__ import annotations


class VarSightError(Exception):
    """Base class for expected failures."""


class InputError(VarSightError):
    """The input could not be understood, or something needed is missing from it."""


class ReferenceMismatch(VarSightError):
    """The stated reference allele is not what the genome or transcript has."""

    def __init__(self, message: str, expected: str = "", found: str = ""):
        super().__init__(message)
        self.expected = expected
        self.found = found


class ResolveError(VarSightError):
    """A c. or p. description could not be placed on the genome."""

    def __init__(self, message: str, tried: list[str] | None = None):
        super().__init__(message)
        self.tried = tried or []


class DataError(VarSightError):
    """Local data is missing, damaged or the wrong version."""


class SetupError(VarSightError):
    """Download or build of the local data failed."""
