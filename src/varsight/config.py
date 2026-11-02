"""Saved settings and where the data folder is."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

from .errors import InputError

_WINDOWS_PATH = re.compile(r"^([A-Za-z]):[\\/]*(.*)$")


def _config_file() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME")
    root = Path(base) if base else Path.home() / ".config"
    return root / "varsight" / "config.json"


def load_settings() -> dict:
    try:
        return json.loads(_config_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def save_setting(name: str, value) -> None:
    settings = load_settings()
    settings[name] = value
    target = _config_file()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")


def default_data_dir() -> Path:
    return Path.home() / "varsight-data"


def resolve_data_dir(explicit: str | None = None) -> Path:
    """The data= option, then $VARSIGHT_DATA, then the saved setting, then a default."""
    if explicit:
        return _given_folder(explicit)
    env = os.environ.get("VARSIGHT_DATA")
    if env:
        return _given_folder(env)
    saved = load_settings().get("data_dir")
    return Path(saved).expanduser() if saved else default_data_dir()


def _given_folder(text: str) -> Path:
    """On Linux, I:\\data would quietly become a new folder of that odd name in the current directory."""
    windows = _WINDOWS_PATH.match(text) if os.name != "nt" else None
    if windows:
        drive, rest = windows.groups()
        raise InputError(
            f"{text} looks like a Windows folder, which only works in PowerShell or Command Prompt. "
            f"In WSL, write it like /mnt/{drive.lower()}/{rest.replace(chr(92), '/')}"
        )
    return Path(text).expanduser().resolve()


def is_windows_mount(path: Path) -> bool:
    """True for a Windows drive seen from WSL. SQLite and memory maps run far slower there."""
    if os.name == "nt":
        return False
    parts = Path(path).resolve().parts
    return len(parts) >= 3 and parts[1] == "mnt" and len(parts[2]) == 1
