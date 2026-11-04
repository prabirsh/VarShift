"""The background half of the update notice. Started by updates.py; asks sources and saves the answer."""

from __future__ import annotations

import sys
from pathlib import Path

from ..errors import VarShiftError
from ..updates import write_status
from .fetch import Fetcher
from .manager import DataManager, load_catalog


def run(data_dir: Path) -> int:
    catalog = load_catalog()
    # Short limits: this runs unattended, and a slow source can wait until tomorrow.
    fetcher = Fetcher(catalog["allowed_hosts"], retry_delay=1.0, attempts=2, timeout=20)
    manager = DataManager(data_dir, catalog, fetcher)
    try:
        items = manager.check_updates()
    except VarShiftError:
        return 1
    waiting = [
        {
            "artifact": item.artifact,
            "title": catalog["artifacts"][item.artifact].get("title") or item.artifact,
            "policy": item.policy,
        }
        for item in items
        if item.changed and item.policy in ("auto", "approval")
    ]
    write_status(data_dir, waiting)
    return 0


if __name__ == "__main__":
    sys.exit(run(Path(sys.argv[1])))
