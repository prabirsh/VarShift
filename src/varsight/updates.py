"""The quick notice that says when newer database versions exist.

Annotation must stay fast and must not depend on the network, so nothing here waits for a connection.
Reading the notice is one small file read. When the saved answer is over a day old, a separate
background process asks each source whether its file changed, and the answer is shown next time. Those
requests are headers only, go to the sources' own servers, and carry nothing about any variant.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

from .config import load_settings

STATUS_FILE = "update_status.json"
STARTED_FILE = "update_check.started"
CHECK_EVERY = timedelta(hours=24)
RETRY_AFTER = timedelta(hours=1)
OLD_DATA = timedelta(days=30)


def status_path(data_dir: Path) -> Path:
    return Path(data_dir) / STATUS_FILE


def read_status(data_dir: Path) -> dict:
    try:
        return json.loads(status_path(data_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def write_status(data_dir: Path, items: list[dict]) -> None:
    body = {"checked_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "items": items}
    target = status_path(data_dir)
    temp = target.with_suffix(".tmp")
    temp.write_text(json.dumps(body, indent=1), encoding="utf-8")
    os.replace(temp, target)


def clear_applied(data_dir: Path, applied: list[str]) -> None:
    """Drop databases that were just updated from the saved notice."""
    status = read_status(data_dir)
    if status:
        write_status(data_dir, [item for item in status.get("items", []) if item["artifact"] not in applied])


def _age(stamp: str | None, now: datetime) -> timedelta | None:
    if not stamp:
        return None
    try:
        return now - datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def notice_lines(data_dir: Path, manifest: dict) -> list[str]:
    """The message to show before annotating. Empty when there is nothing to say."""
    if manifest.get("demo") or load_settings().get("updates") == "off":
        return []
    now = datetime.now(timezone.utc)
    items = read_status(data_dir).get("items", [])
    automatic = sorted({item["title"] for item in items if item.get("policy") == "auto"})
    approval = sorted({item["title"] for item in items if item.get("policy") == "approval"})
    lines = []
    if automatic:
        lines.append(f"Update available: {', '.join(automatic)}. These results use your installed versions.")
        lines.append("  To update: varsight update    To hide this: varsight config updates=off")
    if approval:
        lines.append(f"New version available, your choice: {', '.join(approval)}. It can change c. and p. names.")
        lines.append("  To apply: varsight update approve=transcripts")
    age = _age(manifest.get("updated"), now)
    if not lines and age is not None and age > OLD_DATA:
        lines.append(f"Your databases were last updated {age.days} days ago, and no check has been possible since.")
        lines.append("  To update: varsight update")
    return lines


def _claim_check(data_dir: Path) -> bool:
    """Allow one try an hour, so a script calling varsight in a loop, or a computer offline, starts just one."""
    marker = Path(data_dir) / STARTED_FILE
    try:
        if marker.exists() and time.time() - marker.stat().st_mtime < RETRY_AFTER.total_seconds():
            return False
        marker.touch()
    except OSError:
        return False  # a data folder that can't be written couldn't save the answer either
    return True


def check_in_background(data_dir: Path, manifest: dict) -> bool:
    """Start the update check if it is due, and return at once. Never raises."""
    if manifest.get("demo") or load_settings().get("updates") == "off":
        return False
    checked = _age(read_status(data_dir).get("checked_at"), datetime.now(timezone.utc))
    if checked is not None and checked < CHECK_EVERY or not _claim_check(data_dir):
        return False
    arguments =[sys.executable, "-m", "varsight.pipeline.updatecheck", str(data_dir)]
    environment = dict(os.environ)
    package_root = str(Path(__file__).resolve().parents[1])
    environment["PYTHONPATH"] = package_root
    options: dict = {}
    if os.name == "nt":
        options["creationflags"] = subprocess.DETACHED_PROCESS | subprocess.CREATE_NO_WINDOW
    else:
        options["start_new_session"] = True
    try:
        subprocess.Popen(
            arguments,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=environment,
            close_fds=True,
            **options,
        )
    except OSError:
        return False
    return True
