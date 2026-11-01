"""Build the reference links for one annotation from the templates in data/links.json."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from urllib.parse import quote, urlsplit


@lru_cache(maxsize=1)
def templates() -> tuple[dict, ...]:
    text = resources.files("varsight").joinpath("data/links.json").read_text(encoding="utf-8")
    return tuple(json.loads(text)["links"])


def link_names() -> list[str]:
    return [template["name"] for template in templates()]


def _ready(requires: list[str], values: dict) -> bool:
    return all(values.get(name) not in (None, "") for name in requires)


def build_links(values: dict) -> list[dict]:
    """One entry per database. An entry with url None says why there is no link."""
    build = values.get("build", "")
    # A colon stays readable (HGNC:1097). A slash is quoted, so a value can never add a path segment.
    quoted = {key: quote(str(value), safe=":") for key, value in values.items() if value not in (None, "")}
    links = []
    for template in templates():
        entry = {"name": template["name"], "label": template["label"], "url": None}
        note = template.get("note")
        if _ready(template["requires"], values):
            url = template.get("url") or template[f"url_{build}"]
            if template.get("raw"):
                entry["url"] = _safe_raw_url(values[template["requires"][0]])
            else:
                entry["url"] = url.format_map(quoted)
        elif "fallback" in template and _ready(template["fallback"]["requires"], values):
            entry["url"] = template["fallback"]["url"].format_map(quoted)
            entry["label"] = "built"
            note = template["fallback"].get("note")
        else:
            note = template.get("missing_note") or "Not enough information for this link."
        if note:
            entry["note"] = note
        links.append(entry)
    return links


def _safe_raw_url(url: str) -> str | None:
    """Local data can supply a whole URL. Only plain https links are allowed through."""
    parts = urlsplit(str(url))
    return str(url) if parts.scheme == "https" and parts.netloc else None
