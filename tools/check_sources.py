#!/usr/bin/env python3
"""Ask every source in the catalog whether its file is still there. Header requests only."""

from __future__ import annotations

import json
import sys
from importlib import resources

from varsight.errors import SetupError
from varsight.pipeline.fetch import Fetcher


def main() -> int:
    catalog = json.loads(resources.files("varsight").joinpath("data/sources.json").read_text(encoding="utf-8"))
    fetcher = Fetcher(catalog["allowed_hosts"])
    failed = 0
    for source_id, source in sorted(catalog["sources"].items()):
        try:
            if "list_url" in source:
                url = fetcher.resolve_listing(source["list_url"], source["pattern"])
            else:
                url = source.get("url") or source["url_template"].format(chrom="21")
            remote = fetcher.head(url)
            size = f"{remote.size / 1024**2:,.1f} MB" if remote.size else "size not given"
            limit = source["max_mb"] * 1024**2
            if remote.size and remote.size > limit:
                failed += 1
                print(f"TOO BIG  {source_id:22} {size}, limit {source['max_mb']} MB: {url}")
                continue
            print(f"ok       {source_id:22} {size:>14}  {remote.last_modified}")
        except SetupError as exc:
            failed += 1
            print(f"FAILED   {source_id:22} {exc}")
    print(f"{len(catalog['sources']) - failed} of {len(catalog['sources'])} sources reachable")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
