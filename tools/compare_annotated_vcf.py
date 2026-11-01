#!/usr/bin/env python3
"""Compare VarSight's names with an annotated VCF you already trust (VEP-style CSQ field).

Runs entirely on this computer. It prints counts and writes differences to a local file you choose.
Differences list transcript names only, with no sample IDs or genotypes, but treat that file as
private anyway, because it shows which changes the sample carries.

    python3 tools/compare_annotated_vcf.py sample.annotated.vcf.gz --differences diff.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from varsight.builds import canonical_contig  # noqa: E402
from varsight.config import resolve_data_dir  # noqa: E402
from varsight.consequence import SequenceCache, annotate  # noqa: E402
from varsight.errors import VarSightError  # noqa: E402
from varsight.pipeline.builders import open_text  # noqa: E402
from varsight.pipeline.validate import comparable  # noqa: E402
from varsight.store import DataStore  # noqa: E402
from varsight.vcf import VcfReader  # noqa: E402
from varsight.variant import normalize  # noqa: E402

_CSQ_RE = re.compile(r'^##INFO=<ID=CSQ,.*Format: ([^"]+)"')


def csq_columns(path: Path) -> list[str]:
    with open_text(path) as handle:
        for line in handle:
            match = _CSQ_RE.match(line)
            if match:
                return match.group(1).split("|")
            if not line.startswith("##"):
                break
    raise SystemExit(f"{path} has no CSQ field in its header, so there is nothing to compare")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=Path)
    parser.add_argument("--build", help="hg19 or hg38; read from the VCF header by default")
    parser.add_argument("--data-dir")
    parser.add_argument("--differences", type=Path, help="write each differing name to this file")
    args = parser.parse_args()

    columns = csq_columns(args.vcf)
    for needed in ("HGVSc", "HGVSp"):
        if needed not in columns:
            raise SystemExit(f"the CSQ field has no {needed} column")
    c_at, p_at = columns.index("HGVSc"), columns.index("HGVSp")

    store = DataStore(resolve_data_dir(args.data_dir))
    counts: Counter = Counter()
    differences: list[str] = []
    with VcfReader(args.vcf) as reader:
        build = args.build or reader.header.build or reader.header.named_build
        if build not in store.builds:
            raise SystemExit(f"cannot use build {build}; add --build, or install it with varsight setup")
        genome = store.genome(build)
        cache = SequenceCache(genome)
        for record in reader.records():
            if record.allele_count > 1 or not record.is_plain or len(record.fields) < 8:
                counts["records_skipped"] += 1
                continue
            info = dict(item.split("=", 1) for item in record.fields[7].split(";") if "=" in item)
            if "CSQ" not in info:
                continue
            try:
                variant = normalize(genome, canonical_contig(record.chrom, build), record.pos, record.ref, record.alt)
            except VarSightError:
                counts["records_with_ref_mismatch"] += 1
                continue
            counts["records_compared"] += 1
            start, end, _ = variant.edit()
            nearby = {t.accession: t for t in store.transcripts_at(build, variant.contig, start, max(end, start + 1), flank=5000, source=None)}
            kind = "snv" if len(variant.ref) == len(variant.alt) == 1 else "indel"
            for entry in info["CSQ"].split(","):
                fields = entry.split("|")
                if len(fields) <= max(c_at, p_at) or ":" not in fields[c_at]:
                    continue
                accession, change = fields[c_at].split(":", 1)
                tx = nearby.get(accession)
                if tx is None:
                    counts["transcript_or_version_not_installed"] += 1
                    continue
                effect = annotate(genome, variant, tx, cache)
                ours = effect.hgvs_c.split(":", 1)[1] if effect.hgvs_c else None
                counts[f"c_{kind}_total"] += 1
                if ours == change:
                    counts[f"c_{kind}_same"] += 1
                else:
                    differences.append(f"{accession}\tc.\tother tool {change}\tVarSight {ours}")
                theirs_p = fields[p_at].split(":", 1)[-1].replace("%3D", "=")
                if theirs_p:
                    counts[f"p_{kind}_total"] += 1
                    if comparable(effect.hgvs_p) == comparable(theirs_p):
                        counts[f"p_{kind}_same"] += 1
                    else:
                        differences.append(f"{accession}\tp.\tother tool {theirs_p}\tVarSight {effect.hgvs_p}")
    store.close()
    for name, value in sorted(counts.items()):
        print(f"{name:40} {value:,}")
    for level in ("c_snv", "c_indel", "p_snv", "p_indel"):
        total = counts[f"{level}_total"]
        if total:
            print(f"{level} agreement: {100 * counts[f'{level}_same'] / total:.2f}%")
    if args.differences:
        args.differences.write_text("\n".join(differences) + "\n", encoding="utf-8")
        args.differences.chmod(0o600)
        print(f"{len(differences):,} differences written to {args.differences}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
