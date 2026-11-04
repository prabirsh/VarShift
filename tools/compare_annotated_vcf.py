#!/usr/bin/env python3
"""Compare VarShift's names with an annotated VCF you already trust (VEP-style CSQ field).

Runs entirely on this computer. It prints counts and writes differences to a local file you choose.
Differences list transcript names only, with no sample IDs or genotypes, but treat that file as
private anyway, because it shows which changes the sample carries.

    python3 tools/compare_annotated_vcf.py sample.annotated.vcf.gz differences=diff.txt

Options: build=hg19 or hg38 (read from the VCF header by default), data=FOLDER, differences=FILE.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from varshift.builds import canonical_contig  # noqa: E402
from varshift.config import resolve_data_dir  # noqa: E402
from varshift.consequence import SequenceCache, annotate  # noqa: E402
from varshift.errors import VarShiftError  # noqa: E402
from varshift.pipeline.builders import open_text  # noqa: E402
from varshift.pipeline.validate import comparable  # noqa: E402
from varshift.store import DataStore  # noqa: E402
from varshift.vcf import VcfReader  # noqa: E402
from varshift.variant import normalize  # noqa: E402

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
    options = dict(item.split("=", 1) for item in sys.argv[1:] if "=" in item)
    inputs = [item for item in sys.argv[1:] if "=" not in item]
    if len(inputs) != 1 or set(options) - {"build", "data", "differences"}:
        print(__doc__, file=sys.stderr)
        return 2
    vcf = Path(inputs[0])
    differences_file = Path(options["differences"]) if "differences" in options else None

    columns = csq_columns(vcf)
    for needed in ("HGVSc", "HGVSp"):
        if needed not in columns:
            raise SystemExit(f"the CSQ field has no {needed} column")
    c_at, p_at = columns.index("HGVSc"), columns.index("HGVSp")

    store = DataStore(resolve_data_dir(options.get("data")))
    counts: Counter = Counter()
    differences: list[str] = []
    with VcfReader(vcf) as reader:
        build = options.get("build") or reader.header.build or reader.header.named_build
        if build not in store.builds:
            raise SystemExit(f"cannot use build {build}; add build=hg19 or build=hg38, or install it with varshift setup")
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
            except VarShiftError:
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
                    differences.append(f"{accession}\tc.\tother tool {change}\tVarShift {ours}")
                theirs_p = fields[p_at].split(":", 1)[-1].replace("%3D", "=")
                if theirs_p:
                    counts[f"p_{kind}_total"] += 1
                    if comparable(effect.hgvs_p) == comparable(theirs_p):
                        counts[f"p_{kind}_same"] += 1
                    else:
                        differences.append(f"{accession}\tp.\tother tool {theirs_p}\tVarShift {effect.hgvs_p}")
    store.close()
    for name, value in sorted(counts.items()):
        print(f"{name:40} {value:,}")
    for level in ("c_snv", "c_indel", "p_snv", "p_indel"):
        total = counts[f"{level}_total"]
        if total:
            print(f"{level} agreement: {100 * counts[f'{level}_same'] / total:.2f}%")
    if differences_file:
        differences_file.write_text("\n".join(differences) + "\n", encoding="utf-8")
        differences_file.chmod(0o600)
        print(f"{len(differences):,} differences written to {differences_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
