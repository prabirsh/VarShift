"""Run VarSight end to end on the small real data set kept in demo/, with no download."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

from ..annotate import Annotator
from ..batch import annotate_list, annotate_vcf, vcf_columns
from ..errors import DataError
from ..genome import open_build_genome
from ..output import ResultWriter
from ..store import DB_DIR, GENOME_DIR, MANIFEST, SCHEMA_VERSION, DataStore
from . import builders
from .selftest import run_selftest
from .manager import utc_now

BUILDS = ("hg19", "hg38")


def find_demo_dir(explicit: str | None = None) -> Path:
    candidates = [Path(explicit)] if explicit else []
    candidates += [Path(__file__).resolve().parents[3] / "demo", Path.cwd() / "demo"]
    for candidate in candidates:
        if (candidate / "sources").is_dir():
            return candidate
    raise DataError("cannot find the demo folder; run this from the VarSight repository, or add --demo-dir")


def build_demo_store(demo_dir: Path, data_dir: Path) -> dict:
    """Build a complete data folder from demo/sources using the same builders as a real setup."""
    sources = demo_dir / "sources"
    db_dir = data_dir / DB_DIR
    genome_dir = data_dir / GENOME_DIR
    db_dir.mkdir(parents=True, exist_ok=True)
    genome_dir.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((sources / "provenance.json").read_text(encoding="utf-8"))
    files: dict[str, dict] = {}

    def record(name: str, stats: dict, version_key: str) -> None:
        files[name] = {
            "policy": "demo",
            "kind": name.split(".")[0],
            "built_at": utc_now(),
            "source": {"version": provenance["versions"].get(version_key, provenance["date"]), "name": version_key},
            "stats": stats,
            "files": {},
        }

    for build in BUILDS:
        fasta = genome_dir / f"{build}.fa"
        stats = builders.build_genome(sources / f"{build}.windows.fa.gz", fasta)
        record(f"genome.{build}", stats, f"genome.{build}")
    manifest = {"schema": SCHEMA_VERSION, "builds": list(BUILDS), "created": utc_now(), "demo": True, "files": files}

    for build in BUILDS:
        stats = builders.build_transcripts(
            {"refseq": sources / f"refseq.{build}.gff.gz", "mane": sources / "mane.summary.txt.gz"},
            db_dir / f"transcripts.{build}.sqlite",
            build,
        )
        record(f"transcripts.{build}", stats, f"refseq.{build}")
        genome = open_build_genome(genome_dir, build)
        try:
            stats = builders.build_clinvar(
                sources / f"clinvar.{build}.vcf.gz",
                sources / "variant_summary.txt.gz",
                db_dir / f"clinvar.{build}.sqlite",
                build,
                genome,
            )
        finally:
            genome.close()
        record(f"clinvar.{build}", stats, f"clinvar.{build}")
    for source, target in (("hg19", "hg38"), ("hg38", "hg19")):
        stats = builders.build_liftover(sources / f"chain.{source}_{target}.gz", db_dir / f"liftover.{source}_{target}.sqlite", source, target)
        record(f"liftover.{source}_{target}", stats, f"chain.{source}_{target}")
    stats = builders.build_genes(
        {
            "hgnc": sources / "hgnc.txt.gz",
            "clingen_validity": sources / "clingen_validity.csv.gz",
            "clingen_dosage": sources / "clingen_dosage.tsv.gz",
            "gnomad_constraint": sources / "gnomad_constraint.tsv.gz",
            "hpo_diseases": sources / "hpo_diseases.txt.gz",
        },
        db_dir / "genes.sqlite",
    )
    record("genes", stats, "hgnc")
    with builders.open_text(sources / "pfam_regions.tsv.gz") as pfam:
        stats = builders.build_proteins(sources / "uniprot.dat.gz", pfam, db_dir / "proteins.sqlite")
    record("proteins", stats, "uniprot")
    stats = builders.build_civic(sources / "civic_variants.tsv.gz", sources / "civic_evidence.tsv.gz", db_dir / "civic.sqlite")
    record("civic", stats, "civic")
    for name, entry in files.items():
        if name.startswith("genome."):
            fasta = genome_dir / f"{name.split('.')[1]}.fa"
            paths = [fasta, Path(f"{fasta}.fai")]
        else:
            paths = [db_dir / f"{name}.sqlite"]
        entry["files"] = {path.name: {"size": path.stat().st_size} for path in paths}
    manifest["updated"] = utc_now()
    (data_dir / MANIFEST).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def run_demo(demo_dir_text: str | None, out: Path, force: bool = False) -> int:
    demo_dir = find_demo_dir(demo_dir_text)
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="varsight-demo-"))
    try:
        print(f"Building a small data set from {demo_dir / 'sources'}. Nothing is downloaded.")
        started = time.time()
        build_demo_store(demo_dir, work)
        print(f"  built in {time.time() - started:.1f} s\n")

        store = DataStore(work)
        try:
            annotator = Annotator(store)
            results = run_selftest(store)
            failed = [r for r in results if not r["passed"]]
            print(f"Known variants: {len(results) - len(failed)} of {len(results)} passed")
            for result in failed:
                print(f"  FAIL {result['name']}: {result['detail']}")

            print("\nAnnotating demo/variants.txt (p., c., g. and rs inputs)")
            for suffix in ("json", "tsv", "xlsx"):
                target = out / f"variants.{suffix}"
                writer = ResultWriter(str(target), suffix, force=force)
                counts = annotate_list(annotator, demo_dir / "variants.txt", writer, None)
                writer.close()
                print(f"  wrote {target} ({counts['results']} results, {counts['errors']} not annotated)")

            print("\nAnnotating demo/demo.hg19.vcf")
            vcf = demo_dir / "demo.hg19.vcf"
            for suffix in ("tsv", "xlsx"):
                target = out / f"demo.hg19.annotated.{suffix}"
                writer = ResultWriter(str(target), suffix, vcf_columns(vcf), force)
                counts = annotate_vcf(annotator, vcf, writer, None)
                writer.close()
                print(f"  wrote {target} ({counts['annotated']} annotated, {counts['ref_mismatch']} REF mismatches)")

            print("\nOne result in full:")
            first = annotator.annotate_text("BRAF V600E", "hg19")[0]
            _print_summary(first)
        finally:
            store.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"\nDemo finished. Open {out} to see the files.")
    return 1 if failed else 0


def _print_summary(record: dict) -> None:
    clinvar = record["evidence"]["clinvar"]
    rows = [
        ("input", record["input"]["raw"]),
        ("hg19", record["variant"]["hg19"]["key"]),
        ("hg38", record["variant"]["hg38"]["key"]),
        ("gene", record["gene"]["symbol"]),
        ("HGVS c.", record["hgvs"]["c"]),
        ("HGVS p.", record["hgvs"]["p"]),
        ("effect", record["effect"]["type"]),
        ("ClinVar", f"{clinvar.get('id')}: {clinvar.get('classification')}"),
        ("candidates", record["candidate"]["of"]),
    ]
    for label, value in rows:
        print(f"  {label:11} {value}")
    print("  links:")
    for link in record["links"]:
        print(f"    {link['name']:20} {link['label']:13} {link['url']}")
    sys.stdout.flush()

