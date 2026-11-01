#!/usr/bin/env python3
"""Rebuild demo/sources: real data, cut down to a few genes, in each source's own file format.

Maintainers only. This reads from the same public sources as a real setup, with two exceptions that
keep the transfer small: sequence windows come from the UCSC REST API, and Pfam positions from the
InterPro API. Neither is ever used by VarSight itself.

    python3 tools/build_demo_data.py --out demo/sources
    python3 tools/build_demo_data.py --out /tmp/validation --genes-file genes.txt --no-cap
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import random
import re
import struct
import sys
import time
import urllib.request
import zlib
from datetime import date
from pathlib import Path

DEMO_GENES = ["BRAF", "KRAS", "NRAS", "EGFR", "TP53", "PIK3CA", "IDH1", "BRCA1", "BRCA2", "JAK2", "KIT", "PTEN"]
FLANK = 6000
LINE = 60
UA = {"User-Agent": "VarSight-demo-builder/0.1"}
REFSEQ = {
    "hg19": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz",
    "hg38": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz",
}
CLINVAR = {
    "hg19": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    "hg38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}
SUMMARY = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MANE_LIST = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
HGNC = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
CLINGEN_VALIDITY = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINGEN_DOSAGE = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
CONSTRAINT = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
HPO = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt"
CIVIC_VARIANTS = "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv"
CIVIC_EVIDENCE = "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
CHAINS = {
    ("hg19", "hg38"): "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz",
    ("hg38", "hg19"): "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz",
}

transferred = 0


def get(url: str, headers: dict | None = None, attempts: int = 4) -> bytes:
    global transferred
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(request, timeout=120) as response:
                body = response.read()
            transferred += len(body)
            return body
        except Exception as exc:  # noqa: BLE001 - a maintainer tool; retry anything transient
            if attempt == attempts - 1:
                raise RuntimeError(f"{url}: {exc}") from exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(url)


def stream_gz_lines(url: str):
    """Decompress on the fly and yield text lines, so large files are never stored."""
    global transferred
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=300) as response:
        decoder = zlib.decompressobj(31)
        tail = b""
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            transferred += len(block)
            data = decoder.decompress(block)
            while decoder.unused_data:
                unused = decoder.unused_data
                decoder = zlib.decompressobj(31)
                data += decoder.decompress(unused)
            tail += data
            lines = tail.split(b"\n")
            tail = lines.pop()
            for line in lines:
                yield line.decode("utf-8", "replace") + "\n"
        if tail:
            yield tail.decode("utf-8", "replace") + "\n"


def write_gz(path: Path, lines) -> int:
    count = 0
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"), mtime=0) as raw:
        with io.TextIOWrapper(raw, encoding="utf-8", newline="\n") as handle:
            for line in lines:
                handle.write(line)
                count += 1
    return count


# RefSeq
def filter_refseq(build: str, genes: set[str], out: Path) -> dict[str, tuple[str, int, int]]:
    """Keep every line of the chosen genes, and return their spans."""
    spans: dict[str, tuple[str, int, int]] = {}
    kept = []
    total = 0
    for line in stream_gz_lines(REFSEQ[build]):
        if line[0] == "#":
            if line.startswith(("##gff", "#!")):
                kept.append(line)
            continue
        total += 1
        match = re.search(r"(?:^|;)gene=([^;\n]+)", line)
        if not match or match.group(1) not in genes:
            continue
        fields = line.split("\t")
        if not fields[0].startswith("NC_"):
            continue
        kept.append(line)
        if fields[2] == "gene":
            spans[match.group(1)] = (fields[0], int(fields[3]), int(fields[4]))
    write_gz(out, kept)
    print(f"  refseq {build}: {len(kept):,} of {total:,} lines kept", file=sys.stderr)
    return spans


_ACCESSION_TO_CHROM = {f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)} | {"NC_000023": "chrX", "NC_000024": "chrY"}


def windows_for(spans: dict[str, tuple[str, int, int]]) -> dict[str, list[tuple[int, int]]]:
    per_chrom: dict[str, list[tuple[int, int]]] = {}
    for accession, start, end in spans.values():
        chrom = _ACCESSION_TO_CHROM[accession.split(".")[0]]
        per_chrom.setdefault(chrom, []).append((max(start - FLANK, 1), end + FLANK))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, pairs in per_chrom.items():
        pairs.sort()
        out = [list(pairs[0])]
        for start, end in pairs[1:]:
            if start <= out[-1][1]:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end])
        merged[chrom] = [tuple(p) for p in out]
    return merged


def write_windows(build: str, windows: dict, out: Path) -> None:
    with gzip.GzipFile(filename="", mode="wb", fileobj=open(out, "wb"), mtime=0) as raw:
        for chrom in sorted(windows):
            for start, end in windows[chrom]:
                url = f"https://api.genome.ucsc.edu/getData/sequence?genome={build};chrom={chrom};start={start - 1};end={end}"
                dna = json.loads(get(url))["dna"]
                raw.write(f">{chrom}:{start}-{end}\n".encode())
                for i in range(0, len(dna), LINE):
                    raw.write(dna[i : i + LINE].encode() + b"\n")
                time.sleep(1.1)  # UCSC asks for no more than one request per second
        # The whole mitochondrial sequence, which is also how rCRS reaches hg19.
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome={build};chrom=chrM"
        dna = json.loads(get(url))["dna"]
        raw.write(f">chrM\n".encode())
        for i in range(0, len(dna), LINE):
            raw.write(dna[i : i + LINE].encode() + b"\n")


# ClinVar through its tabix index
def _reg2bins(beg: int, end: int) -> list[int]:
    bins = [0]
    end -= 1
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(offset + (beg >> shift), offset + (end >> shift) + 1))
    return bins


class Tabix:
    def __init__(self, url: str):
        raw = gzip.decompress(get(url + ".tbi"))
        self.url = url
        magic, n_ref, _, _, _, _, _, _, names_len = struct.unpack_from("<4s8i", raw, 0)
        if magic != b"TBI\x01":
            raise RuntimeError("not a tabix index")
        offset = 36
        names = raw[offset : offset + names_len].split(b"\0")[:-1]
        offset += names_len
        self.refs: dict[str, tuple[dict, list]] = {}
        for index in range(n_ref):
            n_bin = struct.unpack_from("<i", raw, offset)[0]
            offset += 4
            bins = {}
            for _ in range(n_bin):
                bin_number, n_chunk = struct.unpack_from("<Ii", raw, offset)
                offset += 8
                chunks = struct.unpack_from(f"<{n_chunk * 2}Q", raw, offset)
                offset += 16 * n_chunk
                bins[bin_number] = list(zip(chunks[::2], chunks[1::2]))
            n_intv = struct.unpack_from("<i", raw, offset)[0]
            offset += 4
            linear = struct.unpack_from(f"<{n_intv}Q", raw, offset)
            offset += 8 * n_intv
            self.refs[names[index].decode()] = (bins, linear)

    def fetch(self, contig: str, beg: int, end: int) -> list[str]:
        """Records overlapping the 1-based closed range."""
        if contig not in self.refs:
            return []
        bins, linear = self.refs[contig]
        floor = linear[min((beg - 1) >> 14, len(linear) - 1)] if linear else 0
        chunks = [c for b in _reg2bins(beg - 1, end) for c in bins.get(b, []) if c[1] > floor]
        if not chunks:
            return []
        start = min(c[0] for c in chunks) >> 16
        stop = (max(c[1] for c in chunks) >> 16) + 65536
        body = get(self.url, {"Range": f"bytes={start}-{stop}"})
        text = bytearray()
        decoder = zlib.decompressobj(31)
        data = body
        while data:
            try:
                text += decoder.decompress(data)
            except zlib.error:
                break
            if decoder.unused_data:
                data = decoder.unused_data
                decoder = zlib.decompressobj(31)
            elif decoder.eof:
                break
            else:
                break
        lines = []
        for line in bytes(text).decode("utf-8", "replace").split("\n"):
            fields = line.split("\t", 2)
            if len(fields) >= 2 and fields[0] == contig and fields[1].isdigit() and beg <= int(fields[1]) <= end:
                lines.append(line + "\n")
        return lines


def filter_clinvar(build: str, windows: dict, out: Path, cap: int | None, keep_ids: set[int] | None = None) -> set[int]:
    """Cut ClinVar down to the windows. keep_ids makes the second build hold the same records as the first."""
    url = CLINVAR[build]
    head = zlib.decompressobj(31).decompress(get(url, {"Range": "bytes=0-65535"}))
    header = [line + "\n" for line in head.decode("utf-8", "replace").split("\n") if line.startswith("#")]
    index = Tabix(url)
    seen: set[str] = set()
    rows: list[str] = []
    for chrom in sorted(windows):
        for start, end in windows[chrom]:
            rows.extend(index.fetch(chrom[3:], start, end))
        rows.extend(index.fetch("MT", 1, 16569) if chrom == sorted(windows)[0] else [])
    unique = []
    for row in rows:
        fields = row.split("\t", 5)
        marker = f"{fields[0]}-{fields[1]}-{fields[3]}-{fields[4]}"
        if marker not in seen:
            seen.add(marker)
            unique.append(row)
    if keep_ids is not None:
        unique = [row for row in unique if int(row.split("\t", 3)[2]) in keep_ids]
    elif cap:
        unique = _cap_rows(unique, cap)
    unique.sort(key=lambda r: (_chrom_order(r.split("\t", 1)[0]), int(r.split("\t", 2)[1])))
    write_gz(out, header + unique)
    print(f"  clinvar {build}: {len(unique):,} records", file=sys.stderr)
    return {int(r.split("\t", 3)[2]) for r in unique}


def _chrom_order(name: str) -> int:
    return {"X": 23, "Y": 24, "MT": 25}.get(name, int(name) if name.isdigit() else 26)


# Well-known records stay in the demo even when a gene is capped: the variants named in demo/variants.txt
# and tests, by ClinVar variation ID, plus every allele of rs113488022.
ALWAYS_KEEP = {
    13961, 40389, 12582, 13900, 16609, 12374, 13652, 156444, 13852, 14662, 7829, 177620, 17662, 266568, 9325, 9589,
}  # fmt: skip
_KEEP_RS = "RS=113488022"


def _cap_rows(rows: list[str], cap: int) -> list[str]:
    rng = random.Random(20261101)
    by_gene: dict[str, list[str]] = {}
    for row in rows:
        match = re.search(r"GENEINFO=([^:;\t]+)", row)
        by_gene.setdefault(match.group(1) if match else "", []).append(row)
    kept = []
    for _, items in sorted(by_gene.items()):
        rng.shuffle(items)
        chosen = items[:cap]
        chosen += [r for r in items[cap:] if int(r.split("\t", 3)[2]) in ALWAYS_KEEP or _KEEP_RS in r]
        kept.extend(chosen)
    return kept


def filter_summary(variation_ids: set[int], out: Path, limit_mb: int) -> int:
    """Read from the start of variant_summary up to limit_mb, keeping names for our records."""
    body = get(SUMMARY, {"Range": f"bytes=0-{limit_mb * 1024 * 1024 - 1}"})
    decoder = zlib.decompressobj(31)
    text = decoder.decompress(body).decode("utf-8", "replace")
    lines = text.split("\n")[:-1]
    header = lines[0].lstrip("#").split("\t")
    at = header.index("VariationID")
    kept = [lines[0] + "\n"]
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) > at and fields[at].isdigit() and int(fields[at]) in variation_ids:
            kept.append(line + "\n")
    write_gz(out, kept)
    print(f"  variant_summary: {len(kept) - 1:,} rows from the first {limit_mb} MB", file=sys.stderr)
    return len(kept) - 1


# gene and protein tables
def filter_table(url: str, genes: set[str], out: Path, column: str, skip_comment: bool = False) -> None:
    body = get(url)
    if body[:2] == b"\x1f\x8b":
        body = _gunzip_all(body)
    lines = body.decode("utf-8", "replace").split("\n")
    kept, at = [], None
    for line in lines:
        if not line:
            continue
        if at is None:
            kept.append(line + "\n")
            fields = line.lstrip("#").split("\t")
            if column in fields:
                at = fields.index(column)
            continue
        fields = line.split("\t")
        if len(fields) > at and fields[at].strip('"') in genes:
            kept.append(line + "\n")
    write_gz(out, kept)
    print(f"  {out.name}: {len(kept)} lines", file=sys.stderr)


def _gunzip_all(body: bytes) -> bytes:
    out = bytearray()
    data = body
    while data:
        decoder = zlib.decompressobj(31)
        out += decoder.decompress(data)
        data = decoder.unused_data
    return bytes(out)


def filter_validity(genes: set[str], out: Path) -> None:
    kept = []
    for line in get(CLINGEN_VALIDITY).decode("utf-8", "replace").splitlines():
        first = line.split('","', 1)[0].strip('"')
        if first in genes or first in ("GENE SYMBOL",) or first.startswith(("+", "CLINGEN", "FILE", "WEBPAGE")):
            kept.append(line + "\n")
    write_gz(out, kept)


def filter_civic(genes: set[str], variants_out: Path, evidence_out: Path) -> None:
    variants = get(CIVIC_VARIANTS).decode("utf-8", "replace").split("\n")
    header = variants[0].split("\t")
    gene_at, profile_at = header.index("gene"), header.index("single_variant_molecular_profile_id")
    kept, profiles = [variants[0] + "\n"], set()
    for line in variants[1:]:
        fields = line.split("\t")
        if len(fields) > gene_at and fields[gene_at] in genes:
            kept.append(line + "\n")
            profiles.add(fields[profile_at])
    write_gz(variants_out, kept)
    evidence = get(CIVIC_EVIDENCE).decode("utf-8", "replace").split("\n")
    header = evidence[0].split("\t")
    at = header.index("molecular_profile_id")
    kept = [evidence[0] + "\n"]
    for line in evidence[1:]:
        fields = line.split("\t")
        if len(fields) > at and fields[at] in profiles:
            # Keep only the columns VarSight reads; the free text is not needed for a demo.
            fields[header.index("evidence_statement")] = ""
            kept.append("\t".join(fields) + "\n")
    write_gz(evidence_out, kept)


def filter_chain(url: str, contigs: set[str], out: Path) -> None:
    kept, active = [], False
    for line in gzip.decompress(get(url)).decode().split("\n"):
        if line.startswith("chain"):
            active = line.split()[2] in contigs
        if active:
            kept.append(line + "\n")
    write_gz(out, kept)


def fetch_uniprot(accessions: list[str], out: Path) -> None:
    parts = []
    for accession in accessions:
        parts.append(get(f"https://rest.uniprot.org/uniprotkb/{accession}.txt").decode("utf-8", "replace"))
        time.sleep(0.2)
    write_gz(out, ["".join(parts)])


def fetch_pfam(accessions: list[str], out: Path) -> None:
    lines = ["pfamseq_acc\tseq_version\tcrc64\tmd5\tpfamA_acc\tseq_start\tseq_end\tali_start\tali_end\n"]
    for accession in accessions:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/?format=json"
        try:
            data = json.loads(get(url, {"Accept": "application/json"}))
        except (RuntimeError, ValueError):
            continue
        for entry in data.get("results", []):
            pfam = entry["metadata"]["accession"].upper().replace("PFAM:", "")
            for protein in entry.get("proteins", []):
                for location in protein.get("entry_protein_locations", []):
                    for fragment in location.get("fragments", []):
                        lines.append(f"{accession}\t1\t-\t-\t{pfam}\t{fragment['start']}\t{fragment['end']}\t{fragment['start']}\t{fragment['end']}\n")
        time.sleep(0.2)
    write_gz(out, lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--genes-file", type=Path, help="one gene per line, instead of the demo genes")
    parser.add_argument("--no-cap", action="store_true", help="keep every ClinVar record in the windows")
    parser.add_argument("--cap", type=int, default=1200, help="ClinVar records per gene when capped")
    parser.add_argument("--summary-mb", type=int, default=30, help="compressed MB of variant_summary to read")
    args = parser.parse_args()

    genes = set(DEMO_GENES)
    if args.genes_file:
        genes = {g.strip() for g in args.genes_file.read_text().splitlines() if g.strip()}
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    cap = None if args.no_cap else args.cap

    spans_by_build = {}
    ids: set[int] | None = None
    for build in ("hg19", "hg38"):
        print(f"{build}: RefSeq, sequence windows, ClinVar", file=sys.stderr)
        spans = filter_refseq(build, genes, out / f"refseq.{build}.gff.gz")
        spans_by_build[build] = spans
        windows = windows_for(spans)
        write_windows(build, windows, out / f"{build}.windows.fa.gz")
        found = filter_clinvar(build, windows, out / f"clinvar.{build}.vcf.gz", cap, ids if cap else None)
        ids = found if ids is None else ids | found

    print("names and gene tables", file=sys.stderr)
    filter_summary(ids, out / "variant_summary.txt.gz", args.summary_mb)
    listing = get(MANE_LIST).decode()
    mane_name = re.search(r'href="(MANE\.GRCh38\.v[0-9.]+\.summary\.txt\.gz)"', listing).group(1)
    filter_table(MANE_LIST + mane_name, genes, out / "mane.summary.txt.gz", "symbol")

    hgnc = get(HGNC).decode("utf-8", "replace").split("\n")
    header = hgnc[0].split("\t")
    symbol_at, uniprot_at = header.index("symbol"), header.index("uniprot_ids")
    kept, accessions = [hgnc[0] + "\n"], []
    for line in hgnc[1:]:
        fields = line.split("\t")
        if len(fields) > uniprot_at and fields[symbol_at] in genes:
            kept.append(line + "\n")
            if fields[uniprot_at]:
                accessions.append(fields[uniprot_at].strip('"').split("|")[0])
    write_gz(out / "hgnc.txt.gz", kept)

    filter_validity(genes, out / "clingen_validity.csv.gz")
    filter_table(CLINGEN_DOSAGE, genes, out / "clingen_dosage.tsv.gz", "Gene Symbol")
    filter_table(CONSTRAINT, genes, out / "gnomad_constraint.tsv.gz", "gene")
    filter_table(HPO, genes, out / "hpo_diseases.txt.gz", "gene_symbol")
    filter_civic(genes, out / "civic_variants.tsv.gz", out / "civic_evidence.tsv.gz")
    fetch_uniprot(accessions, out / "uniprot.dat.gz")
    fetch_pfam(accessions, out / "pfam_regions.tsv.gz")
    contigs = {_ACCESSION_TO_CHROM[s[0].split(".")[0]][3:] for s in spans_by_build["hg19"].values()}
    for pair, url in CHAINS.items():
        filter_chain(url, contigs, out / f"chain.{pair[0]}_{pair[1]}.gz")

    provenance = {
        "date": date.today().isoformat(),
        "genes": sorted(genes),
        "note": "Real public data cut down to these genes. See demo/README.md. Not for clinical use.",
        "versions": {},
        "transferred_mb": round(transferred / 2**20, 1),
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"done; {transferred / 2**20:.1f} MB transferred", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
