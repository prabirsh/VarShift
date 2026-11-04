"""Turn downloaded source files into the local databases.

Each builder streams its input and returns row counts for the report. Source formats change without
notice, so builders check for the columns they need and fail with a clear message instead of
quietly producing empty tables.
"""

from __future__ import annotations

import csv
import functools
import gzip
import io
import logging
import re
import shutil
import sqlite3
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from urllib.parse import unquote

from ..builds import HG19, canonical_contig, is_mitochondrial, primary_contigs
from ..errors import InputError, ReferenceMismatch, SetupError
from ..genome import Genome, write_fai
from ..hgvs import protein_style
from ..store import encode_blocks, encode_exons, ucsc_bin
from ..variant import normalize

log = logging.getLogger("varshift.build")

_BATCH = 20_000
csv.field_size_limit(16 * 1024 * 1024)


def open_text(path) -> io.TextIOBase:
    """Text lines from a plain or gzip file, by content rather than by file name."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="\n")
    return open(path, encoding="utf-8", errors="replace", newline="\n")


_open_builds: list[sqlite3.Connection] = []


def closes_on_failure(builder):
    """Close a half-built database when a builder stops. Windows can't remove a file that is still open."""

    @functools.wraps(builder)
    def run(*args, **kwargs):
        try:
            return builder(*args, **kwargs)
        finally:
            while _open_builds:
                _open_builds.pop().close()

    return run


def _new_db(path: Path, script: str) -> sqlite3.Connection:
    path.unlink(missing_ok=True)
    con = sqlite3.connect(path)
    _open_builds.append(con)
    # Building an index sorts on disk, by default on the system disk, which may have no room. Stay beside the data.
    folder = str(Path(path).resolve().parent).replace("'", "''")
    con.executescript(
        f"PRAGMA temp_store_directory='{folder}'; "
        "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA cache_size=-64000; PRAGMA temp_store=FILE;"
    )
    con.executescript(script)
    return con


def _finish_db(con: sqlite3.Connection, script: str = "") -> None:
    con.commit()
    if script:
        con.executescript(script)
    con.execute("ANALYZE")
    con.commit()
    con.close()


# genome
@closes_on_failure
def build_genome(fasta_gz: Path, out_fa: Path) -> dict:
    """Unpack the reference and index it. Line-based reads keep memory flat for a 3 GB file."""
    temp = out_fa.with_name(out_fa.name + ".tmp")
    needed = fasta_gz.stat().st_size * 4
    if shutil.disk_usage(out_fa.parent).free < needed:
        raise SetupError(f"unpacking {fasta_gz.name} needs about {needed / 1024**3:.1f} GB free in {out_fa.parent}")
    with open_binary(fasta_gz) as source, open(temp, "wb") as out:
        shutil.copyfileobj(source, out, 1024 * 1024)
    temp.replace(out_fa)
    records = write_fai(out_fa)
    if records == 0:
        raise SetupError(f"{out_fa.name} holds no sequences")
    return {"sequences": records, "bytes": out_fa.stat().st_size}


def open_binary(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


# liftover
@closes_on_failure
def build_liftover(chain_gz: Path, out: Path, source_build: str, target_build: str) -> dict:
    con = _new_db(
        out,
        """
        CREATE TABLE blocks(contig TEXT, start INT, end INT, target_contig TEXT, target_start INT,
                            target_end INT, target_strand INT, score INT);
        CREATE TABLE meta(contig TEXT PRIMARY KEY, value INT);
        """,
    )
    primary = set(primary_contigs())
    longest: Counter = Counter()
    batch: list[tuple] = []
    blocks = chains = skipped = 0
    current = None
    with open_text(chain_gz) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                current = None
                continue
            fields = line.split()
            if fields[0] == "chain":
                if len(fields) < 12:
                    raise SetupError(f"{chain_gz}: malformed chain header: {line[:80]}")
                source = _chain_contig(fields[2], source_build)
                target = _chain_contig(fields[7], target_build)
                if source not in primary or target not in primary or fields[4] != "+":
                    current = None
                    skipped += 1
                    continue
                current = {
                    "score": int(fields[1]),
                    "source": source,
                    "target": target,
                    "strand": 1 if fields[9] == "+" else -1,
                    "target_size": int(fields[8]),
                    "s": int(fields[5]),
                    "t": int(fields[10]),
                }
                chains += 1
                continue
            if current is None:
                continue
            size = int(fields[0])
            s_start = current["s"]
            t_start = current["t"]
            if current["strand"] > 0:
                t_low, t_high = t_start, t_start + size
            else:
                t_low, t_high = current["target_size"] - t_start - size, current["target_size"] - t_start
            batch.append(
                (current["source"], s_start, s_start + size, current["target"], t_low, t_high, current["strand"], current["score"])
            )
            longest[current["source"]] = max(longest[current["source"]], size)
            blocks += 1
            if len(fields) == 3:
                current["s"] += size + int(fields[1])
                current["t"] += size + int(fields[2])
            if len(batch) >= _BATCH:
                con.executemany("INSERT INTO blocks VALUES(?,?,?,?,?,?,?,?)", batch)
                batch.clear()
    con.executemany("INSERT INTO blocks VALUES(?,?,?,?,?,?,?,?)", batch)
    con.executemany("INSERT INTO meta VALUES(?,?)", longest.items())
    if blocks == 0:
        con.close()
        raise SetupError(f"{chain_gz} produced no usable alignment blocks")
    _finish_db(con, "CREATE INDEX blocks_pos ON blocks(contig, start);")
    return {"chains": chains, "blocks": blocks, "chains_skipped": skipped}


def _chain_contig(name: str, build: str) -> str:
    try:
        return canonical_contig(name, build)
    except InputError:
        return ""


# transcripts
_TRANSCRIPT_COLUMNS = (
    "accession, base, version, gene, gene_id, hgnc_id, contig, strand, tx_low, tx_high, bin, exons, "
    "cds_start, cds_end, biotype, protein, protein_base, source, mane, mane_accession, refseq_select, "
    "genome_note, translation_exception, alignment"
)


def _attributes(text: str) -> dict[str, str]:
    out = {}
    for item in text.split(";"):
        key, _, value = item.partition("=")
        if value:
            out[key] = value
    return out


def read_mane(path: Path) -> tuple[dict[str, tuple[str, str]], int]:
    """Map base accession to (MANE status, the exact MANE accession) for RefSeq and Ensembl."""
    mapping: dict[str, tuple[str, str]] = {}
    with open_text(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        _require_columns(path, header, ["RefSeq_nuc", "Ensembl_nuc", "MANE_status"])
        refseq, ensembl, status = (header.index(n) for n in ("RefSeq_nuc", "Ensembl_nuc", "MANE_status"))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(refseq, ensembl, status):
                continue
            for accession in (fields[refseq], fields[ensembl]):
                if accession:
                    mapping[accession.split(".")[0]] = (fields[status], accession)
    return mapping, len(mapping)


def read_refseq_alignments(lines: Iterable[str]) -> dict[tuple[str, str], list[tuple[list, list]]]:
    """RefSeq's alignment of each transcript that doesn't match the genome base for base.

    These cDNA_match rows come after the genes on each chromosome, so they get a pass of their own. One row
    is one exon. The exon rows can't stand in for them: where a transcript has an extra base, they skip two
    genome bases instead to keep the reading frame, which looks like an intron and numbers the rest wrongly.
    """
    rows: dict[tuple[str, str, str], list] = {}
    for line in lines:
        if "\tcDNA_match\t" not in line:
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or not fields[0].startswith("NC_"):
            continue
        attrs = _attributes(fields[8])
        target = attrs.get("Target", "").split()
        if len(target) != 4 or target[3] != "+":
            raise SetupError(f"RefSeq alignment: unexpected Target {attrs.get('Target')!r}. The source format has changed, so nothing was installed.")
        key = (fields[0], target[0], attrs.get("ID", ""))
        rows.setdefault(key, []).append((int(target[1]) - 1, int(target[2]), int(fields[3]) - 1, int(fields[4]), fields[6], attrs.get("Gap", "")))
    alignments: dict[tuple[str, str], list[tuple[list, list]]] = {}
    for (seqid, accession, _), group in rows.items():
        group.sort()
        exact = group[0][0] == 0 and all(a[1] == b[0] for a, b in zip(group, group[1:])) and not any(row[5] for row in group)
        if not exact:
            exons = [(row[2], row[3]) for row in group]
            blocks = [stretch for row in group for stretch in _aligned_stretches(*row)]
            alignments.setdefault((seqid, accession), []).append((exons, blocks))
    return alignments


def _aligned_stretches(index: int, last_index: int, start: int, end: int, strand: str, gap: str) -> list[tuple[int, int, int]]:
    """Split one alignment row at its gaps. The Gap operators follow the transcript, on either strand."""
    stretches = []
    position = start if strand == "+" else end
    step = 1 if strand == "+" else -1
    for operator in gap.split() or [f"M{end - start}"]:
        kind, size = operator[0], int(operator[1:])
        if kind == "M":
            stretches.append((position, position + size, index) if strand == "+" else (position - size, position, index))
        elif kind not in "ID":
            raise SetupError(f"RefSeq alignment: unknown Gap operator {operator!r}. The source format has changed, so nothing was installed.")
        if kind in "MI":
            index += size
        if kind in "MD":
            position += step * size
    if index != last_index or position != (end if strand == "+" else start):
        raise SetupError(f"RefSeq alignment: Gap {gap!r} doesn't add up. The source format has changed, so nothing was installed.")
    return stretches


def _alignment_for(record: dict, alignments: dict | None) -> tuple[list, list] | None:
    """The alignment on this part of the chromosome. A few transcripts align to two places on one."""
    low = min(start for start, _ in record["exons"])
    high = max(end for _, end in record["exons"])
    best, best_overlap = None, 0
    for exons, blocks in (alignments or {}).get((record["seqid"], record["accession"]), []):
        overlap = min(high, max(end for _, end in exons)) - max(low, min(start for start, _ in exons))
        if overlap > best_overlap:
            best, best_overlap = (exons, blocks), overlap
    return best


def iter_refseq_transcripts(lines: Iterable[str], build: str) -> Iterator[dict]:
    """Yield transcript records from an NCBI RefSeq GFF3, one gene at a time."""
    records: dict[str, dict] = {}
    for line in lines:
        if not line or line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        seqid, _, feature, start, end, _, strand, _, raw = fields
        if not seqid.startswith("NC_"):
            continue
        if feature in ("gene", "pseudogene"):
            yield from _flush_refseq(records, build)
            records = {}
            continue
        if feature in ("exon", "CDS"):
            attrs = _attributes(raw)
            parent = attrs.get("Parent", "")
            record = records.get(parent)
            if record is None:
                continue
            bounds = (int(start) - 1, int(end))
            if feature == "exon":
                record["exons"].append(bounds)
                note = attrs.get("Note", "")
                if note.startswith("The RefSeq transcript") and not record["note"]:
                    record["note"] = unquote(note)
            else:
                record["cds"].append(bounds)
                record["protein"] = attrs.get("protein_id", record["protein"])
                if "transl_except" in attrs or "exception" in attrs:
                    record["exception"] = True
            continue
        attrs = _attributes(raw)
        transcript_id = attrs.get("transcript_id", "")
        if not transcript_id.startswith(("NM_", "NR_")):
            continue
        if "trans-splicing" in attrs.get("exception", ""):
            continue
        dbxref = attrs.get("Dbxref", "")
        gene_id = re.search(r"GeneID:(\d+)", dbxref)
        hgnc = re.search(r"HGNC:(HGNC:\d+)", dbxref)
        note = attrs.get("Note", "")
        records[attrs.get("ID", transcript_id)] = {
            "accession": transcript_id,
            "gene": attrs.get("gene", ""),
            "gene_id": gene_id.group(1) if gene_id else "",
            "hgnc_id": hgnc.group(1) if hgnc else "",
            "seqid": seqid,
            "strand": 1 if strand == "+" else -1,
            "biotype": feature,
            "exons": [],
            "cds": [],
            "protein": "",
            "source": "RefSeq",
            "note": unquote(note) if note.startswith("The RefSeq transcript") else "",
            "exception": False,
            "tags": unquote(attrs.get("tag", "")),
        }
    yield from _flush_refseq(records, build)


def _flush_refseq(records: dict[str, dict], build: str) -> Iterator[dict]:
    for record in records.values():
        if not record["exons"]:
            continue
        try:
            contig = canonical_contig(record["seqid"], build)
        except InputError:
            continue
        if is_mitochondrial(contig):
            continue
        record["contig"] = contig
        yield record


def iter_gencode_transcripts(lines: Iterable[str], build: str) -> Iterator[dict]:
    """Yield transcript records from a GENCODE GFF3. Its CDS leaves out the stop codon, so add it."""
    records: dict[str, dict] = {}
    for line in lines:
        if not line or line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        seqid, _, feature, start, end, _, strand, _, raw = fields
        if feature == "gene":
            yield from _flush_gencode(records, build)
            records = {}
            continue
        if feature == "transcript":
            attrs = _attributes(raw)
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id.startswith("ENST"):
                continue
            hgnc = attrs.get("hgnc_id", "")
            records[attrs.get("ID", transcript_id)] = {
                "accession": transcript_id,
                "gene": attrs.get("gene_name", ""),
                "gene_id": "",
                "hgnc_id": hgnc,
                "seqid": seqid,
                "strand": 1 if strand == "+" else -1,
                "biotype": attrs.get("transcript_type", ""),
                "exons": [],
                "cds": [],
                "protein": attrs.get("protein_id", ""),
                "source": "GENCODE",
                "note": "",
                "exception": False,
                "tags": attrs.get("tag", ""),
            }
            continue
        if feature in ("exon", "CDS", "stop_codon"):
            attrs = _attributes(raw)
            record = records.get(attrs.get("Parent", ""))
            if record is None:
                continue
            bounds = (int(start) - 1, int(end))
            if feature == "exon":
                record["exons"].append(bounds)
            else:
                record["cds"].append(bounds)
                record["protein"] = attrs.get("protein_id", record["protein"])
    yield from _flush_gencode(records, build)


def _flush_gencode(records: dict[str, dict], build: str) -> Iterator[dict]:
    for record in records.values():
        if not record["exons"] or record["seqid"] == "chrM":
            continue
        try:
            contig = canonical_contig(record["seqid"], build)
        except InputError:
            continue
        record["contig"] = contig
        yield record


def _transcript_row(record: dict, mane: dict, alignments: dict | None = None) -> tuple:
    strand = record["strand"]
    alignment = _alignment_for(record, alignments)
    exons = sorted(alignment[0] if alignment else record["exons"], reverse=strand < 0)
    blocks = alignment[1] if alignment else []
    accession = record["accession"]
    base, _, version = accession.partition(".")
    low = min(s for s, _ in exons)
    high = max(e for _, e in exons)
    cds_start = cds_end = None
    protein = record["protein"]
    if record["cds"]:
        cds_start = min(s for s, _ in record["cds"])
        cds_end = max(e for _, e in record["cds"])
    else:
        protein = ""
    status, mane_accession = mane.get(base, ("", ""))
    tags = record["tags"]
    if not status:
        if "MANE Select" in tags or "MANE_Select" in tags:
            status, mane_accession = "MANE Select", accession
        elif "MANE Plus Clinical" in tags or "MANE_Plus_Clinical" in tags:
            status, mane_accession = "MANE Plus Clinical", accession
    return (
        accession,
        base,
        int(version or 0),
        record["gene"],
        record["gene_id"],
        record["hgnc_id"],
        record["contig"],
        strand,
        low,
        high,
        ucsc_bin(low, high),
        encode_exons(exons),
        cds_start,
        cds_end,
        record["biotype"],
        protein,
        protein.split(".")[0] if protein else None,
        record["source"],
        status,
        mane_accession,
        int("RefSeq Select" in tags or "MANE Select" in tags),
        record["note"],
        int(record["exception"]),
        encode_blocks(blocks),
    )


@closes_on_failure
def build_transcripts(sources: dict[str, Path], out: Path, build: str) -> dict:
    """sources holds refseq, mane and, optionally, gencode paths."""
    mane, mane_rows = read_mane(sources["mane"])
    con = _new_db(
        out,
        """
        CREATE TABLE transcripts(accession TEXT, base TEXT, version INT, gene TEXT, gene_id TEXT, hgnc_id TEXT,
            contig TEXT, strand INT, tx_low INT, tx_high INT, bin INT, exons TEXT, cds_start INT, cds_end INT,
            biotype TEXT, protein TEXT, protein_base TEXT, source TEXT, mane TEXT, mane_accession TEXT,
            refseq_select INT, genome_note TEXT, translation_exception INT, alignment TEXT);
        """,
    )
    insert = f"INSERT INTO transcripts({_TRANSCRIPT_COLUMNS}) VALUES({','.join('?' * 24)})"
    stats = Counter()
    with open_text(sources["refseq"]) as handle:
        alignments = read_refseq_alignments(handle)
    readers = [("refseq", iter_refseq_transcripts)]
    if "gencode" in sources:
        readers.append(("gencode", iter_gencode_transcripts))
    for name, reader in readers:
        batch = []
        with open_text(sources[name]) as handle:
            for record in reader(handle, build):
                row = _transcript_row(record, mane, alignments)
                batch.append(row)
                stats[name] += 1
                if row[-1]:
                    stats[f"{name}_with_alignment"] += 1
                if record["note"]:
                    stats[f"{name}_differs_from_genome"] += 1
                if record["cds"]:
                    stats[f"{name}_coding"] += 1
                if len(batch) >= _BATCH:
                    con.executemany(insert, batch)
                    batch.clear()
        con.executemany(insert, batch)
        if stats[name] == 0:
            con.close()
            raise SetupError(f"{sources[name]} produced no transcripts; the file format may have changed")
    _finish_db(
        con,
        """
        CREATE INDEX tx_location ON transcripts(contig, bin);
        CREATE INDEX tx_base ON transcripts(base);
        CREATE INDEX tx_gene ON transcripts(gene);
        CREATE INDEX tx_protein ON transcripts(protein_base);
        """,
    )
    stats["mane_entries"] = mane_rows
    return dict(stats)


# ClinVar
_CLINVAR_INFO = {"RS", "CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO", "MC", "ORIGIN", "CLNHGVS", "AF_EXAC", "AF_ESP", "AF_TGP", "ONC", "SCI"}
_NAME_RE = re.compile(
    r"^(?P<tx>[A-Z]{2}_\d+\.\d+)(?:\((?P<gene>[^)]+)\))?:(?P<c>[cnm]\.[^ ]+)(?: \((?P<p>p\.[^)]+)\))?"
)


def _float(value: str | None) -> float | None:
    if not value or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


@closes_on_failure
def build_clinvar(vcf_gz: Path, summary: Path | None, out: Path, build: str, genome: Genome) -> dict:
    con = _new_db(
        out,
        """
        CREATE TABLE clinvar(contig TEXT, pos INT, ref TEXT, alt TEXT, variation_id INT, rsid TEXT,
            classification TEXT, review_status TEXT, condition TEXT, gene TEXT, molecular_consequence TEXT,
            origin TEXT, hgvs_g TEXT, af_exac REAL, af_esp REAL, af_tgp REAL, oncogenicity TEXT,
            somatic_impact TEXT, PRIMARY KEY(contig, pos, ref, alt)) WITHOUT ROWID;
        CREATE TABLE clinvar_names(variation_id INT, transcript TEXT, gene TEXT, hgvs_c TEXT, hgvs_p TEXT,
            protein_style TEXT);
        """,
    )
    insert = "INSERT OR IGNORE INTO clinvar VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    stats = Counter()
    batch = []
    contigs: dict[str, str] = {}
    with open_text(vcf_gz) as handle:
        for line in handle:
            if line[0] == "#":
                if line.startswith("##fileDate="):
                    stats["file_date"] = line.strip().split("=", 1)[1]
                continue
            stats["records"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                stats["malformed"] += 1
                continue
            chrom, pos, vid, ref, alt, _, _, info = fields[:8]
            if alt == "." or not alt:
                stats["no_alt"] += 1
                continue
            try:
                contig = contigs.get(chrom) or contigs.setdefault(chrom, canonical_contig(chrom, build))
                if not genome.has(contig):
                    stats["contig_not_installed"] += 1
                    continue
                at, ref, alt = _normalized_quickly(genome, contig, int(pos), ref, alt)
            except ReferenceMismatch:
                stats["ref_mismatch"] += 1
                continue
            except InputError:
                stats["not_normalisable"] += 1
                continue
            data = {}
            for item in info.split(";"):
                key, equals, value = item.partition("=")
                if equals and key in _CLINVAR_INFO:
                    data[key] = value
            rs = data.get("RS", "")
            batch.append(
                (
                    contig,
                    at,
                    ref,
                    alt,
                    int(vid),
                    rs.split("|")[0] if rs else None,
                    data.get("CLNSIG"),
                    data.get("CLNREVSTAT"),
                    data.get("CLNDN"),
                    data.get("GENEINFO", "").split(":")[0] or None,
                    _consequence(data.get("MC", "")),
                    data.get("ORIGIN"),
                    data.get("CLNHGVS"),
                    _float(data.get("AF_EXAC")),
                    _float(data.get("AF_ESP")),
                    _float(data.get("AF_TGP")),
                    data.get("ONC"),
                    data.get("SCI"),
                )
            )
            if len(batch) >= _BATCH:
                con.executemany(insert, batch)
                batch.clear()
    con.executemany(insert, batch)
    con.commit()
    stats["rows"] = con.execute("SELECT count(*) FROM clinvar").fetchone()[0]
    if stats["rows"] == 0:
        con.close()
        raise SetupError(f"{vcf_gz} produced no ClinVar rows; the file format may have changed")
    if summary is not None:
        _load_clinvar_names(con, summary, build, stats)
    _finish_db(
        con,
        """
        CREATE INDEX clinvar_variation ON clinvar(variation_id);
        CREATE INDEX clinvar_rsid ON clinvar(rsid);
        CREATE INDEX names_variation ON clinvar_names(variation_id);
        CREATE INDEX names_protein ON clinvar_names(gene, protein_style);
        """,
    )
    return dict(stats)


def _normalized_quickly(genome: Genome, contig: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """88% of ClinVar is one changed base, where normalizing comes down to comparing REF with the genome.

    Anything else, including a REF that doesn't match, goes to normalize(), which also raises the right error.
    """
    if len(ref) == 1 and len(alt) == 1 and ref != alt and ref in "ACGT" and alt in "ACGT":
        if genome.fetch(contig, pos - 1, pos) == ref:
            return pos, ref, alt
    variant = normalize(genome, contig, pos, ref, alt)
    return variant.pos, variant.ref, variant.alt


def _consequence(mc: str) -> str | None:
    terms = [item.split("|")[-1] for item in mc.split(",") if item]
    return ",".join(terms) or None


def _load_clinvar_names(con, summary: Path, build: str, stats: Counter) -> None:
    assembly = "GRCh37" if build == HG19 else "GRCh38"
    insert = "INSERT INTO clinvar_names VALUES(?,?,?,?,?,?)"
    batch = []
    with open_text(summary) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        _require_columns(summary, header, ["Name", "GeneSymbol", "Assembly", "VariationID"])
        name_at, gene_at, assembly_at, variation_at = (
            header.index(n) for n in ("Name", "GeneSymbol", "Assembly", "VariationID")
        )
        top = max(name_at, gene_at, assembly_at, variation_at)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= top or fields[assembly_at] != assembly:
                continue
            match = _NAME_RE.match(fields[name_at])
            if not match:
                continue
            gene = match.group("gene") or fields[gene_at].split(";")[0]
            protein = match.group("p")
            batch.append((int(fields[variation_at]), match.group("tx"), gene, match.group("c"), protein, protein_style(protein)))
            stats["named"] += 1
            if len(batch) >= _BATCH:
                con.executemany(insert, batch)
                batch.clear()
    con.executemany(insert, batch)


# genes
def _first(value: str) -> str | None:
    value = (value or "").split("|")[0].strip()
    return value or None


@closes_on_failure
def build_genes(sources: dict[str, Path], out: Path) -> dict:
    con = _new_db(
        out,
        """
        CREATE TABLE genes(hgnc_id TEXT PRIMARY KEY, symbol TEXT, name TEXT, locus_type TEXT, location TEXT,
            ncbi_gene TEXT, ensembl_gene TEXT, uniprot TEXT, omim TEXT, mane_select TEXT);
        CREATE TABLE aliases(alias TEXT, hgnc_id TEXT, kind TEXT);
        CREATE TABLE clingen_validity(symbol TEXT, hgnc_id TEXT, disease TEXT, mondo TEXT, moi TEXT,
            classification TEXT, report_url TEXT, classified TEXT);
        CREATE TABLE clingen_dosage(symbol TEXT, ncbi_gene TEXT, haplo_score TEXT, haplo_description TEXT,
            triplo_score TEXT, triplo_description TEXT, evaluated TEXT);
        CREATE TABLE constraint_metrics(symbol TEXT, transcript TEXT, pli REAL, oe_lof REAL, loeuf REAL, mis_z REAL);
        CREATE TABLE gene_diseases(symbol TEXT, ncbi_gene TEXT, association TEXT, disease_id TEXT);
        """,
    )
    stats = Counter()
    with open_text(sources["hgnc"]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(
            sources["hgnc"],
            reader.fieldnames or [],
            ["hgnc_id", "symbol", "name", "alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id", "uniprot_ids", "omim_id"],
        )
        for row in reader:
            if row.get("status") and row["status"] != "Approved":
                continue
            con.execute(
                "INSERT OR IGNORE INTO genes VALUES(?,?,?,?,?,?,?,?,?,?)",
                (
                    row["hgnc_id"],
                    row["symbol"],
                    row["name"],
                    row.get("locus_type"),
                    row.get("location"),
                    _first(row["entrez_id"]),
                    _first(row["ensembl_gene_id"]),
                    _first(row["uniprot_ids"]),
                    _first(row["omim_id"]),
                    row.get("mane_select"),
                ),
            )
            stats["genes"] += 1
            for column, kind in (("alias_symbol", "alias"), ("prev_symbol", "previous")):
                for alias in filter(None, (row.get(column) or "").split("|")):
                    con.execute("INSERT INTO aliases VALUES(?,?,?)", (alias.upper(), row["hgnc_id"], kind))
                    stats["aliases"] += 1
    if stats["genes"] == 0:
        con.close()
        raise SetupError("the HGNC file produced no genes; the format may have changed")

    if "clingen_validity" in sources:
        with open_text(sources["clingen_validity"]) as handle:
            for row in csv.reader(handle):
                if len(row) < 9 or row[0] in ("GENE SYMBOL",) or row[0].startswith("+"):
                    continue
                con.execute(
                    "INSERT INTO clingen_validity VALUES(?,?,?,?,?,?,?,?)",
                    (row[0], row[1], row[2], row[3], row[4], row[6], row[7], row[8][:10]),
                )
                stats["clingen_validity"] += 1

    if "clingen_dosage" in sources:
        with open_text(sources["clingen_dosage"]) as handle:
            header = None
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith("#Gene Symbol"):
                    header = line.lstrip("#").split("\t")
                    _require_columns(
                        sources["clingen_dosage"], header, ["Gene Symbol", "Haploinsufficiency Score", "Triplosensitivity Score"]
                    )
                    continue
                if line.startswith("#") or header is None:
                    continue
                row = dict(zip(header, line.split("\t")))
                con.execute(
                    "INSERT INTO clingen_dosage VALUES(?,?,?,?,?,?,?)",
                    (
                        row.get("Gene Symbol"),
                        row.get("Gene ID"),
                        row.get("Haploinsufficiency Score"),
                        row.get("Haploinsufficiency Description"),
                        row.get("Triplosensitivity Score"),
                        row.get("Triplosensitivity Description"),
                        row.get("Date Last Evaluated"),
                    ),
                )
                stats["clingen_dosage"] += 1

    if "gnomad_constraint" in sources:
        with open_text(sources["gnomad_constraint"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            _require_columns(sources["gnomad_constraint"], reader.fieldnames or [], ["gene", "transcript", "pLI", "oe_lof", "oe_lof_upper"])
            for row in reader:
                con.execute(
                    "INSERT INTO constraint_metrics VALUES(?,?,?,?,?,?)",
                    (row["gene"], row["transcript"], _float(row["pLI"]), _float(row["oe_lof"]), _float(row["oe_lof_upper"]), _float(row.get("mis_z"))),
                )
                stats["constraint"] += 1

    if "hpo_diseases" in sources:
        with open_text(sources["hpo_diseases"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            _require_columns(sources["hpo_diseases"], reader.fieldnames or [], ["ncbi_gene_id", "gene_symbol", "association_type", "disease_id"])
            for row in reader:
                con.execute(
                    "INSERT INTO gene_diseases VALUES(?,?,?,?)",
                    (row["gene_symbol"], row["ncbi_gene_id"].replace("NCBIGene:", ""), row["association_type"], row["disease_id"]),
                )
                stats["gene_diseases"] += 1

    _finish_db(
        con,
        """
        CREATE INDEX genes_symbol ON genes(symbol);
        CREATE INDEX genes_symbol_upper ON genes(upper(symbol));
        CREATE INDEX aliases_alias ON aliases(alias);
        CREATE INDEX validity_symbol ON clingen_validity(symbol);
        CREATE INDEX dosage_symbol ON clingen_dosage(symbol);
        CREATE INDEX constraint_symbol ON constraint_metrics(symbol);
        CREATE INDEX diseases_symbol ON gene_diseases(symbol);
        """,
    )
    return dict(stats)


# proteins
_FEATURES_KEPT = {"DOMAIN", "ACT_SITE", "BINDING", "MOTIF", "DNA_BIND", "ZN_FING", "SITE"}
_FT_RE = re.compile(r"^FT   (\S+)\s+(?:[A-Z0-9]+:)?[<>?]?(\d+)(?:\.\.[<>?]?(\d+))?")
_FT_NOTE_RE = re.compile(r'^FT\s+/note="([^"]*)"?')
_GN_RE = re.compile(r"Name=([^;{\s]+)")


def iter_uniprot(lines: Iterable[str]) -> Iterator[dict]:
    """Yield reviewed entries from a UniProt flat file."""
    entry = None
    feature = None
    sequence: list[str] | None = None
    for line in lines:
        tag = line[:2]
        if tag == "ID":
            entry = {"reviewed": "Reviewed;" in line, "accession": "", "gene": "", "pfam": [], "features": [], "sequence": ""}
            feature = None
            sequence = None
        elif entry is None:
            continue
        elif tag == "AC":
            if not entry["accession"]:
                entry["accession"] = line[5:].split(";")[0].strip()
        elif tag == "GN":
            if not entry["gene"]:
                match = _GN_RE.search(line)
                if match:
                    entry["gene"] = match.group(1)
        elif tag == "DR":
            if line.startswith("DR   Pfam;"):
                parts = [part.strip().rstrip(".") for part in line[5:].split(";")]
                if len(parts) >= 3:
                    entry["pfam"].append((parts[1], parts[2]))
        elif tag == "FT":
            match = _FT_RE.match(line)
            if match:
                feature = None
                if match.group(1) in _FEATURES_KEPT:
                    start = int(match.group(2))
                    feature = [match.group(1), start, int(match.group(3) or start), ""]
                    entry["features"].append(feature)
            elif feature is not None and not feature[3]:
                match = _FT_NOTE_RE.match(line)
                if match:
                    feature[3] = match.group(1)
        elif tag == "SQ":
            sequence = []
        elif tag == "//":
            if sequence is not None:
                entry["sequence"] = "".join(sequence).replace(" ", "")
            if entry["reviewed"] and entry["accession"]:
                yield entry
            entry = None
            sequence = None
        elif sequence is not None:
            sequence.append(line.strip())


@closes_on_failure
def build_proteins(uniprot: Path, pfam_lines: Iterable[str], out: Path) -> dict:
    con = _new_db(
        out,
        """
        CREATE TABLE uniprot(accession TEXT PRIMARY KEY, gene TEXT, reviewed INT, sequence TEXT);
        CREATE TABLE uniprot_pfam(accession TEXT, pfam TEXT, name TEXT);
        CREATE TABLE uniprot_features(accession TEXT, type TEXT, start INT, end INT, note TEXT);
        CREATE TABLE pfam_regions(accession TEXT, pfam TEXT, name TEXT, start INT, end INT);
        """,
    )
    stats = Counter()
    wanted: set[str] = set()
    with open_text(uniprot) as handle:
        for entry in iter_uniprot(handle):
            con.execute(
                "INSERT OR IGNORE INTO uniprot VALUES(?,?,?,?)",
                (entry["accession"], entry["gene"], int(entry["reviewed"]), entry["sequence"]),
            )
            con.executemany(
                "INSERT INTO uniprot_pfam VALUES(?,?,?)", [(entry["accession"], a, n) for a, n in entry["pfam"]]
            )
            con.executemany(
                "INSERT INTO uniprot_features VALUES(?,?,?,?,?)",
                [(entry["accession"], *feature) for feature in entry["features"]],
            )
            wanted.add(entry["accession"])
            stats["reviewed_proteins"] += 1
    if not wanted:
        con.close()
        raise SetupError("the UniProt file produced no reviewed entries; the format may have changed")
    con.executemany("INSERT INTO pfam_regions VALUES(?,?,?,?,?)", _pfam_regions(pfam_lines, wanted, stats))
    if not stats["pfam_regions"]:
        raise SetupError("the Pfam file matched none of the UniProt proteins; the format may have changed")
    _finish_db(
        con,
        """
        CREATE INDEX uniprot_gene ON uniprot(gene);
        CREATE INDEX pfam_acc ON uniprot_pfam(accession);
        CREATE INDEX feature_acc ON uniprot_features(accession);
        CREATE INDEX region_acc ON pfam_regions(accession);
        """,
    )
    return dict(stats)


def _pfam_regions(lines: Iterable[str], wanted: set[str], stats: Counter):
    """Rows of a Pfam proteome file, whose column names sit in a comment: #<seq id> <alignment start> ..."""
    columns: dict[str, int] | None = None
    for line in lines:
        if line.startswith("#"):
            names = re.findall(r"<([^>]+)>", line)
            if names:
                _require_columns("Pfam regions", names, ["seq id", "alignment start", "alignment end", "hmm acc", "hmm name"])
                columns = {name: index for index, name in enumerate(names)}
            continue
        if columns is None:
            _require_columns("Pfam regions", [], ["seq id", "alignment start", "alignment end", "hmm acc", "hmm name"])
        fields = line.rstrip("\n").split("\t")
        if fields[columns["seq id"]] in wanted:
            stats["pfam_regions"] += 1
            yield (
                fields[columns["seq id"]],
                fields[columns["hmm acc"]],
                fields[columns["hmm name"]],
                int(fields[columns["alignment start"]]),
                int(fields[columns["alignment end"]]),
            )


# CIViC
@closes_on_failure
def build_civic(variants: Path, evidence: Path, out: Path) -> dict:
    con = _new_db(
        out,
        """
        CREATE TABLE civic_variants(variant_id INT, gene TEXT, variant TEXT, url TEXT, evidence_items INT,
            chromosome TEXT, start INT, stop INT, ref TEXT, alt TEXT, build TEXT, hgvs TEXT, clinvar_ids TEXT);
        """,
    )
    stats = Counter()
    counts: Counter = Counter()
    with open_text(evidence) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(evidence, reader.fieldnames or [], ["molecular_profile_id", "evidence_status"])
        for row in reader:
            if row["evidence_status"] == "accepted":
                counts[row["molecular_profile_id"]] += 1
                stats["accepted_evidence_items"] += 1
    with open_text(variants) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require_columns(variants, reader.fieldnames or [], ["variant_id", "gene", "variant", "variant_civic_url"])
        for row in reader:
            if not row["gene"]:
                continue
            con.execute(
                "INSERT INTO civic_variants VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    int(row["variant_id"]),
                    row["gene"],
                    row["variant"],
                    row["variant_civic_url"],
                    counts.get(row.get("single_variant_molecular_profile_id", ""), 0),
                    row.get("chromosome"),
                    int(row["start"]) if row.get("start", "").isdigit() else None,
                    int(row["stop"]) if row.get("stop", "").isdigit() else None,
                    row.get("reference_bases"),
                    row.get("variant_bases"),
                    row.get("reference_build"),
                    row.get("hgvs_descriptions"),
                    row.get("clinvar_ids"),
                ),
            )
            stats["variants"] += 1
    if stats["variants"] == 0:
        con.close()
        raise SetupError("the CIViC file produced no variants; the format may have changed")
    _finish_db(con, "CREATE INDEX civic_gene ON civic_variants(gene, variant);")
    return dict(stats)


# optional sources
def _int(value: str | None) -> int | None:
    try:
        return int(value) if value not in (None, ".") else None
    except ValueError:
        return None


def _require_columns(source, header: list[str], needed: list[str]) -> None:
    missing = [name for name in needed if name not in header]
    if missing:
        raise SetupError(
            f"{source}: expected columns are missing ({', '.join(missing)}). The source format has changed, "
            "so nothing was installed. Please report this."
        )
