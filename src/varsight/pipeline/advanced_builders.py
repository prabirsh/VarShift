"""Builders for the advanced databases.

Very large population files are read straight off the network and only the rows inside coding regions
are kept, so a 26 GB file becomes one or two GB on disk. Formats were taken from real samples of each
public file. COSMIC, OMIM, SpliceAI and PrimateAI need an account, so those readers follow the column
names in each provider's documentation and report clearly if a file doesn't match.
"""

from __future__ import annotations

import csv
import io
import json
import re
import zipfile
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

from ..builds import canonical_contig
from ..errors import InputError, ReferenceMismatch, ResolveError, SetupError
from ..genome import Genome
from ..hgvs import parse_dna_change
from ..resolver import _apply_forward
from ..store import decode_exons, open_readonly, ucsc_bin
from ..variant import normalize
from .builders import _finish_db, _float, _int, _new_db, _require_columns, closes_on_failure, open_text

_BATCH = 20_000
_KEY_COLUMNS = "contig TEXT, pos INT, ref TEXT, alt TEXT"
_KEY_PRIMARY = "PRIMARY KEY(contig, pos, ref, alt)"
_DBSNP_STUDIES = {
    "TOPMED": "af_topmed",
    "1000Genomes": "af_1000g",
    "GnomAD_genomes": "af_gnomad_genomes",
    "GnomAD_exomes": "af_gnomad_exomes",
    "ExAC": "af_exac",
    "dbGaP_PopFreq": "af_alfa",
}


class CodingRegions:
    """Exons of every installed transcript, widened a little, for keeping only useful rows."""

    def __init__(self, transcripts_db: Path, pad: int = 50):
        spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
        connection = open_readonly(transcripts_db)
        try:
            for row in connection.execute("SELECT contig, exons FROM transcripts"):
                for start, end in decode_exons(row["exons"]):
                    spans[row["contig"]].append((start + 1 - pad, end + pad))
        finally:
            connection.close()
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        for contig, items in spans.items():
            items.sort()
            starts, ends = [items[0][0]], [items[0][1]]
            for start, end in items[1:]:
                if start <= ends[-1] + 1:
                    ends[-1] = max(ends[-1], end)
                else:
                    starts.append(start)
                    ends.append(end)
            self._starts[contig], self._ends[contig] = starts, ends

    def holds(self, contig: str, pos: int) -> bool:
        starts = self._starts.get(contig)
        if not starts:
            return False
        index = bisect_right(starts, pos) - 1
        return index >= 0 and pos <= self._ends[contig][index]

    def total_bases(self) -> int:
        return sum(end - start + 1 for contig in self._starts for start, end in zip(self._starts[contig], self._ends[contig]))


class _Contigs:
    """Remember the conversion of each source contig name; large files repeat them millions of times."""

    def __init__(self, build: str):
        self._build = build
        self._known: dict[str, str | None] = {}

    def convert(self, name: str) -> str | None:
        if name in self._known:
            return self._known[name]
        try:
            contig = canonical_contig(name, self._build)
        except InputError:
            contig = None
        if contig is not None and not contig.startswith("chr"):
            contig = None
        self._known[name] = contig
        return contig


def _normalised(genome: Genome, contig: str, pos: int, ref: str, alt: str):
    """Left-aligned key. A plain SNV can't move, so it skips the genome read to keep large builds quick."""
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    variant = normalize(genome, contig, pos, ref, alt)
    return variant.pos, variant.ref, variant.alt


def _info(text: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def _allele_value(values: str | None, index: int) -> str | None:
    if values is None:
        return None
    parts = values.split(",")
    return parts[index] if index < len(parts) else None


def _vcf_rows(lines: Iterable[str], build: str, genome: Genome, regions: CodingRegions | None, stats: Counter):
    """Yield (contig, pos, ref, alt, allele index, info text, id, filter), keeping only wanted positions."""
    contigs = _Contigs(build)
    last_name = None
    contig = None
    for line in lines:
        if line[0] == "#":
            continue
        name, pos_text, rest = line.split("\t", 2)
        if name != last_name:
            last_name, contig = name, contigs.convert(name)
        if contig is None:
            stats["other_contigs"] += 1
            continue
        pos = int(pos_text)
        if regions is not None and not regions.holds(contig, pos):
            stats["outside_scope"] += 1
            continue
        fields = rest.rstrip("\r\n").split("\t")
        if len(fields) < 6:
            stats["malformed"] += 1
            continue
        identifier, ref, alts, _, filter_text, info = fields[:6]
        for index, alt in enumerate(alts.split(",")):
            if alt[:1] in ("<", "*", ".") or not alt:
                continue
            try:
                key_pos, key_ref, key_alt = _normalised(genome, contig, pos, ref, alt)
            except (InputError, ReferenceMismatch):
                stats["skipped"] += 1
                continue
            yield contig, key_pos, key_ref, key_alt, index, info, identifier, filter_text


@closes_on_failure
def build_dbsnp(lines: Iterable[str], out: Path, build: str, genome: Genome, regions: CodingRegions | None) -> dict:
    columns = ["rsid"] + list(_DBSNP_STUDIES.values())
    con = _new_db(out, f"CREATE TABLE variants({_KEY_COLUMNS}, rsid TEXT, {', '.join(c + ' REAL' for c in columns[1:])}, {_KEY_PRIMARY}) WITHOUT ROWID;")
    insert = f"INSERT OR IGNORE INTO variants VALUES({','.join('?' * (4 + len(columns)))})"
    stats: Counter = Counter()
    batch = []
    for contig, pos, ref, alt, index, info_text, identifier, _ in _vcf_rows(lines, build, genome, regions, stats):
        frequencies = dict.fromkeys(_DBSNP_STUDIES.values())
        match = re.search(r"(?:^|;)FREQ=([^;]+)", info_text)
        if match:
            for study in match.group(1).split("|"):
                label, _, values = study.partition(":")
                target = _DBSNP_STUDIES.get(label)
                if target:
                    # The list starts with the reference allele, so ALT number n is item n + 1.
                    frequencies[target] = _float(_allele_value(values, index + 1))
        batch.append((contig, pos, ref, alt, identifier.removeprefix("rs"), *frequencies.values()))
        stats["rows"] += 1
        if len(batch) >= _BATCH:
            con.executemany(insert, batch)
            batch.clear()
    con.executemany(insert, batch)
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX variants_rsid ON variants(rsid);")
    return dict(stats)


@closes_on_failure
def build_thousand_genomes(lines: Iterable[str], out: Path, build: str, genome: Genome, regions: CodingRegions | None) -> dict:
    fields = {"af": "AF", "af_afr": "AFR_AF", "af_amr": "AMR_AF", "af_eas": "EAS_AF", "af_eur": "EUR_AF", "af_sas": "SAS_AF"}
    con = _new_db(out, f"CREATE TABLE variants({_KEY_COLUMNS}, {', '.join(c + ' REAL' for c in fields)}, {_KEY_PRIMARY}) WITHOUT ROWID;")
    insert = f"INSERT OR IGNORE INTO variants VALUES({','.join('?' * (4 + len(fields)))})"
    stats: Counter = Counter()
    batch = []
    for contig, pos, ref, alt, index, info_text, _, _ in _vcf_rows(lines, build, genome, regions, stats):
        info = _info(info_text)
        batch.append((contig, pos, ref, alt, *(_float(_allele_value(info.get(tag), index)) for tag in fields.values())))
        stats["rows"] += 1
        if len(batch) >= _BATCH:
            con.executemany(insert, batch)
            batch.clear()
    con.executemany(insert, batch)
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_gnomad(lines_by_chrom: Iterable[tuple[str, Iterable[str]]], out: Path, build: str, genome: Genome, regions: CodingRegions | None) -> dict:
    con = _new_db(
        out,
        f"CREATE TABLE variants({_KEY_COLUMNS}, filter TEXT, ac INT, an INT, af REAL, nhomalt INT, popmax TEXT, af_popmax REAL, {_KEY_PRIMARY}) WITHOUT ROWID;",
    )
    insert = "INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?,?,?,?,?,?,?)"
    stats: Counter = Counter()
    for chrom, lines in lines_by_chrom:
        batch = []
        for contig, pos, ref, alt, index, info_text, _, filter_text in _vcf_rows(lines, build, genome, regions, stats):
            info = _info(info_text)
            batch.append(
                (
                    contig,
                    pos,
                    ref,
                    alt,
                    filter_text,
                    _int(_allele_value(info.get("AC"), index)),
                    _int(info.get("AN")),
                    _float(_allele_value(info.get("AF"), index)),
                    _int(_allele_value(info.get("nhomalt"), index)),
                    _allele_value(info.get("popmax"), index),
                    _float(_allele_value(info.get("AF_popmax"), index)),
                )
            )
            stats["rows"] += 1
            if len(batch) >= _BATCH:
                con.executemany(insert, batch)
                batch.clear()
        con.executemany(insert, batch)
        stats[f"chrom_{chrom}"] = 1
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_gnomad_sv(lines: Iterable[str], out: Path, build: str) -> dict:
    con = _new_db(out, "CREATE TABLE intervals(contig TEXT, start INT, end INT, bin INT, sv_id TEXT, sv_type TEXT, sv_length INT, af REAL, filter TEXT);")
    insert = "INSERT INTO intervals VALUES(?,?,?,?,?,?,?,?,?)"
    contigs = _Contigs(build)
    stats: Counter = Counter()
    batch = []
    for line in lines:
        if line[0] == "#":
            continue
        fields = line.rstrip("\r\n").split("\t", 8)
        if len(fields) < 8:
            stats["malformed"] += 1
            continue
        contig = contigs.convert(fields[0])
        if contig is None:
            continue
        info = _info(fields[7])
        start = int(fields[1])
        end = _int(info.get("END")) or start
        if end < start:
            start, end = end, start
        batch.append((contig, start, end, ucsc_bin(start - 1, end), fields[2], info.get("SVTYPE"), _int(info.get("SVLEN")), _float(_allele_value(info.get("AF"), 0)), fields[6]))
        stats["rows"] += 1
        if len(batch) >= _BATCH:
            con.executemany(insert, batch)
            batch.clear()
    con.executemany(insert, batch)
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX intervals_location ON intervals(contig, bin);")
    return dict(stats)


def _score_table(out: Path, extra: str) -> tuple:
    con = _new_db(out, f"CREATE TABLE variants({_KEY_COLUMNS}, {extra}, {_KEY_PRIMARY}) WITHOUT ROWID;")
    return con, Counter()


@closes_on_failure
def build_revel(zip_path: Path, out: Path, build: str) -> dict:
    position_column = "hg19_pos" if build == "hg19" else "grch38_pos"
    con, stats = _score_table(out, "score REAL")
    with zipfile.ZipFile(zip_path) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as raw:
            reader = csv.DictReader(io.TextIOWrapper(raw, encoding="ascii", newline=""))
            _require_columns(zip_path, reader.fieldnames or [], ["chr", position_column, "ref", "alt", "REVEL"])
            batch = []
            for row in reader:
                if row[position_column] in (".", ""):
                    continue
                batch.append((f"chr{row['chr']}", int(row[position_column]), row["ref"], row["alt"], _float(row["REVEL"])))
                stats["rows"] += 1
                if len(batch) >= _BATCH:
                    con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?)", batch)
                    batch.clear()
            con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?)", batch)
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_alphamissense(tsv_gz: Path, out: Path, build: str) -> dict:
    con, stats = _score_table(out, "am_pathogenicity REAL, am_class TEXT, uniprot TEXT, protein_variant TEXT")
    contigs = _Contigs(build)
    batch = []
    with open_text(tsv_gz) as handle:
        for line in handle:
            if line[0] == "#":
                continue
            fields = line.rstrip("\r\n").split("\t")
            contig = contigs.convert(fields[0]) if len(fields) >= 10 else None
            if contig is None:
                continue
            batch.append((contig, int(fields[1]), fields[2], fields[3], float(fields[8]), fields[9], fields[5], fields[7]))
            stats["rows"] += 1
            if len(batch) >= _BATCH:
                con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?,?,?,?)", batch)
                batch.clear()
    con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?,?,?,?)", batch)
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_spliceai(lines: Iterable[str], out: Path, build: str, genome: Genome, regions: CodingRegions | None) -> dict:
    con, stats = _score_table(out, "max_delta REAL, gene TEXT, details TEXT")
    batch = []
    for contig, pos, ref, alt, _, info_text, _, _ in _vcf_rows(lines, build, genome, regions, stats):
        match = re.search(r"(?:^|;)SpliceAI=([^;]+)", info_text)
        if not match:
            stats["no_score"] += 1
            continue
        best_score, best_gene, best_text = None, None, None
        for item in match.group(1).split(","):
            parts = item.split("|")
            if len(parts) < 6 or parts[0] != alt[-1:] and parts[0] != alt:
                continue
            scores = [_float(value) or 0.0 for value in parts[2:6]]
            if best_score is None or max(scores) > best_score:
                best_score, best_gene, best_text = max(scores), parts[1], item
        if best_score is None:
            continue
        batch.append((contig, pos, ref, alt, best_score, best_gene, best_text))
        stats["rows"] += 1
        if len(batch) >= _BATCH:
            con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?,?,?)", batch)
            batch.clear()
    con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?,?,?)", batch)
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_primateai(path: Path, out: Path, build: str) -> dict:
    con, stats = _score_table(out, "score REAL")
    contigs = _Contigs(build)
    columns: dict[str, int] | None = None
    batch = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if columns is None:
                _require_columns(path, fields, ["chr", "pos", "ref", "alt", "primateDL_score"])
                columns = {name: fields.index(name) for name in ("chr", "pos", "ref", "alt", "primateDL_score")}
                continue
            contig = contigs.convert(fields[columns["chr"]])
            if contig is None:
                continue
            batch.append((contig, int(fields[columns["pos"]]), fields[columns["ref"]], fields[columns["alt"]], _float(fields[columns["primateDL_score"]])))
            stats["rows"] += 1
            if len(batch) >= _BATCH:
                con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?)", batch)
                batch.clear()
    con.executemany("INSERT OR IGNORE INTO variants VALUES(?,?,?,?,?)", batch)
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_hotspots(json_path: Path, out: Path) -> dict:
    with open_text(json_path) as handle:
        try:
            entries = json.load(handle)
        except ValueError as exc:
            raise SetupError(f"{json_path} is not the expected JSON list; the source may have changed") from exc
    con = _new_db(out, "CREATE TABLE hotspots(gene TEXT, residue TEXT, start INT, end INT, tumor_count INT, type TEXT, variants TEXT, q_value TEXT);")
    stats: Counter = Counter()
    for entry in entries:
        position = entry.get("aminoAcidPosition") or {}
        if "hugoSymbol" not in entry or "start" not in position:
            raise SetupError(f"{json_path}: expected fields are missing; the source format has changed")
        variants = entry.get("variantAminoAcid") or {}
        con.execute(
            "INSERT INTO hotspots VALUES(?,?,?,?,?,?,?,?)",
            (
                entry["hugoSymbol"],
                entry.get("residue"),
                position["start"],
                position.get("end", position["start"]),
                entry.get("tumorCount"),
                entry.get("type"),
                ",".join(f"{aa}:{count}" for aa, count in sorted(variants.items())),
                str(entry.get("qValue")),
            ),
        )
        stats["rows"] += 1
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX hotspots_gene ON hotspots(gene, start);")
    return dict(stats)


@closes_on_failure
def build_decipher(hi_bed: Path, ddg2p: Path, out: Path) -> dict:
    con = _new_db(
        out,
        "CREATE TABLE genes(symbol TEXT PRIMARY KEY, hi_score REAL, hi_percent TEXT);"
        "CREATE TABLE ddg2p(symbol TEXT, disease TEXT, inheritance TEXT, confidence TEXT, consequence TEXT);",
    )
    stats: Counter = Counter()
    with open_text(hi_bed) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4 or "|" not in fields[3]:
                continue
            symbol, score, percent = (fields[3].split("|") + ["", ""])[:3]
            con.execute("INSERT OR IGNORE INTO genes VALUES(?,?,?)", (symbol, _float(score), percent))
            stats["hi_genes"] += 1
    if stats["hi_genes"] == 0:
        con.close()
        raise SetupError(f"{hi_bed} produced no genes; the format may have changed")
    with open_text(ddg2p) as handle:
        reader = csv.DictReader(handle)
        _require_columns(ddg2p, reader.fieldnames or [], ["gene symbol", "disease name", "allelic requirement", "confidence"])
        for row in reader:
            con.execute(
                "INSERT INTO ddg2p VALUES(?,?,?,?,?)",
                (row["gene symbol"], row["disease name"], row["allelic requirement"], row["confidence"], row.get("variant consequence")),
            )
            stats["ddg2p"] += 1
    _finish_db(con, "CREATE INDEX ddg2p_symbol ON ddg2p(symbol);")
    return dict(stats)


@closes_on_failure
def build_clingen_variants(tsv: Path, out: Path) -> dict:
    con = _new_db(
        out,
        "CREATE TABLE classifications(variation_id INT, allele_id TEXT, gene TEXT, disease TEXT, inheritance TEXT, assertion TEXT, criteria TEXT, expert_panel TEXT, url TEXT);",
    )
    stats: Counter = Counter()
    with open_text(tsv) as handle:
        header = handle.readline().lstrip("#").rstrip("\r\n").split("\t")
        _require_columns(tsv, header, ["ClinVar Variation Id", "HGNC Gene Symbol", "Assertion"])
        at = {name: index for index, name in enumerate(header)}

        def field(values: list[str], name: str) -> str | None:
            index = at.get(name)
            return values[index] if index is not None and index < len(values) else None

        for line in handle:
            values = line.rstrip("\r\n").split("\t")
            identifier = field(values, "ClinVar Variation Id")
            if not identifier or not identifier.isdigit():
                continue
            con.execute(
                "INSERT INTO classifications VALUES(?,?,?,?,?,?,?,?,?)",
                (
                    int(identifier),
                    field(values, "Allele Registry Id"),
                    field(values, "HGNC Gene Symbol"),
                    field(values, "Disease"),
                    field(values, "Mode of Inheritance"),
                    field(values, "Assertion"),
                    field(values, "Applied Evidence Codes (Met)"),
                    field(values, "Expert Panel"),
                    field(values, "Evidence Repo Link"),
                ),
            )
            stats["rows"] += 1
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX classifications_id ON classifications(variation_id);")
    return dict(stats)


def _table_reader(handle):
    """A dict reader for COSMIC's exports, which have been tab-separated and comma-separated over time."""
    first = handle.readline()
    delimiter = "\t" if "\t" in first else ","
    names = [re.sub(r"[^A-Z0-9]+", "_", name.upper()).strip("_") for name in next(csv.reader([first], delimiter=delimiter))]
    return names, csv.reader(handle, delimiter=delimiter)


@closes_on_failure
def build_cosmic(path: Path, out: Path, build: str, genome: Genome) -> dict:
    con, stats = _score_table(out, "cosmic_id TEXT, gene TEXT, cds TEXT, aa TEXT, samples INT")
    insert = (
        "INSERT INTO variants VALUES(?,?,?,?,?,?,?,?,1) "
        "ON CONFLICT(contig, pos, ref, alt) DO UPDATE SET samples = samples + 1"
    )
    contigs = _Contigs(build)
    with open_text(path) as handle:
        names, reader = _table_reader(handle)
        _require_columns(path, names, ["GENOMIC_MUTATION_ID", "HGVSG"])
        at = {name: index for index, name in enumerate(names)}
        for row in reader:
            if len(row) <= at["HGVSG"] or not row[at["HGVSG"]]:
                stats["no_position"] += 1
                continue
            chrom, _, change = row[at["HGVSG"]].partition(":g.")
            contig = contigs.convert(chrom)
            if contig is None or not change:
                stats["no_position"] += 1
                continue
            try:
                dna = parse_dna_change(change)
                variant = _apply_forward(genome, contig, dna.first.base - 1, dna.last.base, dna)
            except (InputError, ReferenceMismatch, ResolveError):
                stats["skipped"] += 1
                continue
            con.execute(
                insert,
                (
                    variant.contig,
                    variant.pos,
                    variant.ref,
                    variant.alt,
                    row[at["GENOMIC_MUTATION_ID"]],
                    row[at["GENE_SYMBOL"]] if "GENE_SYMBOL" in at else None,
                    row[at["MUTATION_CDS"]] if "MUTATION_CDS" in at else None,
                    row[at["MUTATION_AA"]] if "MUTATION_AA" in at else None,
                ),
            )
            stats["samples_read"] += 1
    stats["rows"] = con.execute("SELECT count(*) FROM variants").fetchone()[0]
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_cosmic_census(path: Path, out: Path) -> dict:
    con = _new_db(out, "CREATE TABLE genes(symbol TEXT PRIMARY KEY, tier TEXT, role TEXT, tumour_types TEXT);")
    stats: Counter = Counter()
    with open_text(path) as handle:
        names, reader = _table_reader(handle)
        symbol_at = _first_of(names, path, "GENE_SYMBOL")
        tier_at = _first_of(names, path, "TIER")
        role_at = _first_of(names, path, "ROLE_IN_CANCER")
        tumour_at = next((names.index(n) for n in ("TUMOUR_TYPES_SOMATIC", "TUMOUR_TYPES") if n in names), None)
        for row in reader:
            if len(row) <= max(symbol_at, tier_at, role_at):
                continue
            tumours = row[tumour_at] if tumour_at is not None and tumour_at < len(row) else None
            con.execute("INSERT OR IGNORE INTO genes VALUES(?,?,?,?)", (row[symbol_at], row[tier_at], row[role_at], tumours))
            stats["rows"] += 1
    _check_rows(con, out, stats)
    _finish_db(con)
    return dict(stats)


@closes_on_failure
def build_cosmic_fusions(path: Path, out: Path) -> dict:
    con = _new_db(out, "CREATE TABLE fusions(symbol TEXT, partner TEXT, samples INT);")
    counts: Counter = Counter()
    with open_text(path) as handle:
        names, reader = _table_reader(handle)
        five_at = _first_of(names, path, "FIVE_PRIME_GENE_SYMBOL")
        three_at = _first_of(names, path, "THREE_PRIME_GENE_SYMBOL")
        for row in reader:
            if len(row) <= max(five_at, three_at):
                continue
            five, three = row[five_at].split("_")[0], row[three_at].split("_")[0]
            if five and three:
                counts[(five, three)] += 1
                counts[(three, five)] += 1
    con.executemany("INSERT INTO fusions VALUES(?,?,?)", [(a, b, n) for (a, b), n in counts.items()])
    stats = Counter(rows=len(counts))
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX fusions_symbol ON fusions(symbol);")
    return dict(stats)


@closes_on_failure
def build_omim(path: Path, out: Path) -> dict:
    con = _new_db(out, "CREATE TABLE genes(symbol TEXT, mim TEXT, phenotypes TEXT);")
    stats: Counter = Counter()
    columns: dict[str, int] | None = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("# Chromosome") or line.startswith("#Chromosome"):
                header = line.lstrip("# ").rstrip("\r\n").split("\t")
                symbol_name = next((n for n in ("Approved Gene Symbol", "Approved Symbol") if n in header), None)
                _require_columns(path, header, ["MIM Number", "Phenotypes"])
                if symbol_name is None:
                    raise SetupError(f"{path}: no approved gene symbol column; is this genemap2.txt?")
                columns = {"symbol": header.index(symbol_name), "mim": header.index("MIM Number"), "phenotypes": header.index("Phenotypes")}
                continue
            if line.startswith("#"):
                continue
            if columns is None:
                raise SetupError(f"{path}: no header line found; is this genemap2.txt?")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= max(columns.values()) or not fields[columns["symbol"]]:
                continue
            con.execute("INSERT INTO genes VALUES(?,?,?)", (fields[columns["symbol"]], fields[columns["mim"]], fields[columns["phenotypes"]]))
            stats["rows"] += 1
    _check_rows(con, out, stats)
    _finish_db(con, "CREATE INDEX genes_symbol ON genes(symbol);")
    return dict(stats)


def _first_of(names: list[str], path: Path, wanted: str) -> int:
    if wanted not in names:
        raise SetupError(f"{path}: expected column {wanted} is missing. Is this the right COSMIC file? Please report it if so")
    return names.index(wanted)


def _check_rows(con, out: Path, stats: Counter) -> None:
    if not stats.get("rows"):
        con.close()
        raise SetupError(f"building {out.name} produced no rows; the file format may have changed. Nothing was installed")
