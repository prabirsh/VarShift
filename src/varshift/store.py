"""Read-only access to the local data folder.

Every database is opened read-only with parameterised queries. No statement here is built from
user text, and SQLite's extension loading and schema triggers stay off.
"""

from __future__ import annotations

import json
import sqlite3
from functools import lru_cache
from pathlib import Path

from .bigwig import BigWig
from .errors import DataError, InputError
from .genome import FastaFile, Genome
from .transcripts import Transcript

MANIFEST = "manifest.json"
DB_DIR = "db"
GENOME_DIR = "genome"
SCORE_DIR = "scores"
SCHEMA_VERSION = 1


def ucsc_bin(start: int, end: int) -> int:
    """Smallest UCSC bin holding a 0-based half-open range (Kent et al., 2002)."""
    start_bin, end_bin = start >> 17, max(end - 1, start) >> 17
    bin_offsets = (512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0)
    for offset in bin_offsets:
        if start_bin == end_bin:
            return offset + start_bin
        start_bin >>= 3
        end_bin >>= 3
    raise ValueError(f"range {start}-{end} is too large to bin")


def bins_for_range(start: int, end: int) -> list[int]:
    bins = []
    start_bin, end_bin = start >> 17, max(end - 1, start) >> 17
    for offset in (512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0):
        bins.extend(range(offset + start_bin, offset + end_bin + 1))
        start_bin >>= 3
        end_bin >>= 3
    return bins


def open_readonly(path: Path) -> sqlite3.Connection:
    uri = f"{Path(path).resolve().as_uri()}?mode=ro&immutable=1"
    try:
        con = sqlite3.connect(uri, uri=True)
        con.execute("PRAGMA query_only=ON")
        con.execute("PRAGMA trusted_schema=OFF")
        con.execute("PRAGMA cache_size=-16000")
        con.execute("PRAGMA mmap_size=268435456")
    except sqlite3.Error as exc:
        raise DataError(f"cannot open {path}: {exc}") from exc
    con.row_factory = sqlite3.Row
    return con


def encode_exons(exons) -> str:
    return ",".join(f"{start}-{end}" for start, end in exons)


def decode_exons(text: str) -> tuple[tuple[int, int], ...]:
    pairs = []
    for item in text.split(","):
        start, end = item.split("-")
        pairs.append((int(start), int(end)))
    return tuple(pairs)


def encode_blocks(blocks) -> str:
    return ",".join(f"{start}-{end}-{index}" for start, end, index in blocks)


def decode_blocks(text: str) -> tuple[tuple[int, int, int], ...]:
    return tuple(tuple(int(number) for number in item.split("-")) for item in text.split(",")) if text else ()


class DataStore:
    def __init__(self, data_dir: Path, overrides: dict[str, Path] | None = None):
        """overrides maps a database name to a file used in its place, for testing a new build."""
        self._overrides = overrides or {}
        self.data_dir = Path(data_dir).expanduser()
        path = self.data_dir / MANIFEST
        if not path.is_file():
            raise DataError(f"no VarShift data in {self.data_dir}; run: varshift setup data={self.data_dir}")
        with open(path, encoding="utf-8") as handle:
            self.manifest = json.load(handle)
        if self.manifest.get("schema") != SCHEMA_VERSION:
            raise DataError(
                f"data in {self.data_dir} uses layout {self.manifest.get('schema')}, "
                f"but this version expects {SCHEMA_VERSION}; run: varshift setup"
            )
        self._connections: dict[str, sqlite3.Connection | None] = {}
        self._genomes: dict[str, Genome] = {}
        self._rcrs: FastaFile | None = None
        self._score_files: dict[str, BigWig] = {}
        self.sql_log: list | None = None
        self._transcript_cache: dict[tuple[str, str, str, int], Transcript] = {}
        self.gene_symbol = lru_cache(maxsize=4096)(self._gene_symbol)

    # plumbing
    @property
    def builds(self) -> list[str]:
        return list(self.manifest.get("builds", []))

    def versions(self) -> dict[str, str]:
        out = {}
        for name, entry in self.manifest.get("files", {}).items():
            source = entry.get("source") or {}
            out[name] = source.get("version") or entry.get("built_at", "")[:10]
        return out

    def _db(self, name: str) -> sqlite3.Connection | None:
        if name not in self._connections:
            path = self._overrides.get(name) or self.data_dir / DB_DIR / f"{name}.sqlite"
            self._connections[name] = open_readonly(path) if Path(path).is_file() else None
            if self.sql_log is not None and self._connections[name] is not None:
                # Used by "varshift explain" to show which query ran against which file.
                self._connections[name].set_trace_callback(lambda statement, db=name: self.sql_log.append((db, statement)))
        return self._connections[name]

    def has(self, name: str) -> bool:
        return self._db(name) is not None

    def has_genome(self, build: str) -> bool:
        return build in self._genomes or (self.data_dir / GENOME_DIR / f"{build}.fa").is_file()

    def genome(self, build: str) -> Genome:
        genome = self._genomes.get(build)
        if genome is None:
            fasta = self.data_dir / GENOME_DIR / f"{build}.fa"
            if not fasta.is_file():
                raise DataError(f"the {build} reference is not installed; run: varshift setup build={build}")
            rcrs = None
            if build == "hg19":
                hg38 = self.data_dir / GENOME_DIR / "hg38.fa"
                if hg38.is_file():
                    rcrs = self._rcrs or FastaFile(hg38)
                    self._rcrs = rcrs
            genome = Genome(build, fasta, rcrs)
            self._genomes[build] = genome
        return genome

    def close(self) -> None:
        for connection in self._connections.values():
            if connection is not None:
                connection.close()
        self._connections.clear()
        for genome in self._genomes.values():
            genome.close()
        self._genomes.clear()
        if self._rcrs is not None:
            self._rcrs.close()
            self._rcrs = None
        for reader in self._score_files.values():
            reader.close()
        self._score_files.clear()

    # transcripts
    def _make_transcript(self, build: str, row: sqlite3.Row) -> Transcript:
        # A few RefSeq transcripts are placed twice on one chromosome, so the locus is part of the key.
        cache_key = (build, row["accession"], row["contig"], row["tx_low"])
        cached = self._transcript_cache.get(cache_key)
        if cached is not None:
            return cached
        tx = Transcript(
            accession=row["accession"],
            gene=row["gene"],
            contig=row["contig"],
            strand=row["strand"],
            exons=decode_exons(row["exons"]),
            cds_start=row["cds_start"],
            cds_end=row["cds_end"],
            biotype=row["biotype"],
            protein=row["protein"] or "",
            source=row["source"],
            gene_id=row["gene_id"] or "",
            hgnc_id=row["hgnc_id"] or "",
            mane=row["mane"] or "",
            mane_accession=row["mane_accession"] or "",
            refseq_select=bool(row["refseq_select"]),
            genome_note=row["genome_note"] or "",
            translation_exception=bool(row["translation_exception"]),
            # A database built before this column existed still reads, with the old numbering.
            blocks=decode_blocks(row["alignment"]) if "alignment" in row.keys() else (),
        )
        if len(self._transcript_cache) > 4096:
            self._transcript_cache.clear()
        self._transcript_cache[cache_key] = tx
        return tx

    def _transcripts_db(self, build: str) -> sqlite3.Connection:
        db = self._db(f"transcripts.{build}")
        if db is None:
            raise DataError(f"transcripts for {build} are not installed; run: varshift setup build={build}")
        return db

    def transcripts_at(
        self, build: str, contig: str, start: int, end: int, flank: int = 0, source: str | None = "RefSeq"
    ) -> list[Transcript]:
        """Transcripts near a span. RefSeq first; GENCODE only where RefSeq has nothing here."""
        db = self._transcripts_db(build)
        low, high = max(start - flank, 0), end + flank
        bins = bins_for_range(low, max(high, low + 1))
        rows = []
        for chunk_start in range(0, len(bins), 400):
            chunk = bins[chunk_start : chunk_start + 400]
            marks = ",".join("?" * len(chunk))
            rows.extend(
                db.execute(
                    f"SELECT * FROM transcripts WHERE contig=? AND bin IN ({marks}) AND tx_low<? AND tx_high>?",
                    [contig, *chunk, max(high, low + 1), low],
                ).fetchall()
            )
        if source is not None:
            preferred = [row for row in rows if row["source"] == source]
            rows = preferred or rows
        return [self._make_transcript(build, row) for row in rows]

    def transcript_versions(self, build: str, accession: str) -> list[Transcript]:
        """Every placed version of an accession. A versioned accession matches that version first."""
        db = self._transcripts_db(build)
        base = accession.split(".")[0]
        rows = db.execute(
            "SELECT * FROM transcripts WHERE base=? ORDER BY version DESC, contig", (base,)
        ).fetchall()
        return [self._make_transcript(build, row) for row in rows]

    def transcripts_for_protein(self, build: str, protein: str) -> list[Transcript]:
        db = self._transcripts_db(build)
        base = protein.split(".")[0]
        rows = db.execute(
            "SELECT * FROM transcripts WHERE protein_base=? ORDER BY version DESC, contig", (base,)
        ).fetchall()
        return [self._make_transcript(build, row) for row in rows]

    def transcripts_for_gene(self, build: str, gene: str) -> list[Transcript]:
        db = self._transcripts_db(build)
        rows = db.execute("SELECT * FROM transcripts WHERE gene=?", (gene,)).fetchall()
        return sorted((self._make_transcript(build, row) for row in rows), key=Transcript.default_rank)

    # genes
    def _gene_symbol(self, text: str) -> str | None:
        """Resolve a symbol, alias or previous symbol to the approved symbol. None if unknown."""
        db = self._db("genes")
        if db is None:
            return None
        upper = text.upper()
        row = db.execute("SELECT symbol FROM genes WHERE upper(symbol)=?", (upper,)).fetchone()
        if row:
            return row["symbol"]
        rows = db.execute(
            "SELECT DISTINCT g.symbol, a.kind FROM aliases a JOIN genes g ON g.hgnc_id=a.hgnc_id WHERE a.alias=?",
            (upper,),
        ).fetchall()
        symbols = sorted({r["symbol"] for r in rows if r["kind"] == "previous"}) or sorted({r["symbol"] for r in rows})
        if not symbols:
            return None
        if len(symbols) > 1:
            raise InputError(f"{text} is an alias of several genes ({', '.join(symbols)}); use the approved symbol")
        return symbols[0]

    def gene_info(self, symbol: str) -> sqlite3.Row | None:
        db = self._db("genes")
        return db.execute("SELECT * FROM genes WHERE symbol=?", (symbol,)).fetchone() if db else None

    def gene_rows(self, table: str, symbol: str) -> list[sqlite3.Row]:
        if table not in ("clingen_validity", "clingen_dosage", "constraint_metrics", "gene_diseases"):
            raise ValueError(f"unknown gene table {table}")
        db = self._db("genes")
        return db.execute(f"SELECT * FROM {table} WHERE symbol=?", (symbol,)).fetchall() if db else []

    # ClinVar
    def clinvar_row(self, build: str, contig: str, pos: int, ref: str, alt: str) -> sqlite3.Row | None:
        db = self._db(f"clinvar.{build}")
        if db is None:
            return None
        return db.execute(
            "SELECT * FROM clinvar WHERE contig=? AND pos=? AND ref=? AND alt=?", (contig, pos, ref, alt)
        ).fetchone()

    def clinvar_by_rsid(self, build: str, rsid: str) -> list[sqlite3.Row]:
        db = self._db(f"clinvar.{build}")
        return db.execute("SELECT * FROM clinvar WHERE rsid=?", (rsid,)).fetchall() if db else []

    def clinvar_by_protein(self, build: str, gene: str, style: str) -> list[sqlite3.Row]:
        db = self._db(f"clinvar.{build}")
        if db is None:
            return []
        return db.execute(
            "SELECT c.* FROM clinvar_names n JOIN clinvar c ON c.variation_id=n.variation_id "
            "WHERE n.gene=? AND n.protein_style=?",
            (gene, style),
        ).fetchall()

    def clinvar_names(self, build: str, variation_id: int) -> list[sqlite3.Row]:
        db = self._db(f"clinvar.{build}")
        if db is None:
            return []
        return db.execute("SELECT * FROM clinvar_names WHERE variation_id=?", (variation_id,)).fetchall()

    # other evidence
    def uniprot_for_gene(self, symbol: str) -> sqlite3.Row | None:
        db = self._db("proteins")
        if db is None:
            return None
        return db.execute("SELECT * FROM uniprot WHERE gene=? ORDER BY reviewed DESC LIMIT 1", (symbol,)).fetchone()

    def uniprot_pfam(self, accession: str) -> list[sqlite3.Row]:
        db = self._db("proteins")
        return db.execute("SELECT * FROM uniprot_pfam WHERE accession=?", (accession,)).fetchall() if db else []

    def pfam_regions_at(self, accession: str, residue: int) -> list[sqlite3.Row]:
        db = self._db("proteins")
        if db is None:
            return []
        return db.execute(
            "SELECT pfam, name, start, end FROM pfam_regions WHERE accession=? AND start<=? AND end>=?",
            (accession, residue, residue),
        ).fetchall()

    def uniprot_features_at(self, accession: str, residue: int) -> list[sqlite3.Row]:
        db = self._db("proteins")
        if db is None:
            return []
        return db.execute(
            "SELECT * FROM uniprot_features WHERE accession=? AND start<=? AND end>=?",
            (accession, residue, residue),
        ).fetchall()

    def civic_variants(self, gene: str, name: str) -> list[sqlite3.Row]:
        db = self._db("civic")
        if db is None or not name:
            return []
        return db.execute("SELECT * FROM civic_variants WHERE gene=? AND upper(variant)=?", (gene, name.upper())).fetchall()

    # advanced databases

    def installed(self) -> list[str]:
        return list(self.manifest.get("files", {}))

    def advanced_variant(self, name: str, build: str, contig: str, pos: int, ref: str, alt: str) -> sqlite3.Row | None:
        db = self._db(f"{name}.{build}")
        if db is None:
            return None
        return db.execute("SELECT * FROM variants WHERE contig=? AND pos=? AND ref=? AND alt=?", (contig, pos, ref, alt)).fetchone()

    def dbsnp_by_rsid(self, build: str, rsid: str) -> list[sqlite3.Row]:
        db = self._db(f"dbsnp.{build}")
        return db.execute("SELECT * FROM variants WHERE rsid=?", (rsid,)).fetchall() if db else []

    def advanced_sv(self, build: str, contig: str, pos: int, limit: int = 5) -> list[sqlite3.Row]:
        db = self._db(f"gnomad_sv.{build}")
        if db is None:
            return []
        bins = bins_for_range(pos - 1, pos)
        marks = ",".join("?" * len(bins))
        return db.execute(
            f"SELECT * FROM intervals WHERE contig=? AND bin IN ({marks}) AND start<=? AND end>=? AND filter='PASS' ORDER BY af DESC",
            [contig, *bins, pos, pos],
        ).fetchmany(limit)

    def advanced_hotspots(self, gene: str, residue: int) -> list[sqlite3.Row]:
        """Every hotspot covering this residue, most tumours first. One residue can be in several."""
        db = self._db("hotspots")
        if db is None:
            return []
        return db.execute("SELECT * FROM hotspots WHERE gene=? AND start<=? AND end>=? ORDER BY tumor_count DESC", (gene, residue, residue)).fetchall()

    def advanced_clingen(self, variation_id: int) -> sqlite3.Row | None:
        db = self._db("clingen_variants")
        return db.execute("SELECT * FROM classifications WHERE variation_id=?", (variation_id,)).fetchone() if db else None

    def advanced_gene(self, name: str, table: str, symbol: str) -> list[sqlite3.Row]:
        if (name, table) not in (("decipher", "genes"), ("decipher", "ddg2p"), ("cosmic_census", "genes"), ("cosmic_fusions", "fusions"), ("omim", "genes")):
            raise ValueError(f"unknown advanced table {name}.{table}")
        db = self._db(name)
        return db.execute(f"SELECT * FROM {table} WHERE symbol=?", (symbol,)).fetchall() if db else []

    def phylop(self, build: str, contig: str, position: int) -> float | None:
        reader = self._score_files.get(build)
        if reader is None:
            path = self.data_dir / SCORE_DIR / f"phylop.{build}.bw"
            if not path.is_file():
                return None
            reader = BigWig(path)
            self._score_files[build] = reader
        return reader.value(contig, position)

    def liftover_db(self, source: str, target: str) -> sqlite3.Connection | None:
        return self._db(f"liftover.{source}_{target}")
