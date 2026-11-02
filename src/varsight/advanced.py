"""The advanced databases: what each one is, how to get it, and what it adds to the results.

Nothing here is installed by default. Each result and report gains these columns only for the databases
you chose to install. Sizes were read from each server on 1 November 2026. "Kept" sizes are estimates,
because the full downloads have not been measured on every computer.
"""

from __future__ import annotations

from dataclasses import dataclass

DOWNLOAD = "download"  # public file; VarSight fetches it
ACCEPT = "accept"  # public file with terms you must accept first
OWN_FILE = "own_file"  # needs your own account or licence, so you supply the file
ALIAS = "alias"  # comes with another database


@dataclass(frozen=True)
class Database:
    name: str
    title: str
    provides: str
    level: str  # variant, gene, residue, region, position, clinvar_id or alias
    access: str
    license: str
    transfer: str
    kept: str
    columns: tuple[tuple[str, str], ...] = ()
    get_it: str = ""
    alias_of: str = ""
    per_build: bool = True
    has_scope: bool = False


DATABASES: dict[str, Database] = {
    db.name: db
    for db in (
        Database(
            "dbsnp",
            "dbSNP",
            "rs numbers, plus the frequencies dbSNP carries from TOPMed, 1000 Genomes, gnomAD, ExAC and ALFA",
            "variant",
            DOWNLOAD,
            "NCBI places no restrictions on reuse",
            "26 GB per build, read without saving",
            "1 to 3 GB per build, coding regions",
            (
                ("dbsnp_rsid", "rsid"),
                ("topmed_af", "af_topmed"),
                ("g1k_af_in_dbsnp", "af_1000g"),
                ("gnomad_genomes_af_in_dbsnp", "af_gnomad_genomes"),
                ("gnomad_exomes_af_in_dbsnp", "af_gnomad_exomes"),
                ("exac_af_in_dbsnp", "af_exac"),
                ("alfa_af", "af_alfa"),
            ),
            has_scope=True,
        ),
        Database(
            "topmed",
            "TOPMed",
            "TOPMed allele frequencies. They come with dbSNP, so installing dbSNP gives you these",
            "alias",
            ALIAS,
            "As dbSNP. TOPMed's own per-sample files need dbGaP access",
            "Included in dbSNP",
            "Included in dbSNP",
            alias_of="dbsnp",
        ),
        Database(
            "thousand_genomes",
            "1000 Genomes Project",
            "Allele frequency overall and in five continental groups",
            "variant",
            DOWNLOAD,
            "Open data from the International Genome Sample Resource",
            "1.4 GB on hg19, 0.9 GB on hg38, read without saving",
            "Under 200 MB per build, coding regions",
            (
                ("g1k_af", "af"),
                ("g1k_afr_af", "af_afr"),
                ("g1k_amr_af", "af_amr"),
                ("g1k_eas_af", "af_eas"),
                ("g1k_eur_af", "af_eur"),
                ("g1k_sas_af", "af_sas"),
            ),
            has_scope=True,
        ),
        Database(
            "gnomad",
            "gnomAD exomes",
            "Allele frequency, popmax and homozygote count from gnomAD v2.1.1 exomes (hg38 uses gnomAD's own liftover)",
            "variant",
            DOWNLOAD,
            "Public release; check gnomad.broadinstitute.org for terms",
            "59 GB on hg19, 85 GB on hg38, read without saving. Limit it with chrom=",
            "2 to 4 GB per build (estimate)",
            (
                ("gnomad_af", "af"),
                ("gnomad_popmax", "popmax"),
                ("gnomad_popmax_af", "af_popmax"),
                ("gnomad_homozygotes", "nhomalt"),
                ("gnomad_filter", "filter"),
            ),
            has_scope=True,
        ),
        Database(
            "gnomad_sv",
            "gnomAD SV",
            "Structural variants from gnomAD that overlap the position",
            "region",
            DOWNLOAD,
            "Public release; check gnomad.broadinstitute.org for terms",
            "134 MB on hg19, 1.6 GB on hg38, read without saving",
            "Under 100 MB per build",
            (("gnomad_sv_overlapping", "count"), ("gnomad_sv_top", "top")),
        ),
        Database(
            "revel",
            "REVEL",
            "Missense pathogenicity score from 0 to 1, for every possible missense change",
            "variant",
            ACCEPT,
            "Free for non-commercial use, per the REVEL site; ask its authors about other uses",
            "About 640 MB",
            "About 3 GB per build (estimate)",
            (("revel", "score"),),
        ),
        Database(
            "alphamissense",
            "AlphaMissense",
            "Missense pathogenicity prediction from Google DeepMind",
            "variant",
            ACCEPT,
            "The file header states CC BY-NC-SA 4.0 (non-commercial); a later Zenodo record says CC BY 4.0",
            "About 600 MB per build",
            "About 3 GB per build (estimate)",
            (("alphamissense_score", "am_pathogenicity"), ("alphamissense_class", "am_class")),
        ),
        Database(
            "spliceai",
            "SpliceAI",
            "Predicted splicing change (highest delta score)",
            "variant",
            OWN_FILE,
            "Free for academic and not-for-profit use; other use needs a licence from Illumina",
            "Your own download",
            "2 to 5 GB per build (estimate)",
            (("spliceai_max", "max_delta"), ("spliceai_gene", "gene")),
            get_it="Sign in at basespace.illumina.com/s/otSPW8hnhaZR and download spliceai_scores.masked.snv.hg19.vcf.gz or the hg38 one",
            has_scope=True,
        ),
        Database(
            "primateai",
            "PrimateAI",
            "Missense pathogenicity score based on variation in primates",
            "variant",
            OWN_FILE,
            "Check Illumina's terms for your use",
            "Your own download",
            "About 3 GB per build (estimate)",
            (("primateai", "score"),),
            get_it="Sign in at basespace.illumina.com/s/yYGFdGih1rXL and download PrimateAI_scores_v0.2.tsv.gz or the hg38 one",
        ),
        Database(
            "phylop",
            "phyloP",
            "Base-level conservation across 100 vertebrates",
            "position",
            DOWNLOAD,
            "UCSC states genome data is free for public and commercial use",
            "8.9 GB on hg19, 9.2 GB on hg38",
            "8.9 GB on hg19, 9.2 GB on hg38",
            (("phylop100way", "value"),),
        ),
        Database(
            "hotspots",
            "Cancer Hotspots",
            "Whether the residue is a recurrent cancer hotspot, and in how many tumours",
            "residue",
            DOWNLOAD,
            "Public website; check cancerhotspots.org for terms",
            "Under 5 MB",
            "Under 5 MB",
            (("cancer_hotspot", "label"), ("cancer_hotspot_tumours", "tumor_count")),
            per_build=False,
        ),
        Database(
            "decipher",
            "DECIPHER",
            "Haploinsufficiency prediction for the gene, and DDG2P developmental-disorder entries. Patient-level DECIPHER data needs a data access agreement, so that stays a link",
            "gene",
            DOWNLOAD,
            "Public downloads; check deciphergenomics.org and Gene2Phenotype for terms",
            "Under 5 MB",
            "Under 10 MB",
            (("decipher_hi_percent", "hi_percent"), ("ddg2p", "ddg2p")),
            per_build=False,
        ),
        Database(
            "clingen_variants",
            "ClinGen variant classifications",
            "Expert panel classification and the evidence codes applied. ClinGen gene curation is already in the core data",
            "clinvar_id",
            DOWNLOAD,
            "CC0, per ClinGen",
            "About 20 MB",
            "About 30 MB",
            (("clingen_classification", "assertion"), ("clingen_evidence_codes", "criteria")),
            per_build=False,
        ),
        Database(
            "cosmic",
            "COSMIC",
            "COSMIC mutation ID and how many tumour samples carry this change",
            "variant",
            OWN_FILE,
            "Free for non-commercial use after registering; commercial users need a licence from QIAGEN",
            "Your own download",
            "About 1 GB per build (estimate)",
            (("cosmic_id", "cosmic_id"), ("cosmic_samples", "samples")),
            get_it="Register at cancer.sanger.ac.uk/cosmic/download, then download Cosmic_GenomeScreensMutant for your genome build",
        ),
        Database(
            "cosmic_census",
            "COSMIC Cancer Gene Census",
            "Whether the gene is a known cancer gene: its tier, role and tumour types",
            "gene",
            OWN_FILE,
            "Free for non-commercial use after registering; commercial users need a licence from QIAGEN",
            "Your own download",
            "Under 5 MB",
            (("cgc_tier", "tier"), ("cgc_role", "role"), ("cgc_tumour_types", "tumour_types")),
            get_it="Register at cancer.sanger.ac.uk/cosmic/download, then download Cosmic_CancerGeneCensus",
            per_build=False,
        ),
        Database(
            "cosmic_fusions",
            "COSMIC gene fusions",
            "Known fusion partners of the gene, with sample counts",
            "gene",
            OWN_FILE,
            "Free for non-commercial use after registering; commercial users need a licence from QIAGEN",
            "Your own download",
            "Under 20 MB",
            (("cosmic_fusion_partners", "partners"),),
            get_it="Register at cancer.sanger.ac.uk/cosmic/download, then download Cosmic_Fusion",
            per_build=False,
        ),
        Database(
            "omim",
            "OMIM",
            "Disease phenotypes and inheritance linked to the gene",
            "gene",
            OWN_FILE,
            "Copyright Johns Hopkins University. Free for research after registering; commercial use needs a licence",
            "Your own download",
            "Under 20 MB",
            (("omim_phenotypes", "phenotypes"),),
            get_it="Register at omim.org/downloads, then download genemap2.txt with your own key",
            per_build=False,
        ),
    )
}


def get_database(name: str) -> Database | None:
    return DATABASES.get(name.strip().lower().replace("-", "_"))


def artifact_names(database: Database, builds: list[str]) -> list[str]:
    """The installed-file names for this database on these builds."""
    if database.access == ALIAS:
        return []
    return [f"{database.name}.{build}" for build in builds] if database.per_build else [database.name]


def output_columns(installed: list[str]) -> list[str]:
    """Extra report columns, only for the databases installed here."""
    columns: list[str] = []
    for database in DATABASES.values():
        if database.access != ALIAS and any(a == database.name or a.startswith(f"{database.name}.") for a in installed):
            columns.extend(column for column, _ in database.columns)
    return columns


def flatten_advanced(evidence: dict) -> dict:
    """The report columns for one result's advanced evidence."""
    row: dict = {}
    for name, values in (evidence.get("advanced") or {}).items():
        database = DATABASES.get(name)
        if database is None or not values:
            continue
        for column, field in database.columns:
            row[column] = values.get(field)
    return row
