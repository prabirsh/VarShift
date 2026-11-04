"""Genome builds: contig names, lengths and RefSeq accessions."""

from __future__ import annotations

import re

from .errors import InputError

HG19 = "hg19"
HG38 = "hg38"
BUILDS = (HG19, HG38)

_ALIASES = {
    "hg19": HG19,
    "grch37": HG19,
    "b37": HG19,
    "37": HG19,
    "hg38": HG38,
    "grch38": HG38,
    "b38": HG38,
    "38": HG38,
}

# Contig: (hg19 length, hg38 length, RefSeq accession suffix on hg19, on hg38)
_CONTIGS = {
    "chr1": (249250621, 248956422, "NC_000001.10", "NC_000001.11"),
    "chr2": (243199373, 242193529, "NC_000002.11", "NC_000002.12"),
    "chr3": (198022430, 198295559, "NC_000003.11", "NC_000003.12"),
    "chr4": (191154276, 190214555, "NC_000004.11", "NC_000004.12"),
    "chr5": (180915260, 181538259, "NC_000005.9", "NC_000005.10"),
    "chr6": (171115067, 170805979, "NC_000006.11", "NC_000006.12"),
    "chr7": (159138663, 159345973, "NC_000007.13", "NC_000007.14"),
    "chr8": (146364022, 145138636, "NC_000008.10", "NC_000008.11"),
    "chr9": (141213431, 138394717, "NC_000009.11", "NC_000009.12"),
    "chr10": (135534747, 133797422, "NC_000010.10", "NC_000010.11"),
    "chr11": (135006516, 135086622, "NC_000011.9", "NC_000011.10"),
    "chr12": (133851895, 133275309, "NC_000012.11", "NC_000012.12"),
    "chr13": (115169878, 114364328, "NC_000013.10", "NC_000013.11"),
    "chr14": (107349540, 107043718, "NC_000014.8", "NC_000014.9"),
    "chr15": (102531392, 101991189, "NC_000015.9", "NC_000015.10"),
    "chr16": (90354753, 90338345, "NC_000016.9", "NC_000016.10"),
    "chr17": (81195210, 83257441, "NC_000017.10", "NC_000017.11"),
    "chr18": (78077248, 80373285, "NC_000018.9", "NC_000018.10"),
    "chr19": (59128983, 58617616, "NC_000019.9", "NC_000019.10"),
    "chr20": (63025520, 64444167, "NC_000020.10", "NC_000020.11"),
    "chr21": (48129895, 46709983, "NC_000021.8", "NC_000021.9"),
    "chr22": (51304566, 50818468, "NC_000022.10", "NC_000022.11"),
    "chrX": (155270560, 156040895, "NC_000023.10", "NC_000023.11"),
    "chrY": (59373566, 57227415, "NC_000024.9", "NC_000024.10"),
}

# hg19 chrM is the old Yoruba sequence (16571 bp). Only rCRS matches ClinVar and RefSeq,
# so hg19 keeps two contigs: chrM (Yoruba) and chrMT (rCRS). On hg38 chrM already is rCRS.
RCRS_ACCESSION = "NC_012920.1"
RCRS_LENGTH = 16569
YORUBA_CHRM_LENGTH = 16571
HG19_RCRS = "chrMT"

_BY_ACCESSION: dict[str, tuple[str, str]] = {RCRS_ACCESSION: ("", "")}
for _name, (_l19, _l38, _ac19, _ac38) in _CONTIGS.items():
    _BY_ACCESSION[_ac19] = (HG19, _name)
    _BY_ACCESSION[_ac38] = (HG38, _name)

_ACCESSION_RE = re.compile(r"^NC_\d{6}\.\d+$")
_CONTIG_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|M|MT)$", re.IGNORECASE)


def parse_build(text: str) -> str:
    build = _ALIASES.get(text.strip().lower())
    if build is None:
        raise InputError(f"unknown genome build {text!r}; use hg19 or hg38")
    return build


def build_of_accession(accession: str) -> str | None:
    """hg19 or hg38 for a chromosome accession, None when it does not decide the build."""
    entry = _BY_ACCESSION.get(accession)
    return entry[0] or None if entry else None


def canonical_contig(name: str, build: str) -> str:
    """Convert any accepted spelling (7, chr7, NC_000007.13, MT) to the name used internally."""
    name = name.strip()
    if _ACCESSION_RE.match(name):
        if name == RCRS_ACCESSION:
            return HG19_RCRS if build == HG19 else "chrM"
        entry = _BY_ACCESSION.get(name)
        if entry is None:
            raise InputError(f"unknown chromosome accession {name}")
        if entry[0] != build:
            raise InputError(f"{name} belongs to {entry[0]}, but the build in use is {build}")
        return entry[1]
    match = _CONTIG_RE.match(name)
    if not match:
        return name  # alt or unplaced contig, passed through as written
    label = match.group(1).upper()
    if label == "MT":
        return HG19_RCRS if build == HG19 else "chrM"
    if label == "M":
        return "chrM"
    if label.isdigit():
        label = str(int(label))
        if not 1 <= int(label) <= 22:
            raise InputError(f"chromosome {name!r} does not exist in human")
    return f"chr{label}"


def contig_length(contig: str, build: str) -> int | None:
    if contig == HG19_RCRS or (contig == "chrM" and build == HG38):
        return RCRS_LENGTH
    if contig == "chrM":
        return YORUBA_CHRM_LENGTH
    entry = _CONTIGS.get(contig)
    if entry is None:
        return None
    return entry[0] if build == HG19 else entry[1]


def contig_accession(contig: str, build: str) -> str | None:
    if contig == HG19_RCRS or (contig == "chrM" and build == HG38):
        return RCRS_ACCESSION
    entry = _CONTIGS.get(contig)
    if entry is None:
        return None
    return entry[2] if build == HG19 else entry[3]


def primary_contigs() -> list[str]:
    return list(_CONTIGS)


def is_mitochondrial(contig: str) -> bool:
    return contig in ("chrM", HG19_RCRS)


def key_contig(contig: str) -> str:
    """Bare contig name used in outside URLs and ClinVar (7, X, MT)."""
    if is_mitochondrial(contig):
        return "MT"
    return contig[3:] if contig.startswith("chr") else contig


def detect_build_from_lengths(lengths: dict[str, int]) -> tuple[str | None, list[str]]:
    """Decide the build from VCF ##contig lengths.

    Returns the build (None when undecided) and the contigs that disagree with it, which
    catches a VCF whose header says one build and whose contigs say another.
    """
    votes = {HG19: 0, HG38: 0}
    for raw_name, length in lengths.items():
        contig = _normalise_vcf_contig(raw_name)
        entry = _CONTIGS.get(contig)
        if entry is None:
            continue
        if entry[0] == length:
            votes[HG19] += 1
        if entry[1] == length:
            votes[HG38] += 1
    if votes[HG19] == votes[HG38]:
        return None, []
    build = HG19 if votes[HG19] > votes[HG38] else HG38
    wrong = []
    for raw_name, length in lengths.items():
        contig = _normalise_vcf_contig(raw_name)
        entry = _CONTIGS.get(contig)
        if entry is not None and (entry[0] if build == HG19 else entry[1]) != length:
            wrong.append(raw_name)
    return build, wrong


def _normalise_vcf_contig(name: str) -> str:
    match = _CONTIG_RE.match(name)
    return f"chr{match.group(1).upper()}" if match else name
