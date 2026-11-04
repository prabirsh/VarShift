"""DNA helpers and the genetic code."""

from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "U": "Sec",
    "O": "Pyl", "*": "Ter", "X": "Xaa",
}  # fmt: skip
AA1 = {three: one for one, three in AA3.items()}

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def translate(cds: str, to_stop: bool = True) -> str:
    """Translate from the first base. With to_stop, the result ends at and includes the first *."""
    protein = []
    for i in range(0, len(cds) - 2, 3):
        aa = CODON_TABLE.get(cds[i : i + 3], "X")
        protein.append(aa)
        if to_stop and aa == "*":
            break
    return "".join(protein)


def aa_three(protein: str) -> str:
    return "".join(AA3[aa] for aa in protein)


def codons_for(amino_acid: str) -> list[str]:
    return [codon for codon, aa in CODON_TABLE.items() if aa == amino_acid]


def is_dna(seq: str) -> bool:
    return bool(seq) and not set(seq) - set("ACGTN")
