"""Effect of a variant on one transcript: Sequence Ontology terms and HGVS c. and p. names."""

from __future__ import annotations

from dataclasses import dataclass, field

from .builds import contig_accession
from .genome import Genome
from .sequence import AA3, reverse_complement, translate
from .transcripts import Transcript
from .variant import Variant

UPSTREAM_DISTANCE = 5000

# Most severe first; the same order Ensembl VEP publishes for these terms.
SEVERITY = [
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
]
_RANK = {term: rank for rank, term in enumerate(SEVERITY)}


@dataclass
class TranscriptEffect:
    transcript: Transcript
    terms: list[str] = field(default_factory=list)
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    protein_short: str | None = None
    exon: int | None = None
    intron: int | None = None
    aa_position: int | None = None
    flags: list[str] = field(default_factory=list)

    @property
    def effect(self) -> str | None:
        return self.terms[0] if self.terms else None


class SequenceCache:
    """Spliced transcript sequences, kept for the transcripts used most recently."""

    def __init__(self, genome: Genome, capacity: int = 128):
        self._genome = genome
        self._capacity = capacity
        self._items: dict[str, tuple[str, str]] = {}

    def get(self, transcript: Transcript) -> tuple[str, str]:
        """The spliced sequence and its translation, computed once for each transcript."""
        # The same accession can sit on chrX and chrY (pseudoautosomal genes), or twice on one chromosome,
        # so the locus is in the key.
        key = f"{transcript.accession}@{transcript.contig}:{transcript.tx_low}"
        item = self._items.pop(key, None)
        if item is None:
            item = _sequences_of(transcript, self._genome)
            if len(self._items) >= self._capacity:
                self._items.pop(next(iter(self._items)))
        self._items[key] = item
        return item


def _sequences_of(transcript: Transcript, genome: Genome) -> tuple[str, str]:
    seq = transcript.spliced(genome)
    first, stop = transcript.cds_tx_range
    return seq, translate(seq[first:stop])


def annotate(genome: Genome, variant: Variant, tx: Transcript, cache: SequenceCache | None = None) -> TranscriptEffect:
    start, end, inserted = variant.edit()
    result = TranscriptEffect(tx)
    if tx.genome_note:
        result.flags.append("transcript_differs_from_genome")
    if tx.mane and tx.mane_accession and tx.mane_accession != tx.accession:
        result.flags.append("mane_version_differs")

    if end <= tx.tx_low or start >= tx.tx_high:
        before = end <= tx.tx_low
        result.terms.append(
            "upstream_gene_variant" if before == (tx.strand > 0) else "downstream_gene_variant"
        )
        result.hgvs_c = _hgvs_c(genome, tx, start, end, inserted)
        return result

    # VCF keys are left-aligned, but HGVS and the effect both use the 3' position on the transcript.
    bounds = _shift_bounds(tx, start, end)
    if bounds is not None:
        start, end, inserted = _shift_three_prime(genome, tx.contig, start, end, inserted, tx.strand, *bounds)

    _classify_regions(tx, start, end, result)
    if bounds is not None:
        # Inside a repeat, every placement is the same DNA change. A splice site touched by the 5'-most
        # placement counts too, whichever strand the gene is on.
        five_start, five_end, _ = _shift_three_prime(genome, tx.contig, start, end, inserted, -tx.strand, *bounds)
        if (five_start, five_end) != (start, end):
            other = TranscriptEffect(tx)
            _classify_regions(tx, five_start, five_end, other)
            result.terms.extend(term for term in other.terms if term.startswith("splice_"))
    result.hgvs_c = _hgvs_c(genome, tx, start, end, inserted)
    if tx.coding and _touches_cds(tx, start, end):
        seq, ref_protein = cache.get(tx) if cache else _sequences_of(tx, genome)
        _protein_effect(seq, ref_protein, tx, start, end, inserted, result)
    elif not tx.coding:
        result.terms.append("non_coding_transcript_variant")
    if not result.terms:
        result.terms.append("coding_sequence_variant")
    result.terms = sorted(set(result.terms), key=_RANK.__getitem__)
    return result


def _touches_cds(tx: Transcript, start: int, end: int) -> bool:
    if start == end:
        return tx.cds_start < start < tx.cds_end
    return start < tx.cds_end and end > tx.cds_start


def _classify_regions(tx: Transcript, start: int, end: int, result: TranscriptEffect) -> None:
    point = start == end
    in_exon = False
    for number, (low, high) in enumerate(tx.exons):
        if (low < start < high) if point else (start < high and end > low):
            in_exon = True
            if result.exon is None:
                result.exon = number + 1
            if not tx.coding:
                result.terms.append("non_coding_transcript_exon_variant")
            else:
                _classify_utr(tx, low, high, start, end, point, result)
    for number, (low, high) in enumerate(tx.introns):
        # An insertion exactly on an exon-intron junction counts as touching the intron.
        hit = (low <= start <= high) if point else (start < high and end > low)
        if not hit:
            continue
        if result.intron is None:
            result.intron = number + 1
        if not (point and in_exon):
            result.terms.append("intron_variant")
        _classify_splice(tx, low, high, start, end, point, result)
    _classify_exon_edges(tx, start, end, point, result)


def _classify_utr(tx, low, high, start, end, point, result) -> None:
    low_side, high_side = "5_prime_UTR_variant", "3_prime_UTR_variant"
    if tx.strand < 0:
        low_side, high_side = high_side, low_side
    if point:
        if start <= tx.cds_start:
            result.terms.append(low_side)
        elif start >= tx.cds_end:
            result.terms.append(high_side)
        return
    if max(start, low) < tx.cds_start:
        result.terms.append(low_side)
    if min(end, high) > tx.cds_end:
        result.terms.append(high_side)


def _overlaps(start: int, end: int, point: bool, low: int, high: int) -> bool:
    if point:
        return low < start < high
    return start < high and end > low


def _classify_splice(tx, low, high, start, end, point, result) -> None:
    donor_end = (low, low + 2) if tx.strand > 0 else (high - 2, high)
    acceptor_end = (high - 2, high) if tx.strand > 0 else (low, low + 2)
    if _overlaps(start, end, point, *donor_end):
        result.terms.append("splice_donor_variant")
    if _overlaps(start, end, point, *acceptor_end):
        result.terms.append("splice_acceptor_variant")
    # Intronic part of the splice region: bases 3 to 8 from either end of the intron.
    if _overlaps(start, end, point, low + 2, low + 8) or _overlaps(start, end, point, high - 8, high - 2):
        result.terms.append("splice_region_variant")
    if point and start in (low, high):
        result.terms.append("splice_region_variant")


def _classify_exon_edges(tx, start, end, point, result) -> None:
    # Exonic part of the splice region: the three bases next to an intron.
    for low, high in tx.introns:
        if _overlaps(start, end, point, low - 3, low) or _overlaps(start, end, point, high, high + 3):
            result.terms.append("splice_region_variant")


# HGVS at DNA level
def _shift_three_prime(genome, contig, start, end, inserted, strand, low, high):
    """Move a pure insertion or deletion as far 3' as its repeat allows, within [low, high)."""
    if start != end and inserted:
        return start, end, inserted
    window = 512
    if strand > 0:
        block_start, block = start, genome.fetch(contig, start, min(start + window, high))
        if start == end:
            while start < high and inserted:
                index = start - block_start
                if index >= len(block):
                    block_start, block = start, genome.fetch(contig, start, min(start + window, high))
                    if not block:
                        break
                    index = 0
                if block[index] != inserted[0]:
                    break
                start, end = start + 1, end + 1
                inserted = inserted[1:] + inserted[0]
            return start, end, inserted
        while end < high:
            index_a, index_b = start - block_start, end - block_start
            if index_b >= len(block):
                block_start, block = start, genome.fetch(contig, start, min(start + (end - start) + window, high))
                index_a, index_b = 0, end - start
                if index_b >= len(block):
                    break
            if block[index_a] != block[index_b]:
                break
            start, end = start + 1, end + 1
        return start, end, inserted

    block_end, block = end, genome.fetch(contig, max(end - window, low), end)
    block_start = block_end - len(block)
    if start == end:
        while start > low and inserted:
            index = start - 1 - block_start
            if index < 0:
                block_end, block = start, genome.fetch(contig, max(start - window, low), start)
                block_start = block_end - len(block)
                if not block:
                    break
                index = start - 1 - block_start
            if block[index] != inserted[-1]:
                break
            start, end = start - 1, end - 1
            inserted = inserted[-1] + inserted[:-1]
        return start, end, inserted
    while start > low:
        index_a, index_b = start - 1 - block_start, end - 1 - block_start
        if index_a < 0:
            block_end, block = end, genome.fetch(contig, max(end - (end - start) - window, low), end)
            block_start = block_end - len(block)
            index_a, index_b = start - 1 - block_start, end - 1 - block_start
            if index_a < 0:
                break
        if block[index_a] != block[index_b]:
            break
        start, end = start - 1, end - 1
    return start, end, inserted


def _shift_bounds(tx: Transcript, start: int, end: int) -> tuple[int, int] | None:
    """How far the 3' rule may move a change. It can cross an exon-intron junction, as ClinVar does."""
    if end <= tx.tx_low:
        return max(0, tx.tx_low - 2 * UPSTREAM_DISTANCE), tx.tx_low
    if start >= tx.tx_high:
        return tx.tx_high, tx.tx_high + 2 * UPSTREAM_DISTANCE
    if tx.tx_low <= start and end <= tx.tx_high:
        return tx.tx_low, tx.tx_high
    return None


def _dup_span(genome, contig, start, inserted, strand, low, high):
    """If inserting here repeats the bases just 5' of the point, return that span."""
    size = len(inserted)
    if strand > 0:
        if start - size >= low and genome.fetch(contig, start - size, start) == inserted:
            return start - size, start
    elif start + size <= high and genome.fetch(contig, start, start + size) == inserted:
        return start, start + size
    return None


def _hgvs_c(genome: Genome, tx: Transcript, start: int, end: int, inserted: str) -> str | None:
    prefix = "c" if tx.coding else "n"
    change = _describe(genome, tx.contig, start, end, inserted, tx.strand, _shift_bounds(tx, start, end), tx.position_of)
    return f"{tx.accession}:{prefix}.{change}" if change else None


def hgvs_g(genome: Genome, variant: Variant) -> str | None:
    """Chromosome-level name, shifted 3' on the forward strand."""
    accession = contig_accession(variant.contig, variant.build)
    if accession is None:
        return None
    start, end, inserted = variant.edit()
    level = "m" if accession.startswith("NC_012920") else "g"
    bounds = (max(0, start - 100_000), end + 100_000)
    change = _describe(genome, variant.contig, start, end, inserted, 1, bounds, lambda g: g + 1)
    return f"{accession}:{level}.{change}"


def _describe(genome, contig, start, end, inserted, strand, bounds, position_of) -> str | None:
    """Name one change. position_of maps a genomic base to its label in the chosen numbering."""
    if bounds is not None:
        start, end, inserted = _shift_three_prime(genome, contig, start, end, inserted, strand, *bounds)
    low, high = bounds if bounds is not None else (start, end)
    strand_seq = (lambda s: s) if strand > 0 else reverse_complement

    def span(first_low: int, first_high: int) -> str | None:
        five = position_of(first_low if strand > 0 else first_high - 1)
        three = position_of(first_high - 1 if strand > 0 else first_low)
        if five is None or three is None:
            return None
        return str(five) if five == three else f"{five}_{three}"

    if start == end:
        dup = _dup_span(genome, contig, start, inserted, strand, low, high)
        if dup is not None:
            label = span(*dup)
            return f"{label}dup" if label else None
        left = position_of(start - 1 if strand > 0 else start)
        right = position_of(start if strand > 0 else start - 1)
        if left is None or right is None:
            return None
        return f"{left}_{right}ins{strand_seq(inserted)}"

    label = span(start, end)
    if label is None:
        return None
    if not inserted:
        return f"{label}del"
    deleted = genome.fetch(contig, start, end)
    if len(deleted) == 1 and len(inserted) == 1:
        return f"{label}{strand_seq(deleted)}>{strand_seq(inserted)}"
    if len(deleted) >= 2 and inserted == reverse_complement(deleted):
        return f"{label}inv"
    return f"{label}delins{strand_seq(inserted)}"


# Protein level
def _protein_effect(seq, ref_protein, tx, start, end, inserted, result) -> None:
    cds_first, cds_stop = tx.cds_tx_range
    span = _tx_span(tx, start, end)
    if span is None:
        if result.exon is not None:
            result.terms.append("coding_sequence_variant")
            one_exon = tx.locate(start) == tx.locate(max(start, end - 1))
            result.flags.append("transcript_gap_here" if tx.blocks and one_exon else "spans_exon_intron_boundary")
        return
    a, b = span
    change = inserted if tx.strand > 0 else reverse_complement(inserted)
    if b <= cds_first or a >= cds_stop:
        return  # entirely in UTR (the classifier already named it)

    cds_length = cds_stop - cds_first
    if cds_length % 3 or not ref_protein.endswith("*") or len(ref_protein) * 3 != cds_length:
        reason = "cds_internal_stop" if tx.translation_exception or "*" in ref_protein[:-1] else "cds_incomplete"
        result.flags.append(reason)
        result.terms.append("coding_sequence_variant")
        return

    a_cds, b_cds = a - cds_first, b - cds_first
    first_codon = max(a_cds, 0) // 3 * 3
    last_codon = (max(b_cds, a_cds + 1) - 1) // 3 * 3
    if tx.blocks and "N" in seq[cds_first + first_codon : cds_first + last_codon + 3]:
        result.flags.append("transcript_gap_here")
        result.terms.append("coding_sequence_variant")
        return

    alt_coding = seq[cds_first:a] + change + seq[b:]
    result.aa_position = max(a_cds, 0) // 3 + 1

    if a_cds < 3:
        if a_cds >= 0 and alt_coding[:3] == "ATG" and b_cds <= 3 and (len(change) - (b - a)) % 3 == 0:
            result.terms.append("start_retained_variant")
            result.hgvs_p = f"{tx.protein}:p.Met1="
            result.protein_short = "M1="
        else:
            result.terms.append("start_lost")
            result.hgvs_p = f"{tx.protein}:p.Met1?"
            result.protein_short = "M1?"
            if a_cds < 0:
                result.flags.append("overlaps_start_codon")
        return

    if b_cds > cds_length:
        result.terms.append("stop_lost")
        result.flags.append("overlaps_stop_codon")
        result.hgvs_p = f"{tx.protein}:p.?"
        return

    if (len(change) - (b - a)) % 3:
        _name_frameshift(tx, ref_protein, alt_coding, a_cds // 3, result)
    else:
        _name_in_frame(tx, ref_protein, seq[cds_first:], a_cds, b_cds, change, alt_coding, result)


def _tx_span(tx: Transcript, start: int, end: int) -> tuple[int, int] | None:
    """Transcript indices [a, b) of the replaced bases, or None if they do not lie in one exon."""
    if start == end:
        before = tx.exon_index_of(start - 1)
        after = tx.exon_index_of(start)
        if before is None or after is None or abs(before - after) != 1:
            return None
        index = max(before, after)
        return index, index
    first = tx.exon_index_of(start)
    last = tx.exon_index_of(end - 1)
    if first is None or last is None or abs(first - last) != end - start - 1:
        return None
    return (first, last + 1) if tx.strand > 0 else (last, first + 1)


def _name_frameshift(tx, ref_protein, alt_coding, near, result) -> None:
    # The protein is unchanged before the codon holding the change, so translate only from there.
    alt_tail = translate(alt_coding[near * 3 :])
    offset = 0
    while (
        offset < len(alt_tail)
        and near + offset < len(ref_protein)
        and ref_protein[near + offset] == alt_tail[offset]
    ):
        offset += 1
    first = near + offset
    result.terms.append("frameshift_variant")
    result.aa_position = first + 1
    if first >= len(ref_protein) or offset >= len(alt_tail):
        result.hgvs_p = f"{tx.protein}:p.?"
        return
    ref_aa, alt_aa = ref_protein[first], alt_tail[offset]
    if ref_aa == "*":
        result.terms.append("stop_lost")
        _name_extension(tx, alt_tail[offset:], first, result)
        return
    if alt_aa == "*":
        result.hgvs_p = f"{tx.protein}:p.{AA3[ref_aa]}{first + 1}Ter"
        result.protein_short = f"{ref_aa}{first + 1}*"
        return
    stop = alt_tail.find("*", offset)
    length = str(stop - offset + 1) if stop >= 0 else "?"
    result.hgvs_p = f"{tx.protein}:p.{AA3[ref_aa]}{first + 1}{AA3[alt_aa]}fsTer{length}"
    result.protein_short = f"{ref_aa}{first + 1}{alt_aa}fs*{length}"


def _name_extension(tx, alt_from_stop, stop_index, result) -> None:
    """alt_from_stop is the changed protein, starting at the residue that replaced the stop."""
    alt_aa = alt_from_stop[0]
    new_stop = alt_from_stop.find("*")
    distance = str(new_stop) if new_stop >= 0 else "?"
    result.hgvs_p = f"{tx.protein}:p.Ter{stop_index + 1}{AA3[alt_aa]}extTer{distance}"
    result.protein_short = f"*{stop_index + 1}{alt_aa}ext*{distance}"


def _name_in_frame(tx, ref_protein, coding, a_cds, b_cds, change, alt_coding, result) -> None:
    last_codon = (max(b_cds, a_cds + 1) - 1) // 3
    first_codon = a_cds // 3
    if a_cds == b_cds and a_cds % 3 == 0:
        first_codon -= 1  # insertion between codons: anchor on the residue before it
        last_codon = first_codon
    ref_codons = coding[first_codon * 3 : (last_codon + 1) * 3]
    alt_codons = (
        coding[first_codon * 3 : a_cds] + change + coding[b_cds : (last_codon + 1) * 3]
    )
    ref_seg = translate(ref_codons, to_stop=False)
    alt_seg = translate(alt_codons)

    begin = 0
    while begin < min(len(ref_seg), len(alt_seg)) and ref_seg[begin] == alt_seg[begin]:
        begin += 1
    ref_tail, alt_tail = len(ref_seg), len(alt_seg)
    if not alt_seg.endswith("*") or ref_seg.endswith("*"):
        while ref_tail > begin and alt_tail > begin and ref_seg[ref_tail - 1] == alt_seg[alt_tail - 1]:
            ref_tail, alt_tail = ref_tail - 1, alt_tail - 1
    removed = ref_seg[begin:ref_tail]
    added = alt_seg[begin:alt_tail]
    position = first_codon + begin  # 0-based residue of the first change

    if not removed and not added:
        aa = ref_seg[0]
        result.aa_position = first_codon + 1
        result.terms.append("stop_retained_variant" if aa == "*" else "synonymous_variant")
        label = f"{AA3[aa]}{first_codon + 1}="
        result.hgvs_p = f"{tx.protein}:p.{label}"
        result.protein_short = f"{aa}{first_codon + 1}="
        return

    result.aa_position = position + 1
    if removed.startswith("*"):
        result.terms.append("stop_lost")
        _name_extension(tx, translate(alt_coding[position * 3 :]), position, result)
        return
    if "*" in removed:
        result.terms.append("stop_lost")
        result.flags.append("overlaps_stop_codon")
        result.hgvs_p = f"{tx.protein}:p.?"
        return

    if not removed:
        _name_insertion(tx, ref_protein, position, added, result)
        return
    if not added:
        _name_deletion(tx, ref_protein, position, len(removed), result)
        return

    if len(removed) == 1 and len(added) == 1:
        ref_aa, alt_aa = removed, added
        if alt_aa == "*":
            result.terms.append("stop_gained")
        else:
            result.terms.append("missense_variant")
        result.hgvs_p = f"{tx.protein}:p.{AA3[ref_aa]}{position + 1}{AA3[alt_aa]}"
        result.protein_short = f"{ref_aa}{position + 1}{alt_aa}"
        return

    result.terms.append("stop_gained" if added.endswith("*") else "protein_altering_variant")
    last = position + len(removed)
    span = f"{AA3[removed[0]]}{position + 1}"
    short = f"{removed[0]}{position + 1}"
    if len(removed) > 1:
        span += f"_{AA3[removed[-1]]}{last}"
        short += f"_{removed[-1]}{last}"
    new = "".join(AA3[aa] for aa in added)
    result.hgvs_p = f"{tx.protein}:p.{span}delins{new}"
    result.protein_short = f"{short}delins{added}"


def _name_deletion(tx, ref_protein, position, size, result) -> None:
    # Protein names follow the 3' rule as well: use the C-terminal copy of a repeated residue.
    while position + size < len(ref_protein) - 1 and ref_protein[position] == ref_protein[position + size]:
        position += 1
    result.terms.append("inframe_deletion")
    result.aa_position = position + 1
    first, last = ref_protein[position], ref_protein[position + size - 1]
    if size == 1:
        result.hgvs_p = f"{tx.protein}:p.{AA3[first]}{position + 1}del"
        result.protein_short = f"{first}{position + 1}del"
    else:
        result.hgvs_p = f"{tx.protein}:p.{AA3[first]}{position + 1}_{AA3[last]}{position + size}del"
        result.protein_short = f"{first}{position + 1}_{last}{position + size}del"


def _name_insertion(tx, ref_protein, position, added, result) -> None:
    if added.endswith("*"):
        result.terms.append("stop_gained")
        before, after = ref_protein[position - 1], ref_protein[position]
        new = "".join(AA3[aa] for aa in added)
        result.hgvs_p = f"{tx.protein}:p.{AA3[before]}{position}_{AA3[after]}{position + 1}ins{new}"
        result.protein_short = f"{before}{position}_{after}{position + 1}ins{added}"
        return
    while position < len(ref_protein) - 1 and ref_protein[position] == added[0]:
        position += 1
        added = added[1:] + added[0]
    size = len(added)
    if position >= size and ref_protein[position - size : position] == added:
        result.terms.append("inframe_insertion")
        result.aa_position = position - size + 1
        first, last = ref_protein[position - size], ref_protein[position - 1]
        if size == 1:
            result.hgvs_p = f"{tx.protein}:p.{AA3[first]}{position}dup"
            result.protein_short = f"{first}{position}dup"
        else:
            result.hgvs_p = f"{tx.protein}:p.{AA3[first]}{position - size + 1}_{AA3[last]}{position}dup"
            result.protein_short = f"{first}{position - size + 1}_{last}{position}dup"
        return
    result.terms.append("inframe_insertion")
    result.aa_position = position
    before, after = ref_protein[position - 1], ref_protein[position]
    new = "".join(AA3[aa] for aa in added)
    result.hgvs_p = f"{tx.protein}:p.{AA3[before]}{position}_{AA3[after]}{position + 1}ins{new}"
    result.protein_short = f"{before}{position}_{after}{position + 1}ins{added}"
