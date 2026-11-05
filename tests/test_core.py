"""Reference access, normalisation, contig names and coordinate arithmetic."""

import tempfile
import unittest
from pathlib import Path

from tests.helpers import CONTIG, FORWARD, ToyGene, write_fasta
from varshift.batch import _vcf_build, split_queries
from varshift.builds import (
    canonical_contig,
    contig_accession,
    detect_build_from_lengths,
    key_contig,
    parse_build,
)
from varshift.errors import DataError, InputError, ReferenceMismatch, ResolveError
from varshift.genome import Genome, write_fai
from varshift.store import DataStore, bins_for_range, decode_exons, encode_exons, ucsc_bin
from varshift.transcripts import TxPosition
from varshift.variant import Variant, normalize
from varshift.vcf import VcfReader


class VcfBuild(unittest.TestCase):
    def _reader(self, header: str):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        path = Path(folder.name) / "s.vcf"
        path.write_text(header + "#CHROM\tPOS\tID\tREF\tALT\n", encoding="utf-8")
        reader = VcfReader(path)
        self.addCleanup(reader.close)
        return reader

    def test_build_contradicting_the_reference_line_stops(self):
        reader = self._reader("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        with self.assertRaises(InputError):
            _vcf_build(reader, "hg19")
        self.assertEqual(_vcf_build(reader, "hg38"), "hg38")
        self.assertEqual(_vcf_build(reader, None), "hg38")


class SplitQueries(unittest.TestCase):
    def test_semicolons_separate_queries(self):
        self.assertEqual(split_queries("a; b;c"), ["a", "b", "c"])

    def test_a_semicolon_inside_square_brackets_is_part_of_the_change(self):
        self.assertEqual(split_queries("NM_1:c.[76A>C;83G>C]; BRAF V600E"), ["NM_1:c.[76A>C;83G>C]", "BRAF V600E"])

    def test_empty_parts_are_dropped(self):
        self.assertEqual(split_queries(";;"), [])
        self.assertEqual(split_queries(" BRAF V600E ;; "), ["BRAF V600E"])

    def test_one_query_stays_one(self):
        self.assertEqual(split_queries("BRAF V600E"), ["BRAF V600E"])

    def test_an_unclosed_bracket_keeps_the_rest_together(self):
        self.assertEqual(split_queries("NM_1:c.[76A>C;83G>C"), ["NM_1:c.[76A>C;83G>C"])


class TranscriptCache(unittest.TestCase):
    def test_one_accession_at_two_places_on_a_chromosome_stays_two_transcripts(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        (Path(folder.name) / "manifest.json").write_text('{"schema": 1, "builds": ["hg19"], "files": {}}', encoding="utf-8")
        store = DataStore(Path(folder.name))
        self.addCleanup(store.close)
        row = {
            "accession": "NM_000001.1", "gene": "G", "contig": "chr1", "strand": 1, "exons": "100-200",
            "cds_start": None, "cds_end": None, "biotype": "mRNA", "protein": "", "source": "RefSeq", "gene_id": "",
            "hgnc_id": "", "mane": "", "mane_accession": "", "refseq_select": 0, "genome_note": "",
            "translation_exception": 0, "alignment": "", "tx_low": 100,
        }  # fmt: skip
        first = store._make_transcript("hg19", row)
        second = store._make_transcript("hg19", {**row, "exons": "900000-900100", "tx_low": 900000})
        self.assertEqual((first.tx_low, second.tx_low), (100, 900000))


class GenomeReading(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)

    def test_fetch_across_line_breaks_matches_python_slice(self):
        sequence = "".join("ACGT"[(i * i) % 4] for i in range(257))
        genome = Genome("hg38", write_fasta(self.root, "a.fa", {"chr1": sequence, "chr2": sequence[::-1]}, width=7))
        self.addCleanup(genome.close)
        for start, end in ((0, 1), (6, 8), (7, 14), (100, 200), (250, 257), (0, 257)):
            self.assertEqual(genome.fetch("chr1", start, end), sequence[start:end])
            self.assertEqual(genome.fetch("chr2", start, end), sequence[::-1][start:end])

    def test_soft_masked_bases_come_back_uppercase(self):
        genome = Genome("hg38", write_fasta(self.root, "b.fa", {"chr1": "acgtACGTnn"}))
        self.addCleanup(genome.close)
        self.assertEqual(genome.fetch("chr1", 0, 10), "ACGTACGTNN")

    def test_fetch_is_clipped_at_the_end(self):
        genome = Genome("hg38", write_fasta(self.root, "c.fa", {"chr1": "ACGT"}))
        self.addCleanup(genome.close)
        self.assertEqual(genome.fetch("chr1", 2, 100), "GT")
        self.assertEqual(genome.fetch("chr1", 100, 200), "")

    def test_window_records_keep_real_coordinates(self):
        fasta = write_fasta(self.root, "d.fa", {"chr7:1001-1010": "ACGTACGTAC", "chr7:2001-2005": "TTTTT"})
        genome = Genome("hg38", fasta)
        self.addCleanup(genome.close)
        self.assertEqual(genome.fetch("chr7", 1000, 1004), "ACGT")
        self.assertEqual(genome.fetch("chr7", 2000, 2002), "TT")
        self.assertTrue(genome.covers("chr7", 1000, 1010))
        self.assertFalse(genome.covers("chr7", 1005, 1500))

    def test_missing_files_give_a_clear_error(self):
        with self.assertRaises(DataError):
            Genome("hg38", self.root / "absent.fa")

    def test_uneven_lines_are_refused(self):
        path = self.root / "bad.fa"
        path.write_text(">chr1\nACGT\nAC\nACGT\n", encoding="ascii")
        with self.assertRaises(DataError):
            write_fai(path)


class Normalisation(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        #            1234567890123456
        sequence = "GGCATATATATCCTTT"
        self.genome = Genome("hg38", write_fasta(Path(self.folder.name), "n.fa", {"chrT": sequence}))
        self.addCleanup(self.genome.close)

    def test_snv_is_unchanged(self):
        self.assertEqual(normalize(self.genome, "chrT", 3, "C", "G").key, "chrT-3-C-G")

    def test_deletion_in_a_repeat_moves_left(self):
        # Deleting one AT from the ATATATAT run is reported at its left end.
        self.assertEqual(normalize(self.genome, "chrT", 9, "TAT", "T").key, "chrT-3-CAT-C")

    def test_insertion_in_a_repeat_moves_left(self):
        self.assertEqual(normalize(self.genome, "chrT", 11, "T", "TAT").key, "chrT-3-C-CAT")

    def test_shared_bases_are_trimmed(self):
        self.assertEqual(normalize(self.genome, "chrT", 3, "CAT", "GAT").key, "chrT-3-C-G")
        self.assertEqual(normalize(self.genome, "chrT", 2, "GCA", "GTA").key, "chrT-3-C-T")

    def test_lowercase_alleles_are_accepted(self):
        self.assertEqual(normalize(self.genome, "chrT", 3, "c", "g").key, "chrT-3-C-G")

    def test_wrong_reference_is_rejected_with_what_was_found(self):
        with self.assertRaises(ReferenceMismatch) as caught:
            normalize(self.genome, "chrT", 3, "A", "G")
        self.assertEqual(caught.exception.found, "C")

    def test_identical_alleles_are_rejected(self):
        with self.assertRaises(InputError):
            normalize(self.genome, "chrT", 3, "C", "C")

    def test_non_dna_is_rejected(self):
        for allele in ("<DEL>", "*", "X", ""):
            with self.assertRaises(InputError):
                normalize(self.genome, "chrT", 3, "C", allele)

    def test_positions_outside_the_sequence_are_rejected(self):
        for position in (0, -1, 17):
            with self.assertRaises(InputError):
                normalize(self.genome, "chrT", position, "C", "G")

    def test_edit_removes_padding(self):
        self.assertEqual(Variant("hg38", "chrT", 3, "CAT", "C").edit(), (3, 5, ""))
        self.assertEqual(Variant("hg38", "chrT", 3, "C", "CAT").edit(), (3, 3, "AT"))
        self.assertEqual(Variant("hg38", "chrT", 3, "C", "G").edit(), (2, 3, "G"))


class ContigNames(unittest.TestCase):
    def test_spellings_of_one_contig(self):
        for name in ("7", "chr7", "CHR7", "NC_000007.13"):
            self.assertEqual(canonical_contig(name, "hg19"), "chr7")
        self.assertEqual(canonical_contig("x", "hg38"), "chrX")

    def test_accession_of_the_wrong_build_is_refused(self):
        with self.assertRaises(InputError):
            canonical_contig("NC_000007.14", "hg19")

    def test_mitochondrial_names(self):
        self.assertEqual(canonical_contig("MT", "hg38"), "chrM")
        self.assertEqual(canonical_contig("chrM", "hg38"), "chrM")
        self.assertEqual(canonical_contig("MT", "hg19"), "chrMT")
        self.assertEqual(canonical_contig("chrM", "hg19"), "chrM")
        self.assertEqual(canonical_contig("NC_012920.1", "hg19"), "chrMT")
        self.assertEqual(key_contig("chrMT"), "MT")

    def test_impossible_chromosome(self):
        with self.assertRaises(InputError):
            canonical_contig("chr99", "hg19")

    def test_accessions(self):
        self.assertEqual(contig_accession("chr7", "hg19"), "NC_000007.13")
        self.assertEqual(contig_accession("chr7", "hg38"), "NC_000007.14")
        self.assertEqual(contig_accession("chrMT", "hg19"), "NC_012920.1")
        self.assertIsNone(contig_accession("chrM", "hg19"))

    def test_build_aliases(self):
        self.assertEqual(parse_build("GRCh37"), "hg19")
        self.assertEqual(parse_build("hg38"), "hg38")
        with self.assertRaises(InputError):
            parse_build("hg18")

    def test_build_from_contig_lengths(self):
        self.assertEqual(detect_build_from_lengths({"chr1": 249250621, "chr2": 243199373})[0], "hg19")
        self.assertEqual(detect_build_from_lengths({"1": 248956422})[0], "hg38")
        self.assertIsNone(detect_build_from_lengths({})[0])
        build, wrong = detect_build_from_lengths({"chr1": 249250621, "chr2": 243199373, "chr3": 5})
        self.assertEqual((build, wrong), ("hg19", ["chr3"]))


class Coordinates(unittest.TestCase):
    def _gene(self, strand):
        gene = ToyGene(strand)
        self.addCleanup(gene.close)
        return gene

    def test_every_transcript_base_round_trips(self):
        for strand in (1, -1):
            tx = self._gene(strand).transcript
            for genomic in range(tx.tx_low - 10, tx.tx_high + 10):
                position = tx.position_of(genomic)
                self.assertEqual(tx.genomic_of(position), genomic, f"{position} strand {strand}")

    def test_spliced_sequence_matches_on_both_strands(self):
        forward = self._gene(1)
        reverse = self._gene(-1)
        self.assertEqual(forward.transcript.spliced(forward.genome), reverse.transcript.spliced(reverse.genome))

    def test_named_positions(self):
        gene = self._gene(1)
        tx = gene.transcript
        expected = {
            25: "1",
            24: "-1",
            39: "15",
            40: "15+1",
            50: "15+11",
            51: "16-10",
            60: "16-1",
            61: "16",
            72: "27",
            73: "*1",
        }
        for genomic, label in expected.items():
            self.assertEqual(str(tx.position_of(genomic)), label)
        self.assertEqual(FORWARD[25:28], "ATG")

    def test_offset_on_a_base_inside_an_exon_is_refused(self):
        tx = self._gene(1).transcript
        with self.assertRaises(ResolveError):
            tx.genomic_of(TxPosition("", 5, 2))
        with self.assertRaises(ResolveError):
            tx.genomic_of(TxPosition("", 15, -2))

    def test_offset_past_the_intron_is_refused(self):
        tx = self._gene(1).transcript
        with self.assertRaises(ResolveError):
            tx.genomic_of(TxPosition("", 15, 500))

    def test_position_past_the_cds_is_refused(self):
        with self.assertRaises(ResolveError):
            self._gene(1).transcript.genomic_of(TxPosition("", 28))

    def test_contig_constant_is_used(self):
        self.assertEqual(self._gene(1).transcript.contig, CONTIG)


class Binning(unittest.TestCase):
    def test_a_range_is_found_through_its_bin(self):
        for start, end in ((0, 1), (131071, 131073), (140453135, 140453136), (1, 5_000_000), (100_000, 100_000_000)):
            self.assertIn(ucsc_bin(start, end), bins_for_range(start, end))
            self.assertIn(ucsc_bin(start, end), bins_for_range(start, start + 1) + bins_for_range(end - 1, end))

    def test_exon_text_round_trips(self):
        exons = ((10, 20), (30, 40))
        self.assertEqual(decode_exons(encode_exons(exons)), exons)


if __name__ == "__main__":
    unittest.main()
