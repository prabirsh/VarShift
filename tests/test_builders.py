"""Source parsers on real-format samples. A format change must fail loudly, never build empty tables."""

import gzip
import sqlite3
import tempfile
import unittest
from pathlib import Path

from tests.helpers import write_fasta
from varshift.errors import InputError, ReferenceMismatch, SetupError
from varshift.genome import Genome
from varshift.pipeline import builders
from varshift.variant import normalize

FIXTURES = Path(__file__).resolve().parent / "fixtures"


class Folder(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)

    def text(self, name, content, gz=False):
        path = self.root / name
        if gz:
            with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
        else:
            path.write_text(content, encoding="utf-8", newline="\n")
        return path


class RefSeq(unittest.TestCase):
    def test_real_sample(self):
        with open(FIXTURES / "refseq.sample.gff", encoding="utf-8") as handle:
            records = list(builders.iter_refseq_transcripts(handle, "hg19"))
        self.assertGreater(len(records), 10)
        by_id = {r["accession"]: r for r in records}
        self.assertTrue(all(a.startswith(("NM_", "NR_")) for a in by_id))
        olfactory = by_id["NM_001005484.2"]
        self.assertEqual((olfactory["gene"], olfactory["contig"], olfactory["strand"]), ("OR4F5", "chr1", 1))
        self.assertEqual(olfactory["protein"], "NP_001005484.2")
        self.assertEqual(olfactory["gene_id"], "79501")
        self.assertEqual(olfactory["hgnc_id"], "HGNC:14825")
        self.assertIn("RefSeq Select", olfactory["tags"])
        noted = by_id["NM_001282585.1"]
        self.assertIn("substitutions", noted["note"])

    def test_alt_contigs_and_model_transcripts_are_left_out(self):
        lines = [
            "NT_113878.1\tBestRefSeq\tgene\t1\t10\t.\t+\t.\tID=gene-A;gene=A\n",
            "NC_000001.10\tGnomon\tgene\t1\t100\t.\t+\t.\tID=gene-B;gene=B\n",
            "NC_000001.10\tGnomon\tmRNA\t1\t100\t.\t+\t.\tID=rna-XM_1.1;Parent=gene-B;gene=B;transcript_id=XM_1.1\n",
            "NC_000001.10\tGnomon\texon\t1\t100\t.\t+\t.\tID=exon-XM_1.1-1;Parent=rna-XM_1.1\n",
        ]
        self.assertEqual(list(builders.iter_refseq_transcripts(lines, "hg19")), [])


class Gencode(unittest.TestCase):
    def test_real_sample(self):
        with open(FIXTURES / "gencode.sample.gff3", encoding="utf-8") as handle:
            records = list(builders.iter_gencode_transcripts(handle, "hg38"))
        self.assertGreater(len(records), 3)
        self.assertTrue(all(r["accession"].startswith("ENST") for r in records))
        self.assertTrue(all(r["exons"] for r in records))
        self.assertEqual(records[0]["contig"], "chr1")

    def test_stop_codon_is_added_to_the_cds(self):
        lines = [
            "chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tID=ENSG1.1;gene_name=TEST\n",
            "chr1\tHAVANA\ttranscript\t1\t100\t.\t+\t.\tID=ENST1.1;Parent=ENSG1.1;transcript_id=ENST1.1;gene_name=TEST;transcript_type=protein_coding;protein_id=ENSP1.1\n",
            "chr1\tHAVANA\texon\t1\t100\t.\t+\t.\tID=exon:1;Parent=ENST1.1\n",
            "chr1\tHAVANA\tCDS\t10\t39\t.\t+\t0\tID=CDS:1;Parent=ENST1.1;protein_id=ENSP1.1\n",
            "chr1\tHAVANA\tstop_codon\t40\t42\t.\t+\t0\tID=stop:1;Parent=ENST1.1\n",
        ]
        record = next(builders.iter_gencode_transcripts(lines, "hg38"))
        row = builders._transcript_row(record, {})
        self.assertEqual(row[12:14], (9, 42))


class Chains(Folder):
    def test_plus_and_minus_blocks(self):
        chain = (
            "chain 1 1 1000 + 100 150 2 2000 + 500 550 1\n"
            "20 5 5\n"
            "25\n"
            "\n"
            "chain 1 1 1000 + 200 230 3 3000 - 100 130 2\n"
            "30\n"
        )
        out = self.root / "c.sqlite"
        stats = builders.build_liftover(self.text("c.chain.gz", chain, gz=True), out, "hg19", "hg38")
        self.assertEqual(stats["blocks"], 3)
        con = sqlite3.connect(out)
        rows = con.execute("SELECT contig,start,end,target_contig,target_start,target_end,target_strand FROM blocks ORDER BY start").fetchall()
        con.close()
        self.assertEqual(rows[0], ("chr1", 100, 120, "chr2", 500, 520, 1))
        self.assertEqual(rows[1], ("chr1", 125, 150, "chr2", 525, 550, 1))
        # On the minus strand the target runs from the far end of the 3000 bp contig.
        self.assertEqual(rows[2], ("chr1", 200, 230, "chr3", 2870, 2900, -1))

    def test_empty_chain_fails(self):
        with self.assertRaises(SetupError):
            builders.build_liftover(self.text("e.chain.gz", "# nothing\n", gz=True), self.root / "e.sqlite", "hg19", "hg38")


class ClinVar(Folder):
    def _genome(self):
        fasta = write_fasta(self.root, "g.fa", {"chr7": "ACGTACGTACGTTTTTACGT"})
        genome = Genome("hg19", fasta)
        self.addCleanup(genome.close)
        return genome

    def test_rows_counts_and_rejects(self):
        vcf = (
            "##fileDate=2026-10-25\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "7\t2\t100\tC\tT\t.\t.\tCLNSIG=Pathogenic;RS=123;GENEINFO=TEST:1;MC=SO:0001583|missense_variant;AF_TGP=0.01\n"
            "7\t3\t101\tA\tT\t.\t.\tCLNSIG=Benign\n"
            "7\t4\t102\tT\t.\t.\t.\tCLNSIG=Benign\n"
            "NT_1.1\t4\t103\tT\tA\t.\t.\tCLNSIG=Benign\n"
            "7\t14\t104\tTT\tT\t.\t.\tCLNSIG=Uncertain_significance\n"
        )
        out = self.root / "c.sqlite"
        stats = builders.build_clinvar(self.text("c.vcf.gz", vcf, gz=True), None, out, "hg19", self._genome())
        self.assertEqual(stats["file_date"], "2026-10-25")
        self.assertEqual(stats["rows"], 2)
        self.assertEqual(stats["ref_mismatch"], 1)
        self.assertEqual(stats["no_alt"], 1)
        self.assertEqual(stats["contig_not_installed"], 1)
        con = sqlite3.connect(out)
        row = con.execute("SELECT variation_id, rsid, classification, gene, molecular_consequence, af_tgp FROM clinvar WHERE pos=2").fetchone()
        # The deletion in the T run is stored left-aligned, so keys match normalised input.
        left = con.execute("SELECT pos, ref, alt FROM clinvar WHERE variation_id=104").fetchone()
        con.close()
        self.assertEqual(row, (100, "123", "Pathogenic", "TEST", "missense_variant", 0.01))
        self.assertEqual(left, (11, "GT", "G"))

    def test_names_come_from_the_right_assembly(self):
        vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n7\t2\t2\tC\tT\t.\t.\tCLNSIG=Benign\n"
        summary = (
            "#AlleleID\tName\tGeneSymbol\tAssembly\tVariationID\n"
            "1\tNM_1.1(TEST):c.5C>T (p.Ala2Val)\tTEST\tGRCh37\t2\n"
            "1\tNM_1.1(TEST):c.5C>T (p.Ala2Val)\tTEST\tGRCh38\t2\n"
            "2\tNM_1.1(TEST):c.9+1G>A\tTEST\tGRCh37\t3\n"
            "3\tsome free text\tTEST\tGRCh37\t4\n"
        )
        out = self.root / "n.sqlite"
        stats = builders.build_clinvar(
            self.text("n.vcf.gz", vcf, gz=True), self.text("s.txt.gz", summary, gz=True), out, "hg19", self._genome()
        )
        self.assertEqual(stats["named"], 2)
        con = sqlite3.connect(out)
        rows = con.execute("SELECT transcript, gene, hgvs_c, hgvs_p FROM clinvar_names ORDER BY variation_id").fetchall()
        con.close()
        self.assertEqual(rows[0], ("NM_1.1", "TEST", "c.5C>T", "p.Ala2Val"))
        self.assertEqual(rows[1], ("NM_1.1", "TEST", "c.9+1G>A", None))

    def test_a_changed_summary_format_fails_loudly(self):
        vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n7\t2\t2\tC\tT\t.\t.\tCLNSIG=Benign\n"
        with self.assertRaises(SetupError) as caught:
            builders.build_clinvar(
                self.text("x.vcf.gz", vcf, gz=True),
                self.text("x.txt", "Name\tGene\n"),
                self.root / "x.sqlite",
                "hg19",
                self._genome(),
            )
        self.assertIn("format has changed", str(caught.exception))

    def test_an_empty_vcf_fails_loudly(self):
        with self.assertRaises(SetupError):
            builders.build_clinvar(self.text("y.vcf", "#CHROM\n"), None, self.root / "y.sqlite", "hg19", self._genome())


class Genes(Folder):
    def test_a_changed_hgnc_format_fails_loudly(self):
        with self.assertRaises(SetupError):
            builders.build_genes({"hgnc": self.text("h.txt", "id\tname\nx\ty\n")}, self.root / "g.sqlite")


class TemporaryFiles(Folder):
    def test_a_build_sorts_beside_its_database_and_not_on_the_system_disk(self):
        folder = self.root / "it's here"
        folder.mkdir()
        con = builders._new_db(folder / "x.sqlite", "CREATE TABLE t(a);")
        self.addCleanup(con.close)
        self.assertEqual(Path(con.execute("PRAGMA temp_store_directory").fetchone()[0]), folder.resolve())
        # The setting belongs to the whole process, so don't leave it pointing at a folder about to be removed.
        con.execute("PRAGMA temp_store_directory=''")
        con.close()
        builders._open_builds.clear()


class ClinvarShortcut(Folder):
    """One changed base skips normalize() during a build, so the two must always agree."""

    SEQUENCE = "ACGTTGCAACGGTTAACCGGTATTTTTTGCA"

    def setUp(self):
        super().setUp()
        self.genome = Genome("hg19", write_fasta(self.root, "g.fa", {"chr1": self.SEQUENCE}))
        self.addCleanup(self.genome.close)

    def test_every_single_base_change_agrees_with_normalize(self):
        for pos, ref in enumerate(self.SEQUENCE, 1):
            for alt in "ACGT".replace(ref, ""):
                variant = normalize(self.genome, "chr1", pos, ref, alt)
                self.assertEqual(builders._normalized_quickly(self.genome, "chr1", pos, ref, alt), (variant.pos, variant.ref, variant.alt))

    def test_a_wrong_reference_base_is_still_refused(self):
        with self.assertRaises(ReferenceMismatch):
            builders._normalized_quickly(self.genome, "chr1", 1, "T", "G")

    def test_a_position_past_the_end_is_still_refused(self):
        with self.assertRaises(InputError):
            builders._normalized_quickly(self.genome, "chr1", len(self.SEQUENCE) + 5, "A", "G")

    def test_an_insertion_is_still_left_aligned(self):
        self.assertEqual(builders._normalized_quickly(self.genome, "chr1", 28, "T", "TT"), (22, "A", "AT"))

    def test_lower_case_and_ambiguous_alleles_take_the_full_route(self):
        self.assertEqual(builders._normalized_quickly(self.genome, "chr1", 1, "a", "g"), (1, "A", "G"))
        with self.assertRaises(InputError):
            builders._normalized_quickly(self.genome, "chr1", 1, "A", "R")


class Uniprot(unittest.TestCase):
    ENTRY = """ID   TEST_HUMAN              Reviewed;         10 AA.
AC   P00001; Q00001;
GN   Name=TEST {ECO:0000303}; Synonyms=OLD;
DR   Pfam; PF00001; 7tm_1; 1.
FT   DOMAIN          2..>8
FT                   /note="Kinase"
FT   ACT_SITE        5
FT   CHAIN           1..10
FT                   /note="ignored"
SQ   SEQUENCE   10 AA;  1000 MW;  ABC CRC64;
     MAEKLGWPAA
//
ID   TREMBL_HUMAN            Unreviewed;         10 AA.
AC   A00001;
SQ   SEQUENCE   3 AA;  1000 MW;  ABC CRC64;
     MAE
//
"""

    def test_flat_file(self):
        entries = list(builders.iter_uniprot(self.ENTRY.splitlines(keepends=True)))
        self.assertEqual(len(entries), 1, "unreviewed entries are dropped")
        entry = entries[0]
        self.assertEqual((entry["accession"], entry["gene"], entry["sequence"]), ("P00001", "TEST", "MAEKLGWPAA"))
        self.assertEqual(entry["pfam"], [("PF00001", "7tm_1")])
        self.assertEqual(entry["features"], [["DOMAIN", 2, 8, "Kinase"], ["ACT_SITE", 5, 5, ""]])

    # The first lines of Pfam's real human file, then a row in its format for each entry above.
    PFAM_HEAD = (
        "#Pfam-A regions from Pfam version 37.0 for ncbi taxid 9606 'Homo sapiens (Human)'\n"
        "#Total number of proteins in proteome: 82677\n"
        "#<seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> "
        "<hmm start> <hmm end> <hmm length> <bit score> <E-value> <clan>\n"
    )
    PFAM_ROWS = (
        "P00001\t2\t8\t1\t9\tPF00001\t7tm_1\tFamily\t1\t7\t7\t10.0\t1e-05\tNo_clan\n"
        "A00001\t1\t3\t1\t3\tPF00002\tOther\tDomain\t1\t3\t3\t10.0\t1e-05\tCL0001\n"
    )

    def build(self, pfam):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        uniprot = Path(folder.name) / "uniprot.dat"
        uniprot.write_text(self.ENTRY, encoding="utf-8")
        out = Path(folder.name) / "proteins.sqlite"
        stats = builders.build_proteins(uniprot, pfam.splitlines(keepends=True), out)
        con = sqlite3.connect(out)
        self.addCleanup(con.close)
        return stats, con.execute("SELECT * FROM pfam_regions").fetchall()

    def test_pfam_rows_are_kept_for_reviewed_proteins_only(self):
        stats, rows = self.build(self.PFAM_HEAD + self.PFAM_ROWS)
        self.assertEqual(rows, [("P00001", "PF00001", "7tm_1", 2, 8)])
        self.assertEqual(stats["pfam_regions"], 1)

    def test_a_changed_pfam_header_stops_the_build(self):
        with self.assertRaises(SetupError):
            self.build(self.PFAM_HEAD.replace("<hmm acc>", "<family>") + self.PFAM_ROWS)

    def test_a_pfam_file_with_no_header_stops_the_build(self):
        with self.assertRaises(SetupError):
            self.build(self.PFAM_ROWS)

    def test_a_pfam_file_matching_no_protein_stops_the_build(self):
        with self.assertRaises(SetupError):
            self.build(self.PFAM_HEAD + "Q99999\t1\t3\t1\t3\tPF00002\tOther\tDomain\t1\t3\t3\t10.0\t1e-05\tCL0001\n")


if __name__ == "__main__":
    unittest.main()
