# Validation

This page says what was tested, how, what came out, and what hasn't been proven.

**VarSight has not been validated for clinical use.** A laboratory using it for patient care must do its own validation.

## Test set

- **Genes:** BRAF, KRAS, NRAS, EGFR, TP53, PIK3CA, IDH1, BRCA1, BRCA2, JAK2, KIT, PTEN, with a mix of forward-strand and reverse-strand genes.
- **Sequence:** real hg19 and hg38 windows around each gene, from UCSC.
- **Transcripts:** the real RefSeq annotation files, filtered to these genes.
- **ClinVar:** every record in those windows (ClinVar version 2026-10-25, about 61,000 per build).
- **Names:** ClinVar's own c. and p. names, from the first 30 MB of `variant_summary`, giving 13,724 usable names per build.
- **Bias:** those names are ClinVar's oldest records, in well-studied genes. Expect the whole genome to score lower.

Date: 1 November 2026. Python 3.14.4 on WSL (Ubuntu 26.04).

## Results

The full report, with examples from every mismatch group, is in [validation-results/clinvar-12-genes.md](validation-results/clinvar-12-genes.md).

| Test | What it checks | hg19 | hg38 |
|---|---|---:|---:|
| `forward_c_snv` | Genome position gives ClinVar's c. name | 99.98% (9,085 of 9,087) | 99.99% (9,086 of 9,087) |
| `forward_c_indel` | The same, for indels | 99.83% (4,022 of 4,029) | 99.83% (4,022 of 4,029) |
| `forward_p_snv` | Genome position gives ClinVar's p. name | 99.99% (6,945 of 6,946) | 99.99% (6,945 of 6,946) |
| `forward_p_indel` | The same, for indels | 99.45% (3,617 of 3,637) | 99.45% (3,617 of 3,637) |
| `reverse_c_snv` | c. name gives ClinVar's position | 99.99% (9,085 of 9,086) | 99.99% (9,086 of 9,087) |
| `reverse_c_indel` | The same, for indels | 99.93% (4,026 of 4,029) | 99.93% (4,026 of 4,029) |
| `reverse_p_snv` | ClinVar's position is among the candidates for its p. name | 100% (6,948 of 6,948) | 100% (6,948 of 6,948) |
| `reverse_p_indel` | The same, for indels | 100% (3,640 of 3,640) | 100% (3,640 of 3,640) |
| `cross_build_snv` | Lifted position equals ClinVar's position on the other build | 100% (9,653 of 9,653) | 100% (9,653 of 9,653) |
| `cross_build_indel` | The same, for indels | 100% (4,070 of 4,070) | 100% (4,070 of 4,070) |
| `consequence_snv` | Effect agrees with ClinVar's molecular consequence | 99.02% (9,072 of 9,162) | 99.02% (9,073 of 9,163) |
| `consequence_indel` | The same, for indels | 99.55% (4,002 of 4,020) | 99.55% (4,002 of 4,020) |

### How names were compared

To avoid counting spelling differences as errors, these are treated as the same:

- Any frameshift form: `p.Glu23fs`, `p.Glu23Valfs`, `p.Glu23ValfsTer17`.
- Any start-codon change: `p.Met1?`, `p.Met1Val`.
- A stop written as an insertion: `p.Glu4_Ser5insTer` and `p.Ser5Ter`.
- An extension with or without its length: `p.Ter9Gln` and `p.Ter9GlnextTer3`.
- Older names listing deleted bases: `c.548delG` and `c.548del`.
- ClinVar's repeat form, such as `c.3901ACT[1]`, counts as a match if it reaches the same genome position (84 names).

Not scored: 4 names using protein repeat notation, and names on other genes' transcripts that the 12-gene set doesn't include.

### What the mismatches are

- **`forward_p_indel`:** mostly a frameshift ending within a few codons, where ClinVar writes `p.Ser205delinsThrTer` and VarSight writes `p.Ser205ThrfsTer2`. Same protein, different spelling. One case, `c.75_80dup`, gets no protein name from VarSight.
- **`forward_c_indel`:** ClinVar names that weren't moved to the 3' end, `delins` where VarSight writes `inv`, and a deletion ClinVar names across two neighbouring exons (`c.8331_8332del`, which VarSight writes as `c.8330_8331del`).
- **`consequence_snv`:** nearly all are PTEN variants upstream of where the hg19 annotation starts that transcript. ClinVar's term comes from a longer, older version. This is a data difference.
- **One SNV in each test** has an ambiguous base in ClinVar (`c.4965C>R`).

### Lower scores on other versions

ClinVar sometimes names an older transcript version than the one installed. Those groups are scored separately and score much lower (31% to 68%), because UTR numbering shifts between versions. **This is a real limit:** `c.-1060C>G` on `NM_000314.4` is a different base from the same name on `NM_000314.8`. VarSight flags this as `transcript_version_substituted`. Check any result carrying that flag.

## Comparison with DRAGEN

Run on local annotated VCFs from the author's samples. Only counts are recorded.

- 16 variants fell inside the 12 genes, all indels, giving 114 names on RefSeq transcripts.
- 113 of 114 c. names were the same. The other was a 796 bp deletion in a 3' UTR, placed 2 bases further 3' by VarSight.
- This is a very small check, because the installed sample only covers 12 genes.

After a full setup, you can run the same comparison on your own files. It stays on your computer:

```bash
python3 tools/compare_annotated_vcf.py sample.annotated.vcf.gz --differences differences.txt
```

## Unit and integration tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -t .
```

The suite covers:

- **HGVS naming on both strands:** a hand-written gene where every expected name was worked out by hand. Every case must give the same answer on the forward and reverse strands.
- **Coordinates:** every base of a transcript round-trips between genome and c. numbering.
- **Parsers:** many input spellings, change types and hostile text.
- **Source parsers:** real-format samples, and failure when a format changes.
- **Downloads:** against a real local HTTP server, including dropped connections, resume, bad checksums, damaged archives, oversized files and blocked redirects.
- **The data lifecycle:** setup, update, rollback, refused updates, approval and tampering.
- **Security rules:** the source code itself is read, and tests fail if forbidden patterns appear.

## Speed and memory

Measured with public ClinVar records from the 12 genes, in a VCF of 3,409 lines, with the data on the WSL Linux disk. WSL had 1.8 GB of memory.

| Task | Time | Peak memory |
|---|---:|---:|
| One variant | About 3 to 5 ms | |
| VCF to TSV | 10 to 15 s | About 100 MB |
| VCF to Excel | About 15 s | About 150 MB |
| VCF to JSON Lines | 10 to 14 s | About 150 MB |

The test databases are much smaller than full ones, so expect slower results on real data, particularly on the first run after a reboot.

**Keep the data on the Linux disk.** SQLite lookups measured 27 microseconds on the WSL Linux disk and 1.3 to 1.8 ms on a Windows drive seen from WSL (`/mnt/i`), about 50 to 65 times slower.

## Not proven yet

| What | Status | How to close it |
|---|---|---|
| Full-genome accuracy | Only 12 genes tested | `varsight validate --sample 20000 --report validation.md` after setup |
| Full setup | Sources and formats checked, and the cycle tested on a local server, but the 8 GB run has not been done | Run `varsight setup` |
| Real setup time and size | Estimates | The report in `reports/` records them |
| GENCODE, AlphaMissense and gnomAD builders | Tested on small real-format samples only | Run them on full files |
| Python 3.10 to 3.13 | Not run yet | CI, once on GitHub |
| Excel itself | Checked with a separate reader (openpyxl) with zero formula cells, but not opened in Microsoft Excel | Open a result in Excel |
| Link formats | See [the table](data-sources.md#reference-links); VarSome, OMIM and UCSC unconfirmed | Click them in a normal browser |
| The weekly Windows task | The script's syntax was checked, but it was not registered or run | Follow [the guide](updates.md) |
