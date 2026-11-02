# Validation

This page says what was tested, how, what came out, and what hasn't been proven.

**VarSight has not been validated for clinical use.** A laboratory using it for patient care must do its own validation.

## Whole-genome sample

On 2 November 2026, `varsight validate sample=20000` ran on a complete hg19 setup. It takes evenly spaced names from all of ClinVar, 20,027 in this run, so it is a much harder test than the 12 genes below. A second run from WSL, on 1,998 names, gave similar rates. hg38 hasn't been measured this way, and no cross-build check was possible with one build installed.

The full report, with examples from every mismatch group, is in [validation-results/clinvar-hg19-whole-genome.md](validation-results/clinvar-hg19-whole-genome.md).

| Test | What it checks | Passed | Rate |
| - | - | -: | -: |
| `forward_c_snv` | Genome position gives ClinVar's c. name | 18,448 of 18,539 | 99.51% |
| `forward_c_indel` | The same, for indels | 1,386 of 1,413 | 98.09% |
| `forward_p_snv` | Genome position gives ClinVar's p. name | 14,917 of 15,010 | 99.38% |
| `forward_p_indel` | The same, for indels | 956 of 996 | 95.98% |
| `reverse_c_snv` | c. name gives ClinVar's position | 18,448 of 18,539 | 99.51% |
| `reverse_c_indel` | The same, for indels | 1,397 of 1,413 | 98.87% |
| `reverse_p_snv` | ClinVar's position is among the candidates for its p. name | 15,022 of 15,048 | 99.83% |
| `reverse_p_indel` | The same, for indels | 991 of 996 | 99.50% |
| `consequence_snv` | Effect agrees with ClinVar's molecular consequence | 18,529 of 18,564 | 99.81% |
| `consequence_indel` | The same, for indels | 1,387 of 1,398 | 99.21% |

### What the mismatches are

Every mismatch in the four `forward` tests, sorted by cause:

| Cause | c. SNV | c. indel | p. SNV | p. indel |
| - | -: | -: | -: | -: |
| The transcript has a gap against hg19 | 84 | 10 | 71 | 10 |
| The same change, spelled or placed differently | 0 | 16 | 0 | 28 |
| No p. name, because the transcript has an in-frame stop on hg19 | 0 | 0 | 19 | 0 |
| Other, including transcripts placed differently on hg19 | 7 | 1 | 3 | 2 |
| **Total** | **91** | **27** | **93** | **40** |

- **Gaps are a real limit in VarSight, and the main one found so far.** About 1,100 of the 82,539 RefSeq transcripts on hg19, 275 of them MANE Select, have an indel or frameshift against the genome. VarSight numbers a transcript from its exons and doesn't yet use RefSeq's alignment of the gap, so everything after the gap is off by its length. In ALMS1 `NM_001378454.1`, `chr2-73675519-A-G` gets `c.1862A>G` and `p.Tyr621Cys`, where ClinVar has `c.1865A>G` and `p.Tyr622Cys`.
- **The same limit affects c. names given as input.** Of the 25 mismatches listed in the report for `reverse_c_snv`, 17 were refused with `no transcript fits`, and 8 were placed 3 bases from ClinVar's position because the wrong base happened to match. Results on these transcripts carry the warning `transcript_differs_from_genome`. **Don't rely on a name or position that carries this warning without checking it against the transcript.**
- **Same change, spelled differently:** ClinVar's `p.Glu43_Ser44insGlu` is VarSight's `p.Glu43dup`, its `p.Ser290delinsArgTer` is VarSight's `p.Ser290ArgfsTer2`, and its `c.1366-23dup` is VarSight's `c.1366-15dup`, moved to the 3' end of the repeat. Counting these as agreement gives 99.22% for c. indels and 98.80% for p. indels.
- **No p. name:** VarSight refuses to translate a transcript with an in-frame stop on hg19 and flags it `cds_internal_stop`. This covers selenoproteins such as SELENON, and genes such as ACTN3.

Names on a different transcript version are scored separately in the report, 53 names here, for the reason given under [lower scores on other versions](#lower-scores-on-other-versions).

## 12-gene test set

- **Genes:** BRAF, KRAS, NRAS, EGFR, TP53, PIK3CA, IDH1, BRCA1, BRCA2, JAK2, KIT, PTEN, with a mix of forward-strand and reverse-strand genes.
- **Sequence:** real hg19 and hg38 windows around each gene, from UCSC.
- **Transcripts:** the real RefSeq annotation files, filtered to these genes.
- **ClinVar:** every record in those windows (ClinVar version 2026-10-25, about 61,000 per build).
- **Names:** ClinVar's own c. and p. names, from the first 30 MB of `variant_summary`, giving 13,724 usable names per build.
- **Bias:** those names are ClinVar's oldest records, in well-studied genes, which is why this set scores higher than the whole-genome sample. It is also the only test so far on hg38 and of mapping between builds.

Date: 1 November 2026. Python 3.14.4 on WSL (Ubuntu 26.04). The whole run was repeated from freshly downloaded sources after the last code change on that date, and gave the same numbers.

## Results on the 12 genes

The full report, with examples from every mismatch group, is in [validation-results/clinvar-12-genes.md](validation-results/clinvar-12-genes.md).

| Test | What it checks | hg19 | hg38 |
| - | - | -: | -: |
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
python3 tools/compare_annotated_vcf.py sample.annotated.vcf.gz differences=differences.txt
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
- **Advanced databases:** samples cut from each real public file, and each provider's documented columns where a file needs an account. Adding, using and removing them against a local server.
- **The command line and setup questions:** every answer to the setup menu, and each command in these docs is parsed by the real parser.
- **The update notice:** when it appears, when it doesn't, and that it never stops annotation.
- **Reports:** the sheets and their order in Excel, sorting, hostile text in the HTML report, and the row limit.
- **House rules:** a test reads every file in the package for double hyphens and long dashes.
- **Downloads:** against a real local HTTP server, including dropped connections, resume, bad checksums, damaged archives, oversized files and blocked redirects.
- **The data lifecycle:** setup, update, rollback, refused updates, approval and tampering.
- **Security rules:** the source code itself is read, and tests fail if forbidden patterns appear.

## Install and command checks

Run on 1 November 2026, from a copy holding only the files git would commit:

| Check | Result |
| - | - |
| `sh install.sh test=yes` on WSL Ubuntu | Installed, all tests passed, demo finished with 11 of 11 known variants |
| `install.ps1 -Test` in PowerShell | The same on Windows, with 4 tests skipped because they cover POSIX file modes |
| `install.bat` in Command Prompt | Installed. `varsight version` and `varsight demo` work there |
| `pip install .` in a new Windows environment | Installed with no dependencies, and the report template is packaged |
| About 50 command lines through the installed launcher, on both systems | Exit code and text checked for every command, including wrong input. `NM_004333.6:c.1799T>A` keeps its `>` from PowerShell and Command Prompt |
| `workers=2` against `workers=1` | Byte-for-byte the same file |
| `tools/check_sources.py` | All 38 download sources answered, within their size limits |

Testing against the real Cancer Hotspots file found a bug the small sample couldn't show: a wide in-frame indel region hid the single-residue BRAF V600 row. It is fixed, and those two rows are now in the test sample.

## Speed and memory

See [benchmarks](benchmarks.md).

## Not proven yet

| What | Status | How to close it |
| - | - | - |
| Transcripts with a gap against hg19 | c. and p. numbers after the gap are off by its length, and a c. name as input can be refused or misplaced. They carry a warning | Read RefSeq's alignment of the gap |
| Whole-genome accuracy on hg38 | hg19 is measured on a 20,027-name sample, but hg38 only on 12 genes | `varsight validate sample=20000 report=validation.md` after an hg38 setup |
| Full setup | hg19 ran from start to the self-test on 2 November 2026, natively on Windows, and all 5 known variants passed. hg38, both builds together, and a full build inside WSL haven't been run | Run `varsight setup` |
| Real setup time and size | [Measured for hg19](benchmarks.md#setup-time). Estimates for hg38 | The report in `reports/` records them |
| Current Pfam positions | Positions come from Pfam 37.0. [What that costs](data-sources.md#core-databases) | A newer human-only file from Pfam |
| GENCODE and the advanced database builders | Tested on samples of real files only. See [how much of this has been tested](data-sources.md#how-much-of-this-has-been-tested) | Run them on full files |
| COSMIC, Census, fusions, OMIM, SpliceAI, PrimateAI | Written from documented columns, without a real file | Add your own file. It stops with a clear message if the columns differ |
| Python 3.10 to 3.13 | Not run yet | CI, once on GitHub |
| Excel itself | Checked with a separate reader (openpyxl): five sheets in order, filters, frozen headers, colours and links, with no warnings. Not opened in Microsoft Excel | Open a result in Excel |
| HTML report | Tested in Chromium: filters, search, sort, details panel, keyboard, TSV export, light and dark, phone width. Not tested in Firefox or Safari | Open one there |
| Link formats | See [the table](data-sources.md#reference-links); VarSome, OMIM and UCSC unconfirmed | Click them in a normal browser |
| The weekly Windows task | A dry run reaches the register step in both forms, and the exact commands it launches were run and worked. The task itself has not been registered, since that changes your computer's schedule | Follow [the guide](updates.md) |
