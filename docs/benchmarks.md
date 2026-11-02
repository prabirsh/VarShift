# Benchmarks and system requirements

Run `varsight benchmark` to measure your own computer. All VarSight figures below come from one laptop: Intel i3-7020U (2 cores, 4 threads, 2.3 GHz), 3.9 GB of memory, SSD. WSL had 1.8 GB of that memory. The test VCFs hold public ClinVar records from the 12 demo genes, on the demo data.

**These are gene-dense records on small databases.** Full databases and the first run after a reboot will be slower, and records outside genes should be faster. Neither has been measured yet.

## Speed

| One variant | WSL, Linux disk | Windows, native |
| - | -: | -: |
| Typical, median of 299 | 1.9 ms | 1.9 ms |
| Slowest of 299 | 12.0 ms | 7.3 ms |

| 3,000 records, one process | WSL, Linux disk | Windows, native |
| - | -: | -: |
| To TSV | 11.4 s | 11.9 s |
| To Excel | 20.3 s | 15.5 s |
| To HTML report | 12.2 s | 10.4 s |
| To JSON Lines | 11.5 s | 10.6 s |

### A larger VCF: 30,000 records

The last row is what you get by default, since a VCF over 1 MB uses up to 4 processes and writes all three reports.

| Output and processes | WSL, Linux disk | Windows, native | Peak memory |
| - | -: | -: | -: |
| TSV, 1 process | 155.9 s | 112.8 s | About 50 MB on Windows |
| TSV, 2 processes | 137.2 s | Not measured | |
| TSV, 4 processes | 83.0 s | 59.8 s | 88 MB on Windows, main process |
| All three reports, 1 process | 276.3 s | Not measured | 472 MB |
| **All three reports, 4 processes** | **169.8 s** | **126.9 s** | 532 MB on WSL, 484 MB on Windows, main process |

This laptop has 2 physical cores, so 4 processes gave about 1.9 times the speed for TSV. Each worker also uses about 75 MB. On a small file, starting the processes uses up the gain, so VarSight uses one process under 1 MB.

The reports for those 30,000 records came to 35.8 MB for `report.html`, 12 to 13 MB for `results.xlsx` and 51.9 MB for `results.tsv`.

## Compared with published figures

These come from other people's papers, on other data and hardware. **They are not head-to-head results.** Tools were not re-run for this page.

### Agreement with ClinVar's names

| Tool | c. names | p. names | Source |
| - | -: | -: | - |
| ANNOVAR | 75.3% | 98.1% | 164,549 two-star ClinVar variants, [Human Genomics 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12181866/) |
| SnpEff | 98.8% | 98.6% | Same |
| VEP | 97.7% | 99.1% | Same |
| VarSight, SNVs | 99.51% | 99.38% | 20,027 names across hg19, [validation](validation.md) |
| VarSight, indels | 98.09% | 95.98% | Same |

VarSight's sample is much smaller than the study's, covers hg19 only and has no review-star filter, so these figures aren't directly comparable. Run `varsight validate` to measure your own.

### Time for a 30x genome on one core

| Tool | Time | Source |
| - | -: | - |
| Nirvana | 7 min | [Nirvana paper, ACM BCB 2017](https://doi.org/10.1145/3107411.3108204) |
| SnpEff | 15 min | Same |
| ANNOVAR | 67 min | Same |
| VEP | 18 h | Same |

VarSight is a Python program and isn't built to win this. At the TSV rates above with 4 processes, about 5 million records would take around 3 to 4 hours on this laptop, if every record were as slow as these gene-dense ones. Its strengths are single variants and exome-sized files on an ordinary laptop.

## Setup time

A complete hg19 setup ran on 2 November 2026, on the laptop above, natively on Windows with the data on a second drive, at about 2.3 MB/s. It took 44 minutes, including a restart that repeated about 3 of them, and ended with all 5 hg19 known variants passing. hg38 hasn't been timed, so its figures are estimates.

| Step | Time | Size kept | Notes |
| - | -: | -: | - |
| hg19 reference | 11 min | 3.0 GB | About 7 of them downloading |
| Genes, CIViC | 2 min | 15 MB | |
| UniProt, Pfam | 2 min | 20 MB | The old 4.7 GB Pfam read was on course for 12 to 14 hours |
| Transcripts | 6 min | 94 MB | 82,539 RefSeq and 196,483 GENCODE transcripts |
| ClinVar | 20 min | 1.6 GB | 5 downloading and 15 building, for 4.6 million records |
| **Total for one build** | **41 min** | **4.7 GB** | |

The build steps peaked at 316 MB of memory, measured in WSL.

### Data on a Windows drive, used from WSL

It works, and setting up from PowerShell and then reading the same folder from WSL is the practical way to do it. Building from WSL is the slow part. Measured on the same data on `I:`:

| Task | PowerShell | WSL through `/mnt/i` |
| - | -: | -: |
| First read of a place in the reference | 0.45 ms | 7 ms, or 54 ms with WSL's memory nearly full |
| The same place again | | 0.007 ms |
| Building one test index of 600,000 rows | 2.4 s | 65 to 96 s |
| One `varsight annotate` command, start to finish | | About 2 s, mostly loading the code through the same link |
| Copying the 3.2 GB reference into WSL | | 13 min |

## System requirements

| Need | Minimum | Notes |
| - | - | - |
| System | Ubuntu or similar Linux, WSL, or Windows 10 and later | PowerShell and Command Prompt both work |
| Python | 3.10 or newer | Nothing else to install |
| Memory | 2 GB free for a large VCF | Measured at under 150 MB per process for one report on 3,000 records, and about 500 MB plus 75 MB per worker for all three reports on 30,000. The HTML report keeps up to 20,000 rows in memory, so this should level off beyond that size. For a very large VCF on a tight machine, add `out=results.tsv`. Setup peaked at 316 MB |
| CPU | Any | More cores help on a large VCF, up to 4 by default |

Disk and download depend on what you choose. Download sizes were read from each server. Disk sizes were measured for hg19 and are estimates for hg38.

| You choose | Download | Kept on disk | Free during setup |
| - | - | -: | -: |
| hg19 | 1.7 GB | About 5 GB | About 15 GB |
| hg38 | 1.9 GB | About 5 GB | About 15 GB |
| Both | 3.0 GB | About 10 GB | About 25 GB |

Advanced databases add to this, from under 5 MB (Cancer Hotspots) to 59 GB read for gnomAD. See [data sources](data-sources.md#advanced-databases).
