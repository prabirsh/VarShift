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
| VarSight, SNVs | 99.98% | 99.99% | 12 genes only, [validation](validation.md) |
| VarSight, indels | 99.83% | 99.45% | 12 genes only |

VarSight's set is much smaller and easier, so expect its whole-genome figures to be lower. Run `varsight validate` to measure your own.

### Time for a 30x genome on one core

| Tool | Time | Source |
| - | -: | - |
| Nirvana | 7 min | [Nirvana paper, ACM BCB 2017](https://doi.org/10.1145/3107411.3108204) |
| SnpEff | 15 min | Same |
| ANNOVAR | 67 min | Same |
| VEP | 18 h | Same |

VarSight is a Python program and isn't built to win this. At the TSV rates above with 4 processes, about 5 million records would take around 3 to 4 hours on this laptop, if every record were as slow as these gene-dense ones. Its strengths are single variants and exome-sized files on an ordinary laptop.

## System requirements

| Need | Minimum | Notes |
| - | - | - |
| System | Ubuntu or similar Linux, WSL, or Windows 10 and later | PowerShell and Command Prompt both work |
| Python | 3.10 or newer | Nothing else to install |
| Memory | 2 GB free for a large VCF | Measured at under 150 MB per process for one report on 3,000 records, and about 500 MB plus 75 MB per worker for all three reports on 30,000. The HTML report keeps up to 20,000 rows in memory, so this should level off beyond that size. For a very large VCF on a tight machine, add `out=results.tsv`. Full setup memory hasn't been measured |
| CPU | Any | More cores help on a large VCF, up to 4 by default |

Disk and download depend on what you choose. Download sizes were read from each server, while disk sizes are estimates until a full setup has been run.

| You choose | Download | Kept on disk | Free during setup |
| - | - | -: | -: |
| hg19 | 2.3 GB, plus 4.7 GB of Pfam read without saving | About 5 GB | About 15 GB |
| hg38 | 2.4 GB, plus 4.7 GB of Pfam read without saving | About 5 GB | About 15 GB |
| Both | 3.3 GB, plus 4.7 GB of Pfam read without saving | About 10 GB | About 25 GB |

Advanced databases add to this, from under 5 MB (Cancer Hotspots) to 59 GB read for gnomAD. See [data sources](data-sources.md#advanced-databases).
