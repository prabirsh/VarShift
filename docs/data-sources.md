# Data sources, licences and links

Every file comes directly from its original source to your computer. The list that controls this is [`src/varsight/data/sources.json`](../src/varsight/data/sources.json). This repository redistributes only a small demo subset.

Sizes and dates were read from each server on 1 November 2026, using header requests.

> [!WARNING]
> The licence notes below summarise what each source says. They aren't legal advice. **Check them against your own use**, particularly if results go into paid reports.

## Core databases

`varsight setup` installs these. You choose hg19, hg38 or both, and only the files for your choice are downloaded.

| Source | Used for | Download | Last changed | Update policy | Licence note |
| - | - | - | -: | - | - |
| hg19 reference (UCSC) | REF check, sequence | 904.8 MB | Aug 2018 | Once | UCSC states genome data is free for public and commercial use |
| hg38 reference (UCSC) | REF check, sequence | 938.1 MB | Jan 2014 | Once | Same |
| Chain files (Ensembl), both directions | Cross-build mapping | 1.0 MB | Jul 2014 | Once | Ensembl data; licence terms for these files not confirmed |
| RefSeq GRCh37 | hg19 transcripts | 28.8 MB | Mar 2022 | Your approval | NCBI places no restrictions on reuse; submitters may hold rights it can't clear |
| RefSeq GRCh38 | hg38 transcripts | 74.2 MB | Aug 2024 | Your approval | Same |
| GENCODE v19 | hg19 Ensembl transcripts | 44.3 MB | Feb 2016 | Your approval | GENCODE describes its data as open access |
| GENCODE latest basic | hg38 Ensembl transcripts | 97.4 MB | Jun 2026 | Your approval | Same |
| MANE summary | Default transcripts | 1.1 MB | Dec 2025 | Your approval | NCBI, as above |
| ClinVar VCF, one per build | Classification, conditions | 191.4 MB each | Weekly | Automatic | NCBI, as above |
| ClinVar variant_summary | HGVS names, p. lookup | 435.4 MB | Weekly | Automatic | NCBI, as above |
| HGNC complete set | Gene names, aliases, IDs | 16.2 MB | Often | Automatic | CC0, per HGNC |
| UniProt human reference proteome | Sequence, Pfam IDs, features | 126.6 MB | Monthly | Automatic | CC BY 4.0: credit UniProt |
| Pfam-A regions | Domain positions | 4.7 GB, read without saving | Jan 2026 | Automatic | CC0, per InterPro |
| CIViC | Cancer evidence | 4.6 MB | Nightly | Automatic | CC0, per CIViC |
| ClinGen validity and dosage | Gene-level curation | 1.3 MB | Often | Automatic | CC0, per ClinGen |
| gnomAD v2.1.1 constraint | pLI, LOEUF | 4.4 MB | Aug 2020 | Automatic | Licence terms not confirmed; see gnomad.broadinstitute.org |
| HPO genes to disease | Gene-disease links | 1.4 MB | Sep 2026 | Automatic | Free with citation of the HPO Consortium; see hpo.jax.org |

Downloads are deleted after each build. If you can't spare the Pfam read, the word `nopfam` on `varsight setup` skips it, but then domain positions are missing.

## Advanced databases

These are optional. A report gains columns only for the ones you install, and the core results never depend on them.

```bash
varsight databases              # list them all, with what is installed
varsight add dbsnp              # add one
varsight remove dbsnp           # remove it and free the disk space
```

### VarSight can fetch these

| Database | Adds | Transfer | Kept on disk (estimate) | Licence note |
| - | - | - | - | - |
| `dbsnp` | rs numbers, plus the frequencies dbSNP carries from TOPMed, 1000 Genomes, gnomAD, ExAC and ALFA | 26 GB per build, read without saving | 1 to 3 GB per build | NCBI places no restrictions on reuse |
| `topmed` | Comes with `dbsnp`, through its `FREQ` field | Included | Included | As dbSNP. TOPMed's own per-sample files need dbGaP access |
| `thousand_genomes` | Allele frequency overall and in five continental groups | 1.4 GB on hg19, 0.9 GB on hg38, read without saving | Under 200 MB per build | Open data from the International Genome Sample Resource |
| `gnomad` | Allele frequency, popmax and homozygote count from v2.1.1 exomes. hg38 uses gnomAD's own liftover | 59 GB on hg19, 85 GB on hg38, read without saving | 2 to 4 GB per build | Public release; check gnomad.broadinstitute.org |
| `gnomad_sv` | Structural variants overlapping the position | 134 MB on hg19, 1.6 GB on hg38 | Under 100 MB per build | Same |
| `phylop` | Conservation across 100 vertebrates | 8.9 GB on hg19, 9.2 GB on hg38 | The same, kept as the original file | UCSC states genome data is free for public and commercial use |
| `hotspots` | Whether the residue is a recurrent cancer hotspot, and in how many tumours | Under 5 MB | Under 5 MB | Public website; check cancerhotspots.org |
| `decipher` | Haploinsufficiency prediction and DDG2P entries for the gene | Under 5 MB | Under 10 MB | Public downloads; check deciphergenomics.org and Gene2Phenotype. Patient-level DECIPHER data needs a data access agreement, so that stays a link |
| `clingen_variants` | Expert panel classification and evidence codes | About 20 MB | About 30 MB | CC0, per ClinGen |

### You must accept the terms first

| Database | Adds | Transfer | Kept (estimate) | Licence note |
| - | - | - | - | - |
| `revel` | Missense pathogenicity score, 0 to 1 | About 650 MB | About 3 GB per build | Check the REVEL site before any commercial use |
| `alphamissense` | Missense prediction from Google DeepMind | About 600 MB per build | About 3 GB per build | The file header states CC BY-NC-SA 4.0, while a later Zenodo record says CC BY 4.0 |

```bash
varsight add revel accept
```

### These need your own account

VarSight can't download these for you, because their licences need your own registration. Download the file yourself, then point VarSight at it. The file stays on your computer.

| Database | How to get it | Licence note |
| - | - | - |
| `cosmic` | Register at cancer.sanger.ac.uk/cosmic/download, and download Cosmic_GenomeScreensMutant for your build | Free for non-commercial use after registering. Commercial users need a licence from QIAGEN, and paid patient reports count as commercial |
| `cosmic_census` | The same site, Cosmic_CancerGeneCensus | Same |
| `cosmic_fusions` | The same site, Cosmic_Fusion | Same |
| `omim` | Register at omim.org/downloads, and download genemap2.txt with your own key | Copyright Johns Hopkins University. Free for research after registering; commercial use needs a licence |
| `spliceai` | Sign in at basespace.illumina.com/s/otSPW8hnhaZR, and download spliceai_scores.masked.snv for your build | Free for academic and not-for-profit use; otherwise a licence from Illumina |
| `primateai` | Sign in at basespace.illumina.com/s/yYGFdGih1rXL, and download the PrimateAI scores | Check Illumina's terms for your use |

```bash
varsight add cosmic_census file=Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz accept
varsight add cosmic file=Cosmic_GenomeScreensMutant_v101_GRCh38.tsv.gz build=hg38 accept
```

If the file name shows the build, you don't need `build=`.

### How much of this has been tested

| Database | Status |
| - | - |
| Cancer Hotspots, DECIPHER, DDG2P, ClinGen classifications | Installed from the real, complete files and read back through annotation on Linux and Windows on 1 November 2026: 1,165 hotspots, 19,124 genes, 2,871 DDG2P entries and 12,207 classifications |
| dbSNP, 1000 Genomes, gnomAD, gnomAD SV, REVEL | Builders tested on samples cut from each real file. The full downloads have not been run |
| AlphaMissense | Builder tested on a small real-format sample |
| phyloP | The reader was run across the network against the real 9 GB files. It gives 9.236 at hg38 chr7:140753336, the same as dbNSFP's value for that site |
| COSMIC, Census, fusions, OMIM, SpliceAI, PrimateAI | **Written from each provider's documented columns, without a real file.** If your file doesn't match, VarSight stops with a clear message and installs nothing. Please report it |

### Scope for very large files

By default, dbSNP, 1000 Genomes, gnomAD and SpliceAI keep only coding exons plus 50 bases either side (500 for SpliceAI). A 26 GB file becomes one or two GB, and most interpretation work is there. To keep everything, add `scope=all`. Limit gnomAD to some chromosomes with `chrom=17,21`. Coverage outside these regions is then reported as `not_found`, so `not_found` for a deep intronic variant means only that it wasn't kept.

## Also available in ClinVar records

Frequencies listed in ClinVar records (ExAC, ESP, 1000 Genomes) are included with the core data, but only for variants ClinVar has.

## Not included

- **OncoKB:** full content needs a licence. Link only.
- **Franklin, VarSome, ProteinPaint:** websites without a free bulk download. Link only.
- **CADD, dbNSFP:** not added. Their terms need checking first, and many restrict commercial use.
- **hg18:** VarSight works with hg19 and hg38 only. hg18 input is refused rather than guessed at.
- **UCSC chain files:** licensed for non-commercial use by independent researchers and nonprofits only, so Ensembl's are used instead.

## Reference links

Each link was opened in an automated browser using BRAF V600E, on the date shown. A person may see something different, and sites change.

| Database | Label | Result | Tested |
| - | - | - | - |
| ClinVar | Confirmed ID | Opens; same classification | 31 Oct 2026 |
| NCBI Gene, RefSeq, dbSNP | Confirmed ID | Open | 31 Oct 2026 |
| HGNC | Confirmed ID | Opens | 1 Nov 2026 |
| Ensembl | Confirmed ID | Opens (GRCh37 site for hg19; the main site redirects for hg38) | 31 Oct 2026 |
| UniProt | Confirmed ID | Opens | 31 Oct 2026 |
| Pfam | Confirmed ID | Opens, hosted on InterPro | 31 Oct 2026 |
| Franklin | Built | Shows this exact variant on both builds | 31 Oct 2026 |
| OncoKB | Built | Shows V600E | 31 Oct 2026 |
| ProteinPaint | Built | **Gene view only**, not one variant | 31 Oct 2026 |
| COSMIC | Built | **Goes to the login page** | 31 Oct 2026 |
| VarSome | Built | **Not confirmed:** the test hit a security page | 31 Oct 2026 |
| gnomAD | Built | Opens on the right variant, on both builds | 1 Nov 2026 |
| ClinGen | Confirmed ID | Opens | 1 Nov 2026 |
| CIViC | Confirmed ID | Opens the V600E page | 1 Nov 2026 |
| cBioPortal | Built | Opens BRAF mutations | 1 Nov 2026 |
| PubMed | Built | Opens a search | 1 Nov 2026 |
| OMIM | Confirmed ID | **Not confirmed:** the test hit a security page | 1 Nov 2026 |
| UCSC Genome Browser | Built | **Not tested:** the test browser refused the site | 1 Nov 2026 |

A format that changes is a one-line fix in [`data/links.json`](../src/varsight/data/links.json).

## Citing

If you publish results, cite the sources you relied on:

- **ClinVar:** Landrum MJ et al., *Nucleic Acids Research*.
- **RefSeq:** O'Leary NA et al., *Nucleic Acids Research*.
- **GENCODE:** Frankish A et al., *Nucleic Acids Research*.
- **MANE:** Morales J et al., *Nature*.
- **HGNC:** Seal RL et al., *Nucleic Acids Research*.
- **UniProt:** The UniProt Consortium, *Nucleic Acids Research*.
- **InterPro and Pfam:** Blum M et al., *Nucleic Acids Research*.
- **CIViC:** Griffith M et al., *Nature Genetics*.
- **ClinGen:** Rehm HL et al., *New England Journal of Medicine*.
- **gnomAD:** Karczewski KJ et al., *Nature*.
- **HPO:** Gargano MA et al., *Nucleic Acids Research*.
- **dbSNP:** Sherry ST et al., *Nucleic Acids Research*.
- **1000 Genomes:** The 1000 Genomes Project Consortium, *Nature*.
- **REVEL:** Ioannidis NM et al., *American Journal of Human Genetics*.
- **AlphaMissense:** Cheng J et al., *Science*.
- **SpliceAI:** Jaganathan K et al., *Cell*.
- **PrimateAI:** Sundaram L et al., *Nature Genetics*.
- **phyloP:** Pollard KS et al., *Genome Research*.
- **COSMIC:** Tate JG et al., *Nucleic Acids Research*.
- **Cancer Hotspots:** Chang MT et al., *Cancer Discovery*.
- **DECIPHER:** Firth HV et al., *American Journal of Human Genetics*.

Check each source's site for its current preferred citation.
