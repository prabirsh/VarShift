# Data sources, licences and links

Every file comes directly from its original source to your computer during `varsight setup`. The list that controls this is [`src/varsight/data/sources.json`](../src/varsight/data/sources.json).

Sizes and dates were read from each server on 1 November 2026, using header requests.

> [!WARNING]
> The licence notes below summarise what each source says. They aren't legal advice. **Check them against your own use**, particularly if results go into paid reports.

## Installed by default

| Source | Used for | Download | Last changed | Update policy | Licence note |
|---|---|---|---:|---|---|
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

**Total:** about 3.3 GB saved and processed, plus 4.7 GB of Pfam read through, so about 8 GB transferred. Downloads are deleted after each build. If you can't spare that, `--skip-pfam-regions` avoids the Pfam read, but then domain positions are missing.

**Final size:** about 10 to 12 GB (estimate). Keep about 25 GB free during setup.

## Added beyond the original list

These weren't in the original request, and were added because they help interpretation and are free to download:

- **ClinGen** gene-disease validity and dosage sensitivity.
- **gnomAD constraint** (pLI, LOEUF), showing how tolerant a gene is of damaging variants.
- **HPO** gene-disease links.
- **CIViC** cancer evidence.
- **MANE**, for choosing the default transcript.
- **GENCODE**, so `ENST` input works.
- **Links** to gnomAD, UCSC Genome Browser, ClinGen, CIViC, OMIM, cBioPortal and PubMed.

## Optional

| Source | Command | Size | Note |
|---|---|---|---|
| AlphaMissense | `varsight setup --with alphamissense --accept-license alphamissense` | About 600 MB per build, and a large database | The file header states **CC BY-NC-SA 4.0, non-commercial**, while a later Zenodo record says CC BY 4.0. You must accept this yourself. Its builder has only been tested on a small real-format sample |
| gnomAD frequencies, hg19 | `varsight setup --with gnomad --gnomad-chrom 21` | About 59 GB read for all chromosomes, without saving the VCF. chr21 alone is 653 MB | Limit it with `--gnomad-chrom`. Its builder has only been tested on a small sample |

Frequencies listed in ClinVar records (ExAC, ESP, 1000 Genomes) are included by default, but only for variants ClinVar has.

## Not included

- **COSMIC:** free for non-commercial use after registering. Its terms count paid patient reports as commercial. Link only.
- **OncoKB:** full content needs a licence. Link only.
- **Franklin, VarSome, ProteinPaint:** websites without a free bulk download. Link only.
- **REVEL, CADD, SpliceAI, dbNSFP:** not added. Their terms need checking first, and many restrict commercial use.
- **UCSC chain files:** licensed for non-commercial use by independent researchers and nonprofits only, so Ensembl's are used instead.

## Reference links

Each link was opened in an automated browser using BRAF V600E, on the date shown. A person may see something different, and sites change.

| Database | Label | Result | Tested |
|---|---|---|---|
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
- **AlphaMissense:** Cheng J et al., *Science*.

Check each source's site for its current preferred citation.
