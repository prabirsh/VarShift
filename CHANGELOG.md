# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- The first whole-genome accuracy measure, on 20,027 ClinVar names across hg19, and the first complete setup from start to self-test. See [docs/validation.md](docs/validation.md).
- The warning `transcript_gap_here`, for a change that touches bases the transcript and the genome don't share.

### Changed

- Setup takes domain positions from Pfam's 2.8 MB human file, release 37.0, instead of reading the 4.7 GB file for every species, which was on course to take 12 to 14 hours on a home connection. The cost is listed in [data sources](docs/data-sources.md). The `nopfam` word is gone, since there is nothing left to skip.
- Long downloads now print how far they have got, their speed and the time left, and long builds say once a minute that they are still going.
- A running setup remembers its data folder, so running it again after an interruption offers the same one.

### Fixed

- `.gitignore` re-included everything under `demo/`, so a sample VCF dropped there to try VarSight, with its results, could have been committed. Only the demo's own files are kept now.
- `varsight add revel` failed with `403 Forbidden`, because REVEL's own server refuses scripted downloads. REVEL now comes from its Zenodo record, which its download page names as the other source, and is checked against Zenodo's MD5.
- c. and p. names on RefSeq transcripts that don't align exactly to the genome, 1,231 on hg19 and 306 of them MANE Select. Counting along RefSeq's exon rows put everything after a gap off by its length, and a c. name as input could land a few bases away. VarSight now numbers through RefSeq's own alignment. On 20,027 ClinVar names, c. SNVs went from 99.51% to 99.94% and p. SNVs from 99.38% to 99.81%. The transcript databases need rebuilding to get this.
- Setup, the demo and the benchmark no longer write to the system disk. SQLite's sort files, which can reach gigabytes while an index is built, and all temporary folders now stay beside the data or the output.
- UCSC now redirects its downloads to its Amazon bucket, `genome-browser.s3.amazonaws.com`, so setup stopped at the first step. That one bucket is now on the allow-list. Redirects still go through the same check, and the reference is still compared with UCSC's published MD5.
- Setup from WSL now refuses a data folder on a Windows drive unless you add `windowsdrive`. It was slow, and with memory tight it stopped with `Cannot allocate memory`. That error now says what to do.
- On Linux and WSL, `data=I:\folder` is now refused with the right spelling. Bash removes the backslash, so it used to create a new folder called `I:folder` in the current directory and download everything again.
- On Windows, setup could stop at the ClinGen servers with `CERTIFICATE_VERIFY_FAILED`, because Windows adds a trusted root only when its own client first needs it.
- A file with a matching MD5 is no longer unpacked once more just to check it.
- ClinVar records are read about 38% faster, with byte-identical results, because one changed base no longer goes through the full normaliser.

## [0.1.0] - 2026-11-01

First release.

### Added

- Annotation of p., c., g., rs, list and VCF input on hg19 and hg38, on Ubuntu, WSL, PowerShell and Command Prompt.
- A setup that asks whether you want hg19, hg38 or both, and installs only that.
- Normalisation, cross-build mapping, HGVS c. and p. names, and Sequence Ontology effects.
- Evidence from ClinVar, CIViC, ClinGen, gnomAD constraint and HPO, and protein context from UniProt and Pfam.
- Labelled reference links for 20 databases.
- Reports: an interactive HTML file, an Excel workbook with a priority sheet and explained warnings, and TSV. A readable card on the screen for one variant, and JSON for programs.
- Review priority, with reasons, as a sorting aid.
- Advanced databases, added one at a time: dbSNP with TOPMed, gnomAD, 1000 Genomes, gnomAD SV, REVEL, AlphaMissense, SpliceAI, PrimateAI, phyloP, Cancer Hotspots, DECIPHER, ClinGen classifications, COSMIC, Cancer Gene Census, gene fusions and OMIM.
- A quick notice when a newer database version exists, which never slows annotation.
- `varsight explain` to show each step, file and query, and `varsight benchmark` to measure a computer.
- Several processes for large VCFs, with identical output.
- Automated setup, weekly update with approval for results-changing sources, checksums and rollback.
- Validation against ClinVar's names, a comparison tool for annotated VCFs, and a bundled demo.

### Known limits

- Accuracy has been measured on 12 genes only. See [docs/validation.md](docs/validation.md).
- The full setup has not yet been run end to end.
- CNVs, structural variants and effects on mitochondrial transcripts are not covered.
