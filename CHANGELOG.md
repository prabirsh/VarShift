# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

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
