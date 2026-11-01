# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-01

First release.

### Added

- Annotation of p., c., g., rs, list and VCF input on hg19 and hg38.
- Normalisation, cross-build mapping, HGVS c. and p. names, and Sequence Ontology effects.
- Evidence from ClinVar, CIViC, ClinGen, gnomAD constraint and HPO, and protein context from UniProt and Pfam.
- Labelled reference links for 20 databases.
- JSON, JSON Lines, TSV and Excel output.
- Automated setup, weekly update with approval for results-changing sources, checksums and rollback.
- Validation against ClinVar's names, a comparison tool for annotated VCFs, and a bundled demo.
- Optional AlphaMissense and gnomAD frequency sources.

### Known limits

- Accuracy has been measured on 12 genes only. See [docs/validation.md](docs/validation.md).
- The full setup has not yet been run end to end.
- CNVs, structural variants and effects on mitochondrial transcripts are not covered.
