# Demo

Run VarSight end to end with no download:

```bash
varsight demo
```

It builds a small data set from `sources/`, using exactly the same code as a real setup, then annotates:

- [`variants.txt`](variants.txt): p., c., g. and rs inputs.
- [`demo.hg19.vcf`](demo.hg19.vcf): a VCF. It includes a two-allele line, a wrong REF and a symbolic allele on purpose, to show how each is handled.

Results go to `demo_output/`, in one folder for each input, each holding `report.html`, `results.xlsx` and `results.tsv`. **Open `demo_output/variants/report.html` in a browser first.**

## What is in `sources/`

Real public data, cut down to 12 genes: BRAF, KRAS, NRAS, EGFR, TP53, PIK3CA, IDH1, BRCA1, BRCA2, JAK2, KIT and PTEN.

| File | From | Cut down to |
| - | - | - |
| `hg19.windows.fa.gz`, `hg38.windows.fa.gz` | UCSC | The sequence around each gene (6 kb either side), at its real coordinates, plus chrM |
| `refseq.*.gff.gz` | RefSeq | Lines for these genes |
| `clinvar.*.vcf.gz` | ClinVar 2026-10-25 | At most 1,200 records per gene, plus the variants used in the demo and tests |
| `variant_summary.txt.gz` | ClinVar | Names for those records, from the start of the file only |
| `mane.summary.txt.gz`, `hgnc.txt.gz` | MANE, HGNC | These genes |
| `clingen_*.gz`, `gnomad_constraint.tsv.gz`, `hpo_diseases.txt.gz` | ClinGen, gnomAD, HPO | These genes |
| `uniprot.dat.gz`, `pfam_regions.tsv.gz` | UniProt, InterPro | These proteins |
| `civic_*.tsv.gz` | CIViC | These genes, without evidence statement text |
| `chain.*.gz` | Ensembl | These chromosomes |

[`tools/build_demo_data.py`](../tools/build_demo_data.py) rebuilds them. [`sources/provenance.json`](sources/provenance.json) records the date and versions.

**Nothing here comes from a patient.** The depth and genotype values in the demo VCF are invented.

**This data is only for trying VarSight out.** Most of the genome is missing, so a variant elsewhere reports `not covered by the installed reference`, and ClinVar is deliberately incomplete. A real `not in ClinVar` result needs the full data.

## Licences

These subsets keep the terms of their sources: see [docs/data-sources.md](../docs/data-sources.md). The repository's MIT licence covers the code, not this data.
