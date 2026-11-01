# Architecture

VarSight is a small Python package with no runtime dependencies. Each module has one job and can be read on its own.

## Modules

| Module | Job |
|---|---|
| `inputs.py` | Decide what kind of description a piece of text is |
| `hgvs.py` | Parse c., g. and p. change descriptions |
| `builds.py` | Contig names, lengths and accessions for hg19 and hg38; work out a VCF's build |
| `genome.py` | Read reference bases through a memory map |
| `variant.py` | Check REF, trim and left-align |
| `transcripts.py` | Transcript model and the arithmetic between genome and c. or n. positions |
| `resolver.py` | Turn c., p. and rs input into normalised variants, keeping every candidate |
| `liftover.py` | Map a variant between builds and repeat the REF check |
| `consequence.py` | Effect, exon, and HGVS c. and p. names for one transcript |
| `store.py` | Read-only access to the local databases |
| `annotate.py` | Run the steps for one variant and assemble the result |
| `links.py`, `data/links.json` | Build labelled reference links from templates |
| `vcf.py`, `batch.py` | Stream a VCF or a list |
| `output.py`, `xlsx.py` | Write JSON, TSV or Excel |
| `pipeline/fetch.py` | Downloads. **The only module that opens a network connection** |
| `pipeline/builders.py` | Turn source files into databases |
| `pipeline/manager.py` | First setup, updates, rollback, checksums |
| `pipeline/selftest.py`, `data/known_variants.json` | Known variants, checked after setup and before each update |
| `pipeline/validate.py` | Measure accuracy against ClinVar's names |
| `pipeline/demo.py` | Build and run the demo with the same builders as a real setup |

## Data folder

```text
varsight-data/
  manifest.json          what is installed: versions, sizes, SHA-256 of each file, sources, licences
  genome/hg19.fa, .fai   uncompressed, for fast random access
  genome/hg38.fa, .fai
  db/
    transcripts.hg19.sqlite, transcripts.hg38.sqlite
    clinvar.hg19.sqlite, clinvar.hg38.sqlite
    liftover.hg19_hg38.sqlite, liftover.hg38_hg19.sqlite
    genes.sqlite         HGNC, ClinGen, gnomAD constraint, HPO
    proteins.sqlite      UniProt, Pfam regions
    civic.sqlite
  snapshots/             the previous version of anything updated
  reports/               a report of each setup and update
  cache/                 the human rows kept from Pfam, so an update need not re-read 5 GB
```

One file per source and build means a weekly update replaces only what changed. A process already running keeps the file it opened, so an update never disturbs it.

## Design decisions

**Why own the HGVS code instead of wrapping VEP?** VEP needs a Perl stack and a large cache per build, and was a poor fit for a laptop with 4 GB of RAM. It also doesn't take a protein change as input on its own. The accuracy risk of writing this ourselves is real, so `validate` measures it against ClinVar on every change.

**Why SQLite?** It comes with Python and needs no server. On a WSL Linux disk, a lookup measured 27 microseconds. Files open read-only, so nothing can change the data during a lookup.

**Why only standard-library code?** A tool that handles patient data benefits from having nothing to install and no dependency chain to trust. It also makes `git clone` and `./install.sh` all that's needed.

**Why keep every candidate?** Picking one DNA change for a protein name would hide real ambiguity, so candidates are ranked and all shown.

**Why no merged verdict?** BRAF V600E shows why: ClinVar's aggregate says "Conflicting classifications of pathogenicity", while cancer sources call it oncogenic. A single merged answer would hide that disagreement.

**Why does the 3' rule cross exon-intron junctions?** To match how ClinVar names variants. The consequence check looks at both ends of a repeat, so a deletion that could sit on a splice site reports that, on either strand.

**Why Ensembl chain files?** UCSC's are licensed for non-commercial use by independent researchers and nonprofits only. Ensembl's gave 100% agreement with ClinVar's hg38 positions on the test set.

## Memory

- VCFs and source files are read a line at a time.
- Reference bases come through a memory map, so only the pages read are loaded.
- The SQLite page cache is capped at 16 MB per database when reading.
- The sequence cache holds 128 transcripts.

## Limits by design

- Only primary-assembly contigs get transcripts. Alt loci and patches are not annotated.
- Only curated RefSeq (`NM_`, `NR_`) and GENCODE transcripts. Model transcripts (`XM_`, `XR_`) are left out.
- Transcript sequence is taken from the genome. Where RefSeq notes that a transcript differs from the genome, the result is flagged.
