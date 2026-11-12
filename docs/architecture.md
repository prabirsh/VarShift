# Architecture

VarShift is a small Python package with no runtime dependencies. Each module has one job and can be read on its own.

## Modules

| Module | Job |
| - | - |
| `inputs.py` | Decide what kind of description a piece of text is |
| `tidy.py` | Fix slips that have only one sensible reading, such as `braf val600glu`, and say what changed |
| `suggest.py` | When a description fails, offer what it probably meant, from the installed data |
| `hgvs.py` | Parse c., g. and p. change descriptions |
| `builds.py` | Contig names, lengths and accessions for hg19 and hg38; work out a VCF's build |
| `genome.py` | Read reference bases through a memory map |
| `variant.py` | Check REF, trim and left-align |
| `transcripts.py` | Transcript model and the arithmetic between genome and c. or n. positions |
| `resolver.py` | Turn c., p. and rs input into normalised variants, keeping every candidate |
| `liftover.py` | Map a variant between builds and repeat the REF check |
| `consequence.py` | Effect, exon, and HGVS c. and p. names for one transcript |
| `store.py` | Read-only access to the local databases |
| `bigwig.py` | Read single values from a bigWig file, for phyloP |
| `annotate.py` | Run the steps for one variant and assemble the result |
| `priority.py` | Sort results for review, with reasons |
| `answer.py` | Sort one description's candidates into one answer (single, clear or choice), and lay out the card, the choice table and the summary line |
| `advanced.py` | What each advanced database is, how to get it, and which columns it adds |
| `links.py`, `data/links.json` | Build labelled reference links from templates |
| `vcf.py`, `batch.py` | Stream a VCF or a list, on several processes for a large file |
| `output.py`, `xlsx.py`, `report.py`, `data/report_template.html` | Write the screen card, Excel, the HTML report, TSV and JSON |
| `flags.py` | Each warning in plain words, with what to do |
| `commandline.py`, `cli.py` | Read the command line, and run each command |
| `term.py` | What kind of terminal this is: colour, width, and whether a person is there to answer. Nothing waits for an answer otherwise |
| `session.py` | The guided start: queries, choices, and saving the answers seen |
| `examples.py` | Every accepted way to write a variant, each with an example that a test annotates |
| `config.py` | Saved settings, and where the data folder is |
| `updates.py` | The quick notice when a newer database exists. **The only module that starts a process of its own for this check**, and it takes no input about variants |
| `explain.py` | Print each step, file and query for one variant |
| `pipeline/fetch.py` | Downloads. **The only module that opens a network connection** |
| `pipeline/builders.py`, `pipeline/advanced_builders.py` | Turn source files into databases |
| `pipeline/updatecheck.py` | The background half of the update notice |
| `pipeline/benchmark.py` | Measure speed and memory on this computer |
| `pipeline/manager.py` | First setup, updates, rollback, checksums |
| `pipeline/selftest.py`, `data/known_variants.json` | Known variants, checked after setup and before each update |
| `pipeline/validate.py` | Measure accuracy against ClinVar's names |
| `pipeline/demo.py` | Build and run the demo with the same builders as a real setup |

## Data folder

```text
varshift-data/
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
    dbsnp.hg19.sqlite, revel.hg38.sqlite, ...   only the advanced databases you added
  scores/                phyloP, if you added it, as the original bigWig files
  update_status.json     what the last update check found
  snapshots/             the previous version of anything updated
  reports/               a report of each setup and update
```

One file per source and build means a weekly update replaces only what changed. A process already running keeps the file it opened, so an update never disturbs it.

## Design decisions

**Why own the HGVS code instead of wrapping VEP?** VEP needs a Perl stack and a large cache per build, which was a poor fit for a laptop with 4 GB of RAM. The accuracy risk of writing this ourselves is real, so `validate` measures it against ClinVar on every change.

**Why SQLite?** It comes with Python and needs no server. On a WSL Linux disk, a lookup measured 27 microseconds. Files open read-only, so nothing can change the data during a lookup. [How lookups work](how-lookups-work.md) goes through it step by step.

**Why only standard-library code?** A tool that handles patient data benefits from having nothing to install and no dependency chain to trust. It also makes `git clone` and one install script all that's needed, on Linux and Windows alike.

**Why keep every candidate?** Picking one DNA change for a protein name would hide real ambiguity, so candidates are ranked and all shown.

**Why word options instead of dashed ones?** `build=hg19 out=results.xlsx` reads the same on every system and needs no explaining, and the same style works in Bash, PowerShell and Command Prompt. Anyone who types dashes out of habit is still understood.

**Why spawn for worker processes?** A forked child would inherit open SQLite connections, which SQLite warns against. Spawned workers open their own, on every system.

**Why an interactive HTML file rather than a local server?** It needs nothing running, and can be emailed or filed with a case. [The options that were weighed](research.md#how-results-should-be-shown).

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
