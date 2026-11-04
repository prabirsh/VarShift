# What people struggle with, and what VarShift does about it

This page records the research behind VarShift's design: what the published evidence says, which tools already exist, and why the reports look the way they do. It was done in November 2026. Other tools were not re-run for this page, so what it says about them comes from their papers and documentation.

## The problems

| Problem | What the evidence says | What VarShift does | Where it stops |
| - | - | - | - |
| The same variant gets different names | 164,549 two-star ClinVar variants run through ANNOVAR, SnpEff and VEP: all three agreed on the HGVS c. name for only 58.52%, on the p. name for 84.04%, and on coding impact for 85.58% ([Human Genomics, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12181866/)). In one exercise, 20 laboratories described one variant 14 different ways (Deans et al., 2016, cited by [Yen et al., 2017](https://doi.org/10.1186/s13073-016-0396-7)) | Measures its own names against ClinVar's, shows ClinVar's names beside its own, and lists every candidate for a p. name | Agreement is measured on a 20,027-name hg19 sample so far, not yet on hg38 ([validation](validation.md)) |
| Indels are where tools differ most | For SNVs, VEP and SnpEff matched ClinVar's syntax over 99.5% of the time. For other variant types it was under 90% (Yen et al., 2017) | Reports SNVs and indels separately, so an average can't hide this | |
| The transcript set changes the answer | Loss-of-function agreement between RefSeq and Ensembl annotations in ANNOVAR was 44%, and between ANNOVAR and VEP it was 65% ([McCarthy et al., 2014](https://doi.org/10.1186/gm543)) | Uses MANE Select first, says which transcript it used, and shows the others | GENCODE is used only where RefSeq has nothing, or when asked |
| Checking one variant means many websites | Each database has its own search and its own naming, so a person repeats the search, copying between tabs | One result with a labelled link per database | Links to COSMIC, OncoKB, Franklin, VarSome and ProteinPaint open those sites. Their data isn't stored locally |
| Patient data leaves the building | A web search sends the variant to that site | Annotation has no network code at all, and a test enforces it | Clicking a link does send that variant to that site |
| hg19 and hg38 get mixed up | A VCF on the wrong build can look fine and give wrong results | Reads the build from the VCF's contig lengths, stops on a contradiction, and always shows both positions when both are installed | |
| Long lists with no starting point | A VCF gives thousands of rows with nothing to say which to read first | Sorts results by review priority, each with its reason | Priority is a sorting aid, not a classification |
| Installing the tools | Common annotators need a language runtime plus a large per-build cache | Needs only Python and works on Ubuntu, WSL, PowerShell and Command Prompt | |
| Not knowing the data is old | Local databases go stale quietly | Tells you when a newer version exists, without slowing the run | Depends on the source reporting a date or version |

**A note on time saved.** Curation time is often quoted, but the evidence is thin. One clinical trial protocol says "hours of skilled person time per variant" ([NCT04196374](https://clinicaltrials.gov/study/NCT04196374)), and most other figures are vendor claims. VarShift therefore makes no claim about hours saved, and you can measure it on your own work.

## Other tools

Each of these is good at what it was built for. This table is about fit for one person on an ordinary laptop, and it does not rank them.

| Tool | Runs on | Input | Output | Notes |
| - | - | - | - | - |
| [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) | Perl, with a large cache per build, or Ensembl's web service | VCF, HGVS and more | VCF, tab, JSON | The widely used reference point. Very broad |
| [ANNOVAR](https://annovar.openbioinformatics.org/) | Perl | Its own format, or converted VCF | Tables | Free for non-commercial use, after registering |
| [SnpEff](https://pcingola.github.io/SnpEff/) | Java | VCF, BED | VCF, HTML summary | Fast, with its own databases |
| [Nirvana](https://illumina.github.io/NirvanaDocumentation/) | .NET | VCF | JSON | Fast, and used in DRAGEN workflows |
| [OpenCRAVAT](https://docs.opencravat.org/en/latest/3.-Viewing-Results.html) | Python, with many modules | VCF and others | A local web interface with Summary, Filter, Variant and Gene tabs, plus Excel, TSV, VCF | The closest in spirit. `oc gui` runs a local server |
| [PCGR](https://pmc.ncbi.nlm.nih.gov/articles/PMC5946881/) | Conda or Docker, with a large data bundle | Tumour VCF | Tiered interactive HTML | Built for cancer genomes |
| Franklin, VarSome | Website | One variant, or uploads | Web pages | The variant goes to their servers, with interpretation support |

**Where VarShift fits:** a person with a p. name, a c. name or a VCF, a laptop, and patient data that shouldn't leave it. It doesn't replace VEP for pipelines, or an interpretation platform for classification.

## How results should be shown

The brief asked for Excel plus something smarter. These were the options:

| Option | Good | Bad | Decision |
| - | - | - | - |
| Excel workbook | Everyone can open it, filter it and send it | Awkward for browsing one variant in detail | **Always produced.** Sheets: Summary, Priority, All results, Warnings explained, Data versions |
| A single HTML file | Interactive. Opens with a double click, works offline, can be emailed or kept with a case | Very large VCFs need a row limit | **Always produced.** Up to 20,000 rows. When there are more, Low results leave first and the page says so |
| TSV | Easy to use from scripts, R or pandas | Not much to look at | **Always produced** |
| A local web server, as OpenCRAVAT does | Handles any size, with rich queries | Something has to keep running, ports and firewalls get in the way, and the result can't be emailed | Not chosen |
| PDF | Good for filing | Can't be filtered, and a VCF gives thousands of pages | Not chosen. The HTML report prints cleanly if one case is needed |
| Web app in the cloud | Nothing to install | Patient data would leave the computer | Ruled out |

Single variants print a readable card in the terminal, because a report folder for one variant would just get in the way. JSON and JSON Lines are available for programs.

### What the interactive report has

- Four cards counting High, Medium, Review and Low results. Clicking one filters the table.
- Bars for effect, genes and ClinVar classes. Clicking one filters the table.
- A search box that checks every field, so `BRAF V600E`, an rs number or a chromosome all work.
- A sortable table, with a panel for each result: why it got its priority, both positions, every name with a copy button, ClinVar, protein domain, gene evidence, advanced databases, warnings explained in plain words, other transcripts, and each link labelled as a confirmed ID or a built address.
- A button that saves the filtered rows as TSV.
- Keyboard use throughout, dark mode, phone width, and print styles.

## Why priority works the way it does

Readers need a starting point, but an automatic pathogenic or benign call would go beyond the evidence and could be mistaken for a classification. So:

- It only sorts. It never says a variant is pathogenic or benign.
- Every result states its reasons, so you can disagree with them.
- Score thresholds are the published calibrated ones: REVEL 0.644 ([Pejaver et al., 2022](https://doi.org/10.1016/j.ajhg.2022.10.013)) and SpliceAI 0.2 ([Walker et al., 2023](https://doi.org/10.1016/j.ajhg.2023.06.002)).
- A variant reported as pathogenic stays High even if common. The common frequency is added to its reasons, so that conflict stays visible.

The rule itself is about a hundred lines in [`priority.py`](../src/varshift/priority.py).
