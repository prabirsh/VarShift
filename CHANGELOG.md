# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- The first whole-genome accuracy measure, on 20,027 ClinVar names across hg19, and the first complete setup from start to self-test. See [docs/validation.md](docs/validation.md).
- The warning `transcript_gap_here`, for a change that touches bases the transcript and the genome don't share.
- Slips with only one sensible reading are fixed, and the reading is recorded. `braf val600glu`, `BRAF c.1799T-A`, `NM_004333.6:1799T>A`, `chr7:140,453,136 A>T`, `BRAF(V600E)`, and text pasted with smart quotes or a fullwidth `>`, now annotate instead of stopping with an error. The input stays as typed, and the new `read_as` field and column hold the text VarShift read, with the gene symbol as HGNC writes it. A gene alias is named there, as in `KRAS G12D (KRAS2 is an alias)`. Text that already reads is never changed.
- Several variants in one command, separated by `;`, as in `varshift annotate "BRAF V600E; KRAS G12D" out=results.xlsx`. They are handled like the lines of a list file, so one that fails doesn't stop the others. A `;` inside square brackets, as in `c.[76A>C;83G>C]`, is left alone.
- Seven columns after `data_versions`: `read_as`, `answer`, `best`, `why`, `also`, `suggestion` and `mode`. Existing columns keep their positions. `suggestion` and `mode` are filled in by later changes.
- Each query now gets one answer instead of a list of equal candidates. `answer` says `single` for one candidate, `clear` when one stands out, and `choice` when nothing decides between them. A candidate stands out when it is the only one in ClinVar, or, with none in ClinVar, the only one that changes the fewest bases. It comes first, with `best` set to yes, `why` saying which rule decided it (and MANE Select when its transcript is), and `also` naming the others, as in `c.1799_1800delinsAA (2 bases, not in ClinVar)`. Each allele of an rs number is a different variant, so an rs number with several is always a `choice`.
- `varshift annotate` on the screen now shows that answer in one short card instead of a full card per candidate: BRAF V600E fits in 20 lines, with the best match, `Why this`, the other candidates on one `Also` line, and five links to check it. Warnings are written as sentences. A `choice` is shown as a numbered table, and in a terminal you pick one, several with `a`, or the first with Enter. Piped or scheduled, it never waits: it prints the table and how to see every candidate. `show=all` gives the old full card for every candidate, and `links=all` every link with its label and note.
- `mode=tumour` and `mode=inherited` on `annotate` choose which evidence and links the card shows first: CIViC, OncoKB and cBioPortal for tumour work, ClinGen, gnomAD constraint and OMIM for inherited disease. The mode is recorded in the new `mode` column for every input, VCFs included. `mode=tumor` is read as `tumour`.
- `varshift config mode=tumour` (or `inherited`) saves the mode, so you don't type it each time, and `mode=none` clears it. `varshift config color=no` (or `yes`, `auto`) saves the colour choice. `mode=` and `color=` on a command still win, and `varshift config` shows both. A saved value that isn't valid says how to clear it. `explain` takes `mode=` too, and shows it.
- Colour on the screen: priority words are coloured, as are the field names and questions. `color=auto` (the default) colours a terminal only, `color=yes` and `color=no` decide, and `NO_COLOR` always turns it off. On Windows the console is switched into colour mode if it can be.
- Several queries separated by `;` with no `out=` now print one numbered line each, then how to see one in full and how to save them all.
- The demo shows the new card.
- The guided start: `varshift` with nothing after it, in a terminal, or `varshift start`, opens a session. You type variants as a colleague wrote them, one or several separated by `;` or pasted together, and each gets its answer: one in full, several a numbered line each. Slips are asked about and choices are offered, as with `annotate`. The data stays open, so each answer after the first is immediate. `save` writes the answers you saw to `varshift-session-DATE-TIME.xlsx`, never over an earlier file, and `save NAME.tsv`, `.html`, `.json` or a folder writes those instead. A number shows that result again, and `all N` shows every candidate. The other words are `list`, `remove N`, `clear`, `mode tumour`, `examples` and `help`. Quitting with results not saved asks first. With no data set up, it offers the 12-gene demo data at once, with no download, or setup. Anything saved from the demo says so. Piped or scheduled, `varshift` alone prints the command list as before, and `varshift start` says it needs a terminal: neither ever waits for an answer.
- A typo that fits a transcript other than the default is no longer passed over in silence. `BRAF A600E` fails on NM_004333.6, the MANE Select transcript, which has Val at residue 600, but fits NM_001378475.1, which has Ala there, so it used to give a different variant with only `other_transcript_used`. Now, at a terminal, VarShift says what the default transcript has and asks whether you meant `BRAF V600E` or the description as typed. Elsewhere the answer stands, with the warning `default_transcript_differs`, a note on the card, and `BRAF V600E` in the `suggestion` column. The same goes for c. names, such as `BRAF c.1799C>G`.
- `varshift examples` shows every way to write a variant, with an example of each that works: protein changes in one- and three-letter form, `NP_` names, c. with a transcript or a gene, n., g. with `NC_` or `chr`, positions in three spellings, m., rs numbers, list files, VCFs and `;` lists. It needs no data. A test annotates every example on the demo data, so none can quietly stop working. When nothing can be suggested for an input, the error now points there.
- When a description fails, VarShift offers what it probably meant, built from the installed data. `BRAF V600` lists the changes ClinVar and CIViC know at that residue, `V600E` the genes that have it, `BRFA V600E` the symbols one or two letters away (BRAF first, ranked by ClinVar records), `BRAF W600E` the residue the protein really has, `BRAF c.1799A>T` the swapped or corrected bases, and a position whose REF fits the other build that build. In a terminal you choose a number and it is annotated, with `Read as` saying what was chosen and the input kept as typed. Otherwise, as with `out=` or `format=`, the error is followed by `Did you mean:` and the exit code stays 1. Nothing is ever applied without being chosen.
- In a list file or a `;` list, a line that fails stays an error row, and its suggestions go in the `suggestion` column, so the file can be fixed and run again. Only the quick suggestions run there (not the search for which genes have a change), and only for the first 500 failing lines, so a broken file can't run for hours. The summary says how many rows have a suggestion. VCF records are never tidied or suggested: they are written by machines, and a change there would hide a real problem.
- Every workbook has a Links sheet, after All results: one row per result and database, with the result number, input, gene, c. and p. names, the database, a clickable link, whether it is a confirmed ID or built from the variant, and its note. A database with no link for that result says why. The rows wait in a private file beside the output while the run goes on, so a large VCF doesn't fill memory. Past Excel's limits, the Summary sheet says what was left out and where to find it.

### Changed

- The tool is renamed from VarSight to VarShift, because another project already uses the name VarSight. *Var* is for variant, and *Shift* is for the move from raw variants to annotation and interpretation. The command, the Python package, the `VARSHIFT_DATA` variable, the default `~/varshift-data` folder and the `_varshift` output folder all use the new name. There is no alias for the old one. Data installed under the old name still works, but the saved setting that remembered its folder is not carried over, so point at it with `data=PATH` or `VARSHIFT_DATA`.
- Setup takes domain positions from Pfam's 2.8 MB human file, release 37.0, instead of reading the 4.7 GB file for every species, which was on course to take 12 to 14 hours on a home connection. The cost is listed in [data sources](docs/data-sources.md). The `nopfam` word is gone, since there is nothing left to skip.
- Long downloads now print how far they have got, their speed and the time left, and long builds say once a minute that they are still going.
- A running setup remembers its data folder, so running it again after an interruption offers the same one.

### Fixed

- `varshift config` showed the saved data folder, or the default when none was saved, even when `VARSHIFT_DATA` named another folder, which is the one every command uses. In WSL, with `VARSHIFT_DATA` set in `~/.bashrc`, it named `~/varshift-data`, a folder that didn't exist. It now shows the folder in use, says when `VARSHIFT_DATA` or the default chose it, and names a saved folder that `VARSHIFT_DATA` overrides. A `VARSHIFT_DATA` or `data=` that can't be used, such as a Windows path in WSL, is now reported before any setting is saved. The help for `data=` names `VARSHIFT_DATA` too.
- On hg19 without hg38, a mitochondrial name such as `NC_012920.1:m.3243A>G` was refused as "not covered by the installed hg19 reference". The rCRS sequence it needs comes with the hg38 reference, and the error now says so, with the command to add it.
- Writing a workbook kept every link of every row in memory, though only the first 65,000 per sheet can be clickable. A large VCF has twenty per row, so this could reach gigabytes. Links past the limit are now plain text and aren't kept.
- `.gitignore` re-included everything under `demo/`, so a sample VCF dropped there to try VarShift, with its results, could have been committed. Only the demo's own files are kept now.
- `varshift add revel` failed with `403 Forbidden`, because REVEL's own server refuses scripted downloads. REVEL now comes from its Zenodo record, which its download page names as the other source, and is checked against Zenodo's MD5.
- c. and p. names on RefSeq transcripts that don't align exactly to the genome, 1,231 on hg19 and 306 of them MANE Select. Counting along RefSeq's exon rows put everything after a gap off by its length, and a c. name as input could land a few bases away. VarShift now numbers through RefSeq's own alignment. On 20,027 ClinVar names, c. SNVs went from 99.51% to 99.94% and p. SNVs from 99.38% to 99.81%. The transcript databases need rebuilding to get this.
- Setup, the demo and the benchmark no longer write to the system disk. SQLite's sort files, which can reach gigabytes while an index is built, and all temporary folders now stay beside the data or the output.
- UCSC now redirects its downloads to its Amazon bucket, `genome-browser.s3.amazonaws.com`, so setup stopped at the first step. That one bucket is now on the allow-list. Redirects still go through the same check, and the reference is still compared with UCSC's published MD5.
- Setup from WSL now refuses a data folder on a Windows drive unless you add `windowsdrive`. It was slow, and with memory tight it stopped with `Cannot allocate memory`. That error now says what to do. `varshift update` and `varshift add` refuse it too, since they build the same databases, and suggest running them from PowerShell. `varshift update check` still works there.
- On Linux and WSL, `data=I:\folder` is now refused with the right spelling. Bash removes the backslash, so it used to create a new folder called `I:folder` in the current directory and download everything again.
- On Windows, setup could stop at the ClinGen servers with `CERTIFICATE_VERIFY_FAILED`, because Windows adds a trusted root only when its own client first needs it.
- A file with a matching MD5 is no longer unpacked once more just to check it.
- ClinVar records are read about 38% faster, with byte-identical results, because one changed base no longer goes through the full normaliser.
- The update check said the gene tables had a new version every time, so each annotation showed the update notice and each `varshift update` rebuilt them for nothing. ClinGen makes its gene-validity file fresh for each request, with no ETag and a date of now. A date like that is now treated as no version, so the file is read again weekly, like other files without version headers.
- A VCF with more results than Excel can open stopped with an error while writing `results.xlsx`. The workbook now stops at Excel's limit, and its Summary sheet says how many results are only in `results.tsv`.
- When the report page was full (20,000 results), a later High result was left out while lower ones stayed. The lowest-priority results now make way first.
- A rollback that failed part way, for example because Windows held a file open, could leave a mix of old and new files. It now puts back everything it moved and says nothing was changed.
- A transcript placed on both chrX and chrY, as in the pseudoautosomal regions, could be read with the other copy's sequence.
- A VCF whose `##reference` line names one build, with no contig lengths to check, is now refused when `build=` names the other, instead of being read on the wrong build.
- An expert panel's benign classification, or a common variant, now gives priority Low. Only a ClinVar benign call did.
- c. and g. names with lowercase bases, such as `c.1799t>a`, are now read.
- Downloads and update checks now retry after a dropped TLS connection or a broken HTTP response, as they already did after a timeout.
- A damaged reference for the other build is now reported as damaged, instead of as not installed.
- In JSON output, `candidate.in_clinvar` was false for a position or c. name that is in ClinVar.

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
- `varshift explain` to show each step, file and query, and `varshift benchmark` to measure a computer.
- Several processes for large VCFs, with identical output.
- Automated setup, weekly update with approval for results-changing sources, checksums and rollback.
- Validation against ClinVar's names, a comparison tool for annotated VCFs, and a bundled demo.

### Known limits

- Accuracy has been measured on 12 genes only. See [docs/validation.md](docs/validation.md).
- The full setup has not yet been run end to end.
- CNVs, structural variants and effects on mitochondrial transcripts are not covered.
