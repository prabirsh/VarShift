# Input and output reference

## Input forms

| Kind | Accepted forms |
| - | - |
| Protein | `BRAF V600E`, `BRAF:V600E`, `BRAF p.V600E`, `BRAF p.Val600Glu`, `BRAF p.(Val600Glu)`, `BRAF (V600E)`, `NP_004324.2:p.Val600Glu` |
| Protein change types | Missense, `R132*` or `R132X` or `Arg132Ter`, `V600=`, `E746_A750del`, `L5dup`, `K10_A11insG`, `K10delinsNG`, `E23fs`, `E23Vfs*17`, `Glu23ValfsTer17`, `Ter110GlnextTer17`, `M1?` |
| cDNA | `NM_004333.6:c.1799T>A`, `NM_004333:c.1799T>A`, `NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)`, `BRAF c.1799T>A`, `ENST00000288602.11:c.1799T>A`, `NR_046018.2:n.100A>G` |
| cDNA change types | `1799T>A`, `68_69del`, `100dup`, `100_101insAT`, `100_102delinsAT`, `100_103inv`, `3901ACT[1]`, `-21G>A`, `*5C>T`, `1234+1G>A`, `1235-2A>G`, `-21-5G>A` |
| Genomic | `NC_000007.13:g.140453136A>T`, `chr7:g.140453136A>T`, `7:g.140453134_140453136del`, `NC_012920.1:m.3243A>G` |
| Position | `chr7-140453136-A-T`, `7-140453136-A-T`, `chr7:140453136:A:T`, `chr7:140453136 A>T`, `7 140453136 A T` |
| dbSNP | `rs113488022`. Through ClinVar, or through dbSNP itself if you added it |
| Several | `BRAF V600E; KRAS G12D; TP53 R175H`, separated by `;`. Each is handled like a line of a list file, so one that fails doesn't stop the others. A `;` inside square brackets, as in `c.[76A>C;83G>C]`, belongs to the change |

Add `gene=BRAF` when a description has no gene. `varshift examples` prints a working example of each form, and a test annotates every one of them.

**Always put a variant in quotes.** A shell can read `>` as "write to a file", so type `varshift annotate "NM_004333.6:c.1799T>A"`. On Windows the installed `varshift.exe` keeps quoted text intact in both PowerShell and Command Prompt.

## When the text is a little off

Every description goes through two steps before it is annotated.

**1. Slips with only one sensible reading are fixed,** and the result says so: `input.raw` keeps the text as typed, `input.read_as` holds the text VarShift read, and `input.read_as_notes` names each change. The card shows them on a `Read as` line. Text that already reads is never changed.

| Typed | Read as | Rule |
| - | - | - |
| `braf val600glu`, `BRAF p.v600e`, `EGFR e746_a750del` | `BRAF Val600Glu`, `BRAF p.V600E`, `EGFR E746_A750del` | Amino acids in standard case, and the gene symbol as HGNC writes it |
| `BRAF(V600E)` | `BRAF V600E` | Brackets straight after the gene |
| `BRAF c.1799T-A`, `c.1799T->A`, `c.1799T/A` | `BRAF c.1799T>A` | Any arrow, dash or slash between the two bases of a substitution |
| `NM_004333.6:1799T>A` | `NM_004333.6:c.1799T>A` | A transcript with no level takes `c.`, or `n.` for `NR_` and `XR_` |
| `chr7:140453136A>T` | `chr7:g.140453136A>T` | A chromosome with no level takes `g.`, or `m.` for the mitochondrial genome |
| `chr7:140,453,136 A>T` | `chr7:140453136 A>T` | Thousands separators |
| Smart quotes, long dashes, a fullwidth `>`, odd spaces, a trailing `.` | Replaced or removed | Text pasted from a document or an email |
| `HER2 V777L` | `ERBB2 V777L (HER2 is an alias)` | A gene alias, resolved through HGNC |

**2. When a description still fails, VarShift offers what it probably meant,** built from the installed data. Nothing is applied without being chosen.

| Typed | Offered |
| - | - |
| `BRAF V600` | The changes ClinVar and CIViC know at that residue, with their classification and evidence count |
| `V600E` | The genes in which ClinVar or CIViC records that change. This search reads the whole ClinVar table, so files skip it |
| `BRFA V600E` | Gene symbols one letter away, or two for longer symbols, a swapped pair counting as one. Ranked by how many ClinVar records each gene has |
| `BRAF W600E` | The residue the protein really has: `BRAF V600E` |
| `BRAF c.1799A>T` | The base the transcript really has: `c.1799T>A`, with REF and ALT swapped, or the REF corrected |
| `chr7-140753336-A-T build=hg19` | The other build, when its REF matches there and it is installed. Otherwise a note that it may be on the other build |

Where the options go depends on who can answer:

- **One description at a terminal,** with `varshift annotate` or in the guided start: a numbered question. Your choice is annotated as if you had typed it, `input.raw` keeps what you typed, and `Read as` says it was chosen.
- **A script, a pipe, `out=` or `format=`:** the error, then `Did you mean: ...`, and exit code 1. Nothing waits for an answer.
- **A list file or a `;` list:** the line stays an error row, and its options go in the `suggestion` column. Only the quick suggestions run there, for the first 500 failing lines, and the summary says how many rows have one.
- **A VCF:** never tidied or suggested. VCFs are written by machines, and a change there would hide a real problem.

When there is nothing to suggest, the error is followed by `Not sure of the format? See: varshift examples`.

### A typo that fits another transcript

A description can fail on its gene's default transcript but fit another one, as `BRAF A600E` does: residue 600 is Val on NM_004333.6, the MANE Select transcript, but Ala on NM_001378475.1. That is more often a typo than a choice of isoform, so it isn't passed over in silence:

- **At a terminal,** VarShift asks: `NM_004333.6, the default transcript, has Val at residue 600, not Ala, so NM_001378475.1 was used. Which did you mean?`, with `BRAF V600E` on the default transcript first, and the description as typed second.
- **Everywhere else,** the answer stands, with the warning `default_transcript_differs`, a `Note` line on the card, `input.default_transcript` in JSON, and the corrected description in the `suggestion` column.

The same holds for a c. name whose base fits another transcript only, such as `BRAF c.1799C>G`. Naming the transcript, as in `NM_001378475.1:c.1799_1800delinsAG`, says which you meant, and nothing is asked.

## What happens with p. input

- **Substitutions:** every codon giving that amino acid is tried, so `BRAF V600E` gives `c.1799T>A` and `c.1799_1800delinsAA`.
- **In-frame deletions and duplications:** the codons named are removed or repeated.
- **Frameshifts, insertions, extensions:** a protein name doesn't define one DNA change, so VarShift looks up ClinVar's names for that gene. If ClinVar has none, it says so and asks for the c. name.
- **Wrong residue on the default transcript:** other transcripts of the gene are tried, and a match carries `other_transcript_used`. That is usually a typo, so VarShift says what the default transcript has there: see [a typo that fits another transcript](#a-typo-that-fits-another-transcript).

## One answer

The candidates for one description are sorted into one kind of answer, recorded in `input.answer`:

| Kind | When | Shown as |
| - | - | - |
| `single` | One candidate | Its card |
| `clear` | Several, and exactly one is in ClinVar. Or none is, and exactly one changes the fewest bases | The best candidate's card, with `Why this` and `Also` lines |
| `choice` | An rs number with several alleles, several candidates in ClinVar, or a tie on bases changed | A numbered table. At a terminal you choose one, `a` for all, or Enter for the first. Otherwise the table and how to see every candidate |

The best candidate comes first, with `candidate.best` true and `candidate.why` naming the rule that chose it, plus `MANE Select transcript` when its transcript is. `also` names up to three of the others. Each allele of an rs number is a different variant, so nothing is preferred among them. `show=all` prints every candidate in full, as before.

## Transcripts shown

The default transcript comes first: MANE Select, then MANE Plus Clinical, RefSeq Select, then the longest coding RefSeq transcript. GENCODE is used only where RefSeq has nothing at that position, or when you name an `ENST` transcript.

Up to five transcripts per gene are computed by default, because some genes have hundreds (BRCA1 has 368 on hg38). `other_transcripts_omitted` says how many were left out. `transcripts=all` shows every one.

## Result fields

| Field | Meaning |
| - | - |
| `input` | `raw` is the text as given, `type` its kind, and `build` the build used. `read_as` is the text read after a slip was fixed, or empty, and `read_as_notes` names each change. `answer` is `single`, `clear` or `choice`. `mode` is `tumour`, `inherited` or null |
| `candidate` | `rank` among candidates and `of` how many there are, `via` how it was found, `in_clinvar`, `bases_changed`, `best` (true on the answer shown), and `why` (the rules that chose it) |
| `also` | On the best candidate of a `clear` answer: the others, as in `c.1799_1800delinsAA (2 bases, not in ClinVar)` |
| `suggestions` | On a list or `;` row that failed: what it probably meant, each a description you can type. Also on a result that fitted only a transcript other than the default: the description written to fit the default |
| `input.default_transcript` | Only then: the default transcript's `accession`, what it `says` it has there, and the `suggestion` |
| `variant.hg19`, `variant.hg38` | `key` is the normalised, left-aligned position in the form `chr7-140453136-A-T`. `ref_check` is `pass`, `fail` or `not_run`. `hgvs_g` is the genomic name |
| `gene` | Symbol, name, HGNC ID, NCBI Gene ID, Ensembl gene, OMIM |
| `transcript` | The default transcript. `genome_note` appears when RefSeq says the transcript differs from the genome |
| `hgvs` | `g`, `c`, `p`, and `protein_short` (such as `V600E`) |
| `effect` | Most severe Sequence Ontology term, all terms, exon, intron, amino-acid position |
| `protein` | UniProt accession, the gene's Pfam domains, and domains and UniProt features at this residue |
| `evidence.clinvar` | `status` is `found`, `not_found` or `not_installed`. When found: ID, classification, review status, condition, oncogenicity, somatic impact, frequencies listed by ClinVar, ClinVar's own names, the build matched on and the version |
| `evidence.dbsnp` | The rs number from the ClinVar record, or from dbSNP if you added it |
| `evidence.civic` | CIViC variant, link, and count of accepted evidence items |
| `evidence.clingen_validity`, `clingen_dosage` | ClinGen curation for the gene |
| `evidence.constraint` | gnomAD v2.1.1 pLI and LOEUF for the gene |
| `evidence.gene_diseases` | OMIM and Orphanet IDs linked to the gene by HPO |
| `evidence.advanced` | One entry per advanced database you installed. Each has a `status` of `found` or `not_found`, plus its values |
| `priority` | `label` is High, Medium, Review or Low, `level` is 1 to 4, and `reasons` lists why |
| `links` | One entry per database: `name`, `label`, `url`, and a `note` when something needs saying. `url` is null when there isn't enough information, and the note says why |
| `versions` | The version of each local database used |
| `other_transcripts` | The same effect and names on other transcripts |
| `flags` | Warnings; see below |

A result whose `evidence.clinvar.status` is `not_found` means only that this variant isn't in that version of ClinVar. It says nothing about whether the variant is harmful.

## What you get

| You give | You get |
| - | - |
| One variant | Its answer on the screen: the best candidate's card, or a choice ([One answer](#one-answer)). `show=all` prints every candidate in full and `links=all` every link with its label and note. Add `format=json`, `format=jsonl` or `format=tsv` for programs, or `out=name.xlsx` for a file |
| Several, separated by `;` | One numbered line each on the screen, or all of them in the file `out=` names |
| A list or a VCF | A folder named after the input, holding `report.html`, `results.xlsx` and `results.tsv`. Name it with `out=FOLDER` |
| A list or a VCF with `out=name.xlsx` | Just that one file. The extension chooses the format: `.xlsx`, `.html`, `.tsv`, `.json` or `.jsonl` |

Existing output is never replaced unless you add the word `replace`.

### The card on the screen

`mode=tumour` puts CIViC first, then ClinVar's oncogenicity and somatic impact, and links to ClinVar, OncoKB, CIViC, COSMIC and cBioPortal. `mode=inherited` puts the ClinVar classification and review status first, then ClinGen gene validity and gnomAD constraint, and links to ClinVar, gnomAD, ClinGen, OMIM and Franklin. With no mode, the links are ClinVar, gnomAD, UCSC Genome Browser, dbSNP and VarSome. `varshift config mode=tumour` saves a mode, and `mode=none` clears it. The mode is recorded in `input.mode` and the `mode` column.

Warnings are written as sentences on `Note` lines, while the codes stay in JSON and TSV. Lines wrap at the terminal's width. Colour is used only in a terminal, and never when `NO_COLOR` is set or `color=no` is given or saved.

### The Excel workbook

| Sheet | Holds |
| - | - |
| Summary | Input, date, counts by priority, effect and gene, and any update notice |
| Priority | High, Medium and Review results, sorted, with reasons and the main links. Start here |
| All results | Every result and every column, with filters |
| Links | One row per result and database: the result number, input, gene, c. and p. names, the database, a clickable link, its label (confirmed ID or built), and its note. A database with no link for that result says why |
| Warnings explained | Each warning found, how often, what it means and what to do |
| Data versions | The version of every database used |

### The HTML report

One file that opens with a double click. It needs no server and no internet, and it can be emailed. It keeps up to 20,000 results. When there are more, Low results leave first and the page says how many, while Excel and TSV still hold everything.

Opening a link in it sends that variant to that website. Nothing else in the report uses the network.

### Columns in TSV and Excel

One row per candidate, on the default transcript. For a VCF, the original columns come first, split per ALT allele, followed by `ALT_ORIGINAL`.

`priority`, `priority_reasons`, `input`, `input_type`, `build`, `candidate_rank`, `candidates_total`, `hg19_key`, `hg38_key`, `gene`, `transcript`, `mane`, `hgvs_g`, `hgvs_c`, `hgvs_p`, `protein_short`, `exon`, `intron`, `effect`, `clinvar_id`, `clinvar_classification`, `clinvar_review_status`, `clinvar_condition`, `clinvar_oncogenicity`, `clinvar_somatic_impact`, `rsid`, `uniprot`, `pfam_at_residue`, `pfam_gene`, `civic_variant`, `civic_evidence_items`, `clingen_validity`, `clingen_dosage`, `gnomad_pli`, `gnomad_loeuf`, then the columns of any advanced database you installed, then `flags`, `error`, one `link_` column per database, and `data_versions`.

The last seven columns come after `data_versions`, so scripts that read columns by position keep working:

| Column | Holds |
| - | - |
| `read_as` | The text VarShift read, when it fixed a slip. Empty otherwise |
| `answer` | `single`, `clear` or `choice` |
| `best` | `yes` on the row shown as the answer, `no` on the others |
| `why` | Why that row is the answer, on the best row |
| `also` | The other candidates, on the best row of a `clear` answer |
| `suggestion` | On a row that failed: what it probably meant, separated by `; `. On a row that fitted only a transcript other than the default: the description written to fit the default |
| `mode` | `tumour`, `inherited` or empty |

Excel can hold at most 1,048,576 rows and 65,000 clickable links per sheet. Past those limits the Summary sheet says what was left out, and `results.tsv` has everything.

### Columns added by advanced databases

| Database | Columns |
| - | - |
| dbSNP | `dbsnp_rsid`, `topmed_af`, `g1k_af_in_dbsnp`, `gnomad_genomes_af_in_dbsnp`, `gnomad_exomes_af_in_dbsnp`, `exac_af_in_dbsnp`, `alfa_af` |
| 1000 Genomes | `g1k_af`, `g1k_afr_af`, `g1k_amr_af`, `g1k_eas_af`, `g1k_eur_af`, `g1k_sas_af` |
| gnomAD | `gnomad_af`, `gnomad_popmax`, `gnomad_popmax_af`, `gnomad_homozygotes`, `gnomad_filter` |
| gnomAD SV | `gnomad_sv_overlapping`, `gnomad_sv_top` |
| REVEL | `revel` |
| AlphaMissense | `alphamissense_score`, `alphamissense_class` |
| SpliceAI | `spliceai_max`, `spliceai_gene` |
| PrimateAI | `primateai` |
| phyloP | `phylop100way` |
| Cancer Hotspots | `cancer_hotspot`, `cancer_hotspot_tumours` |
| DECIPHER | `decipher_hi_percent`, `ddg2p` |
| ClinGen classifications | `clingen_classification`, `clingen_evidence_codes` |
| COSMIC | `cosmic_id`, `cosmic_samples` |
| COSMIC Cancer Gene Census | `cgc_tier`, `cgc_role`, `cgc_tumour_types` |
| COSMIC gene fusions | `cosmic_fusion_partners` |
| OMIM | `omim_phenotypes` |

A file being written is named `results.tsv.incomplete` until the run finishes, so a stopped run can't be mistaken for a finished one.

## Review priority

Priority sorts results so you know where to start. **It is not a classification,** and it never says a variant is pathogenic or benign. Every result states its reasons.

| Level | Given when |
| - | - |
| High | ClinVar or a ClinGen expert panel reports pathogenic or likely pathogenic, ClinVar reports oncogenic or AMP tier I or II, CIViC has evidence items, or the residue is a cancer hotspot |
| Medium | ClinVar has conflicting reports, the change is predicted to disrupt the protein (stop gained, frameshift, splice site, start or stop lost), REVEL is at least 0.644, AlphaMissense says pathogenic, SpliceAI is at least 0.2, or COSMIC has at least 10 samples |
| Review | ClinVar says uncertain significance, or the change alters the protein or nearby splicing |
| Low | ClinVar or ClinGen reports benign, frequency is at least 5% in a population database, or there is no protein change and no record of concern |

High takes precedence. A variant reported as pathogenic stays High even when common, and the frequency is added to its reasons so the conflict stays visible. The rule is in [`priority.py`](../src/varshift/priority.py).

## Large VCFs

For a VCF over 1 MB, VarShift uses up to 4 processes, each with its own read-only connections. The output is identical and in the same order as with one process, and a test checks this. Choose the number with `workers=`, for example `workers=1` to leave the computer free. Each process uses about 75 MB of memory, and writing all three reports for 30,000 records used about 500 MB more. On a very tight machine, `out=results.tsv` needs the least. [Measured figures](benchmarks.md).

## Flags

| Flag | Meaning | What to do |
| - | - | - |
| `multiple_candidates` | More than one DNA change fits this input | Look at all candidates; the first is only the best ranked |
| `other_transcript_used` | The default transcript didn't fit, so another one was used | Check this is the isoform that was meant |
| `default_transcript_differs` | The default transcript has a different residue or base here, so another transcript was used. This is often a typo | See `suggestion`: the change written to fit the default transcript. At a terminal, VarShift asks which you meant |
| `transcript_version_substituted` | The version you named isn't installed, so another version of that transcript was used | Numbering can differ between versions, particularly in UTRs. Check the result |
| `protein_version_substituted` | The same, for a protein accession | Check the result |
| `transcript_differs_from_genome` | RefSeq says this transcript doesn't exactly match the genome here | Numbering follows RefSeq's alignment, but the genome's bases can differ from the transcript's. Check any change near the difference. See [validation](validation.md#what-the-mismatches-are) |
| `transcript_gap_here` | This change touches bases that the transcript and the genome don't share, so no full name can be given | Review against the transcript by hand |
| `mane_version_differs` | The installed transcript is a different version from the MANE one | Common on hg19. Numbering usually agrees, but check |
| `spans_exon_intron_boundary` | The change crosses an exon-intron boundary, so no protein name is given | Review splicing by hand |
| `overlaps_start_codon`, `overlaps_stop_codon` | The change runs over the start or stop codon | Review by hand |
| `cds_internal_stop` | The transcript has a stop inside its coding region, such as a selenoprotein | No protein name is given |
| `cds_incomplete` | The coding region isn't a whole number of codons, or has no stop | No protein name is given |
| `chrM_not_rCRS` | hg19 `chrM` is the old sequence (16571 bp), not the rCRS that ClinVar and RefSeq use | Use `MT` or `NC_012920.1` if your position is on rCRS |
| `mitochondrial_no_transcript_model` | Mitochondrial genes use another genetic code and aren't modelled | Position and ClinVar are still reported |
| `no_liftover` | The position doesn't map to the other build | The other build's key is empty |
| `liftover_span_changed` | The reference span has a different length on the other build | The other build's key is empty |
| `liftover_strand_flipped` | The region is reversed on the other build | Informational; alleles are complemented |
| `ref_mismatch_on_other_build` | The mapped position has a different reference base | The other build's key is empty |
| `ref_alt_swapped_on_other_build` | The other build's reference is this variant's ALT | The "variant" may be the reference allele there |
| `other_build_not_installed`, `liftover_data_not_installed` | Only one build is set up | Install both with `varshift setup build=both` |
| `other_build_region_not_installed` | That region isn't in the installed reference | Only happens with the demo data |
| `uniprot_numbering_differs` | The residue differs between the transcript's protein and UniProt, so domain positions aren't reported | Compare in UniProt |
| `split_from_multiallelic` | This row is one ALT from a multi-allelic VCF line | Informational |
| `symbolic_allele_skipped` | The ALT is `<DEL>`, `*` or similar | Not annotated; CNVs and SVs aren't covered |
| `ref_mismatch` | The VCF's REF isn't what the genome has | Often a wrong build or wrong file. See `error` |
| `not_annotated` | This record couldn't be processed | See `error` |
