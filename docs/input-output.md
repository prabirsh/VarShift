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

Add `gene=BRAF` when a description has no gene.

**Always put a variant in quotes.** A shell can read `>` as "write to a file", so type `varsight annotate "NM_004333.6:c.1799T>A"`. On Windows the installed `varsight.exe` keeps quoted text intact in both PowerShell and Command Prompt.

## What happens with p. input

- **Substitutions:** every codon giving that amino acid is tried, so `BRAF V600E` gives `c.1799T>A` and `c.1799_1800delinsAA`.
- **In-frame deletions and duplications:** the codons named are removed or repeated.
- **Frameshifts, insertions, extensions:** a protein name doesn't define one DNA change, so VarSight looks up ClinVar's names for that gene. If ClinVar has none, it says so and asks for the c. name.
- **Wrong residue on the default transcript:** other transcripts of the gene are tried, and a match carries `other_transcript_used`.

## Transcripts shown

The default transcript comes first: MANE Select, then MANE Plus Clinical, RefSeq Select, then the longest coding RefSeq transcript. GENCODE is used only where RefSeq has nothing at that position, or when you name an `ENST` transcript.

Up to five transcripts per gene are computed by default, because some genes have hundreds (BRCA1 has 368 on hg38). `other_transcripts_omitted` says how many were left out. `transcripts=all` shows every one.

## Result fields

| Field | Meaning |
| - | - |
| `input` | The text given, its type and the build used |
| `candidate` | Rank among candidates, how many there are, how it was found, and whether ClinVar has it |
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
| One variant | A readable card on the screen. Add `format=json`, `format=jsonl` or `format=tsv` for programs, or `out=name.xlsx` for a file |
| A list or a VCF | A folder named after the input, holding `report.html`, `results.xlsx` and `results.tsv`. Name it with `out=FOLDER` |
| A list or a VCF with `out=name.xlsx` | Just that one file. The extension chooses the format: `.xlsx`, `.html`, `.tsv`, `.json` or `.jsonl` |

Existing output is never replaced unless you add the word `replace`.

### The Excel workbook

| Sheet | Holds |
| - | - |
| Summary | Input, date, counts by priority, effect and gene, and any update notice |
| Priority | High, Medium and Review results, sorted, with reasons and the main links. Start here |
| All results | Every result and every column, with filters |
| Warnings explained | Each warning found, how often, what it means and what to do |
| Data versions | The version of every database used |

### The HTML report

One file that opens with a double click. It needs no server and no internet, and it can be emailed. It keeps up to 20,000 results. When there are more, Low results leave first and the page says how many, while Excel and TSV still hold everything.

Opening a link in it sends that variant to that website. Nothing else in the report uses the network.

### Columns in TSV and Excel

One row per candidate, on the default transcript. For a VCF, the original columns come first, split per ALT allele, followed by `ALT_ORIGINAL`.

`priority`, `priority_reasons`, `input`, `input_type`, `build`, `candidate_rank`, `candidates_total`, `hg19_key`, `hg38_key`, `gene`, `transcript`, `mane`, `hgvs_g`, `hgvs_c`, `hgvs_p`, `protein_short`, `exon`, `intron`, `effect`, `clinvar_id`, `clinvar_classification`, `clinvar_review_status`, `clinvar_condition`, `clinvar_oncogenicity`, `clinvar_somatic_impact`, `rsid`, `uniprot`, `pfam_at_residue`, `pfam_gene`, `civic_variant`, `civic_evidence_items`, `clingen_validity`, `clingen_dosage`, `gnomad_pli`, `gnomad_loeuf`, then the columns of any advanced database you installed, then `flags`, `error`, one `link_` column per database, and `data_versions`.

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

High takes precedence. A variant reported as pathogenic stays High even when common, and the frequency is added to its reasons so the conflict stays visible. The rule is in [`priority.py`](../src/varsight/priority.py).

## Large VCFs

For a VCF over 1 MB, VarSight uses up to 4 processes, each with its own read-only connections. The output is identical and in the same order as with one process, and a test checks this. Choose the number with `workers=`, for example `workers=1` to leave the computer free. Each process uses about 75 MB of memory, and writing all three reports for 30,000 records used about 500 MB more. On a very tight machine, `out=results.tsv` needs the least. [Measured figures](benchmarks.md).

## Flags

| Flag | Meaning | What to do |
| - | - | - |
| `multiple_candidates` | More than one DNA change fits this input | Look at all candidates; the first is only the best ranked |
| `other_transcript_used` | The default transcript didn't fit, so another one was used | Check this is the isoform that was meant |
| `transcript_version_substituted` | The version you named isn't installed, so another version of that transcript was used | Numbering can differ between versions, particularly in UTRs. Check the result |
| `protein_version_substituted` | The same, for a protein accession | Check the result |
| `transcript_differs_from_genome` | RefSeq says this transcript doesn't exactly match the genome here | c. and p. are computed from the genome. If the difference is an indel or frameshift, numbers after it are off by its length, and a c. name given as input can land a few bases away. Check them against the transcript. See [validation](validation.md#what-the-mismatches-are) |
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
| `other_build_not_installed`, `liftover_data_not_installed` | Only one build is set up | Install both with `varsight setup build=both` |
| `other_build_region_not_installed` | That region isn't in the installed reference | Only happens with the demo data |
| `uniprot_numbering_differs` | The residue differs between the transcript's protein and UniProt, so domain positions aren't reported | Compare in UniProt |
| `split_from_multiallelic` | This row is one ALT from a multi-allelic VCF line | Informational |
| `symbolic_allele_skipped` | The ALT is `<DEL>`, `*` or similar | Not annotated; CNVs and SVs aren't covered |
| `ref_mismatch` | The VCF's REF isn't what the genome has | Often a wrong build or wrong file. See `error` |
| `not_annotated` | This record couldn't be processed | See `error` |
