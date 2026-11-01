# Input and output reference

## Input forms

| Kind | Accepted forms |
|---|---|
| Protein | `BRAF V600E`, `BRAF:V600E`, `BRAF p.V600E`, `BRAF p.Val600Glu`, `BRAF p.(Val600Glu)`, `BRAF (V600E)`, `NP_004324.2:p.Val600Glu` |
| Protein change types | Missense, `R132*` or `R132X` or `Arg132Ter`, `V600=`, `E746_A750del`, `L5dup`, `K10_A11insG`, `K10delinsNG`, `E23fs`, `E23Vfs*17`, `Glu23ValfsTer17`, `Ter110GlnextTer17`, `M1?` |
| cDNA | `NM_004333.6:c.1799T>A`, `NM_004333:c.1799T>A`, `NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)`, `BRAF c.1799T>A`, `ENST00000288602.11:c.1799T>A`, `NR_046018.2:n.100A>G` |
| cDNA change types | `1799T>A`, `68_69del`, `100dup`, `100_101insAT`, `100_102delinsAT`, `100_103inv`, `3901ACT[1]`, `-21G>A`, `*5C>T`, `1234+1G>A`, `1235-2A>G`, `-21-5G>A` |
| Genomic | `NC_000007.13:g.140453136A>T`, `chr7:g.140453136A>T`, `7:g.140453134_140453136del`, `NC_012920.1:m.3243A>G` |
| Position | `chr7-140453136-A-T`, `7-140453136-A-T`, `chr7:140453136:A:T`, `chr7:140453136 A>T`, `7 140453136 A T` |
| dbSNP | `rs113488022`, through ClinVar only |

Add `--gene BRAF` when a description has no gene.

## What happens with p. input

- **Substitutions:** every codon giving that amino acid is tried, so `BRAF V600E` gives `c.1799T>A` and `c.1799_1800delinsAA`.
- **In-frame deletions and duplications:** the codons named are removed or repeated.
- **Frameshifts, insertions, extensions:** a protein name doesn't define one DNA change, so VarSight looks up ClinVar's names for that gene. If ClinVar has none, it says so and asks for the c. name.
- **Wrong residue on the default transcript:** other transcripts of the gene are tried, and a match carries `other_transcript_used`.

## Transcripts shown

The default transcript comes first: MANE Select, then MANE Plus Clinical, RefSeq Select, then the longest coding RefSeq transcript. GENCODE is used only where RefSeq has nothing at that position, or when you name an `ENST` transcript.

Up to five transcripts per gene are computed by default, because some genes have hundreds (BRCA1 has 368 on hg38). `other_transcripts_omitted` says how many were left out. `--all-transcripts` shows every one.

## Result fields

| Field | Meaning |
|---|---|
| `input` | The text given, its type and the build used |
| `candidate` | Rank among candidates, how many there are, how it was found, and whether ClinVar has it |
| `variant.hg19`, `variant.hg38` | `key` is the normalised, left-aligned position in the form `chr7-140453136-A-T`. `ref_check` is `pass`, `fail` or `not_run`. `hgvs_g` is the genomic name |
| `gene` | Symbol, name, HGNC ID, NCBI Gene ID, Ensembl gene, OMIM |
| `transcript` | The default transcript. `genome_note` appears when RefSeq says the transcript differs from the genome |
| `hgvs` | `g`, `c`, `p`, and `protein_short` (such as `V600E`) |
| `effect` | Most severe Sequence Ontology term, all terms, exon, intron, amino-acid position |
| `protein` | UniProt accession, the gene's Pfam domains, and domains and UniProt features at this residue |
| `evidence.clinvar` | `status` is `found`, `not_found` or `not_installed`. When found: ID, classification, review status, condition, oncogenicity, somatic impact, frequencies listed by ClinVar, ClinVar's own names, the build matched on and the version |
| `evidence.dbsnp` | The rs number from the ClinVar record |
| `evidence.civic` | CIViC variant, link, and count of accepted evidence items |
| `evidence.clingen_validity`, `clingen_dosage` | ClinGen curation for the gene |
| `evidence.constraint` | gnomAD v2.1.1 pLI and LOEUF for the gene |
| `evidence.gene_diseases` | OMIM and Orphanet IDs linked to the gene by HPO |
| `evidence.alphamissense`, `gnomad` | Present only if you installed these sources |
| `links` | One entry per database: `name`, `label`, `url`, and a `note` when something needs saying. `url` is null when there isn't enough information, and the note says why |
| `versions` | The version of each local database used |
| `other_transcripts` | The same effect and names on other transcripts |
| `flags` | Warnings; see below |

A result whose `evidence.clinvar.status` is `not_found` means only that this variant isn't in that version of ClinVar. It says nothing about whether the variant is harmful.

## TSV and Excel columns

One row per candidate, on the default transcript. For a VCF, the original columns come first, split per ALT allele, followed by `ALT_ORIGINAL`.

`input`, `input_type`, `build`, `candidate_rank`, `candidates_total`, `hg19_key`, `hg38_key`, `gene`, `transcript`, `mane`, `hgvs_g`, `hgvs_c`, `hgvs_p`, `protein_short`, `exon`, `intron`, `effect`, `clinvar_id`, `clinvar_classification`, `clinvar_review_status`, `clinvar_condition`, `clinvar_oncogenicity`, `clinvar_somatic_impact`, `rsid`, `uniprot`, `pfam_at_residue`, `pfam_gene`, `civic_variant`, `civic_evidence_items`, `clingen_validity`, `clingen_dosage`, `gnomad_pli`, `gnomad_loeuf`, `gnomad_af`, `alphamissense`, `flags`, `error`, then one `link_` column per database, then `data_versions`.

A file being written is named `results.tsv.incomplete` until the run finishes, so a stopped run can't be mistaken for a finished one.

## Flags

| Flag | Meaning | What to do |
|---|---|---|
| `multiple_candidates` | More than one DNA change fits this input | Look at all candidates; the first is only the best ranked |
| `other_transcript_used` | The default transcript didn't fit, so another one was used | Check this is the isoform that was meant |
| `transcript_version_substituted` | The version you named isn't installed, so another version of that transcript was used | Numbering can differ between versions, particularly in UTRs. Check the result |
| `protein_version_substituted` | The same, for a protein accession | Check the result |
| `transcript_differs_from_genome` | RefSeq says this transcript doesn't exactly match the genome here | c. and p. are computed from the genome, so check them against the transcript |
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
| `other_build_not_installed`, `liftover_data_not_installed` | Only one build is set up | Install both with `varsight setup --build both` |
| `other_build_region_not_installed` | That region isn't in the installed reference | Only happens with the demo data |
| `uniprot_numbering_differs` | The residue differs between the transcript's protein and UniProt, so domain positions aren't reported | Compare in UniProt |
| `split_from_multiallelic` | This row is one ALT from a multi-allelic VCF line | Informational |
| `symbolic_allele_skipped` | The ALT is `<DEL>`, `*` or similar | Not annotated; CNVs and SVs aren't covered |
| `ref_mismatch` | The VCF's REF isn't what the genome has | Often a wrong build or wrong file. See `error` |
| `not_annotated` | This record couldn't be processed | See `error` |
