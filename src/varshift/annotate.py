"""The annotation pipeline: parse, resolve, normalise, cross-build, consequence, evidence, links."""

from __future__ import annotations

from .advanced import ALIAS, DATABASES
from .answer import apply as apply_answer, decide
from .builds import BUILDS, HG19, HG38, is_mitochondrial, key_contig, parse_build
from .consequence import UPSTREAM_DISTANCE, SequenceCache, TranscriptEffect, annotate as effect_on, hgvs_g
from .errors import DataError, InputError
from .inputs import ParsedInput, parse_input
from .links import build_links
from .priority import assign as assign_priority
from .liftover import LiftResult, lift
from .resolver import Candidate, Resolver
from .store import DataStore
from .tidy import Reading, replace_word, tidy
from .variant import Variant


MAX_TRANSCRIPTS_PER_GENE = 5


class Annotator:
    def __init__(self, store: DataStore, all_transcripts: bool = False):
        self.store = store
        self.all_transcripts = all_transcripts
        self.omitted_transcripts = 0
        self.resolver = Resolver(store)
        self._sequences: dict[str, SequenceCache] = {}
        # The manifest doesn't change during a run, so this is read once rather than twice per record.
        self._versions = store.versions()
        installed = store.installed()
        self._installed_advanced = [
            database
            for database in DATABASES.values()
            if database.access != ALIAS and any(a == database.name or a.startswith(f"{database.name}.") for a in installed)
        ]

    def _pick_build(self, parsed: ParsedInput, build: str | None) -> str:
        installed = self.store.builds
        if build:
            build = parse_build(build)
            if parsed.build_hint and parsed.build_hint != build:
                raise InputError(f"{parsed.contig} belongs to {parsed.build_hint}, but build= says {build}")
        elif parsed.build_hint:
            build = parsed.build_hint
        elif parsed.kind in ("position", "genomic"):
            raise InputError("a genomic position needs a build; add build=hg19 or build=hg38")
        else:
            build = HG38 if HG38 in installed else (installed[0] if installed else None)
        if build not in installed:
            raise DataError(f"{build} data is not installed here; run: varshift setup build={build}")
        return build

    def annotate_text(self, text: str, build: str | None = None, gene: str = "", mode: str | None = None) -> list[dict]:
        """Annotate one description. Returns one result per candidate, best candidate first.

        Slips with only one sensible reading are fixed first (tidy.py), and the reading is recorded.
        """
        reading = tidy(text, gene)
        parsed = parse_input(reading.text, gene=gene)
        chosen = self._pick_build(parsed, build)
        candidates = self.resolver.resolve(parsed, chosen)
        read_as, notes = self._read_as(reading)
        input_info = {
            "raw": text.strip(),
            "type": parsed.kind,
            "build": chosen,
            "read_as": read_as,
            "read_as_notes": notes,
            # single, clear or choice: how the candidates add up to an answer, from answer.py.
            "answer": None,
            "mode": mode,
        }
        results = []
        for rank, candidate in enumerate(candidates, 1):
            record = self.annotate_candidate(candidate)
            record["input"] = input_info
            record["candidate"] = {
                "rank": rank,
                "of": len(candidates),
                "via": candidate.via,
                # The resolver only checks ClinVar for some input kinds, so the record's own lookup decides.
                "in_clinvar": candidate.in_clinvar or record["evidence"]["clinvar"].get("status") == "found",
                "bases_changed": candidate.bases_changed,
                "best": False,
                "why": [],
            }
            if len(candidates) > 1 and "multiple_candidates" not in record["flags"]:
                record["flags"].append("multiple_candidates")
            results.append(record)
        # The best candidate, if one stands out, comes first and says why.
        return apply_answer(decide(results))

    def _read_as(self, reading: Reading) -> tuple[str, list[str]]:
        """The text as VarShift read it, with the gene as HGNC writes it, or "" when nothing changed."""
        text, notes = reading.text, list(reading.notes)
        named = self.resolver.last_gene
        if named and named[0] != named[1]:
            typed, approved = named
            # A gene given with gene= isn't in the text, so there is nothing to rewrite.
            rewritten = replace_word(text, typed, approved)
            if rewritten != text and typed.upper() == approved.upper():
                text = rewritten
                notes.append("gene symbol written as HGNC writes it")
            elif rewritten != text:
                text = f"{rewritten} ({typed} is an alias)"
        if text == reading.text and not notes:
            return "", []
        return text, notes

    def annotate_candidate(self, candidate: Candidate) -> dict:
        variant = candidate.variant
        flags = list(candidate.flags)
        if variant.contig == "chrM" and variant.build == HG19:
            flags.append("chrM_not_rCRS")

        effects = self._effects(variant, candidate.transcript)
        chosen = effects[0] if effects else None
        if chosen:
            for flag in chosen.flags:
                if flag not in flags:
                    flags.append(flag)

        positions, lifted = self._both_builds(variant, flags)
        clinvar_variant, clinvar_row = self._clinvar(variant, lifted)
        symbol = chosen.transcript.gene if chosen else ""
        if not symbol and clinvar_row is not None:
            symbol = (clinvar_row["gene"] or "").split("|")[0]
        gene = self._gene_block(symbol)
        protein = self._protein_block(symbol, chosen, flags)
        evidence = self._evidence(variant, lifted, clinvar_row, clinvar_variant, symbol, chosen, flags)

        hgvs = {"g": positions[variant.build].get("hgvs_g"), "c": None, "p": None, "protein_short": None}
        effect = {"type": "intergenic_variant", "all": ["intergenic_variant"], "exon": None, "intron": None, "aa_position": None}
        if is_mitochondrial(variant.contig):
            # Mitochondrial transcripts are not modelled (they use another genetic code), so no effect is claimed.
            effect = {"type": None, "all": [], "exon": None, "intron": None, "aa_position": None}
            flags.append("mitochondrial_no_transcript_model")
        transcript = None
        if chosen:
            hgvs.update(c=chosen.hgvs_c, p=chosen.hgvs_p, protein_short=chosen.protein_short)
            effect = {
                "type": chosen.effect,
                "all": chosen.terms,
                "exon": chosen.exon,
                "intron": chosen.intron,
                "aa_position": chosen.aa_position,
            }
            tx = chosen.transcript
            transcript = {
                "accession": tx.accession,
                "source": tx.source,
                "mane": tx.mane or None,
                "strand": "+" if tx.strand > 0 else "-",
                "protein": tx.protein or None,
                "exon_count": len(tx.exons),
            }
            if tx.genome_note:
                transcript["genome_note"] = tx.genome_note

        record = {
            "input": {},
            "candidate": {"rank": 1, "of": 1},
            "variant": positions,
            "gene": gene,
            "transcript": transcript,
            "hgvs": hgvs,
            "effect": effect,
            "protein": protein,
            "evidence": evidence,
            "links": [],
            "versions": self._versions,
            "other_transcripts": [_other(e) for e in effects[1:]],
            "other_transcripts_omitted": self.omitted_transcripts,
            "flags": flags,
        }
        record["links"] = build_links(_link_values(variant, record, clinvar_row))
        record["priority"] = assign_priority(record)
        return record

    # steps
    def _effects(self, variant: Variant, preferred=None) -> list[TranscriptEffect]:
        start, end, _ = variant.edit()
        if not self.store.has(f"transcripts.{variant.build}"):
            return []
        genome = self.store.genome(variant.build)
        cache = self._sequences.setdefault(variant.build, SequenceCache(genome))
        overlapping = self.store.transcripts_at(
            variant.build, variant.contig, start, max(end, start + 1), flank=UPSTREAM_DISTANCE
        )
        overlapping = self._limit_per_gene(overlapping, preferred)
        effects = [effect_on(genome, variant, tx, cache) for tx in overlapping]

        def order(effect: TranscriptEffect):
            tx = effect.transcript
            is_preferred = preferred is not None and tx.accession == preferred.accession
            outside = effect.effect in ("upstream_gene_variant", "downstream_gene_variant")
            return (not is_preferred, outside, tx.default_rank())

        return sorted(effects, key=order)

    def _limit_per_gene(self, transcripts, preferred):
        """Keep the best few transcripts per gene. BRCA1 alone has 368 RefSeq transcripts on hg38."""
        self.omitted_transcripts = 0
        if self.all_transcripts:
            return transcripts
        groups: dict[str, list] = {}
        for tx in sorted(transcripts, key=lambda t: t.default_rank()):
            groups.setdefault(tx.gene, []).append(tx)
        kept = []
        for items in groups.values():
            chosen = items[:MAX_TRANSCRIPTS_PER_GENE]
            if preferred is not None and preferred.gene == items[0].gene and preferred not in chosen:
                chosen.append(preferred)
            kept.extend(chosen)
            # The preferred transcript may come from outside this group, so count what was really left out.
            self.omitted_transcripts += sum(1 for tx in items if tx not in chosen)
        return kept

    def _both_builds(self, variant: Variant, flags: list[str]) -> tuple[dict, Variant | None]:
        genome = self.store.genome(variant.build)
        positions: dict[str, dict] = {
            variant.build: {"key": variant.key, "ref_check": "pass", "hgvs_g": hgvs_g(genome, variant)}
        }
        other = HG38 if variant.build == HG19 else HG19
        result: LiftResult = lift(self.store, variant, other)
        for flag in result.flags:
            if flag not in flags:
                flags.append(flag)
        if result.variant is not None:
            positions[other] = {
                "key": result.variant.key,
                "ref_check": result.ref_check,
                "hgvs_g": hgvs_g(self.store.genome(other), result.variant),
            }
        else:
            positions[other] = {"key": None, "ref_check": result.ref_check, "hgvs_g": None}
        return {build: positions[build] for build in BUILDS}, result.variant

    def _clinvar(self, variant: Variant, lifted: Variant | None):
        if variant.contig == "chrM" and variant.build == HG19:
            return None, None
        for candidate in (variant, lifted):
            if candidate is None:
                continue
            row = self.store.clinvar_row(candidate.build, candidate.contig, candidate.pos, candidate.ref, candidate.alt)
            if row is not None:
                return candidate, row
        return None, None

    def _gene_block(self, symbol: str) -> dict:
        block = {"symbol": symbol or None, "name": None, "hgnc": None, "ncbi_gene": None, "ensembl": None, "omim": None}
        row = self.store.gene_info(symbol) if symbol else None
        if row is not None:
            block.update(
                name=row["name"],
                hgnc=row["hgnc_id"],
                ncbi_gene=row["ncbi_gene"],
                ensembl=row["ensembl_gene"],
                omim=row["omim"],
            )
        return block

    def _protein_block(self, symbol: str, chosen: TranscriptEffect | None, flags: list[str]) -> dict:
        block = {"uniprot": None, "pfam": [], "domains_at_residue": [], "features_at_residue": []}
        row = self.store.uniprot_for_gene(symbol) if symbol else None
        if row is None:
            return block
        block["uniprot"] = row["accession"]
        block["pfam"] = [{"accession": p["pfam"], "name": p["name"]} for p in self.store.uniprot_pfam(row["accession"])]
        residue = chosen.aa_position if chosen else None
        short = chosen.protein_short if chosen else None
        if not residue or not short:
            return block
        reference = row["sequence"]
        ref_aa = short[0]
        if residue > len(reference) or (ref_aa != "*" and reference[residue - 1] != ref_aa):
            flags.append("uniprot_numbering_differs")
            return block
        block["domains_at_residue"] = [
            {"accession": r["pfam"], "name": r["name"], "start": r["start"], "end": r["end"]}
            for r in self.store.pfam_regions_at(row["accession"], residue)
        ]
        block["features_at_residue"] = [
            {"type": f["type"], "start": f["start"], "end": f["end"], "note": f["note"]}
            for f in self.store.uniprot_features_at(row["accession"], residue)
        ]
        return block

    def _evidence(self, variant, lifted, row, clinvar_variant, symbol, chosen, flags) -> dict:
        versions = self._versions
        version = versions.get(f"clinvar.{variant.build}", "")
        evidence: dict = {"clinvar": {"status": "not_found", "version": version}, "dbsnp": None}
        if not self.store.has(f"clinvar.{variant.build}"):
            evidence["clinvar"] = {"status": "not_installed", "version": None}
        if row is not None:
            names = self.store.clinvar_names(clinvar_variant.build, row["variation_id"])
            evidence["clinvar"] = {
                "status": "found",
                "id": row["variation_id"],
                "classification": _pretty(row["classification"]),
                "review_status": _pretty(row["review_status"]),
                "condition": _pretty(row["condition"]),
                "molecular_consequence": row["molecular_consequence"] or None,
                "origin": row["origin"] or None,
                "oncogenicity": _pretty(row["oncogenicity"]),
                "somatic_impact": _pretty(row["somatic_impact"]),
                "frequencies": {
                    key: row[column]
                    for key, column in (("exac", "af_exac"), ("esp", "af_esp"), ("tgp", "af_tgp"))
                    if row[column] is not None
                },
                "names": [{"transcript": n["transcript"], "c": n["hgvs_c"], "p": n["hgvs_p"]} for n in names],
                "matched_on": clinvar_variant.build,
                "version": versions.get(f"clinvar.{clinvar_variant.build}", ""),
            }
            evidence["dbsnp"] = f"rs{row['rsid']}" if row["rsid"] else None

        short = chosen.protein_short if chosen else None
        evidence["civic"] = [
            {"variant": c["variant"], "url": c["url"], "evidence_items": c["evidence_items"]}
            for c in self.store.civic_variants(symbol, short or "")
        ]
        evidence["clingen_validity"] = [
            {"disease": v["disease"], "moi": v["moi"], "classification": v["classification"], "url": v["report_url"]}
            for v in self.store.gene_rows("clingen_validity", symbol)
        ]
        dosage = self.store.gene_rows("clingen_dosage", symbol)
        evidence["clingen_dosage"] = (
            {"haploinsufficiency": dosage[0]["haplo_description"], "triplosensitivity": dosage[0]["triplo_description"]}
            if dosage
            else None
        )
        constraint = self.store.gene_rows("constraint_metrics", symbol)
        evidence["constraint"] = (
            {"pli": constraint[0]["pli"], "loeuf": constraint[0]["loeuf"], "transcript": constraint[0]["transcript"], "source": "gnomAD v2.1.1"}
            if constraint
            else None
        )
        evidence["gene_diseases"] = sorted({d["disease_id"] for d in self.store.gene_rows("gene_diseases", symbol)})

        evidence["advanced"] = self._advanced(variant, lifted, symbol, chosen, row)
        dbsnp = evidence["advanced"].get("dbsnp") or {}
        if not evidence["dbsnp"] and dbsnp.get("rsid"):
            evidence["dbsnp"] = f"rs{dbsnp['rsid']}"
        return evidence

    def _advanced(self, variant, lifted, symbol, chosen, clinvar_row) -> dict:
        """Look up each advanced database installed here. The core results don't depend on any of them."""
        found: dict = {}
        for database in self._installed_advanced:
            if database.level == "variant":
                found[database.name] = self._advanced_variant(database, variant, lifted)
            elif database.level == "position":
                found[database.name] = self._advanced_position(variant, lifted)
            elif database.level == "region":
                found[database.name] = self._advanced_region(variant)
            elif database.level == "residue":
                found[database.name] = self._advanced_residue(symbol, chosen)
            elif database.level == "clinvar_id":
                hit = self.store.advanced_clingen(clinvar_row["variation_id"]) if clinvar_row is not None else None
                found[database.name] = {key: hit[key] for key in hit.keys()} if hit else {"status": "not_found"}
            else:
                found[database.name] = self._advanced_gene(database.name, symbol)
        return found

    def _advanced_variant(self, database, variant, lifted) -> dict:
        for candidate in (variant, lifted):
            if candidate is None:
                continue
            hit = self.store.advanced_variant(database.name, candidate.build, candidate.contig, candidate.pos, candidate.ref, candidate.alt)
            if hit is not None:
                values = {key: hit[key] for key in hit.keys() if key not in ("contig", "pos", "ref", "alt")}
                return {"status": "found", "matched_on": candidate.build, **values}
        return {"status": "not_found"}

    def _advanced_position(self, variant, lifted) -> dict:
        for candidate in (variant, lifted):
            if candidate is None or candidate.contig == "chrMT":
                continue
            start, _, _ = candidate.edit()
            value = self.store.phylop(candidate.build, candidate.contig, start)
            if value is not None:
                return {"status": "found", "value": round(value, 3), "matched_on": candidate.build}
        return {"status": "not_found"}

    def _advanced_region(self, variant) -> dict:
        rows = self.store.advanced_sv(variant.build, variant.contig, variant.pos)
        if not rows:
            return {"status": "not_found", "count": 0}
        top = rows[0]
        return {"status": "found", "count": len(rows), "top": f"{top['sv_type']} {top['sv_id']} AF {top['af']}"}

    def _advanced_residue(self, symbol, chosen) -> dict:
        if not (symbol and chosen and chosen.aa_position and chosen.protein_short):
            return {"status": "not_found"}
        in_frame = any(term.startswith("inframe_") for term in chosen.terms)
        for hit in self.store.advanced_hotspots(symbol, chosen.aa_position):
            indel_region = hit["type"] == "in-frame indel"
            if indel_region and not in_frame:
                continue
            # The list uses its own transcripts, so a single residue only counts when its letter agrees.
            if not indel_region and (hit["residue"] or "")[:1] != chosen.protein_short[0]:
                continue
            kind = " in-frame indels" if indel_region else ""
            return {
                "status": "found",
                "label": f"{symbol} {hit['residue']}{kind}: {hit['tumor_count']} tumours",
                "tumor_count": hit["tumor_count"],
                "variants": hit["variants"],
            }
        return {"status": "not_found"}

    def _advanced_gene(self, name: str, symbol: str) -> dict:
        if not symbol:
            return {"status": "not_found"}
        if name == "decipher":
            gene = self.store.advanced_gene("decipher", "genes", symbol)
            panel = self.store.advanced_gene("decipher", "ddg2p", symbol)
            if not gene and not panel:
                return {"status": "not_found"}
            return {
                "status": "found",
                "hi_percent": gene[0]["hi_percent"] if gene else None,
                "ddg2p": "; ".join(f"{row['disease']} ({row['inheritance']}, {row['confidence']})" for row in panel) or None,
            }
        if name == "cosmic_census":
            rows = self.store.advanced_gene(name, "genes", symbol)
            if not rows:
                return {"status": "not_found"}
            return {"status": "found", "tier": rows[0]["tier"], "role": rows[0]["role"], "tumour_types": rows[0]["tumour_types"]}
        if name == "cosmic_fusions":
            rows = sorted(self.store.advanced_gene(name, "fusions", symbol), key=lambda row: -row["samples"])
            if not rows:
                return {"status": "not_found"}
            return {"status": "found", "partners": "; ".join(f"{row['partner']} ({row['samples']})" for row in rows[:8])}
        rows = self.store.advanced_gene(name, "genes", symbol)
        if not rows:
            return {"status": "not_found"}
        return {"status": "found", "mim": rows[0]["mim"], "phenotypes": "; ".join(row["phenotypes"] for row in rows if row["phenotypes"]) or None}


def _pretty(value) -> str | None:
    if not value:
        return None
    return str(value).replace("_", " ").replace(" ,", ",")


def _other(effect: TranscriptEffect) -> dict:
    return {
        "accession": effect.transcript.accession,
        "gene": effect.transcript.gene,
        "hgvs_c": effect.hgvs_c,
        "hgvs_p": effect.hgvs_p,
        "effect": effect.effect,
        "flags": effect.flags,
    }


def _link_values(variant: Variant, record: dict, clinvar_row) -> dict:
    gene = record["gene"]
    protein = record["protein"]
    transcript = record["transcript"] or {}
    short = record["hgvs"]["protein_short"]
    if short and any(mark in short for mark in ("=", "?")):
        short = None
    pfam = protein["domains_at_residue"] or protein["pfam"]
    civic = record["evidence"].get("civic") or []
    symbol = gene["symbol"]
    query = f"{symbol} {short}" if symbol and short else (symbol or "")
    return {
        "build": variant.build,
        "contig": "chrM" if variant.contig == "chrMT" else variant.contig,
        "contig_chr": f"chr{key_contig(variant.contig)}",
        "contig_bare": key_contig(variant.contig),
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "gene": symbol,
        "protein_short": short,
        "clinvar_id": clinvar_row["variation_id"] if clinvar_row is not None else None,
        "rsid": (record["evidence"].get("dbsnp") or "").removeprefix("rs") or None,
        "ncbi_gene": gene["ncbi_gene"],
        "hgnc_id": gene["hgnc"],
        "ensembl_gene": gene["ensembl"],
        "omim": gene["omim"],
        "uniprot": protein["uniprot"],
        "pfam": pfam[0]["accession"] if pfam else None,
        "refseq": transcript.get("accession") if transcript.get("source") == "RefSeq" else None,
        "civic_url": civic[0]["url"] if civic else None,
        "cosmic_query": query,
        "pubmed_query": query,
    }

