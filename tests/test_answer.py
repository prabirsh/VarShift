"""Sorting the candidates for one query into one answer, and laying the answer out for a terminal."""

import copy
import re
import unittest

from tests.test_integration import WithData
from varshift.answer import Answer, KEY_LINKS, also_text, choice_table, decide, render_answer, render_choices, summary_line
from varshift.flags import FLAG_HELP
from varshift.links import link_names
from varshift.term import QuietTerminal, Terminal


def record(key, in_clinvar=False, bases=1, kind="protein", mane=None, c="NM_004333.6:c.1799T>A"):
    """A small record with only the fields decide() reads."""
    return {
        "input": {"raw": "BRAF V600E", "type": kind, "build": "hg19", "read_as": ""},
        "candidate": {"in_clinvar": in_clinvar, "bases_changed": bases, "best": False, "why": []},
        "variant": {"hg19": {"key": key}},
        "transcript": {"accession": "NM_004333.6", "mane": mane} if mane is not None else None,
        "hgvs": {"c": c},
    }


class Rules(unittest.TestCase):
    def test_one_candidate_is_single(self):
        answer = decide([record("a")])
        self.assertEqual(answer.kind, "single")
        self.assertEqual(answer.best["variant"]["hg19"]["key"], "a")
        self.assertEqual(answer.why, [])
        self.assertEqual(answer.also, [])

    def test_the_raw_text_and_reading_come_from_the_input(self):
        first = record("a")
        first["input"]["read_as"] = "BRAF Val600Glu"
        answer = decide([first])
        self.assertEqual((answer.raw, answer.read_as), ("BRAF V600E", "BRAF Val600Glu"))

    def test_an_rs_number_with_several_alleles_is_a_choice(self):
        answer = decide([record("a", True, kind="rsid"), record("b", False, kind="rsid")])
        self.assertEqual(answer.kind, "choice")
        self.assertIsNone(answer.best)
        self.assertEqual(answer.why, [])

    def test_the_only_candidate_in_clinvar_is_clear(self):
        answer = decide([record("a", True, 1), record("b", False, 2)])
        self.assertEqual(answer.kind, "clear")
        self.assertEqual(answer.best["variant"]["hg19"]["key"], "a")
        self.assertEqual(answer.why, ["in ClinVar; the others aren't"])
        self.assertEqual([r["variant"]["hg19"]["key"] for r in answer.also], ["b"])

    def test_the_best_comes_first_even_when_it_was_ranked_lower(self):
        answer = decide([record("a", False, 1), record("b", True, 2), record("c", False, 2)])
        self.assertEqual(answer.best["variant"]["hg19"]["key"], "b")
        self.assertEqual([r["variant"]["hg19"]["key"] for r in answer.records], ["b", "a", "c"])
        self.assertEqual([r["variant"]["hg19"]["key"] for r in answer.also], ["a", "c"])

    def test_several_candidates_in_clinvar_are_a_choice(self):
        answer = decide([record("a", True), record("b", True), record("c", False)])
        self.assertEqual(answer.kind, "choice")
        self.assertIsNone(answer.best)

    def test_none_in_clinvar_and_one_needs_fewest_bases_is_clear(self):
        answer = decide([record("a", False, 1), record("b", False, 3), record("c", False, 2)])
        self.assertEqual(answer.kind, "clear")
        self.assertEqual(answer.best["variant"]["hg19"]["key"], "a")
        self.assertEqual(answer.why, ["fewest bases changed (1; the others need 2 or more)"])

    def test_a_tie_on_bases_changed_is_a_choice(self):
        answer = decide([record("a", False, 1), record("b", False, 1), record("c", False, 2)])
        self.assertEqual(answer.kind, "choice")
        self.assertIsNone(answer.best)

    def test_a_mane_select_transcript_is_named_in_why(self):
        answer = decide([record("a", True, mane="MANE Select"), record("b")])
        self.assertEqual(answer.why, ["in ClinVar; the others aren't", "MANE Select transcript"])
        plus_clinical = decide([record("a", False, 1, mane="MANE Plus Clinical"), record("b", False, 2)])
        self.assertNotIn("MANE Select transcript", plus_clinical.why)

    def test_a_single_answer_needs_no_reason(self):
        self.assertEqual(decide([record("a", True, mane="MANE Select")]).why, [])

    def test_nothing_to_decide_is_a_mistake(self):
        with self.assertRaises(ValueError):
            decide([])


class AlsoText(unittest.TestCase):
    def test_one_other_candidate(self):
        answer = decide([record("a", True), record("b", False, 2, c="NM_004333.6:c.1799_1800delinsAA")])
        self.assertEqual(also_text(answer), "c.1799_1800delinsAA (2 bases, not in ClinVar)")

    def test_one_base_and_in_clinvar(self):
        answer = decide([record("a", False, 1, kind="coding"), record("b", False, 3)])
        answer.also[0]["candidate"]["in_clinvar"] = True
        answer.also[0]["candidate"]["bases_changed"] = 1
        self.assertEqual(also_text(answer), "c.1799T>A (1 base, in ClinVar)")

    def test_up_to_three_then_how_many_more(self):
        others = [record(f"k{n}", False, 2 + n, c=f"NM_1:c.{n}A>T") for n in range(5)]
        answer = decide([record("a", True)] + others)
        self.assertEqual(
            also_text(answer),
            "c.0A>T (2 bases, not in ClinVar); c.1A>T (3 bases, not in ClinVar); c.2A>T (4 bases, not in ClinVar); and 2 more",
        )

    def test_a_candidate_with_no_c_name_is_given_by_its_position(self):
        answer = decide([record("a", True), record("chr7-1-A-T", False, 1, c=None)])
        self.assertEqual(also_text(answer), "chr7-1-A-T (1 base, not in ClinVar)")

    def test_nothing_else_gives_no_text(self):
        self.assertEqual(also_text(decide([record("a")])), "")
        self.assertEqual(also_text(decide([record("a", True), record("b", True)])), "")


class OnTheDemoData(WithData):
    def answer_of(self, text):
        records = self.annotator.annotate_text(text, "hg19")
        return records, {r["input"]["answer"] for r in records}

    def test_braf_v600e_is_clear_and_says_why(self):
        records, kinds = self.answer_of("BRAF V600E")
        self.assertEqual(kinds, {"clear"})
        best = [r for r in records if r["candidate"]["best"]]
        self.assertEqual(len(best), 1)
        self.assertIs(best[0], records[0])
        self.assertEqual(best[0]["variant"]["hg19"]["key"], "chr7-140453136-A-T")
        self.assertTrue(best[0]["candidate"]["why"][0].startswith("in ClinVar"))
        self.assertIn("c.1799_1800delinsAA (2 bases, not in ClinVar)", best[0]["also"])
        self.assertEqual([r for r in records[1:] if "also" in r or r["candidate"]["why"]], [])
        self.assertEqual([r["candidate"]["rank"] for r in records], list(range(1, len(records) + 1)))

    def test_an_rs_number_with_several_alleles_is_a_choice(self):
        records, kinds = self.answer_of("rs113488022")
        self.assertEqual(kinds, {"choice"})
        self.assertFalse(any(r["candidate"]["best"] for r in records))
        self.assertFalse(any("also" in r for r in records))

    def test_a_c_name_is_single(self):
        records, kinds = self.answer_of("NM_004333.6:c.1799T>A")
        self.assertEqual(kinds, {"single"})
        self.assertEqual([r["candidate"]["best"] for r in records], [True])
        self.assertNotIn("also", records[0])

    def test_the_answer_reaches_the_columns(self):
        from varshift.output import flatten

        rows = [flatten(r) for r in self.annotator.annotate_text("BRAF V600E", "hg19")]
        self.assertEqual([row["best"] for row in rows[:2]], ["yes", "no"])
        self.assertEqual(rows[0]["answer"], "clear")
        self.assertTrue(rows[0]["why"].startswith("in ClinVar"))
        self.assertTrue(rows[0]["also"])
        self.assertIsNone(rows[1]["also"])


def shown_links(text):
    """The link names on the card's Check it lines, in order."""
    names = sorted(link_names(), key=len, reverse=True)
    shown = []
    for line in text.splitlines():
        if "https://" not in line:
            continue
        rest = line.strip().removeprefix("Check it").strip()
        shown.append(next(name for name in names if rest.startswith(name)))
    return shown


class Card(WithData):
    def card(self, text, mode=None, links_all=False, color=False, width=100):
        records = self.annotator.annotate_text(text, "hg19")
        term = Terminal(interactive=False, color=True, width=width) if color else QuietTerminal(width=width)
        return render_answer(decide(records), term, mode=mode, links_all=links_all)

    def test_braf_v600e_fits_in_twenty_lines_and_says_why(self):
        text = self.card("BRAF V600E")
        lines = text.splitlines()
        self.assertLessEqual(len(lines), 20, text)
        self.assertTrue(lines[0].startswith("BRAF V600E  ->  chr7-140453136-A-T (hg19)"), lines[0])
        for expected in ("Why this", "Also", "c.1799_1800delinsAA", "show=all", "15 more: links=all", "HIGH", "c.1799T>A", "p.Val600Glu"):
            self.assertIn(expected, text)
        for code in FLAG_HELP:
            self.assertNotIn(code, text)

    def test_lines_other_than_links_stay_inside_the_width(self):
        for width in (60, 100):
            with self.subTest(width=width):
                for line in self.card("BRAF V600E", width=width).splitlines():
                    if "https://" not in line:
                        self.assertLessEqual(len(line), width, line)

    def test_links_all_shows_every_link_with_its_label(self):
        text = self.card("BRAF V600E", links_all=True)
        for name in link_names():
            self.assertIn(name, text)
        self.assertIn("confirmed ID", text)
        self.assertIn("COSMIC may ask you to sign in first.", text)
        self.assertNotIn("links=all", text)

    def test_the_mode_chooses_the_links(self):
        self.assertEqual(shown_links(self.card("BRAF V600E")), ["ClinVar", "gnomAD", "UCSC Genome Browser", "dbSNP", "VarSome"])
        tumour = shown_links(self.card("BRAF V600E", mode="tumour"))
        self.assertEqual(len(tumour), 5)
        self.assertTrue({"OncoKB", "CIViC"} <= set(tumour), tumour)
        inherited = shown_links(self.card("BRAF V600E", mode="inherited"))
        self.assertTrue({"ClinGen", "gnomAD"} <= set(inherited), inherited)
        for mode, names in KEY_LINKS.items():
            self.assertTrue(set(names) <= set(link_names()), mode)
            self.assertTrue(set(shown_links(self.card("BRAF V600E", mode=mode))) <= set(names), mode)

    def test_the_mode_recorded_with_the_query_is_used(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19", mode="tumour")
        self.assertIn("OncoKB", shown_links(render_answer(decide(records), QuietTerminal())))

    def test_the_mode_orders_the_evidence(self):
        def order(text):
            labels = [line[2:14].strip() for line in text.splitlines() if line.startswith("  ") and line[2:3] != " "]
            return [label for label in labels if label in ("ClinVar", "CIViC", "Domain", "ClinGen", "Constraint")]

        self.assertEqual(order(self.card("BRAF V600E")), ["ClinVar", "CIViC", "Domain"])
        self.assertEqual(order(self.card("BRAF V600E", mode="tumour"))[:2], ["CIViC", "ClinVar"])
        self.assertEqual(order(self.card("BRAF V600E", mode="inherited"))[0], "ClinVar")
        self.assertEqual(order(self.card("BRAF V600E", mode="inherited"))[-1], "Domain")

    def test_an_unknown_mode_is_a_mistake(self):
        with self.assertRaises(ValueError):
            self.card("BRAF V600E", mode="bogus")

    def test_colour_only_when_asked(self):
        self.assertIn("\x1b[1;31mHIGH", self.card("BRAF V600E", color=True))
        self.assertNotIn("\x1b", self.card("BRAF V600E"))

    def test_the_reading_is_shown(self):
        text = self.card("braf val600glu")
        self.assertIn("Read as", text)
        self.assertIn("BRAF Val600Glu", text)
        self.assertIn("amino acids written in standard case", text)
        self.assertNotIn("Read as", self.card("BRAF V600E"))

    def test_a_single_answer_has_no_why_or_also(self):
        text = self.card("NM_004333.6:c.1799T>A")
        self.assertNotIn("Why this", text)
        self.assertNotIn("Also", text)
        self.assertIn("chr7-140453136-A-T", text.splitlines()[0])

    def test_warnings_are_written_as_sentences(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19")
        records[0]["flags"] += ["other_transcript_used", "liftover_data_not_installed"]
        text = render_answer(decide(records), QuietTerminal())
        self.assertIn(FLAG_HELP["other_transcript_used"][0], text)
        self.assertNotIn(FLAG_HELP["other_transcript_used"][1], text)
        self.assertIn("varshift setup build=both", text)
        self.assertNotIn("other_transcript_used", text)
        self.assertNotIn(FLAG_HELP["multiple_candidates"][0], text)


def fake_choice(record, count):
    records = []
    for number in range(count):
        copied = copy.deepcopy(record)
        copied["variant"]["hg19"]["key"] = f"chr7-{140453000 + number}-A-T"
        records.append(copied)
    return Answer("BRAF V600E", "", "choice", records)


class Choices(WithData):
    def test_an_rs_number_lists_each_allele(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        text = render_choices(decide(records), QuietTerminal())
        lines = text.splitlines()
        self.assertEqual(lines[0], f"rs113488022 covers {len(records)} different variants in ClinVar:")
        rows = [line for line in lines if re.match(r"^\s+\d+\s{2}\S", line)]
        self.assertEqual(len(rows), len(records))
        for record, row in zip(records, rows):
            self.assertIn(record["variant"]["hg19"]["key"], row)
            self.assertIn("BRAF", row)

    def test_a_protein_change_with_several_dna_changes(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19")
        answer = Answer("BRAF V600E", "", "choice", records)
        self.assertEqual(render_choices(answer, QuietTerminal()).splitlines()[0], f"BRAF V600E fits {len(records)} DNA changes:")

    def test_rows_are_cut_at_the_width_and_keep_the_priority(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19")
        for line in render_choices(Answer("BRAF V600E", "", "choice", records), QuietTerminal(width=60)).splitlines()[1:]:
            self.assertLessEqual(len(line), 60, line)
            self.assertTrue(line.rstrip().endswith(("HIGH", "MEDIUM", "REVIEW", "LOW")), line)

    def test_a_narrow_table_drops_columns_before_cutting_the_key(self):
        records = self.annotator.annotate_text("BRAF V600E", "hg19")
        lines = render_choices(Answer("BRAF V600E", "", "choice", records), QuietTerminal(width=60)).splitlines()
        for record, line in zip(records, lines[1:]):
            self.assertIn(record["variant"]["hg19"]["key"], line)
        wide = render_choices(Answer("BRAF V600E", "", "choice", records), QuietTerminal(width=120))
        self.assertIn("Conflicting", wide)

    def test_at_most_twenty_rows(self):
        record = self.annotator.annotate_text("BRAF V600E", "hg19")[0]
        lines = render_choices(fake_choice(record, 25), QuietTerminal()).splitlines()
        self.assertEqual(len([line for line in lines if re.match(r"^\s+\d+\s{2}\S", line)]), 20)
        self.assertEqual(lines[-1].strip(), "5 more: be more specific")

    def test_the_table_can_be_offered_as_options(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        heading, rows = choice_table(decide(records), QuietTerminal())
        self.assertTrue(heading.startswith("rs113488022 covers"))
        self.assertEqual(len(rows), len(records))
        self.assertFalse(any(re.match(r"^\s*\d", row) for row in rows))

    def test_a_card_asked_for_a_choice_shows_the_table(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        self.assertIn("covers", render_answer(decide(records), QuietTerminal()))


class Summary(WithData):
    def test_one_line_per_query(self):
        answer = decide(self.annotator.annotate_text("BRAF V600E", "hg19"))
        line = summary_line(1, answer, QuietTerminal())
        self.assertNotIn("\n", line)
        self.assertLessEqual(len(line), 100)
        for expected in ("1", "BRAF V600E", "chr7-140453136-A-T", "HIGH", "Conflicting", "CIViC"):
            self.assertIn(expected, line)

    def test_a_choice_says_how_many(self):
        records = self.annotator.annotate_text("rs113488022", "hg19")
        line = summary_line(2, decide(records), QuietTerminal())
        self.assertIn(f"{len(records)} choices", line)
        self.assertTrue(line.lstrip().startswith("2"))

    def test_a_narrow_terminal_cuts_the_line(self):
        answer = decide(self.annotator.annotate_text("BRAF V600E", "hg19"))
        self.assertLessEqual(len(summary_line(12, answer, QuietTerminal(width=60))), 60)
        self.assertLessEqual(len(summary_line(12, answer, QuietTerminal(width=60), input_width=40)), 60)

    def test_lines_for_several_queries_line_up(self):
        texts = ("BRAF V600E", "NM_004333.6:c.1799T>A")
        answers = [decide(self.annotator.annotate_text(text, "hg19")) for text in texts]
        width = max(len(text) for text in texts)
        lines = [summary_line(n, answer, QuietTerminal(), input_width=width) for n, answer in enumerate(answers, 1)]
        self.assertEqual(len({line.index("chr7-140453136-A-T") for line in lines}), 1, lines)

    def test_colour_marks_the_priority(self):
        answer = decide(self.annotator.annotate_text("BRAF V600E", "hg19"))
        self.assertIn("\x1b[1;31mHIGH", summary_line(1, answer, Terminal(interactive=False, color=True, width=100)))


if __name__ == "__main__":
    unittest.main()
