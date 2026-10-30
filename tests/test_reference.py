"""
The reference tables, and the one property that makes them worth having.

A reference table typed in by hand can disagree with the app printed next to
it, and when it does the reader has no way of telling which of the two is
wrong. `reference.py` exists to make that impossible: every column is read out
of the module the rest of the app reads it out of. These tests are what holds
that claim up, so most of them are comparisons against the source rather than
against a literal.

The structure drawings are checked here too, lightly: the app only reads them,
and `tools/aminoacids/build.py` is what verifies the chemistry when they are
made. What this file asks is that the twenty files are present, that they are
not empty, and that they still carry the two things the screen depends on,
a viewBox and `currentColor`.
"""

from __future__ import annotations

import math

import pytest

import codons
import reference
import vhh_core

ALL_TWENTY = "ACDEFGHIKLMNPQRSTVWY"


class TestAminoAcids:

    def test_there_are_twenty(self):
        table = reference.aminoacids()
        assert sorted(table["aa"]) == sorted(ALL_TWENTY)

    def test_the_residue_masses_build_the_app_s_molecular_weights(self):
        """The strongest statement this module makes.

        If the residue masses in the table and the masses behind every
        molecular weight the app reports are not the same numbers, one screen
        is lying to the reader. Summing the table across a chain and comparing
        against ProtParam is the cheapest way to say so.
        """
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        masses = reference.aminoacids().set_index("aa")["residue_da"]
        total = sum(masses[residue] for residue in ALL_TWENTY) + 18.01528
        assert total == pytest.approx(
            ProteinAnalysis(ALL_TWENTY).molecular_weight(), abs=0.01)

    def test_hydropathy_is_the_scale_the_profiles_use(self):
        table = reference.aminoacids().set_index("aa")
        for residue in ALL_TWENTY:
            assert table.loc[residue, "hydropathy"] == vhh_core.KD[residue]

    def test_the_charge_column_is_the_profile_s_approximation(self):
        """Including histidine at +0.1, which is the one that looks like a
        typo and is not: it is the app's own compromise between treating
        histidine as protonated and as neutral."""
        table = reference.aminoacids().set_index("aa")
        assert table.loc["H", "charge_profile"] == pytest.approx(0.1)
        for residue in "DEKR":
            assert (table.loc[residue, "charge_profile"]
                    == vhh_core.CHARGE[residue])

    def test_the_drawn_charges_are_the_ones_pH_implies(self):
        """Aspartate and glutamate negative, lysine and arginine positive,
        histidine neutral, everything else neutral. This is what makes the
        pictures specific to a pH rather than decorative."""
        table = reference.aminoacids().set_index("aa")
        expected = {"D": -1, "E": -1, "K": 1, "R": 1}
        for residue in ALL_TWENTY:
            assert table.loc[residue, "charge_at_ph"] == expected.get(residue, 0)

    def test_every_codon_column_agrees_with_the_codon_table(self):
        table = reference.aminoacids().set_index("aa")
        for residue in ALL_TWENTY:
            listed = table.loc[residue, "codons"].split()
            assert listed, f"{residue} has no codons"
            for codon in listed:
                assert reference._table().forward_table[codon] == residue


class TestStructures:

    def test_all_twenty_are_present_and_not_empty(self):
        """They have been zero bytes before, written by a post-processing step
        that opened each file for writing before reading it. Nothing raised;
        the page simply drew nothing."""
        for residue in ALL_TWENTY:
            svg = reference.structure_svg(residue)
            assert svg, f"no drawing for {residue}"
            assert len(svg) > 1000, f"{residue}: {len(svg)} bytes, too short"

    def test_they_follow_the_theme_and_the_container(self):
        for residue in ALL_TWENTY:
            svg = reference.structure_svg(residue)
            assert "currentColor" in svg, \
                f"{residue}: the carbon skeleton is not themeable"
            assert "viewBox" in svg, f"{residue}: no viewBox, it cannot scale"

    def test_the_pH_is_recorded_with_them(self):
        """A drawing of a charged species without its pH says nothing."""
        assert reference.structure_ph() == pytest.approx(7.4)

    def test_a_missing_drawing_is_not_an_exception(self):
        assert reference.structure_svg("X") is None


class TestGeneticCode:

    def test_all_sixty_four_codons(self):
        table = reference.genetic_code()
        assert len(table) == 64
        assert len(set(table["codon"])) == 64

    def test_three_stops(self):
        table = reference.genetic_code()
        assert sorted(table.loc[table["aa"] == "*", "codon"]) == \
            ["TAA", "TAG", "TGA"]

    def test_the_rare_flags_are_the_gene_screen_s_own(self):
        """A codon flagged rare here and not there would be the exact failure
        this module is built to prevent."""
        table = reference.genetic_code().set_index("codon")
        flagged = set(table.index[table["rare_in_ecoli"]])
        assert flagged == set(codons.RARE_CODONS)

    def test_w_is_the_relative_adaptiveness_cai_uses(self):
        table = reference.genetic_code().set_index("codon")
        weights = codons.relative_adaptiveness()
        for codon, weight in weights.items():
            assert table.loc[codon, "w"] == pytest.approx(weight, abs=0.001)

    def test_the_best_synonym_has_w_of_one(self):
        table = reference.genetic_code()
        coding = table[table["aa"] != "*"]
        for amino, group in coding.groupby("aa"):
            assert group["w"].max() == pytest.approx(1.0, abs=0.001), \
                f"{amino}: no codon has w = 1"

    def test_shares_within_an_amino_acid_add_up(self):
        table = reference.genetic_code()
        coding = table[table["aa"] != "*"]
        for amino, group in coding.groupby("aa"):
            assert group["share"].sum() == pytest.approx(1.0, abs=0.02), \
                f"{amino}: shares sum to {group['share'].sum():.3f}"

    def test_stops_have_no_share(self):
        """They do not compete with a synonym for an amino acid, so a share
        would be a number about nothing."""
        table = reference.genetic_code()
        for value in table.loc[table["aa"] == "*", "share"]:
            assert math.isnan(value)


class TestCodonWheel:

    def test_it_is_the_textbook_block(self):
        wheel = reference.codon_wheel()
        assert wheel.shape == (16, 4)
        assert list(wheel.columns) == list("TCAG")

    def test_it_says_the_same_thing_as_the_flat_table(self):
        wheel = reference.codon_wheel()
        flat = reference.genetic_code().set_index("codon")["aa"]
        for pair, row in wheel.iterrows():
            for third, amino in row.items():
                assert flat[pair + third] == amino
