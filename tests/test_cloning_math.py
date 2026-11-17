"""
The ligation arithmetic.

The molar mass of a base pair is derived here rather than typed, so the first
test is that the derivation lands where the chemistry says it should. The
rest is the ratio, which is a count of molecules and not of nanograms: that
distinction is the whole reason this view exists.
"""

from __future__ import annotations

import pytest

import cloning as cl


def test_the_base_pair_mass_falls_out_of_the_formulas():
    """617.9 g/mol, built from base + deoxyribose - water + HPO3 per strand.

    The figure quoted everywhere is 617.96 and circulates without a primary
    source. This is the same number arrived at from the four deoxynucleotide
    monophosphates with the CIAAW atomic weights the buffers section uses, so
    anyone can check it rather than taking it.
    """
    from buffers import molar_mass

    assert cl.base_pair_mass() == pytest.approx(617.9, abs=0.1)

    water = molar_mass("H2O")
    at = (molar_mass("C10H14N5O6P") + molar_mass("C10H15N2O8P") - 2 * water)
    gc = (molar_mass("C10H14N5O7P") + molar_mass("C9H14N3O7P") - 2 * water)
    assert cl.base_pair_mass() == pytest.approx((at + gc) / 2)
    # The spread the screen says is too small to ask the GC content for.
    assert abs(at - gc) / cl.base_pair_mass() < 0.002


def test_a_duplex_adds_one_water_per_strand_for_its_ends():
    from buffers import molar_mass

    assert cl.duplex_mass(1000) == pytest.approx(
        1000 * cl.base_pair_mass() + 2 * molar_mass("H2O"))


def test_picomoles_from_a_mass_and_a_length():
    """100 ng of a 3000 bp fragment is about 0.054 pmol."""
    assert cl.dna_pmol(100.0, 3000) == pytest.approx(0.0539, abs=1e-3)
    # Twice the length, half the moles.
    assert cl.dna_pmol(100.0, 6000) == pytest.approx(
        cl.dna_pmol(100.0, 3000) / 2, rel=1e-3)


def test_the_ratio_counts_molecules_and_not_nanograms():
    """50 ng of a 5000 bp vector at 3:1 with a 500 bp insert is 15 ng.

    Three molecules per molecule, at a tenth of the length, is three tenths
    of the mass. Doing the ratio in nanograms would have put 150 ng in, ten
    times too much, and nothing downstream would have said so.
    """
    plan = cl.ligation_plan(50.0, 5000, 500, 20.0, molar_ratio=3.0)

    assert plan["insert_pmol"] == pytest.approx(3 * plan["vector_pmol"])
    assert plan["insert_ng"] == pytest.approx(15.0, abs=0.1)
    assert plan["mass_ratio"] == pytest.approx(0.3, abs=0.01)
    assert plan["insert_ul"] == pytest.approx(0.75, abs=0.01)


def test_the_mass_ratio_is_spelled_out_when_the_lengths_differ():
    plan = cl.ligation_plan(50.0, 5000, 500, 20.0)
    assert any("per ng of vector" in note for note in plan["notes"])


def test_equal_lengths_make_the_two_ratios_the_same():
    plan = cl.ligation_plan(50.0, 3000, 3000, 20.0, molar_ratio=3.0)
    assert plan["mass_ratio"] == pytest.approx(3.0, rel=1e-6)


def test_an_unpipettable_volume_quotes_the_minimum_from_bench():
    import bench

    plan = cl.ligation_plan(50.0, 5000, 500, 20.0)
    assert any(f"{bench.PIPETTE_MINIMUM_UL:.0f} uL" in note
               for note in plan["notes"])


@pytest.mark.parametrize("kwargs", [
    {"vector_ng": 0.0},
    {"insert_ng_per_ul": 0.0},
    {"molar_ratio": 0.0},
    {"vector_bp": 0},
    {"insert_bp": -5},
])
def test_a_ligation_refuses_what_it_cannot_compute(kwargs):
    args = {"vector_ng": 50.0, "vector_bp": 5000, "insert_bp": 500,
            "insert_ng_per_ul": 20.0}
    args.update(kwargs)
    with pytest.raises(ValueError):
        cl.ligation_plan(**args)
