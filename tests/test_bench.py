"""
Bench tests.

These matter more than most of the suite, because what this module computes
goes into a tube. A wrong CDR annotation is noticed when it disagrees with
ANARCI; a molar excess that is wrong by a factor of ten is noticed weeks later,
in an experiment that did not work, and probably blamed on the chemistry.

The conjugation arithmetic is checked against the worked example in the Thermo
EZ-Link cheat sheet that labs actually use, so the test is a comparison with an
independent calculation rather than with my own.
"""

from __future__ import annotations

import pytest

import bench


# --------------------------------------------------------------------------
# Against the vendor worked example
# --------------------------------------------------------------------------

def test_conjugation_reproduces_the_thermo_cheat_sheet():
    """The EZ-Link NHS-PEG4-Biotin worked example, step by step.

    Its inputs: 0.6 mL of protein at 1.4 mg/mL, molar mass 14949.48 Da, a 4x
    molar excess of biotin, from a 2 mM stock. Its answers: 0.000224 mmol of
    biotin and 112.224 uL of stock.

    The sheet rounds 1/MW to six decimal places, which moves its answer by
    about 0.15%, so the comparison is to a tenth of a percent rather than
    exact. Where the two differ, the sheet is the one carrying the rounding.
    """
    plan = bench.conjugation_plan(
        protein_volume_ul=600, protein_mg_per_ml=1.4,
        protein_mw_da=14949.48, molar_excess=4, reagent_stock_mM=2.0)

    # The sheet carries two roundings: 1/MW to six decimal places, and its
    # displayed result to three significant figures. Its own underlying value
    # is 225.1 nmol and it prints 224; the exact answer is 224.76. A half
    # percent window covers the sheet, its display and the exact value.
    assert plan["reagent_nmol"] == pytest.approx(224.0, rel=5e-3)
    assert plan["reagent_ul"] == pytest.approx(112.224, rel=5e-3)
    assert plan["protein_ug"] == pytest.approx(840.0, rel=1e-6)


def test_stock_preparation_reproduces_the_vendor_protocol():
    """"Add 170 uL of water to the 2 mg vial to make 20 mM."

    That instruction is itself a check on the molar mass: it only gives 20 mM
    if the reagent is 588.7 g/mol, which is how the value in REAGENTS was
    corroborated against a second source.
    """
    prepared = bench.stock_preparation(vial_mg=2.0, reagent_mw=588.7,
                                       target_mM=20.0)
    assert prepared["solvent_ul"] == pytest.approx(170, abs=1.0)
    assert prepared["umol_in_vial"] == pytest.approx(3.4, abs=0.05)


def test_the_listed_reagents_carry_a_source():
    """A molar mass without provenance is a number somebody typed.

    A wrong one does not fail loudly: it silently changes the molar excess, and
    the result gets blamed on the chemistry.
    """
    assert bench.REAGENTS, "the reagent table is empty"
    for name, reagent in bench.REAGENTS.items():
        assert reagent.get("mw", 0) > 0, f"{name} has no molar mass"
        assert reagent.get("source"), f"{name} cites no source for its mass"


# --------------------------------------------------------------------------
# Conjugation, the parts the spreadsheet does not do
# --------------------------------------------------------------------------

def test_an_unpipettable_volume_is_turned_into_a_dilution():
    """The failure this is for: 0.2 uL asked for, 2 uL delivered.

    A 10x error in the molar excess, invisible, in the direction that looks
    like it worked.
    """
    plan = bench.conjugation_plan(20, 0.1, 14000, 2, 20.0)
    assert plan["reagent_ul"] < bench.PIPETTE_MINIMUM_UL
    assert plan["dilution"], "no dilution offered for an unpipettable volume"
    assert plan["dilution"]["factor"] >= 10
    assert plan["dilution"]["volume_ul"] >= bench.PIPETTE_MINIMUM_UL
    # The excess is unchanged by diluting: same nmol, larger volume.
    assert plan["dilution"]["volume_ul"] == pytest.approx(
        plan["reagent_ul"] * plan["dilution"]["factor"], rel=1e-6)
    assert any("1:" in note for note in plan["notes"])


def test_a_large_co_solvent_fraction_is_flagged():
    """DMSO past about 10% of the reaction is a protein problem, not a
    labelling one."""
    plan = bench.conjugation_plan(10, 5.0, 14000, 20, 1.0)
    assert plan["co_solvent_pct"] > 10
    assert any("DMSO" in note for note in plan["notes"])


def test_an_excess_beyond_the_available_amines_is_flagged():
    plan = bench.conjugation_plan(100, 1.0, 14000, 50, 20.0,
                                  n_reactive_amines=4)
    assert any("primary amines" in note for note in plan["notes"])


@pytest.mark.parametrize("kwargs", [
    {"protein_volume_ul": 0, "protein_mg_per_ml": 1, "protein_mw_da": 14000,
     "molar_excess": 4, "reagent_stock_mM": 2},
    {"protein_volume_ul": 100, "protein_mg_per_ml": -1, "protein_mw_da": 14000,
     "molar_excess": 4, "reagent_stock_mM": 2},
    {"protein_volume_ul": 100, "protein_mg_per_ml": 1, "protein_mw_da": 14000,
     "molar_excess": 4, "reagent_stock_mM": 0},
])
def test_conjugation_refuses_impossible_inputs(kwargs):
    with pytest.raises(ValueError):
        bench.conjugation_plan(**kwargs)


def test_desalting_column_follows_the_volume():
    assert "0.5 mL" in bench.desalting_column(100)
    assert "2 mL" in bench.desalting_column(400)
    assert "5 mL" in bench.desalting_column(2000)
    assert "below" in bench.desalting_column(10)


# --------------------------------------------------------------------------
# Quantification
# --------------------------------------------------------------------------

def test_beer_lambert_round_trips():
    """A known concentration, through the calculation and back."""
    mw, epsilon = 14000.0, 28000.0
    molar = 1e-5                                     # 10 uM
    a280 = molar * epsilon                           # 0.28 in a 1 cm cell

    out = bench.concentration_from_a280(a280, epsilon, mw)
    assert out["molar_uM"] == pytest.approx(10.0, rel=1e-6)
    assert out["mg_per_ml"] == pytest.approx(molar * mw, rel=1e-6)


def test_the_dilution_factor_multiplies_and_does_not_divide():
    """The direction that is easy to get backwards, and silently halves or
    doubles every number downstream."""
    neat = bench.concentration_from_a280(0.1, 28000, 14000, dilution_factor=1)
    diluted = bench.concentration_from_a280(0.1, 28000, 14000,
                                            dilution_factor=10)
    assert diluted["mg_per_ml"] == pytest.approx(neat["mg_per_ml"] * 10,
                                                 rel=1e-9)


def test_a_high_260_280_ratio_is_reported_as_an_upper_bound():
    out = bench.concentration_from_a280(0.5, 28000, 14000, a260=0.5)
    assert out["ratio_260_280"] == pytest.approx(1.0)
    assert out["nucleic_acid_warning"]

    clean = bench.concentration_from_a280(0.5, 28000, 14000, a260=0.285)
    assert clean["nucleic_acid_warning"] is None


def test_quantification_refuses_a_protein_with_no_chromophore():
    """A sequence with no Trp, Tyr or Cys has an extinction coefficient of
    zero, and dividing by it would return infinity rather than refusing."""
    with pytest.raises(ValueError):
        bench.concentration_from_a280(0.5, 0.0, 14000)


def test_dilution_cannot_concentrate():
    with pytest.raises(ValueError):
        bench.dilution_plan(stock_mg_per_ml=1.0, target_mg_per_ml=5.0,
                            target_volume_ul=100)


def test_dilution_volumes_add_up():
    plan = bench.dilution_plan(10.0, 2.0, 500)
    assert plan["stock_ul"] + plan["buffer_ul"] == pytest.approx(500)
    assert plan["stock_ul"] == pytest.approx(100)


# --------------------------------------------------------------------------
# Mass
# --------------------------------------------------------------------------

def test_a_mass_that_matches_offers_no_explanation():
    """Explaining a difference that is not there is how a tool teaches someone
    to distrust it."""
    seq = "MQVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYA"
    out = bench.mass_check(seq, measured_da=out_expected(seq))
    assert out["matches"]
    assert not len(out["candidates"])


def out_expected(seq: str) -> float:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return ProteinAnalysis(seq).molecular_weight()


def test_met_excision_is_offered_only_when_the_sequence_starts_with_met():
    with_met = "MQVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGW"
    expected = out_expected(with_met)
    out = bench.mass_check(with_met, measured_da=expected - 131.2)
    hits = out["candidates"]
    assert hits[hits["explains"]]["modification"].str.contains(
        "Met excised").any()

    # Matched on "Met excised", not on "Met": "Oxidation (per Met or Trp)" also
    # contains "Met", and the loose version of this test passed on that.
    without = with_met[1:]
    out2 = bench.mass_check(without, measured_da=out_expected(without) - 131.2)
    assert not out2["candidates"]["modification"].str.contains(
        "Met excised").any()


def test_the_his_tag_artefact_is_offered_only_to_tagged_sequences():
    """+178 on a His-tagged protein from E. coli is gluconoylation. On an
    untagged one it is something else, and saying otherwise would be inventing
    a diagnosis."""
    tagged = "MHHHHHHQVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGW"
    out = bench.mass_check(tagged, measured_da=out_expected(tagged) + 178.05)
    hits = out["candidates"]
    assert hits[hits["explains"]]["modification"].str.contains(
        "Gluconoylation").any()

    plain = "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGW"
    out2 = bench.mass_check(plain, measured_da=out_expected(plain) + 178.05)
    assert not out2["candidates"]["modification"].str.contains(
        "luconoylation").any()


def test_disulfide_formation_costs_two_hydrogens_per_bond():
    seq = "MQVQLQESGGGLVQAGGSLRLSCAASGRTFSSYACMGW"
    expected = out_expected(seq)
    out = bench.mass_check(seq, measured_da=expected - 2.02)
    hits = out["candidates"]
    row = hits[hits["modification"].str.contains("Disulfide")].iloc[0]
    assert row["delta_da"] == pytest.approx(-2.0157, abs=1e-3)
    assert row["explains"]


def test_modification_masses_come_from_their_formulas():
    """Every delta is computed from a chemical formula, so each one can be
    checked by hand. A tabulated number cannot."""
    for modification in bench.MODIFICATIONS:
        assert modification["formula"], f"{modification['name']} has no formula"
        assert modification["why"], f"{modification['name']} has no explanation"
        mass = bench._formula_mass(modification["formula"])
        assert abs(mass) > 1.0, f"{modification['name']} has a null mass delta"
