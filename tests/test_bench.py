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


# --------------------------------------------------------------------------
# Phage
# --------------------------------------------------------------------------

class TestPhageTitre:
    """The arithmetic a panning round is planned with.

    The reference numbers are a real stock measured in the lab: A260 0.459,
    diluted 100-fold, with a 2.77-fold glycerol correction, which the
    spreadsheet this reproduces turned into 3.82e12 particles per millilitre
    and a 26.2 uL aliquot for 1e11 particles. If this test fails, the
    calculator and the bench have stopped agreeing.
    """

    READING = dict(a260=0.459, dilution=100.0, glycerol_factor=2.77)

    def test_reproduces_the_spreadsheet(self):
        out = bench.phage_titre(**self.READING)
        assert out["empirical_per_ml"] == pytest.approx(3.82e12, rel=0.01)

    def test_the_aliquot_matches_too(self):
        out = bench.phage_titre(**self.READING)
        plan = bench.phage_aliquot(out["empirical_per_ml"], 1e11)
        assert plan["volume_ul"] == pytest.approx(26.2, rel=0.01)
        assert not plan["note"]

    def test_glycerol_added_after_the_reading_dilutes_the_tube(self):
        """The correction that was briefly impossible to enter.

        500 uL of phage plus 55.5 uL of pure glycerol is 555.5 uL, 10%
        glycerol, and the tube you pipette from is 1.111 times more dilute
        than the phage whose absorbance you read. So the factor is
        500/555.5 = 0.9, which is *below* one. The first version of the
        screen refused anything below one, which is the only direction this
        protocol ever goes.
        """
        factor = 500.0 / (500.0 + 55.5)
        assert factor == pytest.approx(0.9, abs=0.001)
        neat = bench.phage_titre(a260=0.5, dilution=100.0)["empirical_per_ml"]
        glycerol = bench.phage_titre(a260=0.5, dilution=100.0,
                                     glycerol_factor=factor)["empirical_per_ml"]
        assert glycerol == pytest.approx(neat * factor)
        assert glycerol < neat

    def test_the_two_routes_are_both_reported(self):
        """Neither estimate may quietly replace the other.

        They disagree by more than two orders of magnitude on the same
        reading, which is the single most important thing this screen has to
        say. A version that returned one number would be hiding it.
        """
        out = bench.phage_titre(**self.READING)
        assert out["optical_per_ml"] > 0
        assert out["ratio"] > 100

    def test_the_optical_route_is_days_relation(self):
        """A260 x 6e16 / bases, checked against the arithmetic by hand."""
        out = bench.phage_titre(a260=1.0, genome_bases=7249)
        assert out["optical_per_ml"] == pytest.approx(6e16 / 7249, rel=1e-9)

    def test_a_shorter_genome_gives_a_higher_titre(self):
        """Same absorbance, less DNA per particle, so more particles."""
        short = bench.phage_titre(a260=1.0, genome_bases=4000)
        long = bench.phage_titre(a260=1.0, genome_bases=8000)
        assert short["optical_per_ml"] > long["optical_per_ml"]

    def test_dilution_and_glycerol_both_multiply_up(self):
        plain = bench.phage_titre(a260=0.5)["empirical_per_ml"]
        both = bench.phage_titre(a260=0.5, dilution=10.0,
                              glycerol_factor=2.0)["empirical_per_ml"]
        assert both == pytest.approx(plain * 20)

    @pytest.mark.parametrize("bad", [
        dict(a260=-0.1),
        dict(a260=0.5, dilution=0.0),
        dict(a260=0.5, glycerol_factor=-1.0),
        dict(a260=0.5, per_a260=0.0),
        dict(a260=0.5, genome_bases=0),
    ])
    def test_impossible_inputs_are_refused(self, bad):
        with pytest.raises(ValueError):
            bench.phage_titre(**bad)


class TestPhageAliquot:

    def test_a_volume_too_small_to_pipette_is_called_out(self):
        """The same objection the conjugation planner makes, for the same
        reason: below a couple of microlitres the number is arithmetic rather
        than something you can put in a tube."""
        plan = bench.phage_aliquot(1e15, 1e9)
        assert plan["volume_ul"] < bench.PIPETTE_MINIMUM_UL
        assert "dilute" in plan["note"]

    def test_an_absurdly_large_volume_is_called_out(self):
        plan = bench.phage_aliquot(1e9, 1e13)
        assert plan["volume_ul"] > 1000
        assert plan["note"]

    @pytest.mark.parametrize("titre,target", [(0.0, 1e11), (1e12, 0.0),
                                              (-1.0, 1e11)])
    def test_impossible_inputs_are_refused(self, titre, target):
        with pytest.raises(ValueError):
            bench.phage_aliquot(titre, target)


class TestLiquidPanning:
    """A liquid-phase round takes twice the input of a solid-phase one."""

    def test_the_liquid_volume_is_double(self):
        titre = 3.8e12
        solid = bench.phage_aliquot(titre, 1e11)["volume_ul"]
        liquid = bench.phage_aliquot(
            titre, 1e11 * bench.LIQUID_PANNING_FACTOR)["volume_ul"]
        assert liquid == pytest.approx(solid * 2)


# --------------------------------------------------------------------------
# Batch sheet
# --------------------------------------------------------------------------

VHH = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
       "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDYWGQGTQVTVSS")

PROPERTIES = {"length": len(VHH), "mw_da": 13086.2, "pI": 8.60,
              "ext_coeff_cystines": 35535, "ext_coeff_reduced": 35410,
              "A280_0.1pct": 2.706, "n_Cys": 2, "n_primary_amines": 5,
              "sequence": VHH}


class TestBatchSheet:

    def test_the_blocks_stay_separate(self):
        """The whole point of the sheet.

        A predicted extinction coefficient quoted a year later as a
        measurement is the failure this layout exists to prevent, so nothing
        measured may appear among the sequence-derived values or the other
        way round.
        """
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"a280": 0.512})
        theoretical = {label for label, _value in sheet["theoretical"]}
        measured = {label for label, _value in sheet["measured"]}
        derived = {label for label, _value in sheet["derived"]}
        assert not theoretical & measured
        assert not theoretical & derived
        assert not measured & derived
        assert "A280" in measured
        assert "pI" in theoretical

    def test_a_concentration_from_an_absorbance_is_not_measured(self):
        """It is a reading divided by a predicted extinction coefficient and
        converted with a predicted molar mass, so it belongs in the derived
        block with its provenance, not among the readings."""
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"a280": 0.512})
        measured = {label for label, _value in sheet["measured"]}
        derived = dict(sheet["derived"])
        assert "Concentration" not in measured
        assert "Concentration" in derived
        expected = bench.concentration_from_a280(
            0.512, PROPERTIES["ext_coeff_cystines"], PROPERTIES["mw_da"])
        assert derived["Concentration"] == f"{expected['mg_per_ml']:.3f} mg/mL"
        assert any("predict" in note for note in sheet["derivation_notes"]), \
            "the sheet must say the concentration rests on a prediction"

    def test_a_concentration_determined_another_way_is_measured(self):
        """A BCA does not go through the sequence, so it is a reading, and
        the method is named because a concentration without one cannot be
        interpreted later."""
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"mg_per_ml": 1.84,
                                            "method": "BCA"})
        labels = [label for label, _value in sheet["measured"]]
        assert any("BCA" in label for label in labels)
        assert not any("Concentration" == label
                       for label, _value in sheet["derived"])

    def test_the_dilution_and_path_length_reach_the_concentration(self):
        plain = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"a280": 0.5})
        diluted = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                    measured={"a280": 0.5, "dilution": 10.0})
        # rel rather than exact: the helper rounds, so ten times a rounded
        # number is not the rounding of ten times it.
        assert (diluted["concentration_mg_per_ml"]
                == pytest.approx(plain["concentration_mg_per_ml"] * 10,
                                 rel=1e-3))

    def test_the_extinction_coefficient_choice_changes_the_answer(self):
        formed = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                   measured={"a280": 0.5}, cystines=True)
        reduced = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                    measured={"a280": 0.5}, cystines=False)
        assert (formed["concentration_mg_per_ml"]
                != reduced["concentration_mg_per_ml"])
        assert any("cystines formed" in note
                   for note in formed["derivation_notes"])

    def test_an_empty_measurement_leaves_the_block_empty(self):
        """Zeros are what an untouched number input gives, and a sheet full
        of zeros reads as a measurement of zero."""
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES)
        assert sheet["measured"] == []
        assert sheet["derived"] == []
        assert sheet["mass_difference_da"] is None

    def test_total_protein_is_concentration_times_volume(self):
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"mg_per_ml": 1.84,
                                            "method": "BCA",
                                            "volume_ul": 450.0})
        assert dict(sheet["derived"])["Total protein"] == "0.828 mg"

    def test_the_mass_difference_is_measured_minus_theoretical(self):
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES,
                                  measured={"measured_da": 13245.0})
        assert sheet["mass_difference_da"] == pytest.approx(158.8, abs=0.1)

    def test_the_date_fills_itself_in(self):
        from datetime import date

        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES)
        fields = {key: value for _label, value, key in sheet["fields"]}
        assert fields["date"] == date.today().isoformat()

    def test_every_declared_field_reaches_the_sheet(self):
        """Adding a field to BATCH_FIELDS must be one edit, not two."""
        sheet = bench.batch_sheet("Nb", VHH, PROPERTIES)
        keys = {key for _label, _value, key in sheet["fields"]}
        assert keys == {key for key, _label, _default, _help
                        in bench.BATCH_FIELDS}

    def test_the_pdf_is_one_page_and_says_what_it_is_not(self, tmp_path):
        import pypdf

        sheet = bench.batch_sheet("NbTest", VHH, PROPERTIES,
                                  measured={"a280": 0.512, "mg_per_ml": 1.84},
                                  fields={"lot": "JV-2026-014"})
        path = bench.build_batch_pdf(tmp_path / "batch.pdf", sheet,
                                     version="1.0.0", author="J. Vaneyck")
        assert path.exists() and path.stat().st_size > 1000

        reader = pypdf.PdfReader(str(path))
        assert len(reader.pages) == 1, "a tube label does not need two pages"
        text = reader.pages[0].extract_text()
        assert "JV-2026-014" in text
        assert "not a certificate of analysis" in text, \
            "the sheet must not be mistakable for a regulated document"
        assert "theoretical" in text


# --------------------------------------------------------------------------
# Purification balance
# --------------------------------------------------------------------------
# The arithmetic is deliberately checked by hand here rather than against the
# function's own output: a test that recomputes the implementation is a test
# of nothing. 1.2 AU / 24075 /M/cm * 14000 Da = 0.69782 mg/mL, times 6 mL is
# 4.1869 mg.

EPS = 24075.0
MW = 14000.0


def test_the_balance_totals_match_a_hand_calculation():
    out = bench.purification_balance(
        [{"step": "IMAC eluate", "volume_ml": 6.0, "a280": 1.2},
         {"step": "SEC pool", "volume_ml": 9.0, "a280": 0.55}],
        ext_coeff=EPS, mw_da=MW)

    imac, sec = out["rows"]
    assert imac["concentration_mg_per_ml"] == pytest.approx(0.69782, abs=1e-4)
    assert imac["total_mg"] == pytest.approx(4.1869, abs=1e-3)
    assert sec["total_mg"] == pytest.approx(2.8782, abs=1e-3)

    # 2.8782 / 4.1869 = 68.7%
    assert sec["step_recovery_pct"] == pytest.approx(68.7, abs=0.1)
    assert sec["cumulative_recovery_pct"] == pytest.approx(68.7, abs=0.1)
    # The first counted step is the baseline, so it is 100% of itself and has
    # nothing before it to have survived.
    assert imac["cumulative_recovery_pct"] == 100.0
    assert imac["step_recovery_pct"] is None
    assert out["reference_step"] == "IMAC eluate"
    assert out["final_total_mg"] == pytest.approx(2.8782, abs=1e-3)


def test_a_crude_step_keeps_its_total_and_gets_no_recovery():
    """The whole reason this is not a column of divisions.

    An A280 on a lysate is every protein in it. Divided by the target's
    extinction coefficient it gives a quantity of total protein expressed in
    the target's units, and the ratio between that and a purified pool is not
    a recovery of the target. Reporting one would be the most flattering
    number on the page and the least true.
    """
    out = bench.purification_balance(
        [{"step": "Clarified lysate", "volume_ml": 40.0, "a280": 2.5,
          "crude": True},
         {"step": "IMAC eluate", "volume_ml": 6.0, "a280": 1.2},
         {"step": "SEC pool", "volume_ml": 9.0, "a280": 0.55}],
        ext_coeff=EPS, mw_da=MW)

    lysate, imac, sec = out["rows"]
    assert lysate["total_mg"] is not None, \
        "the crude total is worth recording and should still be computed"
    assert lysate["step_recovery_pct"] is None
    assert lysate["cumulative_recovery_pct"] is None
    assert "crude" in lysate["comment"]

    # The baseline is the first step that is not crude, so the recoveries are
    # from capture onwards and the lysate cannot inflate them.
    assert out["reference_step"] == "IMAC eluate"
    assert imac["cumulative_recovery_pct"] == 100.0
    assert sec["cumulative_recovery_pct"] == pytest.approx(68.7, abs=0.1)


def test_recovery_above_one_hundred_percent_is_called_out():
    """Arithmetically fine, physically not: it is a volume typed wrong."""
    out = bench.purification_balance(
        [{"step": "IMAC", "volume_ml": 5.0, "a280": 1.0},
         {"step": "SEC", "volume_ml": 20.0, "a280": 0.9}],
        ext_coeff=EPS, mw_da=MW)

    assert out["rows"][1]["step_recovery_pct"] > 100.0
    assert any("cannot recover more" in w for w in out["warnings"]), \
        f"no warning about a gain in protein, got {out['warnings']}"


def test_mixing_two_quantification_methods_is_flagged():
    """A BCA against an A280 compares two methods as much as two steps."""
    out = bench.purification_balance(
        [{"step": "IMAC", "volume_ml": 5.0, "mg_per_ml": 0.8, "method": "BCA"},
         {"step": "SEC", "volume_ml": 4.0, "a280": 0.5}],
        ext_coeff=EPS, mw_da=MW)

    assert {row["basis"] for row in out["rows"]} == {"measured", "derived"}
    assert any("not all quantified the same way" in w for w in out["warnings"])
    # Shown anyway, with the basis beside each row, rather than withheld.
    assert out["rows"][1]["step_recovery_pct"] is not None


def test_a_reading_takes_precedence_over_an_a280_on_the_same_row():
    """It does not go through the sequence, so it is the better number."""
    out = bench.purification_balance(
        [{"step": "Pool", "volume_ml": 2.0, "a280": 1.0, "mg_per_ml": 0.5,
          "method": "BCA"}],
        ext_coeff=EPS, mw_da=MW)

    row = out["rows"][0]
    assert row["basis"] == "measured"
    assert row["concentration_mg_per_ml"] == pytest.approx(0.5)
    assert "BCA" in row["basis_detail"]


def test_an_a280_without_a_sequence_is_refused_rather_than_guessed():
    """No extinction coefficient means no concentration, not a default one."""
    out = bench.purification_balance(
        [{"step": "IMAC", "volume_ml": 5.0, "a280": 1.0}])

    assert out["rows"][0]["total_mg"] is None
    assert out["rows"][0]["basis"] == "needs the sequence"
    assert any("no extinction coefficient" in w for w in out["warnings"])


def test_yield_per_litre_divides_by_the_culture_volume():
    out = bench.purification_balance(
        [{"step": "Final", "volume_ml": 1.0, "mg_per_ml": 2.5,
          "method": "BCA"}],
        culture_volume_ml=500.0)

    # 2.5 mg in half a litre of culture is 5 mg/L.
    assert out["yield_mg_per_litre"] == pytest.approx(5.0)


def test_a_step_with_nothing_in_it_is_skipped_not_counted_as_zero():
    out = bench.purification_balance(
        [{"step": "IMAC", "volume_ml": 6.0, "a280": 1.2},
         {"step": "", "volume_ml": 0.0, "a280": 0.0}],
        ext_coeff=EPS, mw_da=MW)

    assert len(out["rows"]) == 1, \
        "an empty row was counted; a zero total would read as a failed step"


@pytest.mark.parametrize("step", [
    {"step": "x", "volume_ml": -1.0},
    {"step": "x", "volume_ml": 1.0, "a280": -0.5},
    {"step": "x", "volume_ml": 1.0, "mg_per_ml": -2.0},
])
def test_the_balance_refuses_negative_inputs(step):
    with pytest.raises(ValueError):
        bench.purification_balance([step], ext_coeff=EPS, mw_da=MW)


def test_an_empty_balance_is_empty_and_not_an_error():
    out = bench.purification_balance([])
    assert out["rows"] == [] and out["final_total_mg"] is None
    assert out["reference_step"] is None


# --------------------------------------------------------------------------
# Batch record, the machine-readable half of a batch sheet
# --------------------------------------------------------------------------

PROPS = {"mw_da": 14000.0, "ext_coeff_cystines": 24075,
         "ext_coeff_reduced": 23950, "pI": 7.2, "length": 117, "n_Cys": 2,
         "n_primary_amines": 5, "A280_0.1pct": 1.72}
SEQ = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
       "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGSTYYDYWGQGTQVTVSS")


def _record(**measured):
    readings = {"a280": 0.84, "dilution": 10.0, "path_cm": 1.0,
                "volume_ul": 600.0}
    readings.update(measured)
    sheet = bench.batch_sheet("NbX", SEQ, PROPS, measured=readings,
                              fields={"lot": "JV-2026-07",
                                      "prepared_by": "JV"})
    return bench.batch_record("NbX", SEQ, {"lot": "JV-2026-07",
                                           "prepared_by": "JV"},
                              readings, True, sheet=sheet,
                              version="1.0.0 (build 20261104.abc1234)")


def test_a_batch_record_round_trips_through_json():
    """The point of the file: the same readings come back out of it."""
    import json

    back = bench.batch_record_inputs(json.loads(json.dumps(_record())))

    assert back["fields"]["lot"] == "JV-2026-07"
    assert back["measured"]["a280"] == pytest.approx(0.84)
    assert back["measured"]["dilution"] == pytest.approx(10.0)
    assert back["measured"]["volume_ul"] == pytest.approx(600.0)
    assert back["cystines"] is True
    assert back["sequence"] == SEQ
    assert back["created_with"] == "1.0.0 (build 20261104.abc1234)"


def test_the_record_keeps_the_inputs_not_only_the_rendered_page():
    """A page of rendered strings cannot be recomputed or recombined.

    Keeping the readings is what makes the file worth having: the derived
    values are regenerated from them on reload, so a reprint says what this
    version of the app makes of those readings rather than repeating a number
    nobody can check.
    """
    record = _record()

    assert record["measured"]["a280"] == pytest.approx(0.84)
    assert record["schema"] == bench.BATCH_RECORD_SCHEMA
    # And the page as it was, so a reprint can be compared against it.
    assert record["rendered"]["derived"], "the rendered blocks were not kept"
    assert any("Concentration" in label
               for label, _value in record["rendered"]["derived"])


def test_the_record_carries_every_batch_field_even_when_blank():
    """A key that disappears when empty reads on reload as a field nobody
    filled in, which is not the same as one that was deliberately left out."""
    record = _record()
    expected = {key for key, _label, _default, _help in bench.BATCH_FIELDS}
    assert set(record["fields"]) == expected


def test_a_mass_left_unmeasured_stays_absent_through_the_round_trip():
    """The checkbox exists so an unmeasured mass is not reported as measured.
    A round trip must not resurrect it as a zero."""
    import json

    back = bench.batch_record_inputs(json.loads(json.dumps(_record())))
    assert back["measured"]["measured_da"] is None


@pytest.mark.parametrize("data, fragment", [
    ("not a dict", "JSON object"),
    ({}, "no `schema` key"),
    ({"schema": "something-else/1"}, "this app reads"),
    ({"schema": "vhh-toolbox-batch-record/9"}, "format"),
    ({"schema": "vhh-toolbox-batch-record/1", "fields": []}, "not an object"),
    ({"schema": "vhh-toolbox-batch-record/1",
      "measured": {"a280": "abc"}}, "not a number"),
])
def test_an_unreadable_record_is_refused_with_a_reason(data, fragment):
    """Refused whole rather than half understood.

    A record from which the unrecognised keys were quietly dropped would
    reprint a sheet missing fields, with nothing on the page saying so.
    """
    with pytest.raises(ValueError, match=fragment):
        bench.batch_record_inputs(data)
