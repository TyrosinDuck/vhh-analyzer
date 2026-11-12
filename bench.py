"""
bench - the arithmetic between a purified protein and the next experiment.

Everything else in this app reasons about a sequence. This module is about the
tube: how much protein is in it, whether it is the protein you think it is, and
how much labelling reagent to add. It is the part that gets done on the back of
a glove box on a Friday afternoon, which is precisely why it is worth writing
down once.

Four calculations, all of them deriving from quantities the rest of the app
already computes from the sequence (molar mass, extinction coefficient, lysine
count), so nothing here needs a parameter that was not already justified
somewhere else.

**No empirical model lives here.** Every number is either a measurement you
supply, a constant with a cited source, or arithmetic on those two. There is no
prediction of a labelling efficiency, of a yield, or of how a conjugate will
behave: those depend on your protein and your hands.
"""

from __future__ import annotations

import pandas as pd

__all__ = [
    "REAGENTS", "SOLVENTS", "concentration_from_a280", "dilution_plan",
    "stock_preparation", "conjugation_plan", "mass_check", "MODIFICATIONS",
    "PIPETTE_MINIMUM_UL", "desalting_column", "PURIFICATION_STEPS",
    "purification_balance", "CULTURE_ROUTES", "culture_for_target",
    "BATCH_RECORD_SCHEMA", "batch_record", "batch_record_inputs",
    "compare_batch_records", "__version__",
]

__version__ = "2026.11.12d"


# --------------------------------------------------------------------------
# Reagents
# --------------------------------------------------------------------------
# Only reagents whose molar mass was checked against a primary source are
# listed, and each carries where the number came from. A conjugation is molar
# arithmetic: a wrong molar mass does not fail loudly, it quietly gives you the
# wrong molar excess and a result you will blame on the chemistry.
#
# Anything not in this table is entered by hand, which is the honest option: a
# catalogue number is not a molar mass, vendors renumber products, and a dye's
# extinction coefficient and 280 nm correction factor are printed on the vial's
# certificate because they vary between lots.

REAGENTS: dict[str, dict] = {
    "NHS-PEG4-Biotin (Thermo 21329 / 21330 / A39259)": {
        "mw": 588.7,
        "formula": "C25H40N4O10S",
        "cas": "459426-22-3",
        "kind": "biotin",
        "solvents": ["water", "DMSO", "DMF"],
        "vial_mg": 2.0,
        "source": "Cayman Chemical product 29114 (Biotin-PEG4-NHS ester); "
                  "consistent with the Thermo No-Weigh protocol, where a 2 mg "
                  "vial in 170 uL gives 20 mM.",
        "note": "Not sulfonated, despite how it is often spoken about. It is "
                "water soluble because of the PEG4 spacer, not a sulfo group, "
                "so it needs no organic co-solvent.",
    },
    "Alexa Fluor 488 NHS ester (Thermo A20000 / A20100)": {
        "mw": 643.0,
        "formula": None,
        "cas": None,
        "kind": "dye",
        "solvents": ["DMSO", "DMF"],
        "vial_mg": 1.0,
        "epsilon": 71000.0,
        "lambda_max": 495,
        "cf280": 0.11,
        "source": "Thermo Fisher 'Amine-Reactive Probes' protocol MP10168, "
                  "Table 4: MW 643, eps 71000 /M/cm at 495 nm, CF280 0.11. "
                  "The A20000 product page now lists MW 643.4 and eps 73000. "
                  "Both are Thermo's own figures; the protocol table is the "
                  "one used here, and the 2.8% difference in eps carries "
                  "straight into a degree of labelling.",
        "note": "The molar mass, extinction coefficient and correction factor "
                "are the dye's published values. Your vial's certificate of "
                "analysis is the authority for the lot you hold.",
    },
}

SOLVENTS = ["water", "DMSO", "DMF", "PBS"]

# Below this, a volume cannot be pipetted accurately enough for the molar
# excess to mean anything, and the stock should be diluted first. The Thermo
# cheat sheet makes the same point with a fixed 1:10 dilution.
PIPETTE_MINIMUM_UL = 2.0


def desalting_column(volume_ul: float) -> str:
    """Which spin desalting column suits a reaction volume.

    The ranges are the ones in the Thermo EZ-Link protocol, kept because they
    are the columns most labs have on the shelf rather than because there is
    anything fundamental about them.
    """
    if volume_ul < 30:
        return ("below the usable range of a 0.5 mL column; make the reaction "
                "up to at least 30 uL, or dialyse")
    if volume_ul <= 170:
        return "0.5 mL spin desalting column"
    if volume_ul <= 700:
        return "2 mL spin desalting column"
    return "5 mL spin desalting column, or dialysis"


# --------------------------------------------------------------------------
# Quantification
# --------------------------------------------------------------------------

def concentration_from_a280(a280: float, ext_coeff: float, mw_da: float,
                            path_length_cm: float = 1.0,
                            dilution_factor: float = 1.0,
                            a260: float | None = None,
                            volume_ul: float | None = None) -> dict:
    """Concentration from a measured A280, by Beer-Lambert.

    ``ext_coeff`` is the molar extinction coefficient at 280 nm in /M/cm, which
    `vhh_core.characterize` reports in two forms: one assuming every cysteine is
    reduced and one assuming they are all paired. For a VHH out of the
    periplasm the canonical disulfide is formed, so the cystine value is usually
    the right one, and the difference between them is small enough that it is
    the choice nobody makes explicitly and everybody gets slightly wrong.

    ``a260`` is optional and buys you the 260/280 ratio, which is the only
    cheap way to notice that what you are quantifying is partly nucleic acid.
    """
    if ext_coeff <= 0:
        raise ValueError("the extinction coefficient must be positive; a "
                         "protein with no Trp, Tyr or Cys cannot be "
                         "quantified at 280 nm")
    if path_length_cm <= 0 or dilution_factor <= 0:
        raise ValueError("path length and dilution factor must be positive")

    absorbance = float(a280) * float(dilution_factor)
    molar = absorbance / (ext_coeff * path_length_cm)      # mol/L
    mg_per_ml = molar * mw_da                              # g/L == mg/mL

    out = {
        "a280_corrected": round(absorbance, 4),
        "molar_uM": round(molar * 1e6, 3),
        "mg_per_ml": round(mg_per_ml, 4),
        "ratio_260_280": None,
        "nucleic_acid_warning": None,
        "total_ug": None,
        "total_nmol": None,
    }

    if volume_ul:
        out["total_ug"] = round(mg_per_ml * float(volume_ul), 2)
        out["total_nmol"] = round(molar * float(volume_ul) * 1e-3 * 1e6, 3)

    if a260 is not None and a280:
        ratio = float(a260) / float(a280)
        out["ratio_260_280"] = round(ratio, 3)
        if ratio > 0.8:
            out["nucleic_acid_warning"] = (
                f"A260/A280 is {ratio:.2f}. A pure protein sits near 0.57; "
                f"above about 0.8 there is enough nucleic acid for the A280 to "
                f"overstate the protein, and this concentration is an upper "
                f"bound rather than a measurement.")
    return out


def dilution_plan(stock_mg_per_ml: float, target_mg_per_ml: float,
                  target_volume_ul: float) -> dict:
    """How much stock and how much buffer, for a wanted concentration.

    Refuses to dilute upwards rather than returning a negative volume of
    buffer, which is the arithmetic a spreadsheet happily does at the bench.
    """
    if stock_mg_per_ml <= 0 or target_mg_per_ml <= 0 or target_volume_ul <= 0:
        raise ValueError("concentrations and volume must be positive")
    if target_mg_per_ml > stock_mg_per_ml:
        raise ValueError(
            f"the target ({target_mg_per_ml:g} mg/mL) is above the stock "
            f"({stock_mg_per_ml:g} mg/mL); dilution cannot concentrate")

    stock_ul = target_volume_ul * target_mg_per_ml / stock_mg_per_ml
    return {
        "stock_ul": round(stock_ul, 2),
        "buffer_ul": round(target_volume_ul - stock_ul, 2),
        "fold": round(stock_mg_per_ml / target_mg_per_ml, 2),
    }


# --------------------------------------------------------------------------
# Conjugation
# --------------------------------------------------------------------------

def stock_preparation(vial_mg: float, reagent_mw: float,
                      target_mM: float) -> dict:
    """Volume of solvent that turns a weighed vial into a stock of known molarity."""
    if min(vial_mg, reagent_mw, target_mM) <= 0:
        raise ValueError("vial mass, molar mass and target molarity must be positive")
    umol = vial_mg / reagent_mw * 1000.0            # mg / (g/mol) = mmol -> umol
    volume_ul = umol / target_mM * 1000.0           # umol / (umol/mL) -> mL -> uL
    return {"umol_in_vial": round(umol, 3),
            "solvent_ul": round(volume_ul, 1),
            "target_mM": target_mM}


def conjugation_plan(protein_volume_ul: float, protein_mg_per_ml: float,
                     protein_mw_da: float, molar_excess: float,
                     reagent_stock_mM: float,
                     n_reactive_amines: int | None = None) -> dict:
    """How much labelling reagent to add, and whether you can pipette it.

    The arithmetic is the one every vendor cheat sheet prints, written out so
    that the units are visible rather than folded into a spreadsheet column:

        nmol protein  = volume (mL) x concentration (mg/mL) / molar mass (kDa)
        nmol reagent  = nmol protein x molar excess
        volume (uL)   = nmol reagent / stock (mM)

    because a millimolar stock is a nanomole per microlitre, which is the one
    unit identity worth remembering at the bench.
    """
    for name, value in (("volume", protein_volume_ul),
                        ("concentration", protein_mg_per_ml),
                        ("molar mass", protein_mw_da),
                        ("molar excess", molar_excess),
                        ("stock", reagent_stock_mM)):
        if value <= 0:
            raise ValueError(f"{name} must be positive")

    protein_mg = protein_volume_ul * 1e-3 * protein_mg_per_ml
    protein_nmol = protein_mg / protein_mw_da * 1e6
    reagent_nmol = protein_nmol * molar_excess
    reagent_ul = reagent_nmol / reagent_stock_mM        # nmol / (nmol/uL)
    # Rounded here rather than only on the way out, so that the diluted volume
    # below is the shown figure times the shown factor. A reader checking the
    # arithmetic on screen should not find it off by a rounding they cannot see;
    # the discrepancy is far below pipetting precision either way.
    reagent_ul = round(reagent_ul, 3)
    final_ul = protein_volume_ul + reagent_ul

    notes: list[str] = []
    dilution = None
    if reagent_ul < PIPETTE_MINIMUM_UL:
        factor = 10
        while reagent_ul * factor < PIPETTE_MINIMUM_UL:
            factor *= 10
        dilution = {
            "factor": factor,
            "diluted_stock_mM": reagent_stock_mM / factor,
            "volume_ul": round(reagent_ul * factor, 2),   # = shown value x factor
        }
        notes.append(
            f"{reagent_ul:.2f} uL cannot be pipetted accurately. Dilute the "
            f"stock 1:{factor} and add {reagent_ul * factor:.1f} uL of that "
            f"instead; the molar excess is unchanged, the error on it is not.")

    co_solvent_pct = reagent_ul / final_ul * 100.0
    if co_solvent_pct > 10:
        notes.append(
            f"The reagent brings {co_solvent_pct:.0f}% of the final volume. If "
            f"its solvent is DMSO or DMF, that is above the ~10% a protein "
            f"usually tolerates; use a more concentrated stock.")

    if n_reactive_amines is not None and molar_excess > n_reactive_amines:
        notes.append(
            f"The excess ({molar_excess:g}x) is higher than the number of "
            f"primary amines ({n_reactive_amines}). Nothing stops you, but past "
            f"this point the extra reagent can only react with amines that are "
            f"already labelled or hydrolyse.")

    return {
        "protein_nmol": round(protein_nmol, 3),
        "protein_ug": round(protein_mg * 1000, 2),
        "reagent_nmol": round(reagent_nmol, 3),
        "reagent_ul": reagent_ul,
        "final_volume_ul": round(final_ul, 2),
        "co_solvent_pct": round(co_solvent_pct, 1),
        "dilution": dilution,
        "desalting": desalting_column(final_ul),
        "notes": notes,
    }


# --------------------------------------------------------------------------
# Identity by mass
# --------------------------------------------------------------------------
# Monoisotopic masses of the elements, so every modification below is computed
# from its chemical formula rather than copied from a table. A tabulated delta
# is a number nobody can check; a formula is one anybody can.

_ATOM = {"H": 1.0078250319, "C": 12.0, "N": 14.0030740052,
         "O": 15.9949146221, "S": 31.97207069, "P": 30.97376151,
         "Na": 22.98976928, "K": 38.96370649}


def _formula_mass(formula: dict[str, int]) -> float:
    return sum(_ATOM[element] * count for element, count in formula.items())


# Each entry: what it is, the formula added (negative counts for a loss), and
# when it applies. The conditions matter as much as the masses: an explanation
# that cannot apply to your sequence is not an explanation.
MODIFICATIONS: list[dict] = [
    {"name": "N-terminal Met excised",
     "formula": {"C": -5, "H": -9, "N": -1, "O": -1, "S": -1},
     "requires": "starts_with_M",
     "why": "Methionine aminopeptidase removes it in E. coli when residue 2 is "
            "small (Ala, Cys, Gly, Pro, Ser, Thr, Val)."},
    {"name": "Disulfide formed (per bond)",
     "formula": {"H": -2},
     "requires": "has_two_cys",
     "why": "Two thiols become one bond and lose two hydrogens. A VHH out of "
            "the periplasm normally has its canonical pair formed."},
    {"name": "Pyroglutamate from N-terminal Gln",
     "formula": {"N": -1, "H": -3},
     "requires": "starts_with_Q",
     "why": "The N-terminal glutamine cyclises and loses ammonia. Spontaneous, "
            "and common enough to be the first thing to check on a -17."},
    {"name": "Pyroglutamate from N-terminal Glu",
     "formula": {"H": -2, "O": -1},
     "requires": "starts_with_E",
     "why": "The same cyclisation from glutamate, losing water instead."},
    {"name": "Gluconoylation of an N-terminal His-tag",
     "formula": {"C": 6, "H": 10, "O": 6},
     "requires": "his_tag",
     "why": "A documented E. coli artefact on N-terminal His-tags: the classic "
            "+178 that makes a clean preparation look wrong."},
    {"name": "Phosphogluconoylation of an N-terminal His-tag",
     "formula": {"C": 6, "H": 11, "O": 9, "P": 1},
     "requires": "his_tag",
     "why": "The phosphorylated form of the same artefact, +258."},
    {"name": "Oxidation (per Met or Trp)",
     "formula": {"O": 1},
     "requires": "has_met_or_trp",
     "why": "A single oxygen, usually on a methionine. Check the liabilities "
            "tab for which residues are exposed to it."},
    {"name": "Sodium adduct",
     "formula": {"Na": 1, "H": -1},
     "requires": None,
     "why": "A salt adduct that survives the ionisation, not a modification of "
            "the protein."},
    {"name": "Potassium adduct",
     "formula": {"K": 1, "H": -1},
     "requires": None,
     "why": "The same, from potassium."},
]


def _applies(requirement: str | None, seq: str) -> bool:
    if requirement is None:
        return True
    if requirement == "starts_with_M":
        return seq.startswith("M")
    if requirement == "starts_with_Q":
        return seq.startswith("Q")
    if requirement == "starts_with_E":
        return seq.startswith("E")
    if requirement == "has_two_cys":
        return seq.count("C") >= 2
    if requirement == "has_met_or_trp":
        return "M" in seq or "W" in seq
    if requirement == "his_tag":
        return "HHHHHH" in seq[:40] or "HHHHHH" in seq[-40:]
    return False


def mass_check(seq: str, measured_da: float | None = None,
               expected_da: float | None = None,
               tolerance_da: float = 2.0) -> dict:
    """Compare a measured mass with the sequence, and explain the difference.

    ``expected_da`` defaults to the average mass of the sequence as written. A
    deconvoluted mass from an intact LC-MS run is an average mass, so that is
    what is compared; the monoisotopic value is reported alongside because it is
    what you need if you are looking at an isotope envelope instead.

    Every candidate explanation is computed from a chemical formula and offered
    only when the sequence allows it. Nothing is ranked or chosen: more than one
    combination can land on the same number, and picking one for you would be
    presenting arithmetic as an identification.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    clean = "".join(c for c in seq.upper() if c.isalpha())
    analysis = ProteinAnalysis(clean)
    average = analysis.molecular_weight()
    monoisotopic = ProteinAnalysis(clean, monoisotopic=True).molecular_weight()

    expected = float(expected_da) if expected_da else average
    out = {
        "average_da": round(average, 2),
        "monoisotopic_da": round(monoisotopic, 2),
        "expected_da": round(expected, 2),
        "measured_da": None,
        "difference_da": None,
        "matches": None,
        "candidates": pd.DataFrame(
            columns=["modification", "delta_da", "n", "explains", "why"]),
    }
    if measured_da is None:
        return out

    difference = float(measured_da) - expected
    out["measured_da"] = round(float(measured_da), 2)
    out["difference_da"] = round(difference, 2)
    out["matches"] = abs(difference) <= tolerance_da
    if out["matches"]:
        return out

    rows = []
    for modification in MODIFICATIONS:
        if not _applies(modification["requires"], clean):
            continue
        delta = _formula_mass(modification["formula"])
        if abs(delta) < 1e-9:
            continue
        # How many of this modification would account for the difference, if a
        # whole number of them does.
        count = round(difference / delta)
        explains = (count >= 1
                    and abs(count * delta - difference) <= tolerance_da)
        rows.append({
            "modification": modification["name"],
            "delta_da": round(delta, 4),
            "n": count if explains else None,
            "explains": explains,
            "why": modification["why"],
        })

    frame = pd.DataFrame(rows, columns=["modification", "delta_da", "n",
                                        "explains", "why"])
    if len(frame):
        frame = frame.sort_values(["explains", "modification"],
                                  ascending=[False, True]).reset_index(drop=True)
    out["candidates"] = frame
    return out


# --------------------------------------------------------------------------
# Phage
# --------------------------------------------------------------------------
# Panning is the step before everything else on this screen: the titre of a
# phage stock decides how much of it goes into the next round, and it is read
# off a spectrophotometer like any other concentration.
#
# There are two ways to turn an A260 into a titre, and they do not agree.
# Both are offered, side by side, because choosing one silently would hide a
# difference of more than two orders of magnitude.

# The empirical route: a titration in which an absorbance of 1 corresponded to
# this many particles per millilitre. It is the number used in the lab
# protocol this calculator reproduces. It is a measured correspondence for one
# phage preparation and one instrument, not a physical constant.
PHAGE_EMPIRICAL_PER_A260 = 3.0e10

# The optical route for filamentous phage,
#   virions/mL = (A269 - A320) x 6e16 / (number of bases in the genome)
# which is the ssDNA absorbance of the particle worked out from its genome
# length. Reading at 260 rather than 269 nm, and against a blank rather than
# with the A320 scatter correction, costs some accuracy but not an order of
# magnitude.
#
# Attribution, honestly: this was carried here as "Day's relation", and that
# attribution has not been confirmed against any primary source. The formula
# in exactly this form is in wide circulation (MIT 20.109 course protocols,
# published filamentous-phage methods, vendor calculators), always without a
# citation of its own. So: the arithmetic is standard and reproducible, the
# constant's provenance is not established, and that is why the empirical
# route below is shown beside it rather than behind it.
PHAGE_OPTICAL_CONSTANT = 6.0e16

# M13mp18 and its phagemid derivatives. A phagemid carrying an insert is
# shorter than the helper phage, which is why this is an input rather than a
# constant.
PHAGE_GENOME_BASES = 7249

# A liquid-phase round is run with twice the input of a solid-phase one: the
# antigen is in solution and the capture step loses more of the library. This
# is a protocol convention, not something derived from anything here.
LIQUID_PANNING_FACTOR = 2.0


def phage_titre(a260: float, dilution: float = 1.0,
                glycerol_factor: float = 1.0,
                per_a260: float = PHAGE_EMPIRICAL_PER_A260,
                genome_bases: int = PHAGE_GENOME_BASES) -> dict:
    """Particles per millilitre in the undiluted stock, both ways.

    ``a260`` is what the instrument read, so ``dilution`` (how many times the
    stock was diluted to get there) and ``glycerol_factor`` (how much the stock
    itself was diluted when glycerol was added to it) multiply back up to the
    stock. A stock made by mixing one volume of phage with one volume of
    glycerol has a glycerol factor of 2.

    Returns both estimates and the ratio between them. Nothing is chosen: the
    two numbers describe different things, one a correspondence measured on a
    particular preparation and one a calculation from the genome, and when they
    disagree the disagreement is the result.
    """
    if a260 < 0:
        raise ValueError("an absorbance cannot be negative")
    if dilution <= 0 or glycerol_factor <= 0:
        raise ValueError("a dilution factor must be greater than zero")
    if per_a260 <= 0:
        raise ValueError("the empirical constant must be greater than zero")
    if genome_bases <= 0:
        raise ValueError("a genome has a positive number of bases")

    corrected = float(a260) * float(dilution) * float(glycerol_factor)
    empirical = corrected * float(per_a260)
    optical = corrected * PHAGE_OPTICAL_CONSTANT / float(genome_bases)
    return {
        "a260_corrected": corrected,
        "empirical_per_ml": empirical,
        "optical_per_ml": optical,
        "ratio": (optical / empirical) if empirical else None,
    }


def phage_aliquot(titre_per_ml: float, target_particles: float) -> dict:
    """How much stock to take to get a given number of particles.

    The volume comes straight out of the titre; what is worth saying is when
    it is too small to pipette, which is the same objection the conjugation
    planner makes, for the same reason.
    """
    if titre_per_ml <= 0:
        raise ValueError("the titre must be greater than zero")
    if target_particles <= 0:
        raise ValueError("the target must be greater than zero")

    volume_ul = target_particles / float(titre_per_ml) * 1000.0
    note = ""
    if volume_ul < PIPETTE_MINIMUM_UL:
        factor = max(2, int(PIPETTE_MINIMUM_UL / volume_ul) + 1)
        note = (f"{volume_ul:.2f} uL is below {PIPETTE_MINIMUM_UL:.0f} uL and "
                f"cannot be pipetted accurately; dilute the stock 1:{factor} "
                f"and take {volume_ul * factor:.1f} uL of that instead")
    elif volume_ul > 1000:
        note = ("more than a millilitre of stock; either the titre is low for "
                "this target or the target is high for this stock")
    return {"volume_ul": volume_ul, "note": note}


# --------------------------------------------------------------------------
# Panning enrichment
# --------------------------------------------------------------------------
# After a round, the eluted phage are titrated on two plates: the output of
# the panning on the antigen (+) and on the negative control without it (-).
# Each plate is counted at the dilution that gave countable colonies, so the
# count times the dilution factor is the output per volume plated, and the
# ratio of the two outputs is the enrichment. Nothing here is a constant: it
# is the laboratory's own counts and dilutions, multiplied and divided.
#
# The only thing added is the counting uncertainty. A colony count is a
# Poisson count, so its relative standard error is 1/sqrt(N), and a ratio of
# two independent counts has a relative error of sqrt(1/N+ + 1/N-). That is
# the floor set by counting alone; pipetting and the dilution series add to
# it, and nothing here knows by how much.
#
# The dilution is the part that is easy to get wrong, because a plate is
# labelled with the exponent ("-3") while the arithmetic needs the factor
# (1000). Typing 3 where 1000 was meant does not fail: it gives an answer
# that is wrong by whatever the two exponents differ by, and 30 colonies at
# "3" against 30 at "2" reads as an enrichment of 1.5 instead of 10. The
# screen therefore asks which of the two conventions is being used instead
# of naming one in a help text, and the conversion lives here.

# A dilution written as a power of ten this small is almost certainly an
# exponent typed where a factor was asked for: a real 1:2 or 1:5 dilution is
# possible, so this is a warning and not a refusal.
DILUTION_LOOKS_LIKE_EXPONENT = 10.0


def dilution_factor(value: float, as_power_of_ten: bool) -> float:
    """The dilution factor, from either convention.

    ``as_power_of_ten``: 3 means a 10^-3 dilution, so a factor of 1000.
    Otherwise the value is the factor itself, which is what a serial
    dilution that is not a power of ten (1:2, 1:5) needs.
    """
    if as_power_of_ten:
        if value < 0:
            raise ValueError("a dilution is 10^-n with n of zero or more")
        return 10.0 ** float(value)
    if value < 1:
        raise ValueError("a dilution factor is 1 or more (1 for neat)")
    return float(value)


def panning_enrichment(plus_colonies: float, plus_dilution: float,
                       minus_colonies: float, minus_dilution: float,
                       plated_ul: float | None = None) -> dict:
    """Output of the + and - plates of one round, and their ratio.

    ``*_dilution`` is the dilution factor of the plate counted, 1000 for a
    10^-3 dilution. ``plated_ul`` is the volume spread on each plate, the same
    on both; without it the outputs are colonies times dilution, which is
    enough for the ratio, and with it they are colony-forming units per mL.

    The enrichment is None when the - plate has no colony: a ratio to zero is
    not a large number, it is a plate that was not diluted little enough to
    say anything, and the note says so.
    """
    for name, value in [("colony count", plus_colonies),
                        ("colony count", minus_colonies)]:
        if value is None or value < 0 or float(value) != int(value):
            raise ValueError(f"a {name} is a whole number of zero or more")
    for value in (plus_dilution, minus_dilution):
        if value is None or value < 1:
            raise ValueError("a dilution factor is 1 or more (1000 for a "
                             "10^-3 dilution)")
    if plated_ul is not None and plated_ul <= 0:
        raise ValueError("the plated volume must be above zero")

    per_ml = 1000.0 / plated_ul if plated_ul else 1.0
    plus = plus_colonies * plus_dilution * per_ml
    minus = minus_colonies * minus_dilution * per_ml
    enrichment, uncertainty, note = None, None, ""
    if minus_colonies == 0:
        note = ("no colony on the - plate: the enrichment cannot be "
                "computed; count a less diluted - plate")
    else:
        enrichment = plus / minus
        if plus_colonies > 0:
            uncertainty = (1 / plus_colonies + 1 / minus_colonies) ** 0.5
        else:
            note = "no colony on the + plate"
    return {"plus_output": plus, "minus_output": minus,
            "per_ml": plated_ul is not None, "enrichment": enrichment,
            "relative_uncertainty": uncertainty, "note": note}


# --------------------------------------------------------------------------
# Purification balance
# --------------------------------------------------------------------------
# The other views in this section are about one tube. This one follows the
# material through a preparation: how much protein was in the pool after each
# step, how much of the previous step survived it, and how much came out per
# litre of culture. It is the page a platform hands to whoever asked for the
# protein, and it is arithmetic on volumes and absorbances, with no model of
# what a step ought to recover.
#
# The reason this is not simply a column of divisions is the crude steps. An
# A280 on a lysate or a clarified supernatant is the absorbance of every
# protein in it; dividing it by the extinction coefficient of the target gives
# a quantity of total protein expressed in the target's units, and a ratio
# between that and a purified pool is not a recovery of the target. Anything
# before the first capture step is therefore marked crude, kept in the table
# because it is worth recording, and left out of the recovery arithmetic.

# The steps of a typical VHH preparation, used to prefill the editor. Only
# names: no volume, no absorbance and no expectation of what any step
# recovers.
PURIFICATION_STEPS = [
    "Clarified lysate",
    "IMAC eluate",
    "Size exclusion pool",
    "Final concentrate",
]


def purification_balance(steps: list[dict], ext_coeff: float | None = None,
                         mw_da: float | None = None,
                         culture_volume_ml: float | None = None) -> dict:
    """Total protein, step recovery and cumulative recovery, step by step.

    Each entry of ``steps`` is a dict with a ``step`` name, a ``volume_ml``,
    and a concentration given one of two ways: ``a280`` (with optional
    ``dilution`` and ``path_cm``), which needs ``ext_coeff`` and ``mw_da`` and
    is *derived*, or ``mg_per_ml`` with ``method``, which is a reading by
    something that does not go through the sequence and is *measured*. A step
    with no volume or no concentration is skipped rather than counted as zero.

    ``crude`` marks a step where the target is not the only thing absorbing:
    a lysate, a clarified supernatant, a flow-through. Those rows carry a
    total, because the number is worth writing down, but no recovery, because
    the ratio would not be the target's.

    Recovery is computed only between two rows that are both non-crude and
    quantified the same way. A step recovery is against the row before it in
    the table, and the cumulative recovery is against the first non-crude row,
    which is named in ``reference_step``.

    Returns plain data: ``rows``, ``warnings``, ``final_total_mg``,
    ``yield_mg_per_litre`` and ``reference_step``. Nothing here predicts a
    yield; every number is a volume or an absorbance you supplied, or
    arithmetic on those.
    """
    rows: list[dict] = []
    warnings: list[str] = []

    for raw in steps or []:
        name = str(raw.get("step") or "").strip()
        volume_ml = raw.get("volume_ml")
        a280 = raw.get("a280")
        direct = raw.get("mg_per_ml")

        volume_ml = float(volume_ml or 0.0)
        if volume_ml < 0:
            raise ValueError(f"{name or 'a step'} has a negative volume")
        if a280 is not None and float(a280 or 0.0) < 0:
            raise ValueError(f"{name or 'a step'} has a negative A280")
        if direct is not None and float(direct or 0.0) < 0:
            raise ValueError(f"{name or 'a step'} has a negative concentration")

        concentration = None
        basis = None
        detail = ""
        if direct:
            concentration = float(direct)
            basis = "measured"
            detail = str(raw.get("method") or "method not stated")
        elif a280:
            if not (ext_coeff and mw_da):
                # Said once, below, rather than per row: a table of eight
                # identical warnings is a table nobody reads.
                basis = "needs the sequence"
            else:
                out = concentration_from_a280(
                    float(a280), float(ext_coeff), float(mw_da),
                    path_length_cm=float(raw.get("path_cm") or 1.0),
                    dilution_factor=float(raw.get("dilution") or 1.0))
                concentration = out["mg_per_ml"]
                basis = "derived"
                detail = f"A280, eps {ext_coeff:g} /M/cm"

        if not name and concentration is None and not volume_ml:
            continue

        total_mg = (round(concentration * volume_ml, 4)
                    if (concentration and volume_ml) else None)
        rows.append({
            "step": name or "unnamed step",
            "volume_ml": volume_ml or None,
            "concentration_mg_per_ml": (round(concentration, 4)
                                        if concentration else None),
            "total_mg": total_mg,
            "basis": basis,
            "basis_detail": detail,
            "crude": bool(raw.get("crude")),
            "step_recovery_pct": None,
            "cumulative_recovery_pct": None,
            "comment": "",
        })

    if any(row["basis"] == "needs the sequence" for row in rows):
        warnings.append(
            "Some steps give an A280 but there is no extinction coefficient to "
            "divide it by. Pick a sequence at the top, or enter those "
            "concentrations as a reading with the method named.")

    countable = [row for row in rows if row["total_mg"] and not row["crude"]]
    reference = countable[0] if countable else None

    if reference:
        bases = {row["basis"] for row in countable}
        if len(bases) > 1:
            warnings.append(
                "The steps counted here were not all quantified the same way, "
                "so a ratio between them compares two methods as much as two "
                "steps. The recoveries below are shown anyway, and the basis "
                "column says which is which, but treat them as approximate.")

        previous = None
        for row in rows:
            if not row["total_mg"] or row["crude"]:
                if row["crude"]:
                    row["comment"] = ("crude, left out of the recovery: the "
                                      "A280 is not only the target")
                continue
            if previous is not None:
                row["step_recovery_pct"] = round(
                    100.0 * row["total_mg"] / previous["total_mg"], 1)
            row["cumulative_recovery_pct"] = round(
                100.0 * row["total_mg"] / reference["total_mg"], 1)
            previous = row

        # A recovery above 100% is arithmetically fine and physically not: the
        # later pool cannot hold more of the target than the one it came from.
        # It is almost always a volume or a dilution factor, so say that
        # rather than printing a number that looks like a good yield.
        gained = [row["step"] for row in rows
                  if row["step_recovery_pct"] and row["step_recovery_pct"] > 100.0]
        if gained:
            warnings.append(
                "More protein came out of " + ", ".join(gained) + " than went "
                "in. A step cannot recover more than it was given, so one of "
                "the volumes or dilution factors above is wrong, or the two "
                "pools were not quantified the same way.")

    final = countable[-1]["total_mg"] if countable else None
    per_litre = None
    if final and culture_volume_ml:
        culture_volume_ml = float(culture_volume_ml)
        if culture_volume_ml <= 0:
            raise ValueError("the culture volume must be positive")
        per_litre = round(final / (culture_volume_ml / 1000.0), 3)

    return {
        "rows": rows,
        "warnings": warnings,
        "final_total_mg": final,
        "yield_mg_per_litre": per_litre,
        "reference_step": reference["step"] if reference else None,
    }


# --------------------------------------------------------------------------
# How much culture, for how much protein
# --------------------------------------------------------------------------
# The purification balance runs forwards: these volumes and absorbances came
# out of a preparation, here is what it recovered. The question asked before a
# preparation runs the other way, and it is the one a logistics role is asked
# most often: somebody needs 5 mg, how much culture do we set up.
#
# There is no model here and there will not be one. A yield cannot be
# predicted from a sequence, and anything claiming to would be the most
# dangerous number in this app: plausible, precise, and wrong in a way nobody
# notices until a purification has already been run. What this does is divide
# a target by a yield *you measured on a previous preparation of the same
# protein in the same host*, which is arithmetic on your own data and nothing
# else. With no previous preparation to divide by, the honest answer is that
# the app cannot help, and it says so rather than offering a default.

# The two ways a measured yield can be stated, which must not be mixed. They
# are listed here because the screen draws them from this and because the
# distinction is the whole point of the function below.
CULTURE_ROUTES = {
    "final": "mg of purified protein per litre of culture, from a previous "
             "preparation taken all the way through",
    "expression": "mg per litre expressed, plus the fraction of it a "
                  "purification recovered",
}


def culture_for_target(target_mg: float | None = None,
                       target_mg_per_ml: float | None = None,
                       target_volume_ml: float | None = None,
                       route: str = "final",
                       yield_mg_per_litre: float | None = None,
                       recovery_pct: float | None = None,
                       margin_pct: float = 0.0,
                       flask_working_ml: float | None = None) -> dict:
    """How much culture a target quantity needs, given a yield you measured.

    The target is given either as a mass in mg, or as a concentration and a
    volume, in which case the mass is their product. Giving both is an error
    rather than a precedence rule: the two disagreeing is a sign that one of
    them is left over from a different plan, and silently preferring one would
    hide that.

    ``route`` picks which measured yield is being divided by, and the two must
    not be combined:

    * ``"final"``: ``yield_mg_per_litre`` is mg of *purified* protein per
      litre, from a preparation taken all the way through. The recovery is
      already inside that figure, so supplying ``recovery_pct`` as well would
      count it twice and is refused.
    * ``"expression"``: ``yield_mg_per_litre`` is how much was *expressed* per
      litre, and ``recovery_pct`` is the fraction a purification recovered.
      Both are needed, because expression alone does not say what reaches the
      tube.

    ``margin_pct`` adds to the target before dividing, for the material a
    preparation loses to assay aliquots and to the volume that stays in the
    column. It is yours to set: nothing here knows what your steps lose.

    ``flask_working_ml`` is the working volume of one flask, which turns a
    culture volume into a number of flasks. It is optional because it is a
    property of your shaker, not of the protein.

    Returns ``required_mg``, ``culture_ml``, ``culture_litres``, ``flasks``,
    ``basis`` naming what was divided by what, and ``warnings``. Every figure
    is arithmetic on what you supplied: there is no yield in this module and
    no default for one.
    """
    if route not in CULTURE_ROUTES:
        raise ValueError(f"the route must be one of {sorted(CULTURE_ROUTES)}, "
                         f"not {route!r}")

    given_mass = target_mg not in (None, 0)
    given_product = (target_mg_per_ml not in (None, 0)
                     and target_volume_ml not in (None, 0))
    if given_mass and given_product:
        raise ValueError(
            "the target is given twice, once as a mass and once as a "
            "concentration times a volume. Give one of the two: if they "
            "disagree, one of them belongs to a different plan")
    if not given_mass and not given_product:
        raise ValueError(
            "no target: give a mass in mg, or a concentration and the volume "
            "you need it in")

    if given_mass:
        target = float(target_mg)
        target_basis = f"{target:g} mg asked for"
    else:
        concentration = float(target_mg_per_ml)
        volume_ml = float(target_volume_ml)
        if concentration <= 0 or volume_ml <= 0:
            raise ValueError("a target concentration and volume must both be "
                             "positive")
        target = concentration * volume_ml
        target_basis = (f"{concentration:g} mg/mL in {volume_ml:g} mL, which "
                        f"is {target:g} mg")
    if target <= 0:
        raise ValueError("the target mass must be positive")

    margin_pct = float(margin_pct or 0.0)
    if margin_pct < 0:
        raise ValueError("the margin cannot be negative")
    required = target * (1.0 + margin_pct / 100.0)

    if yield_mg_per_litre in (None, 0):
        raise ValueError(
            "no measured yield to divide by. This needs mg per litre from a "
            "previous preparation of the same protein in the same host: "
            "nothing here can predict one")
    measured_yield = float(yield_mg_per_litre)
    if measured_yield <= 0:
        raise ValueError("the measured yield must be positive")

    warnings: list[str] = []
    if route == "final":
        if recovery_pct not in (None, 0):
            raise ValueError(
                "a recovery was given alongside a final yield. A final yield "
                "is what came out after purification, so the recovery is "
                "already in it and applying it again would count the losses "
                "twice. Either drop the recovery, or switch to the "
                "expression route and give what was expressed")
        per_litre = measured_yield
        basis = (f"{target_basis}, divided by {measured_yield:g} mg of "
                 f"purified protein per litre measured on a previous "
                 f"preparation")
    else:
        if recovery_pct in (None, 0):
            raise ValueError(
                "the expression route needs the fraction a purification "
                "recovered, as a percentage. Without it, mg expressed per "
                "litre says nothing about what reaches the tube")
        recovery = float(recovery_pct)
        if not 0 < recovery <= 100:
            raise ValueError("the recovery is a percentage above 0 and at "
                             "most 100")
        per_litre = measured_yield * recovery / 100.0
        basis = (f"{target_basis}, divided by {measured_yield:g} mg expressed "
                 f"per litre times a {recovery:g} per cent recovery, which is "
                 f"{per_litre:g} mg per litre reaching the tube")

    litres = required / per_litre
    flasks = None
    if flask_working_ml:
        working = float(flask_working_ml)
        if working <= 0:
            raise ValueError("a flask working volume must be positive")
        # Rounded up: half a flask is a flask, and rounding down would plan a
        # preparation that cannot reach the target it was asked for.
        import math
        flasks = int(math.ceil(litres * 1000.0 / working))

    if margin_pct:
        warnings.append(
            f"The target was raised by {margin_pct:g} per cent before "
            f"dividing, to {required:.3f} mg, which is the margin you asked "
            f"for and not an estimate of anything.")
    warnings.append(
        "This is one division, on a yield you measured. A different host, a "
        "different induction or a different construct makes that yield "
        "inapplicable, and the app cannot tell that it has changed.")

    return {
        "target_mg": round(target, 4),
        "required_mg": round(required, 4),
        "effective_yield_mg_per_litre": round(per_litre, 4),
        "culture_litres": round(litres, 4),
        "culture_ml": round(litres * 1000.0, 1),
        "flasks": flasks,
        "route": route,
        "basis": basis,
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# Batch sheet
# --------------------------------------------------------------------------
# What goes on the label of a tube that leaves your bench. Not a certificate
# of analysis in the regulated sense: there are no specifications here, no
# acceptance criteria and no release signature, and calling it one would
# promise a document this is not. It is the page that stops a freezer box
# turning into a row of unlabelled aliquots in six months.

# The fields the preparer fills in, in the order they go on the page. Kept
# here rather than in the screen so that a notebook can build the same sheet
# without Streamlit, and so that adding a field is one edit.
BATCH_FIELDS = [
    ("lot", "Lot", "", "Whatever your freezer box calls it."),
    ("prepared_by", "Prepared by", "", ""),
    ("date", "Date", "", "Defaults to today."),
    ("host", "Expression host", "", "E. coli WK6, HEK293, and so on."),
    ("purification", "Purification", "",
     "IMAC, size exclusion, the order they were run in."),
    ("buffer", "Buffer", "PBS, pH 7.4", ""),
    ("storage", "Storage", "-80 C", ""),
    ("notes", "Notes", "", "Anything the next person needs to know."),
]


def batch_sheet(name: str, seq: str, properties: dict,
                measured: dict | None = None,
                fields: dict | None = None,
                cystines: bool = True) -> dict:
    """Assemble what goes on a batch sheet, as plain data.

    Four blocks, kept apart on purpose and labelled as such on the page:

    * **Written down**: the free-text fields, which nothing here can check.
    * **Measured**: only what an instrument actually read. An absorbance, a
      volume, a deconvoluted mass. Nothing in this block depends on the
      sequence.
    * **Derived**: what comes from combining the two. A concentration in
      mg/mL is the one everybody calls a measurement and none of them is:
      it is an absorbance divided by a *predicted* extinction coefficient,
      and converted to a molarity with a *predicted* molar mass. Putting it
      in with the readings is how a prediction gets quoted a year later as
      something that was measured, so it sits in its own block, with the
      extinction coefficient it used written next to it.
    * **From the sequence**: molar mass, pI, extinction coefficients, counts.
      Theoretical, every one of them.

    ``measured`` takes ``a280``, ``dilution``, ``path_cm``, ``volume_ul``,
    ``measured_da``, and, when the concentration did not come from an
    absorbance at all, ``mg_per_ml`` with ``method`` naming how it was
    obtained. A concentration given that way is a measurement and is listed
    as one.

    ``cystines`` picks which of the two extinction coefficients is used.
    """
    from datetime import date as _date

    fields = dict(fields or {})
    fields.setdefault("date", _date.today().isoformat())
    measured = dict(measured or {})

    clean = "".join(c for c in (seq or "").upper() if c.isalpha())
    eps_key = "ext_coeff_cystines" if cystines else "ext_coeff_reduced"
    eps = properties.get(eps_key)
    eps_label = "cystines formed" if cystines else "all cysteines reduced"
    mw_da = properties.get("mw_da")

    theoretical = [
        ("Length", f"{properties.get('length', len(clean))} residues"),
        ("Molar mass", f"{mw_da:.1f} Da" if mw_da else ""),
        ("pI", f"{properties['pI']:.2f}" if properties.get("pI") else ""),
        ("eps 280, cystines formed",
         f"{properties['ext_coeff_cystines']} /M/cm"
         if properties.get("ext_coeff_cystines") else ""),
        ("eps 280, reduced", f"{properties['ext_coeff_reduced']} /M/cm"
         if properties.get("ext_coeff_reduced") else ""),
        ("A280 of 1 mg/mL", f"{properties['A280_0.1pct']:.3f}"
         if properties.get("A280_0.1pct") else ""),
        ("Cysteines", str(properties.get("n_Cys", clean.count("C")))),
        ("Primary amines", str(properties.get("n_primary_amines", ""))),
    ]
    theoretical = [(label, value) for label, value in theoretical if value != ""]

    read, derived, notes = [], [], []
    a280 = measured.get("a280")
    dilution = float(measured.get("dilution") or 1.0)
    path_cm = float(measured.get("path_cm") or 1.0)
    volume_ul = measured.get("volume_ul")
    direct = measured.get("mg_per_ml")

    if a280 is not None:
        row = f"{a280:.3f}"
        if dilution != 1.0:
            row += f" (diluted {dilution:g}x)"
        if path_cm != 1.0:
            row += f", {path_cm:g} cm path"
        read.append(("A280", row))
    if volume_ul:
        read.append(("Volume", f"{volume_ul:.0f} uL"))

    concentration = None
    if direct:
        # Determined by something that does not go through the sequence, so
        # it belongs with the readings, named by how it was obtained.
        method = measured.get("method") or "method not stated"
        read.append((f"Concentration ({method})", f"{direct:.3f} mg/mL"))
        concentration = float(direct)
    elif a280 is not None and eps and mw_da:
        out = concentration_from_a280(a280, eps, mw_da, path_length_cm=path_cm,
                                      dilution_factor=dilution)
        concentration = out["mg_per_ml"]
        derived.append(("Concentration", f"{concentration:.3f} mg/mL"))
        derived.append(("Concentration, molar", f"{out['molar_uM']:.1f} uM"))
        notes.append(
            f"The concentration is the A280 divided by the extinction "
            f"coefficient predicted from the sequence, "
            f"{eps} /M/cm with {eps_label}, and converted to a molarity "
            f"with the predicted molar mass of {mw_da:.1f} Da. Both are "
            f"predictions, so this is not a measured concentration.")

    if concentration and volume_ul:
        derived.append(("Total protein",
                        f"{concentration * volume_ul / 1000:.3f} mg"))

    difference = None
    if measured.get("measured_da"):
        read.append(("Measured mass", f"{measured['measured_da']:.1f} Da"))
        if mw_da:
            difference = float(measured["measured_da"]) - float(mw_da)
            derived.append(("Mass minus the sequence", f"{difference:+.1f} Da"))

    return {
        "name": name,
        "sequence": clean,
        "fields": [(label, fields.get(key, ""), key)
                   for key, label, _default, _help in BATCH_FIELDS],
        "theoretical": theoretical,
        "measured": read,
        "derived": derived,
        "derivation_notes": notes,
        "concentration_mg_per_ml": concentration,
        "mass_difference_da": difference,
    }


# The identifier a reader checks before trusting the keys below. The number is
# the format's own, not the app's: a record written today has to stay readable
# by a later version of the app, and the only way to notice that it is not is
# to refuse a schema the code does not know rather than to guess at its keys.
BATCH_RECORD_SCHEMA = "vhh-toolbox-batch-record/1"


def batch_record(name: str, seq: str, fields: dict, measured: dict,
                 cystines: bool, sheet: dict | None = None,
                 version: str = "") -> dict:
    """The batch sheet as data a program can read, for saving beside the PDF.

    The PDF is for a human and for a freezer box. This is the same batch in a
    form something else can ingest: a lab's own database, a spreadsheet, a
    script that collects a year of preparations. The app stores nothing, so the
    file goes to whoever generated it and comes back the same way.

    What is kept is the *inputs*: the free-text fields, the readings, and which
    extinction coefficient was chosen. Those are what cannot be recomputed.
    The rendered blocks are kept alongside them, under ``rendered``, so that a
    sheet reprinted later can be compared against what the original said
    rather than assumed to match it: the derived values are recomputed on
    reload, and ``created_with`` is there because a later version of the app
    could compute them differently.

    This is not a certificate, an audit trail or a signature. It carries no
    release decision, nothing stops it being edited in a text editor, and
    nothing here can tell whether it was.
    """
    from datetime import datetime, timezone

    record = {
        "schema": BATCH_RECORD_SCHEMA,
        "created_with": version,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "name": name or "",
        "sequence": "".join(c for c in (seq or "").upper() if c.isalpha()),
        "cystines_formed": bool(cystines),
        "fields": {key: str((fields or {}).get(key, "") or "")
                   for key, _label, _default, _help in BATCH_FIELDS},
        # Only the keys batch_sheet actually reads, so that a record cannot
        # carry a number the sheet silently ignores.
        "measured": {key: (measured or {}).get(key)
                     for key in ("a280", "dilution", "path_cm", "volume_ul",
                                 "measured_da", "mg_per_ml", "method")},
    }
    if sheet:
        record["rendered"] = {
            "measured": [list(pair) for pair in sheet.get("measured") or []],
            "derived": [list(pair) for pair in sheet.get("derived") or []],
            "theoretical": [list(pair) for pair in sheet.get("theoretical") or []],
            "derivation_notes": list(sheet.get("derivation_notes") or []),
        }
    return record


def batch_record_inputs(data: dict) -> dict:
    """Read a saved record back into the inputs that made it.

    Returns ``name``, ``sequence``, ``fields``, ``measured``, ``cystines``,
    ``created_with`` and ``rendered``. The derived values are deliberately not
    returned: they are recomputed from the readings by `batch_sheet`, so a
    reprint says what this version of the app makes of those readings rather
    than repeating a number whose derivation nobody can check.

    Refuses a schema it does not know instead of picking out the keys it
    recognises. A record half understood is worse than one rejected: the
    fields that did not survive would be missing from the reprinted sheet
    without anything saying so.
    """
    if not isinstance(data, dict):
        raise ValueError("a batch record is a JSON object, and this is not one")
    schema = str(data.get("schema") or "")
    if not schema:
        raise ValueError(
            "this file carries no `schema` key, so it is not a batch record "
            "written by this app")
    family, _, number = schema.rpartition("/")
    expected_family, _, expected_number = BATCH_RECORD_SCHEMA.rpartition("/")
    if family != expected_family:
        raise ValueError(
            f"this file says it is `{schema}`, and this app reads "
            f"`{BATCH_RECORD_SCHEMA}`")
    if number != expected_number:
        raise ValueError(
            f"this record is format {schema} and this app reads "
            f"{BATCH_RECORD_SCHEMA}. A newer app wrote it, or an older one: "
            f"either way the keys are not guaranteed to mean the same thing, "
            f"so it is not being read rather than read wrongly")

    # Typed before the `or {}`, not after. With the fallback first, an empty
    # list is falsy and became an empty dict, so `"fields": []` passed a check
    # written to catch exactly that and reprinted a sheet with every text
    # field blank and nothing saying why.
    fields = data.get("fields")
    if fields is not None and not isinstance(fields, dict):
        raise ValueError("`fields` in this record is not an object")
    fields = fields or {}
    measured = data.get("measured")
    if measured is not None and not isinstance(measured, dict):
        raise ValueError("`measured` in this record is not an object")
    measured = measured or {}

    numeric = {}
    for key in ("a280", "dilution", "path_cm", "volume_ul", "measured_da",
                "mg_per_ml"):
        value = measured.get(key)
        if value in (None, ""):
            numeric[key] = None
            continue
        try:
            numeric[key] = float(value)
        except (TypeError, ValueError):
            raise ValueError(
                f"`measured.{key}` in this record is {value!r}, which is not "
                f"a number") from None
    numeric["method"] = str(measured.get("method") or "")

    return {
        "name": str(data.get("name") or ""),
        "sequence": "".join(c for c in str(data.get("sequence") or "").upper()
                            if c.isalpha()),
        "fields": {key: str(fields.get(key, "") or "")
                   for key, _label, _default, _help in BATCH_FIELDS},
        "measured": numeric,
        "cystines": bool(data.get("cystines_formed", True)),
        "created_with": str(data.get("created_with") or ""),
        "rendered": data.get("rendered") or {},
    }


def compare_batch_records(records: list[tuple[str, dict]],
                          properties_for=None,
                          version: str = "") -> dict:
    """Several saved batch records side by side, one row per lot.

    A batch record is the one thing this app produces that accumulates. One of
    them describes a tube; a year of them describes a preparation, and the
    question they answer together, whether the last three lots came out like
    the first three, is not visible in any single sheet. The app stores
    nothing, so the records come back from wherever the user kept them, which
    is also why this takes them as a list rather than reading a directory.

    ``records`` is a list of ``(label, data)``, where ``label`` is whatever
    the record came from, usually a file name, and ``data`` is the parsed
    JSON. A record that cannot be read is reported in ``rejected`` with the
    reason and left out of the rows, rather than failing the whole comparison:
    one unreadable file among twelve should cost that file, not the other
    eleven.

    ``properties_for`` is called with a sequence and returns the properties
    dict `batch_sheet` needs. It is a parameter rather than an import because
    this module computes nothing from a sequence and is not about to start:
    the caller owns that, and a test can pass a stub. Without it, the rows
    carry what the records state and leave the derived concentration out.

    What is deliberately *recomputed* rather than read back: every derived
    value. A record carries the readings and, under ``rendered``, what the
    sheet said at the time. The comparison recomputes from the readings with
    the version running now, so a column that disagrees with an old printout
    is visible instead of being carried forward. ``created_with`` is in each
    row for exactly that reason.

    What this is not: a trend. There is no fit, no control limit and no
    verdict on whether a lot is acceptable. Those need a specification, and a
    specification is a decision somebody makes, not something arithmetic can
    supply.
    """
    rows: list[dict] = []
    rejected: list[dict] = []
    warnings: list[str] = []

    for label, data in records or []:
        try:
            inputs = batch_record_inputs(data)
        except ValueError as exc:
            rejected.append({"source": str(label), "reason": str(exc)})
            continue

        properties: dict = {}
        if properties_for and inputs["sequence"]:
            try:
                properties = properties_for(inputs["sequence"]) or {}
            except Exception as exc:          # a stub, or a sequence it hates
                rejected.append({
                    "source": str(label),
                    "reason": f"the sequence in it could not be "
                              f"characterised ({exc})"})
                continue

        sheet = batch_sheet(inputs["name"], inputs["sequence"], properties,
                            measured=inputs["measured"],
                            fields=inputs["fields"],
                            cystines=inputs["cystines"])

        measured = inputs["measured"]
        if measured.get("mg_per_ml"):
            basis = "measured"
            basis_detail = measured.get("method") or "method not stated"
        elif measured.get("a280") is not None and sheet[
                "concentration_mg_per_ml"]:
            basis = "derived"
            eps_key = ("ext_coeff_cystines" if inputs["cystines"]
                       else "ext_coeff_reduced")
            eps = properties.get(eps_key)
            basis_detail = (f"A280, eps {eps:g} /M/cm" if eps
                            else "A280")
        else:
            basis = None
            basis_detail = ""

        concentration = sheet["concentration_mg_per_ml"]
        volume_ul = measured.get("volume_ul")
        rows.append({
            "source": str(label),
            "lot": inputs["fields"].get("lot") or "",
            "date": inputs["fields"].get("date") or "",
            "name": inputs["name"],
            "prepared_by": inputs["fields"].get("prepared_by") or "",
            "host": inputs["fields"].get("host") or "",
            "purification": inputs["fields"].get("purification") or "",
            "buffer": inputs["fields"].get("buffer") or "",
            "storage": inputs["fields"].get("storage") or "",
            "concentration_mg_per_ml": concentration,
            "basis": basis,
            "basis_detail": basis_detail,
            "volume_ul": volume_ul,
            "total_mg": (round(concentration * volume_ul / 1000.0, 4)
                         if (concentration and volume_ul) else None),
            "measured_da": measured.get("measured_da"),
            "mass_difference_da": sheet["mass_difference_da"],
            "sequence": inputs["sequence"],
            "created_with": inputs["created_with"],
        })

    # Lots grouped by the sequence they state, in the order the sequences
    # first appear. More than one group means the comparison spans more than
    # one protein, which is worth doing deliberately and not by accident.
    groups: dict[str, list[str]] = {}
    for row in rows:
        if row["sequence"]:
            groups.setdefault(row["sequence"], []).append(
                row["lot"] or row["source"])
    sequence_groups = list(groups.values())
    if len(sequence_groups) > 1:
        warnings.append(
            "These records do not all carry the same sequence, so the rows "
            "below describe " + str(len(sequence_groups)) + " different "
            "proteins. The molar mass and extinction coefficient differ "
            "between them, which means the concentrations are not comparable "
            "even where the readings look similar.")

    bases = {row["basis"] for row in rows if row["basis"]}
    if len(bases) > 1:
        warnings.append(
            "Some concentrations here were read by a method that does not go "
            "through the sequence and others were derived from an A280 with a "
            "predicted extinction coefficient. Comparing the two compares two "
            "methods as much as two lots: the basis column says which is "
            "which.")

    missing = [row["lot"] or row["source"] for row in rows
               if row["concentration_mg_per_ml"] is None]
    if missing:
        warnings.append(
            "No concentration could be put on " + ", ".join(missing) +
            ". Either the record carries no reading, or it carries an A280 "
            "and no sequence to divide it by.")

    stamped = {row["created_with"] for row in rows if row["created_with"]}
    older = sorted(stamp for stamp in stamped if version and stamp != version)
    if older:
        warnings.append(
            "Some of these were written by another version of the app (" +
            ", ".join(older) + f", against {version} running now). The "
            f"derived columns were recomputed here, so a figure that "
            f"disagrees with the printed sheet is a real disagreement and "
            f"worth looking into rather than a rounding difference.")

    lots = [row["lot"] for row in rows if row["lot"]]
    repeated = sorted({lot for lot in lots if lots.count(lot) > 1})
    if repeated:
        warnings.append(
            "More than one record calls itself " + ", ".join(repeated) +
            ". Nothing here can tell which is which, so the source column is "
            "the only thing separating them.")

    return {
        "rows": rows,
        "rejected": rejected,
        "warnings": warnings,
        "sequence_groups": sequence_groups,
    }


def build_batch_pdf(path, sheet: dict, version: str = "", author: str = "",
                    app_name: str = "VHH tool box", doi: str = "") -> "Path":
    """Render a batch sheet to a one-page PDF.

    Deliberately plain: a header, three blocks, the sequence, and a line
    saying where the theoretical numbers came from. A tube label does not
    need a cover page.
    """
    from datetime import datetime
    from pathlib import Path as _Path

    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (Paragraph, SimpleDocTemplate, Spacer, Table,
                                    TableStyle)

    path = _Path(path)
    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["BodyText"], fontSize=9,
                          leading=12)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=7.4,
                          leading=10)
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontSize=15,
                        spaceAfter=2)
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontSize=11,
                        spaceBefore=11, spaceAfter=3)
    small = ParagraphStyle("small", parent=body, fontSize=7.2,
                           textColor=colors.grey, leading=9.5)

    doc = SimpleDocTemplate(
        str(path), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm,
        topMargin=15 * mm, bottomMargin=15 * mm,
        title=f"Batch sheet, {sheet['name']}", author=author or app_name)

    def block(rows):
        table = Table([[label, value] for label, value in rows],
                      colWidths=[55 * mm, 119 * mm], hAlign="LEFT")
        table.setStyle(TableStyle([
            ("FONTSIZE", (0, 0), (-1, -1), 8.5),
            ("LEADING", (0, 0), (-1, -1), 11),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cccccc")),
            ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#eef1f5")),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        return table

    story = [Paragraph(f"Batch sheet: {sheet['name']}", h1)]
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    line = f"Generated {stamp}"
    if version:
        line += f" by {app_name} {version}"
    story += [Paragraph(line, small), Spacer(1, 6)]

    story.append(Paragraph("Preparation", h2))
    story.append(block([(label, value or "-")
                        for label, value, _key in sheet["fields"]]))

    if sheet["measured"]:
        story.append(Paragraph("Measured", h2))
        story.append(block(sheet["measured"]))
        story.append(Spacer(1, 3))
        story.append(Paragraph(
            "Readings only. Nothing in this block depends on the sequence.",
            small))

    if sheet.get("derived"):
        story.append(Paragraph("Derived from both", h2))
        story.append(block(sheet["derived"]))
        story.append(Spacer(1, 3))
        for note in sheet.get("derivation_notes") or []:
            story.append(Paragraph(note, small))

    story.append(Paragraph("From the sequence", h2))
    story.append(block(sheet["theoretical"]))
    story.append(Spacer(1, 3))
    story.append(Paragraph(
        "Every value in this block is theoretical, computed from the sequence "
        "alone. The extinction coefficient and the pI in particular are "
        "predictions; confirm them experimentally before relying on them.",
        small))

    if sheet["sequence"]:
        story.append(Paragraph("Sequence", h2))
        wrapped = "<br/>".join(sheet["sequence"][i:i + 60]
                               for i in range(0, len(sheet["sequence"]), 60))
        story.append(Paragraph(wrapped, mono))

    footer = ("This is a laboratory batch record, not a certificate of "
              "analysis: it carries no specifications, no acceptance criteria "
              "and no release decision.")
    if doi:
        footer += f" Tool: https://doi.org/{doi}"
    story += [Spacer(1, 8), Paragraph(footer, small)]

    doc.build(story)
    return path
