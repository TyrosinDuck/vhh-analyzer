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
    "PIPETTE_MINIMUM_UL", "desalting_column", "__version__",
]

__version__ = "2026.10.30c"


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
                  "Table 4: MW 643, eps 71,000 /M/cm at 495 nm, CF280 0.11.",
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

# The optical route: Day's relation for filamentous phage,
#   virions/mL = A269 x 6e16 / (number of bases in the genome)
# which is the ssDNA absorbance of the particle worked out from its genome
# length. Reading at 260 rather than 269 nm, and against a blank rather than
# with the A320 scatter correction, costs some accuracy but not an order of
# magnitude.
PHAGE_OPTICAL_CONSTANT = 6.0e16

# M13mp18 and its phagemid derivatives. A phagemid carrying an insert is
# shorter than the helper phage, which is why this is an input rather than a
# constant.
PHAGE_GENOME_BASES = 7249


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
