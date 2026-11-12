"""
provenance - every constant the app computes with, and where it came from.

The app already refuses to invent a constant: each one sits beside a comment
saying what it was checked against, and the November 2026 audit went through
them one by one. That work was invisible to a user. It lived in comments, in
captions on three different screens and in two sections of the manual, and
somebody deciding whether to trust a number had to find the right caption or
read the source.

This module is the one place that information is enumerated, and the screen
that shows it (`sections/provenance.py`) reads nothing else.

**The value is never typed here.** Each entry names the module and the
attribute that owns the constant, and the value is read from the live module
at display time. A constant changed in `bench.py` changes on this screen with
no edit here, and an entry naming an attribute that no longer exists fails a
test rather than showing a stale number. `tests/test_provenance.py` parses
this file and fails if a numeric literal appears in any entry, which is what
keeps that property true rather than merely intended.

**What this screen cannot do**, stated because the opposite would be easy to
assume: it cannot verify that the source text still describes the value. The
text is prose, written by a person, and the test suite can check that it
exists and that it is attached to a constant that exists. If somebody changes
a threshold and not the sentence under it, this screen will show the new
number with the old justification. Only reading the two together catches that,
which is why the basis vocabulary below is deliberately blunt about which
constants are conventions rather than measurements: a convention that moves is
a decision, and a published value that moves is a mistake.
"""

from __future__ import annotations

import dataclasses
import importlib

__version__ = "2026.11.12c"

__all__ = [
    "BASES", "Constant", "entries", "resolve", "value_text", "by_basis",
    "__version__",
]


# How much weight a number carries, in six kinds, from the best understood to
# the least. The vocabulary is small on purpose: "theoretical", "derived" and "measured"
# already mean something precise in this app, and they describe *results*. A
# constant is not a result, so it gets its own words, and the point of each is
# to answer one question, what would it take to change this number.
BASES = {
    "published": "From a publication, checked against it. Changing it means "
                 "the reading of the paper was wrong.",
    "vendor": "A manufacturer's own figure for a product. Lot-specific in "
              "principle, so the certificate of analysis for the vial you "
              "hold outranks it.",
    "standard": "A value every table agrees on, with no single paper behind "
                "it in this app's reading.",
    "convention": "A choice, not a finding. Defensible, widely used, and "
                  "nothing in this app or in the literature fixes it at this "
                  "value rather than a nearby one.",
    "observed here": "Set from a measurement made while building this app, "
                     "on a small number of cases named in the source. Not a "
                     "validated parameter.",
    "unverified": "In wide circulation and traceable to no primary source "
                  "this project could find. Used because the arithmetic is "
                  "standard, and shown beside an alternative route wherever "
                  "the app can offer one.",
}


@dataclasses.dataclass(frozen=True)
class Constant:
    """One constant, named by where it lives rather than by its value.

    ``where`` is the import path to the value: a module name followed by an
    attribute, then any number of mapping keys. ``("bench", "REAGENTS",
    "NHS-PEG4-Biotin", "mw")`` and ``("structure", "BURIED_MAX")`` are both
    valid, and `resolve` walks them the same way.

    ``used_for`` names the screen the constant reaches, because "where does
    this bite" is the question a user asks next and the answer is not derivable
    from the module it lives in: `vhh_core.CHOU_FASMAN` feeds both a profile in
    the analysis and the aggregation score beside it.
    """

    key: str
    label: str
    where: tuple
    unit: str
    basis: str
    source: str
    used_for: str
    note: str = ""


def resolve(entry: Constant):
    """The constant's current value, read from the module that owns it."""
    module_name, attribute, *path = entry.where
    value = getattr(importlib.import_module(module_name), attribute)
    for step in path:
        value = value[step]
    return value


def value_text(value) -> str:
    """A constant formatted for a table cell, whatever kind of thing it is.

    A table is reported as its size rather than printed: the Tien accessibility
    values and the Chou-Fasman parameters are twenty rows each, and a cell
    holding twenty numbers is a cell nobody reads. The screen links to the
    reference tables, which draw them properly.
    """
    from pathlib import Path

    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int, float)):
        return f"{value:g}"
    if isinstance(value, Path):
        return value.name
    if isinstance(value, dict):
        return f"{len(value)} entries"
    if isinstance(value, (set, frozenset)):
        return f"{len(value)} values"
    if isinstance(value, (list, tuple)):
        return ", ".join(value_text(item) for item in value)
    return str(value)


# --------------------------------------------------------------------------
# The registry
# --------------------------------------------------------------------------
# Ordered by where a user meets them, not by module. Every `source` is the
# standing of the number as the code beside it records it, and several of them
# say plainly that the attribution could not be confirmed: those sentences are
# the useful ones, and softening them here would undo the audit they came from.

_REGISTRY = (
    # --- sequence analysis -------------------------------------------------
    Constant(
        key="chou_fasman",
        label="Chou-Fasman helix and sheet propensities",
        where=("vhh_core", "CHOU_FASMAN"),
        unit="propensity, scaled by 100",
        basis="published",
        source="Chou and Fasman, Adv. Enzymol. 47:45-148 (1978), as published "
               "by Expasy ProtScale. The 1978 chapter itself was not read. "
               "Nineteen of the twenty helix values and all twenty sheet "
               "values agree across every source checked; lysine's helix "
               "value was corrected from 114 to 116, where Expasy, the MATLAB "
               "Bioinformatics Toolbox and the Chou-Fasman teaching tables "
               "all agree against the one table it came from.",
        used_for="Sequence analysis: the secondary-structure propensity "
                 "profile, and the beta term of the aggregation score.",
    ),
    Constant(
        key="kyte_doolittle",
        label="Kyte and Doolittle hydropathy scale",
        where=("vhh_core", "KD"),
        unit="hydropathy index",
        basis="standard",
        source="The Kyte and Doolittle scale as it appears in every "
               "implementation checked, including Biopython's own ProtParam, "
               "which the app uses for GRAVY.",
        used_for="Sequence analysis: the hydropathy profile and the "
                 "hydrophobic term of the aggregation score.",
    ),
    Constant(
        key="imgt_cdr",
        label="IMGT CDR boundaries",
        where=("vhh_core", "IMGT_CDR"),
        unit="IMGT position ranges",
        basis="published",
        source="The IMGT unique numbering for V domains, as ANARCI applies "
               "it: CDR1 27 to 38, CDR2 56 to 65, CDR3 105 to 117. These are "
               "positions in the numbering, not residue counts, which is why "
               "a CDR1 annotated by ANARCI and one annotated by the regex "
               "anchors have different lengths.",
        used_for="Sequence analysis and Structure: CDR and framework "
                 "annotation whenever ANARCI runs.",
    ),
    Constant(
        key="tnp_thresholds",
        label="Therapeutic Nanobody Profiler flagging thresholds",
        where=("tnp", "METRICS"),
        unit="per metric: residues, residues per angstrom, or patch score",
        basis="published",
        source="Table 2 of Gordon GL, Gervasio J, Souders C and Deane CM, "
               "Commun Biol 9:344 (2026), doi 10.1038/s42003-026-09594-y, "
               "transcribed from the published PDF and checked a second time "
               "in tests/test_tnp.py. The ranges are those spanned by the 36 "
               "clinical-stage nanobody sequences available to the authors, "
               "and the amber bands the top and bottom 5% of each "
               "distribution. The paper states that they will be revised as "
               "more candidates progress through trials, so unlike most "
               "entries here this one is expected to change: "
               "tnp.THRESHOLD_DATE carries the calibration date and the app "
               "shows it beside the flags.",
        used_for="Sequence analysis: the clinical context tab. Structure: the "
                 "flag on CDR3 compactness.",
        note="Four of the six metrics need a structure, and three of those "
             "are not computed anywhere in this app. Their thresholds are "
             "carried anyway, so the screen can show what it is not "
             "answering.",
    ),
    Constant(
        key="sabdab_stale_days",
        label="When the structure snapshot is called old",
        where=("sabdab", "STALE_DAYS"),
        unit="days",
        basis="convention",
        source="Hand set. SAbDab2 adds entries weekly, so this many days is "
               "roughly thirteen updates missed, which was judged enough "
               "that somebody deciding whether to crystallise a candidate "
               "should check the live database rather than trust the shipped "
               "file. Nothing fixes it here rather than at sixty or a "
               "hundred and twenty days.",
        used_for="Sequence analysis: whether the closest-structure search "
                 "warns that its snapshot is out of date.",
        note="The extraction date it is measured against is read from "
             "data/sabdab/meta.json, which the build script writes from the "
             "export's own update dates rather than from the clock.",
    ),
    Constant(
        key="apr_threshold",
        label="Aggregation-prone window threshold",
        where=("vhh_core", "APR_THRESHOLD"),
        unit="smoothed score",
        basis="convention",
        source="Hand set. The score it cuts is this app's own combination of "
               "charge, hydrophobicity and beta propensity, so no published "
               "threshold applies to it, and it is a slider in the app for "
               "that reason.",
        used_for="Sequence analysis: which windows are called "
                 "aggregation-prone.",
        note="Not SAP, not CamSol, not TANGO. A triage heuristic.",
    ),
    Constant(
        key="side_chain_pka",
        label="Side-chain pKa",
        where=("reference", "SIDE_CHAIN_PKA"),
        unit="pKa",
        basis="standard",
        source="The values of the standard free-amino-acid tables. Entered as "
               "Nelson and Cox, Lehninger 7th ed., which was not consulted to "
               "confirm them. Six of the seven agree across every table "
               "checked; cysteine agrees too, but at least one teaching table "
               "gives 8.08 rather than 8.18, so the second decimal is not "
               "something to lean on.",
        used_for="Reference tables: the amino-acid screen. Display only, "
                 "nothing computes with it, and the pI the app reports comes "
                 "from ProtParam's own constants.",
    ),
    # --- DNA ---------------------------------------------------------------
    Constant(
        key="codon_usage",
        label="E. coli codon usage table",
        where=("codons", "USAGE_FILE"),
        unit="relative frequency per codon",
        basis="published",
        source="Kazusa codon usage database, taxid 316407, CC0. This is a "
               "whole-organism table. Sharp and Li (1987) define the CAI "
               "against a reference set of highly expressed genes, so the CAI "
               "this app reports uses their arithmetic and not their "
               "reference set, and is not comparable with a CAI quoted from "
               "elsewhere.",
        used_for="Gene diagnosis: codon frequencies, the CAI, and which "
                 "synonymous swaps are offered.",
    ),
    Constant(
        key="min_usage",
        label="Floor on a synonymous substitute's share",
        where=("codons", "MIN_USAGE"),
        unit="fraction of its codon family",
        basis="observed here",
        source="Chosen by measurement on test genes rather than by feel: a 20 "
               "per cent floor clears every GC-rich stretch while costing "
               "0.04 of CAI, where a 10 per cent floor costs 0.08 by letting "
               "leucine's TTA, at 13 per cent, replace CTG, at 50 per cent.",
        used_for="Gene diagnosis: which synonymous swaps are offered.",
    ),
    Constant(
        key="rare_codons",
        label="Rare codons, and the tRNA each is short of",
        where=("codons", "RARE_CODONS"),
        unit="codons",
        basis="convention",
        source="The codons the Rosetta and CodonPlus strains supplement, "
               "which is a statement about tRNA supplementation and not about "
               "frequency. Cross-checked against the measured table: five are "
               "genuinely scarce, three are not rare by frequency at all, and "
               "two carry no documented limiting tRNA. The measured frequency "
               "is reported beside each so the two can be judged separately.",
        used_for="Gene diagnosis and Reference tables: the rare-codon flags "
                 "and the genetic code table.",
    ),
    # --- structure ---------------------------------------------------------
    Constant(
        key="tien_max_asa",
        label="Maximum accessible surface area per residue",
        where=("structure", "TIEN_MAX_ASA"),
        unit="square angstrom",
        basis="published",
        source="Tien et al. 2013, Table 1, theoretical column, verified "
               "against the published table rather than reproduced from "
               "memory. The theoretical column is the one the paper "
               "recommends as the upper bound; its empirical column is "
               "systematically a few per cent lower, which would raise every "
               "relative accessibility this app reports.",
        used_for="Structure: relative accessibility, and therefore which "
                 "residues are called buried or exposed.",
    ),
    Constant(
        key="buried_max",
        label="Relative accessibility at or below which a residue is buried",
        where=("structure", "BURIED_MAX"),
        unit="fraction of maximum",
        basis="convention",
        source="The usual convention. Nothing in the app or in Tien et al. "
               "fixes it at 0.20, and the band between this and the exposed "
               "cut-off is called intermediate rather than forced into one of "
               "the two.",
        used_for="Structure: the exposure table and the liabilities in 3D.",
    ),
    Constant(
        key="exposed_min",
        label="Relative accessibility above which a residue is exposed",
        where=("structure", "EXPOSED_MIN"),
        unit="fraction of maximum",
        basis="convention",
        source="The usual convention, and the counterpart of the buried "
               "cut-off above.",
        used_for="Structure: the exposure table and the liabilities in 3D.",
    ),
    Constant(
        key="ss_max_distance",
        label="Sulfur to sulfur distance that counts as a disulfide",
        where=("structure", "SS_MAX_DISTANCE"),
        unit="angstrom",
        basis="convention",
        source="A disulfide S-S bond is about 2.05 angstrom. This cut-off is "
               "comfortably past any real bond and still well short of the "
               "3.5 angstrom at which two unbonded sulfurs can sit.",
        used_for="Structure: which cysteines are reported as paired.",
    ),
    Constant(
        key="patch_min_rel",
        label="Exposure a residue needs to join a surface patch",
        where=("structure", "PATCH_MIN_REL"),
        unit="fraction of maximum",
        basis="observed here",
        source="Checked against one complex, the HIV-1 capsid C-terminal "
               "domain, where the known W184/M185 dimerisation patch is found "
               "at 0.15 and missed at 0.20. That cliff is one residue wide, "
               "so it is an observation on a single structure and "
               "emphatically not a validated parameter. It is a slider in the "
               "app, and the app says so when it finds nothing.",
        used_for="Structure: the surface hydrophobic patches.",
    ),
    Constant(
        key="patch_contact",
        label="Distance at which two side chains count as adjacent",
        where=("structure", "PATCH_CONTACT"),
        unit="angstrom",
        basis="convention",
        source="Measured between closest heavy atoms rather than between "
               "centroids, because tryptophan and arginine carry their bulk "
               "far from the backbone and a centroid rule misses exactly the "
               "large residues a patch depends on. On the one complex checked, "
               "the result barely moves between 5.0 and 7.0 angstrom.",
        used_for="Structure: the surface hydrophobic patches.",
    ),
    Constant(
        key="interface_min_delta",
        label="Buried area that makes a residue part of an interface",
        where=("structure", "INTERFACE_MIN_DELTA"),
        unit="square angstrom",
        basis="convention",
        source="Small enough to catch a glancing contact and large enough to "
               "exclude numerical noise in the accessibility calculation.",
        used_for="Structure: the interface residues of a complex.",
    ),
    Constant(
        key="iptm_high",
        label="ipTM above which AlphaFold3 calls an interface confident",
        where=("af3", "IPTM_HIGH"),
        unit="ipTM",
        basis="convention",
        source="The thresholds AlphaFold3's own documentation uses when "
               "describing interface confidence. They band a continuous score, "
               "and a complex either side of one of them is not a different "
               "kind of prediction.",
        used_for="Structure: how an uploaded AlphaFold3 model's confidence is "
                 "described.",
    ),
    # --- bench -------------------------------------------------------------
    Constant(
        key="phage_optical_constant",
        label="Filamentous phage optical constant",
        where=("bench", "PHAGE_OPTICAL_CONSTANT"),
        unit="bases per mL per absorbance unit",
        basis="unverified",
        source="The optical route is virions/mL = A260 times this constant "
               "divided by the number of bases in the genome. It was carried "
               "here as Day's relation, and that attribution could not be "
               "traced to any primary source. The formula circulates in "
               "exactly this form in course protocols, published methods and "
               "vendor calculators, always uncited. The arithmetic is standard "
               "and reproducible; the constant's provenance is not "
               "established, which is why the empirical route is shown beside "
               "it rather than behind it.",
        used_for="Phage and panning: the titre from an A260.",
    ),
    Constant(
        key="phage_empirical_per_a260",
        label="Particles per mL per A260, empirical route",
        where=("bench", "PHAGE_EMPIRICAL_PER_A260"),
        unit="particles per mL per absorbance unit",
        basis="convention",
        source="The correspondence quoted in phage display protocols for a "
               "filamentous phage preparation. It is a correspondence measured "
               "on particular preparations, not a property of the particle, "
               "and it is shown beside the optical route so that a "
               "disagreement between the two is visible rather than resolved "
               "for you.",
        used_for="Phage and panning: the titre from an A260.",
    ),
    Constant(
        key="phage_genome_bases",
        label="Phage genome length used by the optical route",
        where=("bench", "PHAGE_GENOME_BASES"),
        unit="bases",
        basis="standard",
        source="M13mp18 and its phagemid derivatives. A phagemid carrying an "
               "insert is shorter than the helper phage, which is why this is "
               "an input in the app rather than a fixed constant.",
        used_for="Phage and panning: the titre from an A260.",
    ),
    Constant(
        key="liquid_panning_factor",
        label="Input multiplier for a liquid-phase panning round",
        where=("bench", "LIQUID_PANNING_FACTOR"),
        unit="times the solid-phase input",
        basis="convention",
        source="A protocol convention: a liquid-phase round is run with twice "
               "the input of a solid-phase one, because the antigen is in "
               "solution and the capture step loses more of the library. It is "
               "not derived from anything in this app.",
        used_for="Phage and panning: how much stock a round takes.",
    ),
    Constant(
        key="pipette_minimum",
        label="Smallest volume the conjugation plan will ask for",
        where=("bench", "PIPETTE_MINIMUM_UL"),
        unit="microlitre",
        basis="convention",
        source="Below this, a volume cannot be pipetted accurately enough for "
               "the molar excess to mean anything, and the stock should be "
               "diluted first. The Thermo cheat sheet makes the same point "
               "with a fixed one in ten dilution.",
        used_for="Batch and QC: the conjugation plan.",
    ),
    Constant(
        key="atomic_weights",
        label="Standard atomic weights",
        where=("buffers", "ATOMIC_WEIGHTS"),
        unit="dalton",
        basis="published",
        source="",       # filled from buffers.ATOMIC_WEIGHTS_SOURCE, below
        used_for="Buffers and media: every molar mass, which is computed from "
                 "the formula of the form you weigh rather than looked up.",
    ),
)


def entries() -> list[Constant]:
    """Every constant on the register, with the vendor reagents derived.

    The reagent molar masses are not written out above. They live in
    `bench.REAGENTS`, each already beside the source it was checked against,
    so they are read from there: adding a reagent to that table adds a row to
    this screen, and a reagent whose source string somebody deletes shows up
    as a missing source rather than as a row that quietly stops being listed.

    The atomic weights are the same case in miniature: `buffers` already
    publishes its source as a constant, so the entry above carries no source
    text of its own and takes it from there.
    """
    import bench as bn
    import buffers as bf

    out = []
    for entry in _REGISTRY:
        if entry.key == "atomic_weights":
            entry = dataclasses.replace(entry, source=bf.ATOMIC_WEIGHTS_SOURCE)
        out.append(entry)

    for name, reagent in bn.REAGENTS.items():
        out.append(Constant(
            key=f"reagent:{name}",
            label=f"{name}: molar mass"
                  + (", extinction coefficient and A280 correction"
                     if reagent.get("epsilon") else ""),
            where=("bench", "REAGENTS", name, "mw"),
            unit="dalton",
            basis="vendor",
            source=reagent.get("source") or "",
            used_for="Batch and QC: the conjugation plan.",
            note=reagent.get("note") or "",
        ))
    return out


def by_basis() -> dict[str, list[Constant]]:
    """The register grouped by how much weight each constant carries.

    Grouped in `BASES` order, which runs from the best understood kind to the
    least. A basis with no constant is left out rather than drawn empty, and a
    basis that appears in no entry is therefore invisible rather than an empty
    heading nobody can interpret.
    """
    grouped = {basis: [entry for entry in entries() if entry.basis == basis]
               for basis in BASES}
    return {basis: found for basis, found in grouped.items() if found}
