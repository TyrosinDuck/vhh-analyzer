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

**The register is checked for completeness, not merely assembled.** Until
November 2026 it listed the constants somebody had thought to add, and the
manual claimed it listed every constant the app computes with. Those were not
the same set: a dozen numbers that change what a user reads, the lower ipTM
band and the histidine charge among them, were nowhere on it. `NOT_SCIENTIFIC`
below names the constants that carry no scientific claim, each with the reason,
and `tests/test_provenance.py` walks the modules `app.py` requires and fails on
any public constant that is neither registered nor named there. A new threshold
now has to be classified before the suite passes, which is the only way a
claim of completeness stays true.

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

__version__ = "2026.11.13a"

__all__ = [
    "BASES", "Constant", "NOT_SCIENTIFIC", "entries", "resolve", "value_text",
    "by_basis", "scan_module", "unclassified", "__version__",
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
        # A table of records is reported as its size, the same way a mapping is:
        # the mass-difference modifications and the liability rules are a dict
        # or a tuple each, and flattening one into a cell gives a line of
        # fragments nobody reads. A flat list of numbers is still printed.
        if any(isinstance(item, (dict, list, tuple, set, frozenset))
               for item in value):
            return f"{len(value)} entries"
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
        key="hallmark_positions",
        label="FR2 hallmark positions, with the residues read as camelid",
        where=("vhh_core", "HALLMARK_POSITIONS"),
        unit="IMGT positions, and a residue set per position",
        basis="convention",
        source="The positions are published; the residue sets are not. IMGT "
               "42, 49, 50 and 52, equivalent to Kabat 37, 44, 45 and 47, are "
               "the tetrad every description of the VHH framework names, and "
               "the human VH3 consensus carried beside each (V, G, L, W) is "
               "the standard reading. Which residues at those positions count "
               "as camelid is this app's own list and is wider than any one "
               "publication gives: position 49 in particular accepts G, which "
               "VH3 also carries, and the code records it as the least "
               "discriminating of the four for that reason. A sequence is "
               "never called a VHH on the tetrad alone.",
        used_for="Sequence analysis: the hallmark table and the VHH versus VH "
                 "verdict. Library: the same tetrad when a generated "
                 "framework is checked.",
        note="generator.HALLMARK_POSITIONS and generator.VH_LIKE are derived "
             "from this tuple rather than written out again, so the library "
             "generator and the analysis cannot drift apart.",
    ),
    Constant(
        key="region_bounds",
        label="IMGT region boundaries of the generated framework",
        where=("generator", "REGION_BOUNDS"),
        unit="IMGT position ranges",
        basis="standard",
        source="The IMGT region boundaries, frameworks and CDRs alike, over "
               "a domain running to position 128. The three CDR windows are "
               "the same ones vhh_core.IMGT_CDR carries and the frameworks "
               "are the stretches between them, so this is the one place in "
               "the app where those boundaries are written down twice; they "
               "agree today, and nothing but this sentence and a test "
               "enforces that.",
        used_for="Library: the regions a generated sequence is assembled "
                 "from, and the position matrix it is sampled from.",
    ),
    Constant(
        key="domain_coverage",
        label="Fraction of the sequence the V domain must span",
        where=("vhh_core", "DOMAIN_COVERAGE"),
        unit="fraction of the submitted sequence",
        basis="convention",
        source="Hand set. Below this, ANARCI's numbered domain covers too "
               "little of what was submitted for the VHH analyses to be about "
               "the whole sequence, and the app reports a fusion or "
               "multi-domain protein instead of analysing it as a domain. "
               "Nothing fixes it at this value, and the coverage it was "
               "compared against is reported beside the verdict so the "
               "decision is visible.",
        used_for="Sequence analysis: whether a sequence is treated as a "
                 "single V domain or as a fusion.",
    ),
    Constant(
        key="canonical_pair",
        label="The canonical disulfide, by IMGT position",
        where=("vhh_core", "CANONICAL_PAIR"),
        unit="IMGT positions",
        basis="standard",
        source="The intradomain disulfide of an immunoglobulin V domain, at "
               "IMGT 23 and 104. Stated the same way in every IMGT-numbered "
               "description checked; no single paper behind it here.",
        used_for="Sequence analysis: which cysteines are the canonical pair "
                 "and which are candidate free thiols or a second, "
                 "stabilising pair.",
    ),
    Constant(
        key="conserved_imgt",
        label="Positions a repair suggestion must not touch",
        where=("vhh_core", "CONSERVED_IMGT"),
        unit="IMGT positions",
        basis="standard",
        source="The canonical disulfide cysteines at 23 and 104, the FR2 "
               "tryptophan at 41 that packs the domain core, and the FR4 "
               "tryptophan at 118 of the WGxG motif. Each is conserved in "
               "every description of the fold consulted. The list is "
               "deliberately short and is not a complete inventory of "
               "conserved positions: it exists to stop the oxidation rules, "
               "which match Met and Trp anywhere, from offering W41F or "
               "W118F as a repair.",
        used_for="Sequence analysis: which liability repairs are offered and "
                 "which are withheld.",
    ),
    Constant(
        key="family_cutoff",
        label="CDR3 identity at which two sequences are one family",
        where=("vhh_core", "FAMILY_DEFAULTS", "cutoff"),
        unit="per cent identity over CDR3",
        basis="convention",
        source="Hand set, and kept identical to the Sanger screening notebook "
               "so that the app and the notebook put the same reads in the "
               "same families. It is a clustering cut-off, not a measure of "
               "relatedness: a slider in the app, and the dendrogram is drawn "
               "so the effect of moving it is visible.",
        used_for="Sequence analysis: the CDR3 families and the dendrogram.",
        note="The identity behind it divides by the longer of the two CDR3s, "
             "which is the notebook's denominator and not the one "
             "pairwise_identity uses; the two agree for equal lengths and "
             "differ by a few per cent otherwise.",
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
        key="solubility_weights",
        label="Weights of the four terms in the aggregation score",
        where=("vhh_core", "SOLUBILITY_WEIGHTS"),
        unit="relative weight",
        basis="convention",
        source="A judgement call, and the code says so where they are "
               "defined. Each of the four terms (net charge, Kyte and "
               "Doolittle hydropathy, Chou-Fasman beta propensity, and the "
               "gatekeeper residues P, D, E, K and R) is individually well "
               "established as a determinant of aggregation propensity. The "
               "weights combining them are fitted to nothing. CamSol "
               "intrinsic, AGGRESCAN and TANGO fit coefficients to "
               "experimental data and this reproduces none of them, which is "
               "why the profile is for reading where along a sequence the "
               "prone stretches sit and for ranking sequences against each "
               "other, never as a predicted solubility.",
        used_for="Sequence analysis: the aggregation profile and the "
                 "aggregation score.",
    ),
    Constant(
        key="charge_table",
        label="Per-residue charges used by the profiles",
        where=("vhh_core", "CHARGE"),
        unit="elementary charge",
        basis="convention",
        source="Formal charges at neutral pH, hand set: minus one for D and "
               "E, plus one for K and R, and a fractional value for H as a "
               "stand-in for its partial protonation near pH 7.4. The "
               "histidine value has no source and a pKa-based calculation "
               "would not give it; it is small on purpose, so a "
               "histidine-rich stretch is neither ignored nor counted as a "
               "full charge. Nothing quantitative rests on it: the net charge "
               "and the isoelectric point the app reports come from "
               "Biopython's ProtParam and its own constants, not from this "
               "table, which sets the shape of a profile.",
        used_for="Sequence analysis: the charge profile, and the charge term "
                 "of the aggregation score.",
        note="The aggregation score takes the absolute value of these, "
             "because charge of either sign opposes aggregation. The profile "
             "keeps the sign.",
    ),
    Constant(
        key="liability_rules",
        label="Chemical liability motifs and their base scores",
        where=("vhh_core", "LIABILITY_RULES"),
        unit="motif, regular expression and base score",
        basis="convention",
        source="The motifs are standard: the N-x-S/T sequon, NG and the "
               "slower N-S/T/N contexts, DG and the slower D-S/T/D contexts, "
               "methionine and tryptophan oxidation, the acid-labile DP bond "
               "and the RGD integrin motif are described as liabilities "
               "throughout the antibody engineering literature. The numbers "
               "attached to them are not published. The entries marked "
               "notebook in the table carry the base scores of an in-house "
               "liability notebook, where each motif is weighted by its "
               "solvent accessibility and the pLDDT of an AlphaFold model; "
               "this module is sequence-only, so neither weighting is "
               "applied here and the score is the base score alone. The "
               "entries marked added are a judgement call and are meant to be "
               "edited.",
        used_for="Sequence analysis: the liabilities table, the risk score "
                 "and the repair suggestions. Structure: the same motifs "
                 "weighted by exposure, where a structure is supplied.",
        note="Because the exposure weighting is absent, the scores compress "
             "towards the top: an NG scores the same buried or fully "
             "exposed. The ranking between sequences carries more than the "
             "band any single sequence falls in.",
    ),
    Constant(
        key="risk_thresholds",
        label="Score bands for low, medium and high liability risk",
        where=("vhh_core", "RISK_THRESHOLDS"),
        unit="liability score",
        basis="convention",
        source="Hand set, from the same in-house notebook as the base "
               "scores, and inherited together with them. They band a sum of "
               "scores that are themselves a judgement call, so the band is "
               "the least load-bearing thing on that screen: the motifs and "
               "their positions are what to read.",
        used_for="Sequence analysis: the Low, Medium and High label on the "
                 "liabilities table and in triage.",
    ),
    Constant(
        key="protease_rules",
        label="Protease recognition sites and their per-site scores",
        where=("vhh_core", "PROTEASE_RULES"),
        unit="site, regular expression and per-site score",
        basis="convention",
        source="The specificities are textbook (trypsin after K and R but "
               "not before P, chymotrypsin after F, W and Y, GluC after E, "
               "AspN before D, and the IEGR and DDDDK recognition sites). The "
               "per-site scores are from the in-house notebook and are "
               "deliberately small, because a VHH carries many tryptic sites "
               "and their number on its own means very little; the two "
               "engineered sites score higher because finding one in a VHH is "
               "a real finding.",
        used_for="Sequence analysis: the protease site table.",
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
    Constant(
        key="rare_threshold",
        label="Frequency at or below which a residue is called unusual",
        where=("repertoire", "RARE_THRESHOLD"),
        unit="fraction of the reference sequences at that position",
        basis="convention",
        source="Hand set against the bundled PLAbDab-nano position matrix. "
               "It is a frequency in that reference repertoire and not a "
               "probability of anything: a residue below it may be a Sanger "
               "miscall, a PCR artefact or a genuine novelty, and the app "
               "reports the frequency beside the flag rather than deciding "
               "which. A slider in the app.",
        used_for="Repertoire: which positions are flagged as carrying an "
                 "unusual residue.",
    ),
    Constant(
        key="min_support",
        label="Reference sequences a position needs before it is judged",
        where=("repertoire", "MIN_SUPPORT"),
        unit="sequences",
        basis="convention",
        source="Hand set. A position covered by a handful of reference "
               "sequences has no usable frequency, and calling a residue rare "
               "there would be reporting the thinness of the reference as a "
               "property of the sequence. Positions below this are set aside "
               "and counted rather than flagged, which matters most in long "
               "CDR3s, where the reference thins out exactly where the "
               "sequences differ.",
        used_for="Repertoire: which positions are judged at all.",
    ),
    Constant(
        key="min_for_bars",
        label="Set size below which lengths are drawn as markers",
        where=("repertoire", "MIN_FOR_BARS"),
        unit="submitted sequences",
        basis="convention",
        source="A drawing rule with a statistical reason. Below this many "
               "sequences the CDR length comparison is drawn as markers "
               "placed on the reference distribution instead of as a "
               "histogram, because a five-bar histogram invites reading a "
               "shape that five sequences cannot support. Nothing computes "
               "with it.",
        used_for="Repertoire: how the CDR length comparison is drawn.",
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
    Constant(
        key="cluster_window",
        label="Window in which two rare codons count as a cluster",
        where=("codons", "CLUSTER_WINDOW"),
        unit="codons",
        basis="convention",
        source="Hand set. Rare codons close together are held to slow "
               "translation more than the same codons spread out, which is "
               "why a cluster is reported at all, but no published window "
               "defines this one. The individual rare codons are listed "
               "whether or not they fall in a cluster.",
        used_for="Gene diagnosis: the rare-codon cluster flags.",
    ),
    Constant(
        key="gc_window",
        label="Window over which local GC content is measured",
        where=("codons", "GC_WINDOW"),
        unit="nucleotides",
        basis="convention",
        source="Hand set, around the length of a synthesis oligo, so that a "
               "GC excursion the window finds is one a synthesis vendor "
               "would also meet. A different window moves the boundaries of "
               "every region it reports.",
        used_for="Gene diagnosis: the local GC profile and the GC-rich or "
                 "GC-poor regions.",
    ),
    Constant(
        key="gc_low",
        label="Local GC content below which a region is flagged",
        where=("codons", "GC_LOW"),
        unit="fraction",
        basis="convention",
        source="A practical synthesis band, not a biological one, and the "
               "code says so where it is defined. Vendors differ on what "
               "they will quote for, and a sequence inside the band can "
               "still be hard to synthesise for other reasons, such as a "
               "hairpin or a repeat.",
        used_for="Gene diagnosis: the GC-poor regions and the synthesis "
                 "warnings.",
    ),
    Constant(
        key="gc_high",
        label="Local GC content above which a region is flagged",
        where=("codons", "GC_HIGH"),
        unit="fraction",
        basis="convention",
        source="The upper edge of the same practical synthesis band, and the "
               "one the synonymous-swap suggestions work to bring a region "
               "back under.",
        used_for="Gene diagnosis: the GC-rich regions, the synthesis "
                 "warnings and the swaps offered.",
    ),
    Constant(
        key="head_codons",
        label="Length of the translation initiation region",
        where=("codons", "HEAD_CODONS"),
        unit="codons",
        basis="convention",
        source="The first stretch of the gene, taken as the region where "
               "codon choice and secondary structure bear on initiation "
               "rather than on elongation. Hand set: descriptions of the "
               "initiation region put it at a few tens of codons without "
               "agreeing on a number, and this app reports what it finds "
               "there separately rather than scoring it.",
        used_for="Gene diagnosis: the separate report on the start of the "
                 "gene.",
    ),
    Constant(
        key="sanger_defaults",
        label="Read trimming and quality thresholds",
        where=("sanger", "DEFAULTS"),
        unit="error probability, bases, Phred score",
        basis="convention",
        source="The Mott trimming cutoff is the value the phred and Staden "
               "implementations use as their default, about Q13, and the "
               "algorithm is Mott's modified trimming as those tools apply "
               "it. The minimum trimmed length and the mean Phred floor are "
               "hand-set QC limits: a read below either is flagged, never "
               "discarded, and the trimmed window is drawn so the call can "
               "be checked by eye.",
        used_for="Sequencing: which reads are flagged on quality, and where "
                 "each read is trimmed.",
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
        key="patch_min_size",
        label="Residues a surface patch must contain",
        where=("structure", "PATCH_MIN_SIZE"),
        unit="residues",
        basis="convention",
        source="Fewer than three adjacent exposed hydrophobic residues is "
               "not called a patch. Two such residues side by side are "
               "common on any folded domain and say little, so the floor is "
               "there to keep the list of patches readable rather than to "
               "mark a physical boundary. Raising it drops small patches "
               "silently, which is why it sits beside the exposure cut-off "
               "in the app rather than inside the code alone.",
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
    Constant(
        key="iptm_low",
        label="ipTM below which AlphaFold3 calls an interface unreliable",
        where=("af3", "IPTM_LOW"),
        unit="ipTM",
        basis="convention",
        source="The lower of the same pair of thresholds as the entry above, "
               "and from the same source: AlphaFold3's own documentation on "
               "interface confidence. The two band one continuous score into "
               "three, and the middle band is reported as what it is, a "
               "prediction the model itself does not commit to.",
        used_for="Structure: how an uploaded AlphaFold3 model's confidence is "
                 "described.",
    ),
    Constant(
        key="plddt_bands",
        label="pLDDT confidence bands used to colour a model",
        where=("structure", "PLDDT_COLOURS"),
        unit="pLDDT band floors, with the colour of each",
        basis="standard",
        source="The bands the AlphaFold Protein Structure Database uses for "
               "its own colouring, very high above 90, confident from 70 to "
               "90, low from 50 to 70 and very low below 50, with that "
               "database's four colours so that a model shown here looks "
               "like the same model shown there. Taken from the database's "
               "description of its colouring rather than from the AlphaFold "
               "papers. They band a continuous per-residue score and are not "
               "a measurement of anything.",
        used_for="Structure: the colouring of an uploaded predicted model in "
                 "the 3D viewer.",
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
        key="mass_modifications",
        label="Mass differences a modification would explain",
        where=("bench", "MODIFICATIONS"),
        unit="dalton, computed from each formula",
        basis="standard",
        source="No mass delta is tabulated here. Each one is computed from "
               "the chemical formula of the modification, with negative "
               "counts for a loss, so what a user can check is the formula "
               "rather than a number somebody typed: the classic plus 178 of "
               "a gluconoylated His-tag is C6H10O6 added, and it is right or "
               "wrong for a reason anybody can follow. The monoisotopic atom "
               "masses behind them are the values of the standard tables, "
               "which agree across every table consulted and were not "
               "checked here against a primary source. The conditions "
               "attached to each modification carry as much as the masses: "
               "an explanation that cannot apply to the sequence you "
               "submitted is not an explanation, so a modification is only "
               "offered when its condition is met.",
        used_for="Batch and QC: the explanations offered for a difference "
                 "between an observed mass and the sequence mass.",
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


# --------------------------------------------------------------------------
# The other half of the claim: what is deliberately not on the register
# --------------------------------------------------------------------------
# A register that lists what somebody remembered to add cannot claim to be
# complete, and the manual claimed exactly that. So every public constant in
# the modules `app.py` requires is either registered above or named here, with
# the reason it carries no scientific claim, and `tests/test_provenance.py`
# fails on anything that is neither.
#
# The test is what makes this a decision rather than a comment: a new threshold
# cannot reach the app without somebody saying which of the two it is. The cost
# is one line per guard rail, which is the right price for the claim.
#
# What belongs here: guard rails that stop a container being killed, cache
# settings, figure sizes, widget defaults and bounds, unit conversions that are
# exact by definition, and counts that are what a word means rather than
# something measured. What does not: anything a result depends on, however
# arbitrary it looks. An arbitrary threshold that changes a verdict belongs on
# the register, labelled a convention, which is what most of the register is.

NOT_SCIENTIFIC: dict[str, str] = {
    # Guard rails. Each of these has a measurement behind it, recorded beside
    # the constant, but what it bounds is this app's memory and not biology.
    "af3.MAX_MEMBER_BYTES": "ceiling on what is read out of an uploaded "
                            "archive, so a large model is refused rather than "
                            "exhausting the container",
    "compute.PROFILES_COMBINED_MAX": "above this many sequences the profile "
                                     "panels are rendered one at a time; a "
                                     "measured memory limit",
    "structure.MAX_ATOMS": "ceiling on the size of an uploaded structure",
    "vhh_core.MIN_LENGTH": "below this ProtParam's figures stop meaning much, "
                           "so the app declines to report them",
    "vhh_core.IDENTITY_FIG_MAX": "longest side of the identity heatmap, in "
                                 "inches; a measured rendering limit",
    "vhh_core.IDENTITY_LABEL_MAX": "above this many sequences the heatmap "
                                   "cells stop carrying a number",
    "shell.CACHE_TTL": "how long a cached result lives",
    "shell.CACHE_MAX_ENTRIES": "how many cached results are kept",

    # Units and counts. Exact by definition, or definitional in the word.
    "buffers.MOLAR": "molar to millimolar and micromolar; exact by definition",
    "buffers.MASS": "mass concentrations to grams per litre, per cent weight "
                    "per volume included; exact by definition",
    "buffers.VOLUME": "per cent volume per volume to itself; the entry exists "
                      "so a volume unit cannot be silently treated as a mass",
    "constructs.FORMATS": "how many VHH domains monovalent, bivalent and "
                          "biparatopic carry, which is what the words mean",

    # Widget defaults and bounds. A default is recorded with every run rather
    # than assumed, which is why it is not a hidden constant.
    "vhh_core.DEFAULTS": "the starting values of the sliders, including the "
                         "working pH and the profile windows. Every one is "
                         "written into the session file with the run, so a "
                         "figure can be reproduced at the value that made it",
    "session.METHODS_SETTINGS": "the bounds of those sliders, used to refuse "
                                "a session file asking for a value no widget "
                                "could show",
    "sabdab.DEFAULT_TOP": "how many nearest structures the table shows before "
                          "it is expanded",
    "vhh_core.PLOT_STYLE": "figure resolution and spines",

    # Rules about text and positions rather than quantities.
    "vhh_core.LIABILITY_REPAIRS": "the substitution offered for each "
                                  "liability, with its offset inside the "
                                  "motif. Positions and residues, with the "
                                  "reasoning in the table itself",
    "bench.DILUTION_LOOKS_LIKE_EXPONENT": "the point below which a dilution "
                                          "typed as a factor is probably an "
                                          "exponent. A warning about a typo, "
                                          "and never a refusal",
}


def scan_module(name: str) -> dict[str, int]:
    """Public constants of a module that hold a number, and where.

    Parsed from the source rather than read off the imported module, for two
    reasons: an attribute is not a line number, and the point of the scan is to
    name the line somebody has to go and classify. Reads module-level
    assignments to upper-case names, single and tuple targets alike, keeping
    those whose value contains a numeric literal anywhere inside it, which is
    what makes ``GC_LOW, GC_HIGH = 0.30, 0.70`` two findings and
    ``REGIONS = ["FR1", ...]`` none.
    """
    import ast
    from pathlib import Path

    path = Path(__file__).parent / f"{name}.py"
    found: dict[str, int] = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            targets = [node.target]
        else:
            continue

        names = []
        for target in targets:
            if isinstance(target, ast.Name):
                names.append(target.id)
            elif isinstance(target, ast.Tuple):
                names += [element.id for element in target.elts
                          if isinstance(element, ast.Name)]
        if not any(child for child in ast.walk(node.value)
                   if isinstance(child, ast.Constant)
                   and isinstance(child.value, (int, float))
                   and not isinstance(child.value, bool)):
            continue
        for constant in names:
            if constant.isupper() and not constant.startswith("_"):
                found[constant] = node.lineno
    return found


def unclassified(modules) -> dict[str, int]:
    """Constants in those modules that are neither registered nor declared.

    The answer the completeness test asserts is empty, returned as a mapping of
    ``module.NAME`` to its line number so that a failure names the line rather
    than the fact of a failure. Pass the module names `app.py` requires; a
    screen in `sections/` draws and does not hold constants, and passing one
    here would look in the wrong place for them.
    """
    registered = {f"{entry.where[0]}.{entry.where[1]}" for entry in entries()}
    missing: dict[str, int] = {}
    for name in modules:
        for constant, line in scan_module(name).items():
            key = f"{name}.{constant}"
            if key not in registered and key not in NOT_SCIENTIFIC:
                missing[key] = line
    return missing


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
