"""
sections.bench - the arithmetic between a purified protein and the next
experiment.

Every other screen reasons about a sequence. This one is about the tube. The
calculations live in `bench.py`; this is the screen that puts the sequence's own
molar mass, extinction coefficient and lysine count into them, so that none of
the three has to be typed again or looked up wrong.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import tempfile
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

import bench as bn
import vhh_core as vc
from compute import analyze_one, classify_cached
from shell import (AUTHOR, DOI, FIT, SECTIONS, VERSION, Settings,
                   current_sequences, data_download, usage_prompt)

__version__ = "2026.11.23a"

# Named rather than reached by index. The dispatch below used VIEWS[0] to
# VIEWS[3] with an else, so inserting a view in the middle silently moved every
# screen after it to the wrong function, and the only visible symptom would
# have been the wrong calculator under the right label.
V_QUANTIFY = "🧮 Quantify"
V_MASS = "⚖️ Mass check"
V_CONJUGATE = "🏷️ Conjugate"
V_PAIR = "🔗 Couple two proteins"
V_CONVERT = "🔄 Convert"
V_BALANCE = "🧫 Purification"
V_PLAN = "📐 How much culture"
V_BATCH = "📄 Batch sheet"
V_LOTS = "🗂️ Across lots"


@dataclasses.dataclass(frozen=True)
class Tube:
    """What the views of this screen share: one protein, and the settings.

    The same reasoning as `sections.analysis.Run` and `shell.Settings`: every
    view takes this and nothing else, whether or not it reads any of it, so
    the dispatch below is one table rather than a chain of calls with
    different signatures. A chain is what let the earlier version reach a view
    by index.

    Any of the three may be empty. The calculators work on typed numbers, and
    two of the views do not use a sequence at all.
    """

    name: Any
    sequences: Any
    props: Any
    ph: Any
    method: Any


# Every view, the function that draws it, and whether it works on the sequence
# picked at the top of the screen. Listed once: `VIEWS` is the radio's options,
# the dispatch is a lookup in this, and whether to draw the sequence picker is
# read from it rather than decided by an `if` that names views a second time.
# A view added here and nowhere else is reachable and drawn.
#
# The two that need no sequence are the ones that supply their own: the culture
# plan is arithmetic on a target and a yield you measured, and the lot
# comparison takes each sequence from the record it is reading. Drawing the
# picker above them would suggest the view used it. The phage titre, which also
# needs none, is its own section now (sections/phage.py).
def _views() -> dict:
    return {
        V_QUANTIFY: (_quantify, True),
        V_MASS: (_mass, True),
        V_CONJUGATE: (_conjugate, True),
        V_PAIR: (_pair, True),
        V_CONVERT: (_convert, True),
        V_BALANCE: (_balance, True),
        V_PLAN: (_plan, False),
        V_BATCH: (_batch, True),
        V_LOTS: (_lots, False),
    }


def render(settings: Settings) -> None:
    """Draw the Batch & QC screen (key "bench").

    Reads the pH and the annotation method, which decide the net charge and
    the CDR annotation of the sequence the screen pulls properties from. The
    profile window belongs to the analysis profiles and is ignored here.
    """
    ph, method = settings.ph, settings.method
    # From SECTIONS rather than written out, so the title follows a rename.
    icon, label, *_ = SECTIONS["bench"]
    st.title(f"{icon} {label}")
    st.markdown(
        "The arithmetic between a purified protein and the next experiment: "
        "how much is in the tube, whether it is the protein you think it is, "
        "how much labelling reagent to add, how much survived each "
        "purification step, and what goes on the page that "
        "leaves with the tube. They all start from the molar mass, extinction "
        "coefficient and lysine count the analysis section already derives "
        "from the sequence, so none of those has to be typed again.")

    view = st.radio("View", VIEWS, horizontal=True, label_visibility="collapsed",
                    key="bench_view")

    draw, needs_sequence = _views()[view]
    if needs_sequence:
        name, sequences, props = _sequence_source(ph, method)
    else:
        name, sequences, props = None, {}, {}
    st.divider()

    draw(Tube(name=name, sequences=sequences, props=props, ph=ph,
              method=method))

    _disclaimer()


def _disclaimer() -> None:
    st.divider()
    st.caption(
        "Arithmetic on numbers you supply, with constants that cite their "
        "source. Nothing here predicts a labelling efficiency, a yield or how "
        "a conjugate will behave: those depend on your protein and your hands.")


def _sequence_source(ph: float, method: str):
    """Ask for the one sequence the first three views work on.

    Returns ``(name, {name: sequence}, properties)``, any of which may be
    empty: the calculators all work on typed numbers, and refusing to draw
    them until a sequence is supplied would be worse than letting someone
    quantify one tube by hand.
    """
    # Two ways in. The analysis list is the usual one, and keeps the molar
    # mass, extinction coefficient and lysine count consistent with everything
    # else the app says about that sequence. But arriving here with a sequence
    # in the clipboard and no interest in the rest of the app is a real way to
    # use this screen, and sending someone to another section first to quantify
    # one tube would be ceremony.
    loaded = current_sequences()
    sources = ["From the analysis section", "Paste one here"]
    source = st.radio("Sequence", sources, horizontal=True, key="bench_source",
                      index=0 if loaded else 1,
                      label_visibility="collapsed")

    name, seq = None, None
    if source == sources[0]:
        if not loaded:
            st.info(
                "No sequences in the analysis section. Load them there, or "
                "switch to **Paste one here**.")
        else:
            pick = st.columns([2, 3])
            with pick[0]:
                name = st.selectbox("Protein", list(loaded), key="bench_seq")
            seq = loaded[name]
    else:
        pasted = st.text_area(
            "Protein sequence", key="bench_paste", height=110,
            placeholder=">my VHH\nQVQLQESGGGLVQAGGSLRLSCAAS...",
            help="FASTA or a bare sequence. Only the first one is used: this "
                 "screen is about one tube.")
        if pasted.strip():
            try:
                parsed = vc.parse_input(pasted)
                name = next(iter(parsed))
                seq = parsed[name]
            except Exception as exc:
                st.error(f"That could not be read as a protein sequence: {exc}")

    props: dict = {}
    sequences = {name: seq} if (name and seq) else {}
    if seq:
        try:
            props, _annotation, _error = analyze_one(name, seq, ph, method)
            props = props or {}
        except Exception as exc:
            st.warning(f"**{name}** could not be characterised: {exc}")
            props = {}

    if props:
        st.caption(
            f"**{name}** · {props['mw_da']:.1f} Da · "
            f"eps(280) {props['ext_coeff_cystines']} /M/cm with the cystines "
            f"formed, {props['ext_coeff_reduced']} reduced · "
            f"{props['n_primary_amines']} primary amines")
        # A bivalent or biparatopic chain is one molecule in one tube, so every
        # figure on this screen is the whole chain's and that is the right
        # answer here. Said out loud because the analysis screen does the
        # opposite with the same sequence, cutting it into its domains, and two
        # screens disagreeing without explanation is how a number gets used for
        # the wrong thing.
        try:
            n_dom = classify_cached(seq, method).get("n_domains") or 0
        except Exception:
            n_dom = 0
        if n_dom > 1:
            st.caption(
                f"This chain carries **{n_dom} variable domains**. Everything on "
                "this screen, mass, extinction coefficient, concentration, "
                "labelling ratios, is for the whole construct, which is what "
                "goes in the tube. The sequence analysis screen cuts the same "
                "chain into its domains; use that one for the per-domain "
                "readings.")
    elif not seq:
        st.info(
            "The calculators below still work if you type the molar mass and "
            "extinction coefficient yourself. Giving a sequence fills them in "
            "and removes the most common way of getting this wrong.")

    return name, sequences, props


# --------------------------------------------------------------------------
def _quantify(t: Tube) -> None:
    props = t.props

    st.subheader("Concentration from a measured A280")
    st.markdown(
        "Beer-Lambert, with the extinction coefficient computed from the "
        "sequence rather than guessed from a generic 1 mg/mL = 1 AU. For a VHH "
        "that assumption is usually wrong by a third. Computed is not "
        "measured: the result below is derived from your A280 and a predicted "
        "coefficient, and it is labelled as such.")

    left, right = st.columns(2)
    with left:
        a280 = st.number_input("A280 measured", 0.0, 100.0, 0.500, 0.001,
                               format="%.3f", key="q_a280")
        dilution = st.number_input("Dilution factor", 1.0, 1000.0, 1.0, 1.0,
                                   key="q_dil",
                                   help="1 for a neat sample. 10 if you "
                                        "measured 10 uL in 90 uL of buffer.")
        path = st.number_input("Path length (cm)", 0.01, 10.0, 1.0, 0.01,
                               key="q_path",
                               help="1 cm for a cuvette. A NanoDrop pedestal "
                                    "already reports a 1 cm equivalent.")
    with right:
        cys_state = st.radio(
            "Cysteines", ["Cystines formed", "All reduced"], horizontal=True,
            key="q_cys",
            help="The canonical VHH disulfide is formed in material from the "
                 "periplasm. Pairing two cysteines adds exactly 125 /M/cm, so "
                 "the gap is 125 times the number of disulfides: about 0.5% "
                 "for a VHH with the one canonical pair, and 6.6% for the "
                 "most cysteine-rich sequence in this app's own library. It "
                 "is the choice nobody makes explicitly.")
        # No fallback value here. A prefilled eps or molar mass would be
        # arbitrary: it has no source, and the mg/mL it produces looks exactly
        # like one computed from a real sequence. Empty until the user either
        # loads a sequence or types their own.
        default_eps = None
        if props:
            default_eps = float(props["ext_coeff_cystines"]
                                if cys_state == "Cystines formed"
                                else props["ext_coeff_reduced"])
        eps = st.number_input("eps 280 (/M/cm)", 0.0, 1e7, default_eps, 10.0,
                              key="q_eps",
                              placeholder="from a sequence, or type it")
        # Whether the coefficient is a ProtParam prediction or the user's own
        # number decides what the result below may be called.
        _eps_from_sequence = (default_eps is not None
                              and eps == default_eps)
        mw = st.number_input("Molar mass (Da)", 0.0, 1e7,
                             float(props["mw_da"]) if props else None, 1.0,
                             key="q_mw",
                             placeholder="from a sequence, or type it")

    if not eps or not mw:
        st.info(
            "Give an extinction coefficient and a molar mass, by loading a "
            "sequence above or by typing both. Nothing is assumed for you: a "
            "prefilled value would make this read like a result.")
        return

    more = st.columns(2)
    with more[0]:
        volume = st.number_input("Sample volume (uL), optional", 0.0, 1e6,
                                 0.0, 10.0, key="q_vol")
    with more[1]:
        a260 = st.number_input("A260, optional", 0.0, 100.0, 0.0, 0.001,
                               format="%.3f", key="q_a260",
                               help="Gives the 260/280 ratio, the cheap way to "
                                    "notice you are quantifying nucleic acid.")

    try:
        out = bn.concentration_from_a280(
            a280, eps, mw, path_length_cm=path, dilution_factor=dilution,
            a260=a260 or None, volume_ul=volume or None)
    except ValueError as exc:
        st.error(str(exc))
        return

    cols = st.columns(4)
    cols[0].metric("mg/mL (derived)", f"{out['mg_per_ml']:.3f}")
    cols[1].metric("Molar (derived)", f"{out['molar_uM']:.1f} uM")
    cols[2].metric("Total protein",
                   f"{out['total_ug']:.1f} ug" if out["total_ug"] else "-")
    cols[3].metric("Total", f"{out['total_nmol']:.2f} nmol"
                   if out["total_nmol"] else "-")

    # The A280 is measured; the concentration is not. It is Beer-Lambert run
    # on that reading with an extinction coefficient that is predicted from the
    # sequence unless the user typed their own, so the label has to say which.
    st.caption(
        "**Derived, not measured.** The A280 is your measurement. The "
        + ("extinction coefficient above is predicted from the sequence by "
           "ProtParam, so these two values inherit that prediction's error "
           "and are not a measured concentration."
           if _eps_from_sequence else
           "concentration follows from the extinction coefficient you "
           "supplied, and is only as right as that value.")
        + " An amino-acid analysis or a measured A280 of a known standard is "
          "what makes a concentration a measurement.")

    if out["nucleic_acid_warning"]:
        st.warning(out["nucleic_acid_warning"])
    elif out["ratio_260_280"]:
        st.caption(f"A260/A280 = {out['ratio_260_280']:.2f}. A pure protein "
                   f"sits near 0.57.")

    st.divider()
    st.subheader("Dilute it")
    dil = st.columns(3)
    with dil[0]:
        target = st.number_input("Target (mg/mL)", 0.0, 1000.0,
                                 max(0.1, round(out["mg_per_ml"] / 2, 3)),
                                 0.1, key="q_target")
    with dil[1]:
        target_vol = st.number_input("Final volume (uL)", 1.0, 1e6, 100.0, 10.0,
                                     key="q_tvol")
    try:
        plan = bn.dilution_plan(out["mg_per_ml"], target, target_vol)
    except ValueError as exc:
        with dil[2]:
            st.info(str(exc))
        return
    with dil[2]:
        st.success(f"**{plan['stock_ul']:.1f} uL** of sample "
                   f"+ **{plan['buffer_ul']:.1f} uL** of buffer "
                   f"({plan['fold']:g}x)")


# --------------------------------------------------------------------------
def _mass(t: Tube) -> None:
    name, sequences, props = t.name, t.sequences, t.props

    st.subheader("Is it the protein you think it is?")
    if not sequences:
        st.info("Load a sequence in the analysis section to compare a measured "
                "mass against it.")
        return

    seq = vc.clean_sequence(sequences[name])
    st.markdown(
        "A deconvoluted mass from an intact LC-MS run is an **average** mass, "
        "so that is what is compared. The monoisotopic value is given too, for "
        "when you are looking at an isotope envelope instead.")

    cols = st.columns(3)
    with cols[0]:
        measured = st.number_input("Measured mass (Da)", 0.0, 1e7, 0.0, 0.1,
                                   format="%.2f", key="m_measured")
    with cols[1]:
        tolerance = st.number_input("Tolerance (Da)", 0.1, 100.0, 2.0, 0.1,
                                    key="m_tol",
                                    help="How close counts as a match. Intact "
                                         "mass on a Q-TOF is often within 1-2 "
                                         "Da on a 14 kDa protein.")
    with cols[2]:
        st.write("")
        st.caption(f"Expected **{props.get('mw_da', 0):.2f} Da** from the "
                   f"sequence as written.")

    out = bn.mass_check(seq, measured_da=measured or None,
                        tolerance_da=tolerance)

    metric = st.columns(3)
    metric[0].metric("Average", f"{out['average_da']:.2f} Da")
    metric[1].metric("Monoisotopic", f"{out['monoisotopic_da']:.2f} Da")
    if out["measured_da"]:
        metric[2].metric("Difference", f"{out['difference_da']:+.2f} Da")

    if not measured:
        st.info("Enter a measured mass to compare.")
        return

    if out["matches"]:
        st.success(
            f"Within {tolerance:g} Da of the sequence as written. Note that "
            f"this is the mass of the construct you loaded: if the real "
            f"protein carries a tag or has had one cleaved, compare against "
            f"that sequence instead.")
        return

    st.warning(f"**{out['difference_da']:+.2f} Da** from the sequence as "
               f"written.")
    frame = out["candidates"]
    if not len(frame):
        st.info("No candidate modification applies to this sequence.")
        return

    explains = frame[frame["explains"]]
    if len(explains):
        st.markdown("**Differences this size, that this sequence allows:**")
        st.dataframe(
            explains[["modification", "delta_da", "n", "why"]],
            hide_index=True, **FIT,
            column_config={
                "modification": "Modification",
                "delta_da": st.column_config.NumberColumn("Delta (Da)",
                                                          format="%.4f"),
                "n": st.column_config.NumberColumn("How many", format="%d"),
                "why": st.column_config.TextColumn("Why it happens",
                                                   width="large")})
        st.info(
            "Offered, not chosen. More than one combination can land on the "
            "same number, and picking one for you would be presenting "
            "arithmetic as an identification. A peptide map settles it.")
    else:
        st.markdown("**No single modification accounts for it.** The ones this "
                    "sequence allows, for arithmetic of your own:")

    with st.expander("Every modification considered, and why it was or was not "
                     "offered"):
        st.dataframe(frame, hide_index=True, **FIT)
        data_download("Every modification considered", frame,
                      "mass_candidates.csv", key="dl_mass_candidates")
        st.caption(
            "Each delta is computed from a chemical formula rather than copied "
            "from a table, so you can check any of them by hand. A modification "
            "is only offered when the sequence allows it: gluconoylation needs "
            "an N-terminal His-tag, Met excision needs an N-terminal Met.")


# --------------------------------------------------------------------------
def _conjugate(t: Tube) -> None:
    name, sequences, props = t.name, t.sequences, t.props
    method = t.method

    st.subheader("How much labelling reagent to add")
    st.markdown(
        "NHS esters react with the lysine side-chain amines and the "
        "N-terminal alpha-amine. The arithmetic is the one on every vendor "
        "cheat sheet, written out so the units are visible: a millimolar stock "
        "is **one nanomole per microlitre**, which is the only unit identity "
        "worth remembering here.")

    choice = st.selectbox("Reagent", list(bn.REAGENTS) + ["Something else"],
                          key="c_reagent")
    if choice == "Something else":
        custom = st.columns(2)
        with custom[0]:
            reagent_mw = st.number_input("Molar mass (g/mol)", 1.0, 1e6, 500.0,
                                         1.0, key="c_mw")
        with custom[1]:
            vial_mg = st.number_input("Vial contents (mg)", 0.0, 1000.0, 1.0,
                                      0.1, key="c_vial")
        st.caption(
            "Take the molar mass from the vendor's page or the certificate, "
            "not from a catalogue number. Only reagents whose mass was checked "
            "against a primary source are listed above, which is why this entry "
            "exists rather than a longer list of numbers nobody verified.")
        reagent = {"mw": reagent_mw, "vial_mg": vial_mg, "solvents": bn.SOLVENTS}
    else:
        reagent = bn.REAGENTS[choice]
        reagent_mw = reagent["mw"]
        st.caption(f"**{reagent_mw} g/mol**"
                   + (f" · {reagent['formula']}" if reagent.get("formula") else "")
                   + (f" · CAS {reagent['cas']}" if reagent.get("cas") else "")
                   + f" · soluble in {', '.join(reagent['solvents'])}")
        if reagent.get("note"):
            st.info(reagent["note"])
        with st.expander("Where this molar mass comes from"):
            st.markdown(reagent["source"])

    # Reconstitution and the stock you pipette from are two different numbers,
    # and the vendor protocol uses two: a concentrated vial, then a working
    # dilution. Folding them into one makes the app tell you to dissolve a 2 mg
    # vial in 1.7 mL, which is not wrong but is not what anybody does.
    st.markdown("**Reconstitute the vial**")
    stock_cols = st.columns(3)
    with stock_cols[0]:
        vial = st.number_input("Vial contents (mg)", 0.0, 1000.0,
                               float(reagent.get("vial_mg") or 2.0), 0.1,
                               key="c_vialmg")
    with stock_cols[1]:
        neat_mM = st.number_input("Reconstitute to (mM)", 0.001, 1000.0, 20.0,
                                  1.0, key="c_neat",
                                  help="NHS esters hydrolyse in water. Make "
                                       "this concentrated, use it at once, and "
                                       "do not keep it.")
    with stock_cols[2]:
        st.write("")
        if vial > 0:
            prep = bn.stock_preparation(vial, reagent_mw, neat_mM)
            st.success(f"Dissolve in **{prep['solvent_ul']:.0f} uL** "
                       f"({prep['umol_in_vial']:.2f} umol in the vial)")

    work = st.columns(3)
    with work[0]:
        stock_mM = st.number_input("Working stock (mM)", 0.001, 1000.0,
                                   min(2.0, neat_mM), 0.1, key="c_stock",
                                   help="What you actually pipette from. The "
                                        "reaction volume below is computed "
                                        "against this.")
    with work[1]:
        st.write("")
        if stock_mM < neat_mM:
            fold = neat_mM / stock_mM
            st.caption(f"Dilute the reconstituted vial **1:{fold:g}** "
                       f"(e.g. 100 uL + {100 * (fold - 1):.0f} uL solvent).")
        elif stock_mM > neat_mM:
            st.caption("Above the reconstituted concentration; reconstitute "
                       "more concentrated instead.")

    st.markdown("**The reaction**")
    react = st.columns(3)
    with react[0]:
        volume = st.number_input("Protein volume (uL)", 0.1, 1e6, 100.0, 10.0,
                                 key="c_vol")
    with react[1]:
        concentration = st.number_input("Protein (mg/mL)", 0.001, 1000.0, 1.0,
                                        0.1, key="c_conc",
                                        help="The Quantify view works this out "
                                             "from an A280.")
    with react[2]:
        excess = st.number_input("Molar excess of reagent", 0.1, 200.0, 20.0,
                                 0.5, key="c_excess",
                                 help="A starting point, not a recommendation. "
                                      "The excess that gives the labelling you "
                                      "want is something you titrate.")

    amines = props.get("n_primary_amines")
    if props:
        mw = float(props["mw_da"])
    else:
        # No fallback: an arbitrary 14000 Da would silently set the molar
        # excess of reagent, which is the one number this view exists to give.
        mw = st.number_input("Protein molar mass (Da)", 1.0, 1e7, None, 1.0,
                             key="c_pmw",
                             placeholder="from a sequence, or type it")
        if not mw:
            st.info(
                "Give a molar mass, by loading a sequence above or by typing "
                "it. Nothing is assumed for you.")
            return

    try:
        plan = bn.conjugation_plan(volume, concentration, mw, excess, stock_mM,
                                   n_reactive_amines=amines)
    except ValueError as exc:
        st.error(str(exc))
        return

    st.divider()
    result = st.columns(4)
    result[0].metric("Protein", f"{plan['protein_nmol']:.2f} nmol",
                     help=f"{plan['protein_ug']:.1f} ug")
    result[1].metric("Reagent needed", f"{plan['reagent_nmol']:.2f} nmol")
    result[2].metric("Add", f"{plan['reagent_ul']:.2f} uL")
    result[3].metric("Final volume", f"{plan['final_volume_ul']:.1f} uL")

    if plan["dilution"]:
        st.warning(
            f"**{plan['reagent_ul']:.3f} uL cannot be pipetted accurately.** "
            f"Dilute the {stock_mM:g} mM stock 1:{plan['dilution']['factor']} "
            f"to {plan['dilution']['diluted_stock_mM']:g} mM and add "
            f"**{plan['dilution']['volume_ul']:.1f} uL** of that. The molar "
            f"excess is unchanged; the error on it is not.")
    for note in plan["notes"]:
        if "pipetted" not in note:
            st.warning(note)

    st.caption(f"Then desalt: **{plan['desalting']}**.")

    if amines:
        st.divider()
        st.markdown(f"**{amines} primary amines** on this protein: "
                    f"{props.get('n_Lys', 0)} lysine side chains and the "
                    f"N-terminus.")
        try:
            annotation, _ = vc.annotate(sequences[name], method)
            sites = vc.lysine_sites(sequences[name], annotation)
        except Exception:
            sites = None
        if sites is not None and len(sites):
            in_cdr = sites[sites["in_CDR"]]
            if len(in_cdr):
                st.warning(
                    f"**{len(in_cdr)} of them sit in a CDR** "
                    f"({', '.join(in_cdr['site'].astype(str))}). An NHS ester "
                    f"cannot tell them from the others, so random labelling "
                    f"risks the paratope. If that matters, a site-specific "
                    f"route (sortase, an engineered cysteine, a C-terminal "
                    f"tag) is the alternative. Which lysines are actually "
                    f"*accessible* is a structural question: upload a model to "
                    f"the structure section.")
            st.dataframe(sites, hide_index=True, **FIT)
            data_download("Labelling sites", sites, "labelling_sites.csv",
                          key="dl_conjugate_sites")


# --------------------------------------------------------------------------
def _convert(t: Tube) -> None:
    props = t.props

    st.subheader("Convert")
    st.markdown(
        "The two conversions that get done constantly and from a half "
        "remembered formula. Neither has a fitted constant in it.")

    st.markdown("**A concentration, between mass and molarity**")
    st.caption(
        "Within one family it is a decimal shift and no molar mass is read "
        "at all. Crossing between them needs one, because a molar mass is "
        "the only thing that relates a mass to a count of molecules.")

    units = list(bn.MASS_UNITS) + list(bn.MOLAR_UNITS)
    row = st.columns([2, 2, 1, 2])
    with row[0]:
        value = st.number_input("Value", 0.0, 1e12, 1.0, 0.1, format="%.6f",
                                key="cv_value")
    with row[1]:
        unit_in = st.selectbox("From", units, index=units.index("mg/mL"),
                               key="cv_in")
    with row[2]:
        st.write("")
        st.markdown("→")
    with row[3]:
        unit_out = st.selectbox("To", units, index=units.index("uM"),
                                key="cv_out")

    crosses = ((unit_in in bn.MASS_UNITS) != (unit_out in bn.MASS_UNITS))
    mw = None
    if crosses:
        mw = st.number_input(
            "Molar mass (Da)", 0.0, 1e7,
            float(props["mw_da"]) if props else None, 1.0, key="cv_mw",
            placeholder="from a sequence, or type it")
        if not mw:
            st.info(
                "Give a molar mass, by loading a sequence above or by typing "
                "it. Nothing is assumed for you: without one, a mass and a "
                "molarity are not related at all.")
        else:
            try:
                out = bn.convert_concentration(value, unit_in, unit_out, mw)
            except ValueError as exc:
                st.error(str(exc))
            else:
                st.success(f"**{out['value']:.6g} {out['unit']}**")
                st.caption(
                    f"Through a molar mass of {mw:.1f} Da. When that mass "
                    f"came from a sequence it is theoretical, computed rather "
                    f"than measured, and this figure inherits that: it is as "
                    f"good as the mass behind it.")
    else:
        try:
            out = bn.convert_concentration(value, unit_in, unit_out)
        except ValueError as exc:
            st.error(str(exc))
        else:
            st.success(f"**{out['value']:.6g} {out['unit']}**")
            st.caption(
                "Same family, so this is a decimal shift and nothing about "
                "the protein enters into it.")

    st.divider()
    st.markdown("**A centrifuge, between rpm and relative centrifugal field**")
    st.caption(
        "RCF = omega squared times the radius, over the standard "
        "acceleration of free fall, with omega = 2 pi x rpm / 60. That is the "
        "definition of the quantity rather than a formula to look up, which "
        "is why no fitted constant appears: the widely quoted "
        "1.118e-6 x r(mm) x rpm squared is this same expression with the "
        "factors collapsed.")

    spin = st.columns(3)
    with spin[0]:
        radius = st.number_input(
            "Rotor radius (mm)", 0.0, 1000.0, None, 1.0, key="cv_radius",
            placeholder="from your rotor's manual",
            help="The distance from the axis to the point of the tube you "
                 "care about, usually the bottom. It differs between rotors "
                 "and between positions in the same rotor, so the app never "
                 "assumes one.")
    with spin[1]:
        direction = st.radio("Direction", ["rpm → g", "g → rpm"],
                             horizontal=True, key="cv_dir",
                             label_visibility="collapsed")
    with spin[2]:
        if direction == "rpm → g":
            speed = st.number_input("Speed (rpm)", 0.0, 1e6, None, 100.0,
                                    key="cv_rpm", placeholder="what you set")
            field = None
        else:
            field = st.number_input("Field (x g)", 0.0, 1e7, None, 100.0,
                                    key="cv_rcf",
                                    placeholder="what the protocol asks for")
            speed = None

    if not radius or (speed is None and field is None):
        st.info(
            "Give the rotor radius and one of the two. There is no such thing "
            "as a default rotor, and a speed converted with the wrong radius "
            "is wrong in proportion to it.")
        return

    try:
        spun = bn.centrifuge(radius, rpm=speed, rcf=field)
    except ValueError as exc:
        st.error(str(exc))
        return

    cards = st.columns(3)
    with cards[0]:
        st.metric("Speed", f"{spun['rpm']:,.0f} rpm".replace(",", " "))
    with cards[1]:
        st.metric("Field", f"{spun['rcf']:,.0f} x g".replace(",", " "))
    with cards[2]:
        st.metric("Angular velocity", f"{spun['omega_rad_s']:.0f} rad/s")
    st.caption(
        f"At a radius of {spun['radius_mm']:.0f} mm. The field is a ratio to "
        f"the standard acceleration of free fall, a value fixed by convention "
        f"rather than measured where you are, which is what makes it "
        f"comparable between laboratories. The Provenance screen carries the "
        f"declaration it comes from.")

    usage_prompt("Batch & QC, convert")


# --------------------------------------------------------------------------
def _pair(t: Tube) -> None:
    props = t.props

    st.subheader("Couple two proteins")
    st.markdown(
        "For a covalent couple between two proteins rather than between a "
        "protein and a small reagent: a SpyTag/SpyCatcher pair, a sortase "
        "ligation, an enzymatic conjugation. **Conjugate** is the other case, "
        "where the partner comes as a weighed vial at a known molarity; here "
        "it comes as a protein solution in mg/mL, so the arithmetic runs "
        "through two molar masses.")
    st.info(
        "**No yield is predicted.** The ratio says what to add for a "
        "stoichiometry. What fraction actually ligates depends on your "
        "partners and conditions, so the conjugate concentration below is "
        "what would hold if every molecule of the limiting partner reacted: "
        "an upper bound.", icon="🔗")

    left, right = st.columns(2)
    with left:
        st.markdown("**Protein A**, the one in the tube")
        a_vol = st.number_input("Volume of A (uL)", 0.1, 1e6, 100.0, 10.0,
                                key="pr_avol")
        a_conc = st.number_input("Concentration of A (mg/mL)", 0.0001, 1000.0,
                                 1.0, 0.1, format="%.4f", key="pr_aconc")
        a_mw = st.number_input(
            "Molar mass of A (Da)", 0.0, 1e7,
            float(props["mw_da"]) if props else None, 1.0, key="pr_amw",
            placeholder="from a sequence, or type it",
            help="From the sequence picked at the top when there is one, "
                 "which is theoretical: computed from the sequence, not "
                 "measured.")
    with right:
        st.markdown("**Protein B**, the partner you add")
        b_conc = st.number_input("Concentration of B (mg/mL)", 0.0001, 1000.0,
                                 2.0, 0.1, format="%.4f", key="pr_bconc")
        # No fallback here either, and B never comes from the sequence picker:
        # it is the partner, not the protein the screen is pointed at. An
        # invented molar mass would set the volume to add, which is the one
        # number this view exists to give.
        b_mw = st.number_input(
            "Molar mass of B (Da)", 0.0, 1e7, None, 1.0, key="pr_bmw",
            placeholder="from a sequence, or type it",
            help="Typed, because B is usually not the sequence this screen "
                 "is pointed at. Run it through the analysis section to get "
                 "its mass from its sequence.")
        ratio = st.number_input(
            "Molar ratio B:A", 0.01, 100.0, 1.0, 0.1, key="pr_ratio",
            help="1:1 for a stoichiometric couple. Above 1 drives the "
                 "reaction with an excess of B that a one-to-one couple "
                 "cannot consume.")

    # Asked rather than assumed: a condensation releases one water, and which
    # chemistry a given pair uses is a fact about the user's reaction. The
    # water itself is derived from a formula in `buffers`, not typed here.
    loses_water = st.checkbox(
        "The bond is a condensation, so it releases one water", value=False,
        key="pr_water",
        help="True of an isopeptide bond, which is what a SpyTag/SpyCatcher "
             "pair forms. The app does not tabulate which pair does which: "
             "that is your chemistry, and a delta taken from a table nobody "
             "checks is how a wrong mass becomes a wrong identity.")

    if not a_mw or not b_mw:
        st.info(
            "Give both molar masses, by loading a sequence above for A and "
            "typing B. Nothing is assumed for you: a prefilled mass would set "
            "the volume to add and make this read like a result.")
        return

    try:
        out = bn.protein_pair_plan(a_vol, a_conc, a_mw, b_conc, b_mw,
                                   molar_ratio=ratio,
                                   loses_water=loses_water)
    except ValueError as exc:
        st.error(f"That couple cannot be planned: {exc}")
        return

    cards = st.columns(4)
    with cards[0]:
        st.metric("A in the tube", f"{out['a_nmol']:.3f} nmol")
    with cards[1]:
        st.metric("B to add", f"{out['b_nmol']:.3f} nmol")
    with cards[2]:
        st.metric("Volume of B", f"{out['b_volume_ul']:.2f} uL")
    with cards[3]:
        st.metric("Final volume", f"{out['final_volume_ul']:.1f} uL")

    st.caption(
        f"{out['b_ug']:.1f} ug of B into {out['a_ug']:.1f} ug of A. The "
        f"limiting partner is **{out['limiting']}**, so at most "
        f"{out['conjugate_nmol_max']:.3f} nmol of conjugate can exist, which "
        f"in {out['final_volume_ul']:.1f} uL is "
        f"{out['conjugate_mg_per_ml_max']:.3f} mg/mL. That is a ceiling, not "
        f"a prediction.")

    mass = st.columns(2)
    with mass[0]:
        st.metric("Conjugate molar mass", f"{out['conjugate_mw_da']:.1f} Da")
    with mass[1]:
        st.metric("A diluted by", f"{out['dilution_of_a']:.2f}x")
    if out["water_lost_da"]:
        st.caption(
            f"The sum of the two masses minus one water "
            f"({out['water_lost_da']:.3f} Da, derived from its formula with "
            f"the atomic weights in the buffers module). A's mass is "
            f"theoretical when it came from the sequence; B's is whatever you "
            f"typed, which this view cannot classify. Treat the sum as no "
            f"better than the weaker of the two: it is what to compare a "
            f"measured mass against, not a measurement.")
    else:
        st.caption(
            "The plain sum of the two masses, since no condensation was "
            "declared. A's mass is theoretical when it came from the "
            "sequence; B's is whatever you typed.")

    if out["dilution"]:
        d = out["dilution"]
        st.warning(
            f"Dilute B 1:{d['factor']} to {d['diluted_mg_per_ml']:.4f} mg/mL "
            f"and add {d['volume_ul']:.1f} uL of that instead.")
    for note in out["notes"]:
        st.info(note)

    usage_prompt("Batch & QC, couple two proteins")


# --------------------------------------------------------------------------
def _balance(t: Tube) -> None:
    props = t.props

    st.subheader("Purification balance")
    st.markdown(
        "What came out of each step, how much of the previous step survived "
        "it, and how much per litre of culture. The page a platform hands to "
        "whoever asked for the protein. Nothing here says what a step *ought* "
        "to recover: it is arithmetic on the volumes and absorbances you type.")

    left, right = st.columns(2)
    with left:
        # No fallback here either: every row of the balance below is scaled by
        # these two, so an arbitrary default would produce a whole table of
        # plausible-looking yields resting on nothing.
        eps = st.number_input(
            "eps 280 (/M/cm)", 0.0, 1e7,
            float(props["ext_coeff_cystines"]) if props else None, 10.0,
            key="pb_eps",
            placeholder="from a sequence, or type it",
            help="Prefilled from the sequence, with the cystines formed. Used "
                 "for every step quantified by A280.")
        mw = st.number_input("Molar mass (Da)", 0.0, 1e7,
                             float(props["mw_da"]) if props else None, 1.0,
                             key="pb_mw",
                             placeholder="from a sequence, or type it")
    with right:
        culture = st.number_input(
            "Culture volume (mL), optional", 0.0, 1e6, 0.0, 100.0,
            key="pb_culture",
            help="Turns the final pool into mg per litre of culture, which is "
                 "the number a platform is usually asked for.")
        st.caption(
            "Tick **Crude** for any step where the target is not the only "
            "thing absorbing at 280 nm: a lysate, a clarified supernatant, a "
            "flow-through. Those rows keep their total, because it is worth "
            "recording, and stay out of the recovery arithmetic, because the "
            "ratio would not be the target's.")

    if not eps or not mw:
        st.info(
            "Give an extinction coefficient and a molar mass, by loading a "
            "sequence above or by typing both. Every row below is scaled by "
            "them, so nothing is assumed for you.")
        return

    blank = pd.DataFrame({
        "Step": bn.PURIFICATION_STEPS,
        "Volume (mL)": [0.0] * len(bn.PURIFICATION_STEPS),
        "A280": [0.0] * len(bn.PURIFICATION_STEPS),
        "Dilution": [1.0] * len(bn.PURIFICATION_STEPS),
        "mg/mL, other method": [0.0] * len(bn.PURIFICATION_STEPS),
        "Method": [""] * len(bn.PURIFICATION_STEPS),
        # Only the first of the prefilled steps is crude by default, because
        # it is the only one of the four that always is.
        "Crude": [i == 0 for i in range(len(bn.PURIFICATION_STEPS))],
    })
    table = st.data_editor(
        blank, key="pb_steps", num_rows="dynamic", hide_index=True,
        **FIT,
        column_config={
            "Step": st.column_config.TextColumn(width="medium"),
            "Volume (mL)": st.column_config.NumberColumn(
                min_value=0.0, step=0.1, format="%.3f"),
            "A280": st.column_config.NumberColumn(
                min_value=0.0, step=0.01, format="%.3f",
                help="Divided by the extinction coefficient above, so this "
                     "route is derived, not measured."),
            "Dilution": st.column_config.NumberColumn(
                min_value=1.0, step=1.0, format="%.1f"),
            "mg/mL, other method": st.column_config.NumberColumn(
                min_value=0.0, step=0.01, format="%.3f",
                help="A BCA, a Bradford, a gel against a standard. Does not "
                     "go through the sequence, so it is a reading. Takes "
                     "precedence over the A280 on the same row."),
            "Method": st.column_config.TextColumn(
                width="small", help="Named in the table, because a "
                                    "concentration without its method cannot "
                                    "be interpreted later."),
            "Crude": st.column_config.CheckboxColumn(width="small"),
        })

    steps = [{"step": row["Step"],
              "volume_ml": row["Volume (mL)"],
              "a280": row["A280"] or None,
              "dilution": row["Dilution"] or 1.0,
              "mg_per_ml": row["mg/mL, other method"] or None,
              "method": row["Method"],
              "crude": bool(row["Crude"])}
             for _index, row in table.iterrows()]

    try:
        balance = bn.purification_balance(
            steps, ext_coeff=eps or None, mw_da=mw or None,
            culture_volume_ml=culture or None)
    except ValueError as exc:
        st.error(f"That balance could not be computed: {exc}")
        return

    if not any(row["total_mg"] for row in balance["rows"]):
        st.info(
            "Give at least one step a volume and a concentration. A step with "
            "neither is left out rather than counted as zero.")
        return

    shown = pd.DataFrame([{
        "Step": row["step"],
        "Volume (mL)": row["volume_ml"],
        "mg/mL": row["concentration_mg_per_ml"],
        "Total (mg)": row["total_mg"],
        "Step recovery (%)": row["step_recovery_pct"],
        "Cumulative (%)": row["cumulative_recovery_pct"],
        "Basis": row["basis_detail"] or row["basis"] or "",
        "Note": row["comment"],
    } for row in balance["rows"]])
    st.dataframe(shown, **FIT, hide_index=True)
    data_download("The balance as shown", shown, "purification_balance.csv",
                  key="dl_balance_shown",
                  help="The table above. The full record, with every column "
                       "the computation carries, is below.")

    if balance["reference_step"]:
        st.caption(
            f"Cumulative recovery is against **{balance['reference_step']}**, "
            f"the first step counted. The crude rows are excluded, so this is "
            f"recovery from capture onwards and not from the culture.")

    for warning in balance["warnings"]:
        st.warning(warning)

    cards = st.columns(2)
    with cards[0]:
        st.metric("Final pool",
                  f"{balance['final_total_mg']:.3f} mg"
                  if balance["final_total_mg"] else "-")
    with cards[1]:
        st.metric("Per litre of culture",
                  f"{balance['yield_mg_per_litre']:.2f} mg/L"
                  if balance["yield_mg_per_litre"] else "-",
                  help="The final counted pool divided by the culture volume. "
                       "Give a culture volume above to see it.")

    st.caption(
        "A total in mg from an A280 is a reading divided by an extinction "
        "coefficient predicted from the sequence, so it is derived and not "
        "measured, exactly as on the batch sheet. The **Basis** column says "
        "which route each row took.")

    st.download_button(
        "⬇️ Download the balance (CSV)", shown.to_csv(index=False),
        file_name="purification_balance.csv", mime="text/csv", key="pb_csv",
        **FIT)
    st.caption(
        "Nothing is carried from here into the batch sheet. Retype the final "
        "row there: a number you have not retyped is a number you have not "
        "checked.")


# --------------------------------------------------------------------------
def _plan(t: Tube) -> None:
    st.subheader("How much culture, for how much protein")
    st.markdown(
        "The purification balance runs forwards: here is what a preparation "
        "recovered. This runs the other way, which is the question asked "
        "before the preparation: somebody needs a quantity, how much culture "
        "does that take. It is one division, by a yield **you measured**.")
    st.info(
        "**There is no yield model here and there will not be one.** A yield "
        "cannot be predicted from a sequence, and a number that looked like "
        "one would be the most dangerous figure in this app: plausible, "
        "precise, and wrong in a way nobody notices until the purification "
        "has already run. With no previous preparation to divide by, this "
        "view cannot help you.", icon="⚠️")

    st.markdown("**What is needed**")
    ways = ["A mass", "A concentration and a volume"]
    way = st.radio("Target", ways, horizontal=True, key="cp_way",
                   help="The second is the usual way a request arrives: so "
                        "many mL at so many mg/mL.")
    target_mg = target_conc = target_vol = None
    asked = st.columns(3)
    if way == ways[0]:
        with asked[0]:
            target_mg = st.number_input(
                "Protein needed (mg)", 0.0, 1e6, None, 0.1, key="cp_mg",
                placeholder="how many mg")
    else:
        with asked[0]:
            target_conc = st.number_input(
                "Concentration needed (mg/mL)", 0.0, 1000.0, None, 0.1,
                key="cp_conc", placeholder="mg/mL")
        with asked[1]:
            target_vol = st.number_input(
                "Volume needed (mL)", 0.0, 1e5, None, 0.1, key="cp_vol",
                placeholder="mL")
    with asked[2]:
        margin = st.number_input(
            "Margin (%)", 0.0, 500.0, 0.0, 5.0, key="cp_margin",
            help="Added to the target before dividing, for the material that "
                 "goes to assay aliquots and stays in the column. Yours to "
                 "set: nothing here knows what your steps lose.")

    st.markdown("**The yield you measured**")
    # No default, for the same reason the purification balance has none for
    # the extinction coefficient: a plausible mg/L sitting in the box would
    # produce a whole plan resting on nothing, and a plan is acted on.
    labels = {
        "final": "mg of purified protein per litre of culture",
        "expression": "mg per litre expressed, plus the recovery",
    }
    route = st.radio(
        "Route", list(bn.CULTURE_ROUTES), horizontal=True, key="cp_route",
        format_func=lambda key: labels[key],
        help="These must not be mixed. A final yield already has the "
             "purification losses inside it, so applying a recovery to it "
             "would count them twice.")
    st.caption(bn.CULTURE_ROUTES[route].capitalize() + ".")

    measured = st.columns(3)
    with measured[0]:
        observed = st.number_input(
            "Yield (mg/L)", 0.0, 1e5, None, 0.1, key="cp_yield",
            placeholder="from a previous prep",
            help="The **Purification** view computes exactly this figure "
                 "from a preparation you have already run.")
    recovery = None
    with measured[1]:
        if route == "expression":
            recovery = st.number_input(
                "Recovery (%)", 0.0, 100.0, None, 1.0, key="cp_recovery",
                placeholder="what purification kept")
    with measured[2]:
        flask = st.number_input(
            "Working volume per flask (mL), optional", 0.0, 1e5, 0.0, 50.0,
            key="cp_flask",
            help="A property of your shaker, not of the protein. Turns the "
                 "volume into a number of flasks, rounded up.")

    try:
        plan = bn.culture_for_target(
            target_mg=target_mg, target_mg_per_ml=target_conc,
            target_volume_ml=target_vol, route=route,
            yield_mg_per_litre=observed, recovery_pct=recovery,
            margin_pct=margin, flask_working_ml=flask or None)
    except ValueError as exc:
        st.info(f"{exc}.".replace("..", "."))
        return

    st.divider()
    cards = st.columns(3)
    with cards[0]:
        st.metric("Culture needed", f"{plan['culture_litres']:.2f} L",
                  help=f"{plan['culture_ml']:.0f} mL")
    with cards[1]:
        st.metric("Protein to make", f"{plan['required_mg']:.2f} mg",
                  help="The target, plus the margin you set.")
    with cards[2]:
        st.metric("Flasks", plan["flasks"] if plan["flasks"] else "-",
                  help="Rounded up: half a flask is a flask. Give a working "
                       "volume above to see it.")

    st.caption(plan["basis"] + ".")
    for warning in plan["warnings"]:
        st.caption(warning)
    st.caption(
        "Nothing here is a prediction. Change the host, the induction or the "
        "construct and the yield you divided by no longer applies, and the "
        "app cannot tell that it has changed.")


# --------------------------------------------------------------------------
# The widget keys the batch form uses, and the record key each one holds. Kept
# here so that loading a record and drawing the form cannot disagree about
# where a value lives: a reload that filled seven of eight fields would look
# like a sheet somebody had simply left incomplete.
_RECORD_NUMBERS = {
    # widget key: (record key, the widget's own minimum)
    "b_a280": ("a280", 0.0),
    "b_dil": ("dilution", 1.0),
    "b_path": ("path_cm", 0.01),
    "b_vol": ("volume_ul", 0.0),
    "b_conc": ("mg_per_ml", 0.0),
    "b_mass": ("measured_da", 0.0),
}


def _reload_record(props: dict) -> None:
    """Load a saved batch record back into this form.

    The app stores nothing, so this is the only way a sheet printed in March
    can be reprinted in September: the record went to whoever generated it and
    comes back from them. What is restored is the readings and the text, not
    the derived values; those are recomputed below from the readings, which is
    the point. If a later version of the app computes an extinction
    coefficient differently, the reprint will differ from the original, and
    the version each was made with is on the page.
    """
    with st.expander("📂 Reuse a saved record"):
        st.caption(
            "Takes the JSON file saved alongside a previous sheet, and fills "
            "this form with what was typed then. Nothing was kept on the "
            "server; the file is the record.")
        upload = st.file_uploader("Batch record (JSON)", type=["json"],
                                  key="b_record_in",
                                  label_visibility="collapsed")
        if upload is None:
            return

        raw = upload.getvalue()
        digest = hashlib.sha256(raw).hexdigest()
        if st.session_state.get("b_record_digest") == digest:
            # Already applied on an earlier run. Re-applying on every rerun
            # would overwrite whatever the user has edited since.
            st.caption("Loaded. Edit the fields above and download again.")
            return

        try:
            loaded = bn.batch_record_inputs(json.loads(raw.decode("utf-8")))
        except (ValueError, UnicodeDecodeError) as exc:
            st.error(f"That is not a batch record this app can read: {exc}")
            return
        except json.JSONDecodeError as exc:
            st.error(f"That file is not valid JSON: {exc}")
            return

        for key, _label, _default, _help in bn.BATCH_FIELDS:
            st.session_state[f"b_{key}"] = loaded["fields"].get(key, "")

        clamped = []
        for widget_key, (record_key, minimum) in _RECORD_NUMBERS.items():
            value = loaded["measured"].get(record_key)
            if value is None:
                continue
            if value < minimum:
                clamped.append(record_key)
                value = minimum
            st.session_state[widget_key] = float(value)

        st.session_state["b_route"] = ("Determined another way"
                                       if loaded["measured"].get("mg_per_ml")
                                       else "From the A280")
        st.session_state["b_method"] = loaded["measured"].get("method") or "BCA"
        st.session_state["b_cys"] = bool(loaded["cystines"])
        st.session_state["b_ms"] = loaded["measured"].get("measured_da") is not None
        st.session_state["b_record_digest"] = digest
        st.session_state["b_record_from"] = loaded["created_with"]
        st.session_state["b_record_sequence"] = loaded["sequence"]
        if clamped:
            st.session_state["b_record_clamped"] = clamped
        st.rerun()


def _batch(t: Tube) -> None:
    name, props = t.name, t.props

    st.subheader("Batch sheet")
    st.markdown(
        "The page that goes with a tube when it leaves your bench. Three "
        "blocks, kept apart on purpose: what you wrote down, what you "
        "measured, and what the sequence implies. Mixing those is how a "
        "predicted extinction coefficient ends up being quoted as a "
        "measurement a year later.")
    st.info(
        "**Not a certificate of analysis.** There are no specifications here, "
        "no acceptance criteria and no release decision. It is a laboratory "
        "batch record, and the PDF says so at the bottom.", icon="📄")

    if not props:
        st.warning(
            "Pick a sequence at the top. Without one there is nothing to put "
            "in the third block, and a batch sheet with only free text is a "
            "sticky label.")
        return

    _reload_record(props)

    # Said once, after a reload, and only where it changes what the sheet
    # means: which version made the record, whether it was made for this
    # protein, and whether a value had to be moved to fit its input.
    loaded_from = st.session_state.get("b_record_from")
    if loaded_from and loaded_from != VERSION:
        st.caption(
            f"These values came from a record made with {loaded_from}; this "
            f"app is {VERSION}. The readings are what they were. Anything "
            f"derived from them is recomputed here, so it can differ from the "
            f"original sheet.")
    loaded_sequence = st.session_state.get("b_record_sequence")
    if loaded_sequence and loaded_sequence != (props.get("sequence") or ""):
        st.warning(
            "The loaded record was made for a different sequence from the one "
            "selected at the top. The text and the readings are the record's; "
            "everything in the third block is this sequence's. Check which "
            "protein this sheet is for before printing it.")
    clamped = st.session_state.get("b_record_clamped")
    if clamped:
        st.warning(
            "These values in the record are outside what the form accepts and "
            "were raised to its minimum: " + ", ".join(clamped) +
            ". Check them before printing.")

    written = {}
    pairs = st.columns(2)
    for position, (key, label, default, hint) in enumerate(bn.BATCH_FIELDS):
        with pairs[position % 2]:
            if key == "notes":
                continue
            if key == "date":
                written[key] = st.text_input(
                    label, value=date.today().isoformat(), key=f"b_{key}",
                    help=hint or None)
            else:
                written[key] = st.text_input(label, value=default,
                                             key=f"b_{key}",
                                             help=hint or None)
    written["notes"] = st.text_area("Notes", key="b_notes", height=70,
                                    help="Anything the next person needs to "
                                         "know.")

    st.markdown("**Measured**")
    st.caption(
        "Leave a field at zero and it stays off the sheet. Nothing here is "
        "carried over from the other views: a number you have not retyped is "
        "a number you have not checked.")
    routes = ["From the A280", "Determined another way"]
    route = st.radio("Concentration", routes, horizontal=True, key="b_route",
                     help="A concentration in mg/mL from an A280 is not a "
                          "measurement: it is a reading divided by a "
                          "predicted extinction coefficient. The sheet keeps "
                          "the two apart. A BCA or a Bradford does not go "
                          "through the sequence, so it is listed as measured, "
                          "with the method named.")

    got = st.columns(4)
    a280 = mg_per_ml = None
    method = ""
    if route == routes[0]:
        with got[0]:
            a280 = st.number_input("A280", 0.0, 100.0, 0.0, 0.001,
                                   format="%.3f", key="b_a280")
        with got[1]:
            dilution = st.number_input("Dilution factor", 1.0, 1000.0, 1.0,
                                       1.0, key="b_dil")
        with got[2]:
            path_cm = st.number_input("Path length (cm)", 0.01, 10.0, 1.0,
                                      0.01, key="b_path")
    else:
        dilution, path_cm = 1.0, 1.0
        with got[0]:
            mg_per_ml = st.number_input("Concentration (mg/mL)", 0.0, 1000.0,
                                        0.0, 0.01, key="b_conc")
        with got[1]:
            method = st.text_input("How", value="BCA", key="b_method",
                                   help="Named on the sheet, because a "
                                        "concentration without its method is "
                                        "not interpretable later.")
    with got[3]:
        volume = st.number_input("Volume (uL)", 0.0, 1e6, 0.0, 10.0,
                                 key="b_vol")

    cystines = st.checkbox(
        "The canonical disulfide is formed", value=True, key="b_cys",
        help="Which extinction coefficient the concentration uses. True for "
             "material out of the periplasm.")

    # The mass is almost never remeasured, so the field is prefilled with the
    # theoretical value rather than with a zero to retype. It stays behind a
    # checkbox all the same: a sheet that reports a measured mass identical to
    # the sequence, because nobody touched a prefilled box, is exactly the
    # prediction-read-as-measurement this layout exists to prevent.
    theoretical_da = float(props.get("mw_da") or 0.0)
    confirmed = st.checkbox("Mass confirmed by MS", value=False, key="b_ms")
    measured_da = None
    if confirmed:
        measured_da = st.number_input(
            "Measured mass (Da)", 0.0, 1e7, theoretical_da, 0.1,
            format="%.1f", key="b_mass",
            help="Prefilled with the sequence mass, so you only edit the "
                 "digits that differ.")

    sheet = bn.batch_sheet(
        name, props.get("sequence", ""), props,
        measured={"a280": a280 or None, "mg_per_ml": mg_per_ml or None,
                  "method": method, "dilution": dilution, "path_cm": path_cm,
                  "volume_ul": volume or None, "measured_da": measured_da},
        fields=written, cystines=cystines)

    difference = sheet["mass_difference_da"]
    if difference is not None and abs(difference) > 2.0:
        st.warning(
            f"The measured mass is {difference:+.1f} Da from the sequence. "
            f"The **Mass check** view lists the modifications that would "
            f"explain that and that this sequence allows.")

    st.divider()
    blocks = st.columns(3)
    with blocks[0]:
        st.markdown("**Measured**")
        st.caption("Readings. None of these uses the sequence.")
        if sheet["measured"]:
            st.dataframe(pd.DataFrame(sheet["measured"], columns=["", "value"]),
                         **FIT, hide_index=True)
        else:
            st.caption("Nothing entered yet.")
    with blocks[1]:
        st.markdown("**Derived from both**")
        st.caption("A reading and a prediction, combined.")
        if sheet["derived"]:
            st.dataframe(pd.DataFrame(sheet["derived"], columns=["", "value"]),
                         **FIT, hide_index=True)
        else:
            st.caption("Nothing to derive yet.")
    with blocks[2]:
        st.markdown("**From the sequence**")
        st.caption("Theoretical, all of it.")
        st.dataframe(pd.DataFrame(sheet["theoretical"], columns=["", "value"]),
                     **FIT, hide_index=True)

    for note in sheet["derivation_notes"]:
        st.caption(note)

    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "batch.pdf"
        try:
            bn.build_batch_pdf(path, sheet, version=VERSION, author=AUTHOR,
                               doi=DOI)
        except Exception as exc:
            st.error(f"The PDF could not be built: {exc}")
            return
        stamp = (written.get("lot") or name or "batch").replace(" ", "_")
        both = st.columns(2)
        with both[0]:
            st.download_button(
                "📄 Download the batch sheet", path.read_bytes(),
                file_name=f"batch_{stamp}.pdf", mime="application/pdf",
                key="b_pdf", **FIT)
        with both[1]:
            record = bn.batch_record(
                name, props.get("sequence", ""), written,
                {"a280": a280 or None, "mg_per_ml": mg_per_ml or None,
                 "method": method, "dilution": dilution, "path_cm": path_cm,
                 "volume_ul": volume or None, "measured_da": measured_da},
                cystines, sheet=sheet, version=VERSION)
            st.download_button(
                "⬇️ Download the record (JSON)",
                json.dumps(record, indent=1, ensure_ascii=False),
                file_name=f"batch_{stamp}.json", mime="application/json",
                key="b_json", **FIT)

    st.caption(
        "The PDF is for a human and a freezer box. The JSON is the same batch "
        "in a form a database or a script can read, and the file **Reuse a "
        "saved record** above takes back. Keep them together: nothing is "
        "stored here, so the pair of files is the whole record. It is not an "
        "audit trail, and nothing in it can tell whether it has been edited.")

    with st.expander("🏷️ Tube labels for this batch"):
        st.caption(
            "The few lines that fit on a tube, taken from the sheet above: "
            "the lot, the protein, the concentration when there is one, the "
            "buffer, the storage and the date. A label is not a record. What "
            "it cannot carry, including how that concentration was arrived "
            "at, is on the sheet the tube travels with.")
        lines = bn.label_lines(sheet)
        if not lines:
            st.info(
                "Give the batch a lot number above, at least. A label with "
                "nothing on it identifies nothing.")
        else:
            st.code("\n".join(lines), language=None)
            how = st.columns(3)
            with how[0]:
                copies = st.number_input("Labels", 1, 200, 24, 1,
                                         key="b_lab_n")
            with how[1]:
                cols = st.number_input(
                    "Columns", 1, 8, bn.LABEL_SHEET["columns"], 1,
                    key="b_lab_c")
            with how[2]:
                rws = st.number_input("Rows", 1, 20, bn.LABEL_SHEET["rows"],
                                      1, key="b_lab_r")
            # Width and height were parameters of `build_label_sheet` with no
            # widget, while this caption, the manual and the register all said
            # the whole grid was settable. Printing onto real stock is the one
            # case where the millimetres matter more than the counts, so the
            # honest repair is to expose them rather than withdraw the claim.
            size = st.columns(2)
            with size[0]:
                label_w = st.number_input(
                    "Label width (mm)", 10.0, 200.0,
                    float(bn.LABEL_SHEET["width_mm"]), 0.1, key="b_lab_w")
            with size[1]:
                label_h = st.number_input(
                    "Label height (mm)", 10.0, 200.0,
                    float(bn.LABEL_SHEET["height_mm"]), 0.1, key="b_lab_h")
            st.caption(
                f"The grid is {bn.LABEL_SHEET['width_mm']:g} by "
                f"{bn.LABEL_SHEET['height_mm']:g} mm by default, and all "
                f"four values are yours to set. It matches no commercial sheet in "
                f"particular, and says so rather than naming a product it "
                f"might be half a millimetre out on: each label prints its "
                f"own cut outline.")
            with tempfile.TemporaryDirectory() as folder:
                sheet_path = Path(folder) / "labels.pdf"
                try:
                    bn.build_label_sheet(
                        sheet_path, lines, copies=int(copies),
                        columns=int(cols), rows=int(rws),
                        width_mm=float(label_w), height_mm=float(label_h),
                        # Only a derived concentration carries the asterisk,
                        # so only a derived one gets the footnote explaining
                        # it. A measured figure prints bare, with its method
                        # on the batch sheet where there is room for it.
                        footnote=("* derived from an A280, see the batch "
                                  "sheet"
                                  if sheet.get("concentration_basis")
                                  == "derived" else ""))
                except Exception as exc:
                    st.error(f"The labels could not be built: {exc}")
                else:
                    st.download_button(
                        "🏷️ Download the labels (PDF)",
                        sheet_path.read_bytes(),
                        file_name=f"labels_{stamp}.pdf",
                        mime="application/pdf", key="b_labels", **FIT)
    usage_prompt("Batch & QC, batch sheet")


# --------------------------------------------------------------------------
def _lots(t: Tube) -> None:
    st.subheader("Across lots")
    st.markdown(
        "One batch record describes a tube. A year of them describes a "
        "preparation, and the question they answer together, whether the last "
        "three lots came out like the first three, is not visible in any "
        "single sheet. Drop the JSON records here and they line up in one "
        "table.")
    st.caption(
        "Nothing was kept on the server, so these are your files coming back "
        "from wherever you put them. They are read in your browser session "
        "and are not stored here either.")

    uploads = st.file_uploader(
        "Batch records (JSON)", type=["json"], accept_multiple_files=True,
        key="lots_in",
        help="The files the **Batch sheet** view saves beside each PDF.")
    if not uploads:
        st.info(
            "No records yet. Each one is the JSON file the **Batch sheet** "
            "view offers next to the PDF.")
        return

    records: list[tuple[str, dict]] = []
    unreadable: list[tuple[str, str]] = []
    for upload in uploads:
        try:
            records.append((upload.name,
                            json.loads(upload.getvalue().decode("utf-8"))))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            unreadable.append((upload.name, f"not valid JSON ({exc})"))

    # The properties come from the same cached analysis the rest of the app
    # uses, with the settings this screen was given, so a concentration
    # recomputed here is the one the analysis section would report for that
    # sequence. `bench.py` takes this as a parameter because it computes
    # nothing from a sequence itself.
    def properties_for(sequence: str) -> dict:
        props, _annotation, _error = analyze_one("lot", sequence, t.ph,
                                                 t.method)
        return props or {}

    comparison = bn.compare_batch_records(records,
                                          properties_for=properties_for,
                                          version=VERSION)

    for name, reason in unreadable:
        st.error(f"**{name}**: {reason}")
    for bad in comparison["rejected"]:
        st.error(f"**{bad['source']}**: {bad['reason']}")

    rows = comparison["rows"]
    if not rows:
        st.warning("None of those files could be read as a batch record.")
        return

    table = pd.DataFrame([{
        "Lot": row["lot"] or "(none)",
        "Date": row["date"],
        "Protein": row["name"],
        "Host": row["host"],
        "Purification": row["purification"],
        "mg/mL": row["concentration_mg_per_ml"],
        "Basis": row["basis_detail"] or row["basis"] or "",
        "Volume (uL)": row["volume_ul"],
        "Total (mg)": row["total_mg"],
        "Mass minus sequence (Da)": row["mass_difference_da"],
        "Written with": row["created_with"],
        "File": row["source"],
    } for row in rows])
    st.dataframe(table, **FIT, hide_index=True)
    data_download("The lots as shown", table, "across_lots.csv",
                  key="dl_across_lots")

    for warning in comparison["warnings"]:
        st.warning(warning)

    # Drawn only when every row was quantified the same way. A line through a
    # mixture of A280-derived and separately measured concentrations invites
    # exactly the comparison the warning above says not to make, and a chart
    # is more persuasive than the sentence under it.
    bases = {row["basis"] for row in rows if row["basis"]}
    plottable = [row for row in rows if row["total_mg"]]
    if len(bases) == 1 and len(plottable) > 1:
        st.markdown("**Total per lot, in the order given**")
        st.line_chart(
            pd.DataFrame({"Total (mg)": [row["total_mg"]
                                         for row in plottable]},
                         index=[row["lot"] or row["source"]
                                for row in plottable]))
        st.caption(
            "A line between lots, not a trend: there is no fit here, no "
            "control limit and no verdict on whether a lot is acceptable. "
            "Those need a specification, which is a decision somebody makes.")
    elif len(bases) > 1:
        st.caption(
            "No chart, because these lots were not all quantified the same "
            "way. A line through two methods looks like one measurement "
            "drifting.")

    st.download_button(
        "⬇️ Download the comparison (CSV)", table.to_csv(index=False),
        file_name="batch_lots.csv", mime="text/csv", key="lots_csv", **FIT)
    st.caption(
        "Every derived column was recomputed here from the readings in each "
        "record, with this version of the app, rather than read back from "
        "what the record said at the time. A figure that disagrees with an "
        "old printout is a real disagreement, and the **Written with** column "
        "says which version made each record.")
    usage_prompt("Batch & QC, across lots")


# The radio's options, derived from the dispatch table above rather than typed
# beside it. At the bottom of the module because `_views()` names the view
# functions and they have to exist before it is called: a constant written up
# beside the table raised a NameError on import, which is a screen nobody can
# open. The smoke test compares these labels with the radio the app draws.
VIEWS = list(_views())
