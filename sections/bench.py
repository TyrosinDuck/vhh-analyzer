"""
sections.bench - the arithmetic between a purified protein and the next
experiment.

Every other screen reasons about a sequence. This one is about the tube. The
calculations live in `bench.py`; this is the screen that puts the sequence's own
molar mass, extinction coefficient and lysine count into them, so that none of
the three has to be typed again or looked up wrong.
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

import bench as bn
import vhh_core as vc
from compute import analyze_one
from shell import FIT, current_sequences, footer

__version__ = "2026.10.30b"

VIEWS = ["🧮 Quantify", "⚖️ Mass check", "🏷️ Conjugate"]


def render(ph: float, method: str) -> None:
    """Draw the Bench screen."""
    st.title("⚗️ Bench")
    st.markdown(
        "The arithmetic between a purified protein and the next experiment: "
        "how much is in the tube, whether it is the protein you think it is, "
        "and how much labelling reagent to add. Everything here starts from "
        "the molar mass, extinction coefficient and lysine count the analysis "
        "section already derives from the sequence, so none of those has to be "
        "typed again.")

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
            f"eps(280) {props['ext_coeff_cystines']:,} /M/cm with the cystines "
            f"formed, {props['ext_coeff_reduced']:,} reduced · "
            f"{props['n_primary_amines']} primary amines")
    elif not seq:
        st.info(
            "The calculators below still work if you type the molar mass and "
            "extinction coefficient yourself. Giving a sequence fills them in "
            "and removes the most common way of getting this wrong.")

    view = st.radio("View", VIEWS, horizontal=True, label_visibility="collapsed",
                    key="bench_view")
    st.divider()

    if view == VIEWS[0]:
        _quantify(props)
    elif view == VIEWS[1]:
        _mass(name, sequences, props)
    else:
        _conjugate(name, sequences, props, method)

    st.divider()
    st.caption(
        "Arithmetic on numbers you supply, with constants that cite their "
        "source. Nothing here predicts a labelling efficiency, a yield or how "
        "a conjugate will behave: those depend on your protein and your hands.")


# --------------------------------------------------------------------------
def _quantify(props: dict) -> None:
    st.subheader("Concentration from a measured A280")
    st.markdown(
        "Beer-Lambert, with the extinction coefficient computed from the "
        "sequence rather than guessed from a generic 1 mg/mL = 1 AU. For a VHH "
        "that assumption is usually wrong by a third.")

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
        default_eps = float(props.get("ext_coeff_cystines") or 28000)
        cys_state = st.radio(
            "Cysteines", ["Cystines formed", "All reduced"], horizontal=True,
            key="q_cys",
            help="The canonical VHH disulfide is formed in material from the "
                 "periplasm. The two values differ by about 1%, which is the "
                 "choice nobody makes explicitly.")
        if props:
            default_eps = float(props["ext_coeff_cystines"]
                                if cys_state == "Cystines formed"
                                else props["ext_coeff_reduced"])
        eps = st.number_input("eps 280 (/M/cm)", 0.0, 1e7, default_eps, 10.0,
                              key="q_eps")
        mw = st.number_input("Molar mass (Da)", 0.0, 1e7,
                             float(props.get("mw_da") or 14000.0), 1.0,
                             key="q_mw")

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
    cols[0].metric("mg/mL", f"{out['mg_per_ml']:.3f}")
    cols[1].metric("Molar", f"{out['molar_uM']:.1f} uM")
    cols[2].metric("Total protein",
                   f"{out['total_ug']:.1f} ug" if out["total_ug"] else "-")
    cols[3].metric("Total", f"{out['total_nmol']:.2f} nmol"
                   if out["total_nmol"] else "-")

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
def _mass(name, sequences, props: dict) -> None:
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
        st.caption(f"Expected **{props.get('mw_da', 0):,.2f} Da** from the "
                   f"sequence as written.")

    out = bn.mass_check(seq, measured_da=measured or None,
                        tolerance_da=tolerance)

    metric = st.columns(3)
    metric[0].metric("Average", f"{out['average_da']:,.2f} Da")
    metric[1].metric("Monoisotopic", f"{out['monoisotopic_da']:,.2f} Da")
    if out["measured_da"]:
        metric[2].metric("Difference", f"{out['difference_da']:+,.2f} Da")

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

    st.warning(f"**{out['difference_da']:+,.2f} Da** from the sequence as "
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
        st.caption(
            "Each delta is computed from a chemical formula rather than copied "
            "from a table, so you can check any of them by hand. A modification "
            "is only offered when the sequence allows it: gluconoylation needs "
            "an N-terminal His-tag, Met excision needs an N-terminal Met.")


# --------------------------------------------------------------------------
def _conjugate(name, sequences, props: dict, method: str) -> None:
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

    mw = float(props.get("mw_da") or 14000.0)
    amines = props.get("n_primary_amines")
    if not props:
        mw = st.number_input("Protein molar mass (Da)", 1.0, 1e7, 14000.0, 1.0,
                             key="c_pmw")

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
