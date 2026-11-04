"""
sections.bench - the arithmetic between a purified protein and the next
experiment.

Every other screen reasons about a sequence. This one is about the tube. The
calculations live in `bench.py`; this is the screen that puts the sequence's own
molar mass, extinction coefficient and lysine count into them, so that none of
the three has to be typed again or looked up wrong.
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from datetime import date
from pathlib import Path

import pandas as pd
import streamlit as st

import bench as bn
import vhh_core as vc
from compute import analyze_one
from shell import (AUTHOR, DOI, FIT, SECTIONS, VERSION, current_sequences,
                   footer)

__version__ = "2026.11.04e"

# Named rather than reached by index. The dispatch below used VIEWS[0] to
# VIEWS[3] with an else, so inserting a view in the middle silently moved every
# screen after it to the wrong function, and the only visible symptom would
# have been the wrong calculator under the right label.
V_QUANTIFY = "🧮 Quantify"
V_MASS = "⚖️ Mass check"
V_CONJUGATE = "🏷️ Conjugate"
V_PHAGE = "🦠 Phage titre"
V_BALANCE = "🧫 Purification"
V_BATCH = "📄 Batch sheet"

VIEWS = [V_QUANTIFY, V_MASS, V_CONJUGATE, V_PHAGE, V_BALANCE, V_BATCH]

# The views that need no sequence: asking for one would be asking a question
# with no bearing on the answer.
NO_SEQUENCE = {V_PHAGE}


def render(ph: float, method: str) -> None:
    """Draw the Batch & QC screen (key "bench")."""
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
        "from the sequence, so none of those has to be typed again. The "
        "exception is the phage titre, which is about the step before all of "
        "that and needs no sequence at all.")

    # The view is chosen before the sequence is asked for, because one of them
    # does not use a sequence at all and asking for one anyway would be
    # asking a question with no bearing on the answer.
    view = st.radio("View", VIEWS, horizontal=True, label_visibility="collapsed",
                    key="bench_view")

    if view in NO_SEQUENCE:
        st.divider()
        _phage()
        _disclaimer()
        return

    name, sequences, props = _sequence_source(ph, method)
    st.divider()

    if view == V_QUANTIFY:
        _quantify(props)
    elif view == V_MASS:
        _mass(name, sequences, props)
    elif view == V_CONJUGATE:
        _conjugate(name, sequences, props, method)
    elif view == V_BALANCE:
        _balance(props)
    else:
        _batch(name, props)

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
    elif not seq:
        st.info(
            "The calculators below still work if you type the molar mass and "
            "extinction coefficient yourself. Giving a sequence fills them in "
            "and removes the most common way of getting this wrong.")

    return name, sequences, props


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


# --------------------------------------------------------------------------
def _phage() -> None:
    st.subheader("Titre of a phage stock from its A260")
    st.markdown(
        "The step before the rest of this screen. An absorbance at 260 nm "
        "turns into particles per millilitre, and that decides how much stock "
        "goes into the next round of panning. This one needs no sequence.")

    read = st.columns(2)
    with read[0]:
        a260 = st.number_input("A260 measured", 0.0, 100.0, 0.459, 0.001,
                               format="%.3f", key="ph_a260")
    with read[1]:
        dilution = st.number_input(
            "Dilution factor", 1.0, 1e6, 100.0, 1.0, key="ph_dil",
            help="How many times the stock was diluted to make the "
                 "measurement. 1 for a neat sample.")

    # Glycerol is asked for as the two volumes rather than as a factor,
    # because the factor is the part that is easy to get backwards. Measuring
    # the phage and then adding glycerol to it means the tube you pipette from
    # is more dilute than what you measured, so the correction is below 1; a
    # field that refused anything below 1 was asking for the number the wrong
    # way round.
    glycerol_modes = ["No glycerol", "From the volumes", "Known factor"]
    mode = st.radio("Glycerol", glycerol_modes, horizontal=True, key="ph_glymode",
                    label_visibility="collapsed")
    glycerol = 1.0
    if mode == glycerol_modes[1]:
        volumes = st.columns(3)
        with volumes[0]:
            v_phage = st.number_input("Phage volume (uL)", 0.1, 1e6, 500.0, 1.0,
                                      key="ph_vphage")
        with volumes[1]:
            v_glycerol = st.number_input("Glycerol added (uL)", 0.0, 1e6, 55.5,
                                         0.1, key="ph_vgly")
        glycerol = v_phage / (v_phage + v_glycerol)
        with volumes[2]:
            st.write("")
            st.info(f"Factor **{glycerol:.4f}**")
        st.caption(
            f"The absorbance was read on the phage before the glycerol went "
            f"in, so the tube you now pipette from is "
            f"{1 / glycerol:.3f} times more dilute than what you measured. If "
            f"you read the absorbance *after* mixing, the glycerol is already "
            f"in the reading and the factor is 1.")
    elif mode == glycerol_modes[2]:
        glycerol = st.number_input(
            "Glycerol factor", 0.01, 100.0, 0.9, 0.01, format="%.4f",
            key="ph_gly",
            help="Below 1 when glycerol was added after the reading, which "
                 "dilutes the tube you pipette from. Above 1 when the stock "
                 "was already in glycerol and you are correcting back to it.")

    with st.expander("The two constants, and why they disagree"):
        st.markdown(f"""
There are two ways to turn an absorbance into a titre, and on the same
reading they differ by more than two orders of magnitude. Both are shown
below rather than one being picked for you.

**The empirical route** multiplies the absorbance by a number measured
directly: a titration in which an A260 of 1 corresponded to
{bn.PHAGE_EMPIRICAL_PER_A260:.1e} particles per millilitre. It is a
correspondence observed for a particular preparation on a particular
instrument, not a physical constant, which is why it is editable here.

**The optical route** is Day's relation for filamentous phage,

    virions/mL = A269 x {bn.PHAGE_OPTICAL_CONSTANT:.0e} / (bases in the genome)

which works the absorbance out from the ssDNA the particle contains. Reading
at 260 nm rather than 269, and without subtracting the A320 scatter, costs
some accuracy but not an order of magnitude.

A gap this size between them is not rounding. It means one of the two is not
describing your stock: most often the empirical factor was established on a
different preparation, or the reading has a large scattering contribution
that the genome calculation has no way to know about. **Use the one you can
defend, and say which in the notebook.** A plaque or colony titration settles
it outright, and is the only one of the three that counts particles that
actually infect.
""")

    fine = st.columns(2)
    with fine[0]:
        per_a260 = st.number_input(
            "Empirical constant (particles/mL per A260)",
            1e6, 1e18, float(bn.PHAGE_EMPIRICAL_PER_A260), format="%.3e",
            key="ph_const")
    with fine[1]:
        bases = st.number_input(
            "Genome length (bases)", 100, 100_000, bn.PHAGE_GENOME_BASES, 1,
            key="ph_bases",
            help="M13mp18 is 7249. A phagemid carrying an insert is a "
                 "different length, which is why this is not fixed.")

    try:
        out = bn.phage_titre(a260, dilution=dilution, glycerol_factor=glycerol,
                             per_a260=per_a260, genome_bases=int(bases))
    except ValueError as exc:
        st.error(str(exc))
        return

    st.divider()
    shown = st.columns(3)
    shown[0].metric("Empirical", f"{out['empirical_per_ml']:.2e} /mL",
                    help="Absorbance x dilution x glycerol x the constant "
                         "above.")
    shown[1].metric("From the genome", f"{out['optical_per_ml']:.2e} /mL",
                    help="Day's relation, from the genome length above.")
    shown[2].metric("Corrected A260", f"{out['a260_corrected']:.1f}",
                    help="What the undiluted stock would read, if it could be "
                         "read undiluted.")

    ratio = out["ratio"]
    if ratio and (ratio > 3 or ratio < 1 / 3):
        st.warning(
            f"The two estimates differ by a factor of "
            f"{max(ratio, 1 / ratio):.0f}. See the note above: pick the one "
            f"you can defend and record which.")

    st.divider()
    st.subheader("How much to take")
    want = st.columns(2)
    with want[0]:
        basis = st.radio("Titre to use", ["Empirical", "From the genome"],
                         horizontal=True, key="ph_basis")
    titre = (out["empirical_per_ml"] if basis == "Empirical"
             else out["optical_per_ml"])
    with want[1]:
        target = st.number_input("Particles wanted", 1e6, 1e18, 1e11,
                                 format="%.3e", key="ph_target",
                                 help="The input for one round on a solid "
                                      "phase. The liquid-phase volume beside "
                                      "it is twice this.")
    if titre <= 0:
        st.info("Give an absorbance above zero.")
        return

    # Two formats, one titre. A solid-phase round takes the input as given; a
    # liquid-phase round is run with twice as many particles, because the
    # antigen is in solution and the capture step loses more of them. Showing
    # one and leaving the reader to double it is how a round gets set up with
    # half the phage it should have had.
    plans = {"Solid phase (adsorbed)": bn.phage_aliquot(titre, target),
             "Liquid phase": bn.phage_aliquot(
                 titre, target * bn.LIQUID_PANNING_FACTOR)}
    shown = st.columns(2)
    for column, (label, plan) in zip(shown, plans.items()):
        with column:
            st.success(f"**{plan['volume_ul']:.1f} uL** of stock\n\n{label}")
            if plan["note"]:
                st.caption(plan["note"])
    st.caption(
        f"Liquid phase is {bn.LIQUID_PANNING_FACTOR:g}x the solid-phase input, "
        f"which is a protocol convention and not something this app derives. "
        f"Change the target above if yours differs.")


# --------------------------------------------------------------------------
def _balance(props: dict) -> None:
    st.subheader("Purification balance")
    st.markdown(
        "What came out of each step, how much of the previous step survived "
        "it, and how much per litre of culture. The page a platform hands to "
        "whoever asked for the protein. Nothing here says what a step *ought* "
        "to recover: it is arithmetic on the volumes and absorbances you type.")

    left, right = st.columns(2)
    with left:
        eps = st.number_input(
            "eps 280 (/M/cm)", 0.0, 1e7,
            float(props.get("ext_coeff_cystines") or 28000.0), 10.0,
            key="pb_eps",
            help="Prefilled from the sequence, with the cystines formed. Used "
                 "for every step quantified by A280.")
        mw = st.number_input("Molar mass (Da)", 0.0, 1e7,
                             float(props.get("mw_da") or 14000.0), 1.0,
                             key="pb_mw")
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


def _batch(name, props: dict) -> None:
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
