"""
sections.phage - the phage stocks that go into each round of panning.

Moved out of Batch & QC, where the titre was the only view that needed no
sequence. The screen draws; the arithmetic, and the constants with their
sources, live in `bench.py`, unchanged by the move. The widget keys keep their
``ph_`` prefix, so nothing typed before the move is lost within a session.
"""

from __future__ import annotations

import streamlit as st

import bench as bn
from shell import SECTIONS

__version__ = "2026.11.12a"


def render() -> None:
    """Draw the Phage & panning screen (key "phage")."""
    icon, label, *_ = SECTIONS["phage"]
    st.title(f"{icon} {label}")
    st.markdown(
        "The phage stocks that go into each round of panning. Nothing here "
        "needs a sequence: the purified protein that comes out at the end of "
        "a campaign is the business of Batch & QC.")
    st.divider()
    _titre()
    st.divider()
    st.caption(
        "Arithmetic on readings you supply, with constants that cite their "
        "source and can be edited. Nothing here predicts how many phage "
        "infect, bind or come back from a round.")


def _titre() -> None:
    st.subheader("Titre of a phage stock from its A260")
    st.markdown(
        "An absorbance at 260 nm turns into particles per millilitre, and "
        "that decides how much stock goes into the next round of panning.")

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
