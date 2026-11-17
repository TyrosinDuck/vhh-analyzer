"""
sections.plates - dilution series and plate maps.

Its own section rather than a view of Batch & QC, for the reason the phage
titre moved out: neither of these needs a sequence, and a screen that asks for
one before it can draw is a screen asking a question with no bearing on the
answer. The arithmetic lives in `plates.py`, importable from a notebook.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

import plates as pt
from shell import (AUTHOR, DOI, FIT, SECTIONS, VERSION, Settings,
                   usage_prompt)

__version__ = "2026.11.17b"

V_SERIES = "💧 Dilution series"
V_LAYOUT = "🔲 Plate layout"
VIEWS = [V_SERIES, V_LAYOUT]

# How the two protocols are named on screen. The keys are what `plates` takes.
MODE_LABELS = {
    pt.SERIES: "By transfer, each from the one before",
    pt.INDEPENDENT: "Independently, each from the stock",
}


def render(settings: Settings) -> None:
    """Draw the Plates screen (key "plates").

    Reads none of the settings: a dilution is a volume and a factor, and a
    plate map is a grid. See `shell.Settings` for why the argument is taken
    anyway.
    """
    icon, label, *_ = SECTIONS["plates"]
    st.title(f"{icon} {label}")
    st.markdown(
        "What to pipette for a dilution series, and which well holds what. "
        "Neither needs a sequence. Both are arithmetic on the numbers you "
        "give: a layout says what goes where, not what the readings will "
        "mean, and a series says what to pipette, not what the curve will "
        "look like.")

    view = st.radio("View", VIEWS, horizontal=True,
                    label_visibility="collapsed", key="plates_view")
    st.divider()
    if view == V_LAYOUT:
        _layout()
    else:
        _series()
    st.divider()
    st.caption(
        "Arithmetic and combinatorics on what you supply. The only constant "
        f"here is the {pt.PIPETTE_MINIMUM_UL:.0f} uL this app treats as the "
        "smallest volume worth pipetting, which it shares with the "
        "conjugation planner in Batch & QC.")


# --------------------------------------------------------------------------
def _series() -> None:
    st.subheader("Dilution series")
    st.markdown(
        "The calculation that is easy to get wrong without noticing. In a "
        "series made by transfer, a position that must still hold its final "
        "volume **after** passing some on has to be made up with the final "
        "volume plus the transfer. Forget that and the last positions come "
        "out short, or the last dilution is not the one you meant.")

    mode = st.radio(
        "How the series is made", list(MODE_LABELS),
        format_func=lambda key: MODE_LABELS[key], key="pt_mode",
        help="Two protocols, not two estimates of one number. The app does "
             "not choose between them.")

    top = st.columns(4)
    with top[0]:
        start = st.number_input("Starting concentration", 0.0, 1e12, 10.0,
                                0.1, format="%.4f", key="pt_start")
    with top[1]:
        unit = st.text_input("Unit", value="mg/mL", key="pt_unit",
                             help="Carried through to the table. Nothing is "
                                  "converted: the arithmetic is the same "
                                  "whatever the unit.")
    with top[2]:
        fold = st.number_input("Fold per step", 1.01, 1e6, 2.0, 0.5,
                               key="pt_fold")
    with top[3]:
        points = st.number_input("Positions", 2, 384, 8, 1, key="pt_points")

    volume = st.number_input(
        "Volume each position must still hold (uL)", 0.1, 1e6, 100.0, 10.0,
        key="pt_vol",
        help="What stays in the well or tube once the series is finished, "
             "not what you pipette into it.")

    try:
        out = pt.serial_dilution(start, fold, int(points), volume, mode=mode,
                                 unit=unit or "")
    except ValueError as exc:
        st.error(f"That series cannot be made: {exc}")
        return

    if mode == pt.SERIES:
        cards = st.columns(4)
        with cards[0]:
            st.metric("Transfer each step", f"{out['transfer_ul']:.2f} uL")
        with cards[1]:
            st.metric("Diluent per position",
                      f"{out['diluent_per_position_ul']:.1f} uL")
        with cards[2]:
            st.metric("Stock needed", f"{out['stock_needed_ul']:.1f} uL")
        with cards[3]:
            st.metric("Final dilution", f"1:{out['final_dilution']:g}")
        st.caption(
            f"Every position takes the same {out['diluent_per_position_ul']:.1f} "
            f"uL of diluent and the same {out['transfer_ul']:.2f} uL transfer, "
            f"which is what the carry-over buys you: one volume to set on the "
            f"pipette for the whole series. The last position keeps the "
            f"{out['surplus_ul']:.2f} uL it never passes on.")
    else:
        cards = st.columns(3)
        with cards[0]:
            st.metric("Stock needed in total",
                      f"{out['stock_needed_ul']:.1f} uL")
        with cards[1]:
            st.metric("Diluent in total", f"{out['diluent_total_ul']:.1f} uL")
        with cards[2]:
            st.metric("Final dilution", f"1:{out['final_dilution']:g}")
        st.caption(
            "Each position is made from the stock, so a pipetting error in "
            "one does not travel to the next. The cost is more stock and, as "
            "the dilution deepens, volumes too small to pipette.")

    shown = pd.DataFrame([{
        "Position": row["position"],
        "Dilution": f"1:{row['dilution']:g}",
        f"Concentration ({row['unit']})": row["concentration"],
        "Stock (uL)": row["stock_ul"] or None,
        "Diluent (uL)": row["diluent_ul"] or None,
        "From previous (uL)": row["from_previous_ul"] or None,
        "Pass on (uL)": row["pass_on_ul"] or None,
        "Note": row["note"],
    } for row in out["rows"]])
    st.dataframe(shown, hide_index=True, **FIT)

    for warning in out["warnings"]:
        st.warning(warning)

    st.download_button(
        "⬇️ Download the series (CSV)", shown.to_csv(index=False),
        file_name="dilution_series.csv", mime="text/csv", key="pt_csv", **FIT)

    # Kept so the layout view can offer these as well labels. Labels, not
    # measurements: the batch sheet's rule about retyping is about numbers
    # that go on a record, and this is what is written on a plate map.
    st.session_state["pt_series_labels"] = [
        f"{row['concentration']:.4g} {row['unit']}".strip()
        for row in out["rows"]]
    st.caption(
        f"The **{V_LAYOUT}** view can label a column of wells with these "
        f"{len(out['rows'])} concentrations instead of numbering them D1 to "
        f"D{len(out['rows'])}.")
    usage_prompt("Plates, dilution series")


# --------------------------------------------------------------------------
def _layout() -> None:
    st.subheader("Plate layout")
    st.markdown(
        "Which well holds what, for a plate you are about to set up. The "
        "placement is predictable rather than packed: wells fill in order, "
        "and each replicate of each sample starts at the top of a fresh "
        "column, so the printed map matches the way the series is pipetted. "
        "A layout whose positions move when an unrelated sample is added is "
        "one nobody can check against the bench.")

    left, right = st.columns([3, 2])
    with left:
        names = st.text_area(
            "Samples, one per line", value="NbA\nNbB\nNbC", height=120,
            key="pl_samples")
        controls = st.text_area(
            "Controls, one per line", value="blank\nno antigen", height=80,
            key="pl_controls",
            help="Placed after the samples, with the same number of "
                 "replicates.")
    with right:
        plate = st.radio("Format", sorted(pt.PLATE_FORMATS), horizontal=True,
                         key="pl_format",
                         format_func=lambda n: f"{n} wells")
        dil = st.number_input("Dilution points per sample", 1, 48, 8, 1,
                              key="pl_dil")
        reps = st.number_input("Replicates", 1, 48, 2, 1, key="pl_reps")
        by_column = st.toggle(
            "Series down a column", value=True, key="pl_bycol",
            help="Off puts the series along a row instead.")
        group_per_line = st.toggle(
            "One group per column or row", value=True, key="pl_group",
            help="Off fills every usable well in order, which packs the "
                 "plate but stops the map lining up with the pipetting.")
        skip_edges = st.toggle(
            "Leave the outer wells empty", value=False, key="pl_edges",
            help="Many protocols leave the outer ring out, or fill it with "
                 "buffer, because it evaporates differently. The app does "
                 "not quantify that and takes no view on whether it matters "
                 "for your assay.")

    samples = [line.strip() for line in names.splitlines() if line.strip()]
    control_names = [line.strip() for line in controls.splitlines()
                     if line.strip()]

    # The series from the other view, offered only when it has exactly as
    # many points as this layout, because a label per point or none at all is
    # the only pair that cannot slip by one.
    series = st.session_state.get("pt_series_labels") or []
    labels = None
    if series:
        if len(series) == int(dil):
            if st.checkbox(
                    f"Label the dilutions with the series from "
                    f"**{V_SERIES}** ({series[0]} to {series[-1]})",
                    value=False, key="pl_uselabels",
                    help="Writes the concentration in each well instead of "
                         "D1 to Dn. A label on a map, not a number carried "
                         "onto a record."):
                labels = series
        else:
            st.caption(
                f"The series in **{V_SERIES}** has {len(series)} positions "
                f"and this layout has {int(dil)} dilution points. Match them "
                f"to label the wells with their concentrations.")

    try:
        out = pt.plate_layout(samples, dilution_points=int(dil),
                              replicates=int(reps), controls=control_names,
                              plate=int(plate), skip_edges=skip_edges,
                              by_column=by_column,
                              group_per_line=group_per_line,
                              dilution_labels=labels)
    except ValueError as exc:
        st.error(f"That layout cannot be drawn: {exc}")
        return

    counts = st.columns(4)
    with counts[0]:
        st.metric("Wells used", out["used"])
    with counts[1]:
        st.metric("Usable", out["usable"])
    with counts[2]:
        st.metric("Free", out["free"])
    with counts[3]:
        st.metric("Asked for", out["requested"])

    for warning in out["warnings"]:
        st.warning(warning)

    st.dataframe(out["grid"], **FIT)
    st.caption(
        "Empty cells are wells this layout does not use. The map is what goes "
        "where; it carries no volumes, because what each well takes depends "
        "on the plate and the protocol rather than on the count of its wells.")

    long = pd.DataFrame(out["wells"])
    with st.expander("One row per well"):
        st.dataframe(long, hide_index=True, **FIT)
        st.caption(
            "The form a plate reader or a database ingests, rather than the "
            "grid a person reads.")

    title = st.text_input(
        "Title for the printed sheet", value="Plate layout", key="pl_title",
        help="What the sheet beside the plate is called.")

    three = st.columns(3)
    with three[0]:
        with tempfile.TemporaryDirectory() as folder:
            pdf = Path(folder) / "plate.pdf"
            try:
                pt.build_plate_pdf(
                    pdf, out, title=title or "Plate layout",
                    subtitle=f"{plate}-well", version=VERSION, author=AUTHOR,
                    doi=DOI)
            except Exception as exc:
                st.error(f"The sheet could not be built: {exc}")
            else:
                st.download_button(
                    "📄 Download the sheet (PDF)", pdf.read_bytes(),
                    file_name=f"plate_{plate}.pdf", mime="application/pdf",
                    key="pl_pdf", **FIT)
    with three[1]:
        st.download_button(
            "⬇️ The map (CSV)", out["grid"].to_csv(),
            file_name=f"plate_{plate}_map.csv", mime="text/csv",
            key="pl_grid_csv", **FIT)
    with three[2]:
        st.download_button(
            "⬇️ One row per well (CSV)", long.to_csv(index=False),
            file_name=f"plate_{plate}_wells.csv", mime="text/csv",
            key="pl_long_csv", **FIT)
    st.caption(
        "The PDF is the sheet that goes to the bench next to the plate, laid "
        "out the way the plate is held. The two CSVs are for a machine.")
    usage_prompt("Plates, layout")
