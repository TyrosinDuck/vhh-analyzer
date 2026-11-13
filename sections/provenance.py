"""
sections.provenance - where every constant in the app came from, in one place.

This screen holds no constant and no source text of its own. Both come from
`provenance.py`, which names the module and attribute that owns each value and
reads it live, so a constant changed anywhere in the app changes here with no
edit. See that module for what this can and cannot guarantee.

Why it exists as a screen rather than as a section of the manual: the manual is
prose, and a reader checking whether to trust a number wants the number itself,
as the running app holds it. A manual can say what a threshold is; only the app
can show what it is.
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

import provenance as pv
from shell import FIT, REPO_URL, Settings, usage_prompt

__version__ = "2026.11.13a"

VIEWS = ["📋 By constant", "⚖️ By standing"]


def render(settings: Settings) -> None:
    """Draw the Provenance screen.

    The router gives every screen the same argument. This one reads none of
    it: see `shell.Settings` for why that is cheaper than a second list of the
    screens.
    """
    st.title("🔎 Provenance")
    st.markdown(
        "Every constant the app computes with, what it is worth, and where it "
        "came from. The values on this page are read from the modules that own "
        "them while the page is drawn, so this cannot show a number the app is "
        "not using.")
    st.info(
        "The point of this page is the constants that are **not** measurements. "
        "A threshold chosen by convention and a value verified against a "
        "publication both look like a number on a screen, and the difference "
        "decides how much weight a result can carry.", icon="⚖️")

    entries = pv.entries()
    st.caption(
        f"{len(entries)} constants on the register, and "
        f"{len(pv.NOT_SCIENTIFIC)} more declared to carry no scientific "
        f"claim: guard rails, cache settings, figure sizes, widget defaults "
        f"and unit conversions. The test suite walks the app's modules and "
        f"fails on a public constant that is in neither set, so this page "
        f"covers all of them rather than the ones somebody remembered to add.")

    view = st.radio("View", VIEWS, horizontal=True, key="pv_view",
                    label_visibility="collapsed")
    st.divider()

    if view == VIEWS[0]:
        _table(entries)
    else:
        _by_basis()

    st.divider()
    st.caption(
        "What this page cannot do: check that the sentence still describes the "
        "value. The value is read from the code; the source text is prose "
        "written by a person. If a threshold moves and the sentence under it "
        "does not, this page shows the new number with the old justification. "
        f"The code beside each constant is the authority, and it is public: "
        f"[{REPO_URL}]({REPO_URL}).")
    usage_prompt("Provenance")


# --------------------------------------------------------------------------
def _table(entries) -> None:
    st.subheader("Every constant, as the app holds it")
    frame = pd.DataFrame([{
        "Constant": entry.label,
        "Value": pv.value_text(pv.resolve(entry)),
        "Unit": entry.unit,
        "Standing": entry.basis,
        "Lives in": ".".join(str(step) for step in entry.where[:2]),
        "Reaches": entry.used_for,
    } for entry in entries])
    st.dataframe(frame, **FIT, hide_index=True)
    st.caption(
        "A table of twenty values is reported as its size rather than printed "
        "out: the accessibility maxima and the Chou-Fasman parameters are one "
        "row per residue each, and the Reference screen draws them properly.")

    st.download_button(
        "⬇️ Download the register (CSV)", frame.to_csv(index=False),
        file_name="vhh_toolbox_constants.csv", mime="text/csv", key="pv_csv",
        **FIT)
    st.caption(
        "The CSV carries what this table shows, not the source sentences: "
        "those are below, under the standing each constant has, where they "
        "can be read whole.")


# --------------------------------------------------------------------------
def _by_basis() -> None:
    st.subheader("Grouped by what the number is worth")
    st.markdown(
        "Read from the bottom up if you are deciding how far to trust a "
        "result: the last groups are the ones where a number could reasonably "
        "be somewhere else.")

    grouped = pv.by_basis()
    for basis, entries in grouped.items():
        st.markdown(f"#### {basis} · {len(entries)}")
        st.caption(pv.BASES[basis])
        for entry in entries:
            value = pv.value_text(pv.resolve(entry))
            unit = f" {entry.unit}" if entry.unit else ""
            with st.expander(f"{entry.label}  ·  {value}{unit}"):
                st.markdown(entry.source or
                            "_No source recorded, which is itself a defect: "
                            "please report it._")
                if entry.note:
                    st.caption(entry.note)
                st.caption(f"Reaches: {entry.used_for}")
                st.caption("Lives in `" +
                           ".".join(str(step) for step in entry.where) + "`")
