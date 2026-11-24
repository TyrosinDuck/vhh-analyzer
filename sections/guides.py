"""
sections.guides - Step-by-step paths through the app, one per question.

The screen only. The guides themselves are data in `guides.py`; the screen
labels and the analysis families they point at come from `shell`. Nothing is
computed here, so nothing here is cached.

A step can take the reader where it happens. For a step inside the analysis,
the button opens the right family and, where the installed Streamlit can
select a tab (see `shell.tabs_default`), the right tab in it. On a screen
with a view selector of its own, Batch & QC and Structure, it opens the
right view. The first
version opened the family only and left the tab to the reader; read as "Go
there", that landed on the family's first tab, which was the wrong one for
most steps.
"""

from __future__ import annotations

import streamlit as st

import guides as gd
from shell import (ANALYSIS_FAMILIES, FIT, SCREEN_VIEWS, Settings,
                   screen_labels)

__version__ = "2026.11.24c"


def _open(step: gd.Step) -> None:
    """Go to the step's screen, and to its family on the analysis screen.

    Run as a button callback, so before the next run draws anything. The
    family goes through the `pending_` mechanism of `app.py`, the same one the
    sequencing screen uses to fill the analysis input: the radio that owns
    `analysis_family` may already exist in the session, and a widget's state
    is only written before the widget is drawn.
    """
    st.session_state["section"] = step.screen
    if step.family is not None:
        st.session_state["pending_analysis_family"] = gd.family_label(
            step.family, ANALYSIS_FAMILIES)
    if step.tab is not None:
        # Not a widget key, so no `pending_` detour: the analysis screen reads
        # it once, when it draws the family's tabs.
        st.session_state["analysis_open_tab"] = step.tab
    if step.view is not None:
        # A widget key, like the family above, so it goes through the
        # `pending_` detour: the radio that owns it may already exist in the
        # session, and a widget's state is only written before it is drawn.
        state_key, _views = SCREEN_VIEWS[step.screen]
        st.session_state[f"pending_{state_key}"] = gd.view_label(
            step.screen, step.view, SCREEN_VIEWS)


def render(settings: Settings) -> None:
    """Draw the screen.

    The router gives every screen the same argument. This one reads
    none of it: see `shell.Settings` for why that is cheaper than a
    second list of the screens.
    """
    st.title("🧭 Guides")
    st.markdown(
        "The app is organised by tool. A guide starts from a question instead, "
        "and lists the screens that answer it, in order. Nothing is moved: "
        "every step happens where it always does, and **Go there** takes you "
        "to it.")

    keys = list(gd.GUIDES)
    # `keep_` so the question stays chosen while the reader is away on a
    # step's screen: Streamlit discards a widget's state once its screen
    # stops being drawn, and `app.py` holds the keys with this prefix. Coming
    # back to the guide after step 2 and finding the first guide selected is
    # exactly the thing these guides exist to prevent.
    chosen = st.selectbox("Question", keys, key="keep_guide_chosen",
                          format_func=lambda key: gd.GUIDES[key].title)
    guide = gd.GUIDES[chosen]
    labels = screen_labels()

    st.subheader(guide.title)
    st.markdown(f"*{guide.question}*")

    # The screens the path crosses, in the order it first reaches them: the
    # whole route in one line before the detail.
    route = list(dict.fromkeys(labels[step.screen] for step in guide.steps))
    st.caption("Route: " + " → ".join(route))

    for number, step in enumerate(guide.steps, 1):
        with st.container(border=True):
            cols = st.columns([5, 1])
            with cols[0]:
                tag = " *(optional)*" if step.optional else ""
                st.markdown(f"**{number}. {step.action}**{tag}")
                st.caption(gd.location(step, labels, ANALYSIS_FAMILIES, SCREEN_VIEWS))
                if step.why:
                    st.markdown(step.why)
            with cols[1]:
                st.button("Go there", key=f"guide_{chosen}_{number}", **FIT,
                          on_click=_open, args=(step,),
                          help=gd.location(step, labels, ANALYSIS_FAMILIES, SCREEN_VIEWS))

    if guide.note:
        st.info(guide.note, icon="ℹ️")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it.
