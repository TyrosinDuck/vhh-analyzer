"""
sections.manual - The user manual, rendered from MANUAL.md so there is one
source of truth.

One section at a time rather than seventeen hundred lines in one scroll.
Nobody reads a reference manual from the top: you arrive with a question about
one screen, and what you want is that screen's entry.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st

from compute import manual_sections
from shell import FIT, footer, load_manual

__version__ = "2026.10.29a"

WHOLE = "Everything, in one page"


def render() -> None:
    """Draw the screen."""
    manual = load_manual()
    sections = manual_sections(manual)

    st.title("📖 Manual")

    if not sections:
        # No headings parsed: show the file rather than an empty page.
        st.markdown(manual)
        return

    titles = [s["title"] for s in sections]
    top = st.columns([3, 1])
    with top[0]:
        chosen = st.selectbox(
            "Jump to", titles + [WHOLE], index=0, key="manual_section",
            help="The manual is long on purpose: it documents what every number "
                 "means and how far to trust it. This shows one section at a "
                 "time.")
    with top[1]:
        st.write("")
        st.download_button("Download (Markdown)", manual,
                           file_name="VHH_manual.md", mime="text/markdown",
                           **FIT)

    if chosen == WHOLE:
        st.caption(f"{len(sections)} sections, {len(manual.splitlines())} lines.")
        st.divider()
        st.markdown(manual)
        return

    section = sections[titles.index(chosen)]

    # The subheadings of the chosen section, as a local table of contents. Not
    # links: Streamlit's anchors are generated per rendered page and are not
    # reliable to target from here, and a list of what is in the section is
    # most of the value anyway.
    if section["subs"]:
        with st.expander(f"What is in section {chosen.split('.')[0]} "
                         f"({len(section['subs'])} parts)"):
            st.markdown("\n".join(f"- {sub}" for sub in section["subs"]))

    position = titles.index(chosen)
    nav = st.columns([1, 1, 2])
    if position > 0:
        with nav[0]:
            st.button(f"← {titles[position - 1].split('.')[0]}.",
                      key="manual_prev", **FIT,
                      on_click=_go, args=(titles[position - 1],),
                      help=titles[position - 1])
    if position < len(titles) - 1:
        with nav[1]:
            st.button(f"{titles[position + 1].split('.')[0]}. →",
                      key="manual_next", **FIT,
                      on_click=_go, args=(titles[position + 1],),
                      help=titles[position + 1])

    st.divider()
    st.markdown(f"## {chosen}\n\n{section['body']}")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.


def _go(title: str) -> None:
    """Move the selectbox, on the next rerun, from the prev/next buttons."""
    st.session_state["manual_section"] = title
