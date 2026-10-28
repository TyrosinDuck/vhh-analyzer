"""
sections.manual - The user manual, rendered from MANUAL.md so there is one source of truth.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st

from shell import FIT, footer, load_manual

__version__ = "2026.10.28c"


def render() -> None:
    """Draw the screen."""
    manual = load_manual()
    top = st.columns([3, 1])
    with top[1]:
        st.download_button("Download (Markdown)", manual, file_name="VHH_manual.md",
                           mime="text/markdown", **FIT)
    st.markdown(manual)

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
