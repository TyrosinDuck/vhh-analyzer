"""
sections.home - The landing page: one card per section, and the manual.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st

from shell import APP_ICON, APP_NAME, FIT, SECTIONS, current_sequences, footer, go_to

__version__ = "2026.10.28a"


def render() -> None:
    """Draw the screen."""
    st.title(f"{APP_ICON} {APP_NAME}")
    st.markdown(
        "A workbench for single-domain antibodies: characterize them, pull them out "
        "of sequencing reads, and plan the cloning. Pick where to start.")

    # Two rows of two rather than one row of four: four cards side by side are
    # too narrow for their descriptions to be readable.
    entries = list(SECTIONS.items())
    for row_start in range(0, len(entries), 2):
        for col, (key, (icon, label, desc)) in zip(
                st.columns(2), entries[row_start:row_start + 2]):
            with col, st.container(border=True):
                st.markdown(f"### {icon} {label}")
                st.markdown(desc)
                st.button(f"Open {label.lower()}", key=f"home_{key}",
                          type="primary", **FIT, on_click=go_to, args=(key,))

    st.write("")
    hcols = st.columns([1, 2])
    with hcols[0]:
        st.button("📖 Read the manual", key="home_manual", **FIT,
                  on_click=go_to, args=("manual",))
    with hcols[1]:
        st.caption("Every screen explained, what each number means, and how far to "
                   "trust it.")

    staged = current_sequences()
    if staged:
        st.info(f"{len(staged)} sequence(s) are currently loaded in the analysis "
                f"input: {', '.join(list(staged)[:6])}"
                f"{'…' if len(staged) > 6 else ''}")

    with st.expander("How the sections fit together"):
        st.markdown("""
**Sanger sequencing → Sequence analysis.** After extracting VHHs from a set of
chromatograms you select the clones worth keeping and send them straight to the
analysis section, which is where every characterization lives. Duplicates are
deselected by default, so one entry per unique protein sequence.

**Sequence analysis → In-silico cloning.** The cloning section reads the same
sequence list, so whatever you are analysing is available as an insert without
retyping it.

**Sequence analysis → Structure.** Everything the analysis section derives from
sequence alone — which cysteines pair, whether a deamidation site matters, where
the aggregation-prone stretches are — becomes a measurement once you upload a
model. The structure section reads the same sequence list, so it can tell you
which of the liabilities it already flagged are buried and which are on the
surface.

**Sequence library → anywhere.** The library holds proteins and DNA; send the
proteins to the analysis and the DNA to gene diagnosis.

Each section also works on its own: you can go straight to cloning with a pasted
DNA insert, or straight to analysis with a FASTA file.
""")
    footer()
