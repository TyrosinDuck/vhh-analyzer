"""
sections.home - The landing page: one card per section, and the manual.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st

from shell import (APP_ICON, APP_NAME, AUTHOR_URL, DOI, DOI_URL, FIT,
                   LICENCE, REPO_URL, SECTIONS, current_sequences,
                   footer, go_to)

__version__ = "2026.10.28d"


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

    st.write("")
    with st.expander("📑 Citing this tool"):
        st.markdown(
            f"Archived on Zenodo, which mints a DOI and keeps a copy of the "
            f"code. That matters more than the links: a GitHub repository and a "
            f"Streamlit URL both depend on accounts somebody has to keep alive, "
            f"and a DOI does not.\n\n"
            f"The DOI below is the **concept DOI** — it always resolves to the "
            f"latest release. If you need to say exactly which code produced a "
            f"result, use the version DOI of the release you ran, shown on the "
            f"Zenodo record.")
        st.code(
            f"Vaneyck, J. {APP_NAME} [Computer software]. "
            f"https://doi.org/{DOI}", language=None)
        st.caption("BibTeX")
        st.code(
            "@software{vaneyck_vhh_tool_box,\n"
            "  author  = {Vaneyck, Jonathan},\n"
            f"  title   = {{{APP_NAME}}},\n"
            f"  doi     = {{{DOI}}},\n"
            f"  url     = {{{REPO_URL}}},\n"
            "  license = {MIT}\n"
            "}", language="bibtex")
        st.markdown(
            f"[DOI]({DOI_URL}) · [Source]({REPO_URL}) · "
            f"[ORCID]({AUTHOR_URL}) · {LICENCE} licence\n\n"
            f"The {LICENCE} licence covers this app's own code. One dependency, "
            f"`pyfamsa` (FAMSA, used for the multiple alignment), is "
            f"GPL-3.0-or-later — see `LICENSE` in the repository.")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
