"""
sections.home - The landing page: one card per section, and the manual.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

from shell import (APP_ICON, APP_NAME, AUTHOR_URL, DOI, DOI_URL, FIT,
                   ISSUES_URL, LICENCE, REPO_URL, SECTIONS, VERSION,
                   current_sequences, footer, go_to, report_url)

__version__ = "2026.10.29b"

# Drawn as inline SVG rather than with a diagramming library: it is six boxes
# and six arrows, it has to be legible on a light and a dark theme, and a
# dependency for that would be a dependency to keep working.
#
# It is rendered through components.html rather than st.markdown, which was the
# first attempt and does not work: Streamlit's markdown sanitiser keeps the SVG
# elements but strips their positioning attributes, so the shapes vanish and the
# labels reflow as a paragraph of loose words. In a component the markup is
# untouched.
#
# The cost of a component is an iframe, which inherits no colour from the page,
# so the theme is handled inside it: `currentColor` follows the body colour, and
# a prefers-color-scheme rule sets that colour to match a light or dark theme.
FLOW_DIAGRAM = """<svg viewBox="0 0 900 560" xmlns="http://www.w3.org/2000/svg"
     style="width:100%;height:auto;max-width:900px" role="img"
     aria-label="How the six sections feed each other. You supply the input
     yourself to Sanger sequencing (.ab1 files), Sequence analysis (sequences),
     Structure (a model) and Gene diagnosis (coding DNA). Sanger sequencing and
     the sequence library also feed Sequence analysis, which feeds Structure,
     Gene diagnosis and In-silico cloning; the library feeds Gene diagnosis
     directly with DNA.">
  <defs>
    <marker id="vhh-arrow" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor" fill-opacity="0.5"/>
    </marker>
  </defs>

  <g fill="none" stroke="currentColor" stroke-opacity="0.4" stroke-width="1.6"
     marker-end="url(#vhh-arrow)" stroke-linejoin="round">
    <path d="M 224 131 H 277 V 224 H 324"/>
    <path d="M 224 303 H 277 V 254 H 324"/>
    <path d="M 560 214 H 605 V 131 H 644"/>
    <path d="M 560 238 H 605 V 285 H 644"/>
    <path d="M 560 262 H 605 V 433 H 644"/>
    <path d="M 124 340 V 505 H 878 V 297 H 866"/>
  </g>

  <g fill="none" stroke="currentColor" stroke-opacity="0.38" stroke-width="1.6"
     stroke-dasharray="5 4" marker-end="url(#vhh-arrow)">
    <path d="M 124 74 V 94"/>
    <path d="M 443 168 V 188"/>
    <path d="M 755 74 V 94"/>
    <path d="M 755 228 V 248"/>
  </g>

  <g fill="none" stroke="currentColor" stroke-opacity="0.34" stroke-width="1.3"
     stroke-dasharray="5 4">
    <rect x="29" y="30" width="190" height="42" rx="8"/>
    <rect x="348" y="124" width="190" height="42" rx="8"/>
    <rect x="660" y="30" width="190" height="42" rx="8"/>
    <rect x="660" y="184" width="190" height="42" rx="8"/>
  </g>

  <g font-size="11" fill="currentColor" fill-opacity="0.72" text-anchor="middle"
     font-family="system-ui, -apple-system, sans-serif">
    <text x="124" y="48">you upload</text>
    <text x="124" y="63" font-weight="600">.ab1 chromatograms</text>
    <text x="443" y="142">you paste or upload</text>
    <text x="443" y="157" font-weight="600">sequences, FASTA</text>
    <text x="755" y="48">you upload</text>
    <text x="755" y="63" font-weight="600">a PDB or AlphaFold model</text>
    <text x="755" y="202">you paste</text>
    <text x="755" y="217" font-weight="600">coding DNA</text>
  </g>

  <g font-size="11.5" fill="currentColor" fill-opacity="0.72"
     font-family="system-ui, -apple-system, sans-serif" text-anchor="middle">
    <text x="277" y="119">selected clones</text>
    <text x="277" y="326">proteins</text>
    <text x="603" y="122">sequences</text>
    <text x="603" y="276">DNA</text>
    <text x="603" y="424">insert</text>
    <text x="470" y="521">DNA, straight to the diagnosis</text>
  </g>

  <g stroke="#4c72b0" stroke-opacity="0.55" stroke-width="1.4">
    <rect x="24" y="96" width="200" height="70" rx="10"
          fill="#4c72b0" fill-opacity="0.10"/>
    <rect x="24" y="268" width="200" height="70" rx="10"
          fill="#4c72b0" fill-opacity="0.10"/>
    <rect x="324" y="190" width="236" height="96" rx="12"
          fill="#4c72b0" fill-opacity="0.22" stroke-width="2.2"/>
    <rect x="644" y="96" width="220" height="70" rx="10"
          fill="#4c72b0" fill-opacity="0.10"/>
    <rect x="644" y="250" width="220" height="70" rx="10"
          fill="#4c72b0" fill-opacity="0.10"/>
    <rect x="644" y="398" width="220" height="70" rx="10"
          fill="#4c72b0" fill-opacity="0.10"/>
  </g>

  <g fill="currentColor" font-family="system-ui, -apple-system, sans-serif"
     text-anchor="middle">
    <text x="124" y="127" font-size="14" font-weight="600">🧪 Sanger sequencing</text>
    <text x="124" y="147" font-size="11" fill-opacity="0.68">VHHs out of the reads</text>

    <text x="124" y="299" font-size="14" font-weight="600">📚 Sequence library</text>
    <text x="124" y="319" font-size="11" fill-opacity="0.68">22 published sequences</text>

    <text x="442" y="227" font-size="15.5" font-weight="700">🔬 Sequence analysis</text>
    <text x="442" y="248" font-size="11" fill-opacity="0.68">holds the sequence list</text>
    <text x="442" y="265" font-size="11" fill-opacity="0.68">every other section reads</text>

    <text x="754" y="127" font-size="14" font-weight="600">🧊 Structure</text>
    <text x="754" y="147" font-size="11" fill-opacity="0.68">measured, not inferred</text>

    <text x="754" y="281" font-size="14" font-weight="600">🔧 Gene diagnosis</text>
    <text x="754" y="301" font-size="11" fill-opacity="0.68">rare codons, GC stretches</text>

    <text x="754" y="429" font-size="14" font-weight="600">🧬 In-silico cloning</text>
    <text x="754" y="449" font-size="11" fill-opacity="0.68">backbone + insert</text>
  </g>
</svg>"""


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
        components.html(
            "<style>body{margin:0;background:transparent;color:#31333F;"
            "font-family:system-ui,-apple-system,sans-serif}"
            "@media (prefers-color-scheme: dark){body{color:#fafafa}}</style>"
            + FLOW_DIAGRAM, height=560, scrolling=False)
        st.markdown("""
Three things the diagram cannot carry:

- **Duplicates are deselected by default** when the sequencing section hands
  clones to the analysis, so you get one entry per unique protein sequence.
- **The structure section needs the sequence list** to map a model onto the
  full-length construct. A model is routinely missing termini and loops, and
  its numbering may be sequential, IMGT, Kabat or the depositor's own.
- **Each section also works on its own.** You can go straight to cloning with a
  pasted DNA insert, or straight to analysis with a FASTA file.
""")

    st.write("")
    with st.expander("🐛 Found a problem, or something that looks wrong?"):
        st.markdown(
            "**A wrong number is the most useful thing you can report.** Most of "
            "what this app computes looks equally plausible whether it is right "
            "or not, which is exactly what makes a quiet error dangerous. If a "
            "value disagrees with ANARCI, with a published figure, with a gel or "
            "with what you measured at the bench, that is worth more than a "
            "feature request.\n\n"
            "⚠️ **The issue tracker is public.** Reproduce the problem with a "
            "published sequence if you can, the library has 22 of them. Nothing "
            "unpublished or under an NDA should go into an issue. To analyse "
            "confidential sequences at all, run the app on your own machine "
            "rather than here; the manual explains how.")
        st.link_button("Report a problem on GitHub", report_url(),
                       type="primary", **FIT)
        st.caption(
            "The link opens a new issue with the version already filled in, "
            f"this one is **{VERSION}**. Without it there is no way to tell "
            "whether a defect still exists.")
        st.markdown(
            f"Questions are welcome in the same place. For what a number means "
            f"and how far to trust it, the manual is usually faster. "
            f"[All issues]({ISSUES_URL}) · "
            f"[How to contribute]({REPO_URL}/blob/main/CONTRIBUTING.md)")

    with st.expander("📑 Citing this tool"):
        st.markdown(
            f"Archived on Zenodo, which mints a DOI and keeps a copy of the "
            f"code. That matters more than the links: a GitHub repository and a "
            f"Streamlit URL both depend on accounts somebody has to keep alive, "
            f"and a DOI does not.\n\n"
            f"The DOI below is the **concept DOI**, it always resolves to the "
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
            f"GPL-3.0-or-later, see `LICENSE` in the repository.")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
