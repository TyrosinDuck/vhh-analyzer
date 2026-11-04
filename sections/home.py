"""
sections.home - The landing page: one card per section, grouped, and the manual.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import html

import streamlit as st

from shell import (APP_ICON, APP_NAME, AUTHOR_URL, DOI, DOI_URL, FIT, GROUPS,
                   ISSUES_URL, LICENCE, REPO_URL, SECTIONS, VERSION,
                   current_sequences, embed_html, footer, go_to,
                   group_header_html, report_url, sections_by_group)

__version__ = "2026.11.04e"

# Drawn as inline SVG rather than with a diagramming library: it is six boxes
# and six arrows, it has to be legible on a light and a dark theme, and a
# dependency for that would be a dependency to keep working.
#
# It is rendered through an iframe rather than st.markdown, which was the
# first attempt and does not work: Streamlit's markdown sanitiser keeps the SVG
# elements but strips their positioning attributes, so the shapes vanish and the
# labels reflow as a paragraph of loose words. In a component the markup is
# untouched.
#
# The cost of a component is an iframe, which inherits no colour from the page,
# so the theme is handled inside it: `currentColor` follows the body colour, and
# a prefers-color-scheme rule sets that colour to match a light or dark theme.
FLOW_DIAGRAM = """<svg viewBox="0 0 900 760" xmlns="http://www.w3.org/2000/svg"
     style="width:100%;height:auto;max-width:900px" role="img"
     aria-label="How the seven sections feed each other. You supply the input
     yourself to Sanger sequencing (.ab1 files), Sequence analysis (sequences),
     Structure (a model), Gene diagnosis (coding DNA) and Batch &amp; QC (a
     sequence).
     Sanger sequencing and the sequence library feed Sequence analysis, which
     feeds Structure, Gene diagnosis, In-silico cloning and Batch &amp; QC; the
     library
     feeds Gene diagnosis directly with DNA. Boxes are tinted by group:
     Protein, DNA and Bench.">
  <defs>
    <marker id="vhh-arrow" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor" fill-opacity="0.5"/>
    </marker>
  </defs>

  <!-- Into the analysis, and the library straight to the gene diagnosis -->
  <g fill="none" stroke="currentColor" stroke-opacity="0.4" stroke-width="1.6"
     marker-end="url(#vhh-arrow)" stroke-linejoin="round">
    <path d="M 224 131 H 277 V 224 H 324"/>
    <path d="M 224 303 H 277 V 254 H 324"/>
    <path d="M 124 340 V 680 H 878 V 297 H 866"/>
  </g>

  <!-- Out of the analysis, as one trunk with four stubs rather than four
       lines that would have to cross each other -->
  <g fill="none" stroke="currentColor" stroke-opacity="0.4" stroke-width="1.6"
     stroke-linejoin="round">
    <path d="M 560 238 H 592"/>
    <path d="M 592 131 V 593"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-opacity="0.4" stroke-width="1.6"
     marker-end="url(#vhh-arrow)">
    <path d="M 592 131 H 642"/>
    <path d="M 592 285 H 642"/>
    <path d="M 592 439 H 642"/>
    <path d="M 592 593 H 642"/>
  </g>

  <!-- What you provide yourself -->
  <g fill="none" stroke="currentColor" stroke-opacity="0.38" stroke-width="1.6"
     stroke-dasharray="5 4" marker-end="url(#vhh-arrow)">
    <path d="M 124 74 V 94"/>
    <path d="M 443 168 V 188"/>
    <path d="M 755 74 V 94"/>
    <path d="M 755 228 V 248"/>
    <path d="M 755 536 V 556"/>
  </g>
  <g fill="none" stroke="currentColor" stroke-opacity="0.34" stroke-width="1.3"
     stroke-dasharray="5 4">
    <rect x="29" y="30" width="190" height="42" rx="8"/>
    <rect x="348" y="124" width="190" height="42" rx="8"/>
    <rect x="660" y="30" width="190" height="42" rx="8"/>
    <rect x="660" y="184" width="190" height="42" rx="8"/>
    <rect x="660" y="492" width="190" height="42" rx="8"/>
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
    <text x="755" y="510">you paste</text>
    <text x="755" y="525" font-weight="600">one sequence</text>
  </g>

  <!-- Edge labels -->
  <g font-size="11.5" fill="currentColor" fill-opacity="0.72"
     font-family="system-ui, -apple-system, sans-serif" text-anchor="middle">
    <text x="277" y="119">selected clones</text>
    <text x="277" y="326">proteins</text>
    <text x="621" y="122">sequences</text>
    <text x="621" y="276">DNA</text>
    <text x="621" y="430">insert</text>
    <text x="621" y="584">properties</text>
    <text x="470" y="696">DNA, straight to the diagnosis</text>
  </g>

<!--BOXES-->

  <g fill="currentColor" font-family="system-ui, -apple-system, sans-serif"
     text-anchor="middle">
    <text x="124" y="147" font-size="11" fill-opacity="0.68">VHHs out of the reads</text>

    <text x="124" y="319" font-size="11" fill-opacity="0.68">22 published sequences</text>

    <text x="442" y="248" font-size="11" fill-opacity="0.68">holds the sequence list</text>
    <text x="442" y="265" font-size="11" fill-opacity="0.68">every other section reads</text>

    <text x="753" y="147" font-size="11" fill-opacity="0.68">measured, not inferred</text>

    <text x="753" y="301" font-size="11" fill-opacity="0.68">rare codons, GC stretches</text>

    <text x="753" y="455" font-size="11" fill-opacity="0.68">backbone + insert</text>

    <text x="753" y="609" font-size="11" fill-opacity="0.68">A280, mass, conjugation</text>
  </g>
<!--LEGEND-->
</svg>"""


# Where each section's box sits in FLOW_DIAGRAM: x, y, width, height, corner
# radius. The analysis box is larger and drawn with a heavier outline, because
# it is the hub. The fill colour is not here: it is the colour of the section's
# group in shell.GROUPS, looked up when the diagram is drawn, so that the
# diagram, the home page cards and the sidebar always tint a section alike.
_BOXES = {
    "sequencing": (24, 96, 200, 70, 10),
    "library": (24, 268, 200, 70, 10),
    "analysis": (324, 190, 236, 96, 12),
    "structure": (642, 96, 222, 70, 10),
    "genes": (642, 250, 222, 70, 10),
    "cloning": (642, 404, 222, 70, 10),
    "bench": (642, 558, 222, 70, 10),
}
_HUB = "analysis"


def flow_diagram() -> str:
    """FLOW_DIAGRAM with its section boxes and the group legend filled in.

    Raises KeyError if a section has no box, rather than drawing a diagram
    that silently omits it; a test draws it for that reason.
    """
    boxes = []
    for key in SECTIONS:
        x, y, w, h, r = _BOXES[key]
        colour = GROUPS[SECTIONS[key][3]][0]
        hub = key == _HUB
        boxes.append(
            f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" '
            f'fill="{colour}" fill-opacity="{0.24 if hub else 0.13}" '
            f'stroke="{colour}" stroke-opacity="0.6" '
            f'stroke-width="{2.2 if hub else 1.4}"/>')
    # The legend is written from GROUPS, one swatch and one name per group,
    # spaced evenly and centred, so a renamed or added group shows up here.
    names = list(sections_by_group())
    step, left = 150, 450 - 150 * (len(names) - 1) / 2
    legend = ['  <g font-family="system-ui, -apple-system, sans-serif" '
              'font-size="13" fill="currentColor">']
    for i, name in enumerate(names):
        cx, colour = left + i * step, GROUPS[name][0]
        legend.append(
            f'    <rect x="{cx - 45:.0f}" y="728" width="16" height="16" '
            f'rx="4" fill="{colour}" fill-opacity="0.35" stroke="{colour}"/>'
            f'<text x="{cx - 22:.0f}" y="741" font-weight="600">{name}</text>')
    legend.append("  </g>")
    # Box titles come from SECTIONS, so a renamed section is renamed here too.
    # Escaped, because a label may hold an "&" ("Batch & QC").
    for key in SECTIONS:
        x, y, w, _h, _r = _BOXES[key]
        icon, label, *_ = SECTIONS[key]
        hub = key == _HUB
        boxes.append(
            f'  <text x="{x + w // 2}" y="{y + (37 if hub else 31)}" '
            f'text-anchor="middle" fill="currentColor" '
            f'font-family="system-ui, -apple-system, sans-serif" '
            f'font-size="{15.5 if hub else 14}" '
            f'font-weight="{700 if hub else 600}">'
            f'{html.escape(icon + " " + label)}</text>')
    return (FLOW_DIAGRAM.replace("<!--BOXES-->", "\n".join(boxes))
            .replace("<!--LEGEND-->", "\n".join(legend)))


# White on every group colour, by choice: black was tried on the green and the
# orange and rejected for its look. GROUPS is darkened to keep this readable.
BUTTON_TEXT = "#ffffff"


def _open_button_css() -> str:
    """CSS that paints each card's "Open" button in the colour of its group.

    Streamlit gives every keyed widget the CSS class ``st-key-<key>``, and the
    card buttons are keyed ``home_<section>``, so each rule targets exactly
    one button and nothing else on the page. Fill colours come from GROUPS;
    the text is white on all of them, which shell.GROUPS darkens its colours
    to keep readable (at least 4.5:1, checked by a test). Hover and focus
    darken the fill a little rather than switching to the theme colour, which
    is what Streamlit's own primary button style would otherwise do.

    Relies on the ``st-key-`` class, which Streamlit documents for
    containers and applies to widgets too in the version pinned here; if a
    later version drops it, the buttons fall back to the theme colour, they
    do not break. The smoke test checks the rules are emitted, not that the
    browser applies them; that was checked by eye, in a headless browser.
    """
    rules = []
    for key, (*_rest, group) in SECTIONS.items():
        fill = GROUPS[group][0]
        text = BUTTON_TEXT
        sel = f".st-key-home_{key} button"
        rules.append(
            f"{sel}{{background-color:{fill};border-color:{fill};color:{text}}}"
            f"{sel} p{{color:{text}}}"
            f"{sel}:hover,{sel}:focus,{sel}:active{{background-color:{fill};"
            f"border-color:{fill};color:{text};filter:brightness(0.9)}}")
    return "<style>" + "".join(rules) + "</style>"


def _card(key: str, height: str = "content") -> None:
    """One section card: title, description, and the button that opens it.

    `height` is passed to the card's container: "stretch" when the card has
    to share out the height of a group drawn beside a taller one.
    """
    icon, label, desc, _group = SECTIONS[key]
    with st.container(border=True, height=height):
        st.markdown(f"#### {icon} {label}")
        st.markdown(desc)
        # The label as written, not lowercased: lowercasing made "Open sanger
        # sequencing" and would make "Open batch & qc".
        st.button(f"Open {label}", key=f"home_{key}",
                  type="primary", **FIT, on_click=go_to, args=(key,))


def render() -> None:
    """Draw the screen."""
    st.title(f"{APP_ICON} {APP_NAME}")
    st.markdown(
        "A workbench for single-domain antibodies: characterize them, pull them out "
        "of sequencing reads, and plan the cloning. Pick where to start.")

    # Cards grouped by the molecule they work on (see shell.GROUPS), under a
    # group title larger than the card titles, which is the hierarchy. Groups
    # are paired two per row, each group a column with its cards stacked in
    # it; a group left without a partner, Bench today, gets the full width.
    #
    # Stacked per group rather than laid out row by row across both groups,
    # although row by row would line the n-th cards up on a wide screen: on a
    # narrow one Streamlit stacks columns, and row by row then interleaves
    # the groups (Structure, Gene diagnosis, Sequence library...), which
    # undoes the grouping on exactly the screens with the least room to
    # compensate. Seen on a 420 px viewport before this was changed.
    #
    # Two groups side by side are made the same height: each group is a
    # container with height="stretch", which fills the column, and every card
    # in it stretches too, so the slack of the shorter group is shared between
    # its cards instead of leaving a gap under its last one. The bottoms line
    # up; the tops of the n-th cards only line up when the cards above them
    # happen to be the same height, which per-group stacking cannot promise.
    # On a narrow screen the columns stack and "stretch" falls back to the
    # content height, so nothing changes there. Checked on a prototype and on
    # the app, at 1400 and 420 px.
    st.markdown(_open_button_css(), unsafe_allow_html=True)
    groups = list(sections_by_group().items())
    for pair_start in range(0, len(groups), 2):
        pair = groups[pair_start:pair_start + 2]
        height = "stretch" if len(pair) > 1 else "content"
        for col, (name, keys) in zip(st.columns(len(pair)), pair):
            with col, st.container(height=height):
                st.markdown(group_header_html(name, "2.1rem"),
                            unsafe_allow_html=True)
                for key in keys:
                    _card(key, height)
        st.write("")

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
        embed_html(
            "<style>body{margin:0;background:transparent;color:#31333F;"
            "font-family:system-ui,-apple-system,sans-serif}"
            "@media (prefers-color-scheme: dark){body{color:#fafafa}}</style>"
            + flow_diagram(), height=760)
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
