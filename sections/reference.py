"""
sections.reference - the two tables you look up rather than compute.

This screen holds no calculation of its own. Everything on it is read out of
the modules the rest of the app uses, which is the only reason it is worth
having: a reference table that is typed in separately can disagree with the
app printed next to it, and these cannot.

Structures are committed SVG files rather than drawn at run time; see
`tools/aminoacids/README.md` for how they are made and what is checked.
"""

from __future__ import annotations

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import reference as rf
from shell import FIT, go_to

__version__ = "2026.10.30a"

VIEWS = ["🧬 Genetic code", "🔤 Amino acids"]

# Tall enough for the drawing plus its caption, in the iframe components.html
# puts it in. The SVGs have no intrinsic height, so this is what decides it.
CARD_HEIGHT = 250


def render() -> None:
    """Draw the Reference screen."""
    st.title("📐 Reference")
    st.markdown(
        "The two tables you look up rather than work out. Neither is typed in: "
        "the codon table is Biopython's, the usage figures are the measured "
        "*E. coli* table the gene diagnosis screen scores CAI against, and the "
        "residue masses, hydropathy and charges are the same constants the "
        "analysis screen builds its profiles from. So nothing here can "
        "contradict what the rest of the app tells you.")

    view = st.radio("View", VIEWS, horizontal=True, key="ref_view",
                    label_visibility="collapsed")
    st.divider()

    if view == VIEWS[0]:
        _genetic_code()
    else:
        _aminoacids()


# --------------------------------------------------------------------------
def _genetic_code() -> None:
    st.subheader("The standard genetic code")

    wheel = rf.codon_wheel()
    st.markdown(
        "Rows are the first two bases, columns the third. Read across a row "
        "and the amino acid usually does not change: that is the wobble "
        "position, and it is why a synonymous swap is available at all.")
    # Tall enough for all sixteen rows: the point of this layout is reading
    # down a column and across a row, and half of it behind a scrollbar is
    # worse than no block at all.
    st.dataframe(wheel, use_container_width=True, height=16 * 35 + 40)

    st.divider()
    st.subheader("Codon by codon, with how *E. coli* uses it")
    st.markdown(
        "`share` is the codon's share of its own amino acid's codons. `w` is "
        "that share divided by the best synonym's, which is the quantity CAI "
        "is the geometric mean of, so a codon with `w` near 1 is the one the "
        "host prefers. `rare_in_ecoli` marks the eight codons the gene "
        "diagnosis screen flags, with the tRNA that is limiting for each.")

    table = rf.genetic_code()
    picked = st.multiselect(
        "Filter by amino acid", sorted(set(table["name"])), default=[],
        key="ref_code_filter",
        help="Leave empty for all 64.")
    if picked:
        table = table[table["name"].isin(picked)]
    only_rare = st.checkbox("Only the rare ones", key="ref_code_rare")
    if only_rare:
        table = table[table["rare_in_ecoli"]]

    st.dataframe(table, use_container_width=True, hide_index=True,
                 column_config={
                     "share": st.column_config.ProgressColumn(
                         "share", min_value=0.0, max_value=1.0,
                         format="%.2f"),
                     "w": st.column_config.NumberColumn("w", format="%.2f"),
                 })
    st.caption(
        "Usage from `data/ecoli_codon_usage.csv`; stop codons have no share "
        "because they do not compete with a synonym for an amino acid.")

    st.button("Open the gene diagnosis screen", key="ref_to_genes", **FIT,
              on_click=go_to, args=("genes",))


# --------------------------------------------------------------------------
def _aminoacids() -> None:
    table = rf.aminoacids()
    ph = rf.structure_ph()

    st.subheader("The twenty amino acids")
    if ph is not None:
        st.markdown(
            f"Drawn as they exist **at pH {ph}**: the zwitterion, with "
            f"aspartate and glutamate deprotonated and lysine and arginine "
            f"protonated. Histidine is drawn neutral, because its side chain "
            f"pKa is near 6 and at {ph} the neutral imidazole dominates; it is "
            f"the one residue whose picture would be different at a pH you "
            f"might plausibly work at. A drawing of a charged species without "
            f"its pH is a drawing of nothing in particular, which is why the "
            f"number is on the screen and not in a footnote.")
    else:
        st.warning(
            "The structure drawings are missing: `data/aminoacids/` was not "
            "copied with the rest of the app. The table below still works.")

    families = ["All"] + sorted(set(table["family"]))
    family = st.radio("Family", families, horizontal=True, key="ref_aa_family")
    shown = table if family == "All" else table[table["family"] == family]

    if ph is not None:
        columns = st.columns(5)
        for position, row in enumerate(shown.itertuples()):
            svg = rf.structure_svg(row.aa)
            with columns[position % 5]:
                if svg:
                    # components.html rather than st.markdown: Streamlit's
                    # markdown sanitiser keeps SVG elements but strips the
                    # attributes that position them, so an inline drawing
                    # comes out as an empty box.
                    components.html(
                        f"<div style='font:12px system-ui;text-align:center;"
                        f"color:inherit'>"
                        f"<div style='height:200px'>{svg}</div>"
                        f"<b>{row.name}</b> &middot; {row.code} ({row.aa})"
                        f"</div>"
                        f"<style>:root{{color-scheme:light dark}}"
                        f"body{{margin:0;color:#262730}}"
                        f"@media (prefers-color-scheme:dark)"
                        f"{{body{{color:#fafafa}}}}</style>",
                        height=CARD_HEIGHT)
                else:
                    st.caption(f"{row.name} ({row.aa}): no drawing")
        st.caption(
            "Stereo descriptors are the CIP ones: every amino acid here is L, "
            "which is S at the alpha carbon except in cysteine, where the "
            "sulfur outranks the carboxyl and the same spatial arrangement is "
            "called R. Threonine is 2S,3R and isoleucine 2S,3S.")

    st.divider()
    st.subheader("Properties")
    st.markdown(
        "`residue_da` is what the residue adds to a chain, which is the number "
        "the app's molecular weights are built from; `free_da` is the free "
        "amino acid. `hydropathy` is Kyte and Doolittle, the scale behind the "
        "hydrophobicity profile. `charge_profile` is the approximation used in "
        "the charge profile, where histidine counts as +0.1 rather than being "
        "treated as fully protonated or fully neutral.")
    st.dataframe(shown, use_container_width=True, hide_index=True)
    st.caption(
        "Side-chain pKa values are Lehninger's, for the free amino acid. In a "
        "folded protein the local environment moves them, sometimes by "
        "several units, so they say what to expect and not what is the case.")
