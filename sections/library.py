"""
sections.library - The curated sequence library, read-only in the app by design.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import re

import streamlit as st

import library as lib
from compute import load_library_cached
from shell import FIT, footer, go_to

__version__ = "2026.10.28c"


def render() -> None:
    """Draw the screen."""
    st.title("📚 Sequence library")
    st.markdown(
        "The sequences you come back to, kept as a file in the repository so "
        "they are one click away. Pick what you need and send it to the section "
        "that will work on it.")

    table, problems = load_library_cached()

    st.info(
        "**Read-only here, and public.** Streamlit rebuilds this app from the "
        "repository on every reboot and redeploy, so anything written from the "
        "browser would be erased — the library lives in `data/library.csv` "
        "instead, where it is durable and versioned. Edit it on your own "
        "machine, commit and push, and it appears here on the next deploy.\n\n"
        "Because it sits in a public repository it is readable by anyone, "
        "whether or not they open this app. **Nothing unpublished or "
        "confidential belongs in it.**", icon="📚")

    for message in problems:
        st.warning(message, icon="⚠️")

    if not len(table):
        st.info("The library is empty. Add rows to `data/library.csv` — the "
                "columns are `name`, `protein`, `dna`, `tags`, `note`, and at "
                "least one of protein or DNA is required.")
    else:
        filters = st.columns([2, 1, 1])
        with filters[0]:
            search = st.text_input("Search", placeholder="name or note",
                                   key="lib_search")
        with filters[1]:
            tags = lib.all_tags(table)
            chosen_tags = st.multiselect("Tags", tags, key="lib_tags") if tags else []
        with filters[2]:
            needs = st.selectbox("Must carry", ["anything", "a protein", "DNA"],
                                 key="lib_needs")

        view = table
        if search.strip():
            needle = search.strip().lower()
            view = view[view["name"].str.lower().str.contains(needle)
                        | view["note"].str.lower().str.contains(needle)]
        for tag in chosen_tags:
            view = view[view["tags"].str.contains(rf"\b{re.escape(tag)}\b",
                                                  case=False, regex=True)]
        if needs == "a protein":
            view = view[view["has_protein"]]
        elif needs == "DNA":
            view = view[view["has_dna"]]

        st.dataframe(
            view[["name", "has_protein", "has_dna", "protein_length",
                  "dna_length", "tags", "note", "issue"]],
            hide_index=True, **FIT,
            column_config={
                "has_protein": st.column_config.CheckboxColumn("Protein"),
                "has_dna": st.column_config.CheckboxColumn("DNA"),
                "protein_length": st.column_config.NumberColumn("aa"),
                "dna_length": st.column_config.NumberColumn("nt"),
                "note": st.column_config.TextColumn("Note", width="large"),
                "issue": st.column_config.TextColumn("Problem", width="medium")})
        st.caption(f"{len(view)} of {len(table)} entr{'y' if len(table) == 1 else 'ies'} "
                   "shown.")

        st.divider()
        st.subheader("Send to a section")
        picked = st.multiselect("Entries", list(view["name"]), key="lib_picked")

        if not picked:
            st.caption("Pick one or more entries to send them on.")
        else:
            chosen = table[table["name"].isin(picked)]
            n_prot = int(chosen["has_protein"].sum())
            n_dna = int(chosen["has_dna"].sum())
            send = st.columns(2)
            with send[0]:
                if st.button(f"🔬 Send {n_prot} protein(s) to analysis",
                             type="primary", **FIT, disabled=not n_prot):
                    st.session_state["pending_seq_text"] = lib.library_fasta(
                        table, picked, "protein")
                    go_to("analysis")
                    st.rerun()
                if n_prot < len(picked):
                    st.caption(f"{len(picked) - n_prot} of the selected entries "
                               "carry no protein sequence and will be left out.")
            with send[1]:
                if st.button(f"🔧 Send {n_dna} DNA sequence(s) to gene diagnosis",
                             type="primary", **FIT, disabled=not n_dna):
                    st.session_state["pending_gene_dna"] = lib.library_fasta(
                        table, picked, "dna")
                    go_to("genes")
                    st.rerun()
                if n_dna < len(picked):
                    st.caption(f"{len(picked) - n_dna} of the selected entries "
                               "carry no DNA and will be left out.")

            dl = st.columns(2)
            dl[0].download_button("Selected proteins (FASTA)",
                                  lib.library_fasta(table, picked, "protein"),
                                  file_name="library_proteins.fasta",
                                  mime="text/plain", disabled=not n_prot)
            dl[1].download_button("Selected DNA (FASTA)",
                                  lib.library_fasta(table, picked, "dna"),
                                  file_name="library_dna.fasta",
                                  mime="text/plain", disabled=not n_dna)

    with st.expander("How to add or remove entries"):
        st.markdown(
            "The library is the file `data/library.csv` in the repository. Edit "
            "it in a spreadsheet or a text editor on your own machine, commit "
            "it, and push — the app picks it up on its next deploy.\n\n"
            "| Column | |\n|---|---|\n"
            "| `name` | what it is called here; must be unique |\n"
            "| `protein` | the amino-acid sequence, or empty |\n"
            "| `dna` | the coding sequence, or empty |\n"
            "| `tags` | comma-separated, used by the filter above |\n"
            "| `note` | free text |\n\n"
            "At least one of `protein` or `dna` is required. When both are "
            "given, the app checks on load that the DNA translates to the "
            "protein and flags the row if it does not. A malformed row is "
            "shown with its problem rather than dropped, so the library never "
            "loses an entry quietly.")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
