"""
sections.sequencing - Sanger .ab1 chromatograms in, VHH sequences out.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import io
import zipfile

import matplotlib.pyplot as plt
import streamlit as st

import demo
import sanger as sg
from compute import png_quality, process_ab1
from shell import FIT, Settings, go_to

__version__ = "2026.11.16b"


def _zipped(files) -> bytes:
    """The example reads as a .zip, so they can leave the app as files."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.txt", demo.SANGER_NOTE + "\n")
        for name, data in files:
            archive.writestr(name, data)
    return buffer.getvalue()


def render(settings: Settings) -> None:
    """Draw the screen.

    The router gives every screen the same argument. This one reads
    none of it: see `shell.Settings` for why that is cheaper than a
    second list of the screens.
    """
    st.title("🧪 Sanger sequencing")
    st.markdown(
        "Upload the chromatograms from a screening run. Each read is quality-trimmed, "
        "then the PelB–VHH–tag cassette is located in any strand and reading frame, "
        "and the VHH protein sequence is extracted. Pick the clones you want and send "
        "them to the analysis section.")

    ab1_files = st.file_uploader("Chromatograms (.ab1)", type=["ab1"],
                                 accept_multiple_files=True, key="ab1_uploader")

    # An example, for somebody who has no .ab1 to hand. It is parked under a
    # key of our own rather than pushed into the uploader, which Streamlit
    # does not allow a script to fill, and it is dropped the moment a real
    # file is uploaded: two sets of reads on one screen would be a screen
    # that cannot say which it is describing.
    example = st.columns([2, 5])
    with example[0]:
        if st.button("Load an example run", key="ab1_example",
                     disabled=bool(ab1_files), **FIT,
                     help="Six simulated reads, enough to see every status "
                          "the table reports."):
            st.session_state["ab1_example_files"] = demo.sanger_example_files()
            st.rerun()
    if ab1_files:
        st.session_state.pop("ab1_example_files", None)
    example_files = st.session_state.get("ab1_example_files")
    with example[1]:
        if example_files:
            side = st.columns([1, 1])
            with side[0]:
                if st.button("Clear the example", key="ab1_example_clear"):
                    st.session_state.pop("ab1_example_files", None)
                    st.rerun()
            with side[1]:
                st.download_button(
                    "The example files (.zip)", _zipped(example_files),
                    file_name="vhh_toolbox_example_reads.zip",
                    mime="application/zip", key="ab1_example_zip",
                    help="The same six .ab1 files, to open in a chromatogram "
                         "viewer or to run through your own pipeline.")

    opts = st.columns(4)
    with opts[0]:
        # The list is the bundled maps, now five rather than the two that had
        # hand-written entries, so it is sorted by name and index 0 is no longer
        # the default vector. Point at the default by name instead.
        vector_names = list(sg.VECTORS)
        vector = st.selectbox(
            "Vector", vector_names,
            index=(vector_names.index(sg.DEFAULTS["vector"])
                   if sg.DEFAULTS["vector"] in vector_names else 0),
            help="Determines which tags are expected after the VHH. The layout is "
                 "read from the vector's own map.")
    with opts[1]:
        mott_cutoff = st.slider("Trim cutoff", 0.01, 0.20,
                                sg.DEFAULTS["mott_cutoff"], 0.01,
                                help="Mott algorithm error-probability cutoff. "
                                     "Lower is stricter.")
    with opts[2]:
        min_len = st.number_input("Min trimmed length (nt)", 50, 800,
                                  sg.DEFAULTS["min_trimmed_len"], 10)
    with opts[3]:
        min_phred = st.number_input("Min mean Phred", 5, 50,
                                    sg.DEFAULTS["min_mean_phred"], 1)
    st.caption(f"{vector} layout: " + " → ".join(sg.VECTORS[vector]["layout"]))

    if not ab1_files and not example_files:
        st.info("No chromatogram uploaded yet.")
    else:
        if ab1_files:
            payload = tuple((f.name, f.getvalue()) for f in ab1_files)
        else:
            payload = tuple(example_files)
            st.warning("**Example run, not a measurement.** " + demo.SANGER_NOTE
                       + " The six reads are a clean clone read from both "
                         "ends, a second clone, a clone with one nucleotide "
                         "deleted, a read that stops inside the VHH, and an "
                         "empty well: one of each status the table below "
                         "reports.")
        bar = st.progress(0.0, text="Reading chromatograms…")
        try:
            sanger_df, reads = process_ab1(payload, vector, mott_cutoff,
                                           int(min_len), float(min_phred))
        finally:
            bar.empty()

        n_total = len(sanger_df)
        n_usable = int(sanger_df["usable"].sum()) if "usable" in sanger_df else 0
        n_unique = sanger_df[sanger_df["usable"]]["vhh_aa"].nunique() if n_usable else 0
        cols_m = st.columns(4)
        cols_m[0].metric("Reads", n_total)
        cols_m[1].metric("VHH extracted", n_usable)
        cols_m[2].metric("Unique sequences", n_unique)
        cols_m[3].metric("Failed", n_total - n_usable)

        show_cols = [c for c in ["read", "status", "qc_pass", "trimmed_len",
                                 "mean_phred_trim", "strand", "frame", "vhh_len",
                                 "tags_found", "tags_missing", "premature_stops",
                                 "n_copies", "duplicate_of", "warnings"]
                     if c in sanger_df.columns]
        st.dataframe(
            sanger_df[show_cols], hide_index=True, **FIT,
            column_config={
                "warnings": st.column_config.TextColumn("Warnings", width="large"),
                "mean_phred_trim": st.column_config.NumberColumn("Mean Q",
                                                                 format="%.1f"),
                "trimmed_len": st.column_config.NumberColumn("Trimmed (nt)"),
                "vhh_len": st.column_config.NumberColumn("VHH (aa)"),
                "duplicate_of": st.column_config.TextColumn("Duplicate of"),
            })
        st.caption("Drag a column edge to widen it, or use the fullscreen button at "
                   "the top right of the table.")

        flagged = (sanger_df[sanger_df["warnings"].astype(bool)]
                   if "warnings" in sanger_df.columns else sanger_df.iloc[0:0])
        if len(flagged):
            with st.expander(f"Warnings in full ({len(flagged)} read(s))",
                             expanded=True):
                for _, row in flagged.iterrows():
                    st.markdown(f"**{row['read']}**, *{row['status']}*")
                    for w in str(row["warnings"]).split("; "):
                        if w:
                            st.markdown(f"- {w}")

        if n_total - n_usable:
            st.caption("A failed read is not necessarily a failed clone: "
                       "'no FR4' usually means the read stopped before the end of the "
                       "VHH, and a frameshift means PelB and FR4 sit in different "
                       "reading frames, so an indel lies between them.")

        if reads:
            with st.expander("Quality traces"):
                img_q = png_quality(payload, mott_cutoff)
                if img_q:
                    st.image(img_q, **FIT)
                    st.caption("Blue is the window kept after trimming; the dashed "
                               "line is Phred 20.")
                plt.close("all")

        usable_df = sanger_df[sanger_df["usable"]] if n_usable else sanger_df.iloc[0:0]
        if n_usable:
            st.subheader("Send to the analysis section")
            unique_first = usable_df[usable_df["duplicate_of"] == ""]["read"].tolist()
            chosen = st.multiselect(
                "Sequences to analyse", usable_df["read"].tolist(),
                default=unique_first,
                help="Duplicates are left out by default, one entry per unique "
                     "protein sequence.")
            replace = st.checkbox("Replace what is already loaded there", value=True)
            if st.button("Send selected to analysis", type="primary",
                         disabled=not chosen, **FIT):
                picked = usable_df[usable_df["read"].isin(chosen)]
                fasta = sg.reads_to_fasta(picked, only_usable=False)
                current = st.session_state.get("seq_text", "").strip()
                st.session_state["pending_seq_text"] = (
                    fasta if replace or not current else current + "\n" + fasta)
                go_to("analysis")
                st.rerun()

        st.subheader("Download")
        dl = st.columns(3)
        with dl[0]:
            st.download_button("Sequencing results (CSV)",
                               sanger_df.drop(columns=["vhh_nt"], errors="ignore")
                               .to_csv(index=False),
                               file_name="sanger_results.csv", mime="text/csv", **FIT)
        with dl[1]:
            st.download_button("VHH proteins (FASTA)", sg.reads_to_fasta(sanger_df),
                               file_name="vhh_proteins.fasta", mime="text/plain",
                               disabled=not n_usable, **FIT)
        with dl[2]:
            st.download_button("VHH DNA (FASTA)",
                               sg.reads_to_fasta(sanger_df, column="vhh_nt"),
                               file_name="vhh_dna.fasta", mime="text/plain",
                               disabled=not n_usable, **FIT)

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
