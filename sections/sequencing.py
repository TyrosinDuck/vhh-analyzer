"""
sections.sequencing - Sanger .ab1 chromatograms in, VHH sequences out.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import hashlib
import io
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import demo
import sanger as sg
from compute import png_quality, process_ab1, quality_table_for
from shell import FIT, Settings, data_download, go_to

__version__ = "2026.11.23a"


def _zipped(files) -> bytes:
    """The example reads as a .zip, so they can leave the app as files."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.txt", demo.SANGER_NOTE + "\n")
        for name, data in files:
            archive.writestr(name, data)
    return buffer.getvalue()


def _name_reads(stems: list[str]) -> list[str]:
    """The table where the reads get the user's names. Returns the names used.

    What the user typed is the editor's own state, and nothing here writes it
    back into the frame the editor was built from: a data_editor whose data
    changes under it on every rerun loses the edit in progress. The frame is
    rebuilt only on three occasions, each of which bumps the editor's key: a
    table of names is loaded, the names are reset, or the set of files changes.
    On each, the names last typed are carried into the new frame, so loading a
    sheet for half the reads keeps what was typed for the other half, and
    adding a file to the upload keeps the names already given.
    """
    state = st.session_state
    base = state.setdefault("ab1_names", {})
    last = state.get("ab1_names_last", {})
    files_key = hashlib.sha1("\x00".join(stems).encode()).hexdigest()[:10]
    if state.get("ab1_names_files") != files_key:
        base.update(last)
        for key in [k for k, v in base.items() if not v]:
            del base[key]
        state["ab1_names_files"] = files_key
        state.pop("ab1_names_msg", None)
    rev = state.setdefault("ab1_names_rev", 0)

    renamed = sum(1 for s in stems
                  if (last[s] if s in last else base.get(s, "")).strip())
    label = "✏️ Name the reads" + (f" ({renamed} renamed)" if renamed else "")
    with st.expander(label, expanded=False):
        st.caption(
            "Each read is named after its file, which is the provider's "
            "reference. Type the name you want next to it, or paste a column "
            "from a spreadsheet into *Your name*; an empty cell keeps the "
            "reference. The names are used everywhere below: the table, the "
            "quality traces, the FASTA files and the analysis section. The "
            "reference stays in a column of its own, so a renamed read can be "
            "traced back to its chromatogram.")
        edited = st.data_editor(
            pd.DataFrame({"Reference name": stems,
                          "Your name": [base.get(s, "") for s in stems]}),
            key=f"ab1_names_editor_{files_key}_{rev}", hide_index=True,
            num_rows="fixed", disabled=["Reference name"], **FIT,
            height=min(38 + 35 * len(stems), 400),
            column_config={
                "Reference name": st.column_config.TextColumn(width="medium"),
                "Your name": st.column_config.TextColumn(
                    width="medium",
                    help="Spaces become underscores, because the analysis "
                         "section reads a FASTA name up to its first space."),
            })
        requested = ["" if not isinstance(v, str) else v
                     for v in edited["Your name"].tolist()]
        state["ab1_names_last"] = dict(zip(stems, requested))
        resolved = sg.assign_names(stems, requested)

        changed = resolved[resolved["note"].astype(bool)]
        if len(changed):
            st.warning("Some names were adjusted so that every read has a "
                       "name of its own:\n\n" + "\n".join(
                           f"- **{r.file}**: {r.note}"
                           for r in changed.itertuples()))

        side = st.columns([3, 2, 2])
        with side[0]:
            sheet = st.file_uploader(
                "Or load a table of names (.csv, .tsv, .txt)",
                type=["csv", "tsv", "txt"], key="ab1_names_file",
                help="Two columns: the reference name, with or without .ab1, "
                     "then your name. A header row is optional. An Excel "
                     "sheet can be saved as CSV, or its two columns pasted "
                     "into the table above.")
        if sheet is not None:
            raw = sheet.getvalue()
            digest = hashlib.sha1(raw).hexdigest()
            # Applied once per file, not on every rerun: re-applying it would
            # overwrite whatever is typed into the table after loading it.
            if state.get("ab1_names_sheet") != digest:
                state["ab1_names_sheet"] = digest
                try:
                    mapping, unmatched = sg.read_name_map(raw, stems)
                except ValueError as exc:
                    state["ab1_names_msg"] = ("error", str(exc))
                else:
                    state["ab1_names"] = {**{k: v for k, v in
                                             state["ab1_names_last"].items() if v},
                                          **mapping}
                    state["ab1_names_rev"] = rev + 1
                    msg = f"{len(mapping)} name(s) applied from {sheet.name}."
                    if unmatched:
                        msg += (f" {len(unmatched)} reference(s) in the table "
                                "match no loaded read and were not applied: "
                                + ", ".join(unmatched[:10])
                                + (", …" if len(unmatched) > 10 else "") + ".")
                    state["ab1_names_msg"] = ("info", msg)
                st.rerun()
        kind, msg = state.get("ab1_names_msg", (None, None))
        if msg:
            (st.error if kind == "error" else st.info)(msg)
        with side[1]:
            data_download("The names (CSV)",
                          resolved[["file", "name"]].rename(
                              columns={"file": "reference", "name": "name"}),
                          "read_names.csv", key="dl_read_names",
                          help="The names in use, to keep with the run or to "
                               "load again next time.")
        with side[2]:
            if st.button("Reset to the reference names", key="ab1_names_reset",
                         disabled=not renamed, **FIT):
                state["ab1_names"] = {}
                state["ab1_names_last"] = {}
                state.pop("ab1_names_msg", None)
                state["ab1_names_rev"] = rev + 1
                st.rerun()
    return resolved["name"].tolist()


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
    # The stretch after the last annotated tag is read from the vector's
    # sequence, so say so: on this screen everything else in the layout comes
    # from an annotation on the map, and the two are not equally well attested.
    trailer = sg.VECTORS[vector].get("trailer") or {}
    if trailer.get("aa"):
        note = (f"The {len(trailer['aa'])} residues after the last tag "
                f"({trailer['aa']}, {3 * len(trailer['aa'])} nt) are translated "
                "from this vector's sequence, not read from an annotation on its "
                "map")
        if trailer.get("cys"):
            note += (", and they end on a cysteine immediately before the stop. "
                     "A read has to reach that far to show the clone carries it; "
                     "one that stops earlier is reported as ending inside the "
                     "stretch rather than as missing it")
        st.caption(note + ".")

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
        names = tuple(_name_reads([Path(name).stem for name, _ in payload]))
        bar = st.progress(0.0, text="Reading chromatograms…")
        try:
            sanger_df, reads = process_ab1(payload, vector, mott_cutoff,
                                           int(min_len), float(min_phred),
                                           names)
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

        # The reference column only when a read was renamed: otherwise it
        # repeats the first column word for word.
        renamed = ("file" in sanger_df.columns
                   and bool((sanger_df["file"] != sanger_df["read"]).any()))
        show_cols = [c for c in ["read"] + (["file"] if renamed else []) + ["status", "qc_pass", "trimmed_len",
                                 "mean_phred_trim", "strand", "frame", "vhh_len",
                                 "tags_found", "tags_missing", "tags_partial",
                                 "premature_stops",
                                 "n_copies", "duplicate_of", "warnings"]
                     if c in sanger_df.columns]
        data_download("The results as shown", sanger_df[show_cols],
                      "sanger_results_shown.csv", key="dl_sanger_shown",
                      help="The columns above. The full table, with every "
                           "column the run produced, is further down.")
        st.dataframe(
            sanger_df[show_cols], hide_index=True, **FIT,
            column_config={
                "warnings": st.column_config.TextColumn("Warnings", width="large"),
                "mean_phred_trim": st.column_config.NumberColumn("Mean Q",
                                                                 format="%.1f"),
                "trimmed_len": st.column_config.NumberColumn("Trimmed (nt)"),
                "vhh_len": st.column_config.NumberColumn("VHH (aa)"),
                "duplicate_of": st.column_config.TextColumn("Duplicate of"),
                "file": st.column_config.TextColumn("Reference"),
            })
        st.caption("Drag a column edge to widen it, or use the fullscreen button at "
                   "the top right of the table.")

        flagged = (sanger_df[sanger_df["warnings"].astype(bool)]
                   if "warnings" in sanger_df.columns else sanger_df.iloc[0:0])
        if len(flagged):
            with st.expander(f"Warnings in full ({len(flagged)} read(s))",
                             expanded=True):
                for _, row in flagged.iterrows():
                    ref = row.get("file", row["read"])
                    where = f" ({ref})" if ref != row["read"] else ""
                    st.markdown(f"**{row['read']}**{where}, *{row['status']}*")
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
                img_q = png_quality(payload, mott_cutoff, names)
                if img_q:
                    st.image(img_q, **FIT)
                    st.caption("Blue is the window kept after trimming; the dashed "
                               "line is Phred 20.")
                    data_download("Phred per base",
                                  quality_table_for(payload, mott_cutoff, names),
                                  "read_quality.csv", key="dl_read_quality",
                                  help="One row per base of each read drawn, "
                                       "with the bases kept by the trim "
                                       "marked.")
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
