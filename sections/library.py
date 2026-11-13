"""
sections.library - The public sequence library, read-only in the app by design,
and a private one held in the browser session.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import hashlib
import re

import pandas as pd
import streamlit as st

import library as lib
from compute import load_library_cached
from shell import FIT, Settings, go_to

__version__ = "2026.11.13c"

# Everything the private library keeps between reruns lives under keys with
# this prefix in st.session_state, and nowhere else: not on disk, not in
# st.cache_data (which every session on the server shares).
_MINE = "lib_mine_"


# --------------------------------------------------------------------------- #
# The private library's state
# --------------------------------------------------------------------------- #

def _has_mine() -> bool:
    return f"{_MINE}current" in st.session_state


def _editor_key() -> str:
    return f"{_MINE}editor_{st.session_state.get(f'{_MINE}gen', 0)}"


def _replace_mine(frame, file_name: str | None = None,
                  saved: bool = False) -> None:
    """Make ``frame`` the private library and give the editor a fresh key.

    The editor's own state holds the edits as differences from the frame it
    was given. Handing it a different frame under the same key would apply
    those differences to the wrong rows, so every replacement goes through a
    new key instead.
    """
    frame = with_send(frame)
    st.session_state[f"{_MINE}base"] = frame
    st.session_state[f"{_MINE}current"] = frame
    st.session_state[f"{_MINE}gen"] = st.session_state.get(f"{_MINE}gen", 0) + 1
    if file_name is not None:
        st.session_state[f"{_MINE}file"] = file_name
    if saved:
        st.session_state[f"{_MINE}saved"] = lib.library_digest(frame)


SEND = "send"


def with_send(frame=None):
    """The private library as the editor shows it: a ``send`` tick box,
    then the columns of the file.

    The tick belongs to the row, so it follows its entry through an edit.
    It is not part of the library: `library.editable_frame`, and therefore
    the download and the unsaved-changes digest, leave it out, so ticking a
    row neither changes the file nor counts as a change.
    """
    ticks = None
    if frame is not None and SEND in getattr(frame, "columns", []):
        ticks = frame[SEND].fillna(False).astype(bool).reset_index(drop=True)
    out = lib.editable_frame(frame)
    out.insert(0, SEND, ticks if ticks is not None
               else pd.Series(False, index=out.index, dtype=bool))
    return out


def ticked_names(frame) -> list[str]:
    """Names of the rows whose ``send`` box is ticked, as the library check
    will spell them; a ticked row without a name is left out."""
    if frame is None or SEND not in frame.columns:
        return []
    names = frame.loc[frame[SEND].fillna(False).astype(bool), "name"]
    return [n.strip() for n in names.fillna("").astype(str) if n.strip()]


def _close_mine() -> None:
    for key in [k for k in st.session_state if str(k).startswith(_MINE)]:
        del st.session_state[key]


def _mark_saved(digest: str) -> None:
    st.session_state[f"{_MINE}saved"] = digest


def _flash(message: str, kind: str = "success", where: str = "mine") -> None:
    """A message for the next run, shown in the tab where the action was."""
    st.session_state[f"{_MINE}flash_{where}"] = (kind, message)


def _show_flash(where: str) -> None:
    flash = st.session_state.pop(f"{_MINE}flash_{where}", None)
    if flash:
        getattr(st, flash[0])(flash[1])


# --------------------------------------------------------------------------- #
# Shared by both tabs
# --------------------------------------------------------------------------- #

def _table_key(prefix: str, *state) -> str:
    """Key of a selectable table, changed whenever what it shows changes.

    Streamlit records a selection as row *positions*. If the rows under them
    changed, a filter typed or a cell edited, the same positions would
    silently pick other entries; a new key clears the selection instead.
    """
    digest = hashlib.sha256(repr(state).encode()).hexdigest()[:10]
    return f"{prefix}_table_{digest}"


def _selectable_table(view, key: str, column_config: dict) -> list[str]:
    """Draw ``view`` with a tick box on each row; return the ticked names."""
    event = st.dataframe(view, hide_index=True, **FIT, key=key,
                         on_select="rerun", selection_mode="multi-row",
                         column_config=column_config)
    rows = [r for r in event.selection.rows if 0 <= r < len(view)]
    return list(view.iloc[rows]["name"])


_COLUMN_CONFIG = {
    "has_protein": st.column_config.CheckboxColumn("Protein"),
    "has_dna": st.column_config.CheckboxColumn("DNA"),
    "protein_length": st.column_config.NumberColumn("aa"),
    "dna_length": st.column_config.NumberColumn("nt"),
    "note": st.column_config.TextColumn("Note", width="large"),
    "issue": st.column_config.TextColumn("Problem", width="medium")}


def _send_block(table, picked: list[str], prefix: str,
                offer_copy: bool, how: str = "Tick the box at the left of "
                "one or more rows in the table to send them on.") -> None:
    """Send the ticked entries of ``table`` to a section, as FASTA, or,
    from the public library, into the private one."""
    st.subheader("Send to a section")
    if not picked:
        st.caption(how)
        return
    st.caption(f"Selected: {', '.join(picked)}")

    chosen = table[table["name"].isin(picked)]
    n_prot = int(chosen["has_protein"].sum())
    n_dna = int(chosen["has_dna"].sum())
    send = st.columns(2)
    with send[0]:
        if st.button(f"🔬 Send {n_prot} protein(s) to analysis",
                     type="primary", **FIT, disabled=not n_prot,
                     key=f"{prefix}_to_analysis"):
            st.session_state["pending_seq_text"] = lib.library_fasta(
                table, picked, "protein")
            go_to("analysis")
            st.rerun()
        if n_prot < len(picked):
            st.caption(f"{len(picked) - n_prot} of the selected entries "
                       "carry no protein sequence and will be left out.")
    with send[1]:
        if st.button(f"🔧 Send {n_dna} DNA sequence(s) to gene diagnosis",
                     type="primary", **FIT, disabled=not n_dna,
                     key=f"{prefix}_to_genes"):
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
                          file_name=f"{prefix}_proteins.fasta",
                          mime="text/plain", disabled=not n_prot,
                          key=f"{prefix}_dl_protein")
    dl[1].download_button("Selected DNA (FASTA)",
                          lib.library_fasta(table, picked, "dna"),
                          file_name=f"{prefix}_dna.fasta",
                          mime="text/plain", disabled=not n_dna,
                          key=f"{prefix}_dl_dna")

    if offer_copy:
        _show_flash("public")
        label = (f"➕ Copy {len(picked)} to my library" if _has_mine()
                 else f"➕ Start my library with these {len(picked)}")
        if st.button(label, **FIT, key=f"{prefix}_to_mine"):
            mine = (st.session_state[f"{_MINE}current"] if _has_mine()
                    else lib.editable_frame())
            if not _has_mine():
                # A library that has never been saved: nothing matches it.
                st.session_state[f"{_MINE}saved"] = ""
            mine, added, skipped = lib.add_entries(mine, table, picked)
            _replace_mine(mine)
            parts = []
            if added:
                parts.append(f"Copied to your library: {', '.join(added)}.")
            if skipped:
                parts.append("Already in your library under the same name, "
                             f"left as they were: {', '.join(skipped)}.")
            _flash(" ".join(parts), "success" if added else "warning",
                   where="public")
            st.rerun()
        st.caption("Copies the entries into the **My library** tab. The public "
                   "library is not changed; your copy is yours to edit.")


# --------------------------------------------------------------------------- #
# The tabs
# --------------------------------------------------------------------------- #

def _public_tab() -> None:
    table, problems = load_library_cached()

    st.info(
        "**Read-only here, and public.** Streamlit rebuilds this app from the "
        "repository on every reboot and redeploy, so anything written from the "
        "browser would be erased, the library lives in `data/library.csv` "
        "instead, where it is durable and versioned. Edit it on your own "
        "machine, commit and push, and it appears here on the next deploy.\n\n"
        "Because it sits in a public repository it is readable by anyone, "
        "whether or not they open this app. **Nothing unpublished or "
        "confidential belongs in it.** For your own sequences, use the "
        "**My library** tab.", icon="📚")

    for message in problems:
        st.warning(message, icon="⚠️")

    if not len(table):
        st.info("The library is empty. Add rows to `data/library.csv`, the "
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

        picked = _selectable_table(
            view[["name", "has_protein", "has_dna", "protein_length",
                  "dna_length", "tags", "note", "issue"]],
            _table_key("library", search, chosen_tags, needs),
            _COLUMN_CONFIG)
        st.caption(f"{len(view)} of {len(table)} entr{'y' if len(table) == 1 else 'ies'} "
                   "shown. Changing a filter clears the ticks.")

        st.divider()
        _send_block(table, picked, "library", offer_copy=True)

    with st.expander("How to add or remove entries"):
        st.markdown(
            "The library is the file `data/library.csv` in the repository. Edit "
            "it in a spreadsheet or a text editor on your own machine, commit "
            "it, and push, the app picks it up on its next deploy.\n\n"
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


def _mine_tab() -> None:
    st.info(
        "**Yours, and held in this browser session only.** Load a CSV in the "
        "same format as the public library, edit it here, and download it to "
        "keep it. The app never writes it to disk and never keeps it once the "
        "session ends: closing or reloading the tab loses it, and so, on the "
        "public site, does the app going to sleep. **Download before you "
        "leave.**\n\n"
        "On the public site the file passes through the memory of a server "
        "that is not yours, Streamlit Community Cloud. For confidential "
        "sequences, run the app on your own machine, where nothing leaves "
        "your computer.", icon="🔒")

    _show_flash("mine")

    upload = st.file_uploader("Load a library (CSV)", type=["csv"],
                              key=f"{_MINE}upload")
    if upload is not None and upload.file_id != st.session_state.get(
            f"{_MINE}upload_id"):
        dirty = _has_mine() and lib.library_digest(
            st.session_state[f"{_MINE}current"]) != st.session_state.get(
                f"{_MINE}saved")
        if dirty:
            st.warning(f"Loading **{upload.name}** replaces the library open "
                       "here, which has changes you have not downloaded.")
        if not dirty or st.button("Replace it, discarding those changes",
                                  key=f"{_MINE}replace"):
            table, problems = lib.read_library_bytes(upload.getvalue(),
                                                     upload.name)
            st.session_state[f"{_MINE}upload_id"] = upload.file_id
            if problems and not len(table):
                for message in problems:
                    st.error(message)
            else:
                _replace_mine(table, file_name=upload.name, saved=True)
                _flash(f"Loaded {upload.name}: {len(table)} "
                       f"entr{'y' if len(table) == 1 else 'ies'}.")
                st.rerun()

    if not _has_mine():
        st.caption("No library loaded. Load a file above, start an empty one, "
                   "or copy entries from the public library in the other tab.")
        if st.button("Start an empty library", key=f"{_MINE}new"):
            _replace_mine(lib.editable_frame(), file_name=None)
            st.session_state[f"{_MINE}saved"] = ""
            st.rerun()
        return

    # The editor's own state is dropped by Streamlit whenever this screen is
    # not drawn, so coming back from another section would show the frame it
    # was first given, without the edits. What the editor last returned is
    # kept separately, and becomes its starting point again in that case.
    if _editor_key() not in st.session_state:
        st.session_state[f"{_MINE}base"] = st.session_state[f"{_MINE}current"]

    edited = st.data_editor(
        st.session_state[f"{_MINE}base"], key=_editor_key(),
        num_rows="dynamic", hide_index=True, **FIT,
        column_config={
            SEND: st.column_config.CheckboxColumn(
                "Send", default=False, width="small",
                help="Tick to send this entry with the buttons below. "
                     "Not saved in the file."),
            "name": st.column_config.TextColumn("name", width="medium"),
            "protein": st.column_config.TextColumn("protein", width="large"),
            "dna": st.column_config.TextColumn("dna", width="large"),
            "tags": st.column_config.TextColumn("tags"),
            "note": st.column_config.TextColumn("note", width="large")})
    current = with_send(edited)
    st.session_state[f"{_MINE}current"] = current
    picked = ticked_names(current)
    st.caption("Tick **Send** for the entries to send on. Click a cell to "
               "edit it, the last row to add one, and the box at the far "
               "left of a row then the bin icon to delete it.")

    table, problems = lib.validate_library(current, "your library")
    for message in problems:
        st.warning(message, icon="⚠️")

    if len(table):
        st.divider()
        _send_block(table, picked, "my_library", offer_copy=False,
                    how="Tick **Send** on one or more rows of the table above "
                        "to send them on.")

    st.divider()
    st.subheader("Save")
    csv = lib.library_csv(current)
    digest = lib.library_digest(current)
    if digest == st.session_state.get(f"{_MINE}saved"):
        st.success("No changes since the file was loaded or last downloaded.")
    else:
        st.warning("**Changes not downloaded yet.** They exist only in this "
                   "browser tab.")
    file_name = st.session_state.get(f"{_MINE}file") or "my_library.private.csv"
    save = st.columns([2, 1])
    save[0].download_button(
        f"💾 Download {file_name}", csv, file_name=file_name, mime="text/csv",
        type="primary", **FIT, key=f"{_MINE}download",
        on_click=_mark_saved, args=(digest,))
    st.caption("The app counts a click on this button as saved. It cannot see "
               "whether your browser kept the file, or where.")
    with save[1]:
        discard = True
        if digest != st.session_state.get(f"{_MINE}saved"):
            discard = st.checkbox("Discard my unsaved changes",
                                  key=f"{_MINE}discard")
        if st.button("Close my library", disabled=not discard, **FIT,
                     key=f"{_MINE}close"):
            _close_mine()
            st.rerun()



def render(settings: Settings) -> None:
    """Draw the screen.

    The router gives every screen the same argument. This one reads
    none of it: see `shell.Settings` for why that is cheaper than a
    second list of the screens.
    """
    st.title("📚 Sequence library")
    st.markdown(
        "Two libraries. The **public** one is kept in the repository, the "
        "same for everyone. **Yours** is a file you load, edit here and "
        "download again; the app does not keep it. Pick what you need and "
        "send it to the section that will work on it.")

    public, mine = st.tabs(["Public library", "My library"])
    # Mine first: the public tab's copy button reads the private library as
    # the editor last returned it, so the editor has to have run already.
    with mine:
        _mine_tab()
    with public:
        _public_tab()

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
