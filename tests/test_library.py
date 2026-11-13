"""
The sequence library: the public file's checks, and the private library that
is read from an upload, edited in the app and written back as a download.

These exercise `library.py` alone. The screen is driven by the smoke tests,
which cannot hand a file to `st.file_uploader` nor type into `st.data_editor`;
what those two widgets do is checked by hand, and the functions they feed are
checked here.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd
import pytest

import library as lib

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def public():
    table, problems = lib.load_library()
    assert not problems
    return table


def _frame(rows):
    return pd.DataFrame(rows, columns=lib.COLUMNS)


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #

def test_an_upload_reads_like_the_file_in_the_repository(public):
    """The same bytes give the same table whichever door they come in by."""
    data = lib.LIBRARY_FILE.read_bytes()
    table, problems = lib.read_library_bytes(data, "library.csv")
    assert not problems
    pd.testing.assert_frame_equal(table, public)


def test_an_excel_byte_order_mark_is_accepted(public):
    data = "﻿".encode() + lib.LIBRARY_FILE.read_bytes()
    table, problems = lib.read_library_bytes(data, "excel.csv")
    assert not problems
    assert list(table["name"]) == list(public["name"])


@pytest.mark.parametrize("data, says", [
    (b"", "is empty"),
    (b"\xff\xfe\x00n", "not UTF-8"),
    (b"title,sequence\nx,QVQ\n", "missing the column"),
])
def test_an_unreadable_upload_says_why_and_returns_nothing(data, says):
    table, problems = lib.read_library_bytes(data, "mine.csv")
    assert not len(table)
    assert len(problems) == 1 and says in problems[0]
    assert "mine.csv" in problems[0], "the message should name the file"


def test_a_bad_row_in_an_upload_is_kept_and_flagged():
    data = (b"name,protein,dna,tags,note\n"
            b"good,QVQLQESGGG,,,\n"
            b"good,QVQLQ,,,\n"
            b",,,,only a note\n")
    table, problems = lib.read_library_bytes(data, "mine.csv")
    assert len(table) == 3
    assert any("line 3" in p and "duplicate name" in p for p in problems)
    assert any("line 4" in p and "no usable sequence" in p for p in problems)


def test_an_empty_row_left_by_the_editor_is_skipped_without_renumbering():
    """`st.data_editor` adds a blank row before you type into it. It is not a
    problem to report, and the rows after it keep their line numbers."""
    frame = _frame([["a", "QVQ", "", "", ""],
                    ["", "", "", "", ""],
                    ["a", "QVQ", "", "", ""]])
    table, problems = lib.validate_library(frame, "your library")
    assert list(table["name"]) == ["a", "a"]
    assert problems == ["line 4 (a): duplicate name"]


def test_missing_cells_from_the_editor_count_as_empty():
    """A row added in the editor comes back with None, not ''."""
    frame = _frame([["a", "QVQLQESGG", None, None, None]])
    table, problems = lib.validate_library(frame)
    assert not problems
    assert table.loc[0, "tags"] == "" and table.loc[0, "dna"] == ""


# --------------------------------------------------------------------------- #
# Writing back
# --------------------------------------------------------------------------- #

def test_a_download_reads_back_to_the_same_library(public):
    csv = lib.library_csv(public)
    again, problems = lib.read_library_bytes(csv.encode(), "round trip")
    assert not problems
    pd.testing.assert_frame_equal(again, public)


def test_the_download_carries_only_the_columns_a_person_edits(public):
    header = lib.library_csv(public).splitlines()[0]
    assert header == ",".join(lib.COLUMNS)


def test_blank_rows_are_not_written():
    frame = _frame([["a", "QVQ", "", "", ""], [None, "", " ", None, ""]])
    assert lib.library_csv(frame).splitlines() == [",".join(lib.COLUMNS),
                                                   "a,QVQ,,,"]


def test_the_digest_follows_what_would_be_saved_and_nothing_else(public):
    """It decides whether the screen says 'changes not downloaded'. A blank
    row, or a derived column, does not change the file and must not raise
    the warning; a changed sequence must."""
    digest = lib.library_digest(public)
    assert lib.library_digest(lib.editable_frame(public)) == digest
    with_blank = pd.concat([lib.editable_frame(public),
                            _frame([["", "", "", "", ""]])], ignore_index=True)
    assert lib.library_digest(with_blank) == digest

    changed = lib.editable_frame(public)
    changed.loc[0, "protein"] += "K"
    assert lib.library_digest(changed) != digest


def test_editable_frame_of_nothing_has_the_columns():
    frame = lib.editable_frame()
    assert list(frame.columns) == lib.COLUMNS and not len(frame)


# --------------------------------------------------------------------------- #
# Copying from the public library
# --------------------------------------------------------------------------- #

def test_copying_adds_the_entries_with_their_sequences(public):
    names = list(public["name"][:2])
    mine, added, skipped = lib.add_entries(lib.editable_frame(), public, names)
    assert added == names and not skipped
    source = public.set_index("name")
    for _, row in mine.iterrows():
        assert row["protein"] == source.loc[row["name"], "protein"]
        assert row["dna"] == source.loc[row["name"], "dna"]


def test_copying_never_overwrites_or_duplicates_a_name(public):
    name = public["name"].iloc[0]
    mine = _frame([[name, "EDITEDSEQ", "", "", "my version"]])
    mine, added, skipped = lib.add_entries(mine, public,
                                           [name, public["name"].iloc[1]])
    assert skipped == [name]
    assert added == [public["name"].iloc[1]]
    assert mine.loc[mine["name"] == name, "protein"].tolist() == ["EDITEDSEQ"]


def test_copying_leaves_the_public_table_alone(public):
    before = public.copy()
    lib.add_entries(lib.editable_frame(), public, list(public["name"]))
    pd.testing.assert_frame_equal(public, before)


# --------------------------------------------------------------------------- #
# The private library stays out of the shared cache
# --------------------------------------------------------------------------- #

def _calls_in(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            names.add(node.attr)
        elif isinstance(node, ast.Name):
            names.add(node.id)
    return names


PRIVATE_FUNCTIONS = {"read_library_bytes", "validate_library", "library_csv",
                     "library_digest", "add_entries", "editable_frame"}


def test_compute_does_not_cache_the_private_library():
    """`st.cache_data` is shared by every session on the server. compute.py
    is where caching lives, so the private library's functions have no
    business being called there; somebody "harmonising" the two libraries
    by caching both would leak one person's sequences into a cache that
    outlives their session."""
    leaked = PRIVATE_FUNCTIONS & _calls_in(ROOT / "compute.py")
    assert not leaked, f"compute.py calls {leaked}"


def test_the_library_screen_caches_nothing_itself():
    names = _calls_in(ROOT / "sections" / "library.py")
    assert not {"cache_data", "cache_resource"} & names


# --------------------------------------------------------------------------- #
# Comments: whole lines only
# --------------------------------------------------------------------------- #

def test_a_hash_inside_a_row_is_data_not_a_comment():
    """Regression. The reader used `pd.read_csv(comment="#")`, which ends a
    line at any `#`: `NbNAMPT_cluster#1-(18)` lost everything after the `#`,
    sequence included, and two clusters ended up with the same name. Found
    on a real private library converted from a spreadsheet."""
    data = (b"name,protein,dna,tags,note\n"
            b"# a comment line\n"
            b"   # an indented one\n"
            b"Nb_cluster#1-(18),QVQLQESGGGLVQAGG,,,see # 4\n"
            b"Nb_cluster#2-(5),QVQLQESGGGLVQAGK,,,\n")
    table, problems = lib.read_library_bytes(data, "mine.csv")
    assert not problems
    assert list(table["name"]) == ["Nb_cluster#1-(18)", "Nb_cluster#2-(5)"]
    assert list(table["protein"]) == ["QVQLQESGGGLVQAGG", "QVQLQESGGGLVQAGK"]
    assert table.loc[0, "note"] == "see # 4"


def test_a_name_with_a_hash_survives_a_download_and_reload():
    frame = _frame([["Nb#7", "QVQLQESGG", "", "a#b", "c # d"]])
    table, problems = lib.read_library_bytes(lib.library_csv(frame).encode())
    assert not problems
    assert list(table.iloc[0][lib.COLUMNS]) == ["Nb#7", "QVQLQESGG", "",
                                                 "a#b", "c # d"]


@pytest.mark.parametrize("name", ["NA", "None", "null", "NaN", "N/A"])
def test_a_name_pandas_would_read_as_missing_is_kept(name):
    data = f"name,protein,dna,tags,note\n{name},QVQLQESGG,,,\n".encode()
    table, problems = lib.read_library_bytes(data)
    assert not problems and table.loc[0, "name"] == name


def test_the_public_file_still_skips_its_comment_header():
    text = lib.LIBRARY_FILE.read_text(encoding="utf-8")
    assert text.startswith("#"), "the test assumes the file opens with comments"
    table, problems = lib.load_library()
    assert not problems and not any(n.startswith("#") for n in table["name"])


# --------------------------------------------------------------------------- #
# The Send column of the private library's editor
# --------------------------------------------------------------------------- #

def test_a_send_tick_is_neither_saved_nor_a_change():
    """The Send column lives in the editor only: ticking a row must not mark
    the library as changed, nor end up in the downloaded file."""
    import sections.library as screen
    public, _ = lib.load_library()
    plain = screen.with_send(public)
    ticked = plain.copy()
    ticked.loc[[0, 2], "send"] = True
    assert lib.library_digest(ticked) == lib.library_digest(plain)
    assert "send" not in lib.library_csv(ticked).splitlines()[0]
    assert screen.ticked_names(ticked) == [public["name"].iloc[0],
                                           public["name"].iloc[2]]


def test_a_send_tick_follows_its_row_through_an_edit():
    import sections.library as screen
    public, _ = lib.load_library()
    frame = screen.with_send(public.head(3))
    frame.loc[1, "send"] = True
    frame = frame.drop(index=0).reset_index(drop=True)          # delete a row
    frame.loc[0, "protein"] += "K"                              # edit another
    again = screen.with_send(frame)
    assert screen.ticked_names(again) == [public["name"].iloc[1]]


def test_a_row_added_in_the_editor_comes_back_unticked():
    import sections.library as screen
    edited = pd.DataFrame({"send": [True, None], "name": ["a", "b"],
                           "protein": ["QVQ", None], "dna": [None, None],
                           "tags": [None, None], "note": [None, None]})
    frame = screen.with_send(edited)
    assert frame["send"].tolist() == [True, False]
    assert screen.ticked_names(frame) == ["a"]


@pytest.mark.parametrize("ticked, shown, now, expected", [
    ([], ["a", "b"], ["b"], ["b"]),                       # tick one
    (["b"], ["a", "b"], [], []),                          # untick it
    (["x"], ["a", "b"], ["a"], ["x", "a"]),               # hidden tick kept
    (["x", "a"], ["a"], ["a"], ["x", "a"]),               # no duplicates
    (["a"], [], [], ["a"]),                               # all filtered out
])
def test_merge_ticks_touches_only_the_rows_shown(ticked, shown, now, expected):
    import sections.library as screen
    assert screen.merge_ticks(ticked, shown, now) == expected
