"""Reads renamed from the provider's reference to the user's own names.

A sequencing provider names each chromatogram after its own order and well,
and nobody wants `24AB1234_C3_PelB-F` in a figure. The naming is in `sanger`
rather than in the screen so a notebook gets the same rules, and these tests
are about the rules: every read ends up with a name of its own, the reference
is never lost, and nothing a user typed is changed without saying so.

The example reads are simulated (see `demo.SANGER_NOTE`). They are used here
because what is checked is how names travel, not what a sequencer measures.
"""

from __future__ import annotations

import pandas as pd
import pytest

import demo
import sanger as sg
import vhh_core as vc


# --------------------------------------------------------------------------- #
# clean_name
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("given, expected", [
    ("Nb12", "Nb12"),
    ("  Nb12  ", "Nb12"),
    ("Nb 12", "Nb_12"),
    ("Nb \t 12\nx", "Nb_12_x"),
    ("", ""),
    ("   ", ""),
    (None, ""),
    (float("nan"), ""),
])
def test_a_name_is_made_safe_for_a_fasta_header(given, expected):
    assert sg.clean_name(given) == expected


def test_why_spaces_matter_the_analysis_reads_a_header_to_its_first_space():
    """The reason `clean_name` exists, checked on the parser it protects."""
    body = "QVQLVESGGGLVQAGGSLRLSCAASGRTFS"
    parsed = vc.parse_input(f">Nb 12\n{body}\n>Nb 13\n{body}\n")
    assert list(parsed) == ["Nb"], (
        "if this changes, clean_name may no longer be needed, but the two "
        "must be changed together")


# --------------------------------------------------------------------------- #
# assign_names
# --------------------------------------------------------------------------- #

def test_no_names_keeps_the_references():
    out = sg.assign_names(["A01", "B01"])
    assert out["name"].tolist() == ["A01", "B01"]
    assert not out["note"].astype(bool).any()


def test_an_empty_cell_keeps_its_reference():
    out = sg.assign_names(["A01", "B01", "C01"], ["Nb1", "", None])
    assert out["name"].tolist() == ["Nb1", "B01", "C01"]
    assert out["file"].tolist() == ["A01", "B01", "C01"]


def test_a_shorter_list_of_names_leaves_the_rest_alone():
    out = sg.assign_names(["A01", "B01"], ["Nb1"])
    assert out["name"].tolist() == ["Nb1", "B01"]


def test_forward_and_reverse_given_one_name_both_survive():
    """The usual way to arrive at a repeat, and both reads have to be kept."""
    out = sg.assign_names(["A1_fwd", "A1_rev"], ["Nb1", "Nb1"])
    assert out["name"].tolist() == ["Nb1", "Nb1_2"]
    assert out["note"].iloc[0] == ""
    assert "already used by A1_fwd" in out["note"].iloc[1]


def test_an_explicit_name_wins_over_a_reference_kept_by_default():
    """Renaming A to B's reference moves the suffix onto B, which nobody named."""
    out = sg.assign_names(["A", "B"], ["B", ""])
    assert out["name"].tolist() == ["B", "B_2"]
    assert out["note"].iloc[0] == ""
    assert "renamed B_2" in out["note"].iloc[1]


def test_a_suffix_never_lands_on_a_name_already_taken():
    out = sg.assign_names(["a", "b", "X_2"], ["X", "X", ""])
    names = out["name"].tolist()
    assert len(set(names)) == 3
    assert names[:2] == ["X", "X_2"]


def test_the_same_file_uploaded_twice_gives_two_names():
    out = sg.assign_names(["A01", "A01"], [])
    assert out["name"].tolist() == ["A01", "A01_2"]


def test_spaces_are_replaced_and_the_note_says_so():
    out = sg.assign_names(["A01"], ["Nb 12"])
    assert out["name"].iloc[0] == "Nb_12"
    assert "spaces replaced" in out["note"].iloc[0]


def test_outer_whitespace_alone_is_not_worth_a_note():
    out = sg.assign_names(["A01"], ["  Nb12 "])
    assert out["name"].iloc[0] == "Nb12"
    assert out["note"].iloc[0] == ""


def test_names_always_come_out_unique():
    stems = ["a", "b", "c", "d", "e", "f"]
    requested = ["x", "x", "x_2", "", "x", "d"]
    names = sg.assign_names(stems, requested)["name"].tolist()
    assert len(names) == len(set(names))


# --------------------------------------------------------------------------- #
# read_name_map
# --------------------------------------------------------------------------- #

STEMS = ["24AB_A01", "24AB_B01", "24AB_C01"]


@pytest.mark.parametrize("sep", [",", ";", "\t"])
def test_a_sheet_is_read_whatever_its_delimiter(sep):
    text = f"24AB_A01{sep}Nb1\n24AB_B01{sep}Nb2\n"
    mapping, unmatched = sg.read_name_map(text, STEMS)
    assert mapping == {"24AB_A01": "Nb1", "24AB_B01": "Nb2"}
    assert unmatched == []


def test_a_header_row_is_skipped():
    mapping, unmatched = sg.read_name_map(
        b"reference,name\n24AB_A01,Nb1\n", STEMS)
    assert mapping == {"24AB_A01": "Nb1"}
    assert unmatched == []


def test_file_names_with_the_extension_match_their_stems():
    mapping, _ = sg.read_name_map("24AB_A01.ab1,Nb1\n24AB_B01.AB1,Nb2\n", STEMS)
    assert mapping == {"24AB_A01": "Nb1", "24AB_B01": "Nb2"}


def test_a_utf8_bom_from_excel_does_not_hide_the_first_reference():
    mapping, _ = sg.read_name_map("24AB_A01;Nb1\n24AB_B01;Nb2\n"
                                  .encode("utf-8-sig"), STEMS)
    assert "24AB_A01" in mapping


def test_references_that_match_nothing_are_reported_not_dropped():
    mapping, unmatched = sg.read_name_map(
        "24AB_A01,Nb1\n24AB_Z99,Nb9\n", STEMS)
    assert mapping == {"24AB_A01": "Nb1"}
    assert unmatched == ["24AB_Z99"]


def test_extra_columns_are_ignored():
    mapping, _ = sg.read_name_map(
        "24AB_A01,Nb1,plate 3,round 2\n24AB_B01,Nb2,plate 3,round 2\n", STEMS)
    assert mapping == {"24AB_A01": "Nb1", "24AB_B01": "Nb2"}


def test_one_column_is_refused_with_a_readable_message():
    with pytest.raises(ValueError, match="two columns"):
        sg.read_name_map("24AB_A01\n24AB_B01\n", STEMS)


def test_matching_is_case_sensitive():
    mapping, unmatched = sg.read_name_map("24ab_a01,Nb1\n24AB_B01,Nb2\n", STEMS)
    assert mapping == {"24AB_B01": "Nb2"}


# --------------------------------------------------------------------------- #
# Through the pipeline
# --------------------------------------------------------------------------- #

def _example():
    return demo.sanger_example_files()


def test_without_names_the_pipeline_is_unchanged():
    files = _example()
    df, reads = sg.process_reads(files)
    assert df["read"].tolist() == df["file"].tolist()
    assert [r["read"] for r in reads] == [r["file_stem"] for r in reads]


def test_names_reach_the_table_the_duplicates_and_the_fasta():
    files = _example()
    stems = [name.rsplit(".", 1)[0] for name, _ in files]
    names = sg.assign_names(stems, ["Nb_A1", "Nb_A1"])["name"].tolist()
    df, reads = sg.process_reads(files, names=names)

    assert df["read"].tolist()[:2] == ["Nb_A1", "Nb_A1_2"]
    assert df["file"].tolist() == stems, "the reference was lost"
    assert [r["read"] for r in reads][:2] == ["Nb_A1", "Nb_A1_2"]

    # forward and reverse of the same clone: the reverse duplicates the
    # forward, and says so under the forward's new name
    assert df.loc[df["read"] == "Nb_A1_2", "duplicate_of"].iloc[0] == "Nb_A1"

    fasta = sg.reads_to_fasta(df)
    assert ">Nb_A1\n" in fasta
    assert ">clone_A1_fwd\n" not in fasta


def test_renamed_reads_arrive_in_the_analysis_under_their_names():
    """What the user asked for: their names used for the rest of the run."""
    files = _example()
    stems = [name.rsplit(".", 1)[0] for name, _ in files]
    names = sg.assign_names(stems, ["Nb 1", "", "Nb2"])["name"].tolist()
    df, _ = sg.process_reads(files, names=names)
    parsed = vc.parse_input(sg.reads_to_fasta(df))
    assert "Nb_1" in parsed and "Nb2" in parsed
    assert len(parsed) == int(df["usable"].sum())


def test_an_unreadable_file_keeps_both_its_name_and_its_reference():
    df, _ = sg.process_reads([("broken_A01.ab1", b"not a chromatogram")],
                             names=["Nb9"])
    row = df.iloc[0]
    assert row["status"] == "unreadable"
    assert row["read"] == "Nb9" and row["file"] == "broken_A01"


def test_the_quality_figure_uses_the_names():
    data = _example()[0][1]
    read = sg.read_ab1(data, name="clone_A1_fwd.ab1", label="Nb_A1")
    table = sg.quality_table([read])
    assert set(table["read"]) == {"Nb_A1"}
    assert read["file_stem"] == "clone_A1_fwd"


def test_assign_names_returns_a_frame_with_its_documented_columns():
    out = sg.assign_names(["A"], ["B"])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["file", "requested", "name", "note"]
