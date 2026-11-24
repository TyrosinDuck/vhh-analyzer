"""
The count of identical Sanger reads, carried from the Sanger screen into the
analysis in the FASTA header, after the name.

Sending one entry per unique protein, the default, used to drop how many reads
had given it, so a family picked seven times looked like one picked once. The
count is now written as `>name identical_reads=N` and read back from the
analysis input itself: one record of it, the text, and nothing beside it to
fall out of step.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import sanger as sg
import session as ses
import vhh_core as vc

_VHH = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
        "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA{cdr3}WGQGTQVTVSS")
A = _VHH.format(cdr3="DRGYSSSWYDY")
B = _VHH.format(cdr3="DRGYSSSWYDF")


# --------------------------------------------------------------------------- #
# Reading headers: one parser, so names always agree with parse_input
# --------------------------------------------------------------------------- #

def test_descriptions_are_keyed_as_parse_input_names_the_sequences():
    text = (f">Nb1 identical_reads=3 from plate 2\n{A}\n"
            f">\n{B}\n"                          # empty header: Sequence_2
            f">Nb3\n{A}\n")                      # no description
    names = list(vc.parse_input(text))
    descriptions = vc.fasta_descriptions(text)
    assert names == ["Nb1", "Sequence_2", "Nb3"]
    assert descriptions == {"Nb1": "identical_reads=3 from plate 2"}
    assert set(descriptions) <= set(names)


def test_a_repeated_name_keeps_its_last_description_like_its_last_sequence():
    text = f">Nb1 identical_reads=1\n{A}\n>Nb1 identical_reads=4\n{B}\n"
    assert vc.parse_input(text) == {"Nb1": B}
    assert vc.fasta_descriptions(text) == {"Nb1": "identical_reads=4"}


def test_raw_sequences_have_no_descriptions():
    assert vc.fasta_descriptions(A) == {}


def test_parse_input_is_unchanged_by_descriptions():
    plain = f">Nb1\n{A}\n>Nb2\n{B}\n"
    tagged = f">Nb1 identical_reads=2\n{A}\n>Nb2 anything else\n{B}\n"
    assert vc.parse_input(plain) == vc.parse_input(tagged)


# --------------------------------------------------------------------------- #
# The tag itself
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("desc, expected", [
    ("identical_reads=3", 3),
    ("clone from plate 2 identical_reads=12", 12),
    ("identical_reads=0", None),            # zero is not a count of reads
    ("identical_reads=3x", None),           # not exactly the tag
    ("my_identical_reads=3", None),
    ("identical_reads = 3", None),
    ("reads=3", None),
])
def test_only_the_exact_tag_is_read(desc, expected):
    found = sg.identical_reads({"Nb1": desc})
    assert found.get("Nb1") == expected


def test_the_writer_and_the_reader_agree():
    df = pd.DataFrame({"read": ["A1_fwd", "A1_rev", "B2"],
                       "vhh_aa": [A, A, B], "usable": [True] * 3,
                       "n_copies": [2, 2, 1]})
    text = sg.reads_to_fasta(df, counts=True)
    assert text.splitlines()[0] == ">A1_fwd identical_reads=2"
    assert sg.identical_reads(vc.fasta_descriptions(text)) == {
        "A1_fwd": 2, "A1_rev": 2, "B2": 1}
    assert list(vc.parse_input(text)) == ["A1_fwd", "A1_rev", "B2"]


def test_without_counts_the_fasta_is_as_before():
    df = pd.DataFrame({"read": ["A1"], "vhh_aa": [A], "usable": [True],
                       "n_copies": [2]})
    assert sg.reads_to_fasta(df).splitlines()[0] == ">A1"


# --------------------------------------------------------------------------- #
# The tables
# --------------------------------------------------------------------------- #

def test_the_column_sits_after_the_name_and_is_empty_not_zero_when_unknown():
    frame = pd.DataFrame({"name": ["A1", "X"], "family": ["F01", "F02"]})
    out = sg.add_identical_reads(frame, {"A1": 2})
    assert list(out.columns) == ["name", sg.READS_TAG, "family"]
    assert out[sg.READS_TAG].iloc[0] == 2
    assert pd.isna(out[sg.READS_TAG].iloc[1])
    assert list(frame.columns) == ["name", "family"]     # not modified


def test_no_counts_no_column():
    frame = pd.DataFrame({"name": ["A1"]})
    assert sg.add_identical_reads(frame, {}) is frame


def _families(assignment):
    return (assignment.groupby("family").size().rename("n_sequences")
            .reset_index())


def test_a_family_counts_each_distinct_protein_once():
    """The duplicate sent along with its original carries the same count;
    adding both would count those reads twice."""
    assignment = pd.DataFrame({"name": ["A1_fwd", "A1_rev", "B2"],
                               "family": ["F01", "F01", "F01"]})
    seqs = {"A1_fwd": A, "A1_rev": A, "B2": B}
    reads = {"A1_fwd": 2, "A1_rev": 2, "B2": 5}
    out = sg.family_identical_reads(_families(assignment), assignment, seqs,
                                    reads)
    assert out.loc[0, sg.READS_TAG] == 7
    cols = list(out.columns)
    assert cols.index(sg.READS_TAG) == cols.index("n_sequences") + 1


def test_a_family_with_an_uncounted_member_has_no_total():
    assignment = pd.DataFrame({"name": ["A1", "pasted"],
                               "family": ["F01", "F01"]})
    out = sg.family_identical_reads(_families(assignment), assignment,
                                    {"A1": A, "pasted": B}, {"A1": 3})
    assert pd.isna(out.loc[0, sg.READS_TAG])


# --------------------------------------------------------------------------- #
# Session files
# --------------------------------------------------------------------------- #

def test_a_session_keeps_the_counts_through_a_save_and_a_restore():
    record = ses.session_record({"Nb1": A, "Nb2": B}, 7.4, 9, "auto",
                                descriptions={"Nb1": "identical_reads=3",
                                              "gone": "identical_reads=9"})
    assert record["descriptions"] == {"Nb1": "identical_reads=3"}
    back = ses.session_inputs(json.loads(json.dumps(record)))
    text = ses.fasta(back["sequences"], descriptions=back["descriptions"])
    assert sg.identical_reads(vc.fasta_descriptions(text)) == {"Nb1": 3}
    assert vc.parse_input(text) == {"Nb1": A, "Nb2": B}


def test_a_session_file_from_before_reads_with_no_descriptions():
    old = json.loads((Path(__file__).parent / "compat" / "session_1.json")
                     .read_text(encoding="utf-8"))
    assert "descriptions" not in old
    assert ses.session_inputs(old)["descriptions"] == {}


def test_a_description_cannot_break_the_fasta_it_is_written_into():
    record = ses.session_record({"Nb1": A}, 7.4, 9, "auto")
    record["descriptions"] = {"Nb1": "identical_reads=3\nQVQ"}
    back = ses.session_inputs(record)
    text = ses.fasta(back["sequences"], descriptions=back["descriptions"])
    assert vc.parse_input(text) == {"Nb1": A}


# --------------------------------------------------------------------------- #
# The whole path, Sanger example run to clonal families
# --------------------------------------------------------------------------- #

@pytest.mark.smoke
def test_the_count_reaches_the_clonal_families_and_the_triage_table():
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    import shell
    app = str(Path(__file__).resolve().parent.parent / "app.py")
    at = AppTest.from_file(app, default_timeout=900)
    at.run()

    def click(label):
        at.button[[b.label for b in at.button].index(label)].click().run()

    def family(key):
        label = next(lab for lab, (k, *_r) in shell.ANALYSIS_FAMILIES.items()
                     if k == key)
        [r for r in at.radio if r.label == "Section"][0].set_value(label).run()
        assert not at.exception, [e.value for e in at.exception]

    click("🧪 Sanger sequencing")
    click("Load an example run")
    click("Send selected to analysis")
    headers = [line for line in at.session_state["seq_text"].splitlines()
               if line.startswith(">")]
    assert headers and all(f" {sg.READS_TAG}=" in h for h in headers)

    family("compare")
    tables = [d.value for d in at.dataframe if sg.READS_TAG in d.value.columns]
    families = next(t for t in tables if "n_sequences" in t.columns)
    # The example's clean clone is read from both ends: one clone, two reads.
    assert families[sg.READS_TAG].max() == 2
    assert families[sg.READS_TAG].notna().all()

    family("triage")
    assert any(sg.READS_TAG in d.value.columns for d in at.dataframe)
