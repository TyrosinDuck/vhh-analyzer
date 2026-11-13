"""Motif detection on a plasmid, where a repeated residue used to be counted many times.

`cloning.detect_motifs` searches for short protein motifs in all six frames. A
motif made of one repeated residue matches at every offset inside a longer run
of it, so the His6 motif matched five times inside a His10: the feature
inventory and the vector map showed five overlapping tags, at 70-87, 73-90,
76-93, 79-96 and 82-99, where the vector carries one, and every one of them was
labelled His6 on a vector with ten histidines.

These tests pin the run being reported once and named for its real length, and
pin that the five bundled vectors, which all carry six, still read His6.
"""

from __future__ import annotations

from Bio.Data import CodonTable

import cloning as cl

BACK_TABLE = CodonTable.unambiguous_dna_by_id[1].back_table
LEADER = "MKYLLPTAAAGLLLLAAQPAMA"


def _plasmid_encoding(protein: str) -> dict:
    """A linear record whose plus strand encodes `protein`, padded either side."""
    dna = "AAA" + "".join(BACK_TABLE[a] for a in protein) + "AAA"
    return cl.load_plasmid(f">synthetic\n{dna}\n".encode(), name="synthetic.fasta")


def _his_hits(protein: str) -> list[dict]:
    plasmid = _plasmid_encoding(protein)
    return [hit for hit in cl.detect_motifs(plasmid["record"], circular=False)
            if hit["feature"].startswith("His")]


def test_a_six_histidine_run_is_still_his6():
    hits = _his_hits(LEADER + "H" * 6)
    assert len(hits) == 1
    assert hits[0]["feature"] == "His6"
    assert hits[0]["end"] - hits[0]["start"] + 1 == 18


def test_a_longer_run_is_reported_once_and_named_for_its_length():
    for length in (7, 8, 10, 12):
        hits = _his_hits(LEADER + "H" * length)
        assert len(hits) == 1, (
            f"a run of {length} histidines gave {len(hits)} features: "
            f"{[(h['feature'], h['start'], h['end']) for h in hits]}")
        assert hits[0]["feature"] == f"His{length}"
        assert hits[0]["protein"] == "H" * length
        assert hits[0]["end"] - hits[0]["start"] + 1 == 3 * length


def test_two_separate_runs_stay_two_features():
    hits = _his_hits(LEADER + "H" * 8 + "GGSGGS" + "H" * 6)
    assert sorted(h["feature"] for h in hits) == ["His6", "His8"]


def test_a_run_shorter_than_the_motif_is_not_reported():
    assert _his_hits(LEADER + "HHHHH") == []


def test_the_bundled_vectors_all_read_his6():
    for name, path in cl.list_bundled_plasmids().items():
        plasmid = cl.load_plasmid(path)
        hits = [h for h in cl.detect_motifs(plasmid["record"])
                if h["feature"].startswith("His")]
        assert [h["feature"] for h in hits] == ["His6"], (name, hits)


def test_the_feature_table_shows_one_row_for_one_run():
    plasmid = _plasmid_encoding(LEADER + "H" * 10)
    table = cl.feature_table(plasmid)
    rows = table[table["feature"].str.startswith("His")]
    assert len(rows) == 1
    assert rows.iloc[0]["feature"] == "His10"


def test_a_longer_run_does_not_break_the_cloning_window():
    """`insertion_region` matches tag labels by substring, `his` among them.

    Renaming the feature to His10 must not take it out of that match, or the
    window between the leader and the first tag would stop being found on a
    vector with a longer tag.
    """
    plasmid = _plasmid_encoding(LEADER + "QVQLQ" + "GGSGGSGGS" + "H" * 10)
    region = cl.insertion_region(plasmid)
    assert region is not None
