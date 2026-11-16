"""
Every figure's CSV is the figure, checked where it can be.

The app offers the data behind each figure next to it, and the rule that makes
that worth anything is that the download and the drawing come from one call:
the plot functions ask these functions for their numbers rather than computing
their own. A test cannot read a PNG and compare it with a frame, so what is
checked here is the half that is checkable:

* the data function returns the shape and the values it claims, against an
  independent computation where one is cheap (Biopython for the titration
  curves, the profile's own smoothing for the tracks);
* the plot function still draws from it, by handing the plot a set of numbers
  and watching the figure follow, and by patching the data function and
  checking the figure notices.

The second kind is the one that would catch a future edit reintroducing a
second computation inside a plot function, which is the failure this whole
arrangement exists to prevent.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import sanger as sg
import vhh_core as vc

SEQ = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
       "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDYWGQGTQVTVSS")
SEQS = {"Nb1": SEQ, "Nb2": SEQ[:-1] + "A"}


# --------------------------------------------------------------------------- #
# Titration curves
# --------------------------------------------------------------------------- #

def test_the_titration_curves_span_the_grid_the_figure_is_drawn_on():
    curves = vc.titration_curves(SEQS)
    first, last, step = vc.TITRATION_GRID
    per_sequence = len(curves) / len(SEQS)
    assert per_sequence == round((last - first) / step) + 1
    assert curves["pH"].min() == pytest.approx(first)
    assert curves["pH"].max() == pytest.approx(last)
    assert set(curves["name"]) == set(SEQS)


def test_every_charge_is_biopythons_charge_at_that_ph():
    """The claim the caption makes: this is `charge_at_pH`, not a model of it."""
    curves = vc.titration_curves({"Nb1": SEQ})
    pa = ProteinAnalysis(SEQ)
    for row in curves.sample(12, random_state=0).itertuples():
        assert row.net_charge == pytest.approx(pa.charge_at_pH(row.pH))


def test_the_charge_figure_follows_the_curves_it_is_given(monkeypatch):
    """Patch the data, and the drawing has to change with it.

    If `plot_charge` went back to computing its own titration, this test would
    pass a sequence set the patched function never sees and the line drawn
    would be the real curve instead of the flat one.
    """
    flat = pd.DataFrame({"name": ["Nb1"] * 3, "pH": [2.0, 7.0, 12.0],
                         "net_charge": [5.0, 5.0, 5.0]})
    monkeypatch.setattr(vc, "titration_curves", lambda seqs: flat)
    frame = pd.DataFrame({"name": ["Nb1"], "net_charge": [5.0]})
    fig = vc.plot_charge(frame, {"Nb1": SEQ}, ph=7.0)
    line = fig.axes[0].get_lines()[0]
    assert list(line.get_ydata()) == [5.0, 5.0, 5.0]
    assert list(line.get_xdata()) == [2.0, 7.0, 12.0]


# --------------------------------------------------------------------------- #
# Per-residue profiles
# --------------------------------------------------------------------------- #

def test_the_profile_table_has_one_row_per_residue_of_each_sequence():
    table = vc.profile_table(SEQS, None, 9)
    assert len(table) == sum(len(s) for s in SEQS.values())
    assert list(table.columns) == ["name", "position", "residue", "hydropathy",
                                   "charge", "region", "window"]
    one = table[table["name"] == "Nb1"]
    assert list(one["residue"]) == list(SEQ)
    assert one["position"].tolist() == list(range(1, len(SEQ) + 1))


def test_the_profile_tracks_are_the_apps_own_smoothing():
    table = vc.profile_table({"Nb1": SEQ}, None, 9)
    expected = vc.sliding(np.array([vc.KD[a] for a in SEQ], float), 9)
    assert table["hydropathy"].to_numpy() == pytest.approx(expected)


def test_the_window_travels_with_the_profile():
    """A profile without its window is not comparable with another one."""
    assert set(vc.profile_table(SEQS, None, 5)["window"]) == {5}
    assert set(vc.profile_table(SEQS, None, 11)["window"]) == {11}


def test_the_cdr_of_each_residue_is_marked_when_an_annotation_is_given():
    annotation, _method = vc.annotate(SEQ)
    table = vc.profile_table({"Nb1": SEQ}, {"Nb1": annotation}, 9)
    cdr3 = annotation["CDR3"]
    inside = table[(table["position"] >= cdr3["start"])
                   & (table["position"] <= cdr3["end"])]
    assert set(inside["region"]) == {"CDR3"}
    assert (table["region"] == "").any(), "the frameworks carry no region"


def test_the_profile_figure_follows_the_table_it_is_given(monkeypatch):
    table = pd.DataFrame({"name": ["Nb1"] * 3, "position": [1, 2, 3],
                          "residue": list("QVQ"),
                          "hydropathy": [1.0, 2.0, 3.0],
                          "charge": [0.0, 0.0, 0.0],
                          "region": ["", "", ""], "window": [9, 9, 9]})
    monkeypatch.setattr(vc, "profile_table",
                        lambda seqs, anns=None, window=None: table)
    fig = vc.plot_profiles({"Nb1": SEQ}, None, 9)
    assert list(fig.axes[0].get_lines()[0].get_ydata()) == [1.0, 2.0, 3.0]


# --------------------------------------------------------------------------- #
# Chou-Fasman
# --------------------------------------------------------------------------- #

def test_the_chou_fasman_table_is_the_profile_reshaped():
    """A reshaping of the heuristic, not a second copy of it."""
    profile = vc.chou_fasman_profile(SEQ, 5)
    table = vc.chou_fasman_table(SEQ, 5)
    assert len(table) == len(SEQ)
    assert table["helix_propensity"].to_numpy() == pytest.approx(profile["helix"])
    assert table["sheet_propensity"].to_numpy() == pytest.approx(profile["sheet"])
    assert list(table["state"]) == list(profile["state"])
    counts = table["state"].value_counts().to_dict()
    assert counts.get("E", 0) == profile["counts"]["E"]


# --------------------------------------------------------------------------- #
# Liability profile
# --------------------------------------------------------------------------- #

def test_the_liability_matrix_rows_sum_to_the_sequences_own_total():
    liabs = {n: vc.find_liabilities(s, None) for n, s in SEQS.items()}
    matrix = vc.liability_profile_matrix(liabs)
    for name, frame in liabs.items():
        assert matrix.loc[name].sum() == pytest.approx(float(frame["score"].sum()))


def test_the_liability_matrix_keeps_the_order_it_was_given():
    liabs = {n: vc.find_liabilities(s, None) for n, s in SEQS.items()}
    assert list(vc.liability_profile_matrix(liabs).index) == list(SEQS)


def test_a_set_with_no_liability_at_all_has_no_matrix_and_no_figure():
    empty = {"clean": pd.DataFrame()}
    assert vc.liability_profile_matrix(empty) is None
    assert vc.plot_liability_profile(empty) is None


def test_the_liability_figure_follows_the_matrix_it_is_given(monkeypatch):
    matrix = pd.DataFrame({"Met oxidation": [3.0, 7.0]}, index=["Nb1", "Nb2"])
    monkeypatch.setattr(vc, "liability_profile_matrix", lambda liabs: matrix)
    fig = vc.plot_liability_profile({"Nb1": pd.DataFrame(),
                                     "Nb2": pd.DataFrame()})
    heights = [patch.get_height() for patch in fig.axes[0].patches]
    assert heights == [3.0, 7.0]


# --------------------------------------------------------------------------- #
# The dendrogram
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def families() -> dict:
    frame = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
             "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA{cdr3}WGQGTQVTVSS")
    seqs = {f"Nb{i}": frame.format(cdr3="DRGYSSSWYD" + last)
            for i, last in enumerate("YFWL")}
    anns = {n: vc.annotate(s)[0] for n, s in seqs.items()}
    return vc.cluster_families(seqs, anns)


def test_the_merge_table_has_one_row_per_join(families):
    table = vc.linkage_table(families)
    assert len(table) == len(families["identity"]) - 1
    assert list(table["step"]) == list(range(1, len(table) + 1))


def test_every_merge_names_what_joined(families):
    table = vc.linkage_table(families)
    leaves = set(families["identity"].index)
    groups = {f"group_{s}" for s in table["step"]}
    for row in table.itertuples():
        assert row.joined_a in leaves | groups
        assert row.joined_b in leaves | groups


def test_the_merge_identity_is_the_axis_the_dendrogram_draws(families):
    """The figure's x axis is 100 minus the linkage distance, and so is this."""
    table = vc.linkage_table(families)
    distances = [row[2] for row in families["linkage"]]
    assert table["identity_pct"].tolist() == \
        [pytest.approx(100.0 - d * 100.0, abs=0.01) for d in distances]
    assert (table["identity_pct"] <= 100.0).all()


def test_a_set_with_nothing_to_cluster_gives_an_empty_merge_table():
    empty = vc.cluster_families({"only": SEQ}, {"only": vc.annotate(SEQ)[0]})
    table = vc.linkage_table(empty)
    assert not len(table)
    assert list(table.columns) == ["step", "joined_a", "joined_b",
                                   "identity_pct", "n_sequences"]


# --------------------------------------------------------------------------- #
# Sanger quality traces
# --------------------------------------------------------------------------- #

def _read(name: str, quals: list[int], start: int, end: int) -> dict:
    """A read as `sanger.read_ab1` returns it, with only what the figure uses.

    Built here rather than read from an .ab1 so that the trimmed window is a
    known quantity: the point of the test is which bases the table marks as
    kept, and a real chromatogram would make that depend on the trimmer.
    """
    return {"read": name, "quals": quals, "trim_start": start, "trim_end": end,
            "trimmed_len": end - start, "mean_phred_trim": 40.0}


def test_the_quality_table_has_one_row_per_base():
    reads = [_read("r1", [10, 40, 45, 42, 12], 1, 4),
             _read("r2", [20, 30], 0, 2)]
    table = sg.quality_table(reads)
    assert len(table) == 7
    assert list(table.columns) == ["read", "position", "phred", "kept"]
    assert table[table["read"] == "r1"]["position"].tolist() == [1, 2, 3, 4, 5]


def test_kept_marks_the_window_the_figure_draws_in_blue():
    table = sg.quality_table([_read("r1", [10, 40, 45, 42, 12], 1, 4)])
    assert table["kept"].tolist() == [False, True, True, True, False]


def test_the_table_covers_the_reads_the_figure_draws_and_no_more():
    """Both are capped at `max_reads`, and with the same number."""
    reads = [_read(f"r{i}", [30, 30], 0, 2) for i in range(20)]
    table = sg.quality_table(reads, max_reads=3)
    assert set(table["read"]) == {"r0", "r1", "r2"}
    fig = sg.plot_quality(reads, max_reads=3)
    assert len(fig.axes) == 3


def test_no_reads_gives_an_empty_table_and_no_figure():
    table = sg.quality_table([])
    assert not len(table)
    assert list(table.columns) == ["read", "position", "phred", "kept"]
    assert sg.plot_quality([]) is None
