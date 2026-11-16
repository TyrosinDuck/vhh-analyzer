"""
What `ranking.py` promises: ranks within a set, gates that are not points off,
and a Pareto front that is actually a Pareto front.

The tables here are hand-written rather than computed from sequences, which is
deliberate: these tests are about the comparison arithmetic, and a fixture built
by running the analysis would make a failure ambiguous between the two. The
tests that tie the axes back to the real triage table do use the analysis, and
are the ones that catch an axis naming a column that no longer exists.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import ranking as rk
import triage as tg
import vhh_core as vc

# Four sequences with known relationships: B and C tie on the liability score,
# C is the best aggregation window, D is worst at everything, and C has no
# polyreactivity percentile at all.
TABLE = pd.DataFrame({
    "name": ["A", "B", "C", "D"],
    "liability_score": [10.0, 5.0, 5.0, 20.0],
    "n_high_risk": [1, 0, 0, 2],
    "worst_window": [-0.2, 0.1, 0.3, -0.9],
    "humanness": [70.0, 75.0, 80.0, 60.0],
    "polyreactivity_pct": [50.0, 60.0, np.nan, 95.0],
    "has_orphan_cys": [False, True, False, False],
    "has_high_risk": [True, False, False, True],
    "has_apr": [False, False, False, True],
    "has_unusual_residue": [False, False, False, False],
    "liab_in_CDR": [False, False, False, False],
    "polyreactive": [False, False, False, True],
})


# --------------------------------------------------------------------------- #
# The axes, against the table they claim to read
# --------------------------------------------------------------------------- #

def test_every_axis_names_a_column_the_triage_table_produces():
    """An axis naming a column that moved would draw an empty stripe.

    The triage table is the single source of what the app computes per
    sequence, so the axes are checked against its column groups rather than
    against a copy of them here.
    """
    produced = {c for group in tg.COLUMN_GROUPS.values() for c in group}
    missing = [c for c in rk.AXES if c not in produced]
    assert not missing, f"axes naming no triage column: {missing}"


def test_the_defaults_are_axes_and_gates_that_exist():
    assert all(c in rk.AXES for c in rk.DEFAULT_AXES)
    assert all(g in tg.FLAG_COLUMNS for g in rk.DEFAULT_GATES)


def test_every_gate_label_names_a_flag_the_triage_table_raises():
    """The gates come from `triage.FLAG_COLUMNS`; only the wording is local."""
    unknown = [g for g in rk.GATE_LABELS if g not in tg.FLAG_COLUMNS]
    assert not unknown, f"labels for flags that do not exist: {unknown}"


def test_every_axis_declares_a_direction_and_a_reason():
    for column, axis in rk.AXES.items():
        assert axis.better in ("high", "low"), column
        assert axis.basis in ("derived", "theoretical"), (
            f"{column} claims basis {axis.basis!r}; this section never sees a "
            f"measurement, and the app's words for a result are measured, "
            f"derived and theoretical")
        assert len(axis.why) > 40, f"{column} has no stated reason"


def test_no_axis_is_offered_for_a_column_with_no_better_end():
    """pI, mass and CDR3 length are descriptions, not scores.

    They are in the triage table and must stay out of the axes: ranking a set
    by pI would read as an order of merit where none exists.
    """
    for column in ("pI", "mw_kda", "cdr3_len", "length", "n_cys",
                   "hallmark_score"):
        assert column not in rk.AXES


# --------------------------------------------------------------------------- #
# Ranks
# --------------------------------------------------------------------------- #

def test_a_low_is_better_axis_ranks_the_smallest_value_best():
    ranks = rk.within_set_ranks(TABLE, ["liability_score"])
    assert ranks.loc["D", "liability_score"] < ranks.loc["A", "liability_score"]
    assert ranks.loc["B", "liability_score"] == ranks.loc["C", "liability_score"]
    assert ranks["liability_score"].max() == pytest.approx(87.5)


def test_a_high_is_better_axis_ranks_the_largest_value_best():
    ranks = rk.within_set_ranks(TABLE, ["worst_window", "humanness"])
    assert ranks.loc["C", "worst_window"] == pytest.approx(100.0)
    assert ranks.loc["D", "worst_window"] == pytest.approx(25.0)
    assert ranks.loc["C", "humanness"] == pytest.approx(100.0)


def test_a_missing_value_is_not_ranked_last_it_is_not_ranked():
    """C has no percentile. Ranking it against those that do invents an order."""
    ranks = rk.within_set_ranks(TABLE, ["polyreactivity_pct"])
    assert pd.isna(ranks.loc["C", "polyreactivity_pct"])
    assert ranks.loc["A", "polyreactivity_pct"] > ranks.loc["D", "polyreactivity_pct"]


def test_ranks_are_relative_to_the_set_and_move_when_it_changes():
    """The claim the interface makes, checked rather than asserted in prose."""
    small = TABLE[TABLE["name"].isin(["A", "D"])]
    assert rk.within_set_ranks(small, ["liability_score"]).loc["A",
                                                               "liability_score"] \
        == pytest.approx(100.0)
    assert rk.within_set_ranks(TABLE, ["liability_score"]).loc["A",
                                                               "liability_score"] \
        == pytest.approx(50.0)


def test_an_axis_empty_for_every_sequence_is_dropped_rather_than_drawn():
    table = TABLE.assign(polyreactivity_pct=[np.nan] * 4)
    layers = rk.comparison_frame(table, ["liability_score",
                                         "polyreactivity_pct"])
    assert layers["axes"] == ["liability_score"]


# --------------------------------------------------------------------------- #
# Gates
# --------------------------------------------------------------------------- #

def test_a_gate_marks_a_sequence_without_touching_its_ranks():
    """The whole point of a gate: it is not a score.

    A sequence that raises a gate keeps the ranks it would have had, so the
    figure can show both at once and neither can hide the other.
    """
    layers = rk.comparison_frame(TABLE, ["liability_score"],
                                 ["has_orphan_cys"])
    assert bool(layers["gates"].loc["B", "blocked"])
    alone = rk.within_set_ranks(TABLE, ["liability_score"])
    assert layers["ranks"].loc["B", "liability_score"] == \
        alone.loc["B", "liability_score"]


def test_the_reason_names_the_gate_that_was_raised():
    gates = rk.gate_status(TABLE, ["has_orphan_cys", "has_high_risk"])
    assert gates.loc["B", "reasons"] == rk.GATE_LABELS["has_orphan_cys"]
    assert gates.loc["A", "reasons"] == rk.GATE_LABELS["has_high_risk"]
    assert gates.loc["C", "reasons"] == ""


def test_selecting_no_gate_blocks_nothing():
    gates = rk.gate_status(TABLE, [])
    assert not gates["blocked"].any()


# --------------------------------------------------------------------------- #
# The front
# --------------------------------------------------------------------------- #

def test_the_front_holds_what_nothing_else_beats_on_every_axis():
    front = rk.pareto_front(TABLE, ["liability_score", "worst_window"])
    by_name = front.set_index("name")
    # C is best on both, so it dominates everything and is alone on the front.
    assert bool(by_name.loc["C", "front"])
    assert not by_name.loc[["A", "B", "D"], "front"].any()
    assert int(by_name.loc["D", "dominated_by"]) == 3


def test_a_sequence_better_on_one_axis_stays_on_the_front():
    """Two axes that disagree put both sequences on the front.

    This is the case the whole module exists for: nothing says how much a
    liability is worth against an aggregation window, so neither sequence can
    be declared the loser.
    """
    table = pd.DataFrame({"name": ["X", "Y"],
                          "liability_score": [1.0, 9.0],
                          "worst_window": [-0.9, 0.5]})
    front = rk.pareto_front(table, ["liability_score", "worst_window"])
    assert front["front"].all()
    assert (front["dominated_by"] == 0).all()


def test_a_sequence_missing_an_axis_is_reported_rather_than_ranked_or_dropped():
    front = rk.pareto_front(TABLE, ["liability_score", "polyreactivity_pct"])
    by_name = front.set_index("name")
    assert not bool(by_name.loc["C", "comparable"])
    assert not bool(by_name.loc["C", "front"])
    assert len(front) == len(TABLE), "a sequence must not vanish from the table"


def test_excluding_a_sequence_stops_it_keeping_another_off_the_front():
    table = pd.DataFrame({"name": ["X", "Y"],
                          "liability_score": [1.0, 9.0],
                          "worst_window": [0.5, 0.1]})
    assert not bool(rk.pareto_front(table, ["liability_score", "worst_window"])
                    .set_index("name").loc["Y", "front"])
    excluded = rk.pareto_front(table, ["liability_score", "worst_window"],
                               exclude=["X"]).set_index("name")
    assert bool(excluded.loc["Y", "front"])
    assert not bool(excluded.loc["X", "front"])


def test_the_front_is_computed_on_values_not_on_ranks():
    """Ranks are a monotone transform, so the two agree; the test says which.

    If the implementation ever ranked first and compared after, a tie on an
    axis would become a tie in rank and the domination test would still hold,
    so this checks the sign handling instead: an axis whose good end is low
    must not invert the result.
    """
    table = pd.DataFrame({"name": ["X", "Y"], "liability_score": [1.0, 2.0],
                          "n_high_risk": [0, 1]})
    front = rk.pareto_front(table, ["liability_score", "n_high_risk"])
    assert bool(front.set_index("name").loc["X", "front"])
    assert not bool(front.set_index("name").loc["Y", "front"])


def test_no_layer_produces_a_total_score():
    """The refusal, checked.

    Nothing in what the app hands to the figures may be a composite: a column
    that summed the axes would be read as a grade, which is the thing
    `triage.py` refuses in writing and this module repeats.
    """
    layers = rk.comparison_frame(TABLE)
    for frame in (layers["ranks"], layers["values"]):
        for column in frame.columns:
            assert column in rk.AXES
    assert set(layers["pareto"].columns) == {"name", "excluded", "comparable",
                                             "dominated_by", "front"}
    assert "score" not in layers["gates"].columns
    assert "total" not in layers["gates"].columns


# --------------------------------------------------------------------------- #
# Against a real run
# --------------------------------------------------------------------------- #

FRAME = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
         "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA{cdr3}WGQGTQVTVSS")


@pytest.fixture(scope="module")
def real_triage() -> pd.DataFrame:
    seqs = {f"Nb{i}": FRAME.format(cdr3="DRGYSSSWYD" + last)
            for i, last in enumerate("YFWLK")}
    anns = {n: vc.annotate(s)[0] for n, s in seqs.items()}
    kinds = {n: vc.classify(s) for n, s in seqs.items()}
    return tg.triage_table(seqs, anns, kinds)


def test_the_default_axes_all_carry_values_on_a_real_set(real_triage):
    """A default axis that is always empty would show a grey stripe on open."""
    layers = rk.comparison_frame(real_triage)
    empty = [c for c in rk.DEFAULT_AXES if c not in layers["axes"]]
    assert not empty, f"default axes empty on a plain VHH set: {empty}"


def test_the_figures_draw_for_a_real_set(real_triage):
    layers = rk.comparison_frame(real_triage)
    assert rk.plot_rank_heatmap(layers["ranks"], layers["values"],
                                layers["gates"]) is not None
    axes = layers["axes"]
    assert rk.plot_pareto(layers["values"], axes[0], axes[1],
                          layers["pareto"], layers["gates"]) is not None


def test_the_heatmap_draws_something_when_no_axis_survives():
    """An empty figure beats a traceback on a set of non-antibody sequences."""
    empty = pd.DataFrame(index=pd.Index([], name="name"))
    assert rk.plot_rank_heatmap(empty) is not None
