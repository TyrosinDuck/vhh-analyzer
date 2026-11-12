"""
TNP tests: the thresholds, and the one way this can give a wrong answer.

The metric values themselves are counts and a ratio, and neither is hard to
get right. What is easy to get wrong is the *comparison*: TNP's bands are
calibrated on IMGT CDR boundaries, and this app can also annotate with regex
anchors, which give different lengths on the same sequence. A flag drawn
against the wrong boundary definition looks exactly like a correct one. Most
of what is below is about that.

The thresholds are read from `tnp.METRICS` rather than written out again
here, except in the handful of cases where a literal from the paper is the
whole point of the test. Those are marked.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import tnp
import vhh_core as vc

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def cablys3():
    """A real camelid VHH with a long CDR3, from the app's own library."""
    import library as lib
    table, problems = lib.load_library(ROOT / "data" / "library.csv")
    assert not problems
    hits = table[table["name"].str.contains("cAbLys3", case=False)]
    assert len(hits) == 1
    return hits.iloc[0]["protein"]


needs_anarci = pytest.mark.skipif(
    not vc.anarci_available(),
    reason="ANARCI (and its HMMER binary) is not installed here")


# --------------------------------------------------------------------------- #
# Table 2 of Gordon et al. 2026, transcribed once more, on purpose
# --------------------------------------------------------------------------- #

def test_the_published_thresholds_are_what_the_module_carries():
    """A second transcription of Table 2, which is the only way to check one.

    Everywhere else in this app a second list of the same thing is the bug.
    Here it is the test: the values in `tnp.METRICS` came off a published PDF
    by hand, and the only error that catches is reading the table twice and
    getting the same answer. If this test and the module disagree, go back to
    the paper rather than editing either one to match the other.
    """
    published = {
        "total_cdr_length": ((20, 39), ((20, 24), (37, 39))),
        "cdr3_length": ((5, 23), ((5, 8), (21, 23))),
        "cdr3_compactness": ((0.56, 1.61), ((0.56, 0.81), (1.57, 1.61))),
        "psh": ((73.40, 155.47), ((73.40, 79.59), (126.83, 155.47))),
        "ppc": ((0.00, 1.18), ((0.39, 1.18),)),
        "pnc": ((0.00, 1.88), ((1.47, 1.88),)),
    }
    assert set(published) == set(tnp.METRICS), \
        "TNP has six metrics; the module lists a different set"
    for metric, (clinical, amber) in published.items():
        assert tnp.METRICS[metric]["clinical"] == clinical, metric
        assert tnp.METRICS[metric]["amber"] == amber, metric


@pytest.mark.parametrize("value, expected", [
    (19, "red"),        # below the clinical range
    (20, "amber"),      # the bottom edge is inside it, and in the amber band
    (24, "amber"),      # top of the lower amber band, inclusive
    (25, "green"),      # first value clear of both bands
    (36, "green"),
    (37, "amber"),      # bottom of the upper amber band, inclusive
    (39, "amber"),      # top edge of the clinical range
    (40, "red"),
])
def test_total_cdr_length_flags_at_the_published_edges(value, expected):
    """Band edges are inclusive on both sides, the cautious reading.

    Table 2 writes the amber bands as closed intervals and the red regions as
    "above or below", which leaves the edge itself ambiguous. A value sitting
    exactly on 24 is called amber rather than green here: of the two readings
    that is the one that cannot turn a borderline candidate green.
    """
    assert tnp.flag("total_cdr_length", value) == expected


def test_a_one_sided_metric_has_no_lower_amber_band():
    """PPC and PNC start at zero, so only the top tail is flagged.

    Writing them with two bands, by symmetry with the length metrics, would
    paint every well-behaved sequence amber: a patch score of zero is the
    *good* end, not an outlier.
    """
    assert tnp.flag("ppc", 0.0) == "green"
    assert tnp.flag("ppc", 0.38) == "green"
    assert tnp.flag("ppc", 0.39) == "amber"
    assert tnp.flag("ppc", 1.19) == "red"


def test_no_value_is_not_a_flag():
    """A missing metric gives None, never green.

    Four of the six metrics are usually absent, and a default of green would
    read as "checked and fine" for metrics nobody computed.
    """
    assert tnp.flag("cdr3_compactness", None) is None
    assert tnp.flag("cdr3_compactness", float("nan")) is None
    with pytest.raises(KeyError):
        tnp.flag("not_a_metric", 1.0)


# --------------------------------------------------------------------------- #
# The wrong-numbering trap
# --------------------------------------------------------------------------- #

@needs_anarci
def test_imgt_and_regex_disagree_and_only_imgt_is_flagged(cablys3):
    """The bug this guard exists for, demonstrated on a real VHH.

    The two annotators give different CDR lengths for the same sequence,
    because IMGT windows and conserved-anchor boundaries are different
    definitions, not two estimates of one thing. TNP's bands are IMGT, so a
    regex profile reports its lengths and carries no flags at all.
    """
    imgt_ann, imgt_method = vc.annotate(cablys3, method="anarci")
    regex_ann = vc.annotate_regex(cablys3)

    imgt = tnp.profile("x", cablys3, imgt_ann, imgt_method)
    regex = tnp.profile("x", cablys3, regex_ann, "Regex anchors")

    assert imgt["scheme_ok"] and not regex["scheme_ok"]
    assert imgt["total_cdr_length"] != regex["total_cdr_length"], \
        "if the two agreed on this sequence, pick one where they do not"
    assert imgt["total_cdr_length_flag"] is not None
    assert regex["total_cdr_length_flag"] is None
    assert regex["n_flagged_metrics"] == 0


@needs_anarci
def test_a_long_cdr3_vhh_falls_outside_the_clinical_range(cablys3):
    """cAbLys3 is red on both length metrics, and that is the right answer.

    Its IMGT CDR3 is 26 residues against a clinical-stage range that stops at
    23, which puts the total out of range too. Worth pinning for two reasons.
    It is the behaviour the tab exists to produce: a red flag is not a defect,
    it says this molecule does something no trialled nanobody does, and
    cAbLys3 is a crystallography workhorse, not a drug candidate. And the
    numbers themselves catch an off-by-one in the length counting more cheaply
    than reading the code does: CDR1 and CDR2 are the canonical 8 and 8, so
    any drift shows up here first.
    """
    profile = tnp.profile("cAbLys3", cablys3, *vc.annotate(cablys3,
                                                           method="anarci"))
    assert (profile["cdr1_length"], profile["cdr2_length"]) == (8, 8)
    assert profile["cdr3_length"] == 26
    assert profile["total_cdr_length"] == 42
    assert profile["cdr3_length_flag"] == "red"
    assert profile["total_cdr_length_flag"] == "red"
    assert profile["worst"] == "red"


@needs_anarci
def test_the_cdr_lengths_are_the_annotation_lengths(cablys3):
    """Total is the sum of the three, and each is its annotated region.

    Trivial, and worth pinning: the one thing a reader of the table assumes.
    """
    annotation, method = vc.annotate(cablys3, method="anarci")
    lengths = tnp.cdr_lengths(cablys3, annotation, method)
    for region in ("CDR1", "CDR2", "CDR3"):
        assert lengths[region.lower() + "_length"] == \
            len(annotation[region]["seq"])
    assert lengths["total_cdr_length"] == sum(
        lengths[f"cdr{i}_length"] for i in (1, 2, 3))


# --------------------------------------------------------------------------- #
# The screen reads these
# --------------------------------------------------------------------------- #

def test_the_guidelines_table_is_derived_from_the_thresholds():
    """One row per metric, built from `METRICS` and not typed out again."""
    table = tnp.metric_table()
    assert len(table) == len(tnp.METRICS)
    assert set(table["Metric"]) == {m["label"] for m in tnp.METRICS.values()}


def test_the_sequence_metrics_are_the_two_the_app_computes():
    """If a metric moves between sequence and structure, the screens change.

    `SEQUENCE_METRICS` is what the analysis screen flags without a model;
    anything else needs one or is not computed at all.
    """
    assert set(tnp.SEQUENCE_METRICS) == {"total_cdr_length", "cdr3_length"}
    assert "cdr3_compactness" in tnp.STRUCTURAL_METRICS


def test_a_profile_table_survives_a_sequence_it_cannot_annotate():
    """One bad sequence must not take the panel down with it.

    The rest of the app has this property everywhere; the screen builds this
    table from whatever the user pasted, so it needs it too.
    """
    table = tnp.profile_table({"nonsense": "XXXXXXXXXX"})
    assert len(table) == 1
