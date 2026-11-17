"""
What `specificity.py` claims, checked against the authors' own numbers.

The fixture, `fixtures/harvey_pool_sample.csv`, is 150 sequences taken from the
two pools released with

    Harvey EP, Shin JE, Skiba MA, et al. An in silico method to assess antibody
    fragment polyreactivity. Nat Commun 13, 7554 (2022). MIT licence.

by `tools/polyreactivity/build.py`'s reading of their tables: the full-length
sequence rebuilt from their IMGT position columns, their own CDR annotation,
and the decision function they publish for each clone. Picked at even spacing
along their score, so the sample spans the whole range rather than its middle.

That is what makes these tests worth running: the ground truth comes from
outside this repository, so they check the implementation rather than check it
against itself. The one test that cannot work this way is the ANARCI path,
which needs a binary the suite may not have, and it skips.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import specificity as sp
import triage as tg
import vhh_core as vc

FIXTURE = Path(__file__).parent / "fixtures" / "harvey_pool_sample.csv"


@pytest.fixture(scope="module")
def pool() -> pd.DataFrame:
    return pd.read_csv(FIXTURE)


# --------------------------------------------------------------------------- #
# The model
# --------------------------------------------------------------------------- #

def test_the_score_reproduces_the_authors_published_decision_function(pool):
    """The whole point of shipping their weights rather than our own.

    Their decision function is positive for *low* polyreactivity and this
    module negates it, so the comparison is against the negative. A failure
    here means the featurisation order, the L1 normalisation, the intercept or
    the sign changed, and every percentile in the app is wrong.
    """
    cdrs = (pool["cdr1_imgt"] + pool["cdr2_imgt"] + pool["cdr3_imgt"]).tolist()
    score = sp.polyreactivity(cdrs)
    assert np.max(np.abs(score + pool["decision_function"].to_numpy())) < 1e-9


def test_the_feature_space_is_the_one_the_weights_are_aligned_to():
    """20 monomers, then 400 dimers, then 8000 trimers, in that order."""
    assert sp.N_FEATURES == 8420
    matrix = sp.cdrs_to_kmer(["AAA"])
    assert matrix.shape == (1, 8420)
    # "AAA" is one trimer, so the only non-zero column is the first trimer one.
    assert matrix[0, 420] == pytest.approx(1.0)
    assert np.count_nonzero(matrix) == 1


def test_a_sequence_too_short_to_hold_a_trimer_scores_without_raising():
    assert sp.cdrs_to_kmer([""]).sum() == 0.0
    assert np.isfinite(sp.polyreactivity(["", "AC"])).all()


def test_the_percentile_is_bounded_and_rises_with_the_score():
    low, high = sp.polyreactivity_percentile([-6.0, 6.0])
    assert 0.0 <= low < high <= 100.0
    rising = sp.polyreactivity_percentile(np.linspace(-5, 5, 25))
    assert np.all(np.diff(rising) >= 0)


def test_the_reference_report_reads_the_files_rather_than_quoting_them():
    """The screens quote these, so they have to come from the shipped files."""
    report = sp.reference_report()
    assert report["n_features"] == sp.N_FEATURES
    assert 0 < report["n_weights"] < sp.N_FEATURES
    assert report["n_pool"] == sp.polyreactivity_pool_size() > 10_000
    assert report["plastic_modes"] == ["aromatic", "cationic"]
    assert report["n_peptides"] == len(sp.reference_peptides())


# --------------------------------------------------------------------------- #
# Which CDRs the model is fed
# --------------------------------------------------------------------------- #

def test_the_anchors_reproduce_the_authors_cdrs_whenever_they_accept(pool):
    """A rate, not an identity, and the rate is the claim the app displays.

    Writing this as "every accepted sequence is exact" would be a nicer test
    and a false one: the CDR2 and FR3 boundary is genuinely ambiguous for some
    sequences and the anchors get a few of them wrong. The measured figure over
    the authors' whole library is 99.2% of accepted sequences exact; this
    sample is small, so the bar is set lower than that and the point of the
    test is that it cannot quietly collapse.
    """
    accepted = exact = 0
    for row in pool.itertuples(index=False):
        got = sp.cdrs_regex(row.sequence)
        if not got.ok:
            continue
        accepted += 1
        exact += (got.cdr1 == row.cdr1_imgt and got.cdr2 == row.cdr2_imgt
                  and got.cdr3 == row.cdr3_imgt)
    assert accepted >= 0.4 * len(pool), "the anchors have stopped accepting"
    assert exact / accepted >= 0.95


def test_a_rejected_sequence_carries_its_reason_and_no_cdrs():
    got = sp.cdrs_regex("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
    assert not got.ok
    assert got.cdrs == ""
    assert got.note


def test_the_ambiguous_cdr2_boundary_is_refused_rather_than_guessed():
    """The one case where only the numbering can decide, see `cdrs_regex`."""
    ambiguous = ("QVQLVESGGGLVQAGGSLRLSCAASGITFGRNTMGWYRQAPGKERELVAAIDSGGNT"
                 "YADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAVDGGEFSSYGYWGQGTQVTVSS")
    got = sp.cdrs_regex(ambiguous)
    assert not got.ok
    assert "ambiguous" in got.note


def test_the_display_annotation_of_vhh_core_is_never_fed_to_the_model(pool):
    """`model_cdrs` must ignore an annotation it was not told is IMGT.

    This is the mistake the module exists to prevent: the regex annotation of
    `vhh_core` draws boundaries for display, and feeding them to the model
    moves the percentile by a median of 9 points. A caller that passes an
    annotation without the flag gets the anchors, not that annotation.
    """
    sequence = pool["sequence"].iloc[1]
    display = vc.annotate_regex(sequence)
    through_model = sp.model_cdrs(sequence, display, annotation_is_imgt=False)
    assert through_model.annotator == "regex anchors"
    assert through_model.cdrs == sp.cdrs_regex(sequence).cdrs
    display_cdrs = "".join((display.get(k) or {}).get("seq", "")
                           for k in ("CDR1", "CDR2", "CDR3"))
    assert through_model.cdrs != display_cdrs


def test_an_imgt_annotation_is_used_when_the_caller_says_it_is_one(pool):
    """And the model's CDR2 is rebuilt as IMGT 55 plus the app's CDR2.

    The annotation here is the app's shape: its CDR2 is the IMGT loop, 56 to
    65, and IMGT 55 is the last residue of FR2. What comes back must be the
    authors' CDR2, one residue wider.
    """
    row = pool.iloc[1]
    annotation = {"FR2": {"seq": "WYRQAPGKERELVA" + row["cdr2_imgt"][0]},
                  "CDR1": {"seq": row["cdr1_imgt"]},
                  "CDR2": {"seq": row["cdr2_imgt"][1:]},
                  "CDR3": {"seq": row["cdr3_imgt"]}}
    got = sp.model_cdrs(row["sequence"], annotation, annotation_is_imgt=True)
    assert got.ok and got.annotator == "ANARCI / IMGT"
    assert got.cdr2 == row["cdr2_imgt"]
    assert got.cdrs == row["cdr1_imgt"] + row["cdr2_imgt"] + row["cdr3_imgt"]


def test_an_annotation_without_fr2_is_refused_rather_than_scored_short(pool):
    """Without IMGT 55 the model's CDR2 cannot be assembled at all."""
    row = pool.iloc[1]
    got = sp.cdrs_from_annotation({"CDR1": {"seq": row["cdr1_imgt"]},
                                   "CDR2": {"seq": row["cdr2_imgt"][1:]},
                                   "CDR3": {"seq": row["cdr3_imgt"]}})
    assert not got.ok
    assert "IMGT 55" in got.note


def test_anarci_and_the_anchors_agree_where_anarci_runs(pool):
    """The two paths must give one number, or the same file scores two ways.

    Skipped without HMMER, which is the whole reason the anchors exist. On the
    fixture the two agree on every sequence both of them accept.
    """
    if not vc.anarci_available():
        pytest.skip("ANARCI cannot run here, so the two paths cannot be compared")
    agree = compared = 0
    for sequence in pool["sequence"]:
        anchors = sp.cdrs_regex(sequence)
        through_anarci = sp.cdrs_from_annotation(vc.annotate_anarci(sequence))
        if not (anchors.ok and through_anarci.ok):
            continue
        compared += 1
        agree += anchors.cdrs == through_anarci.cdrs
    assert compared > 50, "too few sequences annotated both ways to mean anything"
    assert agree == compared


def test_a_long_cdr1_is_scored_rather_than_called_implausible():
    """An IMGT numbering is not a guess, so its lengths are not checked.

    This VHH is PDB 7dss chain A, from the SAbDab snapshot this repository
    ships: ANARCI gives it a 13-residue CDR1, and the plausibility ceiling
    written for the regex anchors used to refuse it with a note blaming
    anchors that had never run. Ten of the 1 352 single domains in that
    snapshot have a CDR1 longer than 12, and every one of them is a solved
    structure.
    """
    if not vc.anarci_available():
        pytest.skip("ANARCI cannot run here")
    # Read out of the snapshot rather than pasted in here: a sequence copied
    # into a test is a second copy of data the repository already holds, and
    # the one time it was pasted it was pasted wrong.
    snapshot = pd.read_csv(Path(__file__).parent.parent / "data" / "sabdab"
                           / "single_domain.csv")
    row = snapshot[snapshot["structures"].astype(str).str.contains("7dss")]
    assert len(row) == 1, "7dss is no longer in the SAbDab snapshot"
    got = sp.cdrs_from_annotation(vc.annotate_anarci(str(row["sequence"].iloc[0])))
    assert got.ok, got.note
    assert len(got.cdr1) >= 13
    assert 0.0 <= sp.polyreactivity_percentile(sp.polyreactivity([got.cdrs]))[0] <= 100.0


def test_the_anarci_path_reproduces_the_authors_cdrs(pool):
    """Including their CDR2, which is one residue wider than the app's.

    The residual is ANARCI itself: their run and ours place the gaps inside
    CDR2 differently on a few percent of sequences, which no definition here
    can fix. 96% exact is what that looks like; a collapse to the 0% this
    test caught before `cdrs_from_annotation` widened the CDR2 would mean the
    whole ANARCI path had gone back to scoring the wrong input.
    """
    if not vc.anarci_available():
        pytest.skip("ANARCI cannot run here")
    exact = annotated = 0
    for row in pool.itertuples(index=False):
        got = sp.cdrs_from_annotation(vc.annotate_anarci(row.sequence))
        if not got.ok:
            continue
        annotated += 1
        exact += (got.cdr1 == row.cdr1_imgt and got.cdr2 == row.cdr2_imgt
                  and got.cdr3 == row.cdr3_imgt)
    assert annotated > 100
    assert exact / annotated >= 0.90


# --------------------------------------------------------------------------- #
# Plastic resemblance, and what it is not allowed to claim
# --------------------------------------------------------------------------- #

def test_every_plastic_row_says_it_is_not_calibrated():
    frame = sp.plastic_resemblance(["AAAWWWRRR", "GGGG", ""])
    assert list(frame["calibrated"]) == [False, False, False]


def test_the_two_modes_separate_the_two_published_classes():
    """The reason there are two scores: each class is extreme in its own.

    Taken from the shipped reference set itself, so the test moves if the
    peptide table does. A cationic tag scoring higher on the aromatic mode
    than a tryptophan-rich peptide would mean the weights or the features had
    lost their meaning.
    """
    peptides = sp.reference_peptides()
    cationic = peptides.loc[peptides["binding_class"] == "cationic", "peptide"]
    aromatic = peptides.loc[peptides["binding_class"] == "aromatic", "peptide"]
    scored = sp.plastic_resemblance(list(cationic) + list(aromatic))
    n = len(cationic)
    assert scored["cationic_mode"].iloc[:n].min() > scored["cationic_mode"].iloc[n:].max()
    assert scored["aromatic_mode"].iloc[n:].max() > scored["aromatic_mode"].iloc[:n].max()


def test_the_plastic_percentiles_are_bounded(pool):
    cdrs = (pool["cdr1_imgt"] + pool["cdr2_imgt"] + pool["cdr3_imgt"]).tolist()
    frame = sp.plastic_resemblance(cdrs)
    for column in ("cationic_mode_pct", "aromatic_mode_pct"):
        assert frame[column].between(0.0, 100.0).all()


def test_every_reference_peptide_names_a_source():
    """A sequence whose provenance is not in the file is an invented constant."""
    peptides = sp.reference_peptides()
    assert len(peptides) >= 10
    assert peptides["peptide"].str.fullmatch(f"[{sp.AA}]+").all()
    assert peptides["source"].str.len().min() > 30
    assert set(peptides["source_kind"]) <= {"primary", "review"}
    assert set(peptides["binding_class"]) == {"cationic", "aromatic"}


def test_nothing_in_the_module_adds_the_two_scores_together():
    """The composite index the app must not grow.

    Checked on the output rather than on the source: any column that mixes a
    calibrated score with an uncalibrated one would have to appear here.
    """
    frame = sp.table({"a": "QVQLVESGGGLVQAGGSLRLSCAASGITFGRNAMGWYRQAPGKERELVA"
                           "AISGRGGRTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA"
                           "RSRRLPGDYWGQGTQVTVSS"})
    for column in frame.columns:
        name = column.lower()
        assert not any(word in name for word in
                       ("total", "combined", "composite", "overall", "index")), \
            f"{column} looks like a composite of the two scores"


# --------------------------------------------------------------------------- #
# The table the screens draw
# --------------------------------------------------------------------------- #

def test_a_sequence_that_cannot_be_annotated_keeps_its_row_and_gets_no_score():
    """A dropped row is indistinguishable from a sequence never submitted."""
    frame = sp.table({"nonsense": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"})
    assert len(frame) == 1
    assert not frame["annotation_ok"].iloc[0]
    assert frame["note"].iloc[0]
    for column in ("polyreactivity_score", "polyreactivity_pct",
                   "cationic_mode_pct", "aromatic_mode_pct"):
        assert pd.isna(frame[column].iloc[0])
    assert not frame["polyreactive_flag"].iloc[0]


def test_the_flag_follows_the_registered_threshold(pool):
    frame = sp.table({row.id: row.sequence for row in pool.itertuples(index=False)})
    scored = frame[frame["annotation_ok"]]
    assert len(scored)
    expected = scored["polyreactivity_pct"] >= sp.FLAG_PERCENTILE
    assert list(scored["polyreactive_flag"]) == list(expected)


def test_the_triage_table_carries_the_percentile_and_the_flag(pool):
    sequences = {row.id: row.sequence for row in pool.head(5).itertuples(index=False)}
    # As the screen does: an annotation that cannot be produced is a missing
    # annotation, not a failed run. `vhh_core`'s anchors raise on a sequence
    # whose Cys23 they cannot place, and two of these five are such sequences.
    annotations = {}
    for name, seq in sequences.items():
        try:
            annotations[name] = vc.annotate(seq)[0]
        except ValueError:
            annotations[name] = None
    frame = tg.triage_table(sequences, annotations)
    assert "polyreactivity_pct" in frame.columns
    assert "polyreactive" in frame.columns
    assert "Specificity" in tg.COLUMN_GROUPS
    # The uncalibrated score has no business in a table of flags.
    assert not [c for c in frame.columns if "mode" in c]


def test_the_runtime_needs_nothing_but_numpy_and_pandas():
    """The property that keeps this deployable on a 690 MB container.

    scikit-learn opens the authors' pickle in the build script, where it
    belongs. An import of it here, or of torch, or of pickle, would mean the
    weights had stopped being plain numbers in an .npz.

    Stated as what may not be imported rather than as what may. The allow-list
    version of this test failed the day `observed_cdr_bounds` started reading
    the SAbDab snapshot, which is a sibling module of this repository carrying
    no dependency of its own beyond pandas, and the real property, no heavy
    machine-learning runtime, was never in question. So: the standard library
    and the repository's own modules are free, third-party imports are limited
    to numpy and pandas, and the named offenders are refused by name in case
    one of them arrives by a route this does not anticipate.
    """
    root = Path(__file__).parent.parent
    source = (root / "specificity.py").read_text()
    imported = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])

    banned = {"sklearn", "scikit_learn", "torch", "tensorflow", "pickle",
              "joblib", "transformers", "scipy"}
    assert not (imported & banned), imported & banned

    first_party = {path.stem for path in root.glob("*.py")}
    third_party = imported - set(sys.stdlib_module_names) - first_party
    assert third_party <= {"numpy", "pandas"}, third_party
