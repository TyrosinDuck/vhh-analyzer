"""
The SAbDab2 snapshot: that it is what it claims, and that it says how old.

Two kinds of test here, and the second is the one that matters.

The first kind checks the index itself: that every shipped row is readable,
that the sequences are sequences this app can compute with, and that the
CDR columns mean what the rest of the app means by CDR. That last one is not
obvious and was verified rather than assumed, see below.

The second kind checks that the app cannot show a result from this file
without saying when the file was taken. A snapshot of a weekly-updated
database is only honest while its date travels with it, and the easiest way
to lose that is for somebody to simplify the caption.
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import pytest

import sabdab
import vhh_core as vc

ROOT = Path(__file__).resolve().parent.parent

needs_anarci = pytest.mark.skipif(
    not vc.anarci_available(),
    reason="ANARCI (and its HMMER binary) is not installed here")

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


@pytest.fixture(scope="module")
def index():
    return sabdab.load_index()


@pytest.fixture(scope="module")
def cablys3():
    import library as lib
    table, problems = lib.load_library(ROOT / "data" / "library.csv")
    assert not problems
    hits = table[table["name"].str.contains("cAbLys3", case=False)]
    assert len(hits) == 1
    return hits.iloc[0]["protein"]


# --------------------------------------------------------------------------- #
# The index is what it says it is
# --------------------------------------------------------------------------- #

def test_the_index_and_its_metadata_agree(index):
    """The row count in meta.json is the row count in the file.

    They are written by the same script in the same run, so a disagreement
    means somebody edited one of the two by hand, which is the thing this
    directory must never contain.
    """
    meta = sabdab.meta()
    assert meta, "data/sabdab/meta.json is missing"
    assert meta["sequences"] == len(index)
    assert meta["built_by"] == "tools/sabdab/build_index.py"


def test_every_sequence_is_one_this_app_can_compute_with(index):
    """No `X`, no empty rows, nothing the sequence cleaner would refuse.

    `vhh_core.clean_sequence` raises on a non-standard residue, so a single
    such row in the index would turn the whole search into a traceback on
    whichever query reached it. The build script drops them; this is the
    check that it did.
    """
    for sequence in index["sequence"]:
        assert sequence, "an empty sequence reached the index"
        assert not set(sequence) - STANDARD, \
            f"{sequence[:20]}... carries a residue this app cannot handle"


def test_sequences_are_unique(index):
    """One row per distinct sequence is the whole point of the derivation.

    The export has one row per chain instance, about 4700 of them for about
    1400 sequences. If duplicates survive, the hit list fills with the same
    nanobody repeated and the "next best" column stops meaning anything.
    """
    assert index["sequence"].is_unique


@needs_anarci
def test_the_shipped_cdrs_are_the_same_imgt_cdrs_this_app_annotates(index):
    """Checked, not assumed: SAbDab2's CDR columns match ANARCI's windows.

    The whole CDR3-identity column rests on this. SAbDab2 does not have to
    use the same boundaries the app does, and if it did not, a CDR3 identity
    computed between their column and our annotation would be comparing two
    different substrings and reporting a plausible-looking number.

    Twenty-five entries were compared when this was written and all three
    CDRs agreed on every one of them. A sample is checked here rather than
    the whole index, because this costs an ANARCI run per sequence.

    Entries ANARCI refuses to number are skipped rather than failed. A few in
    the index are truncated or divergent enough that it gives up, and that is
    a fact about those depositions, not a disagreement about boundaries. The
    skips are counted, so a change that made ANARCI refuse *most* of the
    index would still fail here instead of passing on an empty comparison.
    """
    sample = index[index["type"] == "SD-H"].head(25)
    compared = skipped = 0
    for entry in sample.itertuples():
        try:
            annotation = vc.annotate_anarci(entry.sequence)
        except Exception:                                 # noqa: BLE001
            skipped += 1
            continue
        compared += 1
        assert annotation["CDR1"]["seq"] == entry.cdr1, entry.structures
        assert annotation["CDR2"]["seq"] == entry.cdr2, entry.structures
        assert annotation["CDR3"]["seq"] == entry.cdr3, entry.structures
    assert compared >= 20, \
        f"only {compared} of {len(sample)} entries could be numbered " \
        f"({skipped} skipped); the comparison is no longer meaningful"


def test_vnar_entries_are_kept_and_labelled(index):
    """A shark domain is not a VHH, and hiding it would be worse than saying so.

    They stay in the file with their type, and the search leaves them out by
    default. Dropping them at build time would mean a user whose closest
    relative really is a VNAR gets a more distant VHH presented as the
    nearest thing, with nothing saying otherwise.
    """
    assert set(index["type"]) <= {"SD-H", "VNAR"}
    assert (index["type"] == "VNAR").any()


def test_the_licence_travels_with_the_data():
    """CC-BY 4.0, not the repository's MIT, and the directory has to say so."""
    licence = (ROOT / "data" / "sabdab" / "LICENSE")
    assert licence.exists(), "the CC-BY licence file is missing"
    text = licence.read_text(encoding="utf-8")
    assert "CC-BY 4.0" in text
    assert "MIT" in text, "the file must say what the rest of the repo is"
    assert "CC-BY 4.0" in sabdab.meta()["licence"]


# --------------------------------------------------------------------------- #
# The snapshot cannot be shown without its date
# --------------------------------------------------------------------------- #

def test_the_extraction_date_is_a_real_date():
    """Parsed, not just present: an unparseable stamp reads as "unknown"."""
    stamp = sabdab.meta()["extracted"]
    assert dt.date.fromisoformat(stamp) <= dt.date.today(), \
        "the snapshot claims to have been extracted in the future"


def test_age_is_counted_from_the_file_not_from_the_clock():
    """`age_days` takes today as an argument so the tests need not wait.

    It also means nothing in the app can read the clock twice inside one
    render and show two different ages.
    """
    extracted = dt.date.fromisoformat(sabdab.meta()["extracted"])
    assert sabdab.age_days(extracted) == 0
    assert sabdab.age_days(extracted + dt.timedelta(days=37)) == 37


def test_an_old_snapshot_says_so_in_words():
    """Past the threshold the message changes, and names the live database.

    This is the one sentence the screen prints verbatim. A user looking at a
    hit list has no other way of knowing that the file is a year old.
    """
    extracted = dt.date.fromisoformat(sabdab.meta()["extracted"])
    fresh = sabdab.staleness(extracted + dt.timedelta(days=1))
    assert not fresh["stale"]

    old = sabdab.staleness(extracted + dt.timedelta(days=sabdab.STALE_DAYS))
    assert old["stale"]
    assert "updates" in old["message"]
    assert str(sabdab.STALE_DAYS // 7) in old["message"]


def test_a_missing_date_is_treated_as_old_not_as_fresh(monkeypatch):
    """No date means no claim of freshness, which is the safe default.

    The opposite, defaulting to "recent" when the stamp is unreadable, would
    silence the warning in exactly the case where it is most needed.
    """
    monkeypatch.setattr(sabdab, "meta", lambda: {})
    state = sabdab.staleness()
    assert state["stale"] and state["days"] is None
    assert "no way to tell" in state["message"]


# --------------------------------------------------------------------------- #
# The search
# --------------------------------------------------------------------------- #

@needs_anarci
def test_cablys3_finds_its_own_crystal_structure(cablys3):
    """End to end, against a case where the right answer is known.

    cAbLys3 is the anti-lysozyme VHH of PDB 1MEL. If the search works, its
    own entry is the top hit at 100% identity, bound to lysozyme. This is
    the test that would catch a mangled index, a column read in the wrong
    order, or an identity computed against the wrong field.
    """
    annotation, _ = vc.annotate(cablys3, method="anarci")
    hits = sabdab.closest(cablys3, top=5, annotation=annotation)
    assert len(hits) == 5
    best = hits.iloc[0]
    assert best["identity"] == 100.0
    assert "1mel" in best["structures"]
    assert "ysozyme" in best["antigen_names"]
    assert best["cdr3_identity"] == 100.0


def test_hits_come_back_in_order(cablys3):
    """Best first, because every caption says "closest" and means the top row."""
    hits = sabdab.closest(cablys3, top=10)
    identities = list(hits["identity"])
    assert identities == sorted(identities, reverse=True)


def test_vnar_entries_are_out_unless_asked_for(cablys3):
    """Default off. A VNAR hit answers a question about a different molecule."""
    without = sabdab.closest(cablys3, top=50)
    assert "VNAR" not in set(without["type"])
    with_them = sabdab.closest(cablys3, top=50, include_vnar=True)
    assert len(with_them) == len(without)        # same count, different members


def test_cdr3_identity_is_empty_rather_than_guessed(cablys3):
    """Without an annotation the column is None, never the full-length number.

    Filling it with the overall identity would be the worst available
    answer: the column exists precisely because the two differ.
    """
    hits = sabdab.closest(cablys3, top=3)
    assert hits["cdr3_identity"].isna().all()


def test_summary_reports_the_gap_to_the_second_hit(cablys3):
    """A top hit clear of the field and one in a cluster are different answers."""
    hits = sabdab.closest(cablys3, top=3)
    found = sabdab.summary(hits)
    assert found["n"] == 3
    assert found["spread"] == pytest.approx(
        hits.iloc[0]["identity"] - hits.iloc[1]["identity"], abs=0.05)
    assert sabdab.summary(hits.head(0))["best"] is None


def test_the_build_script_refuses_an_export_it_does_not_recognise(tmp_path):
    """SAbDab2 is live and its export format can change under us.

    A renamed column must stop the build with a message naming it, not
    produce an index with a silently empty field. The index is rebuilt by
    hand every few months, which is exactly the interval over which nobody
    remembers what the columns were called last time.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "build_index", ROOT / "tools" / "sabdab" / "build_index.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    export = tmp_path / "wrong.csv"
    export.write_text("PDB,VH,CDR-H1\nabcd,EVQL,GRTF\n", encoding="utf-8")
    with pytest.raises(SystemExit) as raised:
        module.build(export)
    assert "missing columns" in str(raised.value)


def test_meta_json_is_valid_json_and_names_its_source():
    """It is read by the screen and by this suite; a typo there breaks both."""
    meta = json.loads((ROOT / "data" / "sabdab" / "meta.json")
                      .read_text(encoding="utf-8"))
    assert meta["source"].startswith("https://")
    assert meta["export_rows"] > meta["sequences"], \
        "the export has one row per chain instance, so it must be larger"
