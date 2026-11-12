"""
Session files: the inputs of one analysis, written down and read back.

A session file exists so that a run can be repeated, which means the only
thing that matters about it is whether what comes back is what went in. These
tests are mostly about the ways it should refuse: a file that restores six of
seven settings puts a screen on the table that looks like the saved one and is
not, and nothing would say so.

The other half is the methods paragraph, where the thing worth testing is not
the prose but what it is not allowed to claim: it must name the annotator that
actually ran rather than the setting that was asked for, and it must say that
the mass and the pI are theoretical.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import session as ses
import vhh_core as vc

ROOT = Path(__file__).resolve().parent.parent

SEQ = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
       "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDWGQGTQVTVSS")


def a_record(**kwargs) -> dict:
    defaults = dict(sequences={"Nb1": SEQ}, ph=7.4, window=9, method="auto",
                    method_used="ANARCI / IMGT", version="1.1.0 (build x.y)",
                    release="1.1.0", build="20261112", code="abc1234",
                    runtime={"Python": "3.12.1"})
    defaults.update(kwargs)
    return ses.session_record(**defaults)


def test_a_record_round_trips():
    back = ses.session_inputs(a_record())
    assert back["sequences"] == {"Nb1": SEQ}
    assert back["ph"] == 7.4
    assert back["window"] == 9
    assert back["method"] == "auto"
    assert back["method_used"] == "ANARCI / IMGT"
    assert back["code"] == "abc1234"
    assert back["runtime"] == {"Python": "3.12.1"}


def test_a_record_survives_json():
    """Written and read as a file, which is how it is actually used."""
    text = json.dumps(a_record())
    assert ses.session_inputs(json.loads(text))["sequences"] == {"Nb1": SEQ}


def test_the_record_holds_no_computed_value():
    """The point of the format: it stores the question, not the answer.

    A stored molar mass would be reprinted a year later by a version that
    computes it differently, and nothing would show the difference. Keeping
    only the inputs means reopening a session recomputes, and
    ``created_with`` is there to say which version produced what you are
    looking at now.
    """
    record = a_record()
    assert set(record) == {"schema", "created_with", "created_utc", "release",
                           "build", "code", "runtime", "settings", "sequences"}
    flat = json.dumps(record)
    for computed in ("mw_da", "pI", "ext_coeff", "gravy", "net_charge"):
        assert computed not in flat


def test_sequences_are_stored_as_residues_only():
    record = a_record(sequences={"Nb1": "qvq lqe\nSGG", "empty": "  "})
    assert record["sequences"] == {"Nb1": "QVQLQESGG"}


def test_the_fasta_it_writes_is_read_back_by_the_app_parser():
    """The text area and the session file have to agree about formatting.

    `fasta` is what gets pushed into the analysis input, so what it writes has
    to survive `vc.parse_input`: a restored session that lands in the box as
    something the parser rejects would be a reload that fails at the next
    click rather than at the upload.
    """
    text = ses.fasta({"Nb1": SEQ, "Nb2": SEQ[:-3]})
    assert vc.parse_input(text) == {"Nb1": SEQ, "Nb2": SEQ[:-3]}


@pytest.mark.parametrize("schema, message", [
    ("", "carries no `schema`"),
    ("something-else/1", "this app reads"),
    ("vhh-toolbox-session/2", "not being read rather than read wrongly"),
])
def test_an_unknown_schema_is_refused(schema, message):
    record = a_record()
    if schema:
        record["schema"] = schema
    else:
        del record["schema"]
    with pytest.raises(ValueError, match=re.escape(message)):
        ses.session_inputs(record)


@pytest.mark.parametrize("key, value", [("ph", 1.0), ("ph", 14.0),
                                        ("window", 1), ("window", 99)])
def test_a_setting_outside_the_widget_is_refused(key, value):
    """A value the sidebar cannot show is not restorable.

    Clamping it would be worse than refusing: the screen would come back with
    a setting the file does not contain, and the numbers under it would not be
    the saved run's.
    """
    record = a_record()
    record["settings"][key] = value
    with pytest.raises(ValueError, match="outside the"):
        ses.session_inputs(record)


def test_a_missing_setting_is_refused():
    record = a_record()
    del record["settings"]["ph"]
    with pytest.raises(ValueError, match="missing from this session file"):
        ses.session_inputs(record)


def test_a_setting_that_is_not_a_number_is_refused():
    record = a_record()
    record["settings"]["ph"] = "neutral"
    with pytest.raises(ValueError, match="not a number"):
        ses.session_inputs(record)


def test_a_fractional_window_is_refused_rather_than_rounded():
    record = a_record()
    record["settings"]["window"] = 9.5
    with pytest.raises(ValueError, match="whole number"):
        ses.session_inputs(record)


def test_an_unknown_annotation_method_is_refused():
    record = a_record(method="auto")
    record["settings"]["annotation_method"] = "imgt-ish"
    with pytest.raises(ValueError, match="this app knows"):
        ses.session_inputs(record)


def test_a_sequence_that_is_not_text_is_refused():
    record = a_record()
    record["sequences"]["Nb2"] = 42
    with pytest.raises(ValueError, match="is not text"):
        ses.session_inputs(record)


def test_an_empty_session_is_allowed_and_comes_back_empty():
    """Saving before analysing anything is a legitimate thing to do.

    The screen decides whether an empty session is worth restoring; the format
    has no reason to refuse one.
    """
    assert ses.session_inputs(a_record(sequences={}))["sequences"] == {}


# --------------------------------------------------------------------------
# The methods paragraph
# --------------------------------------------------------------------------

def test_the_paragraph_names_the_annotator_that_ran_not_the_one_asked_for():
    """The one claim in the paragraph that a user cannot check by looking.

    "auto" resolves to ANARCI where it is installed and to the regex anchors
    where it is not, so a paragraph repeating the setting would describe IMGT
    numbering on a server that never ran it.
    """
    text = ses.methods_paragraph(a_record(method="auto",
                                          method_used="Regex anchors"))
    assert "Regex anchors" in text
    assert "automatic setting resolved to" in text
    assert "ANARCI" not in text


def test_the_paragraph_says_the_sequence_values_are_theoretical():
    text = ses.methods_paragraph(a_record())
    assert "theoretical" in text
    assert "was measured here" in text


def test_the_paragraph_reports_the_settings_that_were_used():
    text = ses.methods_paragraph(a_record(ph=6.5, window=11))
    assert "pH 6.5" in text
    assert "11 residues" in text


def test_the_paragraph_counts_the_sequences_without_naming_them():
    """Clone identifiers are internal, and this text goes into a manuscript."""
    text = ses.methods_paragraph(a_record(sequences={"internal-clone-7": SEQ}))
    assert "One amino-acid sequence was" in text
    assert "internal-clone-7" not in text


def test_the_paragraph_admits_when_the_record_cannot_say_what_ran():
    text = ses.methods_paragraph(a_record(method_used=""))
    assert "does not say which annotator ran" in text


def test_the_paragraph_carries_the_fingerprint_and_the_doi():
    text = ses.methods_paragraph(a_record(), doi="10.5281/zenodo.23742269")
    assert "abc1234" in text
    assert "10.5281/zenodo.23742269" in text


def test_the_paragraph_accepts_a_record_or_its_parsed_form():
    record = a_record()
    assert (ses.methods_paragraph(record)
            == ses.methods_paragraph(ses.session_inputs(record)))


# --------------------------------------------------------------------------
# The bounds, against the widgets that collect them
# --------------------------------------------------------------------------

def test_the_bounds_match_the_sliders_in_the_app():
    """`session.METHODS_SETTINGS` is what a file is validated against.

    If it disagrees with the sliders in `app.py`, one of two things happens:
    a file the app could show is refused, or a file the app cannot show is
    accepted and crashes the sidebar. Both are silent until somebody reopens
    a session, so the two are checked against each other here rather than
    kept in step by hand.
    """
    source = (ROOT / "app.py").read_text(encoding="utf-8")
    ph_slider = re.search(
        r'st\.slider\("pH for net charge",\s*([\d.]+),\s*([\d.]+)', source)
    window_slider = re.search(
        r'st\.slider\("Profile window \(residues\)",\s*(\d+),\s*(\d+)',
        source)
    assert ph_slider and window_slider, \
        "the sidebar sliders were renamed; update this test with them"
    assert ses.METHODS_SETTINGS["ph"] == (float(ph_slider.group(1)),
                                          float(ph_slider.group(2)))
    assert ses.METHODS_SETTINGS["window"] == (int(window_slider.group(1)),
                                              int(window_slider.group(2)))


def test_the_annotation_methods_are_the_ones_the_core_offers():
    """"auto", "anarci" and "regex" are `vhh_core.annotate`'s vocabulary."""
    assert set(ses.ANNOTATION_METHODS) == {"auto", "anarci", "regex"}
    assert vc.DEFAULTS["annotation_method"] in ses.ANNOTATION_METHODS


def test_the_sidebar_widgets_carry_the_keys_the_restore_writes():
    """The restore pushes `pending_set_ph` and the router moves it into place.

    That only works if the widgets are keyed, and the keys are spelt the same
    in two files. A typo here loses the setting silently: the session reloads,
    the sequences come back, and the pH is whatever the slider was on.
    """
    app = (ROOT / "app.py").read_text(encoding="utf-8")
    analysis = (ROOT / "sections" / "analysis.py").read_text(encoding="utf-8")
    for key in ("set_ph", "set_window", "set_method"):
        assert f'key="{key}"' in app, f"{key} is not a widget key in app.py"
        assert f'"pending_{key}"' in analysis, \
            f"nothing in the analysis screen fills pending_{key}"
    assert 'key.startswith("pending_")' in app, \
        "app.py no longer moves pending_ keys into place"
