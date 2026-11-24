"""
The Sanger screen keeps its chromatograms and its reading settings when the
reader leaves it and comes back.

It did not. Streamlit forgets a file uploader's files, and resets an unkeyed
selectbox or slider, as soon as the screen that draws them is not drawn, so
going to the Guides screen between step 1 (upload) and step 2 (choose the
vector) came back to an empty screen with the vector back on its default. The
same happened through the sidebar; the guides only made it the expected path.
Found by hand on the Guides screen, reproduced with these tests before the
fix, which failed on both the files and the vector.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1",
                    reason="this Streamlit is too old for AppTest")
from streamlit.testing.v1 import AppTest        # noqa: E402

import demo  # noqa: E402

pytestmark = pytest.mark.smoke

APP = str(Path(__file__).resolve().parent.parent / "app.py")
FILES = [(name, data, "application/octet-stream")
         for name, data in demo.sanger_example_files()]
SANGER = "🧪 Sanger sequencing"


def _fresh():
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    return at


def _click(at, label):
    at.button[[b.label for b in at.button].index(label)].click().run()


def _reads(at):
    return [m.value for m in at.metric if m.label == "Reads"]


def _upload(at, files=FILES):
    at.file_uploader[0].set_value(files).run()


def test_files_and_vector_survive_a_visit_elsewhere():
    at = _fresh()
    _click(at, SANGER)
    _upload(at)
    other = next(v for v in at.selectbox[0].options
                 if v != at.selectbox[0].value)
    at.selectbox[0].set_value(other).run()

    _click(at, "🏠 Home")
    _click(at, SANGER)
    assert not at.exception
    assert _reads(at) == [str(len(FILES))]
    assert at.selectbox[0].value == other
    assert any("still loaded" in str(i.value) for i in at.info)


def test_the_guide_path_keeps_the_files():
    """The path that found it: step 1, back to the guide, step 2."""
    at = _fresh()
    _click(at, "🧭 Guides")
    next(b for b in at.button
         if b.key == "guide_sanger_screening_1").click().run()
    _upload(at)
    _click(at, "🧭 Guides")
    next(b for b in at.button
         if b.key == "guide_sanger_screening_2").click().run()
    assert any(t.value == SANGER for t in at.title)
    assert _reads(at) == [str(len(FILES))]


def test_the_reading_settings_survive_a_visit_elsewhere():
    at = _fresh()
    _click(at, SANGER)
    at.slider[0].set_value(0.12).run()
    at.number_input[0].set_value(200).run()
    at.number_input[1].set_value(30).run()
    _click(at, "🔬 Sequence analysis")
    _click(at, SANGER)
    assert (at.slider[0].value, at.number_input[0].value,
            at.number_input[1].value) == (0.12, 200, 30)
    # A widget given a default and a Session State value at once warns on
    # screen; the defaults are set in Session State only, so nothing does.
    assert not [w for w in at.warning if "Session State" in str(w.value)]


def test_emptying_the_uploader_on_the_screen_forgets_the_files():
    at = _fresh()
    _click(at, SANGER)
    _upload(at)
    at.file_uploader[0].set_value(None).run()
    assert _reads(at) == []
    assert "ab1_kept" not in at.session_state


def test_forget_these_files_clears_the_kept_files():
    at = _fresh()
    _click(at, SANGER)
    _upload(at)
    _click(at, "🏠 Home")
    _click(at, SANGER)
    _click(at, "Forget these files")
    assert _reads(at) == []
    assert any("No chromatogram" in str(i.value) for i in at.info)


def test_uploading_again_replaces_rather_than_adds():
    at = _fresh()
    _click(at, SANGER)
    _upload(at)
    _click(at, "🏠 Home")
    _click(at, SANGER)
    _upload(at, FILES[:2])
    assert _reads(at) == ["2"]
    _click(at, "🏠 Home")
    _click(at, SANGER)
    assert _reads(at) == ["2"]
