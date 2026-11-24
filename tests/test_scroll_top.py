"""
A change of screen brings the page back to its top; a rerun of the same
screen does not.

Streamlit draws a new screen in the page as it was scrolled, so a navigation
button clicked at the bottom of a long screen landed at the bottom of the
next. The router now adds a small script on the run where the screen changes,
and on no other: a reader who moves a slider at the bottom of a screen has to
stay there.

What AppTest can see is whether the script is in the page. Whether it moves
the page is a browser question, and was measured in headless Chrome when this
was written: from 1947 px down the Guides screen, the Manual, the Sanger
screen through a guide's Go there, and the analysis through Send selected all
landed at 0, where without the script they landed at the bottom (258, 114 and
1074 px). The first version of the script was put in place without a keyed
container, and Chrome showed it never ran after a sidebar click; that is why
`test_each_change_of_screen_is_a_new_element` exists.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1",
                    reason="this Streamlit is too old for AppTest")
from streamlit.testing.v1 import AppTest        # noqa: E402

pytestmark = pytest.mark.smoke

APP = str(Path(__file__).resolve().parent.parent / "app.py")


def _scripts(at):
    return [h for h in at.get("html") if "scrollTop" in str(h.proto.body)]


def _click(at, label):
    at.button[[b.label for b in at.button].index(label)].click().run()


def test_a_change_of_screen_adds_the_script_and_a_rerun_does_not():
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    assert len(_scripts(at)) == 1          # first load: arriving is a change
    at.run()
    assert _scripts(at) == []              # same screen, ordinary rerun
    _click(at, "🧭 Guides")
    assert len(_scripts(at)) == 1
    at.selectbox[0].set_value(at.selectbox[0].value).run()
    assert _scripts(at) == []


def test_each_change_of_screen_is_a_new_element():
    """Swapping the content of one element does not run a script; a new
    element does. The keyed container is what makes it new."""
    import shell
    source = Path(shell.__file__).read_text()
    assert 'st.container(key=f"scroll_top_{screen}")' in source
