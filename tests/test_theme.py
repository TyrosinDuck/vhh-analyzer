"""
Theme tests, light and dark.

Two defects motivated this file, both invisible to every other test because
neither raises anything, and both seen only in a browser:

- A colour under a bare ``[theme]`` in `.streamlit/config.toml` makes
  Streamlit build a single custom theme on the light one. The Light / Dark /
  System choice then disappears from the app menu, and a visitor whose system
  is in dark mode sees the page paint dark and flip to white a moment later.
- The iframes (the home page diagram, the amino acid drawings) chose their
  text colour from ``prefers-color-scheme``, the operating system's setting,
  while the page around them follows the app's theme. When the two disagree
  the text vanishes into the background.

What a browser paints cannot be checked here; what can is that the
configuration and the code are in the shape that was verified in one.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import shell

ROOT = Path(__file__).resolve().parent.parent


def _config() -> dict:
    return tomllib.loads((ROOT / ".streamlit" / "config.toml").read_text())


def test_no_colour_is_set_under_a_bare_theme_table():
    theme = _config().get("theme", {})
    bare = {key: value for key, value in theme.items()
            if not isinstance(value, dict)}
    assert not bare, (
        f"[theme] sets {sorted(bare)}: that turns the theme into a single "
        "light-based custom theme, removes the Light/Dark choice from the "
        "menu and makes the page flash dark then white on a dark system. "
        "Put them under [theme.light] and [theme.dark] instead.")


def test_light_and_dark_themes_share_the_accent_colour():
    theme = _config()["theme"]
    assert theme["light"]["primaryColor"] == theme["dark"]["primaryColor"]


def test_no_section_follows_the_system_colour_scheme():
    """Every iframe goes through shell.follow_app_theme instead."""
    offenders = [p.name for p in (ROOT / "sections").glob("*.py")
                 if "prefers-color-scheme" in p.read_text()]
    assert not offenders, f"still keyed on the system theme: {offenders}"


def test_follow_app_theme_reads_the_page_around_the_iframe():
    snippet = shell.follow_app_theme()
    assert "window.parent.document" in snippet
    assert ".stApp" in snippet
    assert r"match(/\d+/g)" in snippet
    for colour in shell.THEME_TEXT.values():
        assert colour in snippet


def test_follow_app_theme_starts_light_outside_a_running_app():
    """No session, no theme: the light text, which is what Streamlit itself
    shows before it knows better."""
    snippet = shell.follow_app_theme()
    assert f"color:{shell.THEME_TEXT['light']}}}" in snippet
