"""
Tests for the guides: every step points at something that exists, the
checker catches a step that does not, and a step's button lands where it says.

The reference check is the reason the guides are data rather than prose. A
guide written as text would keep naming a tab after it had been renamed, and
nothing would say so; here the rename fails this file.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import guides as gd
import shell

LABELS = shell.screen_labels()
FAMILIES = shell.ANALYSIS_FAMILIES
VIEWS = shell.SCREEN_VIEWS


def test_every_guide_step_resolves():
    assert gd.problems(gd.GUIDES, LABELS, FAMILIES, screen_views=VIEWS) == []


def test_there_is_at_least_the_sanger_guide():
    assert "sanger_screening" in gd.GUIDES
    screens = {step.screen for step in gd.GUIDES["sanger_screening"].steps}
    assert {"sequencing", "analysis"} <= screens


def _one_step_guide(step):
    return {"g": gd.Guide("t", "q", (step,))}


def test_the_guides_cover_the_three_questions_they_claim_to():
    assert list(gd.GUIDES) == ["sanger_screening", "purified_batch",
                               "structure_measured"]
    batch = {s.view for s in gd.GUIDES["purified_batch"].steps}
    assert {"quantify", "balance", "batch"} <= batch
    model = {s.view for s in gd.GUIDES["structure_measured"].steps}
    assert {"disulfides", "liabilities"} <= model


def test_no_step_points_at_a_view_that_is_not_always_drawn():
    """The structure screen draws its interface view only with a binding
    partner, so a step opening it would land on a view that is not there."""
    for key, guide in gd.GUIDES.items():
        for step in guide.steps:
            assert step.view != "interface", key


@pytest.mark.parametrize("step, fragment", [
    (gd.Step("nowhere", "x"), "unknown screen"),
    (gd.Step("bench", "x", view="nope"), "unknown view"),
    (gd.Step("analysis", "x", view="batch"), "has no view selector"),
    (gd.Step("bench", "x", view="batch", family="compare"),
     "a view beside a family"),
    (gd.Step("analysis", "x", family="nope"), "unknown family"),
    (gd.Step("analysis", "x", family="compare", tab="Not a tab"),
     "is not a tab of"),
    (gd.Step("analysis", "x", tab="Clonal families"), "without a family"),
    (gd.Step("sequencing", "x", family="compare"), "outside the analysis"),
])
def test_a_broken_reference_is_reported(step, fragment):
    # One broken step can raise more than one complaint, and a step carrying
    # both a view and a family really is wrong twice over. What matters is
    # that the reason is named and that a sound guide raises nothing, which
    # `test_every_guide_step_resolves` is for.
    found = gd.problems(_one_step_guide(step), LABELS, FAMILIES,
                        screen_views=VIEWS)
    assert any(fragment in line for line in found), found


def test_a_guide_without_steps_is_reported():
    found = gd.problems({"g": gd.Guide("t", "q", ())}, LABELS, FAMILIES,
                        screen_views=VIEWS)
    assert found == ["g: no steps"]


def test_a_renamed_tab_breaks_the_check():
    """What the checker exists for: rename a tab and the guide fails."""
    renamed = {label: (key, [("Families" if t == "Clonal families" else t)
                             for t in tabs], blurb)
               for label, (key, tabs, blurb) in FAMILIES.items()}
    found = gd.problems(gd.GUIDES, LABELS, renamed, screen_views=VIEWS)
    assert any("'Clonal families' is not a tab" in line for line in found)


def test_a_renamed_view_key_breaks_the_check():
    """The same guarantee for the screens that have a view selector."""
    bench_key, bench_views = VIEWS["bench"]
    without = {k: v for k, v in bench_views.items() if k != "batch"}
    found = gd.problems(gd.GUIDES, LABELS, FAMILIES,
                        screen_views={**VIEWS, "bench": (bench_key, without)})
    assert any("unknown view 'batch'" in line for line in found)


def test_the_screens_with_a_view_selector_read_their_labels_from_shell():
    """One copy of each label: the screen and the guides must agree."""
    import sections.bench as bench_screen
    import sections.structure as structure_screen
    assert list(bench_screen.VIEWS) == list(VIEWS["bench"][1].values())
    assert structure_screen.V_DISULFIDES == VIEWS["structure"][1]["disulfides"]
    assert structure_screen.V_INTERFACE == VIEWS["structure"][1]["interface"]


def test_a_view_shows_in_the_location_line():
    step = gd.Step("bench", "x", view="batch")
    assert gd.location(step, LABELS, FAMILIES, VIEWS) == (
        f"{LABELS['bench']} › {VIEWS['bench'][1]['batch']}")


def test_location_is_derived_from_the_labels():
    step = gd.Step("analysis", "x", family="compare", tab="Clonal families")
    compare_label = next(label for label, (key, *_r) in FAMILIES.items()
                         if key == "compare")
    assert gd.location(step, LABELS, FAMILIES) == (
        f"{LABELS['analysis']} › {compare_label} › Clonal families")
    assert gd.location(gd.Step("sequencing", "x"), LABELS, FAMILIES) == \
        LABELS["sequencing"]


def test_screen_labels_cover_sections_and_utilities_only():
    assert set(LABELS) == set(shell.SECTIONS) | set(shell.UTILITIES)
    assert LABELS["guides"] == "🧭 Guides"


def test_the_analysis_screen_uses_the_shared_family_table():
    """The table moved to `shell` so the Guides screen could read it without
    reading another screen's globals. The analysis must use that very object,
    not a copy that could drift from it."""
    import sections.analysis as analysis
    assert analysis._FAMILIES is shell.ANALYSIS_FAMILIES


def test_guide_text_has_no_em_dash():
    """House style for everything the app displays."""
    for guide in gd.GUIDES.values():
        texts = [guide.title, guide.question, guide.note]
        texts += [t for s in guide.steps for t in (s.action, s.why)]
        assert not any("\u2014" in t for t in texts)


# --------------------------------------------------------------------------- #
# The screen, driven through Streamlit's AppTest
# --------------------------------------------------------------------------- #

APP = str(Path(__file__).resolve().parent.parent / "app.py")


@pytest.mark.smoke
def test_go_there_opens_the_analysis_on_the_steps_family():
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    # Something to analyse, so the family selector is drawn at all.
    at.session_state["pending_seq_text"] = (
        ">Nb1\nQVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGST"
        "YYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDYWGQGTQVTVSS\n")
    at.run()

    labels = [b.label for b in at.button]
    at.button[labels.index("🧭 Guides")].click().run()
    assert not at.exception
    assert "🧭 Guides" in [t.value for t in at.title]

    steps = gd.GUIDES["sanger_screening"].steps
    number = next(i for i, s in enumerate(steps, 1) if s.family == "compare")
    key = f"guide_sanger_screening_{number}"
    next(b for b in at.button if b.key == key).click().run()
    assert not at.exception and not at.error

    assert "🔬 Sequence analysis" in [t.value for t in at.title]
    selector = [r for r in at.radio if r.label == "Section"]
    assert selector, "the family selector was not drawn"
    assert selector[0].value == gd.family_label("compare", FAMILIES)
    # The tab request was read when the family drew its tabs, and is gone, so
    # the next visit to the family opens on its first tab as usual. Which tab
    # the browser shows is not visible to AppTest; it was checked in headless
    # Chrome for every step of the Sanger guide when this was written.
    assert "analysis_open_tab" not in at.session_state


def test_a_step_with_a_view_asks_for_it_through_the_pending_key():
    import sections.guides as screen

    step = next(s for s in gd.GUIDES["purified_batch"].steps
                if s.view == "batch")
    state = {}
    original = screen.st.session_state
    try:
        screen.st.session_state = state
        screen._open(step)
    finally:
        screen.st.session_state = original
    state_key = VIEWS["bench"][0]
    assert state["section"] == "bench"
    assert state[f"pending_{state_key}"] == VIEWS["bench"][1]["batch"]


def test_a_step_with_a_tab_asks_for_it():
    import sections.guides as screen

    step = next(s for s in gd.GUIDES["sanger_screening"].steps
                if s.tab == "Clonal families")
    state = {}
    original = screen.st.session_state
    try:
        screen.st.session_state = state
        screen._open(step)
    finally:
        screen.st.session_state = original
    assert state["section"] == "analysis"
    assert state["analysis_open_tab"] == "Clonal families"
    assert state["pending_analysis_family"] == gd.family_label("compare",
                                                                FAMILIES)


def test_tabs_default_only_when_asked_and_supported():
    import inspect
    import streamlit as st
    assert shell.tabs_default(None) == {}
    supported = "default" in inspect.signature(st.tabs).parameters
    assert shell.tabs_default("Session") == (
        {"default": "Session"} if supported else {})


@pytest.mark.smoke
def test_the_chosen_question_survives_a_step():
    """Found in the browser with the second guide: leaving the screen reset
    the selector to the first question, so coming back after step 2 offered
    the wrong guide's steps. Invisible while there was only one guide."""
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    labels = [b.label for b in at.button]
    at.button[labels.index("🧭 Guides")].click().run()
    at.selectbox[0].set_value("purified_batch").run()
    next(b for b in at.button
         if b.key == "guide_purified_batch_2").click().run()
    assert "⚗️ Batch & QC" in [t.value for t in at.title]
    labels = [b.label for b in at.button]
    at.button[labels.index("🧭 Guides")].click().run()
    assert at.selectbox[0].value == "purified_batch"
    assert any(b.key == "guide_purified_batch_2" for b in at.button)


@pytest.mark.smoke
def test_a_view_step_opens_its_view():
    """Which view the radio shows is visible to AppTest, unlike a tab."""
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    labels = [b.label for b in at.button]
    at.button[labels.index("🧭 Guides")].click().run()
    at.selectbox[0].set_value("purified_batch").run()
    step_no, step = next((i, s) for i, s
                         in enumerate(gd.GUIDES["purified_batch"].steps, 1)
                         if s.view == "batch")
    next(b for b in at.button
         if b.key == f"guide_purified_batch_{step_no}").click().run()
    assert not at.exception
    view = [r for r in at.radio if r.key == VIEWS["bench"][0]]
    assert view and view[0].value == VIEWS["bench"][1]["batch"]


@pytest.mark.smoke
def test_go_there_on_a_plain_screen_lands_on_it():
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    at = AppTest.from_file(APP, default_timeout=600)
    at.run()
    labels = [b.label for b in at.button]
    at.button[labels.index("🧭 Follow a guide")].click().run()
    assert "🧭 Guides" in [t.value for t in at.title]
    next(b for b in at.button
         if b.key == "guide_sanger_screening_1").click().run()
    assert not at.exception
    assert "🧪 Sanger sequencing" in [t.value for t in at.title]

