"""
Smoke tests, does every screen of the app actually render?

The other two test files exercise the computation modules. Nothing in them would
notice if `app.py` itself stopped working, and `app.py` is where the damage from
a refactor lands: it is a single script that Streamlit re-executes top to bottom,
so a name that moved to another module without being imported back is a
``NameError`` on one screen and silence everywhere else.

That is not hypothetical. Splitting the cached computation layer out of `app.py`
moved ``analyze_with_progress`` into `compute.py` while the two-line switch it
depends on stayed behind. Every section still loaded; the analysis blew up the
moment you pressed **Analyze**. These tests are the net that catches that class
of break without anyone having to remember to click through six sections and
thirteen panels.

They drive the real app through Streamlit's own ``AppTest`` harness, so they are
slower than the rest of the suite, roughly a minute, against seven seconds for
everything else. Skip them while iterating with:

    pytest -m "not smoke"

and run the whole suite before pushing.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1",
                    reason="this Streamlit is too old for AppTest")
from streamlit.testing.v1 import AppTest        # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
APP = str(ROOT / "app.py")

pytestmark = pytest.mark.smoke

# Seven VHHs differing only in the last CDR3 residue: enough to exercise the
# multi-sequence paths (alignment, logo, identity matrices, clonal families)
# without making the sweep slow.
_FRAME = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
          "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA{cdr3}WGQGTQVTVSS")
FASTA = "".join(f">Nb{i}\n{_FRAME.format(cdr3='DRGYSSSWYD' + last)}\n"
                for i, last in enumerate("YFWLKRD"))

# Derived from shell.SECTIONS rather than written out. The list used to be
# typed here, a second list of the sections, and a section added to SECTIONS
# but not here would have been left out of every smoke test without a sound.
import shell  # noqa: E402

from conftest import HAS_HMMER  # noqa: E402

SECTION_BUTTONS = [f"{icon} {label}"
                   for icon, label, *_ in shell.SECTIONS.values()]
# The manual and the reference tables come from `shell.UTILITIES` for the same
# reason: they were written out here, and a third screen added beside them
# would have been navigated to by nothing.
UTILITY_BUTTONS = [f"{icon} {label}"
                   for icon, label in shell.UTILITIES.values()]
NAV = SECTION_BUTTONS + UTILITY_BUTTONS + ["🏠 Home"]


def _labels(app):
    return [button.label for button in app.button]


def _click(app, label):
    app.button[_labels(app).index(label)].click().run()
    return app


def _goto(app, label):
    """Navigate to a section, whether or not we are already on it.

    The nav button for the current section is disabled, and `AppTest` refuses
    to click a disabled widget, exactly as a browser user could not. The
    module-scoped `app` fixture means a test can inherit whatever section the
    previous one left behind, so a bare `_click` here passes or fails by test
    order.
    """
    if any(title.value == label for title in app.title):
        return app
    return _click(app, label)


def _no_failure(app, where):
    assert not app.exception, \
        f"{where}: {[str(e.value) for e in app.exception]}"
    assert not app.error, f"{where}: {[e.value for e in app.error]}"


@pytest.fixture(scope="module")
def app():
    at = AppTest.from_file(APP, default_timeout=2400)
    at.run()
    _no_failure(at, "initial load")
    return at


# Every screen draws its own title, except the manual, which renders MANUAL.md
# as markdown and has no st.title of its own.
# Every screen draws its own sidebar label as its title, which is what makes
# this derivable from `shell` rather than typed out: a screen added to
# SECTIONS or UTILITIES is covered here the moment it exists, where the
# written-out version raised a KeyError in the middle of a sweep and said
# nothing about which screen was missing a title.
#
# Two exceptions, both deliberate: the home page is titled with the app rather
# than with the button, and the manual draws the document instead of a title
# of its own.
EXPECTED_TITLE = {**{button: button for button in SECTION_BUTTONS},
                  **{button: button for button in UTILITY_BUTTONS},
                  "🏠 Home": "🧰 VHH tool box",
                  "📖 Manual": None}


def _content_count(app):
    return (len(app.markdown) + len(app.title) + len(app.header)
            + len(app.subheader) + len(app.caption))


def test_every_section_renders(app):
    """The version banner calls ``st.stop()``, so a mismatch shows up here too."""
    for label in NAV:
        _click(app, label)
        _no_failure(app, label)


def test_every_section_actually_draws_something(app):
    """Not raising is not the same as working.

    Splitting the screens out of `app.py` put the Gene diagnosis screen's body
    inside its own helper function by accident: `render()` was left holding
    nothing but a docstring. It imported cleanly, raised nothing, and drew an
    empty page, and ``test_every_section_renders`` passed it, because all that
    test asked was whether anything had blown up.

    So this one asks for evidence instead: the screen's own title, and enough
    elements that an empty body cannot pass. The broken version produced no
    title and two elements, both of them the sidebar's.
    """
    for label in NAV:
        _click(app, label)
        _no_failure(app, label)

        expected = EXPECTED_TITLE[label]
        titles = [t.value for t in app.title]
        if expected is not None:
            assert expected in titles, \
                f"{label}: expected the title {expected!r}, found {titles}"

        count = _content_count(app)
        assert count >= 6, \
            f"{label}: only {count} text elements, the screen looks empty"


def test_the_footer_is_drawn_exactly_once(app):
    """Every screen drew it twice after the split, for a week.

    The router calls ``footer()`` once ``render()`` returns. Each screen module
    also ended with its own ``footer()``, left behind when the section bodies
    were lifted out of `app.py`, the extraction stripped the trailing
    ``st.stop()`` and not the ``footer()`` in front of it. Nothing raised and
    nothing looked empty, so neither of the other smoke tests noticed; the only
    symptom was the licence notice and the version line appearing twice at the
    bottom of every page, which a user reported.

    The ``footer()`` calls still inside the screen modules are deliberate: each
    one precedes an early ``st.stop()``, so ``render()`` never returns and the
    router's call is never reached.
    """
    for label in NAV:
        _click(app, label)
        _no_failure(app, label)
        drawn = [c.value for c in app.caption if "Built on Biopython" in str(c.value)]
        assert len(drawn) == 1, \
            f"{label}: the footer was drawn {len(drawn)} times, expected once"


def test_analysis_renders_every_family(app):
    """The break that prompted this file: the sections loaded, Analyze did not."""
    _click(app, "🔬 Sequence analysis")
    app.text_area[0].set_value(FASTA)
    _click(app, "Analyze")
    _no_failure(app, "Analyze")

    selector = [r for r in app.radio if r.label == "Section"]
    assert selector, "the analysis family selector did not render"
    families = selector[0].options
    assert len(families) == 6, f"expected six families, got {families}"

    for family in families:
        [r for r in app.radio if r.label == "Section"][0].set_value(family).run()
        _no_failure(app, f"family {family}")


def test_profiles_draw_one_sequence_unless_all_are_asked_for(app):
    """With 50 sequences, the profile figure alone took 516 MB to draw.

    The Plots tab now draws only the selected sequence, like Chou–Fasman, and
    ``Show all`` draws one image per sequence. Both paths must render.
    """
    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    identity = [o for o in [r for r in app.radio if r.label == "Section"][0].options
                if "Identity" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(identity).run()
    _no_failure(app, "Identity family")

    picks = [s for s in app.selectbox if s.key == "prof_pick"]
    assert picks, "the profile sequence selector did not render"
    picks[0].set_value(picks[0].options[-1]).run()
    _no_failure(app, "profile of the last sequence")

    [c for c in app.checkbox if c.key == "prof_all"][0].check().run()
    _no_failure(app, "all profiles")
    assert not [s for s in app.selectbox if s.key == "prof_pick"], \
        "Show all should replace the selector"
    [c for c in app.checkbox if c.key == "prof_all"][0].uncheck().run()


def test_properties_show_both_extinction_coefficients_and_mass_in_da(app):
    """Batch & QC had both extinction coefficients, Sequence analysis one.

    The Properties table now carries the cystine and the reduced values, the
    A280 factor for each, and the mass in Da to two decimals rather than kDa.
    """
    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    identity = [o for o in [r for r in app.radio if r.label == "Section"][0].options
                if "Identity" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(identity).run()
    _no_failure(app, "Identity family")
    wanted = {"mw_da", "ext_coeff_cystines", "ext_coeff_reduced",
              "A280_0.1pct_cystines", "A280_0.1pct"}
    tables = [d.value for d in app.dataframe]
    tables = [getattr(t, "data", t) for t in tables]     # a Styler or a frame
    match = [t for t in tables if wanted <= set(t.columns)]
    assert match, "the Properties table lost a column"
    assert "mw_kda" not in match[0].columns
    assert (match[0]["ext_coeff_cystines"] >= match[0]["ext_coeff_reduced"]).all()


def test_identity_heatmap_selector_and_cdr3_note(app):
    """One identity heatmap at a time, and both tabs say how CDR3 identity
    differs between them."""
    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    comparison = [o for o in [r for r in app.radio if r.label == "Section"][0].options
                  if "Comparison" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(comparison).run()
    _no_failure(app, "Comparison family")

    view = [r for r in app.radio if r.key == "id_view"]
    assert view and list(view[0].options) == ["Full-length", "CDR3"]
    labels = [b.label for b in app.get("download_button")]
    assert "Full-length identity (CSV)" in labels and "CDR3 identity (CSV)" in labels

    note = "is not computed the same way here and in"
    notes = [str(c.value) for c in app.caption if note in str(c.value)]
    assert not any("**Clonal families**" in n for n in notes), \
        "the Alignment tab carries the note in its CDR3 view only"
    view[0].set_value("CDR3").run()
    _no_failure(app, "CDR3 identity view")
    notes = [str(c.value) for c in app.caption if note in str(c.value)]
    assert any("**Clonal families**" in n for n in notes)
    assert any("**Alignment**" in n for n in notes), \
        "Clonal families should carry the same note"
    [r for r in app.radio if r.key == "id_view"][0].set_value("Full-length").run()


def test_the_reference_comparison_says_what_it_needs_without_hmmer(app):
    """Without HMMER this view used to show the operating system's errno.

    `compare_to_reference` numbers the reference outside the per-sequence
    `try`, so a missing `hmmscan` left the screen reading `Could not compare:
    [Errno 2] No such file or directory: 'hmmscan'`. It is the only view in
    the app with no regex fallback, and the person most likely to meet it is
    the one who followed the privacy notice and cloned the repository, where
    `pip install -r requirements.txt` brings in ANARCI and nothing brings in
    HMMER.

    Both branches assert: with HMMER the comparison must still compute, and
    without it the screen must name what to install. A test that skipped here
    would let that message disappear unnoticed on the only machines that see
    it.
    """
    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    selector = [r for r in app.radio if r.label == "Section"][0]
    comparison = [o for o in selector.options if "Comparison" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(comparison).run()
    _no_failure(app, "Comparison family")

    said = " ".join(str(info.value) for info in app.info)
    if HAS_HMMER:
        assert "ANARCI" not in said, \
            f"nothing is missing, so nothing should be explaining itself: {said!r}"
        tables = [getattr(d.value, "data", d.value) for d in app.dataframe]
        assert any({"name", "n_differences"} <= set(t.columns) for t in tables), \
            "the reference comparison summary did not render"
    else:
        assert "ANARCI" in said and "HMMER" in said, \
            f"the view did not say what it needs, it said: {said!r}"
        assert "install hmmer" in said, \
            "the message should carry the command that fixes it"


@pytest.mark.parametrize("with_cf", [False, True])
def test_the_pdf_report_builds(app, with_cf):
    """The PDF report raised ``NameError`` on every press of its button.

    `VERSION`, `AUTHOR` and `APP_NAME` were used in the Download tab but never
    imported from `shell`, and the Chou–Fasman window was read through a local
    variable of the Plots tab. Every screen still rendered, so nothing but
    pressing the button showed it.
    """
    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    outputs = [o for o in [r for r in app.radio if r.label == "Section"][0].options
               if "Outputs" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(outputs).run()
    app.session_state["pdf_bytes"] = None
    box = [c for c in app.checkbox
           if c.label.startswith("Include the per-sequence secondary")][0]
    (box.check() if with_cf else box.uncheck()).run()
    _click(app, "Build PDF report")
    assert not app.exception, [str(e.value) for e in app.exception]
    pdf = app.session_state["pdf_bytes"]
    assert pdf and pdf[:5] == b"%PDF-", "no PDF was produced"


def test_structure_renders_every_view(app, tmp_path):
    """A synthetic two-chain PDB: enough to reach the interface and 3D views.

    Built here rather than shipped so the suite needs no data files, and
    deliberately tiny, nothing in this test reads the geometry.
    """
    lines, serial = [], 0
    for chain, offset in (("A", 0.0), ("B", 6.0)):
        for index, amino in enumerate("ACDEFGHIKLMNPQRSTVWY" * 2):
            serial += 1
            three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                     "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                     "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
                     "S": "SER", "T": "THR", "V": "VAL", "W": "TRP",
                     "Y": "TYR"}[amino]
            lines.append(
                f"ATOM  {serial:5d}  CA  {three} {chain}{index + 1:4d}    "
                f"{index * 3.8:8.3f}{offset:8.3f}{0.0:8.3f}  1.00 50.00           C")
    pdb = tmp_path / "pair.pdb"
    pdb.write_text("\n".join(lines) + "\nTER\nEND\n")

    app.session_state["struct_raw"] = pdb.read_bytes()
    app.session_state["struct_name"] = "pair.pdb"
    _click(app, "🧊 Structure")
    _no_failure(app, "structure section")

    selector = [r for r in app.radio if r.label == "View"]
    assert selector, "the structure view selector did not render"
    for view in selector[0].options:
        [r for r in app.radio if r.label == "View"][0].set_value(view).run()
        _no_failure(app, f"structure view {view}")


def test_structure_section_reads_an_alphafold3_zip(app):
    """The .zip path: model picker, ipTM panel and every view, with no failure.

    The ipTM here is high on purpose. A low one draws an ``st.error`` by
    design, and ``_no_failure`` would read that as the section breaking.
    """
    from test_af3 import server_zip

    app.session_state["struct_raw"] = server_zip(iptms=(0.9, 0.85))
    app.session_state["struct_name"] = "fold_test_job.zip"
    # The previous test leaves the app on this screen, whose own navigation
    # button is then disabled; going home first makes the order not matter.
    _click(app, "🏠 Home")
    _click(app, "🧊 Structure")
    _no_failure(app, "structure section, AlphaFold3 zip")
    pickers = [s for s in app.selectbox if s.label.startswith("Model from")]
    assert pickers, "the AlphaFold3 model picker did not render"
    selector = [r for r in app.radio if r.label == "View"][0]
    selector.set_value("🧬 Chains & confidence").run()
    _no_failure(app, "AlphaFold3 confidence panel")
    labels = [m.label for m in app.metric]
    assert "ipTM, whole complex" in labels
    assert any("ipTM 0.90" in s.value for s in app.success)
    for view in [r for r in app.radio if r.label == "View"][0].options:
        [r for r in app.radio if r.label == "View"][0].set_value(view).run()
        _no_failure(app, f"AlphaFold3 zip, view {view}")
    pickers = [s for s in app.selectbox if s.label.startswith("Model from")]
    pickers[0].set_value("model 1").run()
    _no_failure(app, "AlphaFold3 zip, second model")


def test_sequence_input_survives_leaving_the_section(app):
    """Streamlit discards a widget's state once the widget stops rendering.

    Leaving a section is exactly that, so the text area is re-asserted on every
    run. If that stops happening, your sequences vanish when you glance at
    another section and come back, which is worse than an error, because
    nothing tells you it happened.
    """
    _click(app, "🔬 Sequence analysis")
    app.text_area[0].set_value(FASTA)
    _click(app, "Analyze")
    _click(app, "🧬 In-silico cloning")
    _click(app, "🔬 Sequence analysis")
    _no_failure(app, "returning to the analysis")
    assert app.text_area[0].value.strip() == FASTA.strip(), \
        "the sequence input was lost on leaving and re-entering the section"


def test_bench_renders_every_view(app):
    """Every view of Batch & QC, with no sequence loaded.

    The expected views are read from the module rather than counted here. A
    hard-coded five was a second list of the same thing: it failed when a view
    was added, which is the right alarm for the wrong reason, and it would
    have stayed quiet if a view had been renamed.
    """
    from sections import bench as bench_screen

    _click(app, "⚗️ Batch & QC")
    _no_failure(app, "bench section")

    selector = [r for r in app.radio if r.label == "View"]
    assert selector, "the bench view selector did not render"
    views = selector[0].options
    assert list(views) == list(bench_screen.VIEWS), (
        f"the selector offers {views}, the module lists "
        f"{bench_screen.VIEWS}")

    for view in views:
        [r for r in app.radio if r.label == "View"][0].set_value(view).run()
        _no_failure(app, f"bench view {view}")
        assert _content_count(app) >= 6, f"bench view {view} looks empty"


def test_the_phage_section_computes_both_titres_and_both_volumes(app):
    """The phage titre in its own section, with no sequence anywhere.

    It used to be a view of Batch & QC reached through an early return;
    moved out, it is checked for what it is for rather than for drawing:
    the two titres and the two volumes. The numbers are the screen's
    defaults worked by hand: 0.459 x 100 = 45.9; 45.9 x 6e16 / 7249 is
    3.80e14 /mL, and 1e11 particles of that is 0.263 uL, below the 2 uL a
    pipette can take, so the volume shown is the one after the dilution
    `bench.phage_aliquot` proposes. Only the presence of both volumes is
    checked here; their arithmetic is pinned in tests/test_bench.py.
    """
    _goto(app, "🦠 Phage & panning")
    _no_failure(app, "phage section")
    assert not [r for r in app.radio if r.label == "Sequence"], \
        "the phage section asks for a sequence it does not use"
    metrics = {m.label: m.value for m in app.metric}
    assert set(metrics) >= {"Empirical", "From the genome", "Corrected A260"}
    assert metrics["Corrected A260"] == "45.9"
    assert metrics["From the genome"] == "3.80e+14 /mL"
    shown = " ".join(s.value for s in app.success)
    assert "Solid phase" in shown and "Liquid phase" in shown


def test_the_buffers_section_prepares_pbs_at_its_usual_stock(app):
    """The screen opens on the first recipe at its usual factor, and changing
    the form of a salt changes the mass shown. PBS is made at 20X: 160 g of
    NaCl per litre, and 35.60 g of the dihydrate phosphate or, as the
    anhydrous salt, 35.6 x 141.958 / 177.988 = 28.39 g."""
    _goto(app, "🧂 Buffers & media")
    _no_failure(app, "buffers section")
    assert app.number_input(key="buf_factor_PBS").value == 20.0
    table = app.dataframe[0].value
    amounts = dict(zip(table["Component"], table["Amount"]))
    assert amounts["Sodium chloride"] == "160.00 g"
    assert amounts["Disodium hydrogen phosphate"] == "35.60 g"
    app.selectbox(key="buf_form_PBS_Disodium hydrogen phosphate").set_value(
        "anhydrous").run()
    _no_failure(app, "buffers section, anhydrous phosphate")
    table = app.dataframe[0].value
    amounts = dict(zip(table["Component"], table["Amount"]))
    assert amounts["Disodium hydrogen phosphate"] == "28.39 g"
    steps = " ".join(m.value for m in app.markdown)
    assert "Do not adjust the pH" in steps


def test_composing_a_buffer_weighs_the_hydrate_chosen(app):
    """The other view of Buffers & media: compounds picked by hand rather
    than read from a recipe. 29.21 mM of phosphate in 1 L is 4.146 g of the
    anhydrous salt and 5.198 g of the dihydrate, the same moles either way.
    29.21 mM of the anhydrous salt is 0.02921 x 141.958 = 4.147 g, and of
    the dihydrate 0.02921 x 177.988 = 5.199 g. `AppTest` has no accessor for
    `st.multiselect`, so the selection is injected through its session-state
    entry, as the purification balance test does for the data editor.
    """
    _goto(app, "🧂 Buffers & media")
    app.radio(key="buf_view").set_value("🧪 Compose one").run()
    _no_failure(app, "buffers, compose")
    assert any("Pick at least one compound" in info.value
               for info in app.info), "the empty state said nothing"

    app.session_state["cus_compounds"] = ["Disodium hydrogen phosphate"]
    app.run()
    app.number_input(key="cus_amount_Disodium hydrogen phosphate").set_value(
        29.21).run()
    _no_failure(app, "buffers, composed buffer")
    amounts = dict(zip(app.dataframe[0].value["Component"],
                       app.dataframe[0].value["Amount"]))
    assert amounts["Disodium hydrogen phosphate"] == "4.147 g"

    app.selectbox(key="cus_form_Disodium hydrogen phosphate").set_value(
        "dihydrate").run()
    _no_failure(app, "buffers, composed buffer, dihydrate")
    amounts = dict(zip(app.dataframe[0].value["Component"],
                       app.dataframe[0].value["Amount"]))
    assert amounts["Disodium hydrogen phosphate"] == "5.199 g"


def test_the_enrichment_view_counts_a_round(app):
    """The second view of Phage & panning, in both dilution conventions.

    120 colonies at 10^-3 against 3 at 10^-1 is an enrichment of 4000,
    whether the dilutions are written as exponents (3 and 1, the default) or
    as factors (1000 and 10). The counts go in through the data editor's
    session-state entry, as the purification balance does.
    """
    from sections import phage as phage_screen

    _goto(app, "🦠 Phage & panning")
    app.radio(key="phage_view").set_value(phage_screen.V_ENRICHMENT).run()
    _no_failure(app, "phage, enrichment, empty")
    assert any("four numbers" in info.value for info in app.info), \
        "the empty state did not say what it needs"
    assert app.radio(key="en_dilution_mode").value == phage_screen.D_POWER, \
        "the exponent is what is written on the plate; it is the default"

    app.session_state["en_rounds"] = {
        "edited_rows": {0: {phage_screen.C_PLUS: 120,
                            phage_screen.C_PLUS_DIL: 3,
                            phage_screen.C_MINUS: 3,
                            phage_screen.C_MINUS_DIL: 1}},
        "added_rows": [], "deleted_rows": []}
    app.run()
    _no_failure(app, "phage, enrichment, exponents")
    # The last dataframe on the page: the first is the editor itself.
    table = app.dataframe[-1].value
    assert list(table["Enrichment (+/-)"]) == ["4 000"]
    assert list(table["+ output (cfu/mL)"]) == ["1.20e+06"]

    app.radio(key="en_dilution_mode").set_value(phage_screen.D_FACTOR).run()
    app.session_state["en_rounds"] = {
        "edited_rows": {0: {phage_screen.C_PLUS: 120,
                            phage_screen.C_PLUS_DIL: 1000,
                            phage_screen.C_MINUS: 3,
                            phage_screen.C_MINUS_DIL: 10}},
        "added_rows": [], "deleted_rows": []}
    app.run()
    _no_failure(app, "phage, enrichment, factors")
    assert list(app.dataframe[-1].value["Enrichment (+/-)"]) == ["4 000"]


def test_an_exponent_typed_into_the_factor_column_is_flagged(app):
    """The silent wrong answer: 30 colonies at "3" against 30 at "2" reads
    as an enrichment of 1.5 where 10 was meant. The number is still shown,
    because a 1:2 dilution is a real thing, with a warning beside it."""
    from sections import phage as phage_screen

    _goto(app, "🦠 Phage & panning")
    app.radio(key="phage_view").set_value(phage_screen.V_ENRICHMENT).run()
    app.radio(key="en_dilution_mode").set_value(phage_screen.D_FACTOR).run()
    app.session_state["en_rounds"] = {
        "edited_rows": {0: {phage_screen.C_PLUS: 30,
                            phage_screen.C_PLUS_DIL: 3,
                            phage_screen.C_MINUS: 30,
                            phage_screen.C_MINUS_DIL: 2}},
        "added_rows": [], "deleted_rows": []}
    app.run()
    _no_failure(app, "phage, enrichment, exponent in the factor column")
    assert list(app.dataframe[-1].value["Enrichment (+/-)"]) == ["1.50"]
    assert any(phage_screen.D_POWER in warning.value
               for warning in app.warning), \
        "nothing pointed at the other convention"


def test_batch_and_qc_no_longer_offers_the_phage_titre(app):
    from sections import bench as bench_screen

    assert not any("Phage" in view for view in bench_screen.VIEWS)
    _goto(app, "⚗️ Batch & QC")
    selector = [r for r in app.radio if r.label == "View"][0]
    assert not any("Phage" in view for view in selector.options)


def test_finished_results_ask_how_the_app_is_used(app):
    """The usage request is repeated under the batch sheet and the analysis
    downloads, with a link that says where it was clicked.

    On the home page it sits in a collapsed panel at the bottom of a page
    people leave as soon as they have picked a section.
    """
    import shell
    from sections import bench as bench_screen

    _goto(app, "⚗️ Batch & QC")
    [r for r in app.radio if r.label == "View"][0].set_value(
        bench_screen.V_BATCH).run()
    _no_failure(app, "batch sheet")
    url = shell.usage_url("Batch & QC, batch sheet")
    assert any(url in str(c.value) for c in app.caption), \
        "no usage request under the batch sheet"

    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    outputs = [o for o in [r for r in app.radio if r.label == "Section"][0].options
               if "Outputs" in o][0]
    [r for r in app.radio if r.label == "Section"][0].set_value(outputs).run()
    url = shell.usage_url("Sequence analysis, Download")
    assert any(url in str(c.value) for c in app.caption), \
        "no usage request under the analysis downloads"


def test_the_purification_balance_computes_in_the_screen(app):
    """The view with the table filled in, not only its empty state.

    `AppTest` has no accessor for `st.data_editor`, so the edit is injected
    through the widget's session-state entry, which is how Streamlit carries
    an edited cell back into the next run. Without this the smoke test only
    ever saw the "give a step a volume" message, and the table, the metrics
    and the recovery column, which is everything this view is for, were
    covered by no test that runs the screen.

    The numbers are checked by hand: 1.2 AU / 28000 /M/cm * 14000 Da is
    0.6 mg/mL, and 6 mL of that is 3.6 mg. The extinction coefficient and the
    molar mass are set here rather than left to the view's defaults, because
    the fixture is shared and an earlier test may have loaded a sequence,
    whose own values would then be the ones used.
    """
    from sections import bench as bench_screen

    _goto(app, "⚗️ Batch & QC")
    [r for r in app.radio if r.label == "View"][0].set_value(
        bench_screen.V_BALANCE).run()
    _no_failure(app, "purification balance, empty")
    assert any("at least one step" in info.value for info in app.info), \
        "the empty balance did not say what it needs"

    app.session_state["pb_eps"] = 28000.0
    app.session_state["pb_mw"] = 14000.0
    app.session_state["pb_steps"] = {
        "edited_rows": {1: {"Volume (mL)": 6.0, "A280": 1.2}},
        "added_rows": [], "deleted_rows": []}
    app.run()
    _no_failure(app, "purification balance, one step")

    totals = [m for m in app.metric if m.label == "Final pool"]
    assert totals, "the final pool metric did not render"
    assert "3.6" in totals[0].value, \
        f"expected 3.600 mg in the final pool, got {totals[0].value}"

    # The fixture is shared, so the edited cell is put back as it was found
    # rather than left for whichever test runs next.
    app.session_state["pb_steps"] = {"edited_rows": {}, "added_rows": [],
                                     "deleted_rows": []}
    app.run()


def test_reference_renders_every_view(app):
    """The structures are files on disk, so a deployment that did not copy
    `data/aminoacids/` draws a table and a warning rather than failing. This
    checks the ordinary case: both views draw, with the files present."""
    _click(app, "📐 Reference")
    _no_failure(app, "reference section")

    selector = [r for r in app.radio if r.label == "View"]
    assert selector, "the reference view selector did not render"
    for view in selector[0].options:
        [r for r in app.radio if r.label == "View"][0].set_value(view).run()
        _no_failure(app, f"reference view {view}")
        assert _content_count(app) >= 6, f"reference view {view} looks empty"


class _BrokenImport:
    """A meta-path finder that fails one import the way a deploy does."""

    def __init__(self, name, error):
        self.name, self.error = name, error

    def find_spec(self, fullname, path=None, target=None):
        if fullname == self.name:
            raise self.error
        return None


@pytest.mark.parametrize("error", [KeyError("generator"),
                                   ImportError("cannot import name 'x'")],
                         ids=["KeyError", "ImportError"])
def test_a_run_mid_deploy_says_the_app_is_updating(error, monkeypatch):
    """``KeyError: 'generator'`` in the console, and a traceback on the page.

    Community Cloud swaps the files under a running server and Streamlit drops
    the changed modules from ``sys.modules``; a run that imports one in that
    window raised out of the first lines of `app.py`. The run must instead end
    on the "updating" notice, with the exception named underneath so that a
    file that is really missing is not mistaken for an update forever.

    Simulated by removing ``generator`` from ``sys.modules`` and making its
    next import raise, which is what the swap does. Not reproduced: the race
    itself, against a real deploy.
    """
    import sys

    monkeypatch.delitem(sys.modules, "generator", raising=False)
    monkeypatch.setattr(sys, "meta_path",
                        [_BrokenImport("generator", error)] + sys.meta_path)
    at = AppTest.from_file(APP, default_timeout=60)
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    assert any("updating" in str(i.value) for i in at.info), \
        "no 'The app is updating' notice"
    assert any(type(error).__name__ in str(c.value) for c in at.caption), \
        "the notice does not name the exception underneath"
    assert not at.title, "the run went on to draw a screen"


def test_the_footer_reports_the_running_versions(app):
    """The fingerprint covers the files, this line covers what runs them."""
    import platform

    lines = [str(c.value) for c in app.caption
             if str(c.value).startswith("Running on ")]
    assert len(lines) == 1, f"expected one 'Running on' line, found {lines}"
    assert f"Python {platform.python_version()}" in lines[0]
    for name in ("Streamlit", "Biopython", "ANARCI"):
        assert name in lines[0], f"{name} is missing from {lines[0]!r}"


def test_sections_are_grouped_on_the_home_page_and_in_the_sidebar(app):
    """The group titles are drawn, once in the sidebar and once on the home
    page, and the section buttons follow the grouping in both places.

    The order of the buttons is what tells a user which group a section is
    in, so it is checked, not only the presence of the titles.
    """
    import shell

    # Not _goto: the home screen's title is the app name, not "🏠 Home",
    # and the Home button is disabled while it is showing.
    if EXPECTED_TITLE["🏠 Home"] not in [t.value for t in app.title]:
        _click(app, "🏠 Home")
    _no_failure(app, "home")
    grouped = shell.sections_by_group()

    for where, tree in (("sidebar", app.sidebar), ("home", app.main)):
        headers = [m.value for m in tree.markdown
                   if "border-left:5px solid" in m.value]
        names = [n for n in grouped if any(f">{n}<" in h for h in headers)]
        assert names == list(grouped), f"{where}: group titles {names}"

    expected = [f"{shell.SECTIONS[k][0]} {shell.SECTIONS[k][1]}"
                for keys in grouped.values() for k in keys]
    sidebar = [b.label for b in app.sidebar.button if b.label in expected]
    assert sidebar == expected, f"sidebar order {sidebar}"
    home = [b.key for b in app.main.button if b.key.startswith("home_")
            and b.key != "home_manual"]
    assert home == [f"home_{k}" for keys in grouped.values() for k in keys], \
        f"home card order {home}"


def test_the_home_page_asks_how_the_app_is_used(app):
    """The invitation is drawn on the home page and offers both routes."""
    import shell

    if EXPECTED_TITLE["🏠 Home"] not in [t.value for t in app.title]:
        _click(app, "🏠 Home")
    _no_failure(app, "home")
    labels = [e.label for e in app.main.expander]
    assert any("Using it in your work" in label for label in labels), labels
    assert any(shell.usage_mailto() in m.value for m in app.main.markdown), \
        "the home page no longer offers the e-mail route"
    footer = [str(c.value) for c in app.caption if "report a problem" in str(c.value)]
    assert len(footer) == 1 and shell.usage_url() in footer[0], \
        "the footer no longer links to the usage report"


def test_the_culture_plan_divides_a_target_by_a_measured_yield(app):
    """The view with its inputs filled, not only its empty state.

    Arithmetic checked by hand: 10 mg at 4 mg per litre is 2.5 litres, and at
    a working volume of 1 L that is three flasks, because half a flask is a
    flask. The empty state is checked first, because "there is no default
    yield" is the claim this view rests on and an empty box is what enforces
    it.
    """
    from sections import bench as bench_screen

    _goto(app, "⚗️ Batch & QC")
    [r for r in app.radio if r.label == "View"][0].set_value(
        bench_screen.V_PLAN).run()
    _no_failure(app, "culture plan, empty")
    assert any("no target" in info.value for info in app.info), \
        "the empty plan did not say that it needs a target"

    # A target and no yield: the state where a tool that guessed would guess.
    app.session_state["cp_mg"] = 10.0
    app.run()
    _no_failure(app, "culture plan, target only")
    assert any("nothing here can predict one" in info.value
               for info in app.info), \
        "with a target and no yield, the plan did not say it cannot predict one"
    assert not [m for m in app.metric if m.label == "Culture needed"], \
        "a culture volume was shown before any yield was given"

    app.session_state["cp_yield"] = 4.0
    app.session_state["cp_flask"] = 1000.0
    app.run()
    _no_failure(app, "culture plan, filled")

    volume = [m for m in app.metric if m.label == "Culture needed"]
    assert volume, "the culture volume metric did not render"
    assert "2.50 L" in volume[0].value, volume[0].value
    flasks = [m for m in app.metric if m.label == "Flasks"]
    assert flasks and "3" in str(flasks[0].value), \
        f"expected three flasks, got {flasks[0].value if flasks else None}"

    # The fixture is shared, so the inputs go back as they were found.
    app.session_state["cp_mg"] = None
    app.session_state["cp_yield"] = None
    app.session_state["cp_flask"] = 0.0
    app.run()


def test_the_lot_comparison_says_what_it_needs_before_any_file(app):
    from sections import bench as bench_screen

    _goto(app, "⚗️ Batch & QC")
    [r for r in app.radio if r.label == "View"][0].set_value(
        bench_screen.V_LOTS).run()
    _no_failure(app, "across lots")
    assert any("Batch sheet" in info.value for info in app.info), \
        "the empty comparison did not say where the records come from"


def test_the_session_tab_offers_the_run_as_a_file_and_a_paragraph(app):
    """The session file and the methods paragraph, on a real run.

    Streamlit executes every tab body of the selected family, so selecting
    Outputs draws the Session tab too. What is checked is what the feature is
    for: a file to download, and a paragraph that names the annotator that
    actually ran rather than the setting, and says the values are
    theoretical.
    """
    import shell

    _goto(app, "🔬 Sequence analysis")
    if not [r for r in app.radio if r.label == "Section"]:
        app.text_area[0].set_value(FASTA)
        _click(app, "Analyze")
    families = [r for r in app.radio if r.label == "Section"][0]
    outputs = [option for option in families.options if "Outputs" in option][0]
    families.set_value(outputs).run()
    _no_failure(app, "analysis outputs")

    labels = [button.label for button in app.button] + [
        element.label for element in app.get("download_button")]
    assert any("Download the session" in label for label in labels), \
        f"no session download button among {labels}"

    paragraphs = [block.value for block in app.code]
    assert paragraphs, "the methods paragraph did not render"
    text = "\n".join(paragraphs)
    assert "theoretical" in text
    assert shell.APP_NAME in text
    # Without HMMER the run falls back to the regex anchors, and the paragraph
    # has to say which one ran either way: that is the claim it exists for.
    assert ("Regex anchors" in text) != HAS_HMMER or "ANARCI" in text


def test_provenance_lists_the_constants_and_what_they_are_worth(app):
    """The register, both views.

    The values are read from the modules while the page is drawn, so this is
    also the test that would catch an entry naming a constant that a refactor
    renamed: the screen would raise instead of drawing.
    """
    import provenance as pv

    _click(app, "🔎 Provenance")
    _no_failure(app, "provenance section")
    assert any(f"{len(pv.entries())} constants" in str(caption.value)
               for caption in app.caption), \
        "the register did not say how many constants it holds"
    assert app.dataframe, "the register table did not render"

    selector = [r for r in app.radio if r.label == "View"][0]
    for view in selector.options:
        [r for r in app.radio if r.label == "View"][0].set_value(view).run()
        _no_failure(app, f"provenance view {view}")
        assert _content_count(app) >= 6, f"provenance view {view} looks empty"

    # The second view groups by standing, and the group that matters most is
    # the one whose source could not be traced.
    body = "\n".join(str(element.value) for element in app.markdown)
    assert "unverified" in body, \
        "the grouped view does not show the unverified constants"


def test_a_pending_setting_reaches_the_sidebar_widget(app):
    """The half of the session restore that can be driven from a test.

    `AppTest` has no way to hand a file to `st.file_uploader`, so the upload
    itself is exercised by hand and not here. What *is* testable is the
    mechanism it depends on, and it is the part that fails silently: the
    screen parks a value under `pending_set_ph`, the router moves it into the
    widget's state before any widget is drawn, and the slider comes up on the
    restored value. If the key were misspelt in either file, the session would
    reload, the sequences would come back, and the pH would quietly be
    whatever the slider was already on.
    """
    _goto(app, "🔬 Sequence analysis")
    before = app.session_state["set_ph"]
    assert before != 6.0, "pick a different pH for this test"

    app.session_state["pending_set_ph"] = 6.0
    app.session_state["pending_set_window"] = 15
    app.run()
    _no_failure(app, "a restored setting")

    assert app.session_state["set_ph"] == 6.0
    assert app.session_state["set_window"] == 15
    assert "pending_set_ph" not in app.session_state, \
        "the pending key was not consumed, so it would be reapplied forever"
    sliders = {slider.label: slider.value for slider in app.slider}
    assert sliders["pH for net charge"] == 6.0
    assert sliders["Profile window (residues)"] == 15

    app.session_state["pending_set_ph"] = before
    app.run()


def _tab(app, label):
    return next(tab for tab in app.tabs if tab.label == label)


def _tick(app, prefix, names):
    """Tick rows of a library table by name.

    `AppTest` cannot click a row of `st.dataframe`; the selection is put in
    the widget's state, where a click in the browser would put it. The rows
    are found by name in the table as drawn, so this goes through the same
    position-to-name mapping as a real click.
    """
    table = next(df for df in app.dataframe
                 if str(df.key).startswith(f"{prefix}_table_"))
    shown = list(table.value["name"])
    rows = [shown.index(name) for name in names]
    app.session_state[table.key] = {"selection": {"rows": rows,
                                                  "columns": [], "cells": []}}
    app.run()
    return app


def test_the_library_has_a_public_and_a_private_tab(app):
    _goto(app, "📚 Sequence library")
    _no_failure(app, "library")
    assert [t.label for t in app.tabs][:2] == ["Public library", "My library"]
    mine = _tab(app, "My library")
    assert "Start an empty library" in [b.label for b in mine.button]
    assert any("held in this browser session only" in i.value
               for i in mine.info)


def test_public_entries_can_be_copied_into_my_library(app):
    """The copy button creates the private library when there is none, and
    a second copy of the same name is refused rather than overwriting."""
    import library as lib
    _goto(app, "📚 Sequence library")
    # Start from no private library, whatever an earlier test left.
    for key in ("lib_mine_current", "lib_mine_base", "lib_mine_saved",
                "lib_mine_file"):
        if key in app.session_state:
            del app.session_state[key]
    app.run()

    public, _ = lib.load_library()
    names = list(public["name"][:2])
    _tick(app, "library", names)
    _click(app, "➕ Start my library with these 2")
    _no_failure(app, "copy to my library")

    mine = app.session_state["lib_mine_current"]
    assert list(mine["name"]) == names
    assert any("Copied to your library" in s.value for s in app.success)
    assert any("Changes not downloaded yet" in w.value for w in app.warning), \
        "a library that was never downloaded must say so"

    _click(app, "➕ Copy 2 to my library")
    assert list(app.session_state["lib_mine_current"]["name"]) == names
    assert any("left as they were" in w.value for w in app.warning)


def test_my_library_survives_a_visit_to_another_section(app):
    """Coming back to the screen must show the library as last edited.

    Two things can happen to the editor's own state while another screen is
    drawn. Streamlit documents that it deletes the state of a widget that is
    not drawn; under `AppTest` with Streamlit 1.66 it was kept instead. If it
    is kept, the edits are replayed onto the frame they were made against and
    nothing is lost. If it is deleted, the screen restarts the editor from
    what it last returned. This test covers the second case, by deleting it
    as Streamlit says it would; `AppTest` cannot type into a data editor, so
    the edit is put where the screen stores the editor's last output.
    """
    import library as lib
    _goto(app, "📚 Sequence library")
    _goto(app, "🔬 Sequence analysis")
    editor = f"lib_mine_editor_{app.session_state['lib_mine_gen']}"
    if editor in app.session_state:
        del app.session_state[editor]
    edited = lib.editable_frame(app.session_state["lib_mine_current"])
    # A real sequence, Nb80's, so the analysis it is sent to accepts it.
    public, _ = lib.load_library()
    nb80 = public.loc[public["name"] == "Nb80", "protein"].iloc[0]
    edited.loc[len(edited)] = ["mine_only", nb80, "", "", ""]
    app.session_state["lib_mine_current"] = edited
    _goto(app, "📚 Sequence library")
    _no_failure(app, "back to the library")
    assert "mine_only" in list(app.session_state["lib_mine_current"]["name"])


def test_my_library_sends_to_analysis_and_closes(app):
    _goto(app, "📚 Sequence library")
    # The ticks are a column of the editor, which `AppTest` cannot click.
    # They are put in what the editor last returned, and the editor is
    # restarted from it, the same path as coming back from another screen.
    mine = app.session_state["lib_mine_current"].copy()
    mine["send"] = mine["name"] == "mine_only"
    app.session_state["lib_mine_current"] = mine
    editor = f"lib_mine_editor_{app.session_state['lib_mine_gen']}"
    if editor in app.session_state:
        del app.session_state[editor]
    app.run()
    assert any(c.value == "Selected: mine_only" for c in app.caption)
    send = [b for b in _tab(app, "My library").button
            if b.label.startswith("🔬 Send 1 protein")]
    assert send, "no send button in the private tab"
    send[0].click().run()
    _no_failure(app, "send from my library")
    assert any(t.value == "🔬 Sequence analysis" for t in app.title)

    _goto(app, "📚 Sequence library")
    closing = _tab(app, "My library")
    assert next(b for b in closing.button
                if b.label == "Close my library").disabled, \
        "closing with unsaved changes must ask first"
    app.checkbox(key="lib_mine_discard").check().run()
    _click(app, "Close my library")
    _no_failure(app, "close my library")
    assert "lib_mine_current" not in app.session_state


def test_a_filter_clears_the_ticks_rather_than_moving_them(app):
    """A selection is stored as row positions. Typing a filter changes which
    entry sits at each position, so the table is drawn under a new key and
    starts unticked, instead of quietly sending other sequences."""
    import library as lib
    _goto(app, "📚 Sequence library")
    public, _ = lib.load_library()
    last = public["name"].iloc[-1]
    _tick(app, "library", [last])
    assert any(f"Selected: {last}" == c.value for c in app.caption)

    app.text_input(key="lib_search").input(public["name"].iloc[0]).run()
    _no_failure(app, "filtering")
    assert not any(c.value.startswith("Selected:") for c in app.caption)
    assert any("Tick the box" in c.value for c in app.caption)
    app.text_input(key="lib_search").input("").run()

