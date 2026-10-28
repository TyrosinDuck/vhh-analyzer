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

NAV = ["🔬 Sequence analysis", "🧊 Structure", "🧪 Sanger sequencing",
       "🔧 Gene diagnosis", "📚 Sequence library", "🧬 In-silico cloning",
       "📖 Manual", "🏠 Home"]


def _labels(app):
    return [button.label for button in app.button]


def _click(app, label):
    app.button[_labels(app).index(label)].click().run()
    return app


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
EXPECTED_TITLE = {
    "🔬 Sequence analysis": "🔬 Sequence analysis",
    "🧊 Structure": "🧊 Structure",
    "🧪 Sanger sequencing": "🧪 Sanger sequencing",
    "🔧 Gene diagnosis": "🔧 Gene diagnosis",
    "📚 Sequence library": "📚 Sequence library",
    "🧬 In-silico cloning": "🧬 In-silico cloning",
    "🏠 Home": "🧰 VHH tool box",
    "📖 Manual": None,
}


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
