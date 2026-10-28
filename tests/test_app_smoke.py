"""
Smoke tests — does every screen of the app actually render?

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
slower than the rest of the suite — roughly a minute, against seven seconds for
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


def test_every_section_renders(app):
    """The version banner calls ``st.stop()``, so a mismatch shows up here too."""
    for label in NAV:
        _click(app, label)
        _no_failure(app, label)


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
    deliberately tiny — nothing in this test reads the geometry.
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
    another section and come back — which is worse than an error, because
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
