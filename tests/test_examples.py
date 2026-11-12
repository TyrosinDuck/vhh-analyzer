"""
The worked example: does the app still give the numbers it gave before, and is
one of them right?

Two different questions, and the file is in two halves because mixing them up
is the usual way a "validation suite" comes to validate nothing.

**Has it changed.** `examples/expected.json` holds descriptors this code
produced for three published sequences. Recomputing them and comparing catches
a refactor that moved a number and a dependency bump that moved one, which
nothing else in the suite does: every other test checks that the code does
what the code says, and a change of constant is perfectly consistent with
that. It is a baseline, not a reference. It cannot tell anyone that a molar
mass is correct, only that it is the same.

**Is it right.** One descriptor can be checked against something outside
Biopython without leaving the repository, and it is the one most results rest
on. The molar mass is recomputed from the per-residue formulas in
`data/aminoacids/meta.json` and the CIAAW atomic weights in `buffers.py`, two
sources that have nothing to do with Biopython's weight table, and compared
with what the app reports. On the three example sequences the two routes agree
to within 14 parts per million, which is the rounding in the published
free-acid masses rather than a disagreement about chemistry.

No equivalent check exists here for the isoelectric point or the extinction
coefficients. The pI comes from an iterative solve over pKa values, the
extinction coefficient from the Pace tabulation of Trp, Tyr and Cys
absorbances, and this repository holds no second source for either. Checking
those means comparing against Expasy ProtParam by hand; `examples/README.md`
records where that has been done and where it has not.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

import buffers as bf                                         # noqa: E402
import library as lib                                        # noqa: E402
import make_examples as mk                                    # noqa: E402
import vhh_core as vc                                         # noqa: E402

EXAMPLES = ROOT / "examples"

# Measured, not chosen: the two mass routes differ by 14 ppm at most on these
# sequences, and the tolerance is set a little above that so a rounding
# difference does not fail the suite while a real divergence does. A hundred
# ppm on a 14 kDa VHH is 1.4 Da, which is smaller than any modification the
# mass check screen looks for, so nothing useful hides under it.
MASS_TOLERANCE_PPM = 100.0


@pytest.fixture(scope="module")
def spec():
    return mk.manifest()


def test_the_manifest_names_sequences_the_library_has(spec):
    table, _warnings = lib.load_library()
    missing = set(spec["library_entries"]) - set(table["name"])
    assert not missing, f"the example names sequences the library lost: {missing}"


def test_the_committed_fasta_is_what_the_library_now_holds(spec):
    """The example sequences are derived from the library, not a copy of it.

    If they were a copy, a correction to `data/library.csv` would leave the
    example quietly describing the old sequence. The FASTA is committed so a
    reviewer can read it without running anything, and this is what keeps the
    committed file honest.
    """
    committed = (EXAMPLES / "sequences.fasta").read_text(encoding="utf-8")
    assert committed == mk._fasta(mk.sequences(spec)), (
        "examples/sequences.fasta no longer matches data/library.csv. Run "
        "`python3 tools/make_examples.py` once you know why it moved.")


def test_the_descriptors_still_match_the_baseline(spec):
    """The regression half: every pinned descriptor, recomputed.

    A failure here is not necessarily a defect. It says a number moved, and
    the two things that move numbers are a change in this code and a change in
    a pinned dependency, Biopython above all. Find out which before updating
    the baseline, and say so in CHANGELOG.md if it was the second: a result
    that silently changed with an upgrade is exactly what a user citing a
    version DOI needs to know about.
    """
    stored = json.loads((EXAMPLES / "expected.json").read_text(encoding="utf-8"))
    fresh = mk.compute(spec)
    assert fresh["settings"] == stored["settings"]
    assert set(fresh["sequences"]) == set(stored["sequences"])
    for name, values in fresh["sequences"].items():
        assert values == stored["sequences"][name], (
            f"{name} differs from examples/expected.json")


def test_the_baseline_does_not_depend_on_anarci_being_installed(spec):
    """The annotator in the example is the regex one, on purpose.

    ANARCI is installed beside the app rather than by it, so a baseline
    annotated with it would pass or fail according to whether HMMER happens to
    be on the machine running the suite. A test that passes on some machines
    and not others teaches nobody anything.
    """
    assert spec["settings"]["annotation_method"] == "regex"
    stored = json.loads((EXAMPLES / "expected.json").read_text(encoding="utf-8"))
    for name, values in stored["sequences"].items():
        assert values["annotator"] == "Regex anchors", name


@pytest.fixture(scope="module")
def residue_masses() -> dict[str, float]:
    """Neutral residue masses, built from the repository's own two sources.

    `data/aminoacids/meta.json` carries each amino acid's formula in the form
    it takes at pH 7.4, so a charged side chain is written with the proton it
    has lost or gained: glutamate is `C5H8NO4-` and lysine `C6H15N2O2+`. The
    neutral free acid is that formula with the charge undone, one hydrogen per
    unit of charge, and its mass comes from the CIAAW atomic weights in
    `buffers.py` through the same formula parser the buffer recipes use.

    Nothing in this path touches Biopython, which is the point.
    """
    meta = json.loads(
        (ROOT / "data" / "aminoacids" / "meta.json").read_text(encoding="utf-8"))
    rows = meta["aminoacids"]
    rows = list(rows.values() if isinstance(rows, dict) else rows)
    hydrogen = bf.ATOMIC_WEIGHTS["H"]
    masses = {}
    for row in rows:
        neutral = re.sub(r"[+-]+$", "", row["formula"])
        masses[row["one"]] = bf.molar_mass(neutral) - row["charge"] * hydrogen
    assert len(masses) == 20
    return masses


def test_the_formula_route_reproduces_the_published_free_acid_masses(
        residue_masses):
    """The independent route, checked against the file it came from.

    `meta.json` publishes a mass per amino acid as well as a formula, and the
    two are independent of each other: one was typed, the other is parsed and
    summed. They agree exactly for eighteen residues. Methionine and cysteine
    differ by 0.007 Da, which is the revision of sulfur's standard atomic
    weight, so the tolerance here is a hundredth of a dalton and the
    disagreement is named rather than hidden.
    """
    meta = json.loads(
        (ROOT / "data" / "aminoacids" / "meta.json").read_text(encoding="utf-8"))
    rows = meta["aminoacids"]
    rows = list(rows.values() if isinstance(rows, dict) else rows)
    for row in rows:
        assert residue_masses[row["one"]] == pytest.approx(
            row["mw_free_acid"], abs=0.01), row["one"]


def test_the_molar_mass_agrees_with_an_independent_calculation(
        spec, residue_masses):
    """The validation half: the app's mass, against formulas and CIAAW weights.

    A protein's average mass is the sum of its residues' free-acid masses less
    one water per peptide bond. The app's figure comes from Biopython's own
    residue weight table; this one comes from the formulas and atomic weights
    in this repository. Two routes, no shared constant, and they have to meet.
    """
    water = bf.molar_mass("H2O")
    for name, sequence in mk.sequences(spec).items():
        independent = (sum(residue_masses[residue] for residue in sequence)
                       - (len(sequence) - 1) * water)
        reported = vc.characterize(name, sequence, 7.4)["mw_da"]
        difference_ppm = abs(independent - reported) / reported * 1e6
        assert difference_ppm < MASS_TOLERANCE_PPM, (
            f"{name}: the app says {reported:.2f} Da, the formulas and CIAAW "
            f"atomic weights say {independent:.2f} Da, which is "
            f"{difference_ppm:.0f} ppm apart")


def test_the_example_is_regenerated_by_its_own_tool(spec, capsys):
    """`tools/make_examples.py --check` is the command a contributor runs.

    Driven here so that the tool cannot rot while the test that duplicates its
    arithmetic keeps passing, which would leave the documented command broken
    and nothing saying so.
    """
    assert mk.check(spec) == 0
    assert "match the baseline" in capsys.readouterr().out
