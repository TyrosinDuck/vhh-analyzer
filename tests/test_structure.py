"""
Structure tests, the residue-number/sequential-index distinction, mostly.

The structure section has one source of subtle error that is worth more than all
the others put together: a model's **residue numbers** and a residue's
**position in the chain** are different things, and they differ the moment a
chain does not start at 1 or has an unmodelled loop, which is most of the time.
Confusing the two does not raise; it quietly puts the CDR bands, the region
column and the 3D highlighting in the wrong place.

The structures here are synthetic and built in-process: CA atoms on a line, with
a deliberate numbering offset and a deliberate gap. Their geometry is meaningless
and no test below asks anything of it. What is being tested is the bookkeeping,
and for that a synthetic chain is better than a real one, because the right
answer is known by construction.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Data import PDBData

import structure as stx
import vhh_core as vc

ROOT = Path(__file__).resolve().parent.parent

THREE_LETTER = {one: three for three, one
                in PDBData.protein_letters_3to1.items()}


def synthetic_pdb(sequence: str, first_resnum: int = 1,
                  missing: set[int] | None = None) -> bytes:
    """A CA-only PDB of ``sequence``, numbered from ``first_resnum``.

    ``missing`` is a set of 0-based indices to leave out, which creates a gap in
    the numbering exactly as an unmodelled loop does. Coordinates are a straight
    line: enough to parse, meaningless to measure.
    """
    missing = missing or set()
    lines, serial = [], 0
    for index, amino in enumerate(sequence):
        if index in missing:
            continue
        serial += 1
        three = THREE_LETTER[amino].upper()
        resnum = first_resnum + index
        lines.append(
            f"ATOM  {serial:5d}  CA  {three:>3s} A{resnum:4d}    "
            f"{index * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           C")
    return ("\n".join(lines) + "\nTER\nEND\n").encode()


@pytest.fixture(scope="module")
def cablys3():
    import library as lib
    table, problems = lib.load_library(ROOT / "data" / "library.csv")
    assert not problems
    hits = table[table["name"].str.contains("cAbLys3", case=False)]
    assert len(hits) == 1
    return hits.iloc[0]["protein"]


# --------------------------------------------------------------------------- #
# The bug this file exists for
# --------------------------------------------------------------------------- #

def test_resnum_map_handles_an_offset_and_a_gap(cablys3):
    """Bug: the mapping was keyed by sequential index, consumed as residue
    numbers.

    Every consumer, the CDR bands, the region column, the viewer's residue
    selection, speaks residue numbers. The mapping spoke positions. On a chain
    starting at residue 7, every annotation was off by six; with a gap, every
    annotation after the gap was off by more.

    Here the chain starts at 7 and eight residues are missing from CDR2, so the
    expected answer is known: residue number *n* maps to sequence index *n - 7*.
    """
    first, missing = 7, set(range(52, 60))
    model = stx.load_structure(
        synthetic_pdb(cablys3, first_resnum=first, missing=missing),
        "synthetic.pdb")

    chain_seq, residues = stx.chain_sequence(model.model["A"])
    assert len(chain_seq) == len(cablys3) - len(missing)
    assert residues[0].id[1] == first

    mapping, identity = stx.resnum_map(model, "A", chain_seq, cablys3)
    assert identity == pytest.approx(1.0), \
        "a chain built from the sequence itself must map at 100% identity"
    assert len(mapping) == len(chain_seq), "some residues went unmapped"

    for resnum, index in mapping.items():
        assert index == resnum - first, (
            f"residue {resnum} mapped to index {index}, expected "
            f"{resnum - first}, the offset is being lost")
        assert cablys3[index] == chain_seq[
            [r.id[1] for r in residues].index(resnum)], \
            f"residue {resnum} maps to the wrong amino acid"


def test_identity_map_is_not_the_identity_function(cablys3):
    """"Annotate against the structure's own sequence" is still a mapping.

    The tempting shortcut is ``{n: n}``. It is wrong: residue number 7 is
    sequence position 0, not 7. This is the same bug wearing a different hat.
    """
    model = stx.load_structure(synthetic_pdb(cablys3, first_resnum=7),
                               "synthetic.pdb")
    mapping = stx.identity_resnum_map(model, "A")
    assert mapping[7] == 0, "the offset was dropped"
    assert max(mapping) == 7 + len(cablys3) - 1
    assert list(mapping.values()) == list(range(len(cablys3)))


def test_regions_land_on_the_right_residue_numbers(cablys3):
    """The consequence of the bug, stated in terms a user would notice.

    The CDR3 of the sequence must be labelled CDR3 on the model residues that
    actually carry it, not on the residues six positions earlier.
    """
    first = 7
    model = stx.load_structure(synthetic_pdb(cablys3, first_resnum=first),
                               "synthetic.pdb")
    chain_seq, residues = stx.chain_sequence(model.model["A"])
    mapping, _ = stx.resnum_map(model, "A", chain_seq, cablys3)
    annotation, _ = vc.annotate(cablys3)

    resnums = [r.id[1] for r in residues]
    regions = stx.residue_regions(resnums, annotation, mapping)
    labelled = dict(zip(resnums, regions))

    start, end = annotation["CDR3"]["start"], annotation["CDR3"]["end"]
    for position in range(start, end + 1):
        resnum = position - 1 + first
        assert labelled[resnum] == "CDR3", (
            f"residue {resnum} carries CDR3 position {position} but is "
            f"labelled {labelled[resnum]!r}")
    assert labelled[start - 2 + first] != "CDR3", \
        "the residue before CDR3 is labelled CDR3, the bands are shifted"


def test_residue_regions_is_empty_without_an_annotation():
    """No annotation is a normal state, not an error."""
    assert stx.residue_regions([1, 2, 3], None, {1: 0}) == [None, None, None]
    assert stx.residue_regions([1, 2, 3], {"CDR3": {}}, None) == [None] * 3


# --------------------------------------------------------------------------- #
# Measurements, on a structure whose answer is known independently
# --------------------------------------------------------------------------- #

def test_liabilities_carry_the_model_residue_number(cablys3):
    """The column the 3D viewer selects on.

    Without it the exposed-liability highlighting had to re-derive the mapping
    in the UI, which is exactly where the index/resnum confusion came from.
    """
    first = 7
    model = stx.load_structure(synthetic_pdb(cablys3, first_resnum=first),
                               "synthetic.pdb")
    annotation, _ = vc.annotate(cablys3)
    liabilities = vc.find_liabilities(cablys3, annotation)
    if not len(liabilities):
        pytest.skip("no liabilities in this sequence")

    joined = stx.annotate_liabilities(liabilities, model, "A", cablys3, cablys3)
    assert "resnum" in joined.columns
    for _, row in joined.iterrows():
        if row["resnum"] is None:
            continue
        assert int(row["resnum"]) == int(row["position"]) - 1 + first, (
            f"liability at sequence position {row['position']} reports model "
            f"residue {row['resnum']}, expected {row['position'] - 1 + first}")


def test_unmodelled_residues_are_reported_as_such_not_as_buried(cablys3):
    """A residue the model does not cover has no accessibility.

    Calling it buried would be the worst possible failure mode here: it would
    silently downgrade a liability on the grounds that it is hidden, when in
    fact nothing is known about it.
    """
    missing = set(range(52, 60))
    model = stx.load_structure(synthetic_pdb(cablys3, missing=missing),
                               "synthetic.pdb")
    chain_seq, _residues = stx.chain_sequence(model.model["A"])
    annotation, _ = vc.annotate(cablys3)
    liabilities = vc.find_liabilities(cablys3, annotation)
    joined = stx.annotate_liabilities(liabilities, model, "A", chain_seq,
                                      cablys3)
    # Whatever lands in the deleted stretch must be "not modelled"
    in_gap = joined[joined["position"].astype(int).sub(1).isin(missing)]
    for _, row in in_gap.iterrows():
        assert row["exposure"] == "not modelled", (
            f"a liability in an unmodelled loop was called {row['exposure']!r}")


def test_annotate_liabilities_refuses_the_wrong_structure_sequence(cablys3):
    """Passing the full-length construct where the chain's sequence belongs.

    Found by this very test file getting it wrong. It did not raise: the mapping
    came back looking plausible and every position after the gap was quietly
    shifted, which is the worst kind of wrong. It is now refused.
    """
    missing = set(range(52, 60))
    model = stx.load_structure(synthetic_pdb(cablys3, missing=missing),
                               "synthetic.pdb")
    annotation, _ = vc.annotate(cablys3)
    liabilities = vc.find_liabilities(cablys3, annotation)
    with pytest.raises(stx.StructureError):
        stx.annotate_liabilities(liabilities, model, "A", cablys3, cablys3)


def test_empty_liabilities_still_gains_the_structure_columns(cablys3):
    """The same empty-table rule as in the sequence-only path."""
    model = stx.load_structure(synthetic_pdb(cablys3), "synthetic.pdb")
    empty = vc.find_liabilities("GSGSGSGSGSGSGSGSGSGSGSGSGSGSGS")
    joined = stx.annotate_liabilities(empty, model, "A", cablys3, cablys3)
    assert joined.empty
    for column in ("rel_sasa", "exposure"):
        assert column in joined.columns


def test_confidence_kind_is_a_guess_and_says_so(cablys3):
    """pLDDT and B-factor live in the same column and mean opposite things.

    Nothing in either file format states which it is, so the only honest
    behaviour is a labelled guess with a reason attached, never a silent
    assumption.
    """
    model = stx.load_structure(synthetic_pdb(cablys3), "synthetic.pdb")
    kind, reason = stx.confidence_kind(model)
    assert kind in ("pLDDT", "B-factor", "unknown")
    assert reason and isinstance(reason, str), \
        "a guess with no stated reason is not a guess, it is an assumption"


def test_a_file_that_is_not_a_structure_raises_cleanly():
    """Garbage in must produce StructureError, not a stack trace in the UI."""
    with pytest.raises(stx.StructureError):
        stx.load_structure(b"this is not a PDB file\n", "nope.pdb")


def test_patch_threshold_is_a_named_constant():
    """The default was tuned on a single complex and must stay visible.

    It is referenced by name in the app and documented in the module, so that
    nobody later reads 0.15 as a validated parameter.
    """
    assert 0.0 < stx.PATCH_MIN_REL < stx.EXPOSED_MIN, (
        "the patch threshold must sit below the 'exposed' cut-off, a residue "
        "inside a patch is half-covered by the patch itself")
