"""
The provenance register: does it show the app's constants, or its own copy?

The register is only worth having if it cannot drift. A page that listed
values typed next to their sources would be a second copy of every constant in
the app, and the first time one of them changed the page would start lying
with no symptom at all: a plausible number beside a real citation is exactly
what the November 2026 audit was cleaning up.

So the register names a module and an attribute, and the value is read live.
These tests are what makes that structural rather than a habit: one of them
parses `provenance.py` and fails if a numeric literal appears inside an entry.

What they cannot check, and it is worth being plain about: whether the source
sentence still describes the value. That needs a person reading the two
together.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

import bench as bn
import buffers as bf
import provenance as pv

ROOT = Path(__file__).resolve().parent.parent


def test_every_entry_resolves_to_a_live_value():
    """An entry naming an attribute that no longer exists fails here.

    This is the failure mode a register has: a constant gets renamed or moved
    during a refactor, and the page either crashes for a user or, worse, is
    quietly built from a stale copy. Resolving every entry is cheap and does
    both jobs.
    """
    for entry in pv.entries():
        value = pv.resolve(entry)
        assert value is not None, entry.key
        assert pv.value_text(value), entry.key


def test_no_entry_types_a_value_of_its_own():
    """The structural guarantee: no numeric literal inside a `Constant(...)`.

    Parsed rather than grepped, so that a digit inside a source sentence
    (there are several, "0.20", "14 ppm", "7249") does not trip it: a string
    is not a numeric literal. What this forbids is `value=6.0e16` or any other
    way of writing a number into an entry, which is the only way this file
    could come to disagree with the app.
    """
    tree = ast.parse((ROOT / "provenance.py").read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree)
             if isinstance(node, ast.Call)
             and getattr(node.func, "id", "") == "Constant"]
    assert len(calls) >= 15, "the registry stopped being built from Constant()"
    for call in calls:
        for node in ast.walk(call):
            if isinstance(node, ast.Constant) and isinstance(
                    node.value, (int, float)) and not isinstance(
                        node.value, bool):
                key = next((kw.value.value for kw in call.keywords
                            if kw.arg == "key"), "?")
                pytest.fail(
                    f"the entry {key!r} carries the numeric literal "
                    f"{node.value!r}. Name the module and attribute instead, "
                    f"so the value is read from the code that uses it")


def test_every_entry_says_what_the_number_is_worth():
    for entry in pv.entries():
        assert entry.basis in pv.BASES, f"{entry.key}: {entry.basis}"
        assert entry.source.strip(), f"{entry.key} has no source"
        assert entry.used_for.strip(), f"{entry.key} says nothing about where "\
                                       f"it is used"
        assert entry.label.strip()


def test_the_keys_are_unique():
    keys = [entry.key for entry in pv.entries()]
    assert len(keys) == len(set(keys))


def test_grouping_by_basis_loses_nothing():
    grouped = pv.by_basis()
    flat = [entry.key for entries in grouped.values() for entry in entries]
    assert sorted(flat) == sorted(entry.key for entry in pv.entries())
    assert list(grouped) == [basis for basis in pv.BASES if basis in grouped]
    assert all(grouped.values()), "a basis with no constant should be left out"


def test_the_reagents_are_derived_from_the_bench_table():
    """Adding a reagent to `bench.REAGENTS` adds a row to the register.

    The reagent molar masses already sit beside the source each was checked
    against, in the table the conjugation plan reads. Writing them out again
    here would be the second list this file exists to avoid, so they are read
    from there, and this is the test that says so.
    """
    rows = {entry.key for entry in pv.entries()}
    for name in bn.REAGENTS:
        assert f"reagent:{name}" in rows
    for entry in pv.entries():
        if entry.key.startswith("reagent:"):
            name = entry.key.split(":", 1)[1]
            assert entry.source == bn.REAGENTS[name]["source"]
            assert pv.resolve(entry) == bn.REAGENTS[name]["mw"]
            assert entry.basis == "vendor"


def test_the_atomic_weights_take_their_source_from_buffers():
    """`buffers` publishes its own source string, so the entry reads it."""
    entry = next(e for e in pv.entries() if e.key == "atomic_weights")
    assert entry.source == bf.ATOMIC_WEIGHTS_SOURCE
    assert "CIAAW" in entry.source


def test_the_unverified_constants_are_still_called_unverified():
    """The phage optical constant is the one with no traceable source.

    It is used because the arithmetic is standard and the app shows a second
    route beside it. The honest label is the whole point, and a later edit
    that promoted it to "published" without finding a source would be a
    regression this catches.
    """
    unverified = {entry.key for entry in pv.entries()
                  if entry.basis == "unverified"}
    assert "phage_optical_constant" in unverified
    entry = next(e for e in pv.entries()
                 if e.key == "phage_optical_constant")
    assert "provenance is not established" in entry.source


def test_the_conventions_do_not_claim_to_be_measurements():
    """A convention is a choice, and the register has to say so.

    Four of the structure thresholds and the panning input factor are
    conventions. If one of them were relabelled, a reader would take a
    threshold somebody picked for a value somebody measured, which is the
    distinction this page exists to draw.
    """
    conventions = {entry.key for entry in pv.entries()
                   if entry.basis == "convention"}
    for key in ("buried_max", "exposed_min", "ss_max_distance",
                "liquid_panning_factor", "apr_threshold"):
        assert key in conventions, key


def test_the_observed_constants_name_what_they_were_observed_on():
    """One structure is not a validation set, and the text has to say which.

    The patch exposure cut-off came from a single complex. The register is the
    place a user meets that number without the comment beside it, so the
    source has to carry the limitation rather than the code alone.
    """
    for entry in pv.entries():
        if entry.basis == "observed here":
            assert any(word in entry.source
                       for word in ("single", "one complex", "test genes")), \
                f"{entry.key} does not say what it was observed on"


def test_the_basis_vocabulary_is_explained_on_the_page():
    """Every kind carries a sentence, because the label alone decides nothing.

    "convention" means something to whoever wrote it and nothing to a reader,
    and the screen draws these sentences under each group.
    """
    for basis, explanation in pv.BASES.items():
        assert explanation.strip().endswith("."), basis
        assert len(explanation) > 40, basis
