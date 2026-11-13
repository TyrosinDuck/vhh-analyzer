"""What a 1.x release promises not to break, checked rather than stated.

The README's *Public API and versioning* section says what is public: the
names in `vhh_core.__all__`, and the three files the app writes for a user to
keep and later reads back (a batch record, a session file, a sequence
library). Breaking any of them is a MAJOR release. These tests are what turns
that sentence into something a push can fail on.

Two kinds of check:

* `vhh_core` against `tests/compat/vhh_core_api.json`, a snapshot written by
  `tools/public_api.py`. A name may be added, a parameter with a default may
  be added; nothing recorded may disappear, move, lose its default or change
  kind. The snapshot is refreshed only at a MAJOR release or to protect what a
  MINOR release added, never to make this test pass.
* every file in `tests/compat/` that an earlier release wrote must still be
  read by this one. The files are found by walking the folder rather than
  listed here, so a new schema number is covered by dropping a file of the
  old one beside it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import bench
import library
import session
import shell
import vhh_core

COMPAT = Path(__file__).resolve().parent / "compat"
SNAPSHOT = json.loads((COMPAT / "vhh_core_api.json").read_text(encoding="utf-8"))

POSITIONAL = {"POSITIONAL_ONLY", "POSITIONAL_OR_KEYWORD"}


def test_the_snapshot_belongs_to_this_major_release():
    """A snapshot from another MAJOR is not a promise this release made."""
    major = int(shell.RELEASE.split(".")[0])
    assert SNAPSHOT["release_major"] == major, (
        f"tests/compat/vhh_core_api.json was taken under release "
        f"{SNAPSHOT['taken_at_release']} and this is {shell.RELEASE}. After a "
        f"MAJOR release, run `python3 tools/public_api.py` to record the new "
        f"API as the one 2.x promises to keep.")


@pytest.mark.parametrize("name", sorted(SNAPSHOT["names"]))
def test_a_public_name_is_still_there(name):
    assert name in vhh_core.__all__, (
        f"vhh_core.{name} was public and is no longer in __all__. Removing it "
        f"breaks code that imports it: that is a MAJOR release, or keep it "
        f"with a DeprecationWarning until the next one.")
    assert hasattr(vhh_core, name)


def _params(obj):
    import inspect
    return [{"name": p.name, "kind": p.kind.name,
             "default": p.default is not inspect.Parameter.empty}
            for p in inspect.signature(obj).parameters.values()]


@pytest.mark.parametrize("name", sorted(n for n, v in SNAPSHOT["names"].items()
                                        if v["callable"]))
def test_a_public_signature_still_accepts_every_call_it_used_to(name):
    old = SNAPSHOT["names"][name]["parameters"]
    new = _params(getattr(vhh_core, name))
    by_name = {p["name"]: p for p in new}
    where = f"vhh_core.{name}"

    for p in old:
        assert p["name"] in by_name, (
            f"{where} lost its parameter `{p['name']}`; a call passing it now "
            f"fails. That is a MAJOR change.")
        q = by_name[p["name"]]
        assert q["kind"] == p["kind"], (
            f"{where}: `{p['name']}` was {p['kind']} and is now {q['kind']}, "
            f"which breaks calls that pass it the old way.")
        if p["default"]:
            assert q["default"], (
                f"{where}: `{p['name']}` had a default and no longer does, so "
                f"a call that left it out now fails.")

    old_positional = [p["name"] for p in old if p["kind"] in POSITIONAL]
    new_positional = [p["name"] for p in new if p["kind"] in POSITIONAL]
    assert new_positional[:len(old_positional)] == old_positional, (
        f"{where}: positional parameters were {old_positional} and are now "
        f"{new_positional}. A call passing them by position would hand each "
        f"value to the wrong parameter, silently.")

    known = {p["name"] for p in old}
    for q in new:
        if q["name"] in known or q["kind"] in ("VAR_POSITIONAL", "VAR_KEYWORD"):
            continue
        assert q["default"], (
            f"{where} gained `{q['name']}` without a default; every existing "
            f"call now fails for want of it. Give it a default (MINOR) or "
            f"make this a MAJOR release.")


# --------------------------------------------------------------------------- #
# Files a user keeps
# --------------------------------------------------------------------------- #

# What reads a file, by the family named in its `schema` key. A family is the
# part before the slash; the number after it is the format version, which is
# the thing allowed to move as long as the old numbers stay readable.
READERS = {
    bench.BATCH_RECORD_SCHEMA.rpartition("/")[0]: bench.batch_record_inputs,
    session.SESSION_SCHEMA.rpartition("/")[0]: session.session_inputs,
}

KEPT_JSON = sorted(COMPAT.glob("*.json"))
KEPT_JSON = [p for p in KEPT_JSON if p.name != "vhh_core_api.json"]


def test_every_kind_of_kept_json_file_has_a_frozen_example():
    """A schema with no frozen file would be a promise nothing checks."""
    families = {json.loads(p.read_text(encoding="utf-8"))["schema"]
                .rpartition("/")[0] for p in KEPT_JSON}
    assert families == set(READERS), (
        f"frozen examples cover {sorted(families)}, readers exist for "
        f"{sorted(READERS)}. Save one file the app writes into tests/compat/.")


@pytest.mark.parametrize("path", KEPT_JSON, ids=lambda p: p.name)
def test_a_file_an_earlier_release_wrote_is_still_read(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    family = data["schema"].rpartition("/")[0]
    back = READERS[family](data)   # raises if this release refuses it
    assert back, f"{path.name} was read back as nothing"


def test_the_frozen_batch_record_still_gives_its_readings():
    data = json.loads((COMPAT / "batch-record_1.json").read_text("utf-8"))
    back = bench.batch_record_inputs(data)
    assert back["fields"]["lot"] == "JV-2026-07"
    assert back["measured"]["a280"] == pytest.approx(0.84)
    assert back["measured"]["measured_da"] == pytest.approx(13650.2)
    assert back["sequence"] == data["sequence"]


def test_the_frozen_session_still_gives_its_settings():
    data = json.loads((COMPAT / "session_1.json").read_text("utf-8"))
    back = session.session_inputs(data)
    assert back["sequences"] == data["sequences"]
    assert back["ph"] == pytest.approx(7.4)
    assert back["window"] == 9


def test_the_frozen_library_is_read_without_a_problem():
    table, problems = library.load_library(COMPAT / "library.csv")
    assert problems == []
    assert list(table["name"]) == ["NbX", "pelB leader"]
    assert table.loc[0, "note"] == "protein only, with a comma in the note"
    assert set(library.COLUMNS) <= set(table.columns)
