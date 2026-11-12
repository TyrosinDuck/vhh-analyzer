#!/usr/bin/env python3
"""
Build, or check, the worked example in `examples/`.

Why this exists. The test suite checks that the code does what the code says:
that a parser raises on malformed input, that a fingerprint moves when a file
moves, that a caption and a constant agree. None of it pins a *number*. A
refactor that quietly changed every molar mass by a tenth of a percent would
pass the whole suite, and a reviewer asking "how do I know this gives the
right answer" has nothing to run.

So: three published sequences, one set of settings, and the descriptors they
produce, committed as a file. The suite recomputes them and compares.

**What this is not.** The expected values were produced by this code. They are
not an independent reference, and this example cannot tell you that a molar
mass is correct, only that it has not changed. Checking correctness means
comparing against something outside this repository, which for the ProtParam
descriptors means Expasy, and that is a separate exercise, done by hand,
recorded in `examples/README.md` where it has been done at all.

What it is good for, and it is worth having for this alone: a dependency bump
that moves a number stops being invisible. Biopython computes the mass, the
pI and the extinction coefficients, and `requirements.txt` pins it; the day
that pin moves, this is what says whether the numbers moved with it.

The annotator is the regex one, deliberately. ANARCI is installed beside the
app rather than by it, so a baseline that depended on it would pass or fail
according to whether HMMER happens to be on the machine, which is the one
thing a baseline must not do.

Usage:

    python3 tools/make_examples.py           # write examples/
    python3 tools/make_examples.py --check   # compare, exit 1 on a difference
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import library as lib                                        # noqa: E402
import vhh_core as vc                                        # noqa: E402

EXAMPLES = ROOT / "examples"
MANIFEST = EXAMPLES / "manifest.json"
FASTA = EXAMPLES / "sequences.fasta"
EXPECTED = EXAMPLES / "expected.json"

# How many decimals a float keeps in the committed file. Four is well past
# anything the app displays and well short of where floating point noise
# lives, so a difference at this precision is a difference in the code or in
# a pinned dependency, not in the hardware.
PLACES = 4

# The descriptors the baseline pins, which is every scalar `characterize`
# returns apart from the ones that cannot drift. `name` and `sequence` are
# inputs, and `seq_id` is a hash of the sequence, so pinning those three would
# pin nothing. Listed rather than derived from a run, so that a descriptor
# which stops being computed fails the check instead of quietly leaving the
# baseline.
PINNED = (
    "length", "mw_da", "mw_kda", "pI", "ext_coeff_reduced",
    "ext_coeff_cystines", "A280_0.1pct", "A280_0.1pct_cystines",
    "instability_index", "aliphatic_index", "gravy", "aromaticity",
    "net_charge", "charge_pH", "n_Cys", "n_Lys", "n_primary_amines",
    "n_Met", "n_Trp", "n_Asn", "helix_frac", "turn_frac", "sheet_frac",
)


def _round(value):
    return round(value, PLACES) if isinstance(value, float) else value


def manifest() -> dict:
    """What the example is made of: which sequences, and which settings."""
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def sequences(spec: dict) -> dict[str, str]:
    """The example sequences, taken from the library shipped with the app.

    Read from `data/library.csv` rather than written out here, so that the
    example cannot carry its own private copy of a sequence that the library
    later corrects. `examples/sequences.fasta` is written from this, and the
    suite checks the committed file still matches, which is what turns the
    derivation into something that fails loudly.
    """
    table, _warnings = lib.load_library()
    known = set(table["name"])
    missing = [name for name in spec["library_entries"] if name not in known]
    if missing:
        raise SystemExit(
            f"the manifest names sequences the library does not have: "
            f"{', '.join(missing)}")
    rows = table.set_index("name")
    return {name: vc.clean_sequence(rows.loc[name, "protein"])
            for name in spec["library_entries"]}


def compute(spec: dict) -> dict:
    """The expected file: descriptors, annotation and classification.

    Everything here goes through the same functions the app calls, with the
    settings in the manifest, so a change in any of them shows up.
    """
    settings = spec["settings"]
    out: dict[str, dict] = {}
    for name, seq in sequences(spec).items():
        properties = vc.characterize(name, seq, settings["ph"])
        annotation, used = vc.annotate(seq, settings["annotation_method"])
        kind = vc.classify(seq, settings["annotation_method"])
        liabilities = vc.find_liabilities(seq, annotation)
        hallmarks = vc.hallmarks(seq, annotation)
        out[name] = {
            "length": len(seq),
            "properties": {key: _round(properties[key]) for key in PINNED
                           if key in properties},
            "annotator": used,
            "regions": {region: [annotation[region]["start"],
                                 annotation[region]["end"]]
                        for region in vc.REGIONS if region in annotation},
            "kind": kind["kind"],
            "camelid_score": hallmarks["camelid_score"],
            "n_liabilities": int(len(liabilities)),
        }
    return {
        "note": ("Generated by tools/make_examples.py from the sequences in "
                 "data/library.csv. A baseline, not an independent reference: "
                 "see examples/README.md."),
        "settings": settings,
        "sequences": out,
    }


def write(spec: dict) -> None:
    FASTA.write_text(_fasta(sequences(spec)), encoding="utf-8")
    EXPECTED.write_text(json.dumps(compute(spec), indent=1, sort_keys=True) +
                        "\n", encoding="utf-8")
    print(f"wrote {FASTA.relative_to(ROOT)} and "
          f"{EXPECTED.relative_to(ROOT)}")


def _fasta(seqs: dict[str, str], width: int = 60) -> str:
    lines = []
    for name, seq in seqs.items():
        lines.append(f">{name}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + "\n"


def check(spec: dict) -> int:
    """Compare what is committed with what the code computes now."""
    problems = []
    if _fasta(sequences(spec)) != FASTA.read_text(encoding="utf-8"):
        problems.append(
            f"{FASTA.relative_to(ROOT)} is not what the library now holds")
    fresh = compute(spec)
    stored = json.loads(EXPECTED.read_text(encoding="utf-8"))
    if fresh["sequences"] != stored["sequences"]:
        for name, values in fresh["sequences"].items():
            was = stored["sequences"].get(name)
            if was != values:
                problems.append(f"{name} differs from the baseline")
    for problem in problems:
        print(f"differs: {problem}", file=sys.stderr)
    if problems:
        print("Run `python3 tools/make_examples.py` to update the baseline, "
              "once you know why it moved.", file=sys.stderr)
        return 1
    print(f"{len(fresh['sequences'])} example sequences match the baseline")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="compare instead of writing")
    args = parser.parse_args()
    spec = manifest()
    if args.check:
        return check(spec)
    write(spec)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
