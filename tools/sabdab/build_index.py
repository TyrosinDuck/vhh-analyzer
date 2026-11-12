"""
Turn a SAbDab2 single-domain export into the small index the app ships.

The export is a 53-column summary with one row per *chain instance*: the same
nanobody appears once for every chain of every structure it is in, so 4759
rows carry about 1390 distinct sequences. Shipping it whole would put three
megabytes of columns the app never reads into a container that may have 690
MB, and would make the file impossible to review in a diff.

So this script derives the index instead, and the derivation is the point:

* **One row per distinct VH sequence.** Every structure that contains it is
  kept, as a list, ordered best-resolution first. That is the shape the
  question has: "has anything like this been solved" wants the sequence once
  and the structures under it, not the sequence repeated.
* **Only the columns the app uses.** Sequence, CDRs, the structures, species,
  method, resolution, antigen. Everything else stays in SAbDab2, where it is
  searchable and current.
* **Rows with no VH are dropped.** The export's `SD-L` entries are
  single-domain *light* chains and carry a VL instead; there is nothing for a
  VHH tool to compare against. The count of what was dropped goes in the
  metadata rather than being silently absent.
* **Sequences carrying `X` are dropped too.** Nine of them do, where a residue
  in the deposited model could not be assigned. The app's own sequence cleaner
  refuses `X` outright, and an identity computed against a sequence with holes
  in it is a number nobody can interpret. Counted in the metadata like the
  rest.
* **`type` is kept, not filtered on.** `VNAR` is a shark domain, not a VHH,
  and a user whose closest hit is a VNAR needs to be told that rather than
  have the hit hidden. The app labels them; this script keeps them.

**The extraction date is read from the data, not from the clock.** The
greatest `SABDABupdate_date` in the export is when SAbDab2 last touched any
row in it, which is what a user needs to know. Taking today's date instead
would stamp a rebuild of an old file as fresh.

Run:  python tools/sabdab/build_index.py <export.csv>

Writes `data/sabdab/single_domain.csv` and `data/sabdab/meta.json`, prints
what it kept and dropped, and exits non-zero if anything about the export is
not what this script expects.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = ROOT / "data" / "sabdab"

# The columns this script reads. Checked against the export's header before
# anything else happens: SAbDab2 is a live service and its summary format can
# change, and a renamed column must stop the build rather than quietly produce
# an index with an empty field in it.
REQUIRED = ("INSTANCE", "PDB", "Hchain", "VH", "CDR-H1", "CDR-H2", "CDR-H3",
            "type", "method", "resolution", "heavy_species", "antigen_type",
            "antigen_name", "antigen_species", "date", "SABDABupdate_date")

# What one row of the shipped index holds.
FIELDS = ("sequence", "cdr1", "cdr2", "cdr3", "type", "species",
          "n_structures", "structures", "best_resolution", "methods",
          "antigen_names", "antigen_types", "antigen_species", "first_date")

LICENCE = ("SAbDab2 data, Oxford Protein Informatics Group, CC-BY 4.0. "
           "Cite Schneider C, Raybould MIJ, Deane CM, Nucleic Acids Res "
           "50(D1):D1368-D1372 (2022), doi 10.1093/nar/gkab1050.")

SOURCE = "https://sabdab.opig.stats.ox.ac.uk/"

# The twenty residues the app can compute with. `vhh_core.clean_sequence`
# raises on anything else, so a sequence carrying one cannot be compared and
# has no business in the index.
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def _pdb_code(raw: str) -> str:
    """`pdb_000010nd` to `10nd`, the code people actually quote.

    SAbDab2 uses the PDB's extended identifiers internally. The four-character
    code is what a user types into the PDB, and what fits in a table, so it is
    what the index carries; the long form is recoverable by prefixing it back.
    """
    code = raw.strip()
    if code.lower().startswith("pdb_") and len(code) > 4:
        code = code[4:].lstrip("0")
    return code.lower()


def _resolution(raw: str):
    """Resolution as a float, or None. Cryo-EM entries sometimes have none."""
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def build(export: Path) -> dict:
    with export.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"the export is missing columns this script needs: {missing}. "
                f"SAbDab2's summary format has changed; update REQUIRED and "
                f"the fields below rather than working around it.")
        rows = list(reader)

    by_sequence: dict[str, list[dict]] = defaultdict(list)
    dropped_no_vh = 0
    dropped_nonstandard = set()
    for row in rows:
        sequence = (row["VH"] or "").strip().upper()
        if not sequence:
            dropped_no_vh += 1
            continue
        if set(sequence) - STANDARD:
            dropped_nonstandard.add(sequence)
            continue
        by_sequence[sequence].append(row)

    out = []
    for sequence, group in by_sequence.items():
        # Best resolution first, then by PDB code, so the order is stable
        # across rebuilds: a set's iteration order is not, and an index whose
        # rows reshuffle on every rebuild makes every diff unreadable.
        ordered = sorted(
            group,
            key=lambda r: (_resolution(r["resolution"]) is None,
                           _resolution(r["resolution"]) or 0.0,
                           _pdb_code(r["PDB"])))
        seen, structures = set(), []
        for row in ordered:
            label = f"{_pdb_code(row['PDB'])}:{row['Hchain']}"
            if label not in seen:
                seen.add(label)
                structures.append(label)
        resolutions = [_resolution(r["resolution"]) for r in ordered]
        resolutions = [r for r in resolutions if r is not None]
        first = ordered[0]

        def unique(column):
            """Distinct non-empty values, in first-seen order."""
            found = []
            for row in ordered:
                value = (row[column] or "").strip()
                if value and value not in found:
                    found.append(value)
            return found

        out.append({
            "sequence": sequence,
            "cdr1": first["CDR-H1"],
            "cdr2": first["CDR-H2"],
            "cdr3": first["CDR-H3"],
            "type": first["type"],
            "species": "; ".join(unique("heavy_species")),
            "n_structures": len(structures),
            "structures": " ".join(structures),
            "best_resolution": min(resolutions) if resolutions else "",
            "methods": "; ".join(unique("method")),
            "antigen_names": " | ".join(unique("antigen_name")),
            "antigen_types": "; ".join(unique("antigen_type")),
            "antigen_species": "; ".join(unique("antigen_species")),
            "first_date": min(r["date"] for r in ordered if r["date"]) or "",
        })

    out.sort(key=lambda r: r["sequence"])

    extracted = max((r["SABDABupdate_date"] for r in rows
                     if r["SABDABupdate_date"]), default="")
    if len(extracted) == 8 and extracted.isdigit():
        extracted = f"{extracted[:4]}-{extracted[4:6]}-{extracted[6:]}"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "single_domain.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(out)

    meta = {
        "source": SOURCE,
        "licence": LICENCE,
        "extracted": extracted,
        "export_rows": len(rows),
        "sequences": len(out),
        "dropped_no_vh": dropped_no_vh,
        "dropped_nonstandard": len(dropped_nonstandard),
        "types": {t: sum(1 for r in out if r["type"] == t)
                  for t in sorted({r["type"] for r in out})},
        "built_by": "tools/sabdab/build_index.py",
    }
    (OUT_DIR / "meta.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    (OUT_DIR / "LICENSE").write_text(
        LICENCE + "\n\nThe full licence text is at "
        "https://creativecommons.org/licenses/by/4.0/\n"
        "This applies to the data in this directory only. The rest of this "
        "repository is MIT; see LICENSE at the root.\n", encoding="utf-8")
    return meta


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    meta = build(Path(sys.argv[1]))
    print(f"{meta['export_rows']} rows in, {meta['sequences']} distinct "
          f"sequences out, {meta['dropped_no_vh']} dropped for having no VH, "
          f"{meta['dropped_nonstandard']} for carrying a non-standard "
          f"residue.")
    print(f"types: {meta['types']}")
    print(f"extracted: {meta['extracted']}")
    print(f"written to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
