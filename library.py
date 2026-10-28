"""
library, the sequence library shipped with the app, read-only.

A curated list of sequences kept as a file in the repository, so that the
sequences you work with repeatedly are one click away instead of being pasted
again every time.

**It is deliberately read-only from the app.** Streamlit Community Cloud rebuilds
the container from the repository on every reboot, redeploy and wake-from-sleep,
so anything the app wrote to disk would be erased, usually on the day you had
forgotten it was fragile. Keeping the library in the repository instead means it
is durable, versioned, and restorable from git history. The cost is that editing
happens on your own machine, not in the browser.

The workflow: edit ``data/library.csv`` locally, commit and push, and the app
serves the new version on its next deploy.

**Everything in this file is public.** It lives in a public repository, so it is
readable by anyone, whether or not they ever open the app. Nothing confidential, unpublished sequences, anything under an NDA, anything from a collaborator who
has not agreed, belongs here.

Format: one row per molecule, columns ``name``, ``protein``, ``dna``, ``tags``,
``note``. At least one of ``protein`` or ``dna`` must be present; an entry may
carry both, in which case the DNA is checked against the protein when the file
is loaded.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

import vhh_core as vc

__all__ = [
    "LIBRARY_FILE", "COLUMNS", "load_library", "library_fasta",
    "blank_library", "__version__",
]

__version__ = "2026.10.28a"

LIBRARY_FILE = Path(__file__).parent / "data" / "library.csv"
COLUMNS = ["name", "protein", "dna", "tags", "note"]


def blank_library() -> pd.DataFrame:
    """An empty library with the right columns, for when the file is missing."""
    return pd.DataFrame(columns=COLUMNS + ["has_protein", "has_dna",
                                           "protein_length", "dna_length",
                                           "issue"])


def _clean_dna(dna: str) -> str:
    dna = re.sub(r"[^A-Za-z]", "", dna or "").upper().replace("U", "T")
    bad = set(dna) - set("ACGTN")
    if bad:
        raise ValueError(f"unexpected characters: {', '.join(sorted(bad))}")
    return dna


def load_library(path: str | Path | None = None) -> tuple[pd.DataFrame, list[str]]:
    """Read the library file. Returns ``(table, problems)``.

    A malformed row is kept and flagged rather than dropped: a library that
    silently loses an entry is worse than one that shows you which line to fix.
    """
    path = Path(path) if path else LIBRARY_FILE
    problems: list[str] = []
    if not path.exists():
        return blank_library(), [f"{path.name} not found, the library is empty."]

    try:
        frame = pd.read_csv(path, comment="#", dtype=str).fillna("")
    except Exception as exc:
        return blank_library(), [f"{path.name} could not be read: {exc}"]

    missing = [c for c in ("name", "protein", "dna") if c not in frame.columns]
    if missing:
        return blank_library(), [
            f"{path.name} is missing the column(s) {', '.join(missing)}. "
            f"Expected: {', '.join(COLUMNS)}."]
    for column in COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    rows = []
    seen: set[str] = set()
    for index, row in frame.iterrows():
        line = index + 2          # +1 for the header, +1 for 1-based counting
        name = (row["name"] or "").strip()
        issue = ""

        if not name:
            name = f"row_{line}"
            issue = "no name, one was generated"
        if name in seen:
            issue = "; ".join(filter(None, [issue, "duplicate name"]))
        seen.add(name)

        protein = (row["protein"] or "").strip()
        dna = (row["dna"] or "").strip()
        if protein:
            try:
                protein = vc.clean_sequence(protein)
            except ValueError as exc:
                issue = "; ".join(filter(None, [issue, f"protein: {exc}"]))
                protein = ""
        if dna:
            try:
                dna = _clean_dna(dna)
            except ValueError as exc:
                issue = "; ".join(filter(None, [issue, f"DNA: {exc}"]))
                dna = ""

        if not protein and not dna:
            issue = "; ".join(filter(None, [issue, "no usable sequence"]))

        # When both are present they should agree; a mismatch is worth knowing
        # about, and this is the only place the file is ever checked.
        if protein and dna:
            try:
                translated = vc.translate_dna(dna).rstrip("*")
                if translated.lstrip("M") != protein.lstrip("M") \
                        and translated != protein:
                    issue = "; ".join(filter(
                        None, [issue, "the DNA does not translate to the protein"]))
            except Exception:
                pass

        if issue:
            problems.append(f"line {line} ({name}): {issue}")

        rows.append({
            "name": name, "protein": protein, "dna": dna,
            "tags": (row.get("tags") or "").strip(),
            "note": (row.get("note") or "").strip(),
            "has_protein": bool(protein), "has_dna": bool(dna),
            "protein_length": len(protein) or None,
            "dna_length": len(dna) or None,
            "issue": issue,
        })

    return pd.DataFrame(rows, columns=blank_library().columns), problems


def library_fasta(table: pd.DataFrame, names: list[str], kind: str) -> str:
    """FASTA of the chosen entries, taking either the protein or the DNA.

    Entries that do not carry the requested kind are skipped, selecting a
    protein-only entry and asking for DNA is a no-op, not an error.
    """
    column = "protein" if kind == "protein" else "dna"
    parts = []
    for _, row in table[table["name"].isin(names)].iterrows():
        sequence = row[column]
        if not sequence:
            continue
        parts.append(f">{row['name']}\n" + "\n".join(
            sequence[i:i + 60] for i in range(0, len(sequence), 60)))
    return "\n".join(parts) + ("\n" if parts else "")


def all_tags(table: pd.DataFrame) -> list[str]:
    """Every distinct tag in the library, for filtering."""
    tags: set[str] = set()
    for value in table.get("tags", []):
        tags.update(t.strip() for t in str(value).split(",") if t.strip())
    return sorted(tags)
