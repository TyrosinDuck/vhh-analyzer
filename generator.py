"""
generator — synthetic VHH sequences built from real position frequencies.

Frameworks are sampled position by position from an IMGT position-frequency
matrix derived from real VHHs, so they look like natural frameworks. CDRs are
drawn at random, with lengths taken from the empirical distribution of the same
reference set.

The bundled matrix comes from **PLAbDab-nano** (OPIG, University of Oxford),
1958 unique VHH sequences numbered with ANARCI against the alpaca germlines.
You can supply your own matrix instead — a panning campaign of your own gives a
generator that mimics your library rather than the literature.

Two honest limitations, both inherent to sampling a position-frequency matrix:

* **Positions are sampled independently.** Real frameworks have correlated
  positions (a residue at one site constrains another), so an occasional
  combination produced here never occurs in nature. The hallmark tetrad is the
  usual casualty, which is why ``generate`` can enforce it.
* **A frequency is only as good as its support.** The ``n_residues`` column says
  how many real sequences back each position; a position seen in 30 sequences
  carries far less information than one seen in 1900.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "AA", "REGION_BOUNDS", "load_matrix", "load_meta", "sample_cdr_length",
    "generate", "generate_many", "matrix_report", "__version__",
]

__version__ = "2026.10.21b"

AA = list("ACDEFGHIKLMNPQRSTVWY")

# IMGT region boundaries, the same definition the rest of the app uses
REGION_BOUNDS = [("FR1", 1, 26), ("CDR1", 27, 38), ("FR2", 39, 55),
                 ("CDR2", 56, 65), ("FR3", 66, 104), ("CDR3", 105, 117),
                 ("FR4", 118, 128)]

DATA_DIR = Path(__file__).parent / "data"
MATRIX_FILE = DATA_DIR / "vhh_position_freq.csv"
META_FILE = DATA_DIR / "vhh_position_meta.json"

# The four FR2 hallmark positions and the residues that make a framework
# camelid rather than VH-like. Used only when the caller asks for enforcement.
HALLMARK_POSITIONS = {42: set("FYLH"), 49: set("EQDG"), 50: set("RCLP"),
                      52: set("FGLSW")}
VH_LIKE = {42: "V", 49: "G", 50: "L", 52: "W"}


def load_matrix(path: str | Path | None = None) -> pd.DataFrame:
    """The position-frequency matrix, indexed by IMGT position label.

    Columns: ``region``, ``n_residues``, ``occupancy``, then one per amino acid.
    ``occupancy`` is the fraction of reference sequences carrying any residue at
    that position, which is what reproduces the natural gap pattern.
    """
    path = Path(path) if path else MATRIX_FILE
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    missing = [a for a in AA if a not in df.columns]
    if missing:
        raise ValueError(f"matrix is missing amino-acid columns: {missing}")
    return df


def load_meta(path: str | Path | None = None) -> dict:
    """Provenance and the empirical CDR length distributions."""
    path = Path(path) if path else META_FILE
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return {}


def matrix_report(matrix: pd.DataFrame, meta: dict | None = None) -> dict:
    """What the matrix is made of — shown to the user before they trust it."""
    meta = meta or {}
    by_region = (matrix.groupby("region")["n_residues"]
                 .agg(["count", "min", "max"]).rename(
                     columns={"count": "positions", "min": "min_support",
                              "max": "max_support"}))
    return {
        "source": meta.get("source", "unknown"),
        "n_sequences": meta.get("n_sequences"),
        "n_positions": int(len(matrix)),
        "median_support": int(matrix["n_residues"].median()),
        "min_support": int(matrix["n_residues"].min()),
        "by_region": by_region,
    }


def _imgt_sort_key(label: str) -> tuple:
    """IMGT ordering, including the insertion codes in CDR3 (111A, 112B...)."""
    m = re.match(r"(\d+)([A-Z]*)", label)
    num, ins = int(m.group(1)), m.group(2)
    if ins and num in (33, 61, 112):
        return (num - 1, 1000 - ord(ins[0]))
    return (num, ord(ins[0]) if ins else 0)


def _region_of(label: str) -> str:
    num = int(re.match(r"\d+", label).group())
    for name, lo, hi in REGION_BOUNDS:
        if lo <= num <= hi:
            return name
    return "?"


def sample_cdr_length(region: str, meta: dict, rng: np.random.Generator) -> int:
    """Draw a CDR length from the reference set's empirical distribution."""
    counts = (meta.get("cdr_length_counts") or {}).get(region)
    if not counts:
        return {"CDR1": 8, "CDR2": 8, "CDR3": 17}.get(region, 8)
    lengths = np.array([int(k) for k in counts])
    weights = np.array([counts[k] for k in counts], dtype=float)
    return int(rng.choice(lengths, p=weights / weights.sum()))


def generate(matrix: pd.DataFrame, meta: dict | None = None,
             rng: np.random.Generator | None = None,
             framework_mode: str = "sample", cdr_alphabet: str = "no_cys",
             cdr_lengths: dict | None = None, enforce_hallmarks: bool = True,
             occupancy_threshold: float = 0.5,
             framework_cys: bool = False) -> dict:
    """Build one synthetic VHH.

    ``framework_mode`` is ``"sample"`` (draw each framework position from its
    frequency column) or ``"consensus"`` (take the most frequent residue, which
    gives the same canonical framework every time).

    ``cdr_alphabet`` is ``"uniform"`` (all 20 residues equally likely),
    ``"no_cys"`` (the same minus cysteine, so no stray disulfide is created) or
    ``"natural"`` (residue frequencies taken from the reference set's own CDR
    composition, which produces loops that look far more plausible).

    ``enforce_hallmarks`` resamples the four FR2 positions until the tetrad is
    camelid. Independent sampling otherwise produces a VH-like residue there now
    and then, which would make the output fail the app's own hallmark check.

    ``framework_cys`` allows cysteines outside the canonical IMGT 23/104 pair.
    They do occur in the reference set — about 2% of sequences carry one at IMGT
    50, where it pairs with a CDR cysteine — but in a generator with random CDRs
    that partner does not exist, so such a residue would be left as an orphan
    thiol. It is therefore off by default.
    """
    meta = meta or load_meta()
    rng = rng or np.random.default_rng()
    cdr_lengths = cdr_lengths or {}

    positions = sorted(matrix.index, key=_imgt_sort_key)
    framework_positions = [p for p in positions
                           if not _region_of(p).startswith("CDR")]

    # ---- CDR alphabet ----
    if cdr_alphabet == "natural":
        cdr_rows = matrix[matrix["region"].astype(str).str.startswith("CDR")]
        weights = cdr_rows[AA].mul(cdr_rows["n_residues"], axis=0).sum()
        probs = (weights / weights.sum()).to_numpy()
        alphabet = np.array(AA)
    else:
        alphabet = np.array([a for a in AA if not (cdr_alphabet == "no_cys"
                                                   and a == "C")])
        probs = np.full(len(alphabet), 1 / len(alphabet))

    # ---- frameworks ----
    residues: dict[str, str] = {}
    for label in framework_positions:
        row = matrix.loc[label]
        if float(row["occupancy"]) < occupancy_threshold:
            continue                      # position absent from most real VHHs
        column = row[AA].to_numpy(dtype=float)
        if column.sum() <= 0:
            continue
        canonical_cys = label in ("23", "104")
        if not framework_cys and not canonical_cys:
            column = column * np.array([0.0 if a == "C" else 1.0 for a in AA])
            if column.sum() <= 0:
                continue
        if framework_mode == "consensus":
            residues[label] = AA[int(np.argmax(column))]
        else:
            residues[label] = str(rng.choice(AA, p=column / column.sum()))

    if enforce_hallmarks and framework_mode == "sample":
        for imgt, allowed in HALLMARK_POSITIONS.items():
            label = str(imgt)
            if label not in matrix.index:
                continue
            column = matrix.loc[label, AA].to_numpy(dtype=float)
            mask = np.array([a in allowed and a != VH_LIKE[imgt]
                             and (framework_cys or a != "C") for a in AA],
                            dtype=float)
            weighted = column * mask
            if weighted.sum() > 0:
                residues[label] = str(rng.choice(AA, p=weighted / weighted.sum()))

    # ---- CDRs ----
    cdrs = {}
    for region in ("CDR1", "CDR2", "CDR3"):
        length = int(cdr_lengths.get(region)
                     or sample_cdr_length(region, meta, rng))
        cdrs[region] = "".join(rng.choice(alphabet, size=length, p=probs))

    # ---- assemble in IMGT order ----
    parts, by_region = [], {}
    for name, _lo, _hi in REGION_BOUNDS:
        if name.startswith("CDR"):
            segment = cdrs[name]
        else:
            segment = "".join(residues.get(p, "")
                              for p in framework_positions
                              if _region_of(p) == name)
        by_region[name] = segment
        parts.append(segment)

    sequence = "".join(parts)
    tetrad = "".join(residues.get(str(p), "?") for p in (42, 49, 50, 52))
    return {"sequence": sequence, "regions": by_region, "length": len(sequence),
            "tetrad": tetrad,
            "cdr_lengths": {k: len(v) for k, v in cdrs.items()}}


def generate_many(n: int, matrix: pd.DataFrame | None = None,
                  meta: dict | None = None, seed: int | None = None,
                  prefix: str = "synthVHH", **kwargs) -> pd.DataFrame:
    """Generate n synthetic VHHs as a table, deterministic for a given seed."""
    matrix = load_matrix() if matrix is None else matrix
    meta = load_meta() if meta is None else meta
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(int(n)):
        made = generate(matrix, meta, rng=rng, **kwargs)
        rows.append({
            "name": f"{prefix}_{i + 1:03d}", "sequence": made["sequence"],
            "length": made["length"], "tetrad": made["tetrad"],
            **{f"{k}_len": v for k, v in made["cdr_lengths"].items()},
            **{k: v for k, v in made["regions"].items()},
        })
    return pd.DataFrame(rows)


def to_fasta(df: pd.DataFrame, width: int = 60) -> str:
    parts = []
    for _, row in df.iterrows():
        seq = row["sequence"]
        parts.append(f">{row['name']}\n" +
                     "\n".join(seq[i:i + width] for i in range(0, len(seq), width)))
    return "\n".join(parts) + ("\n" if parts else "")
