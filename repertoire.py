"""
repertoire — compare sequences against the natural VHH repertoire.

The app already ships an IMGT position-frequency matrix built from 1958 real
VHHs (PLAbDab-nano), used so far only to *generate* synthetic sequences. Turned
around, the same table answers two questions about sequences you actually have:

* **Is any framework residue unusual?** A position carrying a residue seen in
  under a percent or two of the natural repertoire is worth a second look — it
  catches Sanger miscalls, PCR errors and constructs that have drifted without
  anyone noticing.
* **Is my selection output biased?** CDR3 lengths and per-region composition
  compared against the same reference say whether a panning enriched for long
  loops, charged residues, and so on.

Two cautions are built into the output rather than left to the reader.

*Support matters.* A frequency computed from 30 sequences is not a frequency
computed from 1900. Every rarity call carries the ``n_residues`` behind it, and
positions below a minimum support are reported separately instead of being
called rare.

*Rare is not wrong.* A residue absent from the reference set may be a genuine
novelty, a germline this set under-samples, or the very mutation you selected
for. This module says "unusual", never "error".
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import generator as gn
import vhh_core as vc

__all__ = [
    "rare_residues", "rarity_table", "plot_rarity", "cdr_length_comparison",
    "plot_cdr_lengths", "composition_comparison", "plot_composition",
    "RARE_THRESHOLD", "MIN_SUPPORT", "__version__",
]

__version__ = "2026.10.23a"

RARE_THRESHOLD = 0.02      # frequency at or below this counts as unusual
MIN_SUPPORT = 100          # positions backed by fewer reference sequences are set aside


def _matrix():
    return gn.load_matrix()


def _meta():
    return gn.load_meta()


# ---------------------------------------------------------------------------
# Unusual framework residues
# ---------------------------------------------------------------------------

def rare_residues(seq: str, matrix: pd.DataFrame | None = None,
                  threshold: float = RARE_THRESHOLD,
                  min_support: int = MIN_SUPPORT,
                  frameworks_only: bool = True) -> pd.DataFrame:
    """Framework positions whose residue is rare in the reference repertoire.

    Returns one row per compared position — every position, not only the rare
    ones, so the caller can show the distribution as well as the outliers.
    Columns: the IMGT number, the region, the residue, its frequency in the
    reference, the most frequent residue there and its frequency, the support,
    and the ``rare`` and ``low_support`` flags.

    CDRs are excluded by default: a CDR is supposed to be unusual, and flagging
    it would bury the framework signal that actually matters.
    """
    matrix = _matrix() if matrix is None else matrix
    numbering = vc.imgt_numbering(seq)

    rows = []
    for num, ins, aa, idx in numbering:
        label = f"{num}{ins}" if ins else str(num)
        if label not in matrix.index:
            continue
        region = vc.imgt_region(num)
        if frameworks_only and region.startswith("CDR"):
            continue
        row = matrix.loc[label]
        freq = float(row.get(aa, 0.0))
        support = int(row["n_residues"])
        counts = row[gn.AA].astype(float)
        top = counts.idxmax()
        rows.append({
            "imgt": label, "position": idx, "region": region, "residue": aa,
            "frequency": round(freq, 4),
            "consensus": top, "consensus_frequency": round(float(counts[top]), 4),
            "support": support,
            "low_support": support < min_support,
            "rare": bool(freq <= threshold and support >= min_support),
        })
    return pd.DataFrame(rows)


def rarity_table(seqs: dict[str, str], matrix: pd.DataFrame | None = None,
                 threshold: float = RARE_THRESHOLD,
                 min_support: int = MIN_SUPPORT) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(per_sequence_summary, all_rare_positions)`` over a set of sequences."""
    matrix = _matrix() if matrix is None else matrix
    summary, flagged = [], []
    for name, seq in seqs.items():
        try:
            table = rare_residues(seq, matrix, threshold, min_support)
        except Exception as exc:
            summary.append({"name": name, "compared": 0, "n_rare": 0,
                            "n_absent": 0, "rarest": None, "note": str(exc)})
            continue
        rare = table[table["rare"]]
        absent = rare[rare["frequency"] == 0]
        rarest = (rare.sort_values("frequency").iloc[0]
                  if len(rare) else None)
        summary.append({
            "name": name,
            "compared": int(len(table)),
            "n_rare": int(len(rare)),
            "n_absent": int(len(absent)),
            "rarest": (f"{rarest['residue']}{rarest['imgt']} "
                       f"({rarest['frequency']:.1%})") if rarest is not None else "—",
            "low_support_positions": int(table["low_support"].sum()),
            "note": "",
        })
        if len(rare):
            flagged.append(rare.assign(sequence=name))

    return (pd.DataFrame(summary),
            pd.concat(flagged, ignore_index=True) if flagged else pd.DataFrame())


def plot_rarity(seq: str, name: str = "", matrix: pd.DataFrame | None = None,
                threshold: float = RARE_THRESHOLD,
                min_support: int = MIN_SUPPORT):
    """Frequency of each framework residue along the sequence, outliers marked."""
    import matplotlib.pyplot as plt

    table = rare_residues(seq, matrix, threshold, min_support)
    fig, ax = plt.subplots(figsize=(11, 3.2))
    if not len(table):
        ax.text(0.5, 0.5, "No position could be compared", ha="center",
                va="center")
        ax.axis("off")
        return fig

    x = table["position"].to_numpy()
    y = table["frequency"].to_numpy()
    ax.vlines(x, 0, y, color="0.75", lw=1.2)
    ok = ~table["rare"] & ~table["low_support"]
    ax.plot(x[ok.to_numpy()], y[ok.to_numpy()], "o", ms=3.5, color="#4c72b0")
    low = table["low_support"].to_numpy()
    if low.any():
        ax.plot(x[low], y[low], "o", ms=3.5, color="0.6")
    rare = table["rare"].to_numpy()
    if rare.any():
        ax.plot(x[rare], y[rare], "o", ms=6, color="crimson")
        for _, r in table[table["rare"]].iterrows():
            ax.annotate(f"{r['residue']}{r['imgt']}",
                        xy=(r["position"], r["frequency"]),
                        xytext=(0, 9), textcoords="offset points",
                        ha="center", fontsize=7.5, color="crimson")
    ax.axhline(threshold, color="crimson", ls="--", lw=0.9)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Frequency in the\nnatural repertoire", fontsize=9)
    ax.set_title(f"{name} — framework residues against the reference repertoire"
                 if name else "Framework residues against the reference repertoire",
                 fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# The submitted set against the natural repertoire
# ---------------------------------------------------------------------------

def cdr_length_comparison(annotations: dict,
                          meta: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """CDR lengths of the submitted set against the reference distribution.

    Returns ``(per_sequence, distribution)``. The distribution carries, per CDR
    and length, the reference frequency and the submitted frequency, so the two
    can be drawn on the same axes.
    """
    meta = _meta() if meta is None else meta
    counts = meta.get("cdr_length_counts") or {}

    per_seq = []
    for name, ann in annotations.items():
        row = {"name": name}
        for region in ("CDR1", "CDR2", "CDR3"):
            loop = ((ann or {}).get(region) or {}).get("seq") or ""
            row[region] = len(loop)
        per_seq.append(row)
    per_seq = pd.DataFrame(per_seq)

    rows = []
    for region in ("CDR1", "CDR2", "CDR3"):
        ref = counts.get(region) or {}
        ref_total = sum(ref.values()) or 1
        mine = (per_seq[region].value_counts().to_dict()
                if region in per_seq.columns and len(per_seq) else {})
        mine_total = sum(mine.values()) or 1
        lengths = sorted({int(k) for k in ref} | {int(k) for k in mine})
        for length in lengths:
            rows.append({
                "region": region, "length": length,
                "reference_frequency": round(ref.get(str(length),
                                                     ref.get(length, 0)) / ref_total, 4),
                "submitted_frequency": round(mine.get(length, 0) / mine_total, 4),
                "submitted_count": int(mine.get(length, 0)),
            })
    return per_seq, pd.DataFrame(rows)


def plot_cdr_lengths(distribution: pd.DataFrame, n_submitted: int = 0):
    """Reference CDR length distribution with the submitted set overlaid."""
    import matplotlib.pyplot as plt

    regions = ["CDR1", "CDR2", "CDR3"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.1))
    for ax, region in zip(axes, regions):
        block = distribution[distribution["region"] == region]
        if not len(block):
            ax.axis("off")
            continue
        used = block[(block["reference_frequency"] > 0.002)
                     | (block["submitted_count"] > 0)]
        x = used["length"].to_numpy()
        ax.bar(x, used["reference_frequency"], color="0.78", width=0.8,
               label="natural repertoire")
        ax.bar(x, used["submitted_frequency"], color="#c44e52", width=0.42,
               label=f"submitted (n={n_submitted})")
        ax.set_title(region, fontsize=10)
        ax.set_xlabel("Length (aa)")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("Fraction")
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    return fig


def composition_comparison(seqs: dict[str, str], annotations: dict,
                           matrix: pd.DataFrame | None = None) -> pd.DataFrame:
    """Amino-acid composition per region, submitted set against the reference.

    The reference composition is the position-frequency matrix averaged over the
    positions of each region, weighted by their support. ``delta`` is the
    submitted fraction minus the reference fraction — positive means enriched.
    """
    matrix = _matrix() if matrix is None else matrix

    rows = []
    for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
        block = matrix[matrix["region"].astype(str) == region]
        if not len(block):
            continue
        weights = block["n_residues"].astype(float)
        ref = (block[gn.AA].mul(weights, axis=0).sum() / weights.sum())
        ref = ref / ref.sum() if ref.sum() else ref

        mine = pd.Series(0.0, index=gn.AA)
        total = 0
        for name in seqs:
            loop = ((annotations.get(name) or {}).get(region) or {}).get("seq") or ""
            for aa in loop:
                if aa in mine.index:
                    mine[aa] += 1
                    total += 1
        mine = mine / total if total else mine

        for aa in gn.AA:
            rows.append({
                "region": region, "residue": aa,
                "reference": round(float(ref.get(aa, 0.0)), 4),
                "submitted": round(float(mine.get(aa, 0.0)), 4),
                "delta": round(float(mine.get(aa, 0.0) - ref.get(aa, 0.0)), 4),
                "n_residues": total,
            })
    return pd.DataFrame(rows)


def plot_composition(comparison: pd.DataFrame, regions: list[str] | None = None):
    """Enrichment per amino acid, submitted minus reference, by region."""
    import matplotlib.pyplot as plt

    regions = regions or ["CDR1", "CDR2", "CDR3"]
    blocks = [r for r in regions
              if len(comparison[comparison["region"] == r])]
    if not blocks:
        fig, ax = plt.subplots(figsize=(6, 1.6))
        ax.text(0.5, 0.5, "Nothing to compare", ha="center", va="center")
        ax.axis("off")
        return fig

    fig, axes = plt.subplots(len(blocks), 1, figsize=(11, 2.3 * len(blocks)),
                             squeeze=False, sharex=True)
    for ax, region in zip(axes[:, 0], blocks):
        block = comparison[comparison["region"] == region].set_index("residue")
        values = block.loc[gn.AA, "delta"].to_numpy()
        colours = ["#c44e52" if v > 0 else "#4c72b0" for v in values]
        ax.bar(gn.AA, values, color=colours)
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_ylabel("Δ fraction", fontsize=8)
        ax.set_title(f"{region} — enriched above the line, depleted below "
                     f"(n = {int(block['n_residues'].iloc[0])} residues)",
                     fontsize=9, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig
