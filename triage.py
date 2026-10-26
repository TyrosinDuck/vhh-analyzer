"""
triage — one row per sequence, every flag the app can raise, side by side.

The app computes a great deal about each sequence, spread over a dozen tabs.
That is the right shape for studying one candidate and the wrong shape for the
moment you actually have to decide: sixteen clones out of a screening, ten slots
in a synthesis order. This module collapses everything into one table you can
sort by whichever column matters today.

**There is no composite score, deliberately.** Weighting a deamidation site
against a hydrophobic patch against three divergent framework positions means
choosing what matters, and that depends on what the VHH is for — a reagent for
your own ELISA and a candidate you intend to humanise have nothing in common.
A single number would hide that choice behind a decimal point. The columns sit
next to each other; you sort.

Every metric here is computed by the same functions that produce the individual
tabs, so a number in this table and the same number in its own tab cannot
disagree.
"""

from __future__ import annotations

import pandas as pd

import repertoire as rp
import vhh_core as vc

__all__ = ["triage_table", "COLUMN_GROUPS", "FLAG_COLUMNS", "__version__"]

__version__ = "2026.10.26a"

# For the UI, so a column can be shown or hidden by theme rather than one by one
COLUMN_GROUPS = {
    "Identity": ["name", "type", "length", "cdr3_len"],
    "Physico-chemical": ["mw_kda", "pI", "gravy", "instability_index"],
    "Liabilities": ["liability_score", "n_high_risk", "n_liab_in_CDR"],
    "Aggregation": ["worst_window", "n_apr"],
    "Repertoire": ["n_unusual", "rarest"],
    "Cysteines": ["n_cys", "n_orphan_cys", "extra_disulfide"],
    "Antibody": ["hallmark_score", "humanness", "closest_human"],
}

# Boolean columns, so the table can be filtered down to "show me only what has
# something wrong with it"
FLAG_COLUMNS = ["has_high_risk", "has_orphan_cys", "has_unusual_residue",
                "has_apr", "liab_in_CDR"]


def _safe(fn, default=None):
    try:
        return fn()
    except Exception:
        return default


def triage_table(seqs: dict[str, str], annotations: dict,
                 kinds: dict | None = None, ph: float | None = None,
                 include_added: bool = True,
                 rare_threshold: float = rp.RARE_THRESHOLD,
                 min_support: int = rp.MIN_SUPPORT) -> pd.DataFrame:
    """Every per-sequence metric the app computes, one row each.

    ``kinds`` is the classification map; without it the ``type`` column reads
    "unknown". Anything that cannot be computed for a given sequence — a
    non-antibody has no hallmark score — comes back as ``None`` rather than
    raising, so one awkward sequence never empties the table.
    """
    rows = []
    for name, seq in seqs.items():
        ann = annotations.get(name)
        kind = (kinds or {}).get(name) or {}
        row = {"name": name,
               "type": kind.get("kind", "unknown"),
               "length": len(vc.clean_sequence(seq))}

        props = _safe(lambda: vc.characterize(name, seq, ph=ph), {}) or {}
        row.update({
            "mw_kda": props.get("mw_kda"),
            "pI": props.get("pI"),
            "gravy": props.get("gravy"),
            "instability_index": props.get("instability_index"),
            "n_cys": props.get("n_Cys"),
        })

        is_domain = bool(kind.get("full_analysis"))
        row["cdr3_len"] = len(((ann or {}).get("CDR3") or {}).get("seq") or "") or None

        liab = _safe(lambda: vc.find_liabilities(seq, ann,
                                                 include_added=include_added))
        if liab is not None and len(liab):
            row.update({
                "liability_score": round(float(liab["score"].sum()), 1),
                "n_high_risk": int((liab["risk"] == "High").sum()),
                "n_liab_in_CDR": int(liab["in_CDR"].sum()),
            })
        else:
            row.update({"liability_score": 0.0, "n_high_risk": 0,
                        "n_liab_in_CDR": 0})

        def _agg():
            profile = vc.aggregation_profile(seq)
            regions = vc.aggregation_regions(profile, ann)
            summary = vc.aggregation_summary(profile, regions)
            return summary["min_score"], summary["n_regions"]
        worst, n_apr = _safe(_agg, (None, None))
        row.update({"worst_window": worst, "n_apr": n_apr})

        def _rare():
            table = rp.rare_residues(seq, threshold=rare_threshold,
                                     min_support=min_support)
            flagged = table[table["rare"]]
            if not len(flagged):
                return 0, None
            worst_row = flagged.sort_values("frequency").iloc[0]
            return len(flagged), (f"{worst_row['residue']}{worst_row['imgt']} "
                                  f"({worst_row['frequency']:.1%})")
        n_unusual, rarest = _safe(_rare, (None, None)) if is_domain else (None, None)
        row.update({"n_unusual": n_unusual, "rarest": rarest})

        # The disulfide, hallmark and humanness analyses all assume a variable
        # domain. Run on anything else they produce confident nonsense — the
        # canonical-pair logic reports lysozyme's eight perfectly ordinary
        # cysteines as eight orphan thiols — so they are left empty rather than
        # answered wrongly.
        ds = _safe(lambda: vc.disulfides(seq, ann)) if is_domain else None
        row.update({
            "n_orphan_cys": ds["n_orphan"] if ds else None,
            "extra_disulfide": ds["has_extra_pair"] if ds else None,
        })

        row["hallmark_score"] = _safe(
            lambda: vc.hallmarks(seq, ann).get("camelid_score")) if is_domain else None

        hum = _safe(lambda: vc.humanness(seq)) if is_domain else None
        row.update({
            "humanness": (hum or {}).get("identity"),
            "closest_human": (hum or {}).get("gene"),
        })

        row.update({
            "has_high_risk": bool(row["n_high_risk"]),
            "has_orphan_cys": bool(row["n_orphan_cys"]),
            "has_unusual_residue": bool(row["n_unusual"]),
            "has_apr": bool(row["n_apr"]),
            "liab_in_CDR": bool(row["n_liab_in_CDR"]),
        })
        rows.append(row)

    ordered = [c for group in COLUMN_GROUPS.values() for c in group] + FLAG_COLUMNS
    frame = pd.DataFrame(rows)
    return frame[[c for c in ordered if c in frame.columns]]
