"""
Extract the polyreactivity model and both background distributions from the
authors' own release, and check that this repository reproduces their scores.

Why this exists: `specificity.py` must not ship a number nobody can trace. The
three files it reads are all built here, from

    github.com/debbiemarkslab/nanobody-polyreactivity   (MIT, Debora Marks Lab)

which accompanies

    Harvey EP, Shin JE, Skiba MA, et al. An in silico method to assess antibody
    fragment polyreactivity. Nat Commun 13, 7554 (2022).
    doi:10.1038/s41467-022-35276-4

Their repository is 115 MB, almost all of it neural-network checkpoints this
app does not use, so it is not vendored. Clone it next to this script:

    git clone --depth 1 https://github.com/debbiemarkslab/nanobody-polyreactivity
    python build.py nanobody-polyreactivity

What is written, all into `data/specificity/`:

* `harvey_lr3mer_cdrs.npz`, the non-zero weights of their 3-mer logistic
  regression on IMGT CDRs (`logistic_regression_3mer_CDRS.sav`), as index and
  value pairs in float64 plus the intercept. 791 of 8420 weights are non-zero,
  so this is both smaller than a dense array and exact. scikit-learn is needed
  here, to open their pickle, and nowhere in the app.
* `polyreactivity_quantiles.npz`, the quantile grid of the score distribution
  over their low-polyreactivity pool, and the size of that pool.
* `plastic_quantiles.npz`, the quantile grids of the two plastic resemblance
  modes over the same pool, which is the naive synthetic library the app's
  percentiles are read against.

What is checked, and printed, every run:

1. the decision function computed by `specificity.py` against the
   `origFACS lr 3mers` column the authors publish for all 69 702 sequences of
   the clean pool, which is the only real test that the featurisation order,
   the L1 normalisation and the intercept are right;
2. the error introduced by storing quantiles rather than the full sample;
3. the fallback CDR anchors of `specificity.py` against the authors' own IMGT
   CDR columns, over both pools, 141 474 sequences: how many it accepts, and
   whether the CDRs it accepts are identical to theirs.

The authors' CSVs carry one column per IMGT position, so the full-length
sequence of every clone can be rebuilt from them by concatenation; that is what
makes check 3 possible on the whole set rather than on a handful of examples.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUT = REPO / "data" / "specificity"
sys.path.insert(0, str(REPO))

import specificity as sp  # noqa: E402

MODEL = "backend/app/models/logistic_regression_3mer_CDRS.sav"
POOLS = {"low": "backend/app/experiments/low_polyreactivity_high_throughput.csv",
         "high": "backend/app/experiments/high_polyreactivity_high_throughput.csv"}

# 10 001 points, so the grid spans 0 to 100 in steps of 0.01 percentile.
N_QUANTILES = 10_001


def _sequences_from_positions(frame: pd.DataFrame) -> list[str]:
    """Rebuild each full-length sequence from the IMGT position columns.

    Everything from the column after `j_identity` up to the first CDR column is
    one IMGT position each, holding a residue or "-". Taken by position rather
    than by name because the two pools do not carry the same insertion codes,
    so the column names differ between the files while their order does not.
    """
    cols = list(frame.columns)
    start = cols.index("j_identity") + 1
    end = min(cols.index(c) for c in cols if c.startswith("CDR"))
    block = frame.iloc[:, start:end].astype(str)
    joined = block.agg("".join, axis=1)
    return [s.replace("-", "").replace("nan", "").upper() for s in joined]


def main(root: Path) -> int:
    if not (root / MODEL).exists():
        print(f"not found: {root / MODEL}\nclone the repository first, see the "
              f"docstring of this file")
        return 2
    OUT.mkdir(parents=True, exist_ok=True)

    # --- 1. the model ------------------------------------------------------ #
    with open(root / MODEL, "rb") as fh:
        model = pickle.load(fh)
    coef = np.asarray(model.coef_).ravel().astype(np.float64)
    intercept = float(np.asarray(model.intercept_).ravel()[0])
    if coef.size != sp.N_FEATURES:
        print(f"the released model has {coef.size} weights, this featurisation "
              f"produces {sp.N_FEATURES}")
        return 1
    index = np.flatnonzero(coef)
    np.savez_compressed(OUT / "harvey_lr3mer_cdrs.npz",
                        index=index.astype(np.int32),
                        value=coef[index],
                        intercept=np.float64(intercept),
                        n_features=np.int32(coef.size))
    print(f"model: {index.size} non-zero of {coef.size} weights, "
          f"intercept {intercept!r}")

    # --- 2. scoring parity, on the authors' published column --------------- #
    low = pd.read_csv(root / POOLS["low"], low_memory=False)
    published = low["origFACS lr 3mers"].to_numpy(float)
    score = sp.polyreactivity(low["CDRS_nogaps"].astype(str).tolist())
    # `score` is negated by convention, so compare against -published.
    gap = float(np.max(np.abs(score + published)))
    print(f"parity: {len(low)} sequences, max |difference| from their own "
          f"'origFACS lr 3mers' column = {gap:.3g}")
    if gap > 1e-9:
        print("parity check failed, refusing to write the background files")
        return 1

    # --- 3. the polyreactivity background ---------------------------------- #
    probs = np.linspace(0.0, 1.0, N_QUANTILES)
    grid = np.quantile(np.sort(score), probs).astype(np.float32)
    np.savez_compressed(OUT / "polyreactivity_quantiles.npz",
                        quantiles=grid, n_pool=np.int32(len(score)))
    # How much the grid costs, measured rather than asserted: the percentile of
    # every sequence in the pool, read off the grid, against the exact rank.
    sp._poly_q = None
    from_grid = sp.polyreactivity_percentile(score)
    exact = (np.searchsorted(np.sort(score), score, side="right")
             / len(score) * 100.0)
    print(f"background: pool {len(score)}, quantile grid {len(grid)} points, "
          f"max percentile error {np.max(np.abs(from_grid - exact)):.3f}, "
          f"mean {np.mean(np.abs(from_grid - exact)):.4f}")

    # --- 4. the plastic background, from the same pool --------------------- #
    feats = sp.plastic_features(low["CDRS_nogaps"].astype(str).tolist())
    modes = {}
    for mode, weights in sp.PLASTIC_WEIGHTS.items():
        raw = sum(w * feats[k].to_numpy(float) for k, w in weights.items())
        modes[mode] = np.quantile(np.sort(np.asarray(raw)), probs).astype(np.float32)
    np.savez_compressed(OUT / "plastic_quantiles.npz", **modes)
    print(f"plastic background: {', '.join(modes)} over the same {len(low)} "
          f"sequences, {N_QUANTILES} points each")

    # --- 5. the fallback anchors against their IMGT CDRs ------------------- #
    total = accepted = identical = 0
    for name, path in POOLS.items():
        pool = low if name == "low" else pd.read_csv(root / path, low_memory=False)
        seqs = _sequences_from_positions(pool)
        truth = pool[["CDR1_nogaps", "CDR2_nogaps", "CDR3_nogaps"]].astype(str)
        for seq, (c1, c2, c3) in zip(seqs, truth.itertuples(index=False)):
            got = sp.cdrs_regex(seq)
            total += 1
            if got.ok:
                accepted += 1
                identical += (got.cdr1 == c1 and got.cdr2 == c2 and got.cdr3 == c3)
    print(f"anchors: {total} sequences, accepted {accepted} "
          f"({accepted / total:.1%}), of which identical to their IMGT CDRs "
          f"{identical} ({identical / max(accepted, 1):.4%})")

    for file in sorted(OUT.glob("*.npz")):
        print(f"wrote {file.relative_to(REPO)}, {file.stat().st_size / 1024:.0f} kB")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(Path(sys.argv[1]).resolve()))
