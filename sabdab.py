"""
sabdab - has anything like this already been solved?

One question, asked against a snapshot of SAbDab2's single-domain entries:
for a VHH you are holding, which experimentally determined structures come
closest, how close, and what were they bound to.

It is the cheapest thing you can know before deciding to model, crystallise,
or simply trust a prediction. A query that lands at 95% identity to a solved
complex has a template and an epitope to read. One that lands at 60% is a
scaffold nobody has seen, which is interesting for a different reason.

**This is a snapshot, and the app says so everywhere it shows a result.**
SAbDab2 is updated weekly from the PDB; the file in `data/sabdab/` is not.
`meta()` carries the extraction date, `age_days()` turns it into a number,
and `staleness()` turns that into a sentence, so that the warning is derived
from the data rather than from somebody remembering to update a caption. Past
`STALE_DAYS` the app says plainly that the index is old and points at the
live database.

**What the snapshot is.** `tools/sabdab/build_index.py` derives it from
SAbDab2's own single-domain summary export: one row per distinct VH sequence,
with every structure that contains it, ordered best-resolution first. Rows
with no VH (single-domain *light* chains) and sequences carrying `X` are
dropped, both counted in the metadata. Nothing is edited by hand, which is
why the export can be replaced and the index rebuilt in one command.

**The data is CC-BY 4.0**, not MIT like the rest of this repository. See
`data/sabdab/LICENSE`, and cite SAbDab2 if a result from this screen reaches
a paper.

**What this is not.** It is not a structure prediction, not a template
selection, and not a guarantee that a high-identity hit folds the same way:
two VHHs at 90% identity over the framework can have unrelated CDR3 loops,
which is why CDR3 identity is reported separately and should usually be read
first. And a VHH is not the only single domain in there: shark VNAR entries
are kept, labelled, and excluded by default, because a VNAR hit answers a
question you did not ask.
"""

from __future__ import annotations

import datetime as _dt
import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

import vhh_core as vc

__version__ = "2026.11.12a"

__all__ = ["INDEX_FILE", "META_FILE", "STALE_DAYS", "load_index", "meta",
           "age_days", "staleness", "closest", "summary", "__version__"]

DATA_DIR = Path(__file__).parent / "data" / "sabdab"
INDEX_FILE = DATA_DIR / "single_domain.csv"
META_FILE = DATA_DIR / "meta.json"

# When the app starts calling its own index old. A convention, not a finding:
# SAbDab2 adds entries weekly, so ninety days is roughly thirteen updates
# missed, which is enough that a user deciding whether to crystallise
# something should go and look at the live database rather than trust this
# file. Nothing fixes it at ninety rather than sixty or a hundred and twenty.
STALE_DAYS = 90

# How many hits the screen shows unless asked for more. Enough to see whether
# the top hit stands alone or sits in a cluster of near-identical entries,
# which is the thing a single best match cannot tell you.
DEFAULT_TOP = 5


@lru_cache(maxsize=1)
def load_index() -> pd.DataFrame:
    """The shipped snapshot, as a frame. Read once per process."""
    if not INDEX_FILE.exists():
        raise FileNotFoundError(
            f"{INDEX_FILE} is missing. Rebuild it with "
            f"tools/sabdab/build_index.py from a SAbDab2 single-domain "
            f"summary export.")
    frame = pd.read_csv(INDEX_FILE, dtype=str, keep_default_na=False)
    frame["n_structures"] = pd.to_numeric(frame["n_structures"],
                                          errors="coerce")
    frame["best_resolution"] = pd.to_numeric(frame["best_resolution"],
                                             errors="coerce")
    return frame


@lru_cache(maxsize=1)
def meta() -> dict:
    """What the snapshot is and when it was taken."""
    if not META_FILE.exists():
        return {}
    return json.loads(META_FILE.read_text(encoding="utf-8"))


def age_days(today: _dt.date | None = None) -> int | None:
    """Days between the extraction date and ``today``, or None if unknown.

    ``today`` is an argument so that the tests do not have to wait for time to
    pass, and so that nothing in the app reads the clock twice and gets two
    different answers within one render.
    """
    stamp = (meta() or {}).get("extracted")
    if not stamp:
        return None
    try:
        extracted = _dt.date.fromisoformat(stamp)
    except ValueError:
        return None
    return (( today or _dt.date.today()) - extracted).days


def staleness(today: _dt.date | None = None) -> dict:
    """The age of the snapshot, and whether the app should say it is old.

    Returns ``{"days", "extracted", "stale", "message"}``. The message is a
    full sentence because every caller shows it verbatim: a second place that
    turns the number into words is a second place that can disagree with this
    one about what the number means.
    """
    days = age_days(today)
    stamp = (meta() or {}).get("extracted")
    if days is None:
        return {"days": None, "extracted": stamp, "stale": True,
                "message": "This index carries no extraction date, so there "
                           "is no way to tell how old it is. Treat it as old."}
    stale = days >= STALE_DAYS
    when = f"extracted {stamp}, {days} day{'s' if days != 1 else ''} ago"
    if stale:
        message = (f"This index was {when}. SAbDab2 adds entries weekly, so "
                   f"it has missed roughly {days // 7} updates: a structure "
                   f"solved since then cannot appear here. Check the live "
                   f"database before concluding that nothing close exists.")
    else:
        message = (f"Snapshot {when}. SAbDab2 is updated weekly and this file "
                   f"is not, so anything deposited since is missing.")
    return {"days": days, "extracted": stamp, "stale": stale,
            "message": message}


def closest(seq: str, top: int = DEFAULT_TOP, include_vnar: bool = False,
            annotation: dict | None = None) -> pd.DataFrame:
    """The solved single domains closest to ``seq``, best first.

    Identity is full-length and computed the same way as everywhere else in
    this app, through `vhh_core.pairwise_identity`, so a number here and a
    number in the Alignment tab mean the same thing.

    ``cdr3_identity`` is the column to read first, and it is reported
    separately rather than folded into a single score on purpose: two VHHs can
    sit at 90% over the framework and have CDR3 loops with nothing in common,
    which is the case where a high overall identity is most misleading. It is
    computed only when the query's CDR3 is known, which needs ``annotation``;
    without one the column is empty rather than guessed at.

    VNAR entries are shark domains and are left out unless ``include_vnar``,
    because a VNAR hit answers a question about a different molecule.
    """
    query = vc.clean_sequence(seq)
    frame = load_index()
    if not include_vnar:
        frame = frame[frame["type"] != "VNAR"]

    query_cdr3 = ""
    if annotation:
        query_cdr3 = (annotation.get("CDR3") or {}).get("seq") or ""

    rows = []
    for entry in frame.itertuples():
        identity = vc.pairwise_identity(query, entry.sequence)
        rows.append({
            "identity": identity,
            "cdr3_identity": (vc.pairwise_identity(query_cdr3, entry.cdr3)
                              if query_cdr3 and entry.cdr3 else None),
            "type": entry.type,
            "structures": entry.structures,
            "n_structures": entry.n_structures,
            "best_resolution": entry.best_resolution,
            "methods": entry.methods,
            "species": entry.species,
            "antigen_names": entry.antigen_names,
            "antigen_types": entry.antigen_types,
            "cdr3": entry.cdr3,
            "sequence": entry.sequence,
        })

    hits = pd.DataFrame(rows)
    if not len(hits):
        return hits
    hits = hits.sort_values("identity", ascending=False).head(int(top))
    return hits.reset_index(drop=True)


def summary(hits: pd.DataFrame) -> dict:
    """One line's worth of what the hit list says, for the screen's header.

    ``spread`` is the gap between the best hit and the next one. A top hit
    three points clear of the field is a single close relative; a flat list is
    a family, and the difference changes what you would do next.
    """
    if hits is None or not len(hits):
        return {"best": None, "spread": None, "n": 0}
    best = float(hits.iloc[0]["identity"])
    second = float(hits.iloc[1]["identity"]) if len(hits) > 1 else None
    return {
        "best": best,
        "second": second,
        "spread": round(best - second, 1) if second is not None else None,
        "n": len(hits),
        "structures": hits.iloc[0]["structures"],
    }
