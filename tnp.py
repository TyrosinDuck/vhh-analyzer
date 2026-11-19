"""
tnp - the Therapeutic Nanobody Profiler guidelines, applied to your sequences.

What this module does, in one sentence: it says where a VHH falls relative to
the 36 nanobodies that have reached clinical trials, on the metrics the
Therapeutic Nanobody Profiler (TNP) defines.

That is a different question from every other developability number in this
app. The liability score, the aggregation profile and the stability flags all
ask whether a sequence carries a known chemical or biophysical problem. TNP
asks something blunter and, for a candidate you intend to develop, more
useful: has anything like this ever been taken into a clinic. A value outside
the clinical range is not a defect. It is a sequence doing something no
approved or trialled nanobody does, which is worth knowing before you commit
to it.

**Where the numbers come from.** Gordon, Gervasio, Souders and Deane,
"Characterising nanobody developability to improve therapeutic design using
the Therapeutic Nanobody Profiler", Communications Biology 9:344 (2026),
https://doi.org/10.1038/s42003-026-09594-y. The flagging thresholds in
`METRICS` below are Table 2 of that paper, transcribed from the published PDF.
The authors state the thresholds will be revised as more nanobodies progress
through trials, so they are dated here and the app says that date.

**What is computed here and what is not.** TNP has six metrics. Two are pure
sequence and are computed in this module:

    total CDR length      CDR1 + CDR2 + CDR3, IMGT definitions
    CDR3 length           IMGT 105 to 117

One is structural and is computed in `structure.py`, from an uploaded model:

    CDR3 compactness      loop length divided by its reach

Three are structural and are **not** computed anywhere in this app:

    PSH                   patches of surface hydrophobicity
    PPC                   patches of positive charge
    PNC                   patches of negative charge

PSH, PPC and PNC follow the patch methodology of the Therapeutic Antibody
Profiler, which fixes a hydrophobicity scale, a patch definition and a
neighbourhood around the CDRs. Reimplementing that from the description
without the reference code would produce numbers on an unknown scale, and
comparing those to the published thresholds would be worse than reporting
nothing. They are listed with their thresholds, marked as not computed, and
the app points at the authors' own web application for them. See MANUAL.md
section 3.3b.

**IMGT or nothing.** TNP counts CDR residues in IMGT numbering. The app can
annotate by regex anchors when ANARCI is unavailable, and those boundaries are
not the same, so a profile built on a regex annotation carries
``scheme_ok=False`` and the app refuses to flag it. A flag against the wrong
boundary definition is a wrong answer with a colour on it.
"""

from __future__ import annotations

import pandas as pd

import vhh_core as vc

__version__ = "2026.11.19a"

__all__ = [
    "METRICS", "SEQUENCE_METRICS", "STRUCTURAL_METRICS", "FLAG_ORDER",
    "THRESHOLD_DATE", "REFERENCE", "WEBAPP",
    "cdr_lengths", "flag", "profile", "profile_table", "metric_table",
    "__version__",
]


REFERENCE = ("Gordon GL, Gervasio J, Souders C, Deane CM. Characterising "
             "nanobody developability to improve therapeutic design using the "
             "Therapeutic Nanobody Profiler. Commun Biol 9:344 (2026). "
             "https://doi.org/10.1038/s42003-026-09594-y")

WEBAPP = "https://opig.stats.ox.ac.uk/webapps/tnp"

# Table 2 of the paper carries no date of its own. The authors' web application
# labelled the same values "July 2025" when this module was written, and that
# label is the date of the most recent entry in the clinical-stage set behind
# the thresholds, not the day anybody recalculated them. The paper says plainly
# that they will be updated as more candidates progress through trials. Shown
# in the app so that a profile is never read as timeless.
#
# So what this string dates is the newest candidate the thresholds can know
# about. An earlier version of MANUAL.md turned it into "the 36 clinical-stage
# sequences available in July 2025", which claims a calibration date this
# module has never had evidence for.
THRESHOLD_DATE = "July 2025"

# Table 2, Gordon et al. 2026. `clinical` is the full range observed across the
# 36 clinical-stage nanobodies; `amber` are the bands inside it that the paper
# flags (the top and/or bottom 5% of the distribution); anything outside
# `clinical` is red. A metric with one amber band has one amber band: for the
# two charge patch metrics the clinical range starts at zero and only the top
# tail is flagged.
METRICS = {
    "total_cdr_length": {
        "label": "Total CDR length",
        "unit": "residues",
        "clinical": (20, 39),
        "amber": ((20, 24), (37, 39)),
        "computed": "sequence",
        "detail": "CDR1 + CDR2 + CDR3, IMGT boundaries.",
    },
    "cdr3_length": {
        "label": "CDR3 length",
        "unit": "residues",
        "clinical": (5, 23),
        "amber": ((5, 8), (21, 23)),
        "computed": "sequence",
        "detail": "IMGT 105 to 117.",
    },
    "cdr3_compactness": {
        "label": "CDR3 compactness",
        "unit": "residues per angstrom",
        "clinical": (0.56, 1.61),
        "amber": ((0.56, 0.81), (1.57, 1.61)),
        "computed": "structure",
        "detail": "Loop length divided by its reach from the rest of the "
                  "domain. Needs a model; see the Structure screen.",
    },
    "psh": {
        "label": "Patches of CDR surface hydrophobicity",
        "unit": "TAP patch score",
        "clinical": (73.40, 155.47),
        "amber": ((73.40, 79.59), (126.83, 155.47)),
        "computed": "not computed",
        "detail": "TAP patch methodology on a structure. Not implemented "
                  "here; use the TNP web application.",
    },
    "ppc": {
        "label": "Patches of CDR positive charge",
        "unit": "TAP patch score",
        "clinical": (0.00, 1.18),
        "amber": ((0.39, 1.18),),
        "computed": "not computed",
        # A patch score is a count of surface patches and cannot be negative,
        # and the clinical range starts at that floor, so there is no reachable
        # red band below it. `metric_table` reads this rather than guessing
        # from the shape of `amber`.
        "floor": 0.00,
        "detail": "TAP patch methodology on a structure. Not implemented "
                  "here; use the TNP web application.",
    },
    "pnc": {
        "label": "Patches of CDR negative charge",
        "unit": "TAP patch score",
        "clinical": (0.00, 1.88),
        "amber": ((1.47, 1.88),),
        "computed": "not computed",
        "floor": 0.00,
        "detail": "TAP patch methodology on a structure. Not implemented "
                  "here; use the TNP web application.",
    },
}

SEQUENCE_METRICS = tuple(k for k, m in METRICS.items()
                         if m["computed"] == "sequence")
STRUCTURAL_METRICS = tuple(k for k, m in METRICS.items()
                           if m["computed"] != "sequence")

# Worst first, which is the order a panel should be sorted in.
FLAG_ORDER = ("red", "amber", "green")

_FLAG_SYMBOL = {"red": "🔴", "amber": "🟠", "green": "🟢"}


def flag(metric: str, value) -> str | None:
    """Green, amber or red for one metric value, or None if there is no value.

    Red is outside the range spanned by the 36 clinical-stage nanobodies.
    Amber is inside it but in one of the tail bands the paper marks. Green is
    everything else. Band edges are inclusive on both sides, which is how
    Table 2 is written; a value sitting exactly on an edge is therefore amber
    rather than green, the cautious reading of an ambiguous boundary.

    One exception, and it is about a band that cannot be reached rather than
    about a judgement: where a metric declares a `floor` and its clinical range
    starts at that floor, nothing below the range exists to be red. PPC and PNC
    are patch scores bounded below by zero with a clinical range starting at
    0.00, so this used to promise a red verdict for a negative patch score.
    """
    if metric not in METRICS:
        raise KeyError(f"unknown TNP metric: {metric!r}")
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if value != value:                                   # NaN
        return None

    lo, hi = METRICS[metric]["clinical"]
    floor = METRICS[metric].get("floor")
    if value > hi or (value < lo and (floor is None or lo > floor)):
        return "red"
    for band_lo, band_hi in METRICS[metric]["amber"]:
        if band_lo <= value <= band_hi:
            return "amber"
    return "green"


def symbol(flag_name: str | None) -> str:
    """The coloured dot for a flag, empty for no flag."""
    return _FLAG_SYMBOL.get(flag_name or "", "")


def cdr_lengths(seq: str, annotation: dict | None = None,
                method: str | None = None) -> dict:
    """CDR1, CDR2, CDR3 and total lengths, with the annotation method used.

    Pass ``annotation`` when the caller already has one, which is the usual
    case in the app and avoids numbering the sequence twice. ``method`` then
    says how that annotation was made, and is what decides ``scheme_ok``:
    TNP's thresholds are IMGT counts, so a regex annotation gives lengths that
    are real but not comparable to them.
    """
    if annotation is None:
        annotation, method = vc.annotate(seq)
    method = method or "unknown"

    lengths = {"method": method, "scheme_ok": "IMGT" in method.upper(),
               "multi_domain": bool((annotation or {}).get("domains"))}

    # A chain carrying several variable domains has one set of CDRs per domain
    # and no single set of its own, so there is no length to return and no
    # total to add up. The lengths come back as None rather than as zero, or as
    # the first domain's: both of those are numbers, and a caller cannot tell a
    # number that means "nothing here" from one that was measured. The domains
    # are measured individually by passing each of them in turn.
    if lengths["multi_domain"]:
        for key in ("cdr1_length", "cdr2_length", "cdr3_length",
                    "total_cdr_length"):
            lengths[key] = None
        return lengths

    for region in ("CDR1", "CDR2", "CDR3"):
        block = annotation.get(region) or {}
        lengths[region.lower() + "_length"] = len(block.get("seq") or "")

    lengths["total_cdr_length"] = sum(
        lengths[f"cdr{i}_length"] for i in (1, 2, 3))
    return lengths


def profile(name: str, seq: str, annotation: dict | None = None,
            method: str | None = None,
            compactness: float | None = None) -> dict:
    """One sequence's TNP profile: the metrics, their flags, and the caveats.

    ``compactness`` is accepted rather than computed: it needs a structure,
    and this module never sees one. The Structure screen computes it and hands
    it back, so that a user who has uploaded a model gets three of the six
    metrics instead of two.
    """
    lengths = cdr_lengths(seq, annotation, method)
    scheme_ok = lengths["scheme_ok"]

    values = {
        "total_cdr_length": lengths["total_cdr_length"],
        "cdr3_length": lengths["cdr3_length"],
        "cdr3_compactness": compactness,
    }

    row = {"name": name, "method": lengths["method"], "scheme_ok": scheme_ok,
           "multi_domain": lengths.get("multi_domain", False),
           "cdr1_length": lengths["cdr1_length"],
           "cdr2_length": lengths["cdr2_length"]}
    flags = {}
    for metric in METRICS:
        value = values.get(metric)
        row[metric] = value
        # A length counted on regex boundaries is a real number against the
        # wrong definition, so it is reported and not flagged.
        flags[metric] = flag(metric, value) if scheme_ok else None
        row[f"{metric}_flag"] = flags[metric]

    row["n_red"] = sum(f == "red" for f in flags.values())
    row["n_amber"] = sum(f == "amber" for f in flags.values())
    row["n_flagged_metrics"] = sum(f is not None for f in flags.values())
    row["worst"] = next((f for f in FLAG_ORDER
                         if f in flags.values()), None)
    return row


def profile_table(seqs: dict[str, str], annotations: dict | None = None,
                  methods: dict | None = None,
                  compactness: dict | None = None) -> pd.DataFrame:
    """TNP profiles for a set of sequences, one row each.

    ``annotations``, ``methods`` and ``compactness`` are all keyed by sequence
    name and all optional; a name missing from any of them falls back to
    annotating on the spot, an unknown method, and no compactness.
    """
    annotations = annotations or {}
    methods = methods or {}
    compactness = compactness or {}

    rows = []
    for name, seq in seqs.items():
        try:
            rows.append(profile(name, seq,
                                annotation=annotations.get(name),
                                method=methods.get(name),
                                compactness=compactness.get(name)))
        except Exception as exc:                          # noqa: BLE001
            rows.append({"name": name, "error": str(exc)})
    return pd.DataFrame(rows)


def metric_table() -> pd.DataFrame:
    """The guidelines themselves, for the screen that shows them.

    Built from `METRICS` rather than typed out again, so the table a user
    reads and the thresholds the app flags against cannot drift apart.
    """
    rows = []
    for key, spec in METRICS.items():
        lo, hi = spec["clinical"]
        amber = "; ".join(f"{a} to {b}" for a, b in spec["amber"])
        rows.append({
            "Metric": spec["label"],
            "Clinical range": f"{lo} to {hi}",
            "Amber band": amber,
            # The lower half of the red band is printed only when a value can
            # land in it. PPC and PNC are patch scores that cannot go below
            # zero and whose clinical range starts there, so "below 0.00"
            # described a band no sequence can reach. `floor` says so per
            # metric; the condition used to be `len(spec["amber"]) > 1`, which
            # happened to agree for these two and said nothing about why.
            "Red": f"above {hi}" if spec.get("floor") is not None
                   and lo <= spec["floor"]
                   else f"below {lo} or above {hi}",
            "Computed from": spec["computed"],
            "Unit": spec["unit"],
        })
    return pd.DataFrame(rows)
