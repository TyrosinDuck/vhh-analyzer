"""
specificity - two sequence-based readouts of how promiscuous a VHH is likely
to be, and nothing else.

    polyreactivity        a published predictor, reimplemented here
    plastic resemblance   NOT a predictor, and the module says so in its output

Both read the same input: the CDRs of a VHH, concatenated. Everything else in
this file is about getting those CDRs right, scoring them, and refusing to
return a number when the CDRs are not trustworthy.

What this module deliberately does not do
-----------------------------------------
It does not add the two scores together. They answer different questions, one
of them is calibrated and the other is not, and a composite index would hide
exactly the distinction the app is trying to make. `tests/test_specificity.py`
fails if a combined column appears.

It does not import `vhh_core`: numpy and pandas are the whole runtime, so a
notebook can `import specificity` without pulling in Biopython, matplotlib or
ANARCI. That is also why only the fitted coefficients of the polyreactivity
model are shipped, as an .npz of its non-zero weights: no scikit-learn, no
pickle, nothing that breaks on a version bump. The price is that this file
holds its own CDR anchors, see `cdrs_regex`.

Polyreactivity
--------------
Logistic regression on CDR 3-mer frequencies, from

    Harvey EP, Shin JE, Skiba MA, Nemeth GR, Hurley JD, Wellner A, Shaw AY,
    Miranda VG, Min JK, Liu CC, Marks DS, Kruse AC. An in silico method to
    assess antibody fragment polyreactivity. Nat Commun 13, 7554 (2022).
    doi:10.1038/s41467-022-35276-4

Their code and data are MIT-licensed (Copyright (c) 2022 Debora Marks Lab,
github.com/debbiemarkslab/nanobody-polyreactivity); `LICENSE` carries the
notice. The weights in `data/specificity/` were extracted from their published
model by `tools/polyreactivity/build.py`, which also checks that this
implementation reproduces their own published scores. Measured, not assumed:
over the 69 702 sequences of their low-polyreactivity pool, the decision
function computed here matches the `origFACS lr 3mers` column of their released
table to a maximum absolute difference of 1.8e-15.

**Sign convention.** The authors' decision function is positive for LOW
polyreactivity. It is negated here, so that higher means more polyreactive and
the column reads like every other flag in this app. Any comparison with their
web server or their tables has to account for the flip.

Plastic resemblance
-------------------
How much a clone's CDRs look like published plastic-binding peptides, as a
percentile against a naive library background. It is not a prediction of
binding to plastic, because no dataset labels VHHs for adsorption on
polystyrene: nothing here was fitted to, or tested against, the quantity it is
named after. The weights in `PLASTIC_WEIGHTS` are set by hand to
reproduce the physicochemical signature of each published class, and that is
all the authority they have. Every row of `plastic_resemblance` carries
`calibrated = False`, and a test enforces it.

Two scores rather than one, because the published binders fall into two
physicochemically opposite groups: arginine-rich immobilisation tags, whose
composition a polyreactivity screen already catches, and tryptophan-rich
target-unrelated peptides, which it does not. Collapsing them would average
away the only part that adds information. The reference set, with one citation
per peptide, is `data/specificity/plastic_binding_peptides.csv`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

__version__ = "2026.11.16a"

__all__ = ["AA", "CdrSet", "FLAG_PERCENTILE", "N_FEATURES", "cdrs_from_annotation",
           "cdrs_regex", "cdrs_to_kmer", "model_cdrs", "plastic_features",
           "PLASTIC_WEIGHTS", "plastic_resemblance", "polyreactivity",
           "polyreactivity_percentile",
           "polyreactivity_pool_size", "reference_peptides", "reference_report",
           "table", "__version__"]

DATA = Path(__file__).parent / "data" / "specificity"

AA = "ACDEFGHIKLMNPQRSTVWY"

# The percentile against the authors' clean pool at or above which this app
# calls a sequence polyreactive. Their own default; it is a convention, not a
# measurement, and `provenance.py` says so.
FLAG_PERCENTILE = 90.0


# --------------------------------------------------------------------------- #
# Featurisation
# --------------------------------------------------------------------------- #

def _kmer_index(k: int = 3) -> dict[str, int]:
    """The index space of the original implementation: 1-mers, then 2-mers, then 3-mers.

    The order matters and must not be touched: the shipped coefficients are
    aligned to it, column by column. 20 + 400 + 8000 = 8420 columns.
    """
    kmers = list(AA)
    last = list(AA)
    for _ in range(1, k):
        last = [p + a for p in last for a in AA]
        kmers.extend(last)
    return {km: i for i, km in enumerate(kmers)}


_KMER_TO_IDX = _kmer_index(3)
N_FEATURES = len(_KMER_TO_IDX)


def cdrs_to_kmer(seqs: Sequence[str], k: int = 3) -> np.ndarray:
    """L1-normalised k-mer frequency matrix, in the layout of Harvey et al.

    One row per sequence, `N_FEATURES` columns for k = 3. A sequence shorter
    than k, or empty, gives an all-zero row rather than an error: the caller
    decides what to do with a sequence it could not annotate, and the decision
    is not "crash in the middle of a batch".
    """
    idx = _KMER_TO_IDX if k == 3 else _kmer_index(k)
    out = np.zeros((len(seqs), len(idx)), dtype=np.float32)
    for i, raw in enumerate(seqs):
        seq = ("" if raw is None else str(raw)).upper()
        counts: dict[str, int] = {}
        for j in range(len(seq) - k + 1):
            km = seq[j:j + k]
            counts[km] = counts.get(km, 0) + 1
        total = float(sum(counts.values()))
        if total == 0.0:
            continue
        for km, c in counts.items():
            col = idx.get(km)
            if col is not None:
                out[i, col] = c / total
    return out


# --------------------------------------------------------------------------- #
# Polyreactivity
# --------------------------------------------------------------------------- #

_model: tuple[np.ndarray, float] | None = None
_poly_q: np.ndarray | None = None
_poly_n: int | None = None


def _load_model() -> tuple[np.ndarray, float]:
    """The fitted coefficients, expanded from their non-zero entries.

    The model is L1-regularised: 791 of its 8420 weights are non-zero, so the
    file stores index and value pairs in float64 rather than a dense float32
    array. That is both smaller and exact, and "exact" is the point, it is what
    lets the build script assert reproduction of the authors' own scores to
    1e-15 instead of to float32 rounding.
    """
    global _model
    if _model is None:
        z = np.load(DATA / "harvey_lr3mer_cdrs.npz")
        coef = np.zeros(int(z["n_features"]), dtype=np.float64)
        coef[z["index"]] = z["value"]
        if coef.size != N_FEATURES:
            raise ValueError(
                f"the shipped model has {coef.size} features, this "
                f"featurisation produces {N_FEATURES}")
        _model = (coef, float(z["intercept"]))
    return _model


def _poly_quantiles() -> tuple[np.ndarray, int]:
    """The clean pool's score distribution, as a quantile grid.

    Stored as 10 001 quantiles rather than as the 69 702 scores themselves: a
    percentile read off the grid differs from one read off the full sample by
    at most 0.045 and on average by 0.004, measured over the pool itself by
    `tools/polyreactivity/build.py`, which is finer than the one decimal the
    app displays, and the file is a fifth of the size.
    """
    global _poly_q, _poly_n
    if _poly_q is None:
        z = np.load(DATA / "polyreactivity_quantiles.npz")
        _poly_q = z["quantiles"].astype(np.float64)
        _poly_n = int(z["n_pool"])
    return _poly_q, int(_poly_n or 0)


def polyreactivity_pool_size() -> int:
    """How many sequences the percentile is read against, from the file itself."""
    return _poly_quantiles()[1]


def reference_report() -> dict:
    """Open every shipped file once, and say what is in them.

    The screens quote these figures instead of hard-coding them, so the pool
    size on the page cannot drift from the file the percentile is read against.
    Also the right thing to put behind `st.cache_resource`: one call warms the
    model, both backgrounds and the peptide table for the whole process.
    """
    coef, intercept = _load_model()
    grid, pool = _poly_quantiles()
    peptides = reference_peptides()
    modes = _plastic_quantiles() or {}
    return {"n_weights": int(np.count_nonzero(coef)),
            "n_features": int(coef.size),
            "intercept": intercept,
            "n_pool": pool,
            "n_quantiles": int(len(grid)),
            "n_peptides": int(len(peptides)),
            "peptide_classes": {str(k): int(v) for k, v in
                                peptides["binding_class"].value_counts().items()},
            "plastic_modes": sorted(modes),
            "flag_percentile": FLAG_PERCENTILE}


def polyreactivity(cdrs: Iterable[str], batch_size: int = 4000) -> np.ndarray:
    """Polyreactivity score of concatenated CDR strings. Higher = more polyreactive.

    Batched because the dense feature matrix is the memory cost here: 69 702
    sequences at 8420 float32 columns is 2.3 GB, which is how the first version
    of the build script died.
    """
    seqs = ["" if s is None else str(s) for s in cdrs]
    coef, intercept = _load_model()
    out = np.empty(len(seqs), dtype=np.float64)
    for i in range(0, len(seqs), batch_size):
        chunk = seqs[i:i + batch_size]
        out[i:i + len(chunk)] = cdrs_to_kmer(chunk).astype(np.float64) @ coef + intercept
    return -out


def polyreactivity_percentile(scores: Iterable[float]) -> np.ndarray:
    """Where each score falls in the authors' low-polyreactivity pool, 0 to 100.

    A sequence at 95 is more polyreactive than 95% of the sequences that
    survived the authors' selection for cleanliness. The pool is a naive
    synthetic library sorted once; it is not a set of therapeutic VHHs, so this
    is a position in a reference library, not a position among good antibodies.
    """
    q, _n = _poly_quantiles()
    x = np.asarray(list(scores), dtype=np.float64)
    if x.size == 0:
        return np.zeros(0)
    return _percentile_from_grid(q, x)


def _percentile_from_grid(grid: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Position of each value in a quantile grid, as 0 to 100.

    Clipped at the top on purpose: a value above the grid's last point lands on
    index `len(grid)`, one step past the end, which arithmetic alone turns into
    100.01. A percentile over a hundred is a nonsense the screens would print.
    """
    idx = np.searchsorted(grid, values, side="right")
    return np.clip(idx * (100.0 / (len(grid) - 1)), 0.0, 100.0)


# --------------------------------------------------------------------------- #
# CDR annotation for the model's input
# --------------------------------------------------------------------------- #
#
# The app annotates CDRs already, in `vhh_core`, and this is not a second copy
# of that for the sake of it. The two serve different purposes and do not
# agree:
#
#   `vhh_core.annotate_regex` draws boundaries for display. Its CDR1 runs to
#   the start of the FR2 motif and its CDR2 to the FR3 motif, which takes in
#   residues IMGT places in the frameworks. The model was trained on IMGT
#   windows, so feeding it those boundaries shifts the score. Measured against
#   the authors' own IMGT CDR columns over 14 192 of their sequences: not one
#   reproduced, the percentile moves by a median of 9.2 points and by 29 at the
#   99th, and the polyreactive flag of this module would land differently on
#   8.45% of sequences. That is why the display annotation is not used here,
#   and why this is not a duplicate of it.
#
# `vhh_core.annotate_anarci` does produce IMGT windows, and where ANARCI is
# installed its CDRs are what this module should be given, through
# `cdrs_from_annotation`. The anchors below are the fallback for a server
# without HMMER. Their offsets were chosen by measurement against the authors'
# IMGT columns, over the 141 204 sequences of their two pools that carry a
# complete annotation: the anchors accept 66.6% of them, and of those accepted
# 99.2% carry exactly the authors' CDRs. The rest are rejected rather than
# guessed, including an ambiguity no regex can settle, see `cdrs_regex`. On the
# accepted sequences the percentile differs from the IMGT one by 0.14 point at
# the 99th and the flag by 0.05%.

_CYS23 = re.compile(r"C[AVTGSL][AVTGSD]S")                     # ..LSCAAS..
_FR2_WYRQ = re.compile(r"W[YFHL][RKQ][QPGRH]")                 # ..WYRQ..
_FR2_END = re.compile(r"[EQKVMGR].[EQDGKR][LFVIYMW][VAIT]A")   # ..ERELVA..
_FR3_START = re.compile(r"[A-Z][YFHS]A[DENG]S[VLMA]")          # ..YYADSV..
_CYS104 = re.compile(r"[YFHWVL][YFHCSN]C")                     # ..YYC..
_FR4 = re.compile(r"W[GASV].[GAEQ]T[QHLKRPSA][VLIA]T[VLQA]S[SAP]")   # WGQGTQVTVSS
_FR4_LOOSE = re.compile(r"W[GASV].[GAEQ]T")


@dataclass(frozen=True)
class CdrSet:
    """The three CDRs of one sequence, and how they were obtained.

    `ok` is the only field a caller has to respect: when it is False the CDRs
    are not the model's input and no score may be shown for that sequence.
    `note` says why, for the user.
    """

    cdr1: str = ""
    cdr2: str = ""
    cdr3: str = ""
    annotator: str = "none"
    ok: bool = False
    note: str = ""

    @property
    def cdrs(self) -> str:
        """CDR1 + CDR2 + CDR3, which is what the model takes."""
        return f"{self.cdr1}{self.cdr2}{self.cdr3}"


def cdrs_regex(seq: str) -> CdrSet:
    """IMGT-like CDRs from sequence anchors, or `ok = False` and the reason.

    Dependency-free fallback for a server where ANARCI cannot run. It is
    deliberately strict: every anchor must be found, the CDR2 / FR3 boundary
    must not be the ambiguous case below, and the three lengths must be
    plausible, otherwise it returns no CDRs at all. A gap in the table is
    recoverable by the reader; a wrong score is not. It rejects a third of the
    authors' own library on these grounds, and that is the intended trade: of
    what it does accept, 99.2% carries exactly the IMGT CDRs.
    """
    s = re.sub(r"[^A-Za-z]", "", seq or "").upper()
    if len(s) < 60:
        return CdrSet(annotator="regex anchors", note="sequence too short to annotate")

    m23 = _CYS23.search(s, 15, 32)
    cys23 = m23.start() if m23 else s.find("C", 18, 30)
    m_fr2 = _FR2_WYRQ.search(s, 28, 70)
    if cys23 == -1 or m_fr2 is None:
        return CdrSet(annotator="regex anchors",
                      note="no conserved Cys23, or no FR2 motif (WYRQ-like) after it")
    # IMGT CDR1 is 27-38: it starts four residues after Cys23 (IMGT 23) and ends
    # two before the FR2 tryptophan (IMGT 41). Both offsets were chosen by
    # measurement, not by reading the scheme, see the comment above.
    cdr1 = s[cys23 + 4:m_fr2.start() - 2]

    m_fr2e = _FR2_END.search(s, m_fr2.end() + 3, 80)
    m_fr3 = _FR3_START.search(s, m_fr2e.end(), 100) if m_fr2e else None
    if m_fr2e is None or m_fr3 is None:
        return CdrSet(annotator="regex anchors",
                      note="no FR2 end or FR3 start motif (ERELVA / YYADSV-like)")
    # The residue the FR3 anchor opens on belongs to IMGT position 66, in the
    # framework, for every residue but one: where it is a threonine, the
    # authors' numbering assigns it to CDR2 about half the time (630 cases in
    # the first 12 000 of their pool, 49% CDR2). Sequence alone cannot tell the
    # two apart, only the numbering can, so the sequence is not scored rather
    # than scored on a coin flip.
    if s[m_fr3.start()] == "T":
        return CdrSet(annotator="regex anchors",
                      note="the CDR2 / FR3 boundary is ambiguous here, only "
                           "IMGT numbering can place it")
    cdr2 = s[m_fr2e.end():m_fr3.start()]

    m_cys104 = _CYS104.search(s, m_fr3.end())
    if m_cys104 is None:
        return CdrSet(annotator="regex anchors", note="no second conserved Cys (YYC)")
    m_fr4 = _FR4.search(s, m_cys104.end()) or _FR4_LOOSE.search(s, m_cys104.end())
    if m_fr4 is None:
        return CdrSet(annotator="regex anchors",
                      note="no FR4 motif (WGQGTQVTVSS-like) after CDR3")
    cdr3 = s[m_cys104.end():m_fr4.start()]

    bad = _implausible(cdr1, cdr2, cdr3)
    if bad:
        return CdrSet(annotator="regex anchors", note=bad)
    return CdrSet(cdr1, cdr2, cdr3, annotator="regex anchors", ok=True)


def _implausible(cdr1: str, cdr2: str, cdr3: str) -> str:
    """Length checks that catch an anchor that matched in the wrong place.

    Only for the anchors. A length is evidence about a regex match and nothing
    at all about an IMGT numbering: when ANARCI says a CDR1 is nineteen
    residues long it is nineteen residues long, and this function used to
    refuse to score such a sequence with a note blaming anchors that had not
    run. Measured over the 1 352 single domains of the SAbDab snapshot in
    `data/sabdab/`, ANARCI gives CDR1 lengths of 13, 14, 15, 17 and 19 to ten
    of them, every one of them a solved structure.

    The ranges below are therefore the observed range of those real domains,
    not the narrow range of the authors' library, which would reject a long
    CDR1 as nonsense: 4 to 19 for CDR1, 6 to 14 for CDR2, 3 to 64 for CDR3.
    What has to be caught here is a match in the wrong place, which produces a
    CDR of one or two residues or one that swallows a framework, and the
    anchors' agreement with the IMGT CDRs is unchanged by the widening, see
    `tools/polyreactivity/build.py`.
    """
    for name, cdr, lo, hi in (("CDR1", cdr1, 4, 19), ("CDR2", cdr2, 6, 14),
                              ("CDR3", cdr3, 3, 64)):
        if not (lo <= len(cdr) <= hi):
            return (f"{name} came out {len(cdr)} residues long, outside the "
                    f"{lo} to {hi} seen in solved single domains, so the "
                    f"anchors matched in the wrong place")
    return ""


def _unusable(cdr1: str, cdr2: str, cdr3: str) -> str:
    """The only thing that makes an IMGT annotation unusable: a missing CDR.

    No ceiling and no floor beyond emptiness. A sequence numbered by ANARCI
    has the CDRs the numbering gives it, and if one of them is empty there is
    nothing to featurise rather than something implausible.
    """
    for name, cdr in (("CDR1", cdr1), ("CDR2", cdr2), ("CDR3", cdr3)):
        if not cdr:
            return f"{name} is empty in the IMGT numbering, so there is nothing to score"
    return ""


def cdrs_from_annotation(annotation: dict | None) -> CdrSet:
    """The CDRs of a `vhh_core` IMGT annotation, widened to the model's CDR2.

    The caller is responsible for this being an ANARCI annotation: the regex
    annotation of `vhh_core` carries display boundaries and must not be passed
    here. `model_cdrs` is the function that makes that choice once.

    **The model's CDR2 is IMGT 55 to 65, not 56 to 65.** The app's CDR2 is the
    IMGT loop, 56 to 65, which is the right thing to display; the authors
    trained on one residue more at the N-terminus. Read off their own
    numbering, not inferred: their CDR2 column lands on IMGT 55, 56, 57, 58,
    59, 62, 63, 64, 65 for the sequences of the fixture, while their CDR1 is 27
    to 38 and their CDR3 105 to 117, which are the app's. Passing the app's
    CDR2 straight through costs the first residue of every CDR2 and scores
    every sequence slightly wrong, which is what this function existed to do
    until the ANARCI path was compared against their numbering. IMGT 55 is the
    last residue of FR2, so it comes out of the same annotation and needs no
    second call.
    """
    if not annotation:
        return CdrSet(note="no annotation")
    try:
        cdr1 = (annotation.get("CDR1") or {}).get("seq") or ""
        cdr2 = (annotation.get("CDR2") or {}).get("seq") or ""
        cdr3 = (annotation.get("CDR3") or {}).get("seq") or ""
        fr2 = (annotation.get("FR2") or {}).get("seq") or ""
    except AttributeError:
        return CdrSet(note="annotation is not in the expected form")
    if not fr2:
        return CdrSet(annotator="ANARCI / IMGT",
                      note="no FR2, so IMGT 55 is missing and the model's CDR2 "
                           "cannot be assembled")
    cdr2 = fr2[-1] + cdr2
    bad = _unusable(cdr1, cdr2, cdr3)
    if bad:
        return CdrSet(annotator="ANARCI / IMGT", note=bad)
    return CdrSet(cdr1, cdr2, cdr3, annotator="ANARCI / IMGT", ok=True)


def model_cdrs(seq: str, annotation: dict | None = None,
               annotation_is_imgt: bool = False) -> CdrSet:
    """The CDRs to score, from ANARCI when the caller has it, anchors otherwise.

    One place decides, so no screen has to remember that the display
    annotation is the wrong input for the model.
    """
    if annotation_is_imgt:
        got = cdrs_from_annotation(annotation)
        if got.ok:
            return got
    return cdrs_regex(seq)


# --------------------------------------------------------------------------- #
# Plastic resemblance
# --------------------------------------------------------------------------- #

# Kyte-Doolittle hydropathy, the published scale.
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
       "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
       "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
       "Y": -1.3, "V": 4.2}

_PS_AROMATIC_MOTIF = re.compile(r"W..W|FH..W|W..W...W")
_CATIONIC_MOTIF = re.compile(r"RR.R|R..R.R")

# Hand-set weights, public so that `provenance.py` can register them and the
# app can say on its own provenance screen that they are fitted to nothing. One
# dict per mode; the keys are columns of `plastic_features`.
PLASTIC_WEIGHTS = {
    "cationic": {"net_charge": 0.30, "arg_frac": 3.0, "basic_frac": 2.0,
                 "cationic_motif": 0.8, "cationic_kmer_share": 2.0},
    "aromatic": {"aromatic_frac": 3.0, "trp_frac": 4.0, "whfpg_frac": 1.5,
                 "ps_aromatic_motif": 1.0, "aromatic_kmer_share": 2.0,
                 "hydropathy": -0.05},
}

_peptides: pd.DataFrame | None = None
_vocab: dict[str, frozenset[str]] | None = None
_plastic_q: dict[str, np.ndarray] | None = None


def reference_peptides() -> pd.DataFrame:
    """The published plastic-binding peptides, as shipped, with their citations."""
    global _peptides
    if _peptides is None:
        _peptides = pd.read_csv(DATA / "plastic_binding_peptides.csv", comment="#")
    return _peptides.copy()


def _reference_vocab() -> dict[str, frozenset[str]]:
    """The 3-mer vocabulary of each published class, derived from the table."""
    global _vocab
    if _vocab is None:
        d = reference_peptides()
        _vocab = {}
        for cls, group in d.groupby("binding_class"):
            kmers: set[str] = set()
            for pep in group["peptide"].astype(str):
                p = pep.upper()
                kmers.update(p[i:i + 3] for i in range(len(p) - 2))
            _vocab[str(cls)] = frozenset(kmers)
    return _vocab


def _plastic_quantiles() -> dict[str, np.ndarray] | None:
    global _plastic_q
    if _plastic_q is None:
        path = DATA / "plastic_quantiles.npz"
        _plastic_q = {k: v.astype(np.float64) for k, v in np.load(path).items()} \
            if path.exists() else {}
    return _plastic_q or None


def _clean(seq: str | None) -> str:
    return re.sub(r"[^A-Z]", "", ("" if seq is None else str(seq)).upper())


def _net_charge(seq: str) -> float:
    """Net charge at neutral pH, counting histidine at a tenth of a unit.

    The crude count every peptide-composition paper uses, kept crude on
    purpose: it is a composition descriptor here, not a titration. The app's
    real charge calculation is in `vhh_core`, from pKa values, and the two are
    not interchangeable.
    """
    return (sum(c in "KR" for c in seq) + 0.1 * seq.count("H")
            - sum(c in "DE" for c in seq))


def _kmer_share(seq: str, vocab: frozenset[str]) -> float:
    n = len(seq) - 2
    if n <= 0:
        return 0.0
    return sum(seq[i:i + 3] in vocab for i in range(n)) / n


def plastic_features(cdrs: Sequence[str]) -> pd.DataFrame:
    """The per-sequence composition features both modes are built from."""
    voc = _reference_vocab()
    cationic = voc.get("cationic", frozenset())
    aromatic = voc.get("aromatic", frozenset())
    rows = []
    for raw in cdrs:
        s = _clean(raw)
        n = max(len(s), 1)
        rows.append({
            "net_charge": _net_charge(s),
            "arg_frac": s.count("R") / n,
            "basic_frac": sum(c in "KR" for c in s) / n,
            "cationic_motif": int(bool(_CATIONIC_MOTIF.search(s))),
            "cationic_kmer_share": _kmer_share(s, cationic),
            "aromatic_frac": sum(c in "WYF" for c in s) / n,
            "trp_frac": s.count("W") / n,
            "whfpg_frac": sum(c in "WHFPG" for c in s) / n,
            "hydropathy": sum(_KD.get(c, 0.0) for c in s) / n,
            "ps_aromatic_motif": int(bool(_PS_AROMATIC_MOTIF.search(s))),
            "aromatic_kmer_share": _kmer_share(s, aromatic),
        })
    return pd.DataFrame(rows, columns=["net_charge", "arg_frac", "basic_frac",
                                       "cationic_motif", "cationic_kmer_share",
                                       "aromatic_frac", "trp_frac", "whfpg_frac",
                                       "hydropathy", "ps_aromatic_motif",
                                       "aromatic_kmer_share"])


def plastic_resemblance(cdrs: Iterable[str]) -> pd.DataFrame:
    """Resemblance to the published plastic binders, by class. UNCALIBRATED.

    Returns the two raw scores, their percentile against the naive-library
    background, and a `calibrated` column that is always False. The percentile
    says how a clone compares with a naive library, nothing else: 99 means it
    looks more like the published binders of that class than 99% of a naive
    library does, not that it has a 99% chance of sticking to a plate.
    """
    seqs = [_clean(c) for c in cdrs]
    feats = plastic_features(seqs)
    raw = {mode: sum(w * feats[k].to_numpy(float) for k, w in weights.items())
           for mode, weights in PLASTIC_WEIGHTS.items()}
    out = pd.DataFrame({"cationic_mode": np.round(raw["cationic"], 3),
                        "aromatic_mode": np.round(raw["aromatic"], 3)})
    grids = _plastic_quantiles()
    for mode in ("cationic", "aromatic"):
        grid = (grids or {}).get(mode)
        if grid is None:
            out[f"{mode}_mode_pct"] = np.nan
            continue
        pct = _percentile_from_grid(grid, np.asarray(raw[mode], dtype=np.float64))
        out[f"{mode}_mode_pct"] = np.round(pct, 1)
    out["calibrated"] = False
    return pd.concat([out, feats], axis=1)


# --------------------------------------------------------------------------- #
# The table both screens draw
# --------------------------------------------------------------------------- #

def table(sequences: dict[str, str], annotations: dict | None = None,
          annotation_is_imgt: bool = False) -> pd.DataFrame:
    """One row per sequence: the CDRs used, the scores, and why a row is empty.

    A sequence whose CDRs could not be established keeps its row, with
    `annotation_ok` False, the reason in `note`, and every score left as NaN.
    Dropping it would be worse: the reader would have no way of telling a
    sequence that was not scored from one that was never submitted.
    """
    names = list(sequences)
    anns = annotations or {}
    found = [model_cdrs(sequences[n], anns.get(n), annotation_is_imgt) for n in names]

    frame = pd.DataFrame({
        "name": names,
        "annotator": [c.annotator for c in found],
        "annotation_ok": [c.ok for c in found],
        "CDR1": [c.cdr1 for c in found],
        "CDR2": [c.cdr2 for c in found],
        "CDR3": [c.cdr3 for c in found],
        "cdrs": [c.cdrs for c in found],
        "note": [c.note for c in found],
    })

    for col in ("polyreactivity_score", "polyreactivity_pct", "cationic_mode",
                "cationic_mode_pct", "aromatic_mode", "aromatic_mode_pct"):
        frame[col] = np.nan
    frame["polyreactive_flag"] = False
    frame["calibrated"] = False

    ok = frame.index[frame["annotation_ok"]].to_numpy()
    if len(ok):
        cdrs = frame.loc[ok, "cdrs"].tolist()
        score = polyreactivity(cdrs)
        pct = polyreactivity_percentile(score)
        frame.loc[ok, "polyreactivity_score"] = np.round(score, 3)
        frame.loc[ok, "polyreactivity_pct"] = np.round(pct, 1)
        frame.loc[ok, "polyreactive_flag"] = pct >= FLAG_PERCENTILE
        plastic = plastic_resemblance(cdrs)
        for col in ("cationic_mode", "cationic_mode_pct",
                    "aromatic_mode", "aromatic_mode_pct"):
            frame.loc[ok, col] = plastic[col].to_numpy()
    return frame
