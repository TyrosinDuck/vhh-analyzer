"""
vhh_core — sequence-level characterization of VHH / single-domain antibodies.

This module holds all the analysis logic, with no dependency on Streamlit or on a
notebook. It is the single source of truth shared by the Jupyter notebook and the
web application: import it from either side rather than copying functions around.

Public API
----------
clean_sequence, seq_id, parse_input            input handling
characterize, stability_flags                  physico-chemical descriptors
annotate, annotate_anarci, annotate_regex      CDR / FR annotation
register_feature, compute_features             pluggable custom descriptors
plot_secondary_structure, plot_charge,
plot_profiles, plot_msa                        figures (return matplotlib objects)
analyze                                        one-call pipeline
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # safe in a server process with no display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

__all__ = [
    "clean_sequence", "seq_id", "parse_input", "characterize", "stability_flags",
    "annotate", "annotate_anarci", "annotate_regex", "anarci_available",
    "register_feature", "compute_features", "FEATURE_REGISTRY",
    "plot_secondary_structure", "plot_charge", "plot_profiles", "plot_msa",
    "analyze", "REGIONS", "REGION_COLORS", "DEFAULTS",
]

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

DEFAULTS = {
    "ph": 7.4,                 # pH used for net-charge calculations
    "profile_window": 9,       # sliding window for per-residue profiles
    "chou_fasman_window": 5,   # sliding window for the Chou-Fasman profile
    "annotation_method": "auto",  # "auto" | "anarci" | "regex"
    "max_sequences": 50,       # guard rail for the web front-end
    "max_length": 400,         # a VHH is ~110-135 aa; reject obvious non-VHH input
}

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

REGIONS = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
REGION_COLORS = {
    "FR1": "#d9d9d9", "CDR1": "#e41a1c", "FR2": "#d9d9d9",
    "CDR2": "#377eb8", "FR3": "#d9d9d9", "CDR3": "#4daf4a", "FR4": "#d9d9d9",
}

PLOT_STYLE = {
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
}
plt.rcParams.update(PLOT_STYLE)


# --------------------------------------------------------------------------- #
# Input handling
# --------------------------------------------------------------------------- #

def clean_sequence(seq: str) -> str:
    """Normalise a raw sequence string and validate its alphabet.

    Strips whitespace, digits and gap characters, upper-cases the result and
    rejects any residue outside the 20 standard amino acids.
    """
    s = re.sub(r"[^A-Za-z]", "", seq).upper()
    if not s:
        raise ValueError("Empty sequence")
    bad = sorted(set(s) - AA20)
    if bad:
        raise ValueError(
            f"Non-standard residue(s) {', '.join(bad)} — remove or substitute them first"
        )
    return s


def seq_id(seq: str) -> str:
    """Stable 12-character fingerprint of a sequence (same sequence, same id)."""
    return hashlib.sha1(clean_sequence(seq).encode()).hexdigest()[:12]


def parse_input(text: str, max_sequences: int | None = None,
                max_length: int | None = None) -> dict[str, str]:
    """Parse FASTA or a bare sequence into {name: sequence}.

    Accepts a multi-record FASTA block, a single raw sequence (named "Sequence_1"),
    or several raw sequences separated by blank lines. Raises ValueError with a
    human-readable message on any problem, so the caller can show it directly.
    """
    max_sequences = max_sequences or DEFAULTS["max_sequences"]
    max_length = max_length or DEFAULTS["max_length"]
    text = (text or "").strip()
    if not text:
        raise ValueError("No sequence provided")

    seqs: dict[str, str] = {}

    if ">" in text:
        name, buf = None, []
        for line in text.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                raw = line[1:].strip() or f"Sequence_{len(seqs) + 1}"
                name, buf = raw.split()[0], []
            elif line:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    else:
        blocks = [b for b in re.split(r"\n\s*\n", text) if b.strip()]
        for i, block in enumerate(blocks, start=1):
            seqs[f"Sequence_{i}"] = block

    # Validate, de-duplicate names, and enforce the guard rails
    cleaned: dict[str, str] = {}
    for name, raw in seqs.items():
        try:
            s = clean_sequence(raw)
        except ValueError as exc:
            raise ValueError(f"{name}: {exc}") from None
        if len(s) > max_length:
            raise ValueError(
                f"{name}: {len(s)} aa exceeds the {max_length} aa limit "
                "(this tool expects single-domain antibodies)"
            )
        if len(s) < 60:
            raise ValueError(f"{name}: {len(s)} aa is too short to be a VHH domain")
        key, n = name, 2
        while key in cleaned:
            key, n = f"{name}_{n}", n + 1
        cleaned[key] = s

    if not cleaned:
        raise ValueError("No valid sequence found in the input")
    if len(cleaned) > max_sequences:
        raise ValueError(
            f"{len(cleaned)} sequences submitted, limit is {max_sequences} per run"
        )
    return cleaned


# --------------------------------------------------------------------------- #
# Physico-chemical characterization
# --------------------------------------------------------------------------- #

def aliphatic_index(seq: str) -> float:
    """Ikai (1980) aliphatic index: A + 2.9*V + 3.9*(I + L), in mole percent.

    Higher values correlate with thermostability.
    """
    n = len(seq)
    mol = {aa: 100.0 * seq.count(aa) / n for aa in "AVIL"}
    return mol["A"] + 2.9 * mol["V"] + 3.9 * (mol["I"] + mol["L"])


def characterize(name: str, seq: str, ph: float | None = None) -> dict:
    """Full physico-chemical descriptor set for one sequence.

    All values are theoretical, derived from the sequence alone (ProtParam
    algorithms via Biopython). ``A280_0.1pct`` is the absorbance of a 1 mg/mL
    solution, i.e. the factor to convert a measured A280 into a concentration.
    The secondary-structure fractions are amino-acid propensities, not a
    structure prediction.
    """
    ph = DEFAULTS["ph"] if ph is None else ph
    seq = clean_sequence(seq)
    pa = ProteinAnalysis(seq)
    helix, turn, sheet = pa.secondary_structure_fraction()
    eps_red, eps_ox = pa.molar_extinction_coefficient()
    mw = pa.molecular_weight()
    counts = pa.count_amino_acids()
    return {
        "name": name,
        "seq_id": seq_id(seq),
        "length": len(seq),
        "mw_da": round(mw, 2),
        "mw_kda": round(mw / 1000, 2),
        "pI": round(pa.isoelectric_point(), 2),
        "ext_coeff_reduced": eps_red,
        "ext_coeff_cystines": eps_ox,
        "A280_0.1pct": round(eps_red / mw, 3),
        "instability_index": round(pa.instability_index(), 2),
        "aliphatic_index": round(aliphatic_index(seq), 2),
        "gravy": round(pa.gravy(), 3),
        "aromaticity": round(pa.aromaticity(), 3),
        "net_charge": round(pa.charge_at_pH(ph), 2),
        "charge_pH": ph,
        "n_Cys": counts["C"],
        "n_Met": counts["M"],
        "n_Trp": counts["W"],
        "n_Asn": counts["N"],
        "helix_frac": round(helix, 3),
        "turn_frac": round(turn, 3),
        "sheet_frac": round(sheet, 3),
        "sequence": seq,
    }


def stability_flags(row: dict | pd.Series, ph: float | None = None) -> list[str]:
    """Rule-based developability warnings. Thresholds are deliberately simple —
    treat them as triage, not as a verdict."""
    ph = DEFAULTS["ph"] if ph is None else ph
    flags = []
    if row["instability_index"] >= 40:
        flags.append("Instability index ≥ 40 — predicted unstable in vitro")
    if row["gravy"] > -0.2:
        flags.append("GRAVY > −0.2 — relatively hydrophobic, aggregation-prone")
    if abs(row["net_charge"]) < 2:
        flags.append(f"|net charge| < 2 at pH {ph} — low colloidal stability")
    if abs(row["pI"] - ph) < 1:
        flags.append("Theoretical pI close to the working pH — risk of precipitation")
    if row["n_Cys"] > 2:
        flags.append(f"{row['n_Cys']} cysteines — possible unpaired or extra disulfide")
    if row["n_Met"] >= 3:
        flags.append(f"{row['n_Met']} methionines — oxidation liability")
    return flags


# --------------------------------------------------------------------------- #
# CDR / FR annotation
# --------------------------------------------------------------------------- #

# Conserved landmarks of a canonical VHH scaffold
RX_FR2 = re.compile(r"W[A-Z]RQ")                      # VHH hallmark WFRQ / WYRQ
RX_FR3 = re.compile(r"[KR][FYW][TASV][IVLFM]S")       # the RFTIS motif
RX_CYS2 = re.compile(r"[YFHV][YFHVAC]C")              # the YYC motif
RX_FR4 = re.compile(r"WG[A-Z]G[TAS][A-Z]{1,2}VTVSS")  # J-segment


def annotate_regex(seq: str) -> dict:
    """Annotate FR1–CDR3–FR4 from conserved sequence anchors (no dependency).

    Fast and dependency-free, but it assumes a canonical VHH scaffold. A heavily
    engineered or truncated construct may break an anchor; the error then names
    the anchor that could not be located.
    """
    seq = clean_sequence(seq)

    cys1 = None
    for m in re.finditer("C", seq[:30]):              # conserved Cys23
        if m.start() >= 15:
            cys1 = m.start()
    if cys1 is None:
        raise ValueError("anchor not found: conserved Cys23 (residues 16–30)")

    m2 = RX_FR2.search(seq, cys1)
    if not m2:
        raise ValueError("anchor not found: FR2 motif WxRQ")
    fr2_s, fr2_e = m2.start(), m2.start() + 14

    m3 = RX_FR3.search(seq, fr2_e)
    if not m3:
        raise ValueError("anchor not found: FR3 motif (RFTIS-like)")
    fr3_s = m3.start()

    m4 = RX_CYS2.search(seq, fr3_s)
    if not m4:
        raise ValueError("anchor not found: second conserved Cys (YYC motif)")
    cys2 = m4.end() - 1

    m5 = RX_FR4.search(seq, cys2)
    if not m5:
        raise ValueError("anchor not found: FR4 motif WGxGTxVTVSS")
    fr4_s = m5.start()

    bounds = {
        "FR1": (0, cys1 + 4), "CDR1": (cys1 + 4, fr2_s),
        "FR2": (fr2_s, fr2_e), "CDR2": (fr2_e, fr3_s),
        "FR3": (fr3_s, cys2 + 1), "CDR3": (cys2 + 1, fr4_s),
        "FR4": (fr4_s, len(seq)),
    }
    return {r: {"start": s + 1, "end": e, "seq": seq[s:e]} for r, (s, e) in bounds.items()}


IMGT_CDR = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}


def anarci_available() -> bool:
    """True when ANARCI and its HMMER backend can actually number a sequence."""
    try:
        from anarci import anarci
        probe = ("EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEG"
                 "RFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS")
        return anarci([("p", probe)], scheme="imgt", output=False)[0][0] is not None
    except Exception:
        return False


def annotate_anarci(seq: str, scheme: str = "imgt") -> dict:
    """Annotate using ANARCI numbering (reference method; needs anarci + HMMER).

    The sequence is aligned to germline HMMs and each residue receives an IMGT
    position; CDRs are then the fixed IMGT windows 27–38, 56–65 and 105–117.
    Because it is numbering-based it tolerates unusual CDR lengths.
    """
    from anarci import anarci

    seq = clean_sequence(seq)
    numbered, _details, _ = anarci([("query", seq)], scheme=scheme, output=False)
    if numbered[0] is None:
        raise ValueError("ANARCI could not number this sequence")
    numbering, start, _ = numbered[0][0]

    out = {r: {"seq": [], "pos": []} for r in REGIONS}
    idx = start
    for (num, _ins), aa in numbering:
        if aa == "-":
            continue
        region = next((c for c, (lo, hi) in IMGT_CDR.items() if lo <= num <= hi), None)
        if region is None:
            region = "FR1" if num < 27 else "FR2" if num < 56 else "FR3" if num < 105 else "FR4"
        out[region]["seq"].append(aa)
        out[region]["pos"].append(idx)
        idx += 1

    return {r: {"start": (out[r]["pos"][0] + 1) if out[r]["pos"] else None,
                "end": (out[r]["pos"][-1] + 1) if out[r]["pos"] else None,
                "seq": "".join(out[r]["seq"])} for r in REGIONS}


# --------------------------------------------------------------------------- #
# VHH hallmarks and germline assignment
# --------------------------------------------------------------------------- #

# The four FR2 positions that distinguish a camelid VHH from a conventional VH.
# IMGT 42/49/50/52 are Kabat 37/44/45/47. In a human VH3 these are V/G/L/W; the
# camelid substitutions replace the VH–VL interface with a hydrophilic surface.
HALLMARK_POSITIONS = [
    # IMGT, Kabat, camelid-preferred, VH3 consensus, note
    (42, 37, set("FYLH"), "V", "VH–VL interface; F/Y is the classic VHH signature"),
    (49, 44, set("EQDG"), "G", "least discriminating — G also occurs in VH3"),
    (50, 45, set("RCLP"), "L", "R is the classic VHH signature"),
    (52, 47, set("GLFWS"), "W", "W indicates a VH-like interface"),
]


def imgt_position_map(seq: str) -> dict[int, str]:
    """{IMGT position: residue} for one sequence. Requires ANARCI."""
    from anarci import anarci

    numbered, _d, _h = anarci([("q", clean_sequence(seq))], scheme="imgt", output=False)
    if numbered[0] is None:
        raise ValueError("ANARCI could not number this sequence")
    return {num: aa for (num, _ins), aa in numbered[0][0][0] if aa != "-"}


def hallmarks(seq: str, annotation: dict | None = None) -> dict:
    """Read the four FR2 hallmark residues and score them camelid vs VH-like.

    Uses IMGT numbering when ANARCI is available, otherwise reads the positions
    off the regex FR2 block (14 residues starting at the WxRQ motif, Kabat
    36–49). Returns the per-position detail plus a 0–4 ``camelid_score``.

    The score is an indication, not a classification: position 49 in particular
    is shared with VH3, and engineered or humanised constructs sit in between.
    """
    residues: dict[int, str] = {}
    source = "regex FR2 offsets"
    try:
        pos_map = imgt_position_map(seq)
        residues = {imgt: pos_map.get(imgt) for imgt, *_ in HALLMARK_POSITIONS}
        source = "ANARCI / IMGT numbering"
    except Exception:
        ann = annotation or annotate_regex(seq)
        fr2 = ann["FR2"]["seq"]
        for imgt, kabat, *_ in HALLMARK_POSITIONS:
            idx = kabat - 36                       # FR2 starts at Kabat 36
            residues[imgt] = fr2[idx] if 0 <= idx < len(fr2) else None

    rows, score = [], 0
    for imgt, kabat, camelid, vh3, note in HALLMARK_POSITIONS:
        aa = residues.get(imgt)
        is_camelid = bool(aa) and aa in camelid and aa != vh3
        score += int(is_camelid)
        rows.append({
            "IMGT": imgt, "Kabat": kabat, "residue": aa or "?",
            "VH3": vh3, "camelid": "/".join(sorted(camelid)),
            "verdict": "camelid" if is_camelid else ("VH-like" if aa == vh3 else "other"),
            "note": note,
        })
    return {"positions": rows, "camelid_score": score, "source": source,
            "tetrad": "".join(residues.get(p) or "?" for p, *_ in HALLMARK_POSITIONS)}


def germline(seq: str, species: tuple[str, ...] = ("alpaca", "human")) -> dict:
    """Closest germline V and J gene per species, with percent identity.

    Running the assignment once per species gives both the camelid germline the
    sequence derives from and its identity to the closest human germline — the
    latter being the usual first-pass humanness indicator.
    """
    from anarci import anarci

    seq = clean_sequence(seq)
    out = {}
    for sp in species:
        try:
            _n, details, _h = anarci([("q", seq)], scheme="imgt", output=False,
                                     assign_germline=True, allowed_species=[sp])
            d = details[0][0]
            g = d.get("germlines", {})
            (_s, v_gene), v_id = g["v_gene"]
            (_s2, j_gene), j_id = g["j_gene"]
            out[sp] = {"v_gene": v_gene, "v_identity": round(100 * v_id, 1),
                       "j_gene": j_gene, "j_identity": round(100 * j_id, 1)}
        except Exception:
            out[sp] = None
    return out


# --------------------------------------------------------------------------- #
# Sequence liabilities
# --------------------------------------------------------------------------- #

LIABILITY_PATTERNS = [
    ("N-glycosylation", r"N[^P][ST]", "high",
     "N-x-S/T sequon — potential glycosylation, heterogeneity in eukaryotic expression"),
    ("Deamidation", r"N[GS]", "high", "NG/NS — fast asparagine deamidation"),
    ("Deamidation (slow)", r"N[AHNT]", "low", "slower deamidation context"),
    ("Isomerisation", r"D[GS]", "medium", "DG/DS — aspartate isomerisation"),
    ("Fragmentation", r"DP", "medium", "DP — acid-labile peptide bond"),
    ("Oxidation (Met)", r"M", "low", "methionine — oxidation under stress"),
    ("Oxidation (Trp)", r"W", "low", "tryptophan — oxidation, especially in CDRs"),
    ("Integrin binding", r"RGD", "medium", "RGD motif — potential integrin binding"),
    ("Hydrophobic patch", r"[FILVWY]{4,}", "medium",
     "four or more consecutive hydrophobic residues"),
]

# The two cysteines forming the canonical intradomain disulfide
CANONICAL_CYS_REGIONS = {"FR1", "FR3"}


def region_of(position: int, annotation: dict | None) -> str:
    """Region containing a 1-based residue position, or '?' without annotation."""
    if not annotation:
        return "?"
    for region in REGIONS:
        r = annotation[region]
        if r["start"] and r["start"] <= position <= r["end"]:
            return region
    return "?"


def find_liabilities(seq: str, annotation: dict | None = None) -> pd.DataFrame:
    """Locate sequence liabilities, with their position and region.

    A motif in a CDR does not carry the same weight as the same motif in a
    framework: in a CDR it may alter binding, in a framework it is mostly a
    stability or heterogeneity concern. The region column is what makes this
    table actionable rather than a count.
    """
    seq = clean_sequence(seq)
    rows = []
    for label, pattern, severity, note in LIABILITY_PATTERNS:
        for m in re.finditer(f"(?=({pattern}))", seq):      # overlapping matches
            motif = m.group(1)
            start = m.start() + 1
            end = start + len(motif) - 1
            region = region_of(start, annotation)
            rows.append({
                "liability": label, "motif": motif,
                # always a string: a mixed int/str column breaks Arrow
                # serialisation in the Streamlit dataframe widget
                "position": str(start) if start == end else f"{start}–{end}",
                "region": region,
                "severity": severity,
                "in_CDR": region.startswith("CDR"),
                "note": note,
            })

    # Cysteines: the canonical pair is expected, anything else is flagged
    cys = [i + 1 for i, a in enumerate(seq) if a == "C"]
    if len(cys) > 2:
        for pos in cys:
            region = region_of(pos, annotation)
            if region not in CANONICAL_CYS_REGIONS:
                rows.append({
                    "liability": "Extra cysteine", "motif": "C", "position": str(pos),
                    "region": region, "severity": "high", "in_CDR": region.startswith("CDR"),
                    "note": "outside the canonical FR1/FR3 pair — possible free thiol "
                            "or additional disulfide",
                })
    elif len(cys) < 2:
        rows.append({
            "liability": "Missing cysteine", "motif": "-", "position": "-",
            "region": "-", "severity": "high", "in_CDR": False,
            "note": f"only {len(cys)} Cys — the canonical intradomain disulfide is incomplete",
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    order = {"high": 0, "medium": 1, "low": 2}
    return (df.assign(_s=df["severity"].map(order))
              .sort_values(["_s", "liability"]).drop(columns="_s")
              .reset_index(drop=True))


# --------------------------------------------------------------------------- #
# Pairwise identity
# --------------------------------------------------------------------------- #

def pairwise_identity(a: str, b: str) -> float:
    """Percent identity of two sequences, from a global pairwise alignment.

    Identity is expressed over the alignment length excluding positions that are
    a gap in both sequences, which is the usual convention for comparing
    antibody domains of unequal length.
    """
    from Bio import Align

    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    aln = aligner.align(clean_sequence(a), clean_sequence(b))[0]
    top, bottom = str(aln[0]), str(aln[1])
    pairs = [(x, y) for x, y in zip(top, bottom) if not (x == "-" and y == "-")]
    if not pairs:
        return float("nan")
    return round(100 * sum(x == y for x, y in pairs) / len(pairs), 1)


def identity_matrix(seqs: dict[str, str], annotations: dict | None = None,
                    region: str | None = None) -> pd.DataFrame:
    """Square matrix of pairwise percent identity.

    ``region`` restricts the comparison to one annotated region, typically
    ``"CDR3"``: two clones from the same panning can share 95% identity overall
    while having unrelated CDR3s, so the full-length figure alone is misleading.
    """
    names = list(seqs)
    if region:
        if not annotations:
            raise ValueError("Region-restricted identity needs annotations")
        subject = {n: annotations[n][region]["seq"] for n in names
                   if n in annotations and annotations[n][region]["seq"]}
        names = list(subject)
    else:
        subject = {n: seqs[n] for n in names}

    mat = pd.DataFrame(index=names, columns=names, dtype=float)
    for i, a in enumerate(names):
        mat.loc[a, a] = 100.0
        for b in names[i + 1:]:
            val = pairwise_identity(subject[a], subject[b])
            mat.loc[a, b] = mat.loc[b, a] = val
    return mat


def plot_identity(mat: pd.DataFrame, title: str = "Pairwise identity (%)"):
    """Heatmap of an identity matrix, annotated with the values."""
    n = len(mat)
    fig, ax = plt.subplots(figsize=(1.1 * n + 3, 0.9 * n + 2.4))
    data = mat.to_numpy(dtype=float)
    im = ax.imshow(data, cmap="viridis", vmin=min(50, np.nanmin(data)), vmax=100)
    ax.set_xticks(range(n), mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n), mat.index, fontsize=8)
    for i in range(n):
        for j in range(n):
            v = data[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                        color="white" if v < 85 else "black")
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.75, label="% identity")
    fig.tight_layout()
    return fig


def annotate(seq: str, method: str | None = None) -> tuple[dict, str]:
    """Dispatch to ANARCI or the regex annotator.

    ``method`` is "auto" (ANARCI with silent fallback), "anarci" (strict) or
    "regex". Returns ``(annotation, method_used)``.
    """
    method = method or DEFAULTS["annotation_method"]
    if method in ("auto", "anarci"):
        try:
            return annotate_anarci(seq), "ANARCI / IMGT"
        except Exception:
            if method == "anarci":
                raise
    return annotate_regex(seq), "Regex anchors"


# --------------------------------------------------------------------------- #
# Pluggable custom features
# --------------------------------------------------------------------------- #

FEATURE_REGISTRY: dict[str, dict] = {}


def register_feature(name: str, description: str = ""):
    """Register a custom descriptor ``f(seq, ctx) -> scalar``.

    ``ctx`` carries the context already computed for that sequence, notably
    ``ctx["annotation"]``, so a feature can address CDR3 specifically. Every
    registered feature automatically becomes a column of the results table.
    """
    def deco(fn):
        FEATURE_REGISTRY[name] = {"fn": fn, "doc": description or (fn.__doc__ or "").strip()}
        return fn
    return deco


def compute_features(seq: str, ctx: dict | None = None) -> dict:
    """Run every registered feature; a failing feature yields NaN, never a crash."""
    ctx = ctx or {}
    out = {}
    for name, entry in FEATURE_REGISTRY.items():
        try:
            out[name] = entry["fn"](seq, ctx)
        except Exception:
            out[name] = np.nan
    return out


@register_feature("charge_pH5", "Net charge at pH 5 — relevant for cation-exchange capture")
def _charge_ph5(seq, ctx):
    return round(ProteinAnalysis(seq).charge_at_pH(5.0), 2)


@register_feature("frac_KR", "Fraction of Lys + Arg — proxy for positive surface patches")
def _frac_kr(seq, ctx):
    return round((seq.count("K") + seq.count("R")) / len(seq), 3)


@register_feature("cdr3_len", "Length of CDR3, the dominant paratope loop in VHHs")
def _cdr3_len(seq, ctx):
    ann = ctx.get("annotation")
    return len(ann["CDR3"]["seq"]) if ann else np.nan


@register_feature("cdr3_gravy", "GRAVY restricted to CDR3 — hydrophobic loops drive non-specific binding")
def _cdr3_gravy(seq, ctx):
    ann = ctx.get("annotation")
    if not ann or not ann["CDR3"]["seq"]:
        return np.nan
    return round(ProteinAnalysis(ann["CDR3"]["seq"]).gravy(), 3)


@register_feature("n_motif_NG", "NG motifs — deamidation hotspots")
def _ng(seq, ctx):
    return len(re.findall(r"NG", seq))


@register_feature("n_motif_DP_DG", "DP / DG motifs — isomerisation and fragmentation hotspots")
def _dp(seq, ctx):
    return len(re.findall(r"D[PG]", seq))


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #

KD = {  # Kyte & Doolittle hydropathy scale
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}  # approximate at pH 7.4


def sliding(values: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average, edge-padded to keep the sequence length."""
    pad = window // 2
    return np.convolve(np.pad(values, pad, mode="edge"),
                       np.ones(window) / window, mode="valid")[:len(values)]


def plot_secondary_structure(df: pd.DataFrame):
    """Stacked bar chart of the helix / turn / sheet propensity fractions."""
    fig, ax = plt.subplots(figsize=(1.6 * len(df) + 3, 4))
    bottom = np.zeros(len(df))
    for col, label, color in [("sheet_frac", "Sheet", "#4c72b0"),
                              ("helix_frac", "Helix", "#dd8452"),
                              ("turn_frac", "Turn", "#55a868")]:
        vals = df[col].to_numpy()
        ax.bar(df["name"], vals, bottom=bottom, label=label, color=color, width=0.6)
        for x, (b, v) in enumerate(zip(bottom, vals)):
            if v > 0.04:
                ax.text(x, b + v / 2, f"{v:.0%}", ha="center", va="center",
                        color="white", fontsize=9, fontweight="bold")
        bottom += vals
    ax.set_ylabel("Propensity fraction")
    ax.set_title("Secondary-structure propensity (ProtParam)")
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    return fig


def plot_charge(df: pd.DataFrame, seqs: dict[str, str], ph: float | None = None):
    """Titration curves plus the net charge at the working pH."""
    ph = DEFAULTS["ph"] if ph is None else ph
    grid = np.arange(2, 12.01, 0.1)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    for name, seq in seqs.items():
        pa = ProteinAnalysis(seq)
        axes[0].plot(grid, [pa.charge_at_pH(p) for p in grid], lw=1.8, label=name)
    axes[0].axhline(0, color="k", lw=0.8, ls="--")
    axes[0].axvline(ph, color="crimson", lw=0.9, ls=":", label=f"pH {ph}")
    axes[0].set(xlabel="pH", ylabel="Net charge", title="Titration curves")
    axes[0].legend(frameon=False, fontsize=8)

    colors = ["#c44e52" if v < 0 else "#4c72b0" for v in df["net_charge"]]
    axes[1].bar(df["name"], df["net_charge"], color=colors, width=0.6)
    axes[1].axhline(0, color="k", lw=0.8)
    for x, v in enumerate(df["net_charge"]):
        axes[1].text(x, v, f"{v:+.1f}", ha="center",
                     va="bottom" if v >= 0 else "top", fontsize=9)
    axes[1].set(ylabel="Net charge", title=f"Net charge at pH {ph}")
    axes[1].tick_params(axis="x", rotation=30)
    fig.tight_layout()
    return fig


def plot_profiles(seqs: dict[str, str], annotations: dict | None = None,
                  window: int | None = None):
    """Per-residue charge and Kyte–Doolittle hydropathy, CDRs shaded."""
    w = window or DEFAULTS["profile_window"]
    fig, axes = plt.subplots(len(seqs), 1, figsize=(11, 2.4 * len(seqs)),
                             squeeze=False, sharex=True)
    for ax, (name, seq) in zip(axes[:, 0], seqs.items()):
        x = np.arange(1, len(seq) + 1)
        charge = sliding(np.array([CHARGE.get(a, 0) for a in seq], float), w)
        hydro = sliding(np.array([KD[a] for a in seq], float), w)
        ax.plot(x, hydro, color="#dd8452", lw=1.5)
        ax.axhline(0, color="k", lw=0.6)
        ax2 = ax.twinx()
        ax2.fill_between(x, charge, color="#4c72b0", alpha=0.3)
        ax2.set_ylabel("Charge", color="#4c72b0", fontsize=8)
        ax.set_ylabel(f"Hydropathy (w={w})", color="#dd8452", fontsize=8)
        ax.set_title(name, fontsize=10, loc="left")
        for a, v in ((ax, hydro), (ax2, charge)):   # align both zero lines
            lim = float(np.abs(v).max()) * 1.15 or 1.0
            a.set_ylim(-lim, lim)
        if annotations and name in annotations:
            for region in ("CDR1", "CDR2", "CDR3"):
                r = annotations[name].get(region)
                if r and r["start"]:
                    ax.axvspan(r["start"], r["end"], color=REGION_COLORS[region],
                               alpha=0.18, lw=0)
                    ax.text((r["start"] + r["end"]) / 2, ax.get_ylim()[1] * 0.78,
                            region, ha="center", fontsize=7,
                            color=REGION_COLORS[region], fontweight="bold")
    axes[-1, 0].set_xlabel("Residue position")
    fig.tight_layout()
    return fig


def build_alignment(seqs: dict[str, str], path: str | Path) -> Path:
    """Align with FAMSA (pure Python wheel, no external binary) -> aligned FASTA."""
    from pyfamsa import Aligner, Sequence

    path = Path(path)
    msa = Aligner(guide_tree="upgma").align(
        [Sequence(n.encode(), s.encode()) for n, s in seqs.items()])
    path.write_text("".join(f">{s.id.decode()}\n{s.sequence.decode()}\n" for s in msa))
    return path


def aln_columns(aligned: str, start: int, end: int) -> tuple[int, int]:
    """Map a 1-based ungapped residue range onto alignment column indices."""
    mapping, u = {}, 0
    for col, ch in enumerate(aligned, start=1):
        if ch != "-":
            u += 1
            mapping[u] = col
    return mapping[start], mapping[end]


def plot_msa(seqs: dict[str, str], annotations: dict | None = None,
             workdir: str | Path = ".", wrap_length: int = 60):
    """Render the alignment with pymsaviz, marking the CDRs of the first sequence.

    Returns ``(figure, path_to_aligned_fasta)``.
    """
    from pymsaviz import MsaViz

    aln_path = build_alignment(seqs, Path(workdir) / "alignment.fasta")
    mv = MsaViz(str(aln_path), wrap_length=wrap_length, show_consensus=True,
                show_grid=True, color_scheme="Clustal", show_count=True)

    if annotations:
        records = dict(ln.split("\n", 1) for ln in
                       aln_path.read_text().lstrip(">").split("\n>"))
        ref = list(seqs)[0]
        aligned_ref = records[ref].replace("\n", "")
        for region in ("CDR1", "CDR2", "CDR3"):
            r = annotations[ref].get(region)
            if r and r["start"]:
                c0, c1 = aln_columns(aligned_ref, r["start"], r["end"])
                mv.add_text_annotation((c0, c1), region,
                                       text_color=REGION_COLORS[region],
                                       range_color=REGION_COLORS[region])
    return mv.plotfig(), aln_path


def render_annotated_html(name: str, annotation: dict) -> str:
    """Coloured HTML rendering of one annotated sequence."""
    spans = []
    for region in REGIONS:
        s = annotation[region]["seq"]
        if not s:
            continue
        style = (f"background:{REGION_COLORS[region]};color:#fff;font-weight:600;"
                 if region.startswith("CDR") else "background:#eee;color:#333;")
        spans.append(f'<span style="{style}padding:2px 0" title="{region}">{s}</span>')
    return (f'<div style="margin-bottom:14px"><b>{name}</b><br>'
            f'<span style="font-family:ui-monospace,monospace;font-size:13px;'
            f'letter-spacing:.5px;word-break:break-all">{"".join(spans)}</span></div>')


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #

def analyze(seqs: dict[str, str], ph: float | None = None,
            method: str | None = None) -> tuple[pd.DataFrame, dict, dict]:
    """Characterize and annotate a set of sequences.

    Returns ``(results_dataframe, annotations, errors)``. A sequence whose
    annotation fails is still characterized: the reason lands in ``errors`` and
    its region columns stay empty, so one awkward construct never kills the run.
    """
    ph = DEFAULTS["ph"] if ph is None else ph
    rows, annotations, errors = [], {}, {}

    for name, raw in seqs.items():
        seq = clean_sequence(raw)
        row = characterize(name, seq, ph=ph)

        try:
            ann, used = annotate(seq, method=method)
            annotations[name] = ann
            row["annotation_method"] = used
            for region in REGIONS:
                row[f"{region}_seq"] = ann[region]["seq"]
        except Exception as exc:
            errors[name] = str(exc)
            ann = None
            row["annotation_method"] = "failed"

        row.update(compute_features(seq, {"annotation": ann, "name": name}))
        row["flags"] = "; ".join(stability_flags(row, ph=ph)) or "none"
        rows.append(row)

    return pd.DataFrame(rows), annotations, errors


# --------------------------------------------------------------------------- #
# PDF report
# --------------------------------------------------------------------------- #

SEVERITY_COLORS = {"high": "#c44e52", "medium": "#dd8452", "low": "#8c8c8c"}


def _fig_bytes(fig, dpi: int = 200) -> bytes:
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()


def build_pdf_report(path: str | Path, df: pd.DataFrame, seqs: dict[str, str],
                     annotations: dict, figures: dict | None = None,
                     ph: float | None = None, version: str = "",
                     author: str = "", with_hallmarks: bool = True,
                     with_liabilities: bool = True) -> Path:
    """Render a complete analysis report as a PDF.

    ``figures`` maps a caption to a matplotlib figure; they are embedded in the
    order given. Returns the path written.
    """
    import io
    from datetime import datetime

    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (Image, PageBreak, Paragraph, SimpleDocTemplate,
                                    Spacer, Table, TableStyle)

    ph = DEFAULTS["ph"] if ph is None else ph
    path = Path(path)
    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["BodyText"], fontSize=8.5, leading=11)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=7.4,
                          leading=10, alignment=TA_LEFT)
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontSize=16, spaceAfter=4)
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontSize=12, spaceBefore=10)
    h3 = ParagraphStyle("h3", parent=styles["Heading3"], fontSize=9.5, spaceBefore=7)
    small = ParagraphStyle("small", parent=body, fontSize=7, textColor=colors.grey)

    doc = SimpleDocTemplate(
        str(path), pagesize=A4,
        leftMargin=16 * mm, rightMargin=16 * mm,
        topMargin=14 * mm, bottomMargin=16 * mm,
        title="VHH analysis report", author=author or "VHH Analyzer")

    def table(data, widths=None, font=7.2, header=True):
        t = Table(data, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
        style = [
            ("FONTSIZE", (0, 0), (-1, -1), font),
            ("LEADING", (0, 0), (-1, -1), font + 2.5),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cccccc")),
            ("TOPPADDING", (0, 0), (-1, -1), 2.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
        ]
        if header:
            style += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef1f5")),
                      ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold")]
        t.setStyle(TableStyle(style))
        return t

    def df_table(frame, widths=None, font=7.2):
        data = [list(frame.columns)] + frame.astype(str).values.tolist()
        return table(data, widths, font)

    story = []

    # ---- Header -----------------------------------------------------------
    story.append(Paragraph("VHH analysis report", h1))
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    meta = f"Generated {stamp} · {len(seqs)} sequence(s) · net charge at pH {ph}"
    if version:
        meta += f" · VHH Analyzer {version}"
    if author:
        meta += f" · {author}"
    story += [Paragraph(meta, small), Spacer(1, 7)]

    # ---- Physico-chemical summary -----------------------------------------
    story.append(Paragraph("Physico-chemical properties", h2))
    cols = [c for c in ["name", "length", "mw_kda", "pI", "ext_coeff_reduced",
                        "A280_0.1pct", "instability_index", "aliphatic_index",
                        "gravy", "aromaticity", "net_charge", "n_Cys"]
            if c in df.columns]
    header = ["Name", "Len", "MW kDa", "pI", "e280", "A280 1mg/mL", "Instab.",
              "Aliph.", "GRAVY", "Arom.", f"Charge", "Cys"]
    summary = [header[:len(cols)]] + df[cols].astype(str).values.tolist()
    story += [table(summary, font=6.8), Spacer(1, 4),
              Paragraph("All values are theoretical, computed from the sequence alone "
                        "(ProtParam). Instability index below 40 predicts a stable "
                        "protein; GRAVY above -0.2 flags a relatively hydrophobic, "
                        "aggregation-prone candidate.", small)]

    # ---- Per-sequence sections --------------------------------------------
    for name, seq in seqs.items():
        story += [PageBreak(), Paragraph(name, h2)]
        row = df[df["name"] == name]
        row = row.iloc[0] if len(row) else None
        ann = annotations.get(name)

        if row is not None:
            facts = [["Length", f"{row['length']} aa", "MW",
                      f"{row['mw_kda']} kDa"],
                     ["Theoretical pI", f"{row['pI']}", f"Net charge (pH {ph})",
                      f"{row['net_charge']:+}"],
                     ["e280 (reduced)", f"{row['ext_coeff_reduced']} M-1cm-1",
                      "A280 at 1 mg/mL", f"{row['A280_0.1pct']}"],
                     ["Instability index", f"{row['instability_index']}",
                      "Aliphatic index", f"{row['aliphatic_index']}"],
                     ["GRAVY", f"{row['gravy']}", "Aromaticity",
                      f"{row['aromaticity']}"]]
            story.append(table(facts, widths=[34 * mm, 42 * mm, 36 * mm, 42 * mm],
                               header=False))

            flags = stability_flags(row, ph=ph)
            story.append(Paragraph("Developability flags", h3))
            story.append(Paragraph(
                "<br/>".join("• " + f for f in flags) if flags
                else "No flag raised.", body))

        # Sequence with regions
        if ann:
            story.append(Paragraph("Regions", h3))
            reg_rows = [["Region", "Range", "Length", "Sequence"]]
            for region in REGIONS:
                r = ann[region]
                if r["seq"]:
                    reg_rows.append([region, f"{r['start']}–{r['end']}",
                                     str(len(r["seq"])), r["seq"]])
            story.append(table(reg_rows, widths=[16 * mm, 20 * mm, 14 * mm, 104 * mm]))
            if row is not None and "annotation_method" in row:
                story.append(Paragraph(f"Annotation method: {row['annotation_method']}",
                                       small))

        story.append(Paragraph("Sequence", h3))
        wrapped = "<br/>".join(seq[i:i + 70] for i in range(0, len(seq), 70))
        story.append(Paragraph(wrapped, mono))

        # Hallmarks and germline
        if with_hallmarks:
            try:
                hm = hallmarks(seq, ann)
                story.append(Paragraph(
                    f"VHH hallmarks — tetrad {hm['tetrad']}, "
                    f"{hm['camelid_score']}/4 camelid ({hm['source']})", h3))
                hrows = [["IMGT", "Kabat", "Residue", "VH3", "Camelid", "Verdict"]]
                hrows += [[str(p["IMGT"]), str(p["Kabat"]), p["residue"], p["VH3"],
                           p["camelid"], p["verdict"]] for p in hm["positions"]]
                story.append(table(hrows, widths=[16 * mm, 16 * mm, 20 * mm,
                                                  16 * mm, 28 * mm, 26 * mm]))
                story.append(Paragraph(
                    "Indicative only: IMGT 49 is shared with VH3, and engineered or "
                    "humanised constructs sit between the two profiles.", small))
            except Exception as exc:
                story.append(Paragraph(f"Hallmarks unavailable: {exc}", small))

            gl = germline(seq)
            if any(gl.values()):
                grows = [["Species", "Closest V gene", "V identity", "J gene", "J identity"]]
                for sp, g in gl.items():
                    if g:
                        grows.append([sp, g["v_gene"], f"{g['v_identity']} %",
                                      g["j_gene"], f"{g['j_identity']} %"])
                story.append(Paragraph("Closest germline", h3))
                story.append(table(grows, widths=[24 * mm, 36 * mm, 26 * mm,
                                                  30 * mm, 26 * mm]))
                story.append(Paragraph(
                    "Identity to the closest human germline is a first-pass humanness "
                    "indicator, not an immunogenicity prediction.", small))

        # Liabilities
        if with_liabilities:
            liab = find_liabilities(seq, ann)
            story.append(Paragraph("Sequence liabilities", h3))
            if liab.empty:
                story.append(Paragraph("None detected.", body))
            else:
                lrows = [["Liability", "Motif", "Position", "Region", "Severity"]]
                lrows += liab[["liability", "motif", "position", "region",
                               "severity"]].astype(str).values.tolist()
                t = table(lrows, widths=[44 * mm, 20 * mm, 24 * mm, 22 * mm, 22 * mm])
                extra = []
                for i, sev in enumerate(liab["severity"], start=1):
                    extra.append(("TEXTCOLOR", (4, i), (4, i),
                                  colors.HexColor(SEVERITY_COLORS.get(sev, "#333333"))))
                t.setStyle(TableStyle(extra))
                story.append(t)
                story.append(Paragraph(
                    "A motif inside a CDR may affect binding; the same motif in a "
                    "framework is mainly a stability or heterogeneity concern. "
                    "Severity is a heuristic, not an experimental measurement.", small))

    # ---- Figures ----------------------------------------------------------
    if figures:
        from reportlab.lib.pagesizes import A4 as _A4

        avail = _A4[0] - 32 * mm
        for caption, fig in figures.items():
            if fig is None:
                continue
            story.append(PageBreak())
            story.append(Paragraph(caption, h2))
            png = _fig_bytes(fig)
            img = Image(io.BytesIO(png))
            ratio = img.imageHeight / img.imageWidth
            img.drawWidth = avail
            img.drawHeight = avail * ratio
            max_h = _A4[1] - 60 * mm
            if img.drawHeight > max_h:                 # keep tall figures on one page
                img.drawHeight, img.drawWidth = max_h, max_h / ratio
            story.append(img)

    # ---- Caveats ----------------------------------------------------------
    story += [PageBreak(), Paragraph("Notes and caveats", h2), Paragraph(
        "• Every value in this report is theoretical, derived from the amino-acid "
        "sequence alone. Confirm the isoelectric point by IEF or cIEF and the "
        "extinction coefficient against a measured A280 before using them "
        "quantitatively.<br/>"
        "• Secondary-structure fractions are composition-based propensities, not a "
        "structure prediction. Use DSSP on a model for real secondary structure.<br/>"
        "• CDR boundaries depend on the method: IMGT windows via ANARCI are tighter "
        "than the regex anchors, so CDR1 and CDR2 lengths are not comparable between "
        "the two.<br/>"
        "• Liability motifs and developability flags are triage heuristics. They "
        "indicate where to look, not what will happen.<br/>"
        "• Germline identity reflects similarity to the closest germline gene; it is "
        "not an immunogenicity prediction.", body)]

    footer_bits = [b for b in ["VHH Analyzer", version, author] if b]
    story += [Spacer(1, 8), Paragraph(" · ".join(footer_bits), small)]

    doc.build(story)
    return path


# --------------------------------------------------------------------------- #
# Chou-Fasman secondary structure profile
# --------------------------------------------------------------------------- #

# Chou & Fasman (1978) conformational parameters: helix and sheet propensities.
# A value above 100 means the residue is found in that conformation more often
# than expected by chance.
CHOU_FASMAN = {
    "A": (142, 83), "R": (98, 93), "N": (67, 89), "D": (101, 54), "C": (70, 119),
    "Q": (111, 110), "E": (151, 37), "G": (57, 75), "H": (100, 87), "I": (108, 160),
    "L": (121, 130), "K": (114, 74), "M": (145, 105), "F": (113, 138), "P": (57, 55),
    "S": (77, 75), "T": (83, 119), "W": (108, 137), "Y": (69, 147), "V": (106, 170),
}


def chou_fasman_profile(seq: str, window: int = 5) -> dict:
    """Per-residue helix and sheet propensities, plus a predicted state track.

    The state at each position is helix (H) when the smoothed helix propensity
    exceeds both the sheet propensity and 100, sheet (E) in the mirror case, and
    coil (C) otherwise. This is the classic Chou–Fasman heuristic, which is
    roughly 50–60% accurate: it shows where the sequence *tends* towards each
    conformation, it does not predict the fold. For real secondary structure,
    run DSSP on a model.
    """
    seq = clean_sequence(seq)
    helix = sliding(np.array([CHOU_FASMAN[a][0] for a in seq], float), window)
    sheet = sliding(np.array([CHOU_FASMAN[a][1] for a in seq], float), window)
    state = np.where((helix > sheet) & (helix > 100), "H",
                     np.where((sheet > helix) & (sheet > 100), "E", "C"))
    counts = {s: int((state == s).sum()) for s in "HEC"}
    return {"helix": helix, "sheet": sheet, "state": state,
            "counts": counts, "window": window,
            "fractions": {s: round(counts[s] / len(seq), 3) for s in "HEC"}}


STATE_COLORS = {"H": "#c44e52", "E": "#4c72b0", "C": "#d9d9d9"}
STATE_NAMES = {"H": "Helix", "E": "Sheet", "C": "Coil"}


def plot_chou_fasman(seq: str, name: str = "", window: int = 5,
                     annotation: dict | None = None):
    """Chou–Fasman propensity curves with a predicted-state track underneath."""
    prof = chou_fasman_profile(seq, window)
    x = np.arange(1, len(seq) + 1)

    fig, (ax, ax_track) = plt.subplots(
        2, 1, figsize=(12, 4.4), sharex=True,
        gridspec_kw={"height_ratios": [5, 1], "hspace": 0.12})

    ax.plot(x, prof["helix"], color=STATE_COLORS["H"], lw=1.4,
            label=f"Helix propensity (Chou–Fasman, w={window})")
    ax.plot(x, prof["sheet"], color=STATE_COLORS["E"], lw=1.4,
            label=f"Sheet propensity (Chou–Fasman, w={window})")
    ax.axhline(100, color="#555555", lw=0.8, ls="--")
    ax.set_ylabel("Propensity")
    title = "Predicted secondary structure (Chou–Fasman heuristic)"
    ax.set_title(f"{title} — {name}" if name else title, fontsize=11)

    for i, s in enumerate(prof["state"]):
        ax_track.axvspan(i + 0.5, i + 1.5, color=STATE_COLORS[s], lw=0)
    ax_track.set_yticks([])
    ax_track.set_ylabel("Predicted", fontsize=8, rotation=0,
                        ha="right", va="center")
    ax_track.set_xlabel("Residue position")
    ax_track.set_xlim(0.5, len(seq) + 0.5)
    for spine in ax_track.spines.values():
        spine.set_visible(True)

    if annotation:
        for region in ("CDR1", "CDR2", "CDR3"):
            r = annotation.get(region)
            if r and r["start"]:
                ax.axvspan(r["start"], r["end"], color=REGION_COLORS[region],
                           alpha=0.12, lw=0)
                lo, hi = ax.get_ylim()
                ax.text((r["start"] + r["end"]) / 2, hi - 0.03 * (hi - lo),
                        region, ha="center", va="top", fontsize=7,
                        color=REGION_COLORS[region], fontweight="bold")

    # One figure-level legend for the predicted-state track, below everything.
    # tight_layout cannot place it, so the margins are set explicitly.
    # One figure-level legend for both the curves and the track, so nothing
    # overlaps the CDR labels inside the axes.
    line_handles, line_labels = ax.get_legend_handles_labels()
    patches = [plt.Rectangle((0, 0), 1, 1, color=STATE_COLORS[s]) for s in "HEC"]
    fig.legend(line_handles + patches,
               line_labels + [f"Predicted {STATE_NAMES[s].lower()}" for s in "HEC"],
               ncol=5, fontsize=7.5, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(top=0.90, bottom=0.22, left=0.08, right=0.98)
    return fig


# --------------------------------------------------------------------------- #
# Reverse translation (E. coli codon optimisation)
# --------------------------------------------------------------------------- #

# Relative synonymous codon usage in E. coli K-12, highly expressed genes,
# ordered from most to least preferred. Rare codons (AGG/AGA, CTA, ATA, CGA)
# are kept at the end of their list so they are only used as a last resort when
# a forbidden site has to be broken.
ECOLI_CODONS = {
    "A": ["GCG", "GCC", "GCA", "GCT"], "R": ["CGT", "CGC", "CGG", "CGA"],
    "N": ["AAC", "AAT"], "D": ["GAT", "GAC"], "C": ["TGC", "TGT"],
    "Q": ["CAG", "CAA"], "E": ["GAA", "GAG"], "G": ["GGC", "GGT", "GGG", "GGA"],
    "H": ["CAT", "CAC"], "I": ["ATC", "ATT", "ATA"],
    "L": ["CTG", "TTA", "CTT", "CTC", "TTG", "CTA"], "K": ["AAA", "AAG"],
    "M": ["ATG"], "F": ["TTT", "TTC"], "P": ["CCG", "CCA", "CCT", "CCC"],
    "S": ["AGC", "TCT", "TCC", "AGT", "TCG", "TCA"],
    "T": ["ACC", "ACG", "ACT", "ACA"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTG", "GTT", "GTC", "GTA"],
}

# Restriction sites commonly used for cloning, avoided by default
DEFAULT_AVOID_SITES = {
    "NdeI": "CATATG", "XhoI": "CTCGAG", "EcoRI": "GAATTC", "BamHI": "GGATCC",
    "HindIII": "AAGCTT", "NcoI": "CCATGG", "NotI": "GCGGCCGC", "SalI": "GTCGAC",
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "SapI": "GCTCTTC",
}


def _violations(dna: str, sites: dict[str, str], max_run: int) -> list[tuple[int, str]]:
    """(position, label) of every forbidden site or homopolymer run in a sequence."""
    found = []
    for label, site in sites.items():
        for m in re.finditer(f"(?=({site}))", dna):
            found.append((m.start(), label))
        rc = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if rc != site:
            for m in re.finditer(f"(?=({rc}))", dna):
                found.append((m.start(), f"{label} (reverse)"))
    for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % ((max_run + 1,) * 4), dna):
        found.append((m.start(), f"{len(m.group(1))}-nt homopolymer"))
    return found


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def reverse_translate(protein: str, avoid_sites: dict[str, str] | None = None,
                      max_homopolymer: int = 6, add_stop: bool = True,
                      stop_codon: str = "TAA", gc_max: float = 0.60,
                      gc_min: float = 0.40) -> dict:
    """Back-translate a protein into E. coli codon-optimised DNA.

    Each residue first takes the most frequent E. coli codon; the sequence is
    then repaired where it contains a forbidden restriction site or a long
    homopolymer run, by swapping the offending codon for the next synonymous
    one. Returns the DNA plus a report of what was changed and what remains.

    This is a first-pass optimiser: it maximises codon frequency and removes
    cloning-unfriendly motifs. It does not model mRNA secondary structure near
    the ribosome binding site, codon-pair bias or tRNA availability, so for a
    construct that matters, run the output through your synthesis vendor's own
    checker before ordering.
    """
    protein = clean_sequence(protein)
    sites = DEFAULT_AVOID_SITES if avoid_sites is None else avoid_sites

    codons = [ECOLI_CODONS[aa][0] for aa in protein]
    repairs = []

    for _ in range(12):                                   # bounded repair passes
        dna = "".join(codons)
        problems = _violations(dna, sites, max_homopolymer)
        if not problems:
            break
        pos, label = problems[0]
        fixed = False
        # try every codon overlapping the offending position, nearest first
        for idx in range(max(0, pos // 3 - 1), min(len(codons), pos // 3 + 3)):
            aa = protein[idx]
            for alt in ECOLI_CODONS[aa][1:]:
                trial = codons.copy()
                trial[idx] = alt
                if len(_violations("".join(trial), sites, max_homopolymer)) < len(problems):
                    repairs.append({"position_aa": idx + 1, "residue": aa,
                                    "from": codons[idx], "to": alt, "reason": label})
                    codons = trial
                    fixed = True
                    break
            if fixed:
                break
        if not fixed:
            break

    # Codon-frequency optimisation alone pushes GC high, because the preferred
    # E. coli codons are GC-rich (GCG, CTG, CGT...). Synthesis vendors flag
    # sequences above ~60%, so trade a little codon frequency for GC balance by
    # swapping in the next-best synonymous codon where it helps.
    gc_swaps = 0
    for _ in range(len(codons)):
        dna_now = "".join(codons)
        current = _gc(dna_now)
        if gc_min <= current <= gc_max:
            break
        lowering = current > gc_max
        best = None
        for idx, aa in enumerate(protein):
            for alt in ECOLI_CODONS[aa][1:3]:
                delta = _gc(alt) - _gc(codons[idx])
                if (lowering and delta < 0) or (not lowering and delta > 0):
                    trial = codons.copy()
                    trial[idx] = alt
                    if _violations("".join(trial), sites, max_homopolymer):
                        continue
                    gain = abs(delta)
                    if best is None or gain > best[0]:
                        best = (gain, idx, alt)
        if best is None:
            break
        _g, idx, alt = best
        codons[idx] = alt
        gc_swaps += 1

    dna = "".join(codons) + (stop_codon if add_stop else "")
    remaining = _violations(dna, sites, max_homopolymer)
    gc = 100 * (dna.count("G") + dna.count("C")) / len(dna)

    return {
        "dna": dna,
        "length_nt": len(dna),
        "gc_percent": round(gc, 1),
        "repairs": repairs,
        "gc_swaps": gc_swaps,
        "remaining_issues": sorted({label for _p, label in remaining}),
        "sites_checked": sorted(sites),
        "stop_codon": stop_codon if add_stop else None,
    }


def translate_dna(dna: str) -> str:
    """Translate DNA back to protein — used to verify a back-translation."""
    from Bio.Seq import Seq

    return str(Seq(dna).translate(to_stop=True))
