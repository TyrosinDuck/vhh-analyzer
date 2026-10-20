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
