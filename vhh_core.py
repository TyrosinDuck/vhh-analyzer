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
    "plot_logo",
    "analyze", "REGIONS", "REGION_COLORS", "DEFAULTS", "MIN_LENGTH",
    "classify", "classify_all", "VHH_KINDS", "KIND_VHH", "KIND_VH",
    "KIND_VL", "KIND_AMBIGUOUS", "KIND_FUSION", "KIND_OTHER",
    "pairwise_scores", "cluster_families", "plot_family_dendrogram",
    "FAMILY_DEFAULTS", "imgt_numbering", "imgt_region",
    "aggregation_profile", "aggregation_regions", "aggregation_summary",
    "plot_aggregation", "weakest_windows", "SOLUBILITY_WEIGHTS",
    "APR_THRESHOLD",
    "humanness", "human_germlines", "disulfides",
    "suggest_repairs", "apply_repairs", "LIABILITY_REPAIRS",
    "CONSERVED_IMGT",
]

# Bumped whenever the public API changes, so app.py can tell the user that a
# module is out of date instead of crashing on a missing attribute.
__version__ = "2026.10.21g"

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

MIN_LENGTH = 25          # below this ProtParam figures stop meaning much

DEFAULTS = {
    "ph": 7.4,                 # pH used for net-charge calculations
    "profile_window": 9,       # sliding window for per-residue profiles
    "chou_fasman_window": 5,   # sliding window for the Chou-Fasman profile
    "annotation_method": "auto",  # "auto" | "anarci" | "regex"
    "max_sequences": 50,       # guard rail for the web front-end
    "max_length": 1500,        # any protein is accepted, not only VHHs; this is
                               # a guard rail for the shared web instance
    "max_total_residues": 12000,   # ceiling on one submission, all sequences
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
                f"{name}: {len(s)} aa exceeds the {max_length} aa limit"
            )
        if len(s) < MIN_LENGTH:
            raise ValueError(
                f"{name}: {len(s)} aa is below the {MIN_LENGTH} aa minimum")
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
    total = sum(len(v) for v in cleaned.values())
    ceiling = DEFAULTS["max_total_residues"]
    if total > ceiling:
        # one long sequence is cheap, fifty of them is not, and the Streamlit
        # instance is shared with whoever else is using the app right now
        raise ValueError(
            f"{total} residues submitted in total, limit is {ceiling} per run — "
            "split the submission or shorten the longest sequences")
    return cleaned


# ---------------------------------------------------------------------------
# What kind of sequence is this?
# ---------------------------------------------------------------------------
# The app accepts any protein, not only VHHs, and has to decide per sequence
# which analyses apply. There is no structure here, so the decision is made from
# the sequence alone:
#
#   1. ANARCI numbering succeeds  -> it is an immunoglobulin variable domain,
#      and ANARCI also reports whether the chain is heavy (H) or light (K/L).
#   2. For a heavy-chain domain, the FR2 hallmark tetrad separates VHH from a
#      conventional VH.
#   3. ANARCI fails -> not a V domain. Everything generic still applies; the
#      VHH-specific analyses do not.
#
# A V domain that covers only part of a longer sequence is reported as such
# rather than silently treated as a VHH: a 900-residue protein carrying a VHH
# domain is a fusion, and running CDR analysis on the whole thing would be
# meaningless.

KIND_VHH = "VHH"
KIND_VH = "Conventional VH"
KIND_VL = "Light-chain V domain"
KIND_AMBIGUOUS = "V domain, VHH/VH ambiguous"
KIND_FUSION = "Contains a V domain"
KIND_OTHER = "Other protein"

VHH_KINDS = {KIND_VHH, KIND_AMBIGUOUS}
DOMAIN_COVERAGE = 0.80      # fraction of the sequence the V domain must cover


def classify(seq: str, method: str | None = None) -> dict:
    """Decide what a sequence is, and therefore which analyses apply.

    Returns ``kind`` (one of the ``KIND_*`` constants), ``is_vhh``, the
    ``chain`` ANARCI reported, the ``hallmark_score`` out of 4, the ``coverage``
    of the sequence by the V domain, the ``method`` the decision came from, and
    a one-line ``note`` for display.
    """
    seq = clean_sequence(seq)
    out = {"length": len(seq), "chain": None, "hallmark_score": None,
           "coverage": None, "evalue": None, "method": "ANARCI"}

    try:
        from anarci import anarci
        numbered, details, _ = anarci([("q", seq)], scheme="imgt", output=False)
    except Exception:
        # No ANARCI: fall back to the regex annotator, which only ever
        # recognises a VHH-like domain. It cannot tell VH from VHH, nor find a
        # light chain, so the verdict is deliberately coarser.
        out["method"] = "regex anchors"
        try:
            annotate_regex(seq)
            score = hallmarks(seq).get("camelid_score")
            out["hallmark_score"] = score
            out["kind"] = KIND_VHH if (score or 0) >= 2 else KIND_AMBIGUOUS
            out["note"] = ("recognised by regex anchors; without ANARCI this "
                           "cannot be told apart from a conventional VH")
        except Exception:
            out["kind"] = KIND_OTHER
            out["note"] = ("no VHH anchors found; generic analyses only "
                           "(install ANARCI for a reliable verdict)")
        out["is_vhh"] = out["kind"] in VHH_KINDS
        return out

    if numbered[0] is None:
        out.update({"kind": KIND_OTHER, "is_vhh": False,
                    "note": "ANARCI finds no immunoglobulin variable domain — "
                            "generic analyses only"})
        return out

    _numbering, start, end = numbered[0][0]
    detail = details[0][0]
    chain = detail.get("chain_type")
    coverage = (end - start + 1) / len(seq)
    out.update({"chain": chain, "coverage": round(coverage, 3),
                "evalue": detail.get("evalue")})

    if coverage < DOMAIN_COVERAGE:
        out.update({
            "kind": KIND_FUSION, "is_vhh": False,
            "note": f"a {chain}-chain V domain spans residues {start + 1}-{end + 1} "
                    f"({coverage:.0%} of the sequence) — this looks like a fusion "
                    f"or a multi-domain protein, so the VHH analyses are skipped",
        })
        return out

    if chain in ("K", "L"):
        out.update({"kind": KIND_VL, "is_vhh": False,
                    "note": f"ANARCI numbers this as a {chain} light-chain "
                            "variable domain, not a VHH"})
        return out

    score = hallmarks(seq).get("camelid_score")
    out["hallmark_score"] = score
    if (score or 0) >= 2:
        out.update({"kind": KIND_VHH, "is_vhh": True,
                    "note": f"heavy-chain V domain with {score}/4 camelid "
                            "hallmarks"})
    elif score == 1:
        out.update({"kind": KIND_AMBIGUOUS, "is_vhh": True,
                    "note": "heavy-chain V domain with only 1/4 camelid "
                            "hallmarks — humanised VHH or VH-like construct; "
                            "VHH analyses are still run, read them with that in "
                            "mind"})
    else:
        out.update({"kind": KIND_VH, "is_vhh": False,
                    "note": "heavy-chain V domain with a VH-like FR2 tetrad — "
                            "a conventional VH, so the VHH analyses are skipped"})
    return out


def classify_all(seqs: dict[str, str], method: str | None = None) -> pd.DataFrame:
    """``classify`` over a batch, as a table."""
    rows = []
    for name, seq in seqs.items():
        rows.append({"name": name, **classify(seq, method=method)})
    return pd.DataFrame(rows)


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

# Scoring scale for chemical liabilities.
#
# The base scores and the Low/Medium/High thresholds come from Jonathan's
# in-silico liability notebook, where each motif is weighted by its solvent
# exposure (SASA) and the model confidence (pLDDT) taken from an AlphaFold
# structure. This module is sequence-only, so both of those factors are absent
# and the score is the base score alone. Two consequences worth keeping in mind:
# the categories compress towards the top (an NG scores 10 whether it is buried
# or fully exposed), and the ranking *between* sequences is more informative
# than the absolute category of any one motif.
#
# Entries marked source="added" are not from that notebook; their scores are a
# judgement call and are meant to be edited.

RISK_THRESHOLDS = (2, 6)          # < 2 Low, < 6 Medium, >= 6 High

LIABILITY_RULES = [
    # (label, regex, base score, source, note)
    ("N-glycosylation", r"N[^P][ST]", 15, "notebook",
     "N-x-S/T sequon — glycosylation and heterogeneity in eukaryotic expression"),
    ("Deamidation (high)", r"N[G]", 10, "notebook",
     "NG — fast asparagine deamidation"),
    ("Deamidation (moderate)", r"N[STN]", 5, "notebook",
     "NS/NT/NN — slower deamidation context"),
    ("Isomerisation (high)", r"D[G]", 10, "notebook",
     "DG — fast aspartate isomerisation"),
    ("Isomerisation (moderate)", r"D[STD]", 5, "notebook",
     "DS/DT/DD — slower isomerisation context"),
    ("Oxidation (Met)", r"M", 4, "notebook", "methionine — oxidation under stress"),
    ("Oxidation (Trp)", r"W", 3, "notebook",
     "tryptophan — oxidation, especially in CDRs"),
    ("Fragmentation", r"DP", 3, "added", "DP — acid-labile peptide bond"),
    ("Integrin binding", r"RGD", 2, "added",
     "RGD motif — potential integrin binding"),
    ("Hydrophobic patch", r"[FILVWY]{4,}", 2, "added",
     "four or more consecutive hydrophobic residues"),
]

# Proteases, straight from the notebook. Scores are per-site and deliberately
# small: a VHH has many tryptic sites and that on its own means little.
PROTEASE_RULES = [
    ("Trypsin-like", r"[KR](?![P])", 0.5, "cleaves after K/R, not before P"),
    ("Chymotrypsin-like", r"[FWY](?![P])", 0.5, "cleaves after F/W/Y, not before P"),
    ("GluC-like", r"E", 0.25, "cleaves after E"),
    ("AspN-like", r"(?=.D)", 0.25, "cleaves before D"),
    ("Factor Xa-like", r"[IVL][ED]GR", 5, "IEGR/VDGR-type recognition site"),
    ("Enterokinase-like", r"DDDDK", 5, "DDDDK recognition site"),
]

# The canonical intradomain disulfide of a VHH; a Cys anywhere else is a
# candidate free thiol. With a structure you would check the SG-SG distances
# instead, which is what the notebook does.
CANONICAL_CYS_REGIONS = {"FR1", "FR3"}



# Every liabilities table carries these columns, even when it is empty: callers
# do arithmetic on `score` and a bare empty DataFrame has no columns at all,
# which turns "this sequence is clean" into a KeyError.
LIABILITY_COLUMNS = ["liability", "motif", "position", "range", "region",
                     "in_CDR", "score", "risk", "source", "note"]

def risk_category(score: float) -> str:
    """Low / Medium / High, on the notebook's thresholds."""
    low, high = RISK_THRESHOLDS
    return "Low" if score < low else ("Medium" if score < high else "High")


def region_of(position: int, annotation: dict | None) -> str:
    """Region containing a 1-based residue position, or '?' without annotation."""
    if not annotation:
        return "?"
    for region in REGIONS:
        r = annotation[region]
        if r["start"] and r["start"] <= position <= r["end"]:
            return region
    return "?"


def find_liabilities(seq: str, annotation: dict | None = None,
                     include_added: bool = True) -> pd.DataFrame:
    """Locate and score sequence liabilities, with position and region.

    A motif in a CDR does not carry the same weight as the same motif in a
    framework: in a CDR it may alter binding, in a framework it is mostly a
    stability or heterogeneity concern. The region column is what makes this
    table actionable rather than a count.
    """
    seq = clean_sequence(seq)
    rows = []

    for label, pattern, base, source, note in LIABILITY_RULES:
        if source == "added" and not include_added:
            continue
        for m in re.finditer(f"(?=({pattern}))", seq):
            motif = m.group(1)
            start = m.start() + 1
            end = start + max(len(motif), 1) - 1
            region = region_of(start, annotation)
            rows.append({
                "liability": label, "motif": motif, "position": start,
                "range": str(start) if start == end else f"{start}–{end}",
                "region": region, "in_CDR": region.startswith("CDR"),
                "score": float(base), "risk": risk_category(base),
                "source": source, "note": note,
            })

    # Cysteines. The canonical-pair logic below is specific to a VHH: it assumes
    # exactly one intradomain disulfide, between FR1 and FR3. Applied to any other
    # protein it is simply wrong — lysozyme's eight cysteines form four ordinary
    # disulfides and are not eight free thiols. So without an annotation (which is
    # what "not a VHH" looks like here) fall back to the only thing that holds for
    # any protein: an odd number of cysteines guarantees one unpaired thiol, an
    # even number tells you nothing without a structure.
    cys = [i + 1 for i, a in enumerate(seq) if a == "C"]
    if annotation is None:
        if len(cys) % 2 == 1:
            rows.append({
                "liability": "Unpaired cysteine", "motif": "C",
                "position": cys[-1], "range": str(cys[-1]), "region": "-",
                "in_CDR": False, "score": 8.0, "risk": risk_category(8),
                "source": "added",
                "note": f"{len(cys)} cysteines — an odd count means at least one "
                        "free thiol",
            })
    elif len(cys) > 2:
        for pos in cys:
            region = region_of(pos, annotation)
            if region not in CANONICAL_CYS_REGIONS:
                rows.append({
                    "liability": "Free cysteine", "motif": "C", "position": pos,
                    "range": str(pos), "region": region,
                    "in_CDR": region.startswith("CDR"), "score": 8.0,
                    "risk": risk_category(8), "source": "notebook",
                    "note": "outside the canonical FR1/FR3 pair — possible free "
                            "thiol or additional disulfide",
                })
    elif len(cys) < 2:
        rows.append({
            "liability": "Missing cysteine", "motif": "-", "position": 0,
            "range": "-", "region": "-", "in_CDR": False, "score": 8.0,
            "risk": risk_category(8), "source": "added",
            "note": f"only {len(cys)} Cys — the canonical intradomain disulfide is "
                    "incomplete",
        })

    df = pd.DataFrame(rows, columns=LIABILITY_COLUMNS)
    if df.empty:
        return df
    order = {"High": 0, "Medium": 1, "Low": 2}
    return (df.assign(_s=df["risk"].map(order))
              .sort_values(["_s", "position"]).drop(columns="_s")
              .reset_index(drop=True))


def find_protease_sites(seq: str, annotation: dict | None = None) -> pd.DataFrame:
    """Predicted protease cleavage sites, scored on the notebook's scale.

    Every protein has tryptic sites; the count on its own means little. What is
    worth looking at is a high-scoring recognition site, or a cluster of sites
    inside a CDR.
    """
    seq = clean_sequence(seq)
    rows = []
    for label, pattern, base, note in PROTEASE_RULES:
        for m in re.finditer(f"(?=({pattern}))", seq):
            motif = m.group(1) or seq[m.start():m.start() + 2]
            pos = m.start() + 1
            if label == "AspN-like":
                pos = m.start() + 2
                motif = seq[m.start():m.start() + 2]
            region = region_of(pos, annotation)
            rows.append({"protease": label, "position": pos, "motif": motif,
                         "region": region, "in_CDR": region.startswith("CDR"),
                         "score": float(base), "note": note})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["score", "position"],
                          ascending=[False, True]).reset_index(drop=True)


def liability_summary(name: str, liabilities: pd.DataFrame,
                      proteases: pd.DataFrame | None = None) -> dict:
    """One-line summary per sequence: total score and the risk breakdown."""
    if liabilities is None or liabilities.empty:
        base = {"name": name, "total_score": 0.0, "n_liabilities": 0,
                "n_high": 0, "n_medium": 0, "n_low": 0, "n_in_CDR": 0}
    else:
        base = {
            "name": name,
            "total_score": round(float(liabilities["score"].sum()), 2),
            "n_liabilities": int(len(liabilities)),
            "n_high": int((liabilities["risk"] == "High").sum()),
            "n_medium": int((liabilities["risk"] == "Medium").sum()),
            "n_low": int((liabilities["risk"] == "Low").sum()),
            "n_in_CDR": int(liabilities["in_CDR"].sum()),
        }
    if proteases is not None and not proteases.empty:
        base["protease_score"] = round(float(proteases["score"].sum()), 2)
        base["n_protease_sites"] = int(len(proteases))
    return base


LIABILITY_COLORS = {
    "N-glycosylation": "#7b3294", "Deamidation (high)": "#d7191c",
    "Deamidation (moderate)": "#fdae61", "Isomerisation (high)": "#2c7bb6",
    "Isomerisation (moderate)": "#abd9e9", "Oxidation (Met)": "#1a9641",
    "Oxidation (Trp)": "#8c6d31", "Fragmentation": "#e377c2",
    "Integrin binding": "#17becf", "Hydrophobic patch": "#999999",
    "Free cysteine": "#000000", "Missing cysteine": "#555555",
}


def plot_liability_map(liabilities: dict[str, pd.DataFrame],
                       annotations: dict | None = None, lengths: dict | None = None):
    """Liabilities along the sequence, one row per VHH, bubble size by score.

    This is the view that answers the question you actually ask of a panel:
    *where* do the problems sit, and do they fall in a CDR.
    """
    names = list(liabilities)
    if not names:
        return None
    fig, ax = plt.subplots(figsize=(12, 0.72 * len(names) + 2.2))

    # CDR bands per row, so a liability landing in a loop is visible at a glance
    if annotations:
        for y, name in enumerate(names):
            ann = annotations.get(name)
            if not ann:
                continue
            for region in ("CDR1", "CDR2", "CDR3"):
                r = ann.get(region)
                if r and r["start"]:
                    ax.add_patch(plt.Rectangle(
                        (r["start"], y - 0.34), r["end"] - r["start"], 0.68,
                        facecolor=REGION_COLORS[region], alpha=0.16, lw=0, zorder=0))

    seen = []
    for y, name in enumerate(names):
        df = liabilities[name]
        if df is None or df.empty:
            continue
        for label, group in df[df["position"] > 0].groupby("liability", sort=False):
            ax.scatter(group["position"], [y] * len(group),
                       s=group["score"] * 22 + 18,
                       color=LIABILITY_COLORS.get(label, "#777777"),
                       alpha=0.75, edgecolor="white", linewidth=0.6,
                       label=label if label not in seen else None, zorder=3)
            if label not in seen:
                seen.append(label)

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_ylim(len(names) - 0.45, -0.55)
    ax.set_xlim(0, max((lengths or {}).values(), default=140) + 4)
    ax.set_xlabel("Residue position")
    ax.set_title("Chemical liabilities along the sequence "
                 "(bubble size = score; shaded bands = CDRs)", fontsize=11)
    ax.grid(alpha=0.25, axis="x")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8,
              frameon=False, labelspacing=0.9)
    fig.tight_layout()
    return fig


def plot_liability_profile(liabilities: dict[str, pd.DataFrame]):
    """Stacked bar of the score contributed by each liability type, per VHH."""
    names = list(liabilities)
    frames = []
    for name in names:
        df = liabilities[name]
        if df is None or df.empty:
            continue
        grouped = df.groupby("liability")["score"].sum()
        grouped.name = name
        frames.append(grouped)
    if not frames:
        return None
    matrix = pd.concat(frames, axis=1).fillna(0).T.reindex(names).fillna(0)

    fig, ax = plt.subplots(figsize=(1.5 * len(names) + 4, 4.4))
    bottom = np.zeros(len(matrix))
    for column in matrix.columns:
        values = matrix[column].to_numpy()
        ax.bar(matrix.index, values, bottom=bottom, label=column,
               color=LIABILITY_COLORS.get(column, "#777777"), width=0.6)
        bottom += values
    for x, total in enumerate(bottom):
        ax.text(x, total, f"{total:.0f}", ha="center", va="bottom", fontsize=9,
                fontweight="bold")
    ax.set_ylabel("Summed liability score")
    ax.set_title("Liability profile", fontsize=11)
    ax.tick_params(axis="x", rotation=30)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    return fig


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

    # The CDR bands refer to the first sequence, which may have no annotation at
    # all now that the app accepts proteins that are not VHHs. Fall back to the
    # first sequence that does have one, and draw no bands if none does.
    ref = next((n for n in seqs if n in (annotations or {})), None)
    if annotations and ref:
        records = dict(ln.split("\n", 1) for ln in
                       aln_path.read_text().lstrip(">").split("\n>"))
        aligned_ref = records[ref].replace("\n", "")
        for region in ("CDR1", "CDR2", "CDR3"):
            r = annotations[ref].get(region)
            if r and r["start"]:
                c0, c1 = aln_columns(aligned_ref, r["start"], r["end"])
                mv.add_text_annotation((c0, c1), region,
                                       text_color=REGION_COLORS[region],
                                       range_color=REGION_COLORS[region])
    return mv.plotfig(), aln_path


def plot_logo(seqs: dict[str, str], annotations: dict | None = None,
              to_type: str = "probability", wrap_length: int = 65,
              color_scheme: str = "chemistry", workdir: str | Path = ".",
              figwidth: float = 16.0):
    """Sequence logo of the alignment, wrapped over several rows.

    Reproduces the logo from the clinical-VHH alignment notebook: the aligned
    sequences go through ``logomaker.alignment_to_matrix`` with gaps and X
    ignored and no pseudocount, amino-acid columns never observed are dropped,
    and the logo is wrapped at ``wrap_length`` columns per row with ticks every
    ten positions carrying absolute alignment coordinates.

    ``to_type`` is ``"probability"`` (height = frequency, the notebook's choice,
    which reads directly as composition) or ``"information"`` (height in bits,
    the classic Schneider logo, where a column of mixed residues is short even
    if one residue is commonest).

    Returns ``(figure, matrix, path_to_aligned_fasta)``.
    """
    import logomaker

    aln_path = build_alignment(seqs, Path(workdir) / "alignment.fasta")
    records = dict(ln.split("\n", 1) for ln in
                   aln_path.read_text().lstrip(">").split("\n>"))
    records = {k: v.replace("\n", "") for k, v in records.items()}
    aligned = list(records.values())

    matrix = logomaker.alignment_to_matrix(aligned, to_type=to_type,
                                           characters_to_ignore="-Xx.",
                                           pseudocount=0)
    # drop amino acids that never occur, so the colour key stays readable
    matrix = matrix.loc[:, (matrix > 0).any(axis=0)]

    n_pos = matrix.shape[0]
    n_rows = max(1, (n_pos + wrap_length - 1) // wrap_length)
    fig, axes = plt.subplots(n_rows, 1, figsize=(figwidth, 2.0 * n_rows),
                             squeeze=False)

    # CDR spans of the first sequence, in alignment columns, for the shading
    spans = {}
    ref = next((n for n in seqs if n in (annotations or {})), None)
    if annotations and ref:
        ref_aligned = records.get(ref, "")
        for region in ("CDR1", "CDR2", "CDR3"):
            info = (annotations.get(ref) or {}).get(region) or {}
            if ref_aligned and info.get("start") and info.get("end"):
                try:
                    spans[region] = aln_columns(ref_aligned, info["start"],
                                                info["end"])
                except KeyError:
                    pass

    for row in range(n_rows):
        start = row * wrap_length
        end = min((row + 1) * wrap_length, n_pos)
        ax = axes[row, 0]
        block = matrix.iloc[start:end].copy()
        block.index = range(len(block))         # logomaker wants 0-based indices
        logomaker.Logo(block, ax=ax, color_scheme=color_scheme)

        for region, (c0, c1) in spans.items():
            lo, hi = max(c0, start + 1), min(c1, end)
            if lo <= hi:
                ax.axvspan(lo - start - 1.5, hi - start - 0.5,
                           color=REGION_COLORS[region], alpha=0.18, lw=0,
                           zorder=0)
                ax.text((lo + hi) / 2 - start - 1, 1.02, region,
                        transform=ax.get_xaxis_transform(), ha="center",
                        fontsize=7.5, color="0.3")

        ax.set_xlim(-0.5, len(block) - 0.5)
        ticks = range(0, len(block), 10)
        ax.set_xticks(list(ticks))
        ax.set_xticklabels([str(start + p + 1) for p in ticks])
        ax.set_ylabel("Bits" if to_type == "information" else "Freq")
    axes[-1, 0].set_xlabel("Alignment position")
    fig.tight_layout()
    return fig, matrix, aln_path


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

RISK_COLORS = {"High": "#c44e52", "Medium": "#dd8452", "Low": "#8c8c8c"}


def _fig_bytes(fig, dpi: int = 150) -> bytes:
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()



# ---------------------------------------------------------------------------
# Clonal families
# ---------------------------------------------------------------------------
# After a panning you rarely care how good each clone is before you know how
# many *different* clones you have: a plate of 16 picked colonies is routinely
# the same three lineages over and over. Grouping on CDR3 answers that, because
# CDR3 is the region that differs between independent lineages while the
# frameworks stay nearly constant.
#
# The method is the one used in the Sanger screening notebook, kept identical so
# the app and the notebook give the same families for the same reads:
#   * identity and similarity from a global BLOSUM62 alignment (gaps -11/-1),
#     over a denominator of max(len(a), len(b))
#   * distance = 1 - identity, hierarchical clustering, cut at 1 - cutoff
#
# Note the denominator differs from ``pairwise_identity`` above, which divides by
# the alignment length. The two agree for equal-length CDR3s and differ by a few
# percent otherwise; this one is kept here for continuity with the notebook.

FAMILY_DEFAULTS = {"cutoff": 80.0, "linkage": "complete",
                   "require_same_length": False}

_FAMILY_ALIGNER = None
_BLOSUM62 = None


def _family_aligner():
    """Global BLOSUM62 aligner, built once — constructing it is not free."""
    global _FAMILY_ALIGNER, _BLOSUM62
    if _FAMILY_ALIGNER is None:
        from Bio import Align
        from Bio.Align import substitution_matrices
        _BLOSUM62 = substitution_matrices.load("BLOSUM62")
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = _BLOSUM62
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        aligner.mode = "global"
        _FAMILY_ALIGNER = aligner
    return _FAMILY_ALIGNER


def _blosum(x: str, y: str) -> float:
    try:
        return float(_BLOSUM62[x, y])
    except (KeyError, IndexError):
        try:
            return float(_BLOSUM62[y, x])
        except (KeyError, IndexError):
            return -4.0


def pairwise_scores(a: str, b: str) -> tuple[float, float]:
    """``(identity, similarity)`` in percent.

    *Identity* counts strictly identical residues; *similarity* also counts
    conservative substitutions (BLOSUM62 > 0, what BLAST calls positives). Both
    are divided by the length of the longer sequence, so a CDR3 that is a
    truncation of another cannot score 100%.
    """
    a, b = clean_sequence(a), clean_sequence(b)
    if not a or not b:
        return float("nan"), float("nan")
    aligner = _family_aligner()
    aln = aligner.align(a, b)[0]
    top, bottom = str(aln[0]), str(aln[1])
    ident = sim = 0
    for x, y in zip(top, bottom):
        if x == "-" or y == "-":
            continue
        if x == y:
            ident += 1
            sim += 1
        elif _blosum(x, y) > 0:
            sim += 1
    denom = max(len(a), len(b))
    return round(100 * ident / denom, 1), round(100 * sim / denom, 1)


def cluster_families(seqs: dict[str, str], annotations: dict,
                     cutoff: float = 80.0, linkage_method: str = "complete",
                     require_same_length: bool = False,
                     region: str = "CDR3") -> dict:
    """Group sequences into clonal families on their CDR3.

    ``cutoff`` is a percent identity: two sequences at or above it belong to the
    same family. ``linkage_method`` decides what that means for a family of more
    than two members — ``"complete"`` guarantees every pair inside a family
    meets the cutoff, ``"average"`` merges on the mean distance and is more
    permissive, so a family can then contain pairs below the cutoff. The returned
    diagnostic reports each family's worst internal pair, which is how you see
    that happening.

    ``require_same_length`` forces CDR3s of different lengths into different
    families, whatever their alignment. Off by default: a single-residue indel in
    a CDR3 is a plausible sibling, not a new lineage.

    Returns a dict with ``assignment`` (one row per sequence), ``families`` (one
    row per family), ``identity`` and ``similarity`` matrices, the scipy
    ``linkage`` matrix for the dendrogram, and ``skipped`` — sequences left out
    because the region could not be annotated.
    """
    import numpy as np

    subject, skipped = {}, []
    for name in seqs:
        loop = ((annotations.get(name) or {}).get(region) or {}).get("seq") or ""
        if loop:
            subject[name] = loop
        else:
            skipped.append(name)

    names = list(subject)
    n = len(names)
    empty = pd.DataFrame()
    if n == 0:
        return {"assignment": empty, "families": empty, "identity": empty,
                "similarity": empty, "linkage": None, "skipped": skipped,
                "cutoff": cutoff, "linkage_method": linkage_method,
                "region": region, "n_families": 0}

    ident = np.eye(n) * 100.0
    simil = np.eye(n) * 100.0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = subject[names[i]], subject[names[j]]
            v, w = pairwise_scores(a, b)
            if require_same_length and len(a) != len(b):
                v = w = 0.0
            ident[i, j] = ident[j, i] = v
            simil[i, j] = simil[j, i] = w

    ident_df = pd.DataFrame(ident, index=names, columns=names)
    simil_df = pd.DataFrame(simil, index=names, columns=names)

    link = None
    if n >= 2:
        from scipy.cluster.hierarchy import fcluster, linkage as scipy_linkage
        from scipy.spatial.distance import squareform

        dist = 1.0 - ident / 100.0
        np.fill_diagonal(dist, 0.0)
        dist = np.clip((dist + dist.T) / 2, 0.0, None)   # enforce exact symmetry
        link = scipy_linkage(squareform(dist, checks=False), method=linkage_method)
        labels = fcluster(link, t=1.0 - cutoff / 100.0, criterion="distance")
    else:
        labels = np.array([1])

    # number families by decreasing size, so F01 is always the dominant lineage
    order = pd.Series(labels).value_counts().index
    remap = {old: f"F{i + 1:02d}" for i, old in enumerate(order)}
    family_of = {names[i]: remap[labels[i]] for i in range(n)}

    assignment = pd.DataFrame({
        "name": names,
        "family": [family_of[x] for x in names],
        region.lower(): [subject[x] for x in names],
        f"{region.lower()}_len": [len(subject[x]) for x in names],
    }).sort_values(["family", "name"]).reset_index(drop=True)

    rows = []
    for fam, group in assignment.groupby("family"):
        members = list(group["name"])
        idx = [names.index(m) for m in members]
        sub_i = ident[np.ix_(idx, idx)]
        sub_s = simil[np.ix_(idx, idx)]
        loops = [subject[m] for m in members]
        # representative: the most frequent CDR3 of the family, ties by name
        counts = pd.Series(loops).value_counts()
        top_loop = counts.index[0]
        rep = sorted(m for m in members if subject[m] == top_loop)[0]
        rows.append({
            "family": fam,
            "n_sequences": len(members),
            "n_unique_cdr3": int(pd.Series(loops).nunique()),
            "representative": rep,
            "cdr3_len": f"{min(map(len, loops))}-{max(map(len, loops))}"
                        if len(set(map(len, loops))) > 1 else str(len(loops[0])),
            "min_identity": round(float(sub_i.min()), 1),
            "min_similarity": round(float(sub_s.min()), 1),
            "members": ", ".join(members),
            "cdr3_list": " | ".join(sorted(set(loops))),
        })
    families = (pd.DataFrame(rows)
                .sort_values("n_sequences", ascending=False)
                .reset_index(drop=True))

    return {"assignment": assignment, "families": families,
            "identity": ident_df, "similarity": simil_df, "linkage": link,
            "skipped": skipped, "cutoff": cutoff,
            "linkage_method": linkage_method, "region": region,
            "n_families": int(families["family"].nunique())}


def plot_family_dendrogram(result: dict, title: str | None = None):
    """Dendrogram of the family clustering, with the cutoff drawn in."""
    from scipy.cluster.hierarchy import dendrogram

    link = result.get("linkage")
    names = list(result["identity"].index)
    if link is None:
        fig, ax = plt.subplots(figsize=(5, 1.6))
        ax.text(0.5, 0.5, "A dendrogram needs at least two annotated CDR3s",
                ha="center", va="center", fontsize=9)
        ax.axis("off")
        return fig

    cutoff = result["cutoff"]
    height = max(2.4, 0.32 * len(names) + 1.2)
    fig, ax = plt.subplots(figsize=(8.2, height))
    dd = dendrogram(link, labels=names, orientation="right",
                    color_threshold=1.0 - cutoff / 100.0, ax=ax,
                    above_threshold_color="0.6")
    ax.axvline(1.0 - cutoff / 100.0, color="crimson", ls="--", lw=1.2)
    ax.annotate(f"{cutoff:.0f}% identity", xy=(1.0 - cutoff / 100.0, 1.0),
                xycoords=("data", "axes fraction"), xytext=(4, -12),
                textcoords="offset points", color="crimson", fontsize=8,
                va="top")
    ticks = [t / 100 for t in range(0, 101, 20)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{100 - t * 100:.0f}" for t in ticks])
    ax.set_xlabel(f"{result['region']} identity (%)")
    ax.set_title(title or f"Clonal families — {result['region']}, "
                          f"{result['linkage_method']} linkage", fontsize=10)
    ax.tick_params(axis="y", labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig, dd



# ---------------------------------------------------------------------------
# IMGT numbering, exposed for the modules that need positions rather than regions
# ---------------------------------------------------------------------------

def imgt_numbering(seq: str) -> list[tuple[int, str, str, int]]:
    """``[(imgt_number, insertion_code, residue, residue_index_1based), ...]``.

    ``annotate_anarci`` throws the numbering away once it has the regions, but
    humanness and disulfide analysis both need to talk about IMGT positions, so
    this returns it. Needs ANARCI.
    """
    from anarci import anarci

    seq = clean_sequence(seq)
    numbered, _details, _ = anarci([("query", seq)], scheme="imgt", output=False)
    if numbered[0] is None:
        raise ValueError("ANARCI could not number this sequence")
    numbering, start, _end = numbered[0][0]
    out, idx = [], start
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        out.append((int(num), (ins or " ").strip(), aa, idx + 1))
        idx += 1
    return out


def imgt_region(num: int) -> str:
    """Region of an IMGT position number."""
    for region, (lo, hi) in IMGT_CDR.items():
        if lo <= num <= hi:
            return region
    return "FR1" if num < 27 else "FR2" if num < 56 else "FR3" if num < 105 else "FR4"


# ---------------------------------------------------------------------------
# Aggregation-prone regions
# ---------------------------------------------------------------------------
# What this is, and what it is not.
#
# The published sequence-only solubility predictors (CamSol intrinsic,
# AGGRESCAN, TANGO) use coefficients fitted to experimental data. This module
# does NOT reproduce any of them: the numbers here are a transparent composite
# of four quantities that are each individually well established as determinants
# of aggregation propensity, combined with weights that are a judgement call.
#
# So read the profile for WHERE along the sequence the aggregation-prone stretches
# are, and compare sequences against each other. Do not read the absolute value
# as a predicted solubility — for a validated number, run the sequence through
# the CamSol or AGGRESCAN server.
#
# The four components, per residue:
#   charge       |net charge| at the working pH. Charged residues of either sign
#                oppose aggregation, so the absolute value is what counts.
#   hydrophobic  Kyte & Doolittle hydropathy. Exposed hydrophobic stretches are
#                the classic driver of aggregation.
#   beta         Chou-Fasman beta-sheet propensity. An aggregation-prone region is
#                typically a short beta-prone segment.
#   gatekeeper   P, D, E, K, R. Residues that flank and interrupt aggregation-prone
#                segments in natural sequences.

SOLUBILITY_WEIGHTS = {"charge": 1.0, "hydrophobic": 0.6, "beta": 0.8,
                      "gatekeeper": 0.5}
GATEKEEPERS = set("PDEKR")
APR_THRESHOLD = -0.3      # smoothed score below this marks an aggregation-prone window


def aggregation_profile(seq: str, ph: float | None = None,
                        window: int | None = None,
                        weights: dict | None = None) -> pd.DataFrame:
    """Per-residue aggregation-propensity profile, smoothed over a window.

    Returns one row per residue with each component and the combined ``score``:
    **positive leans soluble, negative leans aggregation-prone**. ``raw`` is the
    unsmoothed value, ``score`` the sliding-window mean.
    """
    import numpy as np

    seq = clean_sequence(seq)
    window = DEFAULTS["profile_window"] if window is None else window
    w = {**SOLUBILITY_WEIGHTS, **(weights or {})}
    del ph   # formal charges are used, not pH-dependent ones; see CHARGE above

    # CHARGE is the formal charge at neutral pH used by the profile plots; the
    # absolute value is what matters here, since either sign opposes aggregation.
    charge = np.array([abs(CHARGE.get(a, 0)) for a in seq], float)
    hydro = np.array([KD[a] for a in seq], float) / 4.5          # to roughly [-1, 1]
    # CHOU_FASMAN stores propensities x100 (P-beta 1.09 is held as 109), so
    # rescale before centring on 1.0 = no preference.
    beta = np.array([CHOU_FASMAN[a][1] for a in seq], float) / 100.0 - 1.0
    gate = np.array([1.0 if a in GATEKEEPERS else 0.0 for a in seq], float)

    raw = (w["charge"] * charge - w["hydrophobic"] * hydro
           - w["beta"] * beta + w["gatekeeper"] * gate)
    smooth = sliding(raw, window)
    return pd.DataFrame({
        "position": range(1, len(seq) + 1), "residue": list(seq),
        "charge": charge.round(3), "hydropathy": hydro.round(3),
        "beta_propensity": beta.round(3), "gatekeeper": gate,
        "raw": raw.round(3), "score": np.round(smooth, 3),
    })


def aggregation_regions(profile: pd.DataFrame, annotation: dict | None = None,
                        threshold: float = APR_THRESHOLD,
                        min_length: int = 3) -> pd.DataFrame:
    """Contiguous stretches whose smoothed score stays below ``threshold``."""
    rows, run = [], []
    for _, r in profile.iterrows():
        if r["score"] < threshold:
            run.append(r)
        elif run:
            rows.append(run)
            run = []
    if run:
        rows.append(run)

    out = []
    for run in rows:
        if len(run) < min_length:
            continue
        start, end = int(run[0]["position"]), int(run[-1]["position"])
        regions = sorted({region_of(int(r["position"]), annotation) for r in run})
        out.append({
            "start": start, "end": end, "length": end - start + 1,
            "range": f"{start}-{end}",
            "sequence": "".join(r["residue"] for r in run),
            "region": "/".join(regions),
            "in_CDR": any(x.startswith("CDR") for x in regions),
            "min_score": round(min(float(r["score"]) for r in run), 3),
            "mean_score": round(sum(float(r["score"]) for r in run) / len(run), 3),
        })
    return pd.DataFrame(out)


def aggregation_summary(profile: pd.DataFrame, regions: pd.DataFrame) -> dict:
    """Headline numbers for the profile."""
    worst = profile.loc[profile["score"].idxmin()] if len(profile) else None
    return {
        "mean_score": round(float(profile["score"].mean()), 3) if len(profile) else None,
        "min_score": round(float(profile["score"].min()), 3) if len(profile) else None,
        "worst_position": int(worst["position"]) if worst is not None else None,
        "n_regions": int(len(regions)),
        "residues_in_regions": int(regions["length"].sum()) if len(regions) else 0,
        "fraction_prone": round(float(regions["length"].sum()) / len(profile), 3)
                          if len(profile) and len(regions) else 0.0,
    }


def weakest_windows(profile: pd.DataFrame, annotation: dict | None = None,
                    n: int = 3, width: int = 7) -> pd.DataFrame:
    """The ``n`` lowest-scoring windows, whether or not they cross the threshold.

    A well-behaved VHH often has no region below the threshold at all, which
    would leave the table empty and tell you nothing. These are its weakest
    spots regardless, which is what you actually want to compare between
    candidates.
    """
    if not len(profile):
        return pd.DataFrame()
    scores = profile["score"].to_numpy()
    seq = "".join(profile["residue"])
    taken, rows = set(), []
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    for centre in order:
        if len(rows) >= n:
            break
        if any(abs(centre - t) < width for t in taken):
            continue          # keep the windows from overlapping each other
        taken.add(centre)
        lo = max(0, centre - width // 2)
        hi = min(len(seq), lo + width)
        regions = sorted({region_of(i + 1, annotation) for i in range(lo, hi)})
        rows.append({
            "rank": len(rows) + 1,
            "centre": centre + 1,
            "range": f"{lo + 1}-{hi}",
            "sequence": seq[lo:hi],
            "region": "/".join(regions),
            "in_CDR": any(x.startswith("CDR") for x in regions),
            "score": round(float(scores[centre]), 3),
            "below_threshold": bool(scores[centre] < APR_THRESHOLD),
        })
    return pd.DataFrame(rows)


def plot_aggregation(profile: pd.DataFrame, annotation: dict | None = None,
                     threshold: float = APR_THRESHOLD,
                     title: str = "Aggregation-prone regions"):
    """Profile with aggregation-prone stretches shaded and CDR bands marked."""
    import numpy as np

    fig, ax = plt.subplots(figsize=(10, 3.4))
    x = profile["position"].to_numpy()
    y = profile["score"].to_numpy()

    if annotation:
        for region, colour in REGION_COLORS.items():
            info = annotation.get(region) or {}
            if info.get("start") and info.get("end") and region.startswith("CDR"):
                ax.axvspan(info["start"] - 0.5, info["end"] + 0.5, color=colour,
                           alpha=0.22, lw=0, zorder=0)
                ax.text((info["start"] + info["end"]) / 2, 1.02, region,
                        transform=ax.get_xaxis_transform(), ha="center",
                        fontsize=7.5, color="0.3")

    ax.axhline(0, color="0.5", lw=0.8)
    ax.axhline(threshold, color="crimson", lw=0.9, ls="--")
    ax.plot(x, y, color="#1f3b73", lw=1.5, zorder=3)
    ax.fill_between(x, y, threshold, where=(y < threshold), color="crimson",
                    alpha=0.3, zorder=2, interpolate=True)
    ax.set_xlabel("Residue")
    ax.set_ylabel("← aggregation-prone    soluble →", fontsize=8.5)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Humanness
# ---------------------------------------------------------------------------
# Percent identity to the closest human germline V gene, computed position by
# position in IMGT numbering rather than by alignment, using the IMGT-gapped
# germline sequences that ship with ANARCI.
#
# This is the usual first-pass humanness indicator and nothing more. It is NOT
# an immunogenicity prediction: immunogenicity depends on T-cell epitopes,
# formulation, dose and route, none of which is in the sequence alone.

_HUMAN_V: dict | None = None


def human_germlines(family: str | None = None) -> dict:
    """IMGT-gapped human heavy-chain V germlines from ANARCI's own data."""
    global _HUMAN_V
    if _HUMAN_V is None:
        from anarci.germlines import all_germlines
        _HUMAN_V = dict(all_germlines["V"]["H"]["human"])
    if family:
        return {k: v for k, v in _HUMAN_V.items() if k.startswith(family)}
    return dict(_HUMAN_V)


def humanness(seq: str, family: str | None = "IGHV3",
              annotation: dict | None = None) -> dict:
    """Identity to the closest human germline, with the divergent positions.

    Identity is computed over **framework positions only**, because CDR
    divergence from a human germline is expected and not something you would
    ever "fix". ``family`` restricts the search ("IGHV3" is the family VHHs
    belong to); pass ``None`` to search all human heavy V genes.

    The returned ``divergences`` table is the actionable part: one row per
    framework position where the VHH differs from the closest human germline,
    flagged for whether it is one of the four hallmark positions and whether the
    substitution is conservative.
    """
    numbering = imgt_numbering(seq)
    by_imgt = {num: aa for num, ins, aa, _i in numbering if not ins}
    idx_of = {num: i for num, ins, _aa, i in numbering if not ins}
    hallmark_imgt = {h[0] for h in HALLMARK_POSITIONS}

    best, best_id = None, -1.0
    scores = {}
    for gene, gl in human_germlines(family).items():
        match = total = 0
        for num, aa in by_imgt.items():
            if num > len(gl):
                continue
            g = gl[num - 1]
            if g == "-" or imgt_region(num).startswith("CDR"):
                continue
            total += 1
            match += (g == aa)
        if total:
            pct = 100 * match / total
            scores[gene] = round(pct, 1)
            if pct > best_id:
                best, best_id = gene, pct

    if best is None:
        return {"gene": None, "identity": None, "divergences": pd.DataFrame(),
                "by_region": pd.DataFrame(), "all_scores": scores}

    gl = human_germlines(family)[best]
    rows, per_region = [], {}
    for num, aa in sorted(by_imgt.items()):
        if num > len(gl):
            continue
        g = gl[num - 1]
        if g == "-":
            continue
        region = imgt_region(num)
        stat = per_region.setdefault(region, {"match": 0, "total": 0})
        stat["total"] += 1
        stat["match"] += (g == aa)
        if g != aa and not region.startswith("CDR"):
            rows.append({
                "imgt": num, "region": region,
                "position": idx_of.get(num),
                "vhh": aa, "human_germline": g,
                "hallmark": num in hallmark_imgt,
                "conservative": _blosum(aa, g) > 0 if _blosum_ready() else None,
            })

    divergences = pd.DataFrame(rows)
    by_region = pd.DataFrame([
        {"region": r, "identity": round(100 * v["match"] / v["total"], 1),
         "positions": v["total"]}
        for r, v in sorted(per_region.items())
    ])
    return {"gene": best, "identity": round(best_id, 1),
            "divergences": divergences, "by_region": by_region,
            "all_scores": dict(sorted(scores.items(), key=lambda kv: -kv[1])),
            "n_divergent": int(len(divergences)),
            "n_hallmark": int(divergences["hallmark"].sum()) if len(divergences) else 0}


def _blosum_ready() -> bool:
    """BLOSUM62 is loaded lazily by the family aligner; make sure it is there."""
    _family_aligner()
    return _BLOSUM62 is not None


# ---------------------------------------------------------------------------
# Disulfides
# ---------------------------------------------------------------------------
# The canonical VHH disulfide is IMGT 23-104. Camelid VHHs with long CDR3s
# frequently carry a second, non-canonical pair, most often between FR2
# (IMGT 50) or CDR1 and a cysteine in CDR3. That pair is a stabiliser, not a
# liability, which is the opposite of how a lone extra cysteine should be read -
# hence this analysis, which the liabilities table on its own cannot give you.

CANONICAL_PAIR = (23, 104)
EXTRA_PAIR_REGIONS = {("CDR1", "CDR3"), ("FR2", "CDR3"), ("CDR2", "CDR3"),
                      ("CDR1", "CDR2")}


def disulfides(seq: str, annotation: dict | None = None) -> dict:
    """Cysteine inventory: the canonical pair, plausible extra pairs, orphans.

    Pairing is inferred from position alone - there is no structure here, so a
    "plausible pair" means two cysteines sitting in regions known to pair in
    real VHHs, not a verified bond.
    """
    seq = clean_sequence(seq)
    try:
        numbering = imgt_numbering(seq)
        imgt_of = {i: num for num, ins, _aa, i in numbering if not ins}
    except Exception:
        imgt_of = {}

    cys = []
    for i, aa in enumerate(seq, start=1):
        if aa != "C":
            continue
        num = imgt_of.get(i)
        cys.append({"position": i, "imgt": num,
                    "region": imgt_region(num) if num else region_of(i, annotation),
                    "canonical": num in CANONICAL_PAIR if num else False})

    canonical = [c for c in cys if c["canonical"]]
    extra = [c for c in cys if not c["canonical"]]

    pairs, orphans = [], []
    if len(extra) == 2:
        a, b = extra
        key = tuple(sorted((a["region"], b["region"])))
        plausible = key in {tuple(sorted(k)) for k in EXTRA_PAIR_REGIONS}
        pairs.append({
            "cys_1": a["position"], "cys_2": b["position"],
            "imgt_1": a["imgt"], "imgt_2": b["imgt"],
            "regions": f"{a['region']} / {b['region']}",
            "plausible": plausible,
            "note": "regions that pair in real VHHs - likely an additional, "
                    "stabilising disulfide" if plausible else
                    "an unusual region combination - treat as unverified",
        })
    elif len(extra) == 1:
        orphans = extra
    elif len(extra) > 2:
        orphans = extra

    if len(canonical) == 2:
        verdict = "canonical disulfide present"
    elif len(canonical) == 1:
        verdict = f"canonical disulfide incomplete - only IMGT {canonical[0]['imgt']} found"
    else:
        verdict = "canonical IMGT 23-104 disulfide not found"
    if pairs and pairs[0]["plausible"]:
        verdict += "; plus a plausible additional pair"
    elif len(extra) == 1:
        verdict += "; one extra cysteine with no partner - likely a free thiol"
    elif len(extra) > 2:
        verdict += f"; {len(extra)} extra cysteines - pairing cannot be inferred from sequence"

    return {"cysteines": pd.DataFrame(cys), "n_cys": len(cys),
            "canonical": pd.DataFrame(canonical), "extra": pd.DataFrame(extra),
            "pairs": pd.DataFrame(pairs), "orphans": pd.DataFrame(orphans),
            "verdict": verdict,
            "has_extra_pair": bool(pairs and pairs[0]["plausible"]),
            "n_orphan": len(orphans)}


# ---------------------------------------------------------------------------
# Suggested repairs for liabilities
# ---------------------------------------------------------------------------
# One conservative substitution per liability, so the table stops being merely
# descriptive. Two rules of thumb are built in: a substitution inside a CDR is
# flagged because it risks affinity, and a motif that has no single clean fix
# says so instead of offering one.

# liability label -> list of (offset within the motif, from, to, rationale)
LIABILITY_REPAIRS = {
    "N-glycosylation": [(2, "ST", "A", "breaks the N-x-S/T sequon at the +2 residue"),
                        (0, "N", "Q", "removes the acceptor asparagine")],
    "Deamidation (high)": [(1, "G", "A", "the +1 glycine is what makes NG fast"),
                           (0, "N", "Q", "removes the deamidation site")],
    "Deamidation (moderate)": [(0, "N", "Q", "removes the deamidation site")],
    "Isomerisation (high)": [(1, "G", "A", "the +1 glycine is what makes DG fast"),
                             (0, "D", "E", "aspartate to glutamate blocks succinimide")],
    "Isomerisation (moderate)": [(0, "D", "E", "aspartate to glutamate")],
    "Oxidation (Met)": [(0, "M", "L", "isosteric, not oxidisable")],
    "Oxidation (Trp)": [(0, "W", "F", "keeps the aromatic, loses the indole")],
    "Fragmentation": [(0, "D", "E", "removes the acid-labile Asp-Pro bond")],
    "Integrin binding": [(2, "D", "E", "breaks the RGD motif")],
    "Free cysteine": [(0, "C", "S", "isosteric, no free thiol"),
                      (0, "C", "A", "if serine is not tolerated")],
}

# Structurally conserved IMGT positions. Offering a substitution here would
# "repair" a liability by breaking the fold: Cys23/Cys104 form the canonical
# disulfide, Trp41 packs the FR2 core, and Trp118 opens the FR4 WGxG motif. The
# oxidation rules match Met and Trp anywhere, so without this guard the table
# cheerfully suggests W41F and W118F.
CONSERVED_IMGT = {
    23: "canonical disulfide cysteine",
    41: "conserved FR2 tryptophan, packs the domain core",
    104: "canonical disulfide cysteine",
    118: "conserved FR4 tryptophan of the WGxG motif",
}

NO_CLEAN_FIX = {
    "Hydrophobic patch": "a hydrophobic stretch has no single-residue fix - it needs "
                         "redesign of the segment, or acceptance with a formulation "
                         "strategy",
    "Missing cysteine": "restoring the canonical disulfide means putting a cysteine "
                        "back at IMGT 23 or 104, which is a decision about the "
                        "scaffold rather than a substitution",
}


def suggest_repairs(seq: str, liabilities: pd.DataFrame,
                    annotation: dict | None = None,
                    include_cdr: bool = True) -> pd.DataFrame:
    """One row per suggested substitution, derived from a liabilities table."""
    seq = clean_sequence(seq)
    imgt_of: dict[int, int] = {}
    try:
        imgt_of = {i: num for num, ins, _aa, i in imgt_numbering(seq) if not ins}
    except Exception:
        pass        # no ANARCI: the conserved-position guard is simply unavailable

    rows = []
    for _, lb in liabilities.iterrows():
        label = lb["liability"]
        start = int(lb.get("position") or 0)
        if label in NO_CLEAN_FIX:
            rows.append({"liability": label, "motif": lb.get("motif", ""),
                         "position": start or None,
                         "region": lb.get("region", "?"),
                         "in_CDR": bool(lb.get("in_CDR", False)),
                         "substitution": "-", "repaired_motif": "-",
                         "rationale": NO_CLEAN_FIX[label], "applicable": False})
            continue
        for offset, frm, to, why in LIABILITY_REPAIRS.get(label, []):
            pos = start + offset
            if not (1 <= pos <= len(seq)) or seq[pos - 1] not in frm:
                continue
            region = region_of(pos, annotation)
            in_cdr = region.startswith("CDR")
            if in_cdr and not include_cdr:
                continue
            conserved = CONSERVED_IMGT.get(imgt_of.get(pos, -1))
            if conserved:
                rows.append({
                    "liability": label, "motif": str(lb.get("motif", "") or ""),
                    "position": pos, "region": region, "in_CDR": in_cdr,
                    "substitution": "-", "repaired_motif": "-",
                    "rationale": f"no substitution offered: IMGT "
                                 f"{imgt_of[pos]} is the {conserved}",
                    "applicable": False,
                })
                continue
            motif = str(lb.get("motif", "") or "")
            repaired = (motif[:offset] + to + motif[offset + 1:]) if motif else to
            rows.append({
                "liability": label, "motif": motif, "position": pos,
                "region": region, "in_CDR": in_cdr,
                "substitution": f"{seq[pos - 1]}{pos}{to}",
                "repaired_motif": repaired,
                "rationale": why + (" - inside a CDR, so this may cost affinity"
                                    if in_cdr else ""),
                "applicable": True,
            })
    out = pd.DataFrame(rows)
    if len(out):
        # position must be part of the key: every refused row carries
        # substitution "-", so without it one conserved-position warning
        # silently swallows the others.
        out = out.drop_duplicates(subset=["substitution", "liability", "position"])
        out = out.sort_values(["position", "liability"]).reset_index(drop=True)
    return out


def apply_repairs(seq: str, substitutions: list[str]) -> tuple[str, list[str]]:
    """Apply substitutions written as ``N55Q``; returns ``(sequence, problems)``.

    A substitution whose original residue does not match is refused rather than
    applied silently, and two substitutions at the same position are refused as
    a pair - that conflict is a decision for you, not for this function.
    """
    import re as _re

    seq = clean_sequence(seq)
    chars, problems, seen = list(seq), [], {}
    for sub in substitutions:
        m = _re.fullmatch(r"([A-Z])(\d+)([A-Z])", sub.strip().upper())
        if not m:
            problems.append(f"{sub}: not of the form N55Q")
            continue
        frm, pos, to = m.group(1), int(m.group(2)), m.group(3)
        if not 1 <= pos <= len(chars):
            problems.append(f"{sub}: position outside the sequence")
            continue
        if seq[pos - 1] != frm:
            problems.append(f"{sub}: residue {pos} is {seq[pos - 1]}, not {frm}")
            continue
        if pos in seen and seen[pos] != to:
            problems.append(f"{sub}: conflicts with {seen[pos]} already set at {pos}")
            continue
        chars[pos - 1] = to
        seen[pos] = to
    return "".join(chars), problems



def build_pdf_report(path: str | Path, df: pd.DataFrame, seqs: dict[str, str],
                     annotations: dict, figures: dict | None = None,
                     ph: float | None = None, version: str = "",
                     author: str = "", with_hallmarks: bool = True,
                     with_liabilities: bool = True,
                     app_name: str = "VHH tool box") -> Path:
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
        title="VHH analysis report", author=author or app_name)

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
        meta += f" · {app_name} {version}"
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
                summary = liability_summary(name, liab)
                story.append(Paragraph(
                    f"Total score {summary['total_score']} over "
                    f"{summary['n_liabilities']} motif(s): {summary['n_high']} high, "
                    f"{summary['n_medium']} medium, {summary['n_low']} low; "
                    f"{summary['n_in_CDR']} inside a CDR.", body))
                lrows = [["Liability", "Motif", "Position", "Region", "Score", "Risk"]]
                lrows += liab[["liability", "motif", "range", "region", "score",
                               "risk"]].astype(str).values.tolist()
                t = table(lrows, widths=[42 * mm, 18 * mm, 22 * mm, 20 * mm,
                                         16 * mm, 18 * mm])
                extra = []
                for i, level in enumerate(liab["risk"], start=1):
                    extra.append(("TEXTCOLOR", (5, i), (5, i),
                                  colors.HexColor(RISK_COLORS.get(level, "#333333"))))
                t.setStyle(TableStyle(extra))
                story.append(t)
                story.append(Paragraph(
                    "A motif inside a CDR may affect binding; the same motif in a "
                    "framework is mainly a stability or heterogeneity concern. "
                    "Scores follow the in-house scale without the solvent-exposure "
                    "weighting a structure would provide, so the categories compress "
                    "towards the top and the ranking between sequences is more "
                    "informative than any single category.", small))

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

    footer_bits = [b for b in [app_name, version, author] if b]
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
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "SapI": "GCTCTTC", "PstI": "CTGCAG",
    "XbaI": "TCTAGA", "SpeI": "ACTAGT", "KpnI": "GGTACC", "SacI": "GAGCTC",
    # BstEII is deliberately offered: its site sits in the VHH FR4 and is used
    # on purpose in the classic nanobody cloning strategy, so removing it is a
    # choice, not a default.
    "BstEII": "GGTNACC",
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
