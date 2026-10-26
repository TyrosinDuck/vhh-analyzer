"""
codons — diagnose an existing coding sequence that expresses poorly.

The reverse of the DNA tab: instead of designing a gene, you paste one that is
already in a plasmid and ask what might be wrong with it.

**Read this before reading the output.** Poor VHH expression in *E. coli* has
several more likely causes than codon usage, and they are worth eliminating
first: periplasmic toxicity, inefficient signal-peptide cleavage, failure to form
the intradomain disulfide in the periplasm, inclusion bodies, plasmid loss,
induction conditions. A gene with perfect codon usage expresses badly for all of
those reasons. This module tells you whether the *sequence* is part of the
problem — it cannot tell you that it is the problem.

**There is no CAI here.** A Codon Adaptation Index needs a full codon usage table
for the host, and this module does not ship one; inventing sixty-four relative
adaptiveness values from memory would produce a number that looks precise and
means nothing. What it does instead is report the handful of codons whose scarcity
in *E. coli* is well established and routinely supplemented by the Rosetta and
CodonPlus strains, where they sit, and whether they cluster — tandem rare codons
being the documented problem rather than isolated ones.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import vhh_core as vc

__all__ = [
    "RARE_CODONS", "orf_report", "codon_table", "rare_clusters", "gc_windows",
    "gc_regions",
    "diagnose", "plot_diagnosis", "__version__",
]

__version__ = "2026.10.26a"

# Codons whose scarcity in E. coli is well established — these are the ones the
# Rosetta and BL21-CodonPlus strains supplement with extra tRNA genes. The gene
# name is given only where it is unambiguous.
RARE_CODONS: dict[str, dict] = {
    "AGA": {"aa": "R", "trna": "argU", "note": "rare arginine, the classic one"},
    "AGG": {"aa": "R", "trna": "argU/argW", "note": "rare arginine"},
    "CGA": {"aa": "R", "trna": None, "note": "uncommon arginine"},
    "CGG": {"aa": "R", "trna": None, "note": "uncommon arginine"},
    "ATA": {"aa": "I", "trna": "ileY", "note": "rare isoleucine"},
    "CTA": {"aa": "L", "trna": "leuW", "note": "rare leucine"},
    "CCC": {"aa": "P", "trna": "proL", "note": "rare proline"},
    "GGA": {"aa": "G", "trna": "glyT", "note": "uncommon glycine"},
}

CLUSTER_WINDOW = 6        # codons; two rare codons this close count as a cluster
GC_WINDOW = 30            # nucleotides
GC_LOW, GC_HIGH = 0.30, 0.70   # practical synthesis band, not a biological one
HEAD_CODONS = 15          # the translation initiation region, in codons
STOPS = {"TAA", "TAG", "TGA"}


def _clean_dna(dna: str) -> str:
    dna = re.sub(r"[^A-Za-z]", "", dna or "").upper().replace("U", "T")
    bad = set(dna) - set("ACGTN")
    if bad:
        raise ValueError(f"not a DNA sequence — unexpected characters: "
                         f"{', '.join(sorted(bad))}")
    if not dna:
        raise ValueError("empty sequence")
    return dna


def orf_report(dna: str) -> dict:
    """Frame, start, stop and internal stop codons — the deterministic checks.

    These come first because they are the ones that are definitely wrong when
    they are wrong, as opposed to the heuristics below.
    """
    dna = _clean_dna(dna)
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    internal = [(i + 1, c) for i, c in enumerate(codons[:-1]) if c in STOPS]
    protein = vc.translate_dna(dna)
    return {
        "length_nt": len(dna),
        "n_codons": len(codons),
        "in_frame": len(dna) % 3 == 0,
        "trailing_nt": len(dna) % 3,
        "starts_with_ATG": bool(codons) and codons[0] == "ATG",
        "ends_with_stop": bool(codons) and codons[-1] in STOPS,
        "internal_stops": internal,
        "n_internal_stops": len(internal),
        "protein": protein.rstrip("*"),
        "protein_length": len(protein.rstrip("*")),
        "has_N": "N" in dna,
    }


def codon_table(dna: str) -> pd.DataFrame:
    """One row per codon, with its amino acid and whether it is a rare one."""
    dna = _clean_dna(dna)
    rows = []
    for i in range(0, len(dna) - len(dna) % 3, 3):
        codon = dna[i:i + 3]
        info = RARE_CODONS.get(codon)
        rows.append({
            "codon_index": i // 3 + 1,
            "nt_start": i + 1,
            "codon": codon,
            "aa": vc.translate_dna(codon) or "?",
            "rare": info is not None,
            "trna": (info or {}).get("trna") or "",
            "note": (info or {}).get("note") or "",
        })
    return pd.DataFrame(rows)


def rare_clusters(table: pd.DataFrame, window: int = CLUSTER_WINDOW) -> pd.DataFrame:
    """Groups of rare codons sitting within ``window`` codons of each other.

    Isolated rare codons are tolerated well; it is tandem and clustered ones that
    are associated with ribosome stalling, so this is the part of the output
    worth acting on.
    """
    rare = table[table["rare"]].reset_index(drop=True)
    if len(rare) < 2:
        return pd.DataFrame()

    groups, current = [], [rare.iloc[0]]
    for _, row in rare.iloc[1:].iterrows():
        if row["codon_index"] - current[-1]["codon_index"] <= window:
            current.append(row)
        else:
            groups.append(current)
            current = [row]
    groups.append(current)

    out = []
    for group in groups:
        if len(group) < 2:
            continue
        out.append({
            "start_codon": int(group[0]["codon_index"]),
            "end_codon": int(group[-1]["codon_index"]),
            "n_rare": len(group),
            "span": int(group[-1]["codon_index"] - group[0]["codon_index"] + 1),
            "codons": " ".join(g["codon"] for g in group),
            "tandem": any(group[i + 1]["codon_index"] - group[i]["codon_index"] == 1
                          for i in range(len(group) - 1)),
        })
    return pd.DataFrame(out)


def gc_windows(dna: str, window: int = GC_WINDOW) -> pd.DataFrame:
    """GC fraction in a sliding window, with the extremes flagged."""
    dna = _clean_dna(dna)
    if len(dna) < window:
        window = max(3, len(dna))
    gc = np.array([1.0 if b in "GC" else 0.0 for b in dna])
    kernel = np.ones(window) / window
    smooth = np.convolve(gc, kernel, mode="valid")
    starts = np.arange(1, len(smooth) + 1)
    return pd.DataFrame({
        "start": starts, "end": starts + window - 1,
        "gc": np.round(smooth, 3),
        "extreme": (smooth < GC_LOW) | (smooth > GC_HIGH),
    })


def gc_regions(gc: pd.DataFrame, window: int = GC_WINDOW) -> pd.DataFrame:
    """Contiguous stretches where the GC window stays outside the band.

    Counting *windows* overstates the problem badly: they overlap by a single
    nucleotide, so one awkward stretch of 60 nt shows up as thirty "extreme
    windows". Merging them into regions, and keeping only those that persist for
    at least a window length, reports the handful of places that actually matter.
    """
    if not len(gc):
        return pd.DataFrame()
    flags = gc["extreme"].to_numpy()
    rows, start = [], None
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            rows.append((start, i - 1))
            start = None
    if start is not None:
        rows.append((start, len(flags) - 1))

    spans = []
    for lo, hi in rows:
        span_start = int(gc["start"].iloc[lo])
        span_end = int(gc["end"].iloc[hi])
        if span_end - span_start + 1 < window:
            continue            # a brief excursion, not a stretch
        block = gc["gc"].iloc[lo:hi + 1]
        spans.append([span_start, span_end, float(block.min()),
                      float(block.max()), float(block.mean())])

    # Each window covers `window` nucleotides, so two groups of flagged windows
    # that are separate in window-space can still cover overlapping nucleotides.
    # Merge them, otherwise the table reports the same stretch twice.
    merged = []
    for span in sorted(spans):
        if merged and span[0] <= merged[-1][1]:
            last = merged[-1]
            last[1] = max(last[1], span[1])
            last[2] = min(last[2], span[2])
            last[3] = max(last[3], span[3])
            last[4] = (last[4] + span[4]) / 2
        else:
            merged.append(list(span))

    return pd.DataFrame([{
        "start": int(a), "end": int(b), "length": int(b - a + 1),
        "min_gc": round(lo, 3), "max_gc": round(hi, 3),
        "direction": "GC-rich" if mean > 0.5 else "AT-rich",
    } for a, b, lo, hi, mean in merged])


def diagnose(dna: str, expected_protein: str | None = None,
             cluster_window: int = CLUSTER_WINDOW,
             gc_window: int = GC_WINDOW) -> dict:
    """Everything above, plus a short list of findings in plain words."""
    dna = _clean_dna(dna)
    orf = orf_report(dna)
    table = codon_table(dna)
    clusters = rare_clusters(table, cluster_window)
    gc = gc_windows(dna, gc_window)
    regions = gc_regions(gc, gc_window)

    head = table.head(HEAD_CODONS)
    head_nt = dna[:HEAD_CODONS * 3]
    head_gc = (sum(1 for b in head_nt if b in "GC") / len(head_nt)
               if head_nt else float("nan"))
    homopolymers = [(m.start() + 1, m.group()) for m in
                    re.finditer(r"(A{7,}|C{7,}|G{7,}|T{7,})", dna)]

    findings = []
    if not orf["in_frame"]:
        findings.append(
            ("error", f"The sequence is not a whole number of codons — "
                      f"{orf['trailing_nt']} nucleotide(s) left over. Everything "
                      "below assumes the first base starts a codon, so check the "
                      "boundaries before reading anything else."))
    if orf["n_internal_stops"]:
        first = orf["internal_stops"][0]
        findings.append(
            ("error", f"{orf['n_internal_stops']} internal stop codon(s), the "
                      f"first at codon {first[0]} ({first[1]}). In the frame as "
                      "given this gene cannot produce a full-length protein."))
    if orf["has_N"]:
        findings.append(("warning", "The sequence contains N. Those positions "
                                    "are translated as X and are not analysed."))
    if not orf["starts_with_ATG"]:
        findings.append(("info", "No ATG at the first codon — expected if you "
                                 "pasted a fragment or a sequence starting inside "
                                 "a signal peptide."))
    if not orf["ends_with_stop"]:
        findings.append(("info", "No stop codon at the end — again expected for a "
                                 "fragment."))

    n_rare = int(table["rare"].sum())
    pct_rare = n_rare / max(1, len(table))
    if len(clusters):
        tandem = int(clusters["tandem"].sum())
        findings.append(
            ("warning", f"{len(clusters)} cluster(s) of rare codons"
                        + (f", {tandem} containing codons in tandem" if tandem else "")
                        + ". Clustered rare codons are the ones associated with "
                          "ribosome stalling; isolated ones are generally "
                          "tolerated."))
    if n_rare and not len(clusters):
        findings.append(
            ("info", f"{n_rare} rare codon(s) ({pct_rare:.1%}), none clustered. "
                     "Scattered rare codons are usually harmless."))
    if not n_rare:
        findings.append(("ok", "No rare E. coli codon at all."))

    if len(regions):
        worst = regions.loc[regions["length"].idxmax()]
        findings.append(
            ("warning", f"{len(regions)} stretch(es) where GC stays outside "
                        f"{GC_LOW:.0%}–{GC_HIGH:.0%} for at least {gc_window} nt, "
                        f"the longest {int(worst['length'])} nt at position "
                        f"{int(worst['start'])} ({worst['direction']}). Extreme GC "
                        "complicates synthesis and can affect mRNA structure."))
    else:
        findings.append(("ok", f"GC stays within {GC_LOW:.0%}–{GC_HIGH:.0%} "
                               "throughout, apart from brief excursions."))
    if head_gc == head_gc and head_gc > 0.60:
        findings.append(
            ("info", f"The first {HEAD_CODONS} codons are {head_gc:.0%} GC. A "
                     "GC-rich start can fold the ribosome binding region and "
                     "slow initiation — this is one of the few places where "
                     "silent changes reliably help."))
    if homopolymers:
        findings.append(
            ("warning", f"{len(homopolymers)} homopolymer run(s) of 7 or more: "
                        + ", ".join(f"{run} at {pos}" for pos, run in homopolymers[:4])
                        + ". These complicate synthesis and Sanger reading."))

    protein_match = None
    if expected_protein:
        expected = vc.clean_sequence(expected_protein)
        protein_match = orf["protein"].lstrip("M") == expected.lstrip("M") \
            or orf["protein"] == expected
        findings.append(
            ("ok" if protein_match else "error",
             "The translation matches the expected protein."
             if protein_match else
             "The translation does **not** match the expected protein — compare "
             "them before looking at anything else here."))

    return {"orf": orf, "codons": table, "clusters": clusters, "gc": gc,
            "gc_regions": regions,
            "head_gc": head_gc, "homopolymers": homopolymers,
            "n_rare": n_rare, "fraction_rare": round(pct_rare, 4),
            "gc_overall": round(float(np.mean([1 if b in "GC" else 0
                                               for b in dna])), 4),
            "findings": findings, "protein_match": protein_match}


def plot_diagnosis(result: dict, gc_window: int = GC_WINDOW):
    """Rare codons along the gene, over the GC profile."""
    import matplotlib.pyplot as plt

    table = result["codons"]
    gc = result["gc"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 4.2), sharex=True,
                             height_ratios=[1, 2])

    ax = axes[0]
    rare = table[table["rare"]]
    ax.vlines(rare["nt_start"], 0, 1, color="crimson", lw=1.4)
    for _, cl in (result["clusters"] if len(result["clusters"])
                  else pd.DataFrame()).iterrows() if len(result["clusters"]) else []:
        ax.axvspan((cl["start_codon"] - 1) * 3, cl["end_codon"] * 3,
                   color="crimson", alpha=0.18, lw=0)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_ylabel("Rare\ncodons", fontsize=8)
    ax.set_title(f"{result['n_rare']} rare codon(s), "
                 f"{len(result['clusters'])} cluster(s)", fontsize=10, loc="left")

    ax = axes[1]
    ax.plot(gc["start"], gc["gc"], color="#4c72b0", lw=1.3)
    ax.axhline(GC_LOW, color="0.6", ls="--", lw=0.8)
    ax.axhline(GC_HIGH, color="0.6", ls="--", lw=0.8)
    hot = gc[gc["extreme"]]
    if len(hot):
        ax.plot(hot["start"], hot["gc"], ".", color="crimson", ms=3)
    ax.set_ylim(0, 1)
    ax.set_ylabel(f"GC ({gc_window} nt)", fontsize=9)
    ax.set_xlabel("Position (nt)")
    for a in axes:
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    return fig
