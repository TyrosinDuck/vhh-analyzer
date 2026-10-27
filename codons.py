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
    "gc_regions", "suggest_gc_fix", "parse_dna_input", "summarise",
    "plot_gc_comparison",
    "diagnose", "plot_diagnosis", "__version__",
]

__version__ = "2026.10.27c"

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



# ---------------------------------------------------------------------------
# Lowering GC without touching the protein
# ---------------------------------------------------------------------------
# A GC-rich stretch is one of the few sequence problems with a clean fix: swap
# codons for lower-GC synonyms. The constraints are what make it non-trivial —
# the substitute must not be a rare codon, must not create a restriction site
# you are using, and must not open a homopolymer run. Every suggestion below has
# been checked against all three, and the protein is verified unchanged.
#
# What this cannot promise is that the gene will then express. Lowering local GC
# removes one plausible obstacle; see the warning at the top of this module for
# the others.

_SYNONYMS: dict[str, list[str]] | None = None


def _synonyms() -> dict[str, list[str]]:
    """Amino acid -> every codon, bacterial table."""
    global _SYNONYMS
    if _SYNONYMS is None:
        from Bio.Data import CodonTable
        table = CodonTable.unambiguous_dna_by_id[11]
        out: dict[str, list[str]] = {}
        for codon, aa in table.forward_table.items():
            out.setdefault(aa, []).append(codon)
        _SYNONYMS = out
    return _SYNONYMS


def _gc(seq: str) -> float:
    return sum(1 for b in seq if b in "GC") / len(seq) if seq else 0.0


def suggest_gc_fix(dna: str, regions: pd.DataFrame | None = None,
                   avoid_sites: dict | None = None,
                   gc_window: int = GC_WINDOW,
                   max_homopolymer: int = 6,
                   allow_rare: bool = False,
                   strategy: str = "minimal") -> dict:
    """Synonymous swaps that lower GC inside the flagged stretches.

    Only codons inside a GC-rich region are considered — there is no reason to
    disturb the rest of the gene. For each, the lowest-GC synonym is taken that
    is not a documented rare codon, does not create one of ``avoid_sites``, and
    does not produce a homopolymer run longer than ``max_homopolymer``.

    ``strategy`` is ``"minimal"`` (default) or ``"maximal"``. Minimal applies the
    biggest GC reductions first and **stops as soon as the region is back inside
    the band** — you are paying to synthesise a gene, and a construct that
    differs from the original in twenty places is easier to reason about than one
    that differs in sixty. Maximal swaps every codon it can, which lowers GC
    further and is worth having when minimal is not enough.

    Returns the suggestions, the rewritten DNA, and the before/after GC of each
    region so you can see whether it was worth it.
    """
    dna = _clean_dna(dna)
    avoid_sites = vc.DEFAULT_AVOID_SITES if avoid_sites is None else avoid_sites
    if regions is None:
        regions = gc_regions(gc_windows(dna, gc_window), gc_window)
    rich = (regions[regions["direction"] == "GC-rich"]
            if len(regions) else pd.DataFrame())

    protein_before = vc.translate_dna(dna)
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    working = list(codons)
    syn = _synonyms()
    suggestions = []

    targets: set[int] = set()
    for _, region in rich.iterrows():
        first = max(0, (int(region["start"]) - 1) // 3)
        last = min(len(codons) - 1, (int(region["end"]) - 1) // 3)
        targets.update(range(first, last + 1))

    def region_within_band(sequence: str) -> bool:
        """True when no GC-rich stretch survives in the rewritten sequence."""
        remaining = gc_regions(gc_windows(sequence, gc_window), gc_window)
        return not len(remaining[remaining["direction"] == "GC-rich"]) \
            if len(remaining) else True

    # Minimal first applies the swaps that buy the most GC, so it can stop early
    order = sorted(targets)
    if strategy == "minimal":
        def gain(i: int) -> float:
            aa = vc.translate_dna(codons[i])
            alts = [c for c in syn.get(aa, []) if c != codons[i]]
            if not allow_rare:
                alts = [c for c in alts if c not in RARE_CODONS]
            return -(_gc(codons[i]) - min((_gc(c) for c in alts), default=_gc(codons[i])))
        order = sorted(targets, key=gain)

    for index in order:
        if strategy == "minimal" and suggestions and region_within_band(
                "".join(working)):
            break
        codon = working[index]
        aa = vc.translate_dna(codon)
        if aa in ("", "*", "X"):
            continue
        options = [c for c in syn.get(aa, []) if c != codon]
        if not allow_rare:
            options = [c for c in options if c not in RARE_CODONS]
        options = [c for c in options if _gc(c) < _gc(codon)]
        if not options:
            continue
        options.sort(key=lambda c: (_gc(c), vc.ECOLI_CODONS.get(aa, []).index(c)
                                    if c in vc.ECOLI_CODONS.get(aa, []) else 99))

        for candidate in options:
            trial = list(working)
            trial[index] = candidate
            spliced = "".join(trial)
            if re.search(rf"(A{{{max_homopolymer + 1},}}|C{{{max_homopolymer + 1},}}"
                         rf"|G{{{max_homopolymer + 1},}}|T{{{max_homopolymer + 1},}})",
                         spliced):
                continue
            window = spliced[max(0, index * 3 - 9):index * 3 + 12]
            if any(site in window for site in avoid_sites.values()):
                continue
            working[index] = candidate
            suggestions.append({
                "codon_index": index + 1,
                "nt_start": index * 3 + 1,
                "aa": aa,
                "from": codon, "to": candidate,
                "gc_from": round(_gc(codon), 2), "gc_to": round(_gc(candidate), 2),
                "rare_substitute": candidate in RARE_CODONS,
                "app_preference": (vc.ECOLI_CODONS.get(aa, []).index(candidate) + 1
                                   if candidate in vc.ECOLI_CODONS.get(aa, [])
                                   else None),
            })
            break

    if strategy == "minimal":
        suggestions.sort(key=lambda row: row["codon_index"])
    fixed = "".join(working) + dna[len(dna) - len(dna) % 3:]
    protein_after = vc.translate_dna(fixed)

    before_after = []
    for _, region in rich.iterrows():
        lo, hi = int(region["start"]) - 1, int(region["end"])
        before_after.append({
            "start": int(region["start"]), "end": int(region["end"]),
            "gc_before": round(_gc(dna[lo:hi]), 3),
            "gc_after": round(_gc(fixed[lo:hi]), 3),
            "n_changes": sum(1 for s in suggestions
                             if lo < s["nt_start"] <= hi + 2),
        })

    return {
        "suggestions": pd.DataFrame(suggestions),
        "dna": fixed,
        "protein_unchanged": protein_before == protein_after,
        "regions": pd.DataFrame(before_after),
        "gc_before": round(_gc(dna), 4),
        "gc_after": round(_gc(fixed), 4),
        "n_changes": len(suggestions),
        "remaining_regions": gc_regions(gc_windows(fixed, gc_window), gc_window),
    }




def plot_gc_comparison(before: str, after: str, suggestions=None,
                       gc_window: int = GC_WINDOW,
                       labels: tuple[str, str] = ("before", "after")):
    """The GC profile before and after the swaps, on one axis.

    The point of the figure is to show that the stretches that crossed the band
    no longer do, and at what cost — the changed codons are marked underneath so
    you can see where the edits landed.
    """
    import matplotlib.pyplot as plt

    before, after = _clean_dna(before), _clean_dna(after)
    gc_before = gc_windows(before, gc_window)
    gc_after = gc_windows(after, gc_window)

    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.axhspan(GC_LOW, GC_HIGH, color="0.92", zorder=0)
    ax.axhline(GC_LOW, color="0.65", ls="--", lw=0.8)
    ax.axhline(GC_HIGH, color="0.65", ls="--", lw=0.8)

    ax.plot(gc_before["start"], gc_before["gc"], color="#c44e52", lw=1.3,
            label=labels[0], zorder=2)
    ax.plot(gc_after["start"], gc_after["gc"], color="#4c72b0", lw=1.6,
            label=labels[1], zorder=3)
    ax.fill_between(gc_before["start"], gc_before["gc"], GC_HIGH,
                    where=gc_before["gc"] > GC_HIGH, color="#c44e52",
                    alpha=0.22, zorder=1, interpolate=True)

    if suggestions is not None and len(suggestions):
        positions = suggestions["nt_start"].to_numpy()
        ax.plot(positions, [0.04] * len(positions), "|", color="0.35", ms=9,
                zorder=4)
        ax.text(0.99, 0.11, f"{len(positions)} codon(s) changed",
                transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom",
                color="0.35")

    ax.set_ylim(0, 1)
    ax.set_xlim(gc_before["start"].min(), gc_before["start"].max())
    ax.set_xlabel("Position (nt)")
    ax.set_ylabel(f"GC fraction ({gc_window} nt window)", fontsize=9)
    ax.set_title("GC profile before and after the suggested changes", fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


def parse_dna_input(text: str, max_sequences: int = 20,
                    max_length: int = 6000) -> dict[str, str]:
    """Parse FASTA or bare DNA into ``{name: sequence}``.

    The protein side has its own parser in ``vhh_core``; this is the DNA
    equivalent, with its own limits because a coding sequence is three times the
    length of what it encodes.
    """
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

    cleaned: dict[str, str] = {}
    for name, raw in seqs.items():
        try:
            dna = _clean_dna(raw)
        except ValueError as exc:
            raise ValueError(f"{name}: {exc}") from None
        if len(dna) > max_length:
            raise ValueError(f"{name}: {len(dna)} nt exceeds the {max_length} nt limit")
        if len(dna) < 30:
            raise ValueError(f"{name}: {len(dna)} nt is too short to diagnose")
        key, n = name, 2
        while key in cleaned:
            key, n = f"{name}_{n}", n + 1
        cleaned[key] = dna

    if not cleaned:
        raise ValueError("No valid DNA sequence found")
    if len(cleaned) > max_sequences:
        raise ValueError(f"{len(cleaned)} sequences, limit is {max_sequences} per run")
    return cleaned


def summarise(seqs: dict[str, str], gc_window: int = GC_WINDOW,
              cluster_window: int = CLUSTER_WINDOW) -> pd.DataFrame:
    """One row per sequence — the overview for a batch."""
    rows = []
    for name, dna in seqs.items():
        try:
            d = diagnose(dna, cluster_window=cluster_window, gc_window=gc_window)
        except Exception as exc:
            rows.append({"name": name, "note": str(exc)})
            continue
        gc_rich = (int((d["gc_regions"]["direction"] == "GC-rich").sum())
                   if len(d["gc_regions"]) else 0)
        rows.append({
            "name": name,
            "length_nt": d["orf"]["length_nt"],
            "codons": d["orf"]["n_codons"],
            "in_frame": d["orf"]["in_frame"],
            "internal_stops": d["orf"]["n_internal_stops"],
            "gc": d["gc_overall"],
            "head_gc": round(d["head_gc"], 3),
            "rare_codons": d["n_rare"],
            "rare_clusters": len(d["clusters"]),
            "gc_rich_stretches": gc_rich,
            "homopolymers": len(d["homopolymers"]),
            "note": "",
        })
    return pd.DataFrame(rows)


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
