"""
sanger, extraction of VHH sequences from Sanger .ab1 chromatograms.

Reads .ab1 files, trims them on quality, then locates the expression cassette
(PelB leader → VHH → tags → stop) in any strand and any reading frame, and
returns the protein sequence of each VHH together with the problems found along
the way: frameshifts, premature stops, missing tags, truncated reads.

The module has no Streamlit dependency: the web app and a notebook both import
it the same way.

Main entry points
-----------------
read_ab1(data, name)            parse one .ab1 into a read dict
extract_cassette(nt, vector)    locate PelB / VHH / tags in a nucleotide read
process_reads(files, ...)       the whole pipeline over several files
plot_quality(reads)             per-read Phred trace with the trimmed window
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

__all__ = [
    "VECTORS", "DEFAULTS", "mott_trim", "read_ab1", "extract_cassette",
    "process_reads", "plot_quality", "quality_table", "reads_to_fasta",
]

# Bumped whenever the public API changes, so app.py can tell the user that a
# module is out of date instead of crashing on a missing attribute.
__version__ = "2026.11.17a"

# --------------------------------------------------------------------------- #
# Vector definitions
# --------------------------------------------------------------------------- #

# Error probability Mott's modified trimming algorithm scores each base
# against. Not a choice made here: it is the value Biopython's own ABI
# reader hard-codes, and the register says so.
MOTT_CUTOFF = 0.05

PELB_AA = "MKYLLPTAAAGLLLLAAQPAMA"      # pelB leader (pMES3 / pMECS)
PELB_CORE = "AAQPAMA"                   # short anchor, survives truncated reads

# The HA epitope on record is the nine-residue YPYDVPDYA, and that is what this
# module looked for. pMECS does not carry it: its map annotates `HA-tag` over 24
# nucleotides, which encode the eight residues YPYDVPDY, and the two codons
# after them are GS, not A. The nine-residue form appears nowhere in the pMECS
# map, in any of the six frames, so every correct pMECS clone came back with
# `missing tag(s): HA`, on the sequencing screen and in the cloning screen's
# construct check alike. That is the exact failure the comment under VECTORS
# warned about, and it was sitting in VECTORS.
#
# What is searched for is now the eight residues the vector really encodes. They
# are a prefix of the nine-residue form, so a construct carrying the canonical
# epitope still matches: the change makes the test accept pMECS without making
# it reject anything it used to accept. What it cannot tell apart is a cassette
# that stops after the eighth residue from one that carries the ninth, which for
# these vectors is the distinction that does not exist.
#
# Whether an anti-HA antibody detects YPYDVPDY followed by GS as well as it
# detects YPYDVPDYA is a question about the antibody, and this module does not
# answer it. It reports what the vector encodes.
# Both forms, kept here as named constants because the layouts are no longer
# written out in this module and these two are what `tests/test_derive_layout.py`
# checks `cloning.PROTEIN_MOTIFS` against, so the motif table and this module
# cannot go back to disagreeing about what an HA tag is.
HA_AA = "YPYDVPDYA"                     # the epitope itself, nine residues
HA_SEARCH = "YPYDVPDY"                  # what pMECS and pMES4 encode, eight
HIS_AA = "HHHHHH"

def _bundled_layouts() -> dict:
    """Read every bundled plasmid map and derive its cassette layout.

    These layouts used to be written out here, two of them, pMES3 and pMECS,
    with their tag sequences typed in by hand. That was a second account of
    plasmid maps the app already ships in `data/plasmids/`, and the two
    disagreed: the pMECS map annotates an HA tag whose 24 nucleotides encode
    YPYDVPDY, while the entry here carried the nine-residue YPYDVPDYA, so every
    correct pMECS clone was reported as missing its HA tag. Fixing that entry by
    hand left the two lists separate and free to drift again. Reading the maps
    removes the second list and gives the sequencing screen the three vectors it
    never had, pMES4, pHEN6 and pHEN25.

    Measured at 0.06 s for the five bundled maps on this container, once, at
    import. A map that will not parse is skipped rather than breaking the import:
    the cloning screen loads the same file itself and reports the failure there,
    with the file name.
    """
    import cloning as cl                  # local: cloning imports sanger lazily

    layouts = {}
    for name, path in cl.list_bundled_plasmids().items():
        try:
            derived = cl.derive_layout(cl.load_plasmid(path))
        except Exception:
            continue
        layouts[name] = {"layout": derived["layout"], "tags": derived["tags"],
                         "tag_table": derived["tag_table"]}
    return layouts


VECTORS = _bundled_layouts()

# To add a vector, drop its map in `data/plasmids/`: the layout is read from the
# file, so there is no second list to keep in step with it.
# `tests/test_vector_layouts.py` pins what the five bundled maps derive to,
# including the pMECS HA case, so a change in the derivation that alters a
# layout fails the suite instead of being found at the bench.

# C-terminal FR4 signature delimiting the end of the VHH
FR4_RE = re.compile(r"[WF]G[A-Z]G[TA][A-Z][VLIMQ]TV[STF][SPAL]")

DEFAULTS = {
    "vector": "pMES3",
    "mott_cutoff": 0.05,      # Mott trimming error-probability cutoff (~Q13)
    "min_trimmed_len": 150,   # reads shorter than this after trimming are flagged
    "min_mean_phred": 20,     # QC threshold on the trimmed window
}


# --------------------------------------------------------------------------- #
# Reading and trimming
# --------------------------------------------------------------------------- #

def mott_trim(quals, cutoff: float = MOTT_CUTOFF) -> tuple[int, int]:
    """Mott's modified trimming algorithm, the one used by phred and Staden.

    Each base scores ``cutoff - P(error)``; the highest-scoring contiguous run
    is the good window. Returns ``(start, end)``, end exclusive.
    """
    quals = np.asarray(quals, dtype=float)
    if not len(quals):
        return 0, 0
    scores = cutoff - 10 ** (-quals / 10.0)
    best_sum, best = -1.0, (0, 0)
    run_sum, run_start = 0.0, 0
    for i, s in enumerate(scores):
        run_sum += s
        if run_sum < 0:
            run_sum, run_start = 0.0, i + 1
        elif run_sum > best_sum:
            best_sum, best = run_sum, (run_start, i + 1)
    return best if best_sum > 0 else (0, 0)


def read_ab1(data, name: str | None = None, cutoff: float | None = None) -> dict:
    """Parse one .ab1 chromatogram and apply quality trimming.

    ``data`` is a path, a bytes object or a file-like object, so the same
    function serves a notebook reading from disk and a web upload held in memory.
    """
    cutoff = DEFAULTS["mott_cutoff"] if cutoff is None else cutoff

    if isinstance(data, (str, Path)):
        handle, stem = str(data), Path(data).stem
        record = SeqIO.read(handle, "abi")
    else:
        raw = data if isinstance(data, bytes) else data.read()
        stem = Path(name).stem if name else "read"
        record = SeqIO.read(io.BytesIO(raw), "abi")

    quals = np.array(record.letter_annotations.get("phred_quality", []), dtype=float)
    seq = str(record.seq).upper()
    start, end = mott_trim(quals, cutoff)

    return {
        "read": stem, "file_stem": stem,
        "raw_seq": seq, "raw_len": len(seq),
        "trim_start": int(start), "trim_end": int(end),
        "trimmed_seq": seq[start:end], "trimmed_len": int(end - start),
        "mean_phred_raw": float(np.mean(quals)) if len(quals) else np.nan,
        "mean_phred_trim": float(np.mean(quals[start:end])) if end > start else np.nan,
        "quals": quals,
    }


# --------------------------------------------------------------------------- #
# Cassette extraction
# --------------------------------------------------------------------------- #

def three_frame_translations(nt: str) -> dict:
    """{(strand, frame): (protein, nt_offset, nt_string)} for both strands."""
    out = {}
    for strand, s in (("+", nt), ("-", str(Seq(nt).reverse_complement()))):
        for frame in range(3):
            sub = s[frame:]
            sub = sub[: len(sub) - len(sub) % 3]
            if len(sub) < 3:
                continue
            out[(strand, frame)] = (str(Seq(sub).translate()), frame, s)
    return out


def _find_pelb(protein: str) -> tuple[int, int, str | None]:
    """Locate the PelB leader; returns (pelb_start, vhh_start, mode).

    The VHH start is derived from whichever anchor actually matched, rather than
    from a back-computed PelB start. When a read begins part-way through the
    leader, only the short core anchor matches, and back-computing the full
    leader start would run off the beginning of the read and shift the VHH start
    by a few residues, silently truncating the N-terminus.
    """
    i = protein.find(PELB_AA)
    if i >= 0:
        return i, i + len(PELB_AA), "full"
    i = protein.find(PELB_CORE)
    if i >= 0:
        return max(0, i - (len(PELB_AA) - len(PELB_CORE))), i + len(PELB_CORE), "core"
    return -1, -1, None


def tags_for(vector) -> list[tuple[str, str]]:
    """The (name, protein) tags to look for, from a vector name or a layout.

    The cloning screen can load a vector that ships with no entry here at all,
    an upload, so it passes the layout it derived from that file rather than a
    name. Accepting both keeps one code path for the sequencing screen, which
    picks a bundled vector by name, and the cloning screen, which does not.
    An unknown *name* yields no tags, which reports "nothing to check" instead
    of silently checking pMES3's tags against someone else's vector.
    """
    if isinstance(vector, str):
        return list(VECTORS.get(vector, {}).get("tags", []))
    if isinstance(vector, dict):
        return list(vector.get("tags", []))
    return [(str(n), str(s)) for n, s in (vector or ())]


def extract_cassette(nt: str, vector: str = "pMES3") -> dict:
    """Locate the PelB–VHH–tags cassette in a read, in any strand and frame.

    The read and its reverse complement are translated in three frames each.
    The PelB anchor fixes the start and the reading frame; the FR4 signature
    fixes the end of the VHH. If PelB and FR4 land in different frames, an indel
    sits between them and the VHH cannot be trusted, the read is flagged as a
    frameshift rather than silently returning a wrong sequence.
    """
    res = {
        "strand": None, "frame": None, "pelb_found": False, "pelb_mode": None,
        "fr4_found": False, "frameshift": False, "frameshift_hint": "",
        "n_premature_stops": 0, "premature_stop_positions": "",
        "stop_codon_present": False, "vhh_nt": "", "vhh_aa": "",
        "cassette_aa": "", "cassette_nt": "", "tags_found": "", "tags_missing": "",
        "warnings": [],
    }
    if not nt or len(nt) < 60:
        res["warnings"].append("read too short to analyse")
        return res

    frames = three_frame_translations(nt)

    pelb_hits = {}
    for key, (protein, offset, strand_nt) in frames.items():
        i, vhh_i, mode = _find_pelb(protein)
        if i >= 0:
            pelb_hits[key] = (i, vhh_i, mode, protein, offset, strand_nt)
    if not pelb_hits:
        res["warnings"].append("PelB leader not found in any frame")
        return res

    # prefer a full-length PelB match, then the earliest position
    key = sorted(pelb_hits, key=lambda k: (pelb_hits[k][2] != "full", pelb_hits[k][0]))[0]
    pelb_i, vhh_start, pelb_mode, protein, offset, strand_nt = pelb_hits[key]
    res.update(strand=key[0], frame=key[1], pelb_found=True, pelb_mode=pelb_mode)
    if pelb_mode == "core":
        res["warnings"].append(
            "only the PelB core anchor matched, read starts inside the leader")

    match = FR4_RE.search(protein, pelb_i)
    if match is None:
        for other, (prot2, _off2, _s2) in frames.items():
            if other[0] != key[0] or other == key:
                continue
            if FR4_RE.search(prot2):
                res["frameshift"] = True
                res["frameshift_hint"] = (
                    f"PelB in frame {key[1]}, FR4 in frame {other[1]}, "
                    f"net {(other[1] - key[1]) % 3} nt indel between them")
                res["warnings"].append("frameshift (indel) suspected")
                break
        if not res["frameshift"]:
            res["warnings"].append(
                "FR4 signature (WGxGTxVTVSS) not found, truncated read or heavy mutation")
        return res

    res["fr4_found"] = True
    vhh_end = match.end()
    res["vhh_aa"] = protein[vhh_start:vhh_end]
    res["vhh_nt"] = strand_nt[offset + 3 * vhh_start: offset + 3 * vhh_end]

    downstream = protein[vhh_end:]
    stop_i = downstream.find("*")
    res["stop_codon_present"] = stop_i >= 0
    cassette_aa = protein[pelb_i: vhh_end + (stop_i + 1 if stop_i >= 0 else len(downstream))]
    res["cassette_aa"] = cassette_aa
    res["cassette_nt"] = strand_nt[offset + 3 * pelb_i:
                                   offset + 3 * (pelb_i + len(cassette_aa))]

    body = cassette_aa[:-1] if cassette_aa.endswith("*") else cassette_aa
    positions = [i for i, c in enumerate(body) if c == "*"]
    res["n_premature_stops"] = len(positions)
    res["premature_stop_positions"] = ";".join(str(p + 1) for p in positions)
    if positions:
        res["warnings"].append(f"{len(positions)} premature stop codon(s)")

    tail = downstream[:stop_i] if stop_i >= 0 else downstream
    found, missing = [], []
    for tag_name, tag_seq in tags_for(vector):
        (found if tag_seq in tail else missing).append(tag_name)
    res["tags_found"] = ",".join(found)
    res["tags_missing"] = ",".join(missing)
    if missing:
        res["warnings"].append(f"missing tag(s): {','.join(missing)}")
    if not res["stop_codon_present"]:
        res["warnings"].append("no stop codon after the cassette (read too short?)")

    # A VHH that still contains a stop is not usable even if everything else fit
    if "*" in res["vhh_aa"]:
        res["warnings"].append("stop codon inside the VHH")

    return res


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #

def process_reads(files, vector: str | None = None, mott_cutoff: float | None = None,
                  min_trimmed_len: int | None = None, min_mean_phred: float | None = None,
                  progress=None) -> tuple[pd.DataFrame, list[dict]]:
    """Run the whole pipeline over several .ab1 inputs.

    ``files`` is a list of paths, or of ``(name, bytes)`` pairs for uploads.
    ``progress`` is an optional callable ``(fraction, label)`` for a progress bar.
    Returns ``(results_dataframe, reads)`` where reads keep the quality arrays
    needed for the chromatogram figure.
    """
    vector = vector or DEFAULTS["vector"]
    mott_cutoff = DEFAULTS["mott_cutoff"] if mott_cutoff is None else mott_cutoff
    min_len = DEFAULTS["min_trimmed_len"] if min_trimmed_len is None else min_trimmed_len
    min_phred = DEFAULTS["min_mean_phred"] if min_mean_phred is None else min_mean_phred

    rows, reads = [], []
    total = max(len(files), 1)
    for i, item in enumerate(files):
        if isinstance(item, (tuple, list)):
            name, data = item
        else:
            name, data = Path(item).name, item
        if progress:
            progress(i / total, f"Reading {name} ({i + 1}/{total})…")

        try:
            read = read_ab1(data, name=name, cutoff=mott_cutoff)
        except Exception as exc:
            rows.append({"read": Path(name).stem, "status": "unreadable",
                         "warnings": f"could not parse: {exc}", "qc_pass": False,
                         "vhh_aa": "", "usable": False})
            continue

        reads.append(read)
        # a read too short after trimming is analysed raw rather than discarded
        seq = read["trimmed_seq"] if read["trimmed_len"] >= min_len else read["raw_seq"]
        ex = extract_cassette(seq, vector)

        qc_pass = (read["trimmed_len"] >= min_len
                   and (read["mean_phred_trim"] or 0) >= min_phred)
        usable = bool(ex["vhh_aa"]) and not ex["frameshift"] and "*" not in ex["vhh_aa"]
        status = ("ok" if usable and qc_pass else
                  "low quality" if usable else
                  "frameshift" if ex["frameshift"] else
                  "no PelB" if not ex["pelb_found"] else
                  "no FR4" if not ex["fr4_found"] else "unusable")

        rows.append({
            "read": read["read"], "status": status, "usable": usable,
            "qc_pass": bool(qc_pass),
            "raw_len": read["raw_len"], "trimmed_len": read["trimmed_len"],
            "mean_phred_trim": round(read["mean_phred_trim"], 1)
            if read["mean_phred_trim"] == read["mean_phred_trim"] else np.nan,
            "strand": ex["strand"], "frame": ex["frame"],
            "pelb": ex["pelb_mode"] or "", "fr4": ex["fr4_found"],
            "vhh_len": len(ex["vhh_aa"]), "vhh_aa": ex["vhh_aa"],
            "tags_found": ex["tags_found"], "tags_missing": ex["tags_missing"],
            "premature_stops": ex["n_premature_stops"],
            "stop_codon": ex["stop_codon_present"],
            "frameshift": ex["frameshift"],
            "vhh_nt": ex["vhh_nt"],
            "warnings": "; ".join(ex["warnings"]) or "",
        })

    if progress:
        progress(1.0, "Done")

    df = pd.DataFrame(rows)
    if len(df):
        # identical protein sequences = the same clone picked more than once
        df["duplicate_of"] = ""
        seen: dict[str, str] = {}
        for idx, row in df.iterrows():
            aa = row.get("vhh_aa") or ""
            if not aa:
                continue
            if aa in seen:
                df.at[idx, "duplicate_of"] = seen[aa]
            else:
                seen[aa] = row["read"]
        df["n_copies"] = df["vhh_aa"].map(
            df[df["vhh_aa"].astype(bool)]["vhh_aa"].value_counts()).fillna(0).astype(int)
    return df, reads


def reads_to_fasta(df: pd.DataFrame, only_usable: bool = True,
                   column: str = "vhh_aa") -> str:
    """FASTA of the extracted sequences, protein by default, DNA with 'vhh_nt'."""
    sub = df[df["usable"]] if only_usable and "usable" in df.columns else df
    parts = []
    for _, row in sub.iterrows():
        seq = row.get(column) or ""
        if not seq:
            continue
        wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
        parts.append(f">{row['read']}\n{wrapped}")
    return "\n".join(parts) + ("\n" if parts else "")


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #

def quality_table(reads: list[dict], max_reads: int = 12) -> pd.DataFrame:
    """Phred quality per base for each read: the data of the quality figure.

    One row per base, with ``kept`` marking the bases inside the trimmed
    window, which is the blue part of the figure. The same ``max_reads`` slice
    the figure applies is applied here, so the table covers what was drawn and
    not more; a run with more reads than that is trimmed in both.

    Measured, in the sense that the Phred scores come out of the basecaller's
    own estimate in the .ab1 file. What the basecaller estimates is the
    probability it called the wrong base, which is not a measurement of the
    DNA.
    """
    rows = []
    for read in reads[:max_reads]:
        quals = read["quals"]
        start, end = read["trim_start"], read["trim_end"]
        for i, q in enumerate(quals):
            rows.append({"read": read["read"], "position": i + 1,
                         "phred": int(q), "kept": bool(start <= i < end)})
    return pd.DataFrame(rows, columns=["read", "position", "phred", "kept"])


def plot_quality(reads: list[dict], max_reads: int = 12):
    """Phred quality along each read, with the trimmed window highlighted."""
    reads = reads[:max_reads]
    if not reads:
        return None
    fig, axes = plt.subplots(len(reads), 1, figsize=(11, 1.5 * len(reads)),
                             squeeze=False, sharex=True)
    for ax, read in zip(axes[:, 0], reads):
        quals = read["quals"]
        x = np.arange(1, len(quals) + 1)
        ax.fill_between(x, quals, color="#b0b0b0", lw=0, step="mid")
        start, end = read["trim_start"], read["trim_end"]
        if end > start:
            ax.fill_between(x[start:end], quals[start:end], color="#4c72b0",
                            lw=0, step="mid")
        ax.axhline(20, color="#c44e52", lw=0.8, ls="--")
        ax.set_ylim(0, 62)
        ax.set_ylabel("Phred", fontsize=7)
        ax.set_title(f"{read['read']}, {read['trimmed_len']} nt kept, "
                     f"mean Q {read['mean_phred_trim']:.0f}", fontsize=8, loc="left")
    axes[-1, 0].set_xlabel("Base position")
    fig.tight_layout()
    return fig
