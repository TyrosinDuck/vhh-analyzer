"""
cloning — plasmid handling for in-silico restriction cloning of VHHs.

Phase 1: load a plasmid (SnapGene .dna, GenBank or FASTA), inventory what is in
it, and draw its map.

Two sources of truth, deliberately kept apart:

* **Annotated files** (SnapGene, GenBank) already carry their features. They are
  read as they are — nothing is guessed, and that is by far the reliable path.
* **Bare FASTA** carries only sequence, so features have to be found. This
  module detects *short, unambiguous protein motifs* — the PelB leader, His6,
  HA, Myc, a VHH scaffold — by translating in all six frames. It deliberately
  does **not** try to guess large elements such as gIII, AmpR or the origin:
  doing that properly needs a reference sequence database, and a confident
  wrong annotation on a cloning tool is worse than no annotation at all.

The detected-features table always says where each feature came from, so an
annotation read from the file is never confused with one inferred here.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__all__ = [
    "PROTEIN_MOTIFS", "load_plasmid", "detect_motifs", "feature_table",
    "restriction_profile", "unique_cutters", "plot_map", "sequence_stats",
]

# --------------------------------------------------------------------------- #
# Motif library — short protein elements only, see the module docstring
# --------------------------------------------------------------------------- #

PROTEIN_MOTIFS = {
    "PelB leader": ("MKYLLPTAAAGLLLLAAQPAMA", "signal peptide, periplasmic export"),
    "PelB core": ("AAQPAMA", "C-terminal part of the PelB leader"),
    "His6": ("HHHHHH", "purification tag"),
    "HA tag": ("YPYDVPDYA", "detection tag"),
    "Myc tag": ("EQKLISEEDL", "detection tag"),
    "FLAG tag": ("DYKDDDDK", "detection tag"),
    "Strep-tag II": ("WSHPQFEK", "purification tag"),
    "TEV site": ("ENLYFQ", "protease cleavage site"),
    "Thrombin site": ("LVPRGS", "protease cleavage site"),
    "Amber stop (TAG)": (None, "suppressible stop, phage display vectors"),
}

# VHH scaffold signatures, reused from the Sanger module's logic
VHH_FR1_RE = re.compile(r"[QEDA][VIL][QKE]L[VQKAL]")
VHH_FR2_RE = re.compile(r"W[VYFGIL]R[QKR]")
VHH_FR4_RE = re.compile(r"[WF]G[A-Z]G[TA][A-Z][VLIMQ]TV[STF][SPAL]")

FEATURE_COLORS = {
    "file": "#4c72b0", "motif": "#dd8452", "vhh": "#4daf4a",
    "site": "#c44e52", "other": "#b0b0b0",
}


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

def _guess_format(name: str, head: bytes) -> str:
    suffix = Path(name).suffix.lower()
    if suffix == ".dna":
        return "snapgene"
    if suffix in (".gb", ".gbk", ".genbank", ".ape"):
        return "genbank"
    if suffix in (".fa", ".fasta", ".fna", ".seq", ".txt"):
        return "fasta"
    if head.startswith(b">"):
        return "fasta"
    if head.startswith(b"LOCUS"):
        return "genbank"
    return "snapgene"


def load_plasmid(data, name: str = "plasmid", topology: str | None = None) -> dict:
    """Read a plasmid from SnapGene .dna, GenBank or FASTA.

    ``data`` is a path, bytes or a file-like object. Returns a dict with the
    SeqRecord, the format actually used, the topology and any parse warning.
    SnapGene and GenBank carry their own features; FASTA does not, which the
    caller is expected to surface to the user.
    """
    if isinstance(data, (str, Path)):
        raw = Path(data).read_bytes()
        name = Path(data).name
    elif isinstance(data, bytes):
        raw = data
    else:
        raw = data.read()

    fmt = _guess_format(name, raw[:16])
    warning = ""
    record, last = None, None
    # SnapGene is binary; GenBank and FASTA are text, and Biopython needs the
    # matching kind of handle for each.
    for attempt in dict.fromkeys([fmt, "genbank", "fasta", "snapgene"]):
        try:
            handle = (io.BytesIO(raw) if attempt == "snapgene"
                      else io.StringIO(raw.decode("utf-8", errors="replace")))
            record = SeqIO.read(handle, attempt)
            fmt = attempt
            break
        except Exception as exc:
            last = exc
    if record is None:
        raise ValueError(f"could not parse {name} as SnapGene, GenBank or FASTA: {last}")

    if topology is None:
        topo = str(record.annotations.get("topology", "")).lower()
        topology = "circular" if topo.startswith("circ") else (
            "circular" if fmt == "snapgene" else "unknown")
    if fmt == "fasta":
        warning = ("FASTA carries no features — only short protein motifs can be "
                   "detected, and the topology is assumed circular. Export the map "
                   "from SnapGene as .dna or GenBank for an exact inventory.")
        topology = topology if topology != "unknown" else "circular"

    record.seq = Seq(str(record.seq).upper().replace("U", "T"))
    stem = Path(name).stem
    rec_name = (record.name or "").strip()
    if not rec_name or rec_name.lower().startswith("<unknown"):
        rec_name = stem
    return {"record": record, "name": rec_name,
            "format": fmt, "topology": topology, "length": len(record.seq),
            "warning": warning}


def sequence_stats(record: SeqRecord) -> dict:
    seq = str(record.seq)
    gc = 100 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0
    return {"length": len(seq), "gc_percent": round(gc, 1)}


# --------------------------------------------------------------------------- #
# Feature detection
# --------------------------------------------------------------------------- #

def _six_frames(seq: str):
    """Yield (strand, frame, protein, nt_offset) for all six reading frames."""
    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            sub = s[frame:]
            sub = sub[: len(sub) - len(sub) % 3]
            if len(sub) >= 3:
                yield strand, frame, str(Seq(sub).translate()), frame


def _to_plus(start_nt: int, length_nt: int, strand: str, total: int) -> tuple[int, int]:
    """Map a hit on the reverse strand back to plus-strand coordinates."""
    if strand == "+":
        return start_nt, start_nt + length_nt
    end = total - start_nt
    return end - length_nt, end


def detect_motifs(record: SeqRecord, circular: bool = True) -> list[dict]:
    """Find short protein motifs by translating in all six frames.

    Returns one entry per hit with plus-strand nucleotide coordinates. For a
    circular plasmid the sequence is searched doubled, so a motif spanning the
    origin of numbering is still found.
    """
    seq = str(record.seq)
    total = len(seq)
    search = seq + seq if circular else seq
    hits = []

    for strand, _frame, protein, offset in _six_frames(search):
        for label, (motif, note) in PROTEIN_MOTIFS.items():
            if not motif:
                continue
            for m in re.finditer(f"(?=({re.escape(motif)}))", protein):
                nt_start = offset + 3 * m.start()
                if nt_start >= total:           # duplicate from the doubled copy
                    continue
                s0, s1 = _to_plus(nt_start, 3 * len(motif), strand,
                                  len(search) if strand == "-" else total)
                hits.append({"feature": label, "kind": "motif", "source": "motif scan",
                             "start": s0 % total + 1, "end": (s1 - 1) % total + 1,
                             "strand": strand, "note": note,
                             "protein": motif})

        # a VHH is recognised by its own scaffold, not by a fixed motif
        for m2 in VHH_FR4_RE.finditer(protein):
            window_start = max(0, m2.start() - 130)
            window = protein[window_start:m2.end()]
            if not VHH_FR2_RE.search(window) or "*" in window:
                continue
            # the 130-residue look-back overshoots the N-terminus and swallows
            # whatever precedes the domain (typically the leader), so trim to the
            # FR1 start if it is there
            fr1 = VHH_FR1_RE.search(window)
            if fr1:
                window_start += fr1.start()
                window = protein[window_start:m2.end()]
            nt_start = offset + 3 * window_start
            if nt_start >= total:
                continue
            s0, s1 = _to_plus(nt_start, 3 * len(window), strand,
                              len(search) if strand == "-" else total)
            hits.append({"feature": "VHH (scaffold)", "kind": "vhh",
                         "source": "motif scan", "start": s0 % total + 1,
                         "end": (s1 - 1) % total + 1, "strand": strand,
                         "note": "FR2 and FR4 signatures in frame",
                         "protein": window})

    # drop duplicates (same feature, same coordinates, found in several frames)
    seen, unique = set(), []
    for h in sorted(hits, key=lambda x: (x["start"], x["feature"])):
        key = (h["feature"], h["start"], h["end"], h["strand"])
        if key not in seen:
            seen.add(key)
            unique.append(h)

    # PelB core is redundant when the full leader was found at the same place
    full = {(h["start"], h["strand"]) for h in unique if h["feature"] == "PelB leader"}
    unique = [h for h in unique
              if not (h["feature"] == "PelB core"
                      and any(abs(h["start"] - s) < 60 and h["strand"] == st
                              for s, st in full))]
    return unique


def feature_table(plasmid: dict) -> pd.DataFrame:
    """Inventory of a plasmid: features from the file plus detected motifs.

    The ``source`` column distinguishes what the file declared from what this
    module inferred — they do not carry the same confidence.
    """
    record = plasmid["record"]
    total = len(record.seq)
    rows = []

    for feat in record.features:
        if feat.type in ("source",):
            continue
        qualifiers = feat.qualifiers or {}
        label = (qualifiers.get("label") or qualifiers.get("gene")
                 or qualifiers.get("product") or qualifiers.get("note")
                 or [feat.type])[0]
        rows.append({
            "feature": str(label), "kind": "file", "source": f"{plasmid['format']} file",
            "start": int(feat.location.start) + 1, "end": int(feat.location.end),
            "strand": "-" if feat.location.strand == -1 else "+",
            "type": feat.type,
            "note": "; ".join(qualifiers.get("note", []))[:120],
        })

    for hit in detect_motifs(record, circular=plasmid["topology"] == "circular"):
        rows.append({**hit, "type": "motif"})

    df = pd.DataFrame(rows)
    if len(df):
        df["length_nt"] = (df["end"] - df["start"] + 1) % (total or 1)
        df = df.sort_values(["start", "feature"]).reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
# Restriction analysis
# --------------------------------------------------------------------------- #

def _enzyme_batch(names: list[str] | None):
    from Bio.Restriction import CommOnly, RestrictionBatch

    if names:
        return RestrictionBatch(names)
    return CommOnly          # commercially available enzymes


def restriction_profile(record: SeqRecord, enzymes: list[str] | None = None,
                        circular: bool = True) -> pd.DataFrame:
    """Cut positions of every enzyme that cuts the sequence at least once."""
    batch = _enzyme_batch(enzymes)
    analysis = batch.search(record.seq, linear=not circular)
    rows = []
    for enzyme, sites in analysis.items():
        if not sites:
            continue
        rows.append({"enzyme": str(enzyme), "site": str(enzyme.site),
                     "n_cuts": len(sites),
                     "positions": ", ".join(str(p) for p in sorted(sites)[:12]),
                     "overhang": ("blunt" if enzyme.is_blunt() else
                                  "5'" if enzyme.is_5overhang() else "3'"),
                     "cut_positions": sorted(sites)})
    return (pd.DataFrame(rows).sort_values(["n_cuts", "enzyme"]).reset_index(drop=True)
            if rows else pd.DataFrame(columns=["enzyme", "site", "n_cuts", "positions",
                                               "overhang", "cut_positions"]))


def unique_cutters(record: SeqRecord, enzymes: list[str] | None = None,
                   circular: bool = True, within: tuple[int, int] | None = None
                   ) -> pd.DataFrame:
    """Enzymes cutting exactly once, optionally restricted to a window.

    ``within`` is a 1-based inclusive ``(start, end)`` range, typically the
    multiple cloning site, so the result is the list of enzymes that open the
    backbone exactly where the insert should go.
    """
    profile = restriction_profile(record, enzymes, circular)
    if profile.empty:
        return profile
    single = profile[profile["n_cuts"] == 1].copy()
    if within and len(single):
        lo, hi = within
        single = single[single["cut_positions"].map(lambda p: lo <= p[0] <= hi)]
    return single.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Map
# --------------------------------------------------------------------------- #

def plot_map(plasmid: dict, features: pd.DataFrame | None = None,
             sites: pd.DataFrame | None = None, max_features: int = 40,
             circular: bool | None = None):
    """Draw the plasmid map, optionally with chosen restriction sites marked."""
    from dna_features_viewer import CircularGraphicRecord, GraphicFeature, GraphicRecord

    record = plasmid["record"]
    total = len(record.seq)
    features = feature_table(plasmid) if features is None else features
    circular = (plasmid["topology"] == "circular") if circular is None else circular

    graphic = []
    for _, row in features.head(max_features).iterrows():
        start, end = int(row["start"]) - 1, int(row["end"])
        if end <= start:                      # wraps the origin — clip for drawing
            end = total
        graphic.append(GraphicFeature(
            start=start, end=end, strand=-1 if row["strand"] == "-" else 1,
            label=str(row["feature"])[:28],
            color=FEATURE_COLORS.get(row.get("kind", "other"), FEATURE_COLORS["other"])))

    if sites is not None and len(sites):
        for _, row in sites.iterrows():
            pos = row["cut_positions"][0]
            graphic.append(GraphicFeature(
                start=max(0, pos - 1), end=min(total, pos + len(row["site"])),
                strand=0, label=str(row["enzyme"]), color=FEATURE_COLORS["site"]))

    klass = CircularGraphicRecord if circular else GraphicRecord
    gr = klass(sequence_length=total, features=graphic)
    fig, ax = plt.subplots(figsize=(9, 9) if circular else (13, 3.6))
    gr.plot(ax=ax, with_ruler=not circular)
    ax.set_title(f"{plasmid['name']} — {total} bp, {plasmid['topology']}", fontsize=11)
    fig.tight_layout()
    return fig


def plot_region(plasmid: dict, features: pd.DataFrame | None = None,
                sites: pd.DataFrame | None = None, start: int | None = None,
                end: int | None = None, padding: int = 150):
    """Linear zoom on the cloning region.

    A circular map stacks its labels when every feature sits in a small arc,
    which is exactly the case for an expression cassette. Defaults to the span
    of the annotated features plus some padding.
    """
    from dna_features_viewer import GraphicFeature, GraphicRecord

    record = plasmid["record"]
    total = len(record.seq)
    features = feature_table(plasmid) if features is None else features

    if start is None or end is None:
        window = insertion_region(plasmid, features)
        if window:
            start = max(0, window[0] - padding)
            end = min(total, window[1] + padding)
        elif len(features):
            start = max(0, int(features["start"].min()) - padding)
            end = min(total, int(features["end"].max()) + padding)
        else:
            start, end = 0, total
    start, end = int(start), int(end)
    if end - start < 60:                      # keep very narrow windows readable
        pad = (60 - (end - start)) // 2 + 1
        start, end = max(0, start - pad), min(total, end + pad)
    span = end - start

    graphic = []
    for _, row in features.iterrows():
        f_start, f_end = int(row["start"]) - 1, int(row["end"])
        if f_end < start or f_start > end:
            continue
        graphic.append(GraphicFeature(
            start=max(f_start, start), end=min(f_end, end),
            strand=-1 if row["strand"] == "-" else 1,
            label=str(row["feature"])[:28],
            color=FEATURE_COLORS.get(row.get("kind", "other"), FEATURE_COLORS["other"])))

    if sites is not None and len(sites):
        for _, row in sites.iterrows():
            pos = row["cut_positions"][0]
            if start <= pos <= end:
                graphic.append(GraphicFeature(
                    start=pos - 1, end=pos + len(row["site"]), strand=0,
                    label=str(row["enzyme"]), color=FEATURE_COLORS["site"]))

    gr = GraphicRecord(sequence_length=total, features=graphic,
                       first_index=0).crop((start, end))
    fig, ax = plt.subplots(figsize=(max(9, min(18, span / 90)), 3.4))
    gr.plot(ax=ax, with_ruler=True)
    ax.set_title(f"{plasmid['name']} — cloning region {start + 1}–{end}", fontsize=10)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------- #
# Bundled plasmid and primer library
# --------------------------------------------------------------------------- #

DATA_DIR = Path(__file__).parent / "data"
PLASMID_DIR = DATA_DIR / "plasmids"
PRIMER_FILE = DATA_DIR / "primers.csv"


def list_bundled_plasmids() -> dict[str, Path]:
    """{name: path} of the plasmids shipped with the app."""
    if not PLASMID_DIR.is_dir():
        return {}
    return {p.stem: p for p in sorted(PLASMID_DIR.glob("*.dna"))}


def load_primers(path: str | Path | None = None) -> pd.DataFrame:
    """The bundled primer list, or a user-supplied CSV with name/sequence columns."""
    path = Path(path) if path else PRIMER_FILE
    if not path.exists():
        return pd.DataFrame(columns=["name", "sequence", "note"])
    df = pd.read_csv(path)
    df["sequence"] = df["sequence"].astype(str).str.upper().str.replace(r"[^ACGT]", "",
                                                                        regex=True)
    if "note" not in df.columns:
        df["note"] = ""
    return df


# --------------------------------------------------------------------------- #
# Primer binding
# --------------------------------------------------------------------------- #

def primer_sites(record: SeqRecord, primer: str, circular: bool = True) -> list[dict]:
    """Exact binding sites of a primer, with orientation, position and Tm.

    A primer written as the plus strand anneals to the minus strand and extends
    to the right (forward). A primer matching the reverse complement of the plus
    strand anneals to the plus strand and extends to the left (reverse). Both
    are reported, because which one a given primer is depends on where the
    target sits, not on the sequence alone.
    """
    from Bio.SeqUtils import MeltingTemp as mt

    primer = re.sub(r"[^ACGT]", "", primer.upper())
    if not primer:
        return []
    seq = str(record.seq)
    total = len(seq)
    search = seq + seq if circular else seq
    rc = str(Seq(primer).reverse_complement())
    tm = round(mt.Tm_NN(Seq(primer)), 1)
    gc = round(100 * (primer.count("G") + primer.count("C")) / len(primer), 1)

    out = []
    for pattern, orientation in ((primer, "forward"), (rc, "reverse")):
        for m in re.finditer(f"(?=({re.escape(pattern)}))", search):
            if m.start() >= total:
                continue
            out.append({
                "orientation": orientation,
                "start": m.start() % total + 1,
                "end": (m.start() + len(pattern) - 1) % total + 1,
                "extends": "rightwards (increasing coordinates)" if orientation == "forward"
                           else "leftwards (decreasing coordinates)",
                "tm": tm, "gc_percent": gc, "length": len(primer),
            })
    return out


def primer_report(plasmid: dict, primers: pd.DataFrame) -> pd.DataFrame:
    """Where each primer binds a plasmid, and in which orientation."""
    rows = []
    circular = plasmid["topology"] == "circular"
    for _, p in primers.iterrows():
        hits = primer_sites(plasmid["record"], p["sequence"], circular)
        if not hits:
            rows.append({"primer": p["name"], "sequence": p["sequence"],
                         "binds": False, "orientation": "", "start": None, "end": None,
                         "tm": None, "note": "no exact match in this plasmid"})
            continue
        for h in hits:
            rows.append({"primer": p["name"], "sequence": p["sequence"], "binds": True,
                         "orientation": h["orientation"], "start": h["start"],
                         "end": h["end"], "tm": h["tm"],
                         "note": str(p.get("note", ""))})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Insertion region and enzyme choice
# --------------------------------------------------------------------------- #

TAG_LABELS = ("his", "ha", "myc", "flag", "strep", "nb back", "tag")


def insertion_region(plasmid: dict, features: pd.DataFrame | None = None
                     ) -> tuple[int, int] | None:
    """Best guess at the cloning window: between the PelB leader and the first tag.

    This is a heuristic on the annotation, offered as a starting point — the user
    is expected to confirm or override it, because a vector with an unusual
    layout will fool it.
    """
    features = feature_table(plasmid) if features is None else features
    if features.empty:
        return None
    pelb = features[features["feature"].str.contains("pelb", case=False, na=False)]
    if pelb.empty:
        return None
    start = int(pelb["end"].max())

    downstream = features[features["start"] > start]
    tags = downstream[downstream["feature"].str.lower().str.contains(
        "|".join(TAG_LABELS), na=False)]
    if tags.empty:
        return None
    return start + 1, int(tags["start"].min()) - 1


def recommend_enzymes(backbone: SeqRecord, insert: str | SeqRecord | None = None,
                      region: tuple[int, int] | None = None,
                      enzymes: list[str] | None = None,
                      circular: bool = True, pad: int = 0) -> pd.DataFrame:
    """Rank enzymes for inserting a fragment into a backbone.

    An enzyme is a candidate when it cuts the backbone **exactly once**, that cut
    falls inside the cloning window, and it does **not** cut the insert. An
    enzyme that cuts the insert is still listed, flagged, because knowing that
    your favourite enzyme has a site inside this particular VHH is exactly the
    kind of thing that is worth seeing before the ligation fails.
    """
    from Bio.Seq import Seq as _Seq

    profile = restriction_profile(backbone, enzymes, circular)
    if profile.empty:
        return profile

    insert_seq = None
    if insert is not None:
        insert_seq = _Seq(re.sub(r"[^ACGTacgt]", "", str(
            insert.seq if hasattr(insert, "seq") else insert)).upper())
    insert_cuts = {}
    if insert_seq is not None and len(insert_seq):
        ins_profile = restriction_profile(SeqRecord(insert_seq), enzymes, circular=False)
        insert_cuts = dict(zip(ins_profile["enzyme"], ins_profile["n_cuts"]))

    rows = []
    for _, row in profile.iterrows():
        pos = row["cut_positions"][0]
        in_region = bool(region and region[0] - pad <= pos <= region[1] + pad)
        n_insert = int(insert_cuts.get(row["enzyme"], 0))
        usable = row["n_cuts"] == 1 and in_region and n_insert == 0
        if row["n_cuts"] == 1 and in_region:
            verdict = "usable" if n_insert == 0 else f"cuts the insert ({n_insert}x)"
        elif row["n_cuts"] == 1:
            verdict = "single cutter, outside the cloning window"
        else:
            verdict = f"cuts the backbone {row['n_cuts']}x"
        rows.append({"enzyme": row["enzyme"], "site": row["site"],
                     "backbone_cuts": row["n_cuts"], "position": pos,
                     "in_region": in_region, "insert_cuts": n_insert,
                     "overhang": row["overhang"], "usable": usable,
                     "verdict": verdict, "cut_positions": row["cut_positions"]})

    df = pd.DataFrame(rows)
    return df.sort_values(["usable", "in_region", "backbone_cuts", "enzyme"],
                          ascending=[False, False, True, True]).reset_index(drop=True)


def enzyme_pairs(candidates: pd.DataFrame, max_pairs: int = 12) -> pd.DataFrame:
    """Pairs of usable enzymes suitable for directional cloning.

    Two different enzymes cutting either side of the window give a directional
    insert. Pairs with identical overhang types are flagged, because compatible
    ends let the insert go in backwards.
    """
    from Bio.Restriction import RestrictionBatch

    usable = candidates[candidates["usable"]] if "usable" in candidates else candidates
    if len(usable) < 2:
        return pd.DataFrame(columns=["enzyme_1", "enzyme_2", "span_bp", "directional",
                                     "note"])
    batch = RestrictionBatch(list(usable["enzyme"]))
    by_name = {str(e): e for e in batch}

    rows = []
    records = usable.to_dict("records")
    for i, a in enumerate(records):
        for b in records[i + 1:]:
            ea, eb = by_name[a["enzyme"]], by_name[b["enzyme"]]
            same_overhang = (ea.ovhgseq == eb.ovhgseq) and ea.ovhgseq != ""
            span = abs(a["position"] - b["position"])
            note = ("compatible overhangs — the insert can ligate in either "
                    "orientation" if same_overhang else "")
            if a["overhang"] == "blunt" or b["overhang"] == "blunt":
                note = (note + "; " if note else "") + "blunt end — less efficient ligation"
            rows.append({"enzyme_1": a["enzyme"], "enzyme_2": b["enzyme"],
                         "span_bp": span, "directional": not same_overhang,
                         "note": note or "directional, distinct overhangs"})
    return (pd.DataFrame(rows).sort_values(["directional", "span_bp"],
                                           ascending=[False, False])
            .head(max_pairs).reset_index(drop=True))
