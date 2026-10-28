"""
cloning, plasmid handling for in-silico restriction cloning of VHHs.

Phase 1: load a plasmid (SnapGene .dna, GenBank or FASTA), inventory what is in
it, and draw its map.

Two sources of truth, deliberately kept apart:

* **Annotated files** (SnapGene, GenBank) already carry their features. They are
  read as they are, nothing is guessed, and that is by far the reliable path.
* **Bare FASTA** carries only sequence, so features have to be found. This
  module detects *short, unambiguous protein motifs*, the PelB leader, His6,
  HA, Myc, a VHH scaffold, by translating in all six frames. It deliberately
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

# Bumped whenever the public API changes, so app.py can tell the user that a
# module is out of date instead of crashing on a missing attribute.
__version__ = "2026.10.21b"

# --------------------------------------------------------------------------- #
# Motif library, short protein elements only, see the module docstring
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
        warning = ("FASTA carries no features, only short protein motifs can be "
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
    module inferred, they do not carry the same confidence.
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
        if end <= start:                      # wraps the origin, clip for drawing
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
    ax.set_title(f"{plasmid['name']}, {total} bp, {plasmid['topology']}", fontsize=11)
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
    ax.set_title(f"{plasmid['name']}, cloning region {start + 1}–{end}", fontsize=10)
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

# Labels that mark the start of the tag region, i.e. the 3' end of the cloning
# window. "Nb back" is deliberately absent: in pHEN6 it names a pseudo-nanobody
# fragment used as a cloning stuffer, not a tag, and treating it as one cuts the
# window short and hides the enzymes that sit just after it.
TAG_LABELS = ("his", "ha", "myc", "flag", "strep", "tag")


def insertion_region(plasmid: dict, features: pd.DataFrame | None = None
                     ) -> tuple[int, int] | None:
    """Best guess at the cloning window: between the PelB leader and the first tag.

    This is a heuristic on the annotation, offered as a starting point, the user
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
            note = ("compatible overhangs, the insert can ligate in either "
                    "orientation" if same_overhang else "")
            if a["overhang"] == "blunt" or b["overhang"] == "blunt":
                note = (note + "; " if note else "") + "blunt end, less efficient ligation"
            rows.append({"enzyme_1": a["enzyme"], "enzyme_2": b["enzyme"],
                         "span_bp": span, "directional": not same_overhang,
                         "note": note or "directional, distinct overhangs"})
    return (pd.DataFrame(rows).sort_values(["directional", "span_bp"],
                                           ascending=[False, False])
            .head(max_pairs).reset_index(drop=True))


# --------------------------------------------------------------------------- #
# PCR, digestion and assembly
# --------------------------------------------------------------------------- #

def _anneal(template: str, primer: str, min_anneal: int = 15) -> list[dict]:
    """Find where a primer's 3' end anneals, allowing a 5' tail.

    Only the 3' portion of a primer has to match: the 5' end may carry a tail
    (a restriction site, for instance) that is absent from the template and ends
    up in the product. The longest matching 3' suffix is used.
    """
    primer = re.sub(r"[^ACGT]", "", primer.upper())
    hits = []
    for length in range(len(primer), min_anneal - 1, -1):
        core = primer[len(primer) - length:]
        tail = primer[:len(primer) - length]
        for m in re.finditer(f"(?=({re.escape(core)}))", template):
            hits.append({"start": m.start(), "end": m.start() + length - 1,
                         "anneal": core, "tail": tail, "anneal_len": length})
        if hits:                       # longest match wins; stop at the first length
            break
    return hits


def simulate_pcr(record: SeqRecord, forward: str, reverse: str,
                 circular: bool = True, min_anneal: int = 15,
                 max_product: int = 20000) -> dict:
    """Amplify a template with two primers and return the product.

    The forward primer anneals to the minus strand and extends rightwards; the
    reverse primer anneals to the plus strand and extends leftwards, so what is
    searched for is its reverse complement. Any 5' tail on either primer is
    carried into the product, which is how restriction sites get added.
    """
    seq = str(record.seq).upper()
    total = len(seq)
    search = seq + seq if circular else seq

    fwd = re.sub(r"[^ACGT]", "", forward.upper())
    rev = re.sub(r"[^ACGT]", "", reverse.upper())
    rev_rc = str(Seq(rev).reverse_complement())

    f_hits = _anneal(search, fwd, min_anneal)
    r_hits = _anneal(search, rev_rc, min_anneal)
    # a reverse primer's tail sits at the 3' end of its reverse complement
    for h in r_hits:
        tail = re.sub(r"[^ACGT]", "", h["tail"])
        h["tail"] = ""
        h["tail_rc"] = str(Seq(rev[:len(rev) - h["anneal_len"]]).reverse_complement())

    result = {"product": "", "length": 0, "forward_site": None, "reverse_site": None,
              "forward_tail": "", "reverse_tail": "", "warnings": [], "n_products": 0}

    if not f_hits:
        result["warnings"].append(
            f"forward primer does not anneal (needs {min_anneal} nt of 3' match)")
    if not r_hits:
        result["warnings"].append(
            f"reverse primer does not anneal (needs {min_anneal} nt of 3' match)")
    if not f_hits or not r_hits:
        return result

    products, seen = [], set()
    for f in f_hits:
        if f["start"] >= total:              # duplicate from the doubled copy
            continue
        for r in r_hits:
            if r["end"] <= f["end"]:
                continue
            size = r["end"] - f["start"] + 1
            if size > max_product:
                continue
            # on a circular template the doubled sequence yields the same product
            # twice; key on the wrapped coordinates to count it once
            key = (f["start"] % total, r["end"] % total)
            if key in seen:
                continue
            seen.add(key)
            products.append((f, r, size))
    if not products:
        result["warnings"].append(
            "primers anneal but never in a productive orientation, check which one "
            "is forward and which is reverse")
        return result

    products.sort(key=lambda t: t[2])
    f, r, size = products[0]
    result["n_products"] = len(products)
    if len(products) > 1:
        result["warnings"].append(
            f"{len(products)} possible products, the shortest ({size} bp) is reported")

    core = search[f["start"]: r["end"] + 1]
    product = f["tail"] + core + r["tail_rc"]
    result.update({
        "product": product, "length": len(product),
        "forward_site": (f["start"] % total + 1, f["end"] % total + 1),
        "reverse_site": (r["start"] % total + 1, r["end"] % total + 1),
        "forward_tail": f["tail"], "reverse_tail": r["tail_rc"],
        "forward_anneal_len": f["anneal_len"], "reverse_anneal_len": r["anneal_len"],
    })
    return result


def cut_positions(seq: str, enzyme_name: str, circular: bool = False) -> list[int]:
    """1-based top-strand cut positions of one enzyme in a sequence."""
    from Bio.Restriction import RestrictionBatch

    batch = RestrictionBatch([enzyme_name])
    enzyme = list(batch)[0]
    return sorted(enzyme.search(Seq(re.sub(r"[^ACGT]", "", seq.upper())),
                                linear=not circular))


def assemble(backbone: SeqRecord, insert_dna: str, enzyme_5: str, enzyme_3: str,
             circular: bool = True) -> dict:
    """Cut backbone and insert with two enzymes and join them.

    ``enzyme_5`` is the enzyme on the 5' side of the cloning window (the lower
    coordinate), ``enzyme_3`` the one on the 3' side. Both must cut the backbone
    exactly once and the insert exactly once, which is what makes the join
    unambiguous and directional.

    The assembly is done on the top strand. Because both partners are cut with
    the same pair of enzymes, their overhangs match by construction, so the
    top-strand join is the product you would get at the bench, but it is a
    design check, not a simulation of ligation efficiency.
    """
    vector = str(backbone.seq).upper()
    insert = re.sub(r"[^ACGT]", "", insert_dna.upper())
    report = {"construct": "", "length": 0, "warnings": [], "ok": False,
              "vector_cuts": {}, "insert_cuts": {}}

    for name, seq_, label, circ in ((enzyme_5, vector, "backbone", circular),
                                    (enzyme_3, vector, "backbone", circular)):
        pass  # positions collected below, kept explicit for clarity

    v5 = cut_positions(vector, enzyme_5, circular)
    v3 = cut_positions(vector, enzyme_3, circular)
    i5 = cut_positions(insert, enzyme_5, False)
    i3 = cut_positions(insert, enzyme_3, False)
    report["vector_cuts"] = {enzyme_5: v5, enzyme_3: v3}
    report["insert_cuts"] = {enzyme_5: i5, enzyme_3: i3}

    for label, name, cuts in (("backbone", enzyme_5, v5), ("backbone", enzyme_3, v3),
                              ("insert", enzyme_5, i5), ("insert", enzyme_3, i3)):
        if len(cuts) != 1:
            report["warnings"].append(
                f"{name} cuts the {label} {len(cuts)} time(s), exactly one is needed")
    if report["warnings"]:
        return report

    p5, p3 = v5[0], v3[0]
    q5, q3 = i5[0], i3[0]
    if p5 > p3:
        report["warnings"].append(
            f"{enzyme_5} cuts the backbone at {p5}, downstream of {enzyme_3} at {p3}, "
            "swap the two enzymes")
        return report
    if q5 > q3:
        report["warnings"].append(
            f"{enzyme_5} cuts the insert at {q5}, downstream of {enzyme_3} at {q3}, "
            "the insert would go in backwards")
        return report

    removed = p3 - p5
    piece = insert[q5 - 1: q3 - 1]
    construct = vector[: p5 - 1] + piece + vector[p3 - 1:]

    report.update({
        "construct": construct, "length": len(construct), "ok": True,
        "removed_bp": removed, "inserted_bp": len(piece),
        "insert_start": p5, "insert_end": p5 + len(piece) - 1,
        "net_change": len(piece) - removed,
    })
    # No frame arithmetic here on purpose. An empty phage-display vector carries a
    # stuffer, so its downstream tag is usually *out* of frame with the leader
    # until a correctly sized VHH is cloned in, comparing the length change to
    # the empty vector's frame would flag every correct construct. The real test
    # is translating the finished construct, which verify_construct does.
    return report


def verify_construct(construct: str, vector_layout: str = "pMES3") -> dict:
    """Translate the assembled construct and check the expression cassette.

    Reuses the Sanger cassette extractor: the same code that validates a
    sequencing read validates a designed construct, so a frame error, a
    premature stop or a lost tag is caught here rather than after the ligation.
    """
    from sanger import VECTORS, extract_cassette

    layout = vector_layout if vector_layout in VECTORS else "pMES3"
    ex = extract_cassette(construct, layout)
    ok = bool(ex["vhh_aa"]) and not ex["frameshift"] and "*" not in ex["vhh_aa"]
    return {
        "ok": ok, "vhh_aa": ex["vhh_aa"], "vhh_len": len(ex["vhh_aa"]),
        "cassette_aa": ex["cassette_aa"], "tags_found": ex["tags_found"],
        "tags_missing": ex["tags_missing"],
        "premature_stops": ex["n_premature_stops"],
        "stop_codon": ex["stop_codon_present"], "frameshift": ex["frameshift"],
        "strand": ex["strand"], "frame": ex["frame"],
        "warnings": list(ex["warnings"]), "layout": layout,
    }


def evaluate_pairs(backbone: SeqRecord, insert_dna: str, candidates: pd.DataFrame,
                   layout: str = "pMES3", circular: bool = True,
                   max_pairs: int = 20) -> pd.DataFrame:
    """Try each enzyme pair for real: assemble, translate, and report the verdict.

    Ranking enzymes by where they cut says which pairs are *plausible*. This
    goes one step further and builds each construct, then reads its expression
    cassette back out, so a pair that cuts in the right places but leaves the
    VHH out of frame with the tags is caught here rather than at the bench.
    """
    pool = (candidates[candidates["backbone_cuts"] == 1]
            if "backbone_cuts" in candidates else candidates)
    names = list(dict.fromkeys(pool["enzyme"])) if len(pool) else []
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            for e5, e3 in ((a, b), (b, a)):
                built = assemble(backbone, insert_dna, e5, e3, circular=circular)
                if not built["ok"]:
                    continue
                check = verify_construct(built["construct"], layout)
                rows.append({
                    "enzyme_5": e5, "enzyme_3": e3,
                    "construct_bp": built["length"],
                    "removed_bp": built["removed_bp"],
                    "inserted_bp": built["inserted_bp"],
                    "cassette_ok": check["ok"],
                    "vhh_len": check["vhh_len"],
                    "tags_found": check["tags_found"],
                    "premature_stops": check["premature_stops"],
                    "verdict": ("full cassette in frame" if check["ok"] else
                                "; ".join(check["warnings"]) or "cassette not readable"),
                    "construct": built["construct"],
                })
                break           # the other orientation of the same pair is redundant
    if not rows:
        return pd.DataFrame(columns=["enzyme_5", "enzyme_3", "construct_bp",
                                     "cassette_ok", "verdict"])
    df = pd.DataFrame(rows)
    return (df.sort_values(["cassette_ok", "vhh_len"], ascending=[False, False])
            .head(max_pairs).reset_index(drop=True))


def construct_record(construct: str, name: str, insert_start: int, insert_end: int,
                     layout: str = "pMES3", circular: bool = True) -> SeqRecord:
    """Wrap an assembled construct as an annotated SeqRecord, ready for GenBank.

    The insert is annotated from the assembly coordinates; the cassette parts
    (leader, VHH, tags) are re-detected on the finished sequence, so the record
    describes what the construct actually encodes rather than what was intended.
    """
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    record = SeqRecord(Seq(construct), id=name[:16] or "construct", name=name[:16],
                       description=f"{name}, assembled in silico",
                       annotations={"molecule_type": "DNA",
                                    "topology": "circular" if circular else "linear"})
    record.features.append(SeqFeature(
        SimpleLocation(insert_start - 1, insert_end), type="misc_feature",
        qualifiers={"label": ["inserted fragment"]}))

    for hit in detect_motifs(record, circular=circular):
        record.features.append(SeqFeature(
            SimpleLocation(hit["start"] - 1, hit["end"],
                           strand=-1 if hit["strand"] == "-" else 1),
            type="CDS" if hit["kind"] == "vhh" else "misc_feature",
            qualifiers={"label": [hit["feature"]]}))
    return record


def record_to_genbank(record: SeqRecord) -> str:
    """GenBank text for a record, the format SnapGene and Benchling both read."""
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


# --------------------------------------------------------------------------- #
# One-call cloning planner
# --------------------------------------------------------------------------- #

def extract_vhh_protein(dna: str, circular: bool = False) -> dict:
    """Find the VHH domain encoded by a piece of DNA, plasmid or bare fragment.

    Works without a leader: the domain is recognised by its own FR1/FR2/FR4
    signatures in any of the six frames. A bare PCR product usually starts at the
    5' restriction site, a few codons into FR1, so the protein returned can be
    slightly shorter at the N-terminus than the domain in a full plasmid, which
    is exactly what ``compare_vhh`` accounts for.
    """
    record = SeqRecord(Seq(re.sub(r"[^ACGTacgt]", "", dna.upper())),
                       id="x", annotations={"molecule_type": "DNA"})
    if len(record.seq) < 60:
        return {"protein": "", "start": None, "end": None, "strand": None,
                "note": "sequence too short to contain a VHH"}
    hits = [h for h in detect_motifs(record, circular=circular) if h["kind"] == "vhh"]
    if not hits:
        return {"protein": "", "start": None, "end": None, "strand": None,
                "note": "no VHH scaffold found (FR1/FR2/FR4 signatures absent)"}
    best = max(hits, key=lambda h: len(h["protein"]))
    return {"protein": best["protein"], "start": best["start"], "end": best["end"],
            "strand": best["strand"],
            "note": f"{len(hits)} scaffold hit(s); longest kept"}


def compare_vhh(donor_protein: str, construct_protein: str) -> dict:
    """Compare the VHH before and after cloning, and say precisely how they differ.

    Three outcomes:

    * **identical**: the construct encodes exactly the donor's domain;
    * **donor part preserved**: the construct carries extra N-terminal residues
      that the recipient vector contributes. A bare PCR product starts at the 5'
      restriction site, so the first few codons of FR1 necessarily come from the
      recipient; as long as the construct ends with the donor's sequence and
      starts on a valid FR1, nothing of the donor was altered;
    * **changed**: a real difference, reported position by position.
    """
    d, c = donor_protein or "", construct_protein or ""
    if not d or not c:
        return {"status": "unknown", "ok": False, "differences": [],
                "message": "one of the two sequences could not be read"}
    if d == c:
        return {"status": "identical", "ok": True, "differences": [],
                "from_recipient": "",
                "message": "the construct encodes exactly the donor's VHH"}
    if c.endswith(d):
        prefix = c[: len(c) - len(d)]
        valid_fr1 = bool(VHH_FR1_RE.match(c[:5]))
        return {
            "status": "donor part preserved", "ok": valid_fr1, "differences": [],
            "from_recipient": prefix,
            "message": (
                f"the construct adds {len(prefix)} N-terminal residue(s) "
                f"({prefix}) that come from the recipient vector, upstream of the "
                "5' restriction site; everything the donor encodes is unchanged"
                + ("" if valid_fr1 else
                   ", but the result does not start on a recognisable FR1, so "
                   "check for a spacer between the leader and the domain")),
        }
    diffs = [{"position": i + 1, "donor": a, "construct": b}
             for i, (a, b) in enumerate(zip(d, c)) if a != b][:25]
    return {"status": "changed", "ok": False, "differences": diffs,
            "from_recipient": "",
            "message": (f"the VHH differs: {len(d)} aa in the donor against {len(c)} "
                        f"aa in the construct, {len(diffs)} substitution(s) listed")}


def auto_primer_pair(donor: SeqRecord, primers: pd.DataFrame,
                     vhh_start: int, vhh_end: int, circular: bool = True) -> dict:
    """Pick the primer pair that brackets the VHH on a donor plasmid.

    Each primer is tried as written and as its reverse complement, because a
    primer is often noted as the plus-strand sequence of its target rather than
    as the oligo that primes back towards the insert. The shortest bracketing
    product wins.
    """
    candidates = []
    for _, row in primers.iterrows():
        for seq_, form in ((row["sequence"], "as written"),
                           (str(Seq(row["sequence"]).reverse_complement()),
                            "reverse complement")):
            for hit in primer_sites(donor, seq_, circular):
                candidates.append({"name": row["name"], "form": form,
                                   "sequence": seq_, **hit})

    forwards = [c for c in candidates
                if c["orientation"] == "forward" and c["end"] < vhh_start]
    reverses = [c for c in candidates
                if c["orientation"] == "reverse" and c["start"] > vhh_end]
    best = None
    for f in forwards:
        for r in reverses:
            size = r["end"] - f["start"] + 1
            if size <= 0:
                continue
            if best is None or size < best["length"]:
                best = {"forward": f, "reverse": r, "length": size}
    return {"pair": best, "candidates": candidates,
            "n_forward": len(forwards), "n_reverse": len(reverses)}


DEFAULT_ENZYMES = ("PstI", "NotI", "BstEII", "NcoI", "SfiI", "XbaI", "EcoRI",
                   "HindIII", "BamHI", "XhoI", "SalI", "NdeI", "KpnI", "SacI",
                   "SpeI", "ApaLI", "AscI", "AgeI", "NheI", "EagI", "ApaI", "SmaI")


def plan_cloning(recipient: dict, donor_dna: str, donor_is_plasmid: bool,
                 primers: pd.DataFrame | None = None,
                 enzymes: tuple[str, ...] | None = None,
                 layout: str = "pMES3") -> dict:
    """Plan the whole cloning in one call.

    Reads the VHH out of the donor, amplifies it when the donor is a plasmid,
    tries every enzyme pair that could open the recipient between its leader and
    its tags, and keeps only the pairs whose construct encodes the donor's VHH
    unchanged. Returns everything needed to show the result and the reason for
    every rejection.
    """
    enzymes = enzymes or DEFAULT_ENZYMES
    primers = cl_primers = (primers if primers is not None else load_primers())
    donor_dna = re.sub(r"[^ACGTacgt]", "", (donor_dna or "").upper())

    out = {"ok": False, "stage": "donor", "messages": [], "donor_vhh": "",
           "amplicon": "", "primer_pair": None, "candidates": pd.DataFrame(),
           "plans": [], "region": None}

    if not donor_dna:
        out["messages"].append("no donor DNA provided")
        return out

    # ---- 1. the VHH in the donor ----
    found = extract_vhh_protein(donor_dna, circular=donor_is_plasmid)
    if not found["protein"]:
        out["messages"].append(f"donor: {found['note']}")
        return out
    out["donor_vhh"] = found["protein"]
    out["donor_vhh_position"] = (found["start"], found["end"])

    # ---- 2. the fragment to clone ----
    insert = donor_dna
    if donor_is_plasmid:
        record = SeqRecord(Seq(donor_dna), id="donor",
                           annotations={"molecule_type": "DNA",
                                        "topology": "circular"})
        chosen = auto_primer_pair(record, primers, found["start"], found["end"])
        out["primer_candidates"] = chosen["candidates"]
        if not chosen["pair"]:
            out["stage"] = "primers"
            out["messages"].append(
                "no primer pair brackets the VHH in this donor, "
                f"{chosen['n_forward']} usable forward and {chosen['n_reverse']} "
                "usable reverse site(s) found")
            return out
        pair = chosen["pair"]
        pcr = simulate_pcr(record, pair["forward"]["sequence"],
                           pair["reverse"]["sequence"], circular=True)
        if not pcr["product"]:
            out["stage"] = "pcr"
            out["messages"].extend(pcr["warnings"] or ["the PCR gave no product"])
            return out
        insert = pcr["product"]
        out["primer_pair"] = pair
        out["pcr"] = pcr
    out["amplicon"] = insert

    # ---- 3. where the insert must land ----
    features = feature_table(recipient)
    region = insertion_region(recipient, features)
    if not region:
        out["stage"] = "region"
        out["messages"].append(
            "could not locate the window between the leader and the first tag in "
            "the recipient")
        return out
    out["region"] = region

    # ---- 4. enzymes, assembly, and the protein check ----
    circular = recipient["topology"] == "circular"
    cand = recommend_enzymes(recipient["record"], insert, region,
                             enzymes=list(enzymes), circular=circular)
    out["candidates"] = cand
    # A fragment can only be cloned by restriction if at least two enzymes cut
    # both it and the recipient inside the window. Saying which of the two is
    # missing is what tells the user whether to change primers or change vector.
    cuts_insert = cand[cand["insert_cuts"] > 0]
    opens_vector = cand[(cand["backbone_cuts"] == 1) & cand["in_region"]]
    shared = set(cuts_insert["enzyme"]) & set(opens_vector["enzyme"])
    if len(shared) < 2:
        out["stage"] = "sites"
        if cuts_insert.empty:
            out["messages"].append(
                "the donor fragment carries none of the restriction sites "
                "considered, a fragment amplified with plain primers cannot be "
                "cloned this way, it needs primers whose 5' tails add the sites")
        else:
            out["messages"].append(
                "the donor fragment and the recipient share fewer than two usable "
                "sites (shared: " + (", ".join(sorted(shared)) or "none") + "), "
                "amplify with primers carrying tails that add sites the recipient "
                "also has, or pick another recipient")
        return out

    pairs = evaluate_pairs(recipient["record"], insert, cand, layout=layout,
                           circular=circular)
    plans = []
    for _, row in pairs.iterrows():
        after = verify_construct(row["construct"], layout)
        comparison = compare_vhh(found["protein"], after["vhh_aa"])
        p5 = cand.loc[cand["enzyme"] == row["enzyme_5"], "position"].iloc[0]
        p3 = cand.loc[cand["enzyme"] == row["enzyme_3"], "position"].iloc[0]
        # A pair that cuts far apart also deletes whatever the recipient has
        # between the two sites. The cassette can still translate while the
        # promoter or the ribosome binding site is gone, so name what is lost.
        lost = [str(f["feature"]) for _, f in features.iterrows()
                if f["start"] >= p5 and f["end"] <= p3
                and not str(f["feature"]).lower().startswith(("nb ", "insertion"))]
        plans.append({**row.to_dict(), "construct_vhh": after["vhh_aa"],
                      "comparison": comparison,
                      "protein_ok": comparison["ok"],
                      "deletes": lost,
                      "usable": bool(row["cassette_ok"] and comparison["ok"]
                                     and not lost)})
    # Rank by: works at all, then an exact protein match, then the pair that
    # disturbs the recipient least, which is what a cloner would pick by hand.
    plans.sort(key=lambda p: (not p["usable"],
                              p["comparison"]["status"] != "identical",
                              p["removed_bp"]))
    out["plans"] = plans
    out["ok"] = any(p["usable"] for p in plans)
    out["stage"] = "done"
    if not out["ok"]:
        out["messages"].append(
            "no enzyme pair gives a construct that preserves the donor's VHH")
    return out
