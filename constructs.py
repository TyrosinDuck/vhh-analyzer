"""
constructs, assemble a VHH into a fusion construct and characterise it.

Monovalent, bivalent or biparatopic, with a chosen linker and N- and
C-terminal tags, then the properties and the codon-optimised DNA of the
finished protein. The point is to answer "what will the thing I actually order
weigh, and what is its pI" before ordering it, rather than after.

**Every peptide sequence in this module must be checked against your own
construct before you order anything.** Tag sequences vary between vectors and
between papers, and a silent one-residue difference is both plausible and
expensive. The library below is deliberately limited to tags whose sequence is
unambiguous.

Two things are deliberately absent, because writing them from memory would be
worse than leaving them out:

* **SpyTag2 / SpyTag3.** Only the original SpyTag is listed. Paste the variant
  you work with as a custom piece, you have the exact sequence, this module
  does not.
* **Camelid hinges** (the llama upper long hinge and friends). Same reason: take
  it from the construct or paper you are reproducing.

Both are fully supported as custom entries, which is the honest way to handle a
sequence this module cannot vouch for.
"""

from __future__ import annotations

import pandas as pd

import vhh_core as vc

__all__ = [
    "TAGS", "LINKERS", "FORMATS", "VERIFY_NOTICE",
    "build_construct", "plot_construct", "construct_dna",
    "dissect_construct", "SPACER_KINDS", "__version__",
]

__version__ = "2026.11.19a"

VERIFY_NOTICE = (
    "Check every tag and linker sequence against your own vector or reference "
    "before ordering. Tag sequences differ between constructs, and a "
    "one-residue difference will not show up until the protein does not work."
)

# name -> (sequence, where it is usually placed, note)
TAGS: dict[str, dict] = {
    "His6": {"seq": "HHHHHH", "ends": "both",
             "note": "IMAC purification; the classic C-terminal tag on VHHs"},
    "His8": {"seq": "HHHHHHHH", "ends": "both",
             "note": "stronger IMAC binding than His6"},
    "HA": {"seq": "YPYDVPDYA", "ends": "both",
           "note": "detection tag, present in pMECS between the VHH and His6"},
    "c-Myc": {"seq": "EQKLISEEDL", "ends": "both",
              "note": "detection tag, anti-Myc 9E10"},
    "FLAG": {"seq": "DYKDDDDK", "ends": "both",
             "note": "detection and elution; note it carries four aspartates, "
                     "which will move the pI"},
    "Strep-II": {"seq": "WSHPQFEK", "ends": "both",
                 "note": "Strep-Tactin purification, mild elution"},
    "Avi": {"seq": "GLNDIFEAQKIEWHE", "ends": "both",
            "note": "BirA biotinylation site, useful for streptavidin capture "
                    "in BLI, which avoids anti-tag capture entirely"},
    "SpyTag (original)": {"seq": "AHIVMVDAYKPTK", "ends": "both",
                          "note": "the original SpyTag. For SpyTag2 or SpyTag3, "
                                  "paste your own sequence as a custom piece"},
    "TEV site": {"seq": "ENLYFQG", "ends": "both",
                 "note": "protease site, cleaves between Q and G"},
    "Thrombin site": {"seq": "LVPRGS", "ends": "both",
                      "note": "protease site, cleaves after R"},
    "Sortase A (LPETG)": {"seq": "LPETG", "ends": "C",
                          "note": "sortase-mediated ligation; must be "
                                  "C-terminal and followed by little else"},
}

LINKERS: dict[str, dict] = {
    "none (direct fusion)": {"seq": "", "note": "no linker at all"},
    "GS": {"seq": "GS", "note": "two residues, often a cloning scar"},
    "GGGGS": {"seq": "GGGGS", "note": "one GS unit, short and flexible"},
    "(GGGGS)2": {"seq": "GGGGSGGGGS", "note": "10 residues"},
    "(GGGGS)3": {"seq": "GGGGSGGGGSGGGGS",
                 "note": "15 residues, the usual default for a bivalent VHH, "
                         "long enough for both domains to bind independently"},
    "(GGGGS)4": {"seq": "GGGGSGGGGSGGGGSGGGGS",
                 "note": "20 residues, for distant epitopes"},
    "GGSGGS": {"seq": "GGSGGS", "note": "shorter flexible variant"},
    "EAAAK": {"seq": "EAAAK", "note": "rigid α-helical unit, keeps the domains apart"},
    "(EAAAK)3": {"seq": "EAAAKEAAAKEAAAK",
                 "note": "rigid helix, use when you want a fixed spacing rather "
                         "than a floppy tether"},
}

FORMATS = {
    "monovalent": {"n_vhh": 1,
                   "note": "a single VHH, with whatever tags you add"},
    "bivalent": {"n_vhh": 2,
                 "note": "the same VHH twice, raises apparent affinity through "
                         "avidity, not affinity"},
    "biparatopic": {"n_vhh": 2,
                    "note": "two different VHHs against the same target; needs "
                            "two sequences and non-overlapping epitopes"},
}


def _clean(seq: str) -> str:
    """``clean_sequence`` refuses an empty string; an absent linker or tag is
    legitimately empty, so swallow that one case rather than guarding at every
    call site."""
    seq = (seq or "").strip()
    return vc.clean_sequence(seq) if seq else ""


def _segment(label: str, kind: str, seq: str, start: int) -> dict:
    return {"label": label, "kind": kind, "sequence": seq,
            "start": start + 1, "end": start + len(seq), "length": len(seq)}


def build_construct(vhhs: list[tuple[str, str]],
                    linker: str = "",
                    n_tags: list[tuple[str, str]] | None = None,
                    c_tags: list[tuple[str, str]] | None = None,
                    tag_linker: str = "",
                    ph: float | None = None) -> dict:
    """Assemble a construct and characterise the finished protein.

    ``vhhs`` is a list of ``(name, sequence)``, one for a monovalent construct,
    two for a bivalent or biparatopic one. ``n_tags`` and ``c_tags`` are lists of
    ``(label, sequence)`` applied in the order given, reading N to C.
    ``tag_linker`` is inserted between the VHH body and each tag; ``linker``
    between two VHHs.

    Returns the sequence, a segment map, the characterisation of the whole
    construct and of each VHH alone, and the warnings worth seeing.
    """
    n_tags = list(n_tags or [])
    c_tags = list(c_tags or [])
    pieces: list[tuple[str, str, str]] = []      # (label, kind, sequence)

    for label, seq in n_tags:
        pieces.append((label, "tag", _clean(seq)))
        if tag_linker:
            pieces.append(("linker", "linker", _clean(tag_linker)))

    for i, (name, seq) in enumerate(vhhs):
        if i:
            pieces.append(("linker", "linker", _clean(linker)))
        pieces.append((name, "vhh", _clean(seq)))

    for label, seq in c_tags:
        if tag_linker:
            pieces.append(("linker", "linker", _clean(tag_linker)))
        pieces.append((label, "tag", _clean(seq)))

    pieces = [p for p in pieces if p[2]]

    segments, cursor, parts = [], 0, []
    for label, kind, seq in pieces:
        segments.append(_segment(label, kind, seq, cursor))
        cursor += len(seq)
        parts.append(seq)
    sequence = "".join(parts)

    warnings = [VERIFY_NOTICE]
    if len(vhhs) == 2 and not linker:
        warnings.append(
            "Two VHHs fused with no linker: the domains will interfere with each "
            "other, and this is almost never what you want. (GGGGS)3 is the usual "
            "starting point.")
    if any(k == "tag" for _l, k, _s in pieces[:1]) and not tag_linker:
        warnings.append(
            "An N-terminal tag sits directly on the VHH's first framework "
            "residue; check that your vector's signal peptide cleavage still "
            "leaves the construct you expect.")
    sortase = [l for l, k, _s in pieces if l.startswith("Sortase")]
    if sortase and pieces and pieces[-1][0] != sortase[-1]:
        warnings.append(
            "The sortase motif is not the last element, LPETG is cleaved "
            "between T and G, so anything after it is lost.")
    if len(sequence) > vc.DEFAULTS["max_length"]:
        warnings.append(
            f"The construct is {len(sequence)} residues, above the "
            f"{vc.DEFAULTS['max_length']}-residue limit the analysis tools use; "
            "the characterisation below is still valid but the annotation tools "
            "expect a single VHH domain.")

    whole = pd.DataFrame([vc.characterize("construct", sequence, ph=ph)])
    per_vhh = pd.DataFrame([vc.characterize(n, _clean(s), ph=ph)
                            for n, s in vhhs])

    return {"sequence": sequence, "length": len(sequence),
            "segments": pd.DataFrame(segments),
            "properties": whole, "per_vhh": per_vhh,
            "warnings": warnings,
            "composition": {k: sum(s["length"] for s in segments if s["kind"] == k)
                            for k in ("vhh", "linker", "tag")}}


def construct_dna(construct: dict, avoid: set[str] | None = None,
                  add_stop: bool = True) -> dict:
    """Codon-optimised DNA for the finished construct, with a round trip check."""
    sites = vc.DEFAULT_AVOID_SITES if avoid is None else {
        k: v for k, v in vc.DEFAULT_AVOID_SITES.items() if k in avoid}
    dna = vc.reverse_translate(construct["sequence"], avoid_sites=sites,
                               add_stop=add_stop)
    back = vc.translate_dna(dna["dna"])
    dna["back_translation_ok"] = back.rstrip("*") == construct["sequence"]
    return dna



# --------------------------------------------------------------------------- #
# The other direction: a construct in, its pieces out
# --------------------------------------------------------------------------- #
#
# `build_construct` goes from pieces to a sequence. Working from a construct
# that already exists, ordered or received, the question is the reverse one:
# given one long chain, which parts of it are VHHs, what joins them, and what
# is stuck on the ends. A bivalent or biparatopic VHH submitted whole used to
# be classified as "contains a V domain" and dropped out of every VHH-specific
# analysis, which is correct for the chain as a whole and useless in practice.
#
# The cut is not a guess. `vhh_core.find_domains` returns the bounds of every
# variable domain ANARCI numbers, so the linker is what is left between two of
# them, by subtraction. Measured on caplacizumab fused to cAbLys3 through five
# different linkers, from a direct fusion to (EAAAK)3: both domains found every
# time, and the recovered inter-domain stretch equal to the linker that went in,
# residue for residue. Without ANARCI the regex scanner finds the same domains
# on a canonical scaffold but places their N-terminal boundary by convention,
# which moves a few residues of framework into the linker; `bounds_quality`
# says which of the two produced the segments below.

SPACER_KINDS = ("linker", "tag", "other")


def _known_pieces() -> list[tuple[str, str, str]]:
    """(sequence, label, kind) for every tag and linker this module knows.

    Derived from TAGS and LINKERS rather than written out again, so a tag added
    above is recognised here without a second edit. Sorted longest first so
    that His8 is matched before His6 finds its first six residues.
    """
    pieces = [(d["seq"], name, "tag") for name, d in TAGS.items() if d["seq"]]
    pieces += [(d["seq"], name, "linker") for name, d in LINKERS.items() if d["seq"]]
    return sorted(pieces, key=lambda p: -len(p[0]))


def _label_spacer(seq: str, default_kind: str) -> list[tuple[str, str, str]]:
    """Cut a stretch between domains into recognised tags and the rest.

    Returns (label, kind, sequence) triples covering `seq` end to end. A stretch
    that matches nothing comes back as one "unassigned" piece: an unrecognised
    sequence is reported as unrecognised rather than called a linker because it
    happens to sit where a linker would.
    """
    if not seq:
        return []
    for piece, label, kind in _known_pieces():
        idx = seq.find(piece)
        if idx == -1:
            continue
        out = []
        out += _label_spacer(seq[:idx], default_kind)
        out.append((label, kind, piece))
        out += _label_spacer(seq[idx + len(piece):], default_kind)
        return out
    return [("unassigned", default_kind, seq)]


def dissect_construct(seq: str, name: str = "construct",
                      method: str | None = None,
                      ph: float | None = None) -> dict:
    """Cut a multi-domain construct into its domains, linkers and tags.

    The inverse of `build_construct`, and it returns the same shape of dict, so
    `plot_construct`, `construct_dna` and the segment table work on a dissected
    construct exactly as they do on an assembled one.

    Beyond that shape it adds ``domains``, one entry per variable domain with
    its name, sequence and bounds, ``bounds_quality`` ("exact" from ANARCI,
    "approximate" from the regex scanner), and ``junction_sites``, the
    predicted protease sites that fall outside the domains, which is where a
    linker is cut in practice and where the per-domain analyses would never
    look.

    ``per_vhh`` characterises each domain on its own. The whole chain is
    characterised too, in ``properties``: mass, pI and extinction coefficient
    are additive and stay meaningful across a fusion. Nothing camelid is
    computed on the whole chain, because a hallmark tetrad or a CDR3 length
    read across two domains would be an artefact of the concatenation.
    """
    seq = vc.clean_sequence(seq)
    domains = vc.find_domains(seq, method=method)

    quality = ({d.get("bounds", "exact") for d in domains} or {"exact"})
    bounds_quality = "approximate" if "approximate" in quality else "exact"

    pieces: list[tuple[str, str, str]] = []      # (label, kind, sequence)
    cursor = 0
    named_domains = []
    for i, d in enumerate(domains, start=1):
        if d["start"] > cursor:
            # Before the first domain this is an N-terminal extension; between
            # two domains it is the linker. The default kind differs so that a
            # stretch matching nothing is not called a linker when it sits on
            # the end of the chain.
            pieces += _label_spacer(seq[cursor:d["start"]],
                                    "linker" if i > 1 else "other")
        label = f"{name}_D{i}"
        body = seq[d["start"]:d["end"] + 1]
        pieces.append((label, "vhh", body))
        named_domains.append({"name": label, "sequence": body,
                              "start": d["start"] + 1, "end": d["end"] + 1,
                              "length": len(body), "chain": d.get("chain"),
                              "evalue": d.get("evalue"),
                              "bounds": d.get("bounds", "exact")})
        cursor = d["end"] + 1
    if cursor < len(seq):
        pieces += _label_spacer(seq[cursor:], "other")

    pieces = [p for p in pieces if p[2]]
    segments, pos = [], 0
    for label, kind, piece in pieces:
        segments.append(_segment(label, kind, piece, pos))
        pos += len(piece)
    segments_df = pd.DataFrame(segments)

    # Protease sites outside the domains. Computed on the whole chain so that a
    # motif straddling a junction is seen, then kept only where it falls in a
    # segment that is not a domain.
    sites = vc.find_protease_sites(seq)
    junction = sites.iloc[0:0] if sites.empty else None
    if not sites.empty:
        spans = [(s["start"], s["end"], s["label"], s["kind"])
                 for s in segments if s["kind"] != "vhh"]
        keep, seg_label = [], []
        for _idx, row in sites.iterrows():
            hit = next((lab for a, b, lab, _k in spans if a <= row["position"] <= b),
                       None)
            keep.append(hit is not None)
            seg_label.append(hit)
        junction = sites[pd.Series(keep, index=sites.index)].copy()
        junction["segment"] = [s for s in seg_label if s is not None]

    warnings = []
    if not domains:
        warnings.append(
            "No immunoglobulin variable domain was found, so there is nothing "
            "to cut: the sequence is reported as one piece.")
    if bounds_quality == "approximate":
        warnings.append(
            "Domain boundaries come from the regex scanner, not ANARCI. The "
            "C-terminal end of each domain is the FR4 motif and is located "
            "directly; the N-terminal end is placed by IMGT convention, so the "
            "first residues of a domain can be reported as part of the linker. "
            "Read the linker length as approximate.")
    unassigned = [s for s in segments if s["label"] == "unassigned"]
    if unassigned:
        warnings.append(
            f"{len(unassigned)} stretch(es) between or around the domains match "
            "no tag or linker this module knows, and are left as "
            "'unassigned' rather than guessed at.")
    if len(domains) > 1:
        warnings.append(
            "Mass, pI and extinction coefficient below are for the whole chain "
            "and stay valid across a fusion. Everything camelid, the hallmark "
            "tetrad, CDR lengths, germline and humanness, is reported per "
            "domain only; read across the junction it would describe an "
            "artefact of the concatenation rather than the protein.")
    if len(seq) > vc.DEFAULTS["max_length"]:
        warnings.append(
            f"The chain is {len(seq)} residues, above the "
            f"{vc.DEFAULTS['max_length']}-residue limit the analysis tools use.")

    whole = pd.DataFrame([vc.characterize(name, seq, ph=ph)])
    per_vhh = (pd.DataFrame([vc.characterize(d["name"], d["sequence"], ph=ph)
                             for d in named_domains])
               if named_domains else pd.DataFrame())

    return {"sequence": seq, "length": len(seq), "name": name,
            "segments": segments_df,
            "properties": whole, "per_vhh": per_vhh,
            "domains": named_domains,
            "n_domains": len(named_domains),
            "bounds_quality": bounds_quality,
            "junction_sites": junction,
            "warnings": warnings,
            "composition": {k: sum(s["length"] for s in segments
                                   if s["kind"] == k)
                            for k in ("vhh", "linker", "tag", "other")}}


KIND_COLORS = {"vhh": "#4c72b0", "linker": "#c7c7c7", "tag": "#dd8452",
               "other": "#8c8c8c"}


def plot_construct(construct: dict, title: str = "Construct"):
    """Schematic of the construct, drawn to scale."""
    import matplotlib.pyplot as plt

    segments = construct["segments"]
    fig, ax = plt.subplots(figsize=(10, 1.9))
    for _, s in segments.iterrows():
        colour = KIND_COLORS.get(s["kind"], "0.7")
        ax.barh(0, s["length"], left=s["start"] - 1, height=0.5,
                color=colour, edgecolor="white", lw=1.2)
        if s["length"] >= max(4, 0.03 * construct["length"]):
            ax.text(s["start"] - 1 + s["length"] / 2, 0, s["label"],
                    ha="center", va="center", fontsize=8,
                    color="white" if s["kind"] != "linker" else "0.25")
        else:
            ax.annotate(s["label"], xy=(s["start"] - 1 + s["length"] / 2, 0.27),
                        xytext=(0, 14), textcoords="offset points",
                        ha="center", fontsize=7.5, rotation=0,
                        arrowprops=dict(arrowstyle="-", lw=0.7, color="0.5"))
    ax.set_xlim(0, construct["length"])
    ax.set_ylim(-0.6, 0.9)
    ax.set_yticks([])
    ax.set_xlabel("Residue")
    ax.set_title(f"{title}, {construct['length']} residues", fontsize=10)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig
