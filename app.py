"""
VHH tool box — Streamlit front-end.

The app is organised in three sections, reachable from the home page and from
the sidebar:

* **Sequence analysis** — characterization, CDR annotation, alignment, DNA
* **Sanger sequencing** — extract VHHs from .ab1 chromatograms, then send the
  ones you pick to the analysis section
* **In-silico cloning** — pick a backbone, an insert and get the enzymes and
  primers that fit

Data handling: submitted sequences are never written to a persistent file,
logged or sent to a third party. They are processed in the server's memory. The
result cache is deliberately short-lived (see CACHE_TTL) and the alignment step
writes a temporary FASTA that is deleted as soon as the figure is rendered.

All analysis logic lives in vhh_core.py, sanger.py and cloning.py, which a
notebook can import exactly the same way.

Run locally:   streamlit run app.py
"""

from __future__ import annotations

import io
import re
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import cloning as cl
import generator as gn
import sanger as sg
import vhh_core as vc

EXAMPLE = """>NbExample1
EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEG
RFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS
>NbExample2
QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG
RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRLSITSRRYNYWGQGTQVTVSS
>cAbLys3
DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKG
RFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS
"""

# Version stamp, YYYYMMDD.HHMM. Bump this by hand when you push a change that
# users would notice — it is what someone will quote back to you in a bug report.
VERSION = "20261021.1215"
AUTHOR = "TyrosinDuck"
AUTHOR_URL = "https://github.com/TyrosinDuck"

# The app name and its icon live here so the page title, the sidebar, the home
# page, the footer and the PDF report can never drift apart.
APP_NAME = "VHH tool box"
APP_ICON = "🧰"

st.set_page_config(page_title=APP_NAME, page_icon=APP_ICON, layout="wide")

# The analysis modules live in separate files that are updated together with this
# one. Copying only app.py across used to fail deep inside a call with an opaque
# AttributeError; this turns that into a message naming the stale file.
REQUIRED_MODULES = {"vhh_core": "2026.10.21c", "sanger": "2026.10.21b",
                    "cloning": "2026.10.21b", "generator": "2026.10.21b"}
_stale = [f"`{name}.py` (has {getattr(mod, '__version__', 'no version')}, "
          f"needs {want})"
          for name, want, mod in (("vhh_core", REQUIRED_MODULES["vhh_core"], vc),
                                  ("sanger", REQUIRED_MODULES["sanger"], sg),
                                  ("cloning", REQUIRED_MODULES["cloning"], cl),
                                  ("generator", REQUIRED_MODULES["generator"], gn))
          if getattr(mod, "__version__", None) != want]
if _stale:
    st.error("Some modules are out of date — update " + ", ".join(_stale) +
             " alongside `app.py`, they are released together.")
    st.stop()


def _fit() -> dict:
    """Return the right 'fill the container' kwarg for the installed Streamlit.

    `use_container_width` is deprecated in favour of `width="stretch"` from 1.49
    onwards; this keeps the app working on both old and new versions.
    """
    try:
        from packaging.version import Version
        if Version(st.__version__.split("+")[0]) >= Version("1.49"):
            return {"width": "stretch"}
    except Exception:
        pass
    return {"use_container_width": True}


FIT = _fit()

# How long a result may stay in the in-memory cache, and how many results are
# kept at once. The cache makes moving a sidebar slider instant instead of
# re-running the whole analysis, but it is shared by the whole server process,
# so it is kept deliberately small and short-lived. Set CACHE_TTL = 0 to switch
# caching off entirely: the app then recomputes on every interaction, which
# costs well under a second for a handful of sequences.
CACHE_TTL = 600        # seconds
CACHE_MAX_ENTRIES = 8


# --------------------------------------------------------------------------- #
# Cached workers
# --------------------------------------------------------------------------- #

@st.cache_resource
def anarci_status() -> bool:
    """Probe ANARCI once per server process, not once per request."""
    return vc.anarci_available()


def _analyze(seqs: dict[str, str], ph: float, method: str):
    return vc.analyze(seqs, ph=ph, method=method)


def _analyze_one(name: str, seq: str, ph: float, method: str):
    """Analyse a single sequence — the unit the progress bar advances on.

    Caching per sequence rather than per batch means adding one sequence to a
    list of ten only costs one analysis, not eleven.
    """
    frame, anns, errs = vc.analyze({name: seq}, ph=ph, method=method)
    return frame.iloc[0].to_dict(), anns.get(name), errs.get(name)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def hallmark_data(seq: str, annotation: dict | None = None):
    """Hallmark tetrad for one sequence (cached: ANARCI numbering is not free)."""
    return vc.hallmarks(seq, annotation)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def germline_data(seq: str):
    """Closest germline per species; two ANARCI passes, worth caching."""
    return vc.germline(seq)


if CACHE_TTL:
    run_analysis = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)(_analyze)
    analyze_one = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL,
        max_entries=CACHE_MAX_ENTRIES * 8)(_analyze_one)
else:
    run_analysis = _analyze
    analyze_one = _analyze_one


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=4)
def process_ab1(payload: tuple, vector: str, cutoff: float, min_len: int,
                min_phred: float):
    """Parse and extract a batch of .ab1 uploads (cached on the file contents)."""
    return sg.process_reads(list(payload), vector=vector, mott_cutoff=cutoff,
                            min_trimmed_len=min_len, min_mean_phred=min_phred)


def analyze_with_progress(seqs: dict[str, str], ph: float, method: str):
    """Run the analysis sequence by sequence, showing progress as it goes."""
    bar = st.progress(0.0, text="Preparing…")
    rows, anns, errs = [], {}, {}
    total = len(seqs)
    try:
        for i, (name, seq) in enumerate(seqs.items()):
            bar.progress(i / total, text=f"Analyzing {name} ({i + 1}/{total})…")
            row, ann, err = analyze_one(name, seq, ph, method)
            rows.append(row)
            if ann:
                anns[name] = ann
            if err:
                errs[name] = err
        bar.progress(1.0, text="Done")
    finally:
        bar.empty()
    return pd.DataFrame(rows), anns, errs


def fig_to_png(fig, dpi: int = 200, close: bool = False) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return buf.getvalue()


# Streamlit re-runs the whole script — and every tab body — on each interaction,
# so without caching, moving one slider re-renders every figure. These return
# PNG bytes rather than figure objects because bytes cache cleanly and a static
# image is all these plots need; the PDF builder rebuilds real figures on demand.
_FIG_CACHE = dict(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=24)


@st.cache_data(**_FIG_CACHE)
def png_secondary_structure(frame):
    return fig_to_png(vc.plot_secondary_structure(frame), close=True)


@st.cache_data(**_FIG_CACHE)
def png_charge(frame, seqs, ph):
    return fig_to_png(vc.plot_charge(frame, seqs, ph=ph), close=True)


@st.cache_data(**_FIG_CACHE)
def png_profiles(seqs, anns, window):
    return fig_to_png(vc.plot_profiles(seqs, anns, window=window), close=True)


@st.cache_data(**_FIG_CACHE)
def png_chou_fasman(name, seq, window, ann):
    return fig_to_png(vc.plot_chou_fasman(seq, name, window=window, annotation=ann),
                      close=True)


@st.cache_resource
def generator_matrix():
    """The position-frequency matrix and its provenance, loaded once."""
    return gn.load_matrix(), gn.load_meta()


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=8)
def generate_cached(n, seed, framework_mode, cdr_alphabet, lengths, hallmarks,
                    framework_cys):
    matrix, meta = generator_matrix()
    return gn.generate_many(n, matrix, meta, seed=seed,
                            framework_mode=framework_mode,
                            cdr_alphabet=cdr_alphabet,
                            cdr_lengths=dict(lengths) if lengths else None,
                            enforce_hallmarks=hallmarks,
                            framework_cys=framework_cys)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def liabilities_for(seq, ann, include_added):
    return vc.find_liabilities(seq, ann, include_added=include_added)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def proteases_for(seq, ann):
    return vc.find_protease_sites(seq, ann)


@st.cache_data(**_FIG_CACHE)
def png_liability_map(liabs, anns, lengths):
    fig = vc.plot_liability_map(liabs, anns, lengths)
    return fig_to_png(fig, close=True) if fig else None


@st.cache_data(**_FIG_CACHE)
def png_liability_profile(liabs):
    fig = vc.plot_liability_profile(liabs)
    return fig_to_png(fig, close=True) if fig else None


@st.cache_data(**_FIG_CACHE)
@st.cache_data(ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES, show_spinner=False)
def families_for(seqs, anns, cutoff, linkage_method, same_length):
    """Clonal families — cached, because the pairwise alignments dominate."""
    return vc.cluster_families(seqs, anns, cutoff=cutoff,
                               linkage_method=linkage_method,
                               require_same_length=same_length)


@st.cache_data(ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES, show_spinner=False)
def png_family_dendrogram(seqs, anns, cutoff, linkage_method, same_length):
    result = families_for(seqs, anns, cutoff, linkage_method, same_length)
    if result["linkage"] is None:
        return None
    fig, _ = vc.plot_family_dendrogram(result)
    return fig_to_png(fig, close=True)


def png_identity(seqs, anns, region):
    mat = vc.identity_matrix(seqs, anns, region=region) if region \
        else vc.identity_matrix(seqs)
    if len(mat) < 2:
        return None
    title = f"{region} identity (%)" if region else "Full-length identity (%)"
    return fig_to_png(vc.plot_identity(mat, title), close=True)


@st.cache_data(**_FIG_CACHE)
def png_msa(seqs, anns):
    with tempfile.TemporaryDirectory() as tmp:
        fig, aln_path = vc.plot_msa(seqs, anns, workdir=tmp)
        return fig_to_png(fig, close=True), aln_path.read_text()


@st.cache_data(**_FIG_CACHE)
def png_quality(payload, cutoff):
    reads = [sg.read_ab1(data, name=name, cutoff=cutoff) for name, data in payload]
    fig = sg.plot_quality(reads)
    return fig_to_png(fig, close=True) if fig else None


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def load_plasmid_cached(payload, name):
    """Parse a plasmid once per file. Returns (plasmid, feature table)."""
    plasmid = cl.load_plasmid(payload, name=name)
    return plasmid, cl.feature_table(plasmid)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def png_plasmid_map(payload, name, circular, with_sites, window=None):
    plasmid, features = load_plasmid_cached(payload, name)
    sites = (pd.DataFrame([{"enzyme": e, "site": s, "cut_positions": list(p)}
                           for e, s, p in with_sites]) if with_sites else None)
    if circular:
        fig = cl.plot_map(plasmid, features, sites)
    else:
        lo, hi = window if window else (None, None)
        fig = cl.plot_region(plasmid, features, sites, start=lo, end=hi)
    return fig_to_png(fig, close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=256)
def dna_for(seq, sites, add_stop, gc_max):
    return vc.reverse_translate(seq, avoid_sites=dict(sites), add_stop=add_stop,
                                gc_max=gc_max)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def pcr_cached(payload, name, fwd, rev, circular):
    plasmid, _f = load_plasmid_cached(payload, name)
    return cl.simulate_pcr(plasmid["record"], fwd, rev, circular=circular)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=12)
def evaluate_pairs_cached(payload, name, insert_dna, enzymes, region, layout):
    plasmid, _f = load_plasmid_cached(payload, name)
    cand = cl.recommend_enzymes(plasmid["record"], insert_dna or None, region,
                                enzymes=list(enzymes) if enzymes else None,
                                circular=plasmid["topology"] == "circular")
    pairs = cl.evaluate_pairs(plasmid["record"], insert_dna, cand, layout=layout,
                              circular=plasmid["topology"] == "circular")
    return cand, pairs


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=8)
def plan_cloning_cached(payload, name, donor_dna, donor_is_plasmid, primer_rows,
                        layout):
    """The whole cloning plan, cached on the vector, the donor and the primers."""
    plasmid, _features = load_plasmid_cached(payload, name)
    primers = pd.DataFrame([{"name": n, "sequence": s, "note": note}
                            for n, s, note in primer_rows])
    return cl.plan_cloning(plasmid, donor_dna, donor_is_plasmid, primers=primers,
                           layout=layout)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def primer_report_cached(payload, name, primer_rows):
    """Where each primer binds a plasmid. primer_rows is a tuple of (name, seq)."""
    plasmid, _features = load_plasmid_cached(payload, name)
    frame = pd.DataFrame([{"name": n, "sequence": s, "note": note}
                          for n, s, note in primer_rows])
    return cl.primer_report(plasmid, frame)


# --------------------------------------------------------------------------- #
# Navigation
# --------------------------------------------------------------------------- #

SECTIONS = {
    "analysis": ("🔬", "Sequence analysis",
                 "Characterize VHHs you provide: physico-chemical properties, "
                 "CDR/framework annotation, hallmarks and germline, liabilities, "
                 "alignment, codon-optimised DNA and a PDF report."),
    "sequencing": ("🧪", "Sanger sequencing",
                   "Upload .ab1 chromatograms from a screening run, extract the VHH "
                   "sequences from the expression cassette, spot frameshifts and "
                   "duplicates — then send the clones you want to the analysis."),
    "cloning": ("🧬", "In-silico cloning",
                "Load a backbone and an insert, and get the restriction enzymes that "
                "open the vector in the right place without cutting the VHH, plus "
                "which of your primers actually bind."),
}

st.session_state.setdefault("section", "home")


def go_to(target: str) -> None:
    """Switch section on the next rerun."""
    st.session_state["section"] = target


@st.cache_data(show_spinner=False, ttl=3600, max_entries=2)
def load_manual() -> str:
    """The user manual, read from MANUAL.md so there is one source of truth."""
    path = Path(__file__).parent / "MANUAL.md"
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        return f"The manual file could not be read ({exc})."


def current_sequences() -> dict[str, str]:
    """Sequences currently in the analysis input, parsed — or an empty dict.

    The cloning section needs them without depending on the analysis section
    having been rendered, so the text box is the single source of truth.
    """
    try:
        return vc.parse_input(st.session_state.get("seq_text", ""))
    except Exception:
        return {}


def footer() -> None:
    """Page footer — shown on every screen."""
    st.divider()
    st.caption("Built on Biopython (ProtParam), ANARCI, FAMSA, pymsaviz and "
               "dna_features_viewer. Theoretical predictions — validate "
               "experimentally. Sequences are not stored; for confidential work, "
               "run the app locally.")
    st.caption(f"{APP_NAME} · version {VERSION} · by [{AUTHOR}]({AUTHOR_URL})")


# A widget's state cannot be written after the widget exists, so anything that
# wants to fill the text area (the sequencing section) parks its text here and
# reruns; it is moved into place before the widget is created.
if "pending_seq_text" in st.session_state:
    st.session_state["seq_text"] = st.session_state.pop("pending_seq_text")
st.session_state.setdefault("seq_text", "")

section = st.session_state["section"]


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #

with st.sidebar:
    st.markdown(f"### {APP_ICON} {APP_NAME}")
    st.button("🏠 Home", key="nav_home", **FIT,
              disabled=section == "home", on_click=go_to, args=("home",))
    for key, (icon, label, _desc) in SECTIONS.items():
        st.button(f"{icon} {label}", key=f"nav_{key}", **FIT,
                  disabled=section == key, on_click=go_to, args=(key,))
    st.button("📖 Manual", key="nav_manual", **FIT,
              disabled=section == "manual", on_click=go_to, args=("manual",))

    st.divider()

    if section == "analysis":
        st.markdown("#### Settings")
        ph = st.slider("pH for net charge", 3.0, 10.0, vc.DEFAULTS["ph"], 0.1,
                       help="Net charge and the stability flags use this pH.")
        window = st.slider("Profile window (residues)", 3, 21,
                           vc.DEFAULTS["profile_window"], 2,
                           help="Sliding window for the per-residue profiles.")
        has_anarci = anarci_status()
        options = ["auto", "regex"] + (["anarci"] if has_anarci else [])
        method = st.selectbox(
            "CDR annotation method", options, index=0,
            help="'auto' uses ANARCI when available and falls back to the regex "
                 "anchors otherwise.")
        st.caption(f"ANARCI: "
                   f"{'available' if has_anarci else 'not installed — regex only'}")
        st.divider()
    else:
        ph = vc.DEFAULTS["ph"]
        window = vc.DEFAULTS["profile_window"]
        method = "auto"

    with st.expander("🔒 How your sequences are handled"):
        st.markdown(f"""
Sequences, chromatograms and plasmids are **never written to a persistent file,
logged, or sent to a third party**. They are processed in the server's memory and
the results are kept in a short-lived cache ({CACHE_TTL // 60} min,
{CACHE_MAX_ENTRIES} entries max) that only exists to keep the controls responsive.

This app runs on Streamlit Community Cloud, i.e. on third-party infrastructure.
Traffic is encrypted, nothing is stored — but the data does pass through a server
that is not yours.

**For confidential or unpublished material**, run the app on your own machine
instead. The code is open source: clone the repository, `pip install -r
requirements.txt`, `streamlit run app.py`, and nothing ever leaves your computer.
""")
    st.caption("All values are theoretical, derived from the sequence alone. "
               "Confirm pI and extinction coefficient experimentally before "
               "relying on them.")


# --------------------------------------------------------------------------- #
# Home
# --------------------------------------------------------------------------- #

if section == "home":
    st.title(f"{APP_ICON} {APP_NAME}")
    st.markdown(
        "A workbench for single-domain antibodies: characterize them, pull them out "
        "of sequencing reads, and plan the cloning. Pick where to start.")

    cols = st.columns(3)
    for col, (key, (icon, label, desc)) in zip(cols, SECTIONS.items()):
        with col, st.container(border=True):
            st.markdown(f"### {icon} {label}")
            st.markdown(desc)
            st.button(f"Open {label.lower()}", key=f"home_{key}", type="primary",
                      **FIT, on_click=go_to, args=(key,))

    st.write("")
    hcols = st.columns([1, 2])
    with hcols[0]:
        st.button("📖 Read the manual", key="home_manual", **FIT,
                  on_click=go_to, args=("manual",))
    with hcols[1]:
        st.caption("Every screen explained, what each number means, and how far to "
                   "trust it.")

    staged = current_sequences()
    if staged:
        st.info(f"{len(staged)} sequence(s) are currently loaded in the analysis "
                f"input: {', '.join(list(staged)[:6])}"
                f"{'…' if len(staged) > 6 else ''}")

    with st.expander("How the three sections fit together"):
        st.markdown("""
**Sanger sequencing → Sequence analysis.** After extracting VHHs from a set of
chromatograms you select the clones worth keeping and send them straight to the
analysis section, which is where every characterization lives. Duplicates are
deselected by default, so one entry per unique protein sequence.

**Sequence analysis → In-silico cloning.** The cloning section reads the same
sequence list, so whatever you are analysing is available as an insert without
retyping it.

Each section also works on its own: you can go straight to cloning with a pasted
DNA insert, or straight to analysis with a FASTA file.
""")
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# Manual
# --------------------------------------------------------------------------- #

if section == "manual":
    manual = load_manual()
    top = st.columns([3, 1])
    with top[1]:
        st.download_button("Download (Markdown)", manual, file_name="VHH_manual.md",
                           mime="text/markdown", **FIT)
    st.markdown(manual)
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# Sanger sequencing
# --------------------------------------------------------------------------- #

if section == "sequencing":
    st.title("🧪 Sanger sequencing")
    st.markdown(
        "Upload the chromatograms from a screening run. Each read is quality-trimmed, "
        "then the PelB–VHH–tag cassette is located in any strand and reading frame, "
        "and the VHH protein sequence is extracted. Pick the clones you want and send "
        "them to the analysis section.")

    ab1_files = st.file_uploader("Chromatograms (.ab1)", type=["ab1"],
                                 accept_multiple_files=True, key="ab1_uploader")

    opts = st.columns(4)
    with opts[0]:
        vector = st.selectbox("Vector", list(sg.VECTORS), index=0,
                              help="Determines which tags are expected after the VHH.")
    with opts[1]:
        mott_cutoff = st.slider("Trim cutoff", 0.01, 0.20,
                                sg.DEFAULTS["mott_cutoff"], 0.01,
                                help="Mott algorithm error-probability cutoff. "
                                     "Lower is stricter.")
    with opts[2]:
        min_len = st.number_input("Min trimmed length (nt)", 50, 800,
                                  sg.DEFAULTS["min_trimmed_len"], 10)
    with opts[3]:
        min_phred = st.number_input("Min mean Phred", 5, 50,
                                    sg.DEFAULTS["min_mean_phred"], 1)
    st.caption(f"{vector} layout: " + " → ".join(sg.VECTORS[vector]["layout"]))

    if not ab1_files:
        st.info("No chromatogram uploaded yet.")
    else:
        payload = tuple((f.name, f.getvalue()) for f in ab1_files)
        bar = st.progress(0.0, text="Reading chromatograms…")
        try:
            sanger_df, reads = process_ab1(payload, vector, mott_cutoff,
                                           int(min_len), float(min_phred))
        finally:
            bar.empty()

        n_total = len(sanger_df)
        n_usable = int(sanger_df["usable"].sum()) if "usable" in sanger_df else 0
        n_unique = sanger_df[sanger_df["usable"]]["vhh_aa"].nunique() if n_usable else 0
        cols_m = st.columns(4)
        cols_m[0].metric("Reads", n_total)
        cols_m[1].metric("VHH extracted", n_usable)
        cols_m[2].metric("Unique sequences", n_unique)
        cols_m[3].metric("Failed", n_total - n_usable)

        show_cols = [c for c in ["read", "status", "qc_pass", "trimmed_len",
                                 "mean_phred_trim", "strand", "frame", "vhh_len",
                                 "tags_found", "tags_missing", "premature_stops",
                                 "n_copies", "duplicate_of", "warnings"]
                     if c in sanger_df.columns]
        st.dataframe(
            sanger_df[show_cols], hide_index=True, **FIT,
            column_config={
                "warnings": st.column_config.TextColumn("Warnings", width="large"),
                "mean_phred_trim": st.column_config.NumberColumn("Mean Q",
                                                                 format="%.1f"),
                "trimmed_len": st.column_config.NumberColumn("Trimmed (nt)"),
                "vhh_len": st.column_config.NumberColumn("VHH (aa)"),
                "duplicate_of": st.column_config.TextColumn("Duplicate of"),
            })
        st.caption("Drag a column edge to widen it, or use the fullscreen button at "
                   "the top right of the table.")

        flagged = (sanger_df[sanger_df["warnings"].astype(bool)]
                   if "warnings" in sanger_df.columns else sanger_df.iloc[0:0])
        if len(flagged):
            with st.expander(f"Warnings in full ({len(flagged)} read(s))",
                             expanded=True):
                for _, row in flagged.iterrows():
                    st.markdown(f"**{row['read']}** — *{row['status']}*")
                    for w in str(row["warnings"]).split("; "):
                        if w:
                            st.markdown(f"- {w}")

        if n_total - n_usable:
            st.caption("A failed read is not necessarily a failed clone: "
                       "'no FR4' usually means the read stopped before the end of the "
                       "VHH, and a frameshift means PelB and FR4 sit in different "
                       "reading frames, so an indel lies between them.")

        if reads:
            with st.expander("Quality traces"):
                img_q = png_quality(payload, mott_cutoff)
                if img_q:
                    st.image(img_q, **FIT)
                    st.caption("Blue is the window kept after trimming; the dashed "
                               "line is Phred 20.")
                plt.close("all")

        usable_df = sanger_df[sanger_df["usable"]] if n_usable else sanger_df.iloc[0:0]
        if n_usable:
            st.subheader("Send to the analysis section")
            unique_first = usable_df[usable_df["duplicate_of"] == ""]["read"].tolist()
            chosen = st.multiselect(
                "Sequences to analyse", usable_df["read"].tolist(),
                default=unique_first,
                help="Duplicates are left out by default — one entry per unique "
                     "protein sequence.")
            replace = st.checkbox("Replace what is already loaded there", value=True)
            if st.button("Send selected to analysis", type="primary",
                         disabled=not chosen, **FIT):
                picked = usable_df[usable_df["read"].isin(chosen)]
                fasta = sg.reads_to_fasta(picked, only_usable=False)
                current = st.session_state.get("seq_text", "").strip()
                st.session_state["pending_seq_text"] = (
                    fasta if replace or not current else current + "\n" + fasta)
                go_to("analysis")
                st.rerun()

        st.subheader("Download")
        dl = st.columns(3)
        with dl[0]:
            st.download_button("Sequencing results (CSV)",
                               sanger_df.drop(columns=["vhh_nt"], errors="ignore")
                               .to_csv(index=False),
                               file_name="sanger_results.csv", mime="text/csv", **FIT)
        with dl[1]:
            st.download_button("VHH proteins (FASTA)", sg.reads_to_fasta(sanger_df),
                               file_name="vhh_proteins.fasta", mime="text/plain",
                               disabled=not n_usable, **FIT)
        with dl[2]:
            st.download_button("VHH DNA (FASTA)",
                               sg.reads_to_fasta(sanger_df, column="vhh_nt"),
                               file_name="vhh_dna.fasta", mime="text/plain",
                               disabled=not n_usable, **FIT)

    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# In-silico cloning
# --------------------------------------------------------------------------- #

if section == "cloning":
    st.title("🧬 In-silico cloning")
    st.markdown(
        "Pick the recipient vector, paste the donor DNA, and the app works out which "
        "restriction enzymes put the VHH between the signal peptide and the tags — "
        "then checks that the cloned VHH encodes exactly the protein you started "
        "from.")

    # ---- 1. recipient ----
    st.subheader("1 · Recipient vector")
    bundled = cl.list_bundled_plasmids()
    rcol = st.columns([2, 1])
    with rcol[0]:
        backbone_name = st.selectbox("Vector", list(bundled) + ["Upload my own…"],
                                     index=0)
    with rcol[1]:
        layout = st.selectbox(
            "Expected tags", list(sg.VECTORS), index=0,
            help="Which tags the finished construct should carry — used when the "
                 "cassette is translated back.")

    payload, display_name = None, backbone_name
    if backbone_name == "Upload my own…":
        up = st.file_uploader("Vector (SnapGene .dna, GenBank or FASTA)",
                              type=["dna", "gb", "gbk", "genbank", "fasta", "fa",
                                    "txt"], key="backbone_upload")
        if up is not None:
            payload, display_name = up.getvalue(), up.name
    else:
        payload = bundled[backbone_name].read_bytes()

    if payload is None:
        st.info("Upload a vector to continue.")
        footer()
        st.stop()

    try:
        plasmid, features = load_plasmid_cached(payload, display_name)
    except Exception as exc:
        st.error(f"Could not read this vector — {exc}")
        footer()
        st.stop()

    if plasmid["warning"]:
        st.warning(plasmid["warning"])
    stats = cl.sequence_stats(plasmid["record"])
    auto_region = cl.insertion_region(plasmid, features)
    mcols = st.columns(4)
    mcols[0].metric("Length", f"{stats['length']} bp")
    mcols[1].metric("GC", f"{stats['gc_percent']} %")
    mcols[2].metric("Features", len(features))
    mcols[3].metric("Cloning window",
                    f"{auto_region[0]}–{auto_region[1]}" if auto_region else "—")
    if auto_region:
        st.caption(f"The insert goes between the end of the signal peptide and the "
                   f"first tag, positions {auto_region[0]}–{auto_region[1]}.")
    else:
        st.error("Could not find a leader and a tag in this vector, so there is no "
                 "window to clone into. Check the annotation.")

    # ---- 2. donor ----
    st.subheader("2 · Donor DNA")
    donor_kind = st.radio(
        "What are you pasting",
        ["A plasmid that carries the VHH", "A VHH fragment amplified with primers"],
        horizontal=True)
    donor_is_plasmid = donor_kind.startswith("A plasmid")

    if donor_is_plasmid:
        st.caption("The app locates the VHH, picks the primer pair that brackets it, "
                   "and uses the amplicon as the insert.")
    else:
        st.caption("The fragment must already carry the restriction sites its primers "
                   "added — a bare coding sequence cannot be cloned by restriction.")

    dfile = st.file_uploader(
        "Optional: load a file instead of pasting (SnapGene .dna, GenBank, FASTA)",
        type=["dna", "gb", "gbk", "genbank", "fasta", "fa", "txt"],
        key="donor_file")
    if dfile is not None and dfile.file_id != st.session_state.get("last_donor_file"):
        st.session_state["last_donor_file"] = dfile.file_id
        try:
            loaded = cl.load_plasmid(dfile.getvalue(), name=dfile.name)
            st.session_state["pending_donor_dna"] = str(loaded["record"].seq)
            st.rerun()
        except Exception as exc:
            st.error(f"Could not read that file — {exc}")

    if "pending_donor_dna" in st.session_state:
        st.session_state["donor_dna_box"] = st.session_state.pop("pending_donor_dna")
    donor_text = st.text_area("Donor DNA (5'→3')", height=150, key="donor_dna_box",
                              placeholder="ATGAAATACCTATTGCCTACGGCAGCCGCT…")
    donor_dna = re.sub(r"[^ACGTacgt]", "", donor_text or "").upper()
    if donor_dna:
        st.caption(f"{len(donor_dna)} nt")

    if not donor_dna:
        st.info("Paste the donor DNA, or load a file, to plan the cloning.")
        footer()
        st.stop()

    # ---- 3. the plan ----
    st.subheader("3 · Plan")
    own = st.file_uploader("Optional: your own primer list (CSV with name,sequence)",
                           type=["csv"], key="primer_upload")
    if own is not None:
        primers = pd.read_csv(io.BytesIO(own.getvalue()))
        primers["sequence"] = primers["sequence"].astype(str).str.upper()
        if "note" not in primers.columns:
            primers["note"] = ""
    else:
        primers = cl.load_primers()
    primer_rows = tuple((r["name"], r["sequence"], str(r.get("note", "")))
                        for _, r in primers.iterrows())

    with st.spinner("Locating the VHH, amplifying, cutting and translating…"):
        plan = plan_cloning_cached(payload, display_name, donor_dna,
                                   donor_is_plasmid, primer_rows, layout)

    if not plan["donor_vhh"]:
        st.error("No VHH could be found in that DNA. " +
                 " ".join(plan["messages"]))
        footer()
        st.stop()

    st.success(f"VHH found in the donor — {len(plan['donor_vhh'])} aa.")
    with st.expander("Donor VHH sequence"):
        st.code("\n".join(plan["donor_vhh"][i:i + 60]
                          for i in range(0, len(plan["donor_vhh"]), 60)),
                language="text")

    if plan.get("primer_pair"):
        f = plan["primer_pair"]["forward"]
        r = plan["primer_pair"]["reverse"]
        st.markdown(
            f"**Primers chosen automatically** — `{f['name']}` as the forward primer "
            f"({f['form']}, position {f['start']}) and `{r['name']}` as the reverse "
            f"({r['form']}, position {r['start']}), giving a "
            f"**{len(plan['amplicon'])} bp** amplicon.")
        if r["form"] == "reverse complement":
            st.caption(
                f"`{r['name']}` is listed as the plus-strand sequence of its target, "
                "so the oligo that primes back towards the VHH is its reverse "
                "complement — that is the form used here.")

    for msg in plan["messages"]:
        st.warning(msg)

    if not plan["plans"]:
        footer()
        st.stop()

    usable = [p for p in plan["plans"] if p["usable"]]
    if not usable:
        st.error("No enzyme pair gives a construct that preserves the donor's VHH. "
                 "The table below says why each one was rejected.")
    else:
        best = usable[0]
        st.markdown("### Recommended")
        bc = st.columns(4)
        bc[0].metric("Enzymes", f"{best['enzyme_5']} + {best['enzyme_3']}")
        bc[1].metric("Construct", f"{best['construct_bp']} bp")
        bc[2].metric("VHH", f"{len(best['construct_vhh'])} aa")
        bc[3].metric("Vector removed", f"{best['removed_bp']} bp")

        comparison = best["comparison"]
        if comparison["status"] == "identical":
            st.success("**Protein check passed** — the VHH encoded by the construct "
                       "is exactly the one in the donor.")
        else:
            st.info(f"**Protein check** — {comparison['message']}")

        pcols = st.columns(2)
        with pcols[0]:
            st.markdown("**VHH in the donor**")
            st.code("\n".join(plan["donor_vhh"][i:i + 40]
                              for i in range(0, len(plan["donor_vhh"]), 40)),
                    language="text")
        with pcols[1]:
            st.markdown("**VHH in the construct**")
            st.code("\n".join(best["construct_vhh"][i:i + 40]
                              for i in range(0, len(best["construct_vhh"]), 40)),
                    language="text")

        check = cl.verify_construct(best["construct"], layout)
        st.markdown(f"**Cassette** — signal peptide, then the VHH, then "
                    f"{check['tags_found'] or 'no tag'}; "
                    f"{check['premature_stops']} premature stop(s).")
        st.code("\n".join(check["cassette_aa"][i:i + 60]
                          for i in range(0, len(check["cassette_aa"]), 60)),
                language="text")

        built = cl.assemble(plasmid["record"], plan["amplicon"],
                            best["enzyme_5"], best["enzyme_3"],
                            circular=plasmid["topology"] == "circular")
        cname = f"{plasmid['name']}_VHH".replace(" ", "_")[:24]
        crec = cl.construct_record(best["construct"], cname,
                                   built.get("insert_start", 1),
                                   built.get("insert_end", 1), layout=layout,
                                   circular=plasmid["topology"] == "circular")
        dl = st.columns(3)
        with dl[0]:
            st.download_button("Construct (GenBank)", cl.record_to_genbank(crec),
                               file_name=f"{cname}.gb", mime="text/plain",
                               type="primary", **FIT)
        with dl[1]:
            st.download_button(
                "Construct (FASTA)",
                f">{cname}\n" + "\n".join(best["construct"][i:i + 60]
                                          for i in range(0, len(best["construct"]),
                                                         60)) + "\n",
                file_name=f"{cname}.fasta", mime="text/plain", **FIT)
        with dl[2]:
            if plan["amplicon"] != donor_dna:
                st.download_button(
                    "Amplicon (FASTA)",
                    ">amplicon\n" + "\n".join(
                        plan["amplicon"][i:i + 60]
                        for i in range(0, len(plan["amplicon"]), 60)) + "\n",
                    file_name="amplicon.fasta", mime="text/plain", **FIT)

        st.markdown(f"""
**Protocol** — digest the {len(plan['amplicon'])} bp insert and
**{plasmid['name']}** with **{best['enzyme_5']}** and **{best['enzyme_3']}**, then
ligate. The vector loses {best['removed_bp']} bp between the two sites and gains
{best['inserted_bp']} bp from the insert. Screen the colonies by sequencing — the
sequencing section of this app reads the .ab1 files back and confirms the cassette.
""")

    # ---- details, folded away ----
    with st.expander(f"All enzyme pairs considered ({len(plan['plans'])})"):
        rows = [{"5' enzyme": p["enzyme_5"], "3' enzyme": p["enzyme_3"],
                 "usable": p["usable"], "construct bp": p["construct_bp"],
                 "vector removed": p["removed_bp"],
                 "protein": p["comparison"]["status"],
                 "deletes": ", ".join(p["deletes"]) or "—",
                 "why not": ("" if p["usable"] else
                             (", ".join(p["deletes"]) + " would be deleted"
                              if p["deletes"] else p["comparison"]["message"]))}
                for p in plan["plans"]]
        st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT,
                     column_config={"why not": st.column_config.TextColumn(
                         "Why not", width="large")})
        st.caption("Each pair is cut, joined and translated for real, then the VHH it "
                   "encodes is compared with the donor's. A pair that would delete an "
                   "annotated part of the vector is rejected even when the cassette "
                   "still translates.")

    with st.expander("Restriction sites"):
        cand = plan["candidates"]
        if len(cand):
            st.dataframe(
                cand[["enzyme", "site", "backbone_cuts", "position", "in_region",
                      "insert_cuts", "overhang", "verdict"]].head(30),
                hide_index=True, **FIT,
                column_config={"verdict": st.column_config.TextColumn(
                    "Verdict", width="large")})
            st.caption("An enzyme that cuts the insert is not necessarily wrong: in "
                       "VHH cloning BstEII is used on purpose, its site sitting in "
                       "FR4 at the 3' end of the domain.")

    with st.expander("Primers on this vector"):
        full_report = primer_report_cached(payload, display_name, primer_rows)
        binding = full_report[full_report["binds"]] if len(full_report) else full_report
        if len(binding):
            st.dataframe(binding[["primer", "sequence", "orientation", "start", "end",
                                  "tm"]], hide_index=True, **FIT,
                         column_config={"tm": st.column_config.NumberColumn(
                             "Tm (°C)", format="%.1f")})
        absent = [p for p in primers["name"] if p not in set(binding["primer"])] \
            if len(binding) else list(primers["name"])
        if absent:
            st.caption(f"No binding site in {plasmid['name']}: {', '.join(absent)}. "
                       "That is usually a property of the vector rather than a "
                       "mistake — pHEN6 has no gIII for a gIII primer to sit in.")

    with st.expander("Vector map"):
        view = st.radio("View", ["Cloning region", "Whole plasmid"], horizontal=True,
                        key="map_view")
        sites = (tuple((p["enzyme_5"], "", tuple([p["removed_bp"]])) for p in [])
                 or ())
        cand = plan["candidates"]
        usable_sites = (tuple((r["enzyme"], r["site"], tuple(r["cut_positions"]))
                              for _, r in cand[cand["usable"]].iterrows())
                        if len(cand) else ())
        try:
            img_map = png_plasmid_map(
                payload, display_name, view == "Whole plasmid", usable_sites,
                (max(1, (auto_region or (1, stats["length"]))[0] - 160),
                 min(stats["length"], (auto_region or (1, stats["length"]))[1] + 160)))
            st.image(img_map, **FIT)
            st.caption("Enzyme sites in red; features from the file in blue, motifs "
                       "detected here in orange, a VHH scaffold in green.")
        except Exception as exc:
            st.warning(f"Map unavailable: {exc}")
        plt.close("all")

    with st.expander(f"Vector features ({len(features)})"):
        st.dataframe(
            features[[c for c in ["feature", "kind", "source", "start", "end",
                                  "strand", "type", "note"]
                      if c in features.columns]], hide_index=True, **FIT)
        st.caption("The *source* column says whether a feature was read from the file "
                   "or inferred here by motif scan.")

    footer()
    st.stop()

# --------------------------------------------------------------------------- #
# Sequence analysis
# --------------------------------------------------------------------------- #

st.title("🔬 Sequence analysis")
st.markdown(
    "Physico-chemical characterization, CDR/framework annotation and alignment of "
    "single-domain antibodies. Paste one or more sequences, in FASTA or raw form.")

col_in, col_act = st.columns([3, 1])
with col_act:
    st.write("")
    if st.button("Load example", **FIT):
        st.session_state["seq_text"] = EXAMPLE
        st.rerun()
    uploaded = st.file_uploader("or upload FASTA", type=["fasta", "fa", "txt", "faa"])
    if uploaded is not None and uploaded.file_id != st.session_state.get("last_upload"):
        st.session_state["last_upload"] = uploaded.file_id
        st.session_state["seq_text"] = uploaded.read().decode("utf-8", errors="replace")
        st.rerun()
    st.caption(f"Up to {vc.DEFAULTS['max_sequences']} sequences, "
               f"{vc.DEFAULTS['max_length']} aa each.")
    st.button("🧪 Get them from .ab1 files", key="to_sanger", **FIT,
              on_click=go_to, args=("sequencing",))

with col_in:
    # The widget owns st.session_state["seq_text"]; the buttons above write to it
    # before the widget is instantiated, which is the supported way to prefill it.
    text = st.text_area("Sequences", height=220, key="seq_text",
                        placeholder=">MyNanobody\nEVQLVESGGGLVQ...")
    # Clicking the button moves focus out of the text area, which is what commits
    # what was just typed — the same thing ⌘/Ctrl + Enter does, without the
    # keyboard shortcut nobody discovers on their own.
    st.button("Analyze", type="primary", **FIT)
    st.caption("Click Analyze, or press ⌘ + Enter (Ctrl + Enter on Windows/Linux).")

# --- Synthetic sequence generator ------------------------------------------- #
with st.expander("🎲 Generate synthetic VHH sequences"):
    matrix, gmeta = generator_matrix()
    report = gn.matrix_report(matrix, gmeta)
    st.markdown(
        "Frameworks are sampled position by position from a real IMGT "
        "position-frequency matrix, so they look like natural frameworks. CDRs are "
        "drawn at random, with lengths taken from the same reference set.")
    st.caption(
        f"Reference: **{report['source']}** — {report['n_sequences']} unique VHH "
        f"sequences, {report['n_positions']} IMGT positions, median support "
        f"{report['median_support']} sequences per position.")
    st.caption(
        f"*{gmeta.get('source_title', '')}* — "
        f"[doi:{gmeta.get('source_doi', '')}]"
        f"(https://doi.org/{gmeta.get('source_doi', '')}) · "
        f"[database]({gmeta.get('source_url', '')}). "
        "Please cite it if you publish work based on these sequences. What is "
        "bundled here is an aggregate frequency matrix, not the sequences "
        "themselves.")

    gc1, gc2, gc3 = st.columns(3)
    with gc1:
        n_gen = st.number_input("How many", 1, 200, 10, 1)
        seed = st.number_input("Random seed", 0, 10**6, 0, 1,
                               help="The same seed always gives the same set.")
    with gc2:
        fw_mode = st.radio("Frameworks", ["sample", "consensus"],
                           help="'sample' draws each position from its frequency "
                                "column and gives a different framework every "
                                "time; 'consensus' takes the most frequent "
                                "residue, giving one canonical framework.")
        cdr_alpha = st.radio(
            "CDR residues", ["natural", "no_cys", "uniform"],
            help="'natural' uses the reference set's own CDR composition; "
                 "'no_cys' is uniform over 19 residues; 'uniform' is all 20.")
    with gc3:
        keep_hallmarks = st.checkbox("Enforce camelid hallmarks", value=True)
        allow_fw_cys = st.checkbox("Allow cysteines in frameworks", value=False)
        fixed_lengths = st.checkbox("Fix the CDR lengths", value=False)

    lengths = {}
    if fixed_lengths:
        lc = st.columns(3)
        for col, region, default in zip(lc, ("CDR1", "CDR2", "CDR3"), (8, 8, 17)):
            with col:
                lengths[region] = st.number_input(f"{region} length", 2, 30,
                                                  default, 1)
    else:
        st.caption("CDR lengths are drawn from the reference distribution — "
                   "CDR1 and CDR2 around 8 residues, CDR3 a mean of 17 with a "
                   "range of 2 to 30.")

    if st.button("Generate", type="primary"):
        with st.spinner(f"Generating {int(n_gen)} sequence(s)…"):
            st.session_state["generated"] = generate_cached(
                int(n_gen), int(seed), fw_mode, cdr_alpha,
                tuple(sorted(lengths.items())), keep_hallmarks, allow_fw_cys)

    made = st.session_state.get("generated")
    if made is not None and len(made):
        st.dataframe(made[["name", "length", "tetrad", "CDR1_len", "CDR2_len",
                           "CDR3_len", "sequence"]], hide_index=True, **FIT,
                     column_config={"sequence": st.column_config.TextColumn(
                         "Sequence", width="large")})
        gcols = st.columns(2)
        with gcols[0]:
            st.download_button("Download (FASTA)", gn.to_fasta(made),
                               file_name="synthetic_vhh.fasta", mime="text/plain",
                               **FIT)
        with gcols[1]:
            if st.button("Send to the analysis above", **FIT):
                st.session_state["pending_seq_text"] = gn.to_fasta(made)
                st.rerun()

    st.caption(
        "Two limitations worth knowing. **Positions are sampled independently**, so "
        "a combination that never co-occurs in nature can come out — the hallmark "
        "tetrad is the usual casualty, which is why it is enforced by default. And "
        "**a frequency is only as good as its support**: a position seen in 30 "
        "reference sequences carries far less information than one seen in 1900. "
        "These are synthetic sequences for controls and benchmarks; nothing here "
        "predicts that they fold or bind anything.")

if not text.strip():
    st.info("Paste sequences above, upload a FASTA file, load the example, generate "
            "synthetic ones, or extract them from .ab1 chromatograms in the "
            "sequencing section.")
    footer()
    st.stop()

try:
    sequences = vc.parse_input(text)
except ValueError as exc:
    st.error(f"Input error — {exc}")
    footer()
    st.stop()

df, annotations, errors = analyze_with_progress(sequences, ph, method)

if errors:
    for name, msg in errors.items():
        st.warning(f"**{name}** — annotation failed: {msg}. "
                   "Physico-chemical values are still reported below.")

# Named rather than indexed: inserting a tab in the middle of an indexed list
# silently reassigns every panel after it.
TAB_LABELS = ["Properties", "Plots", "CDR / Framework", "Hallmarks & germline",
              "Liabilities", "Alignment", "Clonal families", "DNA", "Download"]
tabs = dict(zip(TAB_LABELS, st.tabs(TAB_LABELS)))

# --- Properties ------------------------------------------------------------ #
with tabs["Properties"]:
    main_cols = ["name", "length", "mw_kda", "pI", "ext_coeff_reduced", "A280_0.1pct",
                 "instability_index", "aliphatic_index", "gravy", "aromaticity",
                 "net_charge", "n_Cys"]
    st.dataframe(
        df[main_cols].style.format(precision=2).background_gradient(
            subset=["gravy", "instability_index"], cmap="RdYlGn_r"),
        **FIT, hide_index=True,
        column_config={
            "mw_kda": st.column_config.NumberColumn("MW (kDa)"),
            "ext_coeff_reduced": st.column_config.NumberColumn("ε280 (M⁻¹cm⁻¹)"),
            "A280_0.1pct": st.column_config.NumberColumn("A280 (1 mg/mL)"),
            "net_charge": st.column_config.NumberColumn(f"Charge pH {ph}"),
        })

    with st.expander("What these descriptors mean"):
        st.markdown("""
| Descriptor | Interpretation |
|---|---|
| **Theoretical pI** | pH of zero net charge. A pI close to the buffer pH predicts poor colloidal stability |
| **ε280 / A280 (1 mg/mL)** | extinction coefficient and the factor converting a measured A280 into mg/mL |
| **Instability index** | Guruprasad dipeptide score; **< 40 predicts a stable protein** |
| **Aliphatic index** | relative volume of A, V, I, L; higher correlates with thermostability |
| **GRAVY** | mean Kyte–Doolittle hydropathy; negative = hydrophilic. Above ≈ −0.2 flags aggregation risk |
| **Aromaticity** | fraction of F, W, Y; drives the A280 signal |
| **Secondary structure** | propensity estimate from composition, **not** a structure prediction |
""")

    st.subheader("Developability flags")
    for _, row in df.iterrows():
        flags = vc.stability_flags(row, ph=ph)
        with st.container(border=True):
            st.markdown(f"**{row['name']}**")
            if flags:
                for f in flags:
                    st.markdown(f"- ⚠️ {f}")
            else:
                st.markdown("- ✅ No flag raised")

    st.subheader("Custom features")
    feature_cols = ["name"] + [c for c in vc.FEATURE_REGISTRY if c in df.columns]
    st.dataframe(df[feature_cols], hide_index=True, **FIT)
    with st.expander("Feature definitions"):
        for name, entry in vc.FEATURE_REGISTRY.items():
            st.markdown(f"- **`{name}`** — {entry['doc']}")

# --- Plots ----------------------------------------------------------------- #
with tabs["Plots"]:
    st.image(png_secondary_structure(df), **FIT)

    st.image(png_charge(df, sequences, ph), **FIT)
    st.caption("Each titration curve crosses zero at the theoretical pI of that "
               "sequence.")

    st.image(png_profiles(sequences, annotations, window), **FIT)
    st.caption("Sliding-window charge (blue) and Kyte–Doolittle hydropathy (orange). "
               "Shaded bands are the CDRs.")

    st.subheader("Predicted secondary structure")
    cf_window = st.slider("Chou–Fasman window", 3, 11,
                          vc.DEFAULTS["chou_fasman_window"], 2, key="cf_window")
    # With many sequences, drawing every panel on every rerun is the single
    # slowest thing the app does, so only the selected one is shown by default.
    cf_names = list(sequences)
    if len(cf_names) > 1:
        show_all_cf = st.checkbox(f"Show all {len(cf_names)} sequences", value=False,
                                  key="cf_all")
        cf_shown = cf_names if show_all_cf else [
            st.selectbox("Sequence", cf_names, key="cf_pick")]
    else:
        cf_shown = cf_names

    for name in cf_shown:
        seq = sequences[name]
        st.image(png_chou_fasman(name, seq, cf_window, annotations.get(name)), **FIT)
        prof = vc.chou_fasman_profile(seq, cf_window)
        st.caption(f"{name}: {prof['fractions']['E']:.0%} sheet, "
                   f"{prof['fractions']['H']:.0%} helix, "
                   f"{prof['fractions']['C']:.0%} coil by this heuristic.")
    st.caption("The Chou–Fasman heuristic is roughly 50–60% accurate: it shows where "
               "the sequence tends towards each conformation, it does not predict the "
               "fold. A VHH is an immunoglobulin β-sandwich, so a sheet-dominated "
               "profile is the expected result. For real secondary structure, run "
               "DSSP on a model.")
    plt.close("all")

# --- CDR / Framework ------------------------------------------------------- #
with tabs["CDR / Framework"]:
    if not annotations:
        st.warning("No sequence could be annotated.")
    else:
        rows = [{"name": n,
                 "method": df.loc[df["name"] == n, "annotation_method"].iloc[0],
                 "CDR1": annotations[n]["CDR1"]["seq"],
                 "CDR2": annotations[n]["CDR2"]["seq"],
                 "CDR3": annotations[n]["CDR3"]["seq"],
                 "CDR1 len": len(annotations[n]["CDR1"]["seq"]),
                 "CDR2 len": len(annotations[n]["CDR2"]["seq"]),
                 "CDR3 len": len(annotations[n]["CDR3"]["seq"])}
                for n in annotations]
        st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT)

        st.subheader("Annotated sequences")
        st.markdown(
            '<span style="background:#e41a1c;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR1</span> '
            '<span style="background:#377eb8;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR2</span> '
            '<span style="background:#4daf4a;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR3</span>', unsafe_allow_html=True)
        st.markdown("".join(vc.render_annotated_html(n, annotations[n])
                            for n in annotations), unsafe_allow_html=True)

        with st.expander("Framework regions"):
            st.dataframe(pd.DataFrame(
                [{"name": n, **{r: annotations[n][r]["seq"]
                                for r in ("FR1", "FR2", "FR3", "FR4")}}
                 for n in annotations]), hide_index=True, **FIT)

        with st.expander("How the two annotation methods work"):
            st.markdown("""
**ANARCI (IMGT numbering)** — the sequence is aligned to germline HMMs and each residue
receives an IMGT position; CDRs are the fixed windows **27–38**, **56–65** and **105–117**.
Numbering-based, so it tolerates unusual CDR lengths. Requires HMMER.

**Regex anchors** — dependency-free; locates the conserved landmarks of a VHH and cuts between them:

| Anchor | Pattern | Role |
|---|---|---|
| conserved Cys23 | `C` in residues 16–30 | ends FR1 |
| FR2 | `W[A-Z]RQ` (the `WFRQ`/`WYRQ` hallmark) | starts FR2 |
| FR3 | `[KR][FYW][TASV][IVLFM]S` (the `RFTIS` motif) | ends CDR2 |
| second Cys | `[YFHV][YFHVAC]C` (the `YYC` motif) | ends FR3 |
| FR4 | `WG[A-Z]G[TAS][A-Z]{1,2}VTVSS` | ends CDR3 |

Assumes a canonical scaffold. Where both methods run, a disagreement usually means a
non-canonical framework — which is itself worth knowing.
""")

# --- Hallmarks & germline --------------------------------------------------- #
with tabs["Hallmarks & germline"]:
    st.markdown("The four FR2 positions that distinguish a camelid VHH from a "
                "conventional VH, and the closest germline genes.")
    hm_rows, gl_rows = [], []
    for name, seq in sequences.items():
        try:
            hm = hallmark_data(seq, annotations.get(name))
            hm_rows.append({"name": name, "tetrad": hm["tetrad"],
                            "camelid score": f"{hm['camelid_score']}/4",
                            **{f"IMGT {p['IMGT']}": f"{p['residue']} ({p['verdict']})"
                               for p in hm["positions"]},
                            "source": hm["source"]})
        except Exception as exc:
            st.warning(f"**{name}** — hallmarks unavailable: {exc}")
        for sp, g in germline_data(seq).items():
            if g:
                gl_rows.append({"name": name, "species": sp, "V gene": g["v_gene"],
                                "V identity (%)": g["v_identity"],
                                "J gene": g["j_gene"],
                                "J identity (%)": g["j_identity"]})

    if hm_rows:
        st.dataframe(pd.DataFrame(hm_rows), hide_index=True, **FIT)
        st.caption("A score of 4/4 is the canonical camelid profile. IMGT 49 is "
                   "shared with VH3, so 3/4 with a G at 49 is common and still a "
                   "genuine VHH. Engineered or humanised constructs sit between the "
                   "two profiles.")

    if gl_rows:
        st.subheader("Closest germline")
        st.dataframe(pd.DataFrame(gl_rows), hide_index=True, **FIT)
        st.caption("Identity to the closest human germline is a first-pass humanness "
                   "indicator, not an immunogenicity prediction — a sequence can be "
                   "germline-like and still contain a T-cell epitope.")
    elif not hm_rows:
        st.info("Germline assignment needs ANARCI, which is not available here.")

# --- Liabilities ------------------------------------------------------------ #
with tabs["Liabilities"]:
    st.markdown(
        "Chemical and manufacturability liabilities, scored and placed along the "
        "sequence.")

    include_added = st.checkbox(
        "Include the motifs added on top of the in-house scale",
        value=True,
        help="Fragmentation (DP), RGD and hydrophobic patches are not part of the "
             "original scoring scale; untick to reproduce its totals exactly.")

    liab_by_seq, prot_by_seq, summaries = {}, {}, []
    for name, seq in sequences.items():
        liab_by_seq[name] = liabilities_for(seq, annotations.get(name), include_added)
        prot_by_seq[name] = proteases_for(seq, annotations.get(name))
        summaries.append(vc.liability_summary(name, liab_by_seq[name],
                                              prot_by_seq[name]))

    st.dataframe(pd.DataFrame(summaries), hide_index=True, **FIT,
                 column_config={
                     "total_score": st.column_config.NumberColumn("Total score",
                                                                   format="%.1f"),
                     "n_liabilities": st.column_config.NumberColumn("Motifs"),
                     "n_high": st.column_config.NumberColumn("High"),
                     "n_medium": st.column_config.NumberColumn("Medium"),
                     "n_low": st.column_config.NumberColumn("Low"),
                     "n_in_CDR": st.column_config.NumberColumn("In a CDR"),
                     "protease_score": st.column_config.NumberColumn(
                         "Protease score", format="%.2f"),
                     "n_protease_sites": st.column_config.NumberColumn(
                         "Protease sites")})

    img_map = png_liability_map(liab_by_seq, annotations,
                                {n: len(s) for n, s in sequences.items()})
    if img_map:
        st.image(img_map, **FIT)
    img_profile = png_liability_profile(liab_by_seq)
    if img_profile:
        st.image(img_profile, **FIT)
    plt.close("all")

    st.caption(
        "Scores come from the in-house scale: N-glycosylation sequon 15, NG and DG "
        "10, NS/NT/NN and DS/DT/DD 5, free cysteine 8, Met 4, Trp 3; Low below 2, "
        "Medium below 6, High at or above 6. In its original form each motif is also "
        "weighted by its solvent exposure and the model confidence taken from a 3D "
        "structure. **This section is sequence-only, so neither factor is applied** — "
        "an NG scores 10 whether it is buried or fully exposed. The practical "
        "consequence is that the categories compress towards the top, and comparing "
        "sequences with each other tells you more than the category of any single "
        "motif.")

    st.subheader("Per sequence")
    for name in sequences:
        liab = liab_by_seq[name]
        summary = next(s for s in summaries if s["name"] == name)
        with st.container(border=True):
            st.markdown(f"**{name}** — score {summary['total_score']:.0f}, "
                        f"{summary['n_liabilities']} motif(s), "
                        f"{summary['n_high']} high, {summary['n_in_CDR']} in a CDR")
            if liab.empty:
                st.markdown("No liability motif detected.")
            else:
                st.dataframe(
                    liab[["liability", "motif", "range", "region", "in_CDR", "score",
                          "risk", "source", "note"]],
                    hide_index=True, **FIT,
                    column_config={
                        "range": st.column_config.TextColumn("Position"),
                        "score": st.column_config.NumberColumn("Score",
                                                                format="%.1f"),
                        "note": st.column_config.TextColumn("Note", width="large")})

    with st.expander("Predicted protease cleavage sites"):
        st.markdown(
            "Every protein has tryptic sites, so the count on its own means little. "
            "What is worth looking at is a high-scoring recognition site — Factor Xa "
            "or enterokinase — or a cluster of sites inside a CDR.")
        for name in sequences:
            prot = prot_by_seq[name]
            st.markdown(f"**{name}** — {len(prot)} site(s), "
                        f"score {prot['score'].sum():.2f}" if not prot.empty
                        else f"**{name}** — none")
            if not prot.empty:
                in_cdr = prot[prot["in_CDR"]]
                st.caption(f"{len(in_cdr)} inside a CDR")
                st.dataframe(prot.head(40), hide_index=True, **FIT,
                             column_config={"note": st.column_config.TextColumn(
                                 "Note", width="large")})

    st.download_button(
        "Liabilities (CSV)",
        pd.concat([df.assign(sequence=n) for n, df in liab_by_seq.items()
                   if not df.empty], ignore_index=True).to_csv(index=False)
        if any(not d.empty for d in liab_by_seq.values()) else "",
        file_name="liabilities.csv", mime="text/csv")

# --- Alignment ------------------------------------------------------------- #
with tabs["Alignment"]:
    if len(sequences) < 2:
        st.info("Submit at least two sequences to see an alignment.")
    else:
        with st.spinner("Aligning…"):
            # pymsaviz reads from a path, so the aligned FASTA is written to a
            # temporary directory that is removed as soon as the figure exists.
            img_msa, aln_text = png_msa(sequences, annotations)
        st.image(img_msa, **FIT)
        st.caption("FAMSA alignment rendered with pymsaviz. CDR bands refer to the "
                   "first sequence submitted.")
        st.download_button("Download alignment (FASTA)", aln_text,
                           file_name="alignment.fasta", mime="text/plain")

        st.subheader("Pairwise identity")
        cols_id = st.columns(2)
        with cols_id[0]:
            st.image(png_identity(sequences, None, None), **FIT)
        if annotations and any(a["CDR3"]["seq"] for a in annotations.values()):
            img_id3 = png_identity(sequences, annotations, "CDR3")
            if img_id3:
                with cols_id[1]:
                    st.image(img_id3, **FIT)
        st.caption("Two clones from the same panning can share 90%+ full-length "
                   "identity while having unrelated CDR3s — the second matrix is the "
                   "one that tells you whether they are really the same binder.")
        plt.close("all")

# --- Clonal families ------------------------------------------------------- #
with tabs["Clonal families"]:
    st.markdown(
        "After a panning, the first question is not how good each clone is but "
        "**how many different clones you actually have** — a plate of picked "
        "colonies is routinely the same few lineages over and over. Grouping on "
        "CDR3 answers that, because CDR3 is what differs between independent "
        "lineages while the frameworks stay nearly constant.")

    if len(sequences) < 2:
        st.info("Submit at least two sequences to group them into families.")
    elif not any((a.get("CDR3") or {}).get("seq") for a in annotations.values()):
        st.warning("No CDR3 could be annotated, so there is nothing to group on.")
    else:
        fam_cols = st.columns([2, 1, 1])
        with fam_cols[0]:
            cutoff = st.slider(
                "CDR3 identity threshold (%)", min_value=50, max_value=100,
                value=int(vc.FAMILY_DEFAULTS["cutoff"]), step=1,
                help="Two sequences at or above this identity belong to the same "
                     "family. 100% groups only strictly identical CDR3s.")
        with fam_cols[1]:
            linkage_method = st.selectbox(
                "Linkage", ["complete", "average"],
                index=0 if vc.FAMILY_DEFAULTS["linkage"] == "complete" else 1,
                help="complete: every pair inside a family meets the threshold. "
                     "average: merges on the mean distance, more permissive — a "
                     "family may then contain pairs below the threshold.")
        with fam_cols[2]:
            same_length = st.checkbox(
                "Same CDR3 length only",
                value=bool(vc.FAMILY_DEFAULTS["require_same_length"]),
                help="Forces CDR3s of different lengths apart whatever their "
                     "alignment. Off by default: a one-residue indel is a "
                     "plausible sibling, not a new lineage.")

        with st.spinner("Grouping…"):
            fam = families_for(sequences, annotations, float(cutoff),
                               linkage_method, bool(same_length))

        n_seqs = len(fam["assignment"])
        n_fam = fam["n_families"]
        singletons = int((fam["families"]["n_sequences"] == 1).sum())
        head = st.columns(4)
        head[0].metric("Sequences grouped", n_seqs)
        head[1].metric("Families", n_fam)
        head[2].metric("Singletons", singletons)
        head[3].metric("Redundancy", f"{n_seqs / n_fam:.2f}" if n_fam else "—",
                       help="Sequences per family. Above ~2 you are re-picking "
                            "the same lineages and could screen more colonies.")

        if fam["skipped"]:
            st.warning("Left out for want of an annotated CDR3: "
                       + ", ".join(fam["skipped"]))

        st.subheader("Families")
        st.dataframe(fam["families"], hide_index=True, **FIT)
        st.caption(
            "`min_identity` is the **worst** pair inside the family: with "
            "*complete* linkage it can never fall below your threshold, so a "
            "value under it means *average* linkage merged across it. "
            "`representative` is the family's most frequent CDR3 — a reasonable "
            "default pick for the clone to carry forward.")

        img_dendro = png_family_dendrogram(sequences, annotations, float(cutoff),
                                           linkage_method, bool(same_length))
        if img_dendro:
            st.subheader("Dendrogram")
            st.image(img_dendro, **FIT)
            st.caption("The dashed line is your threshold; branches joining to "
                       "its right are separate families. Move the slider and "
                       "watch which groups merge — that is the honest way to see "
                       "how robust a family is.")

        st.subheader("Per-sequence assignment")
        st.dataframe(fam["assignment"], hide_index=True, **FIT)

        dl = st.columns(3)
        dl[0].download_button("Families (CSV)",
                              fam["families"].to_csv(index=False),
                              file_name=f"families_{int(cutoff)}pct.csv",
                              mime="text/csv")
        dl[1].download_button("Assignment (CSV)",
                              fam["assignment"].to_csv(index=False),
                              file_name=f"family_assignment_{int(cutoff)}pct.csv",
                              mime="text/csv")
        dl[2].download_button("CDR3 identity matrix (CSV)",
                              fam["identity"].round(1).to_csv(),
                              file_name="cdr3_identity_matrix.csv",
                              mime="text/csv")

        with st.expander("One representative per family, as FASTA"):
            reps = {r: sequences[r] for r in fam["families"]["representative"]
                    if r in sequences}
            fasta = "".join(f">{k}\n{v}\n" for k, v in reps.items())
            st.code(fasta or "—", language="text")
            st.download_button("Download representatives (FASTA)", fasta,
                               file_name="family_representatives.fasta",
                               mime="text/plain")
        plt.close("all")


# --- DNA -------------------------------------------------------------------- #
with tabs["DNA"]:
    st.markdown("Back-translation into DNA using *E. coli* preferred codons, with "
                "cloning sites and long homopolymer runs removed.")

    col_a, col_b = st.columns([2, 1])
    with col_a:
        avoid = st.multiselect(
            "Restriction sites to avoid", sorted(vc.DEFAULT_AVOID_SITES),
            default=sorted(vc.DEFAULT_AVOID_SITES),
            help="A codon creating one of these sites is swapped for a synonymous "
                 "one.")
    with col_b:
        add_stop = st.checkbox("Append a stop codon (TAA)", value=True)
        gc_max = st.slider("Maximum GC fraction", 0.45, 0.70, 0.60, 0.01)

    sites = {k: vc.DEFAULT_AVOID_SITES[k] for k in avoid}
    fasta_parts = []
    for name, seq in sequences.items():
        res = dna_for(seq, tuple(sorted(sites.items())), add_stop, gc_max)
        ok = vc.translate_dna(res["dna"]) == seq
        fasta_parts.append(f">{name}_optimized\n" +
                           "\n".join(res["dna"][i:i + 60]
                                     for i in range(0, len(res["dna"]), 60)))
        with st.container(border=True):
            st.markdown(f"**{name}** — {res['length_nt']} nt, "
                        f"GC {res['gc_percent']} %, "
                        f"{len(res['repairs'])} codon(s) swapped for sites, "
                        f"{res['gc_swaps']} for GC balance")
            if not ok:
                st.error("Round-trip check failed — do not order this sequence.")
            if res["remaining_issues"]:
                st.warning("Could not remove: " + ", ".join(res["remaining_issues"]))
            st.code(res["dna"], language="text")
            if res["repairs"]:
                with st.expander(f"{len(res['repairs'])} codon swap(s)"):
                    st.dataframe(pd.DataFrame(res["repairs"]), hide_index=True, **FIT)

    st.download_button("Download DNA (FASTA)", "\n".join(fasta_parts) + "\n",
                       file_name="vhh_optimized_dna.fasta", mime="text/plain",
                       type="primary")
    st.caption(
        "Every output is verified by translating it back and comparing to the input. "
        "This optimiser maximises codon frequency, balances GC and removes "
        "cloning-unfriendly motifs — it does not model mRNA secondary structure near "
        "the ribosome binding site, codon-pair bias or tRNA availability. Run the "
        "sequence through your synthesis vendor's checker before ordering.")
    st.button("🧬 Use these in the cloning section", key="to_cloning", **FIT,
              on_click=go_to, args=("cloning",))

# --- Download -------------------------------------------------------------- #
with tabs["Download"]:
    st.markdown("Export the full results table, including every region sequence, "
                "custom feature and flag.")
    st.download_button("Results (CSV)", df.to_csv(index=False),
                       file_name="vhh_results.csv", mime="text/csv", type="primary")
    st.download_button("Secondary structure (PNG)", png_secondary_structure(df),
                       file_name="secondary_structure.png", mime="image/png")
    st.download_button("Charge (PNG)", png_charge(df, sequences, ph),
                       file_name="charge.png", mime="image/png")
    st.download_button("Profiles (PNG)", png_profiles(sequences, annotations, window),
                       file_name="profiles.png", mime="image/png")

    st.divider()
    st.markdown("### Full PDF report")
    st.markdown("Everything in one document: properties, developability flags, "
                "regions, hallmarks, germline, liabilities and all figures. Building "
                "it takes a few seconds, so it is generated on demand.")

    include_cf = st.checkbox(
        "Include the per-sequence secondary-structure plots",
        value=len(sequences) <= 6,
        help="One extra page per sequence. With many sequences this is most of the "
             "file size and most of the build time.")

    if st.button("Build PDF report", type="secondary"):
        # The displayed plots are cached as PNG bytes, so the report rebuilds real
        # figure objects — only when the button is pressed.
        figures = {
            "Secondary structure composition": vc.plot_secondary_structure(df),
            "Charge": vc.plot_charge(df, sequences, ph=ph),
            "Charge and hydropathy profiles": vc.plot_profiles(
                sequences, annotations, window=window),
        }
        if include_cf:
            for cf_name, cf_seq in sequences.items():
                figures[f"Predicted secondary structure — {cf_name}"] = \
                    vc.plot_chou_fasman(cf_seq, cf_name, window=cf_window,
                                        annotation=annotations.get(cf_name))
        if len(sequences) >= 2:
            try:
                with tempfile.TemporaryDirectory() as tmp:
                    fig_msa_pdf, _ = vc.plot_msa(sequences, annotations, workdir=tmp)
                figures["Multiple sequence alignment"] = fig_msa_pdf
                figures["Pairwise identity"] = vc.plot_identity(
                    vc.identity_matrix(sequences), "Full-length identity (%)")
                mat3 = vc.identity_matrix(sequences, annotations, region="CDR3")
                if len(mat3) >= 2:
                    figures["CDR3 identity"] = vc.plot_identity(mat3,
                                                                "CDR3 identity (%)")
            except Exception as exc:
                st.warning(f"Alignment figures skipped: {exc}")

        with st.spinner("Building the report…"):
            with tempfile.TemporaryDirectory() as tmp:
                pdf_path = vc.build_pdf_report(
                    f"{tmp}/vhh_report.pdf", df, sequences, annotations,
                    figures=figures, ph=ph, version=VERSION, author=AUTHOR,
                    app_name=APP_NAME)
                st.session_state["pdf_bytes"] = pdf_path.read_bytes()
        plt.close("all")

    if st.session_state.get("pdf_bytes"):
        st.download_button("Download report (PDF)", st.session_state["pdf_bytes"],
                           file_name="vhh_report.pdf", mime="application/pdf",
                           type="primary")

    st.divider()
    st.dataframe(df.drop(columns=["sequence"]), hide_index=True, **FIT)

footer()
