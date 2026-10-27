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
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import cloning as cl
import generator as gn
import constructs as ctm
import repertoire as rep
import triage as tg
import codons as cdn
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
VERSION = "20261027.1200"

# Wall-clock start of this script run. Streamlit re-executes the whole file on
# every interaction, so this is set afresh each time and measures what the user
# actually waits for — not just the analysis loop.
_RUN_STARTED = time.perf_counter()
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
REQUIRED_MODULES = {"vhh_core": "2026.10.27a", "sanger": "2026.10.21b",
                    "cloning": "2026.10.21b", "generator": "2026.10.21b",
                    "constructs": "2026.10.21d",
                    "repertoire": "2026.10.23a",
                    "triage": "2026.10.26a", "codons": "2026.10.27b"}
_stale = [f"`{name}.py` (has {getattr(mod, '__version__', 'no version')}, "
          f"needs {want})"
          for name, want, mod in (("vhh_core", REQUIRED_MODULES["vhh_core"], vc),
                                  ("sanger", REQUIRED_MODULES["sanger"], sg),
                                  ("cloning", REQUIRED_MODULES["cloning"], cl),
                                  ("generator", REQUIRED_MODULES["generator"], gn),
                                  ("constructs", REQUIRED_MODULES["constructs"], ctm),
                                  ("repertoire", REQUIRED_MODULES["repertoire"], rep),
                                  ("triage", REQUIRED_MODULES["triage"], tg),
                                  ("codons", REQUIRED_MODULES["codons"], cdn))
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


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def classify_cached(seq: str, method: str):
    """What kind of sequence this is — ANARCI numbering, so worth caching."""
    return vc.classify(seq, method=method)


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


def analyze_with_progress(seqs: dict[str, str], ph: float, method: str,
                          slot=None):
    """Run the analysis sequence by sequence, showing progress as it goes.

    ``slot`` is a placeholder created earlier in the script — normally the one
    just under the Analyze button. Without it the bar would be drawn wherever
    this happens to be called, which is below the generator panel and nowhere
    near the control the user just clicked.

    The bar is replaced by a one-line result rather than cleared: on a handful
    of sequences the analysis takes well under a second and a bar that vanishes
    leaves no evidence anything happened at all — and once the cache is warm it
    is instantaneous, so there is nothing to see unless something is left behind.
    """
    target = slot if slot is not None else st.empty()
    st.session_state["_status_slot"] = target
    bar = target.progress(0.0, text="Preparing…")
    rows, anns, errs = [], {}, {}
    total = len(seqs)
    started = time.perf_counter()
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
        elapsed = time.perf_counter() - started
        plural = "s" if total != 1 else ""
        # Named precisely: this covers the characterisation and annotation loop
        # only. Every tab body still has to render afterwards — the alignment,
        # the logo, the identity matrices, the plasmid maps — and that is where
        # most of the wall-clock time goes and why Streamlit's own spinner keeps
        # turning after this bar is full. The final line is written at the end of
        # the script, once the page is actually finished.
        #
        # No "from cache" note either: whether a result came from the cache
        # cannot be read off the elapsed time, and guessing it would mislabel a
        # genuinely fast run.
        target.info(
            f"{total} sequence{plural} characterised in "
            + (f"{elapsed * 1000:.0f} ms" if elapsed < 1 else f"{elapsed:.1f} s")
            + " — drawing the figures…", icon="⏳")
        st.session_state["_status_parts"] = (total, elapsed)
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


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def aggregation_for(seq, ann, window, threshold):
    profile = vc.aggregation_profile(seq, window=window)
    regions = vc.aggregation_regions(profile, ann, threshold=threshold)
    return (profile, regions, vc.aggregation_summary(profile, regions),
            vc.weakest_windows(profile, ann))


@st.cache_data(ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES, show_spinner=False)
def png_aggregation(seq, ann, window, threshold, title):
    profile, _r, _s, _w = aggregation_for(seq, ann, window, threshold)
    return fig_to_png(vc.plot_aggregation(profile, ann, threshold=threshold,
                                          title=title), close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def humanness_for(seq, family):
    return vc.humanness(seq, family=family)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def disulfides_for(seq, ann):
    return vc.disulfides(seq, ann)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def repairs_for(seq, liabilities, ann, include_cdr):
    return vc.suggest_repairs(seq, liabilities, ann, include_cdr=include_cdr)


@st.cache_data(ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES, show_spinner=False)
def png_construct(construct, title):
    return fig_to_png(ctm.plot_construct(construct, title), close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def rarity_for(seq, threshold, min_support):
    return rep.rare_residues(seq, threshold=threshold, min_support=min_support)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def rarity_tables(seqs, threshold, min_support):
    return rep.rarity_table(seqs, threshold=threshold, min_support=min_support)


@st.cache_data(**_FIG_CACHE)
def png_rarity(seq, name, threshold, min_support):
    return fig_to_png(rep.plot_rarity(seq, name, threshold=threshold,
                                      min_support=min_support), close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def repertoire_stats(seqs, anns):
    per_seq, distribution = rep.cdr_length_comparison(anns)
    composition = rep.composition_comparison(seqs, anns)
    return per_seq, distribution, composition


@st.cache_data(**_FIG_CACHE)
def png_cdr_lengths(seqs, anns):
    _p, distribution, _c = repertoire_stats(seqs, anns)
    return fig_to_png(rep.plot_cdr_lengths(distribution, len(seqs)), close=True)


@st.cache_data(**_FIG_CACHE)
def png_composition(seqs, anns, regions):
    _p, _d, composition = repertoire_stats(seqs, anns)
    return fig_to_png(rep.plot_composition(composition, list(regions)), close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def comparison_for(reference, seqs, reference_name):
    return vc.compare_to_reference(reference, seqs, reference_name)


@st.cache_data(**_FIG_CACHE)
def png_differences(reference, seqs, reference_name):
    return fig_to_png(vc.plot_differences(
        comparison_for(reference, seqs, reference_name)), close=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def triage_for(seqs, anns, kinds, ph, include_added):
    return tg.triage_table(seqs, anns, kinds, ph=ph, include_added=include_added)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def diagnose_gene(dna, expected, cluster_window, gc_window):
    return cdn.diagnose(dna, expected_protein=expected or None,
                        cluster_window=cluster_window, gc_window=gc_window)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def gc_fix_for(dna, gc_window, strategy, avoid):
    sites = {k: v for k, v in vc.DEFAULT_AVOID_SITES.items() if k in avoid}
    return cdn.suggest_gc_fix(dna, avoid_sites=sites, gc_window=gc_window,
                              strategy=strategy)


@st.cache_data(**_FIG_CACHE)
def png_gene_diagnosis(dna, expected, cluster_window, gc_window):
    result = diagnose_gene(dna, expected, cluster_window, gc_window)
    return fig_to_png(cdn.plot_diagnosis(result, gc_window), close=True)


def png_identity(seqs, anns, region):
    mat = vc.identity_matrix(seqs, anns, region=region) if region \
        else vc.identity_matrix(seqs)
    if len(mat) < 2:
        return None
    title = f"{region} identity (%)" if region else "Full-length identity (%)"
    return fig_to_png(vc.plot_identity(mat, title), close=True)


@st.cache_data(**_FIG_CACHE)
def png_msa(seqs, anns, wrap_length=60):
    with tempfile.TemporaryDirectory() as tmp:
        fig, aln_path = vc.plot_msa(seqs, anns, workdir=tmp,
                                    wrap_length=wrap_length)
        return fig_to_png(fig, close=True), aln_path.read_text()


@st.cache_data(**_FIG_CACHE)
def png_logo(seqs, anns, to_type, wrap_length, color_scheme):
    with tempfile.TemporaryDirectory() as tmp:
        fig, matrix, _path = vc.plot_logo(seqs, anns, to_type=to_type,
                                          wrap_length=wrap_length,
                                          color_scheme=color_scheme, workdir=tmp)
        return fig_to_png(fig, close=True), matrix


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
    "genes": ("🔧", "Gene diagnosis",
              "Paste coding sequences that are already in a plasmid and ask what "
              "might be wrong with them: frame and stop checks, documented rare "
              "codons and their clusters, GC stretches — plus synonymous swaps "
              "that lower GC without touching the protein."),
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
if "pending_gene_dna" in st.session_state:
    st.session_state["gene_dna"] = st.session_state.pop("pending_gene_dna")
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

    # Two rows of two rather than one row of four: four cards side by side are
    # too narrow for their descriptions to be readable.
    entries = list(SECTIONS.items())
    for row_start in range(0, len(entries), 2):
        for col, (key, (icon, label, desc)) in zip(
                st.columns(2), entries[row_start:row_start + 2]):
            with col, st.container(border=True):
                st.markdown(f"### {icon} {label}")
                st.markdown(desc)
                st.button(f"Open {label.lower()}", key=f"home_{key}",
                          type="primary", **FIT, on_click=go_to, args=(key,))

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

@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def gene_summary(seqs, gc_window, cluster_window):
    return cdn.summarise(seqs, gc_window=gc_window, cluster_window=cluster_window)


def render_gene_detail(name, gene_dna, diag, expected, cluster_window, gc_window):
    """The per-sequence view: findings, figure, tables, and the GC suggestions."""
    orf = diag["orf"]
    gm = st.columns(5)
    gm[0].metric("Length", f"{orf['length_nt']} nt")
    gm[1].metric("Codons", orf["n_codons"])
    gm[2].metric("GC", f"{diag['gc_overall']:.1%}")
    gm[3].metric("Rare codons", diag["n_rare"],
                 help="Codons whose scarcity in E. coli is well established and "
                      "supplemented by the Rosetta and CodonPlus strains.")
    gm[4].metric("Rare clusters", len(diag["clusters"]))

    st.subheader("Findings")
    renderers = {"error": st.error, "warning": st.warning,
                 "ok": st.success, "info": st.info}
    for level, message in diag["findings"]:
        renderers.get(level, st.info)(message)

    st.image(png_gene_diagnosis(gene_dna, expected, cluster_window, gc_window),
             **FIT)

    if len(diag["clusters"]):
        st.subheader("Rare-codon clusters")
        st.dataframe(diag["clusters"], hide_index=True, **FIT)
        st.caption("Clustered and tandem rare codons are the ones associated "
                   "with ribosome stalling. Isolated ones are generally "
                   "tolerated, which is why they are counted but not listed.")

    if len(diag["gc_regions"]):
        st.subheader("Extreme GC stretches")
        st.dataframe(diag["gc_regions"], hide_index=True, **FIT)

    gc_rich = (diag["gc_regions"][diag["gc_regions"]["direction"] == "GC-rich"]
               if len(diag["gc_regions"]) else [])
    if len(gc_rich):
        st.subheader("Suggested codon changes")
        st.markdown(
            "A GC-rich stretch is one of the few sequence problems with a clean "
            "fix: swap codons for lower-GC synonyms. Every suggestion below has "
            "been checked — it is not a documented rare codon, it creates none "
            "of the restriction sites you select, it opens no homopolymer run, "
            "and the protein is verified unchanged.")
        fx = st.columns([1, 2])
        with fx[0]:
            strategy = st.radio(
                "How far to go", ["minimal", "maximal"], key=f"strat_{name}",
                format_func=lambda x: ("Minimal — stop once the band is cleared"
                                       if x == "minimal"
                                       else "Maximal — swap everything possible"),
                help="Minimal applies the biggest reductions first and stops as "
                     "soon as no GC-rich stretch survives. You are paying to "
                     "synthesise this gene, and a construct differing in twenty "
                     "places is easier to reason about than one differing in "
                     "sixty.")
        with fx[1]:
            fix_avoid = st.multiselect(
                "Restriction sites to preserve", sorted(vc.DEFAULT_AVOID_SITES),
                default=sorted(vc.DEFAULT_AVOID_SITES), key=f"gcfix_{name}",
                help="No suggestion will create one of these.")

        fix = gc_fix_for(gene_dna, gc_window, strategy, tuple(sorted(fix_avoid)))
        fm = st.columns(4)
        fm[0].metric("Codons changed", fix["n_changes"])
        fm[1].metric("GC", f"{fix['gc_after']:.1%}",
                     delta=f"{(fix['gc_after'] - fix['gc_before']) * 100:+.1f} pts",
                     delta_color="off")
        fm[2].metric("GC-rich stretches left",
                     int((fix["remaining_regions"]["direction"] == "GC-rich").sum())
                     if len(fix["remaining_regions"]) else 0)
        fm[3].metric("Protein",
                     "unchanged" if fix["protein_unchanged"] else "CHANGED")

        if not fix["protein_unchanged"]:
            st.error("The rewritten sequence does not translate to the same "
                     "protein. Do not use it, and please report it.")
        elif not fix["n_changes"]:
            st.info("No synonymous swap satisfies all the constraints here — "
                    "every lower-GC alternative is either a rare codon, creates "
                    "a site you asked to keep, or opens a homopolymer run.")
        else:
            st.dataframe(fix["regions"], hide_index=True, **FIT)
            st.dataframe(fix["suggestions"], hide_index=True, **FIT)
            st.markdown("**Rewritten sequence**")
            st.code(fix["dna"], language="text")
            st.download_button(
                "Rewritten DNA (FASTA)",
                f">{name}_gc_{strategy}\n{fix['dna']}\n",
                file_name=f"{name}_gc_{strategy}.fasta", mime="text/plain",
                key=f"dl_fix_{name}")
            st.caption("Lowering local GC removes one plausible obstacle to "
                       "expression. It does not promise the gene will then "
                       "express — the causes at the top of this page are all "
                       "still in play.")

    with st.expander("Every codon"):
        st.dataframe(diag["codons"], hide_index=True, **FIT)
        st.download_button("Codon table (CSV)", diag["codons"].to_csv(index=False),
                           file_name=f"{name}_codons.csv", mime="text/csv",
                           key=f"dl_codons_{name}")

    with st.expander("Translation"):
        st.code(orf["protein"] or "—", language="text")
        st.caption(f"{orf['protein_length']} residues.")

    st.caption(
        "**No CAI is reported.** A Codon Adaptation Index needs a full codon "
        "usage table for the host; this app does not ship one, and inventing the "
        "values would give a number that looks precise and means nothing.")


# --------------------------------------------------------------------------- #
# Gene diagnosis
# --------------------------------------------------------------------------- #

if section == "genes":
    st.title("🔧 Gene diagnosis")
    st.markdown(
        "The reverse of designing a gene: paste coding sequences that already "
        "exist — in a plasmid, from a vendor, from a colleague — and ask what "
        "might be wrong with them. Several at a time.")
    st.warning(
        "**Codon usage is rarely the first suspect.** A VHH that expresses badly "
        "in *E. coli* is more often explained by periplasmic toxicity, "
        "inefficient signal-peptide cleavage, failure to form the intradomain "
        "disulfide, inclusion bodies, plasmid loss or induction conditions. A "
        "gene with flawless codon usage expresses badly for all of those "
        "reasons. This section tells you whether the sequence is *part* of the "
        "problem; it cannot tell you that it is the problem.", icon="⚠️")

    gcol_in, gcol_act = st.columns([3, 1])
    with gcol_act:
        st.write("")
        if st.button("Load example", key="gene_example", **FIT):
            st.session_state["pending_gene_dna"] = (
                ">example\n"
                + vc.reverse_translate(
                    list(vc.parse_input(EXAMPLE).values())[0])["dna"]
                + "\n")
            st.rerun()
        if st.button("🗑️ Clear", key="gene_clear", **FIT,
                     disabled=not st.session_state.get("gene_dna", "").strip()):
            st.session_state["gene_dna"] = ""
            st.rerun()
        gene_upload = st.file_uploader("or upload FASTA", type=["fasta", "fa", "txt"],
                                       key="gene_upload")
        if gene_upload is not None and gene_upload.file_id != st.session_state.get(
                "last_gene_upload"):
            st.session_state["last_gene_upload"] = gene_upload.file_id
            st.session_state["pending_gene_dna"] = gene_upload.read().decode(
                "utf-8", errors="replace")
            st.rerun()
        st.caption("FASTA, or one bare sequence, or several separated by blank "
                   "lines. Up to 20 sequences, 6000 nt each.")
        st.button("🔬 Design DNA from protein instead", key="genes_to_analysis",
                  **FIT, on_click=go_to, args=("analysis",),
                  help="The analysis section's DNA tab back-translates a protein "
                       "and can send the result here.")

    with gcol_in:
        gene_text = st.text_area("Coding DNA (5'→3')", height=200, key="gene_dna",
                                 placeholder=">clone_1\nCAGGTGCAGCTGCAGGAAAGCGGC…")
        st.button("Diagnose", type="primary", key="gene_run", **FIT)

    gset = st.columns(3)
    with gset[0]:
        g_cluster = st.slider("Cluster window (codons)", 2, 20,
                              int(cdn.CLUSTER_WINDOW), key="g_cluster")
    with gset[1]:
        g_gcwin = st.slider("GC window (nt)", 15, 120, int(cdn.GC_WINDOW), step=5,
                            key="g_gcwin")
    with gset[2]:
        g_expected = st.text_input(
            "Expected protein (optional)", key="g_expected",
            help="Compared against the translation of the sequence selected "
                 "below. Worth doing: a mismatch makes everything else beside "
                 "the point.")

    if not gene_text.strip():
        st.info("Paste one or more coding sequences to diagnose.")
    else:
        try:
            gene_seqs = cdn.parse_dna_input(gene_text)
        except ValueError as exc:
            st.error(f"Input error — {exc}")
            gene_seqs = {}

        if gene_seqs:
            st.subheader(f"Overview — {len(gene_seqs)} sequence(s)")
            overview = gene_summary(gene_seqs, int(g_gcwin), int(g_cluster))
            st.dataframe(overview, hide_index=True, **FIT)
            st.caption(
                "`in_frame` and `internal_stops` are the deterministic checks — "
                "wrong there means wrong. The rest is triage: `rare_clusters` "
                "matters more than `rare_codons`, and `gc_rich_stretches` is what "
                "the suggestions below act on.")
            st.download_button("Overview (CSV)", overview.to_csv(index=False),
                               file_name="gene_overview.csv", mime="text/csv")

            st.divider()
            picked = (list(gene_seqs)[0] if len(gene_seqs) == 1
                      else st.selectbox("Look at one in detail", list(gene_seqs),
                                        key="gene_pick"))
            gene_dna = gene_seqs[picked]
            try:
                diag = diagnose_gene(gene_dna, g_expected.strip(),
                                     int(g_cluster), int(g_gcwin))
            except ValueError as exc:
                st.error(f"Could not read that sequence — {exc}")
                diag = None

            if diag:
                render_gene_detail(picked, gene_dna, diag, g_expected.strip(),
                                   int(g_cluster), int(g_gcwin))
    footer()
    st.stop()


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
    "Physico-chemical characterization, CDR/framework annotation and alignment. "
    "Built around VHHs, but **any protein is accepted** — each sequence is "
    "classified, and the VHH-specific analyses only run on the ones that are "
    "VHHs. Paste one or more sequences, in FASTA or raw form.")

col_in, col_act = st.columns([3, 1])
with col_act:
    st.write("")
    if st.button("Load example", **FIT):
        st.session_state["seq_text"] = EXAMPLE
        st.rerun()
    if st.button("🗑️ Clear", **FIT,
                 help="Empty the sequence box.",
                 disabled=not st.session_state.get("seq_text", "").strip()):
        st.session_state["seq_text"] = ""
        # forget the last upload too, otherwise re-selecting the same file after
        # clearing does nothing: the uploader still considers it already handled
        st.session_state["last_upload"] = None
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
    # Filled in later by analyze_with_progress: the progress bar while it runs,
    # then the result line. Created here so it appears under the button rather
    # than wherever in the script the analysis happens to be called.
    status_slot = st.empty()

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

df, annotations, errors = analyze_with_progress(sequences, ph, method,
                                                slot=status_slot)

# What is each sequence? The app accepts any protein, so every VHH-specific
# analysis has to know which sequences it legitimately applies to.
kinds = {name: classify_cached(seq, method) for name, seq in sequences.items()}
# Heavy-chain V domains — VHH *and* conventional VH — get the full analysis.
# Almost none of it is camelid-specific: IMGT numbering, the CDR windows,
# germline assignment, humanness and CDR3 grouping are the same operation on a
# VH. The two camelid readings (the hallmark tetrad, the extra CDR1-CDR3
# disulfide) are labelled in place instead of gating the rest.
heavy_names = [n for n in sequences if kinds[n]["full_analysis"]]
vhh_names = [n for n in sequences if kinds[n]["is_vhh"]]
other_names = [n for n in sequences if not kinds[n]["full_analysis"]]
heavy_sequences = {n: sequences[n] for n in heavy_names}
heavy_annotations = {n: annotations[n] for n in heavy_names if n in annotations}

df.insert(1, "type", [kinds[n]["kind"] for n in df["name"]])

_vh_only = [n for n in heavy_names if not kinds[n]["is_vhh"]]
if _vh_only:
    st.info(
        f"**{len(_vh_only)} sequence(s) are conventional VH, not VHH.** They go "
        "through the full heavy-chain analysis all the same — CDR/framework, "
        "germline, humanness, clonal families — and are labelled as VH "
        "throughout. Only two readings are camelid-specific and are flagged "
        "where they appear: the FR2 hallmark tetrad, and the additional "
        "CDR1–CDR3 disulfide.\n\n"
        + "\n".join(f"- **{n}** — {kinds[n]['note']}" for n in _vh_only),
        icon="🧬")

if other_names:
    st.info(
        f"**{len(other_names)} of {len(sequences)} sequence(s) are not "
        "heavy-chain variable domains.** They are characterised like any "
        "protein — mass, pI, extinction, charge, aggregation, liabilities, DNA — "
        "but the antibody-specific tabs skip them.\n\n"
        + "\n".join(f"- **{n}** — {kinds[n]['kind']}: {kinds[n]['note']}"
                     for n in other_names), icon="ℹ️")

# An annotation failure on something that is not an antibody domain is the
# expected outcome, not an error worth alarming anyone about.
real_errors = {n: m for n, m in errors.items() if kinds[n]["full_analysis"]}
if real_errors:
    for name, msg in real_errors.items():
        st.warning(f"**{name}** — annotation failed: {msg}. "
                   "Physico-chemical values are still reported below.")


def heavy_only_notice(what: str) -> bool:
    """Guard for a panel that needs a heavy-chain V domain (VHH or VH)."""
    if not heavy_names:
        st.info(f"{what} applies to heavy-chain variable domains — VHH or "
                "conventional VH — and none of the submitted sequences is one. "
                "The generic tabs — Properties, Plots, Liabilities, Aggregation, "
                "DNA, Constructs — cover them all.")
        return False
    if other_names:
        st.caption(f"Showing the {len(heavy_names)} heavy-chain domain(s); "
                   f"{len(other_names)} other sequence(s) are not analysed here.")
    return True


def type_of(name: str) -> str:
    """Short label for a sequence's kind, for the per-sequence tables."""
    return "VHH" if kinds[name]["is_vhh"] else kinds[name]["kind"]

# Named rather than indexed: inserting a tab in the middle of an indexed list
# silently reassigns every panel after it.
TAB_LABELS = ["Triage", "Properties", "Plots", "CDR / Framework",
              "Hallmarks & germline", "Liabilities", "Aggregation", "Alignment",
              "vs Reference", "Repertoire", "Clonal families", "DNA",
              "Constructs", "Download"]
tabs = dict(zip(TAB_LABELS, st.tabs(TAB_LABELS)))

# --- Triage ----------------------------------------------------------------- #
with tabs["Triage"]:
    st.markdown(
        "Every flag the app can raise, one row per sequence. This is the tab for "
        "the moment you have to **decide** — which clones to carry forward, which "
        "ten of sixteen to order — rather than to study one candidate.")
    st.info(
        "**There is no overall score, on purpose.** Weighting a deamidation site "
        "against a hydrophobic patch against three divergent framework positions "
        "means choosing what matters, and that depends on what the VHH is for. A "
        "single number would hide that choice behind a decimal point. Sort by the "
        "column you care about today.", icon="ℹ️")

    tri_cols = st.columns([2, 1])
    with tri_cols[0]:
        shown_groups = st.multiselect(
            "Column groups", list(tg.COLUMN_GROUPS),
            default=[g for g in tg.COLUMN_GROUPS if g != "Physico-chemical"])
    with tri_cols[1]:
        flagged_only = st.checkbox(
            "Only sequences with a flag", value=False,
            help="Keeps the rows carrying a high-risk liability, an orphan "
                 "cysteine, an unusual framework residue or an "
                 "aggregation-prone region.")

    with st.spinner("Collecting…"):
        table = triage_for(sequences, annotations, kinds, ph, include_added=True)

    view = table
    if flagged_only and len(view):
        mask = view[tg.FLAG_COLUMNS].fillna(False).any(axis=1)
        view = view[mask]

    keep = [c for g in shown_groups for c in tg.COLUMN_GROUPS[g]]
    keep = [c for c in dict.fromkeys(["name"] + keep) if c in view.columns]
    st.dataframe(view[keep], hide_index=True, **FIT)
    st.caption(
        f"{len(view)} of {len(table)} sequence(s) shown. Empty cells are "
        "analyses that do not apply — a sequence that is not an antibody domain "
        "has no hallmark score, and its cysteines are not read against the VHH "
        "canonical pair. Click a column header to sort.")

    if len(table):
        flags = table[["name"] + tg.FLAG_COLUMNS].copy()
        n_flagged = int(flags[tg.FLAG_COLUMNS].fillna(False).any(axis=1).sum())
        fm = st.columns(4)
        fm[0].metric("Sequences", len(table))
        fm[1].metric("With a flag", n_flagged)
        fm[2].metric("High-risk liability",
                     int(flags["has_high_risk"].fillna(False).sum()))
        fm[3].metric("Orphan cysteine",
                     int(flags["has_orphan_cys"].fillna(False).sum()))
        with st.expander("The flag columns on their own"):
            st.dataframe(flags, hide_index=True, **FIT)

    st.download_button("Triage table (CSV)", table.to_csv(index=False),
                       file_name="triage.csv", mime="text/csv")


# --- Properties ------------------------------------------------------------ #
with tabs["Properties"]:
    main_cols = ["name", "type", "length", "mw_kda", "pI", "ext_coeff_reduced",
                 "A280_0.1pct", "instability_index", "aliphatic_index", "gravy",
                 "aromaticity", "net_charge", "n_Cys", "n_Lys",
                 "n_primary_amines"]
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

    with st.expander("Lysines and primary amines — for conjugation"):
        st.markdown(
            "`n_Lys` counts lysine side chains; `n_primary_amines` adds the free "
            "N-terminal α-amine, which an NHS ester labels just as readily. "
            "Together they are the sites a random amine-coupling chemistry — "
            "NHS-biotin, NHS-fluorophore, amine coupling onto a chip — will react "
            "with.")
        lys_rows = []
        for name, seq in sequences.items():
            table = vc.lysine_sites(seq, annotations.get(name))
            lys_rows.append(table.assign(sequence=name))
        lys_all = pd.concat(lys_rows, ignore_index=True) if lys_rows else pd.DataFrame()
        if len(lys_all):
            st.dataframe(
                lys_all[["sequence", "site", "position", "region", "in_CDR", "note"]],
                hide_index=True, **FIT,
                column_config={"note": st.column_config.TextColumn("Note",
                                                                   width="large")})
            in_cdr = lys_all[lys_all["in_CDR"]]
            if len(in_cdr):
                st.warning(
                    f"{len(in_cdr)} amine site(s) sit inside a CDR: "
                    + ", ".join(f"{r['sequence']} {r['site']}"
                                for _, r in in_cdr.iterrows())
                    + ". Random amine coupling will hit these along with the "
                      "rest, which is the usual explanation for a VHH that loses "
                      "binding after labelling. Site-specific chemistry — a "
                      "C-terminal cysteine, sortase, or an Avi tag — avoids the "
                      "problem entirely.", icon="⚠️")
            st.download_button("Amine sites (CSV)", lys_all.to_csv(index=False),
                               file_name="amine_sites.csv", mime="text/csv")
        st.caption(
            "Counted from sequence, so this is how many amines **exist**, not "
            "how many are solvent-accessible — that needs a structure. Treat it "
            "as the upper bound on labelling sites.")

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
    _cdr_ok = heavy_only_notice("CDR and framework annotation")
    if _cdr_ok and not heavy_annotations:
        st.warning("No sequence could be annotated.")
    elif _cdr_ok:
        rows = [{"name": n, "type": type_of(n),
                 "method": df.loc[df["name"] == n, "annotation_method"].iloc[0],
                 "CDR1": heavy_annotations[n]["CDR1"]["seq"],
                 "CDR2": heavy_annotations[n]["CDR2"]["seq"],
                 "CDR3": heavy_annotations[n]["CDR3"]["seq"],
                 "CDR1 len": len(heavy_annotations[n]["CDR1"]["seq"]),
                 "CDR2 len": len(heavy_annotations[n]["CDR2"]["seq"]),
                 "CDR3 len": len(heavy_annotations[n]["CDR3"]["seq"])}
                for n in heavy_annotations]
        st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT)

        st.subheader("Annotated sequences")
        st.markdown(
            '<span style="background:#e41a1c;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR1</span> '
            '<span style="background:#377eb8;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR2</span> '
            '<span style="background:#4daf4a;color:#fff;padding:1px 6px;'
            'border-radius:3px">CDR3</span>', unsafe_allow_html=True)
        st.markdown("".join(vc.render_annotated_html(n, heavy_annotations[n])
                            for n in heavy_annotations), unsafe_allow_html=True)

        with st.expander("Framework regions"):
            st.dataframe(pd.DataFrame(
                [{"name": n, **{r: heavy_annotations[n][r]["seq"]
                                for r in ("FR1", "FR2", "FR3", "FR4")}}
                 for n in heavy_annotations]), hide_index=True, **FIT)

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
    _has_vhh = heavy_only_notice("Hallmark and germline analysis")
    hm_rows, gl_rows = [], []
    for name, seq in heavy_sequences.items():
        try:
            hm = hallmark_data(seq, annotations.get(name))
            hm_rows.append({"name": name, "type": type_of(name),
                            "tetrad": hm["tetrad"],
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
                   "two profiles. **This is the one reading that is camelid-"
                   "specific**: on a conventional VH a score of 0/4 is simply what "
                   "a VH looks like, not a defect — the row is shown so you can "
                   "see the contrast, and the `type` column says which is which.")

    if gl_rows:
        st.subheader("Closest germline")
        st.dataframe(pd.DataFrame(gl_rows), hide_index=True, **FIT)
        st.caption("Identity to the closest human germline is a first-pass humanness "
                   "indicator, not an immunogenicity prediction — a sequence can be "
                   "germline-like and still contain a T-cell epitope.")
    elif not hm_rows:
        st.info("Germline assignment needs ANARCI, which is not available here.")

    st.divider()
    st.subheader("Humanness, position by position")
    st.markdown(
        "Percent identity to the closest human germline V gene, computed over "
        "**framework positions only** — CDR divergence from a human germline is "
        "expected and is not something you would ever try to fix. The divergent "
        "positions are listed because they are the candidate humanising "
        "mutations.")
    hcol = st.columns([1, 2])
    with hcol[0]:
        family = st.selectbox(
            "Germline family", ["IGHV3", "all human heavy V"], index=0,
            help="IGHV3 is the family VHHs belong to, so it is the fair "
                 "comparison. Searching all families rarely changes the answer.")
    fam_arg = None if family == "all human heavy V" else family

    for name, seq in heavy_sequences.items():
        try:
            h = humanness_for(seq, fam_arg)
        except Exception as exc:
            st.warning(f"**{name}** — humanness unavailable: {exc}")
            continue
        if not h["gene"]:
            st.warning(f"**{name}** — no germline could be compared.")
            continue
        with st.expander(f"**{name}** ({type_of(name)}) — "
                         f"{h['identity']}% to {h['gene']} "
                         f"({h['n_divergent']} divergent framework positions, "
                         f"{h['n_hallmark']} of them hallmarks)",
                         expanded=len(sequences) == 1):
            m = st.columns(3)
            m[0].metric("Framework identity", f"{h['identity']}%")
            m[1].metric("Closest human gene", h["gene"])
            m[2].metric("Divergent FR positions", h["n_divergent"])
            if len(h["by_region"]):
                st.dataframe(h["by_region"], hide_index=True, **FIT)
            if len(h["divergences"]):
                st.dataframe(h["divergences"], hide_index=True, **FIT)
                st.download_button(
                    f"Divergences — {name} (CSV)",
                    h["divergences"].to_csv(index=False),
                    file_name=f"humanness_{name}.csv", mime="text/csv",
                    key=f"dl_human_{name}")
            if not kinds[name]["is_vhh"]:
                st.caption(
                    "This is a conventional VH, so a high framework identity is "
                    "expected rather than informative — and the `hallmark` rows "
                    "below are not a humanisation decision here, since a VH "
                    "already carries the VH residues.")
            st.caption(
                "The rows flagged `hallmark` are the four FR2 positions that make "
                "a domain a VHH rather than a VH. They are the real decision: "
                "reverting them raises the humanness number and simultaneously "
                "removes the camelid surface that keeps a VHH soluble without a "
                "light chain. The other rows are the cheap ones.")

    st.info(
        "This is a first-pass humanness indicator and nothing more. "
        "Immunogenicity depends on T-cell epitopes, dose, route and formulation, "
        "none of which is in the sequence — a germline-like sequence can still "
        "carry an epitope, and a divergent one can be silent.", icon="⚠️")

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

    # ---- disulfides -------------------------------------------------------- #
    st.divider()
    st.subheader("Cysteines and disulfides")
    _ = heavy_only_notice("IMGT-numbered disulfide analysis")
    st.markdown(
        "The table above flags every cysteine outside the canonical FR1/FR3 pair "
        "as a possible free thiol. That is the right default, but it is wrong in "
        "one common case: camelid VHHs with long CDR3s often carry a **second, "
        "stabilising** disulfide, typically between CDR1 or FR2 and CDR3. This "
        "section tells the two situations apart. The canonical IMGT 23/104 pair "
        "is the same in a conventional VH, so that part applies to both; only "
        "the extra-pair reading is camelid.")
    for name, seq in heavy_sequences.items():
        try:
            ds = disulfides_for(seq, annotations.get(name))
        except Exception as exc:
            st.warning(f"**{name}** — cysteine analysis unavailable: {exc}")
            continue
        icon = "✅" if ds["n_orphan"] == 0 else "⚠️"
        with st.expander(f"{icon} **{name}** — {ds['n_cys']} cysteine(s): "
                         f"{ds['verdict']}", expanded=len(sequences) == 1):
            if len(ds["cysteines"]):
                st.dataframe(ds["cysteines"], hide_index=True, **FIT)
            if len(ds["pairs"]):
                st.dataframe(ds["pairs"], hide_index=True, **FIT)
            if ds["has_extra_pair"] and not kinds[name]["is_vhh"]:
                st.info(
                    "An extra cysteine pair was found, but the CDR1–CDR3 pattern "
                    "this is matched against is a **camelid** one, and this "
                    "sequence is a conventional VH. Treat the pairing as "
                    "unverified rather than as the known VHH motif.")
            elif ds["has_extra_pair"]:
                st.success(
                    "Read the extra cysteines as a stabilising pair rather than a "
                    "liability — but remember this is inferred from position "
                    "alone. Only a structure or a non-reducing gel confirms a bond.")
            elif ds["n_orphan"]:
                st.warning(
                    "An unpaired cysteine is a real risk: scrambling, dimerisation "
                    "through an intermolecular bond, and heterogeneity on "
                    "purification. The repair table below offers C→S.")

    # ---- suggested repairs ------------------------------------------------- #
    st.divider()
    st.subheader("Suggested repairs")
    st.markdown(
        "One conservative substitution per liability, so the table above becomes "
        "actionable. Substitutions inside a CDR are flagged separately because "
        "they risk affinity — a liability you keep is often the right call.")
    rep_cols = st.columns([1, 1, 2])
    with rep_cols[0]:
        rep_include_cdr = st.checkbox("Include CDR positions", value=True)
    rep_target = (list(sequences)[0] if len(sequences) == 1 else
                  rep_cols[1].selectbox("Sequence", list(sequences),
                                        key="repair_target"))

    rep_seq = sequences[rep_target]
    rep_ann = annotations.get(rep_target)
    repairs = repairs_for(rep_seq, liab_by_seq[rep_target], rep_ann,
                          bool(rep_include_cdr))
    doable = repairs[repairs["applicable"]] if len(repairs) else repairs
    refused = repairs[~repairs["applicable"]] if len(repairs) else repairs

    if not len(repairs):
        st.success("Nothing to repair — no liability with a known substitution.")
    else:
        st.dataframe(doable, hide_index=True, **FIT,
                     column_config={"rationale": st.column_config.TextColumn(
                         "Rationale", width="large")})
        if len(refused):
            st.caption("Offered no substitution, deliberately:")
            st.dataframe(refused[["liability", "position", "region", "rationale"]],
                         hide_index=True, **FIT,
                         column_config={"rationale":
                                        st.column_config.TextColumn("Reason",
                                                                    width="large")})

        # Several liabilities offer alternatives at the SAME position (a free
        # cysteine can go to S or to A). Selecting both by default would make
        # every one of them collide, so the default keeps the first — the
        # preferred — option per position.
        _seen_pos, _default = set(), []
        for _, r in doable.iterrows():
            if r["in_CDR"] or r["position"] in _seen_pos:
                continue
            _seen_pos.add(r["position"])
            _default.append(r["substitution"])
        picked = st.multiselect(
            "Apply substitutions", list(doable["substitution"]),
            default=_default,
            help="Defaults to one substitution per position, frameworks only — "
                 "the safe subset. Alternatives at the same position are listed "
                 "but not pre-selected, and picking two of them is refused "
                 "rather than silently resolved.")
        if picked:
            repaired, problems = vc.apply_repairs(rep_seq, picked)
            for p in problems:
                st.error(p)
            before = liab_by_seq[rep_target]["score"].sum()
            # the original annotation is reused on purpose: substitutions do not
            # move the region boundaries, and re-annotating would hide a repair
            # that accidentally broke an anchor.
            after_df = liabilities_for(repaired, rep_ann, include_added)
            rm = st.columns(3)
            rm[0].metric("Liability score", f"{after_df['score'].sum():.0f}",
                         delta=f"{after_df['score'].sum() - before:+.0f}",
                         delta_color="inverse")
            rm[1].metric("Motifs", len(after_df),
                         delta=len(after_df) - len(liab_by_seq[rep_target]),
                         delta_color="inverse")
            rm[2].metric("Substitutions", len(picked))
            st.code(repaired, language="text")
            st.download_button(
                "Repaired sequence (FASTA)",
                f">{rep_target}_repaired {' '.join(picked)}\n{repaired}\n",
                file_name=f"{rep_target}_repaired.fasta", mime="text/plain")
            st.caption(
                "The score after repair is recomputed with the **original** "
                "region boundaries, which is correct for substitutions but would "
                "not be if you changed the length. Re-run the full analysis on "
                "the repaired sequence before trusting anything else about it.")

# --- Aggregation ----------------------------------------------------------- #
with tabs["Aggregation"]:
    st.markdown(
        "Where along the sequence the aggregation-prone stretches are. A "
        "stretch becomes aggregation-prone when it is hydrophobic, β-prone and "
        "carries no charged residue to interrupt it — the three quantities this "
        "profile combines, plus a gatekeeper term for P/D/E/K/R.")
    st.warning(
        "**This is not CamSol, AGGRESCAN or TANGO.** Those fit their "
        "coefficients to experimental data; the weights here are a transparent "
        "judgement call, listed in `SOLUBILITY_WEIGHTS`. So read **where** the "
        "weak spots are and compare candidates against each other — do not read "
        "the absolute value as a predicted solubility. For a validated number, "
        "run the sequence through one of those servers.", icon="⚠️")

    agg_cols = st.columns(3)
    with agg_cols[0]:
        agg_window = st.slider("Smoothing window", 3, 15,
                               int(vc.DEFAULTS["profile_window"]), step=2,
                               help="Wider smoothing finds longer stretches and "
                                    "hides short ones.")
    with agg_cols[1]:
        agg_threshold = st.slider("Aggregation-prone threshold", -1.0, 0.5,
                                  float(vc.APR_THRESHOLD), step=0.05,
                                  help="A window scoring below this is called "
                                       "aggregation-prone. Negative leans prone, "
                                       "positive leans soluble.")
    with agg_cols[2]:
        agg_show_all = st.checkbox("Plot every sequence", value=len(sequences) <= 4,
                                   help="Untick with many sequences; the summary "
                                        "table is computed either way.")

    agg_summary, agg_regions_all = [], []
    for name, seq in sequences.items():
        profile, regions, summary, weakest = aggregation_for(
            seq, annotations.get(name), int(agg_window), float(agg_threshold))
        agg_summary.append({"name": name, **summary})
        if len(regions):
            agg_regions_all.append(regions.assign(sequence=name))

    st.dataframe(pd.DataFrame(agg_summary), hide_index=True, **FIT,
                 column_config={
                     "mean_score": st.column_config.NumberColumn("Mean", format="%.3f"),
                     "min_score": st.column_config.NumberColumn("Worst", format="%.3f"),
                     "worst_position": st.column_config.NumberColumn("At residue"),
                     "n_regions": st.column_config.NumberColumn("Regions"),
                     "residues_in_regions": st.column_config.NumberColumn("Residues"),
                     "fraction_prone": st.column_config.NumberColumn(
                         "Fraction", format="%.2f")})
    st.caption("`Mean` higher is better. `Worst` is the single lowest-scoring "
               "window, which is usually what decides whether a candidate "
               "behaves — one bad patch is enough.")

    for name, seq in sequences.items():
        profile, regions, summary, weakest = aggregation_for(
            seq, annotations.get(name), int(agg_window), float(agg_threshold))
        with st.expander(f"**{name}** — worst window {summary['min_score']:+.2f} "
                         f"at residue {summary['worst_position']}, "
                         f"{summary['n_regions']} region(s) below threshold",
                         expanded=agg_show_all and len(sequences) <= 4):
            if agg_show_all:
                st.image(png_aggregation(seq, annotations.get(name),
                                         int(agg_window), float(agg_threshold),
                                         f"{name} — aggregation-prone regions"),
                         **FIT)
            st.markdown("**Weakest windows**")
            st.dataframe(weakest, hide_index=True, **FIT)
            if len(regions):
                st.markdown("**Regions below the threshold**")
                st.dataframe(regions, hide_index=True, **FIT)
            else:
                st.success("No window crosses the threshold. The weakest-window "
                           "table above is still the one to compare candidates on.")

    dl_agg = st.columns(2)
    dl_agg[0].download_button(
        "Summary (CSV)", pd.DataFrame(agg_summary).to_csv(index=False),
        file_name="aggregation_summary.csv", mime="text/csv")
    dl_agg[1].download_button(
        "Regions (CSV)",
        pd.concat(agg_regions_all, ignore_index=True).to_csv(index=False)
        if agg_regions_all else "",
        file_name="aggregation_regions.csv", mime="text/csv",
        disabled=not agg_regions_all)
    plt.close("all")


# --- Alignment ------------------------------------------------------------- #
with tabs["Alignment"]:
    # Aligning a 1200-residue protein against 120-residue VHHs produces a mostly
    # empty figure, so the default is the VHH subset — unless there are no VHHs,
    # in which case aligning the others is the only thing that makes sense.
    if other_names and heavy_names:
        align_all = st.checkbox(
            f"Include the {len(other_names)} sequence(s) that are not antibody "
            "domains", value=False,
            help="Off by default: mixing sequences of very different lengths "
                 "gives an alignment that is mostly gaps, and a logo to match. "
                 "VHH and conventional VH align together perfectly well and are "
                 "both included either way.")
    else:
        align_all = not heavy_names
    aln_seqs = sequences if align_all else heavy_sequences
    aln_anns = annotations if align_all else heavy_annotations
    if other_names and not align_all and heavy_names:
        st.caption(f"Aligning the {len(heavy_names)} heavy-chain domain(s) only.")

    # Identifies the exact set being aligned, so the drawn state does not carry
    # over to a different submission — otherwise a new analysis would silently
    # pay for the heavy figures again.
    aln_signature = hash(tuple(sorted(aln_seqs.items())))

    if len(aln_seqs) < 2:
        st.info("Submit at least two sequences to see an alignment. "
                + ("Tick the box above to include the non-VHH sequences."
                   if other_names and not align_all else ""))
    else:
        # The alignment and the logo are by some distance the two most expensive
        # figures in the app, and Streamlit re-executes every tab body on every
        # interaction — so without this they are redrawn each time, whether or
        # not anyone is looking at this tab. Behind a button they cost one click
        # per set of sequences, in the one tab people open in order to look at a
        # figure rather than to check a number in passing.
        drawn_for = st.session_state.get("alignment_drawn_for")
        if drawn_for != aln_signature:
            st.info(
                "The alignment and the sequence logo are the two heaviest "
                "figures here, so they are drawn on request rather than on "
                "every interaction.", icon="🖼️")
            if st.button("Draw the alignment and logo", type="primary",
                         key="draw_alignment"):
                st.session_state["alignment_drawn_for"] = aln_signature
                st.rerun()

        if st.session_state.get("alignment_drawn_for") == aln_signature:
            msa_wrap = st.slider(
                "Positions per row", 20, 140, 60, step=5, key="msa_wrap",
                help="How many alignment columns each row carries. Fewer "
                     "positions per row means more rows and wider letters.")
            with st.spinner("Aligning…"):
                # pymsaviz reads from a path, so the aligned FASTA is written to
                # a temporary directory removed as soon as the figure exists.
                img_msa, aln_text = png_msa(aln_seqs, aln_anns, int(msa_wrap))
            st.image(img_msa, **FIT)
            st.caption("FAMSA alignment rendered with pymsaviz. CDR bands refer "
                       "to the first sequence submitted.")
            st.download_button("Download alignment (FASTA)", aln_text,
                               file_name="alignment.fasta", mime="text/plain")

        st.subheader("Pairwise identity")
        cols_id = st.columns(2)
        with cols_id[0]:
            st.image(png_identity(aln_seqs, None, None), **FIT)
        if aln_anns and any((a.get("CDR3") or {}).get("seq")
                            for a in aln_anns.values()):
            img_id3 = png_identity(aln_seqs, aln_anns, "CDR3")
            if img_id3:
                with cols_id[1]:
                    st.image(img_id3, **FIT)
        st.caption("Two clones from the same panning can share 90%+ full-length "
                   "identity while having unrelated CDR3s — the second matrix is the "
                   "one that tells you whether they are really the same binder.")

    if st.session_state.get("alignment_drawn_for") == aln_signature:
        st.divider()
        st.subheader("Sequence logo")
        st.markdown(
            "The alignment as a logo: letter height shows what is conserved and "
            "what varies, which reads faster than the alignment itself when you "
            "only want to know where the diversity sits.")

        lg = st.columns([1, 1, 1])
        with lg[0]:
            logo_type = st.selectbox(
                "Height", ["probability", "information"], index=0,
                format_func=lambda t: "Frequency" if t == "probability"
                                      else "Information (bits)",
                help="Frequency: letter height is the residue's share of the "
                     "column, so every column sums to 1 and reads as "
                     "composition. Information: height in bits, the classic "
                     "Schneider logo, where a mixed column stays short even if "
                     "one residue is commonest.")
        with lg[1]:
            logo_wrap = st.slider("Positions per row", 20, 140, 65, step=5,
                                  help="How many alignment columns each row "
                                       "carries. Fewer rows, wider letters.")
        with lg[2]:
            logo_scheme = st.selectbox(
                "Colours", ["chemistry", "charge", "hydrophobicity", "skylign_protein"],
                index=0)

        with st.spinner("Drawing the logo…"):
            img_logo, logo_matrix = png_logo(aln_seqs, aln_anns, logo_type,
                                             int(logo_wrap), logo_scheme)
        st.image(img_logo, **FIT)
        st.caption("CDR bands refer to the first sequence submitted, mapped onto "
                   "alignment columns. Amino acids that never occur are dropped "
                   "from the matrix, so the colouring only shows residues that "
                   "are actually present.")
        st.download_button("Logo matrix (CSV)", logo_matrix.to_csv(),
                           file_name=f"logo_{logo_type}.csv", mime="text/csv")
        plt.close("all")

# --- vs Reference ----------------------------------------------------------- #
with tabs["vs Reference"]:
    st.markdown(
        "Which positions differ from a sequence you choose as the reference. "
        "Humanisation variants against the parent, a Sanger clone against the "
        "construct it should be, a panel against its ancestor — the same "
        "question each time.")
    if not heavy_only_notice("Reference comparison"):
        pass
    elif len(heavy_sequences) < 2 and not st.session_state.get("ref_paste", "").strip():
        st.info("Submit at least two sequences, or paste a reference below, to "
                "compare against.")
        st.text_area("Reference sequence (optional)", key="ref_paste", height=90,
                     placeholder="Paste the parent sequence here")
    else:
        ref_cols = st.columns([2, 1])
        with ref_cols[0]:
            ref_choice = st.selectbox(
                "Reference", list(heavy_sequences) + ["— paste one —"],
                help="Everything else is compared against this sequence.")
        with ref_cols[1]:
            show_cdr_only = st.checkbox("CDR positions only", value=False)

        if ref_choice == "— paste one —":
            pasted = st.text_area("Reference sequence", key="ref_paste", height=90,
                                  placeholder="Paste the parent sequence here")
            ref_seq, ref_name = pasted.strip(), "pasted reference"
            others = dict(heavy_sequences)
        else:
            ref_seq, ref_name = heavy_sequences[ref_choice], ref_choice
            others = {k: v for k, v in heavy_sequences.items() if k != ref_choice}

        if not ref_seq:
            st.info("Paste a reference sequence to compare against.")
        elif not others:
            st.info("Nothing to compare — submit another sequence.")
        else:
            try:
                comp = comparison_for(ref_seq, others, ref_name)
            except Exception as exc:
                st.error(f"Could not compare: {exc}")
                comp = None

            if comp:
                for name, msg in comp["errors"].items():
                    st.warning(f"**{name}** — could not be numbered: {msg}")

                st.dataframe(comp["summary"], hide_index=True, **FIT)
                st.caption(
                    "Positions are matched by **IMGT number**, not by raw "
                    "alignment. A one-residue deletion in CDR3 is therefore "
                    "reported as one deletion rather than as every downstream "
                    "position having changed. Note that IMGT places insertions "
                    "at the centre of a loop, so adding a residue to CDR3 can "
                    "also shift a few neighbours — those show up as "
                    "substitutions alongside the insertion.")

                diffs = comp["differences"]
                if show_cdr_only and len(diffs):
                    diffs = diffs[diffs["in_CDR"]]

                if not len(diffs):
                    st.success("No difference from the reference"
                               + (" in the CDRs." if show_cdr_only else "."))
                else:
                    st.image(png_differences(ref_seq, others, ref_name), **FIT)
                    st.dataframe(diffs, hide_index=True, **FIT)
                    st.download_button(
                        "Differences (CSV)", diffs.to_csv(index=False),
                        file_name=f"differences_vs_{ref_name}.csv",
                        mime="text/csv")
    plt.close("all")


# --- Repertoire ------------------------------------------------------------- #
with tabs["Repertoire"]:
    st.markdown(
        "Your sequences against the natural VHH repertoire — the same "
        "position-frequency matrix the synthetic generator samples from, used "
        "here the other way round.")
    if not heavy_only_notice("Repertoire comparison"):
        pass
    else:
        if _vh_only:
            st.warning(
                "The reference is a **VHH** repertoire. A conventional VH "
                "compared against it will legitimately show its FR2 hallmark "
                "positions as rare — that is the difference between a VH and a "
                "VHH, not a defect in the sequence.", icon="🧬")

        rar, nat = st.tabs(["Unusual framework residues",
                            "Your set vs the natural repertoire"])

        with rar:
            st.markdown(
                "Framework positions whose residue is rare in the reference. "
                "Worth a second look: this is what catches a Sanger miscall, a "
                "PCR error, or a construct that has drifted without anyone "
                "noticing.")
            rc = st.columns(3)
            with rc[0]:
                rare_thr = st.slider("Rare at or below", 0.001, 0.10,
                                     float(rep.RARE_THRESHOLD), step=0.001,
                                     format="%.3f")
            with rc[1]:
                min_sup = st.slider("Minimum support", 0, 1900,
                                    int(rep.MIN_SUPPORT), step=50,
                                    help="Positions backed by fewer reference "
                                         "sequences than this are set aside "
                                         "rather than called rare.")
            with rc[2]:
                rare_plot = st.checkbox("Plot each sequence",
                                        value=len(heavy_sequences) <= 4)

            summary, flagged = rarity_tables(heavy_sequences, float(rare_thr),
                                             int(min_sup))
            st.dataframe(summary, hide_index=True, **FIT)
            st.caption(
                "**Rare is not wrong.** A residue absent from the reference may "
                "be a genuine novelty, a germline this set under-samples, or the "
                "very mutation you selected for. Every call carries the number "
                "of reference sequences behind it.")

            if len(flagged):
                st.subheader("Flagged positions")
                st.dataframe(flagged, hide_index=True, **FIT)
                st.download_button("Unusual residues (CSV)",
                                   flagged.to_csv(index=False),
                                   file_name="unusual_framework_residues.csv",
                                   mime="text/csv")
            else:
                st.success("No framework residue falls below the threshold — "
                           "every position matches the natural repertoire.")

            if rare_plot:
                for name, seq in heavy_sequences.items():
                    try:
                        st.image(png_rarity(seq, name, float(rare_thr),
                                            int(min_sup)), **FIT)
                    except Exception as exc:
                        st.warning(f"**{name}** — {exc}")

        with nat:
            st.markdown(
                "Where your set sits against the reference distribution. After "
                "a panning this is what says whether the selection was biased — "
                "towards long CDR3s, towards charged loops, and so on.")
            if not heavy_annotations:
                st.info("Nothing annotated to compare.")
            else:
                per_seq, distribution, composition = repertoire_stats(
                    heavy_sequences, heavy_annotations)
                st.image(png_cdr_lengths(heavy_sequences, heavy_annotations), **FIT)
                st.caption(
                    "Grey is the natural repertoire, red your set. With only a "
                    "handful of sequences the red bars are coarse by "
                    "construction — one sequence in ten is a 10% bar.")
                st.dataframe(per_seq, hide_index=True, **FIT)

                st.subheader("Amino-acid composition")
                regions = st.multiselect(
                    "Regions", ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"],
                    default=["CDR1", "CDR2", "CDR3"])
                if regions:
                    st.image(png_composition(heavy_sequences, heavy_annotations,
                                             tuple(regions)), **FIT)
                    st.caption(
                        "Above the line: enriched in your set relative to the "
                        "natural repertoire. Below: depleted. The reference "
                        "composition is the position-frequency matrix averaged "
                        "over each region, weighted by support.")
                    shown = composition[composition["region"].isin(regions)]
                    st.dataframe(
                        shown.sort_values("delta", ascending=False),
                        hide_index=True, **FIT)
                    st.download_button("Composition (CSV)",
                                       composition.to_csv(index=False),
                                       file_name="composition_vs_repertoire.csv",
                                       mime="text/csv")
    plt.close("all")


# --- Clonal families ------------------------------------------------------- #
with tabs["Clonal families"]:
    st.markdown(
        "After a panning, the first question is not how good each clone is but "
        "**how many different clones you actually have** — a plate of picked "
        "colonies is routinely the same few lineages over and over. Grouping on "
        "CDR3 answers that, because CDR3 is what differs between independent "
        "lineages while the frameworks stay nearly constant.")

    if not heavy_only_notice("Clonal family grouping"):
        pass
    elif len(heavy_sequences) < 2:
        st.info("Submit at least two VHH sequences to group them into families.")
    elif not any((a.get("CDR3") or {}).get("seq")
                 for a in heavy_annotations.values()):
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
            fam = families_for(heavy_sequences, heavy_annotations, float(cutoff),
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

        img_dendro = png_family_dendrogram(heavy_sequences, heavy_annotations,
                                           float(cutoff),
                                           linkage_method, bool(same_length))
        if img_dendro:
            st.subheader("Dendrogram")
            st.image(img_dendro, **FIT)
            st.caption("The dashed line is your threshold; branches joining to "
                       "its right are separate families. Move the slider and "
                       "watch which groups merge — that is the honest way to see "
                       "how robust a family is.")

        st.subheader("Per-sequence assignment")
        st.dataframe(fam["assignment"].assign(
            type=[type_of(x) for x in fam["assignment"]["name"]]),
            hide_index=True, **FIT)

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
            reps = {r: heavy_sequences[r]
                    for r in fam["families"]["representative"]
                    if r in heavy_sequences}
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
    dna_out = st.columns(2)
    with dna_out[0]:
        st.button("🧬 Use these in the cloning section", key="to_cloning", **FIT,
                  on_click=go_to, args=("cloning",))
    with dna_out[1]:
        if st.button("🔧 Send to gene diagnosis", key="to_genes", **FIT,
                     help="Run these sequences through the frame, rare-codon and "
                          "GC checks, and get synonymous swaps where GC runs "
                          "high."):
            st.session_state["pending_gene_dna"] = "\n".join(fasta_parts) + "\n"
            go_to("genes")
            st.rerun()

# --- Constructs ------------------------------------------------------------- #
with tabs["Constructs"]:
    st.markdown(
        "Assemble a VHH into the construct you will actually order — monovalent, "
        "bivalent or biparatopic, with a linker and tags — and get the mass, pI "
        "and DNA of the finished protein rather than of the bare domain.")
    st.warning(ctm.VERIFY_NOTICE, icon="⚠️")

    fmt_cols = st.columns([1, 1, 1])
    with fmt_cols[0]:
        fmt = st.selectbox("Format", list(ctm.FORMATS),
                           help="\n\n".join(f"**{k}** — {v['note']}"
                                             for k, v in ctm.FORMATS.items()))
    n_needed = ctm.FORMATS[fmt]["n_vhh"]
    with fmt_cols[1]:
        first = st.selectbox("VHH", list(sequences), key="ctor_first")
    second = None
    if n_needed == 2:
        with fmt_cols[2]:
            if fmt == "bivalent":
                st.caption("Bivalent uses the same VHH twice.")
                second = first
            else:
                second = st.selectbox(
                    "Second VHH", list(sequences),
                    index=min(1, len(sequences) - 1), key="ctor_second")
                if second == first:
                    st.caption("Two copies of the same VHH is a bivalent "
                               "construct, not a biparatopic one.")

    lk_cols = st.columns([1, 1])
    with lk_cols[0]:
        linker_name = st.selectbox(
            "Linker between VHHs", list(ctm.LINKERS),
            index=list(ctm.LINKERS).index("(GGGGS)3"),
            disabled=n_needed == 1)
        st.caption(ctm.LINKERS[linker_name]["note"])
    with lk_cols[1]:
        tag_linker = st.text_input(
            "Linker before each tag", value="GS",
            help="Leave empty to fuse tags directly. Two or three residues is "
                 "usual and often just the cloning scar.")

    tag_cols = st.columns(2)
    with tag_cols[0]:
        n_tag_names = st.multiselect(
            "N-terminal tags",
            [k for k, v in ctm.TAGS.items() if v["ends"] in ("N", "both")],
            help="Applied in the order you pick them, reading N to C.")
    with tag_cols[1]:
        c_tag_names = st.multiselect(
            "C-terminal tags",
            [k for k, v in ctm.TAGS.items() if v["ends"] in ("C", "both")],
            default=["His6"])

    with st.expander("Custom pieces — SpyTag2/3, a camelid hinge, anything else"):
        st.markdown(
            "Two things are deliberately missing from the tag list above: "
            "**SpyTag2 and SpyTag3**, and the **camelid hinges**. Only the "
            "original SpyTag is included. Writing the others from memory would "
            "risk handing you a sequence that is one residue wrong, which you "
            "would discover only after the gene arrived — so paste the exact "
            "sequence from your own construct here instead.")
        cu = st.columns(2)
        custom_n = cu[0].text_input("Custom N-terminal peptide", value="")
        custom_n_label = cu[0].text_input("Its label", value="custom-N")
        custom_c = cu[1].text_input("Custom C-terminal peptide", value="")
        custom_c_label = cu[1].text_input("Its label", value="custom-C")

    chosen = [(first, sequences[first])]
    if n_needed == 2 and second:
        chosen.append((second, sequences[second]))

    n_tags = [(k, ctm.TAGS[k]["seq"]) for k in n_tag_names]
    c_tags = [(k, ctm.TAGS[k]["seq"]) for k in c_tag_names]
    if custom_n.strip():
        n_tags.append((custom_n_label or "custom-N", custom_n))
    if custom_c.strip():
        c_tags.append((custom_c_label or "custom-C", custom_c))

    try:
        built = ctm.build_construct(
            chosen, linker=ctm.LINKERS[linker_name]["seq"] if n_needed == 2 else "",
            n_tags=n_tags, c_tags=c_tags, tag_linker=tag_linker, ph=ph)
    except Exception as exc:
        st.error(f"Could not assemble the construct: {exc}")
        built = None

    if built:
        for w in built["warnings"][1:]:        # [0] is the verify notice, shown above
            st.warning(w, icon="⚠️")

        st.image(png_construct(built, f"{first} — {fmt}"), **FIT)

        props = built["properties"].iloc[0]
        base = built["per_vhh"].iloc[0]
        pm = st.columns(4)
        pm[0].metric("Length", f"{built['length']} aa",
                     delta=f"{built['length'] - int(base['length']):+d} vs the VHH")
        pm[1].metric("Mass", f"{props['mw_kda']:.2f} kDa",
                     delta=f"{props['mw_kda'] - base['mw_kda']:+.2f}")
        pm[2].metric("pI", f"{props['pI']:.2f}",
                     delta=f"{props['pI'] - base['pI']:+.2f}")
        pm[3].metric("A280 (1 mg/mL)", f"{props['A280_0.1pct']:.3f}",
                     help="Divide a measured A280 by this to get mg/mL — the "
                          "number changes with the tags, which is exactly why "
                          "you want it for the construct and not the bare VHH.")

        st.markdown("**Segments**")
        st.dataframe(built["segments"], hide_index=True, **FIT)

        st.markdown("**Full properties**")
        st.dataframe(built["properties"], hide_index=True, **FIT)

        st.markdown("**Protein sequence**")
        st.code(built["sequence"], language="text")

        with st.expander("Codon-optimised DNA for the whole construct"):
            avoid_ct = st.multiselect(
                "Restriction sites to avoid", sorted(vc.DEFAULT_AVOID_SITES),
                default=sorted(vc.DEFAULT_AVOID_SITES), key="ctor_avoid")
            stop_ct = st.checkbox("Append a stop codon (TAA)", value=True,
                                  key="ctor_stop")
            dna_ct = ctm.construct_dna(built, avoid=set(avoid_ct),
                                       add_stop=bool(stop_ct))
            if dna_ct["back_translation_ok"]:
                st.success("Round trip verified: translating this DNA gives back "
                           "the construct exactly.")
            else:
                st.error("The DNA does not translate back to the construct — do "
                         "not order this. Please report it.")
            st.code(dna_ct["dna"], language="text")
            st.download_button(
                "Construct DNA (FASTA)",
                f">{first}_{fmt}\n{dna_ct['dna']}\n",
                file_name=f"{first}_{fmt}.fasta", mime="text/plain")

        dlc = st.columns(2)
        dlc[0].download_button(
            "Construct protein (FASTA)",
            f">{first}_{fmt}\n{built['sequence']}\n",
            file_name=f"{first}_{fmt}_protein.fasta", mime="text/plain")
        dlc[1].download_button(
            "Segments (CSV)", built["segments"].to_csv(index=False),
            file_name=f"{first}_{fmt}_segments.csv", mime="text/csv")
    plt.close("all")


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

# The page is now fully built, so the status line under the Analyze button can
# be replaced with the figure the user actually experienced. A placeholder can be
# rewritten later in the same run, which is exactly what this relies on.
_slot = st.session_state.pop("_status_slot", None)
_parts = st.session_state.pop("_status_parts", None)
if _slot is not None and _parts is not None:
    _total, _analysis = _parts
    _whole = time.perf_counter() - _RUN_STARTED
    _fmt = lambda t: f"{t * 1000:.0f} ms" if t < 1 else f"{t:.1f} s"
    _slot.success(
        f"{_total} sequence{'s' if _total != 1 else ''} · "
        f"characterised in {_fmt(_analysis)} · page built in {_fmt(_whole)}",
        icon="✅")

footer()
