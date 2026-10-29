"""
compute - every cached computation the app performs, and nothing else.

These are thin wrappers: each one calls into `vhh_core`, `codons`, `structure`
and the rest, and exists only to attach `st.cache_data` to the call. Keeping them
together, and apart from the panels that draw their results, is what makes the
caching policy auditable - you can read the whole of it in one file instead of
finding 56 decorators spread through three thousand lines of UI.

Two rules hold throughout:

- **A cached function takes everything it depends on as an argument.** A cached
  result keyed on less than it depends on is a wrong answer waiting for the
  second call.
- **Nothing here draws.** A function that both computes and writes to the page
  cannot be cached, because the drawing would be skipped on a cache hit.
"""

from __future__ import annotations

import io
import re
import tempfile
import time
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import cloning as cl
import codons as cdn
import constructs as ctm
import generator as gn
import library as lib
import repertoire as rep
import sanger as sg
import structure as stx
import triage as tg
import vhh_core as vc

from shell import CACHE_MAX_ENTRIES, CACHE_TTL, _FIG_CACHE, fig_to_png

__version__ = "2026.10.29b"

def _analyze(seqs: dict[str, str], ph: float, method: str):
    return vc.analyze(seqs, ph=ph, method=method)


def _analyze_one(name: str, seq: str, ph: float, method: str):
    """Analyse a single sequence, the unit the progress bar advances on.

    Caching per sequence rather than per batch means adding one sequence to a
    list of ten only costs one analysis, not eleven.
    """
    frame, anns, errs = vc.analyze({name: seq}, ph=ph, method=method)
    return frame.iloc[0].to_dict(), anns.get(name), errs.get(name)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=128)
def classify_cached(seq: str, method: str):
    """What kind of sequence this is, ANARCI numbering, so worth caching."""
    return vc.classify(seq, method=method)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def hallmark_data(seq: str, annotation: dict | None = None):
    """Hallmark tetrad for one sequence (cached: ANARCI numbering is not free)."""
    return vc.hallmarks(seq, annotation)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def germline_data(seq: str):
    """Closest germline per species; two ANARCI passes, worth caching."""
    return vc.germline(seq)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=4)
def process_ab1(payload: tuple, vector: str, cutoff: float, min_len: int,
                min_phred: float):
    """Parse and extract a batch of .ab1 uploads (cached on the file contents)."""
    return sg.process_reads(list(payload), vector=vector, mott_cutoff=cutoff,
                            min_trimmed_len=min_len, min_mean_phred=min_phred)


# Caching is a policy, not a property of the function, so it is attached here
# rather than with a decorator: setting CACHE_TTL = 0 must give the uncached
# functions themselves, not a cache with a zero lifetime.
if CACHE_TTL:
    run_analysis = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)(_analyze)
    analyze_one = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL,
        max_entries=CACHE_MAX_ENTRIES * 8)(_analyze_one)
else:
    run_analysis = _analyze
    analyze_one = _analyze_one


def analyze_with_progress(seqs: dict[str, str], ph: float, method: str,
                          slot=None):
    """Run the analysis sequence by sequence, showing progress as it goes.

    ``slot`` is a placeholder created earlier in the script, normally the one
    just under the Analyze button. Without it the bar would be drawn wherever
    this happens to be called, which is below the generator panel and nowhere
    near the control the user just clicked.

    The bar is replaced by a one-line result rather than cleared: on a handful
    of sequences the analysis takes well under a second and a bar that vanishes
    leaves no evidence anything happened at all, and once the cache is warm it
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
        # only. Every tab body still has to render afterwards, the alignment,
        # the logo, the identity matrices, the plasmid maps, and that is where
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
            + ", drawing the figures…", icon="⏳")
        st.session_state["_status_parts"] = (total, elapsed)
    return pd.DataFrame(rows), anns, errs


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
    """Clonal families, cached, because the pairwise alignments dominate."""
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
    per_seq, distribution, _c = repertoire_stats(seqs, anns)
    return fig_to_png(rep.plot_cdr_lengths(distribution, len(seqs), per_seq),
                      close=True)


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
def gc_fix_for(dna, gc_window, strategy, avoid, min_usage):
    sites = {k: v for k, v in vc.DEFAULT_AVOID_SITES.items() if k in avoid}
    return cdn.suggest_gc_fix(dna, avoid_sites=sites, gc_window=gc_window,
                              strategy=strategy, min_usage=min_usage)


@st.cache_data(**_FIG_CACHE)
def png_gc_comparison(dna, gc_window, strategy, avoid, min_usage):
    fix = gc_fix_for(dna, gc_window, strategy, avoid, min_usage)
    if not fix["n_changes"]:
        return None
    return fig_to_png(cdn.plot_gc_comparison(dna, fix["dna"],
                                             fix["suggestions"], gc_window),
                      close=True)


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


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_load(raw: bytes, name: str) -> tuple[pd.DataFrame, str, str, int, list]:
    model = stx.load_structure(raw, name)
    kind, reason = stx.confidence_kind(model)
    return (stx.chain_table(model), kind, reason, model.n_models, model.ssbond)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_disulfides(raw: bytes, name: str) -> tuple[pd.DataFrame, list, int, list]:
    model = stx.load_structure(raw, name)
    found = stx.disulfides_3d(model)
    return (found, found.attrs["free_cysteines"], found.attrs["n_cys_with_sg"],
            stx.compare_disulfides(found, model.ssbond))


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_exposure(raw: bytes, name: str, chain: str) -> pd.DataFrame:
    return stx.residue_sasa(stx.load_structure(raw, name), chain, in_context=True)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_confidence(raw: bytes, name: str, chain: str) -> pd.DataFrame:
    return stx.confidence_table(stx.load_structure(raw, name), chain)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_patches(raw: bytes, name: str, chain: str, min_rel: float,
               contact: float, min_size: int) -> pd.DataFrame:
    return stx.surface_patches(stx.load_structure(raw, name), chain,
                               min_rel_sasa=min_rel, contact=contact,
                               min_size=min_size)


@st.cache_data(show_spinner=False, ttl=3600, max_entries=2)
def manual_sections(text: str) -> list[dict]:
    """Split MANUAL.md on its top-level headings.

    The manual is around 1 700 lines. Rendering all of it on every interaction
    is both slow and useless: nobody reads a reference manual from the top. This
    cuts it into its numbered sections so the screen can show one at a time,
    and collects each section's own subheadings for a local table of contents.
    """
    import re as _re

    # Everything before the first "## " heading is the document title and the
    # introduction. Dropping it is easy to do and invisible afterwards, because
    # the screen only ever shows one section and nothing looks missing, so it
    # is kept as a section of its own.
    sections: list[dict] = []
    preamble = {"title": "Introduction", "lines": [], "subs": []}
    current = preamble
    for line in text.split("\n"):
        top = _re.match(r"^##\s+(.+?)\s*$", line)
        sub = _re.match(r"^###\s+(.+?)\s*$", line)
        if top and not line.startswith("###"):
            current = {"title": top.group(1), "lines": [], "subs": []}
            sections.append(current)
            continue
        if sub:
            current["subs"].append(sub.group(1))
        current["lines"].append(line)
    if "".join(preamble["lines"]).strip():
        title = next((l[2:].strip() for l in preamble["lines"]
                      if l.startswith("# ")), "Introduction")
        preamble["title"] = title
        sections.insert(0, preamble)
    for section in sections:
        section["body"] = "\n".join(section["lines"]).strip("\n")
        del section["lines"]
    return sections


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_resnum_map(raw: bytes, name: str, chain: str,
                  target_seq: str) -> tuple[dict, float]:
    model = stx.load_structure(raw, name)
    structure_seq, _residues = stx.chain_sequence(model.model[chain])
    return stx.resnum_map(model, chain, structure_seq, target_seq)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_identity_map(raw: bytes, name: str, chain: str) -> dict:
    return stx.identity_resnum_map(stx.load_structure(raw, name), chain)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_interface(raw: bytes, name: str, chain: str,
                 partners: tuple) -> pd.DataFrame:
    return stx.interface_table(stx.load_structure(raw, name), chain,
                               list(partners))


@st.cache_data(show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)
def st_interface_summary(raw: bytes, name: str, a: str, b: str) -> dict:
    return stx.interface_summary(stx.load_structure(raw, name), a, b)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=2)
def load_library_cached():
    return lib.load_library()


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=16)
def gene_summary(seqs, gc_window, cluster_window):
    return cdn.summarise(seqs, gc_window=gc_window, cluster_window=cluster_window)
