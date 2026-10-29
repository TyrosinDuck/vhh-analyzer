"""
VHH tool box, Streamlit front-end.

The app is organised in three sections, reachable from the home page and from
the sidebar:

* **Sequence analysis**: characterization, CDR annotation, alignment, DNA
* **Sanger sequencing**: extract VHHs from .ab1 chromatograms, then send the
  ones you pick to the analysis section
* **In-silico cloning**: pick a backbone, an insert and get the enzymes and
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
import importlib
import tempfile
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import cloning as cl
import generator as gn
import constructs as ctm
import repertoire as rep
import triage as tg
import codons as cdn
import library as lib
import sanger as sg
import structure as stx
import vhh_core as vc
from sections import (analysis as sections_analysis,
                      cloning as sections_cloning, genes as sections_genes,
                      home as sections_home, library as sections_library,
                      manual as sections_manual,
                      sequencing as sections_sequencing,
                      structure as sections_structure)
import compute as cpt
import shell as shl
# The shared chrome and the cached computation layer. Listed by name
# rather than star-imported, so that this block is the documented
# surface the section bodies below are allowed to rely on.
from shell import (APP_ICON, APP_NAME, AUTHOR, AUTHOR_URL, CACHE_MAX_ENTRIES,
    CACHE_TTL, EXAMPLE, FIT, SECTIONS, VERSION, _FIG_CACHE, _fit,
    anarci_status, current_sequences, fig_to_png, footer, go_to, load_manual)
from compute import (_analyze, _analyze_one, aggregation_for, analyze_one,
    run_analysis,
    analyze_with_progress, classify_cached, comparison_for, diagnose_gene,
    disulfides_for, dna_for, evaluate_pairs_cached, families_for, gc_fix_for,
    gene_summary, generate_cached, generator_matrix, germline_data,
    hallmark_data, humanness_for, liabilities_for, load_library_cached,
    load_plasmid_cached, pcr_cached, plan_cloning_cached, png_aggregation,
    png_cdr_lengths, png_charge, png_chou_fasman, png_composition,
    png_construct, png_differences, png_family_dendrogram, png_gc_comparison,
    png_gene_diagnosis, png_identity, png_liability_map,
    png_liability_profile, png_logo, png_msa, png_plasmid_map, png_profiles,
    png_quality, png_rarity, png_secondary_structure, primer_report_cached,
    process_ab1, proteases_for, rarity_for, rarity_tables, repairs_for,
    repertoire_stats, st_confidence, st_disulfides, st_exposure,
    st_identity_map, st_interface, st_interface_summary, st_load, st_patches,
    st_resnum_map, triage_for)


# Three published sequences rather than invented ones, each verifiable: see
# data/library.csv for their PDB entries. The first of these used to be labelled
# "NbExample1" and presented as a generic VHH, it is in fact caplacizumab, which
# is a poor thing to leave mislabelled in a public app.

# Version stamp, YYYYMMDD.HHMM. Bump this by hand when you push a change that
# users would notice, it is what someone will quote back to you in a bug report.

# Wall-clock start of this script run. Streamlit re-executes the whole file on
# every interaction, so this is set afresh each time and measures what the user
# actually waits for, not just the analysis loop.
_RUN_STARTED = time.perf_counter()

# The app name and its icon live here so the page title, the sidebar, the home
# page, the footer and the PDF report can never drift apart.

st.set_page_config(page_title=APP_NAME, page_icon=APP_ICON, layout="wide")

# The analysis modules live in separate files that are updated together with this
# one. Copying only app.py across used to fail deep inside a call with an opaque
# AttributeError; this turns that into a message naming the stale file.
REQUIRED_MODULES = {"vhh_core": "2026.10.28b", "sanger": "2026.10.21b",
                    "cloning": "2026.10.21b", "generator": "2026.10.21b",
                    "constructs": "2026.10.21d",
                    "repertoire": "2026.10.27a",
                    "triage": "2026.10.26a", "codons": "2026.10.27e",
                    "library": "2026.10.28a", "structure": "2026.10.28e",
                    "shell": "2026.10.29a", "compute": "2026.10.29b",
                    "sections.analysis": "2026.10.28a",
                    "sections.cloning": "2026.10.28c",
                    "sections.genes": "2026.10.28c",
                    "sections.home": "2026.10.29b",
                    "sections.library": "2026.10.28c",
                    "sections.manual": "2026.10.29b",
                    "sections.sequencing": "2026.10.28c",
                    "sections.structure": "2026.10.28c"}
# Resolved from the dict by name rather than from a second hand-written list of
# (name, wanted, module) triples. That list was the duplication the dict was
# supposed to avoid: adding a module meant editing both, and forgetting the
# second one meant the new module was never checked at all.
_stale = [f"`{name}.py` (has {getattr(_mod, '__version__', 'no version')}, "
          f"needs {want})"
          for name, want in REQUIRED_MODULES.items()
          for _mod in (importlib.import_module(name),)
          if getattr(_mod, "__version__", None) != want]
if _stale:
    st.error("Some modules are out of date, update " + ", ".join(_stale) +
             " alongside `app.py`, they are released together.")
    st.stop()





# How long a result may stay in the in-memory cache, and how many results are
# kept at once. The cache makes moving a sidebar slider instant instead of
# re-running the whole analysis, but it is shared by the whole server process,
# so it is kept deliberately small and short-lived. Set CACHE_TTL = 0 to switch
# caching off entirely: the app then recomputes on every interaction, which
# costs well under a second for a handful of sequences.


# --------------------------------------------------------------------------- #
# Cached workers
# --------------------------------------------------------------------------- #





















# Streamlit re-runs the whole script (and every tab body) on each interaction,
# so without caching, moving one slider re-renders every figure. These return
# PNG bytes rather than figure objects because bytes cache cleanly and a static
# image is all these plots need; the PDF builder rebuilds real figures on demand.






















































































# --------------------------------------------------------------------------- #
# Navigation
# --------------------------------------------------------------------------- #


st.session_state.setdefault("section", "home")










# A widget's state cannot be written after the widget exists, so anything that
# wants to fill the text area (the sequencing section) parks its text here and
# reruns; it is moved into place before the widget is created.
if "pending_seq_text" in st.session_state:
    st.session_state["seq_text"] = st.session_state.pop("pending_seq_text")
if "pending_gene_dna" in st.session_state:
    st.session_state["gene_dna"] = st.session_state.pop("pending_gene_dna")
# Streamlit discards a widget's state once the widget stops being rendered, and
# leaving a section is exactly that. Re-asserting these keys on every run, before
# any section draws, is what keeps an input filled when you come back to it.
for _input_key in ("seq_text", "gene_dna"):
    st.session_state.setdefault(_input_key, "")
    st.session_state[_input_key] = st.session_state[_input_key]

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
                   f"{'available' if has_anarci else 'not installed, regex only'}")
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
Traffic is encrypted, nothing is stored, but the data does pass through a server
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
    sections_home.render()
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# Manual
# --------------------------------------------------------------------------- #

if section == "manual":
    sections_manual.render()
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# Sanger sequencing
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Structure
# --------------------------------------------------------------------------- #
# Everything here is cached on the raw file bytes rather than on a parsed model,
# because a DataFrame pickles cleanly and a Bio.PDB entity is awkward to hash.
# Re-parsing inside each cached function costs a few milliseconds; the SASA
# computation it protects costs seconds.



















if section == "structure":
    sections_structure.render(method)
    footer()
    st.stop()


if section == "sequencing":
    sections_sequencing.render()
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# In-silico cloning
# --------------------------------------------------------------------------- #







# --------------------------------------------------------------------------- #
# Gene diagnosis
# --------------------------------------------------------------------------- #

if section == "genes":
    sections_genes.render()
    footer()
    st.stop()


# --------------------------------------------------------------------------- #
# Sequence library
# --------------------------------------------------------------------------- #

if section == "library":
    sections_library.render()
    footer()
    st.stop()


if section == "cloning":
    sections_cloning.render()
    footer()
    st.stop()
# Everything that is not one of the sections above is the analysis screen.
sections_analysis.render(ph, window, method)

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
