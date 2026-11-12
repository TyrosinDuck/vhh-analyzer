"""
VHH tool box, Streamlit front-end: the router, and nothing else.

What this file does, in order: import every module of the project inside one
guard, refuse to run if any of them is older than this release expects, draw
the sidebar, and hand the run to exactly one screen.

What it does not do: compute, draw a screen, or hold a list of the screens.
Every cached computation is in `compute.py`, every screen is a `render` in
`sections/`, and the set of screens is derived from `shell.screen_keys()`,
which reads `shell.SECTIONS` and `shell.UTILITIES`. Adding a section is one
entry in `SECTIONS`, one file in `sections/`, one line in `REQUIRED_MODULES`,
and no edit here. This file was a single 3 911-line script before the screens
were pulled out of it; keeping it a router is what stops it growing back.

Data handling: submitted sequences are never written to a persistent file,
logged or sent to a third party. They are processed in the server's memory. The
result cache is deliberately short-lived (see CACHE_TTL) and the alignment step
writes a temporary FASTA that is deleted as soon as the figure is rendered.

All analysis logic lives in vhh_core.py, sanger.py, cloning.py, structure.py,
bench.py and buffers.py, which a notebook can import exactly the same way.

Run locally:   streamlit run app.py
"""

from __future__ import annotations

import importlib
import time

import streamlit as st


def _updating(exc: BaseException) -> None:
    """Say the app is updating and end this run, instead of a traceback.

    A deploy swaps the source files under a server that is already running,
    and Streamlit drops the changed modules out of ``sys.modules`` so the next
    run picks them up. A page loaded in that window can ask for a module in
    the middle of that swap and get a ``KeyError`` naming it (seen in the
    console as ``KeyError: 'generator'``), or an ``ImportError`` from a file
    caught half written. Either reaches the visitor as a raw traceback on what
    is otherwise a working app, and it is transient: the next run, a second
    later, imports fine.

    The same message would hide a file that is really missing, which is a
    plausible outcome of copying a build over a folder by hand, so the
    exception itself is printed underneath, small: a message that is still
    there after a reload is not an update.
    """
    st.info("The app is updating. Reload the page in a few seconds.",
            icon="🔄")
    st.caption(f"If this is still here after a reload, the import failed "
               f"with {type(exc).__name__}: {exc}")
    st.stop()


# Every module of this project is imported inside this block, because these
# imports are the first thing a run does and therefore the first thing a
# mid-deploy run trips on. The guard used to sit only around the version check
# further down, which runs after them and so never got the chance.
# tests/test_regressions.py fails if a project module is imported outside it.
#
# The modules are imported for their side effect of being loaded, not for what
# this file calls: a router calls one `render` and the version check. The
# screens and the cached layer import what they need themselves. This block
# used to pull in sixty-two names from `compute` and twenty from `shell`,
# left behind when the screens moved into `sections/`, and a dead import list
# is one more thing to keep in step with nothing.
try:
    import af3
    import bench
    import buffers
    import cloning
    import codons
    import compute
    import constructs
    import generator
    import library
    import reference
    import repertoire
    import sanger
    import shell
    import structure
    import triage
    import vhh_core as vc
    # Submodules, not `from sections import analysis, bench, ...`: half of
    # those names are also the name of a root module imported just above
    # (`bench`, `buffers`, `cloning`, `library`, `reference`, `structure`),
    # and the second binding would silently win. The router reaches a screen
    # through `SCREENS` rather than through any of these names.
    import sections.analysis
    import sections.bench
    import sections.buffers
    import sections.cloning
    import sections.genes
    import sections.home
    import sections.library
    import sections.manual
    import sections.phage
    import sections.reference
    import sections.sequencing
    import sections.structure
    from shell import (APP_ICON, APP_NAME, CACHE_MAX_ENTRIES, CACHE_TTL, FIT,
                       HOME, SECTIONS, Settings, UTILITIES, anarci_status,
                       footer, go_to, group_header_html, screen_keys,
                       sections_by_group)
except (KeyError, ImportError) as _exc:
    _updating(_exc)


# Wall-clock start of this script run. Streamlit re-executes the whole file on
# every interaction, so this is set afresh each time and measures what the user
# actually waits for, not just the analysis loop.
_RUN_STARTED = time.perf_counter()

st.set_page_config(page_title=APP_NAME, page_icon=APP_ICON, layout="wide")

# The analysis modules live in separate files that are updated together with this
# one. Copying only app.py across used to fail deep inside a call with an opaque
# AttributeError; this turns that into a message naming the stale file.
REQUIRED_MODULES = {"vhh_core": "2026.11.10b", "sanger": "2026.10.21b",
                    "cloning": "2026.10.21b", "generator": "2026.10.21b",
                    "constructs": "2026.10.21d",
                    "repertoire": "2026.10.27a",
                    "triage": "2026.10.26a", "codons": "2026.10.27e",
                    "library": "2026.10.28a", "structure": "2026.11.04a",
                    "shell": "2026.11.12c", "compute": "2026.11.12b",
                    "af3": "2026.11.04a",
                    "reference": "2026.10.30a",
                    "sections.analysis": "2026.11.12a",
                    "sections.bench": "2026.11.12b",
                    "sections.phage": "2026.11.12d",
                    "sections.buffers": "2026.11.12c",
                    "buffers": "2026.11.12b",
                    "sections.reference": "2026.11.12a",
                    "bench": "2026.11.12b",
                    "sections.cloning": "2026.11.12a",
                    "sections.genes": "2026.11.12a",
                    "sections.home": "2026.11.12c",
                    "sections.library": "2026.11.12a",
                    "sections.manual": "2026.11.12a",
                    "sections.sequencing": "2026.11.12a",
                    "sections.structure": "2026.11.12a"}


def _imported(name: str):
    """Import a module by name, or end the run while it is being replaced.

    The imports above have already loaded every module, so this is normally a
    lookup in ``sys.modules``. It can still race a deploy that drops a module
    between the two, and gets the same treatment, see ``_updating``.
    """
    try:
        return importlib.import_module(name)
    except (KeyError, ImportError) as exc:
        _updating(exc)


_checked = {name: _imported(name) for name in REQUIRED_MODULES}

# Resolved from the dict by name rather than from a second hand-written list of
# (name, wanted, module) triples. That list was the duplication the dict was
# supposed to avoid: adding a module meant editing both, and forgetting the
# second one meant the new module was never checked at all.
_stale = [f"`{name}.py` (has {getattr(_module, '__version__', 'no version')}, "
          f"needs {want})"
          for name, want in REQUIRED_MODULES.items()
          for _module in (_checked[name],)
          if getattr(_module, "__version__", None) != want]
if _stale:
    st.error("Some modules are out of date, update " + ", ".join(_stale) +
             " alongside `app.py`, they are released together.")
    st.stop()

# One module per screen, taken from the dict that was just version-checked
# rather than from a mapping written out here: a screen this router could draw
# without having checked its version is exactly the stale-file case the check
# exists for. `test_metadata.py` fails if a screen is missing from
# `REQUIRED_MODULES`, which is what would make this raise on a click.
SCREENS = {key: _checked[f"sections.{key}"] for key in screen_keys()}


# --------------------------------------------------------------------------- #
# Session state
# --------------------------------------------------------------------------- #

st.session_state.setdefault("section", HOME)

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

# A key that is not a screen can only come from a session older than a rename.
# It lands on the home page, which is reachable, rather than on a KeyError.
section = st.session_state["section"]
if section not in SCREENS:
    section = HOME


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #

with st.sidebar:
    st.markdown(f"### {APP_ICON} {APP_NAME}")
    st.button("🏠 Home", key="nav_home", **FIT,
              disabled=section == HOME, on_click=go_to, args=(HOME,))
    # Grouped as on the home page, under the same coloured titles, drawn by
    # the same function so the two cannot drift apart. The group title is
    # larger than the section buttons beneath it, which is the hierarchy.
    for group, keys in sections_by_group().items():
        st.markdown(group_header_html(group, "1.1rem"), unsafe_allow_html=True)
        for key in keys:
            icon, label, _desc, _group = SECTIONS[key]
            st.button(f"{icon} {label}", key=f"nav_{key}", **FIT,
                      disabled=section == key, on_click=go_to, args=(key,))
    st.divider()
    # The manual and the reference tables sit below the divider rather than in
    # SECTIONS: they are what you consult, not what you run, they belong to no
    # group, and the home page cards are built from SECTIONS. They are listed
    # in `shell.UTILITIES` all the same, so that the router reaches them the
    # same way it reaches a section.
    for key, (icon, label) in UTILITIES.items():
        st.button(f"{icon} {label}", key=f"nav_{key}", **FIT,
                  disabled=section == key, on_click=go_to, args=(key,))

    st.divider()

    # The settings are only drawn where they have a visible effect, and a
    # screen that does not draw them gets the defaults rather than whatever
    # the analysis screen was last left on: see `Settings.defaults`.
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
        settings = Settings(ph=ph, window=window, method=method)
    else:
        settings = Settings.defaults()

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
# The screen
# --------------------------------------------------------------------------- #

SCREENS[section].render(settings)

# The page is now fully built, so the status line under the Analyze button can
# be replaced with the figure the user actually experienced. A placeholder can be
# rewritten later in the same run, which is exactly what this relies on. Only
# the analysis screen sets these; on every other screen both pops return None.
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
