"""
shell - the app's chrome: identity, navigation, caching policy, footer.

Everything here is shared by every section and belongs to none of them. It lives
in its own module so that `app.py` can be a router and a section file can import
what it needs without importing the app.

Nothing in here computes anything about a sequence. That is `compute.py`.
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

__version__ = "2026.10.28c"

EXAMPLE = """>caplacizumab
EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEG
RFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS
>cAbLys3
DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKG
RFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS
>Nb80
QVQLQESGGGLVQAGGSLRLSCAASGSIFSINTMGWYRQAPGKQRELVAAIHSGGSTNYANSVKGR
FTISRDNAANTVYLQMNSLKPEDTAVYYCNVKDYGAVLYEYDYWGQGTQVTVSS
"""


# The release this build belongs to, as semantic versioning: MAJOR.MINOR.PATCH.
# It must match the git tag and the Zenodo record, so that a user reporting a bug
# or citing a version names something those two know about. Bump it when you cut
# a release, not on every change.
#
#   PATCH  a fix, nothing new                     (1.0.0 -> 1.0.1)
#   MINOR  a new feature, nothing broken          (1.0.1 -> 1.1.0)
#   MAJOR  a change that breaks an existing use   (1.1.0 -> 2.0.0)
#
# MAJOR is not hypothetical here: the analysis modules are importable from a
# notebook, so changing the signature of a vhh_core function breaks somebody's
# notebook even though the app looks identical.
RELEASE = "1.0.0"

# The day this particular build was assembled. Two builds of the same release
# are possible - a fix pushed between releases - and this is what tells them
# apart in a bug report.
BUILD = "20261028"

VERSION = f"{RELEASE} (build {BUILD})"


# The author as the work should be attributed, not as a platform handle: the
# citation metadata, the licence and the Zenodo record all carry the real name,
# and a footer that disagrees with them is a small way of making the tool harder
# to credit.
AUTHOR = "Jonathan Vaneyck"


AUTHOR_URL = "https://orcid.org/0000-0003-3095-7420"


REPO_URL = "https://github.com/TyrosinDuck/vhh-analyzer"


# The Zenodo *concept* DOI: it always resolves to the latest archived release,
# which is the one that belongs in a reference list. Each release also has a
# version DOI of its own, for citing the exact code that produced a result.
DOI = "10.5281/zenodo.23742269"


DOI_URL = f"https://doi.org/{DOI}"


LICENCE = "MIT"


APP_NAME = "VHH tool box"


APP_ICON = "🧰"


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


CACHE_TTL = 600        # seconds


CACHE_MAX_ENTRIES = 8


_FIG_CACHE = dict(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=24)


def fig_to_png(fig, dpi: int = 200, close: bool = False) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return buf.getvalue()


@st.cache_resource
def anarci_status() -> bool:
    """Probe ANARCI once per server process, not once per request."""
    return vc.anarci_available()


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
    st.caption(
        f"{APP_NAME} · v{VERSION} · by [{AUTHOR}]({AUTHOR_URL}) · "
        f"{LICENCE} licence · [source]({REPO_URL}) · [doi:{DOI}]({DOI_URL})")


SECTIONS = {
    "analysis": ("🔬", "Sequence analysis",
                 "Characterize VHHs you provide: physico-chemical properties, "
                 "CDR/framework annotation, hallmarks and germline, liabilities, "
                 "alignment, codon-optimised DNA and a PDF report."),
    "structure": ("🧊", "Structure",
                  "Upload a PDB, mmCIF or AlphaFold model and replace the "
                  "sequence-based guesses with measurements: which cysteines "
                  "actually pair, which liabilities are buried and which are "
                  "exposed, where the hydrophobic surface is, and — for a "
                  "complex — which residues contact the antigen."),
    "sequencing": ("🧪", "Sanger sequencing",
                   "Upload .ab1 chromatograms from a screening run, extract the VHH "
                   "sequences from the expression cassette, spot frameshifts and "
                   "duplicates — then send the clones you want to the analysis."),
    "genes": ("🔧", "Gene diagnosis",
              "Paste coding sequences that are already in a plasmid and ask what "
              "might be wrong with them: frame and stop checks, documented rare "
              "codons and their clusters, GC stretches — plus synonymous swaps "
              "that lower GC without touching the protein."),
    "library": ("📚", "Sequence library",
                "The sequences you come back to, kept in the repository and one "
                "click away: send the proteins to the analysis section or the "
                "DNA to gene diagnosis, instead of pasting them again."),
    "cloning": ("🧬", "In-silico cloning",
                "Load a backbone and an insert, and get the restriction enzymes that "
                "open the vector in the right place without cutting the VHH, plus "
                "which of your primers actually bind."),
}
