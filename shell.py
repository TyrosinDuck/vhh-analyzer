"""
shell - the app's chrome: identity, navigation, caching policy, footer.

Everything here is shared by every section and belongs to none of them. It lives
in its own module so that `app.py` can be a router and a section file can import
what it needs without importing the app.

Nothing in here computes anything about a sequence. That is `compute.py`.
"""

from __future__ import annotations

import dataclasses
import functools
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

__version__ = "2026.11.12d"

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
RELEASE = "1.1.0"

# The day this particular build was assembled. Two builds of the same release
# are possible - a fix pushed between releases - and this is what tells them
# apart in a bug report. Set it to the day you build, not the day of the
# release: a build stamped three weeks before the code it contains is a date
# that argues with the fingerprint beside it.
#
# Deliberately a constant and not today's date computed at import. The app is
# re-imported on every reboot, so a computed date would change without the
# code changing, and the one question this line has to answer, which build am
# I looking at, would get a different answer each morning.
BUILD = "20261112"

# Names that can appear inside data/ without being part of the app: caches an
# interpreter writes, and the file macOS leaves in every folder it has looked
# at. Hashing them would make the same commit fingerprint differently on two
# machines, which is the one thing the fingerprint must not do.
_FINGERPRINT_SKIP = {"__pycache__", ".DS_Store", ".ipynb_checkpoints"}


def _source_fingerprint(root: Path | None = None) -> str:
    """Seven hex characters that change whenever any source file changes.

    The build date alone cannot answer the question it is there for. Two pushes
    on the same day carry the same date, so a user looking at the footer after a
    deploy has no way to tell whether the new code is actually live or whether
    Streamlit is still serving the old container. A hash of the sources changes
    exactly when the sources change, needs no bookkeeping, and is the same in
    every environment: locally, in Docker, and on Community Cloud.

    What is hashed is everything the app serves: the Python files, plus
    `MANUAL.md` and the whole of `data/`, which the manual, library, gene
    diagnosis, reference and cloning screens read at runtime. A change to any
    of them is visible to a user, so it has to move the fingerprint.

    `data/` is walked rather than enumerated. It used to be three entries
    written out by hand, `library.csv`, `aminoacids/meta.json` and the SVG
    structures, which left `ecoli_codon_usage.csv`, `vhh_position_freq.csv`,
    `vhh_position_meta.json`, `primers.csv` and the four plasmid maps under
    `data/plasmids/` unhashed although the app reads all of them: a corrected
    codon usage table or a changed plasmid map went out under an unchanged
    footer. Walking the directory is the same rule as the `Dockerfile` copying
    the whole project instead of listing files, and for the same reason, a
    second list of what `data/` contains goes stale the first time something
    is added to it and says nothing when it does.

    The path is hashed along with the bytes, and it is the path relative to
    the root, not the bare file name: with the name alone, moving a file
    between two hashed directories left the fingerprint where it was.

    `requirements.txt` and `packages.txt` are hashed too, because Community
    Cloud installs them: a changed dependency line changes the deployed app.
    `requirements.txt` is a lockfile generated by pip-compile, every Python
    package pinned, transitive ones included, so a new Streamlit or Biopython
    release can no longer slip in under an unchanged fingerprint: upgrading
    means recompiling, and that changes the file. What the fingerprint still
    cannot see is what sits outside that file: the Python interpreter
    Community Cloud picked, and whatever version of HMMER apt installs from
    `packages.txt`, which names a package and no version. The footer line
    from `runtime_versions()` is there for those. `requirements.in`, the
    file pip-compile reads, is not hashed: Community Cloud never reads it,
    and any change to it that matters reaches `requirements.txt`.

    The test for what belongs here is whether a change to the file changes
    what a user sees online. Files that only a developer or another host
    reads, `README.md`, `CHANGELOG.md`, `tools/`, `Dockerfile` (Community
    Cloud ignores it) and `requirements-dev.txt` among them, are deliberately
    left out: including them would make the fingerprint move for deploys with
    nothing to verify, which is the one thing it must not do.
    """
    import hashlib

    # Path() around it because this is called with a plain string whenever
    # someone compares two checkouts, and a str has no .glob: the except below
    # swallowed the AttributeError and answered "unknown", which looks exactly
    # like a checkout with no source files in it.
    root = Path(root) if root else Path(__file__).parent
    digest = hashlib.sha256()
    try:
        # Everything under data/ that is a file, minus the droppings a file
        # manager or an interpreter leaves there, which are not part of the
        # app and would make the fingerprint differ between two checkouts of
        # the same commit.
        data_files = [p for p in (root / "data").rglob("*")
                      if p.is_file()
                      and not any(part in _FINGERPRINT_SKIP
                                  for part in p.relative_to(root).parts)]
        files = (sorted(root.glob("*.py"))
                 + sorted((root / "sections").glob("*.py"))
                 + sorted(data_files, key=lambda p: p.relative_to(root).parts)
                 + [p for p in (root / "MANUAL.md",
                                # Hidden, so the easiest file in the project to
                                # leave behind when copying a folder, and it
                                # sets the upload cap. Hashing it means a
                                # deployment missing it reads differently in
                                # the footer from the copy it was made of.
                                root / ".streamlit" / "config.toml",
                                # Installed by Community Cloud on each
                                # rebuild. See the docstring for what this
                                # does and does not catch.
                                root / "requirements.txt",
                                root / "packages.txt")
                    if p.exists()])
        for path in files:
            # posix() and not str(): the separator differs between Windows and
            # Linux, and a fingerprint that depends on the operating system
            # cannot be compared between a local checkout and the deployed app.
            digest.update(path.relative_to(root).as_posix().encode())
            digest.update(path.read_bytes())
    except Exception:
        return "unknown"
    return digest.hexdigest()[:7] if files else "unknown"


CODE = _source_fingerprint()

VERSION = f"{RELEASE} (build {BUILD}.{CODE})"


# The author as the work should be attributed, not as a platform handle: the
# citation metadata, the licence and the Zenodo record all carry the real name,
# and a footer that disagrees with them is a small way of making the tool harder
# to credit.
AUTHOR = "Jonathan Vaneyck"


AUTHOR_URL = "https://orcid.org/0000-0003-3095-7420"


# Where a user can write privately. It is the address already public on the
# author's GitHub profile, and the same one CITATION.cff carries; a test checks
# the two agree, so the app cannot advertise an address the citation does not.
# Offered beside the issue tracker because an issue needs a GitHub account,
# which many of the people this app is for do not have.
CONTACT_EMAIL = "jonathan.vaneyck@gmail.com"


REPO_URL = "https://github.com/TyrosinDuck/vhh-analyzer"


# The Zenodo *concept* DOI: it always resolves to the latest archived release,
# which is the one that belongs in a reference list. Each release also has a
# version DOI of its own, for citing the exact code that produced a result.
DOI = "10.5281/zenodo.23742269"


DOI_URL = f"https://doi.org/{DOI}"


LICENCE = "MIT"


ISSUES_URL = f"{REPO_URL}/issues"


def report_url(screen: str = "") -> str:
    """A GitHub "new issue" link with the version and the questions pre-filled.

    A user who hits a defect is in the app, not in the repository, and the thing
    they are least likely to supply unprompted is the one thing that makes a
    report actionable: which version they were running. Putting it in the link
    removes that step rather than asking for it.
    """
    from urllib.parse import quote

    body = (
        f"**Version:** {VERSION}\n"
        f"**Screen:** {screen or '(which screen and panel?)'}\n\n"
        "**What I expected**\n\n\n"
        "**What happened**\n\n\n"
        "**How to reproduce it**\n"
        "Please use a published sequence if you can - this issue is public, so "
        "nothing unpublished or confidential should go in it.\n")
    return f"{REPO_URL}/issues/new?title={quote('')}&body={quote(body)}"


def usage_url(screen: str = "") -> str:
    """A GitHub "new issue" link inviting a user to say how they use the app.

    The app records nothing about what is done in it, so the hosting
    platform's count of visitors is the only usage figure there is, and a
    visit is not a use. This is the one channel through which the author can
    learn whether a section is used, by whom and for what. It is the same
    pre-filled mechanism as `report_url`, rather than an issue template, so
    the version comes with it and a template cannot change what the existing
    bug link opens.

    ``screen`` names where the link was clicked, when it was offered under a
    result rather than on the home page, and is written into the first answer.
    """
    from urllib.parse import quote

    title = f"I use {APP_NAME} for ..."
    body = (
        f"**Version:** {VERSION}\n\n"
        "**Which sections do you use?**\n"
        + (f"Asked from: {screen}\n" if screen else "") + "\n\n"
        "**What for?** (kind of project, how often)\n\n\n"
        "**What is missing or gets in the way?** (optional)\n\n\n"
        "This issue is public. No sequence is needed, and nothing unpublished "
        "or confidential should go in it. To say it privately instead, write "
        f"to {CONTACT_EMAIL}.\n")
    return f"{REPO_URL}/issues/new?title={quote(title)}&body={quote(body)}"


def usage_mailto(screen: str = "") -> str:
    """The private counterpart of `usage_url`: an e-mail with the same prompts."""
    from urllib.parse import quote

    subject = f"{APP_NAME}: how I use it"
    body = (f"Version: {VERSION}\n\n"
            "Which sections do you use?\n"
            + (f"Asked from: {screen}\n" if screen else "") + "\n\n"
            "What for?\n\n\n"
            "What is missing or gets in the way? (optional)\n")
    return f"mailto:{CONTACT_EMAIL}?subject={quote(subject)}&body={quote(body)}"


def usage_prompt(screen: str) -> None:
    """One quiet line under a finished result, asking how the app is used.

    The home page panel is collapsed at the bottom of a page people leave as
    soon as they have picked a section, so almost nobody opens it. The moment
    a user has a reason to answer is right after the app gave them something
    they will keep, a report or a batch sheet, so the request is repeated
    there, once per screen and in the same plain style as the captions around
    it, rather than made louder on the home page.
    """
    st.caption(
        "Used this in your work? A couple of lines on what for decide what "
        f"gets improved next: [tell us on GitHub]({usage_url(screen)}), or "
        f"[by e-mail]({usage_mailto(screen)}) to say it privately.")


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


def anarci_requirement() -> str:
    """Which half of the ANARCI requirement is missing, as a clause.

    `anarci_status()` answers whether a sequence can be numbered, which is all
    a screen needs in order to decide what to draw. It cannot say why not, and
    the two halves are installed differently: the Python package comes from
    `requirements.txt`, the HMMER binary it shells out to does not, it comes
    from `packages.txt` on Community Cloud and from the system package manager
    anywhere else. Telling someone "ANARCI is unavailable" when the package is
    sitting there and only `hmmscan` is missing sends them to the wrong file.

    Only worth calling once the probe has already failed. The last case, both
    halves present and the probe still failing, is not hypothetical: HMMER
    needs ANARCI's germline HMM files, which a partial install can leave out.
    """
    import importlib.util
    import shutil

    has_package = importlib.util.find_spec("anarci") is not None
    has_hmmer = shutil.which("hmmscan") is not None
    if not has_package and not has_hmmer:
        return ("neither the `anarci` package nor the HMMER binary it calls "
                "is installed here")
    if not has_package:
        return "the `anarci` package is not installed here"
    if not has_hmmer:
        return ("the `anarci` package is installed but HMMER is not, so "
                "`hmmscan` cannot be found")
    return ("both are installed, but numbering a test sequence still failed, "
             "which usually means ANARCI's germline HMM files are missing")


def go_to(target: str) -> None:
    """Switch section on the next rerun."""
    st.session_state["section"] = target


# Body text colours of Streamlit's own light and dark themes, so text drawn
# inside an iframe matches the text around it.
THEME_TEXT = {"light": "#31333F", "dark": "#fafafa"}


def follow_app_theme() -> str:
    """A <style> and <script> that make an iframe's text follow the app theme.

    An iframe inherits no colour from the page around it, and the app can be
    light or dark. The iframes used to choose with
    ``@media (prefers-color-scheme: dark)``, which follows the *operating
    system*, not the app: with the system in dark mode and the app in light,
    which is what every visitor got at night while the theme was forced to
    light, the text turned near-white on a white page and the home page
    diagram became unreadable; choosing Dark in the menu on a light system
    does the reverse. Seen in a headless browser in both cases.

    Two layers, because neither is enough alone, both measured with Streamlit
    1.66.1:

    - the initial colour comes from ``st.context.theme.type``, which is right
      on a first load but stays at the old value after the user switches
      theme in the menu, until something reruns the script;
    - a script then reads the background of the page around the iframe and
      sets the text colour from it, on load and every half second, which
      follows a switch at once. It needs the iframe to share the page's
      origin, as Streamlit's component iframes do; if that ever fails, the
      initial colour stays, which is right until the next theme switch.

    "Dark" means a background whose weighted RGB mean, with the Rec. 709
    luma weights, is below half scale: a threshold for telling Streamlit's
    two backgrounds apart (#ffffff and #0e1117), not a measure of anything.
    """
    try:
        kind = st.context.theme.type
    except Exception:       # outside a running app, e.g. in a unit test
        kind = None
    initial = THEME_TEXT["dark" if kind == "dark" else "light"]
    return (
        f"<style>body{{margin:0;background:transparent;color:{initial}}}</style>"
        "<script>(function(){"
        "function sync(){try{"
        "var app=window.parent.document.querySelector('.stApp');"
        "if(!app)return;"
        "var c=getComputedStyle(app).backgroundColor.match(/\\d+/g).map(Number);"
        "var dark=(0.2126*c[0]+0.7152*c[1]+0.0722*c[2])<128;"
        f"document.body.style.color=dark?'{THEME_TEXT['dark']}':'{THEME_TEXT['light']}';"
        "}catch(e){}}"
        "if(document.body){sync();}else{addEventListener('DOMContentLoaded',sync);}"
        "setInterval(sync,500);"
        "})();</script>")


def embed_html(html: str, height: int) -> None:
    """Put a self-contained HTML fragment in an iframe.

    Three things in this app need a real iframe rather than ``st.markdown``:
    the flow diagram on the home page, the amino acid structures, and the
    3Dmol viewer. Streamlit's markdown sanitiser keeps SVG elements but
    strips the attributes that position them, so an inline drawing comes out
    as an empty box, and the viewer needs its own script context.

    Which call does that depends on the Streamlit version. ``st.iframe``
    replaced ``st.components.v1.html``, which Streamlit says it will remove,
    and the deployed app installs whatever version is current while a
    development machine may be behind. Choosing here, once, is what keeps
    both working, and keeps the deprecation warning out of the logs on the
    version that has the replacement.

    ``height`` stays an explicit number rather than ``st.iframe``'s
    "content": the structures and the diagram are SVG with no intrinsic
    height, so measuring them gives a flat panel.
    """
    if hasattr(st, "iframe"):
        st.iframe(html, height=height)
    else:                                     # Streamlit older than 1.60-ish
        import streamlit.components.v1 as components
        components.html(html, height=height, scrolling=False)


@st.cache_data(show_spinner=False, ttl=3600, max_entries=2)
def load_manual() -> str:
    """The user manual, read from MANUAL.md so there is one source of truth."""
    path = Path(__file__).parent / "MANUAL.md"
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        return f"The manual file could not be read ({exc})."


def current_sequences() -> dict[str, str]:
    """Sequences currently in the analysis input, parsed, or an empty dict.

    The cloning section needs them without depending on the analysis section
    having been rendered, so the text box is the single source of truth.
    """
    try:
        return vc.parse_input(st.session_state.get("seq_text", ""))
    except Exception:
        return {}


@functools.lru_cache(maxsize=1)
def runtime_versions() -> dict[str, str]:
    """What this server process is actually running, read at run time.

    The fingerprint says which files were deployed. It cannot say which
    interpreter Community Cloud chose or what apt put under ``packages.txt``,
    and before ``requirements.txt`` was pinned it could not say which
    Streamlit or Biopython pip had picked either. This answers that question
    from the process itself, so the footer of the deployed app can be read
    instead of guessed.

    Streamlit and Biopython are read from the modules already in memory
    (``__version__``), which is what this process runs even if a newer copy
    has since been installed on disk. ANARCI is read from the installed
    package metadata instead, because importing it only to ask its version
    costs about 150 ms on a first page that does not otherwise need it. It
    reads ``"absent"`` when the package is not installed. It says nothing
    about HMMER: an installed ANARCI without ``hmmscan`` still shows a
    version here, and the analysis sidebar is where its availability is
    tested.

    Cached for the life of the process: none of these can change without a
    restart, and the footer is drawn on every rerun.
    """
    import platform
    from importlib import metadata

    import Bio

    try:
        anarci_version = metadata.version("anarci")
    except metadata.PackageNotFoundError:
        anarci_version = "absent"
    return {"Python": platform.python_version(),
            "Streamlit": st.__version__,
            "Biopython": Bio.__version__,
            "ANARCI": anarci_version}


def footer() -> None:
    """Page footer, shown on every screen."""
    st.divider()
    st.caption("Built on Biopython (ProtParam), ANARCI, FAMSA, pymsaviz and "
               "dna_features_viewer. Theoretical predictions, validate "
               "experimentally. Sequences are not stored; for confidential work, "
               "run the app locally.")
    st.caption(
        f"{APP_NAME} · v{VERSION} · by [{AUTHOR}]({AUTHOR_URL}) · "
        f"{LICENCE} licence · [source]({REPO_URL}) · [doi:{DOI}]({DOI_URL}) · "
        f"[report a problem]({ISSUES_URL}) · "
        f"[using it? tell us]({usage_url()})")
    st.caption("Running on " + " · ".join(
        f"{name} {version}" for name, version in runtime_versions().items()))


# The sections are grouped by the molecule they work on, which is a criterion a
# user can apply without reading anything: the amino-acid sequence and its fold,
# the DNA that encodes it, and the purified protein in a tube. The groups are
# not stages of a workflow, data crosses them in both directions (the flow
# diagram on the home page shows how), they are only where to look.
#
# Each group carries the colour used for it on the home page, in the sidebar
# and in the flow diagram, so the three cannot disagree. The blue is the one
# the analysis figures and the theme already use. The green and the orange
# started as the next two colours of that palette, #55a868 and #dd8452, and
# were darkened, same hue and saturation, until white text on them reaches the
# WCAG 2.1 AA minimum of 4.5:1 (4.51 and 4.53; white on the blue is 4.84):
# the Open buttons on the home page are filled with these colours and carry
# white text, and at the original lightness that text measured 2.92 and 2.80.
# A test checks the ratio, so a later change of colour cannot quietly make the
# buttons unreadable. Colour is never the only carrier: every place that uses
# it also writes the group name.
#
# A group's members are not listed here. Each entry of SECTIONS names its own
# group, and `sections_by_group()` derives the membership, so adding a section
# means writing one entry, and a section cannot be silently left out of the
# grouping: an entry with no group, or a misspelt one, makes
# `sections_by_group()` raise instead of drawing the section nowhere.
GROUPS = {
    "Protein": ("#4c72b0", "The amino-acid sequence and its fold"),
    "DNA": ("#438452", "Reads, coding sequences and plasmids"),
    "Bench": ("#bb5b25", "The tube, and the readings taken on it"),
}

# Ordered by group, then by the order of use within a group, because the home
# page and the sidebar both list them in this order.
SECTIONS = {
    "analysis": ("🔬", "Sequence analysis",
                 "Characterize VHHs you provide: physico-chemical properties, "
                 "CDR/framework annotation, hallmarks and germline, liabilities, "
                 "alignment, codon-optimised DNA and a PDF report.",
                 "Protein"),
    "structure": ("🧊", "Structure",
                  "Upload a PDB, mmCIF or AlphaFold model and replace the "
                  "sequence-based guesses with measurements: which cysteines "
                  "actually pair, which liabilities are buried and which are "
                  "exposed, where the hydrophobic surface is, and, for a "
                  "complex, which residues contact the antigen.",
                  "Protein"),
    "library": ("📚", "Sequence library",
                "The sequences you come back to, kept in the repository and one "
                "click away: send the proteins to the analysis section or the "
                "DNA to gene diagnosis, instead of pasting them again.",
                "Protein"),
    "sequencing": ("🧪", "Sanger sequencing",
                   "Upload .ab1 chromatograms from a screening run, extract the VHH "
                   "sequences from the expression cassette, spot frameshifts and "
                   "duplicates, then send the clones you want to the analysis.",
                   "DNA"),
    "genes": ("🔧", "Gene diagnosis",
              "Paste coding sequences that are already in a plasmid and ask what "
              "might be wrong with them: frame and stop checks, documented rare "
              "codons and their clusters, GC stretches, plus synonymous swaps "
              "that lower GC without touching the protein.",
              "DNA"),
    "cloning": ("🧬", "In-silico cloning",
                "Load a backbone and an insert, and get the restriction enzymes that "
                "open the vector in the right place without cutting the VHH, plus "
                "which of your primers actually bind.",
                "DNA"),
    # Not called "Bench" although its key is: the group it sits in is, and a
    # section with the same name as its group read as a duplicate on the home
    # page.
    "bench": ("⚗️", "Batch & QC",
              "The arithmetic between a purified protein and the next "
              "experiment: concentration from a measured A280, whether the "
              "measured mass is the protein you think it is, how much "
              "labelling reagent to add, what each purification step "
              "recovered, and the batch sheet that leaves with the tube.",
              "Bench"),
    "buffers": ("🧂", "Buffers & media",
                "What to weigh for a buffer or a medium at any volume and "
                "concentration factor, with the hydrate of each salt you "
                "actually have: the mass follows from its formula. The "
                "recipes are the laboratory's own, kept in the repository.",
                "Bench"),
    "phage": ("🦠", "Phage & panning",
              "The phage stocks that go into each round of panning: the "
              "titre from an A260, by two routes shown side by side, and how "
              "much stock a solid-phase or liquid-phase round takes. No "
              "sequence needed.",
              "Bench"),
}

# The two screens that belong to no group: what you consult rather than what
# you run. They are not in SECTIONS because the home page builds a card per
# section and the grouping covers every section, and these two would have to
# be excluded from both. They are declared here all the same, rather than
# written out in the sidebar, so that the router can reach them the same way
# it reaches a section: `app.py` knows every screen from SECTIONS, UTILITIES
# and the home page, and holds no list of its own.
UTILITIES = {
    "manual": ("📖", "Manual"),
    "reference": ("📐", "Reference"),
    "provenance": ("🔎", "Provenance"),
}

# The key of the landing page. Named rather than spelt out in the three places
# that need it (the default section, the Home button, the fallback for a
# section key that no longer exists).
HOME = "home"


def screen_keys() -> list[str]:
    """Every screen the router can draw, in the order the sidebar lists them.

    This is the whole navigable surface of the app, derived from the three
    declarations above rather than repeated: a screen is the home page, a
    section, or a utility. `app.py` imports one module per key in this list,
    checks its version and calls its `render`, which is why adding a section
    is one entry in `SECTIONS` plus one file in `sections/` and nothing else.
    """
    return [HOME, *SECTIONS, *UTILITIES]


@dataclasses.dataclass(frozen=True)
class Settings:
    """The analysis settings the sidebar collects, as a single argument.

    Every screen's `render` takes this and nothing else, whether or not it
    reads any of it. The alternative, each screen taking the arguments it
    happens to need, is what made the router a second list beside `SECTIONS`:
    the calls differed, so `app.py` carried one `if` per section, and adding a
    section meant editing the router too. Nothing checked the router, and a
    section missing from it is a section that cannot be opened.

    A screen that ignores the argument pays one parameter name for that, which
    is cheaper than the list. `test_navigation.py` checks that every screen
    takes exactly this, so a screen whose signature drifts fails a test
    instead of failing a click.

    Frozen: a screen must not be able to change what the next one sees. These
    belong to the sidebar that collected them.
    """

    ph: float
    window: int
    method: str

    @classmethod
    def defaults(cls) -> Settings:
        """The settings a screen gets when the sidebar is not showing them.

        The sliders are only drawn on the analysis screen, where they matter
        and where their effect is visible. Every other screen that reads a pH
        or an annotation method reads the same defaults `vhh_core` uses, so a
        value never depends on whether the analysis screen was opened first.
        """
        return cls(ph=vc.DEFAULTS["ph"],
                   window=vc.DEFAULTS["profile_window"],
                   method=vc.DEFAULTS["annotation_method"])


def sections_by_group() -> dict[str, list[str]]:
    """Section keys per group, groups in `GROUPS` order, sections in
    `SECTIONS` order.

    Derived on every call rather than stored, so that it cannot disagree with
    `SECTIONS`. A group with no section is left out instead of drawn empty; a
    section naming a group that does not exist raises, because drawing it
    nowhere would make it unreachable without any message.
    """
    unknown = {group for *_, group in SECTIONS.values()} - set(GROUPS)
    if unknown:
        raise KeyError(f"SECTIONS names unknown group(s): {sorted(unknown)}")
    grouped = {name: [key for key, (*_, group) in SECTIONS.items()
                      if group == name]
               for name in GROUPS}
    return {name: keys for name, keys in grouped.items() if keys}


def contrast_ratio(colour_a: str, colour_b: str) -> float:
    """WCAG 2 contrast ratio between two "#rrggbb" colours, from 1 to 21.

    Relative luminance and contrast ratio as defined in the glossary of
    W3C, "Web Content Accessibility Guidelines (WCAG) 2.1",
    https://www.w3.org/TR/WCAG21/ (definitions *relative luminance* and
    *contrast ratio*; Success Criterion 1.4.3 asks for at least 4.5:1 for
    normal-size text). The linearisation threshold is 0.04045, the value
    since the May 2021 erratum; it was 0.03928 before, which the
    specification notes has no practical effect, and indeed no 8-bit channel
    value lies between the two.
    """
    def luminance(colour: str) -> float:
        channels = [int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)]
        r, g, b = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
                   for c in channels)
        return 0.2126 * r + 0.7152 * g + 0.0722 * b

    high, low = sorted((luminance(colour_a), luminance(colour_b)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def group_header_html(name: str, size: str) -> str:
    """The coloured group title used by the home page and the sidebar.

    One function for both so they cannot drift apart. `size` is a CSS font
    size: the home page and the sidebar want different ones, but the same
    shape. The name is written out, the colour only underlines it.
    """
    colour, tagline = GROUPS[name]
    return (f'<div style="border-left:5px solid {colour};padding-left:0.6rem;'
            f'margin:0.4rem 0 0.5rem 0">'
            f'<div style="font-size:{size};font-weight:700;line-height:1.2">'
            f'{name}</div>'
            f'<div style="font-size:0.85rem;opacity:0.7">{tagline}</div></div>')
