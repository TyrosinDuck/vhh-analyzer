"""
Regression tests, every bug this app has actually had, frozen so it cannot
come back.

These are not here to prove the app is correct. They are here because each one
of them was a real defect that shipped, was noticed by accident, and was fixed.
A test is the only thing that stops a fix from quietly coming undone the next
time the function around it is edited.

Run them from the app folder:

    pip install -r requirements-dev.txt
    pytest

Each test names the bug it guards in its docstring. If one fails, read that
docstring first: it tells you what the wrong behaviour looked like.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

import codons as cdn
import library as lib
import repertoire as rp
import triage as tg
import vhh_core as vc

ROOT = Path(__file__).resolve().parent.parent

# ANARCI needs the HMMER binary, which is not installed everywhere. The regex
# annotator covers canonical scaffolds, so the app works without it, and the
# tests that pin an ANARCI-specific answer skip rather than fail: a contributor
# who cannot install HMMER should still be able to run most of the suite and
# see that their change broke nothing.
needs_anarci = pytest.mark.skipif(
    not vc.anarci_available(),
    reason="ANARCI (and its HMMER binary) is not installed here")


# --------------------------------------------------------------------------- #
# Fixtures, real sequences, from the shipped library, not invented ones
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="session")
def library():
    table, problems = lib.load_library(ROOT / "data" / "library.csv")
    assert not problems, f"the shipped library has problems: {problems}"
    return table


@pytest.fixture(scope="session")
def by_name(library):
    return {row["name"]: row for _, row in library.iterrows()}


def _protein(by_name, substring):
    """The one library entry whose name contains ``substring``."""
    hits = [row for name, row in by_name.items()
            if substring.lower() in name.lower() and row["protein"]]
    assert len(hits) == 1, f"{substring!r} matched {len(hits)} entries, expected 1"
    return hits[0]["protein"]


@pytest.fixture(scope="session")
def cablys3(by_name):
    return _protein(by_name, "cAbLys3")


@pytest.fixture(scope="session")
def lysozyme(by_name):
    return _protein(by_name, "lysozyme")


# --------------------------------------------------------------------------- #
# 1. The module-version banner
# --------------------------------------------------------------------------- #

def test_module_versions_match_app():
    """The staleness banner fires when a module is older than app.py expects.

    It is a real safeguard (the modules are released together) but it is also
    the single most common way to ship a broken deploy: bump a module, forget
    the dict in app.py, and the app refuses to start with a message about a file
    you did update. This test reads both and compares them, so that mismatch is
    caught here instead of in the browser.
    """
    source = (ROOT / "app.py").read_text()
    block = re.search(r"REQUIRED_MODULES\s*=\s*\{(.*?)\}", source, re.S)
    assert block, "REQUIRED_MODULES not found in app.py"
    required = dict(re.findall(r'"(\w+)":\s*"([^"]+)"', block.group(1)))
    assert required, "REQUIRED_MODULES parsed as empty"

    import importlib
    for name, want in required.items():
        module = importlib.import_module(name)
        assert module.__version__ == want, (
            f"{name}.py is {module.__version__} but app.py wants {want}, "
            f"update REQUIRED_MODULES in app.py, or the module's __version__")


# --------------------------------------------------------------------------- #
# 2. Generic proteins must not be judged as antibodies
# --------------------------------------------------------------------------- #

def test_lysozyme_cysteines_are_not_orphan_thiols(lysozyme):
    """Bug: hen lysozyme's 8 perfectly ordinary cysteines were reported as 8
    free thiols.

    The canonical-pair logic assumes a VHH: one disulfide between IMGT 23 and
    104, anything else unpaired. Run on a non-antibody it produced confident
    nonsense. The rule for a sequence with no annotation is now parity: an even
    number of cysteines is not evidence of a free thiol.
    """
    assert lysozyme.count("C") == 8, \
        "the fixture is no longer hen lysozyme; update this test"
    liabilities = vc.find_liabilities(lysozyme, annotation=None)
    free = liabilities[liabilities["liability"].str.contains("thiol", case=False)]
    assert free.empty, f"free thiols reported on lysozyme: {free.to_dict('records')}"


def test_lysozyme_is_not_classified_as_an_antibody(lysozyme):
    """A non-antibody must not get the antibody analyses at all."""
    kind = vc.classify(lysozyme)
    assert kind["kind"] == vc.KIND_OTHER
    assert not kind["is_vhh"]
    assert not kind["full_analysis"]


def test_triage_leaves_antibody_columns_empty_for_a_non_antibody(lysozyme,
                                                                 cablys3):
    """Bug: the same canonical-pair leak, a second time, in the triage table.

    Fixing ``find_liabilities`` was not enough: the triage table called the
    disulfide, hallmark and humanness analyses directly. They are now gated on
    ``full_analysis``, and the cells are left empty rather than answered wrongly.
    """
    seqs = {"lysozyme": lysozyme, "cAbLys3": cablys3}
    annotations = {}
    for name, seq in seqs.items():
        try:
            annotations[name], _ = vc.annotate(seq)
        except Exception:
            annotations[name] = None
    kinds = {name: vc.classify(seq) for name, seq in seqs.items()}

    table = tg.triage_table(seqs, annotations, kinds).set_index("name")
    for column in ("n_orphan_cys", "hallmark_score", "humanness"):
        assert pd.isna(table.loc["lysozyme", column]), (
            f"{column} answered for a non-antibody: "
            f"{table.loc['lysozyme', column]!r}")
    assert pd.notna(table.loc["cAbLys3", "hallmark_score"]), \
        "the VHH lost its hallmark score, the gate is too strict"


# --------------------------------------------------------------------------- #
# 3. Repairs must not break the fold
# --------------------------------------------------------------------------- #

@needs_anarci
def test_repairs_never_touch_a_structurally_conserved_position(cablys3):
    """Bug: the repair suggester proposed W36F and W111F.

    Those are the FR2 and FR4 tryptophans (IMGT 41 and 118) which are
    conserved across every immunoglobulin domain. "Repairing" an oxidation
    liability by removing them trades a chemical risk for a protein that does
    not fold. The four conserved positions are now refused outright.
    """
    liabilities = vc.find_liabilities(cablys3)
    annotation, _ = vc.annotate(cablys3)
    repairs = vc.suggest_repairs(cablys3, liabilities, annotation)
    if not len(repairs):
        pytest.skip("no repairs suggested for this sequence")

    numbering = {i: num for num, ins, _aa, i in vc.imgt_numbering(cablys3)
                 if not ins}
    offered = repairs[repairs["substitution"] != "-"]
    for _, row in offered.iterrows():
        position = int(row["position"])
        imgt = numbering.get(position - 1)
        assert imgt not in vc.CONSERVED_IMGT, (
            f"substitution {row['substitution']} offered at IMGT {imgt}, "
            f"which is conserved ({vc.CONSERVED_IMGT.get(imgt)})")


def test_apply_repairs_refuses_two_substitutions_at_one_position(cablys3):
    """Bug: the default selection offered C→S and C→A at the same position.

    Both are reasonable on their own and applying both is meaningless. The
    conflict is a decision for the user, so it is refused rather than resolved.
    """
    first_cys = cablys3.find("C") + 1
    assert first_cys > 0
    sequence, problems = vc.apply_repairs(
        cablys3, [f"C{first_cys}S", f"C{first_cys}A"])
    assert problems, "a conflicting pair of substitutions was applied silently"
    assert sequence == cablys3, (
        "one of two conflicting substitutions was applied anyway, the order "
        "of the list must not decide a question the caller has not answered")


def test_apply_repairs_refuses_a_substitution_whose_residue_does_not_match(cablys3):
    """A substitution is only meaningful if the original residue is there."""
    wrong = "W1Q" if cablys3[0] != "W" else "Q1W"
    sequence, problems = vc.apply_repairs(cablys3, [wrong])
    assert problems and sequence == cablys3


# --------------------------------------------------------------------------- #
# 4. Empty tables must keep their shape
# --------------------------------------------------------------------------- #

def test_liabilities_table_keeps_its_columns_when_empty():
    """Bug: a clean sequence gave a DataFrame with no columns, and the repair
    panel died on ``KeyError: 'score'``.

    An empty result is a normal outcome and must have the same shape as a
    non-empty one, or every consumer needs its own special case.
    """
    clean = "GSGSGSGSGSGSGSGSGSGSGSGSGSGSGS"
    table = vc.find_liabilities(clean)
    assert list(table.columns) == vc.LIABILITY_COLUMNS
    assert table.empty


# --------------------------------------------------------------------------- #
# 5. IMGT numbering semantics
# --------------------------------------------------------------------------- #

@needs_anarci
def test_imgt_position_map_keeps_insertion_coded_positions(cablys3):
    """Regression from a refactor, not from a user report.

    Caching the ANARCI call changed ``imgt_position_map`` so that an insertion
    code (111A, 112B) no longer overwrote its base number. That was arguably an
    improvement, but it silently changed every region boundary downstream, so
    the original behaviour was restored. This test pins it: the map is keyed by
    the bare IMGT number, and an insertion shares its base key.
    """
    mapping = vc.imgt_position_map(cablys3)
    assert mapping, "no numbering produced, is ANARCI installed?"
    assert all(isinstance(key, int) for key in mapping), \
        "the map must be keyed by bare integers, not by '111A'-style strings"
    for position in (23, 104):
        assert mapping.get(position) == "C", \
            f"IMGT {position} should be the canonical cysteine, got " \
            f"{mapping.get(position)!r}"


@needs_anarci
def test_cablys3_cdr3_is_what_the_pdb_entry_says(cablys3):
    """A frozen annotation, checked against PDB 1MEL when the entry was added.

    If this changes, either the annotator changed or the library entry did.
    Both are worth knowing about.
    """
    annotation, method = vc.annotate(cablys3, "anarci")
    assert "ANARCI" in method
    assert annotation["CDR3"]["seq"] == "AADSTIYASYYECGHGLSTGGYGYDS"


# --------------------------------------------------------------------------- #
# 6. Back-translation must be reversible
# --------------------------------------------------------------------------- #

def test_reverse_translation_round_trips_for_every_library_protein(library):
    """The property that must hold for every sequence, not just the tested one.

    The optimiser removes restriction sites, breaks homopolymers and balances
    GC. Any of those steps could in principle change an amino acid. Translating
    the result back and comparing is the only check that actually proves it did
    not.
    """
    for _, row in library.iterrows():
        protein = row["protein"]
        if not protein:
            continue
        dna = vc.reverse_translate(protein)["dna"]
        back = vc.translate_dna(dna).rstrip("*")
        assert back == protein, f"{row['name']}: round trip changed the protein"


def test_library_dna_translates_to_its_protein(library):
    """Where an entry carries both, they must agree.

    ``load_library`` already flags a mismatch, and the ``library`` fixture
    asserts there are no problems, this states the invariant explicitly so the
    reason is obvious when it breaks.
    """
    both = library[(library["has_protein"]) & (library["has_dna"])]
    for _, row in both.iterrows():
        back = vc.translate_dna(row["dna"]).rstrip("*")
        assert back.lstrip("M") == row["protein"].lstrip("M") \
            or back == row["protein"], f"{row['name']}: DNA and protein disagree"


# --------------------------------------------------------------------------- #
# 7. Lowering GC must not create a codon-usage problem
# --------------------------------------------------------------------------- #

def test_gc_fix_preserves_the_protein_and_respects_the_usage_floor(cablys3):
    """The explicit requirement: do not trade a GC problem for a rare codon.

    Every substituted codon must still be one *E. coli* actually uses, at or
    above ``codons.MIN_USAGE`` of its synonym family's relative adaptiveness.
    """
    dna = vc.reverse_translate(cablys3)["dna"]
    gc = cdn.gc_windows(dna)
    regions = cdn.gc_regions(gc)
    result = cdn.suggest_gc_fix(dna, regions)

    assert result["protein_unchanged"], "the GC fix changed the protein"
    assert vc.translate_dna(result["dna"]).rstrip("*") == cablys3

    _usage, relative = cdn.codon_usage()
    suggestions = result["suggestions"]
    if not len(suggestions):
        pytest.skip("this gene needed no GC fix")
    for _, row in suggestions.iterrows():
        new = row["to"] if "to" in row else row.get("new_codon")
        if not new:
            continue
        amino = vc.translate_dna(new)
        weight = relative.get(amino, {}).get(new)
        assert weight is not None, f"{new} is not in the usage table"
        assert weight >= cdn.MIN_USAGE, (
            f"{new} ({amino}) has relative usage {weight:.2f}, below the "
            f"{cdn.MIN_USAGE} floor, a GC fix traded for a rare codon")
        assert new not in cdn.RARE_CODONS, f"{new} is a documented rare codon"


def test_gc_regions_do_not_overlap(cablys3):
    """Bug: 98 "extreme windows" were reported on the app's own optimised gene.

    Overlapping 30-nt windows were each counted as a separate finding. They are
    now merged into regions, and the merge itself had a bug that left
    overlapping spans, which is what this pins.
    """
    dna = vc.reverse_translate(cablys3)["dna"]
    regions = cdn.gc_regions(cdn.gc_windows(dna))
    if len(regions) < 2:
        pytest.skip("fewer than two regions to compare")
    ordered = regions.sort_values("start")
    previous_end = -1
    for _, row in ordered.iterrows():
        assert int(row["start"]) > previous_end, \
            f"region {row['start']}–{row['end']} overlaps the previous one"
        previous_end = int(row["end"])


# --------------------------------------------------------------------------- #
# 8. Single-sequence edge cases
# --------------------------------------------------------------------------- #

def test_cdr_length_comparison_works_with_one_sequence(cablys3):
    """Bug: with one sequence the bar reached 1.0 and flattened the reference.

    The plot now switches to markers below eight sequences, but the underlying
    comparison must produce a usable table either way.
    """
    annotation, _ = vc.annotate(cablys3)
    yours, reference = rp.cdr_length_comparison({"cAbLys3": annotation})
    assert len(yours) == 1, "no comparison produced for a single sequence"
    assert len(reference), "the natural-repertoire reference came back empty"

    # The lengths must be the annotation's own, not recomputed differently
    row = yours.iloc[0]
    for region in ("CDR1", "CDR2", "CDR3"):
        assert int(row[region]) == len(annotation[region]["seq"]), \
            f"{region} length disagrees with the annotation"


def test_clustering_at_100_percent_groups_only_identical_cdr3s(cablys3):
    """At a 100% threshold, two different CDR3s must never share a family."""
    other = cablys3.replace("AADSTIYASYYECGHGLSTGGYGYDS",
                            "AADSTIYASYYECGHGLSTGGYGYDW")
    assert other != cablys3, "the fixture sequence changed; update this test"
    seqs = {"a": cablys3, "b": other}
    annotations = {name: vc.annotate(seq)[0] for name, seq in seqs.items()}
    result = vc.cluster_families(seqs, annotations, cutoff=100.0)
    assignment = result["assignment"]
    families = (assignment.set_index("name")["family"].to_dict()
                if hasattr(assignment, "set_index") else dict(assignment))
    assert families["a"] != families["b"], \
        "two different CDR3s were grouped at a 100% identity threshold"


# --------------------------------------------------------------------------- #
# 9. Input validation
# --------------------------------------------------------------------------- #

def test_the_manual_splits_without_losing_anything():
    """The manual screen shows one section at a time, so a splitter that drops
    a section hides part of the documentation with nothing to see.

    Rendering 1 700 lines on every interaction is slow and nobody reads a
    reference manual from the top, but the price of cutting it up is that a
    silent loss is invisible: the page still looks complete.
    """
    import compute
    import shell

    text = shell.load_manual()
    headings = [l for l in text.split("\n") if l.startswith("## ")]
    sections = compute.manual_sections(text)

    # One per heading, plus the introduction that precedes the first one.
    assert len(sections) == len(headings) + 1, (
        f"{len(headings)} top-level headings in MANUAL.md plus an intro, but "
        f"{len(sections)} sections parsed")
    for section in sections:
        assert section["title"].strip(), "a section came back with no title"
        assert section["body"].strip(), (
            f"section {section['title']!r} came back empty")

    # Compared on non-blank lines: each section is stripped of the blank lines
    # at its edges, which is harmless, while a lost paragraph is not. Counting
    # every line would make the test fail on whitespace and pass on content.
    def substance(body: str) -> int:
        return sum(1 for line in body.split("\n") if line.strip())

    kept = sum(substance(s["body"]) for s in sections)
    total = substance(text) - len(headings)
    assert kept == total, (
        f"{total - kept} non-blank line(s) went missing between MANUAL.md and "
        f"the parsed sections")


@pytest.mark.parametrize("junk", ["", "   ", "BJOUXZ" * 10])
def test_clean_sequence_rejects_junk(junk):
    """Empty input and non-amino-acid letters are refused."""
    with pytest.raises(ValueError):
        vc.clean_sequence(junk)


def test_clean_sequence_strips_punctuation_without_judging_length():
    """It validates the alphabet, not the length.

    ``ACGT`` is four perfectly good amino acids; whether four residues is a
    useful input is a question for the caller, and the app enforces
    ``MIN_LENGTH`` where it matters. Stating this here stops someone later
    "fixing" clean_sequence to reject short peptides and breaking every caller
    that legitimately handles fragments.
    """
    assert vc.clean_sequence("ACGT123!!") == "ACGT"
    assert vc.MIN_LENGTH > 4, "if this changes, the comment above needs revising"


def test_clean_sequence_strips_whitespace_and_case(cablys3):
    messy = "  " + "\n".join(cablys3[i:i + 10].lower()
                             for i in range(0, len(cablys3), 10)) + "  "
    assert vc.clean_sequence(messy) == cablys3


def test_every_project_import_in_app_py_is_guarded():
    """``KeyError: 'generator'`` reached a visitor as a traceback mid-deploy.

    A deploy swaps the files under the running server, and a run that imports
    a module in the middle of the swap gets a ``KeyError``. The guard that
    turns it into "The app is updating" used to wrap only the version check,
    which runs after the imports at the top of `app.py`, so the import that
    actually raised was unprotected. The guard is now a ``try`` around those
    imports; this fails if a module of this project is imported at the top
    level of `app.py` outside it. Project modules are found from the files on
    disk, not listed here.
    """
    import ast

    local = ({p.stem for p in ROOT.glob("*.py")}
             | {p.name for p in ROOT.iterdir()
                if p.is_dir() and (p / "__init__.py").exists()})
    tree = ast.parse((ROOT / "app.py").read_text(encoding="utf-8"))

    def project_imports(nodes):
        for node in nodes:
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and not node.level:
                names = [node.module or ""]
            else:
                continue
            if any(name.split(".")[0] in local for name in names):
                yield node

    unguarded = [node.lineno for node in project_imports(tree.body)]
    assert not unguarded, (
        f"app.py imports project modules outside the update guard, at lines "
        f"{unguarded}; move them into the try block at the top")

    guards = [node for node in tree.body if isinstance(node, ast.Try)
              and list(project_imports(node.body))]
    assert len(guards) == 1, "expected one try block holding the imports"
    caught = ast.unparse(guards[0].handlers[0].type)
    assert "KeyError" in caught and "ImportError" in caught, \
        f"the import guard catches {caught}, it must catch KeyError and ImportError"


def test_no_cached_function_is_cached_twice():
    """`families_for` carried two `st.cache_data` decorators.

    One of them was `_FIG_CACHE`, copied down from the figure helper above it,
    which caps at 24 entries. The result was correct, so nothing failed: the
    clonal families and their linkage matrix, which is the largest object this
    module caches, were simply stored in two caches instead of one, in a
    container that can have as little as 690 MB. Stacked caches are never
    intentional here, so this looks for them structurally rather than for that
    one function.
    """
    import ast

    offenders = []
    for path in [ROOT / "compute.py", ROOT / "shell.py",
                 *sorted((ROOT / "sections").glob("*.py"))]:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            cached = [d for d in node.decorator_list
                      if "cache_data" in ast.unparse(d)
                      or "cache_resource" in ast.unparse(d)]
            if len(cached) > 1:
                offenders.append(f"{path.name}:{node.lineno} {node.name}")
    assert not offenders, (
        "these functions carry more than one Streamlit cache decorator, which "
        "stores every result twice: " + ", ".join(offenders))


def _identity_frame(n):
    names = [f"s{i:02d}" for i in range(n)]
    return pd.DataFrame(60.0, index=names, columns=names).where(
        ~pd.DataFrame(False, index=names, columns=names), 100.0)


def test_identity_heatmap_stays_bounded_at_the_sequence_limit():
    """At 50 sequences the identity heatmap took the deployed app down.

    The figure grew by 1.1 inch a sequence and wrote a number in every cell:
    58 x 47 inches and 2500 labels at the app's limit of 50 sequences, about
    110 million pixels at the 200 dpi the app renders at. Rendering it peaked
    at about 4 GB, measured, on a Community Cloud container that can have
    690 MB, so opening Comparison with 50 sequences restarted the server.
    """
    import matplotlib.pyplot as plt

    n = vc.DEFAULTS["max_sequences"]
    fig = vc.plot_identity(_identity_frame(n))
    try:
        assert max(fig.get_size_inches()) <= vc.IDENTITY_FIG_MAX + 1e-9
        ax = fig.axes[0]
        assert not ax.texts, "cells should carry no number past the label limit"
    finally:
        plt.close(fig)


def test_identity_heatmap_keeps_its_numbers_for_small_sets():
    """The cap must not cost the usual case its values."""
    import matplotlib.pyplot as plt

    for n in (5, vc.IDENTITY_LABEL_MAX):
        fig = vc.plot_identity(_identity_frame(n))
        try:
            assert len(fig.axes[0].texts) == n * n
        finally:
            plt.close(fig)


def test_every_figure_helper_is_cached():
    """`png_identity` was the one figure helper without a cache.

    Every rerun of the Comparison family, a slider moved, a checkbox ticked,
    recomputed the pairwise identities and redrew both heatmaps, which at 50
    sequences was also the most expensive figure in the app. Checked
    structurally so the next `png_` helper cannot be added without one.
    """
    import ast

    tree = ast.parse((ROOT / "compute.py").read_text(encoding="utf-8"))
    missing = [node.name for node in tree.body
               if isinstance(node, ast.FunctionDef)
               and node.name.startswith("png_")
               and not any("cache_data" in ast.unparse(d)
                           for d in node.decorator_list)]
    assert not missing, "figure helpers without a cache: " + ", ".join(missing)


def test_profiles_zip_holds_one_png_per_sequence():
    """Above PROFILES_COMBINED_MAX the profiles download is a ZIP.

    One tall figure of 50 panels took 516 MB to render, measured, on a server
    that can have 690 MB. The ZIP must carry every sequence, once, under a
    file name that survives odd characters in the sequence name.
    """
    import io
    import zipfile

    import compute

    frame = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
             "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYD{}WGQGTQVTVSS")
    seqs = {"Nb 1": frame.format("Y"), "Nb/2": frame.format("F"),
            "Nb:3": frame.format("W")}
    data = compute.profiles_zip(seqs, {}, 7)
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        names = zf.namelist()
        assert len(names) == len(set(names)) == 3
        assert all(n.startswith("profiles/") and n.endswith(".png") for n in names)
        assert all(zf.read(n)[:8] == b"\x89PNG\r\n\x1a\n" for n in names)
    assert compute.PROFILES_COMBINED_MAX < vc.DEFAULTS["max_sequences"]


def test_the_two_cdr3_identity_definitions_note_is_true():
    """The app says the Alignment and Clonal families tabs define CDR3
    identity differently, and shows one pair where they disagree.

    `pairwise_identity` divides by the alignment length, `pairwise_scores` by
    the longer CDR3, with slightly different gap penalties. The note claims
    they agree when one CDR3 lies within the other and differ when each
    overhangs the other; if a change to either function made the example agree,
    the note would show two equal numbers under a sentence saying they differ.
    """
    import compute

    a, b = compute.CDR3_IDENTITY_EXAMPLE
    assert vc.pairwise_identity(a, b) != vc.pairwise_scores(a, b)[0]

    contained = [("ARDGYSSSWYDY", "ARDGYSSWYDY"),      # internal deletion
                 ("AAGRWYTSGDYDY", "AAGRWYTSGDY")]     # truncated end
    for x, y in contained:
        assert vc.pairwise_identity(x, y) == vc.pairwise_scores(x, y)[0]


def test_both_extinction_coefficients_reach_the_a280_factor():
    """Sequence analysis gave only the reduced extinction coefficient.

    Batch & QC offered both, cystines formed and all reduced, and defaulted to
    the cystine value; the Properties table showed the reduced one and the
    A280 factor derived from it, so the two screens could convert the same
    A280 with different coefficients and nothing said so. ``characterize`` now
    returns the cystine A280 factor beside the reduced one, both from the same
    average mass.
    """
    vhh = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
           "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDYWGQGTQVTVSS")
    props = vc.characterize("vhh", vhh)
    assert props["n_Cys"] == 2
    assert props["ext_coeff_cystines"] > props["ext_coeff_reduced"]
    assert props["A280_0.1pct_cystines"] == pytest.approx(
        props["ext_coeff_cystines"] / props["mw_da"], abs=1e-3)
    assert props["A280_0.1pct"] == pytest.approx(
        props["ext_coeff_reduced"] / props["mw_da"], abs=1e-3)
    assert props["A280_0.1pct_cystines"] > props["A280_0.1pct"]

    no_cys = vc.characterize("no_cys", vhh.replace("C", "S"))
    assert no_cys["ext_coeff_cystines"] == no_cys["ext_coeff_reduced"]
    assert no_cys["A280_0.1pct_cystines"] == no_cys["A280_0.1pct"]


def test_the_molar_mass_keeps_two_decimals_in_dalton():
    """The Properties table rounded the mass to 10 Da, as kDa with 2 decimals.

    ``mw_da`` carries the average mass to 0.01 Da, which is what the table and
    the PDF report now show.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    seq = "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
    props = vc.characterize("x", seq)
    assert props["mw_da"] == round(ProteinAnalysis(seq).molecular_weight(), 2)
    assert props["mw_da"] != props["mw_kda"] * 1000


def test_properties_table_shows_two_decimals():
    """pI and the net charge read 8.620000 and 1.480000 in the Properties table.

    A second ``Styler.format`` with a dict reset every other column to pandas'
    six-decimal default. The table now formats in one call: two decimals by
    default, integers for the extinction coefficients, three decimals for the
    A280 factors, and no thousands separator.
    """
    from sections.analysis import properties_styler
    seq = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
           "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRGYSSSWYDYWGQGTQVTVSS")
    cols = ["name", "mw_da", "pI", "net_charge", "gravy", "instability_index",
            "ext_coeff_cystines", "ext_coeff_reduced", "A280_0.1pct"]
    frame = pd.DataFrame([vc.characterize("vhh", seq)])[cols]
    body = properties_styler(frame)._translate(True, True)["body"][0]
    shown = {cols[i]: cell["display_value"]
             for i, cell in enumerate(c for c in body if c.get("type") == "td")}
    for col in ("pI", "net_charge", "gravy", "instability_index", "mw_da"):
        assert re.fullmatch(r"-?\d+\.\d{2}", shown[col]), (col, shown[col])
    assert re.fullmatch(r"\d+", shown["ext_coeff_cystines"])
    assert re.fullmatch(r"\d\.\d{3}", shown["A280_0.1pct"])
    assert "," not in "".join(shown.values())


# --------------------------------------------------------------------------- #
# The November 2026 independent audit. Each of these was a claim the app made
# to its user that the code next to it did not support.
# --------------------------------------------------------------------------- #

def test_the_rare_codon_caption_counts_the_tRNAs_it_promises():
    """The Reference screen said all eight rare codons name a limiting tRNA.

    CGA and CGG carry ``"trna": None``, so two rows of the table showed an
    empty tRNA column under a caption asserting every row had one. The caption
    now says six, and this test fails if the table and the sentence drift
    apart again, in either direction.
    """
    named = [c for c, info in cdn.RARE_CODONS.items() if info["trna"]]
    unnamed = [c for c, info in cdn.RARE_CODONS.items() if not info["trna"]]
    assert len(cdn.RARE_CODONS) == 8
    assert len(named) == 6
    assert sorted(unnamed) == ["CGA", "CGG"]

    source = (ROOT / "sections" / "reference.py").read_text(encoding="utf-8")
    assert "Six of them carry the name of the tRNA" in source
    manual = (ROOT / "MANUAL.md").read_text(encoding="utf-8")
    assert "Six of them carry the" in manual


def test_the_cai_never_claims_the_sharp_and_li_reference_set():
    """``cai`` was documented as Sharp & Li, over a whole-organism table.

    Sharp & Li define the index against highly expressed genes. The shipped
    table is Kazusa's whole-organism one, so the arithmetic is theirs and the
    reference set is not. Every place that names the index has to say so, and
    the reverse-translation comment must stop claiming a highly-expressed set
    it never used.
    """
    doc = cdn.cai.__doc__
    assert "whole-genome reference set" in doc
    assert "highly expressed" in doc

    caption = (ROOT / "sections" / "genes.py").read_text(encoding="utf-8")
    assert "not against a reference set of highly expressed" in caption

    core = (ROOT / "vhh_core.py").read_text(encoding="utf-8")
    block = core[core.index("# Relative synonymous codon usage"):
                 core.index("ECOLI_CODONS = {")]
    assert "whole-organism table" in block
    # The claim that started it survives only as the note saying it used to be
    # there and was wrong, in quotation marks.
    assert block.count("highly expressed") == 1
    assert '"highly expressed genes", which contradicted' in block


def test_the_bench_screen_prefills_no_extinction_coefficient_or_mass():
    """28000 /M/cm and 14000 Da were prefilled, with no source.

    Typing an A280 with no sequence loaded produced an mg/mL indistinguishable
    from one computed from a real sequence. The fields are empty now, and the
    three views refuse to compute until both are given.
    """
    source = (ROOT / "sections" / "bench.py").read_text(encoding="utf-8")
    assert "28000" not in source
    assert "14000.0" not in source
    # Each of the three views that needs them has to stop, not assume.
    assert source.count("from a sequence, or type it") == 5
    assert source.lower().count("nothing is assumed for you") == 3


def test_the_concentration_is_labelled_derived_not_measured():
    """mg/mL and uM rested on a predicted extinction coefficient, unlabelled.

    The A280 is the measurement; the concentration is Beer-Lambert run on it
    with a ProtParam prediction, so the metric labels and the caption have to
    say which of the two the user is looking at.
    """
    source = (ROOT / "sections" / "bench.py").read_text(encoding="utf-8")
    assert '"mg/mL (derived)"' in source
    assert '"Molar (derived)"' in source
    assert "**Derived, not measured.**" in source
    assert "predicted from the sequence by " in source


def test_the_cystine_gap_is_125_per_disulfide_not_about_one_percent():
    """The help text said the two extinction coefficients differ by about 1%.

    True for a VHH with the single canonical disulfide, wrong by a factor of
    six for a cysteine-rich protein: the gap is a fixed 125 /M/cm per pair, so
    the relative gap falls as the rest of the coefficient grows. Measured over
    the shipped library it runs from 0.0% to 6.6%.
    """
    for seq, pairs in [("CCWY", 1), ("CCCCWY", 2), ("WY", 0)]:
        props = vc.characterize("x", seq)
        gap = props["ext_coeff_cystines"] - props["ext_coeff_reduced"]
        assert gap == 125 * pairs, (seq, gap)

    table, _warnings = lib.load_library()
    entries = [e for e in table.to_dict("records") if e["protein"]]
    gaps = []
    for entry in entries:
        props = vc.characterize(entry["name"], entry["protein"])
        gaps.append(100.0 * (props["ext_coeff_cystines"]
                             - props["ext_coeff_reduced"])
                    / props["ext_coeff_reduced"])
    assert max(gaps) > 6.0, max(gaps)

    source = (ROOT / "sections" / "bench.py").read_text(encoding="utf-8")
    assert "differ by about 1%" not in source
    assert "adds exactly 125 /M/cm" in source


def test_the_chou_fasman_lysine_matches_its_own_citation():
    """The lysine helix propensity read 114 against a citation giving 1.16.

    Expasy ProtScale, the MATLAB Bioinformatics Toolbox and the Chou-Fasman
    teaching tables all publish 1.16 citing Adv. Enzymol. 47:45-148 (1978);
    only one web table gives 114. Every other value in the pair table agrees
    across all of them.
    """
    assert vc.CHOU_FASMAN["K"] == (116, 74)
    assert len(vc.CHOU_FASMAN) == 20
    assert all(isinstance(v, tuple) and len(v) == 2
               for v in vc.CHOU_FASMAN.values())


def test_the_analysis_tab_count_in_the_manual_matches_the_code():
    """The manual said fifteen views; the code built fourteen.

    The families are the single source of truth for how many there are, so
    the manual is checked against them rather than against a second list.

    The count itself is deliberately not written here either. It was, as
    `assert view_count() == 14`, which made this test a third place holding
    the number: adding a view meant editing the code, the manual and this
    line, and the test that existed to catch a stale number needed updating
    by hand every time the number legitimately changed. Now it reads the
    count, spells it, and checks the manual says that and no other.
    """
    from sections.analysis import view_count
    manual = (ROOT / "MANUAL.md").read_text(encoding="utf-8")
    words = {12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen",
             16: "sixteen", 17: "seventeen", 18: "eighteen"}
    count = view_count()
    assert count in words, \
        f"{count} views, and this test cannot spell that; add the word"
    assert f"{words[count]} different views" in manual, \
        f"the analysis builds {count} views and MANUAL.md says otherwise"
    for number, word in words.items():
        if number != count:
            assert f"{word} different views" not in manual


def test_no_constant_claims_a_source_the_repository_cannot_show():
    """Three constants cited sources that could not be confirmed.

    They are kept, because the values are standard and the arithmetic around
    them is shown both ways, but the citation now says what was and was not
    checked. This test is here so the hedge is not quietly deleted later.
    """
    bench_src = (ROOT / "bench.py").read_text(encoding="utf-8")
    # The name survives once, inside the note saying it was never confirmed.
    assert bench_src.count("Day's relation") == 1
    assert 'carried here as "Day\'s relation"' in bench_src
    assert "has not been confirmed against any primary source" in bench_src

    ref_src = (ROOT / "reference.py").read_text(encoding="utf-8")
    assert "that edition has not been consulted" in ref_src

    manual = (ROOT / "MANUAL.md").read_text(encoding="utf-8")
    assert manual.count("Day's relation") == 1
    assert 'called it "Day\'s relation"' in manual
    assert "the standard filamentous-phage relation" in manual
