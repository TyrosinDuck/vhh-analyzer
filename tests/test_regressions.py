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
