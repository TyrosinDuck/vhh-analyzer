"""
Metadata tests — the three lines nobody ever looks at.

`CITATION.cff` and `LICENSE` carry the author's name, his ORCID, the concept
DOI and the copyright holder. None of it affects a single thing the app does, so
nothing else in this suite or in the app would notice if it went missing: the
page still renders, every analysis still runs, and GitHub's **Cite this
repository** button quietly offers a citation with no DOI in it. You would find
out months later, from a colleague's reference list.

They are also the lines most likely to be lost by accident, because they are
edited on GitHub rather than in the working copy: a DOI does not exist until the
release is published, so it is typically added to the file in the browser, and
the next copy of the repository overwrites it.

Hence this file. It is not testing behaviour; it is testing that four specific
facts survive, and it exists because they are invisible when they do not.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml", reason="PyYAML is needed to read CITATION.cff")

ROOT = Path(__file__).resolve().parent.parent

AUTHOR_FAMILY_NAME = "Vaneyck"
ORCID = "https://orcid.org/0000-0003-3095-7420"
# The Zenodo *concept* DOI: it resolves to the latest release. The per-release
# version DOIs are not pinned here, since there is a new one every release.
CONCEPT_DOI = "10.5281/zenodo.23742269"


@pytest.fixture(scope="module")
def citation():
    path = ROOT / "CITATION.cff"
    assert path.exists(), "CITATION.cff is missing — the repository is not citable"
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_citation_file_is_valid_cff(citation):
    """Malformed CFF is worse than none: GitHub shows an error instead of a button."""
    assert citation.get("cff-version", "").startswith("1."), \
        "cff-version is missing or not a 1.x version"
    for field in ("title", "message", "type", "authors", "license"):
        assert citation.get(field), f"CITATION.cff has no {field}"
    assert citation["type"] == "software"


def test_the_concept_doi_is_present_and_well_formed(citation):
    """The line most likely to be lost, and the whole point of archiving.

    It must be the bare identifier, not a URL: CFF wants `10.5281/zenodo.N`, and
    a `https://doi.org/...` here produces a citation with a broken DOI field.
    """
    doi = citation.get("doi")
    assert doi, (
        "CITATION.cff has no doi. If a copy of the repository overwrote it, the "
        f"concept DOI is {CONCEPT_DOI} — it resolves to the latest release.")
    assert not doi.startswith("http"), \
        f"doi must be the bare identifier, not a URL: got {doi!r}"
    assert re.fullmatch(r"10\.5281/zenodo\.\d+", doi), \
        f"{doi!r} does not look like a Zenodo DOI"
    assert doi == CONCEPT_DOI, (
        f"the doi is {doi!r} but the concept DOI is {CONCEPT_DOI!r}. If this is a "
        f"deliberate change, update CONCEPT_DOI in this test; if a per-release "
        f"version DOI has been pasted in by mistake, the concept DOI is the one "
        f"Zenodo describes as representing all versions.")


def test_the_author_is_attributed_correctly(citation):
    """One author, spelled the way he spells it, with his ORCID.

    The ORCID is the part that matters beyond this repository: it is what keeps
    the work attached to the person across institutions, which is the entire
    reason for filling in citation metadata rather than relying on a URL.
    """
    authors = citation["authors"]
    assert len(authors) == 1, \
        f"expected one author, found {len(authors)} — has a contributor been added?"
    author = authors[0]
    assert author.get("family-names") == AUTHOR_FAMILY_NAME, \
        f"family name is {author.get('family-names')!r}, expected {AUTHOR_FAMILY_NAME!r}"
    assert author.get("orcid") == ORCID, \
        f"ORCID is {author.get('orcid')!r}, expected {ORCID!r}"


def test_the_licence_exists_and_names_the_right_holder():
    """MIT with the wrong name in the copyright line grants nothing cleanly."""
    path = ROOT / "LICENSE"
    assert path.exists(), "LICENSE is missing — without it the default is all rights reserved"
    text = path.read_text(encoding="utf-8")
    assert "MIT License" in text
    assert AUTHOR_FAMILY_NAME in text, \
        f"the copyright line does not name {AUTHOR_FAMILY_NAME}"
    assert "Permission is hereby granted" in text, \
        "the MIT grant itself is missing from LICENSE"


def test_the_licence_states_the_gpl_dependency():
    """pyfamsa is GPL-3.0-or-later and the only copyleft dependency.

    An MIT licence that says nothing about it implies the whole application can
    be reused under MIT terms, which is not true of a distribution that includes
    FAMSA. Claiming more than you can grant is its own kind of wrong.
    """
    text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "pyfamsa" in text and "GPL-3.0" in text, \
        "LICENSE no longer mentions the GPL-3.0 dependency on pyfamsa"


def test_the_licence_is_declared_in_the_citation_file(citation):
    assert citation["license"] == "MIT", \
        f"CITATION.cff declares {citation['license']!r}, LICENSE is MIT"
