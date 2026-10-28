"""
Metadata tests, the three lines nobody ever looks at.

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
    assert path.exists(), "CITATION.cff is missing, the repository is not citable"
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
        f"concept DOI is {CONCEPT_DOI}, it resolves to the latest release.")
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
        f"expected one author, found {len(authors)}, has a contributor been added?"
    author = authors[0]
    assert author.get("family-names") == AUTHOR_FAMILY_NAME, \
        f"family name is {author.get('family-names')!r}, expected {AUTHOR_FAMILY_NAME!r}"
    assert author.get("orcid") == ORCID, \
        f"ORCID is {author.get('orcid')!r}, expected {ORCID!r}"


def test_the_licence_exists_and_names_the_right_holder():
    """MIT with the wrong name in the copyright line grants nothing cleanly."""
    path = ROOT / "LICENSE"
    assert path.exists(), "LICENSE is missing, without it the default is all rights reserved"
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


def test_the_community_guidelines_exist():
    """JOSS requires guidance for contributing, reporting issues and support.

    Checked here rather than left to memory because the file is documentation
    nobody runs: it can be deleted in a reorganisation without anything failing.
    """
    path = ROOT / "CONTRIBUTING.md"
    assert path.exists(), "CONTRIBUTING.md is missing"
    text = path.read_text(encoding="utf-8").lower()
    for topic, needle in (("reporting bugs", "reporting a bug"),
                          ("asking for help", "asking for help"),
                          ("contributing code", "contributing code")):
        assert needle in text, f"CONTRIBUTING.md says nothing about {topic}"
    assert "public" in text, \
        "CONTRIBUTING.md no longer warns that issues and the library are public"


def test_the_app_tells_users_how_to_report_a_problem():
    """CONTRIBUTING.md lives where only developers look.

    The people who hit defects are using the app, not reading the repository, so
    the route to the issue tracker has to exist in the interface - and the
    pre-filled report has to carry the version, which is the one thing nobody
    supplies unprompted and the one thing that makes a report actionable.
    """
    import shell

    assert shell.ISSUES_URL.endswith("/issues")
    url = shell.report_url("Some screen")
    assert "/issues/new?" in url, "report_url no longer opens a new issue"
    from urllib.parse import unquote

    body = unquote(url)
    assert shell.VERSION in body, \
        "the pre-filled report no longer carries the version"
    assert "public" in body.lower(), \
        "the pre-filled report no longer warns that the issue is public"


def test_the_licence_is_declared_in_the_citation_file(citation):
    assert citation["license"] == "MIT", \
        f"CITATION.cff declares {citation['license']!r}, LICENSE is MIT"


def test_the_app_and_the_citation_file_agree(citation):
    """`shell.py` shows the DOI in the footer; it must be the same one.

    Two copies of an identifier in one repository is one copy too many, and the
    failure is silent: the app would advertise a DOI that cites something else,
    on every page, and nothing would look wrong.
    """
    import shell

    assert shell.DOI == citation["doi"], (
        f"shell.DOI is {shell.DOI!r} but CITATION.cff says {citation['doi']!r}")
    assert shell.AUTHOR_URL == citation["authors"][0]["orcid"], \
        "the ORCID in the footer and in CITATION.cff differ"
    assert citation["authors"][0]["family-names"] in shell.AUTHOR, \
        f"the footer credits {shell.AUTHOR!r}, which does not contain the " \
        f"family name in CITATION.cff"
    assert shell.LICENCE == citation["license"]
    assert shell.RELEASE == citation.get("version"), (
        f"the app reports release {shell.RELEASE!r} but CITATION.cff says "
        f"{citation.get('version')!r}, the footer, the git tag and the Zenodo "
        f"record are supposed to name the same thing")


def test_the_release_number_is_semantic_versioning():
    """A build datestamp in the footer names a version nobody else knows.

    The footer used to read `version 20261028.2600`, which matched no git tag and
    no Zenodo record, so a bug report or a citation referred to something
    unfindable. RELEASE is the tag; BUILD distinguishes two builds of it.
    """
    import shell

    assert re.fullmatch(r"\d+\.\d+\.\d+", shell.RELEASE), \
        f"RELEASE is {shell.RELEASE!r}, expected MAJOR.MINOR.PATCH"
    assert re.fullmatch(r"\d{8}", shell.BUILD), \
        f"BUILD is {shell.BUILD!r}, expected YYYYMMDD"
    assert re.fullmatch(r"[0-9a-f]{7}", shell.CODE), \
        f"CODE is {shell.CODE!r}, expected seven hex characters"
    for part in (shell.RELEASE, shell.BUILD, shell.CODE):
        assert part in shell.VERSION, f"{part!r} is missing from the version string"


def test_the_fingerprint_changes_when_the_source_changes(tmp_path):
    """The point of it: telling a finished deploy from a stale container.

    The build date cannot do this. Two pushes on the same day carry the same
    date, so a footer that reads `build 20261028` after a deploy says nothing
    about whether the new code is live.
    """
    import shell

    (tmp_path / "sections").mkdir()
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "sections" / "b.py").write_text("y = 2\n")
    first = shell._source_fingerprint(tmp_path)

    assert re.fullmatch(r"[0-9a-f]{7}", first)
    assert shell._source_fingerprint(tmp_path) == first, \
        "the same sources must give the same fingerprint"

    (tmp_path / "sections" / "b.py").write_text("y = 3\n")
    assert shell._source_fingerprint(tmp_path) != first, \
        "a changed source file left the fingerprint unchanged"

    assert shell._source_fingerprint(tmp_path / "nothing-here") == "unknown"


def test_the_changelog_documents_the_current_release():
    """A changelog that stops before the current version is a stale changelog."""
    import shell

    path = ROOT / "CHANGELOG.md"
    assert path.exists(), "CHANGELOG.md is missing"
    text = path.read_text(encoding="utf-8")
    assert f"[{shell.RELEASE}]" in text, (
        f"CHANGELOG.md has no entry for {shell.RELEASE}, bump the changelog "
        f"when you bump the release")
