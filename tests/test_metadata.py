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
    (tmp_path / "data").mkdir()
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "sections" / "b.py").write_text("y = 2\n")
    (tmp_path / "MANUAL.md").write_text("# manual\n")
    (tmp_path / "data" / "library.csv").write_text("name\nx\n")
    first = shell._source_fingerprint(tmp_path)

    assert re.fullmatch(r"[0-9a-f]{7}", first)
    assert shell._source_fingerprint(tmp_path) == first, \
        "the same sources must give the same fingerprint"

    (tmp_path / "sections" / "b.py").write_text("y = 3\n")
    second = shell._source_fingerprint(tmp_path)
    assert second != first, \
        "a changed source file left the fingerprint unchanged"

    # The manual and the library are served by the app, so a change to either is
    # visible to a user and has to move the fingerprint too.
    (tmp_path / "MANUAL.md").write_text("# manual, revised\n")
    third = shell._source_fingerprint(tmp_path)
    assert third != second, "a changed MANUAL.md left the fingerprint unchanged"

    (tmp_path / "data" / "library.csv").write_text("name\ny\n")
    fourth = shell._source_fingerprint(tmp_path)
    assert fourth != third, "a changed library left the fingerprint unchanged"

    # So are the amino acid structures, which are files rather than code and
    # would otherwise be the one thing a user can see change without the
    # version following.
    (tmp_path / "data" / "aminoacids").mkdir()
    (tmp_path / "data" / "aminoacids" / "Ala.svg").write_text("<svg/>")
    fifth = shell._source_fingerprint(tmp_path)
    assert fifth != fourth, "a new structure left the fingerprint unchanged"

    (tmp_path / "data" / "aminoacids" / "Ala.svg").write_text("<svg />")
    assert shell._source_fingerprint(tmp_path) != fifth, \
        "a changed structure left the fingerprint unchanged"

    # The config file sets the upload cap and is hidden, which makes it the
    # easiest file in the project to leave behind when copying a folder: a
    # deployment without it silently runs with Streamlit's 200 MB default.
    # Hashing it is what makes that visible in the footer.
    (tmp_path / ".streamlit").mkdir()
    (tmp_path / ".streamlit" / "config.toml").write_text("[server]\n")
    sixth = shell._source_fingerprint(tmp_path)
    assert sixth != fifth, "a missing config.toml left the fingerprint unchanged"

    # Community Cloud installs both, so an edit to either changes the
    # deployed app. (Only an edit to the file: the versions are not pinned.)
    (tmp_path / "requirements.txt").write_text("streamlit>=1.30\n")
    with_requirements = shell._source_fingerprint(tmp_path)
    assert with_requirements != sixth, \
        "a new requirements.txt left the fingerprint unchanged"
    (tmp_path / "requirements.txt").write_text("streamlit>=1.49\n")
    assert shell._source_fingerprint(tmp_path) != with_requirements, \
        "a changed requirements.txt left the fingerprint unchanged"
    (tmp_path / "requirements.txt").write_text("streamlit>=1.30\n")
    (tmp_path / "packages.txt").write_text("hmmer\n")
    with_packages = shell._source_fingerprint(tmp_path)
    assert with_packages != with_requirements, \
        "a new packages.txt left the fingerprint unchanged"

    # And what Community Cloud never reads must not move it: a fingerprint
    # that changes for a README typo can no longer say whether a deploy took.
    for name in ("README.md", "CHANGELOG.md", "Dockerfile", ".dockerignore",
                 "requirements-dev.txt"):
        (tmp_path / name).write_text("anything\n")
    assert shell._source_fingerprint(tmp_path) == with_packages, \
        "a file Community Cloud does not use moved the fingerprint"
    sixth = with_packages

    # Called with a string rather than a Path, which is what comparing two
    # checkouts does. It used to answer "unknown", indistinguishable from a
    # directory with nothing in it.
    assert shell._source_fingerprint(str(tmp_path)) == sixth

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


# --------------------------------------------------------------------------
# Deployment configuration
# --------------------------------------------------------------------------

def test_xsrf_protection_stays_off():
    """It broke every upload on the deployed app, and only there.

    `.streamlit/config.toml` had ``enableXsrfProtection = true``. Locally that
    is harmless: one origin, no proxy, the token validates. On Streamlit
    Community Cloud the app sits behind their proxy, the upload goes to a
    different origin than the one that issued the token, and the upload
    endpoint answers 500 before the request reaches any of this code. So
    nothing appeared in the logs, every test here passed, and the only symptom
    was a user dropping a PDB into the structure section and getting
    "AxiosError: Request failed with status code 500".

    Community Cloud runs every app it hosts with this off. Setting it back on
    turns off file uploads on the deployed app, and nowhere else, which is the
    hardest kind of breakage to find.
    """
    config = (ROOT / ".streamlit" / "config.toml").read_text()
    setting = re.search(r"^\s*enableXsrfProtection\s*=\s*(\w+)", config,
                        re.MULTILINE)
    assert setting, "enableXsrfProtection is not set in .streamlit/config.toml"
    assert setting.group(1) == "false", (
        "enableXsrfProtection is on again; this disables every file upload on "
        "Community Cloud while working perfectly on your machine")


def test_the_upload_cap_is_what_the_manual_says():
    """Two numbers in two files that have to agree, and nothing links them."""
    config = (ROOT / ".streamlit" / "config.toml").read_text()
    size = re.search(r"^\s*maxUploadSize\s*=\s*(\d+)", config, re.MULTILINE)
    assert size, "maxUploadSize is not set in .streamlit/config.toml"
    manual = (ROOT / "MANUAL.md").read_text()
    assert f"{size.group(1)} MB per file" in manual, (
        f"the config caps uploads at {size.group(1)} MB and the manual says "
        f"something else")


def test_nothing_calls_the_deprecated_html_component_directly():
    """One place decides how an iframe is drawn, and it is `shell.embed_html`.

    Streamlit replaced ``st.components.v1.html`` with ``st.iframe`` and says
    it will remove the old one. The deployed app installs whatever Streamlit
    is current, a development machine may be behind, and three screens need a
    real iframe: the home diagram, the amino acid structures and the 3D
    viewer. A second call site that picked the API for itself would be a
    screen that breaks on one of the two versions, silently, because it only
    draws when you open that panel.
    """
    offenders = []
    for path in sorted(ROOT.glob("sections/*.py")) + sorted(ROOT.glob("*.py")):
        if path.name == "shell.py":
            continue
        text = path.read_text(encoding="utf-8")
        for number, line in enumerate(text.splitlines(), start=1):
            code = line.split("#", 1)[0]
            if "components.html(" in code or "components.v1" in code:
                offenders.append(f"{path.name}:{number}")
    assert not offenders, (
        "these call the deprecated component directly instead of "
        "shell.embed_html: " + ", ".join(offenders))


# --------------------------------------------------------------------------
# Dependencies: the lockfile and what the footer reports
# --------------------------------------------------------------------------

def _requirement_names(text: str) -> set[str]:
    """Package names in a requirements file, normalised the way pip does."""
    names = set()
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or line.startswith("-"):
            continue
        name = re.match(r"[A-Za-z0-9._-]+", line).group(0)
        names.add(re.sub(r"[-_.]+", "-", name).lower())
    return names


def _locked_python() -> str:
    header = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    found = re.search(r"autogenerated by pip-compile with Python (\d+\.\d+)",
                      header)
    assert found, ("requirements.txt has no pip-compile header; it is "
                   "generated from requirements.in, never written by hand")
    return found.group(1)


def test_requirements_txt_is_a_complete_lockfile():
    """With ranges, a rebuild installed whatever was current that day.

    `requirements.txt` used to give minimum versions, so Community Cloud
    picked the newest Streamlit and Biopython on each rebuild while the
    footer fingerprint, which hashes the file, stayed where it was. The app
    could change under an unchanged fingerprint. It is now generated by
    pip-compile from `requirements.in`, and this checks the two things that
    make it a lock: every line is an exact pin, and every direct dependency
    in `requirements.in` is among them. A hand edit that adds a range, or a
    package added to `requirements.in` without recompiling, fails here.
    """
    _locked_python()
    lock = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    loose = [line for line in lock.splitlines()
             if line.strip() and not line.lstrip().startswith("#")
             and not re.fullmatch(r"[A-Za-z0-9._-]+==[^\s;]+", line.strip())]
    assert not loose, f"requirements.txt has lines that are not pins: {loose}"

    wanted = _requirement_names((ROOT / "requirements.in").read_text())
    pinned = _requirement_names(lock)
    missing = wanted - pinned
    assert not missing, (
        f"{sorted(missing)} are in requirements.in but not in requirements.txt: "
        f"recompile with pip-compile, see README.md, 'Updating dependencies'")


def test_the_python_version_is_the_same_everywhere_it_is_named():
    """The lock only installs on the Python it was resolved for.

    pip-compile writes the version it ran under into the header of
    `requirements.txt`, and that is the one source. The deploy instructions,
    the manual and the Dockerfile name it too, because a person or a build has
    to choose an interpreter; this checks they still name the same one. A lock
    resolved for 3.14 does not install on 3.11 (contourpy 1.4.0 needs 3.12 or
    later, measured when the lock was introduced), so a Dockerfile left on the
    old version would fail at its pip step.
    """
    version = _locked_python()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert f"select Python {version}" in readme, \
        f"README.md's deploy steps do not say to select Python {version}"
    assert f"python={version}" in (ROOT / "Dockerfile").read_text(), \
        f"the Dockerfile does not install python={version}"
    assert f"Python {version}" in (ROOT / "MANUAL.md").read_text(
        encoding="utf-8"), f"MANUAL.md never names Python {version}"


def test_the_footer_versions_are_the_running_ones():
    """The footer is how you learn what Community Cloud actually installed.

    The fingerprint covers the files; it does not cover the interpreter or
    what apt chose for `packages.txt`. `runtime_versions()` reads them from
    the process, and has to agree with what the process has loaded.
    """
    import platform
    from importlib import metadata

    import Bio
    import streamlit
    import shell

    shell.runtime_versions.cache_clear()
    versions = shell.runtime_versions()
    assert list(versions) == ["Python", "Streamlit", "Biopython", "ANARCI"]
    assert versions["Python"] == platform.python_version()
    assert versions["Streamlit"] == streamlit.__version__
    assert versions["Biopython"] == Bio.__version__
    try:
        expected = metadata.version("anarci")
    except metadata.PackageNotFoundError:
        expected = "absent"
    assert versions["ANARCI"] == expected


def test_the_footer_says_absent_when_anarci_is_not_installed(monkeypatch):
    """An install without ANARCI is legitimate, and the footer must not crash."""
    from importlib import metadata

    import shell

    real = metadata.version

    def without_anarci(name):
        if name.lower() == "anarci":
            raise metadata.PackageNotFoundError(name)
        return real(name)

    monkeypatch.setattr(metadata, "version", without_anarci)
    shell.runtime_versions.cache_clear()
    try:
        assert shell.runtime_versions()["ANARCI"] == "absent"
    finally:
        shell.runtime_versions.cache_clear()
