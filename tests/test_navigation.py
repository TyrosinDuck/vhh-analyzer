"""
Navigation tests, the grouping of the sections.

The home page, the sidebar and the flow diagram all group the sections into
Protein, DNA and Bench, and all three tint a group with the same colour. None
of them stores that grouping: each derives it from `shell.SECTIONS`, where
every section names its own group, through `shell.sections_by_group()`. What
these tests guard is that derivation, because its failure mode is silent: a
section that belongs to no group is drawn in no group, which on the home page
and in the sidebar means it is not drawn at all.
"""

from __future__ import annotations

import re

import pytest

import shell


def test_every_section_belongs_to_exactly_one_known_group():
    grouped = shell.sections_by_group()
    flat = [key for keys in grouped.values() for key in keys]
    assert sorted(flat) == sorted(shell.SECTIONS), \
        f"sections missing or duplicated in the grouping: {grouped}"
    assert len(flat) == len(set(flat))
    for key, entry in shell.SECTIONS.items():
        assert len(entry) == 4, f"{key}: expected (icon, label, desc, group)"
        assert entry[3] in shell.GROUPS, f"{key}: unknown group {entry[3]!r}"


def test_groups_come_in_the_declared_order_and_none_is_empty():
    grouped = shell.sections_by_group()
    assert list(grouped) == [g for g in shell.GROUPS if g in grouped]
    assert all(grouped.values())


def test_the_agreed_grouping():
    """The grouping itself, as decided: by the molecule a section works on.

    Pinned because it is a design decision rather than something derivable,
    so a change to it should be deliberate and show up as a failing test.
    """
    assert shell.sections_by_group() == {
        "Protein": ["analysis", "structure", "library"],
        "DNA": ["sequencing", "genes", "cloning"],
        "Bench": ["bench"],
    }


def test_a_misspelt_group_raises_instead_of_hiding_the_section(monkeypatch):
    sections = dict(shell.SECTIONS)
    icon, label, desc, _group = sections["bench"]
    sections["bench"] = (icon, label, desc, "Benchh")
    monkeypatch.setattr(shell, "SECTIONS", sections)
    with pytest.raises(KeyError, match="Benchh"):
        shell.sections_by_group()


def test_group_colours_are_hex_and_distinct():
    colours = [colour for colour, _tagline in shell.GROUPS.values()]
    assert all(re.fullmatch(r"#[0-9a-fA-F]{6}", c) for c in colours), colours
    assert len(set(c.lower() for c in colours)) == len(colours)


def test_the_group_header_writes_the_name_and_not_only_the_colour():
    for name, (colour, tagline) in shell.GROUPS.items():
        html = shell.group_header_html(name, "1rem")
        assert f">{name}<" in html
        assert tagline in html
        assert colour in html
        assert "font-size:1rem" in html


def test_the_flow_diagram_tints_every_section_by_its_group():
    home = pytest.importorskip("sections.home")
    svg = home.flow_diagram()
    assert "<!--" + "BOXES-->" not in svg and "<!--" + "LEGEND-->" not in svg
    for key, (_icon, _label, _desc, group) in shell.SECTIONS.items():
        x, y, *_ = home._BOXES[key]
        rect = re.search(rf'<rect x="{x}" y="{y}"[^>]*>', svg)
        assert rect, f"{key}: no box in the diagram"
        assert f'fill="{shell.GROUPS[group][0]}"' in rect.group(0), \
            f"{key}: box not tinted with the {group} colour"
    for name in shell.sections_by_group():
        assert f'font-weight="600">{name}</text>' in svg, \
            f"{name} missing from the diagram legend"


def test_every_section_has_a_box_in_the_flow_diagram():
    home = pytest.importorskip("sections.home")
    assert set(home._BOXES) == set(shell.SECTIONS)


def test_no_section_has_the_name_of_its_group():
    """A section called like its group read as a duplicate on the home page,
    the group title "Bench" right above a card titled "Bench", which is why
    that section became "Batch & QC"."""
    for key, (_icon, label, _desc, group) in shell.SECTIONS.items():
        assert label.casefold() != group.casefold(), \
            f"{key}: the section {label!r} has the name of its group"


def test_the_flow_diagram_titles_come_from_sections_and_are_escaped():
    import html

    home = pytest.importorskip("sections.home")
    svg = home.flow_diagram()
    for key, (icon, label, *_rest) in shell.SECTIONS.items():
        assert f">{html.escape(icon + ' ' + label)}</text>" in svg, \
            f"{key}: title missing from the diagram"
    # A bare "&" is not well-formed XML; inside an SVG it can break the
    # drawing. Every ampersand must start an entity.
    assert not re.search(r"&(?![a-zA-Z]+;|#\d+;)", svg)


@pytest.mark.parametrize("a, b, expected", [
    ("#000000", "#ffffff", 21.0),     # the maximum, by definition
    ("#ffffff", "#ffffff", 1.0),      # the minimum
    ("#777777", "#ffffff", 4.48),     # the classic just-fails-AA grey
])
def test_contrast_ratio_matches_known_values(a, b, expected):
    assert shell.contrast_ratio(a, b) == pytest.approx(expected, abs=0.01)
    assert shell.contrast_ratio(b, a) == pytest.approx(expected, abs=0.01)


def test_every_open_button_text_is_readable_on_its_group_colour():
    """WCAG 2.1 SC 1.4.3, AA for normal text: at least 4.5:1, for the white
    text the Open buttons carry on every group colour."""
    home = pytest.importorskip("sections.home")
    for name, (colour, _tagline) in shell.GROUPS.items():
        ratio = shell.contrast_ratio(colour, home.BUTTON_TEXT)
        assert ratio >= 4.5, \
            f"{name}: {home.BUTTON_TEXT} on {colour} is only {ratio:.2f}"


def test_the_open_buttons_are_painted_with_their_group_colour():
    home = pytest.importorskip("sections.home")
    css = home._open_button_css()
    for key, (*_rest, group) in shell.SECTIONS.items():
        fill = shell.GROUPS[group][0]
        rule = re.search(rf"\.st-key-home_{key} button\{{([^}}]*)\}}", css)
        assert rule, f"{key}: no rule for its Open button"
        assert f"background-color:{fill}" in rule.group(1)
        assert f"color:{home.BUTTON_TEXT}" in rule.group(1)
