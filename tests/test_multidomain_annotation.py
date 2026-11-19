"""Annotating a chain that carries more than one variable domain.

These are regression tests for one bug with three faces. `annotate` used to
return the first domain's regions for any sequence, so a bivalent construct was
annotated as though its second VHH did not exist. Nothing failed, nothing was
logged, and three different views were wrong in three different ways: the
liability map and the aggregation profile shaded the first domain's CDRs only,
liabilities in the second domain were filed under region "?" and never counted
as being in a CDR, and the clinical-context table compared the first domain's
CDR lengths to the guidelines under the whole construct's name.

The fix is that an annotation carries its domains, and that the two functions
every consumer goes through, `region_of` and `region_bands`, read them. So the
tests below are mostly about those two and about what the consumers get.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import tnp
import vhh_core as vc
from conftest import HAS_HMMER

pytestmark = pytest.mark.skipif(not HAS_HMMER, reason="needs ANARCI")

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "sequences.fasta"
LINKER = "GGGGSGGGGSGGGGS"


def _examples() -> list[str]:
    out, buf = [], []
    for line in EXAMPLES.read_text().splitlines():
        if line.startswith(">"):
            if buf:
                out.append("".join(buf))
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if buf:
        out.append("".join(buf))
    return [vc.clean_sequence(s) for s in out]


@pytest.fixture(scope="module")
def pair():
    seqs = _examples()
    return seqs[0], seqs[1]


@pytest.fixture(scope="module")
def construct(pair):
    a, b = pair
    return a + LINKER + b + "HHHHHH"


@pytest.fixture(scope="module")
def annotation(construct):
    ann, _method = vc.annotate(construct)
    return ann


def test_a_single_domain_annotation_is_unchanged(pair):
    """The common case must not have moved: a VHH still gets flat regions and
    no `domains` key, which is what every existing caller reads."""
    a, _b = pair
    ann, method = vc.annotate(a)
    assert "domains" not in ann
    assert method == "ANARCI / IMGT"
    assert ann["CDR3"]["seq"]
    assert vc.region_bands(ann) == [
        ("CDR1", ann["CDR1"]["start"], ann["CDR1"]["end"]),
        ("CDR2", ann["CDR2"]["start"], ann["CDR2"]["end"]),
        ("CDR3", ann["CDR3"]["start"], ann["CDR3"]["end"])]


def test_both_domains_are_annotated(annotation, pair):
    a, b = pair
    assert len(annotation["domains"]) == 2
    assert annotation["domains"][0]["CDR3"]["seq"]
    assert annotation["domains"][1]["CDR3"]["seq"]
    # and they are the two different CDR3s, not the same one twice
    assert (annotation["domains"][0]["CDR3"]["seq"]
            != annotation["domains"][1]["CDR3"]["seq"])


def test_domain_coordinates_are_the_chains_not_the_domains(annotation, construct):
    """A band is drawn on a figure of the whole construct, so its coordinates
    have to index the construct. Checked by cutting the sequence at them."""
    for domain in annotation["domains"]:
        for region in vc.REGIONS:
            r = domain[region]
            if not r["start"]:
                continue
            assert construct[r["start"] - 1:r["end"]] == r["seq"], \
                f"{region} does not sit where its coordinates say"


def test_the_flat_regions_are_empty_rather_than_the_first_domains(annotation, pair):
    """The heart of the bug. A caller that reads `annotation["CDR3"]` without
    knowing about domains must get a visible absence, not the first domain's
    CDR3 silently standing in for a construct that has two."""
    a, _b = pair
    first = vc.annotate(a)[0]["CDR3"]["seq"]
    assert annotation["CDR3"]["seq"] == ""
    assert annotation["CDR3"]["start"] is None
    assert first != ""      # the value that used to leak through


def test_bands_cover_every_cdr_of_every_domain(annotation):
    bands = vc.region_bands(annotation)
    assert len(bands) == 6
    assert [b[0] for b in bands] == ["CDR1", "CDR2", "CDR3"] * 2
    # N to C, no overlap between the two domains
    starts = [b[1] for b in bands]
    assert starts == sorted(starts)


def test_region_of_finds_the_second_domains_cdrs(annotation):
    d1, d2 = annotation["domains"]
    assert vc.region_of(d2["CDR3"]["start"], annotation) == "CDR3"
    assert vc.region_of(d2["CDR1"]["start"], annotation) == "CDR1"
    assert vc.region_of(d2["FR2"]["start"], annotation) == "FR2"
    # the linker belongs to no domain and must stay unassigned
    assert vc.region_of(d1["FR4"]["end"] + 2, annotation) == "?"


def test_liabilities_in_the_second_domain_know_where_they_are(annotation,
                                                              construct):
    """`find_liabilities` asks `region_of`, so this is the user-visible
    consequence: a deamidation site in the second VHH's CDR3 used to be filed
    as region '?' with in_CDR False."""
    liabilities = vc.find_liabilities(construct, annotation)
    d2 = annotation["domains"][1]
    in_d2 = liabilities[(liabilities["position"] >= d2["FR1"]["start"])
                        & (liabilities["position"] <= d2["FR4"]["end"])]
    assert not in_d2.empty, "the fixture carries no liability in the second domain"
    assert (in_d2["region"] != "?").all()
    assert in_d2["in_CDR"].any()


def test_the_figures_draw_six_bands(annotation, construct):
    """The figures take their bands from `region_bands`, so drawing them is
    what checks the wiring rather than the enumeration."""
    profile = vc.aggregation_profile(construct, window=9)
    fig = vc.plot_aggregation(profile, annotation)
    shaded = [p for p in fig.axes[0].patches if p.get_alpha() == 0.22]
    assert len(shaded) == 6
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_cdr3_length_is_absent_rather_than_the_first_domains(construct, pair):
    """The Properties table used to show the construct a cdr3_len, which was
    its first domain's. There is no single CDR3, so the cell is empty."""
    import math
    frame, _anns, _errs = vc.analyze({"bip": construct})
    assert math.isnan(float(frame["cdr3_len"].iloc[0]))
    a, _b = pair
    frame_mono, _a2, _e2 = vc.analyze({"mono": a})
    assert float(frame_mono["cdr3_len"].iloc[0]) > 0


def test_tnp_returns_no_lengths_for_a_construct(construct, annotation):
    """Zero is a number and means "measured as nothing"; None means there was
    nothing to measure. The clinical-context thresholds are drawn from
    single-domain nanobodies and a construct is not one."""
    row = tnp.profile("bip", construct, annotation, "ANARCI / IMGT")
    assert row["multi_domain"] is True
    assert row["total_cdr_length"] is None
    assert row["cdr3_length"] is None
    assert row["n_flagged_metrics"] == 0


def test_one_anarci_pass_serves_both_the_domains_and_the_annotation(construct):
    """`find_domains` and `annotate` read the same cached numbering, so
    dissecting a construct and annotating it costs one ANARCI call, not two.
    Checked on the cache rather than by timing."""
    vc._anarci_all.cache_clear()
    vc.find_domains(construct)
    vc.annotate(construct)
    info = vc._anarci_all.cache_info()
    assert info.misses == 1, f"ANARCI ran {info.misses} times"
    assert info.hits >= 1
