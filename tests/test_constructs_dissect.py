"""Cutting a multi-domain construct back into its pieces.

`constructs.dissect_construct` is the inverse of `build_construct`, and the
cheapest way to test an inverse is to run the pair back to back: assemble a
construct from known pieces, dissect the result, and check that the pieces that
come out are the ones that went in. That is what most of this file does, which
means a test here fails if either direction drifts.

The regex path is tested separately and against a weaker claim. It places the
N-terminal boundary of a domain by IMGT convention rather than by measurement,
so it is checked for finding the right *number* of domains in the right order,
not for residue-exact bounds; asserting bounds it does not promise would be a
test of the convention rather than of the code.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import constructs as ctm
import vhh_core as vc
from conftest import HAS_HMMER

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "sequences.fasta"


def _examples() -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in EXAMPLES.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


@pytest.fixture(scope="module")
def two_vhhs() -> tuple[str, str]:
    seqs = _examples()
    assert len(seqs) >= 2, "the example FASTA needs at least two sequences"
    return vc.clean_sequence(seqs[0][1]), vc.clean_sequence(seqs[1][1])


# --------------------------------------------------------------------------- #
# find_domains
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_single_domain_is_found_whole(two_vhhs):
    a, _b = two_vhhs
    domains = vc.find_domains(a)
    assert len(domains) == 1
    assert domains[0]["start"] == 0
    assert domains[0]["end"] == len(a) - 1
    assert domains[0]["bounds"] == "exact"


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
@pytest.mark.parametrize("linker", ["", "GS", "GGGGS", "GGGGSGGGGSGGGGS",
                                    "EAAAKEAAAKEAAAK"])
def test_both_domains_found_and_linker_recovered(two_vhhs, linker):
    """The linker is what is left between two domains, so it is recovered by
    subtraction. This is the measurement the whole feature rests on: if ANARCI
    stopped reporting the second domain, every construct would silently become
    a monovalent VHH again."""
    a, b = two_vhhs
    domains = vc.find_domains(a + linker + b)
    assert len(domains) == 2
    gap = (a + linker + b)[domains[0]["end"] + 1:domains[1]["start"]]
    assert gap == linker


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_non_antibody_has_no_domains():
    assert vc.find_domains("M" + "AGTSEPKLQDRVN" * 8) == []


def test_regex_scanner_finds_both_domains(two_vhhs):
    """Without ANARCI the scanner still separates two domains.

    Bounds are not asserted: `_domains_regex` places the N-terminal end of a
    domain at the conserved cysteine minus 22 residues, which is a convention,
    and the first residues of a domain can land in the linker because of it.
    What must hold is that two domains are seen, in order, and that neither
    swallows the other."""
    a, b = two_vhhs
    domains = vc._domains_regex(a + "GGGGSGGGGSGGGGS" + b)
    assert len(domains) == 2
    assert all(d["bounds"] == "approximate" for d in domains)
    assert domains[0]["end"] < domains[1]["start"]


# --------------------------------------------------------------------------- #
# classify
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_two_domains_classify_as_a_fusion(two_vhhs):
    a, b = two_vhhs
    kind = vc.classify(a + "GGGGSGGGGSGGGGS" + b)
    assert kind["kind"] == vc.KIND_FUSION
    assert kind["n_domains"] == 2
    assert not kind["full_analysis"]


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_a_tagged_monovalent_is_still_a_vhh(two_vhhs):
    """A His6 on the end must not turn a VHH into a fusion: the coverage rule
    has always allowed a short tag, and the domain count must not change
    that."""
    a, _b = two_vhhs
    kind = vc.classify(a + "HHHHHH")
    assert kind["kind"] == vc.KIND_VHH
    assert kind["n_domains"] == 1
    assert kind["full_analysis"]


# --------------------------------------------------------------------------- #
# dissect_construct, against build_construct
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_round_trip_through_build_and_dissect(two_vhhs):
    a, b = two_vhhs
    built = ctm.build_construct(
        [("first", a), ("second", b)],
        linker=ctm.LINKERS["(GGGGS)3"]["seq"],
        c_tags=[("His6", ctm.TAGS["His6"]["seq"])])
    cut = ctm.dissect_construct(built["sequence"], name="X")

    assert cut["sequence"] == built["sequence"]
    assert cut["n_domains"] == 2
    assert [d["sequence"] for d in cut["domains"]] == [a, b]

    kinds = list(cut["segments"]["kind"])
    assert kinds == ["vhh", "linker", "vhh", "tag"]
    labels = list(cut["segments"]["label"])
    assert labels[1] == "(GGGGS)3"
    assert labels[3] == "His6"


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_segments_tile_the_sequence_without_gaps(two_vhhs):
    """Every residue belongs to exactly one segment. A gap here would mean the
    schematic and the composition quietly disagree with the sequence."""
    a, b = two_vhhs
    cut = ctm.dissect_construct("YPYDVPDYA" + a + "GGGGSGGGGS" + b + "HHHHHH",
                                name="X")
    seg = cut["segments"]
    assert seg["start"].iloc[0] == 1
    assert seg["end"].iloc[-1] == cut["length"]
    for prev, nxt in zip(seg.itertuples(), seg.iloc[1:].itertuples()):
        assert nxt.start == prev.end + 1
    assert sum(cut["composition"].values()) == cut["length"]


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_three_domains(two_vhhs):
    a, b = two_vhhs
    cut = ctm.dissect_construct(a + "GGGGSGGGGS" + b + "GGGGSGGGGS" + a,
                                name="tri")
    assert cut["n_domains"] == 3
    assert [d["name"] for d in cut["domains"]] == ["tri_D1", "tri_D2", "tri_D3"]
    assert list(cut["per_vhh"]["name"]) == ["tri_D1", "tri_D2", "tri_D3"]


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_unknown_spacer_is_not_called_a_linker(two_vhhs):
    """An unrecognised stretch is reported as unrecognised. Naming it after the
    nearest known linker would be a guess presented as an identification."""
    a, b = two_vhhs
    cut = ctm.dissect_construct(a + "PQWTVSDKAQ" + b, name="X")
    middle = cut["segments"].iloc[1]
    assert middle["label"] == "unassigned"
    assert middle["sequence"] == "PQWTVSDKAQ"
    assert any("unassigned" in w for w in cut["warnings"])


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_whole_chain_mass_is_the_sum_of_its_parts(two_vhhs):
    """Mass is additive across a fusion minus the water lost at each bond, so
    the whole chain's mass must sit just under the sum of the pieces. This is
    the claim the Bench screen relies on when it turns an A280 into a
    concentration for a biparatopic."""
    a, b = two_vhhs
    linker = "GGGGSGGGGSGGGGS"
    cut = ctm.dissect_construct(a + linker + b, name="X")
    whole = float(cut["properties"]["mw_da"].iloc[0])
    pieces = sum(float(vc.characterize("p", s)["mw_da"])
                 for s in (a, linker, b))
    water = 18.015 * 2
    assert pieces - whole == pytest.approx(water, abs=0.5)


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_junction_sites_are_outside_the_domains(two_vhhs):
    """(EAAAK)3 carries a lysine per repeat, so a trypsin site in the linker is
    a real finding and not a quirk of the fixture."""
    a, b = two_vhhs
    cut = ctm.dissect_construct(a + "EAAAKEAAAKEAAAK" + b, name="X")
    js = cut["junction_sites"]
    assert js is not None and not js.empty
    domain_spans = [(d["start"], d["end"]) for d in cut["domains"]]
    for pos in js["position"]:
        assert not any(s <= pos <= e for s, e in domain_spans)
    assert "Trypsin-like" in set(js["protease"])


@pytest.mark.skipif(not HAS_HMMER, reason="ANARCI needs HMMER")
def test_dissected_construct_works_with_the_builder_helpers(two_vhhs):
    """The dissection returns the shape `build_construct` returns, which is
    what lets `plot_construct` and `construct_dna` take either one."""
    a, b = two_vhhs
    cut = ctm.dissect_construct(a + "GGGGSGGGGSGGGGS" + b, name="X")
    dna = ctm.construct_dna(cut, add_stop=True)
    assert dna["back_translation_ok"]
    fig = ctm.plot_construct(cut, "X")
    assert fig is not None


def test_no_domain_at_all_is_reported_not_guessed():
    cut = ctm.dissect_construct("AGTSEPKLQDRVN" * 8, name="X")
    assert cut["n_domains"] == 0
    assert cut["per_vhh"].empty
    assert any("nothing" in w for w in cut["warnings"])
