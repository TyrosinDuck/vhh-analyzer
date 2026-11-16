"""The cassette layout is read out of each vector's own map, and nowhere else.

`sanger.VECTORS` used to be two entries written out by hand, pMES3 and pMECS,
beside plasmid maps the app already ships in `data/plasmids/`. Two accounts of
the same thing, and they had drifted: the pMECS entry carried the nine-residue
HA epitope YPYDVPDYA while the pMECS map annotates `HA-tag` over 24 nucleotides,
which encode the eight residues YPYDVPDY, so every correct pMECS clone was
reported as missing its HA tag. That entry was corrected by hand, which left the
two lists separate and free to drift again, and left pMES4, pHEN6 and pHEN25
with no entry at all: the cloning screen checked their constructs against
pMES3's tags.

`cloning.derive_layout` reads the layout from the map, and `sanger.VECTORS` is
now what it returns for the five bundled maps. These tests pin what the five
derive to, so a change in the derivation that silently alters a layout fails
here; they pin the pMECS HA case on its own, because that is the mistake that
already happened; and they pin the two traps the derivation has to avoid, a
label matched as a substring and a tag read in the leader's frame.

What they cannot check: whether a map's own annotation is right. If a vendor's
map labels the wrong span `HA-tag`, the app reads the wrong span. Deriving moves
the trust from a list typed here to the file the user supplies, it does not
remove it.
"""

from __future__ import annotations

import pandas as pd
import pytest

import cloning as cl
import sanger as sg

# What the five bundled maps derive to today. Not a specification of these
# vectors: a record of what the app reads out of the files it ships, so a change
# in the reading is visible in the diff of this test rather than in a result.
EXPECTED = {
    "pMES3": (["PelB", "VHH", "His6", "stop"], [("His6", "HHHHHH")]),
    "pMES4": (["PelB", "VHH", "His6", "stop"], [("His6", "HHHHHH")]),
    "pHEN6": (["PelB", "VHH", "His6", "stop"], [("His6", "HHHHHH")]),
    "pHEN25": (["PelB", "VHH", "His6", "stop"], [("His6", "HHHHHH")]),
    "pMECS": (["PelB", "VHH", "HA", "spacer(GS)", "His6", "stop"],
              [("HA", "YPYDVPDY"), ("His6", "HHHHHH")]),
}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_each_bundled_map_derives_to_the_layout_on_record(name):
    paths = cl.list_bundled_plasmids()
    assert name in paths, f"{name}.dna is no longer bundled"
    derived = cl.derive_layout(cl.load_plasmid(paths[name]))
    layout, tags = EXPECTED[name]
    assert derived["layout"] == layout
    assert derived["tags"] == tags


def test_vectors_is_what_the_maps_derive_to():
    """The module-level table and the derivation cannot disagree, being one thing."""
    assert set(sg.VECTORS) == set(cl.list_bundled_plasmids())
    for name, path in cl.list_bundled_plasmids().items():
        derived = cl.derive_layout(cl.load_plasmid(path))
        assert sg.VECTORS[name]["tags"] == derived["tags"]
        assert sg.VECTORS[name]["layout"] == derived["layout"]


def test_every_bundled_vector_now_has_a_layout():
    """Three of the five had none, and were checked against pMES3's tags."""
    for name in ("pMES4", "pHEN6", "pHEN25"):
        assert sg.VECTORS[name]["tags"], f"{name} derived no tag"


def test_pmecs_ha_is_the_eight_residues_the_map_encodes():
    tags = dict(sg.VECTORS["pMECS"]["tags"])
    assert tags["HA"] == "YPYDVPDY"
    assert len(tags["HA"]) == 8
    assert tags["HA"] == sg.HA_SEARCH
    # and it is a prefix of the canonical epitope, so a cassette carrying the
    # nine-residue form still matches
    assert sg.HA_AA.startswith(tags["HA"])


def test_a_label_is_matched_on_word_boundaries_and_not_as_a_substring():
    """pMECS annotates gIII as `Protein 3 M13 phage`, which contains "ha".

    `cloning.TAG_LABELS`, the substring test used to find the cloning window, is
    satisfied by that label. Deriving a layout with it pulls the whole 1227 nt
    gIII in as an HA tag, 409 residues of phage coat protein reported as a tag
    the construct has to carry.
    """
    assert cl.tag_name_for("Protein 3 M13 phage") is None
    assert cl.tag_name_for("HA-tag") == "HA"
    assert cl.tag_name_for("6xHis") == "His"
    assert cl.tag_name_for("f1 ori") is None
    assert cl.tag_name_for("AmpR") is None
    # the substring test the window uses would have taken all of these
    assert any(lbl in "protein 3 m13 phage" for lbl in cl.TAG_LABELS)


def test_gIII_is_kept_out_three_times_over(monkeypatch):
    """Three things stop gIII, and each has to stop it on its own.

    The label is matched on word boundaries, the length is bounded, and the
    cassette ends at the first stop, which on pMECS sits between the His6 and
    gIII. They are redundant on the bundled vectors, which is the point: a map
    whose gIII is labelled with some other phage word gets past the label, and a
    vector carrying a genuinely long tag needs the length bound raised. Lifting
    each in turn shows what the others still catch.
    """
    import re as _re

    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()["pMECS"])
    substring_patterns = [(_re.compile(lbl), lbl.upper()) for lbl in cl.TAG_LABELS]

    derived = cl.derive_layout(plasmid)
    assert [n for n, _aa in derived["tags"]] == ["HA", "His6"]
    assert derived["beyond_stop"] == []

    # label alone: no length bound at all
    monkeypatch.setattr(cl, "MAX_TAG_NT", 10_000)
    assert [n for n, _aa in cl.derive_layout(plasmid)["tags"]] == ["HA", "His6"]

    # length alone: the substring matching that pulls "ha" out of "phage"
    monkeypatch.setattr(cl, "MAX_TAG_NT", 60)
    monkeypatch.setattr(cl, "TAG_PATTERNS", substring_patterns)
    assert all(len(aa) <= 20 for _n, aa in cl.derive_layout(plasmid)["tags"])

    # neither label nor length: the stop still keeps it out of the cassette, and
    # names it as something the vector carries downstream of it
    monkeypatch.setattr(cl, "MAX_TAG_NT", 10_000)
    loose = cl.derive_layout(plasmid)
    assert all(len(aa) <= 20 for _n, aa in loose["tags"]), loose["tags"]
    assert [(t, s) for t, s, _e in loose["beyond_stop"]] == [("HA", 400)], \
        "gIII should be named as carried past the stop, not expected in the cassette"


def test_a_vector_with_no_readable_tag_derives_nothing_rather_than_a_default():
    """The silent fallback is what this replaces, so it must not come back.

    An uploaded vector with no tag annotation used to be checked against pMES3,
    reporting a missing His6 that nobody ever expected it to carry.
    """
    empty = {"record": _bare_record("ATG" * 40), "name": "nothing",
             "topology": "linear", "format": "fasta", "warning": ""}
    derived = cl.derive_layout(empty, pd.DataFrame())
    assert derived["tags"] == []
    assert derived["layout"] == ["PelB", "VHH", "stop"]
    assert sg.tags_for(derived) == []


def test_an_unknown_vector_name_yields_no_tags():
    assert sg.tags_for("no such vector") == []
    assert sg.tags_for("pMECS") == sg.VECTORS["pMECS"]["tags"]
    assert sg.tags_for([("X", "AAA")]) == [("X", "AAA")]
    assert sg.tags_for({"tags": [("Y", "CCC")]}) == [("Y", "CCC")]


def test_verify_construct_takes_a_derived_layout_as_well_as_a_name():
    """The cloning screen passes a layout, the sequencing screen passes a name."""
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()["pMECS"])
    construct = _pmecs_clone(plasmid)
    derived = cl.derive_layout(plasmid)

    by_name = cl.verify_construct(construct, "pMECS")
    by_layout = cl.verify_construct(construct, derived)
    by_tags = cl.verify_construct(construct, tuple(derived["tags"]))

    for check in (by_name, by_layout, by_tags):
        assert check["ok"]
        assert check["tags_found"] == "HA,His6"
        assert check["tags_missing"] == ""


def _bare_record(dna: str):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    return SeqRecord(Seq(dna), id="x", name="x",
                     annotations={"molecule_type": "DNA"})


def _pmecs_clone(plasmid) -> str:
    """A pMECS backbone carrying a real VHH, cloned through PstI and BstEII."""
    from Bio.Data import CodonTable

    vhh = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
           "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGSRFSSWGQGTQVTVSS")
    back = CodonTable.unambiguous_dna_by_id[1].back_table
    seq = str(plasmid["record"].seq)
    fr1, fr4 = "CAGGTGCAGCTGCAG", "ACCCAGGTCACCGTCTCCTCA"
    i5, i3 = seq.find(fr1), seq.find(fr4)
    dna = fr1 + "".join(back[a] for a in vhh)[len(fr1):]
    dna = dna[:-len(fr4)] + fr4
    insert = seq[i5 - 25:i5] + dna + seq[i3 + len(fr4):i3 + len(fr4) + 8]
    built = cl.assemble(plasmid["record"], insert, "PstI", "BstEII", circular=True)
    assert built["ok"]
    return built["construct"]


def test_both_forms_of_the_ha_tag_are_recognised():
    """The epitope is nine residues; these vectors carry eight of them.

    The HA epitope is YPYDVPDYA, and `PROTEIN_MOTIFS` held exactly that, so the
    motif scan found it in none of the bundled vectors: pMECS was saved only by
    its map annotating `HA-tag` by hand, and pMES4's tag was invisible. Both
    forms are a tag. The note says which one is there, because a vector carrying
    eight residues where the epitope has nine is worth knowing about rather than
    worth rounding off.
    """
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    back = CodonTable.unambiguous_dna_by_id[1].back_table

    for epitope, expected_note in (("YPYDVPDYA", "detection tag"),
                                   ("YPYDVPDY", "truncated form, 8 of 9 residues")):
        protein = "MKYLLPTAAAGLLLLAAQPAMA" + epitope + "HHHHHH"
        dna = "".join(back[a] for a in protein)
        record = SeqRecord(Seq(dna), id="syn", name="syn",
                           annotations={"molecule_type": "DNA"})
        hits = [h for h in cl.detect_motifs(record, circular=False)
                if h["feature"] == "HA tag"]
        assert len(hits) == 1, f"{epitope}: {hits}"
        assert hits[0]["protein"] == epitope
        assert expected_note in hits[0]["note"]


def test_the_full_epitope_is_not_also_reported_as_its_own_prefix():
    """Two forms at one position is one feature, the longer one."""
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    back = CodonTable.unambiguous_dna_by_id[1].back_table
    dna = "".join(back[a] for a in "MKYLLPTAAAGLLLLAAQPAMAYPYDVPDYAHHHHHH")
    record = SeqRecord(Seq(dna), id="syn", name="syn",
                       annotations={"molecule_type": "DNA"})
    hits = [h for h in cl.detect_motifs(record, circular=False)
            if h["feature"] == "HA tag"]
    assert len(hits) == 1
    assert hits[0]["protein"] == "YPYDVPDYA"


def test_the_motif_table_and_the_sanger_constant_agree_on_what_ha_is():
    """Two places naming the same sequence, which is how they drifted before.

    Membership has to be tested against the forms as a collection, not against a
    bare string: every prefix of "YPYDVPDYA" is `in` that string, so a single
    nine-residue entry would satisfy a substring test while searching for only
    one of the two forms.
    """
    forms = cl.PROTEIN_MOTIFS["HA tag"][0]
    assert not isinstance(forms, str), \
        "one form only; the truncated form these vectors carry is not searched for"
    assert sg.HA_AA in tuple(forms), "the canonical epitope is not searched for"
    assert sg.HA_SEARCH in tuple(forms), "the form pMECS carries is not either"


def test_pmes4_carries_an_ha_tag_past_the_amber_stop_and_does_not_expect_it():
    """pMES4 encodes His6, an amber stop, then HA, then gIII.

    In a suppressor strain the amber is read through and the HA is on the phage
    fusion; in a non-suppressor strain the protein stops at the His6 and has no
    HA at all. `sanger.extract_cassette` looks for tags before the first stop,
    so expecting that HA reports `missing tag(s): HA` on every correct pMES4
    clone, which is the pMECS failure over again from the other end. It is named
    as carried past the stop instead of being expected or silently dropped.
    """
    derived = cl.derive_layout(cl.load_plasmid(cl.list_bundled_plasmids()["pMES4"]))
    assert [n for n, _aa in derived["tags"]] == ["His6"]
    assert [(t, s, e) for t, s, e in derived["beyond_stop"]] == [("HA", 361, 384)]

    # the HA is real: it is in the feature inventory, just not in the cassette
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()["pMES4"])
    features = cl.feature_table(plasmid)
    assert (features["feature"] == "HA tag").any(), \
        "the tag pMES4 carries should still be on the map"

    # and a correct pMES4 clone comes back clean
    check = cl.verify_construct(_clone_in("pMES4"), derived)
    assert check["ok"]
    assert check["tags_found"] == "His6"
    assert check["tags_missing"] == ""
    assert "missing tag(s): HA" not in check["warnings"]


def test_pmecs_still_expects_its_ha_because_that_one_precedes_the_stop():
    """The same rule, the other way round: pMECS puts HA before the His6."""
    derived = cl.derive_layout(cl.load_plasmid(cl.list_bundled_plasmids()["pMECS"]))
    assert [n for n, _aa in derived["tags"]] == ["HA", "His6"]
    assert derived["beyond_stop"] == []

    check = cl.verify_construct(_clone_in("pMECS"), derived)
    assert check["ok"]
    assert check["tags_found"] == "HA,His6"


def _clone_in(vector: str) -> str:
    """That vector's backbone carrying a real VHH, cut with PstI and BstEII."""
    from Bio.Data import CodonTable

    vhh = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG"
           "RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGSRFSSWGQGTQVTVSS")
    back = CodonTable.unambiguous_dna_by_id[1].back_table
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()[vector])
    seq = str(plasmid["record"].seq)
    fr1, fr4 = "CAGGTGCAGCTGCAG", "ACCCAGGTCACCGTCTCCTCA"
    i5, i3 = seq.find(fr1), seq.find(fr4)
    dna = fr1 + "".join(back[a] for a in vhh)[len(fr1):]
    dna = dna[:-len(fr4)] + fr4
    insert = seq[i5 - 25:i5] + dna + seq[i3 + len(fr4):i3 + len(fr4) + 8]
    built = cl.assemble(plasmid["record"], insert, "PstI", "BstEII", circular=True)
    assert built["ok"]
    return built["construct"]
