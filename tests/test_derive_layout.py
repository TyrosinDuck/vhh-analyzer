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
    # pHEN25 is the one map whose cassette carries something past the last
    # annotated tag: the sequence between the histidines and the first stop
    # codon encodes a proline-rich linker and then a single cysteine, neither of
    # them annotated. Until that was read from the sequence this entry said
    # PelB - VHH - His6 - stop, and a clone missing those fifteen residues
    # passed the construct check.
    "pHEN25": (["PelB", "VHH", "His6", "linker(SPSTPPTPSPSTPP)", "Cys(C-term)",
                "stop"],
               [("His6", "HHHHHH"), ("linker-Cys", "SPSTPPTPSPSTPPC")]),
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


# --------------------------------------------------------------------------- #
# What sits between the last annotated tag and the stop
# --------------------------------------------------------------------------- #

def test_phen25_carries_a_linker_and_a_cysteine_its_map_does_not_annotate():
    """The element that is in the sequence and in no annotation.

    pHEN25's map annotates `6xHis` and stops there, while the nucleotides
    between those histidines and the first stop codon encode SPSTPPTPSPSTPP and
    then one cysteine. A tag scan driven by annotations cannot see any of it, so
    the sequencing screen announced a cassette ending at the histidines and the
    construct check looked for nothing after them. These assertions read the
    stretch out of the shipped file, so they fail if the file is replaced by a
    map that encodes something else there.
    """
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()["pHEN25"])
    features = cl.feature_table(plasmid)
    derived = cl.derive_layout(plasmid, features)
    trailer = derived["trailer"]

    assert trailer["aa"] == "SPSTPPTPSPSTPPC"
    assert trailer["linker"] == "SPSTPPTPSPSTPP"
    assert trailer["cys"] is True
    assert trailer["labels"] == ["linker(SPSTPPTPSPSTPP)", "Cys(C-term)"]
    # the cysteine is the last residue before the stop, which is what makes it
    # a free C-terminal cysteine rather than one inside a linker
    assert trailer["aa"].endswith("C")

    # it is in the sequence at the coordinates reported, and nothing annotates it
    seq = str(plasmid["record"].seq)
    span = seq[trailer["start"] - 1:trailer["end"]]
    from Bio.Seq import Seq as _Seq
    assert str(_Seq(span).translate()) == trailer["aa"]
    annotated = features[(features["start"] >= trailer["start"])
                         & (features["end"] <= trailer["end"])]
    assert annotated.empty, \
        "the map now annotates this stretch; derive it from the annotation instead"


def test_the_other_four_vectors_stop_right_after_their_last_tag():
    """Deriving the trailer must not invent one where the stop follows the tag.

    This is what keeps the change to pHEN25 from moving the other four layouts,
    which `EXPECTED` above pins residue by residue.
    """
    for name in ("pMES3", "pMES4", "pHEN6", "pMECS"):
        derived = cl.derive_layout(cl.load_plasmid(cl.list_bundled_plasmids()[name]))
        assert derived["trailer"]["aa"] == "", f"{name}: {derived['trailer']}"
        assert derived["trailer"]["labels"] == []
        assert derived["layout"][-1] == "stop"
        assert derived["layout"][-2].startswith("His"), derived["layout"]


def test_a_trailer_with_no_stop_in_frame_is_not_guessed_at():
    """No stop downstream means the end of the cassette is not established."""
    assert cl.c_terminal_trailer("CCT" * 40, 0)["aa"] == ""          # no stop
    assert cl.c_terminal_trailer("TAA", 0)["aa"] == ""               # stop first
    assert cl.c_terminal_trailer("", 0)["aa"] == ""
    assert cl.c_terminal_trailer("CC", 0)["aa"] == ""                # no full codon


def test_a_trailer_without_a_cysteine_is_reported_as_a_linker_alone():
    """The cysteine is split out when it is there, not assumed to be."""
    # GGT GGT TCT TAA -> GGS, stop
    plain = cl.c_terminal_trailer("GGTGGTTCTTAA", 0)
    assert plain["aa"] == "GGS"
    assert plain["cys"] is False
    assert plain["linker"] == "GGS"
    assert plain["labels"] == ["linker(GGS)"]
    assert plain["name"] == "linker"

    # a lone cysteine, no linker in front of it
    lone = cl.c_terminal_trailer("TGTTAA", 0)
    assert lone["aa"] == "C"
    assert lone["cys"] is True
    assert lone["linker"] == ""
    assert lone["labels"] == ["Cys(C-term)"]
    assert lone["name"] == "Cys"


def test_a_correct_phen25_clone_reports_the_cysteine_as_present():
    """The construct check has to find it, which is the point of deriving it."""
    derived = cl.derive_layout(cl.load_plasmid(cl.list_bundled_plasmids()["pHEN25"]))
    check = cl.verify_construct(_clone_in("pHEN25"), derived)

    assert check["ok"]
    assert check["tags_found"] == "His6,linker-Cys"
    assert check["tags_missing"] == ""
    assert check["cassette_aa"].rstrip("*").endswith("SPSTPPTPSPSTPPC")


def test_a_phen25_construct_that_lost_the_cysteine_is_not_reported_as_correct():
    """Before the trailer was derived this construct passed as a good clone."""
    construct = _clone_in("pHEN25")
    # cut the cassette short: stop immediately after the histidines, which is
    # the layout the app used to expect for this vector
    i = construct.find("CACCACCATCACCATCAC")
    assert i > 0
    truncated = construct[:i + 18] + "TAATAG" + construct[i + 18 + 6:]

    check = cl.verify_construct(truncated, derived_phen25())
    assert check["tags_found"] == "His6"
    assert check["tags_missing"] == "linker-Cys"
    assert any("missing tag(s): linker-Cys" in w for w in check["warnings"])


def test_a_read_that_stops_inside_the_linker_is_short_not_a_bad_clone():
    """Fifteen residues past the His tag is far enough out that reads die there.

    Reporting that as a missing element reads as a clone that lost its cysteine,
    which is a different verdict and a different action: re-sequence rather than
    discard. The two are separable only when no stop was reached, so that is the
    only case in which the partial match is reported.
    """
    full = _clone_in("pHEN25")
    i = full.find("CACCACCATCACCATCAC")
    # a read ending part-way through the linker, with no stop codon reached
    short = full[:i + 18 + 24]

    inside = sg.extract_cassette(short, "pHEN25")
    assert inside["tags_missing"] == "linker-Cys"
    assert "linker-Cys (8 of 15 residues read)" == inside["tags_partial"]
    assert any("read ends inside" in w for w in inside["warnings"])

    # the whole thing read: nothing partial, nothing missing
    whole = sg.extract_cassette(full[:i + 18 + 45 + 6], "pHEN25")
    assert whole["tags_found"] == "His6,linker-Cys"
    assert whole["tags_partial"] == ""


def test_a_clone_that_really_lost_the_trailer_is_not_called_a_short_read():
    """The stop was reached, so the absence is the clone's and is reported flat."""
    full = _clone_in("pHEN25")
    i = full.find("CACCACCATCACCATCAC")
    lost = full[:i + 18] + "TAATAG"

    ex = sg.extract_cassette(lost, "pHEN25")
    assert ex["stop_codon_present"]
    assert ex["tags_missing"] == "linker-Cys"
    assert ex["tags_partial"] == "", "a clone that stops at the His tag is not a short read"


def derived_phen25():
    return cl.derive_layout(cl.load_plasmid(cl.list_bundled_plasmids()["pHEN25"]))
