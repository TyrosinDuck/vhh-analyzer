"""The hand-written vector layouts, checked against the plasmid maps beside them.

`sanger.VECTORS` says which tags follow the VHH in each vector, and the comment
under it says to take that from the actual plasmid map rather than from memory,
because a wrong layout makes the app report missing tags on clones that are
perfectly fine. Nothing enforced it, and it had already happened: the pMECS
entry carried the nine-residue HA epitope YPYDVPDYA while the pMECS map in
`data/plasmids/` encodes the eight residues YPYDVPDY followed by GS, so every
correct pMECS clone came back with `missing tag(s): HA` on the sequencing screen
and in the cloning screen's construct check.

These tests close that gap without changing how the layouts are written. Every
tag sequence in `VECTORS` has to be found in the map of the vector it names, and
the HA case is pinned on its own so the specific mistake cannot come back
quietly. Adding a vector whose layout disagrees with its map now fails the suite
instead of being discovered at the bench.

What they cannot check: whether the *set* of tags is complete. A map annotating
a tag the layout omits passes here, because the test asks only that what the
layout claims is really there.
"""

from __future__ import annotations

from Bio.Seq import Seq

import cloning as cl
import sanger as sg


def _six_frame_proteins(dna: str) -> list[str]:
    """Every reading frame of both strands, as protein strings.

    Six frames and not one: an empty phage-display vector carries a stuffer, so
    its tags are usually out of frame with the leader until a VHH of the right
    length is cloned in. Asking for the leader's frame would fail on vectors
    that are perfectly correct.
    """
    out = []
    for strand in (dna, str(Seq(dna).reverse_complement())):
        for frame in range(3):
            body = strand[frame:]
            body = body[: len(body) - len(body) % 3]
            if len(body) >= 3:
                out.append(str(Seq(body).translate()))
    return out


def _proteins_of(name: str) -> list[str]:
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()[name])
    return _six_frame_proteins(str(plasmid["record"].seq))


def test_every_layout_has_a_bundled_map_to_check_it_against():
    bundled = set(cl.list_bundled_plasmids())
    unchecked = sorted(set(sg.VECTORS) - bundled)
    assert not unchecked, (
        "these vectors are offered but ship no plasmid map, so their layout "
        f"cannot be checked against anything: {unchecked}. Add the map to "
        "data/plasmids/, or take the vector out of VECTORS.")


def test_every_tag_a_layout_claims_is_really_in_that_vector():
    missing = []
    for name, entry in sg.VECTORS.items():
        if name not in cl.list_bundled_plasmids():
            continue
        proteins = _proteins_of(name)
        for tag_name, tag_seq in entry["tags"]:
            if not any(tag_seq in protein for protein in proteins):
                missing.append(f"{name}: {tag_name} ({tag_seq})")
    assert not missing, (
        "these tags are in sanger.VECTORS and in none of the six frames of the "
        "vector's own map, so the app will report them missing on every correct "
        "clone:\n  " + "\n  ".join(missing)
        + "\nTake the sequence from the map, not from the epitope on record.")


def test_the_leader_anchors_are_in_the_vectors_that_use_them():
    for name in ("pMES3", "pMECS"):
        proteins = _proteins_of(name)
        assert any(sg.PELB_AA in protein for protein in proteins), name
        assert any(sg.PELB_CORE in protein for protein in proteins), name


def test_pmecs_encodes_eight_residues_of_ha_and_not_nine():
    """The specific mistake, pinned.

    If someone puts the nine-residue epitope back into the pMECS entry because
    it is the one on record, this fails and says why.
    """
    proteins = _proteins_of("pMECS")
    assert any("YPYDVPDY" in protein for protein in proteins)
    assert not any(sg.HA_AA in protein for protein in proteins), (
        "pMECS now contains the nine-residue HA epitope; if the map was "
        "replaced, this test and the HA_SEARCH constant need revisiting")
    assert sg.HA_SEARCH == "YPYDVPDY"
    assert dict(sg.VECTORS["pMECS"]["tags"])["HA"] == sg.HA_SEARCH


def test_a_correct_pmecs_clone_is_not_reported_as_missing_its_ha_tag():
    """End to end, through the function both screens use.

    A VHH is cloned into pMECS in silico and the finished cassette read back,
    which is what the sequencing screen does to an .ab1 and what the cloning
    screen does to a construct. Before the fix this returned
    `tags_missing: HA` and a `missing tag(s): HA` warning.
    """
    from Bio.Data import CodonTable

    back = CodonTable.unambiguous_dna_by_id[1].back_table
    protein = ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYA"
               "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGSRFSSWGQGTQVTVSS")
    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()["pMECS"])
    sequence = str(plasmid["record"].seq)
    fr1, fr4 = "CAGGTGCAGCTGCAG", "ACCCAGGTCACCGTCTCCTCA"
    start, end = sequence.find(fr1), sequence.find(fr4)
    coding = fr1 + "".join(back[a] for a in protein)[len(fr1):]
    coding = coding[: -len(fr4)] + fr4
    insert = (sequence[start - 25:start] + coding
              + sequence[end + len(fr4): end + len(fr4) + 8])

    built = cl.assemble(plasmid["record"], insert, "PstI", "BstEII",
                        circular=True)
    assert built["ok"], built["warnings"]

    check = cl.verify_construct(built["construct"], "pMECS")
    assert check["tags_found"] == "HA,His6"
    assert check["tags_missing"] == ""
    assert check["warnings"] == []
    assert check["vhh_aa"] == protein

    # and as a read, which is the sequencing screen's path
    read = built["construct"][150:1100]
    extracted = sg.extract_cassette(read, "pMECS")
    assert extracted["tags_missing"] == ""
    assert "missing tag(s): HA" not in extracted["warnings"]


def test_a_cassette_carrying_the_canonical_epitope_still_matches():
    """The fix is a relaxation, so it must not reject what it used to accept.

    The eight residues searched for are a prefix of the nine on record, so a
    construct carrying YPYDVPDYA matches as well as one carrying YPYDVPDY.
    """
    assert sg.HA_AA.startswith(sg.HA_SEARCH)
    cassette = sg.PELB_AA + "QVQLQ" + "AAA" + sg.HA_AA + "GAA" + "HHHHHH"
    assert sg.HA_SEARCH in cassette
