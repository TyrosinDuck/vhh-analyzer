"""
The bundled examples: the simulated chromatograms and the bundled structure.

An example that is wrong is worse than no example, because it is the first
thing somebody sees and they have nothing of their own to check it against.
So the simulated reads are checked end to end against what the app says about
them, and the sequences they yield are checked against the library entries
they were built from rather than against a value typed in here.

The `.ab1` writer gets its own round trip. It is the one piece of this module
that implements a published file format, and the only parser it has been
checked against is Biopython's, which is the one the app uses: these tests
prove the app can read what the app writes, and say nothing about an
instrument's own software.
"""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO

import demo
import library as lb
import sanger as sg
import structure as stx

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def example_files():
    return demo.sanger_example_files()


@pytest.fixture(scope="module")
def example_table(example_files):
    table, _reads = sg.process_reads(example_files, "pMES3")
    return table.set_index("read")


# --------------------------------------------------------------------------- #
# The writer
# --------------------------------------------------------------------------- #

def test_the_ab1_writer_round_trips_through_biopython():
    """Bases and quality values come back exactly as they went in."""
    sequence = "ACGTACGTTTGGCCAATTNNACGT" * 5
    quals = np.arange(len(sequence)) % 50 + 2
    raw = demo._ab1_bytes(sequence, quals, "round trip",
                          np.random.default_rng(0))

    assert raw[:4] == b"ABIF"
    record = SeqIO.read(io.BytesIO(raw), "abi")
    assert str(record.seq) == sequence
    assert record.letter_annotations["phred_quality"] == list(quals)


def test_every_example_file_says_it_is_simulated(example_files):
    """The claim has to survive the file leaving the app.

    A caption on the screen is gone the moment the file is downloaded, and a
    file that looks like a sequencing run and does not say otherwise is the
    one way this example could mislead somebody.
    """
    for name, raw in example_files:
        record = SeqIO.read(io.BytesIO(raw), "abi")
        assert "SIMULATED" in record.id, f"{name} does not name itself a simulation"
        comment = record.annotations["abif_raw"].get("CMNT1", b"").decode()
        assert "Simulated" in comment, f"{name} carries no comment saying so"


def test_the_example_is_the_same_bytes_every_time():
    """Reproducible, because an example that drifts cannot be checked.

    The quality values and the trace noise come from a seeded generator, and
    the date in the file is fixed rather than read off the clock, so two runs
    of the module produce identical files and a changed example means changed
    code.
    """
    first = demo.sanger_example_files()
    second = demo.sanger_example_files()
    assert [name for name, _ in first] == [name for name, _ in second]
    for (name, a), (_, b) in zip(first, second):
        assert a == b, f"{name} is not reproducible"


# --------------------------------------------------------------------------- #
# What the reads are meant to show
# --------------------------------------------------------------------------- #

def test_the_example_run_shows_every_status_the_table_reports(example_table):
    """One read of each kind, which is what the example is for."""
    assert dict(example_table["status"]) == {
        "clone_A1_fwd": "ok",
        "clone_A1_rev": "ok",
        "clone_B2_fwd": "ok",
        "clone_C3_fwd": "frameshift",
        "clone_D4_fwd": "no FR4",
        "empty_E5_fwd": "no PelB",
    }


def test_the_clean_reads_give_back_the_library_sequences(example_table):
    """The VHH extracted from the simulated read is the one put into it.

    Checked against `data/library.csv` rather than against a sequence written
    out here: a correction to the library then reaches this test instead of
    leaving it pinning the old sequence, which is the same rule the worked
    example in `examples/` follows.
    """
    library, _warnings = lb.load_library()
    library = library.set_index("name")
    for read, entry in (("clone_A1_fwd", "cAbLys3"),
                        ("clone_A1_rev", "cAbLys3"),
                        ("clone_B2_fwd", "caplacizumab")):
        assert example_table.loc[read, "vhh_aa"] == library.loc[entry, "protein"]


def test_the_reverse_read_is_found_on_the_minus_strand(example_table):
    """The point of reading one clone from both ends.

    Nothing tells the app which strand a read is on; it searches all six
    frames. If this ever came back `+` the search would have stopped working
    and the forward read alone would not show it.
    """
    assert example_table.loc["clone_A1_fwd", "strand"] == "+"
    assert example_table.loc["clone_A1_rev", "strand"] == "-"
    assert example_table.loc["clone_A1_rev", "duplicate_of"] == "clone_A1_fwd"


def test_the_clean_reads_carry_the_tag_the_vector_encodes(example_table):
    """pMES3 is PelB-VHH-His6-stop, and the example has to show that.

    This is also what would catch a clone built out of frame with the tag: the
    cassette would still translate and the His6 would not be there.
    """
    for read in ("clone_A1_fwd", "clone_A1_rev", "clone_B2_fwd"):
        assert example_table.loc[read, "tags_found"] == "His6"
        assert example_table.loc[read, "tags_missing"] == ""
        assert example_table.loc[read, "warnings"] == ""


def test_the_frameshifted_read_reports_no_sequence(example_table):
    """A frameshift is deliberately not reconstructed, so none must appear."""
    assert example_table.loc["clone_C3_fwd", "vhh_aa"] in ("", None)
    assert "frameshift" in example_table.loc["clone_C3_fwd", "warnings"]


def test_the_quality_profile_leaves_the_trimmer_something_to_do(example_table):
    """An example whose every base is perfect would not show trimming at all."""
    trimmed = example_table.loc["clone_A1_fwd", "trimmed_len"]
    assert 0 < trimmed < demo._READ_LENGTH, (
        "the trimmed window is the whole read, so the example shows nothing "
        "about quality trimming")


# --------------------------------------------------------------------------- #
# The bundled structure
# --------------------------------------------------------------------------- #

def test_the_structure_example_is_on_disk_and_loads():
    raw, name = demo.structure_example()
    assert name.endswith(".gz") and raw[:2] == b"\x1f\x8b"
    model = stx.load_structure(raw, name)
    chains = stx.chain_table(model).set_index("chain")
    assert list(chains.index) == ["A", "B", "L", "M"]
    assert chains.loc["A", "length"] == 132


def test_the_structure_example_declares_more_than_it_models():
    """1MEL chain A: the C-terminal tag is in SEQRES and not in the coordinates.

    This is the case the sequence panel exists for, so if the entry were ever
    swapped for one whose coordinates are complete, the panel would still be
    drawn and would have nothing to say.
    """
    raw, name = demo.structure_example()
    declared = stx.deposited_sequences(raw, name)
    chains = stx.chain_table(stx.load_structure(raw, name)).set_index("chain")
    coverage = stx.sequence_coverage(chains.loc["A", "sequence"], declared["A"])
    assert coverage["n_deposited"] == 148
    assert coverage["n_modelled"] == 132
    assert coverage["missing_c"] == 15
    assert coverage["mismatches"] == []


def test_the_structure_example_is_cited_where_it_is_offered():
    """A file taken from somebody else's deposition is named, not just used."""
    assert "1MEL" in demo.STRUCTURE_EXAMPLE["title"]
    assert "Desmyter" in demo.STRUCTURE_EXAMPLE["citation"]
    assert "10.1038/nsb0996-803" in demo.STRUCTURE_EXAMPLE["citation"]
    assert "CC0" in demo.STRUCTURE_EXAMPLE["licence"]
    readme = (ROOT / "data" / "examples" / "README.md").read_text(encoding="utf-8")
    assert "CC0" in readme and "Desmyter" in readme
