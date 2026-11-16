"""
demo, the example inputs the app offers when you have no file of your own.

Two screens take a file and can do nothing without one: Sanger sequencing wants
`.ab1` chromatograms, Structure wants a model. Somebody evaluating the app, or
learning it before a screening run, has neither to hand. This module supplies
one example for each.

**The Sanger example is simulated, and the app says so wherever it appears.**
The reads are not measurements. Each one is built here: a VHH from the bundled
sequence library is put into the bundled pMES3 map in place of its stuffer, a
window of the resulting construct is taken as the read, and a quality profile
of the usual shape, poor at the start, flat in the middle, decaying at the end,
is invented for it. The numbers in that profile carry no claim at all; they are
the shape of a trace, chosen so that the trimming step has something to do.
What the example does demonstrate is real: the reads go through
`sanger.process_reads` exactly as an upload does, the cassette is found in the
strand and frame it really sits in, and the tags are the ones the pMES3 map
itself annotates.

Writing `.ab1` rather than handing `sanger.py` a string is deliberate. The
example then exercises the same path as an upload, parser included, and the
files can be downloaded and opened in a chromatogram viewer, which is also how
somebody can check this module against a real one. The writer follows the ABIF
specification and has been checked by reading its own output back through
Biopython, which is the parser the app uses; it has not been checked against a
sequencing instrument's software.

**The structure example is a real file**, PDB entry 1MEL, shipped in
`data/examples/`. It is a crystal structure rather than a model, so there is no
pLDDT in it, and it is a complex, so the interface view has something to
measure. The PDB archive is public domain (CC0 1.0), and the entry is cited
wherever the app offers it.
"""

from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np

__all__ = [
    "SANGER_NOTE", "STRUCTURE_EXAMPLE", "sanger_example_files",
    "structure_example", "write_ab1", "__version__",
]

__version__ = "2026.11.16a"

# Said on every screen that offers the reads, and written into the files
# themselves, so a file that leaves the app still carries it.
SANGER_NOTE = ("Simulated chromatograms, not a sequencing run. The sequences "
               "are real VHHs put into the bundled pMES3 map; the traces and "
               "the quality values are invented.")

DATA_DIR = Path(__file__).parent / "data" / "examples"

STRUCTURE_EXAMPLE = {
    "file": "1mel.pdb.gz",
    "name": "1MEL.pdb",
    "title": "PDB 1MEL, cAbLys3 bound to hen egg-white lysozyme",
    "citation": ("Desmyter A, Transue TR, Ghahroudi MA, Thi MH, Poortmans F, "
                 "Hamers R, Muyldermans S, Wyns L (1996) Crystal structure of "
                 "a camel single-domain VH antibody fragment in complex with "
                 "lysozyme. Nat Struct Biol 3:803-811. "
                 "doi:10.1038/nsb0996-803"),
    "licence": ("PDB archive files are released under CC0 1.0 Universal "
                "(Public Domain Dedication); wwPDB asks that the depositors "
                "be credited, which the citation above does."),
}


# --------------------------------------------------------------------------- #
# The structure example
# --------------------------------------------------------------------------- #

def structure_example() -> tuple[bytes, str]:
    """The bundled PDB entry, as ``(bytes, filename)``.

    Returned gzipped and named `.pdb.gz`, which is what is on disk and what
    `structure.load_structure` reads directly: decompressing here would only
    mean holding four times as much in memory to no end.
    """
    path = DATA_DIR / STRUCTURE_EXAMPLE["file"]
    return path.read_bytes(), STRUCTURE_EXAMPLE["file"]


# --------------------------------------------------------------------------- #
# Writing ABIF
# --------------------------------------------------------------------------- #

_HEADFMT = ">H4sI2H3I"          # after the four-byte ABIF marker
_DIRFMT = ">4sI2H4I"
_DIR_ENTRY_SIZE = struct.calcsize(_DIRFMT)
_DATA_START = 128               # where real instruments put the first block

# Element type codes used here, from the ABIF specification: 2 is char, which
# is how base calls and quality values are stored, 4 is a signed short, which
# is how the trace channels and the peak positions are, 10 and 11 are the
# date and time, and 18 is a Pascal string, length byte first.
_CHAR, _SHORT, _DATE, _TIME, _PSTRING = 2, 4, 10, 11, 18
_ELEM_SIZE = {_CHAR: 1, _SHORT: 2, _DATE: 4, _TIME: 4, _PSTRING: 1}


def _pstring(text: str) -> bytes:
    """A Pascal string: one length byte, then at most 255 bytes of text."""
    raw = text.encode()[:255]
    return bytes([len(raw)]) + raw


def _tag(name: str, number: int, code: int, count: int, data: bytes) -> dict:
    # A Pascal string counts its length byte, so the element count is the
    # whole block: the reader unpacks the block as one string and drops that
    # first byte itself.
    if code == _PSTRING:
        count = len(data)
    return {"name": name.encode(), "number": number, "code": code,
            "count": count, "data": data}


def write_ab1(tags: list[dict]) -> bytes:
    """Assemble an ABIF file from directory entries.

    ``tags`` are the dicts `_tag` builds. A block of four bytes or fewer is
    held inside its own directory entry, in the field that otherwise holds the
    offset, which is what the format says and what Biopython's reader expects.
    """
    blocks, offsets, cursor = [], [], _DATA_START
    for tag in tags:
        if len(tag["data"]) <= 4:
            offsets.append(None)
            continue
        offsets.append(cursor)
        blocks.append(tag["data"])
        cursor += len(tag["data"])
        if cursor % 2:                      # keep blocks on even boundaries
            blocks.append(b"\x00")
            cursor += 1

    directory_offset = cursor
    out = bytearray(b"ABIF")
    out += struct.pack(_HEADFMT, 101, b"tdir", 1, 1023, _DIR_ENTRY_SIZE,
                       len(tags), len(tags) * _DIR_ENTRY_SIZE,
                       directory_offset)
    out += b"\x00" * (_DATA_START - len(out))
    for block in blocks:
        out += block

    for tag, offset in zip(tags, offsets):
        data = tag["data"]
        if offset is None:
            field = int.from_bytes(data.ljust(4, b"\x00"), "big")
        else:
            field = offset
        out += struct.pack(_DIRFMT, tag["name"], tag["number"], tag["code"],
                           _ELEM_SIZE[tag["code"]], tag["count"], len(data),
                           field, 0)
    return bytes(out)


# --------------------------------------------------------------------------- #
# Simulating one read
# --------------------------------------------------------------------------- #

# The shape of the invented quality profile, in Phred units and in bases.
# None of this is a measurement: it is a curve that rises, holds and falls like
# a Sanger read so that quality trimming has something to remove at each end.
_Q_START, _Q_PLATEAU, _Q_TAIL = 4.0, 52.0, 4.0
_Q_RAMP_BASES, _Q_DECAY_FROM = 30, 620
_Q_JITTER = 2.5

# One peak every ten points, which is the order of magnitude a real trace has.
_TRACE_SPACING = 10
_TRACE_WIDTH = 2.6
_TRACE_HEIGHT = 1400.0
_CHANNELS = "GATC"              # the order declared in the FWO_ tag


def _quality(n: int, rng, plateau: float = _Q_PLATEAU) -> np.ndarray:
    quals = np.full(n, plateau, dtype=float)
    ramp = min(_Q_RAMP_BASES, n)
    quals[:ramp] = np.linspace(_Q_START, plateau, ramp)
    if n > _Q_DECAY_FROM:
        quals[_Q_DECAY_FROM:] = np.linspace(plateau, _Q_TAIL,
                                            n - _Q_DECAY_FROM)
    quals += rng.normal(0.0, _Q_JITTER, n)
    return np.clip(quals, 2, 60).astype(int)


def _traces(seq: str, quals: np.ndarray, rng) -> tuple[dict, list[int]]:
    """Four channels of Gaussian peaks, one peak per base called."""
    n_points = (len(seq) + 2) * _TRACE_SPACING
    channels = {base: np.zeros(n_points) for base in _CHANNELS}
    positions = []
    window = np.arange(-3 * _TRACE_SPACING, 3 * _TRACE_SPACING + 1)
    bell = np.exp(-0.5 * (window / _TRACE_WIDTH) ** 2)
    for index, base in enumerate(seq):
        centre = (index + 1) * _TRACE_SPACING
        positions.append(centre)
        if base not in channels:            # an N is called into no channel
            continue
        # A confident call is a tall clean peak; a poor one is short, and the
        # other three channels pick up what it loses.
        confidence = float(np.clip(quals[index] / _Q_PLATEAU, 0.05, 1.0))
        slice_from = centre + window[0]
        slice_to = centre + window[-1] + 1
        if slice_from < 0 or slice_to > n_points:
            continue
        channels[base][slice_from:slice_to] += bell * _TRACE_HEIGHT * confidence
        for other in _CHANNELS:
            if other != base:
                channels[other][slice_from:slice_to] += (
                    bell * _TRACE_HEIGHT * (1.0 - confidence) / 3.0)
    for base in _CHANNELS:
        channels[base] += rng.normal(0.0, 8.0, n_points)
    return ({base: np.clip(values, 0, 32000).astype(np.int16)
             for base, values in channels.items()}, positions)


def _ab1_bytes(seq: str, quals: np.ndarray, sample: str, rng) -> bytes:
    channels, positions = _traces(seq, quals, rng)
    tags = [
        _tag("CMNT", 1, _PSTRING, 0, _pstring(SANGER_NOTE)),
        _tag("SMPL", 1, _PSTRING, 0, _pstring(sample)),
        _tag("FWO_", 1, _CHAR, 4, _CHANNELS.encode()),
        _tag("PBAS", 1, _CHAR, len(seq), seq.encode()),
        _tag("PBAS", 2, _CHAR, len(seq), seq.encode()),
        _tag("PCON", 1, _CHAR, len(quals), bytes(int(q) for q in quals)),
        _tag("PCON", 2, _CHAR, len(quals), bytes(int(q) for q in quals)),
        _tag("PLOC", 1, _SHORT, len(positions),
             struct.pack(f">{len(positions)}h", *positions)),
        _tag("PLOC", 2, _SHORT, len(positions),
             struct.pack(f">{len(positions)}h", *positions)),
    ]
    for offset, base in enumerate(_CHANNELS):
        values = channels[base]
        tags.append(_tag("DATA", 9 + offset, _SHORT, len(values),
                         struct.pack(f">{len(values)}h", *values.tolist())))
    # A date and a time are part of the format; these are fixed rather than
    # taken from the clock, so that two runs of this module produce the same
    # bytes and a changed example is visible as a changed file.
    for number in (1, 2):
        tags.append(_tag("RUND", number, _DATE, 1,
                         struct.pack(">h2B", 2026, 1, 1)))
        tags.append(_tag("RUNT", number, _TIME, 1,
                         struct.pack(">4B", 0, 0, 0, 0)))
    return write_ab1(tags)


# --------------------------------------------------------------------------- #
# The reads themselves
# --------------------------------------------------------------------------- #

_VECTOR = "pMES3"
# Where each read starts and how long it is, in nucleotides of the construct.
# The cassette of a pMES3 clone runs from about 225 to about 700, so a read
# that starts at 165 carries the whole of it inside the flat part of the
# quality profile, which is the ordinary case this example should show first.
_READ_START, _READ_LENGTH = 165, 780
_TRUNCATED_LENGTH = 330          # stops inside the VHH: no FR4 to find
_FRAMESHIFT_AT = 120             # nucleotides into the VHH

# (file name, library entry, what the read is meant to show)
_READS = [
    ("clone_A1_fwd", "cAbLys3", "forward read of a clean clone"),
    ("clone_A1_rev", "cAbLys3", "the same clone read from the other end"),
    ("clone_B2_fwd", "caplacizumab", "a second, different clone"),
    ("clone_C3_fwd", "anti-GFP-Enhancer", "one nucleotide deleted in the VHH"),
    ("clone_D4_fwd", "cAbLys3", "the read stops before the end of the VHH"),
    ("empty_E5_fwd", None, "an empty well"),
]


def _clone(vhh_dna: str, whole_codons: bool = True) -> str:
    """Put a VHH into the bundled pMES3 map, in place of its stuffer.

    The vector's own nucleotides are used for everything except the VHH: the
    leader up to the end of PelB, and from the first codon of the His6 tag
    onwards. What sits between the two in the empty vector is the stuffer and
    the cloning sites, which is exactly what a real clone replaces. Keeping
    the His6 segment whole is what keeps it in frame with the leader, as long
    as the VHH itself is a whole number of codons. ``whole_codons=False`` is
    how the frameshifted example is built, and it is the only reason that
    check can be lifted.
    """
    import cloning as cl

    plasmid = cl.load_plasmid(cl.list_bundled_plasmids()[_VECTOR])
    sequence = str(plasmid["record"].seq)
    features = cl.feature_table(plasmid)
    pelb = features[(features["feature"] == "PelB leader")
                    & (features["source"] == "motif scan")].iloc[0]
    his6 = features[(features["feature"] == "His6")
                    & (features["source"] == "motif scan")].iloc[0]
    if whole_codons and len(vhh_dna) % 3:
        raise ValueError("the VHH coding sequence is not whole codons")
    # `feature_table` reports 1-based inclusive spans, so an end is already a
    # Python slice bound and a start is one past where the slice begins.
    return (sequence[:int(pelb["end"])] + vhh_dna
            + sequence[int(his6["start"]) - 1:])


def _coding(protein_dna: str) -> str:
    """The library's DNA without its stop codon, which a clone does not carry."""
    from Bio.Seq import Seq

    dna = "".join(protein_dna.split()).upper()
    dna = dna[: len(dna) - len(dna) % 3]
    while dna and str(Seq(dna[-3:]).translate()) == "*":
        dna = dna[:-3]
    return dna


def _window(construct: str, start: int, length: int, reverse: bool) -> str:
    from Bio.Seq import Seq

    read = construct[start: start + length]
    return str(Seq(read).reverse_complement()) if reverse else read


def sanger_example_files() -> list[tuple[str, bytes]]:
    """The example chromatograms, as the uploader would hand them over.

    Six reads of a made-up screening plate: two clean clones, one of them read
    twice so that the duplicate grouping has something to group, one clone
    carrying a single-nucleotide deletion, one read that stops before the end
    of the VHH, and one empty well. Each is a simulation, see the module
    docstring and `SANGER_NOTE`.
    """
    import library as lb

    table, _warnings = lb.load_library()
    table = table.set_index("name")
    rng = np.random.default_rng(20261116)
    files = []
    for stem, entry, _purpose in _READS:
        if entry is None:
            # An empty well: no template, so the base caller returns noise.
            seq = "".join(rng.choice(list("ACGT"), 400))
            quals = np.clip(rng.normal(9.0, 4.0, len(seq)), 2, 25).astype(int)
            files.append((f"{stem}.ab1",
                          _ab1_bytes(seq, quals, f"{stem} (SIMULATED)", rng)))
            continue

        coding = _coding(str(table.loc[entry, "dna"]))
        shifted = stem.startswith("clone_C")
        if shifted:
            # One nucleotide gone: the tags downstream fall out of the
            # leader's frame, which is what the screen calls a frameshift.
            coding = coding[:_FRAMESHIFT_AT] + coding[_FRAMESHIFT_AT + 1:]
        construct = _clone(coding, whole_codons=not shifted)
        length = (_TRUNCATED_LENGTH if stem.startswith("clone_D")
                  else _READ_LENGTH)
        seq = _window(construct, _READ_START, length,
                      reverse=stem.endswith("_rev"))
        quals = _quality(len(seq), rng)
        files.append((f"{stem}.ab1",
                      _ab1_bytes(seq, quals, f"{stem} (SIMULATED)", rng)))
    return files
