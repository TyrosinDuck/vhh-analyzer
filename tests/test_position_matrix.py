"""
The bundled position-frequency matrix, checked against itself.

`data/vhh_position_freq.csv` is the reference the Repertoire screen calls a
sequence's residue unusual against, and the distribution the synthetic
generator samples frameworks from. `tools/plabdab/build_matrix.py` builds it
from a PLAbDab-nano download, and running that script with `--check` is what
proves the committed file is what the download produces.

That command is not run here. The download is a 2 MB file that is not in the
repository, for the same reason the sequences are not redistributed, and
numbering two thousand sequences with ANARCI takes about half a minute. So what
this file checks is the other half: that the committed matrix and its metadata
are internally consistent, and that the aggregation in the build script turns
numbering into counts the way the file says it does.

The consistency checks are not decoration. Every one of them corresponds to a
way the file could be wrong while still loading: a row whose frequencies do not
sum to 1 would make a rarity call meaningless at that position, a support
greater than the sequence count would mean a position was counted twice, a
metadata `n_sequences` that disagrees with the occupancy denominator would put
a wrong number in the caption the user reads before trusting the screen, and
rows out of IMGT order would be sampled out of order by the generator.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

import generator as gn
import repertoire as rep

ROOT = Path(__file__).resolve().parent.parent
MATRIX_FILE = ROOT / "data" / "vhh_position_freq.csv"
META_FILE = ROOT / "data" / "vhh_position_meta.json"
SCRIPT = ROOT / "tools" / "plabdab" / "build_matrix.py"

# The rounding the build script writes: frequencies to four places, occupancy
# to three. A row of twenty values rounded to four places can miss 1 by up to
# ten times half of the last place.
FREQ_TOLERANCE = 20 * 0.00005


@pytest.fixture(scope="module")
def matrix():
    return gn.load_matrix()


@pytest.fixture(scope="module")
def meta():
    return gn.load_meta()


def test_every_row_is_a_distribution_over_the_twenty_residues(matrix):
    """Each row's amino-acid columns are fractions of that position's residues.

    They are not fractions of the reference set: a position carried by half the
    sequences still has its own twenty columns summing to 1, and `occupancy` is
    what says how many sequences are behind it. Getting these two confused is
    the one arithmetic error in this file that would change every rarity call
    on the screen without changing anything that loads.
    """
    sums = matrix[gn.AA].sum(axis=1)
    bad = sums[(sums - 1.0).abs() > FREQ_TOLERANCE]
    assert bad.empty, f"these positions do not sum to 1: {bad.to_dict()}"


def test_occupancy_is_support_over_the_metadata_sequence_count(matrix, meta):
    """The denominator the file was built with is the one it reports.

    If these disagree, one of the two numbers was edited by hand, and the
    caption the user reads is then describing a different reference set than
    the frequencies underneath it.
    """
    n = meta["n_sequences"]
    expected = (matrix["n_residues"] / n).round(3)
    off = matrix.index[(expected - matrix["occupancy"]).abs() > 0.0005]
    assert len(off) == 0, f"occupancy disagrees with n_residues/{n} at {list(off)}"


def test_no_position_is_carried_by_more_sequences_than_exist(matrix, meta):
    assert int(matrix["n_residues"].max()) <= meta["n_sequences"]
    assert int(matrix["n_residues"].min()) >= 1, \
        "a position with no residue behind it has no business in the matrix"


def test_rows_are_in_imgt_order(matrix):
    """The generator walks the matrix in order; the file is already in it."""
    labels = list(matrix.index.astype(str))
    assert labels == sorted(labels, key=gn._imgt_sort_key)


def test_every_position_lands_in_a_region(matrix):
    """`region` is derived from the position number, so it must agree with it."""
    for label, region in matrix["region"].items():
        assert region == gn._region_of(str(label)), \
            f"{label} is labelled {region}"
    assert "?" not in set(matrix["region"]), \
        "a position outside every IMGT region boundary"


def test_cdr_length_counts_cover_every_reference_sequence(meta):
    """Each distribution counts all of them, once.

    The generator draws a CDR length from these, so a distribution that counted
    a subset would be a distribution of something other than the reference set.
    """
    for region in ("CDR1", "CDR2", "CDR3"):
        counts = meta["cdr_length_counts"][region]
        assert sum(counts.values()) == meta["n_sequences"], \
            f"{region} lengths do not account for every reference sequence"
        assert all(int(k) > 0 for k in counts), \
            f"{region} has a sequence with no CDR at all"


def test_the_metadata_says_where_the_matrix_came_from(meta):
    """The provenance block, which is the whole reason the build script exists.

    Before November 2026 this file named its source and nothing else: not the
    download, not its date, not what had been filtered out. A matrix nobody can
    trace is one nobody can rebuild or check, so its absence is a failure here
    rather than a missing nicety.
    """
    for field in ("source", "source_doi", "source_url", "numbering"):
        assert meta.get(field), f"the matrix does not record its {field}"
    prov = meta.get("provenance")
    assert prov, "the matrix records no provenance for the download behind it"
    for field in ("built_by", "built_from", "latest_entry", "n_entries",
                  "n_distinct_sequences", "min_length"):
        assert prov.get(field) is not None, f"provenance has no {field}"
    assert prov["built_by"] == "tools/plabdab/build_matrix.py"


def test_the_counts_in_the_provenance_add_up(meta):
    """Kept plus dropped is every distinct sequence in the download.

    This is what makes the filters auditable rather than merely described: a
    filter that drops something it does not report would show up here as a
    shortfall.
    """
    prov = meta["provenance"]
    dropped = (prov["dropped_nonstandard_alphabet"]
               + prov["dropped_shorter_than_min_length"]
               + prov["dropped_not_numbered"])
    assert meta["n_sequences"] + dropped == prov["n_distinct_sequences"], \
        "kept and dropped do not account for every distinct sequence"
    assert prov["n_distinct_sequences"] <= prov["n_entries"]


def test_reference_report_reads_the_file_rather_than_restating_it(matrix, meta):
    """What the screens put in their captions comes from the data.

    The numbers in the Repertoire caption and the generator's were written out
    by hand until November 2026. This is the check that they are derived: if
    `reference_report` ever stopped agreeing with the file, every caption would
    be stating a reference set the app is not using.
    """
    report = rep.reference_report()
    assert report["n_sequences"] == meta["n_sequences"]
    assert report["n_positions"] == len(matrix)
    assert report["max_support"] == int(matrix["n_residues"].max())
    assert report["median_support"] == int(matrix["n_residues"].median())
    assert report["provenance"]["latest_entry"] == \
        meta["provenance"]["latest_entry"]


def test_the_build_script_aggregates_a_numbering_the_way_it_claims():
    """`tabulate` on three hand-written numberings, where the answer is countable.

    ANARCI is not run here. The input is the shape ANARCI returns, a list of
    ``((position, insertion), residue)`` pairs with ``-`` for a gap, so this
    checks the arithmetic the script does with it and not the numbering itself:
    that gaps are excluded from the counts, that an insertion code becomes its
    own position, and that a CDR length is how many residues the numbering puts
    in that region.
    """
    sys.path.insert(0, str(SCRIPT.parent))
    import build_matrix as bm

    def numbering(pairs):
        return [((pos, ins), aa) for pos, ins, aa in pairs]

    one = numbering([(1, " ", "Q"), (2, " ", "V"), (27, " ", "G"),
                     (28, " ", "T"), (105, " ", "A"), (112, "A", "R"),
                     (113, " ", "Y")])
    two = numbering([(1, " ", "E"), (2, " ", "V"), (27, " ", "G"),
                     (28, " ", "-"), (105, " ", "A"), (113, " ", "Y")])
    three = numbering([(1, " ", "Q"), (2, " ", "-"), (27, " ", "S"),
                       (105, " ", "C"), (113, " ", "Y")])

    counts, lengths = bm.tabulate([one, two, three])

    assert dict(counts["1"]) == {"Q": 2, "E": 1}
    assert dict(counts["2"]) == {"V": 2}, "a gap must not be counted"
    assert dict(counts["27"]) == {"G": 2, "S": 1}
    assert dict(counts["112A"]) == {"R": 1}, "an insertion is its own position"
    assert "28" in counts and dict(counts["28"]) == {"T": 1}

    # CDR1 is IMGT 27-38: two residues in the first, one in each of the others.
    assert dict(lengths["CDR1"]) == {2: 1, 1: 2}
    # CDR3 is 105-117, and 112A falls inside it.
    assert dict(lengths["CDR3"]) == {3: 1, 2: 2}
    assert sum(lengths["CDR2"].values()) == 3, \
        "a sequence with no CDR2 residue still counts, at length zero"
    assert dict(lengths["CDR2"]) == {0: 3}


def test_the_build_script_writes_the_csv_the_app_reads():
    """`matrix_csv` output parses back into the shape `load_matrix` expects."""
    sys.path.insert(0, str(SCRIPT.parent))
    import build_matrix as bm
    import io

    import pandas as pd

    counts = {"1": {"Q": 3, "E": 1}, "105": {"A": 2}}
    text = bm.matrix_csv({k: __import__("collections").Counter(v)
                          for k, v in counts.items()}, n_sequences=4)
    df = pd.read_csv(io.StringIO(text), index_col=0)
    df.index = df.index.astype(str)
    assert list(df.columns) == ["region", "n_residues", "occupancy"] + gn.AA
    assert df.loc["1", "n_residues"] == 4
    assert df.loc["1", "Q"] == 0.75
    assert df.loc["105", "occupancy"] == 0.5
    assert df.loc["105", "region"] == "CDR3"


def test_the_build_script_refuses_a_download_missing_a_column(tmp_path):
    """A renamed column stops the build instead of changing the reference.

    PLAbDab-nano is self-updating and its export format can change. The failure
    that matters is not a crash, it is a build that quietly produces a matrix
    from a column that no longer means what it did.
    """
    bad = tmp_path / "wrong.csv"
    bad.write_text("seq,published\nQVQ,27-FEB-2020\n", encoding="utf-8")
    result = subprocess.run([sys.executable, str(SCRIPT), str(bad)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "missing columns" in (result.stderr + result.stdout)


def test_the_build_script_filters_and_counts_what_it_drops(tmp_path):
    """Distinct, standard-alphabet, long-enough, with every drop counted.

    ANARCI is not involved: `read_sequences` is the filtering stage on its own,
    which is what this checks.
    """
    sys.path.insert(0, str(SCRIPT.parent))
    import build_matrix as bm

    long_one = "Q" * bm.MIN_LENGTH
    with_x = "Q" * (bm.MIN_LENGTH - 1) + "X"
    short = "Q" * (bm.MIN_LENGTH - 1)
    rows = [("sequence", "update_date"),
            (long_one, "27-FEB-2020"),
            (long_one, "01-JAN-2021"),      # the same sequence, a second entry
            (with_x, "15-MAR-2019"),
            (short, "10-APR-2024"),
            ("", "10-APR-2024"),
            (long_one + "A", "not a date")]
    path = tmp_path / "download.csv"
    path.write_text("\n".join(",".join(r) for r in rows) + "\n",
                    encoding="utf-8")

    read = bm.read_sequences(path)
    assert read["sequences"] == sorted([long_one, long_one + "A"])
    assert read["n_entries"] == 6
    assert read["n_distinct"] == 4, "the blank row is not a sequence"
    assert read["dropped_nonstandard"] == 1
    assert read["dropped_short"] == 1
    assert read["undated_entries"] == 1
    assert read["latest_entry"] == "2024-04-10"
    assert read["earliest_entry"] == "2019-03-15"


def test_the_committed_metadata_is_the_json_the_script_writes():
    """Formatting, so that a rebuild of the same download is an empty diff.

    A file written with a different indent or key order would show every line
    as changed on the next rebuild, which is how a real change goes unnoticed.
    """
    text = META_FILE.read_text(encoding="utf-8")
    assert text == json.dumps(json.loads(text), indent=1) + "\n", \
        "the metadata is not in the format build_matrix.py writes"
