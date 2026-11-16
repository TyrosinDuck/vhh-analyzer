"""
Turn a PLAbDab-nano VHH download into the position-frequency matrix the app ships.

`generator.py` samples synthetic frameworks from `data/vhh_position_freq.csv`,
and the Repertoire screen reads the same file the other way round, to say which
residues in a submitted sequence are rare in natural VHHs. Until this script
existed, that matrix had no derivation in the repository: the source was named,
the licence was quoted, and nothing recorded which download, which filters or
which cut-off had produced it. It could not be rebuilt, and it could not be
checked.

This is the missing derivation, and the first thing it was used for was to
rebuild the matrix that was already committed. It reproduces
`data/vhh_position_freq.csv` byte for byte and every count in
`data/vhh_position_meta.json`, from the PLAbDab-nano VHH download of April 2024.
That is the only evidence that the filters written below are the ones that were
originally applied; they were recovered by reproducing the output, not read off
a note.

Run:  python tools/plabdab/build_matrix.py <vhh_sequences.csv[.gz]>
      python tools/plabdab/build_matrix.py <vhh_sequences.csv[.gz]> --check

`--check` builds into memory and compares against the committed files instead of
writing, exiting non-zero on any difference. That is the reproducibility claim
above, run as a command. It is not in the test suite: the input is a 2 MB
download that is not in the repository, and numbering two thousand sequences
with ANARCI takes about thirty seconds, so the suite checks the committed files
for internal consistency instead (`tests/test_position_matrix.py`) and this
command is what you run when you rebuild.

**Input.** The "VHH sequences (.csv.gz)" file from
https://opig.stats.ox.ac.uk/webapps/plabdab-nano/, compressed or not. One row
per database entry, so a nanobody deposited in three places appears three
times: the April 2024 download carries 4457 entries and 1975 distinct
sequences. The columns read are `sequence` and `update_date`, and the header is
checked before anything else happens. PLAbDab-nano is self-updating and its
export format can change; a renamed column has to stop the build rather than
quietly change the reference the whole Repertoire screen is read against.

**What is kept, with every count written to the metadata.**

* **Distinct sequences.** The matrix is a census of sequences, not of entries.
* **The twenty standard residues only.** Ten sequences in the April 2024
  download carry an `X`. `vhh_core.clean_sequence` refuses those outright, and a
  frequency computed over a position some sequences cannot fill is not a
  frequency.
* **At least `MIN_LENGTH` residues** (see that constant for how the value was
  chosen, and for what is not known about it).
* **Numbered by ANARCI**, IMGT scheme, alpaca germlines, which is the numbering
  the rest of the app uses. A sequence ANARCI cannot number is dropped and
  counted.

**What a row of the matrix is.** One IMGT position. `n_residues` is how many
kept sequences carry a residue there, `occupancy` is that over the number of
kept sequences, and the twenty amino-acid columns are the fraction of *those*
residues that are each amino acid, so each row's twenty columns sum to 1 and
not to the occupancy. Rows are in IMGT order, insertions included, taken from
`generator._imgt_sort_key` rather than reimplemented here: the file this script
writes is read by that module, and two orderings that agree until they do not
is the failure this repository has already had twice.

**The date is read from the data, not from the clock.** `update_date` is when
each entry was published, and the greatest one in the file is the most recent
sequence the reference can know about. The day you downloaded the file is not in
it, so `downloaded` is yours to pass with `--downloaded`; left out, it is
recorded as unknown rather than as today, which would stamp a rebuild of an old
download as fresh.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import io
import json
import sys
from datetime import date, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

# The order, the region boundaries and the alphabet all come from the module
# that reads the file this script writes. None of the three is restated here.
from generator import AA, _imgt_sort_key, _region_of  # noqa: E402

OUT_DIR = ROOT / "data"
MATRIX_FILE = OUT_DIR / "vhh_position_freq.csv"
META_FILE = OUT_DIR / "vhh_position_meta.json"

# The columns this script reads out of the download.
REQUIRED = ("sequence", "update_date")

# Shortest sequence the matrix counts, in residues before numbering.
#
# Seven sequences in the April 2024 download are fragments of 70 to 98
# residues, each missing FR1 or FR4 or both. ANARCI numbers them, so they are
# not caught by the numbering step, and they would add counts to the middle
# positions only: FR1 would then be a frequency over a different set of
# sequences than FR3, with nothing in the file saying so.
#
# What is honestly known about the value: any cut-off from 99 to 101 residues
# separates those seven from the 1958 that were originally kept, so the
# original cut-off is somewhere in that range and cannot be recovered more
# precisely than that. 99 is the low end, which keeps the most data. The
# shortest kept sequence is 99 residues and the next fragment down is 98, so
# the gap is one residue wide and nothing sits in it.
MIN_LENGTH = 99

# The citation the matrix carries, and the sentence from the paper's
# data-availability statement. Quoted, not paraphrased, because what a
# redistributor may do with the data is the thing a reader needs verbatim.
SOURCE = "PLAbDab-nano (OPIG, University of Oxford)"
SOURCE_URL = "https://opig.stats.ox.ac.uk/webapps/plabdab-nano/"
SOURCE_DOI = "10.1093/nar/gkae881"
SOURCE_TITLE = ("PLAbDab-nano: a database of camelid and shark nanobodies "
                "from patents and literature")
SOURCE_AVAILABILITY = ("PLAbDab-nano is freely available to download and "
                       "query at https://opig.stats.ox.ac.uk/webapps/"
                       "plabdab-nano/")
DERIVED = ("Aggregate position-frequency matrix derived from the VHH subset; "
           "the sequences themselves are not redistributed here.")

NUMBERING = "IMGT via ANARCI, alpaca germlines"
SPECIES = ["alpaca"]
SCHEME = "imgt"

STANDARD = set(AA)


def _open(path: Path):
    """The download is served gzipped, and browsers sometimes unzip it."""
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace", newline="")
    return path.open(newline="", encoding="utf-8", errors="replace")


def _parse_date(raw: str):
    """`27-FEB-2020` to a date, or None.

    GenBank's day-month-year form is what PLAbDab-nano carries. Anything else
    returns None and is counted rather than guessed at: a misparsed date would
    silently move the extraction date the app reports.
    """
    try:
        return datetime.strptime(raw.strip(), "%d-%b-%Y").date()
    except (TypeError, ValueError):
        return None


def read_sequences(export: Path) -> dict:
    """Distinct, standard-alphabet, long-enough sequences, and what was dropped."""
    with _open(export) as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"the download is missing columns this script needs: {missing}. "
                f"PLAbDab-nano's export format has changed; update REQUIRED and "
                f"the fields that depend on it rather than working around it.")
        rows = list(reader)

    seen, kept = set(), []
    dropped_nonstandard, dropped_short = set(), []
    dates, undated = [], 0
    for row in rows:
        when = _parse_date(row.get("update_date", ""))
        if when is None:
            undated += 1
        else:
            dates.append(when)
        sequence = (row["sequence"] or "").strip().upper()
        if not sequence or sequence in seen:
            continue
        seen.add(sequence)
        if set(sequence) - STANDARD:
            dropped_nonstandard.add(sequence)
            continue
        if len(sequence) < MIN_LENGTH:
            dropped_short.append(sequence)
            continue
        kept.append(sequence)

    # Sorted, so that a rebuild of the same download numbers the sequences in
    # the same order. Nothing downstream depends on the order, but a build that
    # is reproducible only up to a shuffle is not reproducible.
    return {
        "sequences": sorted(kept),
        "n_entries": len(rows),
        "n_distinct": len(seen),
        "dropped_nonstandard": len(dropped_nonstandard),
        "dropped_short": len(dropped_short),
        "latest_entry": max(dates).isoformat() if dates else None,
        "earliest_entry": min(dates).isoformat() if dates else None,
        "undated_entries": undated,
    }


def number(sequences: list[str], ncpu: int = 4) -> tuple[list, int]:
    """IMGT numbering for each sequence, and how many ANARCI could not number."""
    from anarci import anarci

    numbered, _details, _hits = anarci(
        [(f"s{i}", s) for i, s in enumerate(sequences)],
        scheme=SCHEME, allowed_species=SPECIES, output=False, ncpu=ncpu)
    kept = [n[0][0] for n in numbered if n is not None]
    return kept, len(sequences) - len(kept)


def tabulate(numberings: list) -> tuple[dict, dict]:
    """Per-position residue counts, and the CDR length of each sequence."""
    counts = collections.defaultdict(collections.Counter)
    lengths = {region: collections.Counter()
               for region in ("CDR1", "CDR2", "CDR3")}
    for numbering in numberings:
        per_region = collections.Counter()
        for (position, insertion), residue in numbering:
            if residue == "-":
                continue
            label = f"{position}{insertion.strip()}"
            counts[label][residue] += 1
            region = _region_of(label)
            if region in lengths:
                per_region[region] += 1
        for region in lengths:
            lengths[region][per_region[region]] += 1
    return counts, lengths


def matrix_csv(counts: dict, n_sequences: int) -> str:
    """The shipped CSV, as text, so that `--check` can compare it as text."""
    out = ["imgt,region,n_residues,occupancy," + ",".join(AA)]
    for label in sorted(counts, key=_imgt_sort_key):
        total = sum(counts[label].values())
        row = [label, _region_of(label), str(total),
               str(round(total / n_sequences, 3))]
        row += [str(round(counts[label][residue] / total, 4)) for residue in AA]
        out.append(",".join(row))
    return "\n".join(out) + "\n"


def meta_json(read: dict, lengths: dict, n_sequences: int, n_positions: int,
              unnumbered: int, export: Path, downloaded: str | None) -> str:
    meta = {
        "source": SOURCE,
        "source_url": SOURCE_URL,
        "n_sequences": n_sequences,
        "numbering": NUMBERING,
        "cdr_length_counts": {
            region: {str(k): v for k, v in sorted(lengths[region].items())}
            for region in ("CDR1", "CDR2", "CDR3")},
        "source_doi": SOURCE_DOI,
        "source_title": SOURCE_TITLE,
        "source_availability": SOURCE_AVAILABILITY,
        "derived": DERIVED,
        # Everything below this line is what the file did not carry before
        # November 2026, and the reason this script exists. Without it, a
        # reader cannot tell which download the matrix came from, how old its
        # newest sequence is, or what was thrown away to get it.
        "provenance": {
            "built_by": "tools/plabdab/build_matrix.py",
            "built_from": export.name,
            "downloaded": downloaded,
            "latest_entry": read["latest_entry"],
            "earliest_entry": read["earliest_entry"],
            "n_entries": read["n_entries"],
            "n_distinct_sequences": read["n_distinct"],
            "n_positions": n_positions,
            "min_length": MIN_LENGTH,
            "dropped_nonstandard_alphabet": read["dropped_nonstandard"],
            "dropped_shorter_than_min_length": read["dropped_short"],
            "dropped_not_numbered": unnumbered,
            "undated_entries": read["undated_entries"],
        },
    }
    return json.dumps(meta, indent=1) + "\n"


def build(export: Path, downloaded: str | None = None) -> tuple[str, str, dict]:
    read = read_sequences(export)
    numberings, unnumbered = number(read["sequences"])
    counts, lengths = tabulate(numberings)
    n = len(numberings)
    csv_text = matrix_csv(counts, n)
    json_text = meta_json(read, lengths, n, len(counts), unnumbered, export,
                          downloaded)
    report = dict(read, n_sequences=n, n_positions=len(counts),
                  dropped_not_numbered=unnumbered)
    return csv_text, json_text, report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("export", type=Path,
                        help="the PLAbDab-nano VHH download, .csv or .csv.gz")
    parser.add_argument("--check", action="store_true",
                        help="compare with the committed files instead of "
                             "writing, and exit non-zero on any difference")
    parser.add_argument("--downloaded", metavar="YYYY-MM-DD", default=None,
                        help="the day you downloaded the file, recorded in the "
                             "metadata; omitted, it is recorded as unknown")
    args = parser.parse_args(argv)

    if args.downloaded:
        try:
            date.fromisoformat(args.downloaded)
        except ValueError:
            raise SystemExit("--downloaded wants a YYYY-MM-DD date")

    csv_text, json_text, report = build(args.export, args.downloaded)

    print(f"{report['n_entries']} entries, "
          f"{report['n_distinct']} distinct sequences")
    print(f"  dropped, non-standard alphabet: {report['dropped_nonstandard']}")
    print(f"  dropped, shorter than {MIN_LENGTH} residues: "
          f"{report['dropped_short']}")
    print(f"  dropped, not numbered by ANARCI: "
          f"{report['dropped_not_numbered']}")
    print(f"  kept: {report['n_sequences']} sequences, "
          f"{report['n_positions']} IMGT positions")
    print(f"  newest entry in the download: {report['latest_entry']}")

    if args.check:
        differences = []
        for path, built in ((MATRIX_FILE, csv_text), (META_FILE, json_text)):
            if not path.exists():
                differences.append(f"{path.name} is not committed")
            elif path.read_text(encoding="utf-8") != built:
                differences.append(f"{path.name} differs from the rebuild")
        if differences:
            for line in differences:
                print(f"DIFFERS: {line}", file=sys.stderr)
            print("The committed matrix is not what this download and these "
                  "filters produce. Either the download is a different one, or "
                  "a filter here has changed.", file=sys.stderr)
            return 1
        print("The committed files are exactly what this download rebuilds.")
        return 0

    MATRIX_FILE.write_text(csv_text, encoding="utf-8")
    META_FILE.write_text(json_text, encoding="utf-8")
    print(f"wrote {MATRIX_FILE.relative_to(ROOT)} and "
          f"{META_FILE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
