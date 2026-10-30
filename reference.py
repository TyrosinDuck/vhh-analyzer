"""
reference - the two tables you look up rather than compute.

The genetic code and the twenty amino acids. Everything here is *derived*
rather than transcribed: the codon table comes out of Biopython, the usage
figures out of the same measured *E. coli* table `codons.py` scores CAI
against, the hydropathy and charge out of the same dictionaries `vhh_core.py`
uses for its profiles, and the masses out of Biopython's residue weights.

That is the whole point of the module. A reference table typed in by hand can
disagree with the app that sits next to it, and when it does, the user has no
way of telling which of the two is wrong. These cannot disagree, because they
are reading the app's own numbers.

The one thing that is not derived at run time is the twenty structure
drawings, in `data/aminoacids/`. They are produced once by
`tools/aminoacids/build.py`, which checks each one's formal charge and
molecular weight before writing it, and committed, so that the app does not
have to ship RDKit for twenty pictures that never change.

Importable from a notebook like the rest:

    import reference
    reference.genetic_code()
    reference.aminoacids()
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

import codons as cdn
import vhh_core as vc

__version__ = "2026.10.30a"

STRUCTURE_DIR = Path(__file__).parent / "data" / "aminoacids"
META_FILE = STRUCTURE_DIR / "meta.json"

# Three-letter code to one-letter, and the chemical family each belongs to.
# The families are a convention rather than a measurement, which is why they
# are here and not among `vhh_core`'s constants: nothing the app computes
# depends on them.
FAMILY = {
    "G": "Aliphatic", "A": "Aliphatic", "V": "Aliphatic", "L": "Aliphatic",
    "I": "Aliphatic", "P": "Aliphatic", "M": "Sulfur", "C": "Sulfur",
    "F": "Aromatic", "W": "Aromatic", "Y": "Aromatic",
    "S": "Polar", "T": "Polar", "N": "Polar", "Q": "Polar",
    "D": "Acidic", "E": "Acidic", "K": "Basic", "R": "Basic", "H": "Basic",
}

# Side-chain pKa, for the three residues whose ionisation state changes over
# the range the app's pH slider covers. These are the values behind the pH
# 7.4 forms the structures are drawn in, and the reason histidine is drawn
# neutral. Nelson & Cox, Lehninger, 7th ed.
SIDE_CHAIN_PKA = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07, "H": 6.00,
                  "K": 10.53, "R": 12.48}


@lru_cache(maxsize=1)
def _meta() -> dict:
    """What `tools/aminoacids/build.py` recorded alongside the drawings."""
    if not META_FILE.exists():
        return {"ph": None, "aminoacids": {}}
    with open(META_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def structure_ph() -> float | None:
    """The pH the committed structures were drawn at, or None if absent.

    Returned rather than hard-coded here so that the figure on the screen is
    the one the drawings were actually made with: if they are rebuilt at a
    different pH, the caption follows without anyone remembering to change it.
    """
    return _meta().get("ph")


def structure_svg(one_letter: str) -> str | None:
    """The committed drawing for one amino acid, as SVG text.

    Returns None when the file is missing, which is a deployment that did not
    copy `data/aminoacids/`, not a programming error: the screen says so and
    draws the rest of the table.
    """
    meta = _meta().get("aminoacids", {})
    code = next((three for three, row in meta.items()
                 if row.get("one") == one_letter.upper()), None)
    if code is None:
        return None
    path = STRUCTURE_DIR / f"{code}.svg"
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8")


@lru_cache(maxsize=1)
def aminoacids() -> pd.DataFrame:
    """One row per amino acid, every column read from somewhere else.

    The residue mass is Biopython's, minus a water, because that is the mass
    a residue contributes inside a chain and the mass the app's own molecular
    weights are built from. The free amino acid mass is the one in the
    structure metadata, which `build.py` checked against Biopython before
    writing it, so the two cannot drift apart.
    """
    from Bio.Data.IUPACData import protein_letters_3to1, protein_weights

    water = 18.01528
    meta = _meta().get("aminoacids", {})
    by_one = {row["one"]: (three, row) for three, row in meta.items()}

    rows = []
    for three, one in sorted(protein_letters_3to1.items(), key=lambda kv: kv[1]):
        drawn = by_one.get(one, (three, {}))[1]
        rows.append({
            "aa": one,
            "code": three,
            "name": _NAMES[one],
            "family": FAMILY[one],
            "residue_da": round(protein_weights[one] - water, 3),
            "free_da": drawn.get("mw_free_acid",
                                 round(protein_weights[one], 3)),
            "formula_at_ph": drawn.get("formula", ""),
            "charge_at_ph": drawn.get("charge"),
            "side_chain_pka": SIDE_CHAIN_PKA.get(one),
            "hydropathy": vc.KD[one],
            "charge_profile": vc.CHARGE.get(one, 0),
            "codons": " ".join(sorted(_codons_for(one))),
        })
    return pd.DataFrame(rows)


# The names are the only thing in this module that is typed out, because a
# full name is not derivable from anything the app already holds.
_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamate", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
}


def _table():
    from Bio.Data import CodonTable
    return CodonTable.unambiguous_dna_by_name["Standard"]


def _codons_for(one_letter: str) -> list[str]:
    return [codon for codon, aa in _table().forward_table.items()
            if aa == one_letter]


@lru_cache(maxsize=1)
def genetic_code() -> pd.DataFrame:
    """All 64 codons, with what each one means and how often *E. coli* uses it.

    The codon-to-amino-acid mapping is Biopython's standard table, the usage
    is the measured table in `data/ecoli_codon_usage.csv` that `codons.py`
    scores CAI against, and the rare-codon flags are `codons.RARE_CODONS`, the
    same ones the gene diagnosis screen raises. So a codon flagged rare here
    is flagged rare there, necessarily.

    ``share`` is the codon's share of its own amino acid's codons, which is
    what matters when choosing a synonym, and ``w`` is its relative
    adaptiveness, its share over the best synonym's, which is the quantity CAI
    is the geometric mean of. A raw per-thousand frequency is not something
    this app ever reasons about and is deliberately absent.
    """
    table = _table()
    _by_codon, by_aa = cdn.codon_usage()
    weights = cdn.relative_adaptiveness()

    rows = []
    for first in "TCAG":
        for second in "TCAG":
            for third in "TCAG":
                codon = first + second + third
                if codon in table.stop_codons:
                    amino, name = "*", "Stop"
                else:
                    amino = table.forward_table[codon]
                    name = _NAMES[amino]
                rare = cdn.RARE_CODONS.get(codon)
                rows.append({
                    "codon": codon,
                    "aa": amino,
                    "name": name,
                    "share": round(by_aa.get(amino, {}).get(codon, float("nan")), 3)
                    if amino != "*" else float("nan"),
                    "w": round(weights[codon], 3) if codon in weights
                    else float("nan"),
                    "rare_in_ecoli": bool(rare),
                    "tRNA": (rare or {}).get("trna") or "",
                    "note": (rare or {}).get("note") or "",
                })
    return pd.DataFrame(rows)


def codon_wheel() -> pd.DataFrame:
    """The genetic code laid out the way the textbook block is, 16 rows by 4.

    Rows are the first two bases, columns the third, which is the layout that
    makes the wobble position obvious: read across a row and the amino acid
    usually does not change.
    """
    code = genetic_code().set_index("codon")
    rows = []
    for first in "TCAG":
        for second in "TCAG":
            row = {"1st-2nd": first + second}
            for third in "TCAG":
                codon = first + second + third
                row[third] = code.loc[codon, "aa"]
            rows.append(row)
    return pd.DataFrame(rows).set_index("1st-2nd")
