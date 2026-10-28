"""
structure, what a 3D model says, as opposed to what the sequence suggests.

Every other module in this app reasons from sequence alone. That is the right
default, because a sequence is what you always have. But several of the questions
the app answers are really structural questions being approximated from sequence,
and when a model exists the approximation can simply be replaced by a measurement:

- ``vhh_core.disulfides`` *guesses* which cysteines pair, from IMGT positions and
  a canonical-pair assumption. Here the SG–SG distances are measured.
- ``vhh_core.find_liabilities`` reports a deamidation or oxidation site wherever
  the motif occurs. A buried motif is not a liability in the same sense as an
  exposed one; relative solvent accessibility tells them apart.
- ``vhh_core.aggregation_profile`` scores hydrophobicity along the chain. What
  actually drives aggregation is hydrophobic surface that is *contiguous in space*,
  which a sequence window cannot see and a structure can.
- Which residues contact the antigen is not inferable from sequence at all.

Accepts `.pdb`, `.cif`/`.mmcif` and gzipped versions of either, including the
outputs of AlphaFold2, AlphaFold3 and ColabFold.

**What this is not.** It measures one static conformer of one model. It does not
predict structure, does not dock, does not estimate affinity, and does not model
flexibility. For a predicted model, side-chain placement is markedly less reliable
than the backbone, so the solvent accessibility of long side chains (Lys, Arg, Glu)
should be read as indicative. pLDDT is the model's confidence in itself; it is
neither accuracy nor a measurement of disorder, though low pLDDT and disorder do
tend to coincide.

Solvent accessibility is computed with Biopython's Shrake–Rupley implementation
(1.4 Å probe) and normalised by the residue maxima of Tien et al. 2013. One
caveat, stated because it is easy to miss: those maxima were derived with DSSP,
i.e. the Lee–Richards definition. Shrake–Rupley is a different numerical
approximation of the same quantity, so relative accessibilities here are close to,
but not identical with, DSSP-normalised values. The exposed/buried call is robust
to that; a relative accessibility quoted to the percent is not.

  Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum Allowed
  Solvent Accessibilites of Residues in Proteins. PLoS ONE 8(11): e80635.
  doi:10.1371/journal.pone.0080635
"""

from __future__ import annotations

import gzip
import io
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio.Data import PDBData
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

import vhh_core as vc

__all__ = [
    "TIEN_MAX_ASA", "BURIED_MAX", "EXPOSED_MIN", "SS_MAX_DISTANCE",
    "load_structure", "chain_table", "chain_sequence", "confidence_kind",
    "confidence_table", "disulfides_3d", "compare_disulfides",
    "residue_sasa", "exposure_table", "interface_table", "interface_summary",
    "surface_patches", "annotate_liabilities", "map_structure_to_sequence",
    "resnum_map", "identity_resnum_map", "residue_regions",
    "StructureError", "PATCH_RESIDUES", "INTERFACE_MIN_DELTA",
    "plot_confidence", "plot_exposure", "plot_interface",
    "viewer_html", "PLDDT_COLOURS", "EXPOSURE_COLOURS",
    "PATCH_CONTACT", "PATCH_MIN_SIZE", "PATCH_MIN_REL", "__version__",
]

__version__ = "2026.10.28e"


def _three_to_one(resname: str) -> str | None:
    """Three-letter residue name to one letter, or None if it is not a residue.

    Biopython 1.80 removed ``three_to_one``; the extended table it kept covers
    the modified residues a real structure actually contains (MSE for
    selenomethionine above all, which is most of the cysteine-free structures in
    the PDB) and maps them to their parent amino acid.
    """
    return PDBData.protein_letters_3to1_extended.get(resname.strip().upper())


class StructureError(ValueError):
    """A structure file that cannot be used, with a reason fit to show a user."""


# Tien et al. 2013, Table 1, "theoretical" column, Å². The theoretical values are
# the ones the paper recommends as the upper bound on observable accessibility;
# the empirical column is systematically a few percent lower. Verified against the
# published table rather than reproduced from memory.
TIEN_MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

# The usual convention: at or below 20% of maximum is buried, above 50% is
# exposed, and the band between is called intermediate rather than forced into
# one of the two.
BURIED_MAX = 0.20
EXPOSED_MIN = 0.50

# A disulfide S–S bond is ~2.05 Å. 2.5 Å is comfortably past any real bond and
# still far short of the ~3.5 Å at which two unbonded sulfurs can sit.
SS_MAX_DISTANCE = 2.5

# Residues counted as hydrophobic when looking for surface patches. Deliberately
# the strongly hydrophobic set only: Ala and Gly are too small and too common to
# make a patch mean anything.
PATCH_RESIDUES = set("VILMFWY")
# Two residues count as adjacent when their side chains come within this
# distance, measured between closest heavy atoms rather than between centroids.
# Centroid distance looks simpler and is worse: Trp and Arg carry their bulk far
# from the backbone, so two side chains in van der Waals contact can have
# centroids 9 Å apart, and a patch built on centroids quietly misses exactly the
# large residues that matter most to it.
PATCH_CONTACT = 5.0
PATCH_MIN_SIZE = 3       # fewer than three residues is not a patch
# How exposed a residue must be to count as part of a surface patch. Deliberately
# well below the 0.50 that `exposure_table` calls "exposed", for a structural
# reason: a residue sitting in the middle of a patch is half-covered by the
# patch's own other residues, so requiring each one to be exposed in its own
# right is a good way of never finding a patch.
#
# On the one complex this was checked against, the HIV-1 capsid C-terminal
# domain, the known W184/M185 dimerisation patch is found at 0.15 and missed at
# 0.20, and the contact distance barely matters over 5.0–7.0 Å. That cliff is one
# residue wide (M185 sits at 0.195), so it is an observation on a single
# structure and emphatically not a validated parameter. It is a slider in the app
# for that reason, and the app says so when it finds nothing.
PATCH_MIN_REL = 0.15

# ΔSASA above which a residue counts as part of an interface. 1 Å² is small
# enough to catch a glancing contact and large enough to exclude numerical noise.
INTERFACE_MIN_DELTA = 1.0

MAX_ATOMS = 200_000


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

@dataclass
class Model:
    """One parsed structure, with the bits of provenance worth showing."""
    structure: object
    name: str
    fmt: str
    n_models: int
    header: dict = field(default_factory=dict)
    ssbond: list = field(default_factory=list)

    @property
    def model(self):
        return next(self.structure.get_models())


def _decompress(data: bytes) -> bytes:
    if data[:2] == b"\x1f\x8b":
        return gzip.decompress(data)
    return data


def _parse_ssbond(text: str) -> list[dict]:
    """SSBOND records, when the file is a PDB that carries them.

    These are the depositor's annotation, not a measurement. They are read so
    that the measured bonds can be checked against them, if the two disagree
    that is worth knowing, and silently trusting either one would hide it.
    """
    bonds = []
    for line in text.splitlines():
        if not line.startswith("SSBOND"):
            continue
        try:
            bonds.append({
                "chain_1": line[15].strip() or "?",
                "resnum_1": int(line[17:21]),
                "chain_2": line[29].strip() or "?",
                "resnum_2": int(line[31:35]),
                "declared_distance": (float(line[73:78])
                                      if line[73:78].strip() else None),
            })
        except (ValueError, IndexError):
            continue
    return bonds


def load_structure(data: bytes | str, filename: str = "structure") -> Model:
    """Parse a structure file. Raises ``StructureError`` with a readable reason."""
    if isinstance(data, str):
        data = data.encode()
    try:
        data = _decompress(data)
    except Exception as exc:
        raise StructureError(f"could not decompress {filename}: {exc}") from exc

    lower = filename.lower()
    for suffix in (".gz", ".bz2"):
        if lower.endswith(suffix):
            lower = lower[: -len(suffix)]
    is_cif = lower.endswith((".cif", ".mmcif"))

    try:
        text = data.decode("utf-8", errors="replace")
    except Exception as exc:
        raise StructureError(f"{filename} is not text: {exc}") from exc

    # Some files arrive with the wrong extension; the content is authoritative.
    if not is_cif and ("_atom_site." in text and "\nATOM  " not in text):
        is_cif = True

    parser = (MMCIFParser(QUIET=True) if is_cif
              else PDBParser(QUIET=True, PERMISSIVE=True))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            structure = parser.get_structure(filename, io.StringIO(text))
        except Exception as exc:
            raise StructureError(
                f"{filename} could not be parsed as "
                f"{'mmCIF' if is_cif else 'PDB'}: {exc}") from exc

    models = list(structure.get_models())
    if not models:
        raise StructureError(f"{filename} contains no model.")
    n_atoms = sum(1 for _ in models[0].get_atoms())
    if not n_atoms:
        raise StructureError(f"{filename} contains no atoms.")
    if n_atoms > MAX_ATOMS:
        raise StructureError(
            f"{filename} has {n_atoms:,} atoms in its first model, over the "
            f"{MAX_ATOMS:,} limit. Strip it to the chains you care about first.")

    return Model(structure=structure, name=filename,
                 fmt="mmCIF" if is_cif else "PDB",
                 n_models=len(models),
                 header=dict(getattr(structure, "header", {}) or {}),
                 ssbond=[] if is_cif else _parse_ssbond(text))


# --------------------------------------------------------------------------
# Chains
# --------------------------------------------------------------------------

def _residues(chain):
    """Amino-acid residues of a chain, in order, waters and ligands excluded.

    Modified residues are kept. Selenomethionine in particular is flagged as a
    heteroatom by the PDB but is unambiguously part of the chain, dropping it
    silently shortens the sequence, puts a hole in the middle of the numbering,
    and shifts every downstream mapping. It is mapped to its parent amino acid
    instead, which is what ``PDBData``'s extended table is for.
    """
    kept = []
    for residue in chain:
        hetflag = residue.id[0]
        if hetflag == "W":                      # water
            continue
        if hetflag == " ":
            if is_aa(residue, standard=False):
                kept.append(residue)
            continue
        # A heteroatom residue counts only if it is a modified amino acid.
        if hetflag.startswith("H_") and _three_to_one(residue.get_resname()):
            kept.append(residue)
    return kept


def _strip_solvent(entity):
    """Remove waters in place, so they do not occlude surface.

    Crystallographic waters are part of the model, not of the protein. Left in,
    Shrake-Rupley treats each one as an obstacle and every residue it sits on
    comes back more buried than it is, on a typical structure that is enough to
    turn exposed residues into intermediate ones and lose surface patches
    entirely. Ligands are deliberately left alone: a bound ligand really does
    occlude the surface underneath it, and saying so is the honest answer.
    """
    for chain in entity:
        for residue in [r for r in chain if r.id[0] == "W"]:
            chain.detach_child(residue.id)
    return entity


def chain_sequence(chain) -> tuple[str, list]:
    """One-letter sequence of a chain, with the residue objects it came from."""
    seq, kept = [], []
    for residue in _residues(chain):
        letter = _three_to_one(residue.get_resname())
        if letter is None:
            continue
        seq.append(letter)
        kept.append(residue)
    return "".join(seq), kept


def _numbering_gaps(residues) -> int:
    """Breaks in the residue numbering, unmodelled loops, usually."""
    gaps = 0
    for previous, current in zip(residues, residues[1:]):
        if current.id[1] - previous.id[1] > 1:
            gaps += 1
    return gaps


def chain_table(model: Model) -> pd.DataFrame:
    """One row per polypeptide chain of the first model."""
    rows = []
    for chain in model.model:
        seq, residues = chain_sequence(chain)
        if not seq:
            continue
        b = [a.get_bfactor() for r in residues for a in r]
        rows.append({
            "chain": chain.id,
            "length": len(seq),
            "first": residues[0].id[1],
            "last": residues[-1].id[1],
            "numbering_gaps": _numbering_gaps(residues),
            "mean_b": round(float(np.mean(b)), 1) if b else None,
            "sequence": seq,
        })
    if not rows:
        raise StructureError(
            "no polypeptide chain found, the file may contain only nucleic "
            "acids, ligands or waters.")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Confidence: pLDDT or crystallographic B-factor
# --------------------------------------------------------------------------

def confidence_kind(model: Model) -> tuple[str, str]:
    """Guess whether the B-factor column holds pLDDT, and say why.

    The B-factor column of a predicted model carries pLDDT instead; the two
    mean opposite things (high pLDDT is good, high B is bad) so reading one
    as the other inverts every conclusion. There is no field that states which
    it is, so this guesses from the values and the header, and the app lets the
    guess be overridden. Returns ``(kind, reason)``.
    """
    values = np.array([a.get_bfactor() for a in model.model.get_atoms()],
                      dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return "unknown", "the file carries no B-factor column."

    method = str(model.header.get("structure_method", "") or "").lower()
    if any(k in method for k in ("x-ray", "diffraction", "nmr", "electron")):
        return "B-factor", f"the header gives an experimental method ({method})."

    in_range = float(((values >= 0) & (values <= 100)).mean())
    if in_range > 0.999 and values.mean() > 40 and values.max() <= 100:
        return "pLDDT", ("every value sits in 0–100 with a high mean, and the "
                         "header names no experimental method, the signature of "
                         "a predicted model.")
    if values.max() > 100:
        return "B-factor", (f"values run up to {values.max():.0f}, above the "
                            "0–100 range pLDDT is confined to.")
    return "unknown", ("the values are in 0–100 but the header names no method, "
                       "so this could be either.")


def confidence_table(model: Model, chain_id: str) -> pd.DataFrame:
    """Per-residue mean of the B-factor column, whatever it means."""
    chain = model.model[chain_id]
    rows = []
    for residue in _residues(chain):
        letter = _three_to_one(residue.get_resname())
        if letter is None:
            continue
        atoms = [a.get_bfactor() for a in residue]
        rows.append({"resnum": residue.id[1], "residue": letter,
                     "value": round(float(np.mean(atoms)), 2) if atoms else None})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Disulfides, measured
# --------------------------------------------------------------------------

def disulfides_3d(model: Model,
                  max_distance: float = SS_MAX_DISTANCE) -> pd.DataFrame:
    """Every SG–SG pair close enough to be a disulfide bond, measured.

    This is the one place in the app where a cysteine pairing is known rather
    than inferred.
    """
    sulfurs = []
    for chain in model.model:
        for residue in _residues(chain):
            if residue.get_resname().upper() not in ("CYS", "CYX"):
                continue
            if "SG" in residue:
                sulfurs.append((chain.id, residue.id[1], residue["SG"]))

    rows = []
    for i, (chain_a, num_a, atom_a) in enumerate(sulfurs):
        for chain_b, num_b, atom_b in sulfurs[i + 1:]:
            distance = float(atom_a - atom_b)
            if distance <= max_distance:
                rows.append({
                    "chain_1": chain_a, "resnum_1": num_a,
                    "chain_2": chain_b, "resnum_2": num_b,
                    "distance": round(distance, 2),
                    "inter_chain": chain_a != chain_b,
                })

    bonded = {(r["chain_1"], r["resnum_1"]) for r in rows} | \
             {(r["chain_2"], r["resnum_2"]) for r in rows}
    free = [{"chain": c, "resnum": n} for c, n, _ in sulfurs
            if (c, n) not in bonded]

    frame = pd.DataFrame(rows, columns=["chain_1", "resnum_1", "chain_2",
                                        "resnum_2", "distance", "inter_chain"])
    frame.attrs["free_cysteines"] = free
    frame.attrs["n_cys_with_sg"] = len(sulfurs)
    return frame


def compare_disulfides(measured: pd.DataFrame, declared: list[dict]) -> list[str]:
    """Measured bonds against the file's own SSBOND records.

    Disagreements are reported rather than resolved: an SSBOND record that no
    measurement supports usually means the side chain was modelled apart, and a
    measured bond with no record usually means the depositor did not annotate it.
    Both are things the user should see, not things this module should decide.
    """
    if not declared:
        return []

    def key(chain_1, num_1, chain_2, num_2):
        return tuple(sorted([(chain_1, num_1), (chain_2, num_2)]))

    found = {key(r.chain_1, r.resnum_1, r.chain_2, r.resnum_2)
             for r in measured.itertuples()}
    stated = {key(b["chain_1"], b["resnum_1"], b["chain_2"], b["resnum_2"])
              for b in declared}

    notes = []
    for pair in sorted(stated - found):
        (c1, n1), (c2, n2) = pair
        notes.append(f"SSBOND {c1}{n1}–{c2}{n2} is declared in the file but the "
                     f"SG atoms are further apart than {SS_MAX_DISTANCE} Å.")
    for pair in sorted(found - stated):
        (c1, n1), (c2, n2) = pair
        notes.append(f"{c1}{n1}–{c2}{n2} is bonded by distance but carries no "
                     f"SSBOND record.")
    return notes


# --------------------------------------------------------------------------
# Solvent accessibility
# --------------------------------------------------------------------------

def _sasa_by_residue(entity) -> dict:
    """Shrake–Rupley SASA, keyed ``(chain, resnum)``."""
    _strip_solvent(entity)
    ShrakeRupley().compute(entity, level="R")
    out = {}
    for chain in entity:
        for residue in _residues(chain):
            out[(chain.id, residue.id[1])] = float(residue.sasa)
    return out


def _context_copy(model: Model):
    """A throwaway copy of the whole first model.

    SASA is computed on copies throughout, because Shrake-Rupley writes a
    ``.sasa`` attribute onto the entity and stripping waters mutates it. Doing
    either to the loaded structure would make the second analysis of the same
    upload disagree with the first.
    """
    import copy
    return copy.deepcopy(model.model)


def _isolated_copy(model: Model, chain_id: str):
    """A copy of one chain on its own, for the alone-versus-complexed comparison."""
    import copy
    chain = copy.deepcopy(model.model[chain_id])
    holder = type(model.model)(0)
    holder.add(chain)
    return holder


def residue_sasa(model: Model, chain_id: str,
                 in_context: bool = True) -> pd.DataFrame:
    """Per-residue accessibility for one chain.

    ``in_context=True`` computes it with every other chain present, which is
    what you want for a complex: it reports the surface the residue actually has.
    ``False`` computes the chain alone. The difference between the two is the
    interface.
    """
    if in_context:
        values = _sasa_by_residue(_context_copy(model))
    else:
        values = _sasa_by_residue(_isolated_copy(model, chain_id))

    rows = []
    for residue in _residues(model.model[chain_id]):
        letter = _three_to_one(residue.get_resname())
        if letter is None:
            continue
        absolute = values.get((chain_id, residue.id[1]))
        maximum = TIEN_MAX_ASA.get(letter)
        relative = (absolute / maximum) if (absolute is not None and maximum) else None
        rows.append({
            "resnum": residue.id[1],
            "residue": letter,
            "sasa": round(absolute, 1) if absolute is not None else None,
            "rel_sasa": round(relative, 3) if relative is not None else None,
            "exposure": _exposure_label(relative),
        })
    return pd.DataFrame(rows)


def _exposure_label(relative) -> str | None:
    if relative is None:
        return None
    if relative <= BURIED_MAX:
        return "buried"
    if relative >= EXPOSED_MIN:
        return "exposed"
    return "intermediate"


def exposure_table(model: Model, chain_id: str,
                   annotation: dict | None = None,
                   offset_map: dict | None = None) -> pd.DataFrame:
    """Accessibility with the CDR/framework region attached, where known."""
    frame = residue_sasa(model, chain_id, in_context=True)
    if annotation and offset_map:
        frame["region"] = [
            _region_of(offset_map.get(n), annotation) for n in frame["resnum"]]
    return frame


def _region_of(seq_index, annotation) -> str | None:
    """Which CDR or framework a 0-based sequence index falls in.

    Delegates to ``vhh_core.region_of`` rather than re-deriving the boundaries,
    so a region named here and the same region named in the analysis section
    cannot drift apart. That function is 1-based.
    """
    if seq_index is None or not annotation:
        return None
    region = vc.region_of(int(seq_index) + 1, annotation)
    return None if region == "?" else region


# --------------------------------------------------------------------------
# Interfaces
# --------------------------------------------------------------------------

def interface_table(model: Model, chain_id: str, partners: list[str],
                    min_delta: float = INTERFACE_MIN_DELTA) -> pd.DataFrame:
    """Residues of ``chain_id`` that lose accessibility when ``partners`` are present.

    ΔSASA (accessible alone, minus accessible in the complex) is the standard
    operational definition of an interface residue, and unlike a distance cutoff
    it weights a residue by how much surface it actually commits.
    """
    import copy
    alone = _sasa_by_residue(_isolated_copy(model, chain_id))

    holder = type(model.model)(0)
    for cid in [chain_id] + list(partners):
        holder.add(copy.deepcopy(model.model[cid]))
    complexed = _sasa_by_residue(holder)

    rows = []
    for residue in _residues(model.model[chain_id]):
        letter = _three_to_one(residue.get_resname())
        if letter is None:
            continue
        key = (chain_id, residue.id[1])
        free, bound = alone.get(key), complexed.get(key)
        if free is None or bound is None:
            continue
        delta = free - bound
        if delta < min_delta:
            continue
        maximum = TIEN_MAX_ASA.get(letter)
        rows.append({
            "resnum": residue.id[1],
            "residue": letter,
            "sasa_free": round(free, 1),
            "sasa_bound": round(bound, 1),
            "delta_sasa": round(delta, 1),
            "fraction_buried": round(delta / free, 2) if free > 0 else None,
            "rel_sasa_free": round(free / maximum, 3) if maximum else None,
        })
    return pd.DataFrame(rows, columns=["resnum", "residue", "sasa_free",
                                       "sasa_bound", "delta_sasa",
                                       "fraction_buried", "rel_sasa_free"])


def interface_summary(model: Model, chain_a: str, chain_b: str) -> dict:
    """Buried surface area between two chains, and how big the contact is.

    ``bsa`` is the total surface buried on both sides. The number usually quoted
    for an antibody–antigen interface is half of it, the area buried per side,
    which is reported separately rather than left to the reader to halve.
    """
    a_side = interface_table(model, chain_a, [chain_b])
    b_side = interface_table(model, chain_b, [chain_a])
    buried_a = float(a_side["delta_sasa"].sum()) if len(a_side) else 0.0
    buried_b = float(b_side["delta_sasa"].sum()) if len(b_side) else 0.0
    return {
        "chain_a": chain_a, "chain_b": chain_b,
        "bsa_total": round(buried_a + buried_b, 1),
        "bsa_per_side": round((buried_a + buried_b) / 2, 1),
        "buried_on_a": round(buried_a, 1),
        "buried_on_b": round(buried_b, 1),
        "n_residues_a": int(len(a_side)),
        "n_residues_b": int(len(b_side)),
    }


# --------------------------------------------------------------------------
# Surface hydrophobic patches
# --------------------------------------------------------------------------

def _side_chain_coords(residue) -> np.ndarray | None:
    """Heavy-atom coordinates of the side chain, backbone and hydrogens excluded.

    Glycine has no side chain, so it falls back to CA, it is not in
    ``PATCH_RESIDUES`` anyway, but the fallback keeps the function total.
    """
    atoms = [a for a in residue
             if a.get_id() not in ("N", "CA", "C", "O", "OXT")
             and (a.element or "C") != "H"]
    if not atoms:
        atoms = [a for a in residue if a.get_id() == "CA"]
    if not atoms:
        return None
    return np.array([a.get_coord() for a in atoms], dtype=float)


def _min_side_chain_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Closest approach between two sets of heavy atoms."""
    return float(np.min(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)))


def surface_patches(model: Model, chain_id: str,
                    min_rel_sasa: float = PATCH_MIN_REL,
                    contact: float = PATCH_CONTACT,
                    min_size: int = PATCH_MIN_SIZE) -> pd.DataFrame:
    """Groups of exposed hydrophobic residues that are adjacent in space.

    **A transparent composite, not an implementation of SAP or CamSol.** Those
    methods are parameterised against experimental data this module does not
    have; this one simply asks which strongly hydrophobic residues are both
    solvent-exposed and clustered together, which is the structural feature they
    are all ultimately looking at. Read it as "here is a hydrophobic area worth
    looking at", never as a predicted aggregation rate.
    """
    exposure = residue_sasa(model, chain_id, in_context=True)
    exposed = {int(r.resnum): r for r in exposure.itertuples()
               if r.residue in PATCH_RESIDUES
               and r.rel_sasa is not None and r.rel_sasa >= min_rel_sasa}
    if not exposed:
        return pd.DataFrame(columns=["patch", "n_residues", "residues",
                                     "exposed_area"])

    coords = {}
    for residue in _residues(model.model[chain_id]):
        if residue.id[1] in exposed:
            points = _side_chain_coords(residue)
            if points is not None:
                coords[residue.id[1]] = points

    numbers = sorted(coords)
    parent = {n: n for n in numbers}

    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    for i, a in enumerate(numbers):
        for b in numbers[i + 1:]:
            if _min_side_chain_distance(coords[a], coords[b]) <= contact:
                parent[find(a)] = find(b)

    groups: dict[int, list[int]] = {}
    for n in numbers:
        groups.setdefault(find(n), []).append(n)

    rows = []
    for members in sorted(groups.values(), key=len, reverse=True):
        if len(members) < min_size:
            continue
        rows.append({
            "n_residues": len(members),
            "residues": ", ".join(f"{exposed[n].residue}{n}"
                                  for n in sorted(members)),
            "exposed_area": round(sum(exposed[n].sasa for n in members), 1),
        })
    frame = pd.DataFrame(rows, columns=["n_residues", "residues", "exposed_area"])
    if len(frame):
        frame = frame.sort_values("exposed_area", ascending=False)
        frame.insert(0, "patch", range(1, len(frame) + 1))
    else:
        frame.insert(0, "patch", [])
    return frame


# --------------------------------------------------------------------------
# Joining a structure to a sequence the rest of the app knows about
# --------------------------------------------------------------------------

def map_structure_to_sequence(structure_seq: str,
                              target_seq: str) -> tuple[dict, float]:
    """Map structure residue indices onto positions of a reference sequence.

    A model is routinely missing termini and loops, and its residue numbering may
    be sequential, IMGT, Kabat or the depositor's own. Rather than trusting any
    of that, the structure's own sequence is aligned to the sequence the rest of
    the app is analysing. Returns ``({structure_index: target_index}, identity)``.
    """
    if not structure_seq or not target_seq:
        return {}, 0.0
    aligner = vc.protein_aligner()
    alignment = aligner.align(structure_seq, target_seq)[0]
    mapping = {}
    matches = 0
    for (s_start, s_end), (t_start, t_end) in zip(*alignment.aligned):
        for offset in range(s_end - s_start):
            s_index, t_index = s_start + offset, t_start + offset
            mapping[s_index] = t_index
            if structure_seq[s_index] == target_seq[t_index]:
                matches += 1
    identity = matches / max(1, min(len(structure_seq), len(target_seq)))
    return mapping, identity


def annotate_liabilities(liabilities: pd.DataFrame, model: Model, chain_id: str,
                         structure_seq: str, target_seq: str) -> pd.DataFrame:
    """Attach measured solvent accessibility to each sequence-derived liability.

    This is the join that makes the structure worth uploading: a deamidation
    motif buried in the core and one sitting on the CDR3 apex are the same row in
    the sequence-only table and very different problems in practice.
    """
    frame = liabilities.copy()
    if not len(frame):
        for column in ("rel_sasa", "exposure"):
            frame[column] = pd.Series(dtype="object")
        return frame

    residues = [r for r in _residues(model.model[chain_id])
                if _three_to_one(r.get_resname()) is not None]

    # ``structure_seq`` must be this chain's own sequence, not the full-length
    # construct. Passing the wrong one does not fail: the mapping comes back
    # looking plausible and every position after a gap is quietly off. Since the
    # two arguments are both sequences and easy to swap, the mismatch is refused
    # rather than discovered later in a table nobody can check by eye.
    if len(structure_seq) != len(residues):
        raise StructureError(
            f"structure_seq has {len(structure_seq)} residues but chain "
            f"{chain_id} has {len(residues)} - pass the chain's own sequence "
            f"from chain_sequence(), not the full-length construct.")

    mapping, _ = map_structure_to_sequence(structure_seq, target_seq)
    target_to_structure = {t: s for s, t in mapping.items()}
    exposure = residue_sasa(model, chain_id, in_context=True)
    by_resnum = {int(r.resnum): r for r in exposure.itertuples()}

    rel, label, resnums = [], [], []
    for position in frame.get("position", []):
        structure_index = target_to_structure.get(int(position) - 1)
        row = resnum = None
        if structure_index is not None and structure_index < len(residues):
            resnum = residues[structure_index].id[1]
            row = by_resnum.get(resnum)
        # The residue number in the *model*, which is what the 3D viewer and the
        # file itself speak. It is reported even when the residue has no
        # accessibility, so a row can always be located in the structure.
        resnums.append(int(resnum) if resnum is not None else None)
        rel.append(row.rel_sasa if row is not None else None)
        label.append(row.exposure if row is not None else "not modelled")
    frame["resnum"] = resnums
    frame["rel_sasa"] = rel
    frame["exposure"] = label
    return frame


def resnum_map(model: Model, chain_id: str, structure_seq: str,
               target_seq: str) -> tuple[dict, float]:
    """``({residue number: 0-based index in target_seq}, identity)``.

    :func:`map_structure_to_sequence` is keyed by the structure's *sequential*
    index, because that is what aligning two strings gives you. Everything that
    consumes the mapping (the figures, the region column, the 3D viewer) speaks
    **residue numbers** instead, since that is what the file and the viewer use
    and the two differ as soon as the chain does not start at 1 or has a gap.
    Converting once, here, is what keeps that distinction from leaking.
    """
    mapping, identity = map_structure_to_sequence(structure_seq, target_seq)
    residues = [r for r in _residues(model.model[chain_id])
                if _three_to_one(r.get_resname()) is not None]
    out = {}
    for structure_index, target_index in mapping.items():
        if structure_index < len(residues):
            out[int(residues[structure_index].id[1])] = target_index
    return out, identity


def identity_resnum_map(model: Model, chain_id: str) -> dict:
    """The mapping for "annotate against the structure's own sequence".

    The chain is then its own reference, so residue number *n* is whatever
    sequential position it occupies, not *n* itself.
    """
    residues = [r for r in _residues(model.model[chain_id])
                if _three_to_one(r.get_resname()) is not None]
    return {int(r.id[1]): i for i, r in enumerate(residues)}


def residue_regions(resnums, annotation: dict | None,
                    resnum_to_index: dict | None) -> list:
    """The CDR/framework region of each model residue number.

    ``resnum_to_index`` is what :func:`resnum_map` or
    :func:`identity_resnum_map` returns.
    """
    resnums = list(resnums)
    if not annotation or not resnum_to_index:
        return [None] * len(resnums)
    return [_region_of(resnum_to_index.get(int(n)), annotation) for n in resnums]


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def _cdr_bands(ax, annotation, offset_map):
    """Shade the CDRs on an axis whose x is the model's residue numbers.

    ``offset_map`` must be keyed by residue number (:func:`resnum_map`), not by
    sequential index: the bands are drawn in x, and x is residue numbers.
    """
    if not (annotation and offset_map):
        return
    reverse = {}
    for resnum, target_index in offset_map.items():
        reverse.setdefault(target_index, resnum)
    for region, colour in vc.REGION_COLORS.items():
        if not region.startswith("CDR"):
            continue
        info = annotation.get(region) or {}
        start, end = info.get("start"), info.get("end")
        if not (start and end):
            continue
        covered = [reverse[i] for i in range(start - 1, end) if i in reverse]
        if not covered:
            continue
        ax.axvspan(min(covered) - 0.5, max(covered) + 0.5, color=colour,
                   alpha=0.22, lw=0, zorder=0)
        ax.text((min(covered) + max(covered)) / 2, 1.02, region,
                transform=ax.get_xaxis_transform(), ha="center",
                fontsize=7.5, color="0.3")


def plot_confidence(table: pd.DataFrame, kind: str = "pLDDT",
                    title: str | None = None):
    """pLDDT (or B-factor) along the chain.

    For pLDDT the usual confidence bands are drawn, because the number is close
    to meaningless without them: 90 and 60 are both "high-ish" until you know
    that 70 is the line below which the backbone should not be trusted.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 3.0))
    x = table["resnum"].to_numpy()
    y = table["value"].to_numpy(dtype=float)

    if kind == "pLDDT":
        for low, high, colour, label in ((90, 100, "#0053d6", "very high"),
                                        (70, 90, "#65cbf3", "confident"),
                                        (50, 70, "#ffdb13", "low"),
                                        (0, 50, "#ff7d45", "very low")):
            ax.axhspan(low, high, color=colour, alpha=0.18, lw=0)
            ax.text(1.004, (low + high) / 2 / 100, label, transform=ax.transAxes,
                    va="center", fontsize=7, color="0.35")
        ax.set_ylim(0, 100)
        ax.set_ylabel("pLDDT")
    else:
        ax.set_ylabel("B-factor (Å²)")

    ax.plot(x, y, color="#1f3b73", lw=1.4, zorder=3)
    ax.set_xlabel("Residue number in the model")
    ax.set_title(title or f"{kind} along the chain", fontsize=10)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


def plot_exposure(table: pd.DataFrame, annotation: dict | None = None,
                  offset_map: dict | None = None,
                  title: str = "Relative solvent accessibility"):
    """Relative accessibility per residue, with the buried/exposed bands marked."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 3.2))
    x = table["resnum"].to_numpy()
    y = table["rel_sasa"].astype(float).to_numpy()

    _cdr_bands(ax, annotation, offset_map)
    ax.axhspan(0, BURIED_MAX, color="0.6", alpha=0.18, lw=0)
    ax.axhline(EXPOSED_MIN, color="crimson", lw=0.9, ls="--")
    ax.fill_between(x, 0, y, color="#1f3b73", alpha=0.35, zorder=2)
    ax.plot(x, y, color="#1f3b73", lw=1.2, zorder=3)
    ax.text(0.004, BURIED_MAX / 2, "buried", transform=ax.get_yaxis_transform(),
            fontsize=7, color="0.35", va="center")
    ax.set_xlabel("Residue number in the model")
    ax.set_ylabel("rel. accessibility")
    ax.set_title(title, fontsize=10)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    ax.set_ylim(0, max(1.0, float(np.nanmax(y)) * 1.05))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


def plot_interface(interface: pd.DataFrame, all_residues: pd.DataFrame,
                   annotation: dict | None = None, offset_map: dict | None = None,
                   title: str = "Surface buried on binding"):
    """ΔSASA per residue, so the epitope footprint is visible along the chain."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 3.2))
    _cdr_bands(ax, annotation, offset_map)
    delta = dict(zip(interface["resnum"], interface["delta_sasa"]))
    x = all_residues["resnum"].to_numpy()
    y = np.array([float(delta.get(int(n), 0.0)) for n in x])
    ax.bar(x, y, width=0.85, color="#b5451b", zorder=3)
    ax.set_xlabel("Residue number in the model")
    ax.set_ylabel("ΔSASA (Å²)")
    ax.set_title(title, fontsize=10)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# 3D viewer
# --------------------------------------------------------------------------
# The viewer is the one thing here that genuinely needs a picture: a list of
# exposed hydrophobic residues is a list, and the same residues shown as one
# contiguous yellow area on the surface is an argument. It is also the one part
# that depends on a CDN (3Dmol.js), so every table and figure above works
# without it and the viewer is allowed to be the thing that fails.

PLDDT_COLOURS = [(90, "#0053d6"), (70, "#65cbf3"), (50, "#ffdb13"), (0, "#ff7d45")]
EXPOSURE_COLOURS = {"buried": "#4a4a4a", "intermediate": "#8fb3d9",
                    "exposed": "#d94f1e"}
CHAIN_PALETTE = ["#1f3b73", "#b5451b", "#2e7d5b", "#7b4fa0", "#b38600",
                 "#1d7a8c", "#8c1d4a"]


def _plddt_colour(value) -> str:
    for threshold, colour in PLDDT_COLOURS:
        if value is not None and value >= threshold:
            return colour
    return "#ff7d45"


def viewer_html(model: Model, raw_text: str, colour_by: str = "chain",
                chain_id: str | None = None,
                exposure: pd.DataFrame | None = None,
                confidence: pd.DataFrame | None = None,
                highlight: list[int] | None = None,
                surface: bool = False,
                width: int = 760, height: int = 520) -> str:
    """An interactive 3Dmol.js view, returned as HTML to embed.

    ``raw_text`` is the original file content rather than a re-serialised copy,
    so what is displayed is exactly what was uploaded.
    """
    import py3Dmol

    fmt = "cif" if model.fmt == "mmCIF" else "pdb"
    view = py3Dmol.view(width=width, height=height)
    view.addModel(raw_text, fmt)
    view.setStyle({}, {"cartoon": {"color": "#c9ccd1"}})

    chains = [c.id for c in model.model]

    if colour_by == "chain":
        for index, cid in enumerate(chains):
            view.setStyle({"chain": cid},
                          {"cartoon": {"color": CHAIN_PALETTE[index % len(CHAIN_PALETTE)]}})
    elif colour_by == "spectrum":
        view.setStyle({}, {"cartoon": {"color": "spectrum"}})
    elif colour_by == "confidence" and confidence is not None and chain_id:
        groups: dict[str, list[int]] = {}
        for row in confidence.itertuples():
            groups.setdefault(_plddt_colour(row.value), []).append(int(row.resnum))
        for colour, residues in groups.items():
            view.setStyle({"chain": chain_id, "resi": residues},
                          {"cartoon": {"color": colour}})
    elif colour_by == "exposure" and exposure is not None and chain_id:
        groups = {}
        for row in exposure.itertuples():
            if row.exposure:
                groups.setdefault(EXPOSURE_COLOURS[row.exposure], []).append(
                    int(row.resnum))
        for colour, residues in groups.items():
            view.setStyle({"chain": chain_id, "resi": residues},
                          {"cartoon": {"color": colour}})

    if highlight and chain_id:
        view.addStyle({"chain": chain_id, "resi": [int(n) for n in highlight]},
                      {"stick": {"colorscheme": "yellowCarbon", "radius": 0.22}})

    if surface:
        view.addSurface("VDW", {"opacity": 0.72, "color": "#dfe3e8"},
                        {"chain": chain_id} if chain_id else {})

    view.zoomTo()
    return view._make_html()
