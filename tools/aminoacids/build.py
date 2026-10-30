"""
Draw the twenty amino acids once, so the app does not need RDKit.

RDKit is 66 MB installed and is needed for exactly one thing here: twenty
pictures that never change. So they are drawn here, checked, and committed to
`data/aminoacids/`; the app ships none of this.

Three things are deliberate:

* **The form drawn is the one that exists at the stated pH**, zwitterion and
  ionised side chains, not the textbook neutral form. The pH is written on
  the screen that shows them, and the neutral SMILES is kept in the metadata
  so the molecular weight below matches the residue tables in `vhh_core.py`.
* **Every picture is turned the same way**: nitrogen left, carboxylate right,
  side chain up. RDKit's layout is kept and only rotated or mirrored as a
  rigid body, so the geometry stays RDKit's and only the orientation is ours.
* **The carbon skeleton is `currentColor`**, so it follows the light or dark
  theme. The heteroatoms keep CPK-like colours, chosen a shade lighter than
  RDKit's defaults, which are pure blue and pure red and go illegible on a
  dark background.

Run:  python build.py      (writes svg/, meta.json, and reports any problem)
"""

from __future__ import annotations

import json
import math
import os
import re
import sys

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Descriptors, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from Bio.Data.IUPACData import protein_letters_3to1, protein_weights

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smiles import EXPECTED_CHARGE, NEUTRAL, PH74      # noqa: E402

PH = 7.4
rdDepictor.SetPreferCoordGen(True)

# The backbone, matched on every amino acid including proline (secondary
# nitrogen, in a ring) and glycine (no side chain).
BACKBONE = Chem.MolFromSmarts("[NX4,NX3][CX4H1,CX4H2][CX3](=O)[O-]")

# A shade lighter than RDKit's defaults, which are pure blue and pure red:
# legible on white and on Streamlit's dark background alike. Carbon stays
# black so that themeable() can turn it into currentColor.
PALETTE = {7: (0.23, 0.51, 0.96),      # N
           8: (0.90, 0.22, 0.21),      # O
           16: (0.79, 0.60, 0.05),     # S
           6: (0.0, 0.0, 0.0)}         # C, becomes currentColor

BLACK = ("#000000", "#000")


def oriented(smiles: str):
    """Lay the molecule out, then turn it so every picture reads the same way.

    Returns (mol, matched). A rigid rotation and an optional mirror: the
    molecule is unchanged, only its orientation is, and the wedge bonds are
    recomputed from the stereo tags afterwards by PrepareAndDrawMolecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    match = mol.GetSubstructMatch(BACKBONE)
    if not match:
        return mol, False
    n_idx, ca_idx, c_idx = match[0], match[1], match[2]
    conf = mol.GetConformer()
    pos = conf.GetAtomPosition

    angle = -math.atan2(pos(c_idx).y - pos(n_idx).y, pos(c_idx).x - pos(n_idx).x)
    cos, sin = math.cos(angle), math.sin(angle)
    for i in range(mol.GetNumAtoms()):
        p = pos(i)
        conf.SetAtomPosition(i, Geometry.Point3D(p.x * cos - p.y * sin,
                                                 p.x * sin + p.y * cos, 0.0))
    side = [nb.GetIdx() for nb in mol.GetAtomWithIdx(ca_idx).GetNeighbors()
            if nb.GetIdx() not in (n_idx, c_idx)]
    if side and pos(side[0]).y < pos(ca_idx).y:
        for i in range(mol.GetNumAtoms()):
            p = pos(i)
            conf.SetAtomPosition(i, Geometry.Point3D(p.x, -p.y, 0.0))
    return mol, True


def themeable(svg: str) -> str:
    """Make the carbon skeleton follow the page's text colour, and scale.

    RDKit writes black in two shapes: ``stroke:#000000`` inside a style
    attribute for the bonds, and ``fill='#000000'`` as an attribute on the
    glyph paths of a label. Only carbon, hydrogen and the stereo annotations
    are black, so swapping the black leaves the heteroatom colours alone.
    """
    for black in BLACK:
        for prop in ("stroke", "fill"):
            svg = svg.replace(f"{prop}:{black}", f"{prop}:currentColor")
            svg = svg.replace(f"{prop}='{black}'", f"{prop}='currentColor'")
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"width='\d+px' height='\d+px'",
                 "width='100%' height='100%' preserveAspectRatio='xMidYMid meet'",
                 svg)
    return svg


def draw(mol) -> str:
    drawer = rdMolDraw2D.MolDraw2DSVG(320, 240)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 1.6
    opts.padding = 0.10
    opts.minFontSize = 11
    opts.maxFontSize = 16
    for number, colour in PALETTE.items():
        opts.updateAtomPalette({number: colour})
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return themeable(drawer.GetDrawingText())


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "svg")
    os.makedirs(out, exist_ok=True)

    meta, problems = {}, []
    for code, smiles in PH74.items():
        mol, matched = oriented(smiles)
        if not matched:
            problems.append(f"{code}: backbone not matched, orientation is "
                            f"RDKit's own")
        with open(os.path.join(out, f"{code}.svg"), "w") as handle:
            handle.write(draw(mol))

        charge = Chem.GetFormalCharge(mol)
        if charge != EXPECTED_CHARGE[code]:
            problems.append(f"{code}: charge {charge:+d} at pH {PH}, "
                            f"expected {EXPECTED_CHARGE[code]:+d}")

        one = protein_letters_3to1[code]
        neutral = Chem.MolFromSmiles(NEUTRAL[code])
        weight = Descriptors.MolWt(neutral)
        # Biopython keys this by one-letter code, and gives the
        # free amino acid, which is why NEUTRAL is weighed here
        # and not the zwitterion that is drawn.
        reference = protein_weights[one]
        if abs(weight - reference) > 0.01:
            problems.append(f"{code}: {weight:.3f} Da against Biopython's "
                            f"{reference:.3f}")

        meta[code] = {
            "one": one,
            "charge": charge,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mw_free_acid": round(weight, 3),
            "smiles_ph": smiles,
            "smiles_neutral": NEUTRAL[code],
        }

    with open(os.path.join(here, "meta.json"), "w") as handle:
        json.dump({"ph": PH, "drawn_with": f"RDKit {Chem.rdBase.rdkitVersion}",
                   "aminoacids": meta}, handle, indent=1)

    total = sum(os.path.getsize(os.path.join(out, f"{c}.svg")) for c in meta)
    print(f"{len(meta)} structures, {total // 1024} kB, pH {PH}")
    print("problems:", "; ".join(problems) if problems else "none")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
