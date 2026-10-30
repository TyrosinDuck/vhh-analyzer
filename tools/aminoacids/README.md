# The twenty amino acid structures

`data/aminoacids/` holds twenty SVG files and a `meta.json`, drawn once by
`build.py` and committed. The app reads the SVGs and never imports RDKit:
RDKit is 66 MB installed and would be needed for exactly one thing, twenty
pictures that never change.

Nothing here is imported by the app. It is kept so that the pictures can be
rebuilt, and so that anyone can check how they were made.

## Rebuilding them

    pip install rdkit
    python tools/aminoacids/build.py
    cp tools/aminoacids/svg/*.svg tools/aminoacids/meta.json data/aminoacids/

`build.py` exits non-zero and prints what went wrong if any structure fails a
check, so the rebuild is not silent about a mistake.

## What is drawn, and what is checked

The form drawn is the one that exists **at pH 7.4**: the zwitterion, with
aspartate and glutamate deprotonated, lysine and arginine protonated, and
histidine neutral (its side chain pKa is near 6.0, so at 7.4 the neutral
imidazole dominates). The pH is written on the screen that shows them, because
a drawing of a charged species without its pH is a drawing of nothing in
particular.

`smiles.py` holds two SMILES per amino acid, the neutral free acid and the
pH 7.4 species, and the charge each one should carry. `build.py` checks all
three things before it writes anything:

- the formal charge of the drawn species against `EXPECTED_CHARGE`;
- the molecular weight of the *neutral* form against Biopython's
  `protein_weights`, to 0.01 Da, so the pictures and the residue masses used
  everywhere else in the app cannot be describing different molecules;
- that the backbone was matched, which is what the orientation depends on.

The CIP descriptors were checked separately when the SMILES were written: all
L, which is S at the alpha carbon for every amino acid except cysteine, where
the sulfur outranks the carboxyl and the same spatial arrangement is called R.
Threonine is 2S,3R and isoleucine 2S,3S.

## Why the pictures look the way they do

**Every one is turned the same way**, nitrogen left, carboxylate right, side
chain up. RDKit's layout is good and its orientation is arbitrary; a reader
comparing two side chains should not have to find the backbone again in each
picture. So the layout is kept and only rotated, and mirrored if the side
chain came out underneath, as a rigid body. The wedge bonds are recomputed
from the stereo tags afterwards, not carried over from the old coordinates.

**The carbon skeleton is `currentColor`**, so it follows the light or dark
theme instead of being black on both. The heteroatoms keep CPK-like colours,
a shade lighter than RDKit's defaults, which are pure blue and pure red and
disappear into a dark background.

**The SVGs have no fixed width or height**, so the page decides how large they
are.
