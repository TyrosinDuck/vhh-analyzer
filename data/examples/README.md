# Example inputs

What the app loads when you press *Load an example* on a screen that needs a
file you may not have to hand. `demo.py` is the module that reads this folder
and builds the rest.

## `1mel.pdb.gz`

PDB entry **1MEL**, the camel VHH cAbLys3 in complex with hen egg-white
lysozyme, as downloaded from `files.rcsb.org`, gzipped and otherwise
unmodified. Two copies of the complex are in the asymmetric unit: chains A and
B are the VHH, chains L and M the lysozyme.

It is here rather than fetched because the app deliberately does not download
structures on a user's behalf, and because an example that depends on a remote
service is an example that stops working.

> Desmyter A, Transue TR, Ghahroudi MA, Thi MH, Poortmans F, Hamers R,
> Muyldermans S, Wyns L (1996) Crystal structure of a camel single-domain VH
> antibody fragment in complex with lysozyme. *Nat Struct Biol* 3:803-811.
> doi:[10.1038/nsb0996-803](https://doi.org/10.1038/nsb0996-803)

Files in the PDB archive are released under the
[CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/);
wwPDB asks that users credit the depositors, which the citation above does.
The licence covers the file; it says nothing about the protein, and nothing
here is a claim about either.

Why this entry. It is a real crystal structure, so the B-factor column holds a
B-factor and not a pLDDT, which is the case the Structure screen has to guess
and the one an AlphaFold model never shows. It is a complex, so the interface
view has two chains to measure between. Its VHH carries the extra CDR1-CDR3
disulfide, so the measured disulfides and the sequence-based prediction have
something to disagree about. And its chains are already in the bundled
sequence library, as `cAbLys3` and `hen-lysozyme`, so *Annotate against* can
be pointed at the full-length sequence without uploading anything else.

Chain A is modelled over 132 of the 148 residues the entry declares: the first
residue and the C-terminal `GRYPYDVPDYGSGRA` are in the SEQRES records and not
in the coordinates. That is not a defect, it is the ordinary case the
*Amino-acid sequence of each chain* panel exists to show.

## The Sanger chromatograms

Not here: there are no `.ab1` files in this repository. The example reads are
built in memory by `demo.sanger_example_files()` from the bundled pMES3 map and
the bundled sequence library, so they cannot drift away from either, and they
weigh nothing in the repository.

**They are simulated.** The sequences are real VHHs in a real vector, the
traces and the quality values are invented. Each file says so in its own
`SMPL` and `CMNT` tags, and every screen that offers them says so on the page.
