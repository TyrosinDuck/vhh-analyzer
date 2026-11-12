# A worked example

Three published VHH sequences, one set of settings, and the numbers this code
produces from them. Run in one command:

```
python3 tools/make_examples.py --check
```

It prints `3 example sequences match the baseline`, or says which sequence
differs. `tests/test_examples.py` runs the same comparison as part of
`python3 -m pytest`.

## What is here

| File | What it is |
| --- | --- |
| `manifest.json` | Which library entries the example uses, and the settings |
| `sequences.fasta` | Those sequences, written out so they can be read without running anything |
| `expected.json` | The descriptors, annotation and classification the code produced |

The sequences are not stored here. `manifest.json` names three entries of
`data/library.csv`, the library the app ships, and both other files are
generated from it, so a correction to a sequence in the library reaches this
example instead of leaving it describing the old one. The committed
`sequences.fasta` is checked against the library by the test suite.

The annotator is the regex one rather than ANARCI, deliberately. ANARCI is
installed beside the app rather than by it, so a baseline annotated with it
would pass or fail according to whether HMMER happens to be on the machine.

## What this does, and what it does not

**It detects change.** The expected values were produced by this code, so they
are a baseline and not an independent reference. What they catch is a number
that moved: a refactor that changed a descriptor, and, just as important, a
dependency upgrade that changed one. Biopython computes the molar mass, the
isoelectric point and the extinction coefficients; `requirements.txt` pins it,
and the day that pin moves this is what says whether the numbers moved with
it. Nothing else in the test suite would notice, because every other test
checks that the code does what the code says, and a changed constant is
perfectly consistent with that.

**It does not establish that a number is right.** For that, a value has to be
compared against something outside this repository.

## What has been checked against a second source

**The molar mass: yes, and it agrees to 14 ppm.** A protein's average mass is
the sum of its residues' free-acid masses less one water per peptide bond. The
app's figure comes from Biopython's residue weight table. The test suite
recomputes it from the per-residue formulas in `data/aminoacids/meta.json` and
the CIAAW 2021 standard atomic weights in `buffers.py`, through the same
formula parser the buffer recipes use. Those two routes share no constant, and
on these three sequences they differ by at most 0.19 Da on about 13.9 kDa,
which is 14 parts per million. The residual is the rounding in the published
free-acid masses: the same comparison on the twenty individual amino acids is
exact for eighteen of them, and differs by 0.007 Da for methionine and
cysteine, which is a revision of sulfur's standard atomic weight.

**The isoelectric point and the extinction coefficients: not checked here.**
The pI comes from an iterative solve over pKa values and the extinction
coefficients from the Pace tabulation of tryptophan, tyrosine and cystine
absorbances. This repository holds no second source for either, so there is
nothing here to check them against, and the test suite does not pretend
otherwise. Checking them means pasting a sequence into Expasy ProtParam and
comparing by hand. That has not been done as part of this example. The app
says on every screen that both are theoretical and should be confirmed
experimentally, which is a different statement and not a substitute for this
one.

**The CDR boundaries: no, and they are not the kind of thing that has a right
answer.** IMGT and the regex anchors disagree about CDR1 and CDR2 by
construction, and the app reports which method it used for exactly that
reason. What the baseline pins is that the regex annotator still finds the
same boundaries it found before.

## Updating the baseline

```
python3 tools/make_examples.py
```

Only after finding out why a number moved. If it moved because a pinned
dependency changed, that belongs in `CHANGELOG.md`: a result that changed
under an upgrade is something a user who cited a version DOI needs to be able
to find.
