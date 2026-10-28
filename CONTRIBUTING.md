# Contributing

This is a small scientific tool with one maintainer. That shapes what is useful
to send and how quickly it is answered, so this page says both plainly rather
than implying a project larger than it is.

Everything happens in the
[issue tracker](https://github.com/TyrosinDuck/vhh-analyzer/issues).

---

## ⚠️ Before you paste a sequence anywhere

**GitHub issues are public, and so is `data/library.csv`.** Anything you put in
either is readable by anyone, for good.

Do not paste unpublished sequences, anything under an NDA, or anything from a
collaborator who has not agreed. A bug can almost always be reproduced with a
published VHH, there are 22 of them in the library, and `cAbLys3` or
`caplacizumab` will usually do. If a defect genuinely only appears on a
sequence
you cannot share, say so in the issue and describe the shape of it: its length,
whether the CDR3 is unusual, whether it has extra cysteines. That is often
enough to find the problem.

If you need to analyse confidential sequences at all, run the app locally
rather
than using the public deployment, the README explains how. Nothing is stored
server-side either way, but local is local.

---

## Reporting a bug

Open an issue. What makes one easy to fix:

- **The version**, from the footer at the bottom of any page, for example
  `v1.0.0 (build 20261028)`. Without it there is no way to know whether the
  defect still exists.
- **Which screen and which panel**, since the same number can be reached from
  several places.
- **What you expected and what happened.** "The CDR3 is wrong" is hard to act
  on; "the CDR3 is reported as one residue shorter than ANARCI gives" is a
  defect that can be chased.
- **A sequence that reproduces it**, subject to the warning above.
- **A screenshot**, if it is about the display rather than a number.

Not every report has to be complete. A half-described problem is still worth
sending; it is the ones nobody mentions that stay broken.

## A wrong number is the most valuable report there is

Most of this app computes things that look plausible whether or not they are
right, which is precisely what makes a quiet error dangerous. If a value
disagrees with ANARCI, with a published figure, with a gel, or with what you
measured at the bench, that report is worth more than a feature request, and it
will be treated that way.

Please include what you compared against. "The pI disagrees with my IEF" is the
beginning of a useful thread; the theoretical value and the measured one are
allowed to differ, and knowing by how much is the point.

## Asking for help

Open an issue and label it as a question. There is no mailing list and no chat.

Read [`MANUAL.md`](MANUAL.md) first if you have not: it documents every screen,
what each number means, and (more usefully) how far to trust it. Most questions
about what a value means are answered there, including the ones about what the
app deliberately does not claim.

## Suggesting a change

Open an issue before writing code, so that effort is not spent on something
that
will be declined. Some things are out of scope on purpose:

- **Predictions parameterised against data this project does not have.** The
  aggregation score is a stated composite, not an implementation of CamSol or
  SAP, and will not be relabelled as one.
- **A single composite "developability score".** Weighting a deamidation site
  against three divergent framework positions depends on what the VHH is for,
  and a single number hides that choice behind a decimal point.
- **Anything requiring server-side storage of sequences.** The app stores
  nothing, and that is a design constraint rather than an unfinished feature.

## Contributing code

1. Fork the repository and work on a branch.
2. Install the dependencies: `pip install -r requirements.txt -r requirements-dev.txt`.
   Run the tests from the repository root.
3. Make the change.
4. Run the tests: `pytest`. They take about half a minute. Everything must pass.
5. **Add a test.** Almost every test in `tests/` is a defect that shipped once;
   a fix without a test invites the same defect back.
6. Open a pull request describing what changed and why.

Two conventions hold throughout, and the code reads strangely if you miss them:

- **A screen draws; it does not compute.** Anything cached belongs in
  `compute.py`, not in a `sections/` module.
- **No screen reads another screen's globals.** Every screen is a `render()`
  function. This is not tidiness: in an earlier version the Structure screen's
  `kinds` assignment silently overwrote the analysis screen's variable of the
  same name.

When you change a module, bump its `__version__` and its entry in
`REQUIRED_MODULES` in `app.py`. A test checks that the two agree.

## Contributing science rather than code

As valuable here, and often more so:

- **A correction to a threshold or a constant**, with a reference. Several
  numbers in this app are stated as thresholds rather than validated
parameters, the surface-patch exposure cut-off was set from a single structure,
and says
  so in `structure.py`. A better-founded value, with the evidence, is welcome.
- **A reference that should be cited** where the code currently cites nothing.
- **A published sequence for the library**, with its PDB entry or accession and
  a note on what was trimmed. Published only, see the warning at the top.
- **A case where the app is confidently wrong**, which is the single most useful
  thing anyone can send.

## Licence

Contributions are accepted under the [MIT licence](LICENSE), the licence this
project is released under.
