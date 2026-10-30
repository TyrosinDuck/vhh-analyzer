# VHH tool box, user manual

A browser-based workbench for single-domain antibodies. Nothing to install, no
code
to write: you paste sequences or drop files, and it gives you characterization,
sequencing read-back and cloning plans.

This manual describes every screen and, just as importantly, what each number
means
and how far you should trust it.

---

## 1. What it does, and what it does not

**It does** compute descriptors from an amino-acid sequence, annotate CDRs and
frameworks, read VHHs out of Sanger chromatograms, plan a restriction cloning
between two plasmids, and (if you upload a model) measure on a structure what
the rest of the app can only infer from sequence.

**It is built around VHHs, but it accepts any protein.** Each submitted sequence is
classified. Conventional VH domains get the same full analysis as VHHs,
labelled
as VH. For everything else, the analyses that are not antibody-specific (mass,
pI, extinction coefficient, charge, hydrophobicity, secondary-structure
propensity, aggregation, chemical liabilities, codon-optimised DNA, fusion
constructs) run on all of them.
Only the VHH-specific tabs are skipped for the rest. A mixed submission works:
see
§3.1c.

**It does not** predict structure, affinity, expression yield or immunogenicity.
Apart from the structure section (which measures a model **you** supply, and
never produces one) every value it reports is derived from sequence alone.
Treat the
output as triage: it tells you where to look, not what will happen at the
bench.

Three things are worth stating plainly before you rely on any number:

- The physico-chemical values are **theoretical**. Confirm the isoelectric point by
  IEF or cIEF and the extinction coefficient against a measured A280 before
using
  them quantitatively.
- The secondary-structure output in the analysis section is a **composition-based
  propensity**, not a structure prediction. For real secondary structure, run
DSSP
  on a model. Uploading that model to §4 gives you the accessibility, the real
  disulfides and the interface, but not a DSSP assignment.
- The cloning section is a **design check**. It verifies that your construct would
  encode what you intend; it says nothing about ligation efficiency or whether
the
  protein will express.

---

## 2. Getting around

The app opens on a home page with seven sections. The sidebar carries the same
seven entries plus a **Home** button, and, below the divider, the **Manual**
you are reading and a **Reference** screen, available from every screen.

| Section | Use it for |
|---|---|
| 🔬 **Sequence analysis** | characterizing VHH sequences you already have |
| 🧊 **Structure** | a PDB or AlphaFold model: what is measured rather than guessed |
| ⚗️ **Bench** | the tube: concentration from an A280, mass identity, how much labelling reagent to add |
| 🧪 **Sanger sequencing** | pulling VHH sequences out of `.ab1` files |
| 🔧 **Gene diagnosis** | checking coding sequences that already exist, and lowering GC where it runs high |
| 📚 **Sequence library** | the sequences you keep coming back to, ready to send to any of the others |
| 🧬 **In-silico cloning** | moving a VHH from one vector into another |
| 📐 **Reference** | the genetic code and the twenty amino acids, looked up rather than worked out |

They are connected rather than separate. The sequencing section writes into the
analysis input; the analysis input is what the cloning section offers as an
insert;
the analysis section's **DNA** tab hands its back-translated sequences to gene
diagnosis in one click; the library feeds proteins to the analysis and DNA to
gene
diagnosis; and the structure section reads the analysis input, so a model can
be
annotated against the full-length sequence it came from. Each section still
works
entirely on its own.

Inside the analysis section, the views are grouped into six families chosen
with a
selector across the top, see §3.0, which also explains why it is a selector and
not another row of tabs.

---

## 3. Sequence analysis

### 3.0 How the analysis is laid out

The analysis produces about fifteen different views of a set of sequences.
Showing
them as fifteen tabs made two problems at once: nobody could find anything, and
**Streamlit executes the body of every tab on every interaction**: including the
tabs you never open. Moving one slider redrew nineteen figures.

So the views are grouped into six **families**, chosen with the selector across
the
top. Each family answers one question:

| Family | The question |
|---|---|
| 🎯 **Triage** | Which of these do I keep? (§3.2b) |
| 🔍 **Identity** | What is each sequence? (§3.3–3.6) |
| ⚠️ **Developability** | Will it behave? (§3.7–3.8) |
| 🔗 **Comparison** | How do they compare with each other? (§3.9, 3.10, 3.12) |
| 🌍 **Repertoire** | How do they compare with real VHHs? (§3.11) |
| 🧪 **Outputs** | What do I make from them? (§3.13–3.15) |

Only the selected family executes. That is the whole point of a selector rather
than nested tabs, Streamlit would run the bodies of nested tabs too. The
effect,
on seven sequences:

| | Time |
|---|---|
| First display after **Analyze** | 1.7 s |
| 🎯 Triage | 0.3 s |
| 🔍 Identity | 5.8 s |
| ⚠️ Developability | 2.1 s |
| 🔗 Comparison | 9.1 s |
| 🌍 Repertoire | 3.3 s |
| 🧪 Outputs | 0.9 s |

Before this, every analysis paid all of it at once: roughly twelve seconds
before
anything appeared, on every single interaction. Comparison is still the
expensive
one, because it carries the alignment, the logo and the identity matrices, but
you
only pay for it when you ask for it, and the cost has stopped being silent.

Repertoire sits in its own family rather than inside Comparison for the same
reason: it is a different question (your sequences against 1958 natural VHHs
rather
than against each other), and keeping it separate takes Comparison from 14 s to
9 s.

Inside a family, the views are still tabs, there are at most four, and they are
cheap relative to the family itself.

### 3.1 Giving it sequences

Three ways in:

- **Paste** into the text box: FASTA, or a bare sequence, or several bare sequences
  separated by blank lines. Whitespace, digits and gap characters are stripped.
- **Upload a FASTA file** with the uploader on the right.
- **Load example** fills the box with three published VHHs (caplacizumab, cAbLys3 and Nb80) so you can see what every tab looks like before committing your own
  data. They come from the library (§8), where their PDB entries are given.

Click **Analyze**, or press ⌘ + Enter (Ctrl + Enter on Windows and Linux). The
button
exists because the keyboard shortcut is not discoverable; both do the same
thing,
which is to commit what you typed.

A progress bar runs directly under the button while the analysis is working,
and is
replaced by a line that stays, which matters more than it sounds, because on a
handful of sequences the analysis itself takes a fraction of a second and a bar
that cleared itself would leave no evidence anything had happened.

That line reports **two** numbers, and the gap between them is the honest
picture:

- *characterised in …* (the physico-chemical analysis and the CDR annotation. Typically well under a second. - *page built in …*) everything, including every figure. Streamlit re-executes
  the whole script on each interaction, so all the tab bodies render whether
you
  open them or not: the alignment, the logo, the identity matrices, the
  aggregation profiles, the plasmid maps. On a fresh set of six sequences this
is
  routinely a hundred times the analysis time.

So when the bar fills and Streamlit's own spinner (top right, fixed position,
not
movable) keeps turning, that second number is what it is still working on. A
second run over the same sequences is far quicker: the figures are cached.

**Limits.** Up to 50 sequences per run, each between 25 and 1500 amino acids, and
12 000 residues in total across the submission. Anything outside that is
rejected
with a message naming the offending sequence.

The total-residue ceiling is the one that may surprise you: one long protein is
cheap to analyse, fifty of them are not, and on Community Cloud the instance is
shared with whoever else is using the app at that moment. Split a large
submission
rather than working around it.

A sequence containing a non-standard residue (B, X, Z, a stop…) is rejected by
name.
Fix or remove it rather than letting it through.

**Clearing the box.** The 🗑️ **Clear** button empties the sequence field in one
click instead of selecting everything by hand. It is greyed out when the box is
already empty, and it also resets the file uploader, so re-selecting the same
FASTA after clearing works, which it would not if the uploader still thought it
had already handled that file.

### 3.1b Generating synthetic VHHs

The panel **Generate synthetic VHH sequences** builds sequences that look like
real
VHHs without being any of them. Useful for negative controls, for benchmarking
a
pipeline, or for testing this app itself.

**How it works.** Frameworks are sampled position by position from an IMGT
position-frequency matrix built from real VHHs. CDRs are drawn at random, with
lengths taken from the same reference set.

**The reference set.** The bundled matrix comes from **PLAbDab-nano** (OPIG,
University of Oxford): 1958 unique VHH sequences, numbered with ANARCI against
the
alpaca germlines, giving 212 IMGT positions with a median support of 1890
sequences
per position. The panel states this, along with the support, because a
frequency is
only as good as the number of sequences behind it.

> *PLAbDab-nano: a database of camelid and shark nanobodies from patents and
> literature*, [doi:10.1093/nar/gkae881](https://doi.org/10.1093/nar/gkae881).
> PLAbDab-nano is freely available to download and query at
> <https://opig.stats.ox.ac.uk/webapps/plabdab-nano/>. **Please cite it** if you
> publish work based on sequences generated here. What this app bundles is an
> aggregate position-frequency matrix; the sequences themselves are not
> redistributed.

**Controls.**

| Control | Effect |
|---|---|
| Frameworks | `sample` draws each position from its frequency column, giving a different framework each time; `consensus` takes the most frequent residue, giving one canonical framework |
| CDR residues | `natural` uses the reference set's own CDR composition; `no_cys` is uniform over 19 residues; `uniform` is all 20 |
| Enforce camelid hallmarks | resamples IMGT 42/49/50/52 until the tetrad is camelid |
| Allow cysteines in frameworks | off by default, see below |
| Fix the CDR lengths | otherwise lengths come from the reference distribution: CDR1 and CDR2 around 8, CDR3 a mean of 17 over a range of 2 to 30 |
| Random seed | the same seed always gives the same set, so a benchmark is reproducible |

Generated sequences download as FASTA, or go straight into the analysis above.

**Three things to know before using the output.**

*Positions are sampled independently.* Real frameworks have correlated positions (a residue at one site constrains another) so a combination that never co-occurs in
nature can come out. The hallmark tetrad is the usual casualty, which is why it
is
enforced by default.

*Cysteines in frameworks are off by default.* About 2% of the reference sequences
carry a cysteine at IMGT 50, where it pairs with a cysteine in a CDR. With
randomly
generated CDRs that partner does not exist, so such a residue would be left as
an
orphan thiol.

*These are synthetic sequences.* Nothing here predicts that they fold, express or bind
anything. They are statistically plausible strings, and that is all.

A check worth running: generate a set and send it through the analysis. On the
bundled
matrix every generated sequence annotates cleanly and scores 4/4 on the
hallmarks,
with molecular weights around 13.5–14 kDa, which is the sanity test that the
matrix
and the assembly are behaving.

### 3.1c What kind of sequence is this?

Every submitted sequence is classified before anything else runs, because the
app
has to know which analyses legitimately apply to it. There is no structure
here, so
the decision is made from the sequence alone:

1. **ANARCI numbering succeeds** → the sequence is an immunoglobulin variable
   domain, and ANARCI also reports whether the chain is heavy (H) or light
(K/L).
2. For a **heavy-chain** domain, the FR2 hallmark tetrad separates a VHH from a
   conventional VH.
3. **ANARCI fails** → not a variable domain. Everything generic still applies.

That gives six verdicts, shown in the `type` column of the Properties table:

| Verdict | What happens |
|---|---|
| **VHH** | everything runs |
| **V domain, VHH/VH ambiguous** | everything runs, with a caution (1/4 hallmarks usually means a humanised VHH or a VH-like construct | | **Conventional VH** | everything runs, labelled VH throughout) see below |
| **Light-chain V domain** | antibody-specific tabs skipped |
| **Contains a V domain** | antibody-specific tabs skipped, a V domain covering less than 80% of the sequence means a fusion or a multi-domain protein, and running CDR analysis on the whole thing would be meaningless |
| **Other protein** | antibody-specific tabs skipped |

**A conventional VH gets the full analysis.** Most of what this app calls "the
VHH analysis" is not camelid-specific at all: IMGT numbering, the CDR windows,
germline assignment, humanness and CDR3 grouping are exactly the same operation
on a VH. Withholding them would have been an arbitrary restriction, so a VH
goes
through all of it, with a `type` column in every per-sequence table saying
which
is which, and a notice at the top naming the VH sequences.

Exactly two readings are genuinely camelid, and both are flagged where they
appear rather than used to gate everything else:

- **The FR2 hallmark tetrad.** On a VH a score of 0/4 is simply what a VH looks
  like, not a defect. The row is still shown, because seeing `VGLW` next to a
  camelid `FERF` is the clearest possible illustration of the difference.
- **The additional CDR1–CDR3 disulfide.** The pattern matched against is a
  camelid one. On a VH, an extra cysteine pair is reported as unverified rather
  than as the known VHH motif. The canonical IMGT 23/104 disulfide is the same
in
  both and applies to either.

Humanness on a VH behaves as you would expect and the app says so: a high
framework identity to a human germline is expected rather than informative, and
the hallmark rows are not a humanisation decision there, since a VH already
carries the VH residues.

When something that is **not** a heavy-chain domain is submitted, a notice
names
each sequence and says why. The tabs that skip them (**CDR / Framework**,
**Hallmarks & germline** including humanness, **Clonal families**, and the
disulfide section of **Liabilities**) say so rather than quietly showing a
shorter list.

**Two behaviours worth knowing.**

*Cysteines are judged differently.* For a VHH the app expects exactly one
intradomain disulfide between FR1 and FR3, and flags anything else as a
possible
free thiol. That rule is wrong for any other protein, lysozyme's eight
cysteines
form four ordinary disulfides and are not eight liabilities. So for a non-VHH
the
app falls back to the only thing that holds generally: an odd number of
cysteines
guarantees at least one unpaired thiol; an even number tells you nothing
without a
structure.

*The alignment defaults to the antibody domains.* Aligning a 1200-residue protein
against 120-residue VHHs produces a figure that is mostly gaps, and a logo to
match. So the **Alignment** tab aligns the heavy-chain subset by default (VHH
and VH together, which align perfectly well) with a checkbox to include the
rest. If
nothing in the submission is an antibody domain, it aligns what you gave it.

### 3.2 Sidebar settings

| Setting | Default | What it changes |
|---|---|---|
| pH for net charge | 7.4 | the pH at which net charge and the developability flags are evaluated |
| Profile window | 9 residues | the sliding window for the per-residue charge and hydropathy profiles |
| CDR annotation method | auto | `auto` uses ANARCI when it is installed and falls back to the regex anchors; `regex` forces the dependency-free path; `anarci` appears only when ANARCI is actually available |

The caption under the selector tells you which annotator is live. If it says
*not
installed*, every CDR boundary you see comes from the regex anchors.

### 3.2b Triage

The first tab, and the one for the moment you have to **decide** rather than
study: sixteen clones out of a screening, ten slots in a synthesis order. Every
flag the app can raise, one row per sequence, sortable by whichever column
matters today.

The columns come in groups you can show or hide (identity, physico-chemical,
liabilities, aggregation, repertoire, cysteines, antibody-specific) and a
*flags* view reduces everything to booleans, so "show me only the sequences with
something wrong" is one checkbox.

Every number here is produced by the same function that fills its own tab, so
the
triage table and the detailed tab can never disagree.

**There is no overall score, deliberately.** Weighting a deamidation site against
a hydrophobic patch against three divergent framework positions means choosing
what matters, and that depends entirely on what the VHH is for, a reagent for
your own ELISA and a candidate you mean to humanise have nothing in common. A
single number would hide that choice behind a decimal point. The columns sit
next
to each other; you sort.

Empty cells are analyses that do not apply rather than failures. A sequence
that
is not an antibody domain has no hallmark score, and its cysteines are not read
against the VHH canonical pair, lysozyme's eight cysteines form four ordinary
disulfides and are not eight orphan thiols.

### 3.3 Properties

The main table, one row per sequence.

| Column | Meaning |
|---|---|
| Length, MW | a VHH is typically 110–135 aa and 12–15 kDa |
| Theoretical pI | pH of zero net charge, Bjellqvist pKa set. A pI close to your buffer pH predicts poor colloidal stability |
| ε280 | molar extinction coefficient at 280 nm, reduced cysteines |
| A280 (1 mg/mL) | ε / MW: **this is the number you divide your measured A280 by to get mg/mL** |
| Instability index | Guruprasad dipeptide score; below 40 predicts a stable protein |
| Aliphatic index | relative volume of A, V, I, L; higher correlates with thermostability |
| GRAVY | mean Kyte–Doolittle hydropathy; negative means hydrophilic. Above about −0.2 flags an aggregation-prone candidate |
| Aromaticity | fraction of F, W, Y; this is what drives the A280 signal |
| Charge | net charge at the working pH |
| n_Cys | cysteine count; two is the canonical intradomain disulfide |

GRAVY and the instability index are shaded green-to-red so the outliers in a
set jump
out.

**Developability flags** below the table turn those numbers into sentences. A flag is
a threshold crossing, not a verdict: instability index ≥ 40, GRAVY above −0.2,
net
charge magnitude below 2, pI within one unit of the working pH, more than two
cysteines, three or more methionines. Tune the thresholds to your own data in
`vhh_core.stability_flags` if they do not match your experience.

**Custom features** is a second table carrying descriptors that are easy to extend:
charge at pH 5 (relevant for cation-exchange capture), the Lys+Arg fraction,
CDR3
length, CDR3-only GRAVY, and counts of NG and DP/DG motifs. Each definition is
listed
in the expander underneath.

**Lysines and primary amines.** `n_Lys` counts lysine side chains and
`n_primary_amines` adds the free N-terminal α-amine, which an NHS ester labels
just as readily. Together they are the sites a random amine-coupling chemistry
(NHS-biotin, NHS-fluorophore, amine coupling onto a chip) will react with.

An expander lists every one of them with its region, and warns when any sits
inside a CDR. That warning is the point of the feature: random amine coupling
hits a CDR lysine along with the rest, and it is the usual explanation for a
VHH
that binds on the bench and stops binding once labelled. Site-specific
chemistry (a C-terminal cysteine, sortase, or an Avi tag) avoids the problem
entirely, and
the Constructs tab (§3.14) builds those.

Counted from sequence, so this is how many amines **exist**, not how many are
solvent-accessible. That needs a structure. Treat it as the upper bound on
labelling sites.

### 3.4 Plots

Four figures:

1. **Secondary-structure propensity** (stacked bars of helix/turn/sheet fractions. 2. **Titration curves and net charge**) each curve crosses zero at that sequence's
   theoretical pI, and the dotted line marks your working pH.
3. **Per-residue profiles**, sliding-window charge and Kyte–Doolittle hydropathy
   along the sequence, with the CDRs shaded. This is where charged patches and
   hydrophobic stretches become visible.
4. **Predicted secondary structure (Chou–Fasman)**, helix and sheet propensity
   curves with a predicted-state track underneath.

The Chou–Fasman heuristic is roughly 50–60% accurate. It shows where the
sequence
*tends* towards each conformation; it does not predict the fold. A VHH is an
immunoglobulin β-sandwich, so a sheet-dominated profile is the expected result,
and
getting one is not evidence that anything was predicted correctly.

With more than one sequence, only the selected Chou–Fasman panel is drawn. Tick
**Show all** to see every one, with many sequences this is the slowest thing the app
does.

### 3.5 CDR / Framework

A table of CDR1, CDR2 and CDR3 with their lengths, the coloured sequences, and
the
framework regions in an expander.

**Two annotation methods, and they do not agree by design.**

*ANARCI (IMGT numbering)* aligns the sequence to germline HMMs, gives every residue an
IMGT position, and takes the CDRs as the fixed windows 27–38, 56–65 and
105–117.
Because it is numbering-based it tolerates unusual CDR lengths. It needs HMMER.

*Regex anchors* needs nothing. It locates the conserved landmarks of a VHH and cuts
between them:

| Anchor | Pattern | Role |
|---|---|---|
| conserved Cys23 | `C` in residues 16–30 | ends FR1 |
| FR2 | `W[A-Z]RQ` (the `WFRQ`/`WYRQ` hallmark | starts FR2 | | FR3 | `[KR][FYW][TASV][IVLFM]S`) the `RFTIS` motif | ends CDR2 |
| second Cys | `[YFHV][YFHVAC]C`, the `YYC` motif | ends FR3 |
| FR4 | `WG[A-Z]G[TAS][A-Z]{1,2}VTVSS` | ends CDR3 |

**The two give different CDR1 and CDR2 lengths.** IMGT windows are tighter than the
anchor-based cuts, which follow a more Kabat-like logic. CDR3 usually agrees.
This
matters when you compare CDR lengths across datasets: say which method produced
them.

Where both run, a disagreement beyond that usually means a non-canonical
framework,
which is itself worth knowing. The regex method assumes a canonical VHH
scaffold; if
an anchor is missing the app names the one it could not find rather than
guessing.

### 3.6 Hallmarks & germline

**The hallmark tetrad** is the four FR2 positions that distinguish a camelid VHH from
a conventional VH. The table gives IMGT and Kabat numbering, the residue found,
the
VH3 consensus, and a verdict.

| IMGT | Kabat | VH3 | Camelid |
|---|---|---|---|
| 42 | 37 | V | F, H, L, Y |
| 49 | 44 | G | D, E, G, Q |
| 50 | 45 | L | C, L, P, R |
| 52 | 47 | W | F, G, L, S, W |

A score of 4/4 is the canonical camelid profile. **IMGT 49 is shared with VH3,
so 3/4
with a G at 49 is common and still a genuine VHH.** Engineered or humanised
constructs
sit between the two profiles. Read the score as an indication, never as a
classification.

**Closest germline** runs the assignment twice, once against alpaca and once against
human, and reports the closest V and J gene with percent identity. The identity
to the
closest human germline is a first-pass humanness indicator, a sequence can be
germline-like and still contain a T-cell epitope, so this is not an
immunogenicity
prediction.

Both need ANARCI. Without it, the hallmarks fall back to reading the positions
off the
regex FR2 block, and the germline table is absent.

**Humanness, position by position.** Underneath the germline table, percent
identity to the closest human germline V gene, computed **over framework
positions only**, because CDR divergence from a human germline is expected and
is
not something you would ever try to fix.

The identity is computed position by position in IMGT numbering against the
IMGT-gapped germline sequences that ship with ANARCI, not by a fresh alignment.
On a typical VHH it agrees with ANARCI's own independent germline assignment to
within 0.1%, which is a useful cross-check that the numbering is right.

What makes it actionable is the divergence table: one row per framework
position
where the VHH differs from the closest human germline, with the IMGT number,
the
region, both residues, whether the substitution is conservative, and whether
the
position is one of the four hallmarks.

That last flag is the whole point. FR2 typically comes out around 65% human
while
FR1 and FR3 sit near 90%, not because FR2 is unusual, but because FR2 is where
the hallmarks live. Reverting them raises the humanness number and
simultaneously
removes the hydrophilic camelid surface that lets a VHH stay soluble without a
light chain. The other divergent positions are the cheap ones; the hallmarks
are
a real trade-off, and the table separates them so you can see which is which.

**A limit worth repeating.** This is a first-pass humanness indicator and nothing
more. It is not an immunogenicity prediction: immunogenicity depends on T-cell
epitopes, dose, route and formulation, none of which is in the sequence. A
germline-like sequence can carry an epitope; a divergent one can be silent.

### 3.7 Liabilities

Chemical and manufacturability liabilities, **scored** and placed along the
sequence.

**The scoring scale.** Each motif carries a base score, and the total per sequence is
the sum:

| Liability | Motif | Score |
|---|---|---|
| N-glycosylation | `N-x-S/T`, x ≠ P | 15 |
| Deamidation (high) | `NG` | 10 |
| Isomerisation (high) | `DG` | 10 |
| Free cysteine | a Cys outside the canonical FR1/FR3 pair | 8 |
| Deamidation (moderate) | `NS`, `NT`, `NN` | 5 |
| Isomerisation (moderate) | `DS`, `DT`, `DD` | 5 |
| Oxidation (Met) | `M` | 4 |
| Oxidation (Trp) | `W` | 3 |

Risk categories follow the same scale: **Low** below 2, **Medium** below 6,
**High**
at 6 or above.

Three further motifs (fragmentation (`DP`), integrin binding (`RGD`) and
hydrophobic patches of four or more consecutive F/I/L/V/W/Y) are **additions**,
marked as such in
the `source` column. Untick *Include the motifs added on top of the in-house
scale* to
reproduce the original totals exactly.

**An important limitation.** In its original form each motif is also weighted by its
**solvent exposure** (SASA) and by the **model confidence** (pLDDT), both taken from a
3D structure. This section is sequence-only, so neither factor is applied: an
`NG`
scores 10 whether it is buried in the core or fully exposed on a loop. Two
consequences
follow, and they matter:

- the categories **compress towards the top**, almost everything lands in Medium or
  High, where the structure-weighted version would push buried motifs down to
Low;
- **comparing sequences with each other tells you more than the category of any single
  motif**. Use the totals to rank a panel, not to declare one motif dangerous.

If you have an AlphaFold model, the structure-weighted analysis is a separate
notebook
and gives a far sharper answer.

**The liability map** places every motif at its residue position, one row per sequence,
with bubble size proportional to the score and the CDRs shaded behind. This is
the view
that answers the question you actually ask of a panel: *where* do the problems
sit, and
do they fall in a binding loop. The **liability profile** underneath stacks
each
sequence's score by motif type, with the total on top.

**The region column is what makes the table actionable.** The same `NG` means different
things in FR3 and in CDR3: in a framework it is mainly a stability and
heterogeneity
concern, in a CDR it may alter binding.

**Free cysteines.** With only the sequence, the app assumes the canonical intradomain
disulfide joins the FR1 and FR3 cysteines, and flags any Cys outside those
regions as a
candidate free thiol. cAbLys3, which carries a genuine second disulfide between
CDR1
and CDR3, is therefore flagged twice, correctly in the sense that those
cysteines are
extra, and misleadingly in the sense that they are paired. With a structure you
would
measure the SG–SG distances instead.

**Protease sites** sit in an expander: trypsin-like (after K/R, not before P),
chymotrypsin-like (after F/W/Y, not before P), GluC-like, AspN-like, and the
high-scoring Factor Xa and enterokinase recognition sites. Every protein has
tryptic
sites, so the count on its own means little. What is worth looking at is a
high-scoring
recognition site, or a cluster of sites inside a CDR.

Everything downloads as a single CSV.

**Cysteines and disulfides.** The liability table flags every cysteine outside
the canonical FR1/FR3 pair as a possible free thiol. That is the right default
and
it is wrong in one common case: camelid VHHs with long CDR3s frequently carry a
**second, stabilising** disulfide, most often between CDR1 or FR2 and CDR3. Read
as a liability, a stabiliser becomes a problem to fix, which would be exactly
backwards.

This section tells the two apart. It numbers every cysteine in IMGT, identifies
the canonical IMGT 23/104 pair, and then looks at what is left: two extra
cysteines in regions that pair in real VHHs are reported as a plausible
additional
disulfide, one extra cysteine as an orphan thiol, three or more as something
that
cannot be resolved from sequence.

Pairing here is inferred from **position alone**. There is no structure in this
app, so "plausible pair" means two cysteines sitting where real VHHs pair them,
not a verified bond. A non-reducing gel or a structure is what confirms it.

**Suggested repairs.** One conservative substitution per liability, so the table
stops being merely descriptive: NG fixed at the glycine rather than the
asparagine, DG likewise, sequons broken at the +2 residue, M→L, free Cys→Ser,
and
so on. Each row carries its rationale.

Three behaviours are worth knowing about.

*Conserved positions are refused.* The oxidation rules match every Met and Trp in
the sequence, which includes the conserved FR2 tryptophan (IMGT 41) and the FR4
tryptophan of the WGxG motif (IMGT 118). Substituting either would "repair" a
liability by breaking the fold, so those rows appear in a separate table with
the
reason, and no substitution is offered. The canonical disulfide cysteines are
treated the same way.

*CDR positions are flagged, not hidden.* A liability in a CDR is often one you
keep, fixing it may cost the affinity you selected for. The rows say so, and
the
default selection is the framework-only subset.

*Substitutions that do not match are refused rather than applied.* Each one is
written `N55Q`, and if residue 55 is not an asparagine the substitution is
rejected with a message instead of silently doing something else. Two
conflicting
substitutions at the same position are refused as a pair, because that conflict
is
your decision.

Applying a selection shows the repaired sequence, the new liability score with
the
delta, and a FASTA download. The recomputed score reuses the original region
boundaries, which is correct for substitutions; re-run the full analysis on the
repaired sequence before trusting anything else about it.

### 3.8 Aggregation-prone regions

Where along the sequence the aggregation-prone stretches are. A stretch becomes
aggregation-prone when it is hydrophobic, β-prone and carries no charged
residue
to interrupt it, and the profile combines exactly those quantities: the
absolute
charge per residue, Kyte–Doolittle hydropathy, Chou–Fasman β-sheet propensity,
and
a gatekeeper term for P/D/E/K/R. Positive leans soluble, negative leans
aggregation-prone.

**This is not CamSol, AGGRESCAN or TANGO.** Those fit their coefficients to
experimental data. The weights here are a transparent judgement call and sit in
`SOLUBILITY_WEIGHTS` in `vhh_core.py` for you to change. The consequence is
specific: read **where** the weak spots are, and compare candidates against
each
other, but do not read the absolute value as a predicted solubility. For a
validated number, run the sequence through one of those servers.

With that said, the profile behaves the way it should on controls. A
poly-isoleucine
stretch scores about −1.0, a poly-lysine stretch about +2.0, and a well-behaved
VHH
comes out around +0.3 with its weakest spots on the known hydrophobic framework
segments. Give it a deliberately hydrophobic CDR3 and CDR3 goes straight to the
top of the list.

**Two controls.** The *smoothing window* (default 9) sets how long a stretch has
to be to register; the *threshold* (default −0.3) sets what counts as
aggregation-prone.

**Two tables.** *Regions below the threshold* is the headline, but a good VHH
often has none at all, which would leave the tab empty and tell you nothing. So
*weakest windows* always lists the three lowest-scoring windows whether or not
they cross the threshold, with their region and whether they fall in a CDR.
That
is the table to compare candidates on.

In the summary, `Mean` is higher-is-better and `Worst` is the single lowest
window. `Worst` is usually what decides whether a candidate behaves: one bad
patch
is enough.

### 3.9 Alignment

With two or more sequences: a FAMSA alignment rendered with `pymsaviz`, with
the CDR
bands of the first sequence marked, downloadable as aligned FASTA.

These two figures are the most expensive in the app, and they used to sit
behind a
**Draw the alignment and logo** button because Streamlit re-executes every tab body
on every interaction, so left to themselves they were redrawn constantly,
whether
or not anyone was looking. The family selector (§3.0) solves that properly: the
Comparison family only executes when you select it, so the button became pure
friction and was removed. The cost of these figures is now paid when you ask
for
the family that contains them, and at no other time.

A *Positions per row* slider sets the wrap, 20 to 140 columns, exactly as the
logo
below has its own. Fewer positions per row gives more rows and wider letters,
which
is what you want for a figure; more keeps the whole domain on one or two lines
for a
quick look. The two sliders are independent, the alignment and the logo can be
wrapped differently.

Underneath, **two identity matrices**: full-length and CDR3-only. The second is
the
one that matters when triaging panning output. Two clones can share 90%+
full-length
identity while having unrelated CDR3s, they are the same scaffold, not the same
binder. To turn that observation into actual groups, use **Clonal families**
(§3.12).

**Sequence logo.** Underneath, the same alignment as a logo: letter height shows
what is conserved and what varies, which reads faster than the alignment itself
when all you want to know is where the diversity sits. It reproduces the logo
from
the clinical-VHH alignment notebook (`logomaker`, gaps and X ignored, no
pseudocount, amino acids that never occur dropped so the colouring only shows
residues actually present) with the CDR bands of the first sequence mapped onto
alignment columns.

Three controls:

*Height* is **Frequency** or **Information (bits)**. Frequency makes every column
sum to 1, so the logo reads directly as composition, the notebook's choice.
Information is the classic Schneider logo in bits, where a mixed column stays
short even when one residue is commonest; it is the better one for spotting
which
positions genuinely carry signal.

*Positions per row* sets the wrap, 20 to 140 columns. Fewer positions per row means
more rows and wider letters, which is what you want for a figure; more
positions
per row keeps the whole domain on one or two lines for a quick look.

*Colours* offers the usual logomaker schemes, chemistry by default.

The underlying matrix downloads as CSV, one row per alignment position and one
column per amino acid, if you want to plot it yourself.

### 3.10 vs Reference

Which positions differ from a sequence you pick as the reference. Humanisation
variants against the parent, a Sanger clone against the construct it was
supposed
to be, a panel against its ancestor, the same question each time.

Choose the reference from the submitted sequences, or paste one that is not in
the
set. Everything else is compared against it.

**Positions are matched by IMGT number, not by raw alignment.** That is the whole
point of doing it this way: delete one residue from a CDR3 and an
alignment-based
diff reports every downstream position as changed, whereas this reports one
deletion. Insertions and deletions are reported as such rather than dropped, in
a
CDR3 they are usually the most interesting difference there is.

One honest caveat, visible in the app: IMGT places insertions at the *centre*
of a
loop, so adding a residue to CDR3 can shift a couple of its neighbours. Those
show
up as substitutions alongside the insertion. It is a property of the numbering,
not
a mistake, but it means an insertion rarely reports as exactly one change.

The summary gives, per sequence, the number of differences split into
substitutions
and indels, how many fall in a CDR, and how many substitutions are conservative
by
BLOSUM62. The map underneath places every difference along the IMGT axis with
the
CDR bands marked (red for non-conservative, blue for conservative, grey for an
indel) so a humanisation variant reads at a glance as "four changes, all in
FR2".
A *CDR positions only* toggle narrows everything to the loops, which is the
view
you want when checking whether a clone kept its binding site.

### 3.11 Repertoire

Your sequences against the natural VHH repertoire. This is the same IMGT
position-frequency matrix the synthetic generator samples from (§3.1b,
PLAbDab-nano,
1958 VHHs), used the other way round.

**Unusual framework residues.** For each sequence, the framework positions whose
residue is rare in the reference. This is what catches a Sanger miscall, a PCR
error, or a construct that has quietly drifted. On a normal VHH it flags
nothing
across ninety-odd compared positions; plant a proline in the middle of FR3 and
it
comes out at 0.7% against a consensus alanine at 83%, with 1952 reference
sequences
behind the call.

CDRs are excluded on purpose. A CDR is *supposed* to be unusual, and flagging
it
would bury the framework signal that actually matters.

Two guards are built in rather than left to you. The *minimum support* slider
sets
aside positions backed by too few reference sequences instead of calling them
rare, a frequency from 30 sequences is not a frequency from 1900, and every row
carries
its support. And the wording stays "unusual", never "error": **rare is not
wrong**.
A residue absent from the reference may be a genuine novelty, a germline this
set
under-samples, or the very mutation you selected for.

**Your set vs the natural repertoire.** CDR length distributions and per-region
amino-acid composition, your set shown against the reference.

How your sequences are drawn depends on how many there are, and the reason is a
real one. Plotted as a frequency, a single sequence is 100% of itself: its bar
reaches 1.0, the reference distribution it is meant to be compared against is
squashed into the bottom tenth of the axis, and the figure tells you nothing.
So
below eight sequences they are drawn as **markers placed on the reference
curve**
instead, each labelled with how common that length actually is among real VHHs.
At eight or more, both are drawn as bars and the comparison of shapes becomes
the
point.

Either way, a **"where each CDR length sits"** table gives the number that
matters: the frequency of that exact length in the reference, and the
percentile, the fraction of natural VHHs with a shorter CDR there. A CDR3 of 13
comes out at
the 15th percentile against a distribution peaking at 17, which is a far more
useful statement than any bar. Lengths below the 5th or above the 95th
percentile
are called unusual, and unusual is not wrong: a long CDR3 is often exactly what
you selected for. After a panning this is
what says whether the selection was biased: a shift towards long CDR3s, an
enrichment in charged residues, a depletion of the tyrosine and glycine that
are the
staples of natural CDR3s. Drop a deliberately charged set in and lysine comes
out
enriched by thirty points over the reference, with tyrosine and glycine
depleted, which is exactly the kind of signature a charged target leaves
behind.

With only a handful of sequences the comparison is coarse by construction: one
sequence in ten is a 10% bar. Read it for the direction, not the number.

**A note if you submitted a conventional VH.** The reference is a *VHH* repertoire.
A VH compared against it may legitimately show FR2 positions as unusual, that
is
the difference between a VH and a VHH, not a defect in your sequence. The app
warns
when VH sequences are present.

### 3.12 Clonal families

This is the tab for triaging panning output. The first question after a
screening is
not how good each clone is, but **how many different clones you actually
have**, a
plate of picked colonies is routinely the same few lineages over and over, and
properties computed on 16 reads that represent 3 lineages will mislead you
three
times over.

Grouping is done on **CDR3**, because CDR3 is what differs between independent
lineages while the frameworks stay nearly constant.

**How it works.** Every pair of CDR3s is aligned (global, BLOSUM62, gaps −11/−1) and
scored for identity, strictly identical residues over the length of the longer
CDR3. Those identities become a distance matrix (`1 − identity`), which is
clustered
hierarchically and cut at your threshold. This is the same procedure as the
Sanger
screening notebook, so the app and the notebook give the same families for the
same
reads.

**The three controls.**

*CDR3 identity threshold*, 50–100%, default 80%. Two sequences at or above it belong
to the same family. At 100% only strictly identical CDR3s group together.

*Linkage.* `complete` guarantees that **every pair** inside a family meets the
threshold, the conservative reading of "same family". `average` merges on the
mean
distance and is more permissive, so a family can then contain pairs that fall
below
your threshold. The family table's `min_identity` column is there to catch
exactly
that: it reports each family's worst internal pair, and a value under your
threshold
tells you `average` merged across it.

*Same CDR3 length only* forces CDR3s of different lengths into different families
whatever their alignment. Off by default, because a one-residue indel in a CDR3
is a
plausible sibling rather than a new lineage. Turn it on if you want strictly
length-matched families.

**What you get.** Four numbers at the top, sequences grouped, families, singletons,
and redundancy (sequences per family). Redundancy above roughly 2 means you are
re-picking the same lineages and could screen more colonies for the same
effort.

Then a family table, with for each family its size, how many distinct CDR3s it
contains, its worst internal identity and similarity, its CDR3 length range,
and a
**representative**: the family's most frequent CDR3, a reasonable default for the
clone to carry forward.

A dendrogram draws your threshold as a dashed line; branches joining to its
right
are separate families. Move the slider and watch which groups merge: that is
the
honest way to judge how robust a family is. A grouping that holds from 70% to
95% is
real; one that appears only at 78% is an artefact of where you put the line.

Finally a per-sequence assignment table, and CSV downloads of the families, the
assignment and the CDR3 identity matrix, plus a FASTA of one representative per
family, which is usually what you actually want to order or subclone.

**Two honest limits.** Families are a *hypothesis* about shared lineage drawn from
CDR3 similarity alone; two unrelated lineages can converge on a similar CDR3,
and
true siblings with a heavily mutated CDR3 can be split. And sequences whose
CDR3
could not be annotated are left out rather than guessed at, they are listed
explicitly when that happens.

### 3.13 DNA

Back-translation into *E. coli* codon-optimised DNA.

Each residue first takes the most frequent *E. coli* codon. The sequence is
then
repaired wherever it contains a restriction site you asked to avoid or a
homopolymer
run longer than six, by swapping the offending codon for the next synonymous
one.
Finally the GC content is balanced, because optimising codon frequency alone
pushes GC
to about 63% and synthesis vendors flag anything above 60%.

Controls: which restriction sites to avoid (17 available, all selected by
default),
whether to append a TAA stop, and the maximum GC fraction.

**Every output is verified by translating it back and comparing to the input.** If the
round trip fails the app says so in red and tells you not to order that
sequence. For
DNA you are going to pay to have synthesised, that check is worth its cost.

The optimiser maximises codon frequency, balances GC and removes
cloning-unfriendly
motifs. It does **not** model mRNA secondary structure near the ribosome
binding site,
codon-pair bias or tRNA availability. Run the output through your synthesis
vendor's
own checker before ordering.

A **Send to gene diagnosis** button passes the sequences straight to §7, where
the
frame, rare-codon and GC checks run on them, including the GC balancing this
optimiser does globally but not window by window.

### 3.14 Constructs

Assemble a VHH into the construct you will actually order, and get the mass, pI
and DNA of the finished protein rather than of the bare domain. The practical
reason: a VHH is 13 kDa with a pI near 8.6, and the same VHH with a FLAG tag is
a
different protein on a gel and on an ion-exchange column. Four aspartates will
move the pI by two units, and you would rather know that before ordering than
after the purification misbehaves.

**Formats.** *Monovalent* is one VHH with whatever tags you add. *Bivalent* is the
same VHH twice, which raises apparent affinity through avidity, not affinity.
*Biparatopic* is two different VHHs against the same target, which needs
non-overlapping epitopes, picking the same VHH twice there makes a bivalent
construct and the app says so.

**Linkers.** A list from a direct fusion up to (GGGGS)4, plus rigid (EAAAK)
helices for when you want fixed spacing rather than a floppy tether. (GGGGS)3
is
the default for a bivalent construct, long enough for both domains to bind
independently. Fusing two VHHs with no linker raises a warning, because the
domains will interfere and that is almost never what you want.

**Tags.** His6, His8, HA, c-Myc, FLAG, Strep-II, Avi, the original SpyTag, and TEV,
thrombin and sortase sites. Each carries a short note. The sortase motif has to
be
last (LPETG is cleaved between T and G, so anything after it is lost) and the
app
warns if you put something after it.

**What is deliberately missing, and why.** SpyTag2 and SpyTag3 are not in the
list, and neither are the camelid hinges. Only the original SpyTag is included.
Writing the others from memory would risk handing you a sequence that is one
residue wrong, which you would discover only after the synthesised gene
arrived.
The *Custom pieces* expander takes any peptide you paste, with a label, at
either
end, that is the honest way to handle a sequence this app cannot vouch for, and
it is where those belong.

**Verify everything anyway.** Tag sequences differ between vectors and between
papers. The app says so on every run, and it means it: a one-residue difference
will not show up until the protein does not work.

**What you get.** A schematic drawn to scale with VHHs, linkers and tags
distinguished; length, mass, pI and A280 each with the delta against the bare
VHH;
a segment table giving every piece's exact coordinates; the full property set;
the
protein sequence; and codon-optimised DNA for the whole construct, with the
same
restriction-site avoidance as the DNA tab and the same round-trip verification,
the DNA is translated back and compared to the construct before it is shown.

### 3.15 Download

The results table as CSV, the figures as 300 dpi PNG, and a full PDF report.

The report is built on demand because it takes a few seconds. It contains the
physico-chemical summary, then one section per sequence (properties,
developability
flags, regions, the sequence, hallmarks, germline, liabilities), all the
figures, and
a page of caveats.

The **Include the per-sequence secondary-structure plots** checkbox adds one
page per
sequence. It is ticked automatically up to six sequences; beyond that those
plots are
most of the file size and most of the build time.

---

## 4. Structure

Every other section reasons from sequence alone, because a sequence is what you
always have. But several of the questions this app answers are structural
questions being approximated from sequence, and once a model exists the
approximation can be replaced by a measurement:

| Asked from sequence | Measured here |
|---|---|
| which cysteines pair: from IMGT positions and a canonical-pair assumption | the SG–SG distances |
| a liability motif, wherever it occurs | whether that motif is buried or on the surface |
| hydrophobicity in a sliding window | hydrophobic area that is contiguous **in space** |
|, | which residues contact the antigen |

The last row is the point: contact residues are not inferable from sequence at
all.

### 4.1 What it accepts

`.pdb`, `.cif`/`.mmcif`, and gzipped versions of either. Crystal structures
from
the PDB, AlphaFold2, AlphaFold3 and ColabFold models are all read directly,
pLDDT included, for AlphaFold3, use the `.cif` it writes.

The app does **not** fetch structures for you, deliberately: downloading on
your
behalf would mean this server making requests that look like they came from
you.
Download from `rcsb.org` or `alphafold.ebi.ac.uk` yourself and upload the file.

For a VHH–antigen complex, download the **whole entry** rather than one chain,
otherwise there is no interface left to compute.

If the file holds several models (an NMR ensemble, or several predictions) only
the first is used, and the app says so.

### 4.2 Choosing the chain, and what to annotate against

**Chain of interest** is the VHH, usually; everything per-residue describes it.
**Binding partner(s)** are the chains treated as the antigen.

**Annotate against** is the one that is easy to skip and worth understanding. A
model is routinely missing termini and disordered loops, and its residue
numbering may be sequential, IMGT, Kabat or whatever the depositor chose. Left
at *(the structure's own sequence)* the section works on exactly what is in the
file. Point it at one of the sequences loaded in the analysis section instead
and
the model is aligned onto that full-length construct, so CDR bands, regions and
liability positions line up with the rest of the app. The app reports how many
residues aligned and at what identity, and warns below 90%, if you see that,
check you have picked the same molecule.

Internally the mapping is keyed by **residue number**, not by position in the
chain. Those differ the moment a chain does not start at 1 or has a gap, which
is
most of the time.

### 4.3 Chains & confidence

A predicted model stores **pLDDT** in the column a crystal structure uses for
the
**B-factor**. The two mean opposite things (high pLDDT is good, high B is bad) and no field in either format states which it is. The app guesses from the value
range and the header, says why, and lets you override it.

For pLDDT the usual confidence bands are drawn, because the number is close to
meaningless without them: 90 and 60 are both "high-ish" until you know that
**70
is the line below which the backbone should not be trusted**. The residues
below
it are listed.

pLDDT is the model's confidence in itself. It is neither accuracy nor a
measurement of flexibility, although low pLDDT and genuine disorder do tend to
coincide. Side-chain placement is less reliable than the backbone even at high
pLDDT, which matters when you read the accessibility of a long side chain.

### 4.4 Disulfides

Measured, not inferred. A pair of SG atoms closer than **2.5 Å** is a bond;
anything else is a free thiol. The counts are compared with the sequence-based
prediction from §3.7 and with the `SSBOND` records in the file, and
disagreements are reported rather than resolved silently.

When the counts disagree, the measurement is the one to believe, but check
first
that the model covers every cysteine: a cysteine in an unmodelled loop cannot
be
seen to pair.

### 4.5 Exposure

Absolute accessibility in Å², divided by the maximum that residue type can
have,
so a Trp and a Gly are comparable. At or below **20%** is called buried, at or
above **50%** exposed, and the band between is left as *intermediate* rather
than
forced into one of the two.

Computed with Biopython's Shrake–Rupley implementation (1.4 Å probe) and
normalised by the residue maxima of Tien et al. 2013 (PLoS ONE 8:e80635). One
caveat worth stating because it is easy to miss: those maxima were derived with
DSSP, i.e. the Lee–Richards definition, and Shrake–Rupley is a different
numerical approximation of the same quantity. The buried/exposed call is robust
to that; a relative accessibility quoted to the percent is not.

Waters are removed first. Ligands are kept, because a bound ligand really does
cover the surface underneath it. When a binding partner is selected, these are
the accessibilities the residues have **in the complex**.

### 4.6 Liabilities in 3D

The join worth uploading a structure for. The analysis section flags a
deamidation, oxidation or glycosylation motif wherever it occurs; it cannot
know
whether the motif is on the surface or packed into the core. A buried motif is
largely inert. An exposed one in a CDR is a real problem.

Each liability from §3.7 gets its measured relative accessibility, its
buried/intermediate/exposed call, and the residue number it has in the model,
so
you can find it in the 3D view. Motifs in unmodelled regions are reported as
*not modelled* rather than as buried, the two are very different claims.

The sequence-based `risk` column is left untouched on purpose. Downgrading it
automatically because a residue looks buried in one static conformer would be
overclaiming. Read the two columns together.

### 4.7 Surface patches

What drives aggregation is hydrophobic area that is **contiguous on the
surface**, which a sliding window over the sequence cannot see: two residues
far
apart in sequence can sit side by side on the fold. This groups solvent-exposed
V, I, L, M, F, W, Y residues whose side chains come within a contact distance
of
each other, measured between closest heavy atoms, not between centroids, since
Trp and Arg carry their bulk far from the backbone.

Three thresholds, all sliders:

*Minimum relative accessibility* (default **0.15**) is deliberately well below
the 0.50 that counts as "exposed" elsewhere, for a structural reason: a residue
in the middle of a patch is half-covered by the patch's own other residues, so
requiring each one to be exposed in its own right is a good way of never
finding
a patch. On the one complex this was checked against, the known dimerisation
patch is found at 0.15 and missed at 0.20, a cliff one residue wide. Treat the
default as a starting point, not a validated parameter.

*Contact distance* (default 5.0 Å) barely changed the answer over 5.0–7.0 Å in
that same check.

*Minimum patch size* (default 3), fewer than three residues is not a patch.

**A transparent composite, not SAP or CamSol.** Those methods are parameterised
against experimental aggregation data this app does not have. This asks a
simpler
and fully stated question, and the thresholds are yours to move. Read it as
*here
is a hydrophobic area worth looking at*, never as a predicted aggregation rate.

### 4.8 Interface

Only shown when you select a binding partner.

A residue is at the interface if it loses solvent accessibility when the
partner
is present. **ΔSASA** (accessible alone minus accessible in the complex) is the
standard operational definition, and unlike a distance cutoff it weights each
residue by how much surface it actually commits to the contact.

With one partner you also get the buried surface area per side. For reference,
an
antibody–antigen interface typically buries 600–900 Å² per side; a VHH, with
one
domain instead of two, often sits at the lower end or below it. A much smaller
number is worth suspicion, in a crystal structure it may be a packing contact
rather than the biological interface.

When the model is annotated against a sequence, the buried surface is also
broken
down **by region**, which answers the question people actually have: how much
of
the binding is done by the CDRs. A framework contribution is normal for a VHH,
the FR2 face that would be the VL interface in a conventional antibody is
solvent-exposed and often takes part in binding.

The same calculation run the other way round gives the **epitope**: the
residues
of the antigen that the VHH covers. That is a structural epitope from one
model,
not a functional one, it says which residues are covered, not which ones the
binding energy depends on. Those are different sets, and telling them apart
needs
mutagenesis.

### 4.9 3D view

3Dmol.js, loaded from a CDN. Colour by chain, N→C rainbow, pLDDT/B-factor or
buried→exposed; optionally show the surface; and highlight the disulfides, the
free thiols, the largest hydrophobic patch or the exposed liabilities as
sticks.

If the panel stays blank your browser could not reach the CDN. Every table and
figure in this section is computed server-side with Biopython and is
unaffected.

### 4.10 What this section is not

It measures **one static conformer of one model**. It does not predict
structure,
does not dock, does not estimate affinity, and does not model flexibility. For
a
predicted model, side-chain accessibility is indicative rather than measured,
because side chains are the part a predictor places least reliably.

---

## 5. Bench

Every other section reasons about a sequence. This one is about the tube: how
much protein is in it, whether it is the protein you think it is, and how much
labelling reagent to add.

It is the arithmetic that gets done on the back of a glove box on a Friday
afternoon, which is exactly why it is worth writing down once. The views start
from the molar mass, extinction coefficient and lysine count that the analysis
section already derives from the sequence, so none of those has to be typed
again or looked up wrong. The exception is the phage titre, which is about the
step before all of that and needs no sequence at all.

Load your sequences in the analysis section first and pick one at the top. The
calculators still work without that, with the molar mass and extinction
coefficient typed by hand, but typing them is the most common way of getting
this wrong.

**No empirical model lives in this section.** Every number is a measurement you
supply, a constant that cites its source, or arithmetic on those two. Nothing
predicts a labelling efficiency, a yield, or how a conjugate will behave: those
depend on your protein and your hands.

### 5.1 Quantify

Concentration from a measured A280, by Beer-Lambert, with the extinction
coefficient computed from the sequence rather than assumed. The generic
"1 mg/mL = 1 AU" is wrong by about a third for a typical VHH.

Inputs: the A280, the dilution factor (1 for a neat sample, 10 if you measured
10 uL in 90 uL of buffer), and the path length (1 cm for a cuvette; a NanoDrop
pedestal already reports a 1 cm equivalent). Outputs: mg/mL and uM, plus the
total in ug and nmol if you give a volume.

**Cystines formed or all reduced.** The extinction coefficient differs by about
1% between the two, and this is the choice nobody makes explicitly. A VHH from
the periplasm has its canonical disulfide formed, so the cystine value is
usually the right one.

**A260, optional,** buys you the 260/280 ratio. A pure protein sits near 0.57.
Above about 0.8 there is enough nucleic acid for the A280 to overstate the
protein, and the concentration becomes an upper bound rather than a
measurement. The app says so rather than reporting a number that looks the same
as a good one.

Underneath, a dilution: give a target concentration and a final volume and it
tells you how much sample and how much buffer. It refuses to dilute upwards,
which is the arithmetic a spreadsheet happily does.

### 5.2 Mass check

Is the thing in the tube the thing you think it is? Give a measured mass and
the app compares it against the sequence.

A deconvoluted mass from an intact LC-MS run is an **average** mass, so that is
what is compared. The monoisotopic value is reported alongside for when you are
looking at an isotope envelope instead.

When the mass does not match, the app lists the modifications that would
account for the difference **and that this sequence allows**. The second half
matters: gluconoylation needs an N-terminal His-tag, methionine excision needs
an N-terminal methionine, pyroglutamate needs an N-terminal Gln or Glu. An
explanation that cannot apply to your protein is not an explanation, so it is
not offered.

Each delta is computed from a chemical formula rather than copied from a table,
so any of them can be checked by hand:

| Difference | What it usually is |
|---|---|
| about -131 | the N-terminal methionine was excised |
| -2 per bond | a disulfide formed |
| -17 | pyroglutamate from an N-terminal Gln |
| -18 | pyroglutamate from an N-terminal Glu |
| +16 | an oxidation, usually on a methionine |
| +178 | gluconoylation of an N-terminal His-tag, a documented *E. coli* artefact |
| +258 | the phosphorylated form of the same artefact |
| +22, +38 | a sodium or potassium adduct, not a modification of the protein |

**Candidates are offered, not chosen.** More than one combination can land on
the same number, and picking one for you would be presenting arithmetic as an
identification. A peptide map settles it.

One thing to check before believing a mismatch: the comparison is against the
sequence you loaded. If the real protein carries a tag, or has had one cleaved,
compare against that sequence instead.

### 5.3 Conjugate

How much labelling reagent to add. NHS esters react with the lysine side-chain
amines and the N-terminal alpha-amine, and the arithmetic is the one on every
vendor cheat sheet, written out so the units are visible:

```
nmol protein = volume (mL) x concentration (mg/mL) / molar mass (kDa)
nmol reagent = nmol protein x molar excess
volume (uL)  = nmol reagent / stock (mM)
```

The last line is the one worth remembering: **a millimolar stock is one
nanomole per microlitre.**

**Reconstitution and the working stock are two steps,** as they are in the
vendor protocol. NHS esters hydrolyse in water, so the vial is reconstituted
concentrated and used at once, and the working stock is a dilution of that. The
app gives the solvent volume for the first and the dilution for the second.

**Only reagents whose molar mass was checked against a primary source are
listed,** and each one shows where its number came from. Anything else is
entered by hand. That is the honest option: a catalogue number is not a molar
mass, vendors renumber products, and a dye's extinction coefficient and 280 nm
correction factor are printed on the vial because they vary between lots. A
wrong molar mass does not fail loudly, it quietly changes the molar excess and
the result gets blamed on the chemistry.

Three things the app checks that a spreadsheet does not:

- **A volume too small to pipette.** Asking for 0.2 uL and delivering 2 uL is a
  tenfold error in the molar excess, invisible, in the direction that looks like
  it worked. The app offers a dilution instead, with the excess unchanged.
- **Co-solvent above about 10%** of the reaction, which is a protein problem
  rather than a labelling one if the reagent is in DMSO or DMF.
- **An excess beyond the number of primary amines,** past which the extra
  reagent can only react with amines already labelled or hydrolyse.

Then it names the spin desalting column that suits the final volume.

**Lysines in a CDR.** The app lists every primary amine with its region and
warns when any sit in a CDR, because an NHS ester cannot tell them from the
others and random labelling then risks the paratope. If that matters, a
site-specific route (sortase, an engineered cysteine, a C-terminal tag) is the
alternative. Which lysines are actually *accessible* is a structural question
this cannot answer: upload a model to §4.


### 5.4 Phage titre

The step before the rest of this section. An absorbance at 260 nm turns into
particles per millilitre, and that is what decides how much stock goes into the
next round of panning. Nothing here needs a sequence.

You give the reading, how many times the stock was diluted to take it, and a
**glycerol correction**. The three multiply together into the absorbance the
undiluted stock would have, which is then converted twice.

The glycerol part is the one that is easy to get backwards, so the app asks for
the two volumes rather than for a factor. If you read the absorbance on the
phage and *then* add glycerol, the tube you pipette from afterwards is more
dilute than what you measured, so the correction is **below one**: 500 uL of
phage plus 55.5 uL of pure glycerol is 555.5 uL, 10% glycerol, and the factor
is 500/555.5 = 0.900. If instead you read the absorbance after mixing, the
glycerol is already in the reading and the factor is 1. A factor above one is
for the other question, what the concentration was before the glycerol went in,
which is not the titre of the tube you are about to aliquot.

**Two constants, and they do not agree.** This is the whole point of the view,
so it is worth being blunt about it.

The **empirical** route multiplies by a number measured directly: a titration
in which an A260 of 1 corresponded to 3 x 10^10 particles per millilitre. That
is a correspondence observed for a particular preparation on a particular
instrument, which is why the app lets you edit it rather than burying it.

The **optical** route is Day's relation for filamentous phage,

    virions/mL = A269 x 6e16 / (bases in the genome)

which works the absorbance out from the single-stranded DNA the particle
carries. The genome length is an input because a phagemid carrying an insert is
not the length of the helper phage; M13mp18 is 7249 bases. Reading at 260 nm
rather than 269, and without subtracting the A320 scatter, costs some accuracy
but not an order of magnitude.

On the same reading these two differ by a factor of roughly 275. A gap that
size is not rounding, and the app says so rather than choosing for you. It
means one of the two is not describing your stock: most often the empirical
factor was established on a different preparation, or the reading carries a
large scattering contribution that a calculation from the genome has no way to
know about. **Use the one you can defend, and record which one in the
notebook.** A plaque or colony titration settles it outright and is the only
one of the three that counts particles that actually infect.

Below that, **how much to take**: pick which titre to work from, give the
number of particles you want for a solid-phase round, and the app gives two
volumes side by side, the solid-phase one and the liquid-phase one, which is
twice it. Showing only the first and leaving you to double it is how a liquid
round gets set up with half the phage it should have had. The factor of two is
a protocol convention and not something the app derives, so if yours differs,
change the target. A round of panning is usually run with 10^11 to 10^13
particles. If a volume comes out below two microlitres it cannot be pipetted
accurately, and the app offers a dilution instead, exactly as the conjugation
planner does and for the same reason.


### 5.5 Batch sheet

The page that goes with a tube when it leaves your bench. You give it a lot
number, who prepared it, the host, the purification, the buffer and the
storage, plus whatever you measured, and it produces a one-page PDF.

Three blocks, kept apart on purpose:

- **Preparation**, the free text. Nothing here can be checked by anything.
- **Measured**, what was actually put on an instrument. A field left at zero
  stays off the sheet, because a page full of zeros reads as a measurement of
  zero. Nothing is carried over from the other views: a number you have not
  retyped is a number you have not checked.
- **From the sequence**, which is theoretical, every value of it, and labelled
  as such on the PDF.

Keeping them apart is the only reason the sheet is worth generating. A
predicted extinction coefficient quoted a year later as though it had been
measured is the failure this layout exists to prevent.

Two numbers cross the line, because they are arithmetic rather than judgement:
total protein, from the concentration and the volume, and the difference
between a measured mass and the sequence. If that difference is more than a
couple of daltons the screen points you at §5.2, which lists the modifications
that would explain it and that your sequence allows.

**This is not a certificate of analysis.** That term means something specific
in regulated quality control: specifications, acceptance criteria, a release
decision and a signature. There is none of that here, and the PDF says so at
the bottom so the document cannot be mistaken for one later.

---

## 6. Sanger sequencing

### 6.1 What it does

You upload the chromatograms from a screening run. For each read the app:

1. quality-trims with **Mott's algorithm**, the one phred and Staden use;
2. translates the read and its reverse complement in three frames each;
3. locates the **PelB leader** to fix the start and the reading frame;
4. locates the **FR4 signature** to fix the end of the VHH;
5. extracts the VHH protein and DNA, and reports what went wrong where it did.

Because all six frames are searched, forward and reverse reads are handled
without
being told which is which.

### 6.2 Settings

| Setting | Default | Effect |
|---|---|---|
| Vector | pMES3 | which tags are expected after the VHH: pMES3 is PelB–VHH–His6–stop, pMECS is PelB–VHH–HA–His6–stop |
| Trim cutoff | 0.05 | Mott error-probability cutoff, roughly Q13. Lower is stricter |
| Min trimmed length | 150 nt | below this the read is flagged, and the raw sequence is analysed instead of the trimmed one |
| Min mean Phred | 20 | QC threshold on the trimmed window |

To add a vector, edit `sanger.VECTORS`, and take the tag layout from the actual
plasmid map rather than from memory. A wrong layout makes the app report
missing tags
on clones that are perfectly fine.

### 6.3 Reading the results table

Four counters at the top: reads, VHHs extracted, unique sequences, failures.

The table carries the status, QC pass, trimmed length, mean Q, the strand and
frame
where the cassette was found, the VHH length, the tags found and missing,
premature
stops, the copy number, and which read a duplicate duplicates.

Every warning is also reprinted in full underneath the table, grouped by read,
because
a long warning is unreadable in a table cell.

**A failed read is not necessarily a failed clone.**

| Status | Usually means |
|---|---|
| `ok` | cassette found, VHH extracted, no frameshift or internal stop |
| `low quality` | the VHH was extracted but the read failed the QC thresholds |
| `no FR4` | the read stopped before the end of the VHH: sequence further, or from the other end |
| `frameshift` | PelB and FR4 sit in **different reading frames**, so an indel lies between them. The VHH cannot be reconstructed reliably and is deliberately not reported |
| `no PelB` | the leader was not found in any frame: wrong vector, very poor read, or empty well |
| `unreadable` | the file could not be parsed as `.ab1` |

One warning deserves a note: *only the PelB core anchor matched, read starts
inside
the leader*. The short anchor `AAQPAMA` is used when the full leader is not
present.
That is normal for a read that begins part-way through the leader, and the VHH
start is
derived from the anchor that actually matched.

**Identical protein sequences are grouped.** The same colony picked twice shows up as
one sequence with a copy count, and the duplicates are marked.

### 6.4 Quality traces

An expander draws the Phred quality along each read, with the window kept after
trimming in blue and a dashed line at Phred 20. Up to twelve reads are drawn.

### 6.5 Sending clones to the analysis

Below the table, pick the clones worth keeping and press **Send selected to
analysis**.
The app switches to the analysis section with those sequences loaded.

**Duplicates are deselected by default**: one entry per unique protein sequence. Tick
them back on if you want to see them separately.

### 6.6 Downloads

The full per-read table as CSV, and the extracted VHHs as protein or DNA FASTA.

---

## 7. Gene diagnosis

Its own section, reached from the home page or from the **DNA** tab of the
analysis section, which hands its back-translated sequences straight over.

**What you paste stays there.** Leaving the section and coming back does not clear
the box, the same as in the analysis section, Streamlit normally discards an
input's contents once its screen stops being drawn, and both inputs are
explicitly
kept alive against that.

**Several sequences at once.** Paste FASTA, a bare sequence, or several separated
by blank lines, up to 20, 6000 nt each. The overview table puts one row per
sequence: frame, internal stops, GC, rare codons and their clusters, GC-rich
stretches, homopolymers. `in_frame` and `internal_stops` are the deterministic
columns (wrong there means wrong) while the rest is triage. Pick one sequence
from the list for the detailed view and the codon suggestions.

The reverse of the DNA tab: paste a coding sequence that is already in a
plasmid
and ask what might be wrong with it. Written for the case you described, VHHs
that express poorly.

**Read the warning at the top of the tab first.** Codon usage is rarely the first
suspect. A VHH that expresses badly in *E. coli* is more often explained by
periplasmic toxicity, inefficient signal-peptide cleavage, failure to form the
intradomain disulfide, inclusion bodies, plasmid loss or induction conditions.
A
gene with flawless codon usage expresses badly for all of those reasons. This
tab
tells you whether the sequence is *part* of the problem; it cannot tell you
that
it is the problem.

**The deterministic checks come first**, because they are the ones that are
definitely wrong when they are wrong: frame, start and stop codons, internal
stops, and (if you supply the expected protein) whether the translation matches
it. That last one is worth doing every time: a mismatch makes everything else
on
the page beside the point.

**Rare codons.** The codons whose scarcity in *E. coli* is well established, the
ones the Rosetta and CodonPlus strains supplement with extra tRNA genes: AGA,
AGG,
CGA, CGG (Arg), ATA (Ile), CTA (Leu), CCC (Pro), GGA (Gly). What the tab
reports
and ranks are **clusters**, not the total: clustered and tandem rare codons are
the ones associated with ribosome stalling, while scattered ones are generally
tolerated. A gene can carry twenty rare codons and be perfectly fine.

**The CAI is real, and it was not always.** An earlier version of this section
refused to report a Codon Adaptation Index, because computing one needs a codon
usage table and inventing sixty-four values from memory would have produced a
precise-looking number that meant nothing. That refusal has been resolved
rather
than relaxed: the app now ships the measured table,
`data/ecoli_codon_usage.csv`,
from the Kazusa Codon Usage Database (taxid 316407), released under CC0 and
redistributed via the `python-codon-tables` package. Every usage figure and the
CAI come from it.

A CAI near 1 means the gene uses the host's preferred codons throughout.
Methionine, tryptophan and stop codons are excluded, having no synonym to
choose
between. It describes codon choice; it does not predict expression level.

The same table also rebuilt the codon preference order the **DNA** tab uses for
back-translation (§3.13). The previous ordering was written from memory and was
right for 18 of the 20 most-used codons, but wrong in the minor ranks of
isoleucine, leucine, arginine and serine.

**GC.** Overall, and as a sliding window. Stretches that stay outside 30–70% for
at least a window length are listed, merged into regions rather than counted as
windows, since windows overlap by one nucleotide and one awkward 60 nt stretch
would otherwise be reported as thirty separate problems. The GC of the first
fifteen codons is called out separately: a GC-rich start can fold the ribosome
binding region and slow initiation, and it is one of the few places where
silent
changes reliably help.

Homopolymer runs of seven or more are flagged too, they complicate both
synthesis
and Sanger reading.

**Suggested codon changes.** A GC-rich stretch is one of the few sequence problems
with a clean fix: swap codons for lower-GC synonyms. When the tab finds one, it
offers the swaps.

Every suggestion is checked against four constraints before it is offered.

The substitute is not one of the documented rare codons. On its own that is not
enough, so there is a second brake, expressed in the only terms that mean
anything: **the substitute must carry at least a given share of its amino acid
in
real *E. coli*.** Without it the algorithm takes the lowest-GC synonym whatever
its usage, leucine's CTG carries 50% of all *E. coli* leucines and TTA 13%, so
it
would happily make that swap. Neither codon is rare; the swap is still a real
drop. That is the trade you did not want.

The floor is 20% by default, chosen by measurement rather than by feel. On a
GC-rich test gene it clears every stretch at a cost of 0.044 of CAI, where a
10%
floor costs 0.080 by allowing exactly the CTG→TTA swap above. Raise it to 40%
and
the CAI stops falling altogether, but four GC-rich stretches survive, which the
*GC-rich stretches left* figure tells you.

The suggestions table shows `usage_from` and `usage_to` for every change:
measured
frequencies, not an ordering. The **CAI** above the table is the same cost
summed
over the whole gene, so the trade-off reads as two numbers rather than as a
judgement call.

The remaining two constraints are unchanged: the substitute It creates none of
the restriction sites you
choose to preserve (the same list the DNA tab uses, and you can shorten it if a
site no longer matters) and it opens no homopolymer run. The protein is then
translated
from the rewritten sequence and compared to the original; if they differ the
output is refused rather than shown.

*Minimal* and *Maximal* decide how far to go. Minimal applies the biggest GC
reductions first and **stops as soon as no GC-rich stretch survives**. Maximal
swaps every codon it can. On a deliberately GC-rich test gene, minimal cleared
all
four stretches with 26 codon changes where maximal used 56, both left the
protein
untouched and introduced no rare codon. Minimal is the default because you are
paying to synthesise this gene, and a construct that differs from the original
in
twenty places is easier to reason about than one that differs in sixty.

You get a **before-and-after GC profile**: red is the sequence as you pasted
it,
blue after the swaps, the grey band is the acceptable range, and the shaded red
area is what used to sit above it. Ticks along the bottom mark the codons that
changed, when they cluster in one place, that is the stretch that was the
problem. Underneath, the before-and-after GC of each stretch as numbers, the
full
list of swaps with the GC each one buys, and the rewritten sequence as FASTA.

**What this does not promise.** Lowering local GC removes one plausible obstacle to
expression. It does not make the gene express, every cause listed at the top of
this tab is still in play.

---

## 8. Sequence library

The sequences you come back to, kept as a file in the repository so they are
one
click away instead of being pasted again every time.

**Read-only from the app, on purpose.** Streamlit Community Cloud rebuilds the
container from the repository on every reboot, redeploy and wake-from-sleep, so
anything the app wrote to disk would be erased, usually on the day you had
forgotten it was fragile. Keeping the library in the repository instead makes
it
durable, versioned, and restorable from git history. The cost is that editing
happens on your own machine rather than in the browser, which is the trade this
app makes deliberately.

**And it is public.** The file lives in a public repository, so it is readable by
anyone whether or not they ever open the app. Nothing unpublished, nothing
under
an NDA, nothing from a collaborator who has not agreed, belongs in it.

### 8.1 What is in it

Twenty-two published sequences, every one taken from the source named in its
note (a PDB entry or KEGG) rather than written from memory. They are chosen
to span the range the app has to cope with:

| | |
|---|---|
| **caplacizumab** (7EOW) | the first nanobody-based medicine approved |
| **cAbLys3** (1MEL) | the structure that first showed a VHH paratope; extra CDR1–CDR3 disulfide, 26-residue CDR3 |
| **Nb80** (3P0G) | stabilises the active state of the β2 adrenergic receptor: the conformational-chaperone use |
| **anti-GFP Enhancer / Minimizer** (3K1K, 3G9A) | the pair that raises and lowers GFP fluorescence; the Minimizer carries an extra disulfide |
| **NbALFA** (6I2G) | the ALFA detection-tag nanobody |
| **Ty1** (6ZXN) | neutralises SARS-CoV-2 |
| **7D12** (4KRL), **2Rs15d** (5MY6) | anti-EGFR and anti-HER2, the latter an imaging tracer |
| **Nb.b201** (5VNV) | from a synthetic library rather than an immunised animal |
| **trastuzumab VH, VL, VH+CH1** (1N8Z) | a conventional VH, a light chain, and a two-domain chain |
| **hen lysozyme** (1MEL) | not an antibody at all: and the antigen cAbLys3 binds |
| **cAbBCII10** (3DWT) | anti-β-lactamase VHH, widely reused as a stable framework scaffold |
| **NbSyn2** (2X6M) | against the α-synuclein C-terminus; extra CDR1–CDR3 disulfide |
| **α-synuclein** (1XQ8) | the antigen NbSyn2 binds: 140 residues, disordered, pI 4.7 |
| **neutrophil elastase, human and mouse** (1HNE, KEGG) | the human entry is the mature chain, the mouse one the full precursor, the notes say so, and they are not comparable until the mouse one is trimmed |
| **serum albumin, human, bovine and mouse** (1AO6, P02769, P07724) | the half-life-extension partner, all three as mature chains so they are directly comparable: 583–585 residues, 66 kDa, 70–76% identical to each other. Human and bovine have 35 cysteines, an odd count, so the app flags the free thiol, which is Cys34 |

Together they make a test set as much as a reference: hallmark scores from 0 to
4,
CDR3 lengths from 8 to 26 residues, three with an extra disulfide, humanness
from
66% to 90%, and one entry for each branch of the classification in §3.1c.

Two are worth looking at closely. **Ty1** is a real published alpaca nanobody
that
scores **1/4** on the hallmark tetrad and carries a VH-like `WVRQ` in FR2, so
the
app calls it ambiguous, it proves on a real case what the hallmark caption
claims,
that the score is indicative and not a classification. **NbALFA** scores 2/4
and is
a widely used reagent. A VHH is not required to be 4/4.

Where a deposited construct carried a tag or expression residues, the entry is
trimmed to the variable domain and the note says exactly what was removed.
Three
entries carry DNA, labelled in the note as **back-translated by this app, not
the
published gene**.

### 8.2 Using it

Filter by free text over names and notes, by tag, or by whether an entry
carries
a protein, DNA, or either. Select the entries you want and send them on: the
proteins go to the **analysis** section, the DNA to **gene diagnosis**. The
buttons say how many of your selection carry each kind, and entries that do not
carry the requested one are left out rather than treated as an error, selecting
a protein-only entry and asking for DNA simply sends nothing for it.

Both kinds also download as FASTA if you want them outside the app.

### 8.3 Editing it

The library is `data/library.csv`. Edit it in a spreadsheet or a text editor,
commit, push, the app serves the new version on its next deploy.

| Column | |
|---|---|
| `name` | what it is called here; must be unique |
| `protein` | the amino-acid sequence, or empty |
| `dna` | the coding sequence, or empty |
| `tags` | comma-separated, used by the filter |
| `note` | free text |

At least one of `protein` or `dna` is required. An entry may carry both, and
when
it does the app checks on load that the DNA translates to the protein and flags
the row if it does not, the only place that check ever happens.

A malformed row is **shown with its problem rather than dropped**: a missing
name, a duplicate, a non-standard residue, a stray character in the DNA. A
library that silently loses an entry is worse than one that tells you which
line
to fix, and the message carries the line number.

---

## 9. In-silico cloning

One screen, three steps.

### 9.1 Recipient vector

Five vectors ship with the app: **pMECS, pMES3, pMES4, pHEN6, pHEN25**. You can
also
upload your own in SnapGene `.dna`, GenBank or FASTA.

SnapGene and GenBank carry their own features, which are read as they are. A
bare
FASTA carries only sequence, so the app detects what it can from short protein
motifs
(PelB, His6, HA, Myc, FLAG, Strep, protease sites, and the VHH scaffold
itself). It
deliberately does **not** guess large elements such as gIII, AmpR or the
origin: doing
that properly needs a reference database, and a confident wrong annotation on a
cloning
tool is worse than none. The **source** column in the feature inventory always
says
whether a feature came from the file or was inferred here.

The **cloning window** is the span between the end of the signal peptide and
the first
tag. That is where the insert will land, and the app shows the coordinates it
found.

The **Expected tags** selector says which cassette layout the finished
construct should
have; it is used when the construct is translated back.

### 9.2 Donor DNA

Two kinds of donor:

**A plasmid that carries the VHH.** The app locates the domain, picks the primer pair
that brackets it, and amplifies it in silico. The amplicon becomes the insert.

**A VHH fragment amplified with primers.** It must already carry the restriction sites
those primers added. A bare coding sequence cannot be cloned by restriction,
and the
app says so rather than failing quietly.

Paste the DNA, or load a file and the sequence is filled in for you.

**Primer selection is automatic.** Each primer in the list is tried as written *and* as
its reverse complement, because a primer is often noted as the plus-strand
sequence of
its target rather than as the oligo that primes back towards the insert. The
shortest
bracketing product wins, and the app tells you which form of which primer it
used.

The bundled list is MP57, GIII, FP and RP, in `data/primers.csv`. You can
upload your
own CSV with `name,sequence` columns; nothing is stored.

### 9.3 The plan

The app cuts, joins and translates **every candidate enzyme pair for real**,
then
compares the VHH the construct encodes with the VHH in the donor. A pair is
only
recommended when both hold:

- the protein survives unchanged, and
- nothing annotated in the vector is deleted.

Among the pairs that pass, the one disturbing the vector least is recommended.
On the
bundled vectors that is usually **PstI + BstEII**, which is the classic
nanobody
strategy.

**The protein check has three outcomes.**

| Status | Meaning |
|---|---|
| **identical** | the construct encodes exactly the donor's VHH |
| **donor part preserved** | the construct carries extra N-terminal residues contributed by the recipient vector, upstream of the 5' restriction site. Everything the donor encodes is unchanged |
| **changed** | a real difference, listed position by position |

The middle case is normal for a bare amplified fragment: it starts at the 5'
restriction
site, a few codons into FR1, so the first residues necessarily come from the
recipient.
The app names the residues added rather than calling it a mismatch.

**Two things worth knowing, because they bite on real vectors.**

*The vectors do not all start the domain the same way.* pMECS begins the VHH with
**QVQLQ**, pHEN25 with **EVQLQ**. Moving a VHH from pMECS into pHEN25 through PstI
therefore changes residue 1 from Q to E, so that pair is rejected and NcoI +
BstEII is
proposed instead, its site sitting upstream of the difference.

*An enzyme pair that translates perfectly can still be wrong.* HindIII + BstEII gives a
readable cassette on pMES3 while removing 120 bp that include the PelB leader.
Pairs
that delete an annotated feature are rejected even when the protein check
passes, and
the table of all pairs says exactly what each rejected pair would have deleted.

**There is deliberately no "insert length must be a multiple of 3" test.** An empty
phage-display vector carries a stuffer, so its downstream tag is usually *out*
of frame
with the leader until a correctly sized VHH is cloned in. Comparing against the
empty
vector's frame would flag every correct construct. Translating the finished
construct
is the only check that means anything.

### 9.4 What you get

The recommended enzyme pair, the construct size, the VHH length, how much
vector is
removed. The donor and construct VHH side by side. The translated cassette
(signal peptide, VHH, tags) so you can read it yourself.

Downloads: the construct as **GenBank**, annotated with the insert and the
cassette
parts so SnapGene or Benchling opens it ready to read; as FASTA; and the
amplicon as
FASTA. Plus a short protocol paragraph.

### 9.5 Folded-away detail

- **All enzyme pairs considered** (every pair with its verdict and, when rejected, why. - **Restriction sites**) every enzyme that cuts the backbone, how often, where, and
  whether it also cuts the insert. An enzyme that cuts the insert is not
necessarily
  wrong: BstEII is used on purpose, its site sitting in FR4 at the 3' end of
the domain.
- **Primers on this vector**: where each binds and in which orientation. Primers with
  no site are listed as not applicable; that is usually a property of the
vector rather
  than a mistake, since pHEN6 has no gIII for a gIII primer to sit in.
- **Vector map** (circular or zoomed on the cloning region, with enzyme sites in red, file features in blue, detected motifs in orange and a VHH scaffold in green. - **Vector features**) the full inventory with its sources.

---


## 10. Reference tables

The two tables you look up rather than work out: the genetic code, and the
twenty amino acids. Reached from **Reference** in the sidebar, under the
manual, because it is something you consult rather than something you run.

Nothing on that screen is typed in, and that is the only reason it is there.
A reference table maintained separately can disagree with the app printed next
to it, and when it does the reader has no way of telling which of the two is
wrong. So the codon table is Biopython's, the usage figures are the measured
*E. coli* table that §7 scores CAI against, the rare-codon flags are the same
eight that §7 raises, and the residue masses, hydropathy and charges are the
same constants §3 builds its profiles from. They cannot contradict each other,
because they are the same numbers.

### 10.1 Genetic code

First the textbook block: rows are the first two bases, columns the third.
Read across a row and the amino acid usually does not change, which is the
wobble position and the reason a synonymous swap is available at all.

Then all 64 codons, filterable by amino acid or down to the rare ones:

- **share** is the codon's share of its own amino acid's codons, which is what
  matters when choosing a synonym.
- **w** is that share divided by the best synonym's, the relative
  adaptiveness; CAI is the geometric mean of these. A codon with w near 1 is
  the one the host prefers.
- **rare_in_ecoli** marks the eight codons §7 flags, each with the tRNA that
  is limiting for it.

Stop codons have no share, because they do not compete with a synonym for an
amino acid, so a share would be a number about nothing.

### 10.2 Amino acids

Twenty structures, drawn as they exist **at pH 7.4**: the zwitterion, with
aspartate and glutamate deprotonated and lysine and arginine protonated.
Histidine is drawn neutral, because its side-chain pKa is near 6 and at 7.4 the
neutral imidazole dominates; it is the one residue whose picture would be
different at a pH you might plausibly work at. A drawing of a charged species
without its pH is a drawing of nothing in particular, which is why the number
is on the screen rather than in a footnote.

The stereo descriptors are the CIP ones. Every amino acid here is L, which is
S at the alpha carbon for all of them except cysteine, where the sulfur
outranks the carboxyl and the same spatial arrangement is called R. Threonine
is 2S,3R and isoleucine 2S,3S.

Below the structures, the properties:

- **residue_da** is what the residue adds to a chain, and is the number every
  molecular weight in the app is built from. **free_da** is the free amino
  acid.
- **hydropathy** is Kyte and Doolittle, the scale behind the hydrophobicity
  profile in §3.
- **charge_profile** is the approximation the charge profile uses, where
  histidine counts as +0.1 rather than being treated as fully protonated or
  fully neutral.
- **side_chain_pka** is for the free amino acid. In a folded protein the local
  environment moves these, sometimes by several units, so they say what to
  expect and not what is the case.

The structures are SVG files committed to `data/aminoacids/`, not drawn when
the page loads. They are produced once by `tools/aminoacids/build.py`, which
checks each one's formal charge and its molecular weight against Biopython's
before writing it and refuses to finish quietly if either is wrong. Doing it
that way keeps RDKit, 66 MB installed, out of the app for the sake of twenty
pictures that never change. The carbon skeleton is drawn in `currentColor` so
it follows the light or dark theme instead of being black in both.

---

## 11. How your data is handled

Sequences, chromatograms and plasmids are **never written to a persistent file,
logged,
or sent to a third party**. They are processed in the server's memory.

Results are held in a short-lived cache (10 minutes, 8 entries for analyses,
with per-item caps elsewhere) which exists only to keep the controls
responsive. It is
cleared when it expires and when the app restarts. Nothing reaches disk; the
one
exception is the alignment step, which writes a temporary FASTA that is deleted
as soon
as the figure is rendered.

**But the app runs on Streamlit Community Cloud, i.e. on third-party infrastructure.**
Traffic is encrypted and nothing is stored, but the data does pass through a
server
that is not yours. For genuinely confidential or unpublished material, run the
app on
your own machine (see the next section) and nothing ever leaves your computer.

Uploads are capped at 10 MB per file, 50 sequences per analysis.

---

## 12. Running it on your own machine

```bash
git clone https://github.com/TyrosinDuck/vhh-analyzer
cd vhh-analyzer
pip install -r requirements.txt
streamlit run app.py
```

It opens on `http://localhost:8501`.

ANARCI is optional. Without it the app silently uses the regex annotator and
says so in
the sidebar; the hallmark and germline features degrade accordingly. To get it:

```bash
conda install -c conda-forge -c bioconda hmmer anarci
```

A `Dockerfile` is included if you would rather run it in a container, or deploy
it
somewhere other than Community Cloud.

---

## 13. Limits and caveats, in one place

- Every physico-chemical value is theoretical. Confirm pI and ε280 experimentally.
- Secondary structure is a propensity estimate, not a prediction.
- CDR boundaries depend on the method; IMGT and the regex anchors give different CDR1
  and CDR2 lengths, so always say which you used.
- The hallmark score is indicative. IMGT 49 is shared with VH3.
- Germline identity is not an immunogenicity prediction.
- Liability motifs and developability flags are triage heuristics with hand-set
  thresholds.
- The codon optimiser ignores mRNA secondary structure, codon-pair bias and tRNA
  availability.
- The cloning assembly is a top-strand design check. It assumes the ends ligate; it
  does not predict efficiency, and it does not simulate the transformation.
- A Sanger read flagged as a frameshift has its VHH withheld on purpose, because a
  reconstructed sequence across an indel would be wrong.
- The structure section measures one static conformer of one model you supply. It
  does not predict structure, dock, estimate affinity or model flexibility.
- On a predicted model, side-chain accessibility is indicative: side chains are the
  part a predictor places least reliably, even at high pLDDT.
- Surface patches use a transparent threshold, not SAP or CamSol, and the default
  exposure cut-off was checked against a single complex.
- A structural epitope says which residues are covered, not which ones the binding
  energy depends on.

---

## 14. Troubleshooting

**Did my change actually deploy?** The footer carries a seven-character
fingerprint of the source files, for example `v1.0.0 (build 20261028.1aed6b5)`.
It changes whenever any `.py` file changes, so if it still reads the same after
a push, Streamlit is serving the old container, not your new code. **Reboot app**
in the Community Cloud dashboard. A build date alone cannot tell you this: two
pushes on the same day carry the same date.

**"Some modules are out of date."** Every file is released together, and each one
carries a version `app.py` checks on startup. The message names the file and
both
versions; copy the whole folder rather than individual files. If it persists
after
a correct update, Streamlit is serving a cached import, **Reboot app** from the
Community Cloud dashboard clears it.

**No CDR annotation, or no germline table.** ANARCI is not installed in this
deployment. The sidebar caption confirms it. The regex annotator still works.

**"anchor not found: …"** on a sequence. The regex annotator could not find the
landmark it names. Either the construct is non-canonical, or the sequence is
truncated.
Check the sequence around that landmark.

**A `.ab1` read comes back as `no PelB`.** Either the wrong vector is selected, the read
is very poor, or the well was empty. Look at the quality trace.

**The cloning section says the fragment carries no sites.** A fragment amplified with
plain primers cannot be cloned by restriction. You need primers whose 5' tails
add the
sites the recipient also has.

**The cloning section rejects every pair.** Widen the enzyme set to *All commercial
enzymes*, check that the insert really carries the sites you expect, or try
another
recipient. The table of all pairs gives the reason for each rejection.

**A result looks wrong.** Please say so, with the sequence and what you expected. That
is the single most useful thing you can send back.

---

## 15. Under the hood

Built on Biopython (ProtParam, Restriction, pairwise alignment), ANARCI for
IMGT
numbering, FAMSA via `pyfamsa` for multiple alignment, `pymsaviz` for the
alignment
figure, `dna_features_viewer` for the plasmid maps, `reportlab` for the PDF,
and
Streamlit for the interface.

The synthetic-sequence generator uses a position-frequency matrix derived from
**PLAbDab-nano** (OPIG, University of Oxford); please cite that resource if you
publish anything built on the generated sequences.

The code is in two layers. The science is in importable modules with no
dependency
on the interface (`vhh_core.py`, `structure.py`, `codons.py`, `cloning.py`,
`sanger.py`, `repertoire.py`, `generator.py`, `constructs.py`, `triage.py`,
`library.py`) so a Jupyter notebook can use them
exactly the same way:

```python
import vhh_core as vc
df, annotations, errors = vc.analyze({"MyNb": "EVQLVESGGG..."})
```

The interface is the second layer: `app.py` routes, `shell.py` holds the
chrome,
`compute.py` holds every cached call, and `sections/` holds one module per
screen, each a `render()` function that draws and computes nothing. `app.py`
was a
single 3 911-line script before those screens were pulled out of it; the README
lists the files and the two rules that keep it from growing back.

There is a test suite: `pip install -r requirements-dev.txt` then `pytest`.
Almost
every test in it is a bug the app actually had, frozen so it cannot return.

The source is at <https://github.com/TyrosinDuck/vhh-analyzer>, under the **MIT
licence**: use it, change it, redistribute it, build on it, keep the copyright
notice. One caveat is spelled out in `LICENSE`: that covers the app's own code
and cannot change its dependencies' terms, and one of them, `pyfamsa` (FAMSA,
used for the multiple alignment), is GPL-3.0-or-later.

If you use the app in published work, cite it by its DOI:
**10.5281/zenodo.23742269**. That is the *concept* DOI, it always resolves to
the latest archived release. The home page has a **Citing this tool** panel
with
the citation and a BibTeX entry ready to copy, and each release has a version
DOI
of its own on Zenodo for saying exactly which code produced a result.
