# VHH tool box — user manual

A browser-based workbench for single-domain antibodies. Nothing to install, no code
to write: you paste sequences or drop files, and it gives you characterization,
sequencing read-back and cloning plans.

This manual describes every screen and, just as importantly, what each number means
and how far you should trust it.

---

## 1. What it does, and what it does not

**It does** compute descriptors from an amino-acid sequence, annotate CDRs and
frameworks, read VHHs out of Sanger chromatograms, and plan a restriction cloning
between two plasmids.

**It does not** predict structure, affinity, expression yield or immunogenicity.
Every value it reports is derived from sequence alone. Treat the output as triage —
it tells you where to look, not what will happen at the bench.

Three things are worth stating plainly before you rely on any number:

- The physico-chemical values are **theoretical**. Confirm the isoelectric point by
  IEF or cIEF and the extinction coefficient against a measured A280 before using
  them quantitatively.
- The secondary-structure output is a **composition-based propensity**, not a
  structure prediction. For real secondary structure, run DSSP on a model.
- The cloning section is a **design check**. It verifies that your construct would
  encode what you intend; it says nothing about ligation efficiency or whether the
  protein will express.

---

## 2. Getting around

The app opens on a home page with three sections. The sidebar carries the same three
entries plus a **Home** button, available from every screen.

| Section | Use it for |
|---|---|
| 🔬 **Sequence analysis** | characterizing VHH sequences you already have |
| 🧪 **Sanger sequencing** | pulling VHH sequences out of `.ab1` files |
| 🧬 **In-silico cloning** | moving a VHH from one vector into another |

The three share one sequence list. Whatever is loaded in the analysis input is also
what the cloning section offers as an insert, and the sequencing section writes into
it. Each section still works on its own.

---

## 3. Sequence analysis

### 3.1 Giving it sequences

Three ways in:

- **Paste** into the text box: FASTA, or a bare sequence, or several bare sequences
  separated by blank lines. Whitespace, digits and gap characters are stripped.
- **Upload a FASTA file** with the uploader on the right.
- **Load example** fills the box with three VHHs (two generic, plus cAbLys3) so you
  can see what every tab looks like before committing your own data.

Click **Analyze**, or press ⌘ + Enter (Ctrl + Enter on Windows and Linux). The button
exists because the keyboard shortcut is not discoverable; both do the same thing,
which is to commit what you typed.

**Limits.** Up to 50 sequences per run, each between 60 and 400 amino acids. Anything
outside that is rejected with a message naming the offending sequence — a sequence
longer than 400 aa is almost certainly not a single domain, and the app would rather
say so than silently produce nonsense.

A sequence containing a non-standard residue (B, X, Z, a stop…) is rejected by name.
Fix or remove it rather than letting it through.

### 3.1b Generating synthetic VHHs

The panel **Generate synthetic VHH sequences** builds sequences that look like real
VHHs without being any of them. Useful for negative controls, for benchmarking a
pipeline, or for testing this app itself.

**How it works.** Frameworks are sampled position by position from an IMGT
position-frequency matrix built from real VHHs. CDRs are drawn at random, with
lengths taken from the same reference set.

**The reference set.** The bundled matrix comes from **PLAbDab-nano** (OPIG,
University of Oxford): 1958 unique VHH sequences, numbered with ANARCI against the
alpaca germlines, giving 212 IMGT positions with a median support of 1890 sequences
per position. The panel states this, along with the support, because a frequency is
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
| Fix the CDR lengths | otherwise lengths come from the reference distribution — CDR1 and CDR2 around 8, CDR3 a mean of 17 over a range of 2 to 30 |
| Random seed | the same seed always gives the same set, so a benchmark is reproducible |

Generated sequences download as FASTA, or go straight into the analysis above.

**Three things to know before using the output.**

*Positions are sampled independently.* Real frameworks have correlated positions — a
residue at one site constrains another — so a combination that never co-occurs in
nature can come out. The hallmark tetrad is the usual casualty, which is why it is
enforced by default.

*Cysteines in frameworks are off by default.* About 2% of the reference sequences
carry a cysteine at IMGT 50, where it pairs with a cysteine in a CDR. With randomly
generated CDRs that partner does not exist, so such a residue would be left as an
orphan thiol.

*These are synthetic sequences.* Nothing here predicts that they fold, express or bind
anything. They are statistically plausible strings, and that is all.

A check worth running: generate a set and send it through the analysis. On the bundled
matrix every generated sequence annotates cleanly and scores 4/4 on the hallmarks,
with molecular weights around 13.5–14 kDa — which is the sanity test that the matrix
and the assembly are behaving.

### 3.2 Sidebar settings

| Setting | Default | What it changes |
|---|---|---|
| pH for net charge | 7.4 | the pH at which net charge and the developability flags are evaluated |
| Profile window | 9 residues | the sliding window for the per-residue charge and hydropathy profiles |
| CDR annotation method | auto | `auto` uses ANARCI when it is installed and falls back to the regex anchors; `regex` forces the dependency-free path; `anarci` appears only when ANARCI is actually available |

The caption under the selector tells you which annotator is live. If it says *not
installed*, every CDR boundary you see comes from the regex anchors.

### 3.3 Properties

The main table, one row per sequence.

| Column | Meaning |
|---|---|
| Length, MW | a VHH is typically 110–135 aa and 12–15 kDa |
| Theoretical pI | pH of zero net charge, Bjellqvist pKa set. A pI close to your buffer pH predicts poor colloidal stability |
| ε280 | molar extinction coefficient at 280 nm, reduced cysteines |
| A280 (1 mg/mL) | ε / MW — **this is the number you divide your measured A280 by to get mg/mL** |
| Instability index | Guruprasad dipeptide score; below 40 predicts a stable protein |
| Aliphatic index | relative volume of A, V, I, L; higher correlates with thermostability |
| GRAVY | mean Kyte–Doolittle hydropathy; negative means hydrophilic. Above about −0.2 flags an aggregation-prone candidate |
| Aromaticity | fraction of F, W, Y; this is what drives the A280 signal |
| Charge | net charge at the working pH |
| n_Cys | cysteine count; two is the canonical intradomain disulfide |

GRAVY and the instability index are shaded green-to-red so the outliers in a set jump
out.

**Developability flags** below the table turn those numbers into sentences. A flag is
a threshold crossing, not a verdict: instability index ≥ 40, GRAVY above −0.2, net
charge magnitude below 2, pI within one unit of the working pH, more than two
cysteines, three or more methionines. Tune the thresholds to your own data in
`vhh_core.stability_flags` if they do not match your experience.

**Custom features** is a second table carrying descriptors that are easy to extend:
charge at pH 5 (relevant for cation-exchange capture), the Lys+Arg fraction, CDR3
length, CDR3-only GRAVY, and counts of NG and DP/DG motifs. Each definition is listed
in the expander underneath.

### 3.4 Plots

Four figures:

1. **Secondary-structure propensity** — stacked bars of helix/turn/sheet fractions.
2. **Titration curves and net charge** — each curve crosses zero at that sequence's
   theoretical pI, and the dotted line marks your working pH.
3. **Per-residue profiles** — sliding-window charge and Kyte–Doolittle hydropathy
   along the sequence, with the CDRs shaded. This is where charged patches and
   hydrophobic stretches become visible.
4. **Predicted secondary structure (Chou–Fasman)** — helix and sheet propensity
   curves with a predicted-state track underneath.

The Chou–Fasman heuristic is roughly 50–60% accurate. It shows where the sequence
*tends* towards each conformation; it does not predict the fold. A VHH is an
immunoglobulin β-sandwich, so a sheet-dominated profile is the expected result, and
getting one is not evidence that anything was predicted correctly.

With more than one sequence, only the selected Chou–Fasman panel is drawn. Tick
**Show all** to see every one — with many sequences this is the slowest thing the app
does.

### 3.5 CDR / Framework

A table of CDR1, CDR2 and CDR3 with their lengths, the coloured sequences, and the
framework regions in an expander.

**Two annotation methods, and they do not agree by design.**

*ANARCI (IMGT numbering)* aligns the sequence to germline HMMs, gives every residue an
IMGT position, and takes the CDRs as the fixed windows 27–38, 56–65 and 105–117.
Because it is numbering-based it tolerates unusual CDR lengths. It needs HMMER.

*Regex anchors* needs nothing. It locates the conserved landmarks of a VHH and cuts
between them:

| Anchor | Pattern | Role |
|---|---|---|
| conserved Cys23 | `C` in residues 16–30 | ends FR1 |
| FR2 | `W[A-Z]RQ` — the `WFRQ`/`WYRQ` hallmark | starts FR2 |
| FR3 | `[KR][FYW][TASV][IVLFM]S` — the `RFTIS` motif | ends CDR2 |
| second Cys | `[YFHV][YFHVAC]C` — the `YYC` motif | ends FR3 |
| FR4 | `WG[A-Z]G[TAS][A-Z]{1,2}VTVSS` | ends CDR3 |

**The two give different CDR1 and CDR2 lengths.** IMGT windows are tighter than the
anchor-based cuts, which follow a more Kabat-like logic. CDR3 usually agrees. This
matters when you compare CDR lengths across datasets: say which method produced them.

Where both run, a disagreement beyond that usually means a non-canonical framework,
which is itself worth knowing. The regex method assumes a canonical VHH scaffold; if
an anchor is missing the app names the one it could not find rather than guessing.

### 3.6 Hallmarks & germline

**The hallmark tetrad** is the four FR2 positions that distinguish a camelid VHH from
a conventional VH. The table gives IMGT and Kabat numbering, the residue found, the
VH3 consensus, and a verdict.

| IMGT | Kabat | VH3 | Camelid |
|---|---|---|---|
| 42 | 37 | V | F, H, L, Y |
| 49 | 44 | G | D, E, G, Q |
| 50 | 45 | L | C, L, P, R |
| 52 | 47 | W | F, G, L, S, W |

A score of 4/4 is the canonical camelid profile. **IMGT 49 is shared with VH3, so 3/4
with a G at 49 is common and still a genuine VHH.** Engineered or humanised constructs
sit between the two profiles. Read the score as an indication, never as a
classification.

**Closest germline** runs the assignment twice, once against alpaca and once against
human, and reports the closest V and J gene with percent identity. The identity to the
closest human germline is a first-pass humanness indicator — a sequence can be
germline-like and still contain a T-cell epitope, so this is not an immunogenicity
prediction.

Both need ANARCI. Without it, the hallmarks fall back to reading the positions off the
regex FR2 block, and the germline table is absent.

**Humanness, position by position.** Underneath the germline table, percent
identity to the closest human germline V gene — computed **over framework
positions only**, because CDR divergence from a human germline is expected and is
not something you would ever try to fix.

The identity is computed position by position in IMGT numbering against the
IMGT-gapped germline sequences that ship with ANARCI, not by a fresh alignment.
On a typical VHH it agrees with ANARCI's own independent germline assignment to
within 0.1%, which is a useful cross-check that the numbering is right.

What makes it actionable is the divergence table: one row per framework position
where the VHH differs from the closest human germline, with the IMGT number, the
region, both residues, whether the substitution is conservative, and whether the
position is one of the four hallmarks.

That last flag is the whole point. FR2 typically comes out around 65% human while
FR1 and FR3 sit near 90% — not because FR2 is unusual, but because FR2 is where
the hallmarks live. Reverting them raises the humanness number and simultaneously
removes the hydrophilic camelid surface that lets a VHH stay soluble without a
light chain. The other divergent positions are the cheap ones; the hallmarks are
a real trade-off, and the table separates them so you can see which is which.

**A limit worth repeating.** This is a first-pass humanness indicator and nothing
more. It is not an immunogenicity prediction: immunogenicity depends on T-cell
epitopes, dose, route and formulation, none of which is in the sequence. A
germline-like sequence can carry an epitope; a divergent one can be silent.

### 3.7 Liabilities

Chemical and manufacturability liabilities, **scored** and placed along the sequence.

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

Risk categories follow the same scale: **Low** below 2, **Medium** below 6, **High**
at 6 or above.

Three further motifs — fragmentation (`DP`), integrin binding (`RGD`) and hydrophobic
patches of four or more consecutive F/I/L/V/W/Y — are **additions**, marked as such in
the `source` column. Untick *Include the motifs added on top of the in-house scale* to
reproduce the original totals exactly.

**An important limitation.** In its original form each motif is also weighted by its
**solvent exposure** (SASA) and by the **model confidence** (pLDDT), both taken from a
3D structure. This section is sequence-only, so neither factor is applied: an `NG`
scores 10 whether it is buried in the core or fully exposed on a loop. Two consequences
follow, and they matter:

- the categories **compress towards the top** — almost everything lands in Medium or
  High, where the structure-weighted version would push buried motifs down to Low;
- **comparing sequences with each other tells you more than the category of any single
  motif**. Use the totals to rank a panel, not to declare one motif dangerous.

If you have an AlphaFold model, the structure-weighted analysis is a separate notebook
and gives a far sharper answer.

**The liability map** places every motif at its residue position, one row per sequence,
with bubble size proportional to the score and the CDRs shaded behind. This is the view
that answers the question you actually ask of a panel: *where* do the problems sit, and
do they fall in a binding loop. The **liability profile** underneath stacks each
sequence's score by motif type, with the total on top.

**The region column is what makes the table actionable.** The same `NG` means different
things in FR3 and in CDR3: in a framework it is mainly a stability and heterogeneity
concern, in a CDR it may alter binding.

**Free cysteines.** With only the sequence, the app assumes the canonical intradomain
disulfide joins the FR1 and FR3 cysteines, and flags any Cys outside those regions as a
candidate free thiol. cAbLys3, which carries a genuine second disulfide between CDR1
and CDR3, is therefore flagged twice — correctly in the sense that those cysteines are
extra, and misleadingly in the sense that they are paired. With a structure you would
measure the SG–SG distances instead.

**Protease sites** sit in an expander: trypsin-like (after K/R, not before P),
chymotrypsin-like (after F/W/Y, not before P), GluC-like, AspN-like, and the
high-scoring Factor Xa and enterokinase recognition sites. Every protein has tryptic
sites, so the count on its own means little. What is worth looking at is a high-scoring
recognition site, or a cluster of sites inside a CDR.

Everything downloads as a single CSV.

**Cysteines and disulfides.** The liability table flags every cysteine outside
the canonical FR1/FR3 pair as a possible free thiol. That is the right default and
it is wrong in one common case: camelid VHHs with long CDR3s frequently carry a
**second, stabilising** disulfide, most often between CDR1 or FR2 and CDR3. Read
as a liability, a stabiliser becomes a problem to fix — which would be exactly
backwards.

This section tells the two apart. It numbers every cysteine in IMGT, identifies
the canonical IMGT 23/104 pair, and then looks at what is left: two extra
cysteines in regions that pair in real VHHs are reported as a plausible additional
disulfide, one extra cysteine as an orphan thiol, three or more as something that
cannot be resolved from sequence.

Pairing here is inferred from **position alone**. There is no structure in this
app, so "plausible pair" means two cysteines sitting where real VHHs pair them,
not a verified bond. A non-reducing gel or a structure is what confirms it.

**Suggested repairs.** One conservative substitution per liability, so the table
stops being merely descriptive: NG fixed at the glycine rather than the
asparagine, DG likewise, sequons broken at the +2 residue, M→L, free Cys→Ser, and
so on. Each row carries its rationale.

Three behaviours are worth knowing about.

*Conserved positions are refused.* The oxidation rules match every Met and Trp in
the sequence, which includes the conserved FR2 tryptophan (IMGT 41) and the FR4
tryptophan of the WGxG motif (IMGT 118). Substituting either would "repair" a
liability by breaking the fold, so those rows appear in a separate table with the
reason, and no substitution is offered. The canonical disulfide cysteines are
treated the same way.

*CDR positions are flagged, not hidden.* A liability in a CDR is often one you
keep — fixing it may cost the affinity you selected for. The rows say so, and the
default selection is the framework-only subset.

*Substitutions that do not match are refused rather than applied.* Each one is
written `N55Q`, and if residue 55 is not an asparagine the substitution is
rejected with a message instead of silently doing something else. Two conflicting
substitutions at the same position are refused as a pair, because that conflict is
your decision.

Applying a selection shows the repaired sequence, the new liability score with the
delta, and a FASTA download. The recomputed score reuses the original region
boundaries, which is correct for substitutions; re-run the full analysis on the
repaired sequence before trusting anything else about it.

### 3.8 Aggregation-prone regions

Where along the sequence the aggregation-prone stretches are. A stretch becomes
aggregation-prone when it is hydrophobic, β-prone and carries no charged residue
to interrupt it, and the profile combines exactly those quantities: the absolute
charge per residue, Kyte–Doolittle hydropathy, Chou–Fasman β-sheet propensity, and
a gatekeeper term for P/D/E/K/R. Positive leans soluble, negative leans
aggregation-prone.

**This is not CamSol, AGGRESCAN or TANGO.** Those fit their coefficients to
experimental data. The weights here are a transparent judgement call and sit in
`SOLUBILITY_WEIGHTS` in `vhh_core.py` for you to change. The consequence is
specific: read **where** the weak spots are, and compare candidates against each
other, but do not read the absolute value as a predicted solubility. For a
validated number, run the sequence through one of those servers.

With that said, the profile behaves the way it should on controls. A poly-isoleucine
stretch scores about −1.0, a poly-lysine stretch about +2.0, and a well-behaved VHH
comes out around +0.3 with its weakest spots on the known hydrophobic framework
segments. Give it a deliberately hydrophobic CDR3 and CDR3 goes straight to the
top of the list.

**Two controls.** The *smoothing window* (default 9) sets how long a stretch has
to be to register; the *threshold* (default −0.3) sets what counts as
aggregation-prone.

**Two tables.** *Regions below the threshold* is the headline, but a good VHH
often has none at all — which would leave the tab empty and tell you nothing. So
*weakest windows* always lists the three lowest-scoring windows whether or not
they cross the threshold, with their region and whether they fall in a CDR. That
is the table to compare candidates on.

In the summary, `Mean` is higher-is-better and `Worst` is the single lowest
window. `Worst` is usually what decides whether a candidate behaves: one bad patch
is enough.

### 3.9 Alignment

With two or more sequences: a FAMSA alignment rendered with `pymsaviz`, with the CDR
bands of the first sequence marked, downloadable as aligned FASTA.

Underneath, **two identity matrices**: full-length and CDR3-only. The second is the
one that matters when triaging panning output. Two clones can share 90%+ full-length
identity while having unrelated CDR3s — they are the same scaffold, not the same
binder. To turn that observation into actual groups, use **Clonal families** (§3.10).

### 3.10 Clonal families

This is the tab for triaging panning output. The first question after a screening is
not how good each clone is, but **how many different clones you actually have** — a
plate of picked colonies is routinely the same few lineages over and over, and
properties computed on 16 reads that represent 3 lineages will mislead you three
times over.

Grouping is done on **CDR3**, because CDR3 is what differs between independent
lineages while the frameworks stay nearly constant.

**How it works.** Every pair of CDR3s is aligned (global, BLOSUM62, gaps −11/−1) and
scored for identity — strictly identical residues over the length of the longer
CDR3. Those identities become a distance matrix (`1 − identity`), which is clustered
hierarchically and cut at your threshold. This is the same procedure as the Sanger
screening notebook, so the app and the notebook give the same families for the same
reads.

**The three controls.**

*CDR3 identity threshold*, 50–100%, default 80%. Two sequences at or above it belong
to the same family. At 100% only strictly identical CDR3s group together.

*Linkage.* `complete` guarantees that **every pair** inside a family meets the
threshold — the conservative reading of "same family". `average` merges on the mean
distance and is more permissive, so a family can then contain pairs that fall below
your threshold. The family table's `min_identity` column is there to catch exactly
that: it reports each family's worst internal pair, and a value under your threshold
tells you `average` merged across it.

*Same CDR3 length only* forces CDR3s of different lengths into different families
whatever their alignment. Off by default, because a one-residue indel in a CDR3 is a
plausible sibling rather than a new lineage. Turn it on if you want strictly
length-matched families.

**What you get.** Four numbers at the top — sequences grouped, families, singletons,
and redundancy (sequences per family). Redundancy above roughly 2 means you are
re-picking the same lineages and could screen more colonies for the same effort.

Then a family table, with for each family its size, how many distinct CDR3s it
contains, its worst internal identity and similarity, its CDR3 length range, and a
**representative** — the family's most frequent CDR3, a reasonable default for the
clone to carry forward.

A dendrogram draws your threshold as a dashed line; branches joining to its right
are separate families. Move the slider and watch which groups merge: that is the
honest way to judge how robust a family is. A grouping that holds from 70% to 95% is
real; one that appears only at 78% is an artefact of where you put the line.

Finally a per-sequence assignment table, and CSV downloads of the families, the
assignment and the CDR3 identity matrix, plus a FASTA of one representative per
family — which is usually what you actually want to order or subclone.

**Two honest limits.** Families are a *hypothesis* about shared lineage drawn from
CDR3 similarity alone; two unrelated lineages can converge on a similar CDR3, and
true siblings with a heavily mutated CDR3 can be split. And sequences whose CDR3
could not be annotated are left out rather than guessed at — they are listed
explicitly when that happens.

### 3.11 DNA

Back-translation into *E. coli* codon-optimised DNA.

Each residue first takes the most frequent *E. coli* codon. The sequence is then
repaired wherever it contains a restriction site you asked to avoid or a homopolymer
run longer than six, by swapping the offending codon for the next synonymous one.
Finally the GC content is balanced, because optimising codon frequency alone pushes GC
to about 63% and synthesis vendors flag anything above 60%.

Controls: which restriction sites to avoid (17 available, all selected by default),
whether to append a TAA stop, and the maximum GC fraction.

**Every output is verified by translating it back and comparing to the input.** If the
round trip fails the app says so in red and tells you not to order that sequence. For
DNA you are going to pay to have synthesised, that check is worth its cost.

The optimiser maximises codon frequency, balances GC and removes cloning-unfriendly
motifs. It does **not** model mRNA secondary structure near the ribosome binding site,
codon-pair bias or tRNA availability. Run the output through your synthesis vendor's
own checker before ordering.

### 3.12 Constructs

Assemble a VHH into the construct you will actually order, and get the mass, pI
and DNA of the finished protein rather than of the bare domain. The practical
reason: a VHH is 13 kDa with a pI near 8.6, and the same VHH with a FLAG tag is a
different protein on a gel and on an ion-exchange column. Four aspartates will
move the pI by two units, and you would rather know that before ordering than
after the purification misbehaves.

**Formats.** *Monovalent* is one VHH with whatever tags you add. *Bivalent* is the
same VHH twice, which raises apparent affinity through avidity, not affinity.
*Biparatopic* is two different VHHs against the same target, which needs
non-overlapping epitopes — picking the same VHH twice there makes a bivalent
construct and the app says so.

**Linkers.** A list from a direct fusion up to (GGGGS)4, plus rigid (EAAAK)
helices for when you want fixed spacing rather than a floppy tether. (GGGGS)3 is
the default for a bivalent construct — long enough for both domains to bind
independently. Fusing two VHHs with no linker raises a warning, because the
domains will interfere and that is almost never what you want.

**Tags.** His6, His8, HA, c-Myc, FLAG, Strep-II, Avi, the original SpyTag, and TEV,
thrombin and sortase sites. Each carries a short note. The sortase motif has to be
last — LPETG is cleaved between T and G, so anything after it is lost — and the app
warns if you put something after it.

**What is deliberately missing, and why.** SpyTag2 and SpyTag3 are not in the
list, and neither are the camelid hinges. Only the original SpyTag is included.
Writing the others from memory would risk handing you a sequence that is one
residue wrong, which you would discover only after the synthesised gene arrived.
The *Custom pieces* expander takes any peptide you paste, with a label, at either
end — that is the honest way to handle a sequence this app cannot vouch for, and
it is where those belong.

**Verify everything anyway.** Tag sequences differ between vectors and between
papers. The app says so on every run, and it means it: a one-residue difference
will not show up until the protein does not work.

**What you get.** A schematic drawn to scale with VHHs, linkers and tags
distinguished; length, mass, pI and A280 each with the delta against the bare VHH;
a segment table giving every piece's exact coordinates; the full property set; the
protein sequence; and codon-optimised DNA for the whole construct, with the same
restriction-site avoidance as the DNA tab and the same round-trip verification —
the DNA is translated back and compared to the construct before it is shown.

### 3.13 Download

The results table as CSV, the figures as 300 dpi PNG, and a full PDF report.

The report is built on demand because it takes a few seconds. It contains the
physico-chemical summary, then one section per sequence (properties, developability
flags, regions, the sequence, hallmarks, germline, liabilities), all the figures, and
a page of caveats.

The **Include the per-sequence secondary-structure plots** checkbox adds one page per
sequence. It is ticked automatically up to six sequences; beyond that those plots are
most of the file size and most of the build time.

---

## 4. Sanger sequencing

### 4.1 What it does

You upload the chromatograms from a screening run. For each read the app:

1. quality-trims with **Mott's algorithm** — the one phred and Staden use;
2. translates the read and its reverse complement in three frames each;
3. locates the **PelB leader** to fix the start and the reading frame;
4. locates the **FR4 signature** to fix the end of the VHH;
5. extracts the VHH protein and DNA, and reports what went wrong where it did.

Because all six frames are searched, forward and reverse reads are handled without
being told which is which.

### 4.2 Settings

| Setting | Default | Effect |
|---|---|---|
| Vector | pMES3 | which tags are expected after the VHH: pMES3 is PelB–VHH–His6–stop, pMECS is PelB–VHH–HA–His6–stop |
| Trim cutoff | 0.05 | Mott error-probability cutoff, roughly Q13. Lower is stricter |
| Min trimmed length | 150 nt | below this the read is flagged, and the raw sequence is analysed instead of the trimmed one |
| Min mean Phred | 20 | QC threshold on the trimmed window |

To add a vector, edit `sanger.VECTORS` — and take the tag layout from the actual
plasmid map rather than from memory. A wrong layout makes the app report missing tags
on clones that are perfectly fine.

### 4.3 Reading the results table

Four counters at the top: reads, VHHs extracted, unique sequences, failures.

The table carries the status, QC pass, trimmed length, mean Q, the strand and frame
where the cassette was found, the VHH length, the tags found and missing, premature
stops, the copy number, and which read a duplicate duplicates.

Every warning is also reprinted in full underneath the table, grouped by read, because
a long warning is unreadable in a table cell.

**A failed read is not necessarily a failed clone.**

| Status | Usually means |
|---|---|
| `ok` | cassette found, VHH extracted, no frameshift or internal stop |
| `low quality` | the VHH was extracted but the read failed the QC thresholds |
| `no FR4` | the read stopped before the end of the VHH — sequence further, or from the other end |
| `frameshift` | PelB and FR4 sit in **different reading frames**, so an indel lies between them. The VHH cannot be reconstructed reliably and is deliberately not reported |
| `no PelB` | the leader was not found in any frame — wrong vector, very poor read, or empty well |
| `unreadable` | the file could not be parsed as `.ab1` |

One warning deserves a note: *only the PelB core anchor matched — read starts inside
the leader*. The short anchor `AAQPAMA` is used when the full leader is not present.
That is normal for a read that begins part-way through the leader, and the VHH start is
derived from the anchor that actually matched.

**Identical protein sequences are grouped.** The same colony picked twice shows up as
one sequence with a copy count, and the duplicates are marked.

### 4.4 Quality traces

An expander draws the Phred quality along each read, with the window kept after
trimming in blue and a dashed line at Phred 20. Up to twelve reads are drawn.

### 4.5 Sending clones to the analysis

Below the table, pick the clones worth keeping and press **Send selected to analysis**.
The app switches to the analysis section with those sequences loaded.

**Duplicates are deselected by default** — one entry per unique protein sequence. Tick
them back on if you want to see them separately.

### 4.6 Downloads

The full per-read table as CSV, and the extracted VHHs as protein or DNA FASTA.

---

## 5. In-silico cloning

One screen, three steps.

### 5.1 Recipient vector

Five vectors ship with the app: **pMECS, pMES3, pMES4, pHEN6, pHEN25**. You can also
upload your own in SnapGene `.dna`, GenBank or FASTA.

SnapGene and GenBank carry their own features, which are read as they are. A bare
FASTA carries only sequence, so the app detects what it can from short protein motifs
(PelB, His6, HA, Myc, FLAG, Strep, protease sites, and the VHH scaffold itself). It
deliberately does **not** guess large elements such as gIII, AmpR or the origin: doing
that properly needs a reference database, and a confident wrong annotation on a cloning
tool is worse than none. The **source** column in the feature inventory always says
whether a feature came from the file or was inferred here.

The **cloning window** is the span between the end of the signal peptide and the first
tag. That is where the insert will land, and the app shows the coordinates it found.

The **Expected tags** selector says which cassette layout the finished construct should
have; it is used when the construct is translated back.

### 5.2 Donor DNA

Two kinds of donor:

**A plasmid that carries the VHH.** The app locates the domain, picks the primer pair
that brackets it, and amplifies it in silico. The amplicon becomes the insert.

**A VHH fragment amplified with primers.** It must already carry the restriction sites
those primers added. A bare coding sequence cannot be cloned by restriction, and the
app says so rather than failing quietly.

Paste the DNA, or load a file and the sequence is filled in for you.

**Primer selection is automatic.** Each primer in the list is tried as written *and* as
its reverse complement, because a primer is often noted as the plus-strand sequence of
its target rather than as the oligo that primes back towards the insert. The shortest
bracketing product wins, and the app tells you which form of which primer it used.

The bundled list is MP57, GIII, FP and RP, in `data/primers.csv`. You can upload your
own CSV with `name,sequence` columns; nothing is stored.

### 5.3 The plan

The app cuts, joins and translates **every candidate enzyme pair for real**, then
compares the VHH the construct encodes with the VHH in the donor. A pair is only
recommended when both hold:

- the protein survives unchanged, and
- nothing annotated in the vector is deleted.

Among the pairs that pass, the one disturbing the vector least is recommended. On the
bundled vectors that is usually **PstI + BstEII**, which is the classic nanobody
strategy.

**The protein check has three outcomes.**

| Status | Meaning |
|---|---|
| **identical** | the construct encodes exactly the donor's VHH |
| **donor part preserved** | the construct carries extra N-terminal residues contributed by the recipient vector, upstream of the 5' restriction site. Everything the donor encodes is unchanged |
| **changed** | a real difference, listed position by position |

The middle case is normal for a bare amplified fragment: it starts at the 5' restriction
site, a few codons into FR1, so the first residues necessarily come from the recipient.
The app names the residues added rather than calling it a mismatch.

**Two things worth knowing, because they bite on real vectors.**

*The vectors do not all start the domain the same way.* pMECS begins the VHH with
**QVQLQ**, pHEN25 with **EVQLQ**. Moving a VHH from pMECS into pHEN25 through PstI
therefore changes residue 1 from Q to E — so that pair is rejected and NcoI + BstEII is
proposed instead, its site sitting upstream of the difference.

*An enzyme pair that translates perfectly can still be wrong.* HindIII + BstEII gives a
readable cassette on pMES3 while removing 120 bp that include the PelB leader. Pairs
that delete an annotated feature are rejected even when the protein check passes, and
the table of all pairs says exactly what each rejected pair would have deleted.

**There is deliberately no "insert length must be a multiple of 3" test.** An empty
phage-display vector carries a stuffer, so its downstream tag is usually *out* of frame
with the leader until a correctly sized VHH is cloned in. Comparing against the empty
vector's frame would flag every correct construct. Translating the finished construct
is the only check that means anything.

### 5.4 What you get

The recommended enzyme pair, the construct size, the VHH length, how much vector is
removed. The donor and construct VHH side by side. The translated cassette — signal
peptide, VHH, tags — so you can read it yourself.

Downloads: the construct as **GenBank**, annotated with the insert and the cassette
parts so SnapGene or Benchling opens it ready to read; as FASTA; and the amplicon as
FASTA. Plus a short protocol paragraph.

### 5.5 Folded-away detail

- **All enzyme pairs considered** — every pair with its verdict and, when rejected, why.
- **Restriction sites** — every enzyme that cuts the backbone, how often, where, and
  whether it also cuts the insert. An enzyme that cuts the insert is not necessarily
  wrong: BstEII is used on purpose, its site sitting in FR4 at the 3' end of the domain.
- **Primers on this vector** — where each binds and in which orientation. Primers with
  no site are listed as not applicable; that is usually a property of the vector rather
  than a mistake, since pHEN6 has no gIII for a gIII primer to sit in.
- **Vector map** — circular or zoomed on the cloning region, with enzyme sites in red,
  file features in blue, detected motifs in orange and a VHH scaffold in green.
- **Vector features** — the full inventory with its sources.

---

## 6. How your data is handled

Sequences, chromatograms and plasmids are **never written to a persistent file, logged,
or sent to a third party**. They are processed in the server's memory.

Results are held in a short-lived cache — 10 minutes, 8 entries for analyses, with
per-item caps elsewhere — which exists only to keep the controls responsive. It is
cleared when it expires and when the app restarts. Nothing reaches disk; the one
exception is the alignment step, which writes a temporary FASTA that is deleted as soon
as the figure is rendered.

**But the app runs on Streamlit Community Cloud, i.e. on third-party infrastructure.**
Traffic is encrypted and nothing is stored, but the data does pass through a server
that is not yours. For genuinely confidential or unpublished material, run the app on
your own machine — see the next section — and nothing ever leaves your computer.

Uploads are capped at 10 MB per file, 50 sequences per analysis.

---

## 7. Running it on your own machine

```bash
git clone https://github.com/TyrosinDuck/vhh-analyzer
cd vhh-analyzer
pip install -r requirements.txt
streamlit run app.py
```

It opens on `http://localhost:8501`.

ANARCI is optional. Without it the app silently uses the regex annotator and says so in
the sidebar; the hallmark and germline features degrade accordingly. To get it:

```bash
conda install -c conda-forge -c bioconda hmmer anarci
```

A `Dockerfile` is included if you would rather run it in a container, or deploy it
somewhere other than Community Cloud.

---

## 8. Limits and caveats, in one place

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

---

## 9. Troubleshooting

**"Some modules are out of date."** `app.py` and the analysis modules are released
together. Update `vhh_core.py`, `sanger.py` and `cloning.py` alongside it.

**No CDR annotation, or no germline table.** ANARCI is not installed in this
deployment. The sidebar caption confirms it. The regex annotator still works.

**"anchor not found: …"** on a sequence. The regex annotator could not find the
landmark it names. Either the construct is non-canonical, or the sequence is truncated.
Check the sequence around that landmark.

**A `.ab1` read comes back as `no PelB`.** Either the wrong vector is selected, the read
is very poor, or the well was empty. Look at the quality trace.

**The cloning section says the fragment carries no sites.** A fragment amplified with
plain primers cannot be cloned by restriction. You need primers whose 5' tails add the
sites the recipient also has.

**The cloning section rejects every pair.** Widen the enzyme set to *All commercial
enzymes*, check that the insert really carries the sites you expect, or try another
recipient. The table of all pairs gives the reason for each rejection.

**A result looks wrong.** Please say so, with the sequence and what you expected. That
is the single most useful thing you can send back.

---

## 10. Under the hood

Built on Biopython (ProtParam, Restriction, pairwise alignment), ANARCI for IMGT
numbering, FAMSA via `pyfamsa` for multiple alignment, `pymsaviz` for the alignment
figure, `dna_features_viewer` for the plasmid maps, `reportlab` for the PDF, and
Streamlit for the interface.

The synthetic-sequence generator uses a position-frequency matrix derived from
**PLAbDab-nano** (OPIG, University of Oxford); please cite that resource if you
publish anything built on the generated sequences.

The analysis logic lives in four importable modules — `vhh_core.py`, `sanger.py`,
`cloning.py` and `generator.py` — with no dependency on the interface, so a Jupyter notebook can use them
exactly the same way:

```python
import vhh_core as vc
df, annotations, errors = vc.analyze({"MyNb": "EVQLVESGGG..."})
```

The source is at <https://github.com/TyrosinDuck/vhh-analyzer>.
