# Changelog

What changed in each release, written for someone who uses the app rather than
for someone reading the commits. Versions follow `MAJOR.MINOR.PATCH`: a patch
fixes something, a minor release adds something without breaking what was
there,
and a major release breaks an existing use, which, since the analysis modules
are importable from a notebook, includes changing the signature of a function
in
`vhh_core.py` even when the app itself looks identical.

Each release is archived on Zenodo and gets a DOI of its own. The DOI in
`CITATION.cff` is the *concept* DOI and always resolves to the latest; to cite
the exact code that produced a result, use the version DOI of the release you
ran, from its Zenodo record.

## [Unreleased]

### Fixed

A second independent audit, run in November 2026 by a reviewer that had not
seen the code being written, over the roughly 19 000 lines added since the
first one. Three of its findings were about who to credit, which is the kind of
defect no amount of testing the arithmetic would have found. Nothing below
changes a number the app computes except where it says so.

What the audit confirmed is worth recording too, because it is the part with
the most at stake. The shipped polyreactivity weights are bit-identical to the
authors' released model, 8420 features and 791 non-zero with a maximum
coefficient difference of 0.0; the featurisation matches their published
method including a quirk of their loop that the L1 normalisation cancels; the
score parity, the quantile-grid cost and the sign convention all reproduce at
the figures already claimed. The duplicate rows in their released pool move a
percentile by at most 0.49 point, below the decimal the app shows.

- **The polyreactivity training-set size was a row count.** The screen and
  MANUAL.md §4.7 said the model was fitted on 141 474 nanobody sequences. That
  is the number of rows in the two pool files the authors released, 69 702 plus
  71 772. Their paper reports 65 147 unique low-polyreactivity and 69 155
  unique high-polyreactivity sequences, carrying 51 308 and 59 623 distinct CDR
  sets, and it is the CDR sets the model was fitted on. The row count is still
  the right denominator for the anchor measurements, which were run over every
  row of those files; it was never the right one for the paper.

- **A denominator that was neither the file's nor the paper's, in three
  places.** `specificity.py`, the annotator note on the screen and MANUAL.md
  all said the anchors were calibrated over "the 141 204 sequences of their two
  pools that carry a complete annotation". `tools/polyreactivity/build.py`
  iterates every row and applies no annotation filter, so the figure is
  141 474 and the qualifier named a filter that does not exist. The 66.6%
  acceptance and the 99.2% agreement are unchanged and were re-measured
  against the authors' real pools; what was wrong is what they are rates of.

- **The 90th-percentile flag was credited to the model's authors.** Three
  places called it "their own default". Harvey et al. give no score cutoff and
  no percentile anywhere in the paper, and their released code carries none
  either. The value stays, classified a convention as it already was, and the
  credit goes: an attribution nobody can check is the failure this repository
  calls worse than no reference at all.

- **The anchor CDR length limits contradicted the file they claimed to measure.**
  They were 4 to 19 for CDR1, 6 to 14 for CDR2 and 3 to 64 for CDR3, described
  as the observed range of the solved single domains. Against
  `data/sabdab/single_domain.csv` that was false: the snapshot holds a
  20-residue CDR1 and CDR2 lengths of 15 and 19, all of which those limits
  would refuse. They are now derived from the file at run time by
  `specificity.observed_cdr_bounds()` rather than written down, so they cannot
  drift from it again and refreshing the snapshot moves them. Two rows carry no
  CDR2 and are left out of that range, because a missing annotation is not a
  zero-length CDR; the one two-residue CDR1 is kept, which costs the check its
  floor, and the docstring says so rather than quietly hand-picking the data.
  With no snapshot on disk the check now passes everything instead of falling
  back on numbers typed into the file.

- **MANUAL.md capped the snapshot's longest CDR1 at seventeen residues.**
  Eleven heavy-chain single domains carry one longer than twelve, at 13, 14,
  15, 17, 19 and 20. The same paragraph gave 1 354 as the size of the whole
  snapshot and said it included the VNARs; 1 354 is the heavy-chain subset and
  the VNARs are the other 27, of 1 381. The anchor acceptance rates were
  re-measured and are now given as both: 441 accepted, 32.6% of the
  heavy-chain single domains and 31.9% of the file, with no VNAR accepted.

- **Two database dates were read as one, in two different modules.** The
  PLAbDab-nano export was downloaded on 2026-11-16 and carries nothing
  published after 2024-04-10, which is what `downloaded` and `latest_entry` in
  `data/vhh_position_meta.json` have always said; `build_matrix.py` called it
  "the April 2024 download" throughout, which reads as the day it was fetched.
  Separately, MANUAL.md turned the TNP thresholds' "July 2025" label into "the
  36 clinical-stage sequences available in July 2025". That label comes from
  the authors' web application, where it dates the most recent entry behind the
  thresholds, not a calibration; table 2 of the paper carries no date at all,
  as `tnp.py` already said.

- **The PS19-6L row credited the wrong people.** The article at
  10.1039/d5ra06185g is Retnadhas, Hegg and Ducat, RSC Adv 15:41115-41128
  (2025), not Dedisch et al. Volume, first page, year and DOI were right, which
  is how the error survived a DOI check; the sequence and the rest of the row
  are confirmed by that same paper.

- **"100 is the best and 0 the worst" promised a floor that never appears.**
  `rank(pct=True)` divides by the number of sequences, so the bottom cell is
  100/n: 12.5 for eight sequences, 50.0 for two. The screen now works the floor
  out from the set in front of it and says what it is.

- **PPC and PNC were described as red below 0.00.** They are patch scores
  bounded below by zero whose clinical range starts there, so that band cannot
  be reached. The guideline table derives the red column from a declared
  `floor` per metric rather than from the shape of the amber bands, which
  happened to agree for these two and said nothing about why. One consequence
  worth naming: `flag()` now returns green rather than red for an impossible
  negative patch score. Neither metric is computed by this app.

- **Two counters in two units printed side by side.** In
  `tools/sabdab/build_index.py` a row with no VH could only be counted in rows
  while a non-standard sequence was counted once however many rows it came
  from. They are `dropped_no_vh_rows` and `dropped_nonstandard_sequences` now,
  in the script and in `data/sabdab/meta.json`.

- **The JOSS manuscript cited neither method it ships.** `paper.md` now cites
  the Therapeutic Nanobody Profiler and SAbDab, with both entries added to
  `paper.bib` and every field confirmed against two independent records.

### Changed

- **The provenance register sees past module level.** Its completeness scan
  read module-level assignments only, so a scientific constant living inside a
  function was invisible to it, and the page claimed to cover every constant
  the app computes with. That is exactly where the CDR length limits above were
  sitting while they contradicted their own stated source. `scan_module` now
  also reports upper-case names assigned inside a function and numeric default
  arguments, which is where a threshold most often hides: `def mott_trim(quals,
  cutoff=0.05)` is a scientific constant with a name and a line number.

  That widening found 46 named numbers nobody had classified. Thirty are
  declared in `NOT_SCIENTIFIC`: identity defaults meaning "nothing applied",
  input caps, figure sizes, batch sizes and row counts. The other sixteen
  decide a result, and they are now on the register, which took a change of
  shape to do: a register entry names a value by its import path, so a default
  argument could not be registered at all. Each was promoted to a module
  constant with the signature reading it, which is what lets `resolve` show
  the value the code actually runs on, and which puts the justification above
  the function that applies it instead of nowhere.

  A new third list, `UNCLASSIFIED`, exists for the case where a number decides
  a result and nobody can say where the value came from. It is empty today and
  the mechanism stays, because the honest answer to that question is sometimes
  "we do not know" and the page should be able to say so rather than quietly
  inventing a source. The completeness test refuses anything on none of the
  three lists, with a message saying not to invent a source to close it.

  What the sixteen turned out to be, from the author in November 2026: ten are
  working conventions of the Center for Protein Engineering, three are
  conventions of the field that the laboratory teaches (the 5 angstrom contact
  cutoff, the 15-base annealing length and the five-residue Chou-Fasman
  window), one is fixed by a published model (the 3-mer length of the
  polyreactivity weights), and one is not this laboratory's at all:
  `sanger.mott_trim` trims at 0.05 because that is the value Biopython's own
  ABI reader hard-codes, citing the phred and phrap documentation, which was
  read out of `SeqIO.AbiIO._abi_trim` to confirm it. It is kept rather than
  chosen, so a read trimmed here and a read trimmed by the common tools are
  trimmed the same way. None of the sixteen is a measurement and none claims
  to be.

  One of them was a duplicate. The longest tolerated single-base run was the
  literal 6 in two signatures, `codons.suggest_gc_fix` scoring the risk on an
  existing gene and `vhh_core.reverse_translate` avoiding it in a designed
  one. One rule, so one constant now, `vhh_core.MAX_HOMOPOLYMER`, read by
  both.

  One was nearly made into a duplicate by mistake, which is worth recording
  because the mistake was tempting. `vhh_core.reverse_translate` balances a
  gene towards 40 to 60% GC and `codons.GC_LOW, GC_HIGH` flags a local stretch
  outside 30 to 70%. They look like the same band written twice and they are
  not: one asks where a whole gene should sit, the other asks whether a
  stretch is extreme. They stay separate, with the register saying why.

  What the scan still cannot see is a bare literal in the middle of an
  expression, which has no name to report and would drown the result in array
  indices. The page says so now instead of claiming otherwise.

- **The import test for `specificity.py` states what may not be imported**
  rather than what may. The allow-list version failed the moment
  `observed_cdr_bounds` started reading the SAbDab snapshot, a sibling module
  with no dependency of its own beyond pandas, while the property that matters,
  no heavy machine-learning runtime in the container, was never in question.

### Added

- **A Plates & dilutions section**, in the Bench group, with two views and no
  sequence asked for: like the phage titre, neither of them uses one.
  **Dilution series** works out what to pipette, by either of the two
  protocols that are actually used, and does not choose between them. Made by
  transfer, one pipette setting does the whole series however deep it goes,
  and the error of every step travels into the ones after it; made
  independently from the stock, no error travels, at the cost of more stock
  and of volumes that stop being pipettable. The carry-over is why the first
  is worth writing down: a position that must still hold its final volume
  after passing some on has to be made up with the final volume plus the
  transfer, so the transfer is V_final / (fold - 1) and the diluent per
  position is the final volume itself. Forget it and nothing complains, the
  last positions simply come out short, or the last dilution is not the one
  that was meant. When a volume is not pipettable the screen says which
  volume would fix it rather than only refusing: at 10-fold steps a 2 uL
  minimum needs 18 uL per position. The last position keeps the transfer it
  never passes on, which is stated rather than folded into the figures.
  **Plate layout** says which well holds what, on 96 or 384. The placement is
  predictable rather than packed: each replicate of each sample starts at the
  top of a fresh column, so the series runs down the column and the replicates
  sit beside it, which is how the plate is pipetted. Dilution-major order
  would interleave the two and put neither on an axis. Nothing is optimised,
  because a layout whose positions move when an unrelated sample is added is
  one nobody can check against the bench. A series longer than a column is
  refused rather than wrapped across two, where it would read as two
  different samples; a plate that overflows is counted, with the number of
  plates the run would take. Leaving the outer ring empty is offered as an
  option and nothing more: the app neither quantifies evaporation at the edges
  nor takes a view on whether it matters for a given assay. The map carries no
  volumes, since what a well holds depends on the plate and the protocol and
  not on the count of its wells. It downloads as a grid and as one row per
  well, which is the form a reader or a database ingests.
  The arithmetic is in `plates.py`, importable from a notebook, and it takes
  the pipettable minimum from `bench` rather than keeping a second one. The
  home page diagram grew a box for it, which is why the canvas is 880 tall
  and the DNA arrow runs back at y=800: the box went between the hub and the
  buffers first, where it sat on top of a dashed label. Checked by rendering
  the diagram and looking at it.
- **Couple two proteins**, as a view of Batch & QC beside Conjugate. That one
  adds a small reagent from a weighed vial at a known molarity; this one adds
  a second protein from a solution in mg/mL, so the arithmetic runs through
  two molar masses: a SpyTag/SpyCatcher pair, a sortase ligation, an
  enzymatic conjugation. The molar ratio decides which partner is limiting,
  the excess of the other is named in nanomoles with the reminder that a
  one-to-one couple cannot consume it, and a volume of B too small to pipette
  becomes a dilution of B, as in Conjugate.
  Whether the bond is a condensation is **asked, not assumed**: a condensation
  releases one water, so the conjugate weighs the sum minus one water, and an
  isopeptide bond is one. The app does not tabulate which pair forms which
  bond, because that is a fact about the user's reaction, and a delta from a
  table nobody checks is how a wrong mass becomes a wrong identity. The water
  is derived from its formula with the atomic weights in `buffers` rather than
  typed. No yield is predicted: the conjugate concentration reported is what
  would hold if every molecule of the limiting partner reacted, and is
  labelled a ceiling.
  Neither molar mass is prefilled. A comes from the sequence picked at the top
  when there is one; B is always typed, and the view refuses until it has
  both. The first build prefilled B with 28000 Da, which the regression test
  guarding exactly that caught: an invented mass there would have set the
  volume of B to add, which is the one number the view exists to give.

- **A comparison view in the Triage family, which still refuses to score.**
  The triage table is the right shape for sorting and the wrong shape for
  "which of these eight differs from the others, and where". A second tab
  answers that in three layers, none of them added to another: a heatmap of
  each sequence's rank **within the submitted set** on each axis separately,
  with the raw value written in every cell; the flags of the same table read as
  gates, go or no-go, marking a sequence without touching its ranks; and the
  Pareto front over the axes chosen, which is the only answer to "which is
  best" that needs no weights.

  There is still **no composite score**, and that is not an omission. Adding a
  percentile from a model fitted on 69 702 sequences to a count of motifs
  weighted in this laboratory produces a number with no unit, and a weighted
  total would let three comfortable axes outvote one unpaired cysteine. A front
  holding every sequence is a real answer: the axes disagree and the choice is
  yours.

  Only columns with a defensible direction are offered as axes. There is no
  good pI, no good molecular weight and no good CDR3 length, and the hallmark
  score is an identity rather than a quality, so those stay in the table. Each
  axis carries the direction that is better, its basis in the app's own three
  words for a result, and the reason that direction is defensible; GRAVY is
  labelled the weakest of them, because a VHH needs hydrophobic residues in its
  core. Nothing here claims a sequence is fit for use in people: humanness is
  identity to the closest human germline and not a prediction of
  immunogenicity.

  The arithmetic lives in a new importable module, `ranking.py`, which reads
  the triage table and nothing else, so a notebook can have the same ranks,
  gates and front. The axes name triage columns and the gates come from
  `triage.FLAG_COLUMNS` rather than from a second list; a test fails if an axis
  names a column the table no longer produces. `MANUAL.md` §3.2c is the long
  form of the three layers and of what the view will not do.

- **The data behind every figure, beside the figure.** Each figure and table in
  the app now offers a CSV of what it was drawn from, either next to it or, in
  a tab drawing several views of one computation, as a single button holding
  exactly those rows. Roughly thirty new downloads: the titration curves, the
  per-residue charge and hydropathy profiles, the Chou-Fasman propensities and
  state track, the secondary-structure fractions, the liability score by
  sequence and type, the per-window aggregation profile, the dendrogram as a
  table of merges, the CDR-length distribution drawn against the reference, the
  GC profile of a gene and the before-and-after pair of the suggested swaps,
  Phred per base for each Sanger read, the pLDDT or B-factor curve, the PAE
  matrix, and the ordinary tables that had none: properties, hallmarks,
  germlines, repairs, patches, vector features, primers, codon usage.

  What makes the download worth having is that **the figure and the CSV come
  from one call**. Five plot functions used to compute their own numbers inline
  and now ask for a frame: `plot_charge` takes `titration_curves`,
  `plot_profiles` takes `profile_table`, `plot_liability_profile` takes
  `liability_profile_matrix`, and `chou_fasman_table` and `linkage_table`
  reshape what their figures already use. A download that recomputed the same
  thing would agree with its figure until one of the two was edited, with
  nothing to say which one the user was holding; `tests/test_figure_data.py`
  patches each data function and checks the figure follows it, which is the
  test that would catch a second computation creeping back in.

  Three figures have no table behind them and say so instead of offering a
  substitute: the alignment and the logo are drawn from the alignment itself,
  offered as FASTA, and the 3D viewer is the structure file. §2.1 of the manual
  states the rule and the caveats that travel with these files, among them that
  a profile value is a window mean and not a property of its residue.

### Fixed

- **The front figure no longer runs its axes backwards.** It inverted any axis
  whose good end is low, so that the front always sat in the top right, which
  is a real convention in multi-objective plots and the wrong trade here: the
  axis labels have to be read anyway when the directions are mixed, so the
  inversion moved the confusion from the corners to the tick marks, where it is
  met once per point instead of once per figure. Both axes now run in their
  natural direction, the better end of each carries a green arrow, and the line
  above the figure names the corner the arrows point to, derived from the axes
  in hand. Reported from the first build, where the two default axes are both
  "lower is better" and the figure read as though the largest values belonged
  on the left.

- **The README's module paragraph had lost four modules**, `specificity.py`,
  `tnp.py`, `sabdab.py` and `elsewhere.py`, each of which had been added to the
  app without being added to the sentence that says what every file at the root
  is for. The same failure as the line count that paragraph already describes,
  and noticed only because `ranking.py` was being added to it.

  The list stays written out by hand, because prose that groups the modules by
  what they are for says more than a generated enumeration. What is new is that
  `tests/test_metadata.py` now fails on a module at the root that neither the
  paragraph nor the generated section table above it names, so the next one
  cannot be added silently.

- **Polyreactivity, as a published predictor, in the Developability family.**
  A new tab scores the CDRs with the logistic regression of Harvey et al., Nat
  Commun 13:7554 (2022), and reports the score, its percentile against their
  low-polyreactivity pool of 69 702 sequences, and a flag at their default
  90th percentile. Only their fitted coefficients are shipped, as the 791
  non-zero weights of 8420 in an `.npz`, so scoring needs numpy alone: no
  scikit-learn, no pickle, nothing that breaks on a version bump, and 90 kB of
  data files in total. The percentile also appears in the triage table under a
  new `Specificity` column group, with a `polyreactive` flag.

  Two things are negated or widened relative to their publication, and both are
  stated in the interface because they change the number. Their decision
  function is positive for *low* polyreactivity, and it is flipped here so
  higher means more polyreactive, like every other flag in this app. And their
  CDR2 is IMGT 55 to 65, one residue wider at the N-terminus than the CDR2 this
  app displays; the model is fed their definition, rebuilt from the same
  annotation.

  `tools/polyreactivity/build.py` is the derivation of all three data files,
  from their MIT-licensed repository, and it is also the check: on every run it
  verifies that this implementation reproduces the scores they published for
  all 69 702 sequences of their clean pool, which it does to a maximum
  difference of 1.8e-15. Their copyright notice is in `LICENSE`. `MANUAL.md`
  §3.8c has the caveats, which are the point: one synthetic library, one
  framework, CDRs only, and a pre-screen ranking rather than a measurement.

- **A plastic resemblance score, labelled uncalibrated everywhere it appears.**
  A second tab reports how much the CDRs look like published
  polystyrene-binding peptides, in two modes, as a percentile against the same
  naive library. It is **not** a predictor and the interface says so before
  showing a number: no dataset labels VHHs for adsorption on polystyrene, the
  weights were set by hand in this project to reproduce the composition of each
  published class, and every row carries `calibrated = False`. The two modes
  are never added to each other or to the polyreactivity score, and the plastic
  columns are deliberately absent from the triage table, which is read as a
  table of flags. `provenance.py` gained a basis word, `unfitted`, for exactly
  this kind of number, and both weight sets are registered under it.

  The reference set is twelve peptides, each row carrying its own citation.
  Nine peptides from the draft this module came from were dropped: their only
  attribution was a review table whose text does not contain them, and the data
  file names them and says why. §3.8d covers the rest.

- **CDR annotation for the model, separate from the app's own, because the two
  are not interchangeable.** `specificity.py` carries its own IMGT anchors for
  servers without HMMER. This looks like a duplicate of `vhh_core`'s regex
  annotator and is not: measured against the model authors' own annotation of
  14 192 of their sequences, the app's display boundaries reproduce none of
  their CDR sets, move the percentile by a median of 9 points and would place
  the flag differently on 8.45% of sequences. The new anchors accept 66.6% of
  their 141 204 annotated sequences and reproduce the IMGT CDRs exactly for
  99.2% of what they accept; where ANARCI runs, the two paths agree on every
  sequence both accept. **A sequence whose CDRs cannot be placed keeps its row
  and gets no score**, with the reason in the row, including one ambiguity no
  regex can settle: where the residue opening the FR3 anchor is a threonine,
  the authors' numbering assigns it to CDR2 about half the time. The length
  limits that catch a misplaced anchor are the range of the solved single
  domains in `data/sabdab/`, and they are **not** applied to an ANARCI
  annotation: a numbering is not a guess, ANARCI gives ten of those domains a
  CDR1 longer than twelve residues, and the first version of this refused to
  score them with a note blaming anchors that had not run. On that diverse set
  the anchors accept only a third of the sequences, against two thirds of the
  model authors' single-framework library, and §3.8c now says both figures:
  without HMMER, expect a third of a varied set to be scored.

- 23 tests in `tests/test_specificity.py`, against a fixture of 150 sequences
  taken from the model authors' released pools, carrying their own CDR
  annotation and their own published score, so the ground truth comes from
  outside this repository. Two of them exist to stop a composite index growing
  and to stop the runtime acquiring scikit-learn.

- **The reference repertoire matrix has a build script, and therefore a
  provenance.** `data/vhh_position_freq.csv` is what the Repertoire screen calls
  a residue unusual against and what the synthetic generator samples frameworks
  from, and until now it was simply present in the repository: PLAbDab-nano was
  named and cited, but nothing recorded which download it came from, when, or
  what had been filtered out to produce it. It could not be rebuilt and it could
  not be checked. `tools/plabdab/build_matrix.py` is the missing derivation, and
  it reproduces the committed matrix byte for byte from the April 2024
  PLAbDab-nano VHH download, which is the only evidence that the filters it
  applies are the ones originally applied: they were recovered by reproducing
  the output rather than read off a note. Run it with `--check` to compare a
  rebuild against the committed files. `data/vhh_position_meta.json` now carries
  a `provenance` block naming the script, the download, the date of its newest
  entry, the length cut-off and the count of every sequence dropped at each
  step, so both screens can say where their reference comes from. One thing is
  not known and is marked as such in the script: any cut-off from 99 to 101
  residues separates the seven fragments in that download from the 1958
  sequences kept, so the original value is somewhere in that range and
  `MIN_LENGTH` is set to 99, the low end. `MANUAL.md` §3.11c covers all of it,
  including why rebuilding this matrix is less pressing than refreshing the
  structure index, and what has to be re-examined when it is rebuilt.

### Changed

- **Both screens read the size of the reference set out of the file instead of
  stating it.** The number 1958 was written into the Repertoire family's
  description, the generator's caption and three explanatory captions, and the
  support slider's ceiling of 1900 was one matrix's maximum support copied by
  hand. Rebuilding the matrix from a newer download would have left every one of
  them describing a reference set the app was no longer using, silently, which
  is the duplicated-list failure this repository has had twice before. They now
  come from `repertoire.reference_report`, over `generator.matrix_report`, which
  gained `max_support` and the metadata's `provenance` for the purpose. The
  Repertoire screen also now states what the reference actually is: nanobodies
  from patents and the literature, biased towards what somebody thought worth
  publishing, and not a naive immune repertoire.

- **The recipient vector is drawn where you choose it.** The map and the feature
  inventory sat at the foot of the in-silico cloning screen, behind the stop
  that waits for donor DNA, so the vector you had just picked could not be
  looked at until you had pasted a donor. Both now sit in step 1, under the
  vector's length, GC, feature count and cloning window. There is one map on the
  screen and not two: it is drawn into a slot when the vector loads, showing the
  file's features in blue, motifs detected here in orange and a VHH scaffold in
  green, and the same slot is written again once the plan exists, adding the
  usable enzyme sites in red. The caption says which of the two you are looking
  at. The view selector, whole plasmid or cloning region, sits outside the slot,
  because a widget created twice in one run is a duplicate-key error rather than
  a second selector.

### Removed

- **The "Expected tags" selector, which asked a question the vector answers.**
  It offered the two vectors `sanger.VECTORS` had hand-written entries for while
  the app ships five, and it was a free choice: picking pMECS for a pMES3 clone
  made the construct check report a missing HA tag. The expected cassette is now
  read from the vector in hand and stated under the picker, as
  *PelB → VHH → HA → spacer(GS) → His6 → stop* for pMECS. A vector whose
  annotation yields no tag says so, and the construct check reports the
  translated cassette without looking for a tag in it, instead of checking
  another vector's tags against it.

### Fixed

- **Identical CDR3s are visible in the clonal-family dendrogram.** Sequences
  with the same CDR3 join at a distance of zero, which is exactly where the
  identity axis ends, so their link was drawn under the axis spine and the
  clearest result in the plot, two reads of the same lineage, was the one merge
  that could not be seen. Those links are now drawn as a thick bar at 100%,
  each sequence that has an identical twin carries a dot there, and the axis is
  extended a little past 100% to leave room. Nothing is moved off 100% to make
  it visible: the extra room is in the axis, not in the data. Found on a plate
  where several picked colonies carried the same clone, which is the case the
  tab exists for.

- **The HA tag is recognised whole and truncated, and only where it counts.**
  The motif library held the HA epitope as the nine residues YPYDVPDYA, which is
  what the epitope is, and found it in none of the five bundled vectors: pMECS
  and pMES4 both encode eight of them, YPYDVPDY followed by the linker GS, which
  is how the vectors' authors describe the tag themselves, the C-terminal
  extension of a pMECS clone being AAAYPYDVPDYGSHHHHHH (Vincke et al.,
  US 2015/0093336 A1, SEQ ID NO 261 and 262). pMECS got away with it because its
  map annotates `HA-tag` by hand; pMES4 carries the same sequence unannotated and
  the app did not see it at all, and a bare FASTA of either vector lost the tag
  entirely. Both forms are now searched for, longest first, so a position
  claimed by the full epitope is not reported a second time as its own prefix,
  and the feature's note says which form is there rather than rounding eight
  residues up to nine.
  This alone would have broken pMES4. Its cassette is His6, the amber stop, then
  that HA tag, then gIII: in a suppressor strain the amber is read through and
  the tag is on the phage fusion, and in a non-suppressor strain the protein
  stops at the His6 and has no HA at all. Since the cassette extractor looks for
  tags in the stretch before the first stop, expecting that tag reported
  *missing tag(s): HA* on every correct pMES4 clone, which is the pMECS failure
  over again from the other end. Measured on a real pMES4 clone before the
  second half of the fix was in.
  So the derived cassette now ends at the first stop, in the tags' own reading
  frame rather than the leader's, and a tag past it is named under the picker as
  something the vector carries but the soluble protein does not. On pMES4 that
  reads *Also carried, past the stop that ends the cassette: HA at 361–384*. It
  stays on the map and in the feature inventory, where it is real.

- **A vector with no entry of its own was checked against pMES3's tags.** Three
  of the five bundled vectors, pMES4, pHEN6 and pHEN25, had no entry in
  `sanger.VECTORS`, and neither did any uploaded vector. `verify_construct`
  resolved an unknown name to pMES3 and carried on, so a pHEN6 construct was
  checked for pMES3's tags and an upload for a cassette nobody had claimed it
  carried. An unknown name now yields no tags rather than another vector's.
  Where the two vectors happen to carry the same His6 this changed no verdict,
  which is why it went unnoticed; on an upload with a different cassette it
  reported a tag as missing that was never expected.

### Added

- **An example run on the sequencing screen, and an example structure.** Both
  screens take a file and can do nothing without one, which left anybody
  without an `.ab1` or a model looking at an empty page. *Load an example run*
  fills the sequencing screen with six reads, one of each status the table
  reports: a clean clone read from both ends, a second clone, a clone with one
  nucleotide deleted inside the VHH, a read that stops before FR4, and an empty
  well. *Load an example structure* brings in PDB entry 1MEL, cAbLys3 bound to
  hen egg-white lysozyme (Desmyter et al. 1996, doi:10.1038/nsb0996-803), which
  ships in `data/examples/` because the app does not download structures on
  your behalf. PDB archive files are CC0 1.0, and the entry is cited on the
  screen and in `data/examples/README.md`.
  **The reads are simulated and say so**, above the table, in the module that
  makes them and inside the files themselves: their sample name and comment
  fields both carry the warning, so it survives a download. What is real in
  them is the VHH sequences, which come from the bundled library, and the
  vector, which is the bundled pMES3 map with each read's VHH put in place of
  its stuffer. What is invented is every trace and every quality value: a curve
  that rises, holds and falls like a Sanger read, with no measurement behind
  it. No number read off that example describes a real run.
  They are built when the button is pressed rather than committed, so they
  cannot drift away from the library or the map they are derived from, and
  nothing binary is added to the repository for them. `demo.py` writes real
  ABIF files so that the example goes through the same parser an upload does;
  that writer has been checked by reading its own output back through
  Biopython, which is the parser the app uses, and not against an instrument's
  own software.

- **The amino-acid sequence of each chain, on the Structure screen.** Under the
  chain table, two sequences per chain where the file allows it: what the
  coordinates hold, and what the file declares the chain to be, from SEQRES in
  a PDB file or `_entity_poly` in an mmCIF. The second is routinely longer, by
  the termini, tags and loops that were too disordered to model, and the
  residues that are declared and not modelled are shown in lower case and
  counted at each end and inside the chain. On the bundled 1MEL, chain A is 132
  residues of the 148 the entry declares, the missing ones being the first
  residue and the C-terminal `GRYPYDVPDYGSGRA`.
  A modelled residue that disagrees with the declared one is reported rather
  than quietly aligned over, a file that declares nothing says so instead of
  showing the same sequence twice, and breaks in the residue numbering are
  deliberately not marked inside the sequence, because an IMGT-numbered or
  Kabat-numbered file skips numbers by design and a jump in the numbers is not
  a missing residue. Both sequences download as FASTA and any chain can be sent
  to the analysis section; sending the declared sequence sends its tag with it,
  which the screen says, because the analysis would then read that tag as part
  of the protein.
  The mmCIF route is written against the mmCIF dictionary and tested on a file
  written here. No mmCIF from the PDB or from AlphaFold has been run through it
  yet, which the function's docstring says as well.

- **The cassette layout is read out of each vector's own map.** `sanger.VECTORS`
  was two entries typed in beside plasmid maps the app already ships, and the
  release above fixed the pMECS HA tag by correcting one of them by hand, which
  left the second list in place and free to drift again. `cloning.derive_layout`
  now reads the layout from the map: it takes each annotation downstream of the
  cloning window that is short enough to be a tag and whose label names one, and
  translates the span **in its own frame**, because an empty phage-display
  vector carries a stuffer and its tags are out of frame with the leader until a
  VHH of the right length is cloned in.
  pMECS derives to the eight residues YPYDVPDY, which is what its map encodes
  and what the release above had to be told by hand, and the three vectors that
  had no entry get one from their own file. The sequencing screen's vector list
  is the five bundled maps rather than the two with entries; it points at the
  default by name, since the list is sorted and index 0 is no longer pMES3.
  A label has to match on word boundaries: pMECS annotates gIII as
  *Protein 3 M13 phage*, whose lower-cased form contains `ha` inside `phage`,
  and the substring test that finds the cloning window takes it, pulling 409
  residues of phage coat protein in as an HA tag. A length bound stops the same
  thing independently. `tests/test_derive_layout.py` pins what the five maps
  derive to, pins the pMECS HA case, and lifts each of those two guard rails in
  turn to show the other still holds, then lifts both to show what they are for.
  Verified by putting the hand-written table back on a copy, which fails seven
  tests across the two layout files with the vectors named.
  What this cannot check is whether a map's own annotation is right. If a
  vendor's map labels the wrong span, the app reads the wrong span: deriving
  moves the trust from a list typed in the repository to the file you supply, it
  does not remove it.

- **A correct pMECS clone is no longer reported as missing its HA tag.** The tag
  sequences in `sanger.VECTORS` are written out by hand beside the plasmid maps
  the app ships, and the two had drifted apart. The pMECS entry carried the HA
  epitope on record, the nine residues YPYDVPDYA. The pMECS map annotates its
  `HA-tag` over 24 nucleotides, which encode the eight residues YPYDVPDY, and
  the two codons after them are GS, not A: the nine-residue form is in none of
  the six frames of that map. So every correct pMECS clone came back with
  *missing tag(s): HA*, on the sequencing screen reading an .ab1 and in the
  cloning screen's construct check alike, which is exactly the failure the
  comment under that list warned about.
  What is searched for is now the eight residues the vector really encodes.
  They are a prefix of the nine, so a cassette carrying the canonical epitope
  still matches: the change makes the check accept pMECS without making it
  reject anything it accepted before. What it can no longer tell apart is a
  cassette that stops after the eighth residue from one that carries a ninth,
  a distinction these vectors do not have.
  Whether an anti-HA antibody detects YPYDVPDY followed by GS as well as it
  detects YPYDVPDYA is a question about the antibody and the app does not answer
  it; it reports what the vector encodes.
  `tests/test_vector_layouts.py` now checks every entry of `VECTORS` against the
  map of the same name in `data/plasmids/`: each tag sequence has to be found in
  some frame of its own vector, the HA case is pinned on its own, and a cloned
  pMECS cassette is read back end to end. Verified by putting the nine-residue
  value back on a copy, which fails three of those tests with the sequence
  named. A layout that disagrees with its vector now fails the suite instead of
  being found at the bench.
- **A polyhistidine tag longer than six no longer appears several times over.**
  The His6 motif matches at every offset inside a longer run, so a vector with
  ten histidines produced five overlapping features, at 70-87, 73-90, 76-93,
  79-96 and 82-99, all five labelled His6: five tags drawn on the vector map and
  listed in the feature inventory where the vector carries one. The run is now
  reported once and named for the length it has, so ten histidines read His10.
  The five bundled vectors all carry six and still read His6, and the cloning
  window, which finds the first tag by matching `his` in its label, is unaffected
  by the new name.

### Added

- **What a 1.x release promises not to break is written down and tested.**
  The README's new *Public API and versioning* section declares public the
  names in `vhh_core.__all__` and the three files the app writes for a user to
  keep and reads back: the batch record, the session file and a sequence
  library CSV. Breaking any of them is a MAJOR release; every other module,
  importable as it is, may change in a MINOR one. Until now the rule named
  only `vhh_core` signatures, while the private library and the session file
  had put two more formats into users' hands. `tests/test_public_api.py`
  checks `vhh_core` against a snapshot of its 57 public names and their
  parameters, written by `tools/public_api.py`, and reads back a frozen batch
  record, session file and library from `tests/compat/`. Verified by breaking
  each promise on purpose on a copy: a default removed, parameters reordered,
  one renamed, a required one added, a name dropped from `__all__`, the batch
  record's format number moved without reading the old one, and a library
  column renamed each fail the suite; adding a parameter with a default does
  not. The manual says the same to users in §13, including what is not
  promised: values recomputed from a kept file can change when a calculation
  is corrected.
- **My library: a private sequence library, held in the browser session.** The
  Sequence library screen now has two tabs. **Public library** is the file in
  the repository, unchanged. **My library** loads a CSV in the same format,
  lets you add, edit and delete entries in the app, and gives the whole
  library back as a download. Entries go on to the analysis and to gene
  diagnosis exactly as from the public one, and entries of the public library
  can be copied into yours; a name already in your library is skipped, not
  overwritten. Nothing is written to disk or put in the shared cache, so the
  library ends with the session: the screen says when there are changes you
  have not downloaded, and asks before loading another file or closing over
  them. The checks are the ones the public file has always had, now in
  `library.validate_library`, shared by the file, the upload and the editor;
  `library.read_library_bytes`, `library_csv`, `library_digest`,
  `editable_frame` and `add_entries` are new. `load_library` keeps its
  signature and its behaviour, with one difference: a row with nothing in any
  column is now skipped instead of reported. `.gitignore` excludes
  `*.private.csv` and `private/`, as a safety net for a private library saved
  inside a clone. Tested: the functions, and the screen through `AppTest`,
  which cannot load a file or type into the editor; those two were not
  driven by a test.

### Fixed

- **A `#` inside a library row no longer cuts the row.** The library was read
  with pandas' `comment="#"`, which ends a line at any `#`, not only at the
  start. A name such as `NbNAMPT_cluster#1-(18)` became `NbNAMPT_cluster`,
  the sequence after it was lost, the row was reported as having no usable
  sequence, and two clusters collided on the same name. Only lines that
  start with `#` are comments now, as the header of `data/library.csv`
  always said. The same reader no longer turns a cell reading `NA`, `None`,
  `null` or `N/A` into an empty one, which would have erased a sequence
  named that way. Found on a real private library; the public file had no
  `#` in its rows and reads the same as before.

### Changed

- **Library entries are picked by ticking a Send column in the table**,
  in both the public and the private library, instead of from a separate
  list of names under it. The send buttons still say how many of the ticked
  entries carry a protein and how many carry DNA. In the public library the
  ticks are kept by entry name, so a filter no longer has any effect on
  them: an entry it hides stays ticked, is sent with the others, and the
  line under the table says how many are hidden; **Untick all** clears them.
  In the private library the tick is a column of the editor, so it belongs
  to its row and survives edits; it is not written to the downloaded file
  and does not count as an unsaved change. Checked in Chrome as well as by
  the tests: an earlier version redrew the public table from the ticks
  after every click, and the browser then ignored every third click on the
  same box. The table is now redrawn only when the filter changes.

- **The provenance register now covers every constant, and says so because it
  is checked.** The register listed 26 constants and the manual said it listed
  every constant the app computes with. Those were not the same set. An audit
  of the modules found a dozen numbers that change what a user reads and were
  on no page: the lower ipTM band (while the upper one was registered), the
  per-residue charges behind the charge profile and the charge term of the
  aggregation score, the weights combining that score's four terms, the
  liability base scores and their Low/Medium/High bands, the protease
  per-site scores, the coverage a V domain must span before a sequence is
  called a fusion, the rarity threshold and minimum reference support in the
  repertoire comparison, the CDR3 identity at which two sequences become one
  family, the minimum size of a surface patch, the GC synthesis band and the
  rare-codon cluster window, the Sanger trimming and quality thresholds, the
  pLDDT colouring bands, and the mass-difference modifications. The register
  now holds 50 entries.

  Most of the new entries are **conventions**, which is the point: a number
  somebody chose, defensible and arbitrary, is exactly what a reader cannot
  tell apart from a measurement by looking at it. Two are worth naming. The
  per-residue charge table hands histidine a fractional charge as a stand-in
  for partial protonation near pH 7.4, which has no source and is now
  labelled as having none; nothing quantitative rests on it, since the net
  charge and pI the app reports come from ProtParam. And the liability scores
  come from an in-house notebook that weights each motif by solvent exposure
  and pLDDT, neither of which this sequence-only module applies, so the bands
  compress upward and the ranking between sequences carries more than the band
  any one sequence lands in. That was already written in the code and is now
  on the page.

  **Completeness is a test, not a promise.** `provenance.NOT_SCIENTIFIC`
  names the constants that carry no scientific claim, each with its reason:
  guard rails, cache settings, figure sizes, widget defaults, and unit
  conversions exact by definition. The suite walks the modules `app.py`
  requires, finds every public constant holding a numeric literal (tuple
  assignment included, which is how both ends of the GC band had stayed
  invisible), and fails on one that is neither registered nor declared. A new
  threshold now has to be classified before the suite passes. What the test
  cannot check is whether a classification is honest; it can only make the
  decision explicit. MANUAL.md §12b.

- **The library generator reads the FR2 hallmark tetrad from `vhh_core`
  instead of keeping its own copy.** `generator.HALLMARK_POSITIONS` and
  `generator.VH_LIKE` were literals listing the same four positions as
  `vhh_core.HALLMARK_POSITIONS`. They agreed, which is why nothing reported
  them, and that is the failure mode this project has already hit twice: a
  generator enforcing one reading of the tetrad while the analysis screen
  applies another would produce sequences the app then calls VH-like, with
  nothing pointing at the cause. They are now derived, the derived values
  were checked against the literals they replaced, and a test in
  `test_regressions.py` compares them to those literals rather than to
  `vhh_core`, which would pass whatever the derivation did.

### Added

- **Closest structure: has anything like this already been solved?** A new tab
  in the Repertoire family compares each sequence against a snapshot of
  SAbDab2's single-domain entries and reports the nearest solved structures,
  with the PDB codes, resolution, method, species and what each one was bound
  to. A close hit is a modelling template and, where the entry is a complex,
  an epitope somebody has already looked at. cAbLys3 finds its own crystal
  structure, 1MEL, at 100% identity, which is the end-to-end check the test
  suite now runs.

  Framework and CDR3 identity are reported **separately and never combined**:
  two VHHs at 90% over the framework can have CDR3 loops with nothing in
  common, which is exactly where one number misleads. Shark VNAR entries are
  kept in the index, labelled, and left out of the search by default.

  **The snapshot carries its own date and the app says it out loud.** SAbDab2
  is updated weekly from the PDB and the shipped file is not, so the
  extraction date is read from the data, turned into a number of days, and
  printed as a sentence every time results are drawn. Past `STALE_DAYS` that
  sentence becomes a warning naming how many weekly updates have been missed
  and pointing at the live database, so an empty result can never read as
  "nothing like this exists". A missing or unparseable date counts as old,
  not as fresh. MANUAL.md §3.11b.

  **The data is CC-BY 4.0, not MIT.** `data/sabdab/LICENSE`, the root
  `LICENSE` and the README now say so, and the screen carries the attribution
  and the citation.

  `tools/sabdab/build_index.py` derives the shipped index from SAbDab2's own
  export: one row per distinct VH sequence, every structure containing it,
  and only the columns the app reads, 4759 export rows down to 1381
  sequences. It refuses an export whose columns it does not recognise rather
  than writing an index with a silently empty field, and it drops rows with
  no VH and sequences carrying `X`, counting both in the metadata. Nothing in
  that directory is edited by hand.

  One thing was checked rather than assumed before any of this was built:
  SAbDab2's CDR columns are the same IMGT windows this app annotates with.
  Twenty-five entries were compared residue for residue, all three CDRs
  agreed on every one, and a test keeps checking a sample. The CDR3 identity
  column would otherwise be comparing two different substrings.

- **Clinical context: where a sequence sits against nanobodies that reached
  the clinic.** A new tab in the Developability family applies the
  Therapeutic Nanobody Profiler guidelines (Gordon et al., Commun Biol 9:344,
  2026) to the total CDR length and the CDR3 length, flagging each green,
  amber or red against the range spanned by the 36 clinical-stage nanobody
  sequences. This answers a question nothing else here asked: not whether a
  sequence carries a known problem, but whether anything like it has ever
  been taken into a clinic. A red flag is not a defect, cAbLys3 is red on
  both. The thresholds are shown with the date they were calibrated, because
  the authors say plainly that they will move as more candidates progress
  through trials. See MANUAL.md §3.8b.

  Two guards are worth knowing about. TNP counts in IMGT numbering, so a
  sequence annotated with the regex anchors has its lengths reported and
  **no flags at all**, with a warning: the two definitions differ by several
  residues on a real VHH, and a flag against the wrong boundaries looks
  exactly like a correct one. And the three patch scores in TNP (surface
  hydrophobicity, positive and negative charge) are **not computed here**,
  deliberately: reimplementing the TAP patch methodology from its description
  without the reference code would give numbers on an uncalibrated scale, and
  flagging those against published thresholds would be worse than reporting
  nothing. The screen points at the authors' web application for them.

- **CDR3 compactness, measured on an uploaded model.** The one TNP metric
  this app can compute, on a new view in the Structure screen: the CDR3
  loop's length divided by how far it reaches, which separates a loop that
  extends away from the domain from one of the same length folded over the
  framework-2 face. Reported with the reach, the residues measured, and a
  flag against the clinical range. It is a reimplementation from the paper's
  prose and has not been checked against the authors' values, which the
  screen states every time it draws a number: two choices were forced that
  the paper does not specify, so the absolute value may sit on a slightly
  different scale from theirs. Rank your own panel with it. MANUAL.md §4.5b.

- **Liability scores weighted by exposure and model confidence.** A tick box
  on Liabilities in 3D multiplies each sequence liability score by the
  residue's measured relative accessibility and by pLDDT/100. This is the
  correction §3.7 has always said the sequence-only section cannot make: an
  `NG` buried in the core and one on the CDR3 apex score the same 10 points
  from sequence alone, which is why almost everything lands in Medium or
  High. The weights are this app's choice and not a published formula, the
  risk category is still not recomputed, and an unmodelled residue gets an
  empty cell rather than a zero, because "not in the model" and "buried" are
  different claims. Read as a B-factor rather than pLDDT, no confidence
  factor is applied at all. MANUAL.md §4.6.

- **Elsewhere: what this app does not compute, and who does.** A third table
  on the Reference screen, and the only transcribed one there: nine published
  tools and databases, each with what it determines, what this app does
  instead, and why. Grouped by whether the method is not implemented
  deliberately, approximated here with stated weights, or out of reach in a
  690 MB container. Every citation was checked against the publication or its
  record and carries a DOI. Importable as `elsewhere.py`. MANUAL.md §12.3.

### Fixed

- **A loop variable shadowed the analysis run object.** In the Developability
  family, the repair-suggestion block iterated with `for _, r in
  doable.iterrows()`, and `r` is also the name of the `Run` the family is
  handed. Python leaks a plain `for` target, so from that point on `r` was a
  pandas row rather than the run, and anything further down the family that
  reached for `r.sequences` raised `AttributeError`. Nothing did until the
  new tab was added at the end of that family, which is the only reason this
  surfaced now. The loop variable is renamed and a comment says why.

These all come from one independent audit of the repository, run in November
2026 by a reviewer that had not seen the code being written. Every item is a
claim the app or the manual made that the code beside it did not support.
Nothing here changes a result the app computes, with one exception noted
below; what changes is what the app says about its own numbers.

- **The Reference screen promised a limiting tRNA for all eight rare codons,
  and two of them have none.** `CGA` and `CGG` carry no tRNA name, so two
  rows of the codon table showed an empty column under a caption asserting
  every row had one. They are flagged because the Rosetta and CodonPlus
  strains supplement them, which is not the same as one tRNA being the
  documented bottleneck, and the caption and MANUAL.md §12.1 now say six. A
  test counts the named tRNAs against the sentence so the two cannot drift
  apart again in either direction.

- **The CAI was presented as Sharp & Li while using a whole-organism table.**
  Sharp & Li (1987) define the index against a reference set of *highly
  expressed* genes. The shipped table, `data/ecoli_codon_usage.csv`, is
  Kazusa's whole-organism one. The arithmetic is theirs and the reference set
  is not, which makes a number from this app not comparable with a CAI quoted
  from elsewhere: a gene scoring 0.9 here uses the codons of the average
  *E. coli* gene, not of its most highly expressed ones. The deviation stays,
  because a highly-expressed reference set would have to be chosen and
  shipped and that choice is itself contested, but it is now named in the
  caption on screen, in the `codons.cai()` docstring and in MANUAL.md §9.
  Related and worse: the comment above `ECOLI_CODONS` in `vhh_core.py` said
  the ordering came from highly expressed genes, four lines above naming the
  whole-organism file it actually reads. The ordering itself was verified
  correct against that file, all twenty amino acids.

- **28000 /M/cm and 14000 Da were prefilled in Batch & QC, with no source.**
  Typing an A280 with no sequence loaded produced an mg/mL that looked
  exactly like one computed from a real sequence. The values are plausible
  for a VHH, which is what made them dangerous: on the shipped library the
  medians are 26372 /M/cm and 13836 Da, so they were neither wrong enough to
  notice nor right enough to defend. The fields in Quantify, Conjugate and
  Purification balance are empty now, and each view says what it needs rather
  than assuming it. The same prefill would have set a reagent molar excess in
  Conjugate and scaled every row of the purification table.

- **The concentration from an A280 was not labelled derived.** The A280 is a
  measurement; the mg/mL and uM that follow from it are Beer-Lambert run on
  that reading with an extinction coefficient ProtParam predicts from the
  sequence. The metrics now read "mg/mL (derived)" and "Molar (derived)", and
  the caption underneath says which of the two the number rests on, and says
  it differently when the coefficient was typed by the user rather than
  predicted. Only the high-A260 branch declined to call this a measurement
  before.

- **"The two extinction coefficients differ by about 1%" was true for a VHH
  and wrong by a factor of six elsewhere.** Pairing two cysteines adds
  exactly 125 /M/cm, so the gap is 125 times the number of disulfides and the
  *relative* gap shrinks as the rest of the coefficient grows. Measured
  across the shipped library: 0.5% to 0.9% for the single-disulfide VHHs,
  1.3% for hen lysozyme, 6.6% for human serum albumin with its seventeen
  pairs. The help text and MANUAL.md §5.1 now state the rule rather than one
  sample of it.

- **The Chou-Fasman lysine helix propensity disagreed with its own citation.**
  It read 114, taken from the thegpm.org table, where Expasy ProtScale, the
  MATLAB Bioinformatics Toolbox and the Chou-Fasman teaching tables all
  publish 1.16 citing the same Adv. Enzymol. 47:45-148 (1978) chapter. Every
  other value in the pair table, nineteen helix and twenty sheet, agrees
  across all of them. It now reads 116. This is the one change that moves a
  computed value: on the five-residue smoothed profile the app draws, the
  difference is not visible, but a constant that contradicts its own citation
  is wrong whether or not it shows. The 1978 chapter itself was not read, and
  the comment in `vhh_core.py` says so.

- **Three constants cited sources that could not be confirmed, and now say
  so.** The filamentous-phage optical constant was attributed to "Day's
  relation"; that attribution could not be traced to any primary source, so
  it is dropped in favour of naming the formula and stating plainly that the
  constant's provenance is not established. It circulates widely in exactly
  this form, in course protocols, published methods and vendor calculators,
  always uncited. The side-chain pKa table was attributed to Lehninger 7th
  ed., which was not consulted; the values are the standard ones, six of the
  seven agree across every table checked, and the comment now records that
  cysteine at 8.18 has at least one published table giving 8.08. Neither
  changes a computed result: the pKa table is display-only, and the phage
  screen already shows both titre routes side by side. The Alexa Fluor 488
  figures match the Thermo protocol table they cite, but the A20000 product
  page now lists MW 643.4 and eps 73000; the source string says so, since a
  user comparing the two would otherwise think the app stale.

- **MANUAL.md §3.0 said fifteen analysis views where the code builds
  fourteen.** The families dictionary moved to module scope in
  `sections/analysis.py` and a `view_count()` derives the number from it, so
  the manual is now checked against the screen rather than against a second
  list typed by hand.

Not fixed, and deliberately: the audit could verify only three of the ten
`paper.bib` entries field by field, because the container it ran in cannot
reach Crossref, doi.org, PubMed or PMC. The other seven were checked in the
preceding entry below, by a different route; neither pass is a substitute for
the other, and the bibliography should be re-checked once a route to those
services exists.

### Changed

- **Batch & QC dispatches its views through one table.** The screen chose
  between five views with a chain of calls that each took their own arguments,
  which is what had previously let a view be reached by index. The views now
  share one signature, taking a frozen `Tube` with the selected protein and
  the settings, exactly as the analysis families share a `Run`, and the
  dispatch is a lookup. The table also carries whether a view works on the
  sequence picked at the top, so the two views that supply their own, the
  culture plan and the lot comparison, do not have a picker drawn above them
  that they ignore; previously that would have been an `if` naming views a
  second time.

- **Two tests stopped holding a copy of what they were checking.** The test
  guarding the analysis view count asserted the number was fourteen, which
  made it a third place where that number lived: adding a view meant editing
  the code, the manual and the assertion, and the test written to catch a
  stale number had to be hand-updated whenever the number legitimately moved.
  It now reads the count, spells it, and checks the manual says that and no
  other. The smoke test's map of expected screen titles was written out per
  utility screen and raised a `KeyError` in the middle of a sweep when a
  screen was added, saying nothing about which one; it is derived from
  `shell.SECTIONS` and `shell.UTILITIES`, with the home page and the manual as
  the two deliberate exceptions.

- **The analysis offers fifteen views rather than fourteen**, the new one
  being Session, and MANUAL.md §3.0 follows the count as before.

- **The JOSS manuscript's bibliography is verified rather than provisional.**
  Its ten entries were written at the author/title/year level and carried a
  header saying to check every one before submitting. All ten were checked
  field by field: journal, volume, issue, pages, year and DOI. The nine that
  already had a DOI were right, and none resolved to the wrong paper; what
  was missing was four issue numbers, now added. The tenth, PLAbDab-nano, had
  the institution as its author and `TODO` for the journal and the year; it is
  the Nucleic Acids Research 2025 Database issue article (volume 53, D1,
  D535 to D542, doi 10.1093/nar/gkae881), with the year taken from the issue
  rather than from the advance-access date of October 2024.
  Each entry was confirmed against at least two independent records, normally
  PubMed or PMC plus a publisher or indexer page, and anything the two did not
  agree on was left out rather than guessed: the authors of that tenth entry
  keep their initials, because only one first name could be confirmed in full.
  Two oddities that look like errors are now commented as deliberate in the
  file, since a careful reader would otherwise correct them and make the
  references wrong: the published title of Tien et al. spells
  "accessibilites" with an i missing, and Nature's own web page spells two of
  the Hamers-Casterman authors "Hammers".
  Three tests guard what a compiler cannot: that every key cited in
  `paper.md` has an entry and that no entry is uncited, that no field is still
  a placeholder, and that every entry carries enough to be followed. They are
  derived from `paper.md` rather than from a list of keys typed into the test,
  and they parse the file with a regex rather than adding a dependency for ten
  entries. Not verified: whether each paper is the right one to cite at the
  place it is cited, which no test can answer.

### Added

- **A session file, and a methods paragraph built from it**, as a fourth tab
  of the analysis Outputs family. The app stores nothing, which is why no
  sequence submitted to it can be read by anyone else and also why an
  analysis run in March cannot be looked up in September. The session file is
  the answer that needs no database: the sequences and the three settings in
  one small JSON file, which the same tab reads back, putting the sequences in
  the input and the settings in the sidebar. It holds no computed value on
  purpose. Reopening it recomputes with whatever version is running then, and
  the version that wrote it is in the file, so a number that comes out
  differently is visible instead of being carried forward.

  It refuses rather than half restores: an unknown format number, a setting
  outside what the sidebar can show, a fractional profile window, a sequence
  that is not text. One case restores partly and says so, a session saved with
  the annotation method set explicitly to ANARCI and reopened where ANARCI is
  not installed.

  The methods paragraph is assembled from the same record, so every sentence
  names a setting that was actually used. It names the annotator that ran
  rather than the one that was asked for, which is the claim a user cannot
  check by looking: "auto" resolves to ANARCI where it is installed and to the
  regex anchors where it is not, and the same file reopened on another machine
  can be annotated differently. It says the mass, the pI and the extinction
  coefficients are ProtParam values computed from the sequence and were not
  measured, and it carries the source fingerprint, because a release number
  does not identify a build. Sequence names are left out, since they are
  usually internal clone identifiers and the text is meant for a manuscript.

  The sidebar sliders are keyed for this, and `app.py` now moves any
  `pending_<key>` into place before drawing a widget rather than naming the
  two keys that needed it. The restore itself cannot be driven by `AppTest`,
  which has no way to hand a file to an uploader, so the upload was exercised
  by hand; the mechanism under it, which is the part that would fail silently,
  is covered by a smoke test.

- **A provenance register**, as a third screen below the sidebar divider,
  beside the manual and the reference tables. Every constant the app computes
  with, its current value, what it is worth and where it came from, grouped by
  standing: published, vendor, standard, convention, observed here,
  unverified. That information was already in the repository and effectively
  invisible, sitting in comments and in captions on three different screens,
  and the audit above is what made it worth collecting.

  No value is typed into the register. Each entry names the module and the
  attribute that owns its constant, the value is read live while the page is
  drawn, and a test parses `provenance.py` and fails if a numeric literal
  appears in an entry. The vendor rows are derived from `bench.REAGENTS` and
  the atomic weights take their source from `buffers.ATOMIC_WEIGHTS_SOURCE`,
  so adding a reagent adds a row and no source sentence exists twice.

  What it cannot do, and the page says so: check that the sentence still
  describes the value. The value comes from the code and the text from a
  person, so a threshold that moves without its justification will be shown
  with the old one.

- **How much culture, for how much protein**, as the fifth view of Batch & QC.
  The purification balance runs forwards; this runs the other way, which is
  the question asked before a preparation rather than after it. It divides a
  target by a yield measured on a previous preparation and does nothing else:
  there is no yield model and the yield box has no default, for the same
  reason the balance has no default extinction coefficient. A prefilled mg/L
  would produce a complete plan resting on nothing, and a plan is acted on.

  Three refusals carry the view. A target given twice, once as a mass and
  once as a concentration times a volume, is refused rather than resolved,
  because two targets that disagree mean one belongs to a different plan. A
  recovery given alongside a final yield is refused, because a final yield
  already contains the purification losses and applying them again would
  double the culture while looking entirely reasonable. No yield at all is
  refused with the reason, rather than answered.

- **Across lots**, as the seventh view of Batch & QC. The batch record is the
  one thing this app produces that accumulates, and the question a year of
  them answers, whether the last three lots came out like the first three, is
  not visible in any single sheet. Drop the JSON records in and they line up
  in one table.

  Every derived column is recomputed from the readings with the version
  running now rather than read back from what each sheet said at the time,
  which is the only way a disagreement with an old printout can be seen; the
  version that wrote each record is a column, and a warning names the records
  written by another version. It flags lots that do not all carry the same
  sequence, concentrations that were not all obtained the same way, a lot with
  no usable reading, and two records calling themselves the same lot. A file
  that cannot be read costs its own row and no other. The chart is drawn only
  when every lot was quantified the same way, because a line through two
  methods looks like one measurement drifting, and it is stated to be a line
  between lots and not a trend: a verdict on a lot needs a specification, and
  that is a decision somebody makes.

- **A worked example, in `examples/`,** and the two questions it answers kept
  apart. Three published sequences, one set of settings, and the descriptors
  this code produces from them, checked by `python3 tools/make_examples.py
  --check` and by the test suite. The sequences are not copied into it: the
  manifest names three entries of `data/library.csv` and both other files are
  generated from it, so a correction to the library reaches the example
  instead of leaving it describing the old sequence.

  *Has it changed*: the baseline catches a refactor that moved a number and,
  just as important, a dependency upgrade that moved one. Nothing else in the
  suite would notice, because every other test checks that the code does what
  the code says, and a changed constant is perfectly consistent with that.

  *Is it right*: one descriptor can be checked without leaving the repository,
  and it is the one most results rest on. The molar mass is recomputed from
  the per-residue formulas in `data/aminoacids/meta.json` and the CIAAW 2021
  atomic weights in `buffers.py`, through the formula parser the buffer
  recipes use, and compared with what Biopython reports. The two routes share
  no constant and agree to within 0.19 Da on about 13.9 kDa, which is 14 parts
  per million; on the twenty individual amino acids the same comparison is
  exact for eighteen and differs by 0.007 Da for methionine and cysteine,
  which is a revision of sulfur's standard atomic weight. The isoelectric
  point and the extinction coefficients have no second source in this
  repository, so nothing here checks them, and `examples/README.md` says so
  rather than letting the example imply otherwise.

  The example is annotated with the regex anchors deliberately: ANARCI is
  installed beside the app rather than by it, so a baseline that used it would
  pass or fail according to whether HMMER happens to be on the machine.

- **An enrichment calculator**, as the second view of Phage & panning. One
  row per round, four numbers each: the colonies counted on the plate panned
  on the antigen and on the negative control, and the dilution factor of the
  plate each count came from. The output of a plate is the count times its
  dilution, and the enrichment is the ratio of the two, which is why the
  dilutions are asked for separately: the two plates are rarely counted at
  the same dilution, and dividing raw counts would be wrong by whatever they
  differ by. Giving the volume plated adds cfu/mL; the ratio is unchanged.
  A choice above the table says whether that column holds the exponent
  written on the plate (3 for a 10^-3 dilution, the default) or the factor
  (1000 for the same), with the factor kept for serial dilutions that are
  not powers of ten. The first version of this view asked for the factor and
  named it only in a help text, which is not enough: an exponent typed into
  it does not fail, it answers with the wrong order of magnitude, and 30
  colonies at "3" against 30 at "2" read as an enrichment of 1.5 where 10
  was meant. In factor mode a dilution between 1 and 10 now also raises a
  warning pointing at the other convention, without refusing it.
  A plate with no colony counted is 0 and a plate not counted is empty: a
  half-filled row is left out with a line saying so, rather than read as
  zeros. A negative plate with no colony gives no enrichment at all, with a
  note to count a less diluted one, rather than a division by zero or a very
  large number. Each ratio carries the Poisson floor of its two counts,
  sqrt(1/N+ + 1/N-), stated as a floor and not as the total error. A line
  compares the last round with the first, saying in the same breath that a
  rising ratio is not a selection having worked, since a background binder
  and a faster-growing clone enrich alike.

- **Composing a buffer in Buffers & media**, as a second view beside the
  recipes. Pick compounds from `compounds.csv`, give each a concentration, a
  unit and the form you weigh, and get the same table and the same steps as
  a recipe gives. It goes through the same arithmetic rather than a second
  one: `buffers.custom_recipe` builds a recipe of the shape `prepare`
  already takes. The pH is one you set yourself with the reagent you name,
  defaulting to none, and nothing composed is written anywhere: a buffer
  worth keeping belongs in `data/buffers/` with its source. Compounds with
  no formula offer only the mass and percentage units.

- **A Buffers & media section**, in the Bench group: what to weigh or
  pipette for a buffer or a medium, at any final volume and concentration
  factor, with the form of each salt that is actually on the shelf. A
  component that exists as several hydrates gets a selector, and the mass
  follows from the formula of the form picked: 1.78 g/L of Na2HPO4.2H2O is
  1.420 g/L of the anhydrous salt. No molar mass is typed anywhere; each is
  derived from the formula with the CIAAW abridged standard atomic weights
  (Atomic Weights 2021), copied into `buffers.py` with the date, and a
  formula with an element outside that table is refused. The screen says the
  masses are derived, and that it knows nothing of solubility.
  A recipe describes the 1X buffer and declares whether its pH is adjusted
  (with what) or only measured; a 20X stock is the same recipe at a factor of
  20, so the steps can say not to adjust the PBS stock, and, for a buffer
  whose pH is adjusted, that a concentrated stock does not generally read the
  pH it gives once diluted. Components can come from a stock solution, which
  turns a mass into a volume.
  The recipes are the laboratory's own, thirteen to start with, in three CSV
  files under `data/buffers/`, read-only in the app like the sequence library
  and each with its source. A faulty row is reported with its file and line,
  and a recipe with a faulty component is left out whole.
  `tools/check_buffers.py` runs the same loader before a push. The tests
  check the molar masses against sums worked by hand and the preparation
  against the laboratory's written protocols, amount by amount: PBS 20X,
  TES, NaPi, TB, TB-SALT, glucose, MgCl2, ampicillin, 2X-TY and PEG-NaCl
  come back out to the gram.

- **A Phage & panning section**, in the Bench group, beside Batch & QC.
  The phage titre was a view of Batch & QC, the only one of six that needed
  no sequence: the screen had to choose the view before asking for one and
  skip the sequence block through an early return. It is now a section of its
  own, with the same arithmetic, the same constants and the same widget keys,
  so a titre typed before an update is still there within the session. It
  is where the next panning calculators will go. Batch & QC keeps five views,
  all on the sequence picked at the top. The Bench group's tagline is now
  "The tube, and the readings taken on it", since a phage stock is not a
  purified protein. The flow diagram has a box for the new section, under
  the analysis with no arrow coming in, because it takes nothing from the
  other sections; checked by rendering it in a headless browser on a light
  and a dark theme. The manual's chapters from Sanger sequencing onwards are
  renumbered by one (Sanger sequencing is now §7), and the references to
  them follow.
  The smoke tests read the list of sections from `shell.SECTIONS` instead of
  a list typed in the test file, which a new section would have escaped.
  A new smoke test opens the phage section and checks the two titres and the
  two volumes on the screen's defaults, and another that Batch & QC no
  longer offers the view.

- **A way to say how you use the app.** The app records nothing about what
  is done in it, so the hosting platform's count of visitors is the only
  usage figure there is, and a visit is not a use. A **Using it in your work?
  Say so** panel on the home page, and a **using it? tell us** link in the
  footer, open a GitHub issue pre-filled with the version and three short
  questions: which sections, for what, and what is missing. For anyone
  without a GitHub account, or who would rather not post publicly, the panel
  also offers an e-mail with the same questions. The address is the one in
  `CITATION.cff`, which now carries it, and a test fails if the two differ.
  The issue is pre-filled through the URL, like the existing bug report,
  rather than through an issue template, so the "Report a problem" link opens
  exactly what it opened before.
  The home page panel is collapsed at the bottom of a page people leave as
  soon as they have picked a section, so the request is also one line, in the
  style of the captions around it, under the analysis downloads and under the
  batch sheet: right after the app has given the user something to keep. The
  links there pre-fill which screen they came from. On the home page the
  **Tell us on GitHub** button now has the same weight as **Report a problem
  on GitHub** next to it, rather than being made louder than the section
  buttons.

### Changed

- **`app.py` is a router again, and holds no list of the screens.** It had
  one: an `if section == ...` block per section, each calling that screen's
  `render` with its own arguments, which made it a second list beside
  `shell.SECTIONS`. Adding a section meant editing both, nothing checked the
  router, and a section missing from it is a section that cannot be opened.
  Every screen is now `render(settings)`, one frozen `shell.Settings` carrying
  the pH, the profile window and the annotation method, so the router looks
  the section up in the same dict it has just version-checked and calls it.
  The set of screens comes from `shell.screen_keys()`, which reads `SECTIONS`
  and the new `shell.UTILITIES` (the manual and the reference tables, which
  belong to no group and used to be written out in the sidebar). Adding a
  section is one entry in `SECTIONS`, one file in `sections/` and one line in
  `REQUIRED_MODULES`; three tests fail if a screen is unreachable, missing
  from the version check, or has drifted from the signature.
  The same pass removed what the screens left behind when they moved out of
  `app.py` in the first place: 89 of its 116 imported names were unused,
  including all 62 pulled from `compute`, and about 150 lines of comments
  describing code that is no longer there, under headings like "Cached
  workers" followed by nothing. The file is 559 lines of which 215 were blank.
  No behaviour changed, with one exception worth naming: a `section` key that
  no longer exists, which can only come from a session older than a rename,
  now lands on the home page instead of falling through to the analysis
  screen.

- **The Structure screen no longer computes in the drawing layer.** It called
  `structure.load_structure` directly three times, once in the liability view
  and twice in the 3D view, and annotated the liabilities itself, while
  `compute.py` exists to hold exactly that and cache it. Two cached helpers,
  `st_liabilities` and `st_viewer_html`, now serve both views, so the
  liability view and the viewer's "exposed liabilities" highlight share one
  table instead of each building it. Measured, on a CA-only 128-residue
  synthetic PDB, one parse is 1.37 ms: the parse was never the expensive part
  and this is a consistency fix, not a speed one. The annotation on the same
  screen is a different matter. It ran `classify` and `annotate` on every
  rerun, which under the default "auto" means ANARCI, which means an HMMER
  subprocess for every widget touched on that screen; both now go through
  `classify_cached` and the new `annotation_for`.

- **The README is a front door rather than a second manual.** It carried its
  own account of the in-silico cloning, the Sanger workflow, the synthetic
  generator, what every tab reports and the caveats, all of which `MANUAL.md`
  documents, and the two had drifted: the tab table still described six tabs
  that had become families and said nothing of the four sections added since,
  and the caveat list named neither the structure section nor the buffers nor
  a panning enrichment. Those sections are now a short table of what to read
  where, keeping the two things that are about the repository rather than
  about using the app: how to add a Sanger vector, and the PLAbDab-nano
  attribution. What stays in the README stays because nothing else carries
  it.
  Its file table also kept a line count per file, typed by hand, of which
  eleven of thirteen were wrong, `shell.py` by 478 lines; the column is gone.
  The section table is now generated from `shell.SECTIONS` by
  `tools/readme_sections.py` between two markers, with a test that fails when
  it is stale, so the descriptions a reader of the README sees are the ones
  the home page and the sidebar show, in the same words.

- **Both extinction coefficients in Sequence analysis, and the mass in Da.**
  The Properties table showed only the extinction coefficient with every
  cysteine reduced, and the A280 factor derived from it, while Batch & QC
  offers both and defaults to the cystine value, so the two screens could
  turn the same A280 into two concentrations with nothing saying why. The
  table now shows ε280 with the cystines formed and with all cysteines
  reduced, and the A280 (1 mg/mL) factor for each. `characterize` returns the
  new factor as `A280_0.1pct_cystines`; `A280_0.1pct` keeps its meaning
  (reduced), so a notebook reading it gets the same number as before. The
  mass is shown in Da with two decimals instead of kDa rounded to 10 Da, in
  the Properties table, the construct builder and the PDF report, which also
  carries both coefficients now. `mw_kda` is still returned.
  The first build of this change showed pI and the net charge with six
  decimals: a second `Styler.format` call with a dict resets the columns it
  does not name to pandas' default. The table is now formatted in one call,
  two decimals unless a column says otherwise, and a test reads the
  displayed values.
- **One identity heatmap at a time.** The Alignment tab drew the full-length
  and the CDR3 identity matrices side by side; an **Identity over** selector
  now shows one, and both stay downloadable as CSV. The caption says which
  matrix answers which question: CDR3 tells whether two clones are the same
  binder, full-length separates variants of one clone that differ only in
  their frameworks. The PDF report still carries both.
- **The two CDR3 identities are now explained in the app.** The Alignment
  tab and Clonal families both compute a CDR3 identity matrix, with different
  gap penalties (−12/−1 and −11/−1) and different denominators (alignment
  length, and length of the longer CDR3, kept from the Sanger notebook), so
  their CSV files could give different values for the same pair with nothing
  in the app saying why. Both tabs now carry the same note, written by one
  function, with the values the two definitions give for one example pair, a
  CDR3 shifted by two residues: 73.3% against 84.6%. The values are computed
  by the two functions when the page is drawn, and a test fails if the
  example stops differing or if two contained CDR3s stop agreeing. On 396
  generated pairs where one CDR3 lies within the other the two agreed every
  time; on 1500 random pairs, in 966. Both definitions are kept.

### Fixed

- **The reference comparison said `[Errno 2] No such file or directory:
  'hmmscan'`.** It is the one view in the app with no fallback: positions are
  matched by IMGT number, only ANARCI produces one per residue, and ANARCI
  shells out to HMMER, which no `pip install` brings in. The reference's own
  numbering happens outside the per-sequence `try` in `compare_to_reference`,
  so a missing binary did not become a per-sequence warning, it reached the
  page as the operating system's message about a file nobody has heard of.
  The view now checks first and says which half is missing, the package or
  the binary, with the command for each platform, and that the rest of the
  app falls back to the regex anchors while a position-by-position comparison
  cannot. `shell.anarci_requirement()` is what tells the two apart, including
  the case where both are installed and numbering still fails, which usually
  means ANARCI's germline HMM files are absent.
  The person most likely to have met this is the one who followed the app's
  own privacy notice: it sends confidential work to a local clone, and
  `pip install -r requirements.txt` installs ANARCI without HMMER. The suite
  now asserts in both directions, from `conftest.HAS_HMMER`: with HMMER the
  comparison must still compute, without it the screen must name what to
  install. Checked both ways, by running the suite with HMMER installed and
  again with `hmmscan` moved out of the path.

- **A broken bullet in the README's data-handling section.** An edit had
  folded the alignment and the uploads entries into one another, parenthesis
  and all, leaving a sentence that read as neither. The cache knobs it names
  were also still said to live in `app.py`; they have been in `shell.py`
  since the chrome moved there.

- **The PDF summary table could put a header over the wrong column.** It
  cut its list of headers to the number of columns present, so a missing
  column shifted every header after it by one. Each column now carries its
  own header. Not seen in a real report, since `characterize` always returns
  every column; found while adding the second extinction coefficient.
- **Opening Comparison with 50 sequences restarted the server.** The
  pairwise identity heatmap grew by 1.1 inch a sequence and wrote a number in
  every cell: at 50 sequences, the app's limit, that is a 58 x 47 inch figure
  of about 110 million pixels at the app's 200 dpi, with 2500 labels. Driving
  the app through `AppTest` with 50 synthetic VHHs, the Comparison family
  alone peaked at 5671 MB of memory and took 53 s, on a Community Cloud
  container that can have 690 MB; the server was killed and started again,
  which in the log shows only as a fresh `Uvicorn server started` line.
  The figure is now capped at 14 inches on its longest side, and above 20
  sequences the cells carry no number, since they would be too small to read;
  the full-length and CDR3 matrices can be downloaded as CSV under the
  heatmaps instead. `png_identity` was also the only figure helper in
  `compute.py` without a cache, so every rerun of the family redrew both
  heatmaps; it is cached now, and a test fails if a `png_` helper is added
  without a cache. Same measurement after the change: 646 MB peak, 23 s.
- **The charge and hydropathy profiles drew every sequence at once.** One
  2.4 inch panel per sequence made the figure grow with the set: rendering it
  took about 13 MB for one sequence, 106 MB for 10 and 516 MB for 50,
  measured, and with 50 sequences the Identity family peaked at 1089 MB. The
  Plots tab now works like Chou–Fasman: a **Sequence** selector draws one
  profile, and **Show all** draws one image per sequence, one after the other.
  The PNG download keeps the single figure up to 10 sequences
  (`compute.PROFILES_COMBINED_MAX`) and becomes a ZIP with one PNG per
  sequence above that; the PDF report follows the same rule. Identity with
  50 sequences, same measurement: 407 MB peak, 6 s instead of 19 s.
- **The PDF report could not be built.** Pressing **Build PDF report** raised
  `NameError: name 'VERSION' is not defined`, whatever the input:
  `VERSION`, `AUTHOR` and `APP_NAME` were used in the Download tab but never
  imported from `shell`, and with the secondary-structure plots included the
  Chou–Fasman window was read through a variable local to the Plots tab,
  which would have raised next. Every screen rendered, and no test pressed
  the button. The window is now read from the session, and a smoke test
  builds the report with and without those plots and checks that a PDF comes
  out. With 50 sequences the report builds, in about 30 s and 7 MB, but the
  build peaked at about 870 MB in the same `AppTest` measurement, mostly
  from the two identity heatmaps and the alignment; that is not addressed
  here.
- The manual said the downloaded figures were 300 dpi PNG; they are
  200 dpi, the resolution `shell.fig_to_png` renders every figure at.
- **The home page diagram was unreadable at night.** The diagram and the
  amino acid drawings sit in iframes, which inherit no colour from the page,
  and chose their text colour from `prefers-color-scheme`, the operating
  system's setting, while the page around them follows the app's theme. The
  app was always light (see the next entry), so on any device whose system
  had switched to dark mode, at night for many, the text turned near-white
  on a white page. Choosing Dark in the menu on a light system, possible
  from this release, would have done the reverse. Both iframes now take
  their text colour from the app instead, through
  `shell.follow_app_theme()`: the initial colour comes from
  `st.context.theme`, which is right on a first load but stale after a
  switch until the next rerun, and a small script then reads the page's
  background and follows a switch at once. Verified in a headless browser
  with the system light and dark, crossed with the app theme light, dark
  and following the system, on the home page diagram and on the twenty
  amino acid drawings.
- **The app was stuck in a light theme, and flashed dark first on a dark
  system.** `.streamlit/config.toml` set the accent colour under a bare
  `[theme]`, which makes Streamlit build one custom theme on the light one:
  the System / Light / Dark choice was missing from the app menu, and with
  the system in dark mode the page painted dark, then turned white about
  0.3 s later. The accent colour is now set under `[theme.light]` and
  `[theme.dark]`; the menu offers the three choices and the page loads in
  the system's theme without the flip. Measured with Streamlit 1.66.1 in a
  headless browser, by sampling the page background every 0.15 s during
  the load, before and after the change. A test fails if a colour is set
  under a bare `[theme]` again.

## [1.1.0], 2026-11-04

A release built around the bench: a new Batch & QC section, which turns
readings into concentrations, titres and recoveries and records a purified batch
on a one-page sheet, kept apart from what the sequence only predicts. The
sections are now grouped on the home page by what they work on, the footer
says which files and which library versions the server is running, and the
dependencies are pinned. No function signature in `vhh_core.py` changed, so
code importing it from a notebook is unaffected.

Several of the fixes below concern things added in this same release, among
them the footer fingerprint, the batch sheet and its record, and the upload
settings, so 1.0.0 never had those defects. They are kept because they record
why the code is the way it is.

### Added

- **A Batch & QC section**, in the Bench group: the arithmetic between a
  purified protein and the next experiment. Concentration from a measured
  A280 using the sequence's own extinction coefficient, with the 260/280
  ratio when you give an A260 and a dilution that refuses to concentrate. A mass check that compares a measured
  mass against the sequence and lists the modifications that would explain a
  difference *and that the sequence allows*, each delta computed from a
  chemical formula rather than tabulated. And a conjugation planner that works
  out how much NHS-ester reagent to add, separates reconstituting the vial from
  the working stock, offers a dilution when the volume is too small to pipette
  accurately, and warns when lysines sit in a CDR.
- **A phage titre calculator**, as a view of Batch & QC. An A260, the dilution
  it was read at and a glycerol correction turn into particles per millilitre,
  and into the volume of stock that holds the number of particles you want,
  for a solid-phase round and for a liquid-phase one beside it, which is twice
  as much. Both routes to a titre are shown side by side: the empirical
  constant from a titration, which is editable because it belongs to a
  preparation and an instrument rather than to physics, and Day's relation
  from the genome length. On the same reading they differ by a factor of about
  275, and the app says so instead of picking one. The glycerol correction is
  asked for as the two volumes, because the factor is below one when the
  absorbance was read before the glycerol went in, and that is the direction
  this protocol always goes.
- **A batch sheet**, as a view of Batch & QC: a one-page PDF for a purified
  batch, with what you wrote down, what you measured and what the sequence
  implies in three blocks that are kept apart, so a predicted extinction
  coefficient cannot be read a year later as a measurement. It is a laboratory
  batch record, not a certificate of analysis, and the PDF says so.
- **The batch sheet saves a machine-readable record** beside the PDF, and
  takes it back. The JSON file holds the inputs, the text and the readings,
  which are the part that cannot be recomputed, plus the page as it was
  rendered so a reprint can be compared against the original. Loading one
  refills the form; the derived values are recomputed from the readings rather
  than restored, and the view says so when the record was made with another
  version of the app or for another sequence. A record whose format number or
  whose readings this version cannot read is refused with the reason, rather
  than partly read into a sheet with fields silently missing. Still no
  server-side storage: the file goes to whoever generated it and comes back
  from them, it is not an audit trail, and nothing in it can tell whether it
  has been edited.
- **A purification balance**, as a view of Batch & QC. One row per step of a
  preparation, with a volume and a concentration given either as an A280,
  divided by the extinction coefficient and therefore derived, or as a reading
  from a method that does not go through the sequence, which takes precedence
  and is named in the table. It reports the total in mg per step, the recovery
  against the step before it, the cumulative recovery and, given a culture
  volume, the yield per litre. Steps marked **Crude**, a lysate or a
  clarified supernatant, keep their total and stay out of the recovery
  arithmetic: an A280 there is every protein in the tube, so a ratio against a
  purified pool is not the target's recovery, which would have been the most
  flattering number on the page. It says so rather than printing a number when
  a step recovers more than it was given, when the steps compared were not
  quantified the same way, and when an A280 arrives with no extinction
  coefficient to divide it by. Nothing predicts what a step ought to recover.
- **AlphaFold3 confidence in the structure section.** The structure upload now
  takes the `.zip` AlphaFold3 produces, as downloaded from AlphaFold Server or
  zipped from a local run. Pick any of its models; the confidence view then
  shows ipTM (whole complex and for the chosen VHH–partner pair), pTM, the
  ranking score, the clash flag, a table comparing every model in the archive,
  and the PAE matrix. The interface view warns when the pair's ipTM says the
  arrangement it is measuring has likely failed. A `.cif` on its own still
  works and gives pLDDT only, which says nothing about where the chains sit
  relative to each other.
- **A Reference screen**, in the sidebar under the manual: the genetic code
  and the twenty amino acids, with the structures drawn as they exist at
  pH 7.4. Nothing on it is transcribed. The codon table is Biopython's, the
  usage figures and rare-codon flags are the ones the gene diagnosis screen
  uses, and the masses, hydropathy and charges are the constants the analysis
  screen builds its profiles from, so the reference and the app cannot
  disagree. The structures are committed SVG files built once by
  `tools/aminoacids/build.py`, which checks each formal charge and molecular
  weight before writing, so that the app need not ship RDKit for twenty
  pictures that never change.
- The home page explains how the sections feed each other with a diagram rather
  than four paragraphs, drawn as SVG so it follows the light and dark themes.
- The manual is browsed one section at a time, with a **Jump to** selector,
  previous/next buttons and a list of what each section contains. Rendering
  seventeen hundred lines on every interaction was slow and nobody reads a
  reference manual from the top.
- **The footer says what the server is running**: a line under it reads, for
  example, `Running on Python 3.14.6 · Streamlit 1.66.1 · Biopython 1.88 ·
  ANARCI 2026.2.13.2`, taken from the running process. The fingerprint says
  which files were deployed; this says which interpreter and which library
  versions are executing them, including the interpreter, which a lockfile
  cannot pin. ANARCI reads `absent` when the
  package is not installed.
- **Continuous integration.** A GitHub Actions workflow runs the full suite on
  every push and pull request, and once a week on a schedule, on Python 3.14
  and from `requirements.txt`, which is the interpreter and the dependency set
  Community Cloud runs. HMMER comes from `packages.txt` read by the workflow
  rather than named again in it, so the system package list stays in one
  place. The weekly run is there because the lockfile pins the Python packages
  and not the interpreter patch release or the apt packages under them. The
  suite had been run only on one machine, which the 1.0.0 changelog listed as
  a known limitation; it no longer is. Verified: the lockfile installs and all
  158 tests pass in a clean virtual environment built from
  `requirements.txt` and `requirements-dev.txt` alone, and `pip check` is
  clean there. Not verified: the workflow itself, which only a push can run,
  and Python 3.14, which was not available where this was checked (the clean
  install and the suite were run on 3.13). Since then the full suite has been
  run on Python 3.14.6 from the lockfile, see *Dependencies are pinned*.
- **A dev container**, `.devcontainer/devcontainer.json`, for GitHub
  Codespaces or VS Code. It installs HMMER from `packages.txt`, then
  `requirements.txt` and `requirements-dev.txt`, and starts the app. The
  file as first committed was Streamlit's template, on Python 3.11, where the
  lockfile cannot install: contourpy 1.4.0 requires Python 3.12 or later,
  which `uv pip install --dry-run` on 3.11 confirms. It now starts from
  Microsoft's `python:3-3.14-bookworm` image, a variant listed in that
  image's manifest, and no longer installs Streamlit a second time outside
  the lock. Not verified: building the container, which could not be done
  where this was checked.
- `pytest` now says on its first line when HMMER is missing, with the command
  to install it. ANARCI shells out to the `hmmscan` binary, which no Python
  requirements file can provide, and without it the smoke test failed with
  `[Errno 2] No such file or directory: 'hmmscan'` inside an assertion about a
  Streamlit element in the Comparison view, which is accurate and says nothing
  about what to install. It warns rather than skips: a suite that skipped the
  annotator tests would pass on the one machine where the app's default
  annotator is never exercised.

### Changed

- **The sections are grouped by the molecule they work on.** The home page
  and the sidebar now sort the seven sections into three groups, **Protein**
  (Sequence analysis, Structure, Sequence library), **DNA** (Sanger
  sequencing, Gene diagnosis, In-silico cloning) and **Bench**, each under a
  title in a larger type than the section names and with a colour of its
  own, which the flow diagram on the home page also uses to tint its boxes,
  with a legend. On the home page, Protein and DNA sit side by side and
  Bench takes the full width below them. The groups are where to look, not a
  workflow: the diagram still shows how data moves between them. Each entry
  of `SECTIONS` in `shell.py` names its group, and the home page, the sidebar
  and the diagram all derive the grouping from there through
  `sections_by_group()`, so a section cannot be in one place and missing from
  another; a section naming no group, or a misspelt one, raises rather than
  disappearing from the navigation. The cards on the home page are stacked
  per group rather than laid out row by row: row by row lined them up on a
  wide screen, but on a phone, where Streamlit stacks the columns, it
  interleaved the two groups. Verified by rendering the app in a headless
  browser at 1400 and 420 px wide, and the diagram on a light and a dark
  background; by new tests in `tests/test_navigation.py` and a smoke test,
  which fails if the group titles are missing or the order of the buttons
  stops following the groups (checked by breaking each on purpose).
- **The bench section is called Batch & QC.** It was called Bench while this
  release was being written; the group now has that name, and a section under
  it with the same name read as a duplicate. Its title, sidebar button, home page card and diagram box all take the new
  name from `SECTIONS`: the screen title and the diagram titles, written out
  by hand until now, are derived from it, so a later rename happens in one
  place. The internal key stays `bench`. The home page buttons keep the
  label's own case, "Open Sanger sequencing" rather than "Open sanger
  sequencing".
- **The Protein and DNA groups have the same height on the home page.** Each
  group fills its column and its cards share out the slack, so the two
  columns end at the same level instead of one stopping short. The bottoms
  line up; the tops of facing cards line up only when the cards above them
  are the same height. Unchanged on a phone, where the columns stack.
  Verified by rendering at 1400 and 420 px.
- **The Open buttons on the home page take the colour of their group,**
  with white text. To keep that text readable, the group green and orange
  are darker than the palette colours they come from, #438452 instead of
  #55a868 and #bb5b25 instead of #dd8452, same hue and saturation: white on
  them now gives 4.51:1 and 4.53:1, against 2.92:1 and 2.80:1 before, and
  WCAG 2.1 AA asks for 4.5:1 for normal text (white on the blue is 4.84:1).
  The darker colours apply everywhere the group colour does, titles, sidebar
  and diagram included. The ratio is computed by `shell.contrast_ratio`,
  following the WCAG 2.1 glossary, whose linearisation threshold moved from
  0.03928 to 0.04045 in May 2021; no 8-bit channel value lies between the
  two, so the results are the same either way. A test fails if any group
  colour drops below 4.5:1 with white. The buttons are styled through the
  `st-key-home_<section>` CSS class Streamlit gives keyed widgets; if a
  future Streamlit drops it, they fall back to the theme colour. Verified in
  a headless browser, hover included; the tests check the CSS rules and the
  contrast, not what the browser paints.
- **Dependencies are pinned.** `requirements.txt` is now a lockfile generated
  by pip-compile from `requirements.in`, which keeps the ranges it used to
  hold: every package Community Cloud installs, transitive ones included, is
  at an exact version, so a rebuild installs what was tested rather than
  whatever is current that day, and a dependency can no longer change under
  an unchanged footer fingerprint. The lock is resolved for Python 3.14, the
  interpreter of the deployed app, and some of its pins (contourpy 1.4.0
  among them) do not install on 3.11, so the `Dockerfile` moves to Python
  3.14 too, and takes ANARCI from the lockfile rather than from bioconda. The
  update routine is in the README, *Updating dependencies*. Verified: the
  lock installs and the full suite passes in a fresh Python 3.14.6, and in
  a micromamba 1.5 environment built like the image (Python 3.14 and HMMER
  from conda, then `requirements.txt`) running on the files the image
  receives. Not verified: the Docker image itself, whose base image could not
  be pulled where this was checked, and a Community Cloud rebuild from the
  lock, which only a push will show.

### Fixed

- **The footer fingerprint did not cover most of `data/`.** It hashed three
  entries written out by hand, `data/library.csv`, `data/aminoacids/meta.json`
  and the SVG structures, while the app also reads `ecoli_codon_usage.csv`
  (the CAI table and the usage figures on the reference screen),
  `vhh_position_freq.csv` with its metadata (the germline comparison),
  `primers.csv` and the four plasmid maps under `data/plasmids/` (the cloning
  screen). A corrected codon usage table, a changed primer or a changed
  plasmid map therefore deployed under a footer identical to the one before
  it, which is the one thing the fingerprint exists to prevent. `data/` is now
  walked rather than listed, so a file added to it is covered without editing
  anything, and caches such as `__pycache__` and `.DS_Store` are skipped so
  that two checkouts of one commit still agree. The fingerprint also hashes
  each file's path relative to the root rather than its bare name, so moving a
  file between two hashed directories no longer leaves it unchanged. Verified
  by editing each of the ten files `data/` ships, one at a time, on a copy:
  before the fix, five of them left the fingerprint at `944ee32`. A test now
  walks the real `data/` directory and fails if any file in it can be changed
  invisibly. The fingerprint moves once because of this change, which is
  expected: it is the same files hashed differently, not a change to the app.
- `batch_record_inputs` accepted `"fields": []` and `"measured": []`. The
  type check sat after an `or {}` fallback, and an empty list is falsy, so it
  became an empty dict and passed the check written to catch exactly that: the
  record reprinted with every text field blank and nothing saying why. Found
  by the test for it, which failed on the first run.
- Three numbers were formatted with a comma as the thousands separator, among
  them the mass difference on the batch sheet, which would have read
  `+1,024.5 Da` for a difference past a thousand daltons. A comma means a
  decimal point to most of the people this app is written in English for.
- Ten calls still passed `use_container_width=True` directly although
  `shell.FIT` exists to resolve that argument for the installed Streamlit.
  Streamlit says it will remove the argument after 2025-12-31, which has
  passed, and 1.66.1 logs a deprecation warning on every one of those calls.
  They now go through `FIT` like every other call site. Verified: the warning
  no longer appears when the app is driven through `AppTest`.
- `families_for` in `compute.py` carried two `st.cache_data` decorators, the
  upper one copied down from the figure helper above it, so the clonal
  families and their linkage matrix, the largest object this module caches,
  were stored in two caches at once on a server that can have as little as
  690 MB. Results were correct, so nothing failed. A test now reads the
  module with `ast` and fails if any function carries more than one cache
  decorator. The memory this frees was not measured.

- **The Docker image could not have started.** The `Dockerfile` copied
  `vhh_core.py` and `app.py` and nothing else, while `app.py` imports a dozen
  other modules and reads `sections/`, `data/`, `MANUAL.md` and
  `.streamlit/config.toml`; the first import failed. It now copies the whole
  project, and a new `.dockerignore` says what stays out (`.git`, caches,
  zips, `results/`, `paper/`, `tools/`), so the file no longer has to be
  edited each time a module is added, which is how it went stale. Verified:
  the build context, by building it with Docker, and the app, by installing
  the image's environment with micromamba 1.5 (Python 3.11, HMMER and ANARCI
  from bioconda, then `requirements.txt`) and running it, health check and full
  test suite, from exactly the files the image receives. Not verified: the
  image itself, since its base image could not be pulled where this was
  checked.
- **The footer fingerprint now covers `requirements.txt` and `packages.txt`.**
  Community Cloud installs both, so editing either changes the deployed app,
  and the fingerprint did not move. It still does not move when a dependency
  releases a new version under an unchanged file, because the versions are
  minimums, not pins (since fixed, see *Changed*). The README and the manual described the fingerprint as
  a hash of the `.py` files only, which was already incomplete; both now list
  what it covers and what it does not.
- `requirements-dev.txt` was missing `pypdf`, which `tests/test_bench.py`
  needs to open the batch sheet it generates. The test passed only where
  `pypdf` happened to be installed already.

- **File uploads on the deployed app returned 500 while working locally.**
  `.streamlit/config.toml` turned XSRF protection on; Streamlit Community
  Cloud runs every app it hosts with it off, and with it on the upload goes to
  a different origin than the one that issued the token, so the upload
  endpoint refused the request before it reached any application code. There
  was nothing in the logs, and the only symptom was "AxiosError: Request
  failed with status code 500" on dropping a structure, a chromatogram or a
  FASTA into any section. A test now fails if the setting comes back.
- `st.components.v1.html`, which Streamlit has deprecated and says it will
  remove, is no longer called directly. The home diagram, the amino acid
  structures and the 3D viewer all go through one helper that uses `st.iframe`
  where it exists and the old call where it does not, so the app works on both
  and the deployed one, which installs whatever Streamlit is current, does not
  break the day the old call goes. A test fails if a second call site appears.
- **A page loaded during a deploy no longer shows a raw traceback.**
  Community Cloud swaps the source files under a running server, and a run
  that asked for a module mid-swap raised `KeyError` (seen in the console as
  `KeyError: 'generator'`). A first fix caught it only around the module
  version check, which runs after the imports at the top of `app.py`, so the
  import that actually raised was still unprotected. Every project import is
  now inside the guard, which says the app is updating and names the error
  underneath in small type, so that a file that is really missing is not
  taken for an update forever. Tested by making the import of `generator`
  raise inside Streamlit's `AppTest`, which reproduced the original traceback
  on the old `app.py`; the race against a real deploy was not reproduced. A
  test also fails if a project module is imported outside the guard.
- `.streamlit/config.toml` is now in the repository. It never was, so the
  deployed app ran with Streamlit's defaults and the 10 MB upload cap the
  manual quotes applied only on a machine that happened to hold the file.
- The upload cap is 20 MB rather than 10, so that the structure section's own
  200 000 atom limit is what refuses a file rather than a byte count: at
  10 MB that limit could never be reached. It stays deliberately tight, a
  parsed structure costs 0.6 to 1.2 kB an atom against a server that has
  between 690 MB and 2.7 GB.
- `.gitignore` no longer excludes `*.csv`. The app ships four CSV files it
  cannot run without; the ones already tracked were safe, and the rule would
  have silently refused the next one added.

## [1.0.0], 2026-10-28

First public release, and the first that is citable: MIT licence, citation
metadata, and a Zenodo archive.

### Added

- **Structure section.** Upload a PDB, mmCIF or AlphaFold model and replace the
  sequence-based guesses with measurements: SG–SG distances instead of a
  canonical-pair assumption, relative solvent accessibility, which flagged
  liabilities are buried and which are exposed, hydrophobic patches grouped by
  spatial contiguity rather than by a sliding window, and (for a complex) the
  ΔSASA interface and epitope. Includes a 3Dmol.js viewer.
- **Sequence library.** A curated, read-only list of 22 published sequences kept
  in `data/library.csv`, each with its provenance: VHHs from the PDB
  (cAbLys3, caplacizumab, Nb80, anti-GFP Enhancer and Minimizer, NbALFA, Ty1,
  7D12, 2Rs15d, Nb.b201, cAbBCII10, NbSyn2) and their antigens.
- **Gene diagnosis section.** Paste coding sequences that already exist: frame
  and stop checks, documented rare codons and their clusters, GC-rich
stretches,
  homopolymers, CAI against a measured *E. coli* usage table, and synonymous
  swaps that lower GC without changing the protein or dropping below a usage
  floor, so a GC problem is not traded for a rare-codon one.
- **Triage table.** Every flag the app raises, one row per sequence, sortable.
  No composite score, deliberately.
- Clonal family grouping of panning output by CDR3 identity, at a threshold you
  set, with a dendrogram and one representative per family.
- Aggregation profile, humanness against human germline V genes, detection of an
  additional CDR1–CDR3 disulfide, a construct builder, and suggested repairs
for
  chemical liabilities.
- Sequence logo reproducing the alignment notebook, with adjustable wrap.
- Any protein is accepted and classified; conventional VH domains receive the
  full heavy-chain analysis, labelled as VH.
- `LICENSE` (MIT), `CITATION.cff`, and an automated test suite.

### Changed

- **The analysis is grouped into six families chosen with a selector** instead
  of fifteen tabs. Only the selected family executes (Streamlit runs the body
of every tab, nested or not) so first display on seven sequences drops from
  roughly 12 s to under 2 s. The expensive views are still expensive; they are
  no longer paid for by everyone who never opens them.
- ANARCI is called once per sequence and cached, instead of seven times.
- `app.py` was split into a router, a chrome module, a cached-computation module
  and one module per screen. It went from 3 911 lines to 459.
- The footer now carries the version, the licence, the source and the DOI.

### Fixed

- Liabilities, hallmarks, humanness and disulfides are no longer computed for
  proteins that are not variable domains, where they produced confident
  nonsense, hen lysozyme's eight ordinary cysteines were reported as eight free
  thiols.
- Suggested repairs no longer propose substituting the structurally conserved
  positions IMGT 23, 41, 104 and 118: removing the FR2 tryptophan trades a
  chemical risk for a protein that does not fold.
- `apply_repairs` refuses two conflicting substitutions at one position outright,
  instead of applying whichever came first in the list.
- The structure-to-sequence mapping is keyed by residue number rather than by
  sequential index, so CDR bands and region labels land on the right residues
in
  a chain that does not start at 1 or that has an unmodelled loop.
- GC-rich windows are merged into regions instead of each being reported
  separately, 98 "extreme windows" on one optimised gene were 98 views of a
  handful of stretches.
- A liabilities table with no rows keeps its columns, instead of failing
  downstream with `KeyError: 'score'`.
- The CDR-length comparison is readable with a single sequence.
- The footer is drawn once per screen rather than twice.

### Known limitations

- `reverse_translate` balances GC across the whole sequence rather than within a
  window, so a back-translated gene can still contain local GC-rich stretches,
  which the gene diagnosis section will then flag.
- Mixed Sanger chromatograms are not detected. A detector was written and
  withdrawn: without a per-position quality condition it flagged read edges,
  where Phred scores drop, and a false positive on a clone identity call costs
  more than no detection.
- The test suite runs locally; there is no continuous integration.
