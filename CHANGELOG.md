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
