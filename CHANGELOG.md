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
