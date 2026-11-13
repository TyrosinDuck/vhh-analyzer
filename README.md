# VHH tool box

[![tests](https://github.com/TyrosinDuck/vhh-analyzer/actions/workflows/tests.yml/badge.svg)](https://github.com/TyrosinDuck/vhh-analyzer/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23742269.svg)](https://doi.org/10.5281/zenodo.23742269)
[![licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)

A web application for working with VHH / single-domain antibodies, from the
sequences a panning campaign ends with to the tube the purified protein ends
up in: physico-chemical descriptors, CDR/framework annotation (ANARCI or regex
anchors), developability liabilities, alignment and clonal families, measured
properties from a structure you upload, codon-optimised DNA and in-silico
cloning, and the bench arithmetic around a batch (concentration from an A280,
buffers, phage titres, a batch sheet).

```
vhh-analyzer/
├── app.py             the router: version check, sidebar, one call per screen
├── shell.py           the chrome shared by every screen, and the cache policy
├── compute.py         every cached computation, and nothing that draws
├── sections/          one module per screen; a screen draws, it does not compute
├── vhh_core.py        the analysis, importable without Streamlit
├── structure.py  bench.py  buffers.py  cloning.py  codons.py  ...  the rest of it
├── data/              the sequence library, codon usage, buffer recipes, reference SVGs
├── tests/             pytest, including smoke tests that drive the whole app
├── paper/             the JOSS manuscript draft and its bibliography
├── requirements.in    pip dependencies, as ranges: the file you edit
├── requirements.txt   generated from it by pip-compile, every version pinned
├── packages.txt       apt dependencies (HMMER, the ANARCI backend)
├── Dockerfile         for a container host, if you ever outgrow Community Cloud
└── .streamlit/
    └── config.toml
```

The analysis modules are the single source of truth. A notebook and the web app
both import them, so a fix to a descriptor propagates to both instead of
drifting apart, and that is also why changing a `vhh_core.py` signature is a
MAJOR release even when the app looks identical.

## Documentation

A full user manual lives in `MANUAL.md` and is rendered inside the app,
reachable from
the home page and from the sidebar. It documents every screen, what each number
means,
and how far to trust it. Edit the Markdown file and the app picks it up, there
is one
source of truth.

The **Provenance** page, in the sidebar beside the manual, lists every constant
the app computes with, what it is worth (published, vendor, standard,
convention, observed here, unverified) and where it came from. The values are
read from the modules that own them while the page is drawn, so it cannot show
a number the app is not using.

`examples/` holds a worked example: three published sequences, one set of
settings, and the descriptors this code produces from them, with
`python3 tools/make_examples.py --check` to compare. It is a baseline rather
than an independent reference, except for the molar mass, which the test suite
recomputes from the per-residue formulas and the CIAAW atomic weights in this
repository and finds to agree within 14 ppm. `examples/README.md` says exactly
which descriptors have been checked against a second source and which have
not.

### Bundled data

`data/library.csv`, the sequence library. Edit it locally, commit and push; the
app serves it read-only. It is public, so nothing confidential goes in it.
For your own sequences, the app's **My library** tab loads a file in the same
format, lets you edit it, and gives it back as a download; it keeps nothing.

`data/ecoli_codon_usage.csv`, *E. coli* codon usage (Kazusa Codon Usage
Database, taxid 316407), redistributed from the `python-codon-tables` package
under **CC0-1.0**, a public-domain dedication. Used for the CAI, the usage
floor on codon substitutions, and the back-translation preference order.

`data/sabdab/`, a snapshot of the single-domain entries of **SAbDab2**
(Oxford Protein Informatics Group), redistributed under **CC-BY 4.0**; see
`data/sabdab/LICENSE`. It drives the closest-structure search. Derived from
SAbDab2's own export by `tools/sabdab/build_index.py`, never edited by hand,
and carrying the extraction date the app displays. Rebuild it by downloading
a fresh export and rerunning that script.

## Files

| File | What it is |
|---|---|
| `app.py` | the router: version check, sidebar, and one call per screen |
| `shell.py` | the chrome: identity, navigation, `footer()`, caching policy |
| `compute.py` | every cached computation, and nothing that draws |
| `sections/analysis.py` | the Sequence analysis screen, six families |
| `sections/structure.py` | the Structure screen |
| `sections/bench.py` | the Batch & QC screen |
| `sections/buffers.py` | the Buffers & media screen |
| `sections/phage.py` | the Phage & panning screen |
| `sections/cloning.py` | the In-silico cloning screen |
| `sections/genes.py` | the Gene diagnosis screen |
| `sections/sequencing.py` | the Sanger screen |
| `sections/library.py` | the Sequence library screen |
| `sections/home.py`, `sections/manual.py`, `sections/reference.py` | the landing page, the manual and the reference tables |
| `sections/provenance.py` | the provenance register: every constant and where it came from |

Plus the science, which the screens only arrange: `vhh_core.py`, `cloning.py`,
`structure.py`, `codons.py`, `bench.py`, `buffers.py`, `repertoire.py`,
`sanger.py`, `generator.py`, `constructs.py`, `triage.py`, `library.py`,
`reference.py`, `af3.py`. Plus two files that describe the app rather than
compute with it: `session.py`, which writes the inputs of a run to a file so
the run can be repeated, and `provenance.py`, which names the module and
attribute that owns each constant so the register reads live values instead of
a second copy.

This table used to carry a line count per file, typed in by hand. It drifted,
which is what a hand-kept count of something the repository already knows does:
eleven of its thirteen numbers were wrong, `shell.py` by 478 lines. `wc -l` is
one command and is never out of date.

`app.py` was a single 3 911-line script until the screens were pulled out of
it.
Three rules keep it from growing back:

- **A screen module draws; it does not compute.** Anything cached belongs in
  `compute.py`, where the whole caching policy can be read in one file instead
of
  found in a few dozen decorators spread through three thousand lines of UI.
- **The router holds no list of the screens.** Every screen is
  `render(settings)`, one signature, so `app.py` dispatches by looking the
  section up rather than by an `if` per section. It used to carry one, which
  made it a second list beside `shell.SECTIONS`: adding a section meant
  editing the router too, and nothing checked the router.
- **Nothing reads a global from another screen.** Every screen is a `render()`
  function, which is not tidiness: in the old module-scope version the
Structure
  screen's `kinds` assignment silently overwrote the analysis screen's variable
  of the same name, and the analysis family selector was shadowed by a
  same-named germline selectbox. Neither ever broke anything, because two
screens
  never render in the same run, a function scope removes the class of accident
  rather than its instances.

Each module carries a `__version__`, and `REQUIRED_MODULES` in `app.py` lists
the
version it expects. A mismatch stops the app with a message naming the file,
and
`test_module_versions_match_app` catches it before a deploy does.

## Layout

The app opens on a home page with nine sections in three groups, Protein,
DNA and Bench, also reachable from the sidebar, which carries a **Home**
button on every screen. The table below is generated from `shell.SECTIONS`,
which is also what the home page cards and the sidebar read, by
`python tools/readme_sections.py`; a test fails when it is stale, so do not
edit it by hand. Every screen is documented at length in `MANUAL.md`, §2 for
getting around and one section each after that.

<!-- generated: sections, by tools/readme_sections.py -->

| Group | Section | What it does |
|---|---|---|
| **Protein** | **🔬 Sequence analysis** | Characterize VHHs you provide: physico-chemical properties, CDR/framework annotation, hallmarks and germline, liabilities, alignment, codon-optimised DNA and a PDF report. |
|  | **🧊 Structure** | Upload a PDB, mmCIF or AlphaFold model and replace the sequence-based guesses with measurements: which cysteines actually pair, which liabilities are buried and which are exposed, where the hydrophobic surface is, and, for a complex, which residues contact the antigen. |
|  | **📚 Sequence library** | The sequences you come back to, one click away: the public ones kept in the repository, and your own from a file you load, edit and download again. Send the proteins to the analysis section or the DNA to gene diagnosis, instead of pasting them again. |
| **DNA** | **🧪 Sanger sequencing** | Upload .ab1 chromatograms from a screening run, extract the VHH sequences from the expression cassette, spot frameshifts and duplicates, then send the clones you want to the analysis. |
|  | **🔧 Gene diagnosis** | Paste coding sequences that are already in a plasmid and ask what might be wrong with them: frame and stop checks, documented rare codons and their clusters, GC stretches, plus synonymous swaps that lower GC without touching the protein. |
|  | **🧬 In-silico cloning** | Load a backbone and an insert, and get the restriction enzymes that open the vector in the right place without cutting the VHH, plus which of your primers actually bind. |
| **Bench** | **⚗️ Batch & QC** | The arithmetic between a purified protein and the next experiment: concentration from a measured A280, whether the measured mass is the protein you think it is, how much labelling reagent to add, what each purification step recovered, and the batch sheet that leaves with the tube. |
|  | **🧂 Buffers & media** | What to weigh for a buffer or a medium at any volume and concentration factor, with the hydrate of each salt you actually have: the mass follows from its formula. The recipes are the laboratory's own, kept in the repository. |
|  | **🦠 Phage & panning** | The phage stocks that go into each round of panning: the titre from an A260, by two routes shown side by side, and how much stock a solid-phase or liquid-phase round takes. No sequence needed. |

<!-- /generated: sections -->

They share one sequence list: what you load in the analysis section is what the
cloning section offers as an insert, the sequencing section writes into it, the
analysis section's DNA tab hands its back-translations to gene diagnosis, and
the
library feeds both. Each section still works on its own.

Inside the analysis section the views are grouped into six **families** chosen
with
a selector, Triage, Identity, Developability, Comparison, Repertoire, Outputs.
Only the selected family executes, which is why it is a selector and not
another
row of tabs: Streamlit runs the body of every tab, nested or not. First display
after **Analyze** is 1.7 s on seven sequences, against roughly 12 s when all
fifteen views were tabs; the expensive family (Comparison, 9 s, alignment,
logo,
identity matrices) is now paid for only when asked for.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

`requirements.txt` is a lockfile resolved for Python 3.14, the version the
deployed app runs, and some of its pins do not install on older interpreters.
On another Python, install the ranges instead, `pip install -r
requirements.in`, and you get the newest versions that fit, as the deployed
app did before the lock.

Opens on <http://localhost:8501>. ANARCI is optional: without it the app
silently uses
the regex annotator and says so in the sidebar. To get ANARCI:

```bash
conda install -c conda-forge -c bioconda hmmer anarci
```

## Run with Docker

```bash
docker build -t vhh-toolbox .
docker run -p 8501:8501 vhh-toolbox
```

The image is built on micromamba specifically so HMMER can come from
bioconda: a pip-only image cannot install HMMER, and the app would be limited
to the regex annotator. Python comes from conda-forge, at the version
`requirements.txt` was compiled for, and everything else, ANARCI included,
from that lockfile.

The image takes the whole project (`COPY . /app/`) minus what `.dockerignore`
lists: `.git`, caches, zips, `results/`, `paper/` and `tools/`. A new module or
data file therefore needs no edit to either file; only something that should
*stay out* of the image does.

## Run in a dev container

`.devcontainer/devcontainer.json` is what GitHub Codespaces, or VS Code with
the Dev Containers extension, builds from. It starts from Microsoft's Python
3.14 image, the interpreter `requirements.txt` is compiled for, installs
HMMER from `packages.txt`, then the lockfile and the test requirements, and
starts the app on port 8501. The image is pinned to its major version
(`3-3.14-bookworm`), so a new major release of Microsoft's images does not
change it silently.

## Deploy on Streamlit Community Cloud

1. Push this folder to a **public GitHub repository on your personal account** (not an
   institutional one, the app should outlive any particular affiliation).
2. Go to <https://share.streamlit.io>, sign in with GitHub, click **New app**.
3. Point it at the repository, branch `main`, main file `app.py`.
4. In **Advanced settings**, select Python 3.14. That is the version
   `requirements.txt` is compiled for (its header says so), and the lockfile
   does not install on older ones. The Python version can only be chosen
   when the app is deployed; changing it later means deleting and
   redeploying the app, and recompiling the lock for the new version.
5. Deploy. The first build takes a few minutes: `packages.txt` installs HMMER through
   apt, then `requirements.txt` installs the Python packages including ANARCI.

Choose the subdomain under **Settings → General** (6–63 characters), which
gives you
`https://<your-subdomain>.streamlit.app`. Community Cloud does not support
pointing
your own domain name at the app, so if you want a stable URL of your own,
register a
domain and have it redirect to the Streamlit URL. The day you move to a
container host,
you repoint the redirect and no published link breaks.

Apps sleep after a period without traffic and wake on the next visit.

### If the build fails on ANARCI

Delete `packages.txt` and the `anarci` line from `requirements.in`, then
recompile `requirements.txt` (see *Updating dependencies*). The app
then runs
with the regex annotator only, which handles canonical VHH scaffolds perfectly
well, the sidebar says which method is active.

## Other hosts

| Target | ANARCI | Notes |
|---|---|---|
| **Google Cloud Run** | ✅ | Dockerfile included; scales to zero, generous monthly free allowance, custom domain supported |
| **Render, Scaleway** | ✅ | accept the Dockerfile directly |
| **Hugging Face Spaces** | ✅ | Docker SDK, but a Docker or Gradio Space now requires a PRO account for personal accounts |

Behind a reverse proxy (nginx, Apache), forward websockets, Streamlit will not
work
without them:

```nginx
location / {
    proxy_pass http://127.0.0.1:8501;
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
    proxy_set_header Host $host;
    proxy_read_timeout 86400;
}
```

## Data handling

Submitted sequences are **not stored in any persistent file, not logged, and
not sent
to a third party**. They are parsed into memory and analysed there. Precisely:

- **Result cache**: `@st.cache_data` keeps results in the server process memory, keyed
  on the input. It is shared by the whole process rather than per session, so
it is
  capped at 8 entries and 10 minutes (`CACHE_TTL` and `CACHE_MAX_ENTRIES` in
  `shell.py`, which every cached function in `compute.py` reads). Set
  `CACHE_TTL = 0` to disable caching entirely; the analysis takes well under a
  second, so the only cost is slightly less snappy sidebar controls.
- **Alignment**: `pymsaviz` reads from a path, so the aligned FASTA is written
  to a `TemporaryDirectory` that is deleted as soon as the figure is rendered.
- **Uploads**: handled in memory by Streamlit, never written to disk by this app.
- **Hosting**: on Community Cloud the app runs on third-party infrastructure. Traffic
  is encrypted and nothing is stored, but the data does transit a server you do
not own.
  The sidebar says so explicitly, and points users at running the app locally
for
  confidential sequences.

This is deliberate. People will paste unpublished VHH sequences into this form,
and
persisting them would mean holding other labs' intellectual property on your
server.
If you ever add persistence, replace the cache with a real database (SQLite at
minimum, a CSV cannot survive two concurrent writers) and say so explicitly in
the interface.

The upload size is capped at 2 MB and the app refuses more than 50 sequences or
any
sequence longer than 400 aa per run (`vhh_core.DEFAULTS`).

## Adding a descriptor

Register it once in `vhh_core.py`; it becomes a column of the results table, of
the
CSV export and of the "Custom features" panel automatically:

```python
@register_feature("my_descriptor", "What it measures and why it matters")
def _my_descriptor(seq, ctx):
    ann = ctx.get("annotation")        # CDR/FR boundaries, if annotation succeeded
    return round(some_calculation(seq), 3)
```

A feature that raises returns `NaN` rather than breaking the run.

## What each screen does

`MANUAL.md` documents every screen, what each number means and how far to
trust it, and the app renders it from the same file, so there is one text and
not two. It used to be summarised here as well, screen by screen, and the two
drifted: by this release the summary still described six tabs that had become
families, and said nothing of the four sections added since.

| For | Read |
|---|---|
| what a screen does, field by field | `MANUAL.md` §3 to §12, one per section |
| what the numbers do not mean | `MANUAL.md` §15, every caveat in one place |
| how your sequences are handled | `MANUAL.md` §13, and "Data handling" above |
| running it on your own machine | `MANUAL.md` §14, and "Run locally" above |
| when something looks wrong | `MANUAL.md` §16, troubleshooting |

Two things that live here rather than in the manual, because they are about
the repository rather than about using the app:

**Adding a Sanger vector.** The Sanger screen expects the PelB-VHH-tag
cassette of pMES3 or pMECS, which differ only in the tags expected after the
VHH. Add one in `sanger.VECTORS`, taking the layout from the actual plasmid
map: a wrong tag layout makes the app report missing tags on clones that are
fine.

**The synthetic generator's reference set.** `generator.py` samples frameworks
from an IMGT position-frequency matrix built from **PLAbDab-nano** (OPIG,
University of Oxford), 1958 unique VHH sequences, bundled as
`data/vhh_position_freq.csv` and `data/vhh_position_meta.json`.

> *PLAbDab-nano: a database of camelid and shark nanobodies from patents and
> literature*, [doi:10.1093/nar/gkae881](https://doi.org/10.1093/nar/gkae881).
> Its data-availability statement reads: "PLAbDab-nano is freely available to
> download and query at https://opig.stats.ox.ac.uk/webapps/plabdab-nano/."

Please cite it if you publish work based on generated sequences. What is
bundled here is an **aggregate frequency matrix**, not the sequences
themselves. To mimic your own library instead, supply a matrix of the same CSV
shape: one row per IMGT position, with `region`, `n_residues`, `occupancy` and
the twenty amino-acid columns.

## Tests

```bash
pip install -r requirements.txt -r requirements-dev.txt
pytest
```

Run it from the repository root. `pytest.ini` puts that directory on the import
path, so `pytest` and `python -m pytest` behave the same; without that line the
bare command finds only `tests/` and every test module fails to import.

The whole suite takes about 30 seconds (`pytest -m "not smoke"` runs the fast
ones in under ten). `requirements-dev.txt` is **not** installed on
Streamlit Community Cloud (the app itself needs `requirements.txt` alone) so
adding a test dependency can never break the deployment.

These are regression tests, not a correctness proof. Almost every one of them
is a bug this app actually had: it shipped, it was noticed by accident, it was
fixed. Each test names the bug it guards in its docstring, so a failure tells
you what the wrong behaviour looked like rather than just which assertion blew
up. Among them:

- hen lysozyme's eight ordinary cysteines being reported as eight free thiols,
  twice, once in `find_liabilities`, once again in the triage table
- repairs proposing to substitute the conserved FR2 and FR4 tryptophans
- a clean sequence producing a liabilities table with no columns, which killed
  the repair panel with `KeyError: 'score'`
- overlapping GC windows counted as 98 separate findings
- the structure mapping keyed by sequential index while every consumer read it
  as residue numbers

Two more were found by writing the tests themselves: `apply_repairs` applied
the
first of two conflicting substitutions instead of refusing the pair as its own
docstring promised, and `annotate_liabilities` accepted the full-length
construct where the chain's own sequence belongs, silently shifting every
position after a gap. Both are fixed and both now have a test.

One test earns its place for a different reason.
`test_module_versions_match_app`
compares every module's `__version__` against `REQUIRED_MODULES` in `app.py`.
That dict drives the staleness banner, and forgetting to update it is the
easiest way to ship a deploy that refuses to start while complaining about a
file you did update. It is caught here instead of in the browser.

Tests that pin an ANARCI-specific answer skip when ANARCI is absent rather than
fail, so the suite stays useful to someone who cannot install HMMER. The smoke
tests need Streamlit; `pytest -m "not smoke"` leaves them out and runs in about
seven seconds.

`tests/test_metadata.py` guards something different: the author's name, his
ORCID, the Zenodo concept DOI and the licence. None of it changes what the app
does, which is exactly why it needs a test, the DOI is added to `CITATION.cff`
in the browser after a release is published, so the next copy of the repository
overwrites it, and nothing visible breaks when it does. GitHub simply offers a
citation with no DOI in it, and you find out from someone else's reference
list.

The structure tests build their own PDB files in-process, CA atoms on a line,
with a deliberate numbering offset and a deliberate gap. The geometry is
meaningless and nothing asks anything of it; what is being tested is the
bookkeeping, and for that a synthetic chain beats a real one because the right
answer is known by construction.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md): how to report a bug, how to ask a
question, what is out of scope and why, and the two conventions the interface
code follows.

One point from it is worth repeating here. **GitHub issues are public and so is
`data/library.csv`**, do not paste unpublished or confidential sequences into
either. A defect can nearly always be reproduced with one of the 22 published
sequences in the library. To analyse confidential sequences at all, run the app
locally rather than using the public deployment.

## Licence

MIT, see `LICENSE`. Use it, change it, redistribute it, build on it; keep the
copyright notice.

One caveat, stated in `LICENSE` too: the MIT licence covers this repository's
code and cannot change the terms of its dependencies. All of them are
permissive
except **pyfamsa** (a wrapper around FAMSA, used for the multiple alignment),
which is **GPL-3.0-or-later**. Anyone redistributing the whole application with
the alignment in place has to honour GPL-3.0 for that component. It is imported
in exactly one function (`vhh_core.build_alignment`) so a derivative that does
not need FAMSA can drop the dependency and is then free of it.

Two bundled data files are not MIT either and keep their own terms:
`data/sabdab/` is **CC-BY 4.0** (SAbDab2, attribution required, see
`data/sabdab/LICENSE`) and `data/ecoli_codon_usage.csv` is **CC0-1.0**.
Neither affects the licence of the code.

## If this project is no longer maintained

It is written and maintained by one person. If the repository has gone quiet,
the deployed app has stopped answering, or you cannot reach the author, **take
it and run it**. You do not need permission, and this section exists to say so
explicitly: most people will not touch an abandoned project even when they are
entitled to, and a tool nobody dares pick up is a tool that is lost.

Everything you need is preserved independently of any account:

- The complete source of every release is archived on Zenodo at
  [10.5281/zenodo.23742269](https://doi.org/10.5281/zenodo.23742269). That
  archive does not depend on this GitHub repository, on the Streamlit
  deployment, or on the `vhh-tools-app.org` domain, any of which may lapse.
- The **MIT licence grants you the right** to use, modify, redistribute and
  redeploy it, commercially or not.
- `README.md` documents how to run it locally, with Docker, and how to deploy it
  on Streamlit Community Cloud. `MANUAL.md` documents what every number means and
  how far to trust it. `CONTRIBUTING.md` documents the conventions the code
  follows.

Two things are asked of you, and they are not the same kind of thing.

**You must keep the copyright notice.** This is the one condition the MIT
licence attaches: the notice in `LICENSE` has to travel with any copy or
substantial portion of the software. It is a legal obligation, not a courtesy.

**Please cite the DOI** if you publish work that uses this, or that uses a
version you derived from it. That one is an academic norm rather than a licence
term, so it is a request, not a condition. `CITATION.cff` holds the metadata and
GitHub turns it into a **Cite this repository** button.

If you fork it and carry it forward, say so in your own README and give your
fork its own DOI. A successor that is clearly a successor serves everyone better
than a silent copy.

## Maintaining this app

Written for whoever picks this up, including its author in two years.

**Where it runs.** Streamlit Community Cloud, deployed from this repository's
`main` branch, main file `app.py`. The custom domain redirects to the Streamlit
URL; Community Cloud cannot host a domain itself.

**To redeploy**, push to `main`, Community Cloud rebuilds automatically. The
container is rebuilt from the repository every time, so nothing the app writes
to
disk survives a restart. That is by design and is why the sequence library
lives
in `data/library.csv` in the repository rather than in a database.

**To check a deploy landed**, read the fingerprint in the footer:
`v1.0.0 (build 20261028.1aed6b5)`. The last part is a hash of what the
deployed app is made of: the `.py` files, `MANUAL.md`, everything under
`data/` (the library, the codon usage and position tables, the primers, the
plasmid maps and the amino acid structures), `.streamlit/config.toml`,
`requirements.txt` and
`packages.txt`. `data/` is walked, not listed, so a file added to it is
covered without editing anything. It changes when any of them changes,
including when one is moved, and the build date cannot
be mistaken for it. Same fingerprint after a push means the old container is
still being served. Files Community Cloud does not use, this README, the
changelog, the `Dockerfile`, are left out on purpose, so a push that touches
only those leaves the fingerprint where it was. Since `requirements.txt` is a
lockfile, a dependency cannot change version without the file, and so the
fingerprint, changing with it. Two things still can: the Python interpreter
Community Cloud runs, and the HMMER apt installs for `packages.txt`, which
names no version. The line under the footer, **Running on Python … ·
Streamlit … · Biopython … · ANARCI …**, is read from the running process for
that reason; ANARCI reads `absent` when the package is not installed.

### Updating dependencies

`requirements.in` holds the app's direct dependencies as ranges, and is the
only one of the two files you edit. `requirements.txt` is generated from it by
[pip-compile](https://pip-tools.readthedocs.io/) with every package pinned,
transitive ones included, and it is the file Community Cloud installs. It
never reads `requirements.in`. Never edit `requirements.txt` by hand: a test
fails if a line in it is not an exact pin, or if a package in
`requirements.in` is missing from it.

The lock has to be compiled **under the Python the deployed app runs**,
3.14, because pip-compile resolves for the interpreter it runs in. From the
repository root, with any Python 3.14:

```bash
python3.14 -m venv .venv-lock
.venv-lock/bin/pip install pip-tools
# Everything to the newest version requirements.in allows:
.venv-lock/bin/pip-compile --upgrade --strip-extras --no-emit-index-url \
    --output-file=requirements.txt requirements.in
# Or one package only, the rest left where they are:
.venv-lock/bin/pip-compile --upgrade-package streamlit --strip-extras \
    --no-emit-index-url --output-file=requirements.txt requirements.in
```

Without `--upgrade`, pip-compile keeps every pin already in `requirements.txt`
that still satisfies `requirements.in`, so adding a dependency does not
upgrade everything else behind your back. Then install the new lock into a
fresh Python 3.14 environment and run the whole suite there before pushing:

```bash
python3.14 -m venv .venv-check
.venv-check/bin/pip install -r requirements.txt -r requirements-dev.txt
.venv-check/bin/python -m pytest
```

The diff of `requirements.txt` is the list of what changed, worth reading
before the commit, a new major version especially. The push changes the footer fingerprint, since
`requirements.txt` is part of it, and the **Running on** line confirms what
the rebuilt container actually loaded.

If Community Cloud moves the app to another Python (only by deleting and
redeploying it), recompile under that version, and change the version named in
the deploy steps above, the `Dockerfile` and the manual: a test checks that
they all name the one in the lockfile header. The lock is also resolved for
Linux, which is what Community Cloud and the Docker image run; on macOS or
Windows a platform-specific dependency may be missing from it, and
`requirements.in` is the fallback.

**If a change does not appear**, or the app reports modules that are out of date
after you have updated them, Streamlit is serving cached imports. **Reboot
app**
in the Community Cloud dashboard clears it. This is the single most common
confusing failure.

**Before pushing**, run the tests:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pytest
```

**When changing a module**, bump its `__version__` *and* its entry in
`REQUIRED_MODULES` in `app.py`. Mismatches stop the app on purpose, and
`test_module_versions_match_app` catches them before a deploy does.

**Two rules for the interface**, so `app.py` does not grow back into the
3 911-line script it was: a screen draws and does not compute (cached work goes
in `compute.py`), and no screen reads another screen's globals (every screen is
a
`render()` function).

**The domain** renews annually and is the one piece that is neither in git nor
redeployable. If it lapses, every published link breaks, including any in a
paper. Keep auto-renew on, and make sure the card behind it is one that will
still exist next year.

## Releases

Versions follow `MAJOR.MINOR.PATCH`. The major number is not hypothetical here:
the analysis modules are importable from a notebook, so changing the signature
of
a `vhh_core` function breaks somebody's notebook even when the app looks
identical.

To cut a release:

1. Move the entries under `[Unreleased]` in `CHANGELOG.md` into a section
   for the new version (Added, Changed, Fixed), written for someone who uses
   the app.
2. Bump `RELEASE` (and `BUILD`) in `shell.py`, and `version` and
   `date-released` in `CITATION.cff`.
3. Run `pytest`: three tests check that the footer, `CITATION.cff` and the
   changelog all name the same version. They exist because a footer reading
   `version 20261028.2600` matched no git tag and no Zenodo record, so a bug
   report or a citation referred to something nobody could find.
4. Commit, push, then **Releases → Draft a new release** on GitHub with the tag
   `vX.Y.Z` and the changelog entry as the notes.
5. Zenodo archives it automatically and mints a version DOI.

Cut one when something meaningful has changed (a new section, a batch of fixes)
not on every push. Three empty releases in a week read worse than one real one
a quarter.

**Before citing the tool in a paper, cut a release and cite that release's
version DOI**, not the concept DOI. A reader then gets the exact code that
produced the result rather than whatever the latest version happens to be.

## Making the tool citable

`CITATION.cff` carries the citation metadata; GitHub turns it into a
**Cite this repository** button. The ORCID and the concept DOI are both in
it; the steps below are how the DOI was obtained, kept for whoever has to
redo them on another account.

For the DOI:

1. Sign in to [Zenodo](https://zenodo.org) **with GitHub** and authorise it.
2. Zenodo → **GitHub** settings, find this repository, flip its toggle to **On**.
   This installs a webhook. **Do this before releasing**, Zenodo only archives
   releases created after the toggle is on, and a release made first is simply
   not seen.
3. On GitHub: **Releases → Draft a new release**, tag `v1.0.0`, publish. That
   fires the webhook and Zenodo archives a snapshot of the repository at that
   tag.
4. Zenodo mints two DOIs: a **version DOI** for that release, and a
   version-independent **concept DOI** that always resolves to the latest. The
   concept DOI is the one that belongs in `CITATION.cff` and in a paper.
5. Put it in `CITATION.cff`, commit, and tag `v1.0.1`. The DOI cannot be in the
   first release, because it does not exist until that release is made.

The snapshot is of the repository as it stands at the tag, so `LICENSE` and
`CITATION.cff` must be committed before step 3 or the archived copy will not
have them.

That archive is the real answer to "what happens to this when I move on": a
GitHub URL and a Streamlit URL both depend on accounts someone has to keep
alive,
and a DOI does not.

## Caveats

Every caveat is kept in one place, `MANUAL.md` §15, because a second list of
them is a list that goes out of date while reading as though it were complete.
This one did: it carried four entries and never mentioned the structure
section, the buffers or a panning enrichment.

The three that decide whether this tool is the right one at all:

- Every sequence-derived value is **theoretical**. Confirm pI by IEF or cIEF
  and the extinction coefficient against a measured A280 before relying on
  either. The app says which of its numbers are measured, derived or
  theoretical, and that distinction is the point of it.
- Nothing here **predicts**. Not structure, not affinity, not immunogenicity,
  not labelling efficiency, not expression. The structure section measures a
  model you supply; the flags are triage heuristics with hand-set thresholds.
- Where two methods disagree, the app shows both and does not choose, which is
  a feature and not an omission: the two extinction coefficients, the two
  phage-titre routes, ANARCI against the regex anchors.

Built on Biopython (ProtParam), ANARCI, FAMSA (`pyfamsa`), `pymsaviz`,
logomaker, `dna_features_viewer`, reportlab and 3Dmol.js.

To tune the developability thresholds to your own experimental data, edit
`vhh_core.stability_flags`.
