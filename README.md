# VHH Analyzer

A web front-end for sequence-level characterization of VHH / single-domain antibodies:
physico-chemical descriptors, CDR/framework annotation (ANARCI or regex anchors),
plots and multiple sequence alignment.

```
vhh_webapp/
├── vhh_core.py        all analysis logic — no Streamlit, no notebook dependency
├── app.py             Streamlit front-end
├── requirements.txt   pip dependencies
├── packages.txt       apt dependencies (HMMER, the ANARCI backend)
├── Dockerfile         for a container host, if you ever outgrow Community Cloud
└── .streamlit/
    └── config.toml
```

`vhh_core.py` is the single source of truth. The notebook and the web app both import
it, so a fix to a descriptor propagates to both instead of drifting apart.

## Layout

The app opens on a home page with three sections, also reachable from the sidebar,
which carries a **Home** button on every screen:

| Section | What it does |
|---|---|
| **Sequence analysis** | characterization, CDR/framework annotation, hallmarks and germline, liabilities, alignment, codon-optimised DNA, PDF report |
| **Sanger sequencing** | extract VHHs from .ab1 chromatograms, then send the selected clones to the analysis section |
| **In-silico cloning** | backbone + insert → restriction enzymes that fit, primers that actually bind, plasmid map |

The three share one sequence list: what you load in the analysis section is what the
cloning section offers as an insert, and the sequencing section writes into it. Each
section still works on its own.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Opens on <http://localhost:8501>. ANARCI is optional: without it the app silently uses
the regex annotator and says so in the sidebar. To get ANARCI:

```bash
conda install -c conda-forge -c bioconda hmmer anarci
```

## Run with Docker

```bash
docker build -t vhh-analyzer .
docker run -p 8501:8501 vhh-analyzer
```

The image is built on micromamba specifically so HMMER and ANARCI can come from
bioconda — a pip-only image cannot install HMMER, and the app would be limited to the
regex annotator.

## Deploy on Streamlit Community Cloud

1. Push this folder to a **public GitHub repository on your personal account** (not an
   institutional one — the app should outlive any particular affiliation).
2. Go to <https://share.streamlit.io>, sign in with GitHub, click **New app**.
3. Point it at the repository, branch `main`, main file `app.py`.
4. In **Advanced settings**, select a Python version of 3.11 or later.
5. Deploy. The first build takes a few minutes: `packages.txt` installs HMMER through
   apt, then `requirements.txt` installs the Python packages including ANARCI.

Choose the subdomain under **Settings → General** (6–63 characters), which gives you
`https://<your-subdomain>.streamlit.app`. Community Cloud does not support pointing
your own domain name at the app, so if you want a stable URL of your own, register a
domain and have it redirect to the Streamlit URL. The day you move to a container host,
you repoint the redirect and no published link breaks.

Apps sleep after a period without traffic and wake on the next visit.

### If the build fails on ANARCI

Delete `packages.txt` and the `anarci` line from `requirements.txt`. The app then runs
with the regex annotator only, which handles canonical VHH scaffolds perfectly well —
the sidebar says which method is active.

## Other hosts

| Target | ANARCI | Notes |
|---|---|---|
| **Google Cloud Run** | ✅ | Dockerfile included; scales to zero, generous monthly free allowance, custom domain supported |
| **Render, Scaleway** | ✅ | accept the Dockerfile directly |
| **Hugging Face Spaces** | ✅ | Docker SDK, but a Docker or Gradio Space now requires a PRO account for personal accounts |

Behind a reverse proxy (nginx, Apache), forward websockets — Streamlit will not work
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

Submitted sequences are **not stored in any persistent file, not logged, and not sent
to a third party**. They are parsed into memory and analysed there. Precisely:

- **Result cache** — `@st.cache_data` keeps results in the server process memory, keyed
  on the input. It is shared by the whole process rather than per session, so it is
  capped at 8 entries and 10 minutes (`CACHE_TTL` / `CACHE_MAX_ENTRIES` in `app.py`).
  Set `CACHE_TTL = 0` to disable caching entirely; the analysis takes well under a
  second, so the only cost is slightly less snappy sidebar controls.
- **Alignment** — `pymsaviz` reads from a path, so the aligned FASTA is written to a
  `TemporaryDirectory` that is deleted as soon as the figure is rendered.
- **Uploads** — handled in memory by Streamlit, never written to disk by this app.
- **Hosting** — on Community Cloud the app runs on third-party infrastructure. Traffic
  is encrypted and nothing is stored, but the data does transit a server you do not own.
  The sidebar says so explicitly, and points users at running the app locally for
  confidential sequences.

This is deliberate. People will paste unpublished VHH sequences into this form, and
persisting them would mean holding other labs' intellectual property on your server.
If you ever add persistence, replace the cache with a real database (SQLite at minimum —
a CSV cannot survive two concurrent writers) and say so explicitly in the interface.

The upload size is capped at 2 MB and the app refuses more than 50 sequences or any
sequence longer than 400 aa per run (`vhh_core.DEFAULTS`).

## Adding a descriptor

Register it once in `vhh_core.py`; it becomes a column of the results table, of the
CSV export and of the "Custom features" panel automatically:

```python
@register_feature("my_descriptor", "What it measures and why it matters")
def _my_descriptor(seq, ctx):
    ann = ctx.get("annotation")        # CDR/FR boundaries, if annotation succeeded
    return round(some_calculation(seq), 3)
```

A feature that raises returns `NaN` rather than breaking the run.

## In-silico cloning

Five vectors ship with the app (pMECS, pMES3, pMES4, pHEN6, pHEN25) as SnapGene
files in `data/plasmids/`, and you can upload your own in SnapGene, GenBank or FASTA.
The app infers the cloning window between the PelB leader and the first tag, then
ranks restriction enzymes: a candidate cuts the backbone exactly once, inside that
window, and does not cut the insert. Enzymes that do cut the insert are still listed
with their position — in VHH cloning that is sometimes by design, BstEII having its
site in FR4.

Primers live in `data/primers.csv` (or upload your own CSV with `name,sequence`).
Only the primers with an exact binding site in the selected vector are offered; the
others are listed as not applicable, because a primer that does not bind is usually a
property of the vector rather than a mistake — pHEN6 has no gIII for a gIII primer to
sit in. Orientation is reported explicitly: a primer downstream of the insert reads
back into it only when used as its reverse complement.

## Sanger (.ab1) workflow

The panel under the input area takes the chromatograms from a screening run and turns
them into sequences you can analyse:

1. each read is quality-trimmed with Mott's algorithm (the one phred and Staden use);
2. the PelB–VHH–tag cassette is located in **any strand and any reading frame**, so
   forward and reverse reads are handled without being told which is which;
3. the VHH protein and DNA sequences are extracted, with frameshifts, premature stops,
   missing tags and truncated reads reported rather than silently returned;
4. identical protein sequences are grouped, so the same colony picked twice shows up as
   one sequence with a copy count;
5. you pick the ones you want and send them straight to the analysis above.

Supported vectors are pMES3 and pMECS — they differ only in which tags are expected
after the VHH. Add one in `sanger.VECTORS`, taking the layout from the actual plasmid
map: a wrong tag layout makes the app report missing tags on clones that are fine.

Everything is downloadable: the full per-read table as CSV, and the extracted VHHs as
protein or DNA FASTA.

## What the app reports

| Tab | Contents |
|---|---|
| Properties | physico-chemical descriptors, developability flags, custom features |
| Plots | secondary-structure composition, titration curves, per-residue charge/hydropathy profiles, and a Chou–Fasman predicted-secondary-structure plot with a per-residue state track |
| CDR / Framework | region table, coloured sequences, both annotation methods explained |
| Hallmarks & germline | the four FR2 positions (IMGT 42/49/50/52) scored camelid vs VH-like, and the closest alpaca and human germline V/J genes with percent identity |
| Liabilities | N-glycosylation sequons, deamidation, isomerisation, fragmentation, oxidation, extra cysteines and hydrophobic patches — each with its position and region |
| Alignment | FAMSA alignment with CDR bands, plus full-length and CDR3 identity matrices |
| DNA | back-translation into *E. coli* codon-optimised DNA, with restriction sites and homopolymer runs removed, GC balanced, and a round-trip verification |
| Download | CSV, figures as PNG, and a complete PDF report |

Three of these deserve a word of caution. The **Chou–Fasman plot** is a
composition-based heuristic, roughly 50–60% accurate: it shows where the sequence
tends towards each conformation, it does not predict the fold. The **DNA optimiser**
maximises codon frequency, balances GC and removes cloning-unfriendly motifs, but it
does not model mRNA secondary structure near the ribosome binding site, codon-pair
bias or tRNA availability — run its output through your synthesis vendor's checker
before ordering. Every back-translation is verified by translating it back and
comparing to the input. The **hallmark score** is indicative: IMGT 49
is shared with VH3, so 3/4 is common in genuine VHHs, and engineered constructs sit
between the two profiles. **Germline identity** measures similarity to the closest
germline gene — it is a first-pass humanness indicator, not an immunogenicity
prediction.

## Making the tool citable

Once it is stable, connect the GitHub repository to [Zenodo](https://zenodo.org) and cut
a release: Zenodo archives the code and mints a DOI, so the tool can be cited in a paper
and stays retrievable independently of GitHub, of any host, and of your affiliation.

## Caveats

- Every value is **theoretical**, derived from the sequence alone. Confirm pI by IEF or
  cIEF and the extinction coefficient against a measured A280 before relying on them.
- Secondary-structure fractions are composition-based propensities, **not** a structure
  prediction. Use DSSP on a model (AlphaFold3, ESMFold) for real secondary structure.
- The regex annotator assumes a canonical VHH scaffold. Where both methods run, a
  disagreement between ANARCI and the regex anchors usually means a non-canonical
  framework — worth investigating rather than ignoring.
- The developability flags are triage thresholds, not a verdict. Tune them to your own
  experimental data in `vhh_core.stability_flags`.

Built on Biopython (ProtParam), ANARCI, FAMSA (`pyfamsa`) and `pymsaviz`.
