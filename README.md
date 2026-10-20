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
