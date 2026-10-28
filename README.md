# VHH tool box

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

## Documentation

A full user manual lives in `MANUAL.md` and is rendered inside the app, reachable from
the home page and from the sidebar. It documents every screen, what each number means,
and how far to trust it. Edit the Markdown file and the app picks it up — there is one
source of truth.

### Bundled data

`data/library.csv` — the sequence library. Edit it locally, commit and push; the app serves it read-only. It is public, so nothing confidential goes in it.

`data/ecoli_codon_usage.csv` — *E. coli* codon usage (Kazusa Codon Usage Database, taxid 316407), redistributed from the `python-codon-tables` package under **CC0-1.0**, a public-domain dedication. Used for the CAI, the usage floor on codon substitutions, and the back-translation preference order.

## Layout

The app opens on a home page with five sections, also reachable from the sidebar,
which carries a **Home** button on every screen:

| Section | What it does |
|---|---|
| **Sequence analysis** | any protein is accepted and classified; generic characterization for all, antibody analyses for VHH and VH — CDR/framework annotation, hallmarks, humanness, liabilities and suggested repairs, disulfides, aggregation-prone regions, alignment and logo, clonal families, codon-optimised DNA, construct builder, PDF report |
| **Sanger sequencing** | extract VHHs from .ab1 chromatograms, then send the selected clones to the analysis section |
| **Gene diagnosis** | paste coding sequences that already exist, several at once: frame and stop checks, rare codons and their clusters, GC stretches, homopolymers, CAI against the measured *E. coli* usage table — plus synonymous swaps that lower GC without touching the protein, floored on real codon usage |
| **Sequence library** | a curated list kept in `data/library.csv`: filter it, then send proteins to the analysis or DNA to gene diagnosis. Read-only in the app and public in the repo, by design |
| **In-silico cloning** | backbone + insert → restriction enzymes that fit, primers that actually bind, plasmid map |

They share one sequence list: what you load in the analysis section is what the
cloning section offers as an insert, the sequencing section writes into it, the
analysis section's DNA tab hands its back-translations to gene diagnosis, and the
library feeds both. Each section still works on its own.

Inside the analysis section the views are grouped into six **families** chosen with
a selector — Triage, Identity, Developability, Comparison, Repertoire, Outputs.
Only the selected family executes, which is why it is a selector and not another
row of tabs: Streamlit runs the body of every tab, nested or not. First display
after **Analyze** is 1.7 s on seven sequences, against roughly 12 s when all
fifteen views were tabs; the expensive family (Comparison, 9 s — alignment, logo,
identity matrices) is now paid for only when asked for.

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
docker build -t vhh-toolbox .
docker run -p 8501:8501 vhh-toolbox
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

One screen, three steps: pick the recipient vector, paste the donor DNA, read the
plan.

Five vectors ship with the app (pMECS, pMES3, pMES4, pHEN6, pHEN25) in
`data/plasmids/`; you can upload your own in SnapGene, GenBank or FASTA. The donor is
either a whole plasmid carrying the VHH — the app then finds the domain, picks the
primer pair that brackets it and amplifies it in silico — or a fragment already
amplified with primers, which must carry the restriction sites those primers added. A
bare coding sequence cannot be cloned by restriction, and the app says so rather than
failing quietly.

The insert always lands between the end of the signal peptide and the first tag.
Every candidate enzyme pair is actually cut, joined and translated, then the VHH the
construct encodes is compared with the VHH in the donor. A pair is only recommended
when the protein survives unchanged *and* nothing annotated in the vector is deleted —
HindIII + BstEII translates perfectly on pMES3 while removing 120 bp that include the
PelB leader, so it is rejected.

Among the pairs that pass, the one disturbing the vector least is recommended, which
on these vectors is PstI + BstEII.

Three things worth knowing about the protein check:

- A fragment amplified with primers starts at the 5' restriction site, a few codons
  into FR1, so the first residues of the cloned VHH necessarily come from the
  recipient. The app reports this as *donor part preserved* and names the residues
  added, rather than calling it a mismatch.
- The vectors do not all start the domain the same way: **pMECS begins QVQLQ and
  pHEN25 begins EVQLQ**. Moving a VHH from pMECS into pHEN25 through PstI therefore
  changes residue 1, and that pair is rejected — NcoI + BstEII is proposed instead,
  because its site sits upstream of the difference.
- There is deliberately no "insert length must be a multiple of 3" test. An empty
  phage-display vector carries a stuffer, so its downstream tag is usually out of
  frame with the leader until a correctly sized VHH is cloned in. Translating the
  finished construct is the only check that means anything.

The verified construct downloads as GenBank, annotated with the insert and the
cassette parts, plus a protocol summary. Primers live in `data/primers.csv` (or upload
your own CSV with `name,sequence`); each is tried as written and as its reverse
complement, because a primer is often noted as the plus-strand sequence of its target
rather than as the oligo that primes back towards the insert.

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

## Synthetic sequence generator

`generator.py` builds VHH sequences whose frameworks are sampled from an IMGT
position-frequency matrix derived from real data, with random CDRs whose lengths come
from the same reference set. The bundled matrix is built from **PLAbDab-nano** (OPIG,
University of Oxford), 1958 unique VHH sequences — `data/vhh_position_freq.csv` and
`data/vhh_position_meta.json`.

> *PLAbDab-nano: a database of camelid and shark nanobodies from patents and
> literature*, [doi:10.1093/nar/gkae881](https://doi.org/10.1093/nar/gkae881).
> Its data-availability statement reads: "PLAbDab-nano is freely available to
> download and query at https://opig.stats.ox.ac.uk/webapps/plabdab-nano/."

Please cite it if you publish work based on the generated sequences. What is bundled
here is an **aggregate frequency matrix**, not the sequences themselves.

Supply your own matrix to mimic your own library instead: same CSV shape, one row per
IMGT position with `region`, `n_residues`, `occupancy` and the 20 amino-acid columns.

## What the app reports

| Family | Tab | Contents |
|---|---|---|
| — | Sequence type | each input classified from its sequence — VHH, conventional VH, light-chain V domain, a V domain inside a longer protein, or something else. VHH **and** conventional VH get the full heavy-chain analysis, labelled throughout; only the two camelid readings (FR2 hallmarks, extra CDR1–CDR3 disulfide) are flagged as such |
| 🎯 Triage | Triage | every flag the app raises, one row per sequence, sortable — no composite score, deliberately |
| 🔍 Identity | Properties | physico-chemical descriptors, developability flags, lysines and primary amines for conjugation planning, custom features |
| 🔍 Identity | Plots | secondary-structure composition, titration curves, per-residue charge/hydropathy profiles, and a Chou–Fasman predicted-secondary-structure plot with a per-residue state track |
| 🔍 Identity | CDR / Framework | region table, coloured sequences, both annotation methods explained |
| 🔍 Identity | Hallmarks & germline | the four FR2 positions (IMGT 42/49/50/52) scored camelid vs VH-like, and the closest alpaca and human germline V/J genes with percent identity |
| 🔍 Identity | Humanness | framework-only identity to the closest human germline V gene, computed position by position in IMGT numbering, with every divergent position listed and the hallmarks flagged |
| ⚠️ Developability | Liabilities | N-glycosylation sequons, deamidation, isomerisation, fragmentation, oxidation, extra cysteines and hydrophobic patches — each with its position and region |
| ⚠️ Developability | Disulfides | the canonical IMGT 23/104 pair, plausible additional CDR1–CDR3 or FR2–CDR3 pairs, and orphan thiols told apart |
| ⚠️ Developability | Suggested repairs | one conservative substitution per liability, refusing structurally conserved positions, flagging CDR ones, with the repaired sequence and its new score |
| ⚠️ Developability | Aggregation | per-residue aggregation-propensity profile and the weakest windows — a transparent composite, **not** CamSol or AGGRESCAN |
| 🔗 Comparison | Alignment | FAMSA alignment with CDR bands, plus full-length and CDR3 identity matrices |
| 🔗 Comparison | Sequence logo | logomaker logo of the alignment, frequency or bits, adjustable wrap, CDR bands, matrix downloadable as CSV |
| 🔗 Comparison | vs Reference | positions differing from a chosen reference, matched by IMGT number so an indel reports as one change rather than shifting everything downstream |
| 🌍 Repertoire | Repertoire | unusual framework residues against 1958 natural VHHs, and the submitted set's CDR lengths and composition against the same reference |
| 🔗 Comparison | Clonal families | CDR3-based grouping of panning output, with an adjustable identity threshold, complete or average linkage, a dendrogram, and one representative per family |
| 🧪 Outputs | DNA | back-translation into *E. coli* codon-optimised DNA, with restriction sites and homopolymer runs removed, GC balanced, and a round-trip verification |
| 🧪 Outputs | Constructs | monovalent, bivalent or biparatopic assemblies with linkers and tags, characterised and back-translated as the finished protein |
| 🧪 Outputs | Download | CSV, figures as PNG, and a complete PDF report |

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
