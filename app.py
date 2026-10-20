"""
VHH Analyzer — Streamlit front-end for vhh_core.

Data handling: submitted sequences are never written to a persistent file, logged
or sent to a third party. They are processed in the server's memory. The result
cache is deliberately short-lived (see CACHE_TTL) and the alignment step writes a
temporary FASTA that is deleted as soon as the figure is rendered.

All analysis logic lives in vhh_core.py, which the Jupyter notebook imports too.

Run locally:   streamlit run app.py
"""

from __future__ import annotations

import io
import tempfile

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import sanger as sg
import vhh_core as vc

EXAMPLE = """>NbExample1
EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEG
RFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS
>NbExample2
QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKG
RFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADRLSITSRRYNYWGQGTQVTVSS
>cAbLys3
DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKG
RFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS
"""

# Version stamp, YYYYMMDD.HHMM. Bump this by hand when you push a change that
# users would notice — it is what someone will quote back to you in a bug report.
VERSION = "20261020.1845"
AUTHOR = "TyrosinDuck"
AUTHOR_URL = "https://github.com/TyrosinDuck"

st.set_page_config(page_title="VHH Analyzer", page_icon="🧬", layout="wide")

def _fit() -> dict:
    """Return the right 'fill the container' kwarg for the installed Streamlit.

    `use_container_width` is deprecated in favour of `width="stretch"` from 1.49
    onwards; this keeps the app working on both old and new versions.
    """
    try:
        from packaging.version import Version
        if Version(st.__version__.split("+")[0]) >= Version("1.49"):
            return {"width": "stretch"}
    except Exception:
        pass
    return {"use_container_width": True}


FIT = _fit()




@st.cache_resource
def anarci_status() -> bool:
    """Probe ANARCI once per server process, not once per request."""
    return vc.anarci_available()


# How long a result may stay in the in-memory cache, and how many results are
# kept at once. The cache makes moving a sidebar slider instant instead of
# re-running the whole analysis, but it is shared by the whole server process,
# so it is kept deliberately small and short-lived. Set CACHE_TTL = 0 to switch
# caching off entirely: the app then recomputes on every interaction, which
# costs well under a second for a handful of sequences.
CACHE_TTL = 600        # seconds
CACHE_MAX_ENTRIES = 8


def _analyze(seqs: dict[str, str], ph: float, method: str):
    return vc.analyze(seqs, ph=ph, method=method)


def _analyze_one(name: str, seq: str, ph: float, method: str):
    """Analyse a single sequence — the unit the progress bar advances on.

    Caching per sequence rather than per batch means adding one sequence to a
    list of ten only costs one analysis, not eleven.
    """
    frame, anns, errs = vc.analyze({name: seq}, ph=ph, method=method)
    return frame.iloc[0].to_dict(), anns.get(name), errs.get(name)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def hallmark_data(seq: str, annotation: dict | None = None):
    """Hallmark tetrad for one sequence (cached: ANARCI numbering is not free)."""
    return vc.hallmarks(seq, annotation)


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=64)
def germline_data(seq: str):
    """Closest germline per species; two ANARCI passes, worth caching."""
    return vc.germline(seq)


if CACHE_TTL:
    run_analysis = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)(_analyze)
    analyze_one = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL,
        max_entries=CACHE_MAX_ENTRIES * 8)(_analyze_one)
else:
    run_analysis = _analyze
    analyze_one = _analyze_one


@st.cache_data(show_spinner=False, ttl=CACHE_TTL or 60, max_entries=4)
def process_ab1(payload: tuple, vector: str, cutoff: float, min_len: int,
                min_phred: float):
    """Parse and extract a batch of .ab1 uploads (cached on the file contents)."""
    return sg.process_reads(list(payload), vector=vector, mott_cutoff=cutoff,
                            min_trimmed_len=min_len, min_mean_phred=min_phred)


def analyze_with_progress(seqs: dict[str, str], ph: float, method: str):
    """Run the analysis sequence by sequence, showing progress as it goes."""
    bar = st.progress(0.0, text="Preparing…")
    rows, anns, errs = [], {}, {}
    total = len(seqs)
    try:
        for i, (name, seq) in enumerate(seqs.items()):
            bar.progress(i / total, text=f"Analyzing {name} ({i + 1}/{total})…")
            row, ann, err = analyze_one(name, seq, ph, method)
            rows.append(row)
            if ann:
                anns[name] = ann
            if err:
                errs[name] = err
        bar.progress(1.0, text="Done")
    finally:
        bar.empty()
    return pd.DataFrame(rows), anns, errs


def footer() -> None:
    """Page footer — shown on every screen, including before any analysis."""
    st.divider()
    st.caption("Built on Biopython (ProtParam), ANARCI, FAMSA and pymsaviz. "
               "Theoretical predictions — validate experimentally. "
               "Sequences are not stored; for confidential work, run the app locally.")
    st.caption(f"VHH Analyzer · version {VERSION} · by [{AUTHOR}]({AUTHOR_URL})")


def fig_to_png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #

with st.sidebar:
    st.header("Settings")
    ph = st.slider("pH for net charge", 3.0, 10.0, vc.DEFAULTS["ph"], 0.1,
                   help="Net charge and the stability flags are evaluated at this pH.")
    window = st.slider("Profile window (residues)", 3, 21,
                       vc.DEFAULTS["profile_window"], 2,
                       help="Sliding window for the per-residue charge and hydropathy profiles.")

    has_anarci = anarci_status()
    options = ["auto", "regex"] + (["anarci"] if has_anarci else [])
    method = st.selectbox(
        "CDR annotation method", options, index=0,
        help=("'auto' uses ANARCI when available and falls back to the regex "
              "anchors otherwise."))
    st.caption(
        f"ANARCI: {'available' if has_anarci else 'not installed — regex anchors only'}")

    st.divider()
    with st.expander("🔒 How your sequences are handled"):
        st.markdown(f"""
Sequences are **never written to a persistent file, logged, or sent to a third
party**. They are processed in the server's memory and the results are kept in a
short-lived cache ({CACHE_TTL // 60} min, {CACHE_MAX_ENTRIES} entries max) that
only exists to make the sidebar controls responsive.

This app runs on Streamlit Community Cloud, i.e. on third-party infrastructure.
Traffic is encrypted, nothing is stored — but the data does pass through a server
that is not yours.

**For confidential or unpublished sequences**, run the app on your own machine
instead. The code is open source: clone the repository, `pip install -r
requirements.txt`, `streamlit run app.py`, and nothing ever leaves your computer.
""")
    st.caption(
        "All values are theoretical, derived from the sequence alone. "
        "Confirm pI and extinction coefficient experimentally before relying on them.")


# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #

st.title("🧬 VHH Analyzer")
st.markdown(
    "Physico-chemical characterization, CDR/framework annotation and alignment "
    "of single-domain antibodies. Paste one or more sequences, in FASTA or raw form.")

# A widget's state cannot be written after the widget exists, so anything that
# wants to fill the text area (the Sanger panel below) parks its text here and
# reruns; it is moved into place before the widget is created.
if "pending_seq_text" in st.session_state:
    st.session_state["seq_text"] = st.session_state.pop("pending_seq_text")

st.session_state.setdefault("seq_text", "")

col_in, col_act = st.columns([3, 1])
with col_act:
    st.write("")
    if st.button("Load example", **FIT):
        st.session_state["seq_text"] = EXAMPLE
        st.rerun()
    uploaded = st.file_uploader("or upload FASTA",
                                type=["fasta", "fa", "txt", "faa"])
    if uploaded is not None and uploaded.file_id != st.session_state.get("last_upload"):
        st.session_state["last_upload"] = uploaded.file_id
        st.session_state["seq_text"] = uploaded.read().decode("utf-8", errors="replace")
        st.rerun()
    st.caption(f"Up to {vc.DEFAULTS['max_sequences']} sequences, "
               f"{vc.DEFAULTS['max_length']} aa each.")

with col_in:
    # The widget owns st.session_state["seq_text"]; the buttons above write to it
    # before the widget is instantiated, which is the supported way to prefill it.
    text = st.text_area("Sequences", height=220, key="seq_text",
                        placeholder=">MyNanobody\nEVQLVESGGGLVQ...")
    # Clicking the button moves focus out of the text area, which is what commits
    # what was just typed — the same thing ⌘/Ctrl + Enter does, without the
    # keyboard shortcut nobody discovers on their own.
    st.button("Analyze", type="primary", **FIT)
    st.caption("Click Analyze, or press ⌘ + Enter (Ctrl + Enter on Windows/Linux).")

# --- Sanger .ab1 panel ------------------------------------------------------ #
with st.expander("🧬 Extract VHH sequences from Sanger .ab1 files", expanded=False):
    st.markdown(
        "Upload the chromatograms from a screening run. Each read is quality-trimmed, "
        "then the PelB–VHH–tag cassette is located in any strand and reading frame, "
        "and the VHH protein sequence is extracted. Pick the ones you want and send "
        "them to the analysis above.")

    ab1_files = st.file_uploader("Chromatograms (.ab1)", type=["ab1"],
                                 accept_multiple_files=True, key="ab1_uploader")

    opts = st.columns(4)
    with opts[0]:
        vector = st.selectbox("Vector", list(sg.VECTORS), index=0,
                              help="Determines which tags are expected after the VHH.")
    with opts[1]:
        mott_cutoff = st.slider("Trim cutoff", 0.01, 0.20,
                                sg.DEFAULTS["mott_cutoff"], 0.01,
                                help="Mott algorithm error-probability cutoff. "
                                     "Lower is stricter.")
    with opts[2]:
        min_len = st.number_input("Min trimmed length (nt)", 50, 800,
                                  sg.DEFAULTS["min_trimmed_len"], 10)
    with opts[3]:
        min_phred = st.number_input("Min mean Phred", 5, 50,
                                    sg.DEFAULTS["min_mean_phred"], 1)
    st.caption(f"{vector} layout: " + " → ".join(sg.VECTORS[vector]["layout"]))

    if ab1_files:
        payload = tuple((f.name, f.getvalue()) for f in ab1_files)
        bar = st.progress(0.0, text="Reading chromatograms…")
        try:
            sanger_df, reads = process_ab1(payload, vector, mott_cutoff,
                                           int(min_len), float(min_phred))
        finally:
            bar.empty()

        n_total = len(sanger_df)
        n_usable = int(sanger_df["usable"].sum()) if "usable" in sanger_df else 0
        n_unique = sanger_df[sanger_df["usable"]]["vhh_aa"].nunique() if n_usable else 0
        cols_m = st.columns(4)
        cols_m[0].metric("Reads", n_total)
        cols_m[1].metric("VHH extracted", n_usable)
        cols_m[2].metric("Unique sequences", n_unique)
        cols_m[3].metric("Failed", n_total - n_usable)

        show_cols = [c for c in ["read", "status", "qc_pass", "trimmed_len",
                                 "mean_phred_trim", "strand", "frame", "vhh_len",
                                 "tags_found", "tags_missing", "premature_stops",
                                 "n_copies", "duplicate_of", "warnings"]
                     if c in sanger_df.columns]
        st.dataframe(sanger_df[show_cols], hide_index=True, **FIT)

        if n_total - n_usable:
            st.caption("A failed read is not necessarily a failed clone: "
                       "'no FR4' usually means the read stopped before the end of the "
                       "VHH, and a frameshift means PelB and FR4 sit in different "
                       "reading frames, so an indel lies between them.")

        if reads:
            with st.expander("Quality traces"):
                fig_q = sg.plot_quality(reads)
                if fig_q:
                    st.pyplot(fig_q, **FIT)
                    st.caption("Blue is the window kept after trimming; the dashed "
                               "line is Phred 20.")
                plt.close("all")

        usable_df = sanger_df[sanger_df["usable"]] if n_usable else sanger_df.iloc[0:0]
        if n_usable:
            st.subheader("Send to analysis")
            unique_first = usable_df[usable_df["duplicate_of"] == ""]["read"].tolist()
            chosen = st.multiselect(
                "Sequences to analyse", usable_df["read"].tolist(),
                default=unique_first,
                help="Duplicates are left out by default — one entry per unique "
                     "protein sequence.")
            replace = st.checkbox("Replace what is currently in the text area",
                                  value=True)
            if st.button("Send selected to analysis", type="primary",
                         disabled=not chosen):
                picked = usable_df[usable_df["read"].isin(chosen)]
                fasta = sg.reads_to_fasta(picked, only_usable=False)
                current = st.session_state.get("seq_text", "").strip()
                st.session_state["pending_seq_text"] = (
                    fasta if replace or not current else current + "\n" + fasta)
                st.rerun()

        st.subheader("Download")
        dl = st.columns(3)
        with dl[0]:
            st.download_button("Sequencing results (CSV)",
                               sanger_df.drop(columns=["vhh_nt"], errors="ignore")
                               .to_csv(index=False),
                               file_name="sanger_results.csv", mime="text/csv")
        with dl[1]:
            st.download_button("VHH proteins (FASTA)",
                               sg.reads_to_fasta(sanger_df),
                               file_name="vhh_proteins.fasta", mime="text/plain",
                               disabled=not n_usable)
        with dl[2]:
            st.download_button("VHH DNA (FASTA)",
                               sg.reads_to_fasta(sanger_df, column="vhh_nt"),
                               file_name="vhh_dna.fasta", mime="text/plain",
                               disabled=not n_usable)
    else:
        st.info("No chromatogram uploaded yet.")

if not text.strip():
    st.info("Paste sequences above, upload a FASTA file, load the example, or extract "
            "them from .ab1 chromatograms in the panel above.")
    footer()
    st.stop()

# --------------------------------------------------------------------------- #
# Analysis
# --------------------------------------------------------------------------- #

try:
    sequences = vc.parse_input(text)
except ValueError as exc:
    st.error(f"Input error — {exc}")
    footer()
    st.stop()

df, annotations, errors = analyze_with_progress(sequences, ph, method)

if errors:
    for name, msg in errors.items():
        st.warning(f"**{name}** — annotation failed: {msg}. "
                   "Physico-chemical values are still reported below.")

tabs = st.tabs(["Properties", "Plots", "CDR / Framework", "Hallmarks & germline",
                "Liabilities", "Alignment", "DNA", "Download"])

# --- Properties ------------------------------------------------------------ #
with tabs[0]:
    main_cols = ["name", "length", "mw_kda", "pI", "ext_coeff_reduced", "A280_0.1pct",
                 "instability_index", "aliphatic_index", "gravy", "aromaticity",
                 "net_charge", "n_Cys"]
    st.dataframe(
        df[main_cols].style.format(precision=2).background_gradient(
            subset=["gravy", "instability_index"], cmap="RdYlGn_r"),
        **FIT, hide_index=True,
        column_config={
            "mw_kda": st.column_config.NumberColumn("MW (kDa)"),
            "ext_coeff_reduced": st.column_config.NumberColumn("ε280 (M⁻¹cm⁻¹)"),
            "A280_0.1pct": st.column_config.NumberColumn("A280 (1 mg/mL)"),
            "net_charge": st.column_config.NumberColumn(f"Charge pH {ph}"),
        })

    with st.expander("What these descriptors mean"):
        st.markdown("""
| Descriptor | Interpretation |
|---|---|
| **Theoretical pI** | pH of zero net charge. A pI close to the buffer pH predicts poor colloidal stability |
| **ε280 / A280 (1 mg/mL)** | extinction coefficient and the factor converting a measured A280 into mg/mL |
| **Instability index** | Guruprasad dipeptide score; **< 40 predicts a stable protein** |
| **Aliphatic index** | relative volume of A, V, I, L; higher correlates with thermostability |
| **GRAVY** | mean Kyte–Doolittle hydropathy; negative = hydrophilic. Above ≈ −0.2 flags aggregation risk |
| **Aromaticity** | fraction of F, W, Y; drives the A280 signal |
| **Secondary structure** | propensity estimate from composition, **not** a structure prediction |
""")

    st.subheader("Developability flags")
    for _, row in df.iterrows():
        flags = vc.stability_flags(row, ph=ph)
        with st.container(border=True):
            st.markdown(f"**{row['name']}**")
            if flags:
                for f in flags:
                    st.markdown(f"- ⚠️ {f}")
            else:
                st.markdown("- ✅ No flag raised")

    st.subheader("Custom features")
    feature_cols = ["name"] + [c for c in vc.FEATURE_REGISTRY if c in df.columns]
    st.dataframe(df[feature_cols], hide_index=True, **FIT)
    with st.expander("Feature definitions"):
        for name, entry in vc.FEATURE_REGISTRY.items():
            st.markdown(f"- **`{name}`** — {entry['doc']}")

# --- Plots ----------------------------------------------------------------- #
with tabs[1]:
    fig_ss = vc.plot_secondary_structure(df)
    st.pyplot(fig_ss, **FIT)

    fig_charge = vc.plot_charge(df, sequences, ph=ph)
    st.pyplot(fig_charge, **FIT)
    st.caption("Each titration curve crosses zero at the theoretical pI of that sequence.")

    fig_prof = vc.plot_profiles(sequences, annotations, window=window)
    st.pyplot(fig_prof, **FIT)
    st.caption("Sliding-window charge (blue) and Kyte–Doolittle hydropathy (orange). "
               "Shaded bands are the CDRs.")

    st.subheader("Predicted secondary structure")
    cf_window = st.slider("Chou–Fasman window", 3, 11,
                          vc.DEFAULTS["chou_fasman_window"], 2, key="cf_window")
    cf_figs = {}
    for name, seq in sequences.items():
        fig_cf = vc.plot_chou_fasman(seq, name, window=cf_window,
                                     annotation=annotations.get(name))
        cf_figs[name] = fig_cf
        st.pyplot(fig_cf, **FIT)
        prof = vc.chou_fasman_profile(seq, cf_window)
        st.caption(
            f"{name}: {prof['fractions']['E']:.0%} sheet, "
            f"{prof['fractions']['H']:.0%} helix, {prof['fractions']['C']:.0%} coil "
            "by this heuristic.")
    st.caption("The Chou–Fasman heuristic is roughly 50–60% accurate: it shows where "
               "the sequence tends towards each conformation, it does not predict the "
               "fold. A VHH is an immunoglobulin β-sandwich, so a sheet-dominated "
               "profile is the expected result. For real secondary structure, run DSSP "
               "on a model.")
    plt.close("all")

# --- CDR / Framework ------------------------------------------------------- #
with tabs[2]:
    if not annotations:
        st.warning("No sequence could be annotated.")
    else:
        rows = [{"name": n, "method": df.loc[df["name"] == n, "annotation_method"].iloc[0],
                 "CDR1": annotations[n]["CDR1"]["seq"],
                 "CDR2": annotations[n]["CDR2"]["seq"],
                 "CDR3": annotations[n]["CDR3"]["seq"],
                 "CDR1 len": len(annotations[n]["CDR1"]["seq"]),
                 "CDR2 len": len(annotations[n]["CDR2"]["seq"]),
                 "CDR3 len": len(annotations[n]["CDR3"]["seq"])}
                for n in annotations]
        st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT)

        st.subheader("Annotated sequences")
        st.markdown(
            '<span style="background:#e41a1c;color:#fff;padding:1px 6px;border-radius:3px">CDR1</span> '
            '<span style="background:#377eb8;color:#fff;padding:1px 6px;border-radius:3px">CDR2</span> '
            '<span style="background:#4daf4a;color:#fff;padding:1px 6px;border-radius:3px">CDR3</span>',
            unsafe_allow_html=True)
        st.markdown("".join(vc.render_annotated_html(n, annotations[n]) for n in annotations),
                    unsafe_allow_html=True)

        with st.expander("Framework regions"):
            st.dataframe(pd.DataFrame(
                [{"name": n, **{r: annotations[n][r]["seq"]
                                for r in ("FR1", "FR2", "FR3", "FR4")}}
                 for n in annotations]), hide_index=True, **FIT)

        with st.expander("How the two annotation methods work"):
            st.markdown("""
**ANARCI (IMGT numbering)** — the sequence is aligned to germline HMMs and each residue
receives an IMGT position; CDRs are the fixed windows **27–38**, **56–65** and **105–117**.
Numbering-based, so it tolerates unusual CDR lengths. Requires HMMER.

**Regex anchors** — dependency-free; locates the conserved landmarks of a VHH and cuts between them:

| Anchor | Pattern | Role |
|---|---|---|
| conserved Cys23 | `C` in residues 16–30 | ends FR1 |
| FR2 | `W[A-Z]RQ` (the `WFRQ`/`WYRQ` hallmark) | starts FR2 |
| FR3 | `[KR][FYW][TASV][IVLFM]S` (the `RFTIS` motif) | ends CDR2 |
| second Cys | `[YFHV][YFHVAC]C` (the `YYC` motif) | ends FR3 |
| FR4 | `WG[A-Z]G[TAS][A-Z]{1,2}VTVSS` | ends CDR3 |

Assumes a canonical scaffold. Where both methods run, a disagreement usually means a
non-canonical framework — which is itself worth knowing.
""")

# --- Hallmarks & germline --------------------------------------------------- #
with tabs[3]:
    st.markdown(
        "The four FR2 positions that distinguish a camelid VHH from a conventional "
        "VH, and the closest germline genes.")
    hm_rows, gl_rows = [], []
    for name, seq in sequences.items():
        try:
            hm = hallmark_data(seq, annotations.get(name))
            hm_rows.append({"name": name, "tetrad": hm["tetrad"],
                            "camelid score": f"{hm['camelid_score']}/4",
                            **{f"IMGT {p['IMGT']}": f"{p['residue']} ({p['verdict']})"
                               for p in hm["positions"]},
                            "source": hm["source"]})
        except Exception as exc:
            st.warning(f"**{name}** — hallmarks unavailable: {exc}")
        for sp, g in germline_data(seq).items():
            if g:
                gl_rows.append({"name": name, "species": sp, "V gene": g["v_gene"],
                                "V identity (%)": g["v_identity"],
                                "J gene": g["j_gene"],
                                "J identity (%)": g["j_identity"]})

    if hm_rows:
        st.dataframe(pd.DataFrame(hm_rows), hide_index=True, **FIT)
        st.caption("A score of 4/4 is the canonical camelid profile. IMGT 49 is shared "
                   "with VH3, so 3/4 with a G at 49 is common and still a genuine VHH. "
                   "Engineered or humanised constructs sit between the two profiles.")

    if gl_rows:
        st.subheader("Closest germline")
        st.dataframe(pd.DataFrame(gl_rows), hide_index=True, **FIT)
        st.caption("Identity to the closest human germline is a first-pass humanness "
                   "indicator, not an immunogenicity prediction — a sequence can be "
                   "germline-like and still contain a T-cell epitope.")
    elif not hm_rows:
        st.info("Germline assignment needs ANARCI, which is not available here.")

# --- Liabilities ------------------------------------------------------------ #
with tabs[4]:
    st.markdown(
        "Motifs that are known chemical or manufacturability liabilities, with their "
        "position and the region they fall in.")
    any_liab = False
    for name, seq in sequences.items():
        liab = vc.find_liabilities(seq, annotations.get(name))
        with st.container(border=True):
            st.markdown(f"**{name}** — {len(liab)} motif(s)")
            if liab.empty:
                st.markdown("No liability motif detected.")
            else:
                any_liab = True
                st.dataframe(
                    liab[["liability", "motif", "position", "region", "severity",
                          "note"]],
                    hide_index=True, **FIT,
                    column_config={"note": st.column_config.TextColumn("Note",
                                                                       width="large")})
    if any_liab:
        st.caption("A motif inside a CDR may affect binding; the same motif in a "
                   "framework is mainly a stability or heterogeneity concern. "
                   "Severity is a heuristic ranking, not an experimental measurement — "
                   "an NG in a loop that never sees stress may never deamidate.")

# --- Alignment ------------------------------------------------------------- #
with tabs[5]:
    if len(sequences) < 2:
        st.info("Submit at least two sequences to see an alignment.")
    else:
        with st.spinner("Aligning…"):
            # pymsaviz reads from a path, so the aligned FASTA is written to a
            # temporary directory that is removed as soon as the figure exists.
            with tempfile.TemporaryDirectory() as tmp:
                fig_msa, aln_path = vc.plot_msa(sequences, annotations, workdir=tmp)
                aln_text = aln_path.read_text()
        st.pyplot(fig_msa, **FIT)
        st.caption("FAMSA alignment rendered with pymsaviz. CDR bands refer to the "
                   "first sequence submitted.")
        st.download_button("Download alignment (FASTA)", aln_text,
                           file_name="alignment.fasta", mime="text/plain")

        st.subheader("Pairwise identity")
        fig_id = vc.plot_identity(vc.identity_matrix(sequences),
                                  "Full-length identity (%)")
        cols_id = st.columns(2)
        with cols_id[0]:
            st.pyplot(fig_id, **FIT)
        fig_id3 = None
        if annotations and any(a["CDR3"]["seq"] for a in annotations.values()):
            mat3 = vc.identity_matrix(sequences, annotations, region="CDR3")
            if len(mat3) >= 2:
                fig_id3 = vc.plot_identity(mat3, "CDR3 identity (%)")
                with cols_id[1]:
                    st.pyplot(fig_id3, **FIT)
        st.caption("Two clones from the same panning can share 90%+ full-length "
                   "identity while having unrelated CDR3s — the second matrix is the "
                   "one that tells you whether they are really the same binder.")
        plt.close("all")

# --- DNA -------------------------------------------------------------------- #
with tabs[6]:
    st.markdown(
        "Back-translation into DNA using *E. coli* preferred codons, with cloning "
        "sites and long homopolymer runs removed.")

    col_a, col_b = st.columns([2, 1])
    with col_a:
        avoid = st.multiselect(
            "Restriction sites to avoid", sorted(vc.DEFAULT_AVOID_SITES),
            default=sorted(vc.DEFAULT_AVOID_SITES),
            help="A codon creating one of these sites is swapped for a synonymous one.")
    with col_b:
        add_stop = st.checkbox("Append a stop codon (TAA)", value=True)
        gc_max = st.slider("Maximum GC fraction", 0.45, 0.70, 0.60, 0.01)

    sites = {k: vc.DEFAULT_AVOID_SITES[k] for k in avoid}
    fasta_parts = []
    for name, seq in sequences.items():
        res = vc.reverse_translate(seq, avoid_sites=sites, add_stop=add_stop,
                                   gc_max=gc_max)
        ok = vc.translate_dna(res["dna"]) == seq
        fasta_parts.append(f">{name}_optimized\n" +
                           "\n".join(res["dna"][i:i + 60]
                                     for i in range(0, len(res["dna"]), 60)))
        with st.container(border=True):
            st.markdown(f"**{name}** — {res['length_nt']} nt, "
                        f"GC {res['gc_percent']} %, "
                        f"{len(res['repairs'])} codon(s) swapped for sites, "
                        f"{res['gc_swaps']} for GC balance")
            if not ok:
                st.error("Round-trip check failed — do not order this sequence.")
            if res["remaining_issues"]:
                st.warning("Could not remove: " + ", ".join(res["remaining_issues"]))
            st.code(res["dna"], language="text")
            if res["repairs"]:
                with st.expander(f"{len(res['repairs'])} codon swap(s)"):
                    st.dataframe(pd.DataFrame(res["repairs"]), hide_index=True, **FIT)

    st.download_button("Download DNA (FASTA)", "\n".join(fasta_parts) + "\n",
                       file_name="vhh_optimized_dna.fasta", mime="text/plain",
                       type="primary")
    st.caption(
        "Every output is verified by translating it back and comparing to the input. "
        "This optimiser maximises codon frequency, balances GC and removes "
        "cloning-unfriendly motifs — it does not model mRNA secondary structure near "
        "the ribosome binding site, codon-pair bias or tRNA availability. Run the "
        "sequence through your synthesis vendor's checker before ordering.")

# --- Download -------------------------------------------------------------- #
with tabs[7]:
    st.markdown("Export the full results table, including every region sequence, "
                "custom feature and flag.")
    st.download_button("Results (CSV)", df.to_csv(index=False),
                       file_name="vhh_results.csv", mime="text/csv",
                       type="primary")
    st.download_button("Secondary structure (PNG)", fig_to_png(fig_ss),
                       file_name="secondary_structure.png", mime="image/png")
    st.download_button("Charge (PNG)", fig_to_png(fig_charge),
                       file_name="charge.png", mime="image/png")
    st.download_button("Profiles (PNG)", fig_to_png(fig_prof),
                       file_name="profiles.png", mime="image/png")

    st.divider()
    st.markdown("### Full PDF report")
    st.markdown(
        "Everything in one document: properties, developability flags, regions, "
        "hallmarks, germline, liabilities and all figures. Building it takes a few "
        "seconds, so it is generated on demand.")

    if st.button("Build PDF report", type="secondary"):
        figures = {"Secondary structure composition": fig_ss, "Charge": fig_charge,
                   "Charge and hydropathy profiles": fig_prof}
        for cf_name, cf_fig in cf_figs.items():
            figures[f"Predicted secondary structure — {cf_name}"] = cf_fig
        if len(sequences) >= 2:
            try:
                with tempfile.TemporaryDirectory() as tmp:
                    fig_msa_pdf, _ = vc.plot_msa(sequences, annotations, workdir=tmp)
                figures["Multiple sequence alignment"] = fig_msa_pdf
                figures["Pairwise identity"] = vc.plot_identity(
                    vc.identity_matrix(sequences), "Full-length identity (%)")
                mat3 = vc.identity_matrix(sequences, annotations, region="CDR3")
                if len(mat3) >= 2:
                    figures["CDR3 identity"] = vc.plot_identity(
                        mat3, "CDR3 identity (%)")
            except Exception as exc:
                st.warning(f"Alignment figures skipped: {exc}")

        with st.spinner("Building the report…"):
            with tempfile.TemporaryDirectory() as tmp:
                pdf_path = vc.build_pdf_report(
                    f"{tmp}/vhh_report.pdf", df, sequences, annotations,
                    figures=figures, ph=ph, version=VERSION, author=AUTHOR)
                st.session_state["pdf_bytes"] = pdf_path.read_bytes()
        plt.close("all")

    if st.session_state.get("pdf_bytes"):
        st.download_button("Download report (PDF)", st.session_state["pdf_bytes"],
                           file_name="vhh_report.pdf", mime="application/pdf",
                           type="primary")

    st.divider()
    st.dataframe(df.drop(columns=["sequence"]), hide_index=True, **FIT)

footer()
