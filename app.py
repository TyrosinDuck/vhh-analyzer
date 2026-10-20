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


if CACHE_TTL:
    run_analysis = st.cache_data(
        show_spinner=False, ttl=CACHE_TTL, max_entries=CACHE_MAX_ENTRIES)(_analyze)
else:
    run_analysis = _analyze


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

if not text.strip():
    st.info("Paste sequences above, upload a FASTA file, or load the example to start.")
    st.stop()

# --------------------------------------------------------------------------- #
# Analysis
# --------------------------------------------------------------------------- #

try:
    sequences = vc.parse_input(text)
except ValueError as exc:
    st.error(f"Input error — {exc}")
    st.stop()

with st.spinner(f"Analyzing {len(sequences)} sequence(s)…"):
    df, annotations, errors = run_analysis(sequences, ph, method)

if errors:
    for name, msg in errors.items():
        st.warning(f"**{name}** — annotation failed: {msg}. "
                   "Physico-chemical values are still reported below.")

tabs = st.tabs(["Properties", "Plots", "CDR / Framework", "Alignment", "Download"])

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

# --- Alignment ------------------------------------------------------------- #
with tabs[3]:
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
        plt.close("all")

# --- Download -------------------------------------------------------------- #
with tabs[4]:
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
    st.dataframe(df.drop(columns=["sequence"]), hide_index=True, **FIT)

st.divider()
st.caption("Built on Biopython (ProtParam), ANARCI, FAMSA and pymsaviz. "
           "Theoretical predictions — validate experimentally. "
           "Sequences are not stored; for confidential work, run the app locally.")
