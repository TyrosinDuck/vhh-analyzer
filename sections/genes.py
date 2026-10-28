"""
sections.genes - Diagnose a coding sequence that already exists, and lower its GC.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import streamlit as st

import codons as cdn
import vhh_core as vc
from compute import diagnose_gene, gc_fix_for, gene_summary, png_gc_comparison, png_gene_diagnosis
from shell import EXAMPLE, FIT, footer, go_to

__version__ = "2026.10.28c"


def render_gene_detail(name, gene_dna, diag, expected, cluster_window, gc_window):
    """The per-sequence view: findings, figure, tables, and the GC suggestions."""
    orf = diag["orf"]
    gm = st.columns(5)
    gm[0].metric("Length", f"{orf['length_nt']} nt")
    gm[1].metric("Codons", orf["n_codons"])
    gm[2].metric("GC", f"{diag['gc_overall']:.1%}")
    gm[3].metric("Rare codons", diag["n_rare"],
                 help="Codons whose scarcity in E. coli is well established and "
                      "supplemented by the Rosetta and CodonPlus strains.")
    gm[4].metric("Rare clusters", len(diag["clusters"]))
    cai_now = cdn.cai(gene_dna).get("cai")
    if cai_now is not None:
        st.caption(
            f"**CAI {cai_now:.3f}** — Codon Adaptation Index against the measured "
            "*E. coli* usage table (Kazusa, CC0). 1.0 would mean the preferred "
            "codon everywhere. It describes codon choice; it does not predict "
            "expression level.")

    st.subheader("Findings")
    renderers = {"error": st.error, "warning": st.warning,
                 "ok": st.success, "info": st.info}
    for level, message in diag["findings"]:
        renderers.get(level, st.info)(message)

    st.image(png_gene_diagnosis(gene_dna, expected, cluster_window, gc_window),
             **FIT)

    if len(diag["clusters"]):
        st.subheader("Rare-codon clusters")
        st.dataframe(diag["clusters"], hide_index=True, **FIT)
        st.caption("Clustered and tandem rare codons are the ones associated "
                   "with ribosome stalling. Isolated ones are generally "
                   "tolerated, which is why they are counted but not listed.")

    if len(diag["gc_regions"]):
        st.subheader("Extreme GC stretches")
        st.dataframe(diag["gc_regions"], hide_index=True, **FIT)

    gc_rich = (diag["gc_regions"][diag["gc_regions"]["direction"] == "GC-rich"]
               if len(diag["gc_regions"]) else [])
    if len(gc_rich):
        st.subheader("Suggested codon changes")
        st.markdown(
            "A GC-rich stretch is one of the few sequence problems with a clean "
            "fix: swap codons for lower-GC synonyms. Every suggestion below has "
            "been checked — it is not a documented rare codon, it creates none "
            "of the restriction sites you select, it opens no homopolymer run, "
            "and the protein is verified unchanged.")
        fx = st.columns([1, 1, 1])
        with fx[2]:
            min_usage = st.slider(
                "Minimum E. coli usage of the substitute (%)", 0, 50,
                int(cdn.MIN_USAGE * 100), step=5, key=f"usage_{name}",
                help="The safeguard against trading a GC problem for a codon "
                     "one. A substitute must carry at least this share of its "
                     "amino acid in E. coli. Leucine's CTG carries 50% and TTA "
                     "13%, so a floor of 20% refuses that swap.") / 100
        with fx[0]:
            strategy = st.radio(
                "How far to go", ["minimal", "maximal"], key=f"strat_{name}",
                format_func=lambda x: ("Minimal — stop once the band is cleared"
                                       if x == "minimal"
                                       else "Maximal — swap everything possible"),
                help="Minimal applies the biggest reductions first and stops as "
                     "soon as no GC-rich stretch survives. You are paying to "
                     "synthesise this gene, and a construct differing in twenty "
                     "places is easier to reason about than one differing in "
                     "sixty.")
        with fx[1]:
            fix_avoid = st.multiselect(
                "Restriction sites to preserve", sorted(vc.DEFAULT_AVOID_SITES),
                default=sorted(vc.DEFAULT_AVOID_SITES), key=f"gcfix_{name}",
                help="No suggestion will create one of these.")

        fix = gc_fix_for(gene_dna, gc_window, strategy, tuple(sorted(fix_avoid)),
                         float(min_usage))
        fm = st.columns(5)
        fm[0].metric("Codons changed", fix["n_changes"])
        fm[1].metric("GC", f"{fix['gc_after']:.1%}",
                     delta=f"{(fix['gc_after'] - fix['gc_before']) * 100:+.1f} pts",
                     delta_color="off")
        fm[2].metric("GC-rich stretches left",
                     int((fix["remaining_regions"]["direction"] == "GC-rich").sum())
                     if len(fix["remaining_regions"]) else 0)
        fm[3].metric("CAI", f"{fix['cai_after']:.3f}"
                     if fix["cai_after"] is not None else "—",
                     delta=(f"{fix['cai_after'] - fix['cai_before']:+.3f}"
                            if None not in (fix["cai_after"], fix["cai_before"])
                            else None),
                     help="Codon Adaptation Index, computed from the measured "
                          "E. coli usage table. It falls as swaps move away "
                          "from preferred codons — this is the price of the GC "
                          "fix, in a number.")
        fm[4].metric("Protein",
                     "unchanged" if fix["protein_unchanged"] else "CHANGED")

        if not fix["protein_unchanged"]:
            st.error("The rewritten sequence does not translate to the same "
                     "protein. Do not use it, and please report it.")
        elif not fix["n_changes"]:
            st.info("No synonymous swap satisfies all the constraints here — "
                    "every lower-GC alternative is either a rare codon, creates "
                    "a site you asked to keep, or opens a homopolymer run.")
        else:
            comparison = png_gc_comparison(gene_dna, gc_window, strategy,
                                           tuple(sorted(fix_avoid)),
                                           float(min_usage))
            if comparison:
                st.image(comparison, **FIT)
                st.caption(
                    "Red is the sequence as you pasted it, blue after the "
                    "swaps; the grey band is the acceptable range and the "
                    "shaded red area is what used to sit above it. The ticks "
                    "along the bottom mark the codons that changed — if they "
                    "cluster in one place, that is the stretch that was the "
                    "problem.")
            st.dataframe(fix["regions"], hide_index=True, **FIT)
            st.dataframe(fix["suggestions"], hide_index=True, **FIT)
            st.caption(
                "`usage_from` and `usage_to` are the measured share each codon "
                "carries of its amino acid in *E. coli* — real frequencies from "
                "the Kazusa table, not an ordering. A swap from 0.50 to 0.13 is "
                "a genuine drop even though neither codon is rare; the CAI above "
                "is the same cost expressed over the whole gene.")
            st.markdown("**Rewritten sequence**")
            st.code(fix["dna"], language="text")
            st.download_button(
                "Rewritten DNA (FASTA)",
                f">{name}_gc_{strategy}\n{fix['dna']}\n",
                file_name=f"{name}_gc_{strategy}.fasta", mime="text/plain",
                key=f"dl_fix_{name}")
            st.caption("Lowering local GC removes one plausible obstacle to "
                       "expression. It does not promise the gene will then "
                       "express — the causes at the top of this page are all "
                       "still in play.")

    with st.expander("Every codon"):
        st.dataframe(diag["codons"], hide_index=True, **FIT)
        st.download_button("Codon table (CSV)", diag["codons"].to_csv(index=False),
                           file_name=f"{name}_codons.csv", mime="text/csv",
                           key=f"dl_codons_{name}")

    with st.expander("Translation"):
        st.code(orf["protein"] or "—", language="text")
        st.caption(f"{orf['protein_length']} residues.")

    st.caption(
        "**No CAI is reported.** A Codon Adaptation Index needs a full codon "
        "usage table for the host; this app does not ship one, and inventing the "
        "values would give a number that looks precise and means nothing.")


def render() -> None:
    """Draw the screen."""
    st.title("🔧 Gene diagnosis")
    st.markdown(
        "The reverse of designing a gene: paste coding sequences that already "
        "exist — in a plasmid, from a vendor, from a colleague — and ask what "
        "might be wrong with them. Several at a time.")
    st.warning(
        "**Codon usage is rarely the first suspect.** A VHH that expresses badly "
        "in *E. coli* is more often explained by periplasmic toxicity, "
        "inefficient signal-peptide cleavage, failure to form the intradomain "
        "disulfide, inclusion bodies, plasmid loss or induction conditions. A "
        "gene with flawless codon usage expresses badly for all of those "
        "reasons. This section tells you whether the sequence is *part* of the "
        "problem; it cannot tell you that it is the problem.", icon="⚠️")

    gcol_in, gcol_act = st.columns([3, 1])
    with gcol_act:
        st.write("")
        if st.button("Load example", key="gene_example", **FIT):
            st.session_state["pending_gene_dna"] = (
                ">example\n"
                + vc.reverse_translate(
                    list(vc.parse_input(EXAMPLE).values())[0])["dna"]
                + "\n")
            st.rerun()
        if st.button("🗑️ Clear", key="gene_clear", **FIT,
                     disabled=not st.session_state.get("gene_dna", "").strip()):
            st.session_state["gene_dna"] = ""
            st.rerun()
        gene_upload = st.file_uploader("or upload FASTA", type=["fasta", "fa", "txt"],
                                       key="gene_upload")
        if gene_upload is not None and gene_upload.file_id != st.session_state.get(
                "last_gene_upload"):
            st.session_state["last_gene_upload"] = gene_upload.file_id
            st.session_state["pending_gene_dna"] = gene_upload.read().decode(
                "utf-8", errors="replace")
            st.rerun()
        st.caption("FASTA, or one bare sequence, or several separated by blank "
                   "lines. Up to 20 sequences, 6000 nt each.")
        st.button("🔬 Design DNA from protein instead", key="genes_to_analysis",
                  **FIT, on_click=go_to, args=("analysis",),
                  help="The analysis section's DNA tab back-translates a protein "
                       "and can send the result here.")

    with gcol_in:
        gene_text = st.text_area("Coding DNA (5'→3')", height=200, key="gene_dna",
                                 placeholder=">clone_1\nCAGGTGCAGCTGCAGGAAAGCGGC…")
        st.button("Diagnose", type="primary", key="gene_run", **FIT)

    gset = st.columns(3)
    with gset[0]:
        g_cluster = st.slider("Cluster window (codons)", 2, 20,
                              int(cdn.CLUSTER_WINDOW), key="g_cluster")
    with gset[1]:
        g_gcwin = st.slider("GC window (nt)", 15, 120, int(cdn.GC_WINDOW), step=5,
                            key="g_gcwin")
    with gset[2]:
        g_expected = st.text_input(
            "Expected protein (optional)", key="g_expected",
            help="Compared against the translation of the sequence selected "
                 "below. Worth doing: a mismatch makes everything else beside "
                 "the point.")

    if not gene_text.strip():
        st.info("Paste one or more coding sequences to diagnose.")
    else:
        try:
            gene_seqs = cdn.parse_dna_input(gene_text)
        except ValueError as exc:
            st.error(f"Input error — {exc}")
            gene_seqs = {}

        if gene_seqs:
            st.subheader(f"Overview — {len(gene_seqs)} sequence(s)")
            overview = gene_summary(gene_seqs, int(g_gcwin), int(g_cluster))
            st.dataframe(overview, hide_index=True, **FIT)
            st.caption(
                "`in_frame` and `internal_stops` are the deterministic checks — "
                "wrong there means wrong. The rest is triage: `rare_clusters` "
                "matters more than `rare_codons`, and `gc_rich_stretches` is what "
                "the suggestions below act on.")
            st.download_button("Overview (CSV)", overview.to_csv(index=False),
                               file_name="gene_overview.csv", mime="text/csv")

            st.divider()
            picked = (list(gene_seqs)[0] if len(gene_seqs) == 1
                      else st.selectbox("Look at one in detail", list(gene_seqs),
                                        key="gene_pick"))
            gene_dna = gene_seqs[picked]
            try:
                diag = diagnose_gene(gene_dna, g_expected.strip(),
                                     int(g_cluster), int(g_gcwin))
            except ValueError as exc:
                st.error(f"Could not read that sequence — {exc}")
                diag = None

            if diag:
                render_gene_detail(picked, gene_dna, diag, g_expected.strip(),
                                   int(g_cluster), int(g_gcwin))

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
