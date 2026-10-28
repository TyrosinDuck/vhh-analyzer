"""
sections.analysis - the Sequence analysis screen.

Fifteen views of a set of sequences, grouped into six families chosen with a
selector. Only the selected family executes: that is the point of a selector
rather than nested tabs, since Streamlit runs the body of every tab whether or
not anyone opens it. First display on seven sequences is under two seconds,
against roughly twelve when all fifteen views were tabs.

This screen is one module where the others are one module per screen, because
its six families share a substantial prologue - the parsed input, the analysis
results, the annotations, the classifications - and splitting them into six
files would mean threading eighteen values six ways for nothing. They share one
`Run` object instead, and each family opens by unpacking exactly what it uses,
so the coupling is written down rather than implied.
"""

from __future__ import annotations

import io
import re
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import codons as cdn
import constructs as ctm
import generator as gn
import repertoire as rep
import triage as tg
import vhh_core as vc
from compute import (aggregation_for, analyze_with_progress, annotations,
    classify_cached, comparison_for, ctm, disulfides_for, dna_for,
    families_for, generate_cached, generator_matrix, germline_data, gn,
    hallmark_data, humanness_for, liabilities_for, pd, plt,
    png_aggregation, png_cdr_lengths, png_charge, png_chou_fasman,
    png_composition, png_construct, png_differences, png_family_dendrogram,
    png_identity, png_liability_map, png_liability_profile, png_logo,
    png_msa, png_profiles, png_rarity, png_secondary_structure,
    proteases_for, rarity_tables, rep, repairs_for, repertoire_stats, st,
    tempfile, tg, triage_for, vc)
from shell import (CACHE_TTL, EXAMPLE, FIT, anarci_status, current_sequences,
                   fig_to_png, footer, go_to)

__version__ = "2026.10.28a"


@dataclass
class Run:
    """Everything the six families share, assembled once per render.

    A dataclass rather than eighteen parameters: the families need overlapping
    but different subsets of it, and a positional signature would have to be
    edited in seven places the first time one of them needs one more value.

    The last two entries are functions, not data. They are closures built in
    `render` over the classifications, and they belong to the run exactly as
    much as the classifications do - a family that wants to know whether to draw
    an antibody-only panel should ask the run, not re-derive the answer.
    """
    sequences: Any
    df: Any
    annotations: Any
    errors: Any
    kinds: Any
    heavy_names: Any
    vhh_names: Any
    other_names: Any
    heavy_sequences: Any
    heavy_annotations: Any
    _vh_only: Any
    ph: Any
    window: Any
    method: Any
    tabs: Any
    chosen: Any
    heavy_only_notice: Callable
    type_of: Callable


def _triage(r: Run) -> None:
    """The triage family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    sequences = r.sequences
    annotations = r.annotations
    kinds = r.kinds
    ph = r.ph
    tabs = r.tabs

    # --- Triage ----------------------------------------------------------------- #
    with tabs["Triage"]:
        st.markdown(
            "Every flag the app can raise, one row per sequence. This is the tab for "
            "the moment you have to **decide**, which clones to carry forward, which "
            "ten of sixteen to order, rather than to study one candidate.")
        st.info(
            "**There is no overall score, on purpose.** Weighting a deamidation site "
            "against a hydrophobic patch against three divergent framework positions "
            "means choosing what matters, and that depends on what the VHH is for. A "
            "single number would hide that choice behind a decimal point. Sort by the "
            "column you care about today.", icon="ℹ️")

        tri_cols = st.columns([2, 1])
        with tri_cols[0]:
            shown_groups = st.multiselect(
                "Column groups", list(tg.COLUMN_GROUPS),
                default=[g for g in tg.COLUMN_GROUPS if g != "Physico-chemical"])
        with tri_cols[1]:
            flagged_only = st.checkbox(
                "Only sequences with a flag", value=False,
                help="Keeps the rows carrying a high-risk liability, an orphan "
                     "cysteine, an unusual framework residue or an "
                     "aggregation-prone region.")

        with st.spinner("Collecting…"):
            table = triage_for(sequences, annotations, kinds, ph, include_added=True)

        view = table
        if flagged_only and len(view):
            mask = view[tg.FLAG_COLUMNS].fillna(False).any(axis=1)
            view = view[mask]

        keep = [c for g in shown_groups for c in tg.COLUMN_GROUPS[g]]
        keep = [c for c in dict.fromkeys(["name"] + keep) if c in view.columns]
        st.dataframe(view[keep], hide_index=True, **FIT)
        st.caption(
            f"{len(view)} of {len(table)} sequence(s) shown. Empty cells are "
            "analyses that do not apply, a sequence that is not an antibody domain "
            "has no hallmark score, and its cysteines are not read against the VHH "
            "canonical pair. Click a column header to sort.")

        if len(table):
            flags = table[["name"] + tg.FLAG_COLUMNS].copy()
            n_flagged = int(flags[tg.FLAG_COLUMNS].fillna(False).any(axis=1).sum())
            fm = st.columns(4)
            fm[0].metric("Sequences", len(table))
            fm[1].metric("With a flag", n_flagged)
            fm[2].metric("High-risk liability",
                         int(flags["has_high_risk"].fillna(False).sum()))
            fm[3].metric("Orphan cysteine",
                         int(flags["has_orphan_cys"].fillna(False).sum()))
            with st.expander("The flag columns on their own"):
                st.dataframe(flags, hide_index=True, **FIT)

        st.download_button("Triage table (CSV)", table.to_csv(index=False),
                           file_name="triage.csv", mime="text/csv")


def _identity(r: Run) -> None:
    """The identity family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    sequences = r.sequences
    df = r.df
    annotations = r.annotations
    kinds = r.kinds
    heavy_sequences = r.heavy_sequences
    heavy_annotations = r.heavy_annotations
    ph = r.ph
    window = r.window
    tabs = r.tabs
    heavy_only_notice = r.heavy_only_notice
    type_of = r.type_of

    # --- Properties ------------------------------------------------------------ #
    with tabs["Properties"]:
        main_cols = ["name", "type", "length", "mw_kda", "pI", "ext_coeff_reduced",
                     "A280_0.1pct", "instability_index", "aliphatic_index", "gravy",
                     "aromaticity", "net_charge", "n_Cys", "n_Lys",
                     "n_primary_amines"]
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
                st.markdown(f"- **`{name}`**, {entry['doc']}")

        with st.expander("Lysines and primary amines, for conjugation"):
            st.markdown(
                "`n_Lys` counts lysine side chains; `n_primary_amines` adds the free "
                "N-terminal α-amine, which an NHS ester labels just as readily. "
                "Together they are the sites a random amine-coupling chemistry, "
                "NHS-biotin, NHS-fluorophore, amine coupling onto a chip, will react "
                "with.")
            lys_rows = []
            for name, seq in sequences.items():
                table = vc.lysine_sites(seq, annotations.get(name))
                lys_rows.append(table.assign(sequence=name))
            lys_all = pd.concat(lys_rows, ignore_index=True) if lys_rows else pd.DataFrame()
            if len(lys_all):
                st.dataframe(
                    lys_all[["sequence", "site", "position", "region", "in_CDR", "note"]],
                    hide_index=True, **FIT,
                    column_config={"note": st.column_config.TextColumn("Note",
                                                                       width="large")})
                in_cdr = lys_all[lys_all["in_CDR"]]
                if len(in_cdr):
                    st.warning(
                        f"{len(in_cdr)} amine site(s) sit inside a CDR: "
                        + ", ".join(f"{r['sequence']} {r['site']}"
                                    for _, r in in_cdr.iterrows())
                        + ". Random amine coupling will hit these along with the "
                          "rest, which is the usual explanation for a VHH that loses "
                          "binding after labelling. Site-specific chemistry, a "
                          "C-terminal cysteine, sortase, or an Avi tag, avoids the "
                          "problem entirely.", icon="⚠️")
                st.download_button("Amine sites (CSV)", lys_all.to_csv(index=False),
                                   file_name="amine_sites.csv", mime="text/csv")
            st.caption(
                "Counted from sequence, so this is how many amines **exist**, not "
                "how many are solvent-accessible, that needs a structure. Treat it "
                "as the upper bound on labelling sites.")

    # --- Plots ----------------------------------------------------------------- #
    with tabs["Plots"]:
        st.image(png_secondary_structure(df), **FIT)

        st.image(png_charge(df, sequences, ph), **FIT)
        st.caption("Each titration curve crosses zero at the theoretical pI of that "
                   "sequence.")

        st.image(png_profiles(sequences, annotations, window), **FIT)
        st.caption("Sliding-window charge (blue) and Kyte–Doolittle hydropathy (orange). "
                   "Shaded bands are the CDRs.")

        st.subheader("Predicted secondary structure")
        cf_window = st.slider("Chou–Fasman window", 3, 11,
                              vc.DEFAULTS["chou_fasman_window"], 2, key="cf_window")
        # With many sequences, drawing every panel on every rerun is the single
        # slowest thing the app does, so only the selected one is shown by default.
        cf_names = list(sequences)
        if len(cf_names) > 1:
            show_all_cf = st.checkbox(f"Show all {len(cf_names)} sequences", value=False,
                                      key="cf_all")
            cf_shown = cf_names if show_all_cf else [
                st.selectbox("Sequence", cf_names, key="cf_pick")]
        else:
            cf_shown = cf_names

        for name in cf_shown:
            seq = sequences[name]
            st.image(png_chou_fasman(name, seq, cf_window, annotations.get(name)), **FIT)
            prof = vc.chou_fasman_profile(seq, cf_window)
            st.caption(f"{name}: {prof['fractions']['E']:.0%} sheet, "
                       f"{prof['fractions']['H']:.0%} helix, "
                       f"{prof['fractions']['C']:.0%} coil by this heuristic.")
        st.caption("The Chou–Fasman heuristic is roughly 50–60% accurate: it shows where "
                   "the sequence tends towards each conformation, it does not predict the "
                   "fold. A VHH is an immunoglobulin β-sandwich, so a sheet-dominated "
                   "profile is the expected result. For real secondary structure, run "
                   "DSSP on a model.")
        plt.close("all")

    # --- CDR / Framework ------------------------------------------------------- #
    with tabs["CDR / Framework"]:
        _cdr_ok = heavy_only_notice("CDR and framework annotation")
        if _cdr_ok and not heavy_annotations:
            st.warning("No sequence could be annotated.")
        elif _cdr_ok:
            rows = [{"name": n, "type": type_of(n),
                     "method": df.loc[df["name"] == n, "annotation_method"].iloc[0],
                     "CDR1": heavy_annotations[n]["CDR1"]["seq"],
                     "CDR2": heavy_annotations[n]["CDR2"]["seq"],
                     "CDR3": heavy_annotations[n]["CDR3"]["seq"],
                     "CDR1 len": len(heavy_annotations[n]["CDR1"]["seq"]),
                     "CDR2 len": len(heavy_annotations[n]["CDR2"]["seq"]),
                     "CDR3 len": len(heavy_annotations[n]["CDR3"]["seq"])}
                    for n in heavy_annotations]
            st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT)

            st.subheader("Annotated sequences")
            st.markdown(
                '<span style="background:#e41a1c;color:#fff;padding:1px 6px;'
                'border-radius:3px">CDR1</span> '
                '<span style="background:#377eb8;color:#fff;padding:1px 6px;'
                'border-radius:3px">CDR2</span> '
                '<span style="background:#4daf4a;color:#fff;padding:1px 6px;'
                'border-radius:3px">CDR3</span>', unsafe_allow_html=True)
            st.markdown("".join(vc.render_annotated_html(n, heavy_annotations[n])
                                for n in heavy_annotations), unsafe_allow_html=True)

            with st.expander("Framework regions"):
                st.dataframe(pd.DataFrame(
                    [{"name": n, **{r: heavy_annotations[n][r]["seq"]
                                    for r in ("FR1", "FR2", "FR3", "FR4")}}
                     for n in heavy_annotations]), hide_index=True, **FIT)

            with st.expander("How the two annotation methods work"):
                st.markdown("""
    **ANARCI (IMGT numbering)**: the sequence is aligned to germline HMMs and each residue
    receives an IMGT position; CDRs are the fixed windows **27–38**, **56–65** and **105–117**.
    Numbering-based, so it tolerates unusual CDR lengths. Requires HMMER.

    **Regex anchors**: dependency-free; locates the conserved landmarks of a VHH and cuts between them:

    | Anchor | Pattern | Role |
    |---|---|---|
    | conserved Cys23 | `C` in residues 16–30 | ends FR1 |
    | FR2 | `W[A-Z]RQ` (the `WFRQ`/`WYRQ` hallmark) | starts FR2 |
    | FR3 | `[KR][FYW][TASV][IVLFM]S` (the `RFTIS` motif) | ends CDR2 |
    | second Cys | `[YFHV][YFHVAC]C` (the `YYC` motif) | ends FR3 |
    | FR4 | `WG[A-Z]G[TAS][A-Z]{1,2}VTVSS` | ends CDR3 |

    Assumes a canonical scaffold. Where both methods run, a disagreement usually means a
    non-canonical framework, which is itself worth knowing.
    """)

    # --- Hallmarks & germline --------------------------------------------------- #
    with tabs["Hallmarks & germline"]:
        st.markdown("The four FR2 positions that distinguish a camelid VHH from a "
                    "conventional VH, and the closest germline genes.")
        _has_vhh = heavy_only_notice("Hallmark and germline analysis")
        hm_rows, gl_rows = [], []
        for name, seq in heavy_sequences.items():
            try:
                hm = hallmark_data(seq, annotations.get(name))
                hm_rows.append({"name": name, "type": type_of(name),
                                "tetrad": hm["tetrad"],
                                "camelid score": f"{hm['camelid_score']}/4",
                                **{f"IMGT {p['IMGT']}": f"{p['residue']} ({p['verdict']})"
                                   for p in hm["positions"]},
                                "source": hm["source"]})
            except Exception as exc:
                st.warning(f"**{name}**, hallmarks unavailable: {exc}")
            for sp, g in germline_data(seq).items():
                if g:
                    gl_rows.append({"name": name, "species": sp, "V gene": g["v_gene"],
                                    "V identity (%)": g["v_identity"],
                                    "J gene": g["j_gene"],
                                    "J identity (%)": g["j_identity"]})

        if hm_rows:
            st.dataframe(pd.DataFrame(hm_rows), hide_index=True, **FIT)
            st.caption("A score of 4/4 is the canonical camelid profile. IMGT 49 is "
                       "shared with VH3, so 3/4 with a G at 49 is common and still a "
                       "genuine VHH. Engineered or humanised constructs sit between the "
                       "two profiles. **This is the one reading that is camelid-"
                       "specific**: on a conventional VH a score of 0/4 is simply what "
                       "a VH looks like, not a defect, the row is shown so you can "
                       "see the contrast, and the `type` column says which is which.")

        if gl_rows:
            st.subheader("Closest germline")
            st.dataframe(pd.DataFrame(gl_rows), hide_index=True, **FIT)
            st.caption("Identity to the closest human germline is a first-pass humanness "
                       "indicator, not an immunogenicity prediction, a sequence can be "
                       "germline-like and still contain a T-cell epitope.")
        elif not hm_rows:
            st.info("Germline assignment needs ANARCI, which is not available here.")

        st.divider()
        st.subheader("Humanness, position by position")
        st.markdown(
            "Percent identity to the closest human germline V gene, computed over "
            "**framework positions only**, CDR divergence from a human germline is "
            "expected and is not something you would ever try to fix. The divergent "
            "positions are listed because they are the candidate humanising "
            "mutations.")
        hcol = st.columns([1, 2])
        with hcol[0]:
            # Deliberately not called `family`: that name holds the active
            # analysis family (the section selector) in this file's scope, and
            # shadowing it here would silently disable every later family block.
            gl_family = st.selectbox(
                "Germline family", ["IGHV3", "all human heavy V"], index=0,
                help="IGHV3 is the family VHHs belong to, so it is the fair "
                     "comparison. Searching all families rarely changes the answer.")
        fam_arg = None if gl_family == "all human heavy V" else gl_family

        for name, seq in heavy_sequences.items():
            try:
                h = humanness_for(seq, fam_arg)
            except Exception as exc:
                st.warning(f"**{name}**, humanness unavailable: {exc}")
                continue
            if not h["gene"]:
                st.warning(f"**{name}**, no germline could be compared.")
                continue
            with st.expander(f"**{name}** ({type_of(name)}), "
                             f"{h['identity']}% to {h['gene']} "
                             f"({h['n_divergent']} divergent framework positions, "
                             f"{h['n_hallmark']} of them hallmarks)",
                             expanded=len(sequences) == 1):
                m = st.columns(3)
                m[0].metric("Framework identity", f"{h['identity']}%")
                m[1].metric("Closest human gene", h["gene"])
                m[2].metric("Divergent FR positions", h["n_divergent"])
                if len(h["by_region"]):
                    st.dataframe(h["by_region"], hide_index=True, **FIT)
                if len(h["divergences"]):
                    st.dataframe(h["divergences"], hide_index=True, **FIT)
                    st.download_button(
                        f"Divergences, {name} (CSV)",
                        h["divergences"].to_csv(index=False),
                        file_name=f"humanness_{name}.csv", mime="text/csv",
                        key=f"dl_human_{name}")
                if not kinds[name]["is_vhh"]:
                    st.caption(
                        "This is a conventional VH, so a high framework identity is "
                        "expected rather than informative, and the `hallmark` rows "
                        "below are not a humanisation decision here, since a VH "
                        "already carries the VH residues.")
                st.caption(
                    "The rows flagged `hallmark` are the four FR2 positions that make "
                    "a domain a VHH rather than a VH. They are the real decision: "
                    "reverting them raises the humanness number and simultaneously "
                    "removes the camelid surface that keeps a VHH soluble without a "
                    "light chain. The other rows are the cheap ones.")

        st.info(
            "This is a first-pass humanness indicator and nothing more. "
            "Immunogenicity depends on T-cell epitopes, dose, route and formulation, "
            "none of which is in the sequence, a germline-like sequence can still "
            "carry an epitope, and a divergent one can be silent.", icon="⚠️")


def _risk(r: Run) -> None:
    """The risk family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    sequences = r.sequences
    df = r.df
    annotations = r.annotations
    kinds = r.kinds
    heavy_sequences = r.heavy_sequences
    tabs = r.tabs
    heavy_only_notice = r.heavy_only_notice

    # --- Liabilities ------------------------------------------------------------ #
    with tabs["Liabilities"]:
        st.markdown(
            "Chemical and manufacturability liabilities, scored and placed along the "
            "sequence.")

        include_added = st.checkbox(
            "Include the motifs added on top of the in-house scale",
            value=True,
            help="Fragmentation (DP), RGD and hydrophobic patches are not part of the "
                 "original scoring scale; untick to reproduce its totals exactly.")

        liab_by_seq, prot_by_seq, summaries = {}, {}, []
        for name, seq in sequences.items():
            liab_by_seq[name] = liabilities_for(seq, annotations.get(name), include_added)
            prot_by_seq[name] = proteases_for(seq, annotations.get(name))
            summaries.append(vc.liability_summary(name, liab_by_seq[name],
                                                  prot_by_seq[name]))

        st.dataframe(pd.DataFrame(summaries), hide_index=True, **FIT,
                     column_config={
                         "total_score": st.column_config.NumberColumn("Total score",
                                                                       format="%.1f"),
                         "n_liabilities": st.column_config.NumberColumn("Motifs"),
                         "n_high": st.column_config.NumberColumn("High"),
                         "n_medium": st.column_config.NumberColumn("Medium"),
                         "n_low": st.column_config.NumberColumn("Low"),
                         "n_in_CDR": st.column_config.NumberColumn("In a CDR"),
                         "protease_score": st.column_config.NumberColumn(
                             "Protease score", format="%.2f"),
                         "n_protease_sites": st.column_config.NumberColumn(
                             "Protease sites")})

        img_map = png_liability_map(liab_by_seq, annotations,
                                    {n: len(s) for n, s in sequences.items()})
        if img_map:
            st.image(img_map, **FIT)
        img_profile = png_liability_profile(liab_by_seq)
        if img_profile:
            st.image(img_profile, **FIT)
        plt.close("all")

        st.caption(
            "Scores come from the in-house scale: N-glycosylation sequon 15, NG and DG "
            "10, NS/NT/NN and DS/DT/DD 5, free cysteine 8, Met 4, Trp 3; Low below 2, "
            "Medium below 6, High at or above 6. In its original form each motif is also "
            "weighted by its solvent exposure and the model confidence taken from a 3D "
            "structure. **This section is sequence-only, so neither factor is applied**, "
            "an NG scores 10 whether it is buried or fully exposed. The practical "
            "consequence is that the categories compress towards the top, and comparing "
            "sequences with each other tells you more than the category of any single "
            "motif.")

        st.subheader("Per sequence")
        for name in sequences:
            liab = liab_by_seq[name]
            summary = next(s for s in summaries if s["name"] == name)
            with st.container(border=True):
                st.markdown(f"**{name}**, score {summary['total_score']:.0f}, "
                            f"{summary['n_liabilities']} motif(s), "
                            f"{summary['n_high']} high, {summary['n_in_CDR']} in a CDR")
                if liab.empty:
                    st.markdown("No liability motif detected.")
                else:
                    st.dataframe(
                        liab[["liability", "motif", "range", "region", "in_CDR", "score",
                              "risk", "source", "note"]],
                        hide_index=True, **FIT,
                        column_config={
                            "range": st.column_config.TextColumn("Position"),
                            "score": st.column_config.NumberColumn("Score",
                                                                    format="%.1f"),
                            "note": st.column_config.TextColumn("Note", width="large")})

        with st.expander("Predicted protease cleavage sites"):
            st.markdown(
                "Every protein has tryptic sites, so the count on its own means little. "
                "What is worth looking at is a high-scoring recognition site, Factor Xa "
                "or enterokinase, or a cluster of sites inside a CDR.")
            for name in sequences:
                prot = prot_by_seq[name]
                st.markdown(f"**{name}**, {len(prot)} site(s), "
                            f"score {prot['score'].sum():.2f}" if not prot.empty
                            else f"**{name}**, none")
                if not prot.empty:
                    in_cdr = prot[prot["in_CDR"]]
                    st.caption(f"{len(in_cdr)} inside a CDR")
                    st.dataframe(prot.head(40), hide_index=True, **FIT,
                                 column_config={"note": st.column_config.TextColumn(
                                     "Note", width="large")})

        st.download_button(
            "Liabilities (CSV)",
            pd.concat([df.assign(sequence=n) for n, df in liab_by_seq.items()
                       if not df.empty], ignore_index=True).to_csv(index=False)
            if any(not d.empty for d in liab_by_seq.values()) else "",
            file_name="liabilities.csv", mime="text/csv")

        # ---- disulfides -------------------------------------------------------- #
        st.divider()
        st.subheader("Cysteines and disulfides")
        _ = heavy_only_notice("IMGT-numbered disulfide analysis")
        st.markdown(
            "The table above flags every cysteine outside the canonical FR1/FR3 pair "
            "as a possible free thiol. That is the right default, but it is wrong in "
            "one common case: camelid VHHs with long CDR3s often carry a **second, "
            "stabilising** disulfide, typically between CDR1 or FR2 and CDR3. This "
            "section tells the two situations apart. The canonical IMGT 23/104 pair "
            "is the same in a conventional VH, so that part applies to both; only "
            "the extra-pair reading is camelid.")
        for name, seq in heavy_sequences.items():
            try:
                ds = disulfides_for(seq, annotations.get(name))
            except Exception as exc:
                st.warning(f"**{name}**, cysteine analysis unavailable: {exc}")
                continue
            icon = "✅" if ds["n_orphan"] == 0 else "⚠️"
            with st.expander(f"{icon} **{name}**, {ds['n_cys']} cysteine(s): "
                             f"{ds['verdict']}", expanded=len(sequences) == 1):
                if len(ds["cysteines"]):
                    st.dataframe(ds["cysteines"], hide_index=True, **FIT)
                if len(ds["pairs"]):
                    st.dataframe(ds["pairs"], hide_index=True, **FIT)
                if ds["has_extra_pair"] and not kinds[name]["is_vhh"]:
                    st.info(
                        "An extra cysteine pair was found, but the CDR1–CDR3 pattern "
                        "this is matched against is a **camelid** one, and this "
                        "sequence is a conventional VH. Treat the pairing as "
                        "unverified rather than as the known VHH motif.")
                elif ds["has_extra_pair"]:
                    st.success(
                        "Read the extra cysteines as a stabilising pair rather than a "
                        "liability, but remember this is inferred from position "
                        "alone. Only a structure or a non-reducing gel confirms a bond.")
                elif ds["n_orphan"]:
                    st.warning(
                        "An unpaired cysteine is a real risk: scrambling, dimerisation "
                        "through an intermolecular bond, and heterogeneity on "
                        "purification. The repair table below offers C→S.")

        # ---- suggested repairs ------------------------------------------------- #
        st.divider()
        st.subheader("Suggested repairs")
        st.markdown(
            "One conservative substitution per liability, so the table above becomes "
            "actionable. Substitutions inside a CDR are flagged separately because "
            "they risk affinity, a liability you keep is often the right call.")
        rep_cols = st.columns([1, 1, 2])
        with rep_cols[0]:
            rep_include_cdr = st.checkbox("Include CDR positions", value=True)
        rep_target = (list(sequences)[0] if len(sequences) == 1 else
                      rep_cols[1].selectbox("Sequence", list(sequences),
                                            key="repair_target"))

        rep_seq = sequences[rep_target]
        rep_ann = annotations.get(rep_target)
        repairs = repairs_for(rep_seq, liab_by_seq[rep_target], rep_ann,
                              bool(rep_include_cdr))
        doable = repairs[repairs["applicable"]] if len(repairs) else repairs
        refused = repairs[~repairs["applicable"]] if len(repairs) else repairs

        if not len(repairs):
            st.success("Nothing to repair, no liability with a known substitution.")
        else:
            st.dataframe(doable, hide_index=True, **FIT,
                         column_config={"rationale": st.column_config.TextColumn(
                             "Rationale", width="large")})
            if len(refused):
                st.caption("Offered no substitution, deliberately:")
                st.dataframe(refused[["liability", "position", "region", "rationale"]],
                             hide_index=True, **FIT,
                             column_config={"rationale":
                                            st.column_config.TextColumn("Reason",
                                                                        width="large")})

            # Several liabilities offer alternatives at the SAME position (a free
            # cysteine can go to S or to A). Selecting both by default would make
            # every one of them collide, so the default keeps the first, the
            # preferred, option per position.
            _seen_pos, _default = set(), []
            for _, r in doable.iterrows():
                if r["in_CDR"] or r["position"] in _seen_pos:
                    continue
                _seen_pos.add(r["position"])
                _default.append(r["substitution"])
            picked = st.multiselect(
                "Apply substitutions", list(doable["substitution"]),
                default=_default,
                help="Defaults to one substitution per position, frameworks only, "
                     "the safe subset. Alternatives at the same position are listed "
                     "but not pre-selected, and picking two of them is refused "
                     "rather than silently resolved.")
            if picked:
                repaired, problems = vc.apply_repairs(rep_seq, picked)
                for p in problems:
                    st.error(p)
                before = liab_by_seq[rep_target]["score"].sum()
                # the original annotation is reused on purpose: substitutions do not
                # move the region boundaries, and re-annotating would hide a repair
                # that accidentally broke an anchor.
                after_df = liabilities_for(repaired, rep_ann, include_added)
                rm = st.columns(3)
                rm[0].metric("Liability score", f"{after_df['score'].sum():.0f}",
                             delta=f"{after_df['score'].sum() - before:+.0f}",
                             delta_color="inverse")
                rm[1].metric("Motifs", len(after_df),
                             delta=len(after_df) - len(liab_by_seq[rep_target]),
                             delta_color="inverse")
                rm[2].metric("Substitutions", len(picked))
                st.code(repaired, language="text")
                st.download_button(
                    "Repaired sequence (FASTA)",
                    f">{rep_target}_repaired {' '.join(picked)}\n{repaired}\n",
                    file_name=f"{rep_target}_repaired.fasta", mime="text/plain")
                st.caption(
                    "The score after repair is recomputed with the **original** "
                    "region boundaries, which is correct for substitutions but would "
                    "not be if you changed the length. Re-run the full analysis on "
                    "the repaired sequence before trusting anything else about it.")

    # --- Aggregation ----------------------------------------------------------- #
    with tabs["Aggregation"]:
        st.markdown(
            "Where along the sequence the aggregation-prone stretches are. A "
            "stretch becomes aggregation-prone when it is hydrophobic, β-prone and "
            "carries no charged residue to interrupt it, the three quantities this "
            "profile combines, plus a gatekeeper term for P/D/E/K/R.")
        st.warning(
            "**This is not CamSol, AGGRESCAN or TANGO.** Those fit their "
            "coefficients to experimental data; the weights here are a transparent "
            "judgement call, listed in `SOLUBILITY_WEIGHTS`. So read **where** the "
            "weak spots are and compare candidates against each other, do not read "
            "the absolute value as a predicted solubility. For a validated number, "
            "run the sequence through one of those servers.", icon="⚠️")

        agg_cols = st.columns(3)
        with agg_cols[0]:
            agg_window = st.slider("Smoothing window", 3, 15,
                                   int(vc.DEFAULTS["profile_window"]), step=2,
                                   help="Wider smoothing finds longer stretches and "
                                        "hides short ones.")
        with agg_cols[1]:
            agg_threshold = st.slider("Aggregation-prone threshold", -1.0, 0.5,
                                      float(vc.APR_THRESHOLD), step=0.05,
                                      help="A window scoring below this is called "
                                           "aggregation-prone. Negative leans prone, "
                                           "positive leans soluble.")
        with agg_cols[2]:
            agg_show_all = st.checkbox("Plot every sequence", value=len(sequences) <= 4,
                                       help="Untick with many sequences; the summary "
                                            "table is computed either way.")

        agg_summary, agg_regions_all = [], []
        for name, seq in sequences.items():
            profile, regions, summary, weakest = aggregation_for(
                seq, annotations.get(name), int(agg_window), float(agg_threshold))
            agg_summary.append({"name": name, **summary})
            if len(regions):
                agg_regions_all.append(regions.assign(sequence=name))

        st.dataframe(pd.DataFrame(agg_summary), hide_index=True, **FIT,
                     column_config={
                         "mean_score": st.column_config.NumberColumn("Mean", format="%.3f"),
                         "min_score": st.column_config.NumberColumn("Worst", format="%.3f"),
                         "worst_position": st.column_config.NumberColumn("At residue"),
                         "n_regions": st.column_config.NumberColumn("Regions"),
                         "residues_in_regions": st.column_config.NumberColumn("Residues"),
                         "fraction_prone": st.column_config.NumberColumn(
                             "Fraction", format="%.2f")})
        st.caption("`Mean` higher is better. `Worst` is the single lowest-scoring "
                   "window, which is usually what decides whether a candidate "
                   "behaves, one bad patch is enough.")

        for name, seq in sequences.items():
            profile, regions, summary, weakest = aggregation_for(
                seq, annotations.get(name), int(agg_window), float(agg_threshold))
            with st.expander(f"**{name}**, worst window {summary['min_score']:+.2f} "
                             f"at residue {summary['worst_position']}, "
                             f"{summary['n_regions']} region(s) below threshold",
                             expanded=agg_show_all and len(sequences) <= 4):
                if agg_show_all:
                    st.image(png_aggregation(seq, annotations.get(name),
                                             int(agg_window), float(agg_threshold),
                                             f"{name}, aggregation-prone regions"),
                             **FIT)
                st.markdown("**Weakest windows**")
                st.dataframe(weakest, hide_index=True, **FIT)
                if len(regions):
                    st.markdown("**Regions below the threshold**")
                    st.dataframe(regions, hide_index=True, **FIT)
                else:
                    st.success("No window crosses the threshold. The weakest-window "
                               "table above is still the one to compare candidates on.")

        dl_agg = st.columns(2)
        dl_agg[0].download_button(
            "Summary (CSV)", pd.DataFrame(agg_summary).to_csv(index=False),
            file_name="aggregation_summary.csv", mime="text/csv")
        dl_agg[1].download_button(
            "Regions (CSV)",
            pd.concat(agg_regions_all, ignore_index=True).to_csv(index=False)
            if agg_regions_all else "",
            file_name="aggregation_regions.csv", mime="text/csv",
            disabled=not agg_regions_all)
        plt.close("all")


def _compare(r: Run) -> None:
    """The compare family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    sequences = r.sequences
    annotations = r.annotations
    heavy_names = r.heavy_names
    other_names = r.other_names
    heavy_sequences = r.heavy_sequences
    heavy_annotations = r.heavy_annotations
    tabs = r.tabs
    heavy_only_notice = r.heavy_only_notice
    type_of = r.type_of

    # --- Alignment ------------------------------------------------------------- #
    with tabs["Alignment"]:
        # Aligning a 1200-residue protein against 120-residue VHHs produces a mostly
        # empty figure, so the default is the VHH subset, unless there are no VHHs,
        # in which case aligning the others is the only thing that makes sense.
        if other_names and heavy_names:
            align_all = st.checkbox(
                f"Include the {len(other_names)} sequence(s) that are not antibody "
                "domains", value=False,
                help="Off by default: mixing sequences of very different lengths "
                     "gives an alignment that is mostly gaps, and a logo to match. "
                     "VHH and conventional VH align together perfectly well and are "
                     "both included either way.")
        else:
            align_all = not heavy_names
        aln_seqs = sequences if align_all else heavy_sequences
        aln_anns = annotations if align_all else heavy_annotations
        if other_names and not align_all and heavy_names:
            st.caption(f"Aligning the {len(heavy_names)} heavy-chain domain(s) only.")

        if len(aln_seqs) < 2:
            st.info("Submit at least two sequences to see an alignment. "
                    + ("Tick the box above to include the non-VHH sequences."
                       if other_names and not align_all else ""))
        else:
            # These two figures used to sit behind a "draw" button, because every
            # tab body ran on every interaction and they are the most expensive in
            # the app. Grouping the tabs into families made that unnecessary: this
            # family only executes when it is the one open, so the button was pure
            # friction and has been removed.
            msa_wrap = st.slider(
                "Positions per row", 20, 140, 60, step=5, key="msa_wrap",
                help="How many alignment columns each row carries. Fewer "
                     "positions per row means more rows and wider letters.")
            with st.spinner("Aligning…"):
                # pymsaviz reads from a path, so the aligned FASTA is written to
                # a temporary directory removed as soon as the figure exists.
                img_msa, aln_text = png_msa(aln_seqs, aln_anns, int(msa_wrap))
            st.image(img_msa, **FIT)
            st.caption("FAMSA alignment rendered with pymsaviz. CDR bands refer "
                       "to the first sequence submitted.")
            st.download_button("Download alignment (FASTA)", aln_text,
                               file_name="alignment.fasta", mime="text/plain")

            st.subheader("Pairwise identity")
            cols_id = st.columns(2)
            with cols_id[0]:
                st.image(png_identity(aln_seqs, None, None), **FIT)
            if aln_anns and any((a.get("CDR3") or {}).get("seq")
                                for a in aln_anns.values()):
                img_id3 = png_identity(aln_seqs, aln_anns, "CDR3")
                if img_id3:
                    with cols_id[1]:
                        st.image(img_id3, **FIT)
            st.caption("Two clones from the same panning can share 90%+ full-length "
                       "identity while having unrelated CDR3s, the second matrix is the "
                       "one that tells you whether they are really the same binder.")

        if len(aln_seqs) >= 2:
            st.divider()
            st.subheader("Sequence logo")
            st.markdown(
                "The alignment as a logo: letter height shows what is conserved and "
                "what varies, which reads faster than the alignment itself when you "
                "only want to know where the diversity sits.")

            lg = st.columns([1, 1, 1])
            with lg[0]:
                logo_type = st.selectbox(
                    "Height", ["probability", "information"], index=0,
                    format_func=lambda t: "Frequency" if t == "probability"
                                          else "Information (bits)",
                    help="Frequency: letter height is the residue's share of the "
                         "column, so every column sums to 1 and reads as "
                         "composition. Information: height in bits, the classic "
                         "Schneider logo, where a mixed column stays short even if "
                         "one residue is commonest.")
            with lg[1]:
                logo_wrap = st.slider("Positions per row", 20, 140, 65, step=5,
                                      help="How many alignment columns each row "
                                           "carries. Fewer rows, wider letters.")
            with lg[2]:
                logo_scheme = st.selectbox(
                    "Colours", ["chemistry", "charge", "hydrophobicity", "skylign_protein"],
                    index=0)

            with st.spinner("Drawing the logo…"):
                img_logo, logo_matrix = png_logo(aln_seqs, aln_anns, logo_type,
                                                 int(logo_wrap), logo_scheme)
            st.image(img_logo, **FIT)
            st.caption("CDR bands refer to the first sequence submitted, mapped onto "
                       "alignment columns. Amino acids that never occur are dropped "
                       "from the matrix, so the colouring only shows residues that "
                       "are actually present.")
            st.download_button("Logo matrix (CSV)", logo_matrix.to_csv(),
                               file_name=f"logo_{logo_type}.csv", mime="text/csv")
            plt.close("all")

    # --- vs Reference ----------------------------------------------------------- #
    with tabs["vs Reference"]:
        st.markdown(
            "Which positions differ from a sequence you choose as the reference. "
            "Humanisation variants against the parent, a Sanger clone against the "
            "construct it should be, a panel against its ancestor, the same "
            "question each time.")
        if not heavy_only_notice("Reference comparison"):
            pass
        elif len(heavy_sequences) < 2 and not st.session_state.get("ref_paste", "").strip():
            st.info("Submit at least two sequences, or paste a reference below, to "
                    "compare against.")
            st.text_area("Reference sequence (optional)", key="ref_paste", height=90,
                         placeholder="Paste the parent sequence here")
        else:
            ref_cols = st.columns([2, 1])
            with ref_cols[0]:
                ref_choice = st.selectbox(
                    "Reference", list(heavy_sequences) + [",  paste one , "],
                    help="Everything else is compared against this sequence.")
            with ref_cols[1]:
                show_cdr_only = st.checkbox("CDR positions only", value=False)

            if ref_choice == ",  paste one , ":
                pasted = st.text_area("Reference sequence", key="ref_paste", height=90,
                                      placeholder="Paste the parent sequence here")
                ref_seq, ref_name = pasted.strip(), "pasted reference"
                others = dict(heavy_sequences)
            else:
                ref_seq, ref_name = heavy_sequences[ref_choice], ref_choice
                others = {k: v for k, v in heavy_sequences.items() if k != ref_choice}

            if not ref_seq:
                st.info("Paste a reference sequence to compare against.")
            elif not others:
                st.info("Nothing to compare, submit another sequence.")
            else:
                try:
                    comp = comparison_for(ref_seq, others, ref_name)
                except Exception as exc:
                    st.error(f"Could not compare: {exc}")
                    comp = None

                if comp:
                    for name, msg in comp["errors"].items():
                        st.warning(f"**{name}**, could not be numbered: {msg}")

                    st.dataframe(comp["summary"], hide_index=True, **FIT)
                    st.caption(
                        "Positions are matched by **IMGT number**, not by raw "
                        "alignment. A one-residue deletion in CDR3 is therefore "
                        "reported as one deletion rather than as every downstream "
                        "position having changed. Note that IMGT places insertions "
                        "at the centre of a loop, so adding a residue to CDR3 can "
                        "also shift a few neighbours, those show up as "
                        "substitutions alongside the insertion.")

                    diffs = comp["differences"]
                    if show_cdr_only and len(diffs):
                        diffs = diffs[diffs["in_CDR"]]

                    if not len(diffs):
                        st.success("No difference from the reference"
                                   + (" in the CDRs." if show_cdr_only else "."))
                    else:
                        st.image(png_differences(ref_seq, others, ref_name), **FIT)
                        st.dataframe(diffs, hide_index=True, **FIT)
                        st.download_button(
                            "Differences (CSV)", diffs.to_csv(index=False),
                            file_name=f"differences_vs_{ref_name}.csv",
                            mime="text/csv")
        plt.close("all")


    # --- Clonal families ------------------------------------------------------- #
    with tabs["Clonal families"]:
        st.markdown(
            "After a panning, the first question is not how good each clone is but "
            "**how many different clones you actually have**, a plate of picked "
            "colonies is routinely the same few lineages over and over. Grouping on "
            "CDR3 answers that, because CDR3 is what differs between independent "
            "lineages while the frameworks stay nearly constant.")

        if not heavy_only_notice("Clonal family grouping"):
            pass
        elif len(heavy_sequences) < 2:
            st.info("Submit at least two VHH sequences to group them into families.")
        elif not any((a.get("CDR3") or {}).get("seq")
                     for a in heavy_annotations.values()):
            st.warning("No CDR3 could be annotated, so there is nothing to group on.")
        else:
            fam_cols = st.columns([2, 1, 1])
            with fam_cols[0]:
                cutoff = st.slider(
                    "CDR3 identity threshold (%)", min_value=50, max_value=100,
                    value=int(vc.FAMILY_DEFAULTS["cutoff"]), step=1,
                    help="Two sequences at or above this identity belong to the same "
                         "family. 100% groups only strictly identical CDR3s.")
            with fam_cols[1]:
                linkage_method = st.selectbox(
                    "Linkage", ["complete", "average"],
                    index=0 if vc.FAMILY_DEFAULTS["linkage"] == "complete" else 1,
                    help="complete: every pair inside a family meets the threshold. "
                         "average: merges on the mean distance, more permissive, a "
                         "family may then contain pairs below the threshold.")
            with fam_cols[2]:
                same_length = st.checkbox(
                    "Same CDR3 length only",
                    value=bool(vc.FAMILY_DEFAULTS["require_same_length"]),
                    help="Forces CDR3s of different lengths apart whatever their "
                         "alignment. Off by default: a one-residue indel is a "
                         "plausible sibling, not a new lineage.")

            with st.spinner("Grouping…"):
                fam = families_for(heavy_sequences, heavy_annotations, float(cutoff),
                                   linkage_method, bool(same_length))

            n_seqs = len(fam["assignment"])
            n_fam = fam["n_families"]
            singletons = int((fam["families"]["n_sequences"] == 1).sum())
            head = st.columns(4)
            head[0].metric("Sequences grouped", n_seqs)
            head[1].metric("Families", n_fam)
            head[2].metric("Singletons", singletons)
            head[3].metric("Redundancy", f"{n_seqs / n_fam:.2f}" if n_fam else ", ",
                           help="Sequences per family. Above ~2 you are re-picking "
                                "the same lineages and could screen more colonies.")

            if fam["skipped"]:
                st.warning("Left out for want of an annotated CDR3: "
                           + ", ".join(fam["skipped"]))

            st.subheader("Families")
            st.dataframe(fam["families"], hide_index=True, **FIT)
            st.caption(
                "`min_identity` is the **worst** pair inside the family: with "
                "*complete* linkage it can never fall below your threshold, so a "
                "value under it means *average* linkage merged across it. "
                "`representative` is the family's most frequent CDR3, a reasonable "
                "default pick for the clone to carry forward.")

            img_dendro = png_family_dendrogram(heavy_sequences, heavy_annotations,
                                               float(cutoff),
                                               linkage_method, bool(same_length))
            if img_dendro:
                st.subheader("Dendrogram")
                st.image(img_dendro, **FIT)
                st.caption("The dashed line is your threshold; branches joining to "
                           "its right are separate families. Move the slider and "
                           "watch which groups merge, that is the honest way to see "
                           "how robust a family is.")

            st.subheader("Per-sequence assignment")
            st.dataframe(fam["assignment"].assign(
                type=[type_of(x) for x in fam["assignment"]["name"]]),
                hide_index=True, **FIT)

            dl = st.columns(3)
            dl[0].download_button("Families (CSV)",
                                  fam["families"].to_csv(index=False),
                                  file_name=f"families_{int(cutoff)}pct.csv",
                                  mime="text/csv")
            dl[1].download_button("Assignment (CSV)",
                                  fam["assignment"].to_csv(index=False),
                                  file_name=f"family_assignment_{int(cutoff)}pct.csv",
                                  mime="text/csv")
            dl[2].download_button("CDR3 identity matrix (CSV)",
                                  fam["identity"].round(1).to_csv(),
                                  file_name="cdr3_identity_matrix.csv",
                                  mime="text/csv")

            with st.expander("One representative per family, as FASTA"):
                reps = {r: heavy_sequences[r]
                        for r in fam["families"]["representative"]
                        if r in heavy_sequences}
                fasta = "".join(f">{k}\n{v}\n" for k, v in reps.items())
                st.code(fasta or ", ", language="text")
                st.download_button("Download representatives (FASTA)", fasta,
                                   file_name="family_representatives.fasta",
                                   mime="text/plain")
            plt.close("all")


def _repertoire(r: Run) -> None:
    """The repertoire family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    heavy_sequences = r.heavy_sequences
    heavy_annotations = r.heavy_annotations
    _vh_only = r._vh_only
    tabs = r.tabs
    heavy_only_notice = r.heavy_only_notice

    # --- Repertoire ------------------------------------------------------------- #
    with tabs["Repertoire"]:
        st.markdown(
            "Your sequences against the natural VHH repertoire, the same "
            "position-frequency matrix the synthetic generator samples from, used "
            "here the other way round.")
        if not heavy_only_notice("Repertoire comparison"):
            pass
        else:
            if _vh_only:
                st.warning(
                    "The reference is a **VHH** repertoire. A conventional VH "
                    "compared against it will legitimately show its FR2 hallmark "
                    "positions as rare, that is the difference between a VH and a "
                    "VHH, not a defect in the sequence.", icon="🧬")

            rar, nat = st.tabs(["Unusual framework residues",
                                "Your set vs the natural repertoire"])

            with rar:
                st.markdown(
                    "Framework positions whose residue is rare in the reference. "
                    "Worth a second look: this is what catches a Sanger miscall, a "
                    "PCR error, or a construct that has drifted without anyone "
                    "noticing.")
                rc = st.columns(3)
                with rc[0]:
                    rare_thr = st.slider("Rare at or below", 0.001, 0.10,
                                         float(rep.RARE_THRESHOLD), step=0.001,
                                         format="%.3f")
                with rc[1]:
                    min_sup = st.slider("Minimum support", 0, 1900,
                                        int(rep.MIN_SUPPORT), step=50,
                                        help="Positions backed by fewer reference "
                                             "sequences than this are set aside "
                                             "rather than called rare.")
                with rc[2]:
                    rare_plot = st.checkbox("Plot each sequence",
                                            value=len(heavy_sequences) <= 4)

                summary, flagged = rarity_tables(heavy_sequences, float(rare_thr),
                                                 int(min_sup))
                st.dataframe(summary, hide_index=True, **FIT)
                st.caption(
                    "**Rare is not wrong.** A residue absent from the reference may "
                    "be a genuine novelty, a germline this set under-samples, or the "
                    "very mutation you selected for. Every call carries the number "
                    "of reference sequences behind it.")

                if len(flagged):
                    st.subheader("Flagged positions")
                    st.dataframe(flagged, hide_index=True, **FIT)
                    st.download_button("Unusual residues (CSV)",
                                       flagged.to_csv(index=False),
                                       file_name="unusual_framework_residues.csv",
                                       mime="text/csv")
                else:
                    st.success("No framework residue falls below the threshold, "
                               "every position matches the natural repertoire.")

                if rare_plot:
                    for name, seq in heavy_sequences.items():
                        try:
                            st.image(png_rarity(seq, name, float(rare_thr),
                                                int(min_sup)), **FIT)
                        except Exception as exc:
                            st.warning(f"**{name}**, {exc}")

            with nat:
                st.markdown(
                    "Where your set sits against the reference distribution. After "
                    "a panning this is what says whether the selection was biased, "
                    "towards long CDR3s, towards charged loops, and so on.")
                if not heavy_annotations:
                    st.info("Nothing annotated to compare.")
                else:
                    per_seq, distribution, composition = repertoire_stats(
                        heavy_sequences, heavy_annotations)
                    st.image(png_cdr_lengths(heavy_sequences, heavy_annotations), **FIT)
                    if len(heavy_sequences) >= rep.MIN_FOR_BARS:
                        st.caption(
                            "Grey is the natural repertoire, red your set, both as "
                            "fractions of their own total.")
                    else:
                        st.caption(
                            f"Grey is the natural repertoire; the red markers are "
                            f"your {len(heavy_sequences)} sequence(s). Below "
                            f"{rep.MIN_FOR_BARS} sequences they are drawn as markers "
                            "rather than bars, because a single sequence is 100% of "
                            "itself: its bar would reach 1.0 and flatten the "
                            "distribution it is meant to be compared against. The "
                            "grey figure next to each marker is how common that "
                            "length is among real VHHs.")

                    context = rep.cdr_length_context(per_seq, distribution)
                    if len(context):
                        st.markdown("**Where each CDR length sits**")
                        st.dataframe(context, hide_index=True, **FIT,
                                     column_config={
                                         "reference_frequency":
                                             st.column_config.NumberColumn(
                                                 "Frequency in nature",
                                                 format="%.3f"),
                                         "percentile": st.column_config.NumberColumn(
                                             "Percentile", format="%.0f"),
                                         "natural_mode":
                                             st.column_config.NumberColumn(
                                                 "Most common length")})
                        st.caption(
                            "`Percentile` is the fraction of natural VHHs with a "
                            "shorter CDR of that region. The verdict calls a length "
                            "unusual below the 5th or above the 95th percentile, "
                            "unusual is not wrong, and a long CDR3 is often exactly "
                            "what you selected for.")
                    st.dataframe(per_seq, hide_index=True, **FIT)

                    st.subheader("Amino-acid composition")
                    regions = st.multiselect(
                        "Regions", ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"],
                        default=["CDR1", "CDR2", "CDR3"])
                    if regions:
                        st.image(png_composition(heavy_sequences, heavy_annotations,
                                                 tuple(regions)), **FIT)
                        st.caption(
                            "Above the line: enriched in your set relative to the "
                            "natural repertoire. Below: depleted. The reference "
                            "composition is the position-frequency matrix averaged "
                            "over each region, weighted by support.")
                        shown = composition[composition["region"].isin(regions)]
                        st.dataframe(
                            shown.sort_values("delta", ascending=False),
                            hide_index=True, **FIT)
                        st.download_button("Composition (CSV)",
                                           composition.to_csv(index=False),
                                           file_name="composition_vs_repertoire.csv",
                                           mime="text/csv")
        plt.close("all")


def _outputs(r: Run) -> None:
    """The outputs family."""
    # One line each rather than one wrapped tuple: this is the whole of
    # what this family needs from the run, and it reads as a checklist.
    sequences = r.sequences
    df = r.df
    annotations = r.annotations
    ph = r.ph
    window = r.window
    tabs = r.tabs
    chosen = r.chosen

    # --- DNA -------------------------------------------------------------------- #
    with tabs["DNA"]:
        st.markdown("Back-translation into DNA using *E. coli* preferred codons, with "
                    "cloning sites and long homopolymer runs removed.")

        col_a, col_b = st.columns([2, 1])
        with col_a:
            avoid = st.multiselect(
                "Restriction sites to avoid", sorted(vc.DEFAULT_AVOID_SITES),
                default=sorted(vc.DEFAULT_AVOID_SITES),
                help="A codon creating one of these sites is swapped for a synonymous "
                     "one.")
        with col_b:
            add_stop = st.checkbox("Append a stop codon (TAA)", value=True)
            gc_max = st.slider("Maximum GC fraction", 0.45, 0.70, 0.60, 0.01)

        sites = {k: vc.DEFAULT_AVOID_SITES[k] for k in avoid}
        fasta_parts = []
        for name, seq in sequences.items():
            res = dna_for(seq, tuple(sorted(sites.items())), add_stop, gc_max)
            ok = vc.translate_dna(res["dna"]) == seq
            fasta_parts.append(f">{name}_optimized\n" +
                               "\n".join(res["dna"][i:i + 60]
                                         for i in range(0, len(res["dna"]), 60)))
            with st.container(border=True):
                st.markdown(f"**{name}**, {res['length_nt']} nt, "
                            f"GC {res['gc_percent']} %, "
                            f"{len(res['repairs'])} codon(s) swapped for sites, "
                            f"{res['gc_swaps']} for GC balance")
                if not ok:
                    st.error("Round-trip check failed, do not order this sequence.")
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
            "cloning-unfriendly motifs, it does not model mRNA secondary structure near "
            "the ribosome binding site, codon-pair bias or tRNA availability. Run the "
            "sequence through your synthesis vendor's checker before ordering.")
        dna_out = st.columns(2)
        with dna_out[0]:
            st.button("🧬 Use these in the cloning section", key="to_cloning", **FIT,
                      on_click=go_to, args=("cloning",))
        with dna_out[1]:
            if st.button("🔧 Send to gene diagnosis", key="to_genes", **FIT,
                         help="Run these sequences through the frame, rare-codon and "
                              "GC checks, and get synonymous swaps where GC runs "
                              "high."):
                st.session_state["pending_gene_dna"] = "\n".join(fasta_parts) + "\n"
                go_to("genes")
                st.rerun()

    # --- Constructs ------------------------------------------------------------- #
    with tabs["Constructs"]:
        st.markdown(
            "Assemble a VHH into the construct you will actually order, monovalent, "
            "bivalent or biparatopic, with a linker and tags, and get the mass, pI "
            "and DNA of the finished protein rather than of the bare domain.")
        st.warning(ctm.VERIFY_NOTICE, icon="⚠️")

        fmt_cols = st.columns([1, 1, 1])
        with fmt_cols[0]:
            fmt = st.selectbox("Format", list(ctm.FORMATS),
                               help="\n\n".join(f"**{k}**, {v['note']}"
                                                 for k, v in ctm.FORMATS.items()))
        n_needed = ctm.FORMATS[fmt]["n_vhh"]
        with fmt_cols[1]:
            first = st.selectbox("VHH", list(sequences), key="ctor_first")
        second = None
        if n_needed == 2:
            with fmt_cols[2]:
                if fmt == "bivalent":
                    st.caption("Bivalent uses the same VHH twice.")
                    second = first
                else:
                    second = st.selectbox(
                        "Second VHH", list(sequences),
                        index=min(1, len(sequences) - 1), key="ctor_second")
                    if second == first:
                        st.caption("Two copies of the same VHH is a bivalent "
                                   "construct, not a biparatopic one.")

        lk_cols = st.columns([1, 1])
        with lk_cols[0]:
            linker_name = st.selectbox(
                "Linker between VHHs", list(ctm.LINKERS),
                index=list(ctm.LINKERS).index("(GGGGS)3"),
                disabled=n_needed == 1)
            st.caption(ctm.LINKERS[linker_name]["note"])
        with lk_cols[1]:
            tag_linker = st.text_input(
                "Linker before each tag", value="GS",
                help="Leave empty to fuse tags directly. Two or three residues is "
                     "usual and often just the cloning scar.")

        tag_cols = st.columns(2)
        with tag_cols[0]:
            n_tag_names = st.multiselect(
                "N-terminal tags",
                [k for k, v in ctm.TAGS.items() if v["ends"] in ("N", "both")],
                help="Applied in the order you pick them, reading N to C.")
        with tag_cols[1]:
            c_tag_names = st.multiselect(
                "C-terminal tags",
                [k for k, v in ctm.TAGS.items() if v["ends"] in ("C", "both")],
                default=["His6"])

        with st.expander("Custom pieces, SpyTag2/3, a camelid hinge, anything else"):
            st.markdown(
                "Two things are deliberately missing from the tag list above: "
                "**SpyTag2 and SpyTag3**, and the **camelid hinges**. Only the "
                "original SpyTag is included. Writing the others from memory would "
                "risk handing you a sequence that is one residue wrong, which you "
                "would discover only after the gene arrived, so paste the exact "
                "sequence from your own construct here instead.")
            cu = st.columns(2)
            custom_n = cu[0].text_input("Custom N-terminal peptide", value="")
            custom_n_label = cu[0].text_input("Its label", value="custom-N")
            custom_c = cu[1].text_input("Custom C-terminal peptide", value="")
            custom_c_label = cu[1].text_input("Its label", value="custom-C")

        chosen = [(first, sequences[first])]
        if n_needed == 2 and second:
            chosen.append((second, sequences[second]))

        n_tags = [(k, ctm.TAGS[k]["seq"]) for k in n_tag_names]
        c_tags = [(k, ctm.TAGS[k]["seq"]) for k in c_tag_names]
        if custom_n.strip():
            n_tags.append((custom_n_label or "custom-N", custom_n))
        if custom_c.strip():
            c_tags.append((custom_c_label or "custom-C", custom_c))

        try:
            built = ctm.build_construct(
                chosen, linker=ctm.LINKERS[linker_name]["seq"] if n_needed == 2 else "",
                n_tags=n_tags, c_tags=c_tags, tag_linker=tag_linker, ph=ph)
        except Exception as exc:
            st.error(f"Could not assemble the construct: {exc}")
            built = None

        if built:
            for w in built["warnings"][1:]:        # [0] is the verify notice, shown above
                st.warning(w, icon="⚠️")

            st.image(png_construct(built, f"{first}, {fmt}"), **FIT)

            props = built["properties"].iloc[0]
            base = built["per_vhh"].iloc[0]
            pm = st.columns(4)
            pm[0].metric("Length", f"{built['length']} aa",
                         delta=f"{built['length'] - int(base['length']):+d} vs the VHH")
            pm[1].metric("Mass", f"{props['mw_kda']:.2f} kDa",
                         delta=f"{props['mw_kda'] - base['mw_kda']:+.2f}")
            pm[2].metric("pI", f"{props['pI']:.2f}",
                         delta=f"{props['pI'] - base['pI']:+.2f}")
            pm[3].metric("A280 (1 mg/mL)", f"{props['A280_0.1pct']:.3f}",
                         help="Divide a measured A280 by this to get mg/mL, the "
                              "number changes with the tags, which is exactly why "
                              "you want it for the construct and not the bare VHH.")

            st.markdown("**Segments**")
            st.dataframe(built["segments"], hide_index=True, **FIT)

            st.markdown("**Full properties**")
            st.dataframe(built["properties"], hide_index=True, **FIT)

            st.markdown("**Protein sequence**")
            st.code(built["sequence"], language="text")

            with st.expander("Codon-optimised DNA for the whole construct"):
                avoid_ct = st.multiselect(
                    "Restriction sites to avoid", sorted(vc.DEFAULT_AVOID_SITES),
                    default=sorted(vc.DEFAULT_AVOID_SITES), key="ctor_avoid")
                stop_ct = st.checkbox("Append a stop codon (TAA)", value=True,
                                      key="ctor_stop")
                dna_ct = ctm.construct_dna(built, avoid=set(avoid_ct),
                                           add_stop=bool(stop_ct))
                if dna_ct["back_translation_ok"]:
                    st.success("Round trip verified: translating this DNA gives back "
                               "the construct exactly.")
                else:
                    st.error("The DNA does not translate back to the construct, do "
                             "not order this. Please report it.")
                st.code(dna_ct["dna"], language="text")
                st.download_button(
                    "Construct DNA (FASTA)",
                    f">{first}_{fmt}\n{dna_ct['dna']}\n",
                    file_name=f"{first}_{fmt}.fasta", mime="text/plain")

            dlc = st.columns(2)
            dlc[0].download_button(
                "Construct protein (FASTA)",
                f">{first}_{fmt}\n{built['sequence']}\n",
                file_name=f"{first}_{fmt}_protein.fasta", mime="text/plain")
            dlc[1].download_button(
                "Segments (CSV)", built["segments"].to_csv(index=False),
                file_name=f"{first}_{fmt}_segments.csv", mime="text/csv")
        plt.close("all")


    # --- Download -------------------------------------------------------------- #
    with tabs["Download"]:
        st.markdown("Export the full results table, including every region sequence, "
                    "custom feature and flag.")
        st.download_button("Results (CSV)", df.to_csv(index=False),
                           file_name="vhh_results.csv", mime="text/csv", type="primary")
        st.download_button("Secondary structure (PNG)", png_secondary_structure(df),
                           file_name="secondary_structure.png", mime="image/png")
        st.download_button("Charge (PNG)", png_charge(df, sequences, ph),
                           file_name="charge.png", mime="image/png")
        st.download_button("Profiles (PNG)", png_profiles(sequences, annotations, window),
                           file_name="profiles.png", mime="image/png")

        st.divider()
        st.markdown("### Full PDF report")
        st.markdown("Everything in one document: properties, developability flags, "
                    "regions, hallmarks, germline, liabilities and all figures. Building "
                    "it takes a few seconds, so it is generated on demand.")

        include_cf = st.checkbox(
            "Include the per-sequence secondary-structure plots",
            value=len(sequences) <= 6,
            help="One extra page per sequence. With many sequences this is most of the "
                 "file size and most of the build time.")

        if st.button("Build PDF report", type="secondary"):
            # The displayed plots are cached as PNG bytes, so the report rebuilds real
            # figure objects, only when the button is pressed.
            figures = {
                "Secondary structure composition": vc.plot_secondary_structure(df),
                "Charge": vc.plot_charge(df, sequences, ph=ph),
                "Charge and hydropathy profiles": vc.plot_profiles(
                    sequences, annotations, window=window),
            }
            if include_cf:
                for cf_name, cf_seq in sequences.items():
                    figures[f"Predicted secondary structure, {cf_name}"] = \
                        vc.plot_chou_fasman(cf_seq, cf_name, window=cf_window,
                                            annotation=annotations.get(cf_name))
            if len(sequences) >= 2:
                try:
                    with tempfile.TemporaryDirectory() as tmp:
                        fig_msa_pdf, _ = vc.plot_msa(sequences, annotations, workdir=tmp)
                    figures["Multiple sequence alignment"] = fig_msa_pdf
                    figures["Pairwise identity"] = vc.plot_identity(
                        vc.identity_matrix(sequences), "Full-length identity (%)")
                    mat3 = vc.identity_matrix(sequences, annotations, region="CDR3")
                    if len(mat3) >= 2:
                        figures["CDR3 identity"] = vc.plot_identity(mat3,
                                                                    "CDR3 identity (%)")
                except Exception as exc:
                    st.warning(f"Alignment figures skipped: {exc}")

            with st.spinner("Building the report…"):
                with tempfile.TemporaryDirectory() as tmp:
                    pdf_path = vc.build_pdf_report(
                        f"{tmp}/vhh_report.pdf", df, sequences, annotations,
                        figures=figures, ph=ph, version=VERSION, author=AUTHOR,
                        app_name=APP_NAME)
                    st.session_state["pdf_bytes"] = pdf_path.read_bytes()
            plt.close("all")

        if st.session_state.get("pdf_bytes"):
            st.download_button("Download report (PDF)", st.session_state["pdf_bytes"],
                               file_name="vhh_report.pdf", mime="application/pdf",
                               type="primary")

        st.divider()
        st.dataframe(df.drop(columns=["sequence"]), hide_index=True, **FIT)


def render(ph: float, window: int, method: str) -> None:
    """Draw the Sequence analysis screen."""

    # --------------------------------------------------------------------------- #
    # Sequence analysis
    # --------------------------------------------------------------------------- #

    st.title("🔬 Sequence analysis")
    st.markdown(
        "Physico-chemical characterization, CDR/framework annotation and alignment. "
        "Built around VHHs, but **any protein is accepted**, each sequence is "
        "classified, and the VHH-specific analyses only run on the ones that are "
        "VHHs. Paste one or more sequences, in FASTA or raw form.")

    col_in, col_act = st.columns([3, 1])
    with col_act:
        st.write("")
        if st.button("Load example", **FIT):
            st.session_state["seq_text"] = EXAMPLE
            st.rerun()
        if st.button("🗑️ Clear", **FIT,
                     help="Empty the sequence box.",
                     disabled=not st.session_state.get("seq_text", "").strip()):
            st.session_state["seq_text"] = ""
            # forget the last upload too, otherwise re-selecting the same file after
            # clearing does nothing: the uploader still considers it already handled
            st.session_state["last_upload"] = None
            st.rerun()
        uploaded = st.file_uploader("or upload FASTA", type=["fasta", "fa", "txt", "faa"])
        if uploaded is not None and uploaded.file_id != st.session_state.get("last_upload"):
            st.session_state["last_upload"] = uploaded.file_id
            st.session_state["seq_text"] = uploaded.read().decode("utf-8", errors="replace")
            st.rerun()
        st.caption(f"Up to {vc.DEFAULTS['max_sequences']} sequences, "
                   f"{vc.DEFAULTS['max_length']} aa each.")
        st.button("🧪 Get them from .ab1 files", key="to_sanger", **FIT,
                  on_click=go_to, args=("sequencing",))

    with col_in:
        # The widget owns st.session_state["seq_text"]; the buttons above write to it
        # before the widget is instantiated, which is the supported way to prefill it.
        text = st.text_area("Sequences", height=220, key="seq_text",
                            placeholder=">MyNanobody\nEVQLVESGGGLVQ...")
        # Clicking the button moves focus out of the text area, which is what commits
        # what was just typed, the same thing ⌘/Ctrl + Enter does, without the
        # keyboard shortcut nobody discovers on their own.
        st.button("Analyze", type="primary", **FIT)
        st.caption("Click Analyze, or press ⌘ + Enter (Ctrl + Enter on Windows/Linux).")
        # Filled in later by analyze_with_progress: the progress bar while it runs,
        # then the result line. Created here so it appears under the button rather
        # than wherever in the script the analysis happens to be called.
        status_slot = st.empty()

    # --- Synthetic sequence generator ------------------------------------------- #
    with st.expander("🎲 Generate synthetic VHH sequences"):
        matrix, gmeta = generator_matrix()
        report = gn.matrix_report(matrix, gmeta)
        st.markdown(
            "Frameworks are sampled position by position from a real IMGT "
            "position-frequency matrix, so they look like natural frameworks. CDRs are "
            "drawn at random, with lengths taken from the same reference set.")
        st.caption(
            f"Reference: **{report['source']}**, {report['n_sequences']} unique VHH "
            f"sequences, {report['n_positions']} IMGT positions, median support "
            f"{report['median_support']} sequences per position.")
        st.caption(
            f"*{gmeta.get('source_title', '')}*, "
            f"[doi:{gmeta.get('source_doi', '')}]"
            f"(https://doi.org/{gmeta.get('source_doi', '')}) · "
            f"[database]({gmeta.get('source_url', '')}). "
            "Please cite it if you publish work based on these sequences. What is "
            "bundled here is an aggregate frequency matrix, not the sequences "
            "themselves.")

        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            n_gen = st.number_input("How many", 1, 200, 10, 1)
            seed = st.number_input("Random seed", 0, 10**6, 0, 1,
                                   help="The same seed always gives the same set.")
        with gc2:
            fw_mode = st.radio("Frameworks", ["sample", "consensus"],
                               help="'sample' draws each position from its frequency "
                                    "column and gives a different framework every "
                                    "time; 'consensus' takes the most frequent "
                                    "residue, giving one canonical framework.")
            cdr_alpha = st.radio(
                "CDR residues", ["natural", "no_cys", "uniform"],
                help="'natural' uses the reference set's own CDR composition; "
                     "'no_cys' is uniform over 19 residues; 'uniform' is all 20.")
        with gc3:
            keep_hallmarks = st.checkbox("Enforce camelid hallmarks", value=True)
            allow_fw_cys = st.checkbox("Allow cysteines in frameworks", value=False)
            fixed_lengths = st.checkbox("Fix the CDR lengths", value=False)

        lengths = {}
        if fixed_lengths:
            lc = st.columns(3)
            for col, region, default in zip(lc, ("CDR1", "CDR2", "CDR3"), (8, 8, 17)):
                with col:
                    lengths[region] = st.number_input(f"{region} length", 2, 30,
                                                      default, 1)
        else:
            st.caption("CDR lengths are drawn from the reference distribution, "
                       "CDR1 and CDR2 around 8 residues, CDR3 a mean of 17 with a "
                       "range of 2 to 30.")

        if st.button("Generate", type="primary"):
            with st.spinner(f"Generating {int(n_gen)} sequence(s)…"):
                st.session_state["generated"] = generate_cached(
                    int(n_gen), int(seed), fw_mode, cdr_alpha,
                    tuple(sorted(lengths.items())), keep_hallmarks, allow_fw_cys)

        made = st.session_state.get("generated")
        if made is not None and len(made):
            st.dataframe(made[["name", "length", "tetrad", "CDR1_len", "CDR2_len",
                               "CDR3_len", "sequence"]], hide_index=True, **FIT,
                         column_config={"sequence": st.column_config.TextColumn(
                             "Sequence", width="large")})
            gcols = st.columns(2)
            with gcols[0]:
                st.download_button("Download (FASTA)", gn.to_fasta(made),
                                   file_name="synthetic_vhh.fasta", mime="text/plain",
                                   **FIT)
            with gcols[1]:
                if st.button("Send to the analysis above", **FIT):
                    st.session_state["pending_seq_text"] = gn.to_fasta(made)
                    st.rerun()

        st.caption(
            "Two limitations worth knowing. **Positions are sampled independently**, so "
            "a combination that never co-occurs in nature can come out, the hallmark "
            "tetrad is the usual casualty, which is why it is enforced by default. And "
            "**a frequency is only as good as its support**: a position seen in 30 "
            "reference sequences carries far less information than one seen in 1900. "
            "These are synthetic sequences for controls and benchmarks; nothing here "
            "predicts that they fold or bind anything.")

    if not text.strip():
        st.info("Paste sequences above, upload a FASTA file, load the example, generate "
                "synthetic ones, or extract them from .ab1 chromatograms in the "
                "sequencing section.")
        footer()
        st.stop()

    try:
        sequences = vc.parse_input(text)
    except ValueError as exc:
        st.error(f"Input error, {exc}")
        footer()
        st.stop()

    df, annotations, errors = analyze_with_progress(sequences, ph, method,
                                                    slot=status_slot)

    # What is each sequence? The app accepts any protein, so every VHH-specific
    # analysis has to know which sequences it legitimately applies to.
    kinds = {name: classify_cached(seq, method) for name, seq in sequences.items()}
    # Heavy-chain V domains (VHH *and* conventional VH) get the full analysis.
    # Almost none of it is camelid-specific: IMGT numbering, the CDR windows,
    # germline assignment, humanness and CDR3 grouping are the same operation on a
    # VH. The two camelid readings (the hallmark tetrad, the extra CDR1-CDR3
    # disulfide) are labelled in place instead of gating the rest.
    heavy_names = [n for n in sequences if kinds[n]["full_analysis"]]
    vhh_names = [n for n in sequences if kinds[n]["is_vhh"]]
    other_names = [n for n in sequences if not kinds[n]["full_analysis"]]
    heavy_sequences = {n: sequences[n] for n in heavy_names}
    heavy_annotations = {n: annotations[n] for n in heavy_names if n in annotations}

    df.insert(1, "type", [kinds[n]["kind"] for n in df["name"]])

    _vh_only = [n for n in heavy_names if not kinds[n]["is_vhh"]]
    if _vh_only:
        st.info(
            f"**{len(_vh_only)} sequence(s) are conventional VH, not VHH.** They go "
            "through the full heavy-chain analysis all the same, CDR/framework, "
            "germline, humanness, clonal families, and are labelled as VH "
            "throughout. Only two readings are camelid-specific and are flagged "
            "where they appear: the FR2 hallmark tetrad, and the additional "
            "CDR1–CDR3 disulfide.\n\n"
            + "\n".join(f"- **{n}**, {kinds[n]['note']}" for n in _vh_only),
            icon="🧬")

    if other_names:
        st.info(
            f"**{len(other_names)} of {len(sequences)} sequence(s) are not "
            "heavy-chain variable domains.** They are characterised like any "
            "protein (mass, pI, extinction, charge, aggregation, liabilities, DNA) "
            "but the antibody-specific tabs skip them.\n\n"
            + "\n".join(f"- **{n}**, {kinds[n]['kind']}: {kinds[n]['note']}"
                         for n in other_names), icon="ℹ️")

    # An annotation failure on something that is not an antibody domain is the
    # expected outcome, not an error worth alarming anyone about.
    real_errors = {n: m for n, m in errors.items() if kinds[n]["full_analysis"]}
    if real_errors:
        for name, msg in real_errors.items():
            st.warning(f"**{name}**, annotation failed: {msg}. "
                       "Physico-chemical values are still reported below.")


    def heavy_only_notice(what: str) -> bool:
        """Guard for a panel that needs a heavy-chain V domain (VHH or VH)."""
        if not heavy_names:
            st.info(f"{what} applies to heavy-chain variable domains, VHH or "
                    "conventional VH, and none of the submitted sequences is one. "
                    "The generic tabs, Properties, Plots, Liabilities, Aggregation, "
                    "DNA, Constructs, cover them all.")
            return False
        if other_names:
            st.caption(f"Showing the {len(heavy_names)} heavy-chain domain(s); "
                       f"{len(other_names)} other sequence(s) are not analysed here.")
        return True


    def type_of(name: str) -> str:
        """Short label for a sequence's kind, for the per-sequence tables."""
        return "VHH" if kinds[name]["is_vhh"] else kinds[name]["kind"]

    # Named rather than indexed: inserting a tab in the middle of an indexed list
    # silently reassigns every panel after it.
    # Fourteen tabs are not fourteen subjects, they are four questions plus the
    # decision view. Grouping them does two things at once: a newcomer reads five
    # labels instead of scanning fourteen, and Streamlit only executes the family
    # that is open. Nested tabs would have fixed the first problem and not the
    # second, because every tab body runs whether or not its tab is selected.
    FAMILIES = {
        "🎯 Triage": ("triage", [],
                      "Every flag in one table, the view for deciding which "
                      "sequences to carry forward."),
        "🔍 Identity": ("identity",
                        ["Properties", "Plots", "CDR / Framework",
                         "Hallmarks & germline"],
                        "What each sequence is: physico-chemical descriptors, "
                        "profiles, CDR and framework annotation, hallmarks and "
                        "humanness."),
        "⚠️ Developability": ("risk", ["Liabilities", "Aggregation"],
                              "Whether it will behave: chemical liabilities, "
                              "cysteines and disulfides, suggested repairs, "
                              "aggregation-prone regions."),
        "🔗 Comparison": ("compare", ["Alignment", "vs Reference",
                                      "Clonal families"],
                          "How your sequences compare with each other: alignment "
                          "and sequence logo, differences from a reference, clonal "
                          "families."),
        "🌍 Repertoire": ("repertoire", ["Repertoire"],
                          "How they compare with the natural VHH repertoire: "
                          "unusual framework residues, CDR lengths and composition "
                          "against 1958 real VHHs."),
        "🧪 Outputs": ("outputs", ["DNA", "Constructs", "Download"],
                       "What you make from them: codon-optimised DNA, fusion "
                       "constructs, CSV, figures and the PDF report."),
    }

    chosen = st.radio("Section", list(FAMILIES), horizontal=True,
                      label_visibility="collapsed", key="analysis_family")
    family, tab_labels, family_blurb = FAMILIES[chosen]
    st.caption(family_blurb)
    tabs = (dict(zip(tab_labels, st.tabs(tab_labels))) if tab_labels
            else {"Triage": st.container()})



    run = Run(
        sequences=sequences,
        df=df,
        annotations=annotations,
        errors=errors,
        kinds=kinds,
        heavy_names=heavy_names,
        vhh_names=vhh_names,
        other_names=other_names,
        heavy_sequences=heavy_sequences,
        heavy_annotations=heavy_annotations,
        _vh_only=_vh_only,
        ph=ph,
        window=window,
        method=method,
        tabs=tabs,
        chosen=chosen,
        heavy_only_notice=heavy_only_notice,
        type_of=type_of)

    # One dispatch rather than six module-level `if` statements: only the
    # selected family runs - the whole performance argument for a selector -
    # and it is one line here instead of six that could drift apart.
    {
        "triage": _triage,
        "identity": _identity,
        "risk": _risk,
        "compare": _compare,
        "repertoire": _repertoire,
        "outputs": _outputs,
    }[family](run)
