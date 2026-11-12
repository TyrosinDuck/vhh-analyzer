"""
sections.structure - the Structure screen.

What a 3D model says, as opposed to what the sequence suggests. The measurements
themselves live in `structure.py`; this module is only the screen that arranges
them. The cached wrappers it calls live in `compute.py`.
"""

from __future__ import annotations

import posixpath
import re

import matplotlib.pyplot as plt
import streamlit as st

import af3
import structure as stx
import vhh_core as vc
from compute import (annotation_for, classify_cached, st_af3_full,
                     st_af3_index, st_af3_member, st_confidence, st_disulfides,
                     st_exposure, st_identity_map, st_interface,
                     st_interface_summary, st_liabilities, st_load, st_patches,
                     st_resnum_map, st_viewer_html)
from shell import FIT, Settings, current_sequences, embed_html, footer

__version__ = "2026.11.12a"


def render(settings: Settings) -> None:
    """Draw the Structure screen.

    ``settings.method`` is the CDR annotation method from the sidebar. It is
    passed in rather than read from a global: this screen used to live in
    `app.py`'s module scope, where its own ``kinds`` assignment silently
    overwrote the analysis screen's variable of the same name. Nothing broke,
    because two screens never render in the same run - but a function scope
    removes the whole class of accident rather than the one instance of it.
    """
    method = settings.method
    st.title("🧊 Structure")
    st.markdown(
        "Every other section of this app reasons from sequence alone, because a "
        "sequence is what you always have. Several of the questions it answers are "
        "really structural questions being approximated, **which cysteines pair, "
        "whether a liability sits on the surface, where the hydrophobic area is**, "
        "and once a model exists the approximation can be replaced by a "
        "measurement. Upload a crystal structure or a predicted model.")

    up = st.file_uploader(
        "Structure (.pdb, .cif, .mmcif, optionally gzipped), or the .zip of "
        "an AlphaFold3 prediction",
        type=["pdb", "cif", "mmcif", "ent", "gz", "zip"], key="struct_uploader",
        help="Crystal structures and AlphaFold2/3 or ColabFold models are all "
             "accepted. A complex in one file lets the interface be computed. "
             "For AlphaFold3, upload the .zip exactly as AlphaFold Server (or "
             "your own run) gives it: ipTM and PAE live in JSON files beside "
             "the .cif, not in it, and are only shown when the .zip is "
             "uploaded.")

    # The uploader is unmounted when you leave the section, which discards what
    # it held. Parking the bytes under a key of our own is what lets you come
    # back to the structure you were looking at.
    if up is not None:
        st.session_state["struct_raw"] = up.getvalue()
        st.session_state["struct_name"] = up.name
    raw = st.session_state.get("struct_raw")
    fname = st.session_state.get("struct_name", "structure")

    if not raw:
        st.info("No structure loaded yet.")
        with st.expander("Where to get one"):
            st.markdown("""
- **RCSB PDB**: `https://www.rcsb.org/`, *Download Files → PDB Format*. For a
  VHH–antigen complex, download the whole entry rather than one chain, otherwise
  there is no interface left to compute.
- **AlphaFold DB**: `https://alphafold.ebi.ac.uk/`, search a UniProt accession.
- **Your own predictions**: AlphaFold2, AlphaFold3 and ColabFold outputs are
  read directly, pLDDT included. For AlphaFold3, upload the **.zip** AlphaFold
  Server gives you, unchanged: the app then also shows ipTM, the scores of all
  five models and the PAE, which say whether the *arrangement* of the chains
  can be trusted. A `.cif` on its own gives pLDDT only. If the zip is over the
  upload limit, delete its `msas/` and `templates/` folders, which the app
  does not read.

This app does not fetch structures for you, deliberately: downloading on your
behalf would mean this server making requests that look like they came from you.
""")
        footer()
        st.stop()

    def fail(message: str) -> None:
        st.error(message)
        if st.button("Clear and try another file", key="struct_clear_err"):
            for key in ("struct_raw", "struct_name", "struct_af3_model"):
                st.session_state.pop(key, None)
            st.rerun()
        footer()
        st.stop()

    # An AlphaFold3 archive is kept whole, and the model picked from it stands
    # in for an uploaded .cif from here on: everything below reads `raw` and
    # `fname` and never needs to know the structure came out of a zip.
    bundle = sample = zip_raw = None
    if af3.is_zip(raw):
        zip_raw = raw
        try:
            bundle = st_af3_index(zip_raw, fname)
        except stx.StructureError as exc:
            fail(str(exc))
        labels = [s.label for s in bundle.samples]
        if st.session_state.get("struct_af3_model") not in labels:
            st.session_state["struct_af3_model"] = labels[0]
        pick_model = st.columns([2, 3])
        with pick_model[0]:
            chosen = st.selectbox(
                f"Model from this {bundle.source} job", labels,
                key="struct_af3_model",
                format_func=lambda label: next(
                    f"{s.label}, ranking score {s.ranking_score:.2f}"
                    if s.ranking_score is not None else s.label
                    for s in bundle.samples if s.label == label),
                help="Sorted by AlphaFold3's own ranking score, best first. "
                     "Everything on this screen describes the model picked "
                     "here.")
        with pick_model[1]:
            st.caption(f"**{bundle.job_name}**, {len(bundle.samples)} model(s) "
                       f"in the archive. Read in memory and dropped, like any "
                       f"upload: nothing is kept.")
        sample = next(s for s in bundle.samples if s.label == chosen)
        try:
            raw = st_af3_member(zip_raw, sample.cif_member)
        except stx.StructureError as exc:
            fail(str(exc))
        fname = posixpath.basename(sample.cif_member)

    try:
        chains, kind_guess, kind_reason, n_models, ssbond = st_load(raw, fname)
    except stx.StructureError as exc:
        fail(str(exc))

    head = st.columns([3, 1])
    with head[0]:
        st.success(f"**{fname}**, {len(chains)} chain(s), "
                   f"{int(chains['length'].sum())} residues.")
    with head[1]:
        if st.button("Clear structure", key="struct_clear", **FIT):
            for key in ("struct_raw", "struct_name", "struct_af3_model"):
                st.session_state.pop(key, None)
            st.rerun()

    if n_models > 1:
        st.info(f"The file holds {n_models} models (an NMR ensemble, or several "
                f"predictions). Only the first is used, everything below "
                f"describes model 1.")

    st.dataframe(chains, hide_index=True, **FIT,
                 column_config={
                     "chain": "Chain", "length": "Residues",
                     "first": "First", "last": "Last",
                     "numbering_gaps": st.column_config.NumberColumn(
                         "Gaps", help="Breaks in the residue numbering, "
                                      "unmodelled loops, usually."),
                     "mean_b": st.column_config.NumberColumn(
                         "Mean B / pLDDT",
                         help="Mean of the B-factor column, whatever it holds."),
                     "sequence": st.column_config.TextColumn(
                         "Sequence", width="medium")})

    # ---- which chain is the molecule of interest, and what it is bound to ----
    pick = st.columns([1, 2, 2])
    chain_ids = list(chains["chain"])
    with pick[0]:
        chain = st.selectbox("Chain of interest", chain_ids, index=0,
                             key="struct_chain",
                             help="The VHH, usually. Everything per-residue "
                                  "below describes this chain.")
    partners_available = [c for c in chain_ids if c != chain]
    with pick[1]:
        partners = st.multiselect(
            "Binding partner(s)", partners_available,
            default=partners_available[:1] if partners_available else [],
            key="struct_partners",
            help="Chains to treat as the antigen when computing the interface.")
    chain_seq = chains.loc[chains["chain"] == chain, "sequence"].iloc[0]

    # ---- tie the structure to a sequence the rest of the app knows -----------
    loaded = current_sequences()
    with pick[2]:
        link_options = ["(the structure's own sequence)"] + list(loaded)
        link_choice = st.selectbox(
            "Annotate against", link_options, index=0, key="struct_link",
            help="A model is routinely missing termini and loops. Picking one of "
                 "the sequences loaded in the analysis section maps the model "
                 "onto the full-length construct, so positions match.")

    # Keyed by residue number throughout, see structure.resnum_map for why.
    if link_choice == link_options[0]:
        target_seq, target_name = chain_seq, f"{fname}:{chain}"
        offset_map = st_identity_map(raw, fname, chain)
        map_identity = 1.0
    else:
        target_name = link_choice
        target_seq = vc.clean_sequence(loaded[link_choice])
        offset_map, map_identity = st_resnum_map(raw, fname, chain, target_seq)
        covered = len(offset_map)
        message = (f"Chain {chain} ({len(chain_seq)} residues) mapped onto "
                   f"**{target_name}** ({len(target_seq)} residues): "
                   f"{covered} residues aligned at {map_identity:.1%} identity.")
        if map_identity < 0.90:
            st.warning(message + " That is low, check this is the same molecule. "
                                 "Positions below may be wrong.")
        else:
            st.caption(message)

    # Cached, both of them: under "auto" these run ANARCI, which shells out to
    # HMMER, and this block runs on every rerun of the screen, which means
    # every widget the user touches on it.
    annotation = None
    try:
        kinds = classify_cached(target_seq, method)
        if kinds.get("full_analysis"):
            annotation, _ = annotation_for(target_seq, method)
    except Exception:
        annotation = None

    views = ["🧬 Chains & confidence", "🔗 Disulfides", "💧 Exposure",
             "⚠️ Liabilities in 3D", "🫧 Surface patches"]
    if partners:
        views.append("🎯 Interface")
    views.append("👁 3D view")
    view = st.radio("View", views, horizontal=True, label_visibility="collapsed",
                    key="struct_view")
    st.divider()

    # ------------------------------------------------------------------ #
    if view == "🧬 Chains & confidence":
        if bundle is not None:
            _af3_panel(bundle, sample, zip_raw, chain, partners,
                       list(chains["chain"]))
            st.divider()
        st.subheader("What the B-factor column holds")
        st.markdown(
            "A predicted model stores **pLDDT** in the column a crystal structure "
            "uses for the **B-factor**. The two mean opposite things, high pLDDT "
            "is good, high B is bad, and no field in the file states which it is, "
            "so this is a guess you can override.")
        guess_label = {"pLDDT": "pLDDT (predicted model)",
                       "B-factor": "B-factor (experimental structure)",
                       "unknown": "could not tell"}[kind_guess]
        st.caption(f"Guess: **{guess_label}**, {kind_reason}")
        kind = st.radio("Read the column as",
                        ["pLDDT", "B-factor"],
                        index=0 if kind_guess == "pLDDT" else 1,
                        horizontal=True, key="struct_kind")
        conf = st_confidence(raw, fname, chain)
        if conf["value"].notna().any():
            st.pyplot(stx.plot_confidence(conf, kind,
                                          title=f"{kind} along chain {chain}"),
                      **FIT)
            if kind == "pLDDT":
                low = conf[conf["value"] < 70]
                st.markdown(
                    f"**{len(low)} of {len(conf)} residues sit below pLDDT 70**, "
                    f"the line below which the backbone should not be trusted"
                    + (f": {', '.join(f'{r.residue}{r.resnum}' for r in low.head(20).itertuples())}"
                       f"{'…' if len(low) > 20 else ''}." if len(low) else "."))
                st.caption(
                    "pLDDT is the model's confidence in itself. It is neither "
                    "accuracy nor a measurement of flexibility, although low "
                    "pLDDT and genuine disorder do tend to coincide. Side-chain "
                    "placement is less reliable than the backbone even at high "
                    "pLDDT, which is worth remembering when reading the solvent "
                    "accessibility of a long side chain.")
        else:
            st.info("The file carries no B-factor column.")

        if any(int(g) for g in chains["numbering_gaps"]):
            st.warning(
                "There are breaks in the residue numbering, unmodelled loops. "
                "Residues that are not in the model have no accessibility and no "
                "exposure call; they are reported as *not modelled* rather than "
                "as buried.")

    # ------------------------------------------------------------------ #
    elif view == "🔗 Disulfides":
        st.subheader("Measured, not inferred")
        st.markdown(
            "The analysis section decides which cysteines pair from their IMGT "
            "positions and a canonical-pair assumption. Here the distances "
            f"between SG atoms are measured: a pair under **{stx.SS_MAX_DISTANCE} "
            "Å** is bonded, and anything else is a free thiol.")
        found, free, n_sg, notes = st_disulfides(raw, fname)
        cols = st.columns(3)
        cols[0].metric("Cysteines with an SG atom", n_sg)
        cols[1].metric("Disulfide bonds", len(found))
        cols[2].metric("Free thiols", len(free))

        if len(found):
            st.dataframe(found, hide_index=True, **FIT,
                         column_config={
                             "chain_1": "Chain", "resnum_1": "Residue",
                             "chain_2": "Chain", "resnum_2": "Residue",
                             "distance": st.column_config.NumberColumn(
                                 "SG–SG (Å)", format="%.2f"),
                             "inter_chain": "Between chains"})
        else:
            st.info("No SG–SG pair is close enough to be a disulfide bond.")

        if free:
            st.warning("Free thiols: " + ", ".join(
                f"{f['chain']}{f['resnum']}" for f in free)
                + ". A free cysteine on the surface is a handle for conjugation "
                  "and a route to disulfide scrambling, depending on what you "
                  "want from it.")
        for note in notes:
            st.warning(note)
        if ssbond and not notes:
            st.success(f"The {len(ssbond)} SSBOND record(s) in the file agree with "
                       f"the measured distances.")

        # The comparison that makes this worth doing
        if annotation is not None:
            st.divider()
            st.subheader("Against the sequence-based prediction")
            try:
                predicted = vc.disulfides(target_seq, annotation)
            except Exception as exc:
                st.info(f"The sequence-based prediction is unavailable: {exc}")
            else:
                st.markdown(
                    f"Sequence-based: **{predicted.get('n_pairs', '?')}** pair(s) "
                    f"predicted, **{predicted.get('n_orphan', '?')}** orphan "
                    f"thiol(s). Measured here: **{len(found)}** bond(s), "
                    f"**{len(free)}** free thiol(s).")
                if predicted.get("n_pairs") == len(found):
                    st.success("The counts agree.")
                else:
                    st.warning(
                        "The counts disagree. The measurement is the one to "
                        "believe, but check first that the model covers every "
                        "cysteine in the sequence: a cysteine in an unmodelled "
                        "loop cannot be seen to pair.")

    # ------------------------------------------------------------------ #
    elif view == "💧 Exposure":
        st.subheader("Relative solvent accessibility")
        st.markdown(
            "Absolute accessibility in Å², divided by the maximum that residue "
            "can have, so a Trp and a Gly are comparable. At or below "
            f"**{stx.BURIED_MAX:.0%}** of maximum is called buried, at or above "
            f"**{stx.EXPOSED_MIN:.0%}** exposed, and the band between is left as "
            "intermediate rather than forced into one of the two.")
        exposure = st_exposure(raw, fname, chain)
        counts = exposure["exposure"].value_counts()
        cols = st.columns(4)
        cols[0].metric("Buried", int(counts.get("buried", 0)))
        cols[1].metric("Intermediate", int(counts.get("intermediate", 0)))
        cols[2].metric("Exposed", int(counts.get("exposed", 0)))
        cols[3].metric("Total SASA", f"{exposure['sasa'].sum():.0f} Å²")

        st.pyplot(stx.plot_exposure(
            exposure, annotation, offset_map,
            title=f"Relative solvent accessibility, chain {chain}"), **FIT)

        if partners:
            st.caption(
                f"Computed with chain(s) {', '.join(partners)} present, so these "
                f"are the accessibilities the residues actually have in the "
                f"complex. The interface view reports what they lose on binding.")

        show = exposure.copy()
        if annotation and offset_map:
            show["region"] = stx.residue_regions(show["resnum"], annotation,
                                                 offset_map)
        st.dataframe(show, hide_index=True, **FIT, height=320,
                     column_config={
                         "resnum": "Residue no.", "residue": "AA",
                         "sasa": st.column_config.NumberColumn("SASA (Å²)",
                                                               format="%.1f"),
                         "rel_sasa": st.column_config.ProgressColumn(
                             "Relative", min_value=0.0, max_value=1.0,
                             format="%.2f"),
                         "exposure": "Call", "region": "Region"})
        st.download_button("Download accessibility as CSV",
                           show.to_csv(index=False).encode(),
                           file_name=f"{fname}_{chain}_accessibility.csv",
                           mime="text/csv", key="struct_dl_exposure")
        st.caption(
            "Shrake–Rupley with a 1.4 Å probe, normalised by the residue maxima "
            "of Tien et al. 2013 (PLoS ONE 8:e80635). Those maxima were derived "
            "with DSSP, i.e. the Lee–Richards definition; Shrake–Rupley is a "
            "different numerical approximation of the same quantity, so the "
            "buried/exposed call is robust but a relative accessibility quoted "
            "to the percent is not. Waters are removed before the calculation; "
            "ligands are kept, because a bound ligand really does cover the "
            "surface underneath it.")

    # ------------------------------------------------------------------ #
    elif view == "⚠️ Liabilities in 3D":
        st.subheader("Which of the flagged liabilities actually matter")
        st.markdown(
            "This is the join worth uploading a structure for. The analysis "
            "section flags a deamidation, oxidation or glycosylation motif "
            "wherever it occurs in the sequence; it cannot know whether the motif "
            "is on the surface or packed into the core. A buried motif is "
            "largely inert. An exposed one in a CDR is a real problem.")
        try:
            joined = st_liabilities(raw, fname, chain, chain_seq, target_seq,
                                    annotation)
        except Exception as exc:
            st.error(f"Liabilities could not be computed: {exc}")
            joined = None
        if joined is not None:
            if not len(joined):
                st.success("No liability motifs were found in the sequence.")
            else:
                order = {"exposed": 0, "intermediate": 1,
                         "buried": 2, "not modelled": 3}
                joined = joined.assign(
                    _rank=[order.get(e, 4) for e in joined["exposure"]]
                ).sort_values(["_rank", "score"],
                              ascending=[True, False]).drop(columns="_rank")
                counts = joined["exposure"].value_counts()
                cols = st.columns(4)
                cols[0].metric("Exposed", int(counts.get("exposed", 0)))
                cols[1].metric("Intermediate", int(counts.get("intermediate", 0)))
                cols[2].metric("Buried", int(counts.get("buried", 0)),
                               help="Largely inert, whatever the sequence-based "
                                    "risk column says.")
                cols[3].metric("Not modelled", int(counts.get("not modelled", 0)),
                               help="In a loop or terminus the model does not "
                                    "cover, so nothing can be said.")
                st.dataframe(joined, hide_index=True, **FIT,
                             column_config={
                                 "rel_sasa": st.column_config.ProgressColumn(
                                     "Relative accessibility", min_value=0.0,
                                     max_value=1.0, format="%.2f"),
                                 "exposure": "Measured"})
                st.download_button(
                    "Download as CSV", joined.to_csv(index=False).encode(),
                    file_name=f"{fname}_{chain}_liabilities_exposure.csv",
                    mime="text/csv", key="struct_dl_liab")
                st.info(
                    "The `risk` column is the sequence-based one and is left "
                    "untouched on purpose: downgrading it automatically because a "
                    "residue looks buried in one static conformer would be "
                    "overclaiming. Read the two columns together.")

    # ------------------------------------------------------------------ #
    elif view == "🫧 Surface patches":
        st.subheader("Hydrophobic surface, grouped in space")
        st.markdown(
            "What drives aggregation is hydrophobic area that is **contiguous on "
            "the surface**, which a sliding window over the sequence cannot see: "
            "two residues far apart in sequence can sit side by side on the fold. "
            "This groups solvent-exposed "
            f"**{''.join(sorted(stx.PATCH_RESIDUES))}** residues whose side "
            "chains come within the contact distance of each other.")
        ctrl = st.columns(3)
        with ctrl[0]:
            min_rel = st.slider("Minimum relative accessibility", 0.05, 0.60,
                                stx.PATCH_MIN_REL, 0.05, key="struct_patch_rel",
                                help="How exposed a residue must be to count. "
                                     "Lower than you might expect on purpose: a "
                                     "residue inside a patch is half-covered by "
                                     "its own neighbours.")
        with ctrl[1]:
            contact = st.slider("Contact distance (Å)", 3.5, 8.0,
                                stx.PATCH_CONTACT, 0.5, key="struct_patch_d",
                                help="Between closest side-chain heavy atoms, "
                                     "not between centroids, Trp and Arg carry "
                                     "their bulk far from the backbone.")
        with ctrl[2]:
            min_size = st.slider("Minimum patch size", 2, 8,
                                 stx.PATCH_MIN_SIZE, 1, key="struct_patch_n")
        patches = st_patches(raw, fname, chain, min_rel, contact, min_size)
        if not len(patches):
            st.success(
                f"No group of {min_size} or more exposed hydrophobic residues. "
                f"Loosen the thresholds to see smaller or less exposed clusters.")
        else:
            st.dataframe(patches, hide_index=True, **FIT,
                         column_config={
                             "patch": "#", "n_residues": "Residues",
                             "residues": st.column_config.TextColumn(
                                 "Which", width="large"),
                             "exposed_area": st.column_config.NumberColumn(
                                 "Exposed area (Å²)", format="%.0f")})
            biggest = patches.iloc[0]
            st.caption(
                f"The largest patch exposes {biggest['exposed_area']:.0f} Å² "
                f"across {biggest['n_residues']} residues. Switch to the 3D view "
                f"and highlight them to see whether they form one face or are "
                f"strung along a groove.")
        st.warning(
            "**A transparent composite, not SAP or CamSol.** Those methods are "
            "parameterised against experimental aggregation data this app does "
            "not have. This asks a simpler and fully stated question, which "
            "strongly hydrophobic residues are both exposed and adjacent, and "
            "the three thresholds above are yours to move. Read it as *here is a "
            "hydrophobic area worth looking at*, never as a predicted "
            "aggregation rate.")

    # ------------------------------------------------------------------ #
    elif view == "🎯 Interface":
        st.subheader(f"Chain {chain} against {', '.join(partners)}")
        st.markdown(
            "A residue is at the interface if it loses solvent accessibility when "
            "the partner is present. ΔSASA, accessible alone minus accessible in "
            "the complex, is the standard operational definition, and unlike a "
            "distance cutoff it weights each residue by how much surface it "
            "actually commits to the contact.")
        if bundle is not None:
            _af3_interface_note(bundle, sample, zip_raw, chain, partners,
                                list(chains["chain"]))
        iface = st_interface(raw, fname, chain, tuple(partners))
        if not len(iface):
            st.info("No residue of this chain loses accessibility, the chains "
                    "are not in contact in this model.")
        else:
            if len(partners) == 1:
                summary = st_interface_summary(raw, fname, chain, partners[0])
                cols = st.columns(4)
                cols[0].metric("Buried per side",
                               f"{summary['bsa_per_side']:.0f} Å²",
                               help="Half the total buried surface, the number "
                                    "usually quoted for an interface.")
                cols[1].metric("Total buried", f"{summary['bsa_total']:.0f} Å²")
                cols[2].metric(f"Residues on {chain}", summary["n_residues_a"])
                cols[3].metric(f"Residues on {partners[0]}",
                               summary["n_residues_b"])
                st.caption(
                    "For reference, an antibody–antigen interface typically "
                    "buries 600–900 Å² per side; a VHH, with one domain instead "
                    "of two, often sits at the lower end or below it.")

            st.pyplot(stx.plot_interface(
                iface, st_exposure(raw, fname, chain), annotation, offset_map,
                title=f"Surface of chain {chain} buried on binding"), **FIT)

            show = iface.copy()
            if annotation and offset_map:
                show["region"] = stx.residue_regions(show["resnum"], annotation,
                                                     offset_map)
                by_region = (show.groupby("region")["delta_sasa"].sum()
                             .sort_values(ascending=False))
                st.markdown("**Buried surface by region**, which parts of the "
                            "domain actually do the binding:")
                st.dataframe(
                    by_region.rename("ΔSASA (Å²)").reset_index(),
                    hide_index=True, **FIT)
                cdr_share = sum(v for k, v in by_region.items()
                                if str(k).startswith("CDR"))
                total = float(by_region.sum()) or 1.0
                st.caption(
                    f"The CDRs account for {cdr_share / total:.0%} of the buried "
                    f"surface. A framework contribution is normal for a VHH, the "
                    f"FR2 face that would be the VL interface in a conventional "
                    f"antibody is solvent-exposed and often takes part in binding.")

            st.dataframe(show.sort_values("delta_sasa", ascending=False),
                         hide_index=True, **FIT, height=320,
                         column_config={
                             "resnum": "Residue no.", "residue": "AA",
                             "sasa_free": st.column_config.NumberColumn(
                                 "Free (Å²)", format="%.1f"),
                             "sasa_bound": st.column_config.NumberColumn(
                                 "Bound (Å²)", format="%.1f"),
                             "delta_sasa": st.column_config.NumberColumn(
                                 "ΔSASA (Å²)", format="%.1f"),
                             "fraction_buried": st.column_config.ProgressColumn(
                                 "Fraction buried", min_value=0.0, max_value=1.0,
                                 format="%.2f"),
                             "rel_sasa_free": None, "region": "Region"})
            st.download_button(
                "Download interface as CSV", show.to_csv(index=False).encode(),
                file_name=f"{fname}_{chain}_interface.csv", mime="text/csv",
                key="struct_dl_iface")

            if len(partners) == 1:
                st.divider()
                st.subheader(f"Epitope, the {partners[0]} side")
                epitope = st_interface(raw, fname, partners[0], (chain,))
                st.markdown(
                    f"The same calculation run the other way round: residues of "
                    f"chain {partners[0]} that chain {chain} covers. This is the "
                    f"**structural epitope**, as far as one model can define it.")
                st.dataframe(epitope.sort_values("delta_sasa", ascending=False),
                             hide_index=True, **FIT, height=280)
                st.download_button(
                    "Download epitope as CSV",
                    epitope.to_csv(index=False).encode(),
                    file_name=f"{fname}_{partners[0]}_epitope.csv",
                    mime="text/csv", key="struct_dl_epitope")
                st.caption(
                    "A structural epitope from one model, not a functional one: "
                    "it says which residues are covered, not which ones the "
                    "binding energy depends on. Those are different sets, and "
                    "telling them apart needs mutagenesis.")

    # ------------------------------------------------------------------ #
    elif view == "👁 3D view":
        st.subheader("Interactive view")
        opts = st.columns([2, 2, 1])
        with opts[0]:
            colour_by = st.selectbox(
                "Colour by",
                ["chain", "spectrum", "confidence", "exposure"], index=0,
                key="struct_colour",
                format_func=lambda v: {
                    "chain": "chain", "spectrum": "N → C rainbow",
                    "confidence": "pLDDT / B-factor",
                    "exposure": "buried → exposed"}[v])
        with opts[1]:
            highlight_what = st.selectbox(
                "Highlight as sticks",
                ["nothing", "disulfides", "free thiols", "largest patch",
                 "exposed liabilities"], index=0, key="struct_highlight")
        with opts[2]:
            surface = st.checkbox("Surface", value=False, key="struct_surface")

        highlight: list[int] = []
        try:
            if highlight_what == "disulfides":
                found, _free, _n, _notes = st_disulfides(raw, fname)
                highlight = [int(n) for n in
                             list(found["resnum_1"]) + list(found["resnum_2"])]
            elif highlight_what == "free thiols":
                _f, free, _n, _notes = st_disulfides(raw, fname)
                highlight = [int(f["resnum"]) for f in free
                             if f["chain"] == chain]
            elif highlight_what == "largest patch":
                patches = st_patches(raw, fname, chain, stx.PATCH_MIN_REL,
                                     stx.PATCH_CONTACT, stx.PATCH_MIN_SIZE)
                if len(patches):
                    highlight = [int(re.sub(r"\D", "", token))
                                 for token in patches.iloc[0]["residues"].split(", ")]
            elif highlight_what == "exposed liabilities":
                # st_liabilities already carries the model's own residue
                # numbers, which is exactly what the viewer selects on, and it
                # is the same cached table the liability view draws.
                joined = st_liabilities(raw, fname, chain, chain_seq,
                                        target_seq, annotation)
                highlight = [int(n) for n, e in zip(joined["resnum"],
                                                    joined["exposure"])
                             if e == "exposed" and n is not None]
        except Exception as exc:
            st.caption(f"Nothing to highlight ({exc}).")

        try:
            html = st_viewer_html(raw, fname, chain, colour_by,
                                  tuple(highlight), surface)
            embed_html(html, height=540)
        except Exception as exc:
            st.error(f"The 3D view could not be rendered: {exc}")

        if colour_by == "confidence":
            st.caption("AlphaFold's own colours: dark blue above 90, light blue "
                       "70–90, yellow 50–70, orange below 50.")
        elif colour_by == "exposure":
            st.caption("Dark grey buried, pale blue intermediate, orange exposed.")
        st.caption(
            "Rendered by 3Dmol.js, loaded from a CDN. If the panel stays blank, "
            "your browser could not reach it, every table and figure in this "
            "section is computed server-side with Biopython and is unaffected.")

    st.divider()
    st.caption(
        "One static conformer of one model. This section measures; it does not "
        "predict structure, dock, estimate affinity or model flexibility. For a "
        "predicted model, treat side-chain accessibility as indicative.")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.


# --------------------------------------------------------------------------- #
# AlphaFold3 confidence, shown only when the .zip was uploaded
# --------------------------------------------------------------------------- #

def _num(value, digits: int = 2) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def _pair_rows(sample, zip_raw, chain, partners, cif_chains):
    summary = sample.summary
    full = st_af3_full(zip_raw, sample.full_member)
    order, order_source = af3.chain_order(summary, full, cif_chains)
    rows = []
    for partner in partners:
        rows.append({
            "partner": partner,
            "iptm": af3.pair_iptm(summary, order, chain, partner),
            "pae_min": af3.pair_pae_min(summary, order, chain, partner),
            "block": af3.interface_pae(full, chain, partner) if full else None,
        })
    return rows, full, order_source


def _af3_panel(bundle, sample, zip_raw, chain, partners, cif_chains) -> None:
    st.subheader("AlphaFold3: can the arrangement of the chains be trusted?")
    st.markdown(
        "pLDDT, further down, says whether each chain is folded right. It says "
        "nothing about where the chains sit **relative to each other**: two "
        "chains can both be confidently folded and docked on the wrong face. "
        "That is what **ipTM** and the **PAE** are for. AlphaFold3 writes them "
        "beside the `.cif`, not in it, which is why they appear only when the "
        "`.zip` is uploaded.")
    summary = sample.summary
    if not summary:
        st.warning(f"{sample.label} has no summary_confidences file in the "
                   f"archive, so there is no ipTM to show for it.")
        return

    cols = st.columns(4)
    cols[0].metric("ipTM, whole complex", _num(summary.get("iptm")),
                   help="Every chain against every other. With more than two "
                        "chains it mixes all the pairs, so the pair value "
                        "below is the one to read.")
    cols[1].metric("pTM", _num(summary.get("ptm")),
                   help="The whole complex, folding and arrangement together.")
    cols[2].metric("Ranking score", _num(summary.get("ranking_score")),
                   help="What AlphaFold3 sorts its models by: 0.8 ipTM + "
                        "0.2 pTM + 0.5 fraction disordered, minus 100 if "
                        "there is a clash.")
    cols[3].metric("Clash", "yes" if summary.get("has_clash") else "no",
                   help="Chains interpenetrating. A model with a clash is "
                        "ranked last whatever its other scores.")
    if summary.get("has_clash"):
        st.error("AlphaFold3 flags a clash in this model: chains overlap in "
                 "space. Nothing measured on it at the interface is meaningful.")

    if not partners:
        st.info("Pick a binding partner above to see the ipTM of that pair, "
                "which is the number that matters for a VHH and its antigen.")
    else:
        rows, full, order_source = _pair_rows(sample, zip_raw, chain, partners,
                                              cif_chains)
        for row in rows:
            level, sentence = af3.iptm_verdict(row["iptm"])
            text = (f"**Chain {chain} against {row['partner']}: ipTM "
                    f"{_num(row['iptm'])}**, {sentence}")
            if row["pae_min"] is not None:
                text += (f" Lowest PAE between the two chains: "
                         f"{row['pae_min']:.1f} Å.")
            {"high": st.success, "grey": st.warning, "low": st.error,
             "none": st.info}[level](text)
        st.caption(
            f"Chain order for the pair values read from {order_source}. "
            f"Thresholds from AlphaFold3: above {af3.IPTM_HIGH} confident, "
            f"below {af3.IPTM_LOW} likely failed. ipTM is a confidence, not an "
            f"accuracy and not an affinity. It is pulled down by disordered "
            f"segments, a His-tag, a linker, a floppy terminus, even when the "
            f"interface itself is right, so a construct with tags is worth "
            f"predicting without them.")

    st.markdown("**All the models in the archive.** Models that agree on a "
                "high ipTM are reassuring. Models that disagree, or agree on a "
                "low one, mean the interface is not settled.")
    table = af3.samples_table(bundle, chain, partners[0] if partners else None,
                              cif_chains)
    config = {
        "model": "Model",
        "ranking_score": st.column_config.NumberColumn("Ranking", format="%.2f"),
        "iptm": st.column_config.NumberColumn("ipTM", format="%.2f"),
        "ptm": st.column_config.NumberColumn("pTM", format="%.2f"),
        "has_clash": "Clash",
        "fraction_disordered": st.column_config.NumberColumn(
            "Disordered", format="%.2f",
            help="Fraction of the structure AlphaFold3 considers disordered.")}
    if partners:
        config["pair_iptm"] = st.column_config.NumberColumn(
            f"ipTM {chain}–{partners[0]}", format="%.2f")
        config["pair_pae_min"] = st.column_config.NumberColumn(
            f"Min PAE {chain}–{partners[0]} (Å)", format="%.1f")
    st.dataframe(table, hide_index=True, column_config=config, **FIT)

    full = st_af3_full(zip_raw, sample.full_member)
    if full is None:
        st.caption("No PAE file for this model in the archive.")
        return
    st.markdown("**Predicted aligned error**")
    left, right = st.columns([3, 2])
    with left:
        st.pyplot(af3.plot_pae(full, title=f"PAE, {sample.label}"), **FIT)
    with right:
        st.markdown(
            "Each cell is the error AlphaFold3 expects in the position of the "
            "column residue when the model is superposed on the row residue. "
            "**Dark is confident.** The blocks on the diagonal are each chain "
            "against itself. The **off-diagonal blocks** are the arrangement: "
            "pale throughout means the chains' relative placement is a guess, "
            "and a dark stripe means some region of one chain, an epitope "
            "perhaps, is placed with confidence relative to the other.")
        if partners:
            rows, _full, _source = _pair_rows(sample, zip_raw, chain,
                                              partners, cif_chains)
            for row in rows:
                block = row["block"]
                if block:
                    st.caption(
                        f"{chain}–{row['partner']}: mean {block['mean']:.1f} Å, "
                        f"min {block['min']:.1f} Å, "
                        f"{block['fraction_close']:.0%} of residue pairs under "
                        f"{block['close']:.0f} Å. The mean is dominated by the "
                        f"parts of each chain far from the other; the minimum "
                        f"and the fraction are what speak to an interface.")
    plt.close("all")


def _af3_interface_note(bundle, sample, zip_raw, chain, partners,
                        cif_chains) -> None:
    rows, _full, _source = _pair_rows(sample, zip_raw, chain, partners,
                                      cif_chains)
    for row in rows:
        level, _sentence = af3.iptm_verdict(row["iptm"])
        if level == "low":
            st.error(f"AlphaFold3 gives chain {chain} against {row['partner']} "
                     f"an ipTM of {_num(row['iptm'])}: the arrangement has "
                     f"likely failed. The interface below is what this model "
                     f"happens to contain, not evidence of an epitope.")
        elif level == "grey":
            st.warning(f"ipTM {chain}–{row['partner']} is {_num(row['iptm'])}, "
                       f"in the grey zone. Read the interface below as a "
                       f"hypothesis, and check whether the other models in the "
                       f"archive put the VHH in the same place.")
        elif level == "high":
            st.caption(f"ipTM {chain}–{row['partner']}: {_num(row['iptm'])}, "
                       f"a confident arrangement.")
