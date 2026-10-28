"""
sections.structure - the Structure screen.

What a 3D model says, as opposed to what the sequence suggests. The measurements
themselves live in `structure.py`; this module is only the screen that arranges
them. The cached wrappers it calls live in `compute.py`.
"""

from __future__ import annotations

import re

import streamlit as st
import streamlit.components.v1 as components

import structure as stx
import vhh_core as vc
from compute import (st_confidence, st_disulfides, st_exposure,
                     st_identity_map, st_interface, st_interface_summary,
                     st_load, st_patches, st_resnum_map)
from shell import FIT, current_sequences, footer

__version__ = "2026.10.28c"


def render(method: str) -> None:
    """Draw the Structure screen.

    ``method`` is the CDR annotation method from the sidebar. It is passed in
    rather than read from a global: this screen used to live in `app.py`'s
    module scope, where its own ``kinds`` assignment silently overwrote the
    analysis screen's variable of the same name. Nothing broke, because two
    screens never render in the same run - but a function scope removes the
    whole class of accident rather than the one instance of it.
    """
    st.title("🧊 Structure")
    st.markdown(
        "Every other section of this app reasons from sequence alone, because a "
        "sequence is what you always have. Several of the questions it answers are "
        "really structural questions being approximated, **which cysteines pair, "
        "whether a liability sits on the surface, where the hydrophobic area is**, "
        "and once a model exists the approximation can be replaced by a "
        "measurement. Upload a crystal structure or a predicted model.")

    up = st.file_uploader(
        "Structure (.pdb, .cif, .mmcif, optionally gzipped)",
        type=["pdb", "cif", "mmcif", "ent", "gz"], key="struct_uploader",
        help="Crystal structures and AlphaFold2/3 or ColabFold models are all "
             "accepted. A complex in one file lets the interface be computed.")

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
  read directly, pLDDT included. For AlphaFold3, use the `.cif` it writes.

This app does not fetch structures for you, deliberately: downloading on your
behalf would mean this server making requests that look like they came from you.
""")
        footer()
        st.stop()

    try:
        chains, kind_guess, kind_reason, n_models, ssbond = st_load(raw, fname)
    except stx.StructureError as exc:
        st.error(str(exc))
        if st.button("Clear and try another file", key="struct_clear_err"):
            for key in ("struct_raw", "struct_name"):
                st.session_state.pop(key, None)
            st.rerun()
        footer()
        st.stop()

    head = st.columns([3, 1])
    with head[0]:
        st.success(f"**{fname}**, {len(chains)} chain(s), "
                   f"{int(chains['length'].sum())} residues.")
    with head[1]:
        if st.button("Clear structure", key="struct_clear", **FIT):
            for key in ("struct_raw", "struct_name"):
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

    annotation = None
    try:
        kinds = vc.classify(target_seq)
        if kinds.get("full_analysis"):
            annotation, _ = vc.annotate(target_seq, method)
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
            liab = vc.find_liabilities(target_seq, annotation)
        except Exception as exc:
            st.error(f"Liabilities could not be computed: {exc}")
            liab = None
        if liab is not None:
            if not len(liab):
                st.success("No liability motifs were found in the sequence.")
            else:
                model = stx.load_structure(raw, fname)
                joined = stx.annotate_liabilities(liab, model, chain,
                                                  chain_seq, target_seq)
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
                # annotate_liabilities already carries the model's own residue
                # numbers, which is exactly what the viewer selects on.
                joined = stx.annotate_liabilities(
                    vc.find_liabilities(target_seq, annotation),
                    stx.load_structure(raw, fname), chain, chain_seq, target_seq)
                highlight = [int(n) for n, e in zip(joined["resnum"],
                                                    joined["exposure"])
                             if e == "exposed" and n is not None]
        except Exception as exc:
            st.caption(f"Nothing to highlight ({exc}).")

        try:
            html = stx.viewer_html(
                stx.load_structure(raw, fname),
                raw.decode("utf-8", errors="replace") if isinstance(raw, bytes)
                else raw,
                colour_by=colour_by, chain_id=chain,
                exposure=st_exposure(raw, fname, chain)
                if colour_by == "exposure" else None,
                confidence=st_confidence(raw, fname, chain)
                if colour_by == "confidence" else None,
                highlight=highlight or None, surface=surface)
            components.html(html, height=540, scrolling=False)
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
