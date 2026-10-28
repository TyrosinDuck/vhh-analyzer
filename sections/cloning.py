"""
sections.cloning - Move a VHH from one vector into another, on paper.

The screen only. Its computations live in `compute.py`; the chrome it uses
(`FIT`, `footer`, navigation) lives in `shell.py`.
"""

from __future__ import annotations

import io
import re

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import cloning as cl
import sanger as sg
from compute import load_plasmid_cached, plan_cloning_cached, png_plasmid_map, primer_report_cached
from shell import FIT, footer

__version__ = "2026.10.28c"


def render() -> None:
    """Draw the screen."""
    st.title("🧬 In-silico cloning")
    st.markdown(
        "Pick the recipient vector, paste the donor DNA, and the app works out which "
        "restriction enzymes put the VHH between the signal peptide and the tags — "
        "then checks that the cloned VHH encodes exactly the protein you started "
        "from.")

    # ---- 1. recipient ----
    st.subheader("1 · Recipient vector")
    bundled = cl.list_bundled_plasmids()
    rcol = st.columns([2, 1])
    with rcol[0]:
        backbone_name = st.selectbox("Vector", list(bundled) + ["Upload my own…"],
                                     index=0)
    with rcol[1]:
        layout = st.selectbox(
            "Expected tags", list(sg.VECTORS), index=0,
            help="Which tags the finished construct should carry — used when the "
                 "cassette is translated back.")

    payload, display_name = None, backbone_name
    if backbone_name == "Upload my own…":
        up = st.file_uploader("Vector (SnapGene .dna, GenBank or FASTA)",
                              type=["dna", "gb", "gbk", "genbank", "fasta", "fa",
                                    "txt"], key="backbone_upload")
        if up is not None:
            payload, display_name = up.getvalue(), up.name
    else:
        payload = bundled[backbone_name].read_bytes()

    if payload is None:
        st.info("Upload a vector to continue.")
        footer()
        st.stop()

    try:
        plasmid, features = load_plasmid_cached(payload, display_name)
    except Exception as exc:
        st.error(f"Could not read this vector — {exc}")
        footer()
        st.stop()

    if plasmid["warning"]:
        st.warning(plasmid["warning"])
    stats = cl.sequence_stats(plasmid["record"])
    auto_region = cl.insertion_region(plasmid, features)
    mcols = st.columns(4)
    mcols[0].metric("Length", f"{stats['length']} bp")
    mcols[1].metric("GC", f"{stats['gc_percent']} %")
    mcols[2].metric("Features", len(features))
    mcols[3].metric("Cloning window",
                    f"{auto_region[0]}–{auto_region[1]}" if auto_region else "—")
    if auto_region:
        st.caption(f"The insert goes between the end of the signal peptide and the "
                   f"first tag, positions {auto_region[0]}–{auto_region[1]}.")
    else:
        st.error("Could not find a leader and a tag in this vector, so there is no "
                 "window to clone into. Check the annotation.")

    # ---- 2. donor ----
    st.subheader("2 · Donor DNA")
    donor_kind = st.radio(
        "What are you pasting",
        ["A plasmid that carries the VHH", "A VHH fragment amplified with primers"],
        horizontal=True)
    donor_is_plasmid = donor_kind.startswith("A plasmid")

    if donor_is_plasmid:
        st.caption("The app locates the VHH, picks the primer pair that brackets it, "
                   "and uses the amplicon as the insert.")
    else:
        st.caption("The fragment must already carry the restriction sites its primers "
                   "added — a bare coding sequence cannot be cloned by restriction.")

    dfile = st.file_uploader(
        "Optional: load a file instead of pasting (SnapGene .dna, GenBank, FASTA)",
        type=["dna", "gb", "gbk", "genbank", "fasta", "fa", "txt"],
        key="donor_file")
    if dfile is not None and dfile.file_id != st.session_state.get("last_donor_file"):
        st.session_state["last_donor_file"] = dfile.file_id
        try:
            loaded = cl.load_plasmid(dfile.getvalue(), name=dfile.name)
            st.session_state["pending_donor_dna"] = str(loaded["record"].seq)
            st.rerun()
        except Exception as exc:
            st.error(f"Could not read that file — {exc}")

    if "pending_donor_dna" in st.session_state:
        st.session_state["donor_dna_box"] = st.session_state.pop("pending_donor_dna")
    donor_text = st.text_area("Donor DNA (5'→3')", height=150, key="donor_dna_box",
                              placeholder="ATGAAATACCTATTGCCTACGGCAGCCGCT…")
    donor_dna = re.sub(r"[^ACGTacgt]", "", donor_text or "").upper()
    if donor_dna:
        st.caption(f"{len(donor_dna)} nt")

    if not donor_dna:
        st.info("Paste the donor DNA, or load a file, to plan the cloning.")
        footer()
        st.stop()

    # ---- 3. the plan ----
    st.subheader("3 · Plan")
    own = st.file_uploader("Optional: your own primer list (CSV with name,sequence)",
                           type=["csv"], key="primer_upload")
    if own is not None:
        primers = pd.read_csv(io.BytesIO(own.getvalue()))
        primers["sequence"] = primers["sequence"].astype(str).str.upper()
        if "note" not in primers.columns:
            primers["note"] = ""
    else:
        primers = cl.load_primers()
    primer_rows = tuple((r["name"], r["sequence"], str(r.get("note", "")))
                        for _, r in primers.iterrows())

    with st.spinner("Locating the VHH, amplifying, cutting and translating…"):
        plan = plan_cloning_cached(payload, display_name, donor_dna,
                                   donor_is_plasmid, primer_rows, layout)

    if not plan["donor_vhh"]:
        st.error("No VHH could be found in that DNA. " +
                 " ".join(plan["messages"]))
        footer()
        st.stop()

    st.success(f"VHH found in the donor — {len(plan['donor_vhh'])} aa.")
    with st.expander("Donor VHH sequence"):
        st.code("\n".join(plan["donor_vhh"][i:i + 60]
                          for i in range(0, len(plan["donor_vhh"]), 60)),
                language="text")

    if plan.get("primer_pair"):
        f = plan["primer_pair"]["forward"]
        r = plan["primer_pair"]["reverse"]
        st.markdown(
            f"**Primers chosen automatically** — `{f['name']}` as the forward primer "
            f"({f['form']}, position {f['start']}) and `{r['name']}` as the reverse "
            f"({r['form']}, position {r['start']}), giving a "
            f"**{len(plan['amplicon'])} bp** amplicon.")
        if r["form"] == "reverse complement":
            st.caption(
                f"`{r['name']}` is listed as the plus-strand sequence of its target, "
                "so the oligo that primes back towards the VHH is its reverse "
                "complement — that is the form used here.")

    for msg in plan["messages"]:
        st.warning(msg)

    if not plan["plans"]:
        footer()
        st.stop()

    usable = [p for p in plan["plans"] if p["usable"]]
    if not usable:
        st.error("No enzyme pair gives a construct that preserves the donor's VHH. "
                 "The table below says why each one was rejected.")
    else:
        best = usable[0]
        st.markdown("### Recommended")
        bc = st.columns(4)
        bc[0].metric("Enzymes", f"{best['enzyme_5']} + {best['enzyme_3']}")
        bc[1].metric("Construct", f"{best['construct_bp']} bp")
        bc[2].metric("VHH", f"{len(best['construct_vhh'])} aa")
        bc[3].metric("Vector removed", f"{best['removed_bp']} bp")

        comparison = best["comparison"]
        if comparison["status"] == "identical":
            st.success("**Protein check passed** — the VHH encoded by the construct "
                       "is exactly the one in the donor.")
        else:
            st.info(f"**Protein check** — {comparison['message']}")

        pcols = st.columns(2)
        with pcols[0]:
            st.markdown("**VHH in the donor**")
            st.code("\n".join(plan["donor_vhh"][i:i + 40]
                              for i in range(0, len(plan["donor_vhh"]), 40)),
                    language="text")
        with pcols[1]:
            st.markdown("**VHH in the construct**")
            st.code("\n".join(best["construct_vhh"][i:i + 40]
                              for i in range(0, len(best["construct_vhh"]), 40)),
                    language="text")

        check = cl.verify_construct(best["construct"], layout)
        st.markdown(f"**Cassette** — signal peptide, then the VHH, then "
                    f"{check['tags_found'] or 'no tag'}; "
                    f"{check['premature_stops']} premature stop(s).")
        st.code("\n".join(check["cassette_aa"][i:i + 60]
                          for i in range(0, len(check["cassette_aa"]), 60)),
                language="text")

        built = cl.assemble(plasmid["record"], plan["amplicon"],
                            best["enzyme_5"], best["enzyme_3"],
                            circular=plasmid["topology"] == "circular")
        cname = f"{plasmid['name']}_VHH".replace(" ", "_")[:24]
        crec = cl.construct_record(best["construct"], cname,
                                   built.get("insert_start", 1),
                                   built.get("insert_end", 1), layout=layout,
                                   circular=plasmid["topology"] == "circular")
        dl = st.columns(3)
        with dl[0]:
            st.download_button("Construct (GenBank)", cl.record_to_genbank(crec),
                               file_name=f"{cname}.gb", mime="text/plain",
                               type="primary", **FIT)
        with dl[1]:
            st.download_button(
                "Construct (FASTA)",
                f">{cname}\n" + "\n".join(best["construct"][i:i + 60]
                                          for i in range(0, len(best["construct"]),
                                                         60)) + "\n",
                file_name=f"{cname}.fasta", mime="text/plain", **FIT)
        with dl[2]:
            if plan["amplicon"] != donor_dna:
                st.download_button(
                    "Amplicon (FASTA)",
                    ">amplicon\n" + "\n".join(
                        plan["amplicon"][i:i + 60]
                        for i in range(0, len(plan["amplicon"]), 60)) + "\n",
                    file_name="amplicon.fasta", mime="text/plain", **FIT)

        st.markdown(f"""
**Protocol** — digest the {len(plan['amplicon'])} bp insert and
**{plasmid['name']}** with **{best['enzyme_5']}** and **{best['enzyme_3']}**, then
ligate. The vector loses {best['removed_bp']} bp between the two sites and gains
{best['inserted_bp']} bp from the insert. Screen the colonies by sequencing — the
sequencing section of this app reads the .ab1 files back and confirms the cassette.
""")

    # ---- details, folded away ----
    with st.expander(f"All enzyme pairs considered ({len(plan['plans'])})"):
        rows = [{"5' enzyme": p["enzyme_5"], "3' enzyme": p["enzyme_3"],
                 "usable": p["usable"], "construct bp": p["construct_bp"],
                 "vector removed": p["removed_bp"],
                 "protein": p["comparison"]["status"],
                 "deletes": ", ".join(p["deletes"]) or "—",
                 "why not": ("" if p["usable"] else
                             (", ".join(p["deletes"]) + " would be deleted"
                              if p["deletes"] else p["comparison"]["message"]))}
                for p in plan["plans"]]
        st.dataframe(pd.DataFrame(rows), hide_index=True, **FIT,
                     column_config={"why not": st.column_config.TextColumn(
                         "Why not", width="large")})
        st.caption("Each pair is cut, joined and translated for real, then the VHH it "
                   "encodes is compared with the donor's. A pair that would delete an "
                   "annotated part of the vector is rejected even when the cassette "
                   "still translates.")

    with st.expander("Restriction sites"):
        cand = plan["candidates"]
        if len(cand):
            st.dataframe(
                cand[["enzyme", "site", "backbone_cuts", "position", "in_region",
                      "insert_cuts", "overhang", "verdict"]].head(30),
                hide_index=True, **FIT,
                column_config={"verdict": st.column_config.TextColumn(
                    "Verdict", width="large")})
            st.caption("An enzyme that cuts the insert is not necessarily wrong: in "
                       "VHH cloning BstEII is used on purpose, its site sitting in "
                       "FR4 at the 3' end of the domain.")

    with st.expander("Primers on this vector"):
        full_report = primer_report_cached(payload, display_name, primer_rows)
        binding = full_report[full_report["binds"]] if len(full_report) else full_report
        if len(binding):
            st.dataframe(binding[["primer", "sequence", "orientation", "start", "end",
                                  "tm"]], hide_index=True, **FIT,
                         column_config={"tm": st.column_config.NumberColumn(
                             "Tm (°C)", format="%.1f")})
        absent = [p for p in primers["name"] if p not in set(binding["primer"])] \
            if len(binding) else list(primers["name"])
        if absent:
            st.caption(f"No binding site in {plasmid['name']}: {', '.join(absent)}. "
                       "That is usually a property of the vector rather than a "
                       "mistake — pHEN6 has no gIII for a gIII primer to sit in.")

    with st.expander("Vector map"):
        view = st.radio("View", ["Cloning region", "Whole plasmid"], horizontal=True,
                        key="map_view")
        sites = (tuple((p["enzyme_5"], "", tuple([p["removed_bp"]])) for p in [])
                 or ())
        cand = plan["candidates"]
        usable_sites = (tuple((r["enzyme"], r["site"], tuple(r["cut_positions"]))
                              for _, r in cand[cand["usable"]].iterrows())
                        if len(cand) else ())
        try:
            img_map = png_plasmid_map(
                payload, display_name, view == "Whole plasmid", usable_sites,
                (max(1, (auto_region or (1, stats["length"]))[0] - 160),
                 min(stats["length"], (auto_region or (1, stats["length"]))[1] + 160)))
            st.image(img_map, **FIT)
            st.caption("Enzyme sites in red; features from the file in blue, motifs "
                       "detected here in orange, a VHH scaffold in green.")
        except Exception as exc:
            st.warning(f"Map unavailable: {exc}")
        plt.close("all")

    with st.expander(f"Vector features ({len(features)})"):
        st.dataframe(
            features[[c for c in ["feature", "kind", "source", "start", "end",
                                  "strand", "type", "note"]
                      if c in features.columns]], hide_index=True, **FIT)
        st.caption("The *source* column says whether a feature was read from the file "
                   "or inferred here by motif scan.")

    # No footer() here: the router draws it after render() returns, so
    # that every screen gets exactly one and a new screen cannot forget
    # it. The footer() calls above are different - they precede an early
    # st.stop(), which means render() never returns and the router's
    # call is never reached.
