"""
elsewhere - the questions this app does not answer, and who answers them.

Every other module here exists to compute something. This one exists to say
what is not computed, and where to go instead. That is not modesty: it is the
same commitment the rest of the app makes in its captions, collected into one
place a user can read before trusting a number rather than after.

Three kinds of entry sit in the table below.

**Not implemented, deliberately.** A published method whose output this app
could approximate and does not, because approximating it would produce a
number on a scale nobody has calibrated. The TAP patch scores are the clear
case: reimplementing them from a description, without the reference code,
would give values that cannot be compared to the published thresholds, and a
flag against an uncalibrated value is worse than no flag.

**Approximated here, validated there.** This app computes something in the
same spirit with weights it states and does not fit to data. The aggregation
profile against CamSol is the example, and the manual already says so in
§3.8. The entry names the validated method so the comparison is one a user
can actually make.

**Out of reach.** A method this app cannot host, almost always because it
needs model weights or a GPU, against a container that may have 690 MB. The
deep-learning nativeness and thermostability predictors are all here.

**Every citation below was checked against the publication or its record.**
The DOI is the identifier, not the URL: a web server moves, a DOI does not.
Where a tool's own server is given as well, it is given as the authors'
current address and may be stale; the DOI is what stays resolvable.
"""

from __future__ import annotations

import pandas as pd

__version__ = "2026.11.12b"

__all__ = ["TOOLS", "KINDS", "table", "for_question", "__version__"]


KINDS = {
    "not implemented": "A published method this app could approximate and "
                       "does not, because the approximation would not be "
                       "comparable to the published values.",
    "approximated here": "This app computes something in the same spirit, "
                         "with weights it states rather than fits. Use the "
                         "named method for a validated number.",
    "out of reach": "Needs model weights, a GPU or a database this app "
                    "cannot host in a Community Cloud container.",
    "data": "A curated collection rather than a predictor. Useful for "
            "putting your own measurement in context.",
}


# Ordered by the question a user arrives with, not by method name.
TOOLS = (
    {
        "question": "Is this nanobody like the ones in the clinic?",
        "tool": "TNP, Therapeutic Nanobody Profiler",
        "determines": "Six developability metrics benchmarked against the 36 "
                      "clinical-stage nanobodies: total CDR length, CDR3 "
                      "length, CDR3 compactness, and patches of surface "
                      "hydrophobicity, positive charge and negative charge.",
        "here": "The two length metrics are computed on the Analysis screen "
                "(Developability, Clinical context) and compactness on the "
                "Structure screen. The three patch scores are not computed.",
        "kind": "not implemented",
        "citation": "Gordon GL, Gervasio J, Souders C, Deane CM. Commun Biol "
                    "9:344 (2026).",
        "doi": "10.1038/s42003-026-09594-y",
        "server": "https://opig.stats.ox.ac.uk/webapps/tnp",
    },
    {
        "question": "How human, and how VHH-like, is this sequence?",
        "tool": "AbNatiV",
        "determines": "Nativeness: how likely a sequence is to belong to the "
                      "distribution of immune-system-derived human antibodies "
                      "or of camelid nanobodies, as an overall score and a "
                      "per-residue profile, with an automated humanisation "
                      "pipeline built on it.",
        "here": "Percent identity to the closest human germline V gene over "
                "framework positions, with the divergent positions listed "
                "(Analysis, Hallmarks and germline). That is a first-pass "
                "indicator and not a learned distribution.",
        "kind": "out of reach",
        "citation": "Ramon A, Ali M, Atkinson M, et al. Nat Mach Intell "
                    "6:74-91 (2024).",
        "doi": "10.1038/s42256-023-00778-3",
        "server": "",
    },
    {
        "question": "Will it stick to things it should not?",
        "tool": "Nanobody polyreactivity predictor",
        "determines": "Polyreactivity from the CDR sequences of a nanobody, "
                      "with substitutions ranked by how much they are "
                      "predicted to reduce it.",
        "here": "Nothing. Polyreactivity is not predicted anywhere in this "
                "app, and no column here should be read as a proxy for it.",
        "kind": "not implemented",
        "citation": "Harvey EP, Shin JE, Skiba MA, et al. Nat Commun 13:7554 "
                    "(2022).",
        "doi": "10.1038/s41467-022-35276-4",
        "server": "",
    },
    {
        "question": "How soluble is it, as a number I can trust?",
        "tool": "CamSol",
        "determines": "Intrinsic solubility profile and score from sequence, "
                      "with a structurally corrected profile when a structure "
                      "is available, parameterised against measured "
                      "solubilities.",
        "here": "The aggregation profile (Analysis, Aggregation) combines "
                "charge, hydropathy and beta propensity with weights this app "
                "chose and states. It ranks candidates; it does not predict a "
                "solubility.",
        "kind": "approximated here",
        "citation": "Sormanni P, Aprile FA, Vendruscolo M. J Mol Biol "
                    "427:478-490 (2015).",
        "doi": "10.1016/j.jmb.2014.09.026",
        "server": "",
    },
    {
        "question": "What melting temperature should I expect?",
        "tool": "NanoMelt",
        "determines": "Apparent melting temperature of a nanobody, from an "
                      "ensemble of regressions over sequence embeddings "
                      "trained on 640 measurements.",
        "here": "Nothing predicts Tm here. The aliphatic index and the "
                "instability index are correlates of stability in general and "
                "are not a temperature.",
        "kind": "out of reach",
        "citation": "Ramon A, Ni M, Predeina O, et al. mAbs 17:2442750 "
                    "(2025).",
        "doi": "10.1080/19420862.2024.2442750",
        "server": "",
    },
    {
        "question": "What melting temperature should I expect?",
        "tool": "TEMPRO",
        "determines": "Nanobody melting temperature from protein language "
                      "model embeddings, trained on 567 sequences with "
                      "measured Tm.",
        "here": "Nothing. Listed beside NanoMelt because two predictors "
                "trained on overlapping data will not always agree, and "
                "seeing both is more informative than picking one.",
        "kind": "out of reach",
        "citation": "Alvarez JAE, Dean SN. Sci Rep (2024).",
        "doi": "10.1038/s41598-024-70101-6",
        "server": "",
    },
    {
        "question": "What melting temperatures have other people measured?",
        "tool": "NbThermo",
        "determines": "Curated melting temperatures for hundreds of "
                      "nanobodies from the literature, subclassified by the "
                      "method used to measure them (nanoDSF, DSF with SYPRO "
                      "Orange, DSC, CD).",
        "here": "Nothing comparable. If you have measured a Tm, this is "
                "where to see whether it is ordinary.",
        "kind": "data",
        "citation": "Valdes-Tresanco MS, Valdes-Tresanco ME, Molina-Abad E, "
                    "Moreno E. Database 2023:baad021 (2023).",
        "doi": "10.1093/database/baad021",
        "server": "",
    },
    {
        "question": "Has a structure of something like this been solved?",
        "tool": "SAbDab2, single-domain entries",
        "determines": "Every nanobody structure in the PDB, annotated and "
                      "numbered, updated weekly as new ones are deposited. "
                      "Data under CC-BY 4.0.",
        "here": "A snapshot of its single-domain entries ships with the app "
                "and drives the closest-structure search (Analysis, "
                "Repertoire, Closest structure). The snapshot is static and "
                "carries its extraction date; the live database is updated "
                "weekly and has a sequence search of its own.",
        "kind": "data",
        "citation": "Schneider C, Raybould MIJ, Deane CM. Nucleic Acids Res "
                    "50(D1):D1368-D1372 (2022).",
        "doi": "10.1093/nar/gkab1050",
        "server": "",
    },
    {
        "question": "What nanobody sequences are already published or patented?",
        "tool": "PLAbDab-nano",
        "determines": "Camelid and shark single-domain antibody sequences "
                      "mined from patents and the literature, paired with "
                      "their targets where those are stated.",
        "here": "The repertoire comparison uses a fixed set of natural VHHs "
                "for position frequencies; it is not a searchable corpus of "
                "named binders.",
        "kind": "data",
        "citation": "Gordon GL, Raybould MIJ, Wong A, Deane CM. Nucleic Acids "
                    "Res (2024).",
        "doi": "10.1093/nar/gkae881",
        "server": "",
    },
)


def table(kind: str | None = None) -> pd.DataFrame:
    """The tools, as a frame, optionally filtered to one kind.

    Built from `TOOLS` rather than written out again, which is the same rule
    the rest of the app follows: a second copy of a list is a second thing to
    keep right.
    """
    rows = [t for t in TOOLS if kind is None or t["kind"] == kind]
    return pd.DataFrame([{
        "Question": t["question"],
        "Tool": t["tool"],
        "What it determines": t["determines"],
        "What this app does instead": t["here"],
        "Status here": t["kind"],
        "Citation": t["citation"],
        "DOI": t["doi"],
    } for t in rows])


def for_question(fragment: str) -> list[dict]:
    """Entries whose question mentions ``fragment``, case-insensitively.

    For a screen that wants to point at the right line from a caption rather
    than repeat it.
    """
    needle = fragment.lower()
    return [dict(t) for t in TOOLS if needle in t["question"].lower()]
