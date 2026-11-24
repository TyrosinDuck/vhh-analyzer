"""
guides - From a question to the screens that answer it.

The app is organised by tool: one section reads chromatograms, another
annotates sequences, a third compares them. Somebody who arrives with a
question, "I have the .ab1 files of a screening, what did I get?", has to know
which tools to chain and in which order. A guide writes that chain down, step
by step, without moving anything: every step points at a screen that already
exists.

Nothing here computes and nothing here draws. A guide is data, and the
Guides screen (`sections/guides.py`) renders it.

Every step names its screen by key and, inside the analysis, its family by key
and its tab by name, never by a label typed out. `problems` checks each of
them against `shell.SECTIONS`, `shell.UTILITIES` and `shell.ANALYSIS_FAMILIES`,
and a test fails when one no longer exists. A guide that says "Comparison >
Clonal families" after the tab has been renamed would mislead exactly the
reader it was written for, so a rename has to break the build, not the guide.
"""

from __future__ import annotations

from dataclasses import dataclass

__version__ = "2026.11.24c"


@dataclass(frozen=True)
class Step:
    """One thing to do, on one screen.

    `screen` is a key of `shell.SECTIONS` or `shell.UTILITIES`. `family` and
    `tab` are only meaningful on the analysis screen: `family` is the internal
    key of an entry of `shell.ANALYSIS_FAMILIES` ("compare", not
    "🔗 Comparison"), `tab` one of that family's tab names. `view` is for a
    screen that has a selector of its own, and is the internal key of an
    entry of that screen's dict in `shell.SCREEN_VIEWS` ("batch", not
    "📄 Batch sheet"). A step has either a family, with or without a tab, or
    a view, never both: they belong to different screens.
    """
    screen: str
    action: str
    why: str = ""
    family: str | None = None
    tab: str | None = None
    view: str | None = None
    optional: bool = False


@dataclass(frozen=True)
class Guide:
    """A question, and the steps that answer it, in order."""
    title: str
    question: str
    steps: tuple[Step, ...]
    note: str = ""


GUIDES: dict[str, Guide] = {
    "sanger_screening": Guide(
        title="Screening clones by Sanger sequencing",
        question=("You picked colonies after a selection, sequenced them, and "
                  "have the .ab1 files back. Which distinct VHHs did you get, "
                  "and what are they like?"),
        steps=(
            Step("sequencing",
                 "Upload the .ab1 files. Without files to hand, *Load an "
                 "example run* gives six simulated reads to try the path with.",
                 "Each read is quality-trimmed and the VHH is extracted from "
                 "the expression cassette, in whichever strand and frame it "
                 "sits."),
            Step("sequencing",
                 "Choose the vector the clones are in.",
                 "It sets which tags the app expects after the VHH. The wrong "
                 "vector is the usual reason for a read reported as having no "
                 "leader."),
            Step("sequencing",
                 "Give the reads your own names in *Name the reads*.",
                 "The provider's well references are kept otherwise. A name "
                 "given here follows the clone into every table and figure "
                 "of the analysis.",
                 optional=True),
            Step("sequencing",
                 "Read the results table: the status of each read, and which "
                 "clones are duplicates of another.",
                 "A failed read is not necessarily a failed clone. A frameshift "
                 "read has its VHH withheld on purpose, and a read that stops "
                 "early may only need sequencing from the other end."),
            Step("sequencing",
                 "Under *Send to the analysis section*, keep the default "
                 "selection and click *Send selected to analysis*.",
                 "Duplicates are left out by default, one entry per unique "
                 "protein. The analysis starts on its own once the sequences "
                 "arrive."),
            Step("analysis",
                 "See how many distinct clones there are.",
                 "Sequences are grouped into families on their CDR3, so two "
                 "clones differing elsewhere still count as one lineage. The "
                 "identical_reads column says how many reads gave each "
                 "family, which the one-entry-per-protein selection would "
                 "otherwise hide. It counts reads, not clones: a clone read "
                 "from both ends is two.",
                 family="compare", tab="Clonal families"),
            Step("analysis",
                 "Check where the CDRs and frameworks start and end.",
                 family="identity", tab="CDR / Framework"),
            Step("analysis",
                 "Check the hallmark residues and the closest germline.",
                 "The hallmarks tell a VHH from a conventional VH.",
                 family="identity", tab="Hallmarks & germline"),
            Step("analysis",
                 "Look at the physico-chemical properties of each clone.",
                 "Mass, pI and extinction coefficient are computed from the "
                 "sequence. They are theoretical values, not measurements.",
                 family="identity", tab="Properties"),
            Step("analysis",
                 "Put every flag side by side to choose which clones to carry "
                 "forward.",
                 "There is no overall score: sort by the column that matters "
                 "for what the VHH is for.",
                 family="triage", tab="Table"),
            Step("analysis",
                 "Save the session file, and download the tables or the PDF "
                 "report if you need them.",
                 "The app stores nothing. The session file is what lets "
                 "you, or somebody else, reopen this analysis later.",
                 family="outputs", tab="Session"),
        ),
        note=("The settings of the analysis, pH, profile window and CDR "
              "annotation method, are left at their defaults by this path. "
              "They are in the sidebar of the analysis screen."),
    ),
    "purified_batch": Guide(
        title="From a purified tube to its batch sheet",
        question=("You have purified a VHH and the tube is in front of you. "
                  "How much is in it, is it the right protein, what did the "
                  "purification recover, and what leaves with the tube?"),
        steps=(
            Step("analysis",
                 "Load the sequence of the protein in the tube and press "
                 "*Analyze*.",
                 "Every view below starts from the molar mass, the extinction "
                 "coefficient and the lysine count derived from the sequence, "
                 "so none of them has to be typed. The bench screen can also "
                 "take a pasted sequence, and typing those numbers by hand is "
                 "the most common way of getting them wrong."),
            Step("bench",
                 "Pick that protein at the top of the screen, then enter your "
                 "measured A280, the dilution factor and the path length.",
                 "The concentration is **derived**, not measured: it is "
                 "Beer-Lambert on your reading with an extinction coefficient "
                 "predicted from the sequence. The screen labels it as such, "
                 "and so does everything downstream of it. Choose cystines "
                 "formed or all reduced; a VHH from the periplasm has its "
                 "canonical disulfide formed.",
                 view="quantify"),
            Step("bench",
                 "If you measured a mass by MS, check it against the sequence.",
                 "It answers whether the thing in the tube is the thing you "
                 "think it is, and when it is not, it lists the "
                 "modifications your sequence allows that would explain the "
                 "difference.",
                 view="mass", optional=True),
            Step("bench",
                 "Enter one row per purification step, with its volume and "
                 "its concentration, and tick *Crude* on the lysate, the "
                 "supernatant and the flow-through.",
                 "An A280 on a crude step is the absorbance of every protein "
                 "in it, so a recovery computed from it is the most "
                 "flattering number on the page and the least true. Crude "
                 "rows keep their total and stay out of the recovery "
                 "arithmetic.",
                 view="balance"),
            Step("bench",
                 "If you are labelling the protein, work out how much reagent "
                 "to add.",
                 "The excess is arithmetic on the molar mass and the lysine "
                 "count. Nothing here predicts a labelling efficiency.",
                 view="conjugate", optional=True),
            Step("bench",
                 "Fill in the lot number, the host, the purification, the "
                 "buffer and the storage, retype what you measured, and "
                 "download the PDF.",
                 "The sheet keeps four blocks apart: the free text, what an "
                 "instrument read, what is derived from both, and what comes "
                 "from the sequence. Nothing is carried over from the views "
                 "above, because a number you have not retyped is a number "
                 "you have not checked. Save the JSON beside the PDF if the "
                 "batch is going anywhere other than a freezer box.",
                 view="batch"),
            Step("bench",
                 "Once you have several lots, read them side by side.",
                 "It is where a preparation that drifted shows up.",
                 view="lots", optional=True),
        ),
        note=("This is not a certificate of analysis, and the PDF says so at "
              "the bottom: there are no specifications here, no acceptance "
              "criteria and no release decision."),
    ),
    "structure_measured": Guide(
        title="Replacing sequence guesses with measurements from a model",
        question=("You have a crystal structure or a predicted model of your "
                  "VHH. Which cysteines really pair, and which of the "
                  "liabilities the analysis flagged are actually exposed?"),
        steps=(
            Step("analysis",
                 "Load the full-length sequence of the construct and press "
                 "*Analyze*.",
                 "A model is routinely missing termini and disordered loops, "
                 "and its numbering may be sequential, IMGT, Kabat or the "
                 "depositor's own. Having the sequence here lets the next "
                 "step align the model onto it.",
                 optional=True),
            Step("structure",
                 "Upload the .pdb, .cif or AlphaFold3 .zip. Without a file to "
                 "hand, *Load an example structure* brings in PDB 1MEL, a VHH "
                 "bound to lysozyme.",
                 "The app does not download structures for you, on purpose: "
                 "it would mean this server making requests that look like "
                 "they came from you. For a complex, download the whole entry "
                 "rather than one chain."),
            Step("structure",
                 "Set the chain of interest, the binding partners, and point "
                 "*Annotate against* at the sequence you loaded.",
                 "Left on the structure's own sequence, everything is read "
                 "off what is in the file. Pointed at your construct, the "
                 "CDR bands and liability positions line up with the rest of "
                 "the app. The alignment identity is reported, and a warning "
                 "below 90% usually means the wrong molecule."),
            Step("structure",
                 "Check what the B-factor column actually holds.",
                 "A predicted model stores pLDDT where a crystal structure "
                 "stores a B-factor. High pLDDT is good and high B is bad, "
                 "and no field in either format says which it is, so the app "
                 "guesses from the values and lets you override it.",
                 view="chains"),
            Step("structure",
                 "Read which cysteines are actually bonded.",
                 "Measured from the SG-SG distances, not inferred from IMGT "
                 "positions. The count is compared with the sequence-based "
                 "prediction and with the file's own records, and "
                 "disagreements are reported rather than resolved. Check the "
                 "model covers every cysteine: one in an unmodelled loop "
                 "cannot be seen to pair.",
                 view="disulfides"),
            Step("structure",
                 "Find out which liabilities are buried and which are "
                 "exposed.",
                 "The join worth uploading a structure for. A buried "
                 "deamidation motif is largely inert; an exposed one in a CDR "
                 "is a real problem. The sequence-based risk column is left "
                 "untouched, so read the two together. A motif in an "
                 "unmodelled region is reported as not modelled, which is not "
                 "the same as buried.",
                 view="liabilities"),
            Step("structure",
                 "Look at the hydrophobic surface as area that is contiguous "
                 "in space.",
                 "The aggregation profile in the analysis is a sliding window "
                 "over the sequence, which cannot know that two distant "
                 "stretches fold into one patch.",
                 view="patches", optional=True),
            Step("structure",
                 "With a complex, read which residues contact the antigen.",
                 "The one thing on this screen that cannot be inferred from "
                 "sequence at all. The view appears only when a binding "
                 "partner is selected.",
                 optional=True),
        ),
        note=("Everything here is measured on one static conformer, which is "
              "a strong claim about a molecule that moves. Two models of the "
              "same VHH can disagree, and a prediction is not a structure."),
    ),
}


def view_label(screen: str, key: str, screen_views: dict) -> str:
    """The label a view is shown under, from its screen and its internal key.

    `screen_views` is `shell.SCREEN_VIEWS`. Raises KeyError for a screen with
    no views or an unknown key, which is what `problems` reports.
    """
    _state_key, views = screen_views[screen]
    return views[key]


def family_label(key: str, families: dict) -> str:
    """The label a family is shown under, from its internal key.

    `families` is `shell.ANALYSIS_FAMILIES`, passed in rather than imported so
    that this module stays importable without Streamlit. Raises KeyError for
    an unknown key, which is what `problems` reports.
    """
    for label, (family_key, _tabs, _blurb) in families.items():
        if family_key == key:
            return label
    raise KeyError(key)


def location(step: Step, screens: dict[str, str], families: dict,
             screen_views: dict | None = None) -> str:
    """Where a step happens, as a reader finds it: "Screen › Family › Tab".

    `screens` maps a screen key to the label of its sidebar button,
    `screen_views` is `shell.SCREEN_VIEWS`. Every part is looked up, none is
    stored with the step.
    """
    parts = [screens[step.screen]]
    if step.family is not None:
        parts.append(family_label(step.family, families))
    if step.tab is not None:
        parts.append(step.tab)
    if step.view is not None:
        parts.append(view_label(step.screen, step.view, screen_views or {}))
    return " › ".join(parts)


def problems(guides: dict[str, Guide], screens: dict[str, str],
             families: dict, analysis_screen: str = "analysis",
             screen_views: dict | None = None) -> list[str]:
    """Every reference in `guides` that does not resolve, as readable lines.

    An empty list means every step points at something that exists. The
    checks: the screen is a known screen; a family or a tab only appear on the
    analysis screen; a tab always comes with its family; the family exists;
    the tab is one of that family's tabs; a view is on a screen that has one
    and is one of that screen's views, and never sits beside a family or a
    tab. A guide with no step is reported too, since it would draw a title
    and nothing to do.
    """
    screen_views = screen_views or {}
    found = []
    for key, guide in guides.items():
        if not guide.steps:
            found.append(f"{key}: no steps")
        for number, step in enumerate(guide.steps, 1):
            where = f"{key}, step {number}"
            if step.screen not in screens:
                found.append(f"{where}: unknown screen {step.screen!r}")
            if step.view is not None:
                if step.family or step.tab:
                    found.append(f"{where}: a view beside a family or a tab")
                elif step.screen not in screen_views:
                    found.append(f"{where}: {screens.get(step.screen, step.screen)} "
                                 f"has no view selector")
                elif step.view not in screen_views[step.screen][1]:
                    found.append(
                        f"{where}: unknown view {step.view!r}, "
                        f"{step.screen} has "
                        f"{sorted(screen_views[step.screen][1])}")
            if step.screen != analysis_screen and (step.family or step.tab):
                found.append(f"{where}: a family or a tab outside the "
                             f"analysis screen")
                continue
            if step.tab is not None and step.family is None:
                found.append(f"{where}: tab {step.tab!r} without a family")
                continue
            if step.family is None:
                continue
            try:
                label = family_label(step.family, families)
            except KeyError:
                found.append(f"{where}: unknown family {step.family!r}")
                continue
            tabs = families[label][1]
            if step.tab is not None and step.tab not in tabs:
                found.append(f"{where}: {step.tab!r} is not a tab of "
                             f"{label}, which has {tabs}")
    return found
