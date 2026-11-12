"""
session - the inputs of one analysis, written to a file, so the run repeats.

The app stores nothing. That is a deliberate constraint and it has a cost: an
analysis you ran in March cannot be looked up in September, and the figure you
put in a manuscript rests on a page nobody can rebuild. "I used the VHH tool
box" names a tool, not a run: the pH, the profile window and the annotator are
all choices, and the annotator in particular can differ between two servers
running the same code, because ANARCI is installed beside the app rather than
by it.

This module turns one analysis into two things, both of which the user keeps:

* a **session file**, the inputs in JSON, which the app reads back to put the
  same sequences and the same settings on the screen;
* a **methods paragraph**, prose describing what was computed, assembled from
  the same record, so it cannot describe a setting that was not used.

What is *not* in the record: any computed value. A session file is not a
result, it is the question that produced one. Reloading it recomputes
everything with the version that is running at that moment, and
``created_with`` is in the file so that a run which disagrees with the original
can be told from one that does not. Keeping the numbers would hide exactly that
case, which is the one worth catching.

What a session file cannot do: pin the environment. It records which Python,
Streamlit, Biopython and ANARCI the run saw, because those are what the numbers
came out of, but recording is not pinning. Reading it on a machine with a
different ANARCI tells you the two runs are not comparable; it does not make
them comparable.

Nothing here imports Streamlit. A notebook can write and read these files, and
the tests do.
"""

from __future__ import annotations

__version__ = "2026.11.12a"

__all__ = [
    "SESSION_SCHEMA", "METHODS_SETTINGS", "session_record", "session_inputs",
    "methods_paragraph", "__version__",
]

# The identifier a reader checks before trusting the keys below, with the same
# discipline as `bench.BATCH_RECORD_SCHEMA`: a record whose format number this
# code does not know is refused rather than half read, because the settings
# that did not survive would be missing from the restored screen with nothing
# saying so, and a restored screen that quietly differs from the saved one is
# worse than no restore at all.
SESSION_SCHEMA = "vhh-toolbox-session/1"

# The analysis settings a session carries, each with the bounds the sidebar
# widget enforces. This is the one place they are listed: `session_inputs`
# validates against it and `methods_paragraph` describes what it finds in it,
# so a setting added to a session cannot be restored without being described
# or described without being restored.
#
# The bounds are the widgets', not physics: the pH slider runs 3 to 10 and the
# window slider 3 to 21, and a file asking for a value outside them could not
# be put back on the screen.
METHODS_SETTINGS = {
    "ph": (3.0, 10.0),
    "window": (3, 21),
}

ANNOTATION_METHODS = ("auto", "anarci", "regex")


def session_record(sequences: dict[str, str], ph: float, window: int,
                   method: str, method_used: str = "",
                   version: str = "", release: str = "", build: str = "",
                   code: str = "", runtime: dict | None = None) -> dict:
    """One analysis as data: the sequences, the settings, and what ran them.

    ``method`` is what the sidebar was set to, which may be ``"auto"``.
    ``method_used`` is what actually annotated the sequences, as the annotator
    reports itself ("ANARCI", "Regex anchors"). Both are kept, because they
    answer different questions: the first is reproducible by setting a widget,
    the second is not reproducible at all on a machine where ANARCI is
    missing, and a record holding only ``"auto"`` would look reproducible while
    hiding that.

    ``code`` is the source fingerprint and ``build`` the build date, the two
    halves of the footer. Together with ``release`` they identify the code that
    ran more precisely than a release number does: two builds of 1.1.0 are
    possible, and have happened.

    ``runtime`` is what the server process was running, from
    ``shell.runtime_versions()``. It is recorded rather than checked: see the
    module docstring.
    """
    from datetime import datetime, timezone

    clean: dict[str, str] = {}
    for name, seq in (sequences or {}).items():
        residues = "".join(c for c in str(seq or "").upper() if c.isalpha())
        if residues:
            clean[str(name)] = residues

    return {
        "schema": SESSION_SCHEMA,
        "created_with": version,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "release": release,
        "build": build,
        "code": code,
        "runtime": dict(runtime or {}),
        "settings": {
            "ph": float(ph),
            "window": int(window),
            "annotation_method": str(method),
            "annotation_used": str(method_used or ""),
        },
        "sequences": clean,
    }


def session_inputs(data: dict) -> dict:
    """Read a session file back into the inputs that made it.

    Returns ``sequences``, ``ph``, ``window``, ``method``, ``method_used``,
    ``created_with``, ``release``, ``build``, ``code`` and ``runtime``.

    Refuses rather than guesses, in three cases, each of which would otherwise
    restore a screen that is not the one that was saved:

    * a schema this code does not read, as above;
    * a settings value outside the bounds of the widget that collects it,
      which could not be put back on the screen;
    * a sequence that is not a string of letters.

    An empty ``sequences`` is allowed and returned empty: a session saved
    before anything was analysed is a legitimate thing to have saved, and the
    caller decides whether it is worth restoring.
    """
    if not isinstance(data, dict):
        raise ValueError("a session file is a JSON object, and this is not one")
    schema = str(data.get("schema") or "")
    if not schema:
        raise ValueError(
            "this file carries no `schema` key, so it is not a session file "
            "written by this app")
    family, _, number = schema.rpartition("/")
    expected_family, _, expected_number = SESSION_SCHEMA.rpartition("/")
    if family != expected_family:
        raise ValueError(
            f"this file says it is `{schema}`, and this app reads "
            f"`{SESSION_SCHEMA}`")
    if number != expected_number:
        raise ValueError(
            f"this session is format {schema} and this app reads "
            f"{SESSION_SCHEMA}. A newer app wrote it, or an older one: either "
            f"way the keys are not guaranteed to mean the same thing, so it "
            f"is not being read rather than read wrongly")

    settings = data.get("settings")
    if settings is not None and not isinstance(settings, dict):
        raise ValueError("`settings` in this session file is not an object")
    settings = settings or {}

    numbers: dict[str, float] = {}
    for key, (low, high) in METHODS_SETTINGS.items():
        raw = settings.get(key)
        if raw in (None, ""):
            raise ValueError(
                f"`settings.{key}` is missing from this session file, so the "
                f"run cannot be set up as it was")
        try:
            value = float(raw)
        except (TypeError, ValueError):
            raise ValueError(
                f"`settings.{key}` is {raw!r}, which is not a number") from None
        if not low <= value <= high:
            raise ValueError(
                f"`settings.{key}` is {value:g}, outside the {low:g} to "
                f"{high:g} this app offers, so it could not be put back on "
                f"the screen")
        numbers[key] = value

    method = str(settings.get("annotation_method") or "")
    if method not in ANNOTATION_METHODS:
        raise ValueError(
            f"`settings.annotation_method` is {method!r}, and this app knows "
            f"{', '.join(ANNOTATION_METHODS)}")

    raw_sequences = data.get("sequences")
    if raw_sequences is not None and not isinstance(raw_sequences, dict):
        raise ValueError("`sequences` in this session file is not an object")
    sequences: dict[str, str] = {}
    for name, seq in (raw_sequences or {}).items():
        if not isinstance(seq, str):
            raise ValueError(
                f"the sequence stored under {name!r} is not text")
        residues = "".join(c for c in seq.upper() if c.isalpha())
        if not residues:
            raise ValueError(f"the sequence stored under {name!r} is empty")
        sequences[str(name)] = residues

    runtime = data.get("runtime")
    if runtime is not None and not isinstance(runtime, dict):
        raise ValueError("`runtime` in this session file is not an object")

    return {
        "sequences": sequences,
        "ph": numbers["ph"],
        # int() because it indexes a sliding window. A file saying 9.0 is
        # fine, one saying 9.5 is not, and the bounds check above has already
        # let it through, so the rounding has to be refused rather than done.
        "window": _whole(numbers["window"], "settings.window"),
        "method": method,
        "method_used": str(settings.get("annotation_used") or ""),
        "created_with": str(data.get("created_with") or ""),
        "created_utc": str(data.get("created_utc") or ""),
        "release": str(data.get("release") or ""),
        "build": str(data.get("build") or ""),
        "code": str(data.get("code") or ""),
        "runtime": {str(k): str(v) for k, v in (runtime or {}).items()},
    }


def _whole(value: float, where: str) -> int:
    if float(value) != int(value):
        raise ValueError(
            f"`{where}` is {value:g}, and it counts residues, so it has to be "
            f"a whole number")
    return int(value)


def fasta(sequences: dict[str, str], width: int = 60) -> str:
    """Sequences as FASTA, wrapped, for a session file to be pasted back in.

    Here rather than in the screen because the session file and the text area
    should round-trip through the same formatting, and because a notebook
    reading a session file wants the FASTA without importing a screen.
    """
    out = []
    for name, seq in (sequences or {}).items():
        out.append(f">{name}")
        out.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(out) + ("\n" if out else "")


def methods_paragraph(record: dict, app_name: str = "VHH tool box",
                      doi: str = "", repo_url: str = "") -> str:
    """A Material and Methods paragraph for exactly the run in ``record``.

    Written from the record and from nothing else, which is the whole point:
    every sentence names a setting that was actually used, and a setting the
    app did not use produces no sentence. The alternative, a template with
    blanks, is how a methods section comes to describe ANARCI numbering on a
    server that fell back to the regex anchors.

    Three things it says that a user would not think to write, and that a
    reviewer is entitled to:

    * which annotator ran, not which one was requested;
    * that the mass, pI and extinction coefficients are ProtParam values
      computed from the sequence, so they are theoretical and were not
      measured;
    * the version and the source fingerprint of the code, because the release
      number alone does not identify a build.

    Sequence names are deliberately left out. The paragraph is meant to be
    pasted into a manuscript, and the names in an analysis are often internal
    clone identifiers.

    It describes the sequence characterisation and the annotation, which is
    what a session covers. It says nothing about the bench calculators, the
    cloning design or a structure, because a session file records none of
    those: a paragraph claiming them would be describing work this record
    cannot vouch for.
    """
    inputs = record if "sequences" in record and "ph" in record \
        else session_inputs(record)

    count = len(inputs["sequences"])
    subject = ("One amino-acid sequence was" if count == 1
               else f"{count} amino-acid sequences were")
    if count == 0:
        subject = "Amino-acid sequences were"

    version = inputs.get("created_with") or inputs.get("release") or ""
    stamp = f" ({version})" if version else ""
    code = inputs.get("code")
    fingerprint = (f", source fingerprint {code}"
                   if code and code != "unknown" else "")
    identity = f"{app_name}{stamp}{fingerprint}"

    used = inputs.get("method_used") or ""
    asked = inputs.get("method")
    if used:
        annotator = f"{used.rstrip('.')}"
        if asked == "auto":
            annotator += (", which is what the automatic setting resolved to "
                          "on the machine that ran it")
        annotation = (f"Complementarity-determining regions and frameworks "
                      f"were delimited with {annotator}.")
    else:
        annotation = ("The record does not say which annotator ran, so the "
                      "CDR and framework boundaries in it cannot be "
                      "attributed to one.")

    sentences = [
        f"{subject} characterised with {identity}"
        + (f", archived at doi:{doi}" if doi else "")
        + (f" and available at {repo_url}" if repo_url else "") + ".",
        f"Net charge was computed at pH {inputs['ph']:g}, and the "
        f"per-residue charge and hydropathy profiles with a sliding window of "
        f"{inputs['window']} residues.",
        annotation,
        "Molar mass, isoelectric point and the two extinction coefficients at "
        "280 nm are ProtParam values computed from the sequence, with "
        "Biopython's implementation: they are theoretical, and none of them "
        "was measured here.",
    ]

    runtime = inputs.get("runtime") or {}
    if runtime:
        sentences.append(
            "The run used " + ", ".join(f"{name} {value}" for name, value
                                        in runtime.items()) + ".")

    return " ".join(sentences)
