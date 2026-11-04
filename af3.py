"""
af3, the confidence AlphaFold3 writes next to its models, read from the zip it
hands you.

A predicted structure carries one confidence score inside the file itself:
pLDDT, in the B-factor column. pLDDT is *local*. It says the fold around each
residue is probably right, and nothing about where one chain sits relative to
another. A VHH and its antigen can both be dark blue and still be docked on the
wrong face. Whether the arrangement of the chains can be trusted is what ipTM
and the predicted aligned error (PAE) are for, and AlphaFold3 writes those in
JSON files beside the `.cif`, not in it.

This module reads them from the archive AlphaFold3 produces, as downloaded:

- **AlphaFold Server**: ``fold_<job>_model_<n>.cif`` with
  ``fold_<job>_summary_confidences_<n>.json`` and ``fold_<job>_full_data_<n>.json``
  for n = 0..4. Model 0 is the top-ranked one.
- **AlphaFold3 run locally**, zipped: ``<job>_model.cif`` with
  ``<job>_summary_confidences.json`` and ``<job>_confidences.json`` at the top,
  and the same three names without the prefix in each ``seed-*_sample-*/``
  folder.

The pairing rule covers both: a model's JSON files are the ones in the same
folder whose name is the model's with the last ``model`` replaced. Only those
files are read. ``msas/`` and ``templates/`` are never opened, and nothing is
written to disk.

**Chain order.** ``chain_pair_iptm`` is a matrix with no labels. Its order is
the order of the chains in the job, which is the order of the chains in the
`.cif`. The AlphaFold Server summary also lists ``chain_ids``, and the full PAE
file lists ``token_chain_ids``; when either is present it is used, and the
`.cif` order is the fallback. The source is reported, so a mismatch is visible
rather than silent.

Thresholds are the ones given with AlphaFold3 (Abramson et al. 2024, and the
AlphaFold Server documentation): ipTM above 0.8 is a confident prediction of the
arrangement, below 0.6 the arrangement has likely failed, and the band between
is a grey zone that the rest of the evidence has to decide.

  Abramson J, Adler J, Dunger J, et al. (2024) Accurate structure prediction of
  biomolecular interactions with AlphaFold 3. Nature 630: 493-500.
  doi:10.1038/s41586-024-07487-w
"""

from __future__ import annotations

import io
import json
import posixpath
import re
import zipfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from structure import StructureError

__all__ = [
    "AF3Error", "AF3Sample", "AF3Bundle", "is_zip", "read_af3_zip",
    "read_member", "load_full_data", "chain_order", "pair_iptm", "pair_pae_min",
    "interface_pae", "iptm_verdict", "samples_table", "plot_pae",
    "IPTM_HIGH", "IPTM_LOW", "MAX_MEMBER_BYTES", "__version__",
]

__version__ = "2026.11.04a"

IPTM_HIGH = 0.8
IPTM_LOW = 0.6

# Per file read out of the archive. A full PAE file is about 1 MB for a VHH and
# its antigen and grows with the square of the token count, so 1 000 residues
# is roughly 15 MB. This refuses anything far past what the app is for, and
# refuses it before decompressing: the size is read from the archive's
# directory, so a zip that claims a huge member is turned away unopened.
MAX_MEMBER_BYTES = 60 * 1024 * 1024


class AF3Error(StructureError):
    """An archive that cannot be read as AlphaFold3 output, with a reason."""


@dataclass
class AF3Sample:
    """One model of the job and where its confidence files sit in the zip."""
    label: str
    cif_member: str
    summary_member: str | None
    full_member: str | None
    summary: dict = field(default_factory=dict)

    @property
    def iptm(self):
        return self.summary.get("iptm")

    @property
    def ranking_score(self):
        return self.summary.get("ranking_score")


@dataclass
class AF3Bundle:
    """What the archive holds: the job name and its models, best first."""
    job_name: str
    source: str
    samples: list


def is_zip(data: bytes) -> bool:
    return data[:4] == b"PK\x03\x04"


def _skip(member: str) -> bool:
    parts = member.lower().split("/")
    return (member.endswith("/") or "__macosx" in parts
            or any(p in ("msas", "templates") for p in parts[:-1])
            or parts[-1].startswith("."))


def _companion(cif_member: str, replacement: str) -> str:
    folder, base = posixpath.split(cif_member)
    stem = base[: -len(".cif")]
    cut = stem.rfind("model")
    return posixpath.join(folder,
                          stem[:cut] + replacement + stem[cut + len("model"):]
                          + ".json")


def _read(zf: zipfile.ZipFile, member: str) -> bytes:
    info = zf.getinfo(member)
    if info.file_size > MAX_MEMBER_BYTES:
        raise AF3Error(
            f"{member} is {info.file_size / 1e6:.0f} MB once unpacked, over the "
            f"{MAX_MEMBER_BYTES // (1024 * 1024)} MB this app reads per file. "
            f"Upload the model's .cif on its own instead.")
    return zf.read(member)


def _read_json(zf: zipfile.ZipFile, member: str) -> dict:
    try:
        return json.loads(_read(zf, member))
    except AF3Error:
        raise
    except Exception as exc:
        raise AF3Error(f"{member} is not readable JSON: {exc}") from exc


def _label(cif_member: str) -> str:
    folder, base = posixpath.split(cif_member)
    found = re.search(r"model_(\d+)\.cif$", base)
    if found:
        return f"model {found.group(1)}"
    # A local run writes `model.cif` bare inside each seed-*_sample-*/ folder,
    # and `<job>_model.cif` for the top-ranked one beside them.
    if base.lower() == "model.cif" and folder:
        return posixpath.basename(folder)
    return "top-ranked model"


def read_af3_zip(data: bytes, filename: str = "archive.zip") -> AF3Bundle:
    """Index an AlphaFold3 archive and read its small summary files.

    The models' coordinates and full PAE files are read later, one model at a
    time, by whoever needs them. Samples come back sorted by ranking score,
    best first, which for AlphaFold Server is also model 0 first.
    """
    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise AF3Error(f"{filename} is not a readable zip archive: {exc}") from exc

    with zf:
        names = [n for n in zf.namelist() if not _skip(n)]
        present = set(names)
        cifs = [n for n in names
                if n.lower().endswith(".cif")
                and "model" in posixpath.basename(n).lower()]
        if not cifs:
            raise AF3Error(
                f"{filename} holds no AlphaFold3 model (no `*model*.cif` outside "
                f"msas/ and templates/). Is it the archive AlphaFold3 wrote?")

        samples = []
        for cif in cifs:
            summary_member = _companion(cif, "summary_confidences")
            full_member = next((m for m in (_companion(cif, "full_data"),
                                            _companion(cif, "confidences"))
                                if m in present), None)
            if summary_member not in present:
                summary_member = None
            summary = (_read_json(zf, summary_member) if summary_member else {})
            samples.append(AF3Sample(label=_label(cif), cif_member=cif,
                                     summary_member=summary_member,
                                     full_member=full_member, summary=summary))

        if not any(s.summary for s in samples):
            raise AF3Error(
                f"{filename} holds models but no `summary_confidences` file "
                f"next to any of them, so there is no ipTM to read. Upload the "
                f"model's .cif directly for pLDDT alone.")

        job_name, source = posixpath.splitext(posixpath.basename(filename))[0], \
            "AlphaFold3"
        for name in names:
            if name.lower().endswith("job_request.json"):
                try:
                    request = json.loads(_read(zf, name))
                    request = request[0] if isinstance(request, list) else request
                    job_name = request.get("name") or job_name
                    if request.get("dialect") == "alphafoldserver":
                        source = "AlphaFold Server"
                except Exception:
                    pass
                break
        else:
            if any(posixpath.basename(c).startswith("fold_") for c in cifs):
                source = "AlphaFold Server"

    def rank(sample):
        score = sample.ranking_score
        return (score is None, -(score or 0.0), sample.label)

    samples.sort(key=rank)
    return AF3Bundle(job_name=job_name, source=source, samples=samples)


def read_member(data: bytes, member: str) -> bytes:
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        return _read(zf, member)


def load_full_data(data: bytes, member: str | None) -> dict | None:
    """The PAE matrix and token chain ids for one model, or None if absent.

    ``member`` is a sample's ``full_member``. It is passed rather than the
    sample so that the cached wrapper in `compute.py` is keyed on a string.
    """
    if not member:
        return None
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        full = _read_json(zf, member)
    pae = full.get("pae")
    chains = full.get("token_chain_ids")
    if pae is None or chains is None:
        return None
    pae = np.asarray(pae, dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1] or pae.shape[0] != len(chains):
        raise AF3Error(f"{member}: the PAE matrix is {pae.shape} "
                       f"for {len(chains)} tokens, which do not match.")
    return {"pae": pae, "token_chain_ids": list(chains),
            "token_res_ids": list(full.get("token_res_ids") or [])}


# --------------------------------------------------------------------------
# Reading the numbers
# --------------------------------------------------------------------------

def _unique(seq) -> list:
    seen = []
    for item in seq:
        if item not in seen:
            seen.append(item)
    return seen


def chain_order(summary: dict, full: dict | None = None,
                cif_chains: list | None = None) -> tuple[list, str]:
    """The chain ids, in the order ``chain_pair_iptm`` uses, and where from."""
    n = len(summary.get("chain_pair_iptm") or summary.get("chain_ptm") or [])
    candidates = [
        (_unique(summary.get("chain_ids") or []), "the summary file's chain_ids"),
        (_unique((full or {}).get("token_chain_ids") or []),
         "the PAE file's token_chain_ids"),
        (list(cif_chains or []), "the order of chains in the .cif"),
    ]
    for order, source in candidates:
        if order and (not n or len(order) == n):
            return order, source
    return [], "nowhere: the chain count does not match the summary"


def _pair(matrix, order, chain_a, chain_b):
    if not matrix or chain_a not in order or chain_b not in order:
        return None
    i, j = order.index(chain_a), order.index(chain_b)
    try:
        value = matrix[i][j]
    except (IndexError, TypeError):
        return None
    return None if value is None else float(value)


def pair_iptm(summary: dict, order: list, chain_a: str, chain_b: str):
    """ipTM restricted to one pair of chains. Symmetric in AlphaFold3 output."""
    return _pair(summary.get("chain_pair_iptm"), order, chain_a, chain_b)


def pair_pae_min(summary: dict, order: list, chain_a: str, chain_b: str):
    """The smallest PAE between the two chains, both directions, in Å.

    PAE is not symmetric: the row is the residue the model is aligned on, the
    column the residue whose position is judged. The smaller of the two
    directions is reported, since either says the model has placed one chain
    with confidence relative to the other somewhere.
    """
    values = [_pair(summary.get("chain_pair_pae_min"), order, a, b)
              for a, b in ((chain_a, chain_b), (chain_b, chain_a))]
    values = [v for v in values if v is not None]
    return min(values) if values else None


def interface_pae(full: dict, chain_a: str, chain_b: str,
                  close: float = 5.0) -> dict | None:
    """Summary of the PAE block between two chains.

    ``mean`` is over every token pair, both directions, and is dominated by the
    parts of each chain that are nowhere near the other. ``min`` and the
    fraction of pairs under ``close`` Å are what say whether *some* region is
    placed with confidence, which is what an epitope is.
    """
    chains = np.asarray(full["token_chain_ids"])
    a, b = chains == chain_a, chains == chain_b
    if not a.any() or not b.any():
        return None
    pae = full["pae"]
    block = np.concatenate([pae[np.ix_(a, b)].ravel(), pae[np.ix_(b, a)].ravel()])
    return {"mean": float(block.mean()), "min": float(block.min()),
            "fraction_close": float((block < close).mean()), "close": close}


def iptm_verdict(value) -> tuple[str, str]:
    """('high' | 'grey' | 'low' | 'none', sentence) for an ipTM value."""
    if value is None:
        return "none", "No ipTM: the job has a single chain."
    if value > IPTM_HIGH:
        return "high", (f"above {IPTM_HIGH}: the relative placement of the chains "
                        f"is predicted with confidence.")
    if value < IPTM_LOW:
        return "low", (f"below {IPTM_LOW}: the arrangement has likely failed. "
                       f"Interface residues measured on this model are not "
                       f"evidence of an epitope.")
    return "grey", (f"between {IPTM_LOW} and {IPTM_HIGH}: a grey zone. Correct "
                    f"and incorrect arrangements both land here, so the PAE "
                    f"and agreement between models have to decide.")


def samples_table(bundle: AF3Bundle, chain_a: str | None = None,
                  chain_b: str | None = None,
                  cif_chains: list | None = None) -> pd.DataFrame:
    """One row per model: the global scores and, if a pair is given, its ipTM."""
    rows = []
    for sample in bundle.samples:
        s = sample.summary
        row = {"model": sample.label, "ranking_score": s.get("ranking_score"),
               "iptm": s.get("iptm"), "ptm": s.get("ptm")}
        if chain_a and chain_b:
            order, _ = chain_order(s, None, cif_chains)
            row["pair_iptm"] = pair_iptm(s, order, chain_a, chain_b)
            row["pair_pae_min"] = pair_pae_min(s, order, chain_a, chain_b)
        row["has_clash"] = bool(s.get("has_clash")) if s else None
        row["fraction_disordered"] = s.get("fraction_disordered")
        rows.append(row)
    return pd.DataFrame(rows)


def plot_pae(full: dict, title: str = "Predicted aligned error"):
    """The PAE matrix, with the chain boundaries drawn and labelled."""
    import matplotlib.pyplot as plt

    pae, chains = full["pae"], list(full["token_chain_ids"])
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    image = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=31.75,
                      interpolation="nearest")
    bounds, starts, start = [], [], 0
    for index in range(1, len(chains) + 1):
        if index == len(chains) or chains[index] != chains[index - 1]:
            starts.append((start, index, chains[index - 1]))
            if index < len(chains):
                bounds.append(index - 0.5)
            start = index
    for edge in bounds:
        ax.axhline(edge, color="#c44e52", lw=1.3, zorder=3)
        ax.axvline(edge, color="#c44e52", lw=1.3, zorder=3)
    centres = [(lo + hi - 1) / 2 for lo, hi, _ in starts]
    labels = [chain for _, _, chain in starts]
    ax.set_xticks(centres, labels)
    ax.set_yticks(centres, labels)
    ax.set_xlabel("Scored token (its position is judged)")
    ax.set_ylabel("Aligned token")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Expected position error (Å)",
                 fraction=0.046, pad=0.04)
    fig.tight_layout()
    return fig
