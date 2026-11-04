"""
Reading AlphaFold3 archives: ipTM, PAE and which model is which.

Every archive here is built in the test from synthetic coordinates and invented
scores. A real AlphaFold Server zip was used to check the layout these mimic,
and is deliberately not part of the repository: it holds the sequences of the
job, and the repository is public.
"""

from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBParser

import af3
import structure as stx

THREE = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU",
         "K": "LYS", "E": "GLU", "Y": "TYR"}


def two_chain_cif(len_a: int = 12, len_b: int = 8) -> bytes:
    """A CA-only two-chain model, written as mmCIF the way AlphaFold3 does."""
    lines, serial = [], 0
    for chain, length, offset in (("A", len_a, 0.0), ("B", len_b, 6.0)):
        for index in range(length):
            serial += 1
            lines.append(
                f"ATOM  {serial:5d}  CA  {THREE['AGSVLKEY'[index % 8]]} "
                f"{chain}{index + 1:4d}    {index * 3.8:8.3f}{offset:8.3f}"
                f"{0.0:8.3f}  1.00 90.00           C")
    pdb = "\n".join(lines) + "\nTER\nEND\n"
    structure = PDBParser(QUIET=True).get_structure("m", io.StringIO(pdb))
    out = io.StringIO()
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(out)
    return out.getvalue().encode()


def summary(iptm: float, ptm: float = 0.8, ranking: float | None = None,
            with_chain_ids: bool = True, clash: bool = False) -> dict:
    data = {
        "chain_iptm": [iptm, iptm], "chain_pair_iptm": [[0.9, iptm], [iptm, 0.85]],
        "chain_pair_pae_min": [[0.8, 3.1], [2.7, 0.8]], "chain_ptm": [0.9, 0.85],
        "fraction_disordered": 0.0, "has_clash": 1.0 if clash else 0.0,
        "iptm": iptm, "ptm": ptm,
        "ranking_score": (0.8 * iptm + 0.2 * ptm) if ranking is None else ranking,
    }
    if with_chain_ids:
        data["chain_ids"] = ["A"] * 12 + ["B"] * 8
    return data


def full_data(len_a: int = 12, len_b: int = 8) -> dict:
    n = len_a + len_b
    pae = np.full((n, n), 20.0)
    pae[:len_a, :len_a] = 1.0
    pae[len_a:, len_a:] = 1.5
    pae[2, len_a + 1] = 3.0           # one confident inter-chain pair
    return {"pae": pae.tolist(), "token_chain_ids": ["A"] * len_a + ["B"] * len_b,
            "token_res_ids": list(range(1, len_a + 1)) + list(range(1, len_b + 1)),
            "atom_plddts": [], "contact_probs": []}


def server_zip(iptms=(0.85, 0.62, 0.4), extra: dict | None = None) -> bytes:
    """The AlphaFold Server layout: fold_<job>_model_<n>.cif and friends."""
    cif = two_chain_cif()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        prefix = "fold_test_job"
        zf.writestr(f"{prefix}_job_request.json", json.dumps(
            [{"name": "test job", "dialect": "alphafoldserver"}]))
        # Written out of rank order on purpose: the reader has to sort.
        for n, value in reversed(list(enumerate(iptms))):
            zf.writestr(f"{prefix}_model_{n}.cif", cif)
            zf.writestr(f"{prefix}_summary_confidences_{n}.json",
                        json.dumps(summary(value)))
            zf.writestr(f"{prefix}_full_data_{n}.json", json.dumps(full_data()))
        zf.writestr("templates/fold_test_job_template_hit_0_chains_a.cif", cif)
        zf.writestr("msas/fold_test_job_paired_msa_chains_a.a3m", ">q\nAAAA\n")
        zf.writestr("terms_of_use.md", "terms")
        for name, content in (extra or {}).items():
            zf.writestr(name, content)
    return buffer.getvalue()


def local_zip() -> bytes:
    """AlphaFold3 run locally, then zipped: top model plus seed/sample folders."""
    cif = two_chain_cif()
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        zf.writestr("myjob/myjob_model.cif", cif)
        zf.writestr("myjob/myjob_summary_confidences.json",
                    json.dumps(summary(0.7, with_chain_ids=False)))
        zf.writestr("myjob/myjob_confidences.json", json.dumps(full_data()))
        zf.writestr("myjob/seed-1_sample-0/model.cif", cif)
        zf.writestr("myjob/seed-1_sample-0/summary_confidences.json",
                    json.dumps(summary(0.5, with_chain_ids=False)))
        zf.writestr("myjob/seed-1_sample-0/confidences.json",
                    json.dumps(full_data()))
    return buffer.getvalue()


# --------------------------------------------------------------------------- #

def test_server_archive_is_indexed_and_ranked_best_first():
    bundle = af3.read_af3_zip(server_zip(), "fold_test_job.zip")
    assert bundle.source == "AlphaFold Server"
    assert bundle.job_name == "test job"
    assert [s.label for s in bundle.samples] == ["model 0", "model 1", "model 2"]
    assert [s.iptm for s in bundle.samples] == [0.85, 0.62, 0.4]
    assert all(s.full_member.endswith(f"full_data_{s.label[-1]}.json")
               for s in bundle.samples)


def test_templates_and_msas_are_never_taken_for_models():
    bundle = af3.read_af3_zip(server_zip(), "x.zip")
    members = [s.cif_member for s in bundle.samples]
    assert not any(m.startswith(("templates/", "msas/")) for m in members)


def test_local_layout_pairs_files_by_folder():
    bundle = af3.read_af3_zip(local_zip(), "myjob.zip")
    labels = {s.label: s for s in bundle.samples}
    assert set(labels) == {"top-ranked model", "seed-1_sample-0"}
    top = labels["top-ranked model"]
    assert top.summary_member == "myjob/myjob_summary_confidences.json"
    assert top.full_member == "myjob/myjob_confidences.json"
    seed = labels["seed-1_sample-0"]
    assert seed.full_member == "myjob/seed-1_sample-0/confidences.json"
    assert bundle.samples[0] is top        # 0.7 ranks above 0.5


def test_the_model_read_from_the_zip_is_a_usable_structure():
    raw = server_zip()
    sample = af3.read_af3_zip(raw, "x.zip").samples[0]
    model = stx.load_structure(af3.read_member(raw, sample.cif_member),
                               "model_0.cif")
    assert list(stx.chain_table(model)["chain"]) == ["A", "B"]


def test_pair_iptm_is_read_for_the_named_pair_in_either_order():
    s = summary(0.55)
    order, source = af3.chain_order(s)
    assert order == ["A", "B"] and "chain_ids" in source
    assert af3.pair_iptm(s, order, "B", "A") == 0.55
    assert af3.pair_iptm(s, order, "A", "B") == 0.55
    assert af3.pair_iptm(s, order, "A", "Z") is None
    # The minimum over both directions, since PAE is not symmetric.
    assert af3.pair_pae_min(s, order, "A", "B") == 2.7


def test_chain_order_falls_back_and_says_where_it_came_from():
    s = summary(0.5, with_chain_ids=False)
    order, source = af3.chain_order(s, full_data())
    assert order == ["A", "B"] and "token_chain_ids" in source
    order, source = af3.chain_order(s, None, ["A", "B"])
    assert order == ["A", "B"] and ".cif" in source
    # A count that does not match the matrix is refused, not used.
    order, source = af3.chain_order(s, None, ["A", "B", "C"])
    assert order == [] and "does not match" in source


@pytest.mark.parametrize("value, level", [
    (0.81, "high"), (0.8, "grey"), (0.6, "grey"), (0.59, "low"), (None, "none")])
def test_iptm_thresholds(value, level):
    assert af3.iptm_verdict(value)[0] == level


def test_interface_pae_sees_the_one_confident_pair():
    full = af3.load_full_data(server_zip(), "fold_test_job_full_data_0.json")
    block = af3.interface_pae(full, "A", "B")
    assert block["min"] == 3.0
    assert block["mean"] > 19
    assert 0 < block["fraction_close"] < 0.01
    assert af3.interface_pae(full, "A", "Z") is None


def test_samples_table_has_one_row_per_model_with_the_pair():
    bundle = af3.read_af3_zip(server_zip(), "x.zip")
    table = af3.samples_table(bundle, "A", "B")
    assert list(table["model"]) == ["model 0", "model 1", "model 2"]
    assert list(table["pair_iptm"]) == [0.85, 0.62, 0.4]


def test_plot_pae_draws():
    import matplotlib.pyplot as plt
    full = af3.load_full_data(server_zip(), "fold_test_job_full_data_0.json")
    fig = af3.plot_pae(full)
    assert fig.axes
    plt.close(fig)


# ---- refusals ---------------------------------------------------------------

def test_a_zip_without_models_is_refused_with_a_reason():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        zf.writestr("notes.txt", "hello")
    with pytest.raises(af3.AF3Error, match="no AlphaFold3 model"):
        af3.read_af3_zip(buffer.getvalue(), "notes.zip")


def test_models_without_summaries_are_refused_with_a_reason():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        zf.writestr("x_model_0.cif", two_chain_cif())
    with pytest.raises(af3.AF3Error, match="no `summary_confidences`"):
        af3.read_af3_zip(buffer.getvalue(), "x.zip")


def test_not_a_zip_is_refused_and_is_a_structure_error():
    with pytest.raises(stx.StructureError):
        af3.read_af3_zip(b"PK\x03\x04garbage", "bad.zip")


def test_an_oversized_member_is_refused_before_it_is_read(monkeypatch):
    monkeypatch.setattr(af3, "MAX_MEMBER_BYTES", 100)
    with pytest.raises(af3.AF3Error, match="once unpacked"):
        af3.read_af3_zip(server_zip(), "x.zip")


def test_a_pae_matrix_that_does_not_match_its_tokens_is_refused():
    bad = full_data()
    bad["token_chain_ids"] = bad["token_chain_ids"][:-1]
    raw = server_zip(extra={"broken_full_data.json": json.dumps(bad)})
    with pytest.raises(af3.AF3Error, match="do not match"):
        af3.load_full_data(raw, "broken_full_data.json")


def test_is_zip():
    assert af3.is_zip(server_zip())
    assert not af3.is_zip(two_chain_cif())
