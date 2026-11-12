"""
Buffer recipe tests.

Like the bench tests, these guard a number that ends up on a balance. The
independent check is the laboratory's own PBS 20X protocol, written before
this module existed: 160 g NaCl, 4 g KCl, 4.8 g KH2PO4 and 35.6 g
Na2HPO4.2H2O per litre. The recipe file stores the 1X, a twentieth of that,
and preparing 1 L at 20X must give the protocol's masses back.

The molar masses are checked against sums worked by hand from the CIAAW
abridged weights in the module, written out in each test.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

import buffers as bf

ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------
# Formulas and molar masses
# --------------------------------------------------------------------------

def test_formula_parsing_handles_hydrates_parentheses_and_adducts():
    assert bf.parse_formula("NaCl") == {"Na": 1, "Cl": 1}
    assert bf.parse_formula("Na2HPO4.2H2O") == {"Na": 2, "H": 5, "P": 1,
                                                "O": 6}
    assert bf.parse_formula("Ca(OH)2") == {"Ca": 1, "O": 2, "H": 2}
    assert bf.parse_formula("C4H11NO3.HCl") == {"C": 4, "H": 12, "N": 1,
                                                "O": 3, "Cl": 1}
    assert bf.parse_formula("Na2CO3.10H2O")["H"] == 20


@pytest.mark.parametrize("formula", ["", "Xy2", "NaCl)", "Na(Cl", "nacl",
                                     "UO2"])
def test_a_formula_that_cannot_be_weighed_is_refused(formula):
    """Including a real element with no weight in the table (U): it is
    refused rather than given a weight from somewhere else."""
    with pytest.raises(ValueError):
        bf.parse_formula(formula)


def test_molar_masses_by_hand():
    # NaCl: 22.990 + 35.45
    assert bf.molar_mass("NaCl") == pytest.approx(58.440, abs=1e-9)
    # Na2HPO4: 2 x 22.990 + 1.0080 + 30.974 + 4 x 15.999 = 141.958
    assert bf.molar_mass("Na2HPO4") == pytest.approx(141.958, abs=1e-9)
    # Na2HPO4.2H2O: 141.958 + 2 x (2 x 1.0080 + 15.999) = 177.988
    assert bf.molar_mass("Na2HPO4.2H2O") == pytest.approx(177.988, abs=1e-9)
    # Tris base C4H11NO3: 48.044 + 11.088 + 14.007 + 47.997
    assert bf.molar_mass("C4H11NO3") == pytest.approx(121.136, abs=1e-9)


# --------------------------------------------------------------------------
# The shipped recipes
# --------------------------------------------------------------------------

def test_the_shipped_files_load_without_a_problem():
    book, problems = bf.load()
    assert problems == []
    assert book["recipes"], "no recipe in data/buffers"


def test_every_shipped_formula_has_a_molar_mass():
    book, _ = bf.load()
    for name, forms in book["compounds"].items():
        for form, formula in forms.items():
            if formula:
                assert bf.molar_mass(formula) > 0, (name, form)


def test_pbs_20x_gives_back_the_laboratory_protocol():
    book, _ = bf.load()
    plan = bf.prepare(book["recipes"]["PBS"], book["compounds"], 1000, 20)
    grams = {row["compound"]: row["amount"] for row in plan["rows"]}
    assert grams == pytest.approx({
        "Sodium chloride": 160.0,
        "Potassium chloride": 4.0,
        "Potassium dihydrogen phosphate": 4.8,
        "Disodium hydrogen phosphate": 35.6,
    })
    assert all(row["amount_unit"] == "g" for row in plan["rows"])


@pytest.mark.parametrize("recipe, volume_ml, expected", [
    # The laboratory's own document of 2026-11-12, amount by amount.
    ("TES (Tris-EDTA-sucrose)", 1000,
     {"Tris base": 24.22, "EDTA disodium": 0.19, "Sucrose": 171.15}),
    ("NaPi 50 mM pH 7", 1000,
     {"Sodium dihydrogen phosphate": 3.244,
      "Disodium hydrogen phosphate": 4.146}),
    ("NaPi 50 mM pH 8", 1000,
     {"Sodium dihydrogen phosphate": 3.244,
      "Disodium hydrogen phosphate": 4.146}),
    ("TB (base)", 4000,
     {"Yeast extract": 96.0, "Tryptone": 48.0, "Glycerol": 16.0}),
    ("TB-SALT", 1000, {"Potassium dihydrogen phosphate": 23.1,
                       "Dipotassium hydrogen phosphate": 125.4}),
    ("Glucose 20%", 1000, {"D-Glucose": 200.0}),
    ("Ampicillin 100 mg/mL", 26.3, {"Ampicillin sodium salt": 2.63}),
    ("2X-TY", 1000,
     {"Tryptone": 16.0, "Yeast extract": 10.0, "Sodium chloride": 5.0}),
    ("PEG-NaCl", 1000, {"Sodium chloride": 146.0, "PEG 6000": 200.0}),
])
def test_the_laboratory_protocols_come_back_out(recipe, volume_ml, expected):
    book, _ = bf.load()
    plan = bf.prepare(book["recipes"][recipe], book["compounds"], volume_ml)
    got = {row["compound"]: row["amount"] for row in plan["rows"]}
    assert got == pytest.approx(expected)


def test_a_molar_stock_is_weighed_in_the_form_the_recipe_names():
    """MgCl2 1 M is made from the hexahydrate: 0.05 L x 1 mol/L x
    203.295 g/mol = 10.165 g, the protocol's own figure. The IPTG protocol
    says 4.32 g for 18 mL, which is 1.007 M; at 1 M in 18 mL final volume
    the app says 4.289 g, and the 1 M is what the recipe declares."""
    book, _ = bf.load()
    mg = bf.prepare(book["recipes"]["MgCl2 1 M"], book["compounds"],
                    50)["rows"][0]
    assert mg["form"] == "hexahydrate"
    assert mg["amount"] == pytest.approx(10.165, abs=1e-3)
    iptg = bf.prepare(book["recipes"]["IPTG 1 M"], book["compounds"],
                      18)["rows"][0]
    assert iptg["amount"] == pytest.approx(4.289, abs=1e-3)


def test_another_hydrate_carries_the_same_phosphate():
    """1.78 g/L of the dihydrate, weighed as the anhydrous salt, is
    1.78 x 141.958 / 177.988 = 1.41968 g/L: the same moles of phosphate,
    without the water."""
    book, _ = bf.load()
    plan = bf.prepare(book["recipes"]["PBS"], book["compounds"], 1000, 1,
                      {"Disodium hydrogen phosphate": "anhydrous"})
    row = plan["rows"][3]
    assert row["form"] == "anhydrous"
    assert row["amount"] == pytest.approx(1.78 * 141.958 / 177.988)
    assert "dihydrate" in row["source"]


def test_the_ph_is_adjusted_with_what_the_recipe_names():
    """NaPi pH 8 starts near pH 7 and is raised with NaOH; the imidazole
    buffer and TES start above their target and are lowered with HCl."""
    book, _ = bf.load()
    for name, reagent in [("NaPi 50 mM pH 8", "NaOH"),
                          ("NaPi 50 mM imidazole 500 mM pH 8", "HCl"),
                          ("TES (Tris-EDTA-sucrose)", "HCl")]:
        assert f"Adjust the pH to 8 with {reagent}." in bf.steps(
            book["recipes"][name], 1000, 1)


def test_a_measured_ph_is_never_turned_into_an_adjustment():
    book, _ = bf.load()
    text = " ".join(bf.prepare(book["recipes"]["PBS"], book["compounds"],
                               500, 20)["steps"])
    assert "Do not adjust the pH" in text
    assert "Adjust the pH" not in text
    assert "dilute 20-fold" in text


# --------------------------------------------------------------------------
# Arithmetic on recipes written for the test
# --------------------------------------------------------------------------

COMPOUNDS = """compound,form,formula,cas,note
Sodium chloride,anhydrous,NaCl,,
Tris base,anhydrous,C4H11NO3,,
Glycerol,liquid,C3H8O3,,
Tryptone,powder,,,
"""


def _write(tmp_path, recipes, components, compounds=COMPOUNDS):
    (tmp_path / "compounds.csv").write_text(compounds)
    (tmp_path / "recipes.csv").write_text(
        ",".join(bf.RECIPE_COLUMNS) + "\n" + recipes)
    (tmp_path / "recipe_components.csv").write_text(
        ",".join(bf.COMPONENT_COLUMNS) + "\n" + components)
    return bf.load(tmp_path)


def test_molar_stock_and_volume_arithmetic(tmp_path):
    book, problems = _write(
        tmp_path,
        "T,buffer,8.0,adjust,HCl,autoclave,RT,,test,\n",
        "T,Sodium chloride,150,mM,,,\n"
        "T,Tris base,20,mM,,1 M,\n"
        "T,Glycerol,10,% v/v,,,\n"
        "T,Tryptone,10,g/L,,,\n")
    assert problems == []
    rows = {r["compound"]: r for r in bf.prepare(
        book["recipes"]["T"], book["compounds"], 500)["rows"]}
    # 0.150 mol/L x 0.5 L x 58.44 g/mol
    assert rows["Sodium chloride"]["amount"] == pytest.approx(4.383)
    # 20 mM from 1 M in 500 mL: 10 mL
    assert rows["Tris base"]["amount"] == pytest.approx(10.0)
    assert rows["Tris base"]["amount_unit"] == "mL"
    assert rows["Glycerol"]["amount"] == pytest.approx(50.0)
    assert rows["Tryptone"]["amount"] == pytest.approx(5.0)


def test_a_stock_too_dilute_for_the_factor_is_an_error_not_a_volume(tmp_path):
    book, _ = _write(tmp_path, "T,buffer,,,,none,,,test,\n",
                     "T,Tris base,20,mM,,1 M,\n")
    row = bf.prepare(book["recipes"]["T"], book["compounds"], 100,
                     60)["rows"][0]
    assert row["amount"] is None
    assert "stock is only" in row["error"]


def test_adjusting_a_concentrated_stock_says_where_the_ph_belongs(tmp_path):
    book, _ = _write(tmp_path, "T,buffer,8.0,adjust,HCl,none,,,test,\n",
                     "T,Sodium chloride,150,mM,,,\n")
    one = bf.steps(book["recipes"]["T"], 100, 1)
    ten = bf.steps(book["recipes"]["T"], 100, 10)
    assert "Adjust the pH to 8 with HCl." in one
    assert not any("1X dilution" in s for s in one)
    assert any("1X dilution" in s for s in ten)


@pytest.mark.parametrize("recipes, components, expected", [
    ("T,buffer,,,,none,,,EXAMPLE,\n", "T,Sodium chloride,1,mM,,,\n",
     "EXAMPLE"),
    ("T,buffer,,,,none,,,,\n", "T,Sodium chloride,1,mM,,,\n",
     "source is required"),
    ("T,buffer,7.4,,,none,,,s,\n", "T,Sodium chloride,1,mM,,,\n",
     "ph_mode"),
    ("T,buffer,7.4,adjust,,none,,,s,\n", "T,Sodium chloride,1,mM,,,\n",
     "ph_adjust_with is empty"),
    ("T,buffer,7.4,measured,HCl,none,,,s,\n", "T,Sodium chloride,1,mM,,,\n",
     "ph_adjust_with is filled"),
    ("T,buffer,7,4,measured,,none,,,s,\n", "T,Sodium chloride,1,mM,,,\n",
     "fields"),
    ("T,buffer,,,,none,,,s,\n", "T,Tryptone,5,mM,,,\n", "no formula"),
    ("T,buffer,,,,none,,,s,\n", "T,Natrium,5,mM,,,\n", "not in compounds"),
    ("T,buffer,,,,none,,,s,\n", "T,Sodium chloride,0,mM,,,\n",
     "positive number"),
    ("T,buffer,,,,none,,,s,\n", "T,Sodium chloride,5,mM,,50 % v/v,\n",
     "same kind of unit"),
    ("T,buffer,,,,none,,,s,\n", "T,Sodium chloride,5,mM,,1 mM,\n",
     "not more concentrated"),
])
def test_a_faulty_recipe_is_left_out_with_its_line(tmp_path, recipes,
                                                  components, expected):
    book, problems = _write(tmp_path, recipes, components)
    assert "T" not in book["recipes"]
    assert any(expected in p for p in problems), problems
    assert any(".csv:" in p for p in problems), "no file and line given"


def test_a_mass_of_a_salt_with_several_forms_needs_its_form(tmp_path):
    compounds = COMPOUNDS + ("Disodium hydrogen phosphate,anhydrous,Na2HPO4,,\n"
                             "Disodium hydrogen phosphate,dihydrate,"
                             "Na2HPO4.2H2O,,\n")
    book, problems = _write(tmp_path, "T,buffer,,,,none,,,s,\n",
                            "T,Disodium hydrogen phosphate,1.78,g/L,,,\n",
                            compounds)
    assert "T" not in book["recipes"]
    assert any("depends on the form" in p for p in problems)


def test_a_spreadsheet_semicolon_file_is_named_as_such(tmp_path):
    (tmp_path / "compounds.csv").write_text(
        COMPOUNDS.replace(",", ";"))
    (tmp_path / "recipes.csv").write_text(",".join(bf.RECIPE_COLUMNS) + "\n")
    (tmp_path / "recipe_components.csv").write_text(
        ",".join(bf.COMPONENT_COLUMNS) + "\n")
    _book, problems = bf.load(tmp_path)
    assert any("separated by ';'" in p for p in problems)


def test_the_checker_script_runs_the_app_loader():
    """`tools/check_buffers.py` must accept the shipped files and refuse a
    broken copy, through the same `buffers.load` the app uses."""
    script = ROOT / "tools" / "check_buffers.py"
    ok = subprocess.run([sys.executable, str(script)], capture_output=True,
                        text=True, cwd=ROOT)
    assert ok.returncode == 0, ok.stdout + ok.stderr


def test_the_checker_script_fails_on_a_problem(tmp_path):
    _write(tmp_path, "T,buffer,,,,none,,,,\n", "T,Sodium chloride,1,mM,,,\n")
    script = ROOT / "tools" / "check_buffers.py"
    bad = subprocess.run([sys.executable, str(script), str(tmp_path)],
                         capture_output=True, text=True, cwd=ROOT)
    assert bad.returncode == 1
    assert "source is required" in bad.stdout


def test_amounts_are_shown_in_a_unit_a_balance_reads():
    assert bf.format_amount(35.6, "g") == "35.60 g"
    assert bf.format_amount(4.8, "g") == "4.800 g"
    assert bf.format_amount(0.24, "g") == "240.0 mg"
    assert bf.format_amount(0.0004, "g") == "400.0 ug"
    assert bf.format_amount(10, "mL") == "10.00 mL"
    assert bf.format_amount(0.002, "mL") == "2.0 uL"


# --------------------------------------------------------------------------
# Composing a buffer in the app
# --------------------------------------------------------------------------

def test_a_composed_buffer_goes_through_the_same_arithmetic():
    """In 500 mL: 50 mM Tris is 0.025 mol, 0.025 x 121.136 = 3.0284 g, and
    150 mM NaCl is 0.075 mol, 0.075 x 58.44 = 4.383 g. The same numbers a
    recipe file would give, because `custom_recipe` builds the shape
    `prepare` already takes."""
    book, _ = bf.load()
    recipe = bf.custom_recipe(
        [{"compound": "Tris base", "amount": 50, "unit": "mM"},
         {"compound": "Sodium chloride", "amount": 150, "unit": "mM"}],
        ph=8.0, ph_adjust_with="HCl", sterilization="filter")
    plan = bf.prepare(recipe, book["compounds"], 500)
    got = {row["compound"]: row["amount"] for row in plan["rows"]}
    assert got == pytest.approx({"Tris base": 0.025 * 121.136,
                                 "Sodium chloride": 0.075 * 58.44})
    assert "Adjust the pH to 8 with HCl." in plan["steps"]
    assert "Filter-sterilize." in plan["steps"]


def test_a_composed_buffer_with_no_ph_is_never_told_to_adjust_one():
    book, _ = bf.load()
    recipe = bf.custom_recipe(
        [{"compound": "Sodium chloride", "amount": 1, "unit": "M"}])
    steps = bf.prepare(recipe, book["compounds"], 100)["steps"]
    assert not any("pH" in step for step in steps)
    assert recipe["ph"] is None and recipe["ph_mode"] == ""


def test_a_composed_buffer_picks_the_hydrate_it_is_given():
    """The whole point of the form selector: the same 50 mM of phosphate is
    4.146 g/L as Na2HPO4 and 5.198 g/L as Na2HPO4.2H2O, in one litre."""
    book, _ = bf.load()
    recipe = bf.custom_recipe(
        [{"compound": "Disodium hydrogen phosphate", "amount": 29.21,
          "unit": "mM"}])
    anhydrous = bf.prepare(recipe, book["compounds"], 1000, 1,
                           {"Disodium hydrogen phosphate": "anhydrous"})
    hydrate = bf.prepare(recipe, book["compounds"], 1000, 1,
                         {"Disodium hydrogen phosphate": "dihydrate"})
    assert anhydrous["rows"][0]["amount"] == pytest.approx(0.02921 * 141.958)
    assert hydrate["rows"][0]["amount"] == pytest.approx(0.02921 * 177.988)


def test_a_composed_buffer_is_not_written_anywhere():
    """It carries a source saying where it came from, and the files on disk
    are untouched: the recipe book is edited in the repository, not here."""
    before = (bf.DATA_DIR / "recipes.csv").read_text()
    recipe = bf.custom_recipe(
        [{"compound": "Sodium chloride", "amount": 1, "unit": "M"}])
    assert recipe["source"] == "composed in the app"
    assert (bf.DATA_DIR / "recipes.csv").read_text() == before
