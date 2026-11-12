"""
buffers - buffer and medium recipes, and the arithmetic of preparing them.

The recipes live in three CSV files under ``data/buffers/``, edited by hand
and pushed, like the sequence library: ``compounds.csv`` (one row per form of
a compound, with its formula), ``recipes.csv`` (one row per recipe) and
``recipe_components.csv`` (one row per component, at its 1X concentration).

No molar mass is stored anywhere. Each is derived from the formula with the
abridged standard atomic weights below, so the anhydrous and the dihydrate
forms of a salt cannot disagree with their own formulas, and a hydrate added
to the file needs nothing but its formula.

Importable from a notebook; Streamlit is not imported here.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

__version__ = "2026.11.12a"

DATA_DIR = Path(__file__).parent / "data" / "buffers"

# Abridged standard atomic weights, five significant figures, as published by
# the IUPAC Commission on Isotopic Abundances and Atomic Weights (CIAAW):
# https://www.ciaaw.org/abridged-atomic-weights.htm, based on T. Prohaska et
# al., "Standard atomic weights of the elements 2021 (IUPAC Technical
# Report)", Pure Appl. Chem. 94 (2022), doi:10.1515/pac-2019-0603. Values
# copied from that page on 2026-11-12, elements 1 to 56 except Tc (no
# standard atomic weight) and the noble gases, which no buffer contains.
# Only elements in this table can appear in a formula: a symbol outside it is
# reported as a problem rather than given a weight from somewhere else.
ATOMIC_WEIGHTS = {
    "H": 1.0080, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
    "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693,
    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922,
    "Se": 78.971, "Br": 79.904, "Rb": 85.468, "Sr": 87.62, "Y": 88.906,
    "Zr": 91.222, "Nb": 92.906, "Mo": 95.95, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Cs": 132.91, "Ba": 137.33,
}
ATOMIC_WEIGHTS_SOURCE = (
    "CIAAW abridged standard atomic weights, based on Atomic Weights 2021 "
    "(Prohaska et al., Pure Appl. Chem. 2022, doi:10.1515/pac-2019-0603)")

MOLAR = {"M": 1.0, "mM": 1e-3, "uM": 1e-6}            # to mol/L
MASS = {"g/L": 1.0, "mg/mL": 1.0, "ug/mL": 1e-3, "% w/v": 10.0}  # to g/L
VOLUME = {"% v/v": 1.0}                                # to % v/v
UNITS = {**MOLAR, **MASS, **VOLUME}
KINDS = ("buffer", "medium", "solution")
STERILIZATION = ("filter", "autoclave", "none")
PH_MODES = ("adjust", "measured")

COMPOUND_COLUMNS = ["compound", "form", "formula", "cas", "note"]
RECIPE_COLUMNS = ["recipe", "kind", "ph", "ph_mode", "ph_adjust_with",
                  "sterilization", "storage", "usual_stock", "source",
                  "notes"]
COMPONENT_COLUMNS = ["recipe", "compound", "amount", "unit", "form", "stock",
                     "note"]


# --------------------------------------------------------------------------
# Formulas
def parse_formula(formula: str) -> dict[str, int]:
    """Element counts of ``NaH2PO4.2H2O``, ``Ca(OH)2`` or ``C4H11NO3.HCl``.

    Parts separated by a dot are added, each with its leading multiplier.
    Raises ValueError on anything unreadable, and on any element that has no
    weight in ATOMIC_WEIGHTS, so that a formula is either fully weighable or
    refused.
    """
    if not formula or not formula.strip():
        raise ValueError("empty formula")
    total: dict[str, int] = {}
    for part in formula.strip().split("."):
        match = re.fullmatch(r"(\d*)(.+)", part.strip())
        if not match:
            raise ValueError(f"empty part in {formula!r}")
        factor = int(match.group(1) or 1)
        for element, count in _parse_group(match.group(2), formula).items():
            total[element] = total.get(element, 0) + factor * count
    return total


def _parse_group(text: str, whole: str) -> dict[str, int]:
    stack: list[dict[str, int]] = [{}]
    i = 0
    while i < len(text):
        char = text[i]
        if char == "(":
            stack.append({})
            i += 1
        elif char == ")":
            if len(stack) == 1:
                raise ValueError(f"unbalanced ')' in {whole!r}")
            i += 1
            digits = re.match(r"\d*", text[i:]).group()
            i += len(digits)
            group = stack.pop()
            for element, count in group.items():
                stack[-1][element] = (stack[-1].get(element, 0)
                                      + count * int(digits or 1))
        else:
            match = re.match(r"([A-Z][a-z]?)(\d*)", text[i:])
            if not match:
                raise ValueError(f"cannot read {text[i:]!r} in {whole!r}")
            element = match.group(1)
            if element not in ATOMIC_WEIGHTS:
                raise ValueError(
                    f"{element!r} in {whole!r} is not an element with a "
                    f"standard atomic weight in this app's table")
            stack[-1][element] = (stack[-1].get(element, 0)
                                  + int(match.group(2) or 1))
            i += len(match.group())
    if len(stack) != 1:
        raise ValueError(f"unbalanced '(' in {whole!r}")
    return stack[0]


def molar_mass(formula: str) -> float:
    """g/mol from the formula and ATOMIC_WEIGHTS. Derived, not measured."""
    return sum(ATOMIC_WEIGHTS[el] * n for el, n in parse_formula(formula).items())


# --------------------------------------------------------------------------
# Units
def dimension(unit: str) -> str:
    return "molar" if unit in MOLAR else "mass" if unit in MASS else "volume"


def parse_quantity(text: str) -> tuple[float, str]:
    """``"1 M"`` -> ``(1.0, "M")``. Raises ValueError."""
    match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*(.+)", (text or "").strip())
    if not match or match.group(2).strip() not in UNITS:
        raise ValueError(f"cannot read {text!r} as a number and one of "
                         f"{', '.join(UNITS)}")
    return float(match.group(1)), match.group(2).strip()


def _base(value: float, unit: str) -> float:
    """mol/L, g/L or % v/v, so values of one dimension compare."""
    return value * UNITS[unit]


# --------------------------------------------------------------------------
# Loading
def _rows(path: Path, columns: list[str], problems: list[str]):
    """(line number, row dict) for each data line; header checked."""
    if not path.exists():
        problems.append(f"{path.name}: file missing")
        return []
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    numbered = [(n, line) for n, line in enumerate(lines, 1)
                if line.strip() and not line.lstrip().startswith("#")]
    if not numbered:
        problems.append(f"{path.name}: no header line")
        return []
    n0, header_line = numbered[0]
    if ";" in header_line and "," not in header_line:
        problems.append(
            f"{path.name}:{n0}: separated by ';', probably saved by a "
            "spreadsheet in a French or Belgian locale. Save it as CSV with ',' "
            "between fields and '.' as decimal mark.")
        return []
    header = [h.strip() for h in next(csv.reader([header_line]))]
    if header != columns:
        problems.append(f"{path.name}:{n0}: header is {header}, expected "
                        f"{columns}")
        return []
    out = []
    for n, line in numbered[1:]:
        values = next(csv.reader([line]))
        if len(values) != len(columns):
            problems.append(f"{path.name}:{n}: {len(values)} fields, "
                            f"expected {len(columns)}; row left out")
            continue
        out.append((n, {k: v.strip() for k, v in zip(columns, values)}))
    return out


def load(folder: str | Path | None = None) -> tuple[dict, list[str]]:
    """Read and check the three files.

    Returns ``(book, problems)``. ``book`` has ``"compounds"``,
    ``{name: {form: formula or ""}}`` in file order, and ``"recipes"``,
    ``{name: recipe dict with a "components" list}``. A row with a problem
    is left out and the problem names its file and line, like the sequence
    library; a recipe with any faulty component is left out whole, because
    a buffer missing one salt is a different buffer.
    """
    folder = Path(folder) if folder else DATA_DIR
    problems: list[str] = []

    compounds: dict[str, dict[str, str]] = {}
    for n, row in _rows(folder / "compounds.csv", COMPOUND_COLUMNS, problems):
        where = f"compounds.csv:{n}"
        if not row["compound"] or not row["form"]:
            problems.append(f"{where}: compound and form are both required")
            continue
        forms = compounds.setdefault(row["compound"], {})
        if row["form"] in forms:
            problems.append(f"{where}: {row['compound']} / {row['form']} is "
                            "listed twice")
            continue
        if row["formula"]:
            try:
                parse_formula(row["formula"])
            except ValueError as exc:
                problems.append(f"{where}: {exc}; row left out")
                continue
        forms[row["form"]] = row["formula"]
    for name, forms in list(compounds.items()):
        if not forms:
            del compounds[name]
            continue
        with_formula = [f for f, formula in forms.items() if formula]
        if with_formula and len(with_formula) != len(forms):
            problems.append(f"compounds.csv: {name!r} has forms with and "
                            "without a formula; give all or none. Left out.")
            del compounds[name]

    recipes: dict[str, dict] = {}
    lines: dict[str, int] = {}
    for n, row in _rows(folder / "recipes.csv", RECIPE_COLUMNS, problems):
        where = f"recipes.csv:{n}"
        name = row["recipe"]
        if not name:
            problems.append(f"{where}: recipe name missing")
            continue
        if name in lines:
            problems.append(f"{where}: {name!r} already on line {lines[name]}")
            continue
        lines[name] = n
        bad = []
        if row["kind"] not in KINDS:
            bad.append(f"kind {row['kind']!r}, expected one of {KINDS}")
        if row["sterilization"] not in STERILIZATION:
            bad.append(f"sterilization {row['sterilization']!r}, expected one "
                       f"of {STERILIZATION}")
        ph = None
        if row["ph"]:
            try:
                ph = float(row["ph"])
                if not 0 <= ph <= 14:
                    raise ValueError
            except ValueError:
                bad.append(f"pH {row['ph']!r} is not a number between 0 and "
                           "14 (decimal mark is '.')")
            if row["ph_mode"] not in PH_MODES:
                bad.append(f"ph_mode {row['ph_mode']!r}, expected one of "
                           f"{PH_MODES}")
            elif row["ph_mode"] == "adjust" and not row["ph_adjust_with"]:
                bad.append("ph_mode is adjust but ph_adjust_with is empty")
            elif row["ph_mode"] == "measured" and row["ph_adjust_with"]:
                bad.append("ph_mode is measured, so nothing is adjusted, but "
                           "ph_adjust_with is filled")
        elif row["ph_mode"] or row["ph_adjust_with"]:
            bad.append("ph_mode or ph_adjust_with given without a pH")
        usual = 1.0
        if row["usual_stock"]:
            try:
                usual = float(row["usual_stock"])
                if usual < 1:
                    raise ValueError
            except ValueError:
                bad.append(f"usual_stock {row['usual_stock']!r} is not a "
                           "number of 1 or more")
        if not row["source"]:
            bad.append("source is required")
        elif row["source"].upper() == "EXAMPLE":
            bad.append("the EXAMPLE recipe is a format illustration")
        if bad:
            problems.append(f"{where}: {name!r} left out: " + "; ".join(bad))
            continue
        recipes[name] = {**row, "ph": ph, "usual_stock": usual,
                         "components": [], "line": n}

    broken: set[str] = set()
    seen: dict[tuple[str, str], int] = {}
    for n, row in _rows(folder / "recipe_components.csv", COMPONENT_COLUMNS,
                        problems):
        where = f"recipe_components.csv:{n}"
        recipe, compound = row["recipe"], row["compound"]
        if recipe not in lines:
            problems.append(f"{where}: recipe {recipe!r} is not in "
                            "recipes.csv")
            continue
        if recipe not in recipes:
            continue        # the recipe itself was refused, already reported
        key = (recipe, compound)
        if key in seen:
            problems.append(f"{where}: {compound!r} already in {recipe!r} on "
                            f"line {seen[key]}")
            broken.add(recipe)
            continue
        seen[key] = n
        bad = _component_problem(row, compounds)
        if bad:
            problems.append(f"{where}: {bad}")
            broken.add(recipe)
            continue
        amount = float(row["amount"])
        stock = parse_quantity(row["stock"]) if row["stock"] else None
        recipes[recipe]["components"].append(
            {"compound": compound, "amount": amount, "unit": row["unit"],
             "form": row["form"], "stock": stock, "note": row["note"],
             "line": n})

    for name in list(recipes):
        if name in broken:
            problems.append(f"recipes.csv:{recipes[name]['line']}: {name!r} "
                            "left out because a component has a problem")
            del recipes[name]
        elif not recipes[name]["components"]:
            problems.append(f"recipes.csv:{recipes[name]['line']}: {name!r} "
                            "has no component in recipe_components.csv")
            del recipes[name]

    return {"compounds": compounds, "recipes": recipes}, problems


def _component_problem(row: dict, compounds: dict) -> str | None:
    compound = row["compound"]
    if compound not in compounds:
        return (f"compound {compound!r} is not in compounds.csv (names must "
                "match exactly)")
    forms = compounds[compound]
    has_formula = any(forms.values())
    try:
        amount = float(row["amount"])
        if amount <= 0:
            raise ValueError
    except ValueError:
        return (f"amount {row['amount']!r} is not a positive number "
                "(decimal mark is '.')")
    unit = row["unit"]
    if unit not in UNITS:
        return f"unit {unit!r}, expected one of {', '.join(UNITS)}"
    if unit in MOLAR and not has_formula:
        return f"{compound!r} has no formula, so it cannot be given in {unit}"
    if row["form"] and row["form"] not in forms:
        return (f"form {row['form']!r} is not listed for {compound!r} "
                f"(listed: {', '.join(forms)})")
    if unit in MASS and len(forms) > 1 and not row["form"]:
        return (f"{amount:g} {unit} of {compound!r} depends on the form, "
                "and it has several: say which one the amount refers to")
    if row["stock"]:
        try:
            stock, stock_unit = parse_quantity(row["stock"])
        except ValueError as exc:
            return f"stock: {exc}"
        if dimension(stock_unit) != dimension(unit):
            return (f"stock in {stock_unit} but amount in {unit}; use the "
                    "same kind of unit")
        if _base(stock, stock_unit) <= _base(amount, unit):
            return (f"stock {row['stock']} is not more concentrated than "
                    f"{amount:g} {unit}")
    return None


def forms_offered(compounds: dict, component: dict) -> list[str]:
    """The forms one may weigh for a component, recipe form first.

    A compound without formulas cannot be converted between forms, so only
    the form the recipe names (or its only one) is offered.
    """
    forms = compounds[component["compound"]]
    if not any(forms.values()):
        return [component["form"]] if component["form"] else list(forms)[:1]
    names = list(forms)
    if component["form"]:
        names.remove(component["form"])
        names.insert(0, component["form"])
    return names


# --------------------------------------------------------------------------
# Preparing
def prepare(recipe: dict, compounds: dict, volume_ml: float,
            factor: float = 1.0, forms: dict[str, str] | None = None) -> dict:
    """What to weigh or pipette for ``volume_ml`` of the recipe at ``factor``X.

    ``forms`` maps a compound to the form actually weighed; a compound left
    out gets the first of ``forms_offered``. Returns ``{"rows": [...],
    "steps": [...]}``. Each row carries ``amount`` and ``amount_unit`` ("g",
    "mL") or an ``error`` when the row cannot be made (a stock too dilute for
    the factor). Arithmetic only: whether that much dissolves is not known
    here.
    """
    if volume_ml <= 0:
        raise ValueError("the volume must be above zero")
    if factor <= 0:
        raise ValueError("the concentration factor must be above zero")
    forms = forms or {}
    litres = volume_ml / 1000.0
    rows = []
    for comp in recipe["components"]:
        offered = forms_offered(compounds, comp)
        form = forms.get(comp["compound"], offered[0])
        if form not in offered:
            raise ValueError(f"{form!r} is not a form offered for "
                             f"{comp['compound']!r}")
        formula = compounds[comp["compound"]][form]
        mw = molar_mass(formula) if formula else None
        unit = comp["unit"]
        target = comp["amount"] * factor
        row = {"compound": comp["compound"], "form": form,
               "formula": formula, "molar_mass": mw,
               "amount_1x": comp["amount"], "unit": unit,
               "amount_here": target, "source": "", "amount": None,
               "amount_unit": "", "error": "", "note": comp["note"]}
        if comp["stock"]:
            stock, stock_unit = comp["stock"]
            row["source"] = f"from a {stock:g} {stock_unit} stock"
            need, have = _base(target, unit), _base(stock, stock_unit)
            if need >= have:
                row["error"] = (f"needs {target:g} {unit}, the stock is only "
                                f"{stock:g} {stock_unit}: weigh it instead, "
                                "or lower the factor")
            else:
                row["amount"], row["amount_unit"] = need / have * volume_ml, "mL"
        elif unit in MOLAR:
            row["source"] = "weigh"
            row["amount"] = _base(target, unit) * litres * mw
            row["amount_unit"] = "g"
        elif unit in MASS:
            row["source"] = "weigh"
            grams = _base(target, unit) * litres
            if comp["form"] and form != comp["form"]:
                # Same moles of the compound, carried by a different form.
                grams *= mw / molar_mass(compounds[comp["compound"]][comp["form"]])
                row["source"] = f"weigh; recipe gives it as {comp['form']}"
            row["amount"], row["amount_unit"] = grams, "g"
        else:
            row["source"] = "pipette"
            row["amount"] = target / 100.0 * volume_ml
            row["amount_unit"] = "mL"
        rows.append(row)
    return {"rows": rows, "steps": steps(recipe, volume_ml, factor)}


def steps(recipe: dict, volume_ml: float, factor: float) -> list[str]:
    """The order of work, written from what the recipe declares."""
    out = ["Weigh or pipette the components into less water than the final "
           "volume and dissolve them."]
    ph, mode = recipe["ph"], recipe["ph_mode"]
    if ph is not None and mode == "measured":
        out.append(
            f"Do not adjust the pH. The 1X buffer reads about pH {ph:g} "
            f"without adjustment, according to the recipe's source"
            + (f"; this {factor:g}X stock may read "
               f"differently until diluted." if factor != 1 else "."))
    elif ph is not None:
        out.append(f"Adjust the pH to {ph:g} with {recipe['ph_adjust_with']}.")
        if factor != 1:
            out.append(
                f"That pH is the 1X buffer's. A concentrated stock generally "
                f"does not read the pH it gives once diluted, so check the pH "
                f"on a 1X dilution of this {factor:g}X stock rather than "
                f"assuming it.")
    out.append(f"Bring to {volume_ml:g} mL with water.")
    if recipe["sterilization"] == "filter":
        out.append("Filter-sterilize.")
    elif recipe["sterilization"] == "autoclave":
        out.append("Autoclave.")
    if recipe["storage"]:
        out.append(f"Store: {recipe['storage']}.")
    if factor != 1:
        out.append(f"To use, dilute {factor:g}-fold: 1 volume of this stock "
                   f"made up to {factor:g} volumes with water.")
    return out


def format_amount(amount: float | None, unit: str) -> str:
    """Grams and millilitres in a unit a balance or a pipette reads."""
    if amount is None:
        return ""
    if unit == "g":
        if amount >= 10:
            return f"{amount:.2f} g"
        if amount >= 1:
            return f"{amount:.3f} g"
        if amount >= 1e-3:
            return f"{amount * 1e3:.1f} mg"
        return f"{amount * 1e6:.1f} ug"
    if amount >= 1:
        return f"{amount:.2f} mL"
    return f"{amount * 1e3:.1f} uL"
