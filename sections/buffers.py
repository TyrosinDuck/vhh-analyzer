"""
sections.buffers - what to weigh for a buffer or a medium, at any volume and
concentration factor, with the form of each salt you actually have.

The screen draws. The recipes are read by `compute.load_buffers_cached`, and
the arithmetic, the formulas and the atomic weights live in `buffers.py`. A
buffer composed on this screen is turned into a recipe of the same shape and
goes through the same `prepare`, so there is one arithmetic and not two.
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

import buffers as bf
from compute import load_buffers_cached
from shell import FIT, SECTIONS, Settings

__version__ = "2026.11.12c"

V_RECIPE = "📖 From a recipe"
V_CUSTOM = "🧪 Compose one"
VIEWS = [V_RECIPE, V_CUSTOM]

# The units offered when composing a buffer, molar first: a molar
# concentration survives a change of hydrate, a mass does not.
CUSTOM_UNITS = ["mM", "M", "uM", "g/L", "mg/mL", "ug/mL", "% w/v", "% v/v"]


def render(settings: Settings) -> None:
    """Draw the Buffers & media screen (key "buffers").

    The router gives every screen the same argument. This one reads
    none of it: see `shell.Settings` for why that is cheaper than a
    second list of the screens.
    """
    icon, label, *_ = SECTIONS["buffers"]
    st.title(f"{icon} {label}")
    st.markdown(
        "What to weigh for a buffer or a medium, at the volume and the "
        "concentration you need, with the form of each salt that is on your "
        "shelf: the mass follows from the formula of that form. Take one of "
        "the laboratory's recipes, or compose a buffer from the compounds.")

    book, problems = load_buffers_cached()
    for message in problems:
        st.warning(message, icon="⚠️")

    view = st.radio("View", VIEWS, horizontal=True,
                    label_visibility="collapsed", key="buf_view")
    st.divider()

    if view == V_CUSTOM:
        _custom(book)
    else:
        _from_recipe(book)


# --------------------------------------------------------------------------
def _from_recipe(book: dict) -> None:
    st.info(
        "**Read-only here, and public.** The recipes are the files in "
        "`data/buffers/` of the repository, edited there and pushed, like the "
        "sequence library. They are one laboratory's recipes, each with the "
        "source it gives, not a reference collection.", icon="🧂")

    recipes = book["recipes"]
    if not recipes:
        st.info("No recipe could be loaded. Add rows to the files in "
                "`data/buffers/` and check them with "
                "`python3 tools/check_buffers.py`. You can still compose a "
                "buffer in the other view.")
        return

    top = st.columns([2, 1, 1])
    with top[0]:
        name = st.selectbox("Recipe", list(recipes), key="buf_recipe")
    recipe = recipes[name]
    with top[1]:
        volume = st.number_input("Final volume (mL)", 0.1, 1e6, 1000.0, 50.0,
                                 key=f"buf_volume_{name}")
    with top[2]:
        factor = st.number_input(
            "Concentration (X)", 0.1, 1000.0, float(recipe["usual_stock"]),
            1.0, key=f"buf_factor_{name}",
            help="1 for the buffer as used. The recipe's usual stock factor "
                 "is the starting value.")

    st.caption(
        f"{recipe['kind'].capitalize()} · source: {recipe['source']}"
        + (f" · usually made at {recipe['usual_stock']:g}X"
           if recipe["usual_stock"] != 1 else ""))
    if recipe["notes"]:
        st.markdown(f"*{recipe['notes']}*")

    # One choice per compound that exists in more than one form. Keyed on the
    # recipe and the compound, so a form chosen for PBS is not carried into
    # another recipe that uses the same salt as a different form.
    forms: dict[str, str] = {}
    choices = [(c, bf.forms_offered(book["compounds"], c))
               for c in recipe["components"]]
    choices = [(c, offered) for c, offered in choices if len(offered) > 1]
    if choices:
        st.markdown("**The form you weigh**")
        cols = st.columns(3)
        for i, (comp, offered) in enumerate(choices):
            with cols[i % len(cols)]:
                forms[comp["compound"]] = st.selectbox(
                    comp["compound"], offered,
                    format_func=lambda f, c=comp: (
                        f"{f} ({book['compounds'][c['compound']][f]})"),
                    key=f"buf_form_{name}_{comp['compound']}")

    _results(book, recipe, volume, factor, forms,
             f"{name.replace(' ', '_')}_{factor:g}X_{volume:g}mL.csv",
             "buf_csv")


# --------------------------------------------------------------------------
def _custom(book: dict) -> None:
    st.markdown(
        "Pick the compounds, give each a concentration, and say how much you "
        "are making. Nothing composed here is saved: the app writes nothing, "
        "and a buffer you come back to belongs in `data/buffers/` with its "
        "source, where the other view reads it.")

    compounds = book["compounds"]
    if not compounds:
        st.info("No compound could be loaded from `data/buffers/"
                "compounds.csv`.")
        return

    picked = st.multiselect(
        "Compounds", list(compounds), key="cus_compounds",
        help="Only compounds listed in data/buffers/compounds.csv, because "
             "a molar mass is only as good as the formula behind it. Add "
             "what is missing to that file.")
    if not picked:
        st.info("Pick at least one compound.")
        return

    head = st.columns([1, 1, 1])
    with head[0]:
        volume = st.number_input("Final volume (mL)", 0.1, 1e6, 1000.0, 50.0,
                                 key="cus_volume")
    with head[1]:
        ph = st.number_input(
            "pH to adjust to (0 = none)", 0.0, 14.0, 0.0, 0.1, format="%.2f",
            key="cus_ph",
            help="Only a pH you set yourself with the reagent beside it. "
                 "This screen never computes a pH or an acid/base ratio.")
    with head[2]:
        reagent = st.text_input("Adjust with", "HCl", key="cus_reagent",
                                disabled=ph == 0)

    sterilization = st.radio(
        "Sterilization", list(bf.STERILIZATION), horizontal=True,
        key="cus_sterilization",
        format_func=lambda s: {"filter": "Filter", "autoclave": "Autoclave",
                               "none": "None"}[s])

    st.markdown("**Each compound**")
    components = []
    for compound in picked:
        forms = list(compounds[compound])
        weighable = any(compounds[compound].values())
        row = st.columns([2, 1, 1, 2])
        with row[0]:
            st.markdown(f"**{compound}**")
            if not weighable:
                st.caption("No formula: mass and % units only.")
        with row[1]:
            amount = st.number_input(
                "Concentration", 0.0, 1e6, 0.0, 1.0, format="%.4f",
                key=f"cus_amount_{compound}", label_visibility="collapsed")
        with row[2]:
            units = (CUSTOM_UNITS if weighable
                     else [u for u in CUSTOM_UNITS if u not in bf.MOLAR])
            unit = st.selectbox("Unit", units, key=f"cus_unit_{compound}",
                                label_visibility="collapsed")
        with row[3]:
            form = st.selectbox(
                "Form", forms, key=f"cus_form_{compound}",
                label_visibility="collapsed",
                format_func=lambda f, c=compound: (
                    f"{f} ({compounds[c][f]})" if compounds[c][f] else f))
        if amount > 0:
            components.append({"compound": compound, "amount": amount,
                               "unit": unit, "form": form})

    if not components:
        st.info("Give at least one compound a concentration above zero.")
        return

    recipe = bf.custom_recipe(
        components, name="Composed here", ph=ph or None,
        ph_adjust_with=reagent, sterilization=sterilization)
    _results(book, recipe, volume, 1.0,
             {c["compound"]: c["form"] for c in components},
             f"custom_buffer_{volume:g}mL.csv", "cus_csv")


# --------------------------------------------------------------------------
def _results(book: dict, recipe: dict, volume: float, factor: float,
             forms: dict[str, str], filename: str, key: str) -> None:
    """The table, the warnings and the steps, for a recipe of either origin."""
    try:
        plan = bf.prepare(recipe, book["compounds"], volume, factor, forms)
    except ValueError as exc:
        st.error(str(exc))
        return

    columns = []
    for row in plan["rows"]:
        entry = {
            "Component": row["compound"],
            "Form": row["form"],
            "Formula": row["formula"] or "",
            "Molar mass (g/mol)": (f"{row['molar_mass']:.2f}"
                                   if row["molar_mass"] else ""),
        }
        if factor == 1:
            entry["Concentration"] = f"{row['amount_1x']:g} {row['unit']}"
        else:
            entry["1X"] = f"{row['amount_1x']:g} {row['unit']}"
            entry[f"At {factor:g}X"] = f"{row['amount_here']:g} {row['unit']}"
        entry["Amount"] = (bf.format_amount(row["amount"], row["amount_unit"])
                           or "not possible")
        entry["How"] = row["source"]
        columns.append(entry)
    table = pd.DataFrame(columns)
    st.dataframe(table, hide_index=True, **FIT)

    for row in plan["rows"]:
        if row["error"]:
            st.error(f"**{row['compound']}**: {row['error']}.")
        if row["note"]:
            st.caption(f"**{row['compound']}**: {row['note']}")

    st.caption(
        "Molar masses are derived from the formulas with the "
        f"{bf.ATOMIC_WEIGHTS_SOURCE}; they are not measured, and they are "
        "only as right as the formula on the line. The amounts are "
        "arithmetic. Nothing here knows solubilities: a high concentration "
        "can ask for more than will dissolve.")

    st.subheader("Steps")
    st.markdown("\n".join(f"{i}. {step}"
                          for i, step in enumerate(plan["steps"], 1)))

    st.download_button("Download the table (CSV)",
                       table.to_csv(index=False).encode(), mime="text/csv",
                       file_name=filename, key=key)
