"""
sections.buffers - what to weigh for a buffer or a medium, at any volume and
concentration factor, with the form of each salt you actually have.

The screen draws. The recipes are read by `compute.load_buffers_cached`, and
the arithmetic, the formulas and the atomic weights live in `buffers.py`.
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

import buffers as bf
from compute import load_buffers_cached
from shell import FIT, SECTIONS

__version__ = "2026.11.12a"


def render() -> None:
    """Draw the Buffers & media screen (key "buffers")."""
    icon, label, *_ = SECTIONS["buffers"]
    st.title(f"{icon} {label}")
    st.markdown(
        "Pick a recipe, say how much you need and at what concentration, and "
        "say which form of each salt is on your shelf: the mass to weigh "
        "follows from the formula of that form. A 20X stock and the 1X it "
        "gives are one recipe.")

    book, problems = load_buffers_cached()

    st.info(
        "**Read-only here, and public.** The recipes are the files in "
        "`data/buffers/` of the repository, edited there and pushed, like the "
        "sequence library. They are one laboratory's recipes, each with the "
        "source it gives, not a reference collection.", icon="🧂")
    for message in problems:
        st.warning(message, icon="⚠️")

    recipes = book["recipes"]
    if not recipes:
        st.info("No recipe could be loaded. Add rows to the files in "
                "`data/buffers/` and check them with "
                "`python3 tools/check_buffers.py`.")
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

    try:
        plan = bf.prepare(recipe, book["compounds"], volume, factor, forms)
    except ValueError as exc:
        st.error(str(exc))
        return

    table = pd.DataFrame([{
        "Component": row["compound"],
        "Form": row["form"],
        "Formula": row["formula"] or "",
        "Molar mass (g/mol)": (f"{row['molar_mass']:.2f}"
                               if row["molar_mass"] else ""),
        "1X": f"{row['amount_1x']:g} {row['unit']}",
        f"At {factor:g}X": f"{row['amount_here']:g} {row['unit']}",
        "Amount": (bf.format_amount(row["amount"], row["amount_unit"])
                   or "not possible"),
        "How": row["source"],
    } for row in plan["rows"]])
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
        "factor can ask for more than will dissolve.")

    st.subheader("Steps")
    st.markdown("\n".join(f"{i}. {step}"
                          for i, step in enumerate(plan["steps"], 1)))

    csv = table.to_csv(index=False).encode()
    st.download_button(
        "Download the table (CSV)", csv, mime="text/csv",
        file_name=f"{name.replace(' ', '_')}_{factor:g}X_{volume:g}mL.csv",
        key="buf_csv")
