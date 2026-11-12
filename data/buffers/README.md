# Buffer and medium recipes

The data behind the **Buffers & media** section. Read-only in the app, like the
sequence library: edit the files here, check them, commit and push. They are in
a public repository: nothing confidential belongs in them.

| File | One row per | Holds |
|---|---|---|
| `compounds.csv` | form of a compound | name, form, formula |
| `recipes.csv` | recipe | pH and whether it is adjusted, sterilization, storage, usual stock factor, source |
| `recipe_components.csv` | component of a recipe | 1X concentration, unit, optional form and stock |

Each file starts with comment lines (`#`) describing its columns.

## Filling them in

1. **Compounds.** Keep the forms your lab stocks, add what is missing. Copy
   each formula from the label of the bottle, hydrate water after a dot:
   `Na2HPO4.7H2O`. No molar mass is typed anywhere: the app derives it from
   the formula with the CIAAW abridged standard atomic weights.
2. **One row per recipe** in `recipes.csv`, describing the **1X** buffer.
   A stock you make at 20X is the same recipe with `usual_stock` set to 20.
   `ph_mode` says whether the pH is adjusted to that value (`adjust`) or is
   what the 1X reads without adjustment (`measured`). `source` is required.
3. **One row per component** in `recipe_components.csv`, at its 1X
   concentration. Prefer molar units when the protocol gives them: the
   concentration then holds whatever hydrate is weighed, and `form` only says
   which one the app offers first. When the protocol gives a mass (g/L) of a
   compound with several forms, `form` is required: it says which form the
   mass refers to.
4. **Check:** `python3 tools/check_buffers.py` from the repository root. It
   runs the same loader as the app, so what it accepts is what the app shows.

## Spreadsheets

Excel and Numbers in a French or Belgian locale save CSV with `;` between
fields and `,` as decimal mark, which these files do not accept. Use a text
editor, or LibreOffice with "Field delimiter: ," on save. The check recognises
a `;` file and says so.

## What the check cannot see

It checks structure: formulas parse into elements with a standard atomic
weight, names match across files, units fit the compound, stocks are more
concentrated than what they are diluted to. It cannot know that a formula is
the one on your bottle or that a concentration is your protocol's. A typo that
yields valid symbols passes: `NaCI`, with a capital I, reads as sodium, carbon
and iodine.
