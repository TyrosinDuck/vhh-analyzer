"""
plates - serial dilutions and plate maps, the arithmetic of a plate.

Two things the rest of the app did not do. `bench.dilution_plan` works out one
dilution; a series of them is a different calculation, and the one most often
got wrong. And nothing in the app knew what a well was, although the assays
these sequences end up in are run in plates.

Both are arithmetic and combinatorics on numbers you supply. Nothing here
models an assay: a layout says what goes where, not what the readings will be,
and a dilution series says what to pipette, not what the curve will look like.

The carry-over is the whole reason the series is worth writing down. Pipetting
V_t from one position into the next, which holds V_d of diluent, dilutes by

    f = (V_t + V_d) / V_t

so a position that must still hold V_final *after* giving V_t away has to be
made up with V_final + V_t, and

    V_t * f = V_final + V_t   =>   V_t = V_final / (f - 1)

from which the diluent per position falls out as V_t * (f - 1) = V_final: the
same volume of diluent in every position, which is what makes the series easy
to pipette once the transfer is known. Forgetting the carry-over is invisible,
and it either leaves the last positions short or quietly changes the last
dilution.
"""

from __future__ import annotations

import pandas as pd

# The minimum pipettable volume comes from `bench`, where the conjugation
# planner already refuses anything below it. One definition of what a hand can
# pipette, not one per screen.
from bench import PIPETTE_MINIMUM_UL

__all__ = [
    "PLATE_FORMATS", "ROW_LETTERS", "SERIES", "INDEPENDENT", "DILUTION_MODES",
    "serial_dilution", "plate_layout", "well_name", "__version__",
]

__version__ = "2026.11.17a"


# Rows and columns per format, and nothing else: a plate is a grid, and the
# volumes it holds depend on the plate and the protocol rather than on the
# count of its wells, so no working volume is tabulated here.
PLATE_FORMATS: dict[int, tuple[int, int]] = {96: (8, 12), 384: (16, 24)}

ROW_LETTERS = "ABCDEFGHIJKLMNOP"

# The two ways a dilution series is actually made, named because the choice
# changes the numbers and belongs to whoever is at the bench.
SERIES = "series"
INDEPENDENT = "independent"
DILUTION_MODES = (SERIES, INDEPENDENT)


def well_name(row: int, column: int) -> str:
    """``A1`` for (0, 0). Rows are letters down the plate, columns numbers."""
    return f"{ROW_LETTERS[row]}{column + 1}"


def serial_dilution(start: float, fold: float, points: int,
                    volume_ul: float, mode: str = SERIES,
                    unit: str = "mg/mL",
                    pipette_minimum_ul: float = PIPETTE_MINIMUM_UL) -> dict:
    """What to pipette for a dilution series, by either of the two protocols.

    ``mode`` is `SERIES`, where each position is made from the one before it,
    or `INDEPENDENT`, where each is made from the stock. They are not two
    estimates of one number, they are two protocols, and the app does not
    choose between them: a series carries the pipetting error of every step
    before it into every step after, and independent dilutions do not, at the
    cost of more stock and of transfers that become unpipettable as the
    dilution deepens.

    ``volume_ul`` is the volume each position must still hold when the series
    is finished. In `SERIES` mode the positions are therefore made up with
    ``volume_ul + transfer_ul``, and the last one keeps the transfer it never
    passes on, which is listed as ``surplus_ul`` rather than silently left in
    the figures.

    Returns plain data: ``rows``, one per position, plus ``transfer_ul``,
    ``diluent_per_position_ul``, ``stock_needed_ul``, ``diluent_total_ul``,
    ``surplus_ul`` and ``warnings``.
    """
    if mode not in DILUTION_MODES:
        raise ValueError(f"mode must be one of {DILUTION_MODES}, not {mode!r}")
    if start <= 0:
        raise ValueError("the starting concentration must be positive")
    if fold <= 1:
        raise ValueError(
            f"the fold per step must be above 1; {fold:g} would not dilute")
    points = int(points)
    if points < 2:
        raise ValueError("a series needs at least two positions")
    if volume_ul <= 0:
        raise ValueError("the volume per position must be positive")

    warnings: list[str] = []
    rows: list[dict] = []

    if mode == SERIES:
        transfer_ul = volume_ul / (fold - 1.0)
        diluent_each = volume_ul
        for index in range(points):
            dilution = fold ** index
            rows.append({
                "position": index + 1,
                "dilution": round(dilution, 6),
                "concentration": start / dilution,
                "unit": unit,
                # The first position is the undiluted stock, so it is made up
                # from the stock itself rather than from a transfer.
                "stock_ul": round(volume_ul + transfer_ul, 3) if index == 0 else 0.0,
                "diluent_ul": 0.0 if index == 0 else round(diluent_each, 3),
                "from_previous_ul": 0.0 if index == 0 else round(transfer_ul, 3),
                "pass_on_ul": 0.0 if index == points - 1 else round(transfer_ul, 3),
                "final_volume_ul": round(volume_ul, 3),
                "note": "",
            })
        rows[-1]["note"] = (
            "keeps the transfer it never passes on; discard it or accept the "
            "extra volume")
        stock_needed = volume_ul + transfer_ul
        diluent_total = diluent_each * (points - 1)
        surplus = transfer_ul

        if transfer_ul < pipette_minimum_ul:
            # Said as the volume that would fix it, because that is the knob:
            # the transfer follows from the final volume and the fold.
            needed = pipette_minimum_ul * (fold - 1.0)
            warnings.append(
                f"The transfer is {transfer_ul:.2f} uL, below the "
                f"{pipette_minimum_ul:.0f} uL this app treats as pipettable. "
                f"At {fold:g}-fold steps a final volume of at least "
                f"{needed:.0f} uL per position would be needed, or make the "
                f"dilutions independently instead.")
    else:
        transfer_ul = 0.0
        diluent_each = 0.0
        unpipettable = []
        stock_needed = 0.0
        diluent_total = 0.0
        for index in range(points):
            dilution = fold ** index
            stock_ul = volume_ul / dilution
            diluent_ul = volume_ul - stock_ul
            note = ""
            if stock_ul < pipette_minimum_ul:
                note = "stock volume below the pipettable minimum"
                unpipettable.append(index + 1)
            rows.append({
                "position": index + 1,
                "dilution": round(dilution, 6),
                "concentration": start / dilution,
                "unit": unit,
                "stock_ul": round(stock_ul, 3),
                "diluent_ul": round(diluent_ul, 3),
                "from_previous_ul": 0.0,
                "pass_on_ul": 0.0,
                "final_volume_ul": round(volume_ul, 3),
                "note": note,
            })
            stock_needed += stock_ul
            diluent_total += diluent_ul
        surplus = 0.0
        if unpipettable:
            deepest = volume_ul / (fold ** (points - 1))
            warnings.append(
                "Positions " + ", ".join(str(p) for p in unpipettable) +
                f" need less stock than the {pipette_minimum_ul:.0f} uL this "
                f"app treats as pipettable (the last one, {deepest:.3f} uL). "
                f"Make an intermediate dilution first, or run the series by "
                f"transfer, where every step pipettes the same volume.")

    return {
        "mode": mode,
        "rows": rows,
        "transfer_ul": round(transfer_ul, 3),
        "diluent_per_position_ul": round(diluent_each, 3),
        "stock_needed_ul": round(stock_needed, 3),
        "diluent_total_ul": round(diluent_total, 3),
        "surplus_ul": round(surplus, 3),
        "final_dilution": round(fold ** (points - 1), 6),
        "final_concentration": start / (fold ** (points - 1)),
        "warnings": warnings,
    }


def _blocks(samples, controls, dilution_points: int, replicates: int):
    """The wells to place, grouped so that one group is one column of a plate.

    A group is one replicate of one sample across its dilutions, which is the
    order an assay is normally laid out in: the series runs down a column and
    the replicates sit side by side. Dilution-major order would interleave the
    replicates with the dilutions and put neither on an axis of the plate.
    """
    groups = []
    for name in samples:
        for replicate in range(1, replicates + 1):
            wells = []
            for point in range(1, dilution_points + 1):
                label = str(name)
                if dilution_points > 1:
                    label += f" D{point}"
                if replicates > 1:
                    label += f" R{replicate}"
                wells.append({"content": label, "sample": str(name),
                              "dilution": point if dilution_points > 1 else None,
                              "replicate": replicate, "kind": "sample"})
            groups.append(wells)
    for name in controls:
        wells = []
        for replicate in range(1, replicates + 1):
            label = str(name)
            if replicates > 1:
                label += f" R{replicate}"
            wells.append({"content": label, "sample": str(name),
                          "dilution": None, "replicate": replicate,
                          "kind": "control"})
        groups.append(wells)
    return groups


def plate_layout(samples, dilution_points: int = 1, replicates: int = 1,
                 controls=(), plate: int = 96, skip_edges: bool = False,
                 by_column: bool = True,
                 group_per_line: bool = True) -> dict:
    """A plate map: which well holds what, and how much of the plate is left.

    The placement is deliberately predictable rather than optimal. Wells are
    filled in order, down the columns by default, and with ``group_per_line``
    each replicate of each sample starts at the top of a fresh column, which
    is what makes the printed map match the way the series is pipetted. No
    attempt is made to pack the plate: a layout whose positions move when an
    unrelated sample is added is a layout nobody can check against the bench.

    ``skip_edges`` leaves the outer ring of wells out of the usable positions.
    Many protocols do that, or fill the ring with buffer, because the outer
    wells evaporate differently; the app neither quantifies that effect nor
    takes a view on whether it matters for a given assay. It only places
    wells where it is told it may.

    Returns ``grid`` (a DataFrame of labels, rows by columns), ``wells`` (one
    record per filled well, which is the form a reader or a database ingests),
    and the counts: ``capacity``, ``usable``, ``used``, ``free``, plus
    ``warnings``.
    """
    if plate not in PLATE_FORMATS:
        raise ValueError(
            f"plate must be one of {sorted(PLATE_FORMATS)}, not {plate!r}")
    dilution_points, replicates = int(dilution_points), int(replicates)
    if dilution_points < 1:
        raise ValueError("a layout needs at least one dilution point")
    if replicates < 1:
        raise ValueError("a layout needs at least one replicate")
    samples = [s for s in (samples or []) if str(s).strip()]
    controls = [c for c in (controls or []) if str(c).strip()]
    if not samples and not controls:
        raise ValueError("give at least one sample or one control to place")

    n_rows, n_cols = PLATE_FORMATS[plate]
    rows_available = range(1, n_rows - 1) if skip_edges else range(n_rows)
    cols_available = range(1, n_cols - 1) if skip_edges else range(n_cols)
    rows_available, cols_available = list(rows_available), list(cols_available)

    # Fill order, and the lines the groups are aligned to.
    if by_column:
        lines = [[(r, c) for r in rows_available] for c in cols_available]
    else:
        lines = [[(r, c) for c in cols_available] for r in rows_available]

    groups = _blocks(samples, controls, dilution_points, replicates)
    warnings: list[str] = []
    placed: list[dict] = []
    line_index, offset = 0, 0

    line_length = len(lines[0]) if lines else 0
    out_of_room = False

    for group in groups:
        if group_per_line:
            if len(group) > line_length:
                raise ValueError(
                    f"a group of {len(group)} wells does not fit in one "
                    f"{'column' if by_column else 'row'} of {line_length} "
                    f"usable wells; turn off one group per "
                    f"{'column' if by_column else 'row'}, or use fewer "
                    f"dilution points")
            # Start the group at the top of a fresh line, unless the current
            # one is still untouched.
            if offset:
                line_index, offset = line_index + 1, 0
        for well in group:
            while line_index < len(lines) and offset >= len(lines[line_index]):
                line_index, offset = line_index + 1, 0
            if line_index >= len(lines):
                out_of_room = True
                break
            row, column = lines[line_index][offset]
            placed.append({"well": well_name(row, column),
                           "row": ROW_LETTERS[row], "column": column + 1,
                           **well})
            offset += 1
        if out_of_room:
            break

    requested = sum(len(g) for g in groups)
    usable = sum(len(line) for line in lines)
    if len(placed) < requested:
        short = requested - len(placed)
        plates_needed = -(-requested // usable) if usable else 0
        warnings.append(
            f"{requested} wells were asked for and {usable} are usable on a "
            f"{plate}-well plate with these settings, so {short} did not fit. "
            f"Split the run over {plates_needed} plates, or reduce the "
            f"replicates or the dilution points.")

    grid = pd.DataFrame("", index=list(ROW_LETTERS[:n_rows]),
                        columns=[str(c + 1) for c in range(n_cols)])
    for record in placed:
        grid.loc[record["row"], str(record["column"])] = record["content"]

    return {
        "grid": grid,
        "wells": placed,
        "capacity": n_rows * n_cols,
        "usable": usable,
        "used": len(placed),
        "free": usable - len(placed),
        "requested": requested,
        "warnings": warnings,
    }
