"""
ranking - comparing the sequences of one submission with each other, without
inventing a score for them.

`triage.py` puts every flag the app can raise in one table and refuses, in
writing, to collapse it into a single number: weighting a deamidation site
against a hydrophobic patch means choosing what matters, and that depends on
what the VHH is for. That refusal stands, and this module does not undo it.
What it adds is the part the table is bad at: seeing, at a glance, which
sequence of a set differs from the others and where.

Three layers, deliberately separate, none of them summed into another:

* **Ranks.** Each axis is ranked within the submitted set, on its own. A cell
  says "best of these eight on this axis", never "good". Nothing is added
  across axes, because adding a percentile of a polyreactivity model fitted on
  69 702 sequences to a count of motifs weighted by this laboratory produces a
  number with no unit and no meaning.

* **Gates.** A flag is a go or no-go, not points off. An unpaired cysteine does
  not make a sequence *slightly* worse, and a composite score would let three
  good axes outvote it. A gated sequence is marked, and stays in the table with
  its ranks intact, because whether a given flag disqualifies a candidate is
  the user's call and depends on what the protein is for.

* **Pareto front.** The only honest answer to "which is best" without weights:
  a sequence is on the front when no other sequence in the set is at least as
  good on every chosen axis and better on one. That names the few candidates
  worth arguing about without ever saying how much an axis is worth. A set
  where everything is on the front is a real answer too, it says the axes
  disagree and the choice is yours.

What is deliberately absent, and why:

* **No weights, no total, no stars, no grade.** Not even behind a disclaimer.
  A number gets screenshotted and the disclaimer does not travel with it.

* **No axis without a defensible direction.** `AXES` holds only the columns
  where one end is better than the other for a stated reason. There is no good
  pI, no good molecular weight and no good CDR3 length, so those stay in the
  triage table where they are read as descriptions rather than as scores. The
  hallmark score says how camelid-like a framework is, which is an identity,
  not a quality.

* **No claim about clinical suitability.** Humanness here is identity to the
  closest human germline, which is not a prediction of immunogenicity in a
  patient; polyreactivity is a pre-screen ranking from one synthetic library.
  Both are useful for putting eight clones in order and neither says a
  sequence is fit for use in people.

Every number comes from the triage table, which comes from the same functions
that fill the individual tabs, so nothing here can disagree with a tab.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import triage as tg

__all__ = ["Axis", "AXES", "BETTER_COLOUR", "DEFAULT_AXES", "DEFAULT_GATES",
           "GATE_LABELS",
           "axis_labels", "within_set_ranks", "gate_status", "pareto_front",
           "comparison_frame", "plot_rank_heatmap", "plot_pareto",
           "MIN_FOR_RANKS", "__version__"]

__version__ = "2026.11.17a"

# Below this many sequences a rank says nothing a side-by-side reading of two
# numbers does not say better: with two sequences every axis has a winner and a
# loser by construction, and the heatmap would show a wall of best and worst
# whatever the actual difference between them. The figure still draws, the
# interface says what the ranks are worth.
MIN_FOR_RANKS = 3


@dataclass(frozen=True)
class Axis:
    """One comparable axis: a triage column plus the direction that is better.

    ``better`` is ``"high"`` or ``"low"``, and ``why`` is the reason that
    direction is defensible, which is shown in the interface. ``basis`` says
    what kind of number it is, in the app's three words for a result:
    ``derived`` for something computed from the sequence by a published or
    in-house rule, ``theoretical`` for an index with no measurement behind it.
    No axis here is ``measured``, the third word, and none can be: this
    section reads sequences, it never sees the protein. Whether the rule
    behind a derived number was fitted to data, and to what, belongs in
    ``why``, where it can be stated rather than compressed into a word.
    """

    column: str
    label: str
    better: str
    basis: str
    why: str

    @property
    def ascending(self) -> bool:
        """Whether sorting ascending puts the better sequence first."""
        return self.better == "low"


# Ordered as they are drawn. Keyed by the triage column, which is checked
# against `triage.COLUMN_GROUPS` by a test: an axis naming a column the triage
# table does not produce is a typo that would otherwise show up as an empty
# column in the figure.
AXES: dict[str, Axis] = {
    a.column: a for a in (
        Axis("liability_score", "Liability score", "low", "derived",
             "The summed in-house weighting of the chemical liabilities found. "
             "Fewer and milder motifs is less to go wrong in storage."),
        Axis("n_high_risk", "High-risk liabilities", "low", "derived",
             "Motifs scoring at or above 6 on that same scale, which is where "
             "a single site is worth looking at on its own."),
        Axis("n_liab_in_CDR", "Liabilities in a CDR", "low", "derived",
             "A liability in a CDR risks the binding site itself, so the same "
             "motif matters more there than in a framework."),
        Axis("worst_window", "Worst aggregation window", "high", "derived",
             "The lowest solubility score along the sequence. Higher is "
             "further from the aggregation-prone end of the scale."),
        Axis("n_apr", "Aggregation-prone regions", "low", "derived",
             "How many windows fall below the threshold."),
        Axis("n_unusual", "Unusual framework residues", "low", "derived",
             "Framework positions rarer than the threshold in the reference "
             "repertoire. Common is not the same as correct, but a position "
             "almost no natural VHH carries is worth a reason."),
        Axis("n_orphan_cys", "Orphan cysteines", "low", "derived",
             "Thiols with no partner in the canonical pair. Also available as "
             "a gate, which is usually the better way to read it."),
        Axis("instability_index", "Instability index", "low", "theoretical",
             "Guruprasad's index, from dipeptide composition. An index with "
             "no measurement behind it, and its 40 cutoff was drawn on a "
             "small set of proteins; useful between candidates, weak as an "
             "absolute."),
        Axis("gravy", "GRAVY", "low", "theoretical",
             "Mean Kyte-Doolittle hydropathy over the whole chain. Lower is "
             "more hydrophilic, which usually helps solubility, but a VHH "
             "needs hydrophobic residues in its core and a paratope is often "
             "hydrophobic: this is the weakest direction claimed here."),
        Axis("humanness", "Humanness", "high", "derived",
             "Identity to the closest human germline V gene. Relevant when "
             "the plan is to humanise, and not a prediction of immunogenicity "
             "in a patient."),
        Axis("polyreactivity_pct", "Polyreactivity percentile", "low", "derived",
             "Percentile against the low-polyreactivity pool of Harvey et al. "
             "A pre-screen ranking from one synthetic library and one "
             "framework, not a measurement of this clone."),
    )
}

# What a first look shows: the axes that bear on whether a protein survives
# production and storage, which is the question this view is usually opened
# for. Everything else in AXES is one click away. A default is a choice, so it
# is made here, once, where it can be read.
DEFAULT_AXES = ("liability_score", "n_high_risk", "n_liab_in_CDR",
                "worst_window", "n_apr", "polyreactivity_pct")

# The flags offered as gates come from `triage.FLAG_COLUMNS` rather than from a
# list here, so a flag added to the triage table appears as a gate without a
# second edit. Only the labels are local, and a flag with no label is still
# offered, under its column name.
GATE_LABELS = {
    "has_high_risk": "A high-risk liability",
    "has_orphan_cys": "An orphan cysteine",
    "has_unusual_residue": "An unusual framework residue",
    "has_apr": "An aggregation-prone region",
    "liab_in_CDR": "A liability in a CDR",
    "polyreactive": "Flagged polyreactive",
}

# Two gates on by default, the two closest to a go or no-go when the protein is
# meant to be produced and kept: a thiol with no partner, which scrambles, and
# a liability the in-house scale puts at high risk. Neither is a verdict, both
# can be switched off, and the interface says so.
DEFAULT_GATES = ("has_orphan_cys", "has_high_risk")


def axis_labels(columns=None) -> dict[str, str]:
    """``{column: label}`` for the axes asked for, in the order of `AXES`."""
    wanted = set(AXES if columns is None else columns)
    return {c: a.label for c, a in AXES.items() if c in wanted}


def _usable(table: pd.DataFrame, columns) -> list[str]:
    """The asked-for axes that exist in this table and carry any number."""
    return [c for c in columns
            if c in AXES and c in table.columns
            and pd.to_numeric(table[c], errors="coerce").notna().any()]


def within_set_ranks(table: pd.DataFrame, columns=DEFAULT_AXES) -> pd.DataFrame:
    """Percentile rank of each sequence on each axis, within this set only.

    100 is the best of the submitted sequences on that axis, with ties sharing
    the mean of their positions. Indexed by sequence name, one column per axis.

    The worst is not 0 and cannot be. `rank(pct=True)` divides positions by the
    number of sequences, so the bottom cell is 100/n: 12.5 for eight sequences,
    33.3 for three, and 50.0 when two are compared. An earlier version of this
    docstring and of MANUAL.md said "0 the worst", which would have a reader
    expect a floor that never appears on screen.

    **The ranks are relative to what you submitted and to nothing else.** Add a
    ninth sequence and all eight others move. A sequence ranked 100 is the best
    of these, which may still be bad: the raw value is what says whether it is,
    and it is drawn in the cell next to the colour for exactly that reason.

    An empty cell stays empty. A sequence that is not a variable domain has no
    humanness, and one whose CDRs could not be placed has no polyreactivity
    percentile; ranking it against those that do would invent an order.
    """
    columns = _usable(table, columns)
    out = pd.DataFrame(index=pd.Index(table["name"], name="name"),
                       columns=columns, dtype=float)
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce")
        values.index = out.index
        # rank(pct=True) gives 0 < r <= 1 ascending; flip it for an axis whose
        # low end is the good one, so 100 is always the better sequence.
        ranked = values.rank(pct=True, ascending=not AXES[column].ascending)
        out[column] = (ranked * 100.0).round(1)
    return out


def gate_status(table: pd.DataFrame, gates=DEFAULT_GATES) -> pd.DataFrame:
    """Which of the chosen flags each sequence raises, and whether any does.

    Columns are the chosen gates plus ``blocked``, which is true when the
    sequence raises at least one of them, and ``reasons``, the labels of those
    it raises. A gate is not a score: nothing here is counted, weighted or
    subtracted from a rank.

    A flag reading False can mean "not flagged" or "not assessed": a sequence
    whose CDRs could not be placed is not flagged polyreactive because it has
    no percentile, which is not the same as being clean. The tab that owns each
    flag is where the reason sits, and `triage.py` says the same thing about
    its own empty cells.
    """
    gates = [g for g in gates if g in tg.FLAG_COLUMNS and g in table.columns]
    out = pd.DataFrame(index=pd.Index(table["name"], name="name"))
    for gate in gates:
        values = table[gate].fillna(False).astype(bool)
        values.index = out.index
        out[gate] = values
    out["blocked"] = out[gates].any(axis=1) if gates else False
    out["reasons"] = ["; ".join(GATE_LABELS.get(g, g) for g in gates
                                if bool(row[g]))
                      for _name, row in out.iterrows()] if gates else ""
    return out


def pareto_front(table: pd.DataFrame, columns=DEFAULT_AXES,
                 exclude=()) -> pd.DataFrame:
    """Which sequences no other sequence beats on every chosen axis at once.

    A sequence is **dominated** when some other sequence is at least as good as
    it on every axis and strictly better on at least one; the front is what is
    left. No axis is weighted against another, which is the whole point: the
    front is the set of defensible choices, and picking one of them is where
    judgement starts.

    Returns one row per sequence with ``comparable``, ``front`` and
    ``dominated_by``, the number of sequences that beat it outright.

    **A sequence missing a value on any chosen axis is not comparable** and is
    reported as such rather than dropped silently or treated as worst: it could
    dominate or be dominated and there is no way to tell. Drop the axis or the
    sequence. Names in ``exclude``, typically the gated ones, are left out of
    the comparison entirely, so a gated sequence neither sits on the front nor
    keeps another sequence off it.
    """
    columns = _usable(table, columns)
    names = list(table["name"])
    out = pd.DataFrame({"name": names})
    out["excluded"] = [n in set(exclude) for n in names]

    if not columns:
        out["comparable"] = False
        out["front"] = False
        out["dominated_by"] = 0
        return out

    values = table[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    # One sign per axis, so domination is a single "greater is better" test.
    signs = np.array([-1.0 if AXES[c].ascending else 1.0 for c in columns])
    values = values * signs

    complete = ~np.isnan(values).any(axis=1)
    considered = complete & ~out["excluded"].to_numpy()

    dominated_by = np.zeros(len(names), dtype=int)
    for i in range(len(names)):
        if not considered[i]:
            continue
        for j in range(len(names)):
            if i == j or not considered[j]:
                continue
            if np.all(values[j] >= values[i]) and np.any(values[j] > values[i]):
                dominated_by[i] += 1

    out["comparable"] = complete
    out["dominated_by"] = dominated_by
    out["front"] = considered & (dominated_by == 0)
    return out


def comparison_frame(table: pd.DataFrame, columns=DEFAULT_AXES,
                     gates=DEFAULT_GATES,
                     exclude_gated: bool = False) -> dict:
    """The three layers for one set of sequences, computed once.

    Returns ``{"ranks", "values", "gates", "pareto", "axes"}``: the ranks, the
    raw values they were ranked from, the gate table, the front, and the axes
    actually used, which can be fewer than asked for when a column is empty for
    every sequence.

    This is the single entry point the app uses, so the heatmap and the front
    are always drawn from one set of numbers. Nothing is cached here; the cache
    belongs in `compute.py`, which is also where the app's copy of this call
    lives.
    """
    columns = _usable(table, columns)
    gate_table = gate_status(table, gates)
    blocked = [n for n, row in gate_table.iterrows() if bool(row["blocked"])]
    values = table.set_index(pd.Index(table["name"], name="name"))[columns] \
                  .apply(pd.to_numeric, errors="coerce")
    return {
        "ranks": within_set_ranks(table, columns),
        "values": values,
        "gates": gate_table,
        "pareto": pareto_front(table, columns,
                               exclude=blocked if exclude_gated else ()),
        "axes": columns,
    }


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #

# Same reasoning as the identity heatmap in `vhh_core`: a cell that has shrunk
# past readability carries colour only, and the numbers are in the CSV.
HEATMAP_FIG_MAX = 14.0         # inches, longest side
HEATMAP_LABEL_MAX = 24         # sequences beyond which cells carry no number

# The direction arrows of the front figure: how far outside the axes they sit,
# in axes fractions, and the one colour this figure uses for nothing else, so
# an arrow cannot be read as a data point.
ARROW_OFFSET = -0.085
BETTER_COLOUR = "#55a868"
# The axis labels are pushed past the arrows rather than the arrows past the
# labels: an arrow set far enough out to clear any label would end up detached
# from the axis it describes, and on a figure with wide tick labels it would
# drift further still.
ARROW_LABEL_PAD = 26


def plot_rank_heatmap(ranks: pd.DataFrame, values: pd.DataFrame | None = None,
                      gates: pd.DataFrame | None = None,
                      title: str = "Within-set comparison"):
    """One row per sequence, one column per axis, coloured by rank.

    The colour is the rank **within this submission**, from worst to best on
    that axis alone, and the number written in the cell is the raw value, which
    is what says whether the best of the set is any good. Colour is never
    compared across columns: two cells of the same shade on two axes are both
    "best of these", of two different things.

    A gated sequence keeps its row and its ranks, and is marked in the margin.
    It is not greyed out or moved to the bottom: it may well be the sequence
    you take, with the flag on your list of things to watch.

    Needs no colourbar label beyond rank, and deliberately carries no overall
    column: there is no row total here because there is no defensible way to
    compute one.
    """
    import matplotlib.pyplot as plt

    names = list(ranks.index)
    columns = list(ranks.columns)
    n, m = len(names), len(columns)
    if not n or not m:
        fig, ax = plt.subplots(figsize=(5, 1.6))
        ax.text(0.5, 0.5, "Nothing to compare on the chosen axes",
                ha="center", va="center", fontsize=9)
        ax.axis("off")
        return fig

    width, height = 1.45 * m + 3.4, 0.52 * n + 2.0
    scale = min(1.0, HEATMAP_FIG_MAX / max(width, height))
    fig, ax = plt.subplots(figsize=(width * scale, height * scale))
    data = ranks.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(data)
    # cividis reads the same to a colour-blind eye and in greyscale, which the
    # red-to-green a ranking invites does not.
    cmap = plt.get_cmap("cividis").with_extremes(bad="#e8e8e8")
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=100, aspect="auto")

    blocked = set()
    if gates is not None and "blocked" in gates:
        blocked = {str(i) for i, row in gates.iterrows() if bool(row["blocked"])}
    # The marker is U+26A0, which the font matplotlib falls back to actually
    # carries. A nicer glyph that renders as an empty box says nothing.
    labels = [f"⚠ {name}" if str(name) in blocked else str(name)
              for name in names]
    ax.set_yticks(range(n), labels, fontsize=8 if n <= HEATMAP_LABEL_MAX else 6)
    ax.set_xticks(range(m),
                  [AXES[c].label if c in AXES else c for c in columns],
                  rotation=35, ha="right", fontsize=8)

    if n <= HEATMAP_LABEL_MAX and values is not None:
        for i, name in enumerate(names):
            for j, column in enumerate(columns):
                rank = data[i, j]
                raw = values.loc[name, column] if column in values else np.nan
                if np.isnan(rank):
                    text, colour = "n/a", "#555555"
                else:
                    text = "n/a" if pd.isna(raw) else (
                        f"{raw:g}" if abs(float(raw)) < 1000 else f"{raw:.3g}")
                    colour = "white" if rank < 55 else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=7.5,
                        color=colour)

    ax.set_title(title, fontsize=11)
    bar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    bar.set_label("Rank within this set (100 = best of these)", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


def plot_pareto(values: pd.DataFrame, x: str, y: str,
                front: pd.DataFrame | None = None,
                gates: pd.DataFrame | None = None):
    """Two axes against each other, with the non-dominated sequences marked.

    Drawn for two axes at a time because that is the number a scatter can show
    honestly; the front itself is computed on every axis chosen, in
    `pareto_front`, and a sequence can sit on that front while looking
    unremarkable here. The figure says so rather than letting the picture pass
    for the computation.

    **Both axes run in their natural direction**, values increasing to the
    right and upwards, and the end that is better is marked by an arrow along
    each axis with the corner named above the figure. An earlier version
    inverted any axis whose good end is low, so that the front always sat in
    the top right. That is a real convention in multi-objective plots and it
    was the wrong trade here: it spent the strongest expectation a reader
    brings to a chart, that numbers grow rightwards and upwards, to save them
    reading two axis labels. With a mix of directions the labels have to be
    read anyway, so the inversion did not remove the confusion, it moved it
    from the corners to the tick marks, where it is met once per point instead
    of once per figure.

    The cost of the choice, stated because it is real: the preferable corner
    moves with the axes picked. Two liability axes put it at the bottom left,
    humanness against polyreactivity at the bottom right. The arrows and the
    line above the figure are derived from the axes in hand rather than
    written down, so they cannot disagree with what is drawn.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    xs = pd.to_numeric(values[x], errors="coerce")
    ys = pd.to_numeric(values[y], errors="coerce")
    on_front, blocked = set(), set()
    if front is not None:
        on_front = set(front.loc[front["front"], "name"])
    if gates is not None and "blocked" in gates:
        blocked = {str(i) for i, row in gates.iterrows() if bool(row["blocked"])}

    for name in values.index:
        xv, yv = xs.loc[name], ys.loc[name]
        if pd.isna(xv) or pd.isna(yv):
            continue
        in_front = str(name) in on_front
        ax.scatter(xv, yv, s=120 if in_front else 70,
                   facecolor="#4c72b0" if in_front else "none",
                   edgecolor="#c44e52" if str(name) in blocked else "#4c72b0",
                   linewidths=1.8 if str(name) in blocked else 1.2, zorder=3)
        ax.annotate(str(name), (xv, yv), textcoords="offset points",
                    xytext=(7, 4), fontsize=7.5)

    if on_front:
        pts = sorted(((xs.loc[n], ys.loc[n]) for n in values.index
                      if str(n) in on_front
                      and not pd.isna(xs.loc[n]) and not pd.isna(ys.loc[n])),
                     key=lambda p: p[0])
        if len(pts) > 1:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#4c72b0",
                    lw=1.0, ls="--", alpha=0.7, zorder=2)

    axis_x, axis_y = AXES.get(x), AXES.get(y)
    ax.set_xlabel(axis_x.label if axis_x else x, fontsize=9,
                  labelpad=ARROW_LABEL_PAD)
    ax.set_ylabel(axis_y.label if axis_y else y, fontsize=9,
                  labelpad=ARROW_LABEL_PAD)
    ax.grid(alpha=0.25)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # The better direction as an arrow outside the axes, in the one colour this
    # figure uses for nothing else, so it cannot be read as a data point. Drawn
    # in axes fractions, which is why it survives any data range.
    for axis, near, far in ((axis_x, (0.0, ARROW_OFFSET), (0.32, ARROW_OFFSET)),
                            (axis_y, (ARROW_OFFSET, 0.0), (ARROW_OFFSET, 0.32))):
        if axis is None:
            continue
        start, end = (far, near) if axis.ascending else (near, far)
        ax.annotate("", xy=end, xytext=start, xycoords="axes fraction",
                    textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=BETTER_COLOUR,
                                    lw=1.4), annotation_clip=False)
    corner = (f"{'bottom' if axis_y is None or axis_y.ascending else 'top'} "
              f"{'left' if axis_x is None or axis_x.ascending else 'right'}")
    ax.text(0.5, 1.02, f"Better is towards the {corner}, where the green "
                       "arrows point.", transform=ax.transAxes, ha="center",
            fontsize=8, color=BETTER_COLOUR)
    ax.set_title("Filled: on the front over every chosen axis. "
                 "Red outline: gated.", fontsize=9, pad=22)
    fig.tight_layout()
    return fig
