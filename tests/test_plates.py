"""
Dilution series and plate maps.

The arithmetic is checked against hand calculations rather than against the
functions' own output, because a test that recomputes the implementation tests
nothing. The series by transfer is the one worth being careful about: the
carry-over is invisible when it is wrong, and the numbers below were worked
out on paper first.
"""

from __future__ import annotations

import pytest

import plates as pt


# --------------------------------------------------------------------------
# The series, by transfer
# --------------------------------------------------------------------------

def test_a_two_fold_series_transfers_the_whole_final_volume():
    """At 2-fold the transfer equals the final volume: V/(2-1) = V.

    The easiest case to check in the head, and the one that shows the
    carry-over plainly. Each position is made up with 100 uL of diluent,
    receives 100 uL from the one before, and passes 100 uL on, which leaves
    the 100 uL it was asked to hold.
    """
    out = pt.serial_dilution(10.0, 2, 8, 100.0)

    assert out["transfer_ul"] == pytest.approx(100.0)
    assert out["diluent_per_position_ul"] == pytest.approx(100.0)
    assert out["final_dilution"] == pytest.approx(128)        # 2 ** 7
    assert out["final_concentration"] == pytest.approx(10.0 / 128)
    assert [row["final_volume_ul"] for row in out["rows"]] == [100.0] * 8
    assert not out["warnings"]


def test_the_carry_over_is_what_the_first_position_is_made_up_with():
    """V_final + V_t, not V_final.

    The whole point. A first position holding only the final volume cannot
    pass any on without going short, and the error does not show until the
    last positions come out light.
    """
    out = pt.serial_dilution(1.0, 4, 5, 60.0)

    transfer = 60.0 / 3                                        # V/(f-1) = 20
    assert out["transfer_ul"] == pytest.approx(transfer)
    assert out["rows"][0]["stock_ul"] == pytest.approx(60.0 + transfer)
    assert out["stock_needed_ul"] == pytest.approx(80.0)
    # Every later position takes the same diluent and the same transfer, which
    # is what makes one pipette setting do the whole series.
    for row in out["rows"][1:]:
        assert row["diluent_ul"] == pytest.approx(60.0)
        assert row["from_previous_ul"] == pytest.approx(transfer)


def test_the_last_position_keeps_the_transfer_it_never_passes_on():
    """Stated rather than quietly folded into the figures."""
    out = pt.serial_dilution(1.0, 10, 4, 90.0)

    assert out["rows"][-1]["pass_on_ul"] == 0.0
    assert out["surplus_ul"] == pytest.approx(10.0)            # 90/(10-1)
    assert "discard" in out["rows"][-1]["note"]


def test_the_dilutions_are_powers_of_the_fold():
    out = pt.serial_dilution(1000.0, 3, 5, 100.0)

    assert [row["dilution"] for row in out["rows"]] == [1, 3, 9, 27, 81]
    assert [row["concentration"] for row in out["rows"]] == [
        pytest.approx(1000.0 / d) for d in (1, 3, 9, 27, 81)]


def test_an_unpipettable_transfer_says_which_volume_would_fix_it():
    """The transfer follows from the final volume and the fold, so the final
    volume is the knob, and naming it is more use than refusing."""
    out = pt.serial_dilution(10.0, 10, 4, 15.0)

    assert out["transfer_ul"] == pytest.approx(15.0 / 9, abs=1e-3)
    assert out["warnings"], "an unpipettable transfer passed unremarked"
    # 2 uL minimum at 10-fold steps needs 2 x (10-1) = 18 uL per position.
    assert "18" in out["warnings"][0]


# --------------------------------------------------------------------------
# The series, independently
# --------------------------------------------------------------------------

def test_independent_dilutions_take_the_stock_straight():
    """Each position is V/d of stock made up to V, and nothing is transferred."""
    out = pt.serial_dilution(10.0, 10, 4, 200.0, mode=pt.INDEPENDENT)

    assert [row["stock_ul"] for row in out["rows"]] == [
        pytest.approx(v) for v in (200.0, 20.0, 2.0, 0.2)]
    assert [row["diluent_ul"] for row in out["rows"]] == [
        pytest.approx(v) for v in (0.0, 180.0, 198.0, 199.8)]
    assert all(row["from_previous_ul"] == 0.0 for row in out["rows"])
    assert out["transfer_ul"] == 0.0
    assert out["stock_needed_ul"] == pytest.approx(222.2, abs=0.1)


def test_a_deep_independent_dilution_is_flagged_position_by_position():
    """Where the two protocols differ, and why the app does not choose.

    A series by transfer pipettes the same volume at every step however deep
    it goes; independent dilutions stop being pipettable.
    """
    out = pt.serial_dilution(10.0, 10, 5, 200.0, mode=pt.INDEPENDENT)

    flagged = [row["position"] for row in out["rows"] if row["note"]]
    assert flagged == [4, 5]
    assert out["warnings"] and "intermediate dilution" in out["warnings"][0]


@pytest.mark.parametrize("kwargs", [
    {"fold": 1.0},          # would not dilute
    {"fold": 0.5},          # would concentrate
    {"points": 1},          # not a series
    {"volume_ul": 0.0},
    {"start": 0.0},
])
def test_a_series_refuses_inputs_that_are_not_a_series(kwargs):
    args = {"start": 10.0, "fold": 2.0, "points": 4, "volume_ul": 100.0}
    args.update(kwargs)
    with pytest.raises(ValueError):
        pt.serial_dilution(**args)


def test_an_unknown_mode_is_refused_rather_than_defaulted():
    with pytest.raises(ValueError, match="mode"):
        pt.serial_dilution(10.0, 2, 4, 100.0, mode="whatever")


# --------------------------------------------------------------------------
# The plate map
# --------------------------------------------------------------------------

def test_a_series_runs_down_a_column_and_replicates_sit_side_by_side():
    """The layout an assay is normally set up in.

    Dilution-major order would interleave the replicates with the dilutions
    and put neither on an axis of the plate, which is the arrangement nobody
    can check against the bench.
    """
    out = pt.plate_layout(["NbA", "NbB"], dilution_points=8, replicates=2)
    grid = out["grid"]

    assert grid.loc["A", "1"] == "NbA D1 R1"
    assert grid.loc["H", "1"] == "NbA D8 R1"       # the series went down
    assert grid.loc["A", "2"] == "NbA D1 R2"       # the replicate beside it
    assert grid.loc["A", "3"] == "NbB D1 R1"
    assert out["used"] == 32 and out["free"] == 64


def test_controls_are_placed_after_the_samples_with_the_same_replicates():
    out = pt.plate_layout(["NbA"], dilution_points=4, replicates=2,
                          controls=["blank"])

    wells = {w["well"]: w["content"] for w in out["wells"]}
    assert wells["A3"] == "blank R1" and wells["B3"] == "blank R2"
    assert [w["kind"] for w in out["wells"]].count("control") == 2


def test_leaving_the_outer_wells_empty_shrinks_what_is_usable():
    """An option, not a recommendation: the app places wells where it is told
    it may, and quantifies nothing about the edges."""
    out = pt.plate_layout(["X"], dilution_points=6, skip_edges=True)

    assert out["usable"] == 60                      # 6 rows x 10 columns
    assert out["capacity"] == 96
    assert all(w["row"] not in ("A", "H") for w in out["wells"])
    assert all(w["column"] not in (1, 12) for w in out["wells"])


def test_a_layout_that_does_not_fit_says_how_many_plates_it_needs():
    out = pt.plate_layout([f"S{i}" for i in range(1, 11)],
                          dilution_points=8, replicates=2)

    assert out["requested"] == 160 and out["used"] == 96
    assert out["warnings"] and "2 plates" in out["warnings"][0]


def test_a_group_longer_than_a_column_is_refused_with_the_reason():
    """Rather than silently wrapping a series across two columns, which would
    read on the printed map as two different samples."""
    with pytest.raises(ValueError, match="does not fit in one column"):
        pt.plate_layout(["Z"], dilution_points=20)


def test_a_384_plate_is_sixteen_rows_by_twenty_four():
    out = pt.plate_layout(["X", "Y"], dilution_points=12, replicates=3,
                          plate=384)

    assert out["grid"].shape == (16, 24)
    assert out["capacity"] == 384
    assert out["used"] == 72
    assert out["grid"].loc["L", "1"] == "X D12 R1"


def test_along_a_row_puts_the_series_on_the_other_axis():
    out = pt.plate_layout(["X"], dilution_points=12, by_column=False)

    assert out["grid"].loc["A", "1"] == "X D1"
    assert out["grid"].loc["A", "12"] == "X D12"


def test_well_names_are_the_ones_printed_on_the_plate():
    assert pt.well_name(0, 0) == "A1"
    assert pt.well_name(7, 11) == "H12"
    assert pt.well_name(15, 23) == "P24"


@pytest.mark.parametrize("kwargs, fragment", [
    ({"plate": 48}, "plate must be"),
    ({"dilution_points": 0}, "dilution point"),
    ({"replicates": 0}, "replicate"),
    ({"samples": []}, "at least one sample"),
])
def test_a_layout_refuses_what_it_cannot_draw(kwargs, fragment):
    args = {"samples": ["X"], "dilution_points": 2, "replicates": 1}
    args.update(kwargs)
    with pytest.raises(ValueError, match=fragment):
        pt.plate_layout(**args)


# --------------------------------------------------------------------------
# Labelling a layout with a series, and the printed sheet
# --------------------------------------------------------------------------

def test_a_series_can_label_the_wells_with_its_concentrations():
    """What the series view hands over: labels, not measurements."""
    out = pt.plate_layout(["NbA"], dilution_points=3,
                          dilution_labels=["10 mg/mL", "5 mg/mL", "2.5 mg/mL"])

    assert [w["content"] for w in out["wells"]] == [
        "NbA 10 mg/mL", "NbA 5 mg/mL", "NbA 2.5 mg/mL"]


def test_labels_that_do_not_match_the_points_are_refused():
    """A map whose labels have slipped by one is worse than one that only
    numbers them."""
    with pytest.raises(ValueError, match="dilution labels"):
        pt.plate_layout(["NbA"], dilution_points=3,
                        dilution_labels=["a", "b"])


def test_the_printed_sheet_carries_the_wells_and_no_volumes(tmp_path):
    from pypdf import PdfReader

    layout = pt.plate_layout(["NbA", "NbB"], dilution_points=4, replicates=2,
                             controls=["blank"])
    path = pt.build_plate_pdf(tmp_path / "p.pdf", layout,
                              title="Screening plate", version="x")
    text = PdfReader(str(path)).pages[0].extract_text()

    assert "Screening plate" in text
    assert "NbA D1 R1" in text and "blank" in text
    assert "wells used" in text
    assert "carries no volumes" in text


def test_a_384_sheet_builds_too(tmp_path):
    layout = pt.plate_layout(["X"], dilution_points=12, replicates=4,
                             plate=384)
    path = pt.build_plate_pdf(tmp_path / "p384.pdf", layout)

    assert path.exists() and path.stat().st_size > 1000
