"""The campaign's dataset x site matrix, on the rules that need no database.

Three rules decide what the overview shows in a dataset's cell for a site:
which PanDDA event (if any) is "at" the site, whether the dataset's frame can
be compared with the parent's at all, and which job holds the model a click
opens. Each is a pure function, pinned here on plain values.
"""
from pathlib import Path
from types import SimpleNamespace

import pytest

from ....lib import campaign_matrix as matrix
from ....lib.campaign_matrix import (
    DIMPLE_TASKS,
    REFINEMENT_TASKS,
    choose_current_model,
    frame_mismatch,
    site_reach,
    match_events_to_sites,
    read_model_cell,
)

CELL = [60.0, 60.0, 60.0, 90.0, 90.0, 90.0]


def site(uuid, origin, radius=8.0):
    return {"uuid": uuid, "origin": list(origin), "radius": radius}


def event(idx, centroid, probability=0.5, pose=True):
    return {"event_idx": idx, "centroid": list(centroid) if centroid else None,
            "hit_probability": probability, "has_pose": pose}


# --------------------------------------------------------------------------
# Matching events to sites
# --------------------------------------------------------------------------

def test_every_site_gets_a_cell_even_with_no_events():
    cells = match_events_to_sites([], [site("a", (0, 0, 0)), site("b", (9, 9, 9))])
    assert cells == {"a": None, "b": None}


def test_the_nearest_event_within_the_radius_wins():
    cells = match_events_to_sites(
        [event(1, (5, 0, 0), 0.9), event(2, (2, 0, 0), 0.1), event(3, (30, 0, 0))],
        [site("a", (0, 0, 0))],
    )
    assert cells["a"] == {"event_idx": 2, "hit_probability": 0.1,
                          "distance": 2.0, "has_pose": True}


def test_the_radius_boundary_is_inclusive_and_just_beyond_is_out():
    on_edge = match_events_to_sites([event(1, (8.0, 0, 0))], [site("a", (0, 0, 0))])
    assert on_edge["a"]["distance"] == 8.0
    beyond = match_events_to_sites([event(1, (8.001, 0, 0))], [site("a", (0, 0, 0))])
    assert beyond["a"] is None


def test_radius_is_per_site():
    events = [event(1, (5, 0, 0))]
    cells = match_events_to_sites(
        events, [site("tight", (0, 0, 0), radius=4.0), site("loose", (0, 0, 0), radius=6.0)])
    assert cells["tight"] is None
    assert cells["loose"]["event_idx"] == 1


def test_a_site_without_a_radius_uses_the_default():
    cells = match_events_to_sites([event(1, (7.9, 0, 0))],
                                  [{"uuid": "a", "origin": [0, 0, 0], "radius": None}])
    assert cells["a"]["event_idx"] == 1
    assert matrix.DEFAULT_SITE_RADIUS == 8.0


def test_one_event_can_count_for_two_overlapping_sites():
    cells = match_events_to_sites(
        [event(4, (3, 0, 0), pose=False)],
        [site("left", (0, 0, 0)), site("right", (6, 0, 0))],
    )
    assert cells["left"]["event_idx"] == 4 and cells["left"]["distance"] == 3.0
    assert cells["right"]["event_idx"] == 4 and cells["right"]["distance"] == 3.0
    assert cells["left"]["has_pose"] is False


def test_ties_break_on_probability_then_event_number_not_input_order():
    a = event(7, (0, 3, 0), probability=0.2)
    b = event(5, (3, 0, 0), probability=0.8)
    c = event(2, (0, 0, 3), probability=0.8)
    for order in ([a, b, c], [c, b, a], [b, a, c]):
        cells = match_events_to_sites(order, [site("s", (0, 0, 0))])
        assert cells["s"]["event_idx"] == 2


def test_an_event_without_a_centroid_never_matches():
    cells = match_events_to_sites([event(1, None)], [site("a", (0, 0, 0))])
    assert cells["a"] is None


def test_a_frame_mismatch_matches_nothing():
    cells = match_events_to_sites([event(1, (0, 0, 0))], [site("a", (0, 0, 0))],
                                  mismatch="cell a off")
    assert cells == {"a": None}


# --------------------------------------------------------------------------
# The frame check
# --------------------------------------------------------------------------

def test_cells_within_tolerance_are_the_same_frame():
    assert frame_mismatch([61.1, 59.0, 60.5, 90.0, 91.9, 90.0], CELL) is None


def test_an_edge_more_than_two_percent_off_is_flagged():
    reason = frame_mismatch([60.0, 61.3, 60.0, 90, 90, 90], CELL)
    assert reason is not None and "cell b" in reason and "2.2%" in reason


def test_an_angle_more_than_two_degrees_off_is_flagged():
    reason = frame_mismatch([60.0, 60.0, 60.0, 90, 92.5, 90], CELL)
    assert reason is not None and "beta" in reason


def test_an_unknown_cell_is_not_a_mismatch():
    assert frame_mismatch(None, CELL) is None
    assert frame_mismatch(CELL, None) is None
    assert frame_mismatch([60, 60], CELL) is None


# --------------------------------------------------------------------------
# The current model
# --------------------------------------------------------------------------

FINISHED = matrix.FINISHED


def job(id, task, status=FINISHED, parent_id=None):
    return SimpleNamespace(id=id, task_name=task, status=status, parent_id=parent_id)


def test_finished_matches_the_model_status():
    pytest.importorskip("django")
    from ....db import models
    assert matrix.FINISHED == models.Job.Status.FINISHED


def test_the_latest_top_level_refinement_wins():
    jobs = [job(1, "i2Dimple"), job(2, "refmac"), job(3, "servalcat_pipe"),
            job(4, "coot_rebuild")]
    assert choose_current_model(jobs).id == 3


def test_a_refinement_beats_a_newer_dimple_run():
    jobs = [job(2, "servalcat_pipe"), job(9, "i2Dimple")]
    assert choose_current_model(jobs).id == 2


def test_a_refinement_subjob_does_not_count_but_its_pipeline_does():
    jobs = [job(5, "SubstituteLigand"), job(6, "i2Dimple", parent_id=5),
            job(7, "servalcat_pipe", parent_id=5)]
    assert choose_current_model(jobs).id == 5


def test_unfinished_refinements_fall_back_to_dimple_at_any_level():
    jobs = [job(5, "SubstituteLigand", status=5), job(6, "i2Dimple", parent_id=5),
            job(8, "refmac", status=3)]
    assert choose_current_model(jobs).id == 6


def test_no_model_job_is_none():
    assert choose_current_model([job(1, "acedrg"), job(2, "refmac", status=5)]) is None
    assert choose_current_model([]) is None


def test_dicts_are_accepted_as_well_as_rows():
    jobs = [{"id": 1, "task_name": "dimple", "status": FINISHED, "parent_id": None}]
    assert choose_current_model(jobs)["id"] == 1


def test_the_task_sets_are_disjoint_and_name_real_tasks():
    assert not set(REFINEMENT_TASKS) & set(DIMPLE_TASKS)
    pytest.importorskip("django")
    from ....core.tasks import TASKS
    # i2Refmac and dimple are legacy names kept for imported jobs; the rest must
    # be a registered task, or a typo would silently never match.
    missing = (set(REFINEMENT_TASKS) | set(DIMPLE_TASKS)) - set(TASKS) - {"i2Refmac", "dimple"}
    assert not missing


# --------------------------------------------------------------------------
# Reading a cell
# --------------------------------------------------------------------------

def test_read_cell_from_a_pdb_header(tmp_path: Path):
    path = tmp_path / "m.pdb"
    path.write_text(
        "REMARK hello\n"
        "CRYST1   61.200   62.300   63.400  90.00 100.50  90.00 P 1 21 1      2\n"
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n"
    )
    assert read_model_cell(path) == [61.2, 62.3, 63.4, 90.0, 100.5, 90.0]


def test_read_cell_from_an_mmcif_header(tmp_path: Path):
    path = tmp_path / "m.cif"
    path.write_text(
        "data_x\n"
        "_cell.entry_id x\n"
        "_cell.length_a 61.2\n"
        "_cell.length_b 62.3(2)\n"
        "_cell.length_c   63.4\n"
        "_cell.angle_alpha 90\n"
        "_cell.angle_beta 100.5\n"
        "_cell.angle_gamma 90\n"
        "_symmetry.space_group_name_H-M 'P 1 21 1'\n"
    )
    assert read_model_cell(path) == [61.2, 62.3, 63.4, 90.0, 100.5, 90.0]


def test_no_cell_and_placeholder_cells_read_as_none(tmp_path: Path):
    none = tmp_path / "n.pdb"
    none.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000\n")
    assert read_model_cell(none) is None
    unit = tmp_path / "u.pdb"
    unit.write_text("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1\n")
    assert read_model_cell(unit) is None
    assert read_model_cell(tmp_path / "missing.pdb") is None


# The frame check measures how far a cell difference moves the sites. BAZ2B's
# 5e9l (a = 80.92 against the parent's 83.03, 2.5% off) was hidden by a flat
# 2% rule, though its pocket ~27 A from the origin moves by only ~0.7 A.
BAZ2B_PARENT = [83.03, 96.38, 57.81, 90.0, 90.0, 90.0]
BAZ2B_5E9L = [80.919, 96.38, 57.81, 90.0, 90.0, 90.0]


def test_a_small_cell_difference_near_the_origin_is_compared_directly():
    reach = site_reach([{"origin": (-25.7, -9.3, 0.3)}])
    assert 27 < reach < 28
    assert frame_mismatch(BAZ2B_5E9L, BAZ2B_PARENT, reach) is None


def test_the_same_difference_far_from_the_origin_is_flagged():
    reason = frame_mismatch(BAZ2B_5E9L, BAZ2B_PARENT, 80.0)
    assert reason is not None and "at the sites" in reason


def test_another_crystal_form_is_flagged_whatever_the_reach():
    other = [70.0, 96.38, 57.81, 90.0, 90.0, 90.0]  # 16% off
    assert frame_mismatch(other, BAZ2B_PARENT, 1.0) is not None
