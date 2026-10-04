"""The site index, on the parts of it that need no database.

A receipt is the record of what one run found in one dataset, and these are
the assertions that it can be read back: the events in list order with their
positions, a file-valued parameter distinguished from an unset one, and the
grouping of events into the run's sites -- including the two cases that a
naive grouping gets wrong, an event with no site number and a site the run
listed that no receipt reached.
"""
from pathlib import Path
from types import SimpleNamespace

from ....lib.utils.jobs import pandda_site_index as index


def write_receipt(directory: Path, dtag: str, run_uuid: str, events, *,
                  tree="/runs/pandda2_out", apo=True, name="params.xml") -> object:
    """A receipt job whose parameter file says what its arguments say."""
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    for event in events:
        parts = [f"<EVENT_IDX>{event['idx']}</EVENT_IDX>"]
        if event.get("site") is not None:
            parts.append(f"<SITE_IDX>{event['site']}</SITE_IDX>")
        if event.get("score") is not None:
            parts.append(f"<SCORE>{event['score']}</SCORE>")
        if event.get("probability") is not None:
            parts.append(f"<HIT_PROBABILITY>{event['probability']}</HIT_PROBABILITY>")
        if event.get("centroid"):
            x, y, z = event["centroid"]
            parts.append(f"<CENTROID><x>{x}</x><y>{y}</y><z>{z}</z></CENTROID>")
        # An unset CDataFile is written as an empty element, so "the element
        # is there" must not be read as "there is a file".
        parts.append("<POSE>%s</POSE>" % (
            "<baseName>pose.pdb</baseName>" if event.get("pose") else ""))
        parts.append("<EVENT_MAP>%s</EVENT_MAP>" % (
            "<baseName>event.ccp4</baseName>" if event.get("map", True) else ""))
        rows.append("<CPanddaEvent>%s</CPanddaEvent>" % "".join(parts))

    (directory / name).write_text(
        "<ccp4i2><ccp4i2_body id='pandda_events'>"
        "<inputData>"
        f"<PANDDA_OUT_DIR>{tree}</PANDDA_OUT_DIR>"
        f"<DTAG>{dtag}</DTAG>"
        f"<RUN_JOB_UUID>{run_uuid}</RUN_JOB_UUID>"
        "</inputData>"
        "<outputData>"
        "<XYZIN_APO>%s</XYZIN_APO>" % ("<baseName>apo.pdb</baseName>" if apo else "")
        + "<EVENTS>%s</EVENTS>" % "".join(rows)
        + "</outputData></ccp4i2_body></ccp4i2>"
    )
    return SimpleNamespace(directory=str(directory))


def test_events_come_back_in_list_order_with_their_positions(tmp_path):
    job = write_receipt(tmp_path / "j1", "xtal-0000", "run-1", [
        {"idx": 3, "site": 1, "score": 0.5, "pose": True},
        {"idx": 7, "site": 2, "score": 0.9, "pose": True},
    ])
    record = index.read_receipt(job)
    assert record["dtag"] == "xtal-0000"
    assert record["tree"] == "/runs/pandda2_out"
    assert [e["position"] for e in record["events"]] == [0, 1]
    # The position is what a scene reference needs and the event number is
    # what a person reads; PanDDA's numbering need not start at zero, so the
    # two must not be conflated.
    assert [e["event_idx"] for e in record["events"]] == [3, 7]


def test_an_empty_file_element_is_not_a_file(tmp_path):
    job = write_receipt(tmp_path / "j1", "xtal-0000", "run-1", [
        {"idx": 1, "site": 1, "pose": False, "map": True},
    ])
    event = index.read_receipt(job)["events"][0]
    assert event["has_map"] is True
    assert event["has_pose"] is False


def test_centroid_needs_all_three_axes(tmp_path):
    job = write_receipt(tmp_path / "j1", "xtal-0000", "run-1", [
        {"idx": 1, "site": 1, "centroid": (1.5, -2.0, 3.25)},
        {"idx": 2, "site": 1},
    ])
    events = index.read_receipt(job)["events"]
    assert events[0]["centroid"] == [1.5, -2.0, 3.25]
    assert events[1]["centroid"] is None


def test_a_receipt_with_no_parameter_file_is_skipped(tmp_path):
    (tmp_path / "empty").mkdir()
    assert index.read_receipt(SimpleNamespace(directory=str(tmp_path / "empty"))) is None


def test_input_params_alone_is_readable(tmp_path):
    """A receipt created but never run has no outputs, and is still a receipt:
    the panel must be able to say it is there and empty."""
    job = write_receipt(tmp_path / "j1", "xtal-0000", "run-1", [],
                        name="input_params.xml")
    record = index.read_receipt(job)
    assert record["dtag"] == "xtal-0000" and record["events"] == []


# -- grouping ---------------------------------------------------------------

def receipt(dtag, events, project="P", number=1):
    return {
        "dtag": dtag,
        "events": events,
        "project": {"id": 1, "uuid": "u", "name": project},
        "receipt": {"uuid": "r", "number": number, "id": number, "status": 6},
    }


def event(idx, site, score=None, pose=True, probability=None):
    return {"position": idx - 1, "event_idx": idx, "site_idx": site, "score": score,
            "hit_probability": probability, "has_pose": pose, "centroid": None}


def test_members_group_by_site_best_first():
    grouped = index.group_events_by_site([
        receipt("x0001", [event(1, 2, score=0.4), event(2, 5, score=0.9)], project="A"),
        receipt("x0002", [event(1, 2, score=0.8)], project="B", number=2),
    ], centroids={2: [1.0, 2.0, 3.0], 5: [4.0, 5.0, 6.0]})

    sites = {s["site_idx"]: s for s in grouped["sites"]}
    assert sorted(sites) == [2, 5]
    assert sites[2]["centroid"] == [1.0, 2.0, 3.0]
    assert sites[2]["n_events"] == 2 and sites[2]["n_datasets"] == 2
    assert sites[2]["best_score"] == 0.8
    # Best first: the reason to open a site is its strongest event.
    assert [m["dtag"] for m in sites[2]["members"]] == ["x0002", "x0001"]


def test_an_event_with_no_site_is_not_lost():
    """A run whose sites table was never written puts every event here, and a
    caller that cannot see them would show an empty campaign."""
    grouped = index.group_events_by_site(
        [receipt("x0001", [event(1, None, score=0.7)])], centroids={})
    assert grouped["sites"] == []
    assert [m["event_idx"] for m in grouped["unsited"]] == [1]


def test_a_site_no_receipt_reached_is_still_a_site():
    grouped = index.group_events_by_site(
        [receipt("x0001", [event(1, 1)])], centroids={1: None, 4: [0.0, 0.0, 0.0]})
    sites = {s["site_idx"]: s for s in grouped["sites"]}
    assert sorted(sites) == [1, 4]
    assert sites[4]["n_events"] == 0 and sites[4]["members"] == []
    assert sites[4]["best_score"] is None


def test_poses_are_counted_separately_from_events():
    """An event with no autobuilt pose is drawable in no panoptic scene, so a
    site's pose count is what says whether opening it is worth anything."""
    grouped = index.group_events_by_site([
        receipt("x0001", [event(1, 1, pose=True), event(2, 1, pose=False)]),
    ], centroids={})
    site = grouped["sites"][0]
    assert site["n_events"] == 2 and site["n_poses"] == 1


def test_a_site_centroid_comes_from_its_events_not_the_sites_table():
    """Design note section 9: the ``pandda_analyse_sites.csv`` centroid is
    frequently (0, 0, 0), so it is kept for comparison and not navigated by."""
    with_centroids = [dict(event(1, 1), centroid=[10.0, 0.0, 0.0]),
                      dict(event(2, 1), centroid=[12.0, 2.0, 0.0])]
    grouped = index.group_events_by_site(
        [receipt("x0001", with_centroids)], centroids={1: [0.0, 0.0, 0.0]})
    site = grouped["sites"][0]
    assert site["centroid"] == [11.0, 1.0, 0.0]
    assert site["table_centroid"] == [0.0, 0.0, 0.0]


def test_an_all_zero_table_centroid_is_no_position_at_all():
    grouped = index.group_events_by_site(
        [receipt("x0001", [event(1, 1)])], centroids={1: [0.0, 0.0, 0.0]})
    assert grouped["sites"][0]["centroid"] is None


def test_a_dispersed_site_gets_no_centroid():
    """PanDDA's clustering can sweep scattered weak events into one site. The
    mean of that cloud is a confident-looking coordinate pointing at nothing,
    so it is withheld and the spread is reported in its place."""
    scattered = [dict(event(i, 1), centroid=[i * 20.0, 0.0, 0.0]) for i in range(1, 6)]
    grouped = index.group_events_by_site(
        [receipt("x0001", scattered)], centroids={})
    site = grouped["sites"][0]
    assert site["dispersed"] is True
    assert site["centroid"] is None
    assert site["spread"]["rms"] > index.SITE_DISPERSION_LIMIT


def test_a_real_pocket_keeps_its_centroid_and_reports_a_small_spread():
    tight = [dict(event(i, 1), centroid=[10.0 + i * 0.5, 2.0, 3.0]) for i in range(1, 5)]
    site = index.group_events_by_site(
        [receipt("x0001", tight)], centroids={})["sites"][0]
    assert site["dispersed"] is False
    assert site["centroid"] is not None
    assert site["spread"]["rms"] < 2.0


def test_the_rollup_says_what_a_site_is_made_of():
    """A single 1.0 on top of a pile of noise read as a discovery. The median
    and the count of convincing events are what distinguish the two."""
    members = ([dict(event(1, 1), score=1.0)]
               + [dict(event(i, 1), score=0.2) for i in range(2, 12)])
    site = index.group_events_by_site(
        [receipt("x0001", members)], centroids={})["sites"][0]
    assert site["best_score"] == 1.0
    assert site["median_score"] == 0.2
    assert site["n_convincing"] == 1
