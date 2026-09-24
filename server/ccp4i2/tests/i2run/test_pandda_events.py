"""The receipt task end to end (design note sections 7.2 and 15.2): a short
receipt returns UNSATISFACTORY and still gleans what arrived; every nested
file inside the EVENTS list gets a File row; a zero-event dataset is a clean
receipt with no events. The synthetic tree needs no volume; the last test
reads a real tree when it is mounted."""
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from ccp4i2.db import models
from ccp4i2.tests.unit.pandda.synthetic_tree import event_record, make_tree
from .utils import i2run

REAL_TREE = Path("/Volumes/LocalStore/pandda/seeded_ab_sub/A1_stock")


def _job_row(job_dir):
    job_id = ET.parse(job_dir / "params.xml").find(".//jobId").text
    return models.Job.objects.get(uuid=job_id)


def _param_names(db_job):
    return sorted(f.job_param_name for f in models.File.objects.filter(job=db_job))


def test_complete_receipt_gleans_every_nested_file(tmp_path):
    tree = make_tree(tmp_path / "pandda2_out", {
        "xtal-0004": [event_record(1, bdc=0.80), event_record(2, bdc=0.95, build=False),
                      {"pandda_model": True}],
    })
    with i2run(["pandda_events", "--PANDDA_OUT_DIR", str(tree), "--DTAG", "xtal-0004"]) as job:
        db_job = _job_row(job)
        assert db_job.status == models.Job.Status.FINISHED
        names = _param_names(db_job)
        assert names == sorted([
            "XYZIN_APO", "ZMAP", "PANDDA_MODEL",
            "EVENTS[0].EVENT_MAP", "EVENTS[0].POSE", "EVENTS[1].EVENT_MAP",
        ]), names
        for name in ("XYZIN_APO.pdb", "ZMAP.map", "event_1_map.map", "event_1_pose.pdb", "event_2_map.map"):
            assert (job / name).is_file(), name
        assert not (job / "XYZIN_APO.pdb").is_symlink()
        # PanDDA wrote the Z-map as P1; the copy carries the crystal's space
        # group from the apo model, so coot treats it as the X-ray map it is.
        # The boxed event map is left as written; the tree is untouched.
        import gemmi
        assert gemmi.read_ccp4_map(str(job / "ZMAP.map")).grid.spacegroup.number == 20
        assert gemmi.read_ccp4_map(str(job / "event_1_map.map")).grid.spacegroup.number == 1
        assert gemmi.read_ccp4_map(str(tree / "processed_datasets" / "xtal-0004" / "xtal-0004-z_map.native.ccp4")).grid.spacegroup.number == 1
        # PanDDA wrote LIG; the copies carry the true component and the tree is untouched
        assert " MZ0 " in (job / "event_1_pose.pdb").read_text()
        assert " MZ0 " in (job / "PANDDA_MODEL.pdb").read_text()
        assert " LIG " in (tree / "processed_datasets" / "xtal-0004" / "xtal-0004_event_1_best_autobuild.pdb").read_text()
        params = ET.parse(job / "params.xml")
        assert params.find(".//EVENTS/CPanddaEvent/LIGAND_ID").text == "MZ0"
        kpis = {v.key.name: v.value for v in models.JobFloatValue.objects.select_related("key").filter(job=db_job)}
        assert kpis["nEventsExpected"] == 2 and kpis["nEventsDelivered"] == 2
        assert kpis["nPosesExpected"] == 1 and kpis["nPosesDelivered"] == 1
        assert kpis["bestBuildScore"] == pytest.approx(0.57)
        xml = ET.parse(job / "program.xml")
        assert xml.findtext("counts/events_expected") == "2"
        assert not xml.findall("shortfalls/item")


def test_short_receipt_is_unsatisfactory_and_still_gleans(tmp_path):
    tree = make_tree(tmp_path / "pandda2_out", {
        "xtal-0007": [event_record(1), (event_record(2), {"event_map": False}),
                      (event_record(3), {"pose": False})],
    })
    with i2run(["pandda_events", "--PANDDA_OUT_DIR", str(tree), "--DTAG", "xtal-0007",
                "--RUN_INCOMPLETE", "True"], allow_errors=True) as job:
        db_job = _job_row(job)
        assert db_job.status == models.Job.Status.UNSATISFACTORY
        names = _param_names(db_job)
        # what arrived is published: event 1 whole, event 2's pose, event 3's map
        assert names == sorted([
            "XYZIN_APO", "ZMAP",
            "EVENTS[0].EVENT_MAP", "EVENTS[0].POSE", "EVENTS[1].POSE", "EVENTS[2].EVENT_MAP",
        ]), names
        kpis = {v.key.name: v.value for v in models.JobFloatValue.objects.select_related("key").filter(job=db_job)}
        assert (kpis["nEventsExpected"], kpis["nEventsDelivered"]) == (3, 2)
        assert (kpis["nPosesExpected"], kpis["nPosesDelivered"]) == (3, 2)
        xml = ET.parse(job / "program.xml")
        assert [n.text for n in xml.findall("shortfalls/item")] == [
            "event 2: event map", "event 3: candidate pose"]
        assert xml.findtext("run_incomplete") == "True"
        diagnostic = ET.parse(job / "diagnostic.xml")
        codes = {r.findtext("code") for r in diagnostic.findall(".//errorReport")}
        assert {"202", "204"} <= codes


def test_zero_event_dataset_is_a_clean_empty_receipt(tmp_path):
    tree = make_tree(tmp_path / "pandda2_out", {"xtal-0000": []})
    with i2run(["pandda_events", "--PANDDA_OUT_DIR", str(tree), "--DTAG", "xtal-0000"]) as job:
        db_job = _job_row(job)
        assert db_job.status == models.Job.Status.FINISHED
        assert _param_names(db_job) == ["XYZIN_APO", "ZMAP"]
        kpis = {v.key.name: v.value for v in models.JobFloatValue.objects.select_related("key").filter(job=db_job)}
        assert kpis["nEventsExpected"] == 0


def test_missing_dataset_fails(tmp_path):
    tree = make_tree(tmp_path / "pandda2_out", {"xtal-0000": []})
    with i2run(["pandda_events", "--PANDDA_OUT_DIR", str(tree), "--DTAG", "xtal-4242"],
               allow_errors=True) as job:
        assert _job_row(job).status == models.Job.Status.FAILED


@pytest.mark.skipif(not REAL_TREE.is_dir(), reason="real PanDDA tree not mounted")
def test_real_tree_dataset_with_one_event():
    with i2run(["pandda_events", "--PANDDA_OUT_DIR", str(REAL_TREE), "--DTAG", "xtal-0004"]) as job:
        db_job = _job_row(job)
        names = _param_names(db_job)
        assert "EVENTS[0].EVENT_MAP" in names and "EVENTS[0].POSE" in names
        assert "XYZIN_APO" in names and "ZMAP" in names
        assert not (job / "XYZIN_APO.pdb").is_symlink()
        xml = ET.parse(job / "program.xml")
        event = xml.find("events/event")
        assert float(event.findtext("bdc")) == pytest.approx(0.802, abs=1e-3)
        assert float(event.findtext("optimal_contour")) == pytest.approx(0.975, abs=1e-3)
