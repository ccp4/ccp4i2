"""The orchestrator under i2run against a hand-specified dataset list: no
campaign, no database beyond the job itself (design note 15.2, v1 item 4).
If this needed a campaign, section 3.1 would have been violated.

Stage-only runs everywhere: it stages the contract's tree, writes the
manifest and provenance, and finishes without PanDDA. A real run needs the
volume, pandda2.analyse and at least 25 datasets, and is opt-in."""
import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from ccp4i2.db import models
from ccp4i2.tests.unit.pandda.conftest import LABELS, source_files
from .programs import requires_program
from .utils import i2run

REAL_TREE = Path("/Volumes/LocalStore/pandda/BAZ2B/datasets")


def _dataset_args(labels):
    args = []
    for label in labels:
        pdb, mtz, cif = source_files(label)
        args += ["--DATASETS", f"DTAG={label}", f"XYZIN={pdb}", f"HKLIN={mtz}", f"DICT={cif}"]
    return args


def test_stage_only_needs_no_campaign():
    assert models.ProjectGroup.objects.count() == 0, "the architectural test: no campaign"
    args = ["pandda_campaign", "--RUN_MODE", "stage_only", "--LOCAL_CPUS", "3"] + _dataset_args(LABELS)
    with i2run(args) as job:
        staging = job / "staging"
        for i, label in enumerate(LABELS):
            d = staging / "datasets" / f"xtal-{i:04d}"
            assert (d / "final.pdb").is_file() and (d / "final.mtz").is_file() and (d / "dict.cif").is_file(), label
        assert (staging / "Projects.csv").read_text().splitlines()[1] == f"xtal-0000, {LABELS[0]}"
        manifest = json.loads((staging / "manifest.json").read_text())
        assert [e["label"] for e in manifest["datasets"]] == LABELS
        assert manifest["provenance"]["run_mode"] == "stage_only"
        assert manifest["provenance"]["contract"].startswith("CCP4I2_PANDDA_INVOCATION_CONTRACT")
        assert manifest["provenance"]["sizing_hint"] == {"datasets": 3, "cell_volume_class": "small"}

        job_id = ET.parse(job / "params.xml").find(".//jobId").text
        db_job = models.Job.objects.get(uuid=job_id)
        assert db_job.status == models.Job.Status.FINISHED
        rows = models.File.objects.filter(job=db_job)
        assert rows.filter(job_param_name="MANIFEST").count() == 1
        # the nine inputs were imported into the project and registered
        # (imports are named by objectName, not path: XYZIN, not DATASETS[0].XYZIN)
        assert rows.filter(directory=models.File.Directory.IMPORT_DIR).count() == 9
        assert models.ProjectGroup.objects.count() == 0

        params = ET.parse(job / "params.xml")
        argv = params.find(".//outputData/PROVENANCE_ARGV").text
        assert "--dataset_range 0-999999999" in argv and "--ligand_pdb_regex ligand.pdb" in argv
        assert "--local_cpus 3" in argv
        assert params.find(".//outputData/CONTRACT_VERSION").text
        # the submitted list is the record
        items = params.findall(".//inputData/DATASETS/CPanddaDataset")
        assert [i.findtext("DTAG") for i in items] == LABELS
        assert all(i.find("XYZIN/baseName") is not None for i in items)
        diagnostic = ET.parse(job / "diagnostic.xml")
        codes = {r.findtext("code") for r in diagnostic.findall(".//errorReport")}
        assert {"205", "222"} <= codes, codes   # too few datasets to run; stage-only said so


def test_an_empty_list_is_refused_before_anything_runs():
    with i2run(["pandda_campaign", "--RUN_MODE", "stage_only"], allow_errors=True) as job:
        job_id = ET.parse(job / "params.xml").find(".//jobId").text
        assert models.Job.objects.get(uuid=job_id).status == models.Job.Status.FAILED
        assert not (job / "staging").exists()


@pytest.mark.skipif(not REAL_TREE.is_dir(), reason="BAZ2B datasets not mounted")
@pytest.mark.skipif(not os.environ.get("CCP4I2_PANDDA_E2E"), reason="set CCP4I2_PANDDA_E2E=1 to run PanDDA (an afternoon)")
@requires_program("pandda2.analyse")
def test_local_run_over_baz2b():
    n = int(os.environ.get("CCP4I2_PANDDA_E2E_DATASETS", "30"))
    args = ["pandda_campaign", "--LOCAL_CPUS", os.environ.get("CCP4I2_PANDDA_E2E_CPUS", "4")]
    for d in sorted(REAL_TREE.iterdir())[:n]:
        args += ["--DATASETS", f"DTAG={d.name}", f"XYZIN={d / 'final.pdb'}", f"HKLIN={d / 'final.mtz'}"]
        if (d / "dict.cif").exists():
            args += [f"DICT={d / 'dict.cif'}"]
    with i2run(args) as job:
        assert (job / "pandda2_out" / "analyses" / "pandda_analyse_events.csv").is_file()
        assert (job / "pandda2_out" / "processed_datasets").is_dir()
        job_id = ET.parse(job / "params.xml").find(".//jobId").text
        assert models.Job.objects.get(uuid=job_id).status == models.Job.Status.FINISHED
