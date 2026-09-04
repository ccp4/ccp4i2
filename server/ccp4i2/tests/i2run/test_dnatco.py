"""End-to-end tests for the DNATCO wrapper and pipeline.

Gated on DNATCO being launchable the way a job launches it: through the
CCP4 launcher script, or ``node`` driving the bundled dnatco.js.
"""
import json
from pathlib import Path

import gemmi
from pytest import fixture, mark

from ccp4i2.config.program_discovery import discover_program
from ccp4i2.wrappers.dnatco.script.dnatco import find_dnatco_command
from .urls import rcsb_mmcif
from .utils import download, i2run


requires_dnatco = mark.skipif(
    discover_program(find_dnatco_command())["path"] is None,
    reason="DNATCO not installed (no dnatco launcher, and no $CCP4/dnatco/bin/dnatco.js for node)",
)

NTC_OVERALL_ITEMS = [
    "_ndb_struct_ntc_overall.confal_score",
    "_ndb_struct_ntc_overall.confal_percentile",
    "_ndb_struct_ntc_overall.num_classified",
    "_ndb_struct_ntc_overall.num_unclassified",
    "_ndb_struct_ntc_overall.num_unclassified_rmsd_close",
]


@fixture(scope="module")
def cif1q93():
    with download(rcsb_mmcif("1q93")) as path:
        yield path


@fixture(scope="module")
def cif1q96():
    with download(rcsb_mmcif("1q96")) as path:
        yield path


def check_extended_cif(path: Path):
    assert path.is_file(), f"Extended mmCIF not found: {path}"
    gemmi.read_structure(str(path), format=gemmi.CoorFormat.Mmcif)
    block = gemmi.cif.read_file(str(path))[0]
    for item in NTC_OVERALL_ITEMS:
        assert block.find_value(item), f"Missing {item} in {path}"


def check_naval_json(path: Path):
    assert path.is_file(), f"NAVAL JSON not found: {path}"
    with open(path) as handle:
        data = json.load(handle)
    for key in ("navalLengthsStats", "navalAnglesStats", "lengths", "angles"):
        assert key in data, f"NAVAL JSON lacks {key}"
    assert data["lengths"] and data["lengths"][0]["details"]


@requires_dnatco
def test_dnatco_wrapper(cif1q93):
    args = ["dnatco"]
    args += ["--XYZIN", cif1q93]
    args += ["--GENERATE_RESTRAINTS", "True"]
    with i2run(args) as job:
        # The wrapper keeps DNATCO's own file names (--prefix dnatco).
        check_extended_cif(job / "dnatco_extended.cif")
        check_naval_json(job / "dnatco_angles_lengths_by_residue.json")
        restraints = job / "dnatco_restraints_refmac.txt"
        assert restraints.is_file()
        assert "external" in restraints.read_text()[:2000]
        summary = (job / "program.xml").read_text()
        assert "<confal_score>" in summary


@requires_dnatco
def test_dnatco_pipe_one_model(cif1q93):
    args = ["dnatco_pipe"]
    args += ["--XYZIN1", cif1q93]
    args += ["--GENERATE_RESTRAINTS", "True"]
    with i2run(args) as job:
        check_extended_cif(job / "CIFOUT1.cif")
        check_naval_json(job / "JSONOUT1.json")
        assert not (job / "CIFOUT2.cif").exists()
        assert (job / "RESTRAINTS.txt").is_file()


@requires_dnatco
def test_dnatco_pipe_two_models(cif1q93, cif1q96):
    args = ["dnatco_pipe"]
    args += ["--XYZIN1", cif1q93]
    args += ["--TOGGLE_XYZIN2", "True"]
    args += ["--XYZIN2", cif1q96]
    args += ["--GENERATE_RESTRAINTS", "True"]
    args += ["--MAX_RMSD", "0.49"]
    args += ["--RESTRAINTS_SIGMA", "0.9"]
    with i2run(args) as job:
        for index in (1, 2):
            check_extended_cif(job / f"CIFOUT{index}.cif")
            check_naval_json(job / f"JSONOUT{index}.json")
        assert (job / "RESTRAINTS.txt").is_file()
        # Restraints come from the second model, and the parameters reach DNATCO.
        command = (job / "job_2" / "com.txt").read_text() if (job / "job_2" / "com.txt").exists() else ""
        if command:
            assert "--restraintsRmsd 0.49" in command
            assert "--restraintsSigmaFactor 0.9" in command
        summary = (job / "program.xml").read_text()
        assert summary.count("<DNATCO ") == 2
