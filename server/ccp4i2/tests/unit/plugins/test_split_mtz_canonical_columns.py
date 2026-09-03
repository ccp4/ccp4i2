"""
Every mini-MTZ carries canonical column labels.

A mini-MTZ is identified downstream by its column labels against the class's
CONTENT_SIGNATURE_LIST, not by the contentFlag stamped on the output object.
splitMtz kept the source labels (``Fobs,Sigma``) for everything but Phs, and
produced a file no consumer could read: phaser_MR_AUTO's makeHklin asked
ObsDataConverter for FMEAN, setContentFlag() found no signature match, and the
job failed with "Cannot determine contentFlag from input file".

The one relabelling lives on ``CMiniMtzDataFile.canonical_column_mapping``;
the base class then checks every mini-MTZ a task wrote before it is gleaned.
"""

import pytest

from ccp4i2 import I2_TOP
from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.wrappers.import_common import canonical_mapping


# --- the mapping -------------------------------------------------------------

@pytest.mark.parametrize(
    "cls, content_flag, labels, expected",
    [
        (CCP4XtalData.CObsDataFile, 4, ["Fobs", "Sigma"], ["F", "SIGF"]),
        (CCP4XtalData.CObsDataFile, 3, ["IMEAN", "SIGIMEAN"], ["I", "SIGI"]),
        (CCP4XtalData.CObsDataFile, 1,
         ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"],
         ["Iplus", "SIGIplus", "Iminus", "SIGIminus"]),
        (CCP4XtalData.CFreeRDataFile, 1, ["FreeR_flag"], ["FREER"]),
        (CCP4XtalData.CFreeRDataFile, None, ["FreeR_flag"], ["FREER"]),
        (CCP4XtalData.CMapCoeffsDataFile, 1, ["FWT", "PHWT"], ["F", "PHI"]),
        (CCP4XtalData.CPhsDataFile, 1,
         ["HLA", "HLB", "HLC", "HLD"], ["HLA", "HLB", "HLC", "HLD"]),
        (CCP4XtalData.CPhsDataFile, 2, ["PHIB", "FOM"], ["PHI", "FOM"]),
    ],
)
def test_every_type_maps_to_its_signature(cls, content_flag, labels, expected):
    result = cls.canonical_column_mapping(content_flag, labels)
    assert list(result.keys()) == labels
    assert list(result.values()) == expected


def test_column_count_mismatch_keeps_original_names():
    labels = ["F", "SIGF", "EXTRA"]
    assert CCP4XtalData.CObsDataFile.canonical_column_mapping(4, labels) == {c: c for c in labels}


def test_unknown_content_flag_keeps_original_names():
    labels = ["A", "B"]
    cls = CCP4XtalData.CObsDataFile
    assert cls.canonical_column_mapping(0, labels) == {c: c for c in labels}
    assert cls.canonical_column_mapping(9, labels) == {c: c for c in labels}
    assert cls.canonical_column_mapping("junk", labels) == {c: c for c in labels}


def test_import_tasks_share_the_one_mapping():
    """import_common.canonical_mapping is the same answer under the old name."""
    cls = CCP4XtalData.CObsDataFile
    assert canonical_mapping(cls, 4, ["FP", "SIGFP"]) == cls.canonical_column_mapping(4, ["FP", "SIGFP"])


# --- the file that failed in the wild ---------------------------------------

def _split_beta_blip(tmp_path, mapping):
    gemmi = pytest.importorskip("gemmi")
    from ccp4i2.core.CCP4Utils import split_mtz_file

    src = I2_TOP / "demo_data" / "beta_blip" / "beta_blip_P3221.mtz"
    labels = [c.label for c in gemmi.read_mtz_file(str(src)).columns
              if c.label not in ("H", "K", "L")]
    assert labels == ["Fobs", "Sigma"], "fixture no longer exercises the bug"
    out = tmp_path / "obs.mtz"
    split_mtz_file(src, out, mapping(labels))
    return out


def test_split_beta_blip_is_readable_as_fmean(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    cls = CCP4XtalData.CObsDataFile
    out = _split_beta_blip(tmp_path, lambda labels: cls.canonical_column_mapping(4, labels))

    assert [c.label for c in gemmi.read_mtz_file(str(out)).columns] == \
        ["H", "K", "L", "F", "SIGF"]

    obs = cls(fullPath=str(out))
    obs.setContentFlag()
    assert int(obs.contentFlag) == cls.CONTENT_FLAG_FMEAN
    assert obs.column_contract_violation() is None


def test_source_labels_are_a_contract_violation(tmp_path):
    pytest.importorskip("gemmi")
    obs = CCP4XtalData.CObsDataFile(
        fullPath=str(_split_beta_blip(tmp_path, lambda labels: {c: c for c in labels})))
    obs.contentFlag.set(4)
    problem = obs.column_contract_violation()
    assert problem is not None
    assert "Fobs,Sigma" in problem and "F,SIGF" in problem


def test_unset_content_flag_accepts_any_signature(tmp_path):
    pytest.importorskip("gemmi")
    cls = CCP4XtalData.CObsDataFile
    obs = cls(fullPath=str(_split_beta_blip(tmp_path, lambda labels: cls.canonical_column_mapping(4, labels))))
    obs.contentFlag.unSet()
    assert not obs.contentFlag.isSet()
    assert obs.column_contract_violation() is None


# --- the guard in the base class --------------------------------------------

def test_a_task_that_writes_a_non_canonical_mini_mtz_fails(tmp_path):
    """The base class catches the producer, not the next consumer."""
    pytest.importorskip("gemmi")
    task = get_plugin_class("splitMtz")()
    bad = _split_beta_blip(tmp_path, lambda labels: {c: c for c in labels})
    out = task.container.outputData.MINIMTZOUTLIST
    out.append(CCP4XtalData.CObsDataFile(fullPath=str(bad)))
    out[-1].contentFlag.set(4)

    assert task.checkMiniMtzOutputs() == CPluginScript.FAILED
    (err,) = [e for e in task.errorReport.getErrors() if e["code"] == 990]
    assert "Fobs,Sigma" in err["details"]


def test_a_canonical_mini_mtz_passes_the_guard(tmp_path):
    pytest.importorskip("gemmi")
    cls = CCP4XtalData.CObsDataFile
    task = get_plugin_class("splitMtz")()
    good = _split_beta_blip(tmp_path, lambda labels: cls.canonical_column_mapping(4, labels))
    task.container.outputData.MINIMTZOUTLIST.append(cls(fullPath=str(good)))
    assert task.checkMiniMtzOutputs() == CPluginScript.SUCCEEDED
    assert not [e for e in task.errorReport.getErrors() if e["code"] == 990]


def test_guard_runs_between_process_output_files_and_gleaning(tmp_path):
    """postProcess must not glean a file the guard rejected."""
    from pathlib import Path
    from ccp4i2.core.base_object.error_reporting import CErrorReport

    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = Path(tmp_path)
    plugin.errorReport = CErrorReport()
    plugin.LOG_FAILURES = ()
    plugin.TASKNAME = 'testtask'
    plugin.TASKCOMMAND = 'testprog'
    plugin._exitCode = 0
    plugin._exitStatus = 0
    plugin.reported, plugin.gleaned = [], []
    plugin.processOutputFiles = lambda: CPluginScript.SUCCEEDED
    plugin.checkMiniMtzOutputs = lambda: CPluginScript.FAILED
    plugin.reportStatus = lambda status: plugin.reported.append(status) or status
    plugin._glean_output_files_sync = lambda: plugin.gleaned.append(True)

    assert plugin.postProcess() == CPluginScript.FAILED
    assert plugin.gleaned == []
    assert plugin.reported == [CPluginScript.FAILED]
