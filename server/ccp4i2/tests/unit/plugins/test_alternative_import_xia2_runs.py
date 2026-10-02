"""The xia2 import read its runs from runSummaries, which the Qt interface
filled by scanning the directory and the new one leaves empty; and it chose
what to look for by the run's name ("3d..", "2d..", "dials.."), so any other
name gave no pattern and a TypeError. It now finds the runs itself."""
from ccp4i2.core.tasks import get_plugin_class


def _plugin(tmp_path, directory):
    plugin = get_plugin_class("AlternativeImportXIA2")(workDirectory=str(tmp_path), name="imp")
    plugin.container.inputData.XIA2_DIRECTORY.setFullPath(str(directory))
    return plugin


def test_one_run_directory(tmp_path):
    run = tmp_path / "job_1"
    (run / "DataFiles").mkdir(parents=True)
    assert _plugin(tmp_path, run).xia2Runs() == [("job_1", str(run))]


def test_parent_of_several_runs(tmp_path):
    for name in ("dials-run", "3dii-run"):
        (tmp_path / "xia2" / name / "DataFiles").mkdir(parents=True)
    (tmp_path / "xia2" / "notarun").mkdir()
    runs = _plugin(tmp_path, tmp_path / "xia2").xia2Runs()
    assert [name for name, _ in runs] == ["3dii-run", "dials-run"]
