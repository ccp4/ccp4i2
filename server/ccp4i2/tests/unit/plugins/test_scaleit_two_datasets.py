"""SCALEIT needs two datasets, and says so before the job runs (#783).

The task could be run with one dataset and failed inside processInputFiles,
an error shown only under "other errors" because it named no field. The
list now declares listMinLength 2, so a new job starts with two rows (a
CList is made with its minimum), each row must be set (named at the row),
and fewer than two rows is refused at the list.
"""
from pathlib import Path

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _blocking(error):
    return sorted(str(r.get("name", "")).split(".container.")[-1] for r in error._reports
                  if r["severity"] >= CCP4ErrorHandling.SEVERITY_ERROR)


def test_a_new_job_has_two_rows_that_must_both_be_set(tmp_path):
    plugin = get_plugin_class("scaleit")(workDirectory=str(tmp_path), name="scaleit")
    files = plugin.container.inputData.MERGEDFILES
    assert files.get_qualifier("listMinLength") == 2 and len(files) == 2
    assert _blocking(plugin.validity()) == ["inputData.MERGEDFILES[0]", "inputData.MERGEDFILES[1]"]
    files[0].setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    assert _blocking(plugin.validity()) == ["inputData.MERGEDFILES[1]"]
    files[1].setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    assert _blocking(plugin.validity()) == []


def test_fewer_than_two_rows_is_refused_at_the_list(tmp_path):
    plugin = get_plugin_class("scaleit")(workDirectory=str(tmp_path), name="scaleit")
    files = plugin.container.inputData.MERGEDFILES
    files.clear()
    files.append(files.makeItem())
    files[0].setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    assert _blocking(plugin.validity()) == ["inputData.MERGEDFILES"]
