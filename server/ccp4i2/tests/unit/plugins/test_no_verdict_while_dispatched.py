"""
A dispatched job's diagnostics draw no verdict.

recordCauses read "not SUCCEEDED" as "failed", so a plugin that returned
DISPATCHED -- the program handed to a run target, the job parked in
RUNNING_REMOTELY -- got an ERROR-severity "The job failed" (991) quoting its
own dispatch notice as the warning in question. The status said running
remotely, the diagnostics said failed, and an operator on Batch resubmitted a
healthy run (Materia, 2026-09-27). Work still going on is not a verdict.
"""
import xml.etree.ElementTree as ET

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING


def _plugin(tmp_path):
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = tmp_path
    plugin.errorReport = CErrorReport()
    plugin.TASKNAME = 'testtask'
    plugin._pendingCauses = []
    return plugin


def test_dispatched_leaves_no_error_and_keeps_the_notice_as_a_warning(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=228, details="PanDDA dispatched to run target 'batch'",
                              severity=SEVERITY_WARNING)
    plugin.recordCauses(CPluginScript.DISPATCHED)
    assert plugin.errorReport.maxSeverity() < SEVERITY_ERROR
    assert [e['code'] for e in plugin.errorReport.entries()] == [228]
    # What was recorded is on disk, as it stands.
    tree = ET.parse(tmp_path / "diagnostic.xml")
    codes = {r.findtext("code") for r in tree.iter("errorReport")}
    assert codes == {"228"}


def test_running_is_no_verdict_either(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.recordCauses(CPluginScript.RUNNING)
    assert len(plugin.errorReport) == 0


def test_the_verdict_still_records_when_it_comes(tmp_path):
    """A no-verdict call must not consume the once-only guard."""
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=228, details='dispatched', severity=SEVERITY_WARNING)
    plugin.recordCauses(CPluginScript.DISPATCHED)
    plugin.recordCauses(CPluginScript.FAILED)
    assert plugin.errorReport.maxSeverity() >= SEVERITY_ERROR
    assert plugin.errorReport.entries()[-1]['code'] == 991
