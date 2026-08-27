"""
Unit tests for C3 of docs/error-handling-remediation.md -- a subjob's cause
must reach the job the user actually ran.

``SubstituteLigand`` recorded ``i2Dimple pipeline failed`` twice and nothing
else; the child's actual ``RuntimeError: Incorrect file format (perhaps it is
cif not pdb?)`` existed in no artefact on disk. Two things fix that: every
failing job writes its own ``diagnostic.xml`` in its own directory, and the
parent absorbs the child's causes, each tagged with the subjob it came from.

Pure Python -- no CCP4 binaries needed.
"""
import weakref
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import (
    CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING, SEVERITY_OK,
    write_diagnostic_xml,
)


def _plugin(directory, parent=None, name=None):
    """A CPluginScript far enough constructed to report its own failure."""
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = Path(directory)
    plugin.errorReport = CErrorReport()
    plugin.TASKNAME = name or 'testtask'
    plugin._parent_ref = weakref.ref(parent) if parent is not None else None
    plugin._objectName = name or 'testtask'
    return plugin


def _names(report):
    return [e['name'] for e in report.entries()]


# --- absorb(): the flattening contract --------------------------------------

def test_absorb_tags_the_cause_with_the_subjob():
    child = CErrorReport()
    child.append(klass='refmac', code=32, details='no XMLOUT written', name='XYZOUT')
    parent = CErrorReport()

    assert parent.absorb(child, 'job_2') == 1
    assert _names(parent) == ['job_2/XYZOUT']
    assert parent.entries()[0]['details'] == 'no XMLOUT written'
    assert parent.entries()[0]['class'] == 'refmac'


def test_absorb_keeps_the_label_alone_when_the_cause_has_no_name():
    child = CErrorReport()
    child.append(klass='dimple', code=101, details='segmentation fault')
    parent = CErrorReport()
    parent.absorb(child, 'job_1')
    assert _names(parent) == ['job_1']


def test_paths_accumulate_rather_than_nest():
    """Three deep still yields a flat list, each entry naming its whole path."""
    grandchild = CErrorReport()
    grandchild.append(klass='refmac', code=32, details='exited 11', name='XYZIN')
    child = CErrorReport()
    child.absorb(grandchild, 'job_4')
    parent = CErrorReport()
    parent.absorb(child, 'job_1')

    assert _names(parent) == ['job_1/job_4/XYZIN']
    assert len(parent) == 1


def test_absorb_leaves_behind_entries_below_the_threshold():
    child = CErrorReport()
    child.append(klass='t', code=1, details='chatter', severity=SEVERITY_OK)
    child.append(klass='t', code=2, details='real', severity=SEVERITY_ERROR)
    parent = CErrorReport()
    assert parent.absorb(child, 'job_2') == 1
    assert parent.entries()[0]['details'] == 'real'


def test_warnings_are_carried_up():
    child = CErrorReport()
    child.append(klass='t', code=1, details='low completeness', severity=SEVERITY_WARNING)
    parent = CErrorReport()
    assert parent.absorb(child, 'job_2') == 1


def test_absorbing_twice_does_not_duplicate():
    """A child can be reported by both its own reportStatus and the parent's check."""
    child = CErrorReport()
    child.append(klass='refmac', code=32, details='exited 11', name='XYZIN')
    parent = CErrorReport()
    parent.absorb(child, 'job_2')
    parent.absorb(child, 'job_2')
    assert len(parent) == 1


def test_the_same_cause_from_two_different_subjobs_is_kept_twice():
    child = CErrorReport()
    child.append(klass='refmac', code=32, details='exited 11', name='XYZIN')
    parent = CErrorReport()
    parent.absorb(child, 'job_2')
    parent.absorb(child, 'job_3')
    assert sorted(_names(parent)) == ['job_2/XYZIN', 'job_3/XYZIN']


def test_absorb_preserves_description_and_stack():
    child = CErrorReport()
    child.append(klass='t', code=1, details='d', description='what it means',
                 stack='Traceback...')
    parent = CErrorReport()
    parent.absorb(child, 'job_2')
    assert parent.entries()[0]['description'] == 'what it means'
    assert parent.entries()[0]['stack'] == 'Traceback...'


def test_absorb_of_a_non_report_is_a_no_op():
    parent = CErrorReport()
    assert parent.absorb(None, 'job_2') == 0
    assert parent.absorb(CErrorReport(), '') == 0


# --- write_diagnostic_xml() -------------------------------------------------

def test_diagnostic_xml_is_written_where_asked(tmp_path):
    report = CErrorReport()
    report.append(klass='dimple', code=101, details='no solution', name='XYZOUT')
    path = write_diagnostic_xml(report, tmp_path)

    assert path == str(tmp_path / 'diagnostic.xml')
    root = ET.parse(path).getroot()
    assert 'no solution' in ET.tostring(root, encoding='unicode')


def test_writing_a_diagnostic_never_raises(tmp_path):
    """A job that is already failing must not fail again over its diagnostics."""
    assert write_diagnostic_xml(CErrorReport(), tmp_path / 'no' / 'such' / 'dir') is None


# --- the plugin hook --------------------------------------------------------

def test_a_failing_subjob_writes_its_own_diagnostic(tmp_path):
    child_dir = tmp_path / 'job_2'
    child_dir.mkdir()
    child = _plugin(child_dir)
    child.errorReport.append(klass='dimple', code=101, details='cif not pdb')

    child.recordCauses(CPluginScript.FAILED)

    assert (child_dir / 'diagnostic.xml').exists()


def test_the_parent_learns_why_the_child_failed(tmp_path):
    parent = _plugin(tmp_path, name='SubstituteLigand')
    child_dir = tmp_path / 'job_2'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent, name='i2Dimple')
    child.errorReport.append(klass='i2Dimple', code=101,
                             details='Incorrect file format (perhaps it is cif not pdb?)')

    child.recordCauses(CPluginScript.FAILED)
    parent.recordCauses(CPluginScript.FAILED)

    assert _names(parent.errorReport) == ['job_2']
    assert 'cif not pdb' in parent.errorReport.entries()[0]['details']


def test_a_pipeline_that_carried_on_still_shows_what_it_carried_on_past(tmp_path):
    """Which pipelines *should* survive a failed subjob is not decided here.

    The failure is shown either way; what changes is whether it is fatal. A
    job that finished keeps its finished status and gains a warning, so the
    situation is visible and countable rather than silently either way.
    """
    parent = _plugin(tmp_path)
    child_dir = tmp_path / 'job_2'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent)
    child.errorReport.append(klass='t', code=1, details='first attempt failed',
                             severity=SEVERITY_ERROR)

    child.recordCauses(CPluginScript.FAILED)
    parent.recordCauses(CPluginScript.SUCCEEDED)

    assert (child_dir / 'diagnostic.xml').exists()
    entries = parent.errorReport.entries()
    assert [e['name'] for e in entries] == ['job_2']
    assert entries[0]['severity'] == SEVERITY_WARNING, (
        'a job that finished must not be marked failed by a step it survived')
    assert 'first attempt failed' in entries[0]['details']


def test_severity_only_ever_moves_down_as_a_cause_travels(tmp_path):
    top = _plugin(tmp_path)
    middle_dir = tmp_path / 'job_1'
    middle_dir.mkdir()
    middle = _plugin(middle_dir, parent=top)
    bottom_dir = middle_dir / 'job_2'
    bottom_dir.mkdir()
    bottom = _plugin(bottom_dir, parent=middle)
    bottom.errorReport.append(klass='t', code=1, details='boom')

    bottom.recordCauses(CPluginScript.FAILED)
    middle.recordCauses(CPluginScript.SUCCEEDED)   # middle survived it
    top.recordCauses(CPluginScript.FAILED)         # top did not

    entries = top.errorReport.entries()
    assert [e['name'] for e in entries] == ['job_1/job_2']
    assert entries[0]['severity'] == SEVERITY_WARNING, (
        'downgraded on the way past a job that survived, and not re-promoted')


def test_a_successful_job_does_not_repeat_its_own_advisories_upward(tmp_path):
    parent = _plugin(tmp_path)
    child_dir = tmp_path / 'job_2'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent)
    child.errorReport.append(klass='t', code=1, details='low completeness',
                             severity=SEVERITY_WARNING)

    child.recordCauses(CPluginScript.SUCCEEDED)

    assert len(parent.errorReport) == 0


def test_causes_are_recorded_once_however_often_asked(tmp_path):
    parent = _plugin(tmp_path)
    child_dir = tmp_path / 'job_2'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent)
    child.errorReport.append(klass='t', code=1, details='boom')

    child.recordCauses(CPluginScript.FAILED)
    child.recordCauses(CPluginScript.FAILED)
    parent.recordCauses(CPluginScript.FAILED)

    assert len(parent.errorReport) == 1


def test_a_top_level_job_still_writes_its_diagnostic(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=1, details='boom')
    plugin.recordCauses(CPluginScript.FAILED)
    assert (tmp_path / 'diagnostic.xml').exists()


def test_the_job_label_is_the_working_directory_name(tmp_path):
    child_dir = tmp_path / 'job_7'
    child_dir.mkdir()
    assert _plugin(child_dir).jobLabel() == 'job_7'


def test_parent_plugin_is_found_through_the_hierarchy(tmp_path):
    parent = _plugin(tmp_path)
    child = _plugin(tmp_path / 'job_1', parent=parent)
    assert child.parentPlugin() is parent
    assert parent.parentPlugin() is None


def test_report_status_records_causes_on_failure(tmp_path):
    """The funnel: every termination path goes through reportStatus."""
    parent = _plugin(tmp_path)
    child_dir = tmp_path / 'job_3'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent)
    child.errorReport.append(klass='t', code=1, details='boom')
    child.saveParams = lambda: None
    child.finished = type('S', (), {'emit': staticmethod(lambda *a: None)})()

    child.reportStatus(CPluginScript.FAILED)
    parent.recordCauses(CPluginScript.FAILED)

    assert (child_dir / 'diagnostic.xml').exists()
    assert _names(parent.errorReport) == ['job_3']


def test_report_status_stays_quiet_on_success(tmp_path):
    parent = _plugin(tmp_path)
    child_dir = tmp_path / 'job_3'
    child_dir.mkdir()
    child = _plugin(child_dir, parent=parent)
    child.errorReport.append(klass='t', code=1, details='an advisory',
                             severity=SEVERITY_WARNING)
    child.saveParams = lambda: None
    child.finished = type('S', (), {'emit': staticmethod(lambda *a: None)})()

    child.reportStatus(CPluginScript.SUCCEEDED)

    assert not (child_dir / 'diagnostic.xml').exists()
    assert len(parent.errorReport) == 0
