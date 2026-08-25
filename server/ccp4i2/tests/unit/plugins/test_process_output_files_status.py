"""
Unit tests for C1 of docs/error-handling-remediation.md -- a wrapper's
``processOutputFiles()`` verdict must decide the job's status.

The checks a wrapper makes in ``processOutputFiles()`` are the ones the
framework cannot make generically: *refmac exited 0 but wrote no XMLOUT*, *the
ligand had no dictionary*, *no solution was found*. Thirty-five wrappers make
such a check and return ``FAILED``. The synchronous path used to extend the
error report with the int (a silent no-op, since ``extend`` only accepts a
``CErrorReport``), print a warning to the server console, and record the job as
*Finished*.

Pure Python -- no CCP4 binaries needed.
"""
from pathlib import Path

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import (
    CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING,
)


def _plugin(tmp_path, outputResult=None, raises=None, exit_code=0):
    """A CPluginScript far enough constructed to run postProcess().

    ``processOutputFiles`` returns ``outputResult`` (or raises ``raises``);
    ``reportStatus`` and the gleaner are recorded rather than performed.
    """
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = Path(tmp_path)
    plugin.errorReport = CErrorReport()
    plugin.LOG_FAILURES = ()
    plugin.TASKNAME = 'testtask'
    plugin.TASKCOMMAND = 'testprog'
    plugin._exitCode = exit_code
    plugin._exitStatus = 0 if exit_code == 0 else 1

    plugin.reported = []
    plugin.gleaned = []

    def processOutputFiles():
        if raises is not None:
            raise raises
        return outputResult

    plugin.processOutputFiles = processOutputFiles
    plugin.reportStatus = lambda status: plugin.reported.append(status) or status
    plugin._glean_output_files_sync = lambda: plugin.gleaned.append(True)
    return plugin


def _codes(plugin):
    return [e['code'] for e in plugin.errorReport.getErrors()]


# --- the legacy int convention ---------------------------------------------

def test_failed_from_process_output_files_fails_the_job(tmp_path):
    """C1 itself: the program exited 0, the wrapper says its output is unusable."""
    plugin = _plugin(tmp_path, outputResult=CPluginScript.FAILED)
    assert plugin.postProcess() == CPluginScript.FAILED
    assert plugin.reported == [CPluginScript.FAILED]


def test_a_failure_with_no_message_still_says_which_hook_failed(tmp_path):
    """128 sites return FAILED without appending anything. Say something."""
    plugin = _plugin(tmp_path, outputResult=CPluginScript.FAILED)
    plugin.postProcess()
    assert 992 in _codes(plugin)
    details = plugin.errorReport.getErrors()[0]['details']
    assert 'processOutputFiles' in details
    assert plugin.errorReport.maxSeverity() >= SEVERITY_ERROR


def test_the_wrappers_own_message_is_not_replaced(tmp_path):
    """A wrapper that explained itself must not be given a generic message."""
    plugin = _plugin(tmp_path, outputResult=CPluginScript.FAILED)
    original = plugin.processOutputFiles

    def processOutputFiles():
        plugin.errorReport.append(
            klass='testtask', code=201, details='No XMLOUT was written',
            name='XMLOUT', severity=SEVERITY_ERROR,
        )
        return original()

    plugin.processOutputFiles = processOutputFiles
    plugin.postProcess()
    assert _codes(plugin) == [201]
    assert 992 not in _codes(plugin)


def test_succeeded_is_success(tmp_path):
    plugin = _plugin(tmp_path, outputResult=CPluginScript.SUCCEEDED)
    assert plugin.postProcess() == CPluginScript.SUCCEEDED
    assert _codes(plugin) == []


def test_returning_nothing_is_success(tmp_path):
    """Most implementations fall off the end. That has always meant success."""
    plugin = _plugin(tmp_path, outputResult=None)
    assert plugin.postProcess() == CPluginScript.SUCCEEDED


def test_unsatisfactory_is_neither_success_nor_failure(tmp_path):
    """Three wrappers return it; until C1 it could not be reached at all."""
    plugin = _plugin(tmp_path, outputResult=CPluginScript.UNSATISFACTORY)
    assert plugin.postProcess() == CPluginScript.UNSATISFACTORY
    assert plugin.reported == [CPluginScript.UNSATISFACTORY]


# --- the modern CErrorReport convention ------------------------------------

def test_an_error_report_fails_the_job_and_keeps_its_messages(tmp_path):
    error = CErrorReport()
    error.append(klass='testtask', code=42, details='FPHIOUT is empty',
                 name='FPHIOUT', severity=SEVERITY_ERROR)
    plugin = _plugin(tmp_path, outputResult=error)
    assert plugin.postProcess() == CPluginScript.FAILED
    assert _codes(plugin) == [42]


def test_a_report_of_warnings_is_reporting_not_refusing(tmp_path):
    """The threshold is ERROR (C5). Warnings are collected and shown."""
    error = CErrorReport()
    error.append(klass='testtask', code=43, details='PSOUT was not produced',
                 name='PSOUT', severity=SEVERITY_WARNING)
    plugin = _plugin(tmp_path, outputResult=error)
    assert plugin.postProcess() == CPluginScript.SUCCEEDED
    assert _codes(plugin) == [43]


# --- exceptions -------------------------------------------------------------

def test_an_exception_fails_the_job_and_records_the_traceback(tmp_path):
    """Symmetry with the other three hooks, which all fail on an exception."""
    plugin = _plugin(tmp_path, raises=KeyError('XMLOUT'))
    assert plugin.postProcess() == CPluginScript.FAILED
    assert 993 in _codes(plugin)
    details = plugin.errorReport.getErrors()[0]['details']
    assert 'KeyError' in details
    assert 'Traceback' in details


# --- what follows from the status ------------------------------------------

def test_a_failed_job_is_not_gleaned(tmp_path):
    """The outputs of a failed job must not become inputs to another."""
    plugin = _plugin(tmp_path, outputResult=CPluginScript.FAILED)
    plugin.postProcess()
    assert plugin.gleaned == []


def test_a_successful_job_is_gleaned(tmp_path):
    plugin = _plugin(tmp_path, outputResult=CPluginScript.SUCCEEDED)
    plugin.postProcess()
    assert plugin.gleaned == [True]


def test_a_program_that_failed_never_reaches_process_output_files(tmp_path):
    """A non-zero exit code is decided before the outputs are looked at."""
    reached = []
    plugin = _plugin(tmp_path, exit_code=1)
    plugin.processOutputFiles = lambda: reached.append(True)
    assert plugin.postProcess() == CPluginScript.FAILED
    assert reached == []


# --- the two conventions, in one place --------------------------------------

@pytest.mark.parametrize('result,expected', [
    (None, CPluginScript.SUCCEEDED),
    (CPluginScript.SUCCEEDED, CPluginScript.SUCCEEDED),
    (CPluginScript.FAILED, CPluginScript.FAILED),
    (CPluginScript.UNSATISFACTORY, CPluginScript.UNSATISFACTORY),
])
def test_absorb_hook_status_honours_the_int_convention(tmp_path, result, expected):
    plugin = _plugin(tmp_path)
    assert plugin.absorbHookStatus(result, 'processOutputFiles') == expected


def test_absorb_hook_status_does_not_consult_severity_for_an_int(tmp_path):
    """appendErrorReport defaults to WARNING (C2). FAILED still means failed."""
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='testtask', code=201, details='Exit code: 5',
                              name=None, severity=SEVERITY_WARNING)
    assert plugin.absorbHookStatus(CPluginScript.FAILED, 'processOutputFiles') \
        == CPluginScript.FAILED


# --- the synchronous path, which is where C1 actually bit --------------------

def _synchronous_plugin(tmp_path, outputResult):
    """A plugin whose pre-program hooks all succeed, so process() reaches its tail."""
    plugin = _plugin(tmp_path, outputResult=outputResult)
    plugin.runTimeValidity = lambda: CErrorReport()
    plugin.checkOutputData = lambda: CErrorReport()
    plugin.processInputFiles = lambda: CPluginScript.SUCCEEDED
    plugin.makeCommandAndScript = lambda: CErrorReport()
    plugin.startProcess = lambda: CPluginScript.SUCCEEDED
    plugin.get_db_job_id = lambda: None
    return plugin


def test_the_synchronous_path_fails_the_job_too(tmp_path):
    """The C1 defect proper.

    ``process()`` used to extend its error report with the int (a no-op),
    print a warning to the server console, and return SUCCEEDED -- the line
    that would have failed the job was commented out in the source.
    """
    plugin = _synchronous_plugin(tmp_path, CPluginScript.FAILED)
    assert plugin.process() == CPluginScript.FAILED
    assert plugin.reported == [CPluginScript.FAILED]
    assert plugin.gleaned == []
    assert 992 in _codes(plugin)


def test_the_synchronous_path_still_succeeds_when_the_wrapper_is_happy(tmp_path):
    plugin = _synchronous_plugin(tmp_path, CPluginScript.SUCCEEDED)
    assert plugin.process() == CPluginScript.SUCCEEDED
    assert plugin.gleaned == [True]


@pytest.mark.parametrize('result', [
    None,
    CPluginScript.SUCCEEDED,
    CPluginScript.FAILED,
    CPluginScript.UNSATISFACTORY,
])
def test_both_paths_agree(tmp_path, result):
    """process() and _onProcessFinished() run the same code, so they must agree."""
    synchronous = _synchronous_plugin(tmp_path, result)
    asynchronous = _plugin(tmp_path, outputResult=result)
    assert synchronous.process() == asynchronous.postProcess()
    assert synchronous.gleaned == asynchronous.gleaned
    assert _codes(synchronous) == _codes(asynchronous)
