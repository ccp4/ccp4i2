"""
Unit tests for LOG_FAILURES -- declarative log scanning (mechanism M2 of
docs/error-handling-remediation.md).

Many CCP4 programs report a fatal error to their log and then exit 0, so a
clean exit code is not evidence of success. The motivating case is ARCIMBOLDO:
``ARCIMBOLDO_LITE.main()`` wraps the whole run in ``except SystemExit: pass``,
so a ``sys.exit(1)`` deep in SELSLIB2 (missing shelxe, missing phaser) reaches
the wrapper as exit code 0 with an empty stderr. Before this mechanism such a
job was recorded as Finished with no outputs and no explanation.

Pure Python -- no CCP4 binaries needed.
"""
from pathlib import Path

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR


# The real message ARCIMBOLDO wrote to log.txt, ANSI colour codes and all.
ARCIMBOLDO_FATAL = (
    "\x1b[31mFATAL:\x1b[0m The path given for shelxe: /opt/ccp4/bin/shelxe "
    "does not exist or it is not accessible by the user: nmemn"
)
ARCIMBOLDO_PATTERN = r'^\s*FATAL(?:\s+ERROR)?\b[:\s]*(.*)$'


def _plugin(tmp_path, log_text=None, log_failures=(), exit_code=0):
    """A CPluginScript far enough constructed to run postProcessCheck()."""
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = Path(tmp_path)
    plugin.errorReport = CErrorReport()
    plugin.LOG_FAILURES = log_failures
    plugin.TASKNAME = 'testtask'
    plugin.TASKCOMMAND = 'testprog'
    plugin._exitCode = exit_code
    plugin._exitStatus = 0 if exit_code == 0 else 1
    if log_text is not None:
        (Path(tmp_path) / 'log.txt').write_text(log_text, encoding='utf-8')
    return plugin


def test_no_log_failures_declared_is_a_no_op(tmp_path):
    """A task that declares nothing keeps the old behaviour exactly."""
    plugin = _plugin(tmp_path, log_text=ARCIMBOLDO_FATAL, log_failures=())
    status, exit_status, exit_code = plugin.postProcessCheck()
    assert status == CPluginScript.SUCCEEDED
    assert plugin.errorReport.maxSeverity() < SEVERITY_ERROR


def test_fatal_in_log_fails_a_job_that_exited_zero(tmp_path):
    """The whole point: exit code 0, but the log says the run died."""
    plugin = _plugin(
        tmp_path,
        log_text="some output\n" + ARCIMBOLDO_FATAL + "\n",
        log_failures=((ARCIMBOLDO_PATTERN, 301, None),),
        exit_code=0,
    )
    status, exit_status, exit_code = plugin.postProcessCheck()
    assert status == CPluginScript.FAILED
    assert exit_status == 1
    assert exit_code == 0                       # the program really did exit 0
    assert plugin.errorReport.maxSeverity() >= SEVERITY_ERROR
    errors = plugin.errorReport.getErrors()
    assert len(errors) == 1
    assert errors[0]['code'] == 301
    # details=None means "use what the pattern captured", ANSI stripped
    assert errors[0]['details'].startswith('The path given for shelxe:')
    assert '\x1b' not in errors[0]['details']


def test_clean_log_is_not_flagged(tmp_path):
    """Near-misses must not fail a good run."""
    plugin = _plugin(
        tmp_path,
        log_text=(
            "EXIT STATUS: SUCCESS\n"
            "FATALITY rates were tabulated\n"        # FATAL is not a whole word
            "this line only mentions FATAL midway\n"  # pattern is anchored
        ),
        log_failures=((ARCIMBOLDO_PATTERN, 301, None),),
    )
    status, _, _ = plugin.postProcessCheck()
    assert status == CPluginScript.SUCCEEDED
    assert plugin.scanLogForFailures() == []


def test_missing_log_file_is_not_a_failure(tmp_path):
    plugin = _plugin(tmp_path, log_text=None,
                     log_failures=((ARCIMBOLDO_PATTERN, 301, None),))
    status, _, _ = plugin.postProcessCheck()
    assert status == CPluginScript.SUCCEEDED


def test_explicit_details_override_the_matched_text(tmp_path):
    """The (pattern, code, message) form gives the user actionable advice."""
    plugin = _plugin(
        tmp_path,
        log_text="minimum or no description for ligand LIG\n",
        log_failures=((r'minimum or no description', 302,
                       'No geometry dictionary for a ligand -- run Make Ligand first'),),
    )
    status, _, _ = plugin.postProcessCheck()
    assert status == CPluginScript.FAILED
    assert plugin.errorReport.getErrors()[0]['details'] == (
        'No geometry dictionary for a ligand -- run Make Ligand first')


def test_repeated_message_is_reported_once(tmp_path):
    plugin = _plugin(
        tmp_path,
        log_text=(ARCIMBOLDO_FATAL + "\n") * 5,
        log_failures=((ARCIMBOLDO_PATTERN, 301, None),),
    )
    plugin.postProcessCheck()
    assert len(plugin.errorReport.getErrors()) == 1


def test_malformed_declaration_does_not_break_the_run(tmp_path):
    """A bad pattern must not replace the program's result with our crash."""
    plugin = _plugin(
        tmp_path,
        log_text=ARCIMBOLDO_FATAL + "\n",
        log_failures=(
            ('(unclosed group', 300, None),   # invalid regex
            ('too', 'few'),                   # malformed tuple
            (ARCIMBOLDO_PATTERN, 301, None),  # the good one still applies
        ),
    )
    status, _, _ = plugin.postProcessCheck()
    assert status == CPluginScript.FAILED
    assert [e['code'] for e in plugin.errorReport.getErrors()] == [301]


def test_nonzero_exit_still_reported_without_scanning(tmp_path):
    """An honest non-zero exit is unaffected by the mechanism."""
    plugin = _plugin(tmp_path, log_text="all fine\n",
                     log_failures=((ARCIMBOLDO_PATTERN, 301, None),),
                     exit_code=1)
    status, exit_status, exit_code = plugin.postProcessCheck()
    assert status == CPluginScript.FAILED
    assert exit_code == 1


def test_report_hook_receives_the_failures(tmp_path):
    """reportLogFailures() is how a wrapper gets them into its report XML."""
    seen = []

    plugin = _plugin(tmp_path, log_text=ARCIMBOLDO_FATAL + "\n",
                     log_failures=((ARCIMBOLDO_PATTERN, 301, None),))
    plugin.reportLogFailures = seen.extend
    plugin.postProcessCheck()
    assert len(seen) == 1
    assert seen[0]['code'] == 301
    assert seen[0]['line'] == 1


def test_report_hook_raising_does_not_mask_the_failure(tmp_path):
    def boom(failures):
        raise RuntimeError('report generation broke')

    plugin = _plugin(tmp_path, log_text=ARCIMBOLDO_FATAL + "\n",
                     log_failures=((ARCIMBOLDO_PATTERN, 301, None),))
    plugin.reportLogFailures = boom
    status, _, _ = plugin.postProcessCheck()
    assert status == CPluginScript.FAILED


def test_arcimboldo_declares_the_mechanism():
    """The motivating wrapper is wired up (and to the right helper binaries)."""
    from ccp4i2.wrappers.arcimboldo.script.arcimboldo import arcimboldo
    assert arcimboldo.LOG_FAILURES
    assert set(arcimboldo.AUXILIARY_PROGRAMS) == {'phaser', 'shelxe'}
