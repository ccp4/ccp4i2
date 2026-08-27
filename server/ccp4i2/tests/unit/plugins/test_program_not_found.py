"""
A program that is not installed should say so, once, in the right place.

Reconstructed from a real failure: ShelxCD, job 63, with SHELX installed but
outside the CCP4 tree in use. The user saw a single orange WARNING reading
"ShelxCD failed in mtz2various" whose details were a shelxc command line ---
and mtz2various had in fact succeeded, having written a 520 KB .sca file.

Five things went wrong between "shelxc is not on PATH" and that report, and
each has a test here.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript, SHARED_ERROR_CODES
from ccp4i2.core.base_object.error_reporting import (
    CErrorReport, CException, SEVERITY_ERROR, SEVERITY_WARNING,
)


def _plugin(tmp_path):
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.workDirectory = tmp_path
    plugin.errorReport = CErrorReport()
    plugin.TASKNAME = 'testtask'
    return plugin


# --- 1. getCommand honours program preferences ------------------------------

def test_get_command_returns_the_resolved_path(tmp_path, monkeypatch):
    """It used to return the bare name, whatever the preferences said."""
    monkeypatch.setattr('ccp4i2.config.program_discovery.resolve_program',
                        lambda name: f'/somewhere/else/bin/{name}')
    assert _plugin(tmp_path).getCommand('shelxc') == '/somewhere/else/bin/shelxc'


def test_get_command_refuses_rather_than_returning_something_unusable(tmp_path, monkeypatch):
    monkeypatch.setattr('ccp4i2.config.program_discovery.resolve_program',
                        lambda name: None)
    plugin = _plugin(tmp_path)
    with pytest.raises(CException):
        plugin.getCommand('shelxc')


def test_the_refusal_names_the_program_and_where_to_set_it(tmp_path, monkeypatch):
    monkeypatch.setattr('ccp4i2.config.program_discovery.resolve_program',
                        lambda name: None)
    plugin = _plugin(tmp_path)
    with pytest.raises(CException):
        plugin.getCommand('shelxc')
    entry = plugin.errorReport.entries()[-1]
    assert entry['code'] == 299
    assert 'shelxc' in entry['details']
    assert 'Preferences' in entry['details']
    assert entry['severity'] == SEVERITY_ERROR


# --- 2. the process manager will not launch what it cannot name -------------

def test_the_process_manager_refuses_a_missing_program():
    """None used to reach subprocess and fail there, talking about types."""
    from ccp4i2.core.CCP4ProcessManager import CProcessManager

    with pytest.raises(ValueError) as excinfo:
        CProcessManager().startProcess(None, ['result'])
    assert 'not installed' in str(excinfo.value)


def test_the_process_manager_refuses_an_empty_program_name():
    from ccp4i2.core.CCP4ProcessManager import CProcessManager

    with pytest.raises(ValueError):
        CProcessManager().startProcess('   ', ['result'])


# --- 4. a failed job cannot show only warnings ------------------------------

def test_a_failed_job_gains_an_error_when_it_had_only_warnings(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=201, details='something odd',
                              severity=SEVERITY_WARNING)

    plugin.recordThatItFailed()

    assert plugin.errorReport.maxSeverity() >= SEVERITY_ERROR
    added = plugin.errorReport.entries()[-1]
    assert added['code'] == 991
    assert 'something odd' in added['details'], 'quote what was actually recorded'


def test_what_was_already_recorded_is_not_altered(tmp_path):
    """Promoting someone else's advisory would be a lie in the other direction."""
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=200, details='free R recommended',
                              severity=SEVERITY_WARNING)

    plugin.recordThatItFailed()

    assert plugin.errorReport.entries()[0]['severity'] == SEVERITY_WARNING


def test_a_job_that_already_explained_itself_gains_nothing(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.errorReport.append(klass='t', code=213, details='shelxc exited 1',
                              severity=SEVERITY_ERROR)

    plugin.recordThatItFailed()

    assert len(plugin.errorReport) == 1


def test_a_silent_failure_says_where_to_look(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.recordThatItFailed()
    details = plugin.errorReport.entries()[-1]['details']
    assert 'recorded no reason' in details
    assert 'log' in details.lower()


# --- the codes describe themselves ------------------------------------------

@pytest.mark.parametrize('code', [299, 991])
def test_the_new_codes_have_descriptions(code):
    """Otherwise the panel shows a number and the specifics, with no sentence."""
    assert SHARED_ERROR_CODES[code]['description']


# --- 3 & 5. the wrapper declares what it runs, and names the right step ------

def test_shelxcd_declares_the_program_it_drives_itself():
    from ccp4i2.core.tasks import get_plugin_class

    cls = get_plugin_class('ShelxCD')
    if cls is None:
        pytest.skip('ShelxCD is not available in this environment')
    assert 'shelxc' in (cls.AUXILIARY_PROGRAMS or ()), \
        'nothing declared shelxc, so the pre-run check never looked for it'


def test_shelxcd_has_a_code_for_the_step_that_actually_failed():
    from ccp4i2.core.tasks import get_plugin_class

    cls = get_plugin_class('ShelxCD')
    if cls is None:
        pytest.skip('ShelxCD is not available in this environment')
    assert 'shelxc' in cls.ERROR_CODES[213]['description']
    assert 'mtz2various' in cls.ERROR_CODES[201]['description'], \
        '201 keeps its own meaning; 213 stops it being used for shelxc'


def test_get_command_defaults_to_taskcommand_as_it_always_did(tmp_path, monkeypatch):
    """`self.getCommand()` with no argument is the Qt-era spelling."""
    monkeypatch.setattr('ccp4i2.config.program_discovery.resolve_program',
                        lambda name: f'/opt/{name}')
    plugin = _plugin(tmp_path)
    plugin.TASKCOMMAND = 'shelxd'
    assert plugin.getCommand() == '/opt/shelxd'


def test_get_command_with_nothing_to_go_on_says_so(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.TASKCOMMAND = None
    with pytest.raises(CException) as excinfo:
        plugin.getCommand()
    assert 'TASKCOMMAND' in str(excinfo.value)
