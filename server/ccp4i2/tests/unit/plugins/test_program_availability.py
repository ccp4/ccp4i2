"""
Unit tests for the pre-run program-availability check.

Two axes decide whether a missing binary blocks submission or is merely
flagged:

* is this process the one that will run the job?  Desktop/Electron and i2run
  execute locally against a mounted CCP4, so "not found" is a certain failure.
  In Azure mode the job is queued for a worker whose filesystem we cannot see,
  and on the slim CCP4-free API server there is nothing to look in — in both
  the absence is our ignorance, not the user's misconfiguration.
* can we be sure the task launches the program?  AUXILIARY_PROGRAMS is opt-in,
  so it counts. TASKCOMMAND only counts for a plugin that leaves process() to
  the base class; one that overrides process() may never spawn it (crank2
  declares 'crank2.py' but runs in-process; buster declares 'refine' but
  sources $BUSTERDIR/setup.sh to find it).

Pure Python — no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import (
    CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING)


class _PlainWrapper(CPluginScript):
    """Leaves process() to the base class, so it really launches TASKCOMMAND."""
    TASKNAME = 'plainwrapper'
    TASKCOMMAND = 'no-such-program-xyzzy'


class _SelfDrivingPipeline(CPluginScript):
    """Overrides process(), so TASKCOMMAND may never be spawned (cf. crank2)."""
    TASKNAME = 'selfdriving'
    TASKCOMMAND = 'no-such-program-xyzzy'

    def process(self, container=None):
        return self.SUCCEEDED


class _NeedsHelpers(CPluginScript):
    """Declares a helper it drives from inside TASKCOMMAND (cf. arcimboldo)."""
    TASKNAME = 'needshelpers'
    TASKCOMMAND = 'no-such-program-xyzzy'
    AUXILIARY_PROGRAMS = ('no-such-helper-xyzzy',)

    def process(self, container=None):
        return self.SUCCEEDED


@pytest.fixture
def authoritative(monkeypatch):
    """Pretend we are a desktop launch with CCP4 mounted."""
    monkeypatch.setattr(
        'ccp4i2.lib.utils.jobs.context_run.program_checks_are_authoritative',
        lambda: True)


@pytest.fixture
def not_authoritative(monkeypatch):
    """Pretend the job will be run somewhere we cannot see."""
    monkeypatch.setattr(
        'ccp4i2.lib.utils.jobs.context_run.program_checks_are_authoritative',
        lambda: False)


def _check(cls):
    error = CErrorReport()
    cls.__new__(cls)._checkProgramAvailable(error)
    return error


def test_missing_taskcommand_blocks_when_we_will_run_the_job(authoritative):
    error = _check(_PlainWrapper)
    assert error.maxSeverity() >= SEVERITY_ERROR
    assert 'no-such-program-xyzzy' in error.getErrors()[0]['details']


def test_missing_taskcommand_is_advisory_when_someone_else_will(not_authoritative):
    """Azure worker / slim CCP4-free server: we cannot see the run host."""
    error = _check(_PlainWrapper)
    assert len(error.getErrors()) == 1
    assert error.maxSeverity() == SEVERITY_WARNING


def test_overriding_process_keeps_taskcommand_advisory(authoritative):
    """crank2/buster: TASKCOMMAND is not necessarily what gets spawned."""
    error = _check(_SelfDrivingPipeline)
    assert len(error.getErrors()) == 1
    assert error.maxSeverity() == SEVERITY_WARNING


def test_auxiliary_programs_block_even_when_process_is_overridden(authoritative):
    """AUXILIARY_PROGRAMS is opt-in, so declaring one asserts it is required."""
    error = _check(_NeedsHelpers)
    by_severity = {
        e['details'].split("'")[1]: e['severity'] for e in error.getErrors()}
    assert by_severity['no-such-helper-xyzzy'] == SEVERITY_ERROR
    assert by_severity['no-such-program-xyzzy'] == SEVERITY_WARNING


def test_auxiliary_programs_are_advisory_off_box(not_authoritative):
    error = _check(_NeedsHelpers)
    assert error.maxSeverity() == SEVERITY_WARNING


def test_program_that_exists_is_not_reported(authoritative):
    class _Found(CPluginScript):
        TASKNAME = 'found'
        TASKCOMMAND = 'ls'
    assert _check(_Found).getErrors() == []


def test_no_declared_programs_is_a_no_op(authoritative):
    class _Pure(CPluginScript):
        TASKNAME = 'pure'
    assert _check(_Pure).getErrors() == []


def test_check_never_raises_if_discovery_breaks(authoritative, monkeypatch):
    monkeypatch.setattr(
        'ccp4i2.config.program_discovery.resolve_program',
        lambda name: (_ for _ in ()).throw(RuntimeError('discovery broke')))
    assert _check(_PlainWrapper).getErrors() == []


def test_authority_helper_reflects_execution_mode(monkeypatch):
    from ccp4i2.lib.utils.jobs import context_run

    monkeypatch.setattr(context_run, 'ccp4_available', lambda: True)
    monkeypatch.setattr(context_run, 'get_execution_mode', lambda: 'local')
    assert context_run.program_checks_are_authoritative() is True

    monkeypatch.setattr(context_run, 'get_execution_mode', lambda: 'azure')
    assert context_run.program_checks_are_authoritative() is False

    # slim CCP4-free server: local mode, but nothing to look in
    monkeypatch.setattr(context_run, 'get_execution_mode', lambda: 'local')
    monkeypatch.setattr(context_run, 'ccp4_available', lambda: False)
    assert context_run.program_checks_are_authoritative() is False


# --- OPTIONAL_PROGRAMS: listed in Preferences, never blocking ---------------


class _MayUseHelpers(CPluginScript):
    """Declares programs it *may* run (cf. crank2's phasing routes)."""
    TASKNAME = 'mayusehelpers'
    OPTIONAL_PROGRAMS = ('no-such-optional-xyzzy',)

    def process(self, container=None):
        return self.SUCCEEDED


class _NeedsAndMayUse(CPluginScript):
    TASKNAME = 'needsandmayuse'
    AUXILIARY_PROGRAMS = ('no-such-helper-xyzzy',)
    OPTIONAL_PROGRAMS = ('no-such-optional-xyzzy',)

    def process(self, container=None):
        return self.SUCCEEDED


def test_a_missing_optional_program_never_blocks(authoritative):
    """Even where we are the ones who will run the job.

    crank2 drives shelxc, shelxd, shelxe, prasa and parrot, but which of them
    it needs depends on the phasing route chosen at run time. Declaring them
    as auxiliary would refuse runs that never touch SHELX.
    """
    error = _check(_MayUseHelpers)
    assert len(error) == 1
    assert error.maxSeverity() == SEVERITY_WARNING


def test_a_missing_optional_program_is_still_reported(authoritative):
    """Silence would leave the user guessing when the route does need it."""
    details = ' '.join(e['details'] for e in _check(_MayUseHelpers).entries())
    assert 'no-such-optional-xyzzy' in details
    assert 'Preferences' in details


def test_required_and_optional_keep_their_own_severities(authoritative):
    error = _check(_NeedsAndMayUse)
    by_name = {}
    for entry in error.entries():
        for program in ('no-such-helper-xyzzy', 'no-such-optional-xyzzy'):
            if program in entry['details']:
                by_name[program] = entry['severity']
    assert by_name['no-such-helper-xyzzy'] == SEVERITY_ERROR
    assert by_name['no-such-optional-xyzzy'] == SEVERITY_WARNING


def test_optional_programs_reach_the_preferences_page():
    """Membership of that page is about relocatability, not requirement."""
    from ccp4i2.core.tasks import task_commands

    declared = task_commands()
    for program in ('shelxc', 'shelxd', 'shelxe', 'prasa', 'parrot'):
        assert 'crank2' in declared.get(program, []), \
            f'{program} is not offered on the Program locations page'


def test_crank2_no_longer_names_a_program_that_does_not_exist():
    from ccp4i2.core.tasks import task_commands

    assert 'crank2.py' not in task_commands(), \
        'the phantom is back; twelve tasks inherit this declaration'
