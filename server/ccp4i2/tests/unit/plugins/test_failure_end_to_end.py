"""
Deliberate failures, run for real: does the cause reach the user?

Every other i2run and unit test in this tree asserts on a *successful* run, so
the error-surfacing machinery -- C1's status, C6's descriptions, C7's evidence,
C3's attribution -- was only ever exercised by accident, when something broke.
These tests break things on purpose.

They are fast by construction: a job that fails, fails early. The "program" is
this interpreter exiting non-zero, so there is no crystallography, no CCP4
binary and nothing to download; each test is a few hundred milliseconds. The
point is the plumbing between a failure and the panel, which is identical
whether the process that died was ``python -c`` or refmac.
"""
import sys
from pathlib import Path

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript


def _script(body):
    """Command-line arguments that make this interpreter behave like a program."""
    return ['-c', body]


class _ExitsNonZero(CPluginScript):
    TASKNAME = 'test_exits_non_zero'
    TASKCOMMAND = sys.executable

    def makeCommandAndScript(self):
        self.appendCommandLine(_script(
            'import sys; sys.stderr.write("FATAL ERROR: cell is nonsense\\n"); sys.exit(3)'))
        return CPluginScript.SUCCEEDED


class _RefusesItsInput(CPluginScript):
    TASKNAME = 'test_refuses_its_input'
    TASKCOMMAND = sys.executable

    def processInputFiles(self):
        self.appendErrorReport(101, 'XYZIN is a cif, and this program reads pdb')
        return CPluginScript.FAILED

    def makeCommandAndScript(self):
        self.appendCommandLine(_script('pass'))
        return CPluginScript.SUCCEEDED


class _RunsButWritesNothing(CPluginScript):
    """Exits 0 and produces no output --- the case C1 exists for."""
    TASKNAME = 'test_runs_but_writes_nothing'
    TASKCOMMAND = sys.executable

    def makeCommandAndScript(self):
        self.appendCommandLine(_script('pass'))
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.appendErrorReport(101, 'the program exited 0 but wrote no XMLOUT')
        return CPluginScript.FAILED


class _RaisesInAHook(CPluginScript):
    TASKNAME = 'test_raises_in_a_hook'
    TASKCOMMAND = sys.executable

    def makeCommandAndScript(self):
        self.appendCommandLine(_script('pass'))
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        raise RuntimeError('Incorrect file format (perhaps it is cif not pdb?)')


def _run_as_subjob(tmp_path, child_class, label='job_1'):
    """Run *child_class* as job_N of a parent pipeline; return (parent, child)."""
    parent = CPluginScript(workDirectory=str(tmp_path), name='parent')
    child_dir = Path(tmp_path) / label
    child_dir.mkdir(exist_ok=True)
    child = child_class(parent=parent, workDirectory=str(child_dir), name='child')
    child.process()
    return parent, child


def _absorbed(parent):
    return [e for e in parent.errorReport.entries()
            if str(e['name']).startswith('job_')]


@pytest.mark.parametrize('child_class', [
    _ExitsNonZero, _RefusesItsInput, _RunsButWritesNothing, _RaisesInAHook,
])
def test_every_way_of_failing_reaches_the_parent(tmp_path, child_class):
    parent, child = _run_as_subjob(tmp_path, child_class)

    assert child._status == CPluginScript.FAILED, 'the child should know it failed'
    absorbed = _absorbed(parent)
    assert absorbed, f'{child_class.__name__} failed silently as far as its parent knew'
    assert absorbed[0]['name'].startswith('job_1')


@pytest.mark.parametrize('child_class', [
    _ExitsNonZero, _RefusesItsInput, _RunsButWritesNothing, _RaisesInAHook,
])
def test_every_way_of_failing_leaves_a_diagnostic_on_disk(tmp_path, child_class):
    _parent, child = _run_as_subjob(tmp_path, child_class)
    assert (Path(child.workDirectory) / 'diagnostic.xml').exists()


def test_the_programs_own_words_come_up_with_the_failure(tmp_path):
    """C7: what stderr said, not just that the exit code was non-zero."""
    parent, _child = _run_as_subjob(tmp_path, _ExitsNonZero)
    details = ' '.join(e['details'] for e in _absorbed(parent))
    assert 'exited with code 3' in details
    assert 'cell is nonsense' in details


def test_a_wrappers_own_verdict_survives_the_trip(tmp_path):
    """C1: processOutputFiles said FAILED, and said why."""
    parent, _child = _run_as_subjob(tmp_path, _RunsButWritesNothing)
    details = ' '.join(e['details'] for e in _absorbed(parent))
    assert 'wrote no XMLOUT' in details


def test_an_exception_arrives_as_its_message_not_as_a_status(tmp_path):
    parent, _child = _run_as_subjob(tmp_path, _RaisesInAHook)
    entries = _absorbed(parent)
    joined = ' '.join(e['details'] + e.get('stack', '') for e in entries)
    assert 'cif not pdb' in joined


def test_a_cause_two_levels_down_still_arrives_flattened(tmp_path):
    """The shape of the SubstituteLigand failure: pipeline -> pipeline -> program."""
    top = CPluginScript(workDirectory=str(tmp_path), name='top')
    middle_dir = Path(tmp_path) / 'job_1'
    middle_dir.mkdir()
    middle = CPluginScript(parent=top, workDirectory=str(middle_dir), name='middle')

    bottom_dir = middle_dir / 'job_2'
    bottom_dir.mkdir()
    bottom = _ExitsNonZero(parent=middle, workDirectory=str(bottom_dir), name='bottom')
    bottom.process()
    middle.reportStatus(CPluginScript.FAILED)

    names = [e['name'] for e in top.errorReport.entries()]
    assert any(n.startswith('job_1/job_2') for n in names), names
    assert all('\n' not in n for n in names), 'the path should be flat, not nested'
