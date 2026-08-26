"""
Every registered task must actually be usable.

G1 of docs/error-handling-remediation.md. `core/tasks.py` is a single
enumerable list, so these questions can be asked of all of it at once, in a
second, without running any program:

- does the class import?
- does it declare the TASKNAME it is registered under?
- is the def.xml it points at really there?

Three defects were found by asking, none of which had any way to report itself
(C8: a failure before a job exists has nowhere to write a diagnostic):

- ``nucleofind`` declared no ``TASKNAME``, so ``CPluginScript.__init__`` never
  looked for its def.xml. Its container came up empty, i2run built a parser
  with no task arguments, and every run ended in ``unrecognized arguments:
  --FPHIIN``. The task had never run.
- ``SIMBAD`` guarded an optional import with ``except ImportError`` while the
  library raises ``RuntimeError`` when ``$CCP4`` is unset, so the task failed
  to load under any ``ccp4-python`` that had not sourced ``ccp4.setup-sh``.
- ``pisa`` has no ``pluginPath`` at all --- it names a def.xml for one of
  pisapipe's components --- and asking for its class reported an import error
  reading ``'NoneType' object has no attribute 'split'``.

No CCP4 binaries needed: all of it must hold on the slim server too.
"""
import os

import pytest

from ccp4i2.core.tasks import (
    TASKS,
    get_plugin_class,
    get_import_error,
    get_report_class,
    locate_def_xml,
)

WITH_PLUGIN = sorted(name for name, task in TASKS.items() if task.pluginPath)
WITH_DEF_XML = sorted(name for name, task in TASKS.items() if task.defXmlPath)
WITH_REPORT = sorted(name for name, task in TASKS.items() if task.reportPath)


def test_the_registry_is_not_empty():
    """A guard on the guard: parametrised tests over an empty list pass."""
    assert len(TASKS) > 150
    assert len(WITH_PLUGIN) > 150


@pytest.mark.parametrize("name", sorted(TASKS))
def test_every_task_declares_something_to_load(name):
    task = TASKS[name]
    assert task.pluginPath or task.defXmlPath, (
        f"{name} is registered but declares neither a plugin nor a def.xml"
    )


@pytest.mark.parametrize("name", WITH_PLUGIN)
def test_every_plugin_class_imports(name):
    assert get_plugin_class(name) is not None, (
        f"{name} failed to import: {get_import_error(name)}"
    )


@pytest.mark.parametrize("name", WITH_PLUGIN)
def test_every_plugin_declares_the_taskname_it_is_registered_under(name):
    """TASKNAME is how __init__ finds the def.xml. Get it wrong and the task
    loads no parameters at all, silently."""
    cls = get_plugin_class(name)
    if cls is None:
        pytest.skip("import failure is reported by its own test")
    assert getattr(cls, "TASKNAME", None) == name


@pytest.mark.parametrize("name", WITH_DEF_XML)
def test_every_declared_def_xml_is_on_disk(name):
    path = locate_def_xml(name)
    assert path and os.path.exists(str(path)), (
        f"{name} declares {TASKS[name].defXmlPath}, which is not there"
    )


@pytest.mark.parametrize("name", WITH_REPORT)
def test_every_declared_report_class_imports(name):
    assert get_report_class(name) is not None, (
        f"{name} declares a report class that will not import"
    )


def test_a_task_without_a_plugin_path_is_not_an_import_failure():
    """`pisa` names a def.xml for a pisapipe component and has no class. That
    is an answer, not an error."""
    no_plugin = [n for n, t in TASKS.items() if not t.pluginPath]
    for name in no_plugin:
        assert get_plugin_class(name) is None
        assert get_import_error(name) is None, (
            f"{name} declares no plugin, so asking for one is not a failure"
        )
