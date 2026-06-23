"""
Behavioural guard for the `ccp4_free` task flag.

A task flagged `ccp4_free=True` in the registry claims it can run inline on a
CCP4-free server (the env-aware run_local feasibility check trusts that flag —
see lib/utils/jobs/context_run.can_run_local). This guard makes the flag
un-lie-able: if a flagged task would actually touch a CCP4 binary, CI goes red.

Two complementary checks:

1. **Static** — every `ccp4_free` task imports with no CCP4 present and declares
   no `TASKCOMMAND`. The plugin framework executes `TASKCOMMAND` via subprocess,
   so a flagged task naming an external binary is a contradiction. This catches
   the "binary wrapper" class (e.g. freerflag, were it ever mis-flagged).

2. **Behavioural** — a subprocess/`shutil.which` sentinel that raises if invoked,
   used to prove the direct-call class (a helper that shells out without a
   TASKCOMMAND — the matthewsCoeff/411 shape) really is binary-free. The sentinel
   has a self-test so the guard itself is trustworthy.

All CCP4-free; runs in CI on every platform.
"""
import contextlib
import os
import subprocess
import shutil
import sys

import pytest

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
import django  # noqa: E402

django.setup()

from ccp4i2.core.tasks import TASKS, get_plugin_class  # noqa: E402

CCP4_FREE_TASKS = sorted(n for n, t in TASKS.items() if getattr(t, "ccp4_free", False))


class CCP4BinaryInvoked(AssertionError):
    """Raised by the sentinel when guarded code tries to launch a process."""


@contextlib.contextmanager
def forbid_binaries():
    """Make any subprocess launch / binary probe raise CCP4BinaryInvoked.

    A genuinely CCP4-free code path uses gemmi/python only and never reaches
    here; anything that shells out trips the sentinel.
    """
    def boom(*_a, **_k):
        raise CCP4BinaryInvoked("guarded code attempted to launch a subprocess / probe a binary")

    saved = {}
    targets = [
        (subprocess, "run"), (subprocess, "Popen"), (subprocess, "call"),
        (subprocess, "check_call"), (subprocess, "check_output"),
        (os, "system"), (shutil, "which"),
    ]
    if hasattr(os, "posix_spawn"):
        targets.append((os, "posix_spawn"))
    for mod, attr in targets:
        saved[(mod, attr)] = getattr(mod, attr)
        setattr(mod, attr, boom)
    try:
        yield
    finally:
        for (mod, attr), orig in saved.items():
            setattr(mod, attr, orig)


def test_there_are_ccp4_free_tasks():
    # If this fails the parametrised guards below would vacuously pass.
    assert CCP4_FREE_TASKS, "no ccp4_free tasks registered — guard would be vacuous"


def test_sentinel_self_check():
    """The sentinel must actually fire — otherwise the behavioural guard is worthless."""
    with pytest.raises(CCP4BinaryInvoked):
        with forbid_binaries():
            subprocess.run(["true"])
    with pytest.raises(CCP4BinaryInvoked):
        with forbid_binaries():
            shutil.which("ls")
    # ... and it must be cleanly removed afterwards
    assert shutil.which("python3") or shutil.which("python") or sys.executable


@pytest.mark.parametrize("task_name", CCP4_FREE_TASKS)
def test_ccp4_free_task_imports_without_ccp4(task_name, monkeypatch):
    monkeypatch.delenv("CCP4", raising=False)
    cls = get_plugin_class(task_name)
    assert cls is not None, f"{task_name} is ccp4_free but failed to import CCP4-free"


@pytest.mark.parametrize("task_name", CCP4_FREE_TASKS)
def test_ccp4_free_task_declares_no_taskcommand(task_name):
    cls = get_plugin_class(task_name)
    taskcommand = getattr(cls, "TASKCOMMAND", None)
    assert not taskcommand, (
        f"{task_name} is flagged ccp4_free but declares TASKCOMMAND={taskcommand!r}; "
        "the framework would exec that binary, so it is not CCP4-free"
    )


def test_matthews_compute_path_is_binary_free():
    """Behavioural: the ProvideAsuContents compute path (matthewsCoeff) must run
    with no subprocess at all — proving it is gemmi-native, not silently finding
    a binary. This is the direct-call class the TASKCOMMAND check cannot see."""
    import glob
    from ccp4i2 import I2_TOP
    from ccp4i2.core.CCP4XtalData import CMtzDataFile

    mtz = glob.glob(str(I2_TOP / "demo_data" / "gamma" / "*.mtz"))[0]
    f = CMtzDataFile()
    f.setFullPath(mtz)
    f.loadFile()
    with forbid_binaries():
        rv = f.fileContent.matthewsCoeff(molWt=8000.0, polymerMode="P")
    assert rv["results"] and rv["cell_volume"] > 0
