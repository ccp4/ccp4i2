"""
Run targets resolve from settings; with none, everything is local.

The desktop invariant first: an Electron dev tree or a packaged app sets
nothing and must get the shipped subprocess runner. Then the hook a
deployment uses to add a target CCP4i2 has never heard of, and the errors a
mis-registration produces. CCP4-free; no job is started.
"""
import types

import pytest
from django.conf import settings

from ccp4i2.lib import dispatch
from ccp4i2.lib.dispatch import base, local
from ccp4i2.lib.utils.jobs import context_run


class FakeQueueTarget:
    """A deployment's target: runs jobs by queueing them; no program dispatch."""
    calls = []

    def run_job(self, job, *, synchronous=False):
        self.calls.append((job, synchronous))
        return {"success": True, "data": job, "status": 200}


class FakeBatchTarget:
    """A program target only: axis B, no run_job."""
    def submit(self, tree, argv, out_dir, sizing_hint): return "h1"
    def poll(self, handle): return "queued"
    def cancel(self, handle): return None
    def logs(self, handle): return None


class Unconstructable:
    def __init__(self):
        raise RuntimeError("needs a connection string")


def _path(cls):
    return f"{__name__}.{cls.__name__}"


@pytest.fixture
def no_settings(monkeypatch):
    """The desktop: neither setting, no env."""
    monkeypatch.delattr(settings, "CCP4I2_RUN_TARGETS", raising=False)
    monkeypatch.delattr(settings, "CCP4I2_JOB_TARGET", raising=False)
    monkeypatch.delenv("CCP4I2_JOB_TARGET", raising=False)


@pytest.fixture
def fake_job():
    return types.SimpleNamespace(id=1, uuid="u-1", task_name="pointless")


def test_with_no_settings_everything_is_local(no_settings):
    assert dispatch.job_target_name() == "local"
    assert dispatch.run_target_paths() == dispatch.DEFAULT_RUN_TARGETS
    target = dispatch.get_target("local")
    assert isinstance(target, local.LocalTarget)
    assert isinstance(target, base.JobTarget)
    assert not isinstance(target, base.ProgramTarget)
    assert dispatch.available_targets() == [
        {"name": "local", "path": dispatch.DEFAULT_RUN_TARGETS["local"],
         "runs_jobs": True, "runs_programs": False}]


def test_the_shipped_local_target_finds_its_wrapper_scripts():
    """run_job_local launches through scripts/; the move must not lose them."""
    assert (local.SCRIPTS_DIR / "run_job_safe.sh").is_file()
    assert (local.SCRIPTS_DIR / "run_job_safe.cmd").is_file()


def test_a_deployment_registers_a_target_by_dotted_path(monkeypatch, no_settings):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS,
                         "queue": _path(FakeQueueTarget),
                         "batch": _path(FakeBatchTarget)}, raising=False)
    assert isinstance(dispatch.get_target("queue"), FakeQueueTarget)
    assert isinstance(dispatch.get_target("QUEUE"), FakeQueueTarget)  # names are case-insensitive
    caps = {t["name"]: (t["runs_jobs"], t["runs_programs"]) for t in dispatch.available_targets()}
    assert caps == {"local": (True, False), "queue": (True, False), "batch": (False, True)}


def test_env_names_the_job_target_when_settings_do_not(monkeypatch, no_settings):
    monkeypatch.setenv("CCP4I2_JOB_TARGET", "Queue")
    assert dispatch.job_target_name() == "queue"
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "local", raising=False)
    assert dispatch.job_target_name() == "local"  # settings win over env


def test_an_unknown_name_says_what_is_registered(no_settings):
    with pytest.raises(dispatch.UnknownRunTarget, match="no run target named 'azure'; registered: local"):
        dispatch.get_target("azure")


def test_a_bad_registration_is_reported_not_hidden(monkeypatch, no_settings):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {"local": dispatch.DEFAULT_RUN_TARGETS["local"],
                         "gone": "no_such_package.dispatch.Target",
                         "bare": "NotDotted",
                         "broken": _path(Unconstructable)}, raising=False)
    with pytest.raises(dispatch.RunTargetError, match="run target 'gone' could not be loaded"):
        dispatch.get_target("gone")
    with pytest.raises(dispatch.RunTargetError, match="not a dotted"):
        dispatch.get_target("bare")
    with pytest.raises(dispatch.RunTargetError, match="could not be constructed.*connection string"):
        dispatch.get_target("broken")
    reported = {t["name"]: t for t in dispatch.available_targets()}
    assert reported["local"]["runs_jobs"] is True
    assert "could not be loaded" in reported["gone"]["error"]
    assert reported["broken"]["runs_jobs"] is False


def test_run_job_context_aware_uses_the_deployment_target(monkeypatch, no_settings, fake_job):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS, "queue": _path(FakeQueueTarget)}, raising=False)
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "queue", raising=False)
    FakeQueueTarget.calls.clear()
    result = context_run.run_job_context_aware(fake_job, force_dispatch=True, synchronous=True)
    assert result == {"success": True, "data": fake_job, "status": 200}
    assert FakeQueueTarget.calls == [(fake_job, True)]


def test_force_local_ignores_the_deployment_target(monkeypatch, no_settings, fake_job):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS, "queue": _path(FakeQueueTarget)}, raising=False)
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "queue", raising=False)
    seen = []
    monkeypatch.setattr(local, "run_job_local",
                        lambda job, synchronous=False: seen.append((job, synchronous)) or
                        {"success": True, "data": job, "status": 200})
    FakeQueueTarget.calls.clear()
    result = context_run.run_job_context_aware(fake_job, force_local=True, force_dispatch=True)
    assert result["success"] is True
    assert seen == [(fake_job, False)]
    assert FakeQueueTarget.calls == []


def test_a_misconfigured_target_is_a_result_not_an_exception(monkeypatch, no_settings, fake_job):
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "azure", raising=False)
    result = context_run.run_job_context_aware(fake_job, force_dispatch=True)
    assert result["success"] is False and result["status"] == 500
    assert "registered: local" in result["error"]
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS, "batch": _path(FakeBatchTarget)}, raising=False)
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "batch", raising=False)
    result = context_run.run_job_context_aware(fake_job, force_dispatch=True)
    assert result["success"] is False and "does not run jobs" in result["error"]


def test_program_checks_are_authoritative_only_on_the_local_target(monkeypatch, no_settings):
    monkeypatch.setattr(context_run, "ccp4_available", lambda: True)
    assert context_run.program_checks_are_authoritative() is True
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "queue", raising=False)
    assert context_run.program_checks_are_authoritative() is False
