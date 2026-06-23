"""
Unit tests for the environment-aware local-execution feasibility logic
(ccp4i2.lib.utils.jobs.context_run): ccp4_available / task_requires_ccp4 /
can_run_local. These are the server-side decision behind the run_local endpoint
— the client may *request* local execution, but the server decides whether it is
actually possible in this environment.

Pure functions, no database; run on a CCP4-free interpreter.
"""
import pytest

from ccp4i2.lib.utils.jobs import context_run


class TestTaskRequiresCCP4:
    def test_unflagged_task_requires_ccp4(self):
        # freerflag shells out to the freerflag binary -> not ccp4_free
        assert context_run.task_requires_ccp4("freerflag") is True

    def test_unknown_task_is_conservative(self):
        assert context_run.task_requires_ccp4("no_such_task_xyz") is True

    @pytest.mark.parametrize(
        "task", ["ProvideAsuContents", "coordinate_selector", "splitMtz", "mtzheader"]
    )
    def test_ccp4_free_tasks_do_not_require_ccp4(self, task):
        assert context_run.task_requires_ccp4(task) is False


class TestCanRunLocal:
    def test_no_ccp4_ccp4free_task_allowed(self, monkeypatch):
        monkeypatch.setattr(context_run, "ccp4_available", lambda: False)
        assert context_run.can_run_local("ProvideAsuContents") is True

    def test_no_ccp4_binary_task_refused(self, monkeypatch):
        monkeypatch.setattr(context_run, "ccp4_available", lambda: False)
        assert context_run.can_run_local("freerflag") is False

    def test_with_ccp4_anything_allowed(self, monkeypatch):
        monkeypatch.setattr(context_run, "ccp4_available", lambda: True)
        assert context_run.can_run_local("freerflag") is True
        assert context_run.can_run_local("ProvideAsuContents") is True


class TestCCP4Available:
    def test_false_without_ccp4_env(self, monkeypatch):
        context_run.ccp4_available.cache_clear()
        monkeypatch.delenv("CCP4", raising=False)
        try:
            assert context_run.ccp4_available() is False
        finally:
            context_run.ccp4_available.cache_clear()

    def test_false_when_ccp4_points_nowhere(self, monkeypatch):
        context_run.ccp4_available.cache_clear()
        monkeypatch.setenv("CCP4", "/nonexistent/ccp4/path/zzz")
        try:
            assert context_run.ccp4_available() is False
        finally:
            context_run.ccp4_available.cache_clear()
