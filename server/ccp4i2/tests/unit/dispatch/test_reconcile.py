"""
The dispatch record and the reconcile, with a fake program target.

The reconcile is idempotent and caller-agnostic: it reads dispatch.json, asks
the target, and acts only on a terminal answer, by starting the job again
through whatever runner it is given. CCP4-free; the runner is a stub and the
job is a stub with a directory and a status.
"""
import types

import pytest
from django.conf import settings

from ccp4i2.db import models
from ccp4i2.lib import dispatch
from ccp4i2.lib.utils.jobs import dispatch_record as dr


class FakeBatch:
    """A program target whose answers a test sets."""
    state = "running"
    log = None
    polled = []

    def submit(self, tree, argv, out_dir, sizing_hint): return "batch-42"
    def poll(self, handle):
        FakeBatch.polled.append(handle)
        return FakeBatch.state
    def cancel(self, handle): return None
    def logs(self, handle): return FakeBatch.log


@pytest.fixture
def registered(monkeypatch):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS, "batch": f"{__name__}.FakeBatch"}, raising=False)
    FakeBatch.state, FakeBatch.log, FakeBatch.polled = "running", None, []


class FakeJob:
    def __init__(self, directory, status):
        self.directory = str(directory)
        self.status = status
        self.saved = []
    def save(self, update_fields=None):
        self.saved.append((self.status, tuple(update_fields or ())))
    def get_status_display(self):
        return models.Job.Status(self.status).label


def test_the_record_round_trips_and_marks_a_harvest(tmp_path):
    assert dr.read_record(tmp_path) is None and not dr.is_harvest(tmp_path)
    rec = dr.new_record(tmp_path, target="batch", handle="batch-42", out_dir="/x")
    assert rec["state"] == "submitted" and rec["submitted_at"] and rec["handle"] == "batch-42"
    assert not dr.is_harvest(tmp_path)
    dr.write_record(tmp_path, state="succeeded")
    assert dr.is_harvest(tmp_path)
    assert dr.read_record(tmp_path)["target"] == "batch"      # merged, not replaced


def test_a_job_that_never_dispatched_or_is_not_remote_is_left_alone(tmp_path, registered):
    runs = []
    job = FakeJob(tmp_path, models.Job.Status.RUNNING_REMOTELY)
    assert dr.reconcile(job, run=runs.append)["action"] == "none"
    dr.new_record(tmp_path, target="batch", handle="batch-42")
    job.status = models.Job.Status.FINISHED
    out = dr.reconcile(job, run=runs.append)
    assert out["action"] == "none" and "finished" in out["reason"]
    assert runs == [] and FakeBatch.polled == []


def test_a_run_still_going_changes_nothing_but_the_record(tmp_path, registered):
    dr.new_record(tmp_path, target="batch", handle="batch-42")
    job = FakeJob(tmp_path, models.Job.Status.RUNNING_REMOTELY)
    runs = []
    for state in ("queued", "running", "unknown", "Running", "weird"):
        FakeBatch.state = state
        out = dr.reconcile(job, run=runs.append)
        assert out["action"] == "none", out
    assert dr.read_record(tmp_path)["state"] == "unknown"     # "weird" is not a state we know
    assert dr.read_record(tmp_path)["polled_at"]
    assert runs == [] and job.status == models.Job.Status.RUNNING_REMOTELY


@pytest.mark.parametrize("state", ["succeeded", "failed", "cancelled"])
def test_a_terminal_run_starts_the_harvest_once(tmp_path, registered, state):
    dr.new_record(tmp_path, target="batch", handle="batch-42")
    job = FakeJob(tmp_path, models.Job.Status.RUNNING_REMOTELY)
    FakeBatch.state, FakeBatch.log = state, str(tmp_path / "stderr.txt")
    started = []
    def run(j):
        started.append(j.status)
        return {"success": True, "data": j, "status": 200}
    out = dr.reconcile(job, run=run)
    assert out["action"] == "harvest_started" and out["state"] == state
    rec = dr.read_record(tmp_path)
    assert rec["state"] == state and rec["stderr"] == str(tmp_path / "stderr.txt") and rec["harvest_started_at"]
    assert started == [models.Job.Status.PENDING]              # started again, from pending
    # Idempotent: a second call does not start it twice.
    job.status = models.Job.Status.RUNNING_REMOTELY
    out = dr.reconcile(job, run=run)
    assert out["action"] == "none" and "already started" in out["reason"]
    assert len(started) == 1


def test_a_harvest_that_cannot_start_is_reported_and_the_job_stays_remote(tmp_path, registered):
    dr.new_record(tmp_path, target="batch", handle="batch-42")
    job = FakeJob(tmp_path, models.Job.Status.RUNNING_REMOTELY)
    FakeBatch.state = "succeeded"
    out = dr.reconcile(job, run=lambda j: {"success": False, "error": "no interpreter", "status": 500})
    assert out["action"] == "harvest_failed" and "no interpreter" in out["reason"]
    assert job.status == models.Job.Status.RUNNING_REMOTELY
    rec = dr.read_record(tmp_path)
    assert rec["harvest_started_at"] is None and rec["harvest_error"] == "no interpreter"


def test_a_target_the_deployment_no_longer_registers_is_an_error_result(tmp_path, registered, monkeypatch):
    dr.new_record(tmp_path, target="gone", handle="h")
    job = FakeJob(tmp_path, models.Job.Status.RUNNING_REMOTELY)
    out = dr.reconcile(job, run=lambda j: None)
    assert out["action"] == "error" and "no run target named 'gone'" in out["reason"]
    # A job-only target cannot be polled for a program.
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "local", raising=False)
    dr.write_record(tmp_path, target="local")
    out = dr.reconcile(job, run=lambda j: None)
    assert out["action"] == "error" and "does not run programs" in out["reason"]
