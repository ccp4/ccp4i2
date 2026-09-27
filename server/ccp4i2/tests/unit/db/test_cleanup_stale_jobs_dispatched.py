"""The stale sweep must not judge a dispatched job by its age.

RUNNING_REMOTELY means no worker holds the job: its program runs on a program
run target and the job waits for a reconcile. So a worker restart says nothing
about it, and neither does the clock. On 2026-09-27 a deployment restarted the
worker 2.7 hours into a healthy Azure Batch run; the sweep marked the job
FAILED while its task kept running, and nothing could harvest the output
afterwards, because the reconcile only acts on jobs in RUNNING_REMOTELY.

The one dispatched job this command can still judge is one with no dispatch
record: no handle to poll, no target to ask, so no reconcile will ever move it.
"""
from datetime import timedelta
from io import StringIO

import pytest
from django.core.management import call_command
from django.utils import timezone

from ccp4i2.db import models
from ccp4i2.lib.utils.jobs import dispatch_record


def _job(tmp_path, status, *, hours_old, with_record=False, number=1):
    """A job whose directory exists. Job.directory is derived from the
    project's directory and the job number, so the project is what places it."""
    project = models.Project.objects.create(
        name=f"p{number}", directory=str(tmp_path)
    )
    job = models.Job.objects.create(
        project=project,
        number=str(number),
        status=status,
        title=f"job {number}",
    )
    job.directory.mkdir(parents=True, exist_ok=True)
    if with_record:
        dispatch_record.new_record(job.directory, target="batch", handle="batch-42")
    # creation_time is auto_now_add, so it has to be pushed back afterwards.
    models.Job.objects.filter(pk=job.pk).update(
        creation_time=timezone.now() - timedelta(hours=hours_old)
    )
    job.refresh_from_db()
    return job


def _sweep():
    out = StringIO()
    call_command("cleanup_stale_jobs", "--hours", "2", stdout=out)
    return out.getvalue()


@pytest.mark.django_db
def test_an_old_dispatched_run_is_left_alone(tmp_path):
    """The regression. Six hours is normal for a PanDDA campaign."""
    job = _job(tmp_path, models.Job.Status.RUNNING_REMOTELY, hours_old=6,
               with_record=True)
    _sweep()
    job.refresh_from_db()
    assert job.status == models.Job.Status.RUNNING_REMOTELY


@pytest.mark.django_db
def test_a_dispatched_job_with_no_record_is_swept(tmp_path):
    """Nothing can ever reconcile it, so nothing would ever finish it."""
    job = _job(tmp_path, models.Job.Status.RUNNING_REMOTELY, hours_old=6,
               with_record=False)
    _sweep()
    job.refresh_from_db()
    assert job.status == models.Job.Status.FAILED


@pytest.mark.django_db
def test_a_recent_dispatched_job_with_no_record_is_left_alone(tmp_path):
    """It may be between the status change and the record being written."""
    job = _job(tmp_path, models.Job.Status.RUNNING_REMOTELY, hours_old=0.1,
               with_record=False)
    _sweep()
    job.refresh_from_db()
    assert job.status == models.Job.Status.RUNNING_REMOTELY


@pytest.mark.django_db
def test_an_old_running_job_is_still_swept(tmp_path):
    """The behaviour this command exists for is unchanged."""
    job = _job(tmp_path, models.Job.Status.RUNNING, hours_old=6)
    _sweep()
    job.refresh_from_db()
    assert job.status == models.Job.Status.FAILED


@pytest.mark.django_db
def test_an_unreadable_record_counts_as_present(tmp_path):
    """A partial write is a reason to ask the target, not to declare it dead."""
    job = _job(tmp_path, models.Job.Status.RUNNING_REMOTELY, hours_old=6)
    (job.directory / "dispatch.json").write_text("{ truncated")
    _sweep()
    job.refresh_from_db()
    assert job.status == models.Job.Status.RUNNING_REMOTELY
