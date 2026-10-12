"""Cancelling a local job kills what it started.

Local dispatch runs a job in run_job's own process, which recorded no PID:
cancel marked the job Interrupted and killed nothing, and a MrBUMP spun on
at full CPU for two hours (2026-10-12). run_job now records its own PID.
"""
import os
import subprocess
import sys
import time

import psutil
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.db.management.commands.run_job import record_own_process

API = "/api/ccp4i2"


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def job(tmp_path):
    directory = tmp_path / "cancelkills"
    (directory / "CCP4_JOBS" / "job_1").mkdir(parents=True)
    project = models.Project.objects.create(name="cancelkills", directory=str(directory))
    return models.Job.objects.create(project=project, number="1", task_name="mrbump_basic",
                                     status=models.Job.Status.RUNNING)


def test_run_job_records_its_own_process(job):
    record_own_process(job)
    job.refresh_from_db()
    assert job.process_id == os.getpid()


@pytest.mark.skipif(sys.platform == "win32", reason="a POSIX process tree")
def test_cancel_kills_the_job_and_what_it_started(client, job):
    # A stand-in for run_job and the program it launched: a shell with a child
    runner = subprocess.Popen(["/bin/sh", "-c", "sleep 300 & wait"], start_new_session=True)
    deadline = time.time() + 5
    while not psutil.Process(runner.pid).children() and time.time() < deadline:
        time.sleep(0.05)
    program = psutil.Process(runner.pid).children()[0]
    job.process_id = runner.pid
    job.save()

    r = client.post(f"{API}/jobs/{job.id}/cancel/")
    assert r.status_code == 200, r.content
    runner.wait(timeout=5)
    gone, alive = psutil.wait_procs([program], timeout=5)
    assert not alive
    job.refresh_from_db()
    assert job.status == models.Job.Status.INTERRUPTED
