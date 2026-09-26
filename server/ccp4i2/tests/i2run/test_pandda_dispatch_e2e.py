"""
Dispatch end to end: submit, reconcile from the command, harvest in a job process.

Drives manage.py in child processes against a file database (the tier's own
database is in memory, invisible to a child; see test_run_target_local.py).
The program target is a fake registered through CCP4I2_RUN_TARGETS in the
environment; the test steers "the run" by leaving a tree where the argv said
and telling the fake what to answer.
"""
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from ccp4i2.tests.unit.pandda.conftest import LABELS, source_files

pytest.importorskip("yaml", reason="needs PyYAML")

SERVER_DIR = Path(__file__).resolve().parents[3]
FAKE = "ccp4i2.tests.i2run.fake_targets.FakeBatch"
LOCAL = "ccp4i2.lib.dispatch.local.LocalTarget"

I2RUN = """
import django, sys; django.setup()
from django.core.management import call_command
sys.argv = ["manage.py", "i2run"] + sys.argv[1:]
call_command("i2run")
"""


def _env(tmp_path, **extra):
    env = os.environ.copy()
    for k in ("CCP4I2_JOB_TARGET", "EXECUTION_MODE", "SERVICE_BUS_CONNECTION_STRING"):
        env.pop(k, None)
    env.update({
        "DJANGO_SETTINGS_MODULE": "ccp4i2.config.test_settings",
        "CCP4I2_DB_FILE": str(tmp_path / "db.sqlite"),
        "CCP4I2_PROJECTS_DIR": str(tmp_path / "projects"),
        "CCP4I2_HOME": str(tmp_path / "home"),
        "CCP4I2_RUN_TARGETS": json.dumps({"local": LOCAL, "batch": FAKE}),
        "CCP4I2_TEST_FAKE_BATCH_RECORD": str(tmp_path / "submission.json"),
    })
    env.update(extra)
    return env


def _last_json(text):
    """The last line of a child's stdout that is JSON (settings banners and
    loggers print around it)."""
    for line in reversed(text.strip().splitlines()):
        try:
            return json.loads(line)
        except ValueError:
            continue
    raise AssertionError(f"no JSON line in:\n{text[-1500:]}")


def _run(env, *argv, timeout=600, **overrides):
    """Run a child with ``env`` plus ``overrides`` (how a test steers the fake target)."""
    env = {**env, **overrides}
    r = subprocess.run([sys.executable, *argv], cwd=SERVER_DIR, env=env,
                       capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, (argv[:3], r.stdout[-1500:], r.stderr[-2500:])
    return r


def _dataset_args(_tmp_path):
    args = []
    for label in LABELS:
        pdb, mtz, cif = source_files(label)
        args += ["--DATASETS", f"DTAG={label}", f"XYZIN={pdb}", f"HKLIN={mtz}", f"DICT={cif}"]
    return args


def _job(db_path):
    db = sqlite3.connect(db_path)
    row = db.execute("select j.uuid, j.status, p.directory, j.number from ccp4i2_job j "
                     "join ccp4i2_project p on p.id = j.project_id where j.task_name = 'pandda_campaign'").fetchone()
    db.close()
    return row


@pytest.fixture
def project_name():
    """Unique per run: the child's settings put projects under a fixed home
    (test_settings ignores CCP4I2_PROJECTS_DIR), so a reused name would find
    the last run's staging tree. Whatever the run left is removed after."""
    import shutil, uuid
    name = f"disp_{uuid.uuid4().hex[:8]}"
    yield name
    for root in (Path.home() / ".ccp4i2_test" / "test_projects",):
        shutil.rmtree(root / name, ignore_errors=True)


def _setup(tmp_path):
    env = _env(tmp_path)
    (tmp_path / "projects").mkdir()
    (tmp_path / "home").mkdir()
    _run(env, "manage.py", "migrate", "--run-syncdb", "-v", "0")
    return env


def _submit(tmp_path, env, project_name):
    _run(env, "-c", I2RUN, "pandda_campaign", "--project_name", project_name, "--RUN_MODE", "dispatch",
         "--DISPATCH_TARGET", "batch", *_dataset_args(tmp_path))
    uuid, status, project_dir, number = _job(tmp_path / "db.sqlite")
    job_dir = Path(project_dir) / "CCP4_JOBS" / f"job_{number}"
    diagnostic = (job_dir / "diagnostic.xml").read_text()[-3000:] if (job_dir / "diagnostic.xml").is_file() else "(no diagnostic.xml)"
    assert status == 7, f"expected RUNNING_REMOTELY (7), got {status}\n{diagnostic}"     # waits remotely
    submission = json.loads((tmp_path / "submission.json").read_text())
    record = json.loads((job_dir / "dispatch.json").read_text())
    assert record["target"] == "batch" and record["state"] == "submitted"
    return uuid, job_dir, submission


def test_submit_reconcile_harvest_succeeded(tmp_path, project_name):
    from ccp4i2.tests.unit.pandda.synthetic_tree import event_record, make_tree
    env = _setup(tmp_path)
    uuid, job_dir, submission = _submit(tmp_path, env, project_name)

    # Still queued: the command reports and changes nothing.
    r = _run(env, "manage.py", "reconcile_dispatch", "--job", uuid, "--json")
    assert _last_json(r.stdout)["action"] == "none"
    assert _job(tmp_path / "db.sqlite")[1] == 7

    # "The run" finishes: a complete tree where the argv said.
    make_tree(Path(submission["out_dir"]),
              {f"xtal-{i:04d}": ([event_record(1)] if i == 0 else []) for i in range(len(LABELS))})
    r = _run(env, "manage.py", "reconcile_dispatch", "--all", "--json",
             CCP4I2_TEST_FAKE_BATCH_STATE="succeeded")
    out = _last_json(r.stdout)
    assert out["action"] == "harvest_started", out

    # The harvest was a real job process on the local target. Wait for it.
    import time
    deadline = time.time() + 300
    while time.time() < deadline and _job(tmp_path / "db.sqlite")[1] in (1, 2, 3, 7):
        time.sleep(2)
    status = _job(tmp_path / "db.sqlite")[1]
    assert status in (6, 10), f"harvest ended with status {status}"       # FINISHED / UNSATISFACTORY
    program = (job_dir / "program.xml").read_text()
    assert "<state>finished</state>" in program or "<state>partial</state>" in program or "<state>empty</state>" in program
    assert "<dispatch>" in program and "<state>succeeded</state>" in program.split("<dispatch>")[1]
    db = sqlite3.connect(tmp_path / "db.sqlite")
    n_out = db.execute("select count(*) from ccp4i2_file f join ccp4i2_job j on j.id = f.job_id "
                       "where j.uuid = ? and f.directory = 1", (uuid,)).fetchone()[0]
    db.close()
    assert n_out > 0, "a harvested run publishes its outputs"


def test_submit_reconcile_harvest_failed(tmp_path, project_name):
    env = _setup(tmp_path)
    uuid, job_dir, _submission = _submit(tmp_path, env, project_name)
    stderr = job_dir / "remote_stderr.txt"
    stderr.write_text("ray.exceptions.OutOfMemoryError: Task was killed due to the node running low on memory\n")
    r = _run(env, "manage.py", "reconcile_dispatch", "--job", uuid, "--json",
             CCP4I2_TEST_FAKE_BATCH_STATE="failed", CCP4I2_TEST_FAKE_BATCH_LOG=str(stderr))
    assert _last_json(r.stdout)["action"] == "harvest_started"
    import time
    deadline = time.time() + 300
    while time.time() < deadline and _job(tmp_path / "db.sqlite")[1] in (1, 2, 3, 7):
        time.sleep(2)
    assert _job(tmp_path / "db.sqlite")[1] == 5, "FAILED"
    program = (job_dir / "program.xml").read_text()
    assert "<state>failed</state>" in program
    failure = program.split("<failure>")[1].split("</failure>")[0]
    assert failure not in ("", "unclassified_crash"), failure      # classified by the catalogue
