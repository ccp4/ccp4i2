"""
A job started through the run-target registry really runs, as a subprocess.

The other i2run tests share pytest-django's in-memory database, which a
child process can never see, so none of them can prove that
run_job_context_aware -> LocalTarget -> `ccp4-python -m django run_job`
reaches a job and finishes it. This one drives manage.py in subprocesses
against its own file database, with no run-target settings at all: the
desktop's situation. The fan-out task is the vehicle because it starts each
receipt through run_job_context_aware(synchronous=True).
"""
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("yaml", reason="needs PyYAML")

SERVER_DIR = Path(__file__).resolve().parents[3]          # .../server
SYNTHETIC = "ccp4i2.tests.unit.pandda.synthetic_tree"


def _env(tmp_path):
    env = os.environ.copy()
    for k in ("CCP4I2_JOB_TARGET", "EXECUTION_MODE", "SERVICE_BUS_CONNECTION_STRING"):
        env.pop(k, None)
    env.update({
        "DJANGO_SETTINGS_MODULE": "ccp4i2.config.test_settings",
        "CCP4I2_DB_FILE": str(tmp_path / "db.sqlite"),
        "CCP4I2_PROJECTS_DIR": str(tmp_path / "projects"),
        "CCP4I2_HOME": str(tmp_path / "home"),
    })
    return env


def _manage(env, *args, **kw):
    return subprocess.run([sys.executable, "manage.py", *args], cwd=SERVER_DIR, env=env,
                          capture_output=True, text=True, timeout=600, **kw)


@pytest.fixture
def tag():
    """Unique per run, and whatever the run left under the fixed test home
    (test_settings ignores CCP4I2_PROJECTS_DIR) is removed after."""
    import shutil, uuid
    t = uuid.uuid4().hex[:8]
    yield t
    root = Path.home() / ".ccp4i2_test" / "test_projects"
    for d in root.glob(f"rt_*{t}*"):
        shutil.rmtree(d, ignore_errors=True)


def test_a_dispatched_job_runs_to_completion_on_the_local_target(tmp_path, tag):
    env = _env(tmp_path)
    (tmp_path / "projects").mkdir()
    (tmp_path / "home").mkdir()

    r = _manage(env, "migrate", "--run-syncdb", "-v", "0")
    assert r.returncode == 0, r.stderr[-2000:]

    # Two member projects and a two-dataset tree, written by a child process
    # so the rows are committed in the file the jobs will read.
    setup = f"""
import django, json, sys; django.setup()
from pathlib import Path
from django.conf import settings
from ccp4i2.db import models
from {SYNTHETIC} import event_record, make_tree
members = []
for name in ("rt_m1_" + sys.argv[2], "rt_m2_" + sys.argv[2]):
    d = Path(settings.CCP4I2_PROJECTS_DIR) / name
    (d / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    members.append(models.Project.objects.create(name=name, directory=str(d)))
tree = make_tree(Path(sys.argv[1]) / "pandda2_out", {{"xtal-0000": [event_record(1)], "xtal-0001": []}})
manifest = Path(sys.argv[1]) / "manifest.json"
manifest.write_text(json.dumps({{
    "manifest_version": 1, "created": "now", "datasets_dir": "datasets", "projects_csv": "Projects.csv",
    "provenance": {{}},
    "datasets": [
        {{"xtal": "xtal-0000", "label": members[0].name, "project_uuid": str(members[0].uuid), "files": {{}}}},
        {{"xtal": "xtal-0001", "label": members[1].name, "project_uuid": str(members[1].uuid), "files": {{}}}},
    ]}}))
print(tree); print(manifest)
"""
    r = subprocess.run([sys.executable, "-c", setup, str(tmp_path), tag], cwd=SERVER_DIR, env=env,
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    tree, manifest = r.stdout.strip().splitlines()[-2:]

    # i2run reads its task arguments from sys.argv (as tests/i2run/utils.i2run
    # does), so drive it the same way, in its own process.
    i2run = """
import django, sys; django.setup()
from django.core.management import call_command
sys.argv = ["manage.py", "i2run"] + sys.argv[1:]
call_command("i2run")
"""
    r = subprocess.run([sys.executable, "-c", i2run, "pandda_fanout", "--project_name", f"rt_parent_{tag}",
                        "--MANIFEST", manifest, "--PANDDA_OUT_DIR", tree, "--RUN_RECEIPTS", "True"],
                       cwd=SERVER_DIR, env=env, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, (r.stdout[-1500:], r.stderr[-2000:])

    db = sqlite3.connect(tmp_path / "db.sqlite")
    rows = db.execute(
        "select p.name, j.status from ccp4i2_job j join ccp4i2_project p on p.id = j.project_id "
        "where j.task_name = 'pandda_events' order by p.name").fetchall()
    db.close()
    # 6 = FINISHED. A receipt still QUEUED (2) means the local target never
    # started a process, or the process never found the job.
    assert rows == [(f"rt_m1_{tag}", 6), (f"rt_m2_{tag}", 6)], rows


def test_an_unknown_job_target_from_the_environment_is_reported_not_worked_around(tmp_path, tag):
    """CCP4I2_JOB_TARGET=nowhere: the receipt is created but not started, and
    the fan-out's record says so, naming what is registered."""
    env = _env(tmp_path)
    env["CCP4I2_JOB_TARGET"] = "nowhere"
    (tmp_path / "projects").mkdir()
    (tmp_path / "home").mkdir()
    assert _manage(env, "migrate", "--run-syncdb", "-v", "0").returncode == 0

    setup = f"""
import django, json, sys; django.setup()
from pathlib import Path
from django.conf import settings
from ccp4i2.db import models
from {SYNTHETIC} import event_record, make_tree
name = "rt_x_" + sys.argv[2]
d = Path(settings.CCP4I2_PROJECTS_DIR) / name; (d / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
m = models.Project.objects.create(name=name, directory=str(d))
tree = make_tree(Path(sys.argv[1]) / "pandda2_out", {{"xtal-0000": [event_record(1)]}})
manifest = Path(sys.argv[1]) / "manifest.json"
manifest.write_text(json.dumps({{"manifest_version": 1, "created": "now", "datasets_dir": "datasets",
    "projects_csv": "Projects.csv", "provenance": {{}},
    "datasets": [{{"xtal": "xtal-0000", "label": name, "project_uuid": str(m.uuid), "files": {{}}}}]}}))
print(tree); print(manifest)
"""
    r = subprocess.run([sys.executable, "-c", setup, str(tmp_path), tag], cwd=SERVER_DIR, env=env,
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    tree, manifest = r.stdout.strip().splitlines()[-2:]

    i2run = """
import django, sys; django.setup()
from django.core.management import call_command
sys.argv = ["manage.py", "i2run"] + sys.argv[1:]
call_command("i2run")
"""
    r = subprocess.run([sys.executable, "-c", i2run, "pandda_fanout", "--project_name", f"rt_xparent_{tag}",
                        "--MANIFEST", manifest, "--PANDDA_OUT_DIR", tree, "--RUN_RECEIPTS", "True"],
                       cwd=SERVER_DIR, env=env, capture_output=True, text=True, timeout=600)
    db = sqlite3.connect(tmp_path / "db.sqlite")
    receipt_status = db.execute(
        "select status from ccp4i2_job where task_name = 'pandda_events'").fetchone()[0]
    fanout_dir = db.execute(
        "select p.directory, j.number from ccp4i2_job j join ccp4i2_project p on p.id = j.project_id "
        "where j.task_name = 'pandda_fanout'").fetchone()
    db.close()
    assert receipt_status in (1, 2), receipt_status        # PENDING or QUEUED: never ran
    program_xml = Path(fanout_dir[0]) / "CCP4_JOBS" / f"job_{fanout_dir[1]}" / "program.xml"
    text = program_xml.read_text()
    assert "did not start" in text and "no run target named 'nowhere'" in text and "registered: local" in text, text[-800:]
