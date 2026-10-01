"""What every scenario needs: a scratch home, data fetched once, i2run, a clone.

A scenario builds the project a help page is illustrated from. Run it from
server/ with ccp4-python and CCP4I2_HOME set to a scratch directory.
"""
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

DJANGO = {"DJANGO_SETTINGS_MODULE": "ccp4i2.config.settings"}


def scratch_home() -> Path:
    home = os.environ.get("CCP4I2_HOME")
    if not home:
        sys.exit("Set CCP4I2_HOME to a scratch directory: a scenario "
                 "creates a project and must not touch a live database.")
    for live in (".ccp4i2", ".ccp4i2-django", ".ccp4i2x"):
        if Path(home).resolve() == (Path.home() / live).resolve():
            sys.exit(f"CCP4I2_HOME is the live home {home}; use a scratch one.")
    return Path(home)


def inputs_dir() -> Path:
    work = scratch_home() / "scenario_inputs"
    work.mkdir(parents=True, exist_ok=True)
    return work


def fetch(url: str, name: str = None) -> Path:
    """Download once into the scratch home's scenario_inputs."""
    path = inputs_dir() / (name or Path(url).name)
    if not path.exists():
        print(f"Fetching {url}")
        urllib.request.urlretrieve(url, path)
    return path


def i2run(project: str, *args: str):
    print("i2run", args[0])
    subprocess.run([sys.executable, "-m", "ccp4i2.cli.i2run", *args,
                    "--project_name", project],
                   check=True, env={**os.environ, **DJANGO})


def clone_last(project: str, task: str):
    """An unrun clone of the project's last top-level <task> job: the input
    figures show a job being set up, not one that has run. Top-level only: a
    pipeline's sub-job of the same task is newer, and is not what the user
    set up (the Phaser EP pipeline runs phaser_ep_auto_phil inside it)."""
    subprocess.run([sys.executable, "manage.py", "shell", "-c", (
        "from ccp4i2.db.models import Job\n"
        "from ccp4i2.lib.utils.jobs.clone import clone_job\n"
        f"job = Job.objects.filter(project__name='{project}', "
        f"task_name='{task}', parent__isnull=True).order_by('-id').first()\n"
        "clone_job(str(job.uuid))\n")],
        check=True, env={**os.environ, **DJANGO})



def output_file_id(project: str, task: str, param: str) -> str:
    """The database id of the named output of the project's last top-level
    <task> job that has one. For a file inside a list item (a Phaser
    ensemble's structure), where the fileOut= syntax does not reach: given
    as ".../dbFileId=<id>" it is that job's file, recorded as used, where a
    path would be imported again as a new file of no known origin."""
    out = subprocess.run([sys.executable, "manage.py", "shell", "-c", (
        "from ccp4i2.db.models import File\n"
        f"f = File.objects.filter(job__project__name='{project}', "
        f"job__task_name='{task}', job__parent__isnull=True, "
        f"job_param_name='{param}').order_by('-job__id').first()\n"
        "print('ID=' + str(f.uuid))\n")],
        check=True, env={**os.environ, **DJANGO}, capture_output=True, text=True).stdout
    return next(line[3:] for line in out.splitlines() if line.startswith("ID="))
