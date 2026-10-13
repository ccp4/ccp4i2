"""How "fill from campaign" scales with the size of the campaign.

An XChem-scale campaign has 500-1000+ members, and the fill is one HTTP
request (the object_method endpoint), so its cost matters twice: as wall
time on the desktop, and as the number of database queries, which is what a
hosted deployment pays a network round trip for each of.

Opt-in, because seeding a thousand members takes a while and the numbers are
for reading rather than asserting:

    CCP4I2_BENCH=1 ccp4-python -m pytest ccp4i2/tests/db/test_pandda_fanin_scaling.py -s

Each size reports, for the candidate preview, the fill, and the fill reached
the way the interface reaches it (the plugin method through object_method):
seconds, queries, and queries per member.
"""
import os
import time
from pathlib import Path
from shutil import rmtree

import pytest
from django.db import connection
from django.test import TestCase, override_settings

from ...db import models
from ...lib.utils.helpers.object_method import object_method
from ...lib.utils.jobs.create import create_job
from ...lib.utils.jobs.pandda_fanin import campaign_candidates, fill_datasets_from_campaign

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not os.environ.get("CCP4I2_BENCH"),
                       reason="benchmark: set CCP4I2_BENCH=1 to run"),
]

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_FANIN_SCALING_DIR"
PARENT = models.ProjectGroupMembership.MembershipType.PARENT
MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER
FINISHED = models.Job.Status.FINISHED
SIZES = [int(n) for n in os.environ.get("CCP4I2_BENCH_SIZES", "10,100,1000").split(",")]


def seed_campaign(size):
    """A parent and ``size`` members, each with a finished dimple (model and
    reflections registered) and a finished acedrg (dictionary registered):
    the shape of a campaign ready for PanDDA. Bulk-created, so seeding a
    thousand does not dominate the run."""
    tag = f"c{size}"
    parent = models.Project.objects.create(name=f"{tag}-ref", directory=str(PROJECTS_DIR / f"{tag}-ref"))
    members = models.Project.objects.bulk_create([
        models.Project(name=f"{tag}-x{i:04d}", directory=str(PROJECTS_DIR / f"{tag}-x{i:04d}"))
        for i in range(size)])
    group = models.ProjectGroup.objects.create(name=tag, type="fragment_set")
    models.ProjectGroupMembership.objects.bulk_create(
        [models.ProjectGroupMembership(group=group, project=parent, type=PARENT)]
        + [models.ProjectGroupMembership(group=group, project=m, type=MEMBER) for m in members])

    jobs = models.Job.objects.bulk_create(
        [models.Job(project=m, number="1", task_name="i2Dimple", title="dimple", status=FINISHED)
         for m in members]
        + [models.Job(project=m, number="2", task_name="LidiaAcedrgNew", title="acedrg", status=FINISHED)
           for m in members])
    types = {mime: models.FileType.objects.get_or_create(name=mime, defaults={"description": mime})[0]
             for mime in ("chemical/x-pdb", "application/CCP4-mtz", "application/refmac-dictionary")}
    files = []
    for job in models.Job.objects.filter(project__in=members).select_related("project"):
        job.directory.mkdir(parents=True, exist_ok=True)
        outputs = ((("XYZOUT", "final.pdb", "chemical/x-pdb"),
                    ("COMPLETE_MTZ", "final.mtz", "application/CCP4-mtz"))
                   if job.task_name == "i2Dimple"
                   else (("DICTOUT_LIST[0]", "LIG.cif", "application/refmac-dictionary"),))
        for param, name, mime in outputs:
            (job.directory / name).write_text("x\n")
            files.append(models.File(name=name, directory=models.File.Directory.JOB_DIR,
                                     type=types[mime], job=job, job_param_name=param))
    models.File.objects.bulk_create(files)
    assert len(jobs) == 2 * size

    (Path(parent.directory) / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return parent


def pending_job(parent):
    return models.Job.objects.get(uuid=create_job(
        projectId=str(parent.uuid), taskName="pandda_campaign", title="PanDDA"))


def measure(fn):
    """Seconds and queries for one call. Counted with an execute wrapper
    rather than CaptureQueriesContext, whose log stops at 9000 entries."""
    count = [0]

    def counter(execute, sql, params, many, context):
        count[0] += 1
        return execute(sql, params, many, context)

    with connection.execute_wrapper(counter):
        start = time.perf_counter()
        result = fn()
        seconds = time.perf_counter() - start
    return result, seconds, count[0]


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class FanInScalingTest(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)

    def test_scaling(self):
        rows = []
        for size in SIZES:
            parent = seed_campaign(size)

            job = pending_job(parent)
            preview, t, q = measure(lambda: campaign_candidates(job))
            self.assertEqual(len(preview["candidates"]), size, preview["reason"])
            rows.append((size, "candidates", t, q))

            filled, t, q = measure(lambda: fill_datasets_from_campaign(job))
            self.assertTrue(filled.success, filled.error)
            self.assertEqual(len(filled.data["added"]), size)
            rows.append((size, "fill", t, q))

            # The interface's route: a fresh job, the plugin loaded from its
            # params, then the method. Preview then fill, as the page does.
            job = pending_job(parent)
            _, t1, q1 = measure(lambda: object_method(job, "pandda_campaign", "campaignCandidates"))
            filled, t2, q2 = measure(lambda: object_method(job, "pandda_campaign", "fillDatasetsFromCampaign"))
            self.assertTrue(filled["success"], filled)
            rows.append((size, "object_method preview", t1, q1))
            rows.append((size, "object_method fill", t2, q2))

        print()
        print(f"{'members':>8}  {'step':<22} {'seconds':>8} {'queries':>8} {'q/member':>9}")
        for size, step, seconds, queries in rows:
            print(f"{size:>8}  {step:<22} {seconds:>8.2f} {queries:>8} {queries / size:>9.1f}")
