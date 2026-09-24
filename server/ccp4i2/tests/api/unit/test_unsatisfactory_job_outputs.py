"""A job that ran to the end keeps its verdict, and publishes what it made.

aimless_pipe whose ctruncate step crashes reports UNSATISFACTORY: it scaled
the data and it could not convert the intensities. Two faults conspired to
throw all of that away.

The runner stamped FINISHED over the top of the plugin's verdict ("belt and
braces" for legacy pipelines that never set a status at all), so the failure
was invisible. And gleaning was gated on success, so the scaled unmerged
data, the statistics and the pipeline XML -- all declared in outputData, all
sitting on disk -- were never registered. What the user saw was a completed
data-reduction pipeline with no outputs whatsoever.

INTERRUPTED is the deliberate exception: an interrupted job publishes
nothing, so Resume starts from a clean slate (docs/interrupt-and-resume.md).
"""

import pytest
from asgiref.sync import async_to_sync

from ccp4i2.core.CCP4File import CDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.db import models
from ccp4i2.lib.async_run_job import run_job_async


class _StubPlugin(CPluginScript):
    """Writes one output file, then reports whatever it was told to report.

    Stands in for a pipeline whose late step died after its early steps had
    already produced real files.
    """

    TASKNAME = "unsatisfactory_stub"

    #: What process() reports. Set per test before the runner builds it.
    verdict = CPluginScript.UNSATISFACTORY

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.container.outputData.addContent(CDataFile, "SALVAGED")

    def process(self):
        salvaged = self.container.outputData.SALVAGED
        salvaged.annotation = "what the early steps managed"
        path = salvaged.fullPath.get()
        with open(path, "w") as handle:
            handle.write("the data this job did produce\n")

        if self.verdict is None:
            # A legacy pipeline that never reports anything at all.
            return CPluginScript.SUCCEEDED
        self.reportStatus(self.verdict)
        return self.verdict


@pytest.fixture
def project(tmp_path):
    directory = tmp_path / "unsatisfactory"
    directory.mkdir()
    return models.Project.objects.create(
        name="unsatisfactory", directory=str(directory)
    )


@pytest.fixture
def run_stub(project, monkeypatch):
    """Run one stub job reporting `verdict`, and hand back the Job row."""

    def _run(verdict):
        job = models.Job.objects.create(
            project=project,
            number="1",
            title="stub",
            task_name=_StubPlugin.TASKNAME,
            status=models.Job.Status.PENDING,
        )
        job.directory.mkdir(parents=True, exist_ok=True)

        plugin_class = type("_Stub", (_StubPlugin,), {"verdict": verdict})
        monkeypatch.setattr(
            "ccp4i2.lib.async_run_job.get_plugin_class", lambda _name: plugin_class
        )
        async_to_sync(run_job_async)(job.uuid)
        job.refresh_from_db()
        return job

    return _run


def _outputs(job):
    return models.File.objects.filter(
        job=job, directory=models.File.Directory.JOB_DIR
    )


def test_unsatisfactory_verdict_survives_the_runner(run_stub):
    """The runner must not relabel a verdict the plugin gave."""
    job = run_stub(CPluginScript.UNSATISFACTORY)
    assert job.status == models.Job.Status.UNSATISFACTORY, job.get_status_display()


def test_unsatisfactory_job_publishes_the_outputs_it_did_make(run_stub):
    """The whole point: partial outputs reach the project."""
    job = run_stub(CPluginScript.UNSATISFACTORY)
    assert [f.job_param_name for f in _outputs(job)] == ["SALVAGED"]
    assert _outputs(job)[0].annotation == "what the early steps managed"


def test_an_unsatisfactory_job_has_stopped_running(run_stub):
    """A terminal status carries a finish time, or the job reads as ongoing."""
    job = run_stub(CPluginScript.UNSATISFACTORY)
    assert job.finish_time is not None


def test_a_plugin_that_reports_nothing_still_finishes(run_stub):
    """The belt-and-braces net the conditional stamp had to keep: a legacy
    pipeline that never sets a status must not be left at RUNNING."""
    job = run_stub(None)
    assert job.status == models.Job.Status.FINISHED, job.get_status_display()
    assert [f.job_param_name for f in _outputs(job)] == ["SALVAGED"]


def test_a_succeeding_plugin_is_unaffected(run_stub):
    job = run_stub(CPluginScript.SUCCEEDED)
    assert job.status == models.Job.Status.FINISHED, job.get_status_display()
    assert [f.job_param_name for f in _outputs(job)] == ["SALVAGED"]


def test_an_interrupted_job_keeps_its_status_and_publishes_nothing(run_stub):
    """Resume needs a clean slate, so INTERRUPTED stays out of the glean."""
    job = run_stub(CPluginScript.INTERRUPTED)
    assert job.status == models.Job.Status.INTERRUPTED, job.get_status_display()
    assert list(_outputs(job)) == []
