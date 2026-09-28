"""The backfill_complete_mtz management command: registers final.mtz as the
COMPLETE_MTZ output of finished dimple jobs that ran before i2Dimple tracked
it, against the dimple job itself, dry run by default."""
import pytest
from django.core.management import call_command
from django.core.management.base import CommandError

from ccp4i2.db import models

MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER


def _project(tmp_path, name):
    directory = tmp_path / name
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name=name, directory=str(directory))


def _dimple(project, number="1.3", status=models.Job.Status.FINISHED, on_disk=True):
    parent = models.Job.objects.create(project=project, number=number.split(".")[0],
                                       task_name="SubstituteLigand", title="s", status=status)
    job = models.Job.objects.create(project=project, number=number, task_name="i2Dimple",
                                    title="d", status=status, parent=parent)
    job.directory.mkdir(parents=True)
    (job.directory / "final.pdb").write_text("END\n")
    if on_disk:
        (job.directory / "final.mtz").write_bytes(b"MTZ ")
    return job


def _rows(job):
    return models.File.objects.filter(job=job, job_param_name="COMPLETE_MTZ")


@pytest.fixture
def legacy(tmp_path):
    """Three members: one legacy (file on disk, unregistered), one whose
    dimple outputs are gone, one already registered by the modern pipeline;
    plus a project outside the campaign."""
    group = models.ProjectGroup.objects.create(name="camp", type="fragment_set")
    x1, x2, x3 = (_project(tmp_path, n) for n in ("x1", "x2", "x3"))
    for project in (x1, x2, x3):
        models.ProjectGroupMembership.objects.create(group=group, project=project, type=MEMBER)
    unregistered = _dimple(x1)
    gone = _dimple(x2, on_disk=False)
    modern = _dimple(x3)
    file_type, _ = models.FileType.objects.get_or_create(
        name="application/CCP4-mtz", defaults={"description": "MTZ"})
    models.File.objects.create(name="final.mtz", directory=models.File.Directory.JOB_DIR,
                               type=file_type, job=modern, job_param_name="COMPLETE_MTZ")
    outside = _dimple(_project(tmp_path, "elsewhere"))
    unfinished = _dimple(_project(tmp_path, "running"), status=models.Job.Status.RUNNING)
    return {"unregistered": unregistered, "gone": gone, "modern": modern,
            "outside": outside, "unfinished": unfinished, "parent": unregistered.parent}


def test_dry_run_reports_and_writes_nothing(legacy, capsys):
    call_command("backfill_complete_mtz")
    out = capsys.readouterr().out
    assert "DRY RUN" in out
    assert "Would register 2 COMPLETE_MTZ row(s)" in out
    assert "Skipped (already registered): 1" in out
    assert "Skipped (no final.mtz on disk): 1" in out
    assert not _rows(legacy["unregistered"]).exists()


def test_commit_registers_the_row_the_gleaner_would_have(legacy):
    call_command("backfill_complete_mtz", "--commit")
    row = _rows(legacy["unregistered"]).get()
    assert row.name == "final.mtz"
    assert row.directory == models.File.Directory.JOB_DIR
    assert row.type.name == "application/CCP4-mtz"
    assert (row.sub_type, row.content) == (0, 0)
    assert row.annotation == "Complete unsplit reflection file from dimple"
    assert row.path == legacy["unregistered"].directory / "final.mtz"
    # against the dimple job, never the pipeline that ran it
    assert not models.File.objects.filter(job=legacy["parent"]).exists()
    assert _rows(legacy["outside"]).count() == 1
    assert not _rows(legacy["gone"]).exists()
    assert not _rows(legacy["unfinished"]).exists()
    assert _rows(legacy["modern"]).count() == 1


def test_a_second_run_adds_nothing(legacy, capsys):
    call_command("backfill_complete_mtz", "--commit")
    call_command("backfill_complete_mtz", "--commit")
    assert "Registered 0 COMPLETE_MTZ row(s)" in capsys.readouterr().out
    assert models.File.objects.filter(job_param_name="COMPLETE_MTZ").count() == 3


def test_campaign_scopes_to_its_members(legacy):
    call_command("backfill_complete_mtz", "--commit", "--campaign", "camp")
    assert _rows(legacy["unregistered"]).exists()
    assert not _rows(legacy["outside"]).exists()


def test_limit_stops_after_a_batch(legacy):
    call_command("backfill_complete_mtz", "--commit", "--limit", "1")
    assert models.File.objects.filter(job_param_name="COMPLETE_MTZ").count() == 2


def test_unknown_scope_is_an_error(legacy):
    with pytest.raises(CommandError):
        call_command("backfill_complete_mtz", "--campaign", "nope")
    with pytest.raises(CommandError):
        call_command("backfill_complete_mtz", "--projectname", "nope")
