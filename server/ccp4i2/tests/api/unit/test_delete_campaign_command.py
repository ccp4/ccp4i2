"""The delete_campaign management command."""

import pytest
from django.core.management import call_command
from django.core.management.base import CommandError

from ccp4i2.db import models

PARENT = models.ProjectGroupMembership.MembershipType.PARENT
MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER


def _project(tmp_path, name):
    directory = tmp_path / name
    (directory / "CCP4_JOBS").mkdir(parents=True)
    return models.Project.objects.create(name=name, directory=str(directory))


@pytest.fixture
def campaign(tmp_path):
    group = models.ProjectGroup.objects.create(
        name="camp", type=models.ProjectGroup.GroupType.FRAGMENT_SET
    )
    parent = _project(tmp_path, "camp_ref")
    member = _project(tmp_path, "camp_x1")
    models.ProjectGroupMembership.objects.create(group=group, project=parent, type=PARENT)
    models.ProjectGroupMembership.objects.create(group=group, project=member, type=MEMBER)
    site = models.CampaignSite.objects.create(
        group=group, name="pocket", origin_x=0, origin_y=0, origin_z=0, order=0
    )
    job = models.Job.objects.create(
        project=member, number="1", task_name="SubstituteLigand", title="t",
        status=models.Job.Status.FINISHED,
    )
    return {"group": group, "parent": parent, "member": member, "site": site, "job": job}


def test_deletes_group_sites_projects_and_directories(campaign, tmp_path):
    call_command("delete_campaign", "camp", "--yes")

    assert not models.ProjectGroup.objects.exists()
    assert not models.CampaignSite.objects.exists()
    assert not models.Project.objects.exists()
    assert not models.Job.objects.exists()
    assert not (tmp_path / "camp_ref").exists()
    assert not (tmp_path / "camp_x1").exists()


def test_dry_run_deletes_nothing(campaign, tmp_path, capsys):
    call_command("delete_campaign", "camp", "--dry-run")

    assert "nothing deleted" in capsys.readouterr().out
    assert models.ProjectGroup.objects.count() == 1
    assert models.Project.objects.count() == 2
    assert (tmp_path / "camp_x1").exists()


def test_keep_files_leaves_directories(campaign, tmp_path):
    call_command("delete_campaign", "camp", "--yes", "--keep-files")

    assert not models.Project.objects.exists()
    assert (tmp_path / "camp_ref" / "CCP4_JOBS").exists()
    assert (tmp_path / "camp_x1" / "CCP4_JOBS").exists()


def test_project_shared_with_another_group_is_kept(campaign, tmp_path, capsys):
    other = models.ProjectGroup.objects.create(name="other")
    models.ProjectGroupMembership.objects.create(
        group=other, project=campaign["member"], type=MEMBER
    )

    call_command("delete_campaign", "camp", "--yes")

    assert "KEPT" in capsys.readouterr().out
    assert list(models.Project.objects.values_list("name", flat=True)) == ["camp_x1"]
    assert (tmp_path / "camp_x1").exists()
    assert not (tmp_path / "camp_ref").exists()
    # It lost this campaign's membership and kept the other's.
    assert list(
        campaign["member"].group_memberships.values_list("group__name", flat=True)
    ) == ["other"]


def test_refuses_while_a_job_is_running(campaign, tmp_path):
    campaign["job"].status = models.Job.Status.RUNNING
    campaign["job"].save()

    with pytest.raises(CommandError, match="still active"):
        call_command("delete_campaign", "camp", "--yes")
    assert models.Project.objects.count() == 2
    assert (tmp_path / "camp_x1").exists()

    call_command("delete_campaign", "camp", "--yes", "--force")
    assert not models.Project.objects.exists()


def test_wrong_confirmation_deletes_nothing(campaign, monkeypatch):
    monkeypatch.setattr("builtins.input", lambda prompt="": "not the name")

    with pytest.raises(CommandError, match="did not match"):
        call_command("delete_campaign", "camp")
    assert models.ProjectGroup.objects.count() == 1


def test_directory_that_is_not_a_project_is_left(campaign, tmp_path, capsys):
    # Project.directory is a plain text column: a path that does not look like
    # a project must survive even though the records go.
    stray = tmp_path / "not_a_project"
    stray.mkdir()
    (stray / "thesis.tex").write_text("irreplaceable")
    models.Project.objects.filter(pk=campaign["member"].pk).update(directory=str(stray))

    call_command("delete_campaign", "camp", "--yes")

    assert "left on disk" in capsys.readouterr().out
    assert (stray / "thesis.tex").exists()
    assert not models.Project.objects.exists()


def test_unknown_name_lists_campaigns(campaign):
    with pytest.raises(CommandError, match="Campaigns: camp"):
        call_command("delete_campaign", "nonesuch")
