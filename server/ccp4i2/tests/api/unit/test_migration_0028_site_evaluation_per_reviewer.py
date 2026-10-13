"""Migration 0028 on a database that actually holds verdicts.

0028 moves SiteEvaluation's uniqueness from (project, site) to
(project, site, evaluator). A fresh test database has no verdicts, so it
would pass whatever the migration did to them; these tests go back to 0027,
record verdicts the way a campaign has them (one per cell, most with no
reviewer), and migrate forward.
"""

from django.db import IntegrityError, connection, transaction
from django.db.migrations.executor import MigrationExecutor
import pytest

APP = "ccp4i2"
BEFORE = [(APP, "0027_campaign_events_and_site_radius")]
AFTER = [(APP, "0028_site_evaluation_per_reviewer")]


def _migrate(targets):
    executor = MigrationExecutor(connection)
    executor.migrate(targets)
    return MigrationExecutor(connection).loader.project_state(targets).apps


@pytest.fixture
def apps_at_0027():
    yield _migrate(BEFORE)
    executor = MigrationExecutor(connection)
    executor.migrate(executor.loader.graph.leaf_nodes(APP))


def _campaign(apps):
    Project = apps.get_model(APP, "Project")
    group = apps.get_model(APP, "ProjectGroup").objects.create(name="Mpro campaign")
    site = apps.get_model(APP, "CampaignSite").objects.create(
        group=group, name="Active site", origin_x=1, origin_y=2, origin_z=3)
    soaks = [Project.objects.create(name=f"soak_{i}", directory=f"/nowhere/soak_{i}")
             for i in range(2)]
    return site, soaks


def test_existing_verdicts_survive(apps_at_0027):
    site, (soak_1, soak_2) = _campaign(apps_at_0027)
    SiteEvaluation = apps_at_0027.get_model(APP, "SiteEvaluation")
    SiteEvaluation.objects.create(project=soak_1, site=site, verdict="hit", evaluator="")
    SiteEvaluation.objects.create(project=soak_2, site=site, verdict="empty",
                                  evaluator="mn@example.org", note="nothing at 1 sigma")

    apps = _migrate(AFTER)

    rows = apps.get_model(APP, "SiteEvaluation").objects.order_by("project__name")
    assert [(r.project.name, r.verdict, r.evaluator, r.note) for r in rows] == [
        ("soak_0", "hit", "", ""),
        ("soak_1", "empty", "mn@example.org", "nothing at 1 sigma"),
    ]


def test_a_second_reviewer_fits_beside_the_first_afterwards(apps_at_0027):
    site, (soak, _) = _campaign(apps_at_0027)
    apps_at_0027.get_model(APP, "SiteEvaluation").objects.create(
        project=soak, site=site, verdict="hit", evaluator="")

    apps = _migrate(AFTER)
    SiteEvaluation = apps.get_model(APP, "SiteEvaluation")

    SiteEvaluation.objects.create(project_id=soak.pk, site_id=site.pk, verdict="empty",
                                  evaluator="alice@example.org")
    SiteEvaluation.objects.create(project_id=soak.pk, site_id=site.pk, verdict="unclear",
                                  evaluator="bob@example.org")
    assert SiteEvaluation.objects.filter(project_id=soak.pk, site_id=site.pk).count() == 3

    # Still one per reviewer.
    with pytest.raises(IntegrityError), transaction.atomic():
        SiteEvaluation.objects.create(project_id=soak.pk, site_id=site.pk, verdict="hit",
                                      evaluator="alice@example.org")
