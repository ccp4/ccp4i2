"""Migration 0024 must work on a database that actually has campaign sites.

As released in 3.1.0a71 it did not. Its data step created ``CampaignSite``
rows with ``group=...`` before the operation that adds the ``group`` column,
so it raised ``TypeError: CampaignSite() got unexpected keyword arguments:
'group'`` on any database holding a campaign with sites -- and passed on every
database without one, because the loop body never ran. Every test database is
of the second kind, which is how it shipped.

So these tests do what no fresh test database does: go back to 0023, put
sites in the JSON column the way a real campaign has them, and migrate forward.
"""

import pytest
from django.db import connection
from django.db.migrations.executor import MigrationExecutor

APP = "ccp4i2"
BEFORE = [(APP, "0023_stagedupload")]
AFTER = [(APP, "0024_campaign_sites_and_evaluations")]


def _migrate(targets):
    executor = MigrationExecutor(connection)
    executor.migrate(targets)
    # A fresh executor, so the returned state is read from what was applied.
    return MigrationExecutor(connection).loader.project_state(targets).apps


@pytest.fixture
def apps_at_0023():
    """The schema as it was before 0024, on this test's throwaway database."""
    yield _migrate(BEFORE)
    # Leave the database at head: the conftest teardown expects nothing of it,
    # but a half-migrated file is a confusing thing to find after a failure.
    executor = MigrationExecutor(connection)
    executor.migrate(executor.loader.graph.leaf_nodes(APP))


def test_sites_in_json_become_rows(apps_at_0023):
    ProjectGroup = apps_at_0023.get_model(APP, "ProjectGroup")
    group = ProjectGroup.objects.create(
        name="Mpro campaign",
        sites=[
            {"name": "Active site", "origin": [1.0, 2.0, 3.0],
             "quat": [0, 0, 0, 1], "zoom": 0.5},
            {"name": "Allosteric", "origin": [4, 5, 6]},
        ],
    )
    ProjectGroup.objects.create(name="No sites yet", sites=[])

    apps = _migrate(AFTER)

    CampaignSite = apps.get_model(APP, "CampaignSite")
    rows = list(CampaignSite.objects.order_by("order"))
    assert [(r.group_id, r.name, r.order) for r in rows] == [
        (group.pk, "Active site", 0),
        (group.pk, "Allosteric", 1),
    ]
    assert (rows[0].origin_x, rows[0].origin_y, rows[0].origin_z) == (1.0, 2.0, 3.0)
    assert rows[0].quat == [0, 0, 0, 1]
    assert rows[0].zoom == 0.5
    assert rows[1].quat is None and rows[1].zoom is None


def test_duplicate_and_unusable_sites(apps_at_0023):
    """The JSON list never constrained names; the table does."""
    ProjectGroup = apps_at_0023.get_model(APP, "ProjectGroup")
    ProjectGroup.objects.create(
        name="Messy campaign",
        sites=[
            {"name": "Pocket", "origin": [1, 1, 1]},
            {"name": "Pocket", "origin": [2, 2, 2]},
            {"name": "No origin"},
            "not even a dict",
        ],
    )

    apps = _migrate(AFTER)

    CampaignSite = apps.get_model(APP, "CampaignSite")
    assert list(
        CampaignSite.objects.order_by("order").values_list("name", flat=True)
    ) == ["Pocket", "Pocket (2)"]


def test_reverse_puts_the_sites_back(apps_at_0023):
    ProjectGroup = apps_at_0023.get_model(APP, "ProjectGroup")
    group = ProjectGroup.objects.create(
        name="Round trip",
        sites=[{"name": "Active site", "origin": [1.0, 2.0, 3.0], "zoom": 0.5}],
    )

    _migrate(AFTER)
    apps = _migrate(BEFORE)

    restored = apps.get_model(APP, "ProjectGroup").objects.get(pk=group.pk)
    assert restored.sites == [
        {"name": "Active site", "origin": [1.0, 2.0, 3.0], "zoom": 0.5}
    ]


def test_site_tags_become_hits_all_the_way_to_head(apps_at_0023):
    """0024 then 0025, with data: the path a real campaign database takes.

    0025 reads the rows 0024 writes, and like 0024 it does nothing at all on a
    database with no sites, so nothing else exercises the pair together.
    """
    Project = apps_at_0023.get_model(APP, "Project")
    ProjectGroup = apps_at_0023.get_model(APP, "ProjectGroup")
    Membership = apps_at_0023.get_model(APP, "ProjectGroupMembership")
    ProjectTag = apps_at_0023.get_model(APP, "ProjectTag")

    group = ProjectGroup.objects.create(
        name="Mpro campaign",
        sites=[{"name": "Active site", "origin": [1, 2, 3]}],
    )
    soak = Project.objects.create(name="soak_1", directory="/nowhere/soak_1")
    Membership.objects.create(group=group, project=soak, type="member")

    site_tag = ProjectTag.objects.create(text="Active site", path="Active site")
    site_tag.projects.add(soak)
    other_tag = ProjectTag.objects.create(text="to revisit", path="to revisit")
    other_tag.projects.add(soak)

    executor = MigrationExecutor(connection)
    head = executor.loader.graph.leaf_nodes(APP)
    apps = _migrate(head)

    evaluation = apps.get_model(APP, "SiteEvaluation").objects.get()
    assert (evaluation.project_id, evaluation.site.name) == (soak.pk, "Active site")
    assert (evaluation.verdict, evaluation.evaluator) == ("hit", "migrated")

    # The site tag is spent; the ordinary tag is none of 0025's business.
    remaining = apps.get_model(APP, "ProjectTag").objects
    assert list(remaining.values_list("text", flat=True)) == ["to revisit"]
