"""Migration 0026 gives every existing campaign and site its own uuid.

Seeded on purpose: a RunPython step passes on every empty database, and the
failure this guards -- AddField(default=uuid4) stamping one value onto every
existing row, then the unique constraint refusing it -- only shows with rows
present. Run against the migration state before 0026, seeded, then migrated.
"""

from django.db import connection
from django.db.migrations.executor import MigrationExecutor
from django.test import TransactionTestCase

APP = "ccp4i2"
BEFORE = (APP, "0025_site_tags_become_evaluations")
AFTER = (APP, "0026_campaign_uuids")


class CampaignUuidMigrationTest(TransactionTestCase):
    def setUp(self):
        self.executor = MigrationExecutor(connection)
        self.executor.migrate([BEFORE])
        self.apps_before = self.executor.loader.project_state([BEFORE]).apps

    def tearDown(self):
        executor = MigrationExecutor(connection)
        executor.migrate(executor.loader.graph.leaf_nodes())

    def test_existing_rows_get_distinct_uuids(self):
        ProjectGroup = self.apps_before.get_model(APP, "ProjectGroup")
        CampaignSite = self.apps_before.get_model(APP, "CampaignSite")
        groups = [ProjectGroup.objects.create(name=f"campaign {i}") for i in range(2)]
        for group in groups:
            for j in range(2):
                CampaignSite.objects.create(
                    group=group, name=f"site {j}", origin_x=j, origin_y=0, origin_z=0)

        executor = MigrationExecutor(connection)
        executor.loader.build_graph()
        executor.migrate([AFTER])
        apps_after = executor.loader.project_state([AFTER]).apps

        group_uuids = list(apps_after.get_model(APP, "ProjectGroup").objects.values_list("uuid", flat=True))
        site_uuids = list(apps_after.get_model(APP, "CampaignSite").objects.values_list("uuid", flat=True))
        self.assertEqual(len(group_uuids), 2)
        self.assertEqual(len(site_uuids), 4)
        self.assertTrue(all(group_uuids) and all(site_uuids))
        self.assertEqual(len(set(group_uuids + site_uuids)), 6, "every row its own uuid")
