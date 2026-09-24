"""Campaign state survives losing the database.

Until this landed, ProjectGroup, CampaignSite and SiteEvaluation appeared in
none of the snapshot, export, restore or signal code: lose db.sqlite3 and
every verdict in every campaign was gone with the project directories intact.
These tests pin the ownership split (the campaign and its sites with the
parent, verdicts with the dataset they are about), the round trip, the
ordering that makes the round trip work, and that the serialiser is generic
enough for a column to be added later without touching it.
"""

from datetime import timedelta
from pathlib import Path
from shutil import rmtree
from unittest.mock import patch
from xml.etree import ElementTree as ET

from django.test import TestCase, override_settings
from django.utils import timezone

from ...db import project_snapshot
from ...db.models import (
    CampaignSite,
    File,
    Job,
    Project,
    ProjectGroup,
    ProjectGroupMembership,
    SiteEvaluation,
)
from ...db.project_snapshot import SNAPSHOT_NAME, write_snapshot
from ...db.restore_project import inspect, restore_all, restore_from_directory

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_CAMPAIGN_SNAPSHOT_TEST_DIR"
HOME_DIR = PROJECTS_DIR / "home"

PARENT = ProjectGroupMembership.MembershipType.PARENT
MEMBER = ProjectGroupMembership.MembershipType.MEMBER


def make_project(name: str) -> Project:
    directory = PROJECTS_DIR / name
    (directory / "CCP4_JOBS").mkdir(parents=True, exist_ok=True)
    return Project.objects.create(name=name, directory=str(directory))


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class CampaignSnapshotBase(TestCase):
    def setUp(self):
        PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
        HOME_DIR.mkdir(parents=True, exist_ok=True)
        registry = patch.object(
            project_snapshot, "registry_path",
            return_value=HOME_DIR / project_snapshot.REGISTRY_NAME,
        )
        registry.start()
        self.addCleanup(registry.stop)
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)

        # Build the fixture with snapshot scheduling suspended: a TestCase runs
        # inside one transaction, and a snapshot queued here would sit in the
        # on-commit list for the whole test and make schedule_snapshot treat
        # every later change to the same project as already covered.
        with project_snapshot.suspended():
            self._build_campaign()
        for project in (self.parent, self.m1, self.m2):
            write_snapshot(project)

    def _build_campaign(self):
        self.parent = make_project("BAZ2B-ref")
        self.m1 = make_project("BAZ2B-x425")
        self.m2 = make_project("BAZ2B-x427")
        self.group = ProjectGroup.objects.create(
            name="BAZ2B", type=ProjectGroup.GroupType.FRAGMENT_SET)
        for project, kind in ((self.parent, PARENT), (self.m1, MEMBER), (self.m2, MEMBER)):
            ProjectGroupMembership.objects.create(group=self.group, project=project, type=kind)
        self.site1 = CampaignSite.objects.create(
            group=self.group, name="Acetyl-lysine pocket", origin_x=15.5, origin_y=40.6,
            origin_z=30.0, quat=[0.0, 0.7, 0.0, 0.7], zoom=0.8, order=0)
        self.site2 = CampaignSite.objects.create(
            group=self.group, name="ZA channel", origin_x=1.0, origin_y=2.0, origin_z=3.0, order=1)
        self.e1 = SiteEvaluation.objects.create(
            project=self.m1, site=self.site1, verdict=SiteEvaluation.Verdict.HIT,
            evaluator="martin", note="clear density for the fragment")
        self.e2 = SiteEvaluation.objects.create(
            project=self.m1, site=self.site2, verdict=SiteEvaluation.Verdict.EMPTY, evaluator="martin")
        self.e3 = SiteEvaluation.objects.create(
            project=self.m2, site=self.site1, verdict=SiteEvaluation.Verdict.UNCLEAR, note="partial")
        # A verdict recorded a while ago must come back with its own time.
        SiteEvaluation.objects.filter(pk=self.e1.pk).update(
            evaluated_at=timezone.now() - timedelta(days=40))
        self.e1.refresh_from_db()

        self.uuids = {
            "group": self.group.uuid, "site1": self.site1.uuid, "site2": self.site2.uuid,
            "parent": self.parent.uuid, "m1": self.m1.uuid, "m2": self.m2.uuid,
        }

    def lose_the_database(self):
        with project_snapshot.suspended():
            SiteEvaluation.objects.all().delete()
            CampaignSite.objects.all().delete()
            ProjectGroupMembership.objects.all().delete()
            ProjectGroup.objects.all().delete()
            File.objects.all().delete()
            Job.objects.all().delete()
            Project.objects.all().delete()

    @staticmethod
    def snapshot_root(project):
        return ET.parse(Path(project.directory) / SNAPSHOT_NAME).getroot()


class OwnershipTest(CampaignSnapshotBase):
    """What goes in whose snapshot."""

    def test_parent_snapshot_carries_the_campaign_its_roster_and_its_sites(self):
        root = self.snapshot_root(self.parent)
        campaign = root.find("ccp4i2_body/campaignTable/campaign")
        self.assertEqual(campaign.get("uuid"), self.group.uuid.hex)
        self.assertEqual(campaign.get("name"), "BAZ2B")
        self.assertEqual(campaign.get("type"), "fragment_set")
        roster = {(m.get("projectid"), m.get("type")) for m in campaign.findall("membership")}
        self.assertEqual(roster, {
            (self.parent.uuid.hex, "parent"), (self.m1.uuid.hex, "member"), (self.m2.uuid.hex, "member")})
        sites = campaign.findall("site")
        self.assertEqual([s.get("name") for s in sites], ["Acetyl-lysine pocket", "ZA channel"])
        self.assertEqual(sites[0].get("uuid"), self.site1.uuid.hex)
        self.assertEqual(sites[0].get("quat"), "[0.0, 0.7, 0.0, 0.7]")
        self.assertIsNone(sites[1].get("quat"), "an unset nullable column is omitted")
        self.assertIsNone(root.find("ccp4i2_body/siteevaluationTable"),
                          "verdicts on members are not the parent's to hold")

    def test_member_snapshot_carries_its_membership_and_its_own_verdicts_only(self):
        root = self.snapshot_root(self.m1)
        self.assertIsNone(root.find("ccp4i2_body/campaignTable"))
        membership = root.find("ccp4i2_body/campaignmembershipTable/campaignmembership")
        self.assertEqual(membership.get("campaignuuid"), self.group.uuid.hex)
        self.assertEqual(membership.get("type"), "member")
        verdicts = root.findall("ccp4i2_body/siteevaluationTable/siteevaluation")
        self.assertEqual({v.get("siteuuid") for v in verdicts}, {self.site1.uuid.hex, self.site2.uuid.hex})
        hit = next(v for v in verdicts if v.get("siteuuid") == self.site1.uuid.hex)
        self.assertEqual(hit.get("verdict"), "hit")
        self.assertEqual(hit.get("note"), "clear density for the fragment")
        self.assertEqual(hit.get("evaluated_at"), self.e1.evaluated_at.isoformat())

    def test_inspect_counts_campaign_state(self):
        self.assertEqual(inspect(Path(self.parent.directory)).campaigns, 1)
        self.assertEqual(inspect(Path(self.m1.directory)).evaluations, 2)


class RoundTripTest(CampaignSnapshotBase):
    """Destroy, restore, and the verdicts are still there."""

    def test_everything_comes_back_from_the_directories(self):
        self.lose_the_database()
        self.assertEqual(SiteEvaluation.objects.count(), 0)

        # Members listed first, to prove the ordering is the restore's job.
        result = restore_all([Path(self.m1.directory), Path(self.m2.directory),
                              Path(self.parent.directory)])
        self.assertEqual(len(result["restored"]), 3, result)
        self.assertEqual([], [w for r in result["restored"] for w in r["warnings"]])

        group = ProjectGroup.objects.get(uuid=self.uuids["group"])
        self.assertEqual((group.name, group.type), ("BAZ2B", "fragment_set"))
        self.assertEqual(group.parent_project.uuid, self.uuids["parent"])
        self.assertEqual(
            {m.project.uuid for m in group.memberships.filter(type=MEMBER)},
            {self.uuids["m1"], self.uuids["m2"]})

        site1 = CampaignSite.objects.get(uuid=self.uuids["site1"])
        self.assertEqual(site1.group, group)
        self.assertEqual((site1.name, site1.origin, site1.quat, site1.zoom, site1.order),
                         ("Acetyl-lysine pocket", [15.5, 40.6, 30.0], [0.0, 0.7, 0.0, 0.7], 0.8, 0))
        site2 = CampaignSite.objects.get(uuid=self.uuids["site2"])
        self.assertEqual((site2.quat, site2.zoom, site2.order), (None, None, 1))

        m1 = Project.objects.get(uuid=self.uuids["m1"])
        hit = SiteEvaluation.objects.get(project=m1, site=site1)
        self.assertEqual((hit.verdict, hit.evaluator, hit.note),
                         ("hit", "martin", "clear density for the fragment"))
        self.assertEqual(hit.evaluated_at, self.e1.evaluated_at, "the recorded time, not now")
        self.assertEqual(SiteEvaluation.objects.get(project=m1, site=site2).verdict, "empty")
        m2 = Project.objects.get(uuid=self.uuids["m2"])
        self.assertEqual(SiteEvaluation.objects.get(project=m2, site=site1).note, "partial")
        self.assertEqual(SiteEvaluation.objects.count(), 3)

    def test_a_member_restored_before_its_parent_says_so_and_recovers_later(self):
        self.lose_the_database()
        report = restore_from_directory(Path(self.m1.directory))
        self.assertTrue(report.restored)
        self.assertEqual(SiteEvaluation.objects.count(), 0)
        self.assertTrue(any("parent project" in w for w in report.warnings), report.warnings)

        restore_from_directory(Path(self.parent.directory))
        self.assertEqual(ProjectGroupMembership.objects.filter(type=MEMBER).count(), 1,
                         "the parent's roster rejoined the member that is already back")

        report = restore_from_directory(Path(self.m1.directory), replace=True)
        self.assertEqual(report.warnings, [])
        m1 = Project.objects.get(uuid=self.uuids["m1"])
        self.assertEqual(SiteEvaluation.objects.filter(project=m1).count(), 2)

    def test_a_snapshot_never_drops_verdicts_the_database_cannot_hold_yet(self):
        """The hazard the ordering rule exists for, closed at the source: a
        member restored ahead of its parent keeps carrying its verdicts in its
        snapshot until they can be placed."""
        self.lose_the_database()
        restore_from_directory(Path(self.m1.directory))    # rewrote m1's snapshot
        root = self.snapshot_root(Project.objects.get(uuid=self.uuids["m1"]))
        verdicts = root.findall("ccp4i2_body/siteevaluationTable/siteevaluation")
        self.assertEqual({v.get("siteuuid") for v in verdicts},
                         {self.uuids["site1"].hex, self.uuids["site2"].hex})
        self.assertEqual(
            root.find("ccp4i2_body/campaignmembershipTable/campaignmembership").get("campaignuuid"),
            self.uuids["group"].hex)
        # and a later, unrelated rewrite keeps carrying them
        m1 = Project.objects.get(uuid=self.uuids["m1"])
        m1.description = "edited while the parent was still missing"
        m1.save()
        write_snapshot(m1)
        root = self.snapshot_root(m1)
        self.assertEqual(len(root.findall("ccp4i2_body/siteevaluationTable/siteevaluation")), 2)

    def test_restoring_twice_changes_nothing(self):
        self.lose_the_database()
        dirs = [Path(p.directory) for p in (self.parent, self.m1, self.m2)]
        restore_all(dirs)
        restore_all(dirs, replace=True)
        self.assertEqual(ProjectGroup.objects.count(), 1)
        self.assertEqual(CampaignSite.objects.count(), 2)
        self.assertEqual(SiteEvaluation.objects.count(), 3)
        self.assertEqual(ProjectGroupMembership.objects.count(), 3)

    def test_a_campaign_name_held_by_a_different_campaign_is_relabelled_not_merged(self):
        self.lose_the_database()
        other = ProjectGroup.objects.create(name="BAZ2B")   # a different uuid
        restore_all([Path(self.parent.directory)])
        restored = ProjectGroup.objects.get(uuid=self.uuids["group"])
        self.assertEqual(restored.name, "BAZ2B (restored)")
        self.assertEqual(ProjectGroup.objects.get(pk=other.pk).site_set.count(), 0)
        self.assertEqual(restored.site_set.count(), 2)

    def test_the_serialiser_is_generic(self):
        """A column added later rides through; an attribute the model lacks is
        ignored; a field the snapshot predates keeps its default."""
        path = Path(self.m1.directory) / SNAPSHOT_NAME
        root = ET.parse(path).getroot()
        for node in root.findall("ccp4i2_body/siteevaluationTable/siteevaluation"):
            node.set("evidence_file", "some-future-column")
            node.attrib.pop("note", None)
        ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)

        self.lose_the_database()
        result = restore_all([Path(self.parent.directory), Path(self.m1.directory)])
        self.assertEqual(len(result["restored"]), 2)
        m1 = Project.objects.get(uuid=self.uuids["m1"])
        hit = SiteEvaluation.objects.get(project=m1, site__uuid=self.uuids["site1"])
        self.assertEqual((hit.verdict, hit.note), ("hit", ""))


class SignalTest(CampaignSnapshotBase):
    """A change to campaign state rewrites the right project's snapshot."""

    def scheduled_by(self, action):
        written = []
        with patch.object(project_snapshot, "write_snapshot",
                          side_effect=lambda project: written.append(project.name)):
            with self.captureOnCommitCallbacks(execute=True):
                action()
        return sorted(set(written))

    def test_a_verdict_snapshots_the_dataset(self):
        names = self.scheduled_by(lambda: SiteEvaluation.objects.create(
            project=self.m2, site=self.site2, verdict=SiteEvaluation.Verdict.HIT))
        self.assertEqual(names, ["BAZ2B-x427"])

    def test_moving_a_site_snapshots_the_parent(self):
        def move():
            self.site1.origin_x = 16.0
            self.site1.save()
        self.assertEqual(self.scheduled_by(move), ["BAZ2B-ref"])

    def test_deleting_a_site_snapshots_the_parent_and_the_evaluated_datasets(self):
        names = self.scheduled_by(lambda: self.site2.delete())
        self.assertEqual(names, ["BAZ2B-ref", "BAZ2B-x425"])

    def test_adding_a_member_snapshots_both_sides(self):
        with project_snapshot.suspended():   # its own creation is not the point
            newcomer = make_project("BAZ2B-x428")
        names = self.scheduled_by(lambda: ProjectGroupMembership.objects.create(
            group=self.group, project=newcomer, type=MEMBER))
        self.assertEqual(names, ["BAZ2B-ref", "BAZ2B-x428"])

    def test_renaming_the_campaign_snapshots_the_parent(self):
        def rename():
            self.group.name = "BAZ2B bromodomain"
            self.group.save()
        self.assertEqual(self.scheduled_by(rename), ["BAZ2B-ref"])
