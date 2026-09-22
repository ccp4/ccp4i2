"""API tests for campaign sites and what was found at them.

Two endpoints and one payload addition:

* ``sites/<id>/evaluation/<project>/`` records, changes or withdraws a verdict.
* ``member_projects`` carries each dataset's findings, so the campaign
  overview renders a row per dataset without a request per row.

The distinction these exist to preserve: no row means nobody has looked, while
"empty" asserts that somebody looked and found nothing. A tag could not tell
those apart, and several of the tests below pin that they stay apart.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models


@pytest.fixture
def campaign(tmp_path):
    """A fragment campaign: one parent, two members, two sites."""
    group = models.ProjectGroup.objects.create(
        name="TestCampaign", type=models.ProjectGroup.GroupType.FRAGMENT_SET
    )

    def project(name):
        directory = tmp_path / name
        directory.mkdir(parents=True, exist_ok=True)
        return models.Project.objects.create(name=name, directory=str(directory))

    parent = project("ref")
    member_a = project("ds_a")
    member_b = project("ds_b")
    outsider = project("outsider")

    models.ProjectGroupMembership.objects.create(
        group=group, project=parent,
        type=models.ProjectGroupMembership.MembershipType.PARENT,
    )
    for member in (member_a, member_b):
        models.ProjectGroupMembership.objects.create(
            group=group, project=member,
            type=models.ProjectGroupMembership.MembershipType.MEMBER,
        )

    site_a = models.CampaignSite.objects.create(
        group=group, name="Pocket A", origin_x=1, origin_y=2, origin_z=3, order=0
    )
    site_b = models.CampaignSite.objects.create(
        group=group, name="Pocket B", origin_x=4, origin_y=5, origin_z=6, order=1
    )
    return {
        "group": group, "parent": parent, "a": member_a, "b": member_b,
        "outsider": outsider, "site_a": site_a, "site_b": site_b,
    }


def evaluation_url(group, site, project):
    return (
        f"/api/ccp4i2/projectgroups/{group.id}/sites/{site.id}"
        f"/evaluation/{project.id}/"
    )


@pytest.mark.django_db
class TestRecordingAVerdict:
    def test_records_a_hit(self, campaign):
        client = APIClient()
        response = client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit", "evaluator": "alice"},
            format="json",
        )
        assert response.status_code == 201, response.data
        assert response.data["verdict"] == "hit"

        stored = models.SiteEvaluation.objects.get(
            project=campaign["a"], site=campaign["site_a"]
        )
        assert stored.verdict == "hit"
        assert stored.evaluator == "alice"

    def test_changing_a_verdict_replaces_it(self, campaign):
        """One verdict per dataset per site, not a history of them."""
        client = APIClient()
        url = evaluation_url(campaign["group"], campaign["site_a"], campaign["a"])
        client.put(url, {"verdict": "unclear"}, format="json")
        response = client.put(url, {"verdict": "hit"}, format="json")

        assert response.status_code == 200, "second PUT updates rather than creates"
        assert (
            models.SiteEvaluation.objects.filter(
                project=campaign["a"], site=campaign["site_a"]
            ).count()
            == 1
        )

    def test_rejects_a_verdict_outside_the_vocabulary(self, campaign):
        client = APIClient()
        response = client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "probably"},
            format="json",
        )
        assert response.status_code == 400
        assert not models.SiteEvaluation.objects.exists()

    def test_refuses_a_project_outside_the_campaign(self, campaign):
        client = APIClient()
        response = client.put(
            evaluation_url(
                campaign["group"], campaign["site_a"], campaign["outsider"]
            ),
            {"verdict": "hit"},
            format="json",
        )
        assert response.status_code == 404
        assert not models.SiteEvaluation.objects.exists()

    def test_the_parent_can_be_evaluated_too(self, campaign):
        """The reference is a dataset like any other; it can hold a fragment."""
        client = APIClient()
        response = client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["parent"]),
            {"verdict": "hit"},
            format="json",
        )
        assert response.status_code == 201, response.data


@pytest.mark.django_db
class TestWithdrawingAVerdict:
    def test_delete_removes_the_row_rather_than_writing_empty(self, campaign):
        """Withdrawing is not the same as finding nothing.

        No row means nobody has looked; "empty" means somebody looked and
        found nothing. DELETE must restore the first state, not the second.
        """
        client = APIClient()
        url = evaluation_url(campaign["group"], campaign["site_a"], campaign["a"])
        client.put(url, {"verdict": "empty"}, format="json")

        response = client.delete(url)
        assert response.status_code == 204
        assert not models.SiteEvaluation.objects.filter(
            project=campaign["a"], site=campaign["site_a"]
        ).exists()

    def test_deleting_what_was_never_recorded_is_a_404(self, campaign):
        client = APIClient()
        response = client.delete(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"])
        )
        assert response.status_code == 404


@pytest.mark.django_db
class TestOverviewPayload:
    """member_projects carries the findings the campaign overview renders."""

    def _members(self, campaign):
        client = APIClient()
        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}/member_projects/"
        )
        assert response.status_code == 200, response.data
        return {row["name"]: row for row in response.data}

    def test_reports_progress_against_the_campaign_site_count(self, campaign):
        rows = self._members(campaign)
        assert rows["ds_a"]["sites_total"] == 2
        assert rows["ds_a"]["sites_evaluated"] == 0
        assert rows["ds_a"]["site_evaluations"] == []

    def test_lists_hits_and_unclears_but_not_empties(self, campaign):
        """A rich campaign has 30-40 sites, most empty for most datasets.

        Listing every verdict would put a 40-entry list on every row to render
        a cell showing two chips; "empty" stays recoverable from the counts.
        """
        client = APIClient()
        for site, verdict in (
            (campaign["site_a"], "hit"),
            (campaign["site_b"], "empty"),
        ):
            client.put(
                evaluation_url(campaign["group"], site, campaign["a"]),
                {"verdict": verdict},
                format="json",
            )

        row = self._members(campaign)["ds_a"]
        assert row["sites_evaluated"] == 2, "the empty still counts as looked at"
        listed = [e["site_name"] for e in row["site_evaluations"]]
        assert listed == ["Pocket A"], listed

    def test_hits_sort_before_unclears(self, campaign):
        """So a dataset with many unclears never has its hits pushed out of a
        capped chip list."""
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "unclear"}, format="json",
        )
        client.put(
            evaluation_url(campaign["group"], campaign["site_b"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )

        row = self._members(campaign)["ds_a"]
        assert [e["verdict"] for e in row["site_evaluations"]] == ["hit", "unclear"]

    def test_carries_the_site_id_for_navigation(self, campaign):
        """The chip links to that dataset at that site, so it needs the id --
        not the name, which is exactly what used to break on a rename."""
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )
        row = self._members(campaign)["ds_a"]
        assert row["site_evaluations"][0]["site_id"] == campaign["site_a"].id

    def test_one_dataset_s_findings_do_not_leak_into_another(self, campaign):
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )
        rows = self._members(campaign)
        assert len(rows["ds_a"]["site_evaluations"]) == 1
        assert rows["ds_b"]["site_evaluations"] == []
        assert rows["ds_b"]["sites_evaluated"] == 0


@pytest.mark.django_db
class TestRenamingASite:
    def test_a_rename_keeps_the_verdicts_attached(self, campaign):
        """The whole reason sites became rows with ids.

        Under the JSON list the association was the site's NAME, so renaming
        orphaned every reference to it.
        """
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )

        response = client.patch(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/sites/{campaign['site_a'].id}/",
            {"name": "Acetyl-lysine pocket"},
            format="json",
        )
        assert response.status_code == 200, response.data

        stored = models.SiteEvaluation.objects.get(
            project=campaign["a"], site=campaign["site_a"]
        )
        assert stored.verdict == "hit"
        assert stored.site.name == "Acetyl-lysine pocket"

    def test_deleting_a_site_takes_its_verdicts_with_it(self, campaign):
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )
        response = client.delete(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/sites/{campaign['site_a'].id}/"
        )
        assert response.status_code == 204
        assert not models.SiteEvaluation.objects.exists()


@pytest.mark.django_db
class TestReadingOneDatasetsVerdicts:
    """``evaluations/<project>/`` — the per-dataset view, empties included.

    The overview omits empties on purpose; a control that records verdicts
    cannot. Without them it would show "not looked at yet" for a site somebody
    had looked at and found nothing in, and write the wrong thing back.
    """

    def test_lists_every_verdict_including_empty(self, campaign):
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "hit"}, format="json",
        )
        client.put(
            evaluation_url(campaign["group"], campaign["site_b"], campaign["a"]),
            {"verdict": "empty"}, format="json",
        )

        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['a'].id}/"
        )
        assert response.status_code == 200, response.data
        verdicts = {row["site_id"]: row["verdict"] for row in response.data}
        assert verdicts == {
            campaign["site_a"].id: "hit",
            campaign["site_b"].id: "empty",
        }

    def test_a_site_nobody_looked_at_is_simply_absent(self, campaign):
        """The distinction the rows exist for, from the reading side."""
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "empty"}, format="json",
        )

        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['a'].id}/"
        )
        assert [row["site_id"] for row in response.data] == [campaign["site_a"].id]
        assert campaign["site_b"].id not in {r["site_id"] for r in response.data}

    def test_carries_the_note_and_evaluator(self, campaign):
        client = APIClient()
        client.put(
            evaluation_url(campaign["group"], campaign["site_a"], campaign["a"]),
            {"verdict": "unclear", "evaluator": "mn", "note": "density is weak"},
            format="json",
        )

        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['a'].id}/"
        )
        row = response.data[0]
        assert row["evaluator"] == "mn"
        assert row["note"] == "density is weak"
        assert row["project_id"] == campaign["a"].id

    def test_ordered_by_site(self, campaign):
        client = APIClient()
        for site in (campaign["site_b"], campaign["site_a"]):
            client.put(
                evaluation_url(campaign["group"], site, campaign["a"]),
                {"verdict": "hit"}, format="json",
            )

        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['a'].id}/"
        )
        assert [row["site_id"] for row in response.data] == [
            campaign["site_a"].id, campaign["site_b"].id
        ]

    def test_another_campaigns_verdicts_are_not_included(self, campaign, tmp_path):
        """A project can belong to more than one campaign."""
        other_group = models.ProjectGroup.objects.create(
            name="OtherCampaign", type=models.ProjectGroup.GroupType.FRAGMENT_SET
        )
        models.ProjectGroupMembership.objects.create(
            group=other_group, project=campaign["a"],
            type=models.ProjectGroupMembership.MembershipType.MEMBER,
        )
        other_site = models.CampaignSite.objects.create(
            group=other_group, name="Elsewhere",
            origin_x=0, origin_y=0, origin_z=0, order=0,
        )
        models.SiteEvaluation.objects.create(
            project=campaign["a"], site=other_site, verdict="hit"
        )
        models.SiteEvaluation.objects.create(
            project=campaign["a"], site=campaign["site_a"], verdict="empty"
        )

        response = APIClient().get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['a'].id}/"
        )
        assert [row["site_id"] for row in response.data] == [campaign["site_a"].id]

    def test_a_project_outside_the_campaign_is_refused(self, campaign):
        response = APIClient().get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}"
            f"/evaluations/{campaign['outsider'].id}/"
        )
        assert response.status_code == 404


@pytest.mark.django_db
class TestWhatDeletingASiteWouldCost:
    """The site list carries how much work each site is holding.

    Deleting a site cascades to every verdict recorded there, in every
    dataset, and cannot be undone. The confirmation is only worth reading if
    it can tell a site nobody has looked at from one carrying a campaign's
    worth of verdicts, so the count comes down with the list.
    """

    def sites(self, campaign):
        client = APIClient()
        response = client.get(
            f"/api/ccp4i2/projectgroups/{campaign['group'].id}/sites/"
        )
        assert response.status_code == 200
        return {site["name"]: site for site in response.json()}

    def test_counts_every_dataset_evaluated_at_that_site(self, campaign):
        for project in (campaign["a"], campaign["b"]):
            models.SiteEvaluation.objects.create(
                project=project, site=campaign["site_a"], verdict="hit"
            )
        models.SiteEvaluation.objects.create(
            project=campaign["a"], site=campaign["site_b"], verdict="empty"
        )

        sites = self.sites(campaign)
        assert sites["Pocket A"]["evaluation_count"] == 2
        assert sites["Pocket B"]["evaluation_count"] == 1

    def test_says_zero_rather_than_nothing_for_an_unexamined_site(self, campaign):
        # The dialog distinguishes "nothing recorded here" from "40 verdicts",
        # so an absent key and a zero must not look the same to it.
        sites = self.sites(campaign)
        assert sites["Pocket A"]["evaluation_count"] == 0
        assert sites["Pocket B"]["evaluation_count"] == 0

    def test_counts_empties_too(self, campaign):
        # An "empty" verdict is somebody's work: they looked and found
        # nothing. Deleting the site discards that finding like any other.
        models.SiteEvaluation.objects.create(
            project=campaign["a"], site=campaign["site_a"], verdict="empty"
        )
        assert self.sites(campaign)["Pocket A"]["evaluation_count"] == 1
