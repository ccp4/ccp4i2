"""Tests for the tag-tree API: the forest endpoint, bulk membership, and
filtering the project list by a node of the tree.

The counts are the part worth testing hard. An ancestor node reports distinct
projects in its subtree, not the sum of its descendants' counts, because a
project filed on both a tag and its child must not be counted twice.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db.models import Project, ProjectTag

API_PREFIX = "/api/ccp4i2"


class TestProjectTagTreeAPI:
    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions):
        self.client = APIClient()

        self.root = ProjectTag.objects.create(text="SARS-CoV-2")
        self.mpro = ProjectTag.objects.create(text="Mpro", parent=self.root)
        self.soaks = ProjectTag.objects.create(text="soaks", parent=self.mpro)
        self.other = ProjectTag.objects.create(text="in-house")

        self.top = Project.objects.create(name="top", directory="/tmp/top")
        self.deep = Project.objects.create(name="deep", directory="/tmp/deep")
        self.both = Project.objects.create(name="both", directory="/tmp/both")
        self.loose = Project.objects.create(name="loose", directory="/tmp/loose")

        self.root.projects.add(self.top)
        self.soaks.projects.add(self.deep)
        self.mpro.projects.add(self.both)
        self.soaks.projects.add(self.both)

        self.tree_url = f"{API_PREFIX}/projecttags/tree/"

    # -- tree ------------------------------------------------------------

    def test_tree_is_nested_from_the_roots(self):
        response = self.client.get(self.tree_url)
        assert response.status_code == 200
        payload = response.json()

        roots = {node["text"]: node for node in payload["tags"]}
        assert set(roots) == {"SARS-CoV-2", "in-house"}

        mpro = roots["SARS-CoV-2"]["children"][0]
        assert mpro["text"] == "Mpro"
        assert mpro["children"][0]["text"] == "soaks"
        assert mpro["children"][0]["display_path"] == "SARS-CoV-2/Mpro/soaks"

    def test_direct_count_is_what_is_filed_on_the_node(self):
        response = self.client.get(self.tree_url)
        roots = {node["text"]: node for node in response.json()["tags"]}
        assert roots["SARS-CoV-2"]["project_count"] == 1

    def test_rolled_up_count_deduplicates_across_the_subtree(self):
        """`both` is on Mpro and on Mpro/soaks; the root must count it once."""
        response = self.client.get(self.tree_url)
        roots = {node["text"]: node for node in response.json()["tags"]}
        assert roots["SARS-CoV-2"]["total_project_count"] == 3

        mpro = roots["SARS-CoV-2"]["children"][0]
        assert mpro["total_project_count"] == 2

    def test_untagged_projects_are_reported(self):
        response = self.client.get(self.tree_url)
        assert response.json()["untagged_project_count"] == 1

    # -- bulk membership -------------------------------------------------

    def test_add_projects_in_bulk(self):
        response = self.client.post(
            f"{API_PREFIX}/projecttags/{self.other.id}/add_projects/",
            {"project_ids": [self.top.id, self.loose.id]},
            format="json",
        )
        assert response.status_code == 200
        assert response.json()["changed"] == 2
        assert set(self.other.projects.all()) == {self.top, self.loose}

    def test_add_projects_reports_unknown_ids(self):
        response = self.client.post(
            f"{API_PREFIX}/projecttags/{self.other.id}/add_projects/",
            {"project_ids": [self.top.id, 999999]},
            format="json",
        )
        assert response.status_code == 200
        assert response.json()["missing"] == [999999]
        assert response.json()["changed"] == 1

    def test_add_projects_rejects_an_empty_list(self):
        response = self.client.post(
            f"{API_PREFIX}/projecttags/{self.other.id}/add_projects/",
            {"project_ids": []},
            format="json",
        )
        assert response.status_code == 400

    def test_remove_projects_in_bulk(self):
        response = self.client.post(
            f"{API_PREFIX}/projecttags/{self.soaks.id}/remove_projects/",
            {"project_ids": [self.deep.id, self.both.id]},
            format="json",
        )
        assert response.status_code == 200
        assert list(self.soaks.projects.all()) == []
        # Removing from one tag leaves the project's other tags alone.
        assert set(self.both.tags.all()) == {self.mpro}

    # -- filtering the project list --------------------------------------

    def _project_names(self, response):
        payload = response.json()
        results = payload["results"] if isinstance(payload, dict) else payload
        return {project["name"] for project in results}

    def test_filter_by_tag_rolls_up_by_default(self):
        response = self.client.get(f"{API_PREFIX}/projects/?tag={self.root.id}")
        assert response.status_code == 200
        assert self._project_names(response) == {"top", "deep", "both"}

    def test_filter_by_tag_can_be_restricted_to_direct_members(self):
        response = self.client.get(
            f"{API_PREFIX}/projects/?tag={self.root.id}&descendants=false"
        )
        assert self._project_names(response) == {"top"}

    def test_filter_by_untagged(self):
        response = self.client.get(f"{API_PREFIX}/projects/?untagged=true")
        assert self._project_names(response) == {"loose"}

    def test_filter_by_unknown_tag_returns_nothing(self):
        response = self.client.get(f"{API_PREFIX}/projects/?tag=999999")
        assert self._project_names(response) == set()

    def test_rolled_up_filter_does_not_duplicate_a_multiply_tagged_project(self):
        response = self.client.get(f"{API_PREFIX}/projects/?tag={self.mpro.id}")
        payload = response.json()
        results = payload["results"] if isinstance(payload, dict) else payload
        names = [project["name"] for project in results]
        assert sorted(names) == ["both", "deep"]


class TestProjectTagEditingAPI:
    """Renaming, re-parenting and deleting nodes of the tag forest.

    These all arrive as partial updates, which is where the uniqueness check
    used to fall through: with only `parent` in the payload it compared against
    a null `text` and let a collision reach the unique index as a 500.
    """

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions):
        self.client = APIClient()
        self.a = ProjectTag.objects.create(text="A")
        self.b = ProjectTag.objects.create(text="B", parent=self.a)
        self.c = ProjectTag.objects.create(text="C", parent=self.b)
        self.other = ProjectTag.objects.create(text="X")

    def _url(self, tag):
        return f"{API_PREFIX}/projecttags/{tag.id}/"

    def test_create_a_child_tag(self):
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            {"text": "D", "parent": self.c.id, "projects": []},
            format="json",
        )
        assert response.status_code == 201
        assert response.json()["display_path"] == "A/B/C/D"

    def test_rename_rewrites_the_subtree(self):
        response = self.client.patch(self._url(self.a), {"text": "Alpha"}, format="json")
        assert response.status_code == 200

        self.c.refresh_from_db()
        assert self.c.display_path == "Alpha/B/C"

    def test_rename_onto_an_existing_sibling_is_rejected(self):
        ProjectTag.objects.create(text="B2", parent=self.a)
        response = self.client.patch(
            f"{API_PREFIX}/projecttags/{self.b.id}/", {"text": "B2"}, format="json"
        )
        assert response.status_code == 400

    def test_reparent_moves_the_subtree(self):
        response = self.client.patch(
            self._url(self.b), {"parent": self.other.id}, format="json"
        )
        assert response.status_code == 200

        self.c.refresh_from_db()
        assert self.c.display_path == "X/B/C"

    def test_promote_to_root(self):
        response = self.client.patch(self._url(self.b), {"parent": None}, format="json")
        assert response.status_code == 200

        self.c.refresh_from_db()
        assert self.c.display_path == "B/C"

    def test_reparent_onto_a_descendant_is_rejected_not_a_500(self):
        response = self.client.patch(
            self._url(self.a), {"parent": self.c.id}, format="json"
        )
        assert response.status_code == 400
        assert "parent" in response.json()

    def test_reparent_onto_itself_is_rejected(self):
        response = self.client.patch(
            self._url(self.a), {"parent": self.a.id}, format="json"
        )
        assert response.status_code == 400

    def test_reparent_into_a_name_collision_is_rejected_not_a_500(self):
        """`other` already has a `B`, so `A/B` cannot move under it."""
        ProjectTag.objects.create(text="B", parent=self.other)
        response = self.client.patch(
            self._url(self.b), {"parent": self.other.id}, format="json"
        )
        assert response.status_code == 400

    def test_delete_removes_the_whole_subtree(self):
        project = Project.objects.create(name="p", directory="/tmp/p")
        self.c.projects.add(project)

        response = self.client.delete(self._url(self.a))
        assert response.status_code == 204

        assert ProjectTag.objects.filter(text__in=["A", "B", "C"]).count() == 0
        # The projects themselves survive; only the labels went.
        assert Project.objects.filter(name="p").exists()
