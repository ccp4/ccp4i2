"""
Tests for ProjectTag API endpoints.

Converted to pytest fixture-based approach for compatibility with pytest-xdist
parallel test execution. Uses the isolated_test_db fixture from conftest.py.
"""
import json
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.db.models import Project, ProjectTag

# API URL prefix - all API endpoints are under /api/ccp4i2/
API_PREFIX = "/api/ccp4i2"


class TestProjectTagAPI:
    """Tests for ProjectTag API endpoints using pytest fixtures."""

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions):
        """Set up test fixtures for each test."""
        self.client = APIClient()

        # Create a test project
        self.project = Project.objects.create(
            name="Test Project",
            directory="/tmp/test_project",
        )

        # Create a test tag
        self.tag = ProjectTag.objects.create(
            text="test-tag",
        )

        # Base URLs for API endpoints
        self.projecttags_url = f"{API_PREFIX}/projecttags/"
        self.project_tags_url = f"{API_PREFIX}/projects/{self.project.uuid}/tags/"

    def test_create_project_tag(self):
        """Test creating a new ProjectTag via API"""
        # Test data for creating a new tag
        tag_data = {"text": "Test Tag", "parent": None, "projects": []}

        # Make POST request to create the tag
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(tag_data),
        )

        # Verify the response
        assert response.status_code == 201  # HTTP 201 Created

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        assert created_tag["text"] == "Test Tag"
        assert created_tag["parent"] is None
        assert created_tag["projects"] == []
        assert created_tag["id"] is not None

        # Verify tag exists in database
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        assert db_tag.text == "Test Tag"
        assert db_tag.parent is None

    def test_create_project_tag_with_parent(self):
        """Test creating a hierarchical ProjectTag with parent"""
        # First create a parent tag
        parent_tag = models.ProjectTag.objects.create(text="Parent Tag", parent=None)

        # Test data for creating a child tag
        child_tag_data = {"text": "Child Tag", "parent": parent_tag.id, "projects": []}

        # Make POST request to create the child tag
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(child_tag_data),
        )

        # Verify the response
        assert response.status_code == 201

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        assert created_tag["text"] == "Child Tag"
        assert created_tag["parent"] == parent_tag.id

        # Verify tag exists in database with correct parent
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        assert db_tag.text == "Child Tag"
        assert db_tag.parent.id == parent_tag.id

    def test_create_project_tag_with_projects(self):
        """Test creating a ProjectTag and associating it with projects"""
        # Create a second test project
        second_project = models.Project.objects.create(
            name="Second Test Project",
            description="Another test project",
            directory="/tmp/second_test_project",
        )

        # Test data for creating a tag with associated projects
        tag_data = {
            "text": "Multi-Project Tag",
            "parent": None,
            "projects": [self.project.id, second_project.id],
        }

        # Make POST request to create the tag
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(tag_data),
        )

        # Verify the response
        assert response.status_code == 201

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        assert created_tag["text"] == "Multi-Project Tag"
        assert set(created_tag["projects"]) == {self.project.id, second_project.id}

        # Verify tag exists in database with correct project associations
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        assert db_tag.text == "Multi-Project Tag"
        assert set(db_tag.projects.values_list("id", flat=True)) == {
            self.project.id,
            second_project.id,
        }

    def test_create_duplicate_project_tag_fails(self):
        """Test that creating a duplicate tag (same text and parent) fails"""
        # Create initial tag
        models.ProjectTag.objects.create(text="Duplicate Tag", parent=None)

        # Try to create duplicate tag
        duplicate_tag_data = {"text": "Duplicate Tag", "parent": None, "projects": []}

        # Make POST request to create the duplicate tag
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(duplicate_tag_data),
        )

        # Verify the response fails due to unique constraint
        assert response.status_code == 400  # HTTP 400 Bad Request

        # Verify only one tag exists with this text
        tag_count = models.ProjectTag.objects.filter(
            text="Duplicate Tag", parent=None
        ).count()
        assert tag_count == 1

    def test_create_project_tag_missing_text_fails(self):
        """Test that creating a tag without required 'text' field fails"""
        # Test data missing required 'text' field
        invalid_tag_data = {"parent": None, "projects": []}

        # Make POST request with invalid data
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        assert response.status_code == 400  # HTTP 400 Bad Request

    def test_create_project_tag_empty_text_fails(self):
        """Test that creating a tag with empty 'text' field fails"""
        # Test data with empty text
        invalid_tag_data = {"text": "", "parent": None, "projects": []}

        # Make POST request with invalid data
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        assert response.status_code == 400  # HTTP 400 Bad Request

    def test_create_project_tag_too_long_text_fails(self):
        """Test that creating a tag with text longer than max_length fails"""
        # Test data with text exceeding CharField max_length=50
        invalid_tag_data = {
            "text": "A" * 51,  # 51 characters, exceeds max_length=50
            "parent": None,
            "projects": [],
        }

        # Make POST request with invalid data
        response = self.client.post(
            f"{API_PREFIX}/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        assert response.status_code == 400  # HTTP 400 Bad Request

    def test_get_project_tags(self):
        """Test retrieving all ProjectTags via GET API"""
        # Create some test tags
        tag1 = models.ProjectTag.objects.create(text="Tag One", parent=None)
        tag2 = models.ProjectTag.objects.create(text="Tag Two", parent=None)
        tag3 = models.ProjectTag.objects.create(text="Child Tag", parent=tag1)

        # Make GET request to retrieve all tags
        response = self.client.get(f"{API_PREFIX}/projecttags/")

        # Verify the response
        assert response.status_code == 200

        # Parse response data
        tags_list = response.json()

        # Verify we get all created tags (plus the one from setUp)
        assert len(tags_list) >= 3

        # Verify tag data is correct
        tag_texts = [tag["text"] for tag in tags_list]
        assert tag1.text in tag_texts
        assert tag2.text in tag_texts
        assert tag3.text in tag_texts

    def test_get_single_project_tag(self):
        """Test retrieving a single ProjectTag by ID"""
        # Create a test tag
        tag = models.ProjectTag.objects.create(text="Single Tag", parent=None)

        # Make GET request to retrieve specific tag
        response = self.client.get(f"{API_PREFIX}/projecttags/{tag.id}/")

        # Verify the response
        assert response.status_code == 200

        # Parse response data
        tag_data = response.json()

        # Verify tag data is correct
        assert tag_data["id"] == tag.id
        assert tag_data["text"] == "Single Tag"
        assert tag_data["parent"] is None

    def test_add_tag_to_project(self):
        """Test adding an existing tag to a project via the project/tags endpoint"""
        # Create a test tag
        tag = models.ProjectTag.objects.create(
            text="Project Association Tag", parent=None
        )

        # Verify tag is not associated with project initially
        assert self.project.tags.count() == 0

        # Test data for associating tag with project
        association_data = {"tag_id": tag.id}

        # Make POST request to associate tag with project
        response = self.client.post(
            f"{API_PREFIX}/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response
        assert response.status_code == 200

        # Parse response data
        result = response.json()
        assert result["status"] == "success"
        assert "added to project" in result["message"]

        # Verify tag is now associated with project
        assert self.project.tags.count() == 1
        assert self.project.tags.first().id == tag.id

    def test_get_project_tags_endpoint(self):
        """Test retrieving tags for a specific project via project/tags endpoint"""
        # Create test tags and associate with project
        tag1 = models.ProjectTag.objects.create(text="Project Tag 1", parent=None)
        tag2 = models.ProjectTag.objects.create(text="Project Tag 2", parent=None)
        # Create an unassociated tag to verify it's not returned
        models.ProjectTag.objects.create(text="Unassociated Tag", parent=None)

        # Associate only tag1 and tag2 with the project
        self.project.tags.add(tag1, tag2)

        # Make GET request to retrieve project tags
        response = self.client.get(f"{API_PREFIX}/projects/{self.project.id}/tags/")

        # Verify the response
        assert response.status_code == 200

        # Parse response data
        project_tags = response.json()

        # Verify we get only the associated tags
        assert len(project_tags) == 2

        tag_texts = [tag["text"] for tag in project_tags]
        assert "Project Tag 1" in tag_texts
        assert "Project Tag 2" in tag_texts
        assert "Unassociated Tag" not in tag_texts

    def test_remove_tag_from_project(self):
        """Test removing a tag from a project via DELETE endpoint"""
        # Create a test tag and associate with project
        tag = models.ProjectTag.objects.create(text="Tag to Remove", parent=None)
        self.project.tags.add(tag)

        # Verify tag is associated with project
        assert self.project.tags.count() == 1

        # Make DELETE request to remove tag from project
        response = self.client.delete(
            f"{API_PREFIX}/projects/{self.project.id}/tags/{tag.id}/"
        )

        # Verify the response
        assert response.status_code == 200

        # Parse response data
        result = response.json()
        assert result["status"] == "success"
        assert "removed from project" in result["message"]

        # Verify tag is no longer associated with project
        assert self.project.tags.count() == 0

        # Verify the tag itself still exists (not deleted)
        assert models.ProjectTag.objects.filter(id=tag.id).exists()

    def test_add_nonexistent_tag_to_project_fails(self):
        """Test that adding a non-existent tag to a project fails"""
        # Test data with non-existent tag ID
        association_data = {"tag_id": 99999}  # Non-existent ID

        # Make POST request with invalid tag ID
        response = self.client.post(
            f"{API_PREFIX}/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response fails
        assert response.status_code == 404  # HTTP 404 Not Found

        # Parse response data
        result = response.json()
        assert result["error"] == "Tag not found"

    def test_add_tag_without_tag_id_fails(self):
        """Test that adding a tag without providing tag_id fails"""
        # Test data missing required tag_id
        association_data = {}

        # Make POST request with missing tag_id
        response = self.client.post(
            f"{API_PREFIX}/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response fails
        assert response.status_code == 400  # HTTP 400 Bad Request

        # Parse response data
        result = response.json()
        assert result["error"] == "tag_id is required"
