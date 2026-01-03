import json
import logging
from pathlib import Path
from shutil import rmtree
from django.test import TestCase, override_settings
from django.conf import settings
from ...db import models
from ...db.models import Project, ProjectTag
from ...db.import_i2xml import import_ccp4_project_zip

logger = logging.getLogger(f"ccp4i2::{__name__}")


@override_settings(
    CCP4I2_PROJECTS_DIR=Path(__file__).parent.parent
    / "CCP4I2_TEST_PROJECT_TAG_DIRECTORY"
)
class ProjectTagAPITestCase(TestCase):
    def setUp(self):
        """Set up test fixtures."""
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
        self.projecttags_url = "/projecttags/"
        self.project_tags_url = f"/projects/{self.project.uuid}/tags/"

    def tearDown(self):
        """Clean up test environment"""
        if Path(settings.CCP4I2_PROJECTS_DIR).exists():
            rmtree(settings.CCP4I2_PROJECTS_DIR)
        return super().tearDown()

    def test_create_project_tag(self):
        """Test creating a new ProjectTag via API"""
        # Test data for creating a new tag
        tag_data = {"text": "Test Tag", "parent": None, "projects": []}

        # Make POST request to create the tag
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(tag_data),
        )

        # Verify the response
        self.assertEqual(response.status_code, 201)  # HTTP 201 Created

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        self.assertEqual(created_tag["text"], "Test Tag")
        self.assertIsNone(created_tag["parent"])
        self.assertEqual(created_tag["projects"], [])
        self.assertIsNotNone(created_tag["id"])

        # Verify tag exists in database
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        self.assertEqual(db_tag.text, "Test Tag")
        self.assertIsNone(db_tag.parent)

    def test_create_project_tag_with_parent(self):
        """Test creating a hierarchical ProjectTag with parent"""
        # First create a parent tag
        parent_tag = models.ProjectTag.objects.create(text="Parent Tag", parent=None)

        # Test data for creating a child tag
        child_tag_data = {"text": "Child Tag", "parent": parent_tag.id, "projects": []}

        # Make POST request to create the child tag
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(child_tag_data),
        )

        # Verify the response
        self.assertEqual(response.status_code, 201)

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        self.assertEqual(created_tag["text"], "Child Tag")
        self.assertEqual(created_tag["parent"], parent_tag.id)

        # Verify tag exists in database with correct parent
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        self.assertEqual(db_tag.text, "Child Tag")
        self.assertEqual(db_tag.parent.id, parent_tag.id)

    def test_create_project_tag_with_projects(self):
        """Test creating a ProjectTag and associating it with projects"""
        # Create a second test project
        second_project = models.Project.objects.create(
            name="Second Test Project",
            description="Another test project",
            directory=str(settings.CCP4I2_PROJECTS_DIR / "second_test_project"),
        )

        # Test data for creating a tag with associated projects
        tag_data = {
            "text": "Multi-Project Tag",
            "parent": None,
            "projects": [self.project.id, second_project.id],
        }

        # Make POST request to create the tag
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(tag_data),
        )

        # Verify the response
        self.assertEqual(response.status_code, 201)

        # Parse response data
        created_tag = response.json()

        # Verify the created tag has correct data
        self.assertEqual(created_tag["text"], "Multi-Project Tag")
        self.assertCountEqual(
            created_tag["projects"], [self.project.id, second_project.id]
        )

        # Verify tag exists in database with correct project associations
        db_tag = models.ProjectTag.objects.get(id=created_tag["id"])
        self.assertEqual(db_tag.text, "Multi-Project Tag")
        self.assertCountEqual(
            list(db_tag.projects.values_list("id", flat=True)),
            [self.project.id, second_project.id],
        )

    def test_create_duplicate_project_tag_fails(self):
        """Test that creating a duplicate tag (same text and parent) fails"""
        # Create initial tag
        models.ProjectTag.objects.create(text="Duplicate Tag", parent=None)

        # Try to create duplicate tag
        duplicate_tag_data = {"text": "Duplicate Tag", "parent": None, "projects": []}

        # Make POST request to create the duplicate tag
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(duplicate_tag_data),
        )

        # Verify the response fails due to unique constraint
        self.assertEqual(response.status_code, 400)  # HTTP 400 Bad Request

        # Verify only one tag exists with this text
        tag_count = models.ProjectTag.objects.filter(
            text="Duplicate Tag", parent=None
        ).count()
        self.assertEqual(tag_count, 1)

    def test_create_project_tag_missing_text_fails(self):
        """Test that creating a tag without required 'text' field fails"""
        # Test data missing required 'text' field
        invalid_tag_data = {"parent": None, "projects": []}

        # Make POST request with invalid data
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        self.assertEqual(response.status_code, 400)  # HTTP 400 Bad Request

    def test_create_project_tag_empty_text_fails(self):
        """Test that creating a tag with empty 'text' field fails"""
        # Test data with empty text
        invalid_tag_data = {"text": "", "parent": None, "projects": []}

        # Make POST request with invalid data
        response = self.client.post(
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        self.assertEqual(response.status_code, 400)  # HTTP 400 Bad Request

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
            "/projecttags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(invalid_tag_data),
        )

        # Verify the response fails
        self.assertEqual(response.status_code, 400)  # HTTP 400 Bad Request

    def test_get_project_tags(self):
        """Test retrieving all ProjectTags via GET API"""
        # Create some test tags
        tag1 = models.ProjectTag.objects.create(text="Tag One", parent=None)
        tag2 = models.ProjectTag.objects.create(text="Tag Two", parent=None)
        tag3 = models.ProjectTag.objects.create(text="Child Tag", parent=tag1)

        # Make GET request to retrieve all tags
        response = self.client.get("/projecttags/")

        # Verify the response
        self.assertEqual(response.status_code, 200)

        # Parse response data
        tags_list = response.json()

        # Verify we get all created tags (plus the one from setUp)
        self.assertGreaterEqual(len(tags_list), 3)

        # Verify tag data is correct
        tag_texts = [tag["text"] for tag in tags_list]
        self.assertIn(tag1.text, tag_texts)
        self.assertIn(tag2.text, tag_texts)
        self.assertIn(tag3.text, tag_texts)

    def test_get_single_project_tag(self):
        """Test retrieving a single ProjectTag by ID"""
        # Create a test tag
        tag = models.ProjectTag.objects.create(text="Single Tag", parent=None)

        # Make GET request to retrieve specific tag
        response = self.client.get(f"/projecttags/{tag.id}/")

        # Verify the response
        self.assertEqual(response.status_code, 200)

        # Parse response data
        tag_data = response.json()

        # Verify tag data is correct
        self.assertEqual(tag_data["id"], tag.id)
        self.assertEqual(tag_data["text"], "Single Tag")
        self.assertIsNone(tag_data["parent"])

    def test_add_tag_to_project(self):
        """Test adding an existing tag to a project via the project/tags endpoint"""
        # Create a test tag
        tag = models.ProjectTag.objects.create(
            text="Project Association Tag", parent=None
        )

        # Verify tag is not associated with project initially
        self.assertEqual(self.project.tags.count(), 0)

        # Test data for associating tag with project
        association_data = {"tag_id": tag.id}

        # Make POST request to associate tag with project
        response = self.client.post(
            f"/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response
        self.assertEqual(response.status_code, 200)

        # Parse response data
        result = response.json()
        self.assertEqual(result["status"], "success")
        self.assertIn("added to project", result["message"])

        # Verify tag is now associated with project
        self.assertEqual(self.project.tags.count(), 1)
        self.assertEqual(self.project.tags.first().id, tag.id)

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
        response = self.client.get(f"/projects/{self.project.id}/tags/")

        # Verify the response
        self.assertEqual(response.status_code, 200)

        # Parse response data
        project_tags = response.json()

        # Verify we get only the associated tags
        self.assertEqual(len(project_tags), 2)

        tag_texts = [tag["text"] for tag in project_tags]
        self.assertIn("Project Tag 1", tag_texts)
        self.assertIn("Project Tag 2", tag_texts)
        self.assertNotIn("Unassociated Tag", tag_texts)

    def test_remove_tag_from_project(self):
        """Test removing a tag from a project via DELETE endpoint"""
        # Create a test tag and associate with project
        tag = models.ProjectTag.objects.create(text="Tag to Remove", parent=None)
        self.project.tags.add(tag)

        # Verify tag is associated with project
        self.assertEqual(self.project.tags.count(), 1)

        # Make DELETE request to remove tag from project
        response = self.client.delete(f"/projects/{self.project.id}/tags/{tag.id}/")

        # Verify the response
        self.assertEqual(response.status_code, 200)

        # Parse response data
        result = response.json()
        self.assertEqual(result["status"], "success")
        self.assertIn("removed from project", result["message"])

        # Verify tag is no longer associated with project
        self.assertEqual(self.project.tags.count(), 0)

        # Verify the tag itself still exists (not deleted)
        self.assertTrue(models.ProjectTag.objects.filter(id=tag.id).exists())

    def test_add_nonexistent_tag_to_project_fails(self):
        """Test that adding a non-existent tag to a project fails"""
        # Test data with non-existent tag ID
        association_data = {"tag_id": 99999}  # Non-existent ID

        # Make POST request with invalid tag ID
        response = self.client.post(
            f"/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response fails
        self.assertEqual(response.status_code, 404)  # HTTP 404 Not Found

        # Parse response data
        result = response.json()
        self.assertEqual(result["error"], "Tag not found")

    def test_add_tag_without_tag_id_fails(self):
        """Test that adding a tag without providing tag_id fails"""
        # Test data missing required tag_id
        association_data = {}

        # Make POST request with missing tag_id
        response = self.client.post(
            f"/projects/{self.project.id}/tags/",
            content_type="application/json; charset=utf-8",
            data=json.dumps(association_data),
        )

        # Verify the response fails
        self.assertEqual(response.status_code, 400)  # HTTP 400 Bad Request

        # Parse response data
        result = response.json()
        self.assertEqual(result["error"], "tag_id is required")
