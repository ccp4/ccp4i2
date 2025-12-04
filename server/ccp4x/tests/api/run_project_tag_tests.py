#!/usr/bin/env python3
"""
Test runner script for ProjectTag API tests

This script demonstrates how to run the ProjectTag API tests.
Run from the server directory with:

    python -m pytest ccp4x/tests/api/test_project_tag_api.py -v

Or using Django's test runner:

    python manage.py test ccp4x.tests.api.test_project_tag_api

"""


def main():
    """Main function to run tests"""
    print("ProjectTag API Test Suite")
    print("=" * 50)
    print()
    print("To run these tests, use one of the following commands:")
    print()
    print("1. Using pytest (recommended):")
    print("   cd /path/to/server")
    print("   python -m pytest ccp4x/tests/api/test_project_tag_api.py -v")
    print()
    print("2. Using Django test runner:")
    print("   cd /path/to/server")
    print("   python manage.py test ccp4x.tests.api.test_project_tag_api")
    print()
    print("3. Run specific test:")
    print(
        "   python manage.py test ccp4x.tests.api.test_project_tag_api.ProjectTagAPITestCase.test_create_project_tag"
    )
    print()
    print("Test Coverage:")
    print("- ✅ Create new ProjectTag")
    print("- ✅ Create ProjectTag with parent (hierarchical)")
    print("- ✅ Create ProjectTag with project associations")
    print("- ✅ Duplicate tag validation (unique constraint)")
    print("- ✅ Missing/empty text validation")
    print("- ✅ Text length validation (max 50 chars)")
    print("- ✅ GET all ProjectTags")
    print("- ✅ GET single ProjectTag by ID")
    print("- ✅ Add existing tag to project")
    print("- ✅ GET project-specific tags")
    print("- ✅ Remove tag from project")
    print("- ✅ Error handling (non-existent tag, missing tag_id)")
    print()


if __name__ == "__main__":
    main()
