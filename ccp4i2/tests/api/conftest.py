"""
Pytest configuration for API tests that require subprocess job execution.

This conftest configures a file-based database that both the test
and spawned subprocesses can access.
"""
import os
from pathlib import Path
from shutil import rmtree

import pytest


TEST_DIR = Path(__file__).parent / "TEST_AIMLESS_API_PROJECTS"
TEST_DB = TEST_DIR / "test_aimless.db"


@pytest.fixture(scope="session")
def file_based_db(django_db_blocker):
    """Set up a file-based database manually, bypassing pytest-django.

    This is needed because pytest-django creates its own test database
    which subprocesses cannot access.
    """
    from django.conf import settings
    from django.core.management import call_command
    from django import db

    # Clean slate
    if TEST_DIR.exists():
        rmtree(TEST_DIR, ignore_errors=True)
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    # Override settings to use file-based database
    settings.DATABASES["default"]["NAME"] = str(TEST_DB)

    # Set env var for subprocess
    os.environ["CCP4I2_DB_FILE"] = str(TEST_DB)

    # Unblock database access for setup and tests
    with django_db_blocker.unblock():
        db.connections.close_all()
        call_command('migrate', '--run-syncdb', verbosity=0)

        yield TEST_DB

    # Teardown
    db.connections.close_all()
    os.environ.pop("CCP4I2_DB_FILE", None)
    if TEST_DIR.exists():
        rmtree(TEST_DIR, ignore_errors=True)
