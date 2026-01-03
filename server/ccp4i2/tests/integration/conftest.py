import os
import django
import pytest
from pathlib import Path
from pytest import fixture

# Configure Django settings
# Force use of test_settings even if run_test.sh set ccp4i2.config.settings
os.environ["DJANGO_SETTINGS_MODULE"] = "ccp4i2.config.test_settings"

# Set up test projects directory
TEST_PROJECTS_DIR = Path(__file__).parent / "test_projects"
os.environ["CCP4I2_PROJECTS_DIR"] = str(TEST_PROJECTS_DIR)

# Note: CCP4 environment should be set up by run_test.sh before running tests
# Initialize Django
django.setup()

from ccp4i2.tests.i2run.urls import pdbe_fasta, redo_cif, redo_mtz, rcsb_mmcif
from ccp4i2.tests.i2run.utils import download


def pytest_collection_modifyitems(items):
    """Automatically add django_db marker to all test items."""
    for item in items:
        item.add_marker(pytest.mark.django_db(transaction=True))


@fixture(scope="session")
def test_projects_dir():
    """Create and return test projects directory."""
    TEST_PROJECTS_DIR.mkdir(exist_ok=True)
    yield TEST_PROJECTS_DIR
    # Don't clean up - leave for inspection


@fixture(scope="session", autouse=True)
def django_db_setup(django_db_blocker):
    """Set up test database for Django. Auto-use ensures it runs for all tests."""
    from django.core.management import call_command
    from django.conf import settings

    # Use single test database for sequential test execution
    test_db_path = TEST_PROJECTS_DIR / "test_ccp4i2.sqlite"
    settings.DATABASES['default']['NAME'] = str(test_db_path)

    # Also set CCP4I2_PROJECTS_DIR in settings to match our test directory
    settings.CCP4I2_PROJECTS_DIR = TEST_PROJECTS_DIR

    # Unblock database access for migrations
    with django_db_blocker.unblock():
        # Create tables
        call_command('migrate', '--run-syncdb', verbosity=0)

    yield

    # Don't clean up database - leave for inspection


@fixture(autouse=True)
def cleanup_after_test(django_db_blocker):
    """Clean up database and temp directories after each test for proper isolation."""
    import shutil
    from django.core.management import call_command
    from django.conf import settings

    # Run the test
    yield

    # Clean up temp project directories created during the test
    # Keep only the most recent 5 for debugging
    temp_dirs = sorted(TEST_PROJECTS_DIR.glob("tmp_*"), key=lambda p: p.stat().st_mtime, reverse=True)
    for temp_dir in temp_dirs[5:]:  # Keep newest 5, delete older ones
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            print(f"Warning: Failed to clean up {temp_dir}: {e}")

    # Clean up database: flush all data but keep schema
    # This is faster than recreating the database from scratch
    with django_db_blocker.unblock():
        try:
            # Flush database (delete all data, keep schema)
            call_command('flush', '--no-input', verbosity=0)
        except Exception as e:
            print(f"Warning: Failed to flush database: {e}")


@fixture(scope="session")
def cif7beq():
    with download(rcsb_mmcif("7beq")) as path:
        yield path


@fixture(scope="session")
def cif8xfm():
    with download(redo_cif("8xfm")) as path:
        yield path


@fixture(scope="session")
def mtz8xfm():
    with download(redo_mtz("8xfm")) as path:
        yield path


@fixture(scope="session")
def mtz7beq():
    with download(redo_mtz("7beq")) as path:
        yield path


@fixture(scope="session")
def seq8xfm():
    with download(pdbe_fasta("8xfm")) as path:
        yield path
