"""
Pytest configuration for API tests.

This conftest integrates with the same test infrastructure as tests/i2run:
- Uses test_config.py for project directory management
- Uses the same CCP4 environment setup via run_test.sh
- Configures file-based database for subprocess job execution

Usage:
    ./run_test.sh tests/api/ -v
    ./run_test.sh tests/api/test_data_reduction_api.py -v
"""
import os
import gc
from pathlib import Path

import pytest
import django

from ccp4i2 import I2_TOP

# Configure Django settings
# Force use of test_settings even if run_test.sh set ccp4i2.config.settings
os.environ["DJANGO_SETTINGS_MODULE"] = "ccp4i2.config.test_settings"

# Import test configuration utilities (shared with i2run tests)
from ccp4i2.tests.i2run.test_config import get_test_projects_dir, make_test_project_name, get_cleanup_message

# Set up test projects directory (now defaults to ~/.cache/ccp4i2-tests/)
TEST_PROJECTS_DIR = get_test_projects_dir()
os.environ["CCP4I2_PROJECTS_DIR"] = str(TEST_PROJECTS_DIR)

# Initialize Django
django.setup()

# === Module-level permission bypass for Django TestCase tests ===
# This runs at import time to ensure permissions are bypassed for
# Django TestCase classes which don't use pytest fixtures
def _bypass_permissions_globally():
    """Patch all viewsets to allow unauthenticated access."""
    from rest_framework.permissions import AllowAny
    from ccp4i2.api.ProjectViewSet import ProjectViewSet
    from ccp4i2.api.JobViewSet import JobViewSet
    from ccp4i2.api.FileViewSet import FileViewSet
    from ccp4i2.api.ProjectTagViewSet import ProjectTagViewSet
    from ccp4i2.api.ProjectExportViewSet import ProjectExportViewSet
    from ccp4i2.api.ProjectGroupViewSet import ProjectGroupViewSet
    from ccp4i2.api.ProjectGroupMembershipViewSet import ProjectGroupMembershipViewSet
    from ccp4i2.api.FileTypeViewSet import FileTypeViewSet
    from ccp4i2.api.FileUseViewSet import FileUseViewSet
    from ccp4i2.api.FileImportViewSet import FileImportViewSet

    ProjectViewSet.permission_classes = [AllowAny]
    JobViewSet.permission_classes = [AllowAny]
    FileViewSet.permission_classes = [AllowAny]
    ProjectTagViewSet.permission_classes = [AllowAny]
    ProjectExportViewSet.permission_classes = [AllowAny]
    ProjectGroupViewSet.permission_classes = [AllowAny]
    ProjectGroupMembershipViewSet.permission_classes = [AllowAny]
    FileTypeViewSet.permission_classes = [AllowAny]
    FileUseViewSet.permission_classes = [AllowAny]
    FileImportViewSet.permission_classes = [AllowAny]

_bypass_permissions_globally()

# Import URL helpers from i2run
from ccp4i2.tests.i2run.urls import pdbe_fasta, redo_cif, redo_mtz, rcsb_mmcif
from ccp4i2.tests.i2run.utils import download


def pytest_collection_modifyitems(items):
    """Automatically add django_db marker to all test items."""
    for item in items:
        item.add_marker(pytest.mark.django_db(transaction=True))


@pytest.fixture(autouse=True)
def bypass_api_permissions(monkeypatch):
    """
    Bypass API permissions for tests.

    This fixture patches the permission_classes on all viewsets to allow
    unauthenticated access during tests.
    """
    from rest_framework.permissions import AllowAny

    # Import viewsets after Django is set up
    from ccp4i2.api.ProjectViewSet import ProjectViewSet
    from ccp4i2.api.JobViewSet import JobViewSet
    from ccp4i2.api.FileViewSet import FileViewSet
    from ccp4i2.api.ProjectTagViewSet import ProjectTagViewSet
    from ccp4i2.api.ProjectExportViewSet import ProjectExportViewSet
    from ccp4i2.api.ProjectGroupViewSet import ProjectGroupViewSet
    from ccp4i2.api.ProjectGroupMembershipViewSet import ProjectGroupMembershipViewSet
    from ccp4i2.api.FileTypeViewSet import FileTypeViewSet
    from ccp4i2.api.FileUseViewSet import FileUseViewSet
    from ccp4i2.api.FileImportViewSet import FileImportViewSet

    # Patch permission_classes on all viewsets
    monkeypatch.setattr(ProjectViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(JobViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(FileViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(ProjectTagViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(ProjectExportViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(ProjectGroupViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(ProjectGroupMembershipViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(FileTypeViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(FileUseViewSet, 'permission_classes', [AllowAny])
    monkeypatch.setattr(FileImportViewSet, 'permission_classes', [AllowAny])


@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item):
    """Hook to capture test results for cleanup decisions."""
    outcome = yield
    rep = outcome.get_result()
    setattr(item, f"rep_{rep.when}", rep)


def pytest_sessionfinish(session, exitstatus):
    """Called after whole test run finished.

    Prints cleanup instructions for the test cache directory.
    """
    # Only print if we're the main process (not a worker in parallel runs)
    if hasattr(session.config, 'workerinput'):
        return  # Skip in xdist worker processes

    # Print cleanup message
    print(get_cleanup_message())


@pytest.fixture(scope="session")
def test_projects_dir():
    """Create and return test projects directory."""
    TEST_PROJECTS_DIR.mkdir(exist_ok=True)
    yield TEST_PROJECTS_DIR
    # Don't clean up - leave for inspection


@pytest.fixture(scope="session", autouse=True)
def django_db_setup():
    """Set up test database directory for Django."""
    from django.conf import settings

    # Ensure test projects directory exists
    TEST_PROJECTS_DIR.mkdir(exist_ok=True)

    # Set CCP4I2_PROJECTS_DIR in settings to match our test directory
    settings.CCP4I2_PROJECTS_DIR = TEST_PROJECTS_DIR

    yield

    # Session cleanup: remove old test databases (keep most recent 10)
    old_dbs = sorted(TEST_PROJECTS_DIR.glob("test_*.sqlite*"), key=lambda p: p.stat().st_mtime, reverse=True)
    for old_db in old_dbs[10:]:  # Keep newest 10
        try:
            old_db.unlink()
        except Exception as e:
            print(f"Warning: Failed to clean up old database {old_db}: {e}")


@pytest.fixture
def test_project_path(request):
    """
    Provide the test project directory path to test functions.

    Project naming format: YYYYMMDD_HHMMSS_XXXX_api_<test_file>_<test_function>
    """
    # Extract test function name from node
    test_function = request.node.name

    # Extract test file name from node's fspath
    test_file_path = Path(str(request.node.fspath))
    test_file = None
    if test_file_path.stem.startswith('test_'):
        test_file = test_file_path.stem[5:]  # Remove 'test_' prefix

    # Build base name for the project
    if test_function and test_file:
        base_name = f"api_{test_file}_{test_function}"
    elif test_function:
        base_name = f"api_{test_function}"
    else:
        base_name = request.node.nodeid.replace('/', '_').replace('::', '_')

    # Use make_test_project_name to add timestamp and random suffix
    project_name = make_test_project_name(base_name)
    test_project_dir = TEST_PROJECTS_DIR / project_name

    # Store on request for use by isolated_test_db
    request.node._test_project_dir = test_project_dir

    # Also set as environment variable
    os.environ['_PYTEST_PROJECT_DIR'] = str(test_project_dir)

    return test_project_dir


@pytest.fixture(autouse=True)
def isolated_test_db(request, django_db_blocker, monkeypatch, test_project_path):
    """Create isolated database for each test with proper Django connection management."""
    import shutil
    from django.core.management import call_command
    from django.conf import settings
    from django.db import connections

    # Monitor file descriptors at test start (for debugging resource leaks)
    try:
        fd_start = len(os.listdir('/dev/fd'))
    except Exception:
        fd_start = None

    # Use the project directory from test_project_path fixture
    test_project_dir = test_project_path

    # Create project directory
    test_project_dir.mkdir(exist_ok=True)

    # Place SQLite database inside the project directory
    test_db_path = test_project_dir / "project.sqlite"

    # Store original database configuration
    original_db_settings = settings.DATABASES['default'].copy()

    # Update database configuration for this test
    # IMPORTANT: Update settings BEFORE closing connections so Django reconnects to new path
    monkeypatch.setitem(settings.DATABASES['default'], 'NAME', str(test_db_path))

    # Set environment variables for subprocesses
    os.environ["CCP4I2_DB_FILE"] = str(test_db_path)
    os.environ["CCP4I2_PROJECTS_DIR"] = str(test_project_dir)

    # Also update Django settings directly (env var may not be picked up after startup)
    monkeypatch.setattr(settings, 'CCP4I2_PROJECTS_DIR', test_project_dir)

    # Force Django to use the new database by reconfiguring the connection
    # First, close all existing connections
    connections.close_all()

    # Reset the connection handler to force it to re-read settings
    # This is more aggressive than close_all() which just closes the connection
    # but keeps the DatabaseWrapper object with its cached settings_dict
    connections._databases = settings.DATABASES

    # Also update the default connection's settings_dict directly
    # Django caches this when the connection is first created
    if hasattr(connections._connections, 'default'):
        del connections._connections.default

    # Unblock database access for migrations
    with django_db_blocker.unblock():
        # Ensure the database file exists and is ready
        test_db_path.touch()

        # Create and migrate the test database using explicit database option
        call_command('migrate', '--run-syncdb', '--database=default', verbosity=0)

        # Enable WAL mode for better concurrent access between test process and subprocess
        # WAL mode allows multiple readers and one writer simultaneously
        from django.db import connection
        with connection.cursor() as cursor:
            cursor.execute("PRAGMA journal_mode=WAL;")
            cursor.execute("PRAGMA synchronous=NORMAL;")
            # Disable timeout to avoid lock errors (default is 5 seconds)
            cursor.execute("PRAGMA busy_timeout=30000;")  # 30 seconds

            # Verify tables were created
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='ccp4i2_job';")
            result = cursor.fetchone()
            if not result:
                raise RuntimeError(f"Migration failed - ccp4i2_job table not created in {test_db_path}")

    # Run the test
    yield

    # === AGGRESSIVE RESOURCE CLEANUP ===
    gc.collect()
    gc.collect()  # Run twice to catch circular references

    # Close all Django database connections
    with django_db_blocker.unblock():
        connections.close_all()

    # Restore database configuration
    settings.DATABASES['default'] = original_db_settings

    # Clean up environment variables
    os.environ.pop("CCP4I2_DB_FILE", None)

    # Monitor file descriptor leaks
    if fd_start is not None:
        try:
            fd_end = len(os.listdir('/dev/fd'))
            fd_leaked = fd_end - fd_start
            if fd_leaked > 5:  # Allow small variations
                test_name = request.node.name
                print(f"⚠️  File descriptor leak: {fd_leaked} FDs in {test_name}")
        except Exception:
            pass

    # Only remove the project directory if the test passed
    test_failed = request.node.rep_call.failed if hasattr(request.node, 'rep_call') else False

    if not test_failed:
        try:
            shutil.rmtree(test_project_dir, ignore_errors=True)
        except Exception as e:
            print(f"Warning: Failed to clean up test project directory {test_project_dir}: {e}")
    else:
        print(f"Test failed - preserving project directory: {test_project_dir}")


# === Fixture: file_based_db (backwards compatibility) ===
# This fixture is used by APITestBase for subprocess-accessible database

@pytest.fixture(scope="function")
def file_based_db(isolated_test_db, test_project_path):
    """
    Provide the test database path for API tests.

    This fixture wraps isolated_test_db to provide backwards compatibility
    with APITestBase which expects a file_based_db fixture.
    """
    return test_project_path / "project.sqlite"


# === Session-scoped fixtures for demo data ===

@pytest.fixture(scope="session")
def demo_data_dir():
    """Return path to demo_data directory."""
    return I2_TOP / 'demo_data'


@pytest.fixture(scope="session")
def gamma_mtz(demo_data_dir):
    """Path to gamma merged intensities MTZ."""
    path = demo_data_dir / 'gamma' / 'merged_intensities_Xe.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_native_mtz(demo_data_dir):
    """Path to gamma native merged intensities MTZ."""
    path = demo_data_dir / 'gamma' / 'merged_intensities_native.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_unmerged_mtz(demo_data_dir):
    """Path to gamma native unmerged MTZ."""
    path = demo_data_dir / 'gamma' / 'gamma_native.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_hklout_unmerged(demo_data_dir):
    """Path to gamma HKLOUT unmerged MTZ."""
    path = demo_data_dir / 'gamma' / 'HKLOUT_unmerged.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_model_pdb(demo_data_dir):
    """Path to gamma model PDB."""
    path = demo_data_dir / 'gamma' / 'gamma_model.pdb'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_freer_mtz(demo_data_dir):
    """Path to gamma freeR MTZ."""
    path = demo_data_dir / 'gamma' / 'freeR.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_phases_mtz(demo_data_dir):
    """Path to gamma initial phases MTZ."""
    path = demo_data_dir / 'gamma' / 'initial_phases.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_asu_xml(demo_data_dir):
    """Path to gamma ASU XML."""
    path = demo_data_dir / 'gamma' / 'gamma.asu.xml'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def gamma_heavy_atoms_pdb(demo_data_dir):
    """Path to gamma heavy atoms PDB."""
    path = demo_data_dir / 'gamma' / 'heavy_atoms.pdb'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def mdm2_unmerged_mtz(demo_data_dir):
    """Path to mdm2 unmerged MTZ."""
    path = demo_data_dir / 'mdm2' / 'mdm2_unmerged.mtz'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


@pytest.fixture(scope="session")
def rnase_model_pdb(demo_data_dir):
    """Path to rnase model PDB."""
    path = demo_data_dir / 'rnase' / 'rnase_model.pdb'
    if not path.exists():
        pytest.skip(f"Demo data not found: {path}")
    return str(path)


# === Downloaded data fixtures ===
# These use the download helper from i2run.utils

@pytest.fixture(scope="session")
def cif8xfm():
    """Download 8xfm CIF from PDB-REDO."""
    with download(redo_cif("8xfm")) as path:
        yield path


@pytest.fixture(scope="session")
def mtz8xfm():
    """Download 8xfm MTZ from PDB-REDO."""
    with download(redo_mtz("8xfm")) as path:
        yield path


@pytest.fixture(scope="session")
def seq8xfm():
    """Download 8xfm sequence from PDBe."""
    with download(pdbe_fasta("8xfm")) as path:
        yield path


@pytest.fixture(scope="session")
def cif7beq():
    """Download 7beq CIF from RCSB."""
    with download(rcsb_mmcif("7beq")) as path:
        yield path


@pytest.fixture(scope="session")
def mtz7beq():
    """Download 7beq MTZ from PDB-REDO."""
    with download(redo_mtz("7beq")) as path:
        yield path
