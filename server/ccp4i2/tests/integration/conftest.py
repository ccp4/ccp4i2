import os
import sys
import django
import pytest
from pathlib import Path
from pytest import fixture

# Add server directory and project root to Python path
server_path = Path(__file__).parent.parent.parent / "server"
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(server_path))
sys.path.insert(0, str(project_root))

# Configure Django settings
# Force use of test_settings even if run_test.sh set ccp4i2.config.settings
os.environ["DJANGO_SETTINGS_MODULE"] = "ccp4i2.config.test_settings"

# Set up CCP4I2_ROOT for plugin discovery (.def.xml files)
if "CCP4I2_ROOT" not in os.environ:
    os.environ["CCP4I2_ROOT"] = str(project_root)

# Set up test projects directory
TEST_PROJECTS_DIR = Path(__file__).parent / "test_projects"
os.environ["CCP4I2_PROJECTS_DIR"] = str(TEST_PROJECTS_DIR)

# Source CCP4 environment if available
CCP4_SETUP_SCRIPT = "/Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh"
if Path(CCP4_SETUP_SCRIPT).exists():
    import subprocess
    # Source the setup script and export environment variables
    # Use bash to source the script and print all environment variables
    result = subprocess.run(
        f'source {CCP4_SETUP_SCRIPT} && env',
        shell=True,
        executable='/bin/bash',
        capture_output=True,
        text=True
    )
    if result.returncode == 0:
        # Parse the environment variables and add to os.environ
        ccp4_vars = {}
        for line in result.stdout.splitlines():
            if '=' in line:
                key, _, value = line.partition('=')
                # Only set CCP4-related variables and PATH to avoid polluting environment
                if key.startswith('CCP4') or key in ['CBIN', 'CLIB', 'CCP4_SCR', 'PATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']:
                    ccp4_vars[key] = value

        # Update os.environ with CCP4 variables
        for key, value in ccp4_vars.items():
            os.environ[key] = value

        print(f"CCP4 environment loaded from {CCP4_SETUP_SCRIPT}")
        print(f"  CCP4={os.environ.get('CCP4', 'NOT SET')}")
        print(f"  CBIN={os.environ.get('CBIN', 'NOT SET')}")
        print(f"  PATH includes CBIN: {os.environ.get('CBIN', '') in os.environ.get('PATH', '')}")
    else:
        print(f"Warning: Failed to source CCP4 setup script: {result.stderr}")
else:
    print(f"Warning: CCP4 setup script not found at {CCP4_SETUP_SCRIPT}")

# Initialize Django
django.setup()

from ..i2run.urls import pdbe_fasta, redo_cif, redo_mtz, rcsb_mmcif
from ..i2run.utils import download


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
