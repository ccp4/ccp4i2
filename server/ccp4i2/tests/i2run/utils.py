from contextlib import contextmanager
from multiprocessing import Process
from os.path import basename, join
from os import environ
from pathlib import Path
from random import choice
from shutil import rmtree
from string import ascii_letters, digits
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
from urllib.request import urlopen
from xml.etree import ElementTree as ET
from subprocess import Popen, run, PIPE, CalledProcessError
import gemmi
import sys

# Import test configuration utilities for project naming
from ccp4i2.tests.i2run.test_config import make_test_project_name

# Try to import from legacy ccp4i2 for getCCP4I2Dir, but make it optional
try:
    from ccp4i2.core.CCP4Utils import getCCP4I2Dir
except ImportError:
    # Fallback: use demo_data from current project
    def getCCP4I2Dir():
        return str(Path(__file__).parent.parent.parent)


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a string path to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    urlName = unquote(basename(urlparse(url).path))
    with urlopen(url, timeout=30) as response:
        name = response.headers.get_filename() or urlName
        name = name.strip().replace(" ", "_")
        name = "".join(c for c in name if c.isalnum() or c in "-_.")
        with NamedTemporaryFile(suffix=f"_{name}", delete=False) as temp:
            while chunk := response.read(1_000_000):
                temp.write(chunk)
        path = Path(temp.name).resolve()
        try:
            yield str(path)
        finally:
            path.unlink(missing_ok=True)


@contextmanager
def i2run(args: list[str], project_name: str = None, project_path: Path = None, allow_errors: bool = False):
    """
    Run a task by calling Django management command directly (in-process).

    This runs in the same process as the test, preserving:
    - Django database setup and transactions
    - Test fixtures
    - Environment configuration

    Yields a Path object to the job directory (CCP4_JOBS/job_1).
    Use in a with statement and the project will be cleaned up
    as long as an error is not raised.

    Args:
        args: List of arguments starting with task name, followed by task parameters
        project_name: Optional project name (defaults to derived from test name).
                     Should reflect the test name for better organization.
        project_path: Optional explicit project directory path (overrides project_name).
                     When provided, CCP4_JOBS and other project files will be created here.
                     This allows conftest.py to pass a random-suffixed directory for test isolation.
        allow_errors: If True, skip the assertion that diagnostic.xml has no errors.
                     Use for tests that expect task failures (e.g., no solution found).
    """
    # Use explicit project_path if provided, otherwise try to get from environment
    if project_path is None:
        # Check if conftest.py set the project directory via environment variable
        # This allows i2run() to use the same random-suffixed directory created by the fixture
        pytest_project_dir = environ.get('_PYTEST_PROJECT_DIR')
        if pytest_project_dir:
            project_path = Path(pytest_project_dir)

    # If still no project_path, derive from project_name
    if project_path is None:
        # Generate project name if not provided
        if project_name is None:
            # Try to get the current test name from pytest's request context
            # Include test file name to avoid collisions between tests with same name in different files
            base_name = None
            try:
                import inspect
                # Walk up the stack to find the test function and file
                test_function = None
                test_file = None
                for frame_info in inspect.stack():
                    # Look for a frame with 'test_' in the function name
                    if frame_info.function.startswith('test_'):
                        test_function = frame_info.function
                        # Extract just the test file name without path or extension
                        # e.g., /path/to/test_phaser_simple.py -> phaser_simple
                        test_file_path = Path(frame_info.filename)
                        if test_file_path.stem.startswith('test_'):
                            test_file = test_file_path.stem[5:]  # Remove 'test_' prefix
                        break

                if test_function and test_file:
                    # Format: {file}_{function}
                    # e.g., phaser_simple_test_gamma
                    base_name = f"{test_file}_{test_function}"
                elif test_function:
                    # Fallback: just function name
                    base_name = test_function
            except Exception:
                pass

            # Final fallback: random base name
            if base_name is None:
                chars = ascii_letters + digits
                base_name = "unknown_" + "".join(choice(chars) for _ in range(6))
                base_name = base_name.lower()

            # Use make_test_project_name to add timestamp and random suffix
            # Format: YYYYMMDD_HHMMSS_XXXX_{base_name}
            project_name = make_test_project_name(base_name)

        # Check if we're running in test environment with CCP4I2_PROJECTS_DIR set
        projects_dir = environ.get("CCP4I2_PROJECTS_DIR")
        if projects_dir:
            # Test mode: place project in test_projects directory
            project_path = Path(projects_dir) / project_name
        else:
            # Normal mode: place project in current directory
            project_path = Path(project_name)
    else:
        # Explicit project_path provided - derive project_name from it if needed
        if project_name is None:
            project_name = project_path.name

    # Build command-line arguments for i2run management command
    # Format: manage.py i2run task_name --project_name foo --param1 val1 --param2 val2
    i2run_argv = ['manage.py', 'i2run', args[0]]  # args[0] is task name
    i2run_argv.extend(['--project_name', project_name])
    i2run_argv.extend(args[1:])  # Add plugin-specific parameters

    # Save original sys.argv
    original_sys_argv = sys.argv

    try:
        # Update sys.argv for the management command parser
        # The i2run command uses sys.argv[2:] directly, so we need to set it properly
        sys.argv = i2run_argv

        # Import and run the i2run command directly (same process, preserves test database)
        from django.core.management import call_command
        from io import StringIO

        # Capture stdout/stderr
        stdout_capture = StringIO()
        stderr_capture = StringIO()

        try:
            # Call the management command directly without additional arguments
            # The command will read sys.argv[2:] which we've already set
            call_command(
                'i2run',
                stdout=stdout_capture,
                stderr=stderr_capture
            )

            # Print captured output
            stdout_val = stdout_capture.getvalue()
            stderr_val = stderr_capture.getvalue()
            if stdout_val:
                print("=== i2run stdout ===")
                print(stdout_val)
            if stderr_val:
                print("=== i2run stderr ===")
                print(stderr_val)

        except Exception as e:
            print(f"Error running i2run: {e}")
            print(f"=== stdout ===")
            print(stdout_capture.getvalue())
            print(f"=== stderr ===")
            print(stderr_capture.getvalue())
            raise
        finally:
            # Explicitly close StringIO objects to free file descriptors
            stdout_capture.close()
            stderr_capture.close()

    finally:
        # Restore original sys.argv
        sys.argv = original_sys_argv

    # Find the actual job directory created (may be job_1, job_2, etc. due to database state)
    # In comprehensive tests, Django may reuse database with existing jobs, causing auto-increment
    jobs_dir = project_path / "CCP4_JOBS"
    if jobs_dir.exists():
        # Find the highest-numbered job directory (most recent)
        job_dirs = sorted([d for d in jobs_dir.iterdir() if d.is_dir() and d.name.startswith("job_")],
                         key=lambda d: int(d.name.split("_")[1]))
        if job_dirs:
            directory = job_dirs[-1]  # Use the most recent job directory
        else:
            # Fallback: assume job_1
            directory = jobs_dir / "job_1"
    else:
        directory = project_path / "CCP4_JOBS" / "job_1"

    # Debug: Show what was actually created
    if directory.exists():
        print(f"=== Job directory contents: {directory} ===")
        for item in directory.iterdir():
            print(f"  {item.name}")
    else:
        print(f"WARNING: Job directory does not exist: {directory}")

    # Check diagnostic.xml for errors (unless allow_errors is set)
    # Only fail on actual errors (severity >= 4), not warnings (severity 2)
    xml_path = directory / "diagnostic.xml"
    if xml_path.exists():
        all_reports = ET.parse(xml_path).findall(".//errorReport")
        # Filter to only errors (severity 4) not warnings (severity 2)
        errors = [r for r in all_reports if int(r.find("severity").text or "0") >= 4]
        warnings = [r for r in all_reports if int(r.find("severity").text or "0") < 4]

        if warnings:
            print(f"Note: {len(warnings)} warning(s) in diagnostic.xml (non-fatal)")

        if not allow_errors:
            if errors:
                error_details = []
                for e in errors:
                    name = e.find("name")
                    details = e.find("details")
                    error_details.append(f"  - {name.text if name is not None else 'unknown'}: {details.text if details is not None else 'no details'}")
                assert len(errors) == 0, f"Error reports found in diagnostic.xml:\n" + "\n".join(error_details)
        elif errors:
            print(f"Note: {len(errors)} error report(s) found in diagnostic.xml (expected, allow_errors=True)")
    else:
        print(f"Warning: diagnostic.xml not found at {xml_path}")

    # Use try-finally to ensure cleanup happens, but store exception
    error_occurred = False
    try:
        yield directory
    except Exception as e:
        error_occurred = True
        print(f"\n!!! Test failed - preserving project directory for inspection: {project_path}")
        raise
    finally:
        # Force garbage collection to release gemmi file handles
        # This helps prevent resource exhaustion in long test runs
        import gc
        gc.collect()
        gc.collect()  # Run twice to catch circular references

        # Only clean up if no error occurred
        # Note: Failed test directories are preserved in test_projects/tmp_{file}_{test}/
        # With consolidated structure, both CCP4_JOBS and database are in the same directory
        if not error_occurred:
            # Clean up project directory (includes both CCP4_JOBS and project.sqlite)
            rmtree(str(project_path), ignore_errors=True)
        else:
            # Test failed - directory already logged above
            # Directory contains both CCP4_JOBS/job_1 and project.sqlite for inspection
            pass


def demoData(*paths):
    return join(getCCP4I2Dir(), "demo_data", *paths)


def hasLongLigandName(path):
    "Does the structure contains a residue with a name longer than 3 characters?"
    structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Mmcif)
    for model in structure:
        for chain in model:
            for residue in chain:
                if len(residue.name) > 3:
                    return True
    return False
