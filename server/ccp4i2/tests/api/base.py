"""
Base classes and utilities for API-driven pipeline tests.

This module provides:
- APITestBase: A base class for API-driven tests with common methods
- Helper functions for file uploads, parameter setting, and job execution
- Data fetching utilities (demo_data and web downloads)

Re-exports from i2run for convenience:
- download: Context manager for downloading files from URLs
- URLs: Class with static URL builder methods
"""
import json
import os
import sys
import time
from pathlib import Path

import pytest
from django.test import Client


# Add tests/ directory to path so i2run package is importable
# (tests/i2run/ has __init__.py making it a package named 'i2run')
_tests_dir = Path(__file__).parent.parent
if str(_tests_dir) not in sys.path:
    sys.path.insert(0, str(_tests_dir))

# Also add tests/i2run/ to path for direct imports (like test_config)
# This is needed because i2run/utils.py uses 'from test_config import ...'
_i2run_dir = _tests_dir / "i2run"
if str(_i2run_dir) not in sys.path:
    sys.path.insert(0, str(_i2run_dir))

# Re-export download from i2run.utils
from i2run.utils import download

# Re-export URL functions from i2run.urls and wrap in URLs class for backward compatibility
from i2run import urls as _urls


class URLs:
    """URL builders for common crystallographic data sources.

    This class wraps functions from i2run.urls for backward compatibility
    with tests that use URLs.pdbe_fasta() style calls.
    """

    @staticmethod
    def pdbe_fasta(code: str) -> str:
        return _urls.pdbe_fasta(code)

    @staticmethod
    def pdbe_mmcif(code: str) -> str:
        return _urls.pdbe_mmcif(code)

    @staticmethod
    def pdbe_pdb(code: str) -> str:
        return _urls.pdbe_pdb(code)

    @staticmethod
    def pdbe_sfcif(code: str) -> str:
        return _urls.pdbe_sfcif(code)

    @staticmethod
    def rcsb_ligand_cif(code: str) -> str:
        return _urls.rcsb_ligand_cif(code)

    @staticmethod
    def rcsb_ligand_sdf(code: str) -> str:
        return _urls.rcsb_ligand_sdf(code)

    @staticmethod
    def rcsb_mmcif(code: str) -> str:
        return _urls.rcsb_mmcif(code)

    @staticmethod
    def rcsb_pdb(code: str) -> str:
        return _urls.rcsb_pdb(code)

    @staticmethod
    def redo_cif(code: str) -> str:
        return _urls.redo_cif(code)

    @staticmethod
    def redo_mtz(code: str) -> str:
        return _urls.redo_mtz(code)


# Try to import gemmi for validation (optional)
try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False


def get_ccp4i2_root():
    """Get CCP4I2_ROOT directory."""
    return Path(os.environ.get('CCP4I2_ROOT', Path(__file__).parent.parent.parent))


def demo_data(*paths):
    """Construct path to demo_data file."""
    return str(get_ccp4i2_root() / "demo_data" / Path(*paths))


# Job status constants
class JobStatus:
    """Job status enum matching Django model."""
    PENDING = 0
    QUEUED = 1
    RUNNING = 2
    PAUSED = 3
    HELD = 4
    FINISHING = 5
    FINISHED = 6
    FAILED = 7
    UNSATISFACTORY = 8

    @classmethod
    def is_terminal(cls, status):
        return status in (cls.FINISHED, cls.FAILED, cls.UNSATISFACTORY)

    @classmethod
    def is_success(cls, status):
        return status == cls.FINISHED


class APITestBase:
    """
    Base class for API-driven pipeline tests.

    Provides helper methods for:
    - Creating projects and jobs
    - Setting parameters and uploading files
    - Running jobs and waiting for completion
    - Validating outputs

    Usage:
        @pytest.mark.usefixtures("file_based_db")
        class TestMyPipeline(APITestBase):
            task_name = "my_pipeline"

            def test_basic(self):
                self.create_project("test_my_pipeline")
                self.create_job()
                self.set_param("inputData.F_SIGF", demo_data("gamma", "merged.mtz"))
                self.upload_file("inputData.XYZIN", demo_data("gamma", "model.pdb"))
                self.run_and_wait()
                self.assert_file_exists("XYZOUT.pdb")
    """

    # Override in subclasses
    task_name = None
    timeout = 300  # Default job timeout in seconds
    poll_interval = 2  # Seconds between status checks

    @pytest.fixture(autouse=True)
    def setup_client(self):
        """Set up test client for each test."""
        self.client = Client()
        self.project_id = None
        self.job_id = None
        self.job_number = None
        self.project_dir = None

    # --- Project & Job Creation ---

    def create_project(self, name: str, directory: str = "__default__") -> dict:
        """Create a new project via API."""
        response = self.client.post(
            '/projects/',
            data=json.dumps({'name': name, 'directory': directory}),
            content_type='application/json'
        )
        assert response.status_code == 201, f"Failed to create project: {response.content}"
        data = response.json()
        self.project_id = data['id']
        self.project_dir = data.get('directory')
        return data

    def create_job(self, task_name: str = None, title: str = None, auto_context: bool = False) -> dict:
        """Create a job for the task under test."""
        task_name = task_name or self.task_name
        assert task_name, "task_name must be set on class or passed to create_job"
        assert self.project_id, "Must create project first"

        title = title or f"API Test: {task_name}"

        response = self.client.post(
            f'/projects/{self.project_id}/create_task/',
            data=json.dumps({
                'task_name': task_name,
                'title': title,
                'auto_context': auto_context
            }),
            content_type='application/json'
        )
        assert response.status_code == 200, f"Failed to create job: {response.content}"

        data = response.json()
        assert data.get('success'), f"Job creation failed: {data}"

        job_data = data['data']['new_job']
        self.job_id = job_data['id']
        self.job_number = job_data.get('number')
        return job_data

    # --- Parameter Setting ---

    def set_param(self, object_path: str, value, prefix: str = None) -> dict:
        """
        Set a parameter value on the job.

        Args:
            object_path: Path like "inputData.F_SIGF" (task prefix added automatically)
            value: Value to set (dict, list, string, etc.)
            prefix: Override task name prefix (default: self.task_name.container)
        """
        assert self.job_id, "Must create job first"

        # Build full object path
        if prefix is None:
            prefix = f"{self.task_name}.container"
        full_path = f"{prefix}.{object_path}" if prefix else object_path

        response = self.client.post(
            f'/jobs/{self.job_id}/set_parameter/',
            data=json.dumps({'object_path': full_path, 'value': value}),
            content_type='application/json'
        )
        assert response.status_code == 200, f"Failed to set parameter: {response.content}"
        return response.json()

    def add_list_item(self, object_path: str, item_data: dict = None) -> dict:
        """
        Add an item to a list parameter (like UNMERGEDFILES).

        Args:
            object_path: Path to the list (e.g., "inputData.UNMERGEDFILES")
            item_data: Optional initial data for the item (default: {})
        """
        item_data = item_data or {}
        return self.set_param(object_path, [item_data])

    # --- File Upload ---

    def upload_file(self, object_path: str, file_path: str, prefix: str = None) -> dict:
        """
        Upload a file to a parameter.

        Args:
            object_path: Path like "inputData.HKLIN" or "inputData.UNMERGEDFILES[0].file"
            file_path: Local path to the file to upload
            prefix: Override task name prefix
        """
        assert self.job_id, "Must create job first"

        if prefix is None:
            prefix = f"{self.task_name}.container"
        full_path = f"{prefix}.{object_path}" if prefix else object_path

        with open(file_path, 'rb') as f:
            response = self.client.post(
                f'/jobs/{self.job_id}/upload_file_param/',
                {'objectPath': full_path, 'file': f}
            )

        assert response.status_code == 200, f"Failed to upload file: {response.content}"
        return response.json()

    def upload_file_with_columns(self, object_path: str, file_path: str,
                                  column_labels: str = None, prefix: str = None) -> dict:
        """
        Upload an MTZ file with column selection.

        Args:
            object_path: Path like "inputData.F_SIGF"
            file_path: Path to MTZ file
            column_labels: Column selection like "/*/*/[FP,SIGFP]"
            prefix: Override task name prefix
        """
        assert self.job_id, "Must create job first"

        if prefix is None:
            prefix = f"{self.task_name}.container"
        full_path = f"{prefix}.{object_path}" if prefix else object_path

        with open(file_path, 'rb') as f:
            post_data = {'objectPath': full_path, 'file': f}
            if column_labels:
                post_data['column_selector'] = column_labels
            response = self.client.post(
                f'/jobs/{self.job_id}/upload_file_param/',
                post_data
            )

        assert response.status_code == 200, f"Failed to upload file: {response.content}"
        return response.json()

    # --- Validation ---

    def get_validation(self) -> dict:
        """Get job validation status."""
        assert self.job_id, "Must create job first"

        response = self.client.get(f'/jobs/{self.job_id}/validation/')
        assert response.status_code == 200, f"Failed to get validation: {response.content}"
        return response.json()

    def check_validation(self, allow_warnings: bool = True) -> bool:
        """
        Check that job passes validation.

        Args:
            allow_warnings: If False, also fail on warnings

        Returns:
            True if validation passes
        """
        validation = self.get_validation()

        has_errors = False
        has_warnings = False

        for item in validation.get('validation', []):
            level = item.get('level', '').upper()
            if level == 'ERROR':
                has_errors = True
                print(f"Validation ERROR: {item.get('message')}")
            elif level == 'WARNING':
                has_warnings = True
                if not allow_warnings:
                    print(f"Validation WARNING: {item.get('message')}")

        if has_errors:
            return False
        if has_warnings and not allow_warnings:
            return False
        return True

    # --- Job Execution ---

    def run_job(self) -> dict:
        """Start job execution."""
        assert self.job_id, "Must create job first"

        response = self.client.post(
            f'/jobs/{self.job_id}/run/',
            content_type='application/json'
        )
        assert response.status_code == 200, f"Failed to run job: {response.content}"
        return response.json()

    def get_job_status(self) -> dict:
        """Get current job status."""
        assert self.job_id, "Must create job first"

        response = self.client.get(f'/jobs/{self.job_id}/')
        assert response.status_code == 200, f"Failed to get job status: {response.content}"
        return response.json()

    def wait_for_completion(self, timeout: int = None) -> int:
        """
        Wait for job to complete.

        Args:
            timeout: Maximum seconds to wait (default: self.timeout)

        Returns:
            Final job status code
        """
        timeout = timeout or self.timeout
        start = time.time()

        while time.time() - start < timeout:
            job_data = self.get_job_status()
            status = job_data.get('status')

            if JobStatus.is_terminal(status):
                return status

            time.sleep(self.poll_interval)

        raise TimeoutError(f"Job did not complete within {timeout} seconds")

    def run_and_wait(self, timeout: int = None, expect_success: bool = True) -> int:
        """
        Run job and wait for completion.

        Args:
            timeout: Maximum seconds to wait
            expect_success: If True, assert job completed successfully

        Returns:
            Final job status code
        """
        self.run_job()
        status = self.wait_for_completion(timeout)

        if expect_success:
            assert JobStatus.is_success(status), \
                f"Job failed with status {status}"

        return status

    # --- Output Validation ---

    def get_job_files(self) -> dict:
        """Get list of job output files."""
        assert self.job_id, "Must create job first"

        response = self.client.get(f'/jobs/{self.job_id}/files/')
        assert response.status_code == 200, f"Failed to get job files: {response.content}"
        return response.json()

    def get_job_directory(self) -> Path:
        """Get path to job directory."""
        assert self.project_dir and self.job_number, "Must create and run job first"
        return Path(self.project_dir) / "CCP4_JOBS" / f"job_{self.job_number}"

    def assert_file_exists(self, filename: str) -> Path:
        """Assert that an output file exists and return its path."""
        job_dir = self.get_job_directory()
        file_path = job_dir / filename
        assert file_path.exists(), f"Expected output file not found: {file_path}"
        return file_path

    def read_output_file(self, filename: str) -> str:
        """Read content of an output file."""
        file_path = self.assert_file_exists(filename)
        return file_path.read_text()

    def read_program_xml(self, filename: str = "program.xml"):
        """Parse program.xml output file."""
        import xml.etree.ElementTree as ET
        file_path = self.assert_file_exists(filename)
        return ET.parse(file_path)

    # --- Gemmi-based Validation ---

    def validate_mtz(self, filename: str):
        """Validate that an MTZ file can be read."""
        if not HAS_GEMMI:
            pytest.skip("gemmi not available for MTZ validation")
        file_path = self.assert_file_exists(filename)
        return gemmi.read_mtz_file(str(file_path))

    def validate_pdb(self, filename: str):
        """Validate that a PDB file can be read."""
        if not HAS_GEMMI:
            pytest.skip("gemmi not available for PDB validation")
        file_path = self.assert_file_exists(filename)
        return gemmi.read_pdb(str(file_path))

    def validate_structure(self, filename: str, format=None):
        """Validate that a coordinate file can be read."""
        if not HAS_GEMMI:
            pytest.skip("gemmi not available for structure validation")
        file_path = self.assert_file_exists(filename)
        if format:
            return gemmi.read_structure(str(file_path), format=format)
        return gemmi.read_structure(str(file_path))

    # --- Report Access ---

    def get_job_report(self) -> dict:
        """Get job report data."""
        assert self.job_id, "Must create job first"

        response = self.client.get(f'/jobs/{self.job_id}/report/')
        assert response.status_code == 200, f"Failed to get report: {response.content}"
        return response.json()

    # --- API-based File Digest Validation ---
    # These methods use the /files/{id}/digest/ endpoint to validate
    # crystallographic files without needing local gemmi

    def get_output_file_id(self, param_name: str) -> int:
        """
        Get the database file ID for an output parameter.

        Args:
            param_name: Output parameter name like "XYZOUT" or "HKLOUT"

        Returns:
            File database ID
        """
        assert self.job_id, "Must create job first"

        # Get job files - API returns a list directly
        files = self.get_job_files()

        # Find file by param_name
        for file_info in files:
            if file_info.get('job_param_name') == param_name:
                return file_info.get('id')

        raise ValueError(f"Output file not found for param: {param_name}")

    def get_file_digest(self, file_id: int) -> dict:
        """
        Fetch the digest for a file via API.

        The digest endpoint returns crystallographic metadata:
        - MTZ: cell, spaceGroup, lowRes, highRes, columns, etc.
        - PDB/mmCIF: cell, spaceGroup, nModels, chains, ligands, etc.

        Args:
            file_id: Database file ID

        Returns:
            Digest dictionary with file metadata
        """
        response = self.client.get(f'/files/{file_id}/digest/')
        assert response.status_code == 200, f"Failed to get file digest: {response.content}"
        result = response.json()
        # Unwrap the data if wrapped in {success: True, data: {...}}
        if isinstance(result, dict) and 'data' in result:
            return result['data']
        return result

    def get_output_digest(self, param_name: str) -> dict:
        """
        Get the digest for an output file by parameter name.

        Args:
            param_name: Output parameter name like "XYZOUT" or "HKLOUT"

        Returns:
            Digest dictionary
        """
        file_id = self.get_output_file_id(param_name)
        return self.get_file_digest(file_id)

    def assert_valid_mtz_output(self, param_name: str) -> dict:
        """
        Assert that an MTZ output file exists and can be parsed via API.

        Successfully fetching the digest validates that:
        1. The file was registered in the database
        2. The file exists on disk
        3. The file is a valid MTZ that can be read by gemmi (server-side)

        Args:
            param_name: Output parameter name (e.g., "HKLOUT", "FPHIOUT")

        Returns:
            The digest dictionary (with cell, spaceGroup, lowRes, highRes, columns)
        """
        digest = self.get_output_digest(param_name)
        # A valid MTZ digest should have spaceGroup and resolution info
        assert 'spaceGroup' in digest or 'lowRes' in digest, \
            f"Invalid MTZ digest for {param_name}: {digest}"
        return digest

    def assert_valid_coords_output(self, param_name: str) -> dict:
        """
        Assert that a coordinate output file exists and can be parsed via API.

        Successfully fetching the digest validates that:
        1. The file was registered in the database
        2. The file exists on disk
        3. The file is a valid PDB/mmCIF that can be read by gemmi (server-side)

        Args:
            param_name: Output parameter name (e.g., "XYZOUT")

        Returns:
            The digest dictionary (with cell, spaceGroup, chains, etc.)
        """
        digest = self.get_output_digest(param_name)
        # A valid coordinate digest should have structural info
        assert 'spaceGroup' in digest or 'chains' in digest or 'nModels' in digest, \
            f"Invalid coordinate digest for {param_name}: {digest}"
        return digest

    def validate_mtz_via_api(self, param_name: str,
                              expected_spacegroup: str = None,
                              expected_resolution_low: float = None,
                              expected_resolution_high: float = None,
                              required_columns: list = None) -> dict:
        """
        Validate an MTZ output file using the API digest endpoint.

        This is an alternative to validate_mtz() that doesn't require gemmi
        to be installed locally and exercises the server-side validation.

        Args:
            param_name: Output parameter name (e.g., "HKLOUT", "FPHIOUT")
            expected_spacegroup: Expected spacegroup name (optional)
            expected_resolution_low: Expected low resolution limit (optional)
            expected_resolution_high: Expected high resolution limit (optional)
            required_columns: List of required column names (optional)

        Returns:
            The digest dictionary
        """
        digest = self.get_output_digest(param_name)

        if expected_spacegroup:
            actual_sg = digest.get('spaceGroup')
            assert actual_sg == expected_spacegroup, \
                f"Spacegroup mismatch: expected {expected_spacegroup}, got {actual_sg}"

        if expected_resolution_low is not None:
            low_res = digest.get('lowRes')
            assert low_res is not None, "No low resolution in digest"
            assert abs(low_res - expected_resolution_low) < 1.0, \
                f"Low resolution mismatch: expected ~{expected_resolution_low}, got {low_res}"

        if expected_resolution_high is not None:
            high_res = digest.get('highRes')
            assert high_res is not None, "No high resolution in digest"
            assert abs(high_res - expected_resolution_high) < 0.1, \
                f"High resolution mismatch: expected ~{expected_resolution_high}, got {high_res}"

        if required_columns:
            columns = digest.get('columns', [])
            column_names = [c.get('label') for c in columns]
            for col in required_columns:
                assert col in column_names, \
                    f"Required column {col} not found. Available: {column_names}"

        return digest

    def validate_coords_via_api(self, param_name: str,
                                 expected_spacegroup: str = None,
                                 expected_chains: list = None,
                                 min_residues: int = None) -> dict:
        """
        Validate a coordinate output file using the API digest endpoint.

        Args:
            param_name: Output parameter name (e.g., "XYZOUT")
            expected_spacegroup: Expected spacegroup name (optional)
            expected_chains: Expected chain IDs (optional)
            min_residues: Minimum expected residue count (optional)

        Returns:
            The digest dictionary
        """
        digest = self.get_output_digest(param_name)

        if expected_spacegroup:
            actual_sg = digest.get('spaceGroup')
            assert actual_sg == expected_spacegroup, \
                f"Spacegroup mismatch: expected {expected_spacegroup}, got {actual_sg}"

        if expected_chains:
            chains = digest.get('chains', [])
            chain_ids = [c.get('chainId') for c in chains]
            for chain in expected_chains:
                assert chain in chain_ids, \
                    f"Expected chain {chain} not found. Available: {chain_ids}"

        if min_residues:
            total_residues = sum(c.get('residueCount', 0) for c in digest.get('chains', []))
            assert total_residues >= min_residues, \
                f"Too few residues: expected >= {min_residues}, got {total_residues}"

        return digest

    def validate_cell(self, digest: dict,
                       expected_a: float = None,
                       expected_b: float = None,
                       expected_c: float = None,
                       tolerance: float = 0.5) -> None:
        """
        Validate unit cell parameters from a digest.

        Args:
            digest: File digest dictionary
            expected_a, expected_b, expected_c: Expected cell dimensions
            tolerance: Allowed deviation in Angstroms
        """
        cell = digest.get('cell', {})

        if expected_a is not None:
            actual_a = cell.get('a')
            assert actual_a is not None, "No cell dimension 'a' in digest"
            assert abs(actual_a - expected_a) < tolerance, \
                f"Cell a mismatch: expected ~{expected_a}, got {actual_a}"

        if expected_b is not None:
            actual_b = cell.get('b')
            assert actual_b is not None, "No cell dimension 'b' in digest"
            assert abs(actual_b - expected_b) < tolerance, \
                f"Cell b mismatch: expected ~{expected_b}, got {actual_b}"

        if expected_c is not None:
            actual_c = cell.get('c')
            assert actual_c is not None, "No cell dimension 'c' in digest"
            assert abs(actual_c - expected_c) < tolerance, \
                f"Cell c mismatch: expected ~{expected_c}, got {actual_c}"
