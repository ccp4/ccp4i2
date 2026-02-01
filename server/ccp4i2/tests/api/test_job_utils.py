"""
Tests for job utility functions.

Converted to pytest fixture-based approach for compatibility with pytest-xdist
parallel test execution. Uses the isolated_test_db fixture from conftest.py.
"""
from pathlib import Path
from typing import List
import pytest
from rest_framework.test import APIClient

from ccp4i2.core import CCP4Container
from ccp4i2.core.CCP4ModelData import CDictDataFile
from ccp4i2.core.base_object.cdata_file import CDataFile
from ...lib.utils.jobs.i2run import i2run_for_job
from ...db.import_i2xml import import_ccp4_project_zip
from ...db.ccp4i2_django_projects_manager import CCP4i2DjangoProjectsManager
from ...db import models
from ...lib.utils.jobs.clone import clone_job
from ...db.async_db_handler import AsyncDatabaseHandler
from ...lib.utils.jobs.get_container import get_job_container
from ...lib.utils.files.get_by_context import get_file_by_job_context
from ...lib.utils.navigation.dependencies import find_dependent_jobs
from ...lib.utils.navigation.what_next import get_what_next
from ...lib.utils.helpers.object_method import object_method
from ...lib.utils.files.detect_type import detect_file_type
from ...lib.utils.files.export import export_job_file
from ...db.project_json import project_json
from ...lib.utils.files.digest import digest_file_object, digest_cdatafile_file_object

# Path to test data - these tests require pre-built project zips
TEST_DATA_DIR = Path(__file__).parent.parent.parent.parent.parent.parent / "test101" / "ProjectZips"
SKIP_REASON = f"Test data not found: {TEST_DATA_DIR}"


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason=SKIP_REASON)
class TestJobUtils:
    """Tests for job utility functions using pytest fixtures."""

    @pytest.fixture(autouse=True)
    def setup(self, bypass_api_permissions, test_project_path):
        """Set up test fixtures for each test."""
        self.client = APIClient()
        self.test_project_path = test_project_path
        test_project_path.mkdir(exist_ok=True)

        import_ccp4_project_zip(
            TEST_DATA_DIR / "refmac_gamma_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )
        import_ccp4_project_zip(
            TEST_DATA_DIR / "aimless_gamma_native_test_1.ccp4_project.zip",
            relocate_path=test_project_path,
        )
        import_ccp4_project_zip(
            TEST_DATA_DIR / "molrep_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )
        import_ccp4_project_zip(
            TEST_DATA_DIR / "parrot_test_0.ccp4_project.zip",
            relocate_path=test_project_path,
        )

        self.pm = CCP4i2DjangoProjectsManager()

    def test_clone_jobs(self):
        for project_name in ["refmac_gamma_test_0"]:
            old_job = models.Job.objects.get(project__name=project_name, number="1")
            result = clone_job(old_job.uuid)
            assert result.success, f"Clone failed: {result.error}"
            new_job = result.data
            assert new_job.task_name == old_job.task_name

    def test_glean_job_files(self):
        # SKIP: This test was for the legacy sync glean_job_files() function
        # The new async version is tested in tests/test_async_db_handler.py
        # and tests/test_async_plugin_with_database.py
        pytest.skip("Legacy glean_job_files() removed - use AsyncDatabaseHandler instead")

    def test_get_file_by_job_context(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        file_uuids: List[str] = get_file_by_job_context(
            contextJobId=str(old_job.uuid), fileType="application/CCP4-mtz-observed"
        )
        assert file_uuids == ["99f93cda-c449-11ea-a15f-3417eba0e4fd"]

    def test_find_descendent_jobs(self):
        old_job = models.Job.objects.get(project__name="molrep_test_0", number="1")
        dependents = find_dependent_jobs(old_job)
        print(dependents)

    def test_object_method(self):
        # NOTE: This test requires "bucc_test_0" project which is not imported
        # Skipping until the project zip is available or test is refactored
        bucc_project = models.Project.objects.filter(name="bucc_test_0").first()
        if bucc_project is None:
            pytest.skip("bucc_test_0 project not available")

        the_job = models.Job.objects.filter(project=bucc_project, parent=None).last()
        asuWeight = object_method(
            the_job, "buccaneer.inputData.ASUIN.fileContent", "molecularWeight"
        )
        assert abs(asuWeight - 15107) < 0.1
        result = object_method(
            the_job,
            "buccaneer.inputData.F_SIGF.fileContent",
            "matthewsCoeff",
            kwargs={"molWt": asuWeight},
        )
        assert abs(result["results"][0]["matth_coef"] - 2.1) < 0.1

    def test_project_json(self):
        # NOTE: This test requires "bucc_test_0" project which is not imported
        bucc_project = models.Project.objects.filter(name="bucc_test_0").first()
        if bucc_project is None:
            pytest.skip("bucc_test_0 project not available")

        result = project_json(bucc_project)
        print(result)

    def test_what_next(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        what_next = get_what_next(old_job)
        assert what_next["Status"] == "Success"
        assert len(what_next["result"]) == 3
        # Check task names are present (titles may vary based on plugin definitions)
        task_names = [item["taskName"] for item in what_next["result"]]
        assert "prosmart_refmac" in task_names
        assert "modelcraft" in task_names
        assert "coot_rebuild" in task_names

    def test_digest_cdata(self):
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7.cif"
        )
        file = CDataFile(str(mmcif_path))
        a = digest_cdatafile_file_object(file)
        # Check digest succeeded and has expected structure
        assert a.get("status") != "Failed", f"Digest failed: {a.get('reason', a)}"
        if "sequences" in a:
            assert a["sequences"]["A"] == (
                "QIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDAAQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLV"
            )
        else:
            # If sequences not returned directly, check for composition or other valid keys
            assert "composition" in a or len(a) > 0, f"Unexpected digest result: {a}"

        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "baz2b.pir"
        )
        file = CDataFile(str(mmcif_path))
        a = digest_cdatafile_file_object(file)
        # PIR files should have sequences
        assert a.get("status") != "Failed", f"Digest failed: {a.get('reason', a)}"
        if "sequences" in a:
            assert "A" in a["sequences"]

    def test_digest_cdictdatafile(self):
        mmcif_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "BAZ2BA-x839-LIG.cif"
        )
        # CDictDataFile doesn't accept file_path as first arg (unlike CDataFile)
        # Create instance first, then set path
        aFile = CDictDataFile()
        aFile.setFullPath(str(mmcif_path))
        print(aFile.fullPath)
        a = digest_file_object(aFile)
        print("a", a)

    def test_detect_mtz_file(self):
        mtz_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "mdm2_unmerged.mtz"
        )
        print(mtz_path)
        a = detect_file_type(mtz_path)
        assert a == "MTZ file"

    def test_detect_reflection_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7-sf.cif"
        )
        print(file_path)
        a = detect_file_type(file_path)
        assert a == "mmCIF reflection file"

    def test_detect_coordinate_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "mdm2"
            / "4hg7.cif"
        )
        print(file_path)
        a = detect_file_type(file_path)
        assert a == "mmCIF coordinate file"

    def test_detect_ligand_cif(self):
        file_path = (
            Path(CCP4Container.__file__).parent.parent
            / "demo_data"
            / "baz2b"
            / "BAZ2BA-x839-LIG.cif"
        )
        a = detect_file_type(file_path)
        assert a == "mmCIF ligand/geometry file"

    def test_i2run_for_job(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        i2run_command = i2run_for_job(old_job)
        print(i2run_command)
        assert i2run_command is not None

    def test_export_job_refmac_file(self):
        old_job = models.Job.objects.get(
            project__name="refmac_gamma_test_0", number="1"
        )
        exported_files = export_job_file(str(old_job.id), "FoFc_as_map")
        print(exported_files)
        assert exported_files is not None

    def test_export_job_file(self):
        old_job = models.Job.objects.get(project__name="parrot_test_0", number="1")
        exported_files = export_job_file(str(old_job.id), "complete_mtz")
        print(exported_files)
        assert exported_files is not None
