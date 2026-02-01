"""
API-driven tests for data reduction pipelines.

These tests verify the API workflow for:
- aimless_pipe: Data scaling and merging
- import_merged: Import merged reflections
- freerflag: Generate free R flags

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs via API digest endpoints.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline
from pytest import approx

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestAimlessPipeAPI(APITestBase):
    """API tests for aimless_pipe data scaling pipeline."""

    task_name = "aimless_pipe"
    timeout = 180  # Aimless can take a while, especially for larger datasets

    def test_gamma(self, gamma_unmerged_mtz):
        """Test aimless_pipe with gamma native unmerged data."""
        # Create project and job
        self.create_project("test_aimless_gamma")
        self.create_job()

        # Add list item and upload unmerged MTZ
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", gamma_unmerged_mtz)

        # Run and wait
        self.run_and_wait()

        # Validate outputs via API digest - confirms files are valid and readable
        # FREEROUT is the FreeR flag file output
        # HKLOUT is a list of merged data files, so use HKLOUT[0] for the first item
        self.assert_valid_mtz_output("FREEROUT")
        digest = self.assert_valid_mtz_output("HKLOUT[0]")

        # Check crystallographic properties from digest
        assert digest.get('spaceGroup') == "P 21 21 21"
        assert digest.get('highRes') == approx(1.81, rel=0.05)

        # Also verify via program.xml for additional metrics
        xml = self.read_program_xml()
        rmeas = float(xml.find(".//Rmeas/Overall").text)
        assert rmeas == approx(0.061, rel=0.1)

    def test_mdm2(self, mdm2_unmerged_mtz):
        """Test aimless_pipe with mdm2 unmerged data."""
        self.create_project("test_aimless_mdm2")
        self.create_job()

        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", mdm2_unmerged_mtz)

        self.run_and_wait()

        # Validate outputs via API digest
        self.assert_valid_mtz_output("FREEROUT")
        digest = self.assert_valid_mtz_output("HKLOUT[0]")

        # Check spacegroup and resolution from digest
        assert digest.get('spaceGroup') == "P 61 2 2"
        assert digest.get('highRes') == approx(1.25, rel=0.05)

        # Verify statistics from program.xml
        xml = self.read_program_xml()
        rmeas = float(xml.find(".//Rmeas/Overall").text)
        assert rmeas == approx(0.068, rel=0.1)


@pytest.mark.usefixtures("file_based_db")
class TestImportMergedAPI(APITestBase):
    """API tests for import_merged pipeline."""

    task_name = "import_merged"
    timeout = 60

    def test_gamma_mtz(self, gamma_mtz, gamma_freer_mtz):
        """Test importing merged MTZ with existing FreeR."""
        self.create_project("test_import_merged_gamma")
        self.create_job()

        # Upload input files with column specifications
        self.upload_file_with_columns(
            "inputData.HKLIN", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", gamma_freer_mtz,
            column_labels="/*/*/[FREER]"
        )

        self.run_and_wait()

        # Validate outputs via API digest
        obs_digest = self.assert_valid_mtz_output("OBSOUT")
        free_digest = self.assert_valid_mtz_output("FREEOUT")

        # Check that OBSOUT has expected structure
        assert obs_digest.get('spaceGroup') is not None
        # FREEOUT should have hasFreeR=True or be a valid digest
        assert free_digest.get('hasFreeR', True) or free_digest.get('spaceGroup') is not None, \
            f"Invalid FreeR digest: {free_digest}"

    def test_from_sf_cif(self):
        """Test importing merged data from structure factor CIF."""
        from .base import download, URLs

        with download(URLs.pdbe_sfcif("2ceu")) as cif_path:
            self.create_project("test_import_merged_cif")
            self.create_job()

            self.upload_file("inputData.HKLIN", cif_path)
            self.set_param("inputData.SPACEGROUP", "I 2 2 2")

            self.run_and_wait()

            # Validate outputs via API digest
            obs_digest = self.assert_valid_mtz_output("OBSOUT")
            self.assert_valid_mtz_output("FREEOUT")

            # Spacegroup should match what we specified
            assert obs_digest.get('spaceGroup') == "I 2 2 2"


@pytest.mark.usefixtures("file_based_db")
class TestFreeRFlagAPI(APITestBase):
    """API tests for freerflag generation."""

    task_name = "freerflag"
    timeout = 30

    def test_generate_freer(self, gamma_mtz):
        """Test generating free R flags."""
        self.create_project("test_freerflag")
        self.create_job()

        # Upload with column specification for anomalous intensities
        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.set_param("controlParameters.COMPLETE", False)
        self.set_param("controlParameters.FRAC", 0.2)

        self.run_and_wait()

        # Validate output via API digest
        digest = self.assert_valid_mtz_output("FREEROUT")

        # The digest confirms the file is valid
        # For detailed FreeR flag distribution, would need to read file directly
        # or have the server compute it in the digest
        assert digest.get('spaceGroup') is not None
