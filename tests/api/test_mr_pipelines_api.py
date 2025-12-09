"""
API-driven tests for molecular replacement pipelines.

These tests verify the API workflow for:
- phaser_simple: Simple MR pipeline
- molrep_pipe: Molrep MR pipeline
- i2Dimple: Dimple pipeline
- mrbump_basic: MrBUMP automated MR

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestPhaserSimpleAPI(APITestBase):
    """API tests for phaser_simple MR pipeline."""

    task_name = "phaser_simple"
    timeout = 300  # Phaser can take a while

    def test_gamma_basic(self, gamma_mtz, gamma_model_pdb):
        """Test phaser_simple without post-processing."""
        self.create_project("test_phaser_simple_basic")
        self.create_job()

        # Upload input files
        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # Disable post-processing for faster test
        self.set_param("controlParameters.RUNREFMAC", False)
        self.set_param("controlParameters.RUNSHEETBEND", False)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("PHASER.1.pdb")
        self.validate_mtz("DIFMAPOUT_1.mtz")
        self.validate_mtz("MAPOUT_1.mtz")
        self.validate_mtz("PHASEOUT_1.mtz")

        # Check LLG from program.xml
        xml = self.read_program_xml()
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        assert max(llgs) > 1000

    def test_gamma_with_sheetbend(self, gamma_mtz, gamma_model_pdb):
        """Test phaser_simple with sheetbend post-processing."""
        self.create_project("test_phaser_sheetbend")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.set_param("controlParameters.RUNREFMAC", False)
        self.set_param("controlParameters.RUNSHEETBEND", True)

        self.run_and_wait()

        # Validate outputs including sheetbend
        self.validate_pdb("PHASER.1.pdb")
        self.validate_pdb("XYZOUT_SHEETBEND.pdb")
        self.validate_mtz("PHASEOUT_1.mtz")

    def test_gamma_full(self, gamma_mtz, gamma_model_pdb):
        """Test phaser_simple with full post-processing (refmac + sheetbend)."""
        self.create_project("test_phaser_full")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)
        # Can add multiple search models
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # Full pipeline with refmac and sheetbend
        self.set_param("controlParameters.RUNREFMAC", True)
        self.set_param("controlParameters.RUNSHEETBEND", True)

        self.run_and_wait(timeout=600)  # Longer timeout for refinement

        # Validate outputs
        self.validate_pdb("PHASER.1.pdb")
        self.validate_pdb("XYZOUT_REFMAC.pdb")
        self.validate_pdb("XYZOUT_SHEETBEND.pdb")

        # Check R-factor improvement
        xml = self.read_program_xml()
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        assert min(rworks) < 0.27

    def test_no_solution(self, gamma_mtz, rnase_model_pdb):
        """Test phaser with mismatched model (no solution expected)."""
        self.create_project("test_phaser_no_solution")
        self.create_job()

        # Use gamma data with rnase model - shouldn't find solution
        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.XYZIN", rnase_model_pdb)

        # Lower resolution to speed up
        self.set_param("controlParameters.RESOLUTION_HIGH", 5.0)

        # Expect failure - no solution
        self.run_and_wait(expect_success=False)

        # Should not produce output
        job_dir = self.get_job_directory()
        assert not (job_dir / "PHASER.1.pdb").exists()


@pytest.mark.usefixtures("file_based_db")
class TestMolrepPipeAPI(APITestBase):
    """API tests for molrep_pipe MR pipeline."""

    task_name = "molrep_pipe"
    timeout = 300

    def test_gamma(self, gamma_mtz, gamma_freer_mtz, gamma_model_pdb):
        """Test molrep_pipe with gamma data."""
        self.create_project("test_molrep")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.FREERFLAG", gamma_freer_mtz)
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("XYZOUT_MOLREP.pdb")
        self.validate_pdb("XYZOUT_SHEETBEND.pdb")
        self.validate_pdb("XYZOUT.pdb")

        # Check R-factor improvement
        xml = self.read_program_xml()
        rworks = [float(e.text) for e in xml.findall(".//Cycle/r_factor")]
        rfrees = [float(e.text) for e in xml.findall(".//Cycle/r_free")]

        assert rworks[-1] < rworks[0]
        assert rfrees[-1] < rfrees[0]
        assert rworks[-1] < 0.26
        assert rfrees[-1] < 0.28


@pytest.mark.usefixtures("file_based_db")
class TestDimpleAPI(APITestBase):
    """API tests for i2Dimple pipeline."""

    task_name = "i2Dimple"
    timeout = 300

    def test_gamma(self, gamma_native_mtz, gamma_model_pdb):
        """Test dimple with gamma native data."""
        self.create_project("test_dimple")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_native_mtz)
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.run_and_wait()

        # Validate outputs
        self.validate_mtz("FPHIOUT.mtz")
        self.assert_file_exists("program.xml")


@pytest.mark.usefixtures("file_based_db")
class TestMrBumpAPI(APITestBase):
    """API tests for mrbump_basic automated MR."""

    task_name = "mrbump_basic"
    timeout = 600  # MrBUMP searches databases, can be slow

    def test_gamma(self, gamma_mtz, gamma_freer_mtz, gamma_asu_xml):
        """Test MrBUMP with gamma data (sequence-based search)."""
        self.create_project("test_mrbump")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.FREERFLAG", gamma_freer_mtz)
        self.upload_file("inputData.ASUIN", gamma_asu_xml)

        # Limit search for faster testing
        self.set_param("controlParameters.MRMAX", 5)
        self.set_param("controlParameters.REDUNDANCYLEVEL", 95)
        self.set_param("controlParameters.PJOBS", 2)
        self.set_param("controlParameters.NCYC", 10)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("output_mrbump_1.pdb")
        self.validate_mtz("output_mrbump_1.mtz")

        # Check results file
        results_path = self.get_job_directory() / "search_mrbump_1" / "results" / "results.txt"
        if results_path.exists():
            lines = results_path.read_text().splitlines()
            best_rfree = 1.0
            for i, line in enumerate(lines):
                if line.endswith("Rfree"):
                    try:
                        rfree = float(lines[i + 1].split()[-1])
                        best_rfree = min(best_rfree, rfree)
                    except (ValueError, IndexError):
                        pass
            assert best_rfree < 0.35
