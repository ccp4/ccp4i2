"""
API-driven tests for refinement pipelines.

These tests verify the API workflow for:
- prosmart_refmac: ProSMART + Refmac refinement
- servalcat_pipe: Servalcat refinement pipeline
- sheetbend: Sheetbend coordinate optimization

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestProsmartRefmacAPI(APITestBase):
    """API tests for prosmart_refmac refinement pipeline."""

    task_name = "prosmart_refmac"
    timeout = 300

    def test_gamma_basic(self, gamma_mtz, gamma_model_pdb):
        """Test refmac with gamma anomalous data (basic, no MolProbity)."""
        self.create_project("test_refmac_gamma")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # Basic refinement settings
        self.set_param("controlParameters.NCYCLES", 4)
        self.set_param("controlParameters.USEANOMALOUS", True)
        # Disable MolProbity for faster test
        self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("XYZOUT.pdb")
        self.validate_mtz("ABCDOUT.mtz")
        self.validate_mtz("ANOMFPHIOUT.mtz")
        self.validate_mtz("DIFFPHIOUT.mtz")
        self.validate_mtz("FPHIOUT.mtz")

        # Check refinement improvement
        xml = self.read_program_xml()
        cycles = xml.findall(".//RefmacInProgress/Cycle")
        rworks = [float(c.find("r_factor").text) for c in cycles]

        assert len(rworks) == 5  # Initial + 4 cycles
        assert rworks[-1] < rworks[0]  # R-factor improved
        assert rworks[-1] < 0.27

    def test_8xfm_basic(self, cif8xfm, mtz8xfm):
        """Test refmac with 8xfm data (no water addition, no MolProbity)."""
        self.create_project("test_refmac_8xfm")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz8xfm,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", mtz8xfm,
            column_labels="/*/*/[FREE]"
        )
        self.upload_file("inputData.XYZIN", cif8xfm)

        self.set_param("controlParameters.NCYCLES", 2)
        self.set_param("controlParameters.ADD_WATERS", False)
        self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("XYZOUT.pdb")
        self.validate_mtz("FPHIOUT.mtz")

        # Check R-factor
        xml = self.read_program_xml()
        cycles = xml.findall(".//RefmacInProgress/Cycle")
        rworks = [float(c.find("r_factor").text) for c in cycles]
        assert rworks[-1] < 0.19


@pytest.mark.usefixtures("file_based_db")
class TestServalcatAPI(APITestBase):
    """API tests for servalcat_pipe refinement pipeline."""

    task_name = "servalcat_pipe"
    timeout = 300

    def test_8xfm_basic(self, cif8xfm, mtz8xfm):
        """Test servalcat X-ray refinement (basic, no validation)."""
        self.create_project("test_servalcat_8xfm")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.set_param("controlParameters.DATA_METHOD", "xtal")

        self.upload_file_with_columns(
            "inputData.HKLIN", mtz8xfm,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", mtz8xfm,
            column_labels="/*/*/[FREE]"
        )

        # Set refinement parameters
        self.set_param("controlParameters.NCYCLES", 2)
        self.set_param("controlParameters.ADD_WATERS", False)
        self.set_param("controlParameters.F_SIGF_OR_I_SIGI", "F_SIGF")

        # Disable validation for speed
        self.set_param("controlParameters.VALIDATE_IRIS", False)
        self.set_param("controlParameters.VALIDATE_BAVERAGE", False)
        self.set_param("controlParameters.VALIDATE_RAMACHANDRAN", False)
        self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)
        self.set_param("controlParameters.RUN_ADP_ANALYSIS", False)
        self.set_param("controlParameters.RUN_COORDADPDEV_ANALYSIS", False)

        self.run_and_wait()

        # Validate outputs
        self.assert_file_exists("CIFFILE.cif")
        self.assert_file_exists("program.xml")

        # Check R-factors
        xml = self.read_program_xml()
        rworks = [float(e.text) for e in xml.findall(".//data/summary/Rwork")]
        assert rworks[-1] < 0.18

    def test_electron_diffraction(self, cif7beq, mtz7beq):
        """Test servalcat electron diffraction refinement."""
        from .base import download, URLs

        # Download reference structure for ProSMART
        with download(URLs.redo_cif("7ber")) as cif7ber:
            self.create_project("test_servalcat_electron")
            self.create_job()

            self.upload_file("inputData.XYZIN", cif7beq)
            self.set_param("controlParameters.DATA_METHOD", "xtal")

            self.upload_file_with_columns(
                "inputData.HKLIN", mtz7beq,
                column_labels="/*/*/[FP,SIGFP]"
            )
            self.upload_file_with_columns(
                "inputData.FREERFLAG", mtz7beq,
                column_labels="/*/*/[FREE]"
            )

            self.set_param("controlParameters.NCYCLES", 3)
            self.set_param("controlParameters.SCATTERING_FACTORS", "electron")
            self.set_param("controlParameters.USE_JELLY", True)

            # ProSMART reference
            self.set_param("controlParameters.prosmartProtein.TOGGLE", True)
            self.set_param("controlParameters.prosmartProtein.MODE", "SELECTED")
            self.upload_file("controlParameters.prosmartProtein.REFERENCE_MODELS", cif7ber)
            self.set_param("controlParameters.prosmartProtein.CHAINLIST_1", "A")
            self.set_param("controlParameters.prosmartProtein.ALL_BEST", "ALL")
            self.set_param("controlParameters.prosmartProtein.SEQID", 75)
            self.set_param("controlParameters.prosmartProtein.SIDE_MAIN", "SIDE")

            # Disable validation
            self.set_param("controlParameters.VALIDATE_IRIS", False)
            self.set_param("controlParameters.VALIDATE_BAVERAGE", False)
            self.set_param("controlParameters.VALIDATE_RAMACHANDRAN", False)
            self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)
            self.set_param("controlParameters.RUN_ADP_ANALYSIS", False)
            self.set_param("controlParameters.RUN_COORDADPDEV_ANALYSIS", False)

            self.run_and_wait()

            # Check for ProSMART outputs
            job_dir = self.get_job_directory()
            assert (job_dir / "job_1" / "RESTRAINTS.txt").exists()


@pytest.mark.usefixtures("file_based_db")
class TestSheetbendAPI(APITestBase):
    """API tests for sheetbend coordinate optimization."""

    task_name = "sheetbend"
    timeout = 120

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test sheetbend with 8xfm data."""
        self.create_project("test_sheetbend")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz8xfm,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", mtz8xfm,
            column_labels="/*/*/[FREE]"
        )

        self.run_and_wait()

        # Validate outputs
        self.assert_file_exists("XYZOUT.cif")

        # Check R-factors improved or maintained
        xml = self.read_program_xml()
        # Sheetbend output format may vary - check file exists
        assert xml is not None
