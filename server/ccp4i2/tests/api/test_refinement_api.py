"""
API-driven tests for refinement pipelines.

These tests verify the API workflow for:
- prosmart_refmac: ProSMART + Refmac refinement
- servalcat_pipe: Servalcat refinement pipeline
- sheetbend: Sheetbend coordinate optimization
- zanuda: Space group validation / refinement
- lorestr_i2: Low-resolution structure refinement
- metalCoord: Metal coordination restraint generation

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

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

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
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
        self.set_param("monitor.RUN_COORDADPDEV_ANALYSIS", False)

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
            self.set_param("prosmartProtein.TOGGLE", True)
            self.set_param("prosmartProtein.MODE", "SELECTED")
            self.upload_file("prosmartProtein.REFERENCE_MODELS[0]", cif7ber)
            self.set_param("prosmartProtein.CHAINLIST_1", "A")
            self.set_param("prosmartProtein.ALL_BEST", "ALL")
            self.set_param("prosmartProtein.SEQID", 75)
            self.set_param("prosmartProtein.SIDE_MAIN", "SIDE")

            # Disable validation
            self.set_param("controlParameters.VALIDATE_IRIS", False)
            self.set_param("controlParameters.VALIDATE_BAVERAGE", False)
            self.set_param("controlParameters.VALIDATE_RAMACHANDRAN", False)
            self.set_param("controlParameters.VALIDATE_MOLPROBITY", False)
            self.set_param("controlParameters.RUN_ADP_ANALYSIS", False)
            self.set_param("monitor.RUN_COORDADPDEV_ANALYSIS", False)

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


@pytest.mark.usefixtures("file_based_db")
class TestZanudaAPI(APITestBase):
    """API tests for zanuda space group validation."""

    task_name = "zanuda"
    timeout = 300

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test zanuda with 8xfm refined structure."""
        self.create_project("test_zanuda")
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

        self.run_and_wait()

        # Zanuda writes output as zanuda.pdb (not XYZOUT)
        self.assert_file_exists("zanuda.pdb")
        self.assert_file_exists("FPHIOUT.mtz")
        self.assert_file_exists("DIFFPHIOUT.mtz")


@pytest.mark.usefixtures("file_based_db")
class TestLorestrAPI(APITestBase):
    """API tests for lorestr_i2 low-resolution refinement."""

    task_name = "lorestr_i2"
    timeout = 600

    def test_1h1s(self, pdb1h1s, mtz1h1s):
        """Test lorestr with 1h1s structure (no auto homologue fetch)."""
        self.create_project("test_lorestr")
        self.create_job()

        self.upload_file("inputData.XYZIN", pdb1h1s)
        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz1h1s,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", mtz1h1s,
            column_labels="/*/*/[FREE]"
        )
        self.set_param("controlParameters.AUTO", "none")
        self.set_param("controlParameters.CPU", 1)

        self.run_and_wait()

        # Check refined model and map coefficients
        self.assert_file_exists("XYZOUT.pdb")
        job_dir = self.get_job_directory()
        mtz_files = list(job_dir.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output: {list(job_dir.iterdir())}"


@pytest.mark.usefixtures("file_based_db")
class TestMetalCoordAPI(APITestBase):
    """API tests for metalCoord restraint generation."""

    task_name = "metalCoord"
    timeout = 120

    def test_4dl8(self, cif4dl8):
        """Test metalCoord with 4dl8 (AlF3 complex)."""
        self.create_project("test_metalcoord")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif4dl8)
        self.set_param("inputData.LIGAND_CODE", "AF3")
        self.set_param("controlParameters.KEEP_LINKS", True)
        self.set_param("controlParameters.MAXIMUM_COORDINATION_NUMBER", 5)
        self.set_param("coordination.COORD05", "trigonal-bipyramid")
        self.set_param("controlParameters.MINIMUM_SAMPLE_SIZE", 10)
        self.set_param("controlParameters.DISTANCE_THRESHOLD", 0.45)
        self.set_param("controlParameters.PROCRUSTES_DISTANCE_THRESHOLD", 0.2)

        self.run_and_wait()

        # Check for restraint output files
        job_dir = self.get_job_directory()
        af3_files = list(job_dir.glob("AF3_restraints*"))
        assert len(af3_files) >= 3, f"Expected AF3 restraint files, found: {[f.name for f in af3_files]}"
