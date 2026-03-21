"""
API-driven tests for structure analysis tasks.

These tests verify the API workflow for:
- gesamt: Structure superposition
- privateer: Glycan validation
- edstats: Electron density statistics

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from .base import APITestBase, download, URLs


@pytest.mark.usefixtures("file_based_db")
class TestGesamtAPI(APITestBase):
    """API tests for gesamt structure superposition."""

    task_name = "gesamt"
    timeout = 60

    def test_self_superposition(self, cif8xfm):
        """Test gesamt pairwise superposition (identity)."""
        self.create_project("test_gesamt")
        self.create_job()

        # Superpose structure onto itself — should give RMSD ~ 0
        self.upload_file("inputData.XYZIN_QUERY", cif8xfm)
        self.upload_file("inputData.XYZIN_TARGET", cif8xfm)

        self.run_and_wait()

        # Check output coordinates and program.xml
        self.assert_file_exists("XYZOUT.cif")
        xml = self.read_program_xml()
        t = xml.find(".//Transformation")
        if t is not None:
            assert float(t.get("rms", "1")) < 0.1


@pytest.mark.usefixtures("file_based_db")
class TestPrivateerAPI(APITestBase):
    """API tests for privateer glycan validation."""

    task_name = "privateer"
    timeout = 180

    def test_glyco_4iid(self, glyco_4iid_pdb, glyco_4iid_mtz):
        """Test privateer with 4iid glycoprotein."""
        self.create_project("test_privateer")
        self.create_job()

        self.upload_file("inputData.XYZIN", glyco_4iid_pdb)
        self.upload_file_with_columns(
            "inputData.F_SIGF", glyco_4iid_mtz,
            column_labels="/*/*/[FP,SIGFP]"
        )

        # Use single-threaded for test stability (wrapper uses NUMTHREADS, not SINGLETHREADED)
        self.set_param("controlParameters.NUMTHREADS", 1)

        self.run_and_wait()

        # Privateer produces XML output and map coefficients
        self.assert_file_exists("program.xml")


@pytest.mark.usefixtures("file_based_db")
class TestEdstatsAPI(APITestBase):
    """API tests for edstats electron density statistics."""

    task_name = "edstats"
    timeout = 180

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test edstats with 8xfm structure and map coefficients."""
        self.create_project("test_edstats")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        # 2mFo-DFc map coefficients (FWT, PHWT)
        self.upload_file_with_columns(
            "inputData.FPHIIN1", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )
        # mFo-DFc difference map coefficients (DELFWT, PHDELWT)
        self.upload_file_with_columns(
            "inputData.FPHIIN2", mtz8xfm,
            column_labels="/*/*/[DELFWT,PHDELWT]"
        )

        # Set resolution limits (8xfm is ~1.3A resolution)
        self.set_param("inputData.RES_LOW", 50.0)
        self.set_param("inputData.RES_HIGH", 1.3)

        self.run_and_wait()

        # Edstats produces per-residue density statistics XML
        self.assert_file_exists("program.xml")
