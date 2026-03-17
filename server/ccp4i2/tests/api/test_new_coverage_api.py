"""
API-driven tests for tasks that previously had NO test coverage.

These tests are inferred from wrapper .def.xml and Python code:
- zanuda: Space group validation / refinement
- gesamt: Structure superposition
- ctruncate: Intensity truncation to amplitudes
- privateer: Glycan validation
- edstats: Electron density statistics
- pointless_reindexToMatch: Reindex reflections to match reference

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from .base import APITestBase, download, URLs


@pytest.mark.usefixtures("file_based_db")
class TestZanudaAPI(APITestBase):
    """API tests for zanuda space group validation."""

    task_name = "zanuda"
    timeout = 300

    @pytest.mark.skip(reason="Zanuda's internal refmac5 fails to parse 8xfm mmCIF model")
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

        # Zanuda produces a refined model and map coefficients
        self.assert_file_exists("XYZOUT.cif")
        self.assert_file_exists("FPHIOUT.mtz")
        self.assert_file_exists("DIFFPHIOUT.mtz")


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
class TestCtruncateAPI(APITestBase):
    """API tests for ctruncate intensity truncation."""

    task_name = "ctruncate"
    timeout = 120

    def test_gamma_intensities(self, gamma_mtz, demo_data_dir):
        """Test ctruncate with gamma anomalous intensities."""
        # gamma_mtz has Iplus/SIGIplus/Iminus/SIGIminus
        gamma_seq = demo_data_dir / "gamma" / "gamma.pir"
        if not gamma_seq.exists():
            pytest.skip(f"Sequence file not found: {gamma_seq}")

        self.create_project("test_ctruncate")
        self.create_job()

        self.upload_file("inputData.OBSIN", gamma_mtz)
        self.upload_file("inputData.SEQIN", str(gamma_seq))

        self.run_and_wait()

        # Check for truncated output
        self.assert_file_exists("program.xml")
        # ctruncate produces OBSOUT mini-MTZ with amplitudes
        job_dir = self.get_job_directory()
        mtz_files = list(job_dir.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from ctruncate: {job_dir}"


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


@pytest.mark.usefixtures("file_based_db")
class TestPointlessReindexAPI(APITestBase):
    """API tests for pointless_reindexToMatch."""

    task_name = "pointless_reindexToMatch"
    timeout = 120

    def test_reindex_to_coords(self, gamma_mtz, gamma_model_pdb):
        """Test reindexing gamma data to match coordinate reference."""
        self.create_project("test_pointless_reindex")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN_REF", gamma_model_pdb)
        self.set_param("controlParameters.REFERENCE", "XYZIN_REF")

        self.run_and_wait()

        # Check reindexed output
        job_dir = self.get_job_directory()
        mtz_files = list(job_dir.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from pointless: {job_dir}"
