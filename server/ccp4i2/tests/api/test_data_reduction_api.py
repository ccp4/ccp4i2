"""
API-driven tests for data reduction pipelines.

These tests verify the API workflow for:
- aimless_pipe: Data scaling and merging
- import_merged: Import merged reflections
- freerflag: Generate free R flags
- ctruncate: Intensity truncation to amplitudes
- pointless_reindexToMatch: Reindex reflections to match reference
- chltofom: HL coefficients to PHI/FOM conversion
- splitMtz: MTZ column splitting
- scaleit: Derivative vs native scaling
- matthews: Matthews coefficient estimation
- AUSPEX: Data quality analysis plots

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

    def test_mdm2_mmcif_reference(self, mdm2_unmerged_mtz, mdm2_model_cif):
        """Test aimless_pipe with mmCIF coordinate reference for space group."""
        self.create_project("test_aimless_mdm2_mmcif_ref")
        self.create_job()

        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", mdm2_unmerged_mtz)
        self.upload_file("inputData.XYZIN_REF", mdm2_model_cif)

        self.set_param("controlParameters.MODE", "MATCH")
        self.set_param("controlParameters.REFERENCE_DATASET", "XYZ")

        self.run_and_wait()

        # Validate outputs
        self.assert_valid_mtz_output("FREEROUT")
        digest = self.assert_valid_mtz_output("HKLOUT[0]")

        # Space group should match the 4hg7 reference (deposited as P 65 2 2)
        assert digest.get('spaceGroup') == "P 65 2 2"
        assert digest.get('highRes') == approx(1.25, rel=0.05)


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


@pytest.mark.usefixtures("file_based_db")
class TestChltofomAPI(APITestBase):
    """API tests for chltofom HL to PHI/FOM conversion."""

    task_name = "chltofom"
    timeout = 60

    def test_hl_to_phifom(self):
        """Test HL to PHI/FOM conversion with gamma initial phases."""
        self.create_project("test_chltofom")
        self.create_job()

        from ccp4i2.tests.i2run.utils import demoData
        self.upload_file("inputData.HKLIN", demoData("gamma", "initial_phases.mtz"))

        self.run_and_wait()

        # Check output MTZ exists with PHI/FOM columns
        import gemmi
        job_dir = self.get_job_directory()
        hklout = job_dir / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job_dir.iterdir())}"
        mtz = gemmi.read_mtz_file(str(hklout))
        labels = [c.label for c in mtz.columns]
        assert "PHI" in labels, f"No PHI column: {labels}"
        assert "FOM" in labels, f"No FOM column: {labels}"


@pytest.mark.usefixtures("file_based_db")
class TestSplitMtzAPI(APITestBase):
    """API tests for splitMtz column splitting."""

    task_name = "splitMtz"
    timeout = 60

    def test_gamma_abcd(self, gamma_initial_phases_mtz):
        """Test splitting ABCD columns from gamma initial phases."""
        self.create_project("test_split_mtz_abcd")
        self.create_job()

        self.upload_file("inputData.HKLIN", gamma_initial_phases_mtz)

        # Set up column group via set_param with structured data
        self.set_param("inputData.USERCOLUMNGROUP", {
            "columnGroupType": "Phs",
            "contentFlag": 1,
            "dataset": "ds1",
        })
        # Set the column list entries
        columns = [
            {"columnLabel": "HLA", "columnType": "A", "dataset": "ds1"},
            {"columnLabel": "HLB", "columnType": "A", "dataset": "ds1"},
            {"columnLabel": "HLC", "columnType": "A", "dataset": "ds1"},
            {"columnLabel": "HLD", "columnType": "A", "dataset": "ds1"},
        ]
        self.set_param("inputData.USERCOLUMNGROUP.columnList", columns)

        self.run_and_wait()

        # Check output MTZ exists
        job_dir = self.get_job_directory()
        output_files = list(job_dir.glob("*.mtz"))
        assert len(output_files) >= 1, f"No MTZ output files found in {job_dir}"


@pytest.mark.usefixtures("file_based_db")
class TestScaleitAPI(APITestBase):
    """API tests for scaleit derivative scaling."""

    task_name = "scaleit"
    timeout = 120

    def test_gamma_native_vs_xe(self, gamma_native_mtz, gamma_mtz):
        """Test scaling Xe derivative against native."""
        self.create_project("test_scaleit")
        self.create_job()

        # Create list with two empty items, then upload to each
        self.set_param("inputData.MERGEDFILES", [{}, {}])
        self.upload_file("inputData.MERGEDFILES[0]", gamma_native_mtz)
        self.upload_file("inputData.MERGEDFILES[1]", gamma_mtz)

        self.run_and_wait()

        # Validate output
        self.validate_mtz("hklout.mtz")

        # Check program.xml for scaleit results
        xml = self.read_program_xml()
        nderiv = xml.find(".//SCALEITLOG/Nderivatives")
        assert nderiv is not None, "Nderivatives missing from program.xml"
        assert int(nderiv.text) == 1


@pytest.mark.usefixtures("file_based_db")
class TestMatthewsAPI(APITestBase):
    """API tests for matthews coefficient estimation."""

    task_name = "matthews"
    timeout = 60

    def test_gamma_nres(self, gamma_mtz):
        """Test Matthews coefficient from number of residues."""
        self.create_project("test_matthews_nres")
        self.create_job()

        self.upload_file("inputData.HKLIN", gamma_mtz)
        self.set_param("inputData.MODE", "nres")
        self.set_param("inputData.NRES", 120)

        self.run_and_wait()

        xml = self.read_program_xml()
        compositions = xml.findall(".//matthewsCompositions/composition")
        assert len(compositions) > 0, "No matthews compositions in output"

        summary = xml.find(".//summary")
        assert summary is not None, "No summary element"
        solvent = float(summary.find("solventContent").text)
        assert 20 < solvent < 90, f"Unreasonable solvent content: {solvent}%"


@pytest.mark.usefixtures("file_based_db")
class TestAuspexAPI(APITestBase):
    """API tests for AUSPEX data quality analysis."""

    task_name = "AUSPEX"
    timeout = 120

    def test_gamma(self, gamma_mtz):
        """Test AUSPEX with gamma anomalous intensities."""
        self.create_project("test_auspex")
        self.create_job()

        self.upload_file("inputData.F_SIGF", gamma_mtz)

        self.run_and_wait()

        for plot in ["F", "FSigF", "I", "ISigI", "SigF", "SigI"]:
            self.assert_file_exists(f"{plot}_plot.png")
