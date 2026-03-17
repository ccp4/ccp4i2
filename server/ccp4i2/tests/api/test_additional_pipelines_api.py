"""
API-driven tests for tasks that previously only had i2run coverage.

These tests translate existing i2run tests to the API workflow for:
- comit: Omit map calculation
- metalCoord: Metal coordination restraint generation
- SubtractNative: Native subtraction map
- editbfac: AlphaFold B-factor editing
- splitMtz: MTZ column splitting
- scaleit: Derivative vs native scaling
- phaser_pipeline: Expert Phaser MR with multiple ensembles
- phaser_rnp_pipeline: Phaser MR with chain selections

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import json
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from .base import APITestBase, download, URLs


@pytest.mark.usefixtures("file_based_db")
class TestComitAPI(APITestBase):
    """API tests for comit omit map calculation."""

    task_name = "comit"
    timeout = 120

    def test_8xfm(self, mtz8xfm):
        """Test comit with 8xfm data."""
        self.create_project("test_comit")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz8xfm,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.F_PHI_IN", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )

        self.run_and_wait()

        self.validate_mtz("F_PHI_OUT.mtz")


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


@pytest.mark.usefixtures("file_based_db")
class TestSubtractNativeAPI(APITestBase):
    """API tests for SubtractNative map calculation."""

    task_name = "SubtractNative"
    timeout = 120

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test native subtraction with 8xfm data."""
        self.create_project("test_subtract_native")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.upload_file_with_columns(
            "inputData.MAPIN", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )

        self.run_and_wait()

        self.assert_file_exists("MAPOUT.map")


@pytest.mark.usefixtures("file_based_db")
class TestEditbfacAPI(APITestBase):
    """API tests for editbfac B-factor editing."""

    task_name = "editbfac"
    timeout = 120

    def test_alphafold_pdb(self):
        """Test editbfac with AlphaFold model."""
        # Use AlphaFold API to get the dynamic PDB URL (like the i2run test)
        api_url = "https://alphafold.ebi.ac.uk/api/prediction/Q8W3K0"
        with download(api_url) as api_path:
            with open(api_path, encoding="utf-8") as f:
                info = json.load(f)[0]
            with download(info["pdbUrl"]) as pdb_path:
                self.create_project("test_editbfac")
                self.create_job()

                self.upload_file("inputData.XYZIN", pdb_path)

                self.run_and_wait()

                self.assert_file_exists("converted_model.pdb")
                self.assert_file_exists("converted_model_chainA1.pdb")


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
class TestPhaserPipelineAPI(APITestBase):
    """API tests for phaser_pipeline expert MR."""

    task_name = "phaser_pipeline"
    timeout = 600

    def test_beta_blip(self, beta_blip_mtz, beta_pdb, blip_pdb):
        """Test phaser_pipeline with beta+blip two-ensemble MR."""
        self.create_project("test_phaser_pipeline")
        self.create_job()

        # Upload F/SIGF data with correct column names for beta_blip MTZ
        self.upload_file_with_columns(
            "inputData.F_SIGF", beta_blip_mtz,
            column_labels="/*/*/[Fobs,Sigma]"
        )

        # Set F or I mode
        self.set_param("inputData.F_OR_I", "F")
        self.set_param("inputData.COMP_BY", "DEFAULT")

        # Add two ensembles with correct nested structure
        # ENSEMBLES is a CEnsembleList: each item has use, pdbItemList
        self.set_param("inputData.ENSEMBLES", [{}, {}])
        self.set_param("inputData.ENSEMBLES[0].use", True)
        self.set_param("inputData.ENSEMBLES[0].pdbItemList", [{}])
        self.set_param("inputData.ENSEMBLES[0].pdbItemList[0].identity_to_target", 0.9)
        self.upload_file("inputData.ENSEMBLES[0].pdbItemList[0].structure", beta_pdb)

        self.set_param("inputData.ENSEMBLES[1].use", True)
        self.set_param("inputData.ENSEMBLES[1].pdbItemList", [{}])
        self.set_param("inputData.ENSEMBLES[1].pdbItemList[0].identity_to_target", 0.9)
        self.upload_file("inputData.ENSEMBLES[1].pdbItemList[0].structure", blip_pdb)

        # Disable post-processing for speed
        self.set_param("inputData.RUNREFMAC", False)
        self.set_param("inputData.RUNSHEETBEND", False)

        self.run_and_wait()

        # Check for Phaser output
        job_dir = self.get_job_directory()
        phaser_pdbs = list(job_dir.glob("PHASER.*.pdb")) + list(job_dir.glob("PHASER.*.cif"))
        assert len(phaser_pdbs) >= 1, f"No Phaser output found in {job_dir}"

        # Check LLG
        xml = self.read_program_xml()
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        if llgs:
            assert max(llgs) > 1000


@pytest.mark.usefixtures("file_based_db")
class TestPhaserRnpPipelineAPI(APITestBase):
    """API tests for phaser_rnp_pipeline MR with chain selections."""

    task_name = "phaser_rnp_pipeline"
    timeout = 600

    def test_1h1s_chain_selections(self, pdb1h1s, mtz1h1s):
        """Test phaser_rnp_pipeline with chain selections from 1h1s."""
        self.create_project("test_phaser_rnp")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz1h1s,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file("inputData.XYZIN_PARENT", pdb1h1s)

        # Add chain selections
        self.set_param("inputData.SELECTIONS", [
            {"text": "A/"},
            {"text": "B/"},
            {"text": "C/"},
            {"text": "D/"},
        ])

        # Disable post-processing for speed
        self.set_param("inputData.RUNREFMAC", False)
        self.set_param("inputData.RUNSHEETBEND", False)

        self.run_and_wait()

        # Check for Phaser output
        job_dir = self.get_job_directory()
        phaser_outputs = list(job_dir.glob("PHASER.*"))
        assert len(phaser_outputs) >= 1, f"No Phaser output found in {job_dir}"
