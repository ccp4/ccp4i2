# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
API-driven tests for molecular replacement pipelines.

These tests verify the API workflow for:
- phaser_simple: Simple MR pipeline
- molrep_pipe: Molrep MR pipeline
- i2Dimple: Dimple pipeline
- mrbump_basic: MrBUMP automated MR
- molrep_selfrot: Self-rotation function analysis
- phaser_pipeline: Expert Phaser MR with multiple ensembles
- phaser_rnp_pipeline: Phaser MR with chain selections
- chainsaw: Model editing with alignment
- sculptor: Model trimming with alignment
- phaser_ensembler: Ensemble superposition
- findmyseq: Sequence identification from map
- molrep_den: Molrep density/Patterson-mode MR

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from ..base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestPhaserSimpleAPI(APITestBase):
    """API tests for phaser_simple MR pipeline."""

    task_name = "phaser_simple"
    timeout = 300  # Phaser can take a while

    def test_gamma_basic(self, gamma_mtz, gamma_model_pdb):
        """Test phaser_simple without post-processing."""
        self.create_project("test_phaser_simple_basic")
        self.create_job()

        # Upload input files with column specification for anomalous intensities
        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # Disable post-processing for faster test
        self.set_param("inputData.RUNREFMAC", False)
        self.set_param("inputData.RUNSHEETBEND", False)

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

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.set_param("inputData.RUNREFMAC", False)
        self.set_param("inputData.RUNSHEETBEND", True)

        self.run_and_wait()

        # Validate outputs including sheetbend
        self.validate_pdb("PHASER.1.pdb")
        self.validate_pdb("XYZOUT_SHEETBEND.pdb")
        self.validate_mtz("PHASEOUT_1.mtz")

    def test_gamma_full(self, gamma_mtz, gamma_model_pdb):
        """Test phaser_simple with full post-processing (refmac + sheetbend)."""
        self.create_project("test_phaser_full")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        # Can add multiple search models
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        # Full pipeline with refmac and sheetbend
        self.set_param("inputData.RUNREFMAC", True)
        self.set_param("inputData.RUNSHEETBEND", True)

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
        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN", rnase_model_pdb)

        # Lower resolution to speed up
        self.set_param("inputData.RESOLUTION_HIGH", 5.0)

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

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", gamma_freer_mtz,
            column_labels="/*/*/[FREER]"
        )
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

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_native_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
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

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", gamma_freer_mtz,
            column_labels="/*/*/[FREER]"
        )
        self.upload_file("inputData.ASUIN", gamma_asu_xml)

        # Limit search for faster testing
        self.set_param("modelParameters.MRMAX", 5)
        self.set_param("modelParameters.REDUNDANCYLEVEL", 95)
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


@pytest.mark.usefixtures("file_based_db")
class TestMolrepSelfrotAPI(APITestBase):
    """API tests for molrep_selfrot self-rotation function."""

    task_name = "molrep_selfrot"
    timeout = 120

    def test_1h1s(self, mtz1h1s):
        """Test self-rotation function with 1h1s reflection data."""
        self.create_project("test_molrep_selfrot")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz1h1s,
            column_labels="/*/*/[FP,SIGFP]"
        )

        self.run_and_wait()

        # Check PostScript output
        job_dir = self.get_job_directory()
        ps_files = list(job_dir.glob("*.ps"))
        assert len(ps_files) >= 1, f"No PS output: {list(job_dir.iterdir())}"

        # Check Patterson peaks in program.xml
        xml = self.read_program_xml()
        peaks = xml.findall(".//Patterson/Peak")
        assert len(peaks) >= 10, f"Expected >= 10 Patterson peaks, got {len(peaks)}"

        # Verify peak structure
        first_peak = peaks[0]
        for field in ("No", "Xfrac", "Yfrac", "Zfrac", "Dens", "Dens_sigma"):
            elem = first_peak.find(field)
            assert elem is not None, f"Missing '{field}' in Patterson peak"


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


@pytest.mark.usefixtures("file_based_db")
class TestChainsawAPI(APITestBase):
    """API tests for chainsaw model editing."""

    task_name = "chainsaw"
    timeout = 120

    def test_gamma(self, gamma_model_pdb, demo_data_dir):
        """Test chainsaw with self-alignment."""
        alignment = demo_data_dir / "gamma" / "gamma_self_alignment.fas"
        if not alignment.exists():
            pytest.skip(f"Alignment file not found: {alignment}")

        self.create_project("test_chainsaw")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.upload_file("inputData.ALIGNIN", str(alignment))

        self.run_and_wait()

        self.validate_pdb("XYZOUT.pdb")


@pytest.mark.usefixtures("file_based_db")
class TestSculptorAPI(APITestBase):
    """API tests for sculptor model trimming."""

    task_name = "sculptor"
    timeout = 120

    def test_gamma(self, gamma_model_pdb, demo_data_dir):
        """Test sculptor with gamma FASTA alignment."""
        aln_path = demo_data_dir / "gamma" / "gamma_self_alignment.fas"
        if not aln_path.exists():
            pytest.skip(f"Alignment file not found: {aln_path}")

        self.create_project("test_sculptor")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.set_param("inputData.ALIGNMENTORSEQUENCEIN", "ALIGNMENT")
        self.upload_file("inputData.ALIGNIN", str(aln_path))
        self.set_param("controlParameters.TARGETINDEX", 0)

        self.run_and_wait()

        # sculptor outputs a list of PDB files
        job_dir = self.get_job_directory()
        pdb_files = list(job_dir.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job_dir.iterdir())}"


@pytest.mark.usefixtures("file_based_db")
class TestPhaserEnsemblerAPI(APITestBase):
    """API tests for phaser_ensembler ensemble superposition."""

    task_name = "phaser_ensembler"
    timeout = 120

    def test_gamma_two_models(self, gamma_model_pdb):
        """Test phaser.ensembler with two copies of gamma model."""
        self.create_project("test_phaser_ensembler")
        self.create_job()

        # XYZIN_LIST is a CList of CPdbDataFile
        self.set_param("inputData.XYZIN_LIST", [{}, {}])
        self.upload_file("inputData.XYZIN_LIST[0]", gamma_model_pdb)
        self.upload_file("inputData.XYZIN_LIST[1]", gamma_model_pdb)

        self.run_and_wait()

        job_dir = self.get_job_directory()
        pdb_files = list(job_dir.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job_dir.iterdir())}"


@pytest.mark.usefixtures("file_based_db")
class TestFindmyseqAPI(APITestBase):
    """API tests for findmyseq sequence identification."""

    task_name = "findmyseq"
    timeout = 300

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test findmyseq with 8xfm data."""
        self.create_project("test_findmyseq")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.upload_file_with_columns(
            "inputData.F_SIGF", mtz8xfm,
            column_labels="/*/*/[FP,SIGFP]"
        )
        self.upload_file_with_columns(
            "inputData.FPHI", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )

        self.run_and_wait()

        self.assert_file_exists("SEQOUT.fasta")
        content = self.read_output_file("SEQOUT.fasta")
        assert len(content) > 0, "Empty sequence output"


@pytest.mark.usefixtures("file_based_db")
class TestMolrepDenAPI(APITestBase):
    """API tests for molrep_den density-mode MR."""

    task_name = "molrep_den"
    timeout = 300

    @pytest.mark.skip(reason="Needs a different model from the map — self-search fails in molrep")
    def test_8xfm_density(self, cif8xfm, mtz8xfm):
        """Test molrep density search with 8xfm data."""
        self.create_project("test_molrep_den")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_PHI_MAP", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )
        self.upload_file("inputData.XYZIN", cif8xfm)

        self.run_and_wait()

        # Molrep density search may not find a solution (XYZOUT not always
        # produced) but should complete and write program.xml
        self.assert_file_exists("program.xml")
