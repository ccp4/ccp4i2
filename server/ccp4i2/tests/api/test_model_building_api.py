"""
API-driven tests for model building pipelines.

These tests verify the API workflow for:
- modelcraft: Modelcraft automated model building
- parrot: Parrot density modification

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import json
import pytest

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestModelcraftAPI(APITestBase):
    """API tests for Modelcraft automated model building."""

    task_name = "modelcraft"
    timeout = 900  # Model building can be slow

    def test_gamma_ep(self, gamma_mtz, gamma_freer_mtz, gamma_phases_mtz, gamma_asu_xml):
        """Test Modelcraft from experimental phases."""
        self.create_project("test_modelcraft_ep")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file_with_columns(
            "inputData.FREERFLAG", gamma_freer_mtz,
            column_labels="/*/*/[FREER]"
        )
        self.upload_file_with_columns(
            "inputData.PHASES", gamma_phases_mtz,
            column_labels="/*/*/[HLA,HLB,HLC,HLD]"
        )
        self.upload_file("inputData.ASUIN", gamma_asu_xml)

        self.set_param("controlParameters.USE_MODEL_PHASES", False)
        self.set_param("controlParameters.CYCLES", 2)

        self.run_and_wait()

        # Validate outputs
        self.assert_file_exists("XYZOUT.cif")
        self.validate_mtz("ABCDOUT.mtz")
        self.validate_mtz("DIFFPHIOUT.mtz")
        self.validate_mtz("FPHIOUT.mtz")

        # Check modelcraft results
        results_path = self.get_job_directory() / "modelcraft" / "modelcraft.json"
        if results_path.exists():
            with results_path.open() as f:
                results = json.load(f)
                assert results["final"]["r_free"] < 0.35

    def test_8xfm(self, cif8xfm, mtz8xfm, seq8xfm):
        """Test Modelcraft from structure + data."""
        pytest.skip(
            "ASUIN is a required input for modelcraft. "
            "Cannot upload sequence to ASUIN.seqFile via parameter API - "
            "ASUIN must be a complete ASU XML file. "
            "Use ProvideAsuContents task to create ASUIN from sequence first."
        )


@pytest.mark.usefixtures("file_based_db")
class TestParrotAPI(APITestBase):
    """API tests for Parrot density modification."""

    task_name = "parrot"
    timeout = 180

    def test_gamma(self, gamma_mtz, gamma_phases_mtz, gamma_asu_xml):
        """Test Parrot density modification."""
        self.create_project("test_parrot")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file_with_columns(
            "inputData.ABCD", gamma_phases_mtz,
            column_labels="/*/*/[HLA,HLB,HLC,HLD]"
        )
        self.upload_file("inputData.ASUIN", gamma_asu_xml)

        self.run_and_wait()

        # Validate outputs
        self.validate_mtz("ABCDOUT.mtz")
        self.validate_mtz("FPHIOUT.mtz")

        # Check FOM improvement
        xml = self.read_program_xml()
        foms = [float(e.text) for e in xml.findall(".//MeanFOM")]
        assert max(foms) > 0.794


@pytest.mark.usefixtures("file_based_db")
class TestShelxeMRAPI(APITestBase):
    """API tests for SHELXE molecular replacement extension."""

    task_name = "shelxeMR"
    timeout = 300

    def test_gamma(self, gamma_mtz, gamma_freer_mtz, gamma_model_pdb):
        """Test SHELXE MR with gamma data."""
        self.create_project("test_shelxe_mr")
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

        self.set_param("inputData.NTCYCLES", 5)
        self.set_param("inputData.NMCYCLES", 5)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("shelxrun.pdb")
        self.validate_mtz("FPHIOUT.mtz")

        # Check CC from program.xml
        xml = self.read_program_xml()
        best_cc_elem = xml.find(".//BestCC")
        if best_cc_elem is not None:
            best_cc = float(best_cc_elem.text)
            assert best_cc > 43


@pytest.mark.usefixtures("file_based_db")
class TestAcornAPI(APITestBase):
    """API tests for ACORN density improvement."""

    task_name = "acorn"
    timeout = 180

    def test_gamma(self, gamma_mtz, gamma_model_pdb):
        """Test ACORN with gamma data."""
        self.create_project("test_acorn")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.run_and_wait()

        # Validate outputs
        self.validate_mtz("PHSOUT.mtz")

        # Check correlation coefficient
        xml = self.read_program_xml()
        final_cc = xml.find(".//FinalCC")
        if final_cc is not None:
            assert float(final_cc.text) > 0.1
