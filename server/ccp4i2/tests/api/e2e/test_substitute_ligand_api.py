"""
API-driven tests for SubstituteLigand pipeline.

These tests verify the API workflow for:
- SubstituteLigand with LIGANDAS=NONE (refinement only, no ligand fitting)
- SubstituteLigand with SMILES input (full pipeline with ligand fitting)

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestSubstituteLigandAPI(APITestBase):
    """API tests for SubstituteLigand pipeline."""

    task_name = "SubstituteLigand"
    timeout = 300

    def test_no_ligand(self, mdm2_model_cif, mdm2_unmerged_mtz):
        """Test SubstituteLigand with LIGANDAS=NONE (refinement only)."""
        self.create_project("test_substitute_ligand_no_ligand")
        self.create_job()

        # Upload input files
        self.upload_file("inputData.XYZIN", mdm2_model_cif)
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", mdm2_unmerged_mtz)

        # Set parameters
        self.set_param("controlParameters.LIGANDAS", "NONE")
        self.set_param("inputData.PIPELINE", "DIMPLE")

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("XYZOUT.pdb")
        self.validate_mtz("F_SIGF_OUT.mtz")
        self.validate_mtz("FREERFLAG_OUT.mtz")
        self.validate_mtz("DIFFPHIOUT.mtz")
        self.validate_mtz("ANOMFPHIOUT.mtz")
        self.assert_file_exists("selected_atoms.pdb")

        # Check R-factors
        xml = self.read_program_xml()
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23, f"Final R-work {rworks[-1]} >= 0.23"
        assert rfrees[-1] < 0.25, f"Final R-free {rfrees[-1]} >= 0.25"

    def test_with_smiles(self, mdm2_model_cif, mdm2_unmerged_mtz):
        """Test SubstituteLigand with SMILES input (full pipeline)."""
        self.create_project("test_substitute_ligand_smiles")
        self.create_job()

        # Upload input files
        self.upload_file("inputData.XYZIN", mdm2_model_cif)
        self.add_list_item("inputData.UNMERGEDFILES")
        self.upload_file("inputData.UNMERGEDFILES[0].file", mdm2_unmerged_mtz)

        # Set parameters - use MDM2 ligand SMILES
        mdm2_ligand_smiles = "CC(C)OC1=C(C=CC(=C1)OC)C2=NC(C(N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)C"
        self.set_param("controlParameters.LIGANDAS", "SMILES")
        self.set_param("inputData.SMILESIN", mdm2_ligand_smiles)
        self.set_param("inputData.PIPELINE", "DIMPLE")

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("XYZOUT.pdb")
        self.assert_file_exists("DICTOUT.cif")
        self.validate_mtz("F_SIGF_OUT.mtz")
        self.validate_mtz("FREERFLAG_OUT.mtz")
        self.validate_mtz("DIFFPHIOUT.mtz")
        self.assert_file_exists("selected_atoms.pdb")

        # Check R-factors
        xml = self.read_program_xml()
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23, f"Final R-work {rworks[-1]} >= 0.23"
        assert rfrees[-1] < 0.25, f"Final R-free {rfrees[-1]} >= 0.25"
