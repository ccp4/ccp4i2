"""
API-driven tests for experimental phasing (EP) pipelines.

These tests verify the API workflow for:
- shelx: SHELX phasing pipeline
- crank2: Crank2 automated SAD/MAD phasing
- phaser_EP: Phaser experimental phasing

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import re
import pytest

from .base import APITestBase


@pytest.mark.usefixtures("file_based_db")
class TestShelxAPI(APITestBase):
    """API tests for SHELX phasing pipeline."""

    task_name = "shelx"
    timeout = 600  # EP pipelines can be slow

    def test_substructure_determination(self, gamma_mtz, gamma_asu_xml):
        """Test SHELX substructure determination only."""
        self.create_project("test_shelx_substrdet")
        self.create_job()

        # Upload anomalous data with column specification
        self.upload_file_with_columns(
            "inputData.F_SIGFanom", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.SEQIN", gamma_asu_xml)

        # Set phasing parameters (these are in inputData for shelx)
        self.set_param("inputData.ATOM_TYPE", "Xe")
        self.set_param("inputData.NUMBER_SUBSTRUCTURE", 2)
        self.set_param("inputData.WAVELENGTH", 1.54179)
        self.set_param("inputData.FPRIME", -0.79)
        self.set_param("inputData.FDPRIME", 7.36)
        self.set_param("inputData.END_PIPELINE", "substrdet")

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("PDBCUR.pdb")

    def test_gamma_sad(self, gamma_mtz, gamma_asu_xml):
        """Test full SHELX SAD phasing."""
        self.create_project("test_shelx_sad")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGFanom", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.SEQIN", gamma_asu_xml)

        self.set_param("inputData.ATOM_TYPE", "Xe")
        self.set_param("inputData.NUMBER_SUBSTRUCTURE", 2)
        self.set_param("inputData.WAVELENGTH", 1.54179)
        self.set_param("inputData.FPRIME", -0.79)
        self.set_param("inputData.FDPRIME", 7.36)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("n_part.pdb")
        self.validate_pdb("n_PDBCUR.pdb")
        self.validate_mtz("FPHOUT_2FOFC.mtz")
        self.validate_mtz("FPHOUT_DIFF.mtz")
        self.validate_mtz("FPHOUT_HL.mtz")
        self.validate_mtz("FREEROUT.mtz")

        # Check quality metrics from log
        log = self.read_output_file("log.txt")
        foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
        rworks = [float(x) for x in re.findall(r"R factor .* (0\.\d+)", log)]
        rfrees = [float(x) for x in re.findall(r"R-free .* (0\.\d+)", log)]

        assert max(foms) > 0.7
        assert min(rworks) < 0.22
        assert min(rfrees) < 0.25

    def test_gamma_siras(self, gamma_mtz, gamma_native_mtz, gamma_asu_xml):
        """Test SHELX SIRAS phasing with native data."""
        self.create_project("test_shelx_siras")
        self.create_job()

        # Anomalous derivative
        self.upload_file_with_columns(
            "inputData.F_SIGFanom", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        # Native data
        self.upload_file_with_columns(
            "inputData.F_SIGFnative", gamma_native_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.SEQIN", gamma_asu_xml)

        self.set_param("inputData.NATIVE", True)
        self.set_param("inputData.EXPTYPE", "SIRAS")
        self.set_param("inputData.ATOM_TYPE", "Xe")
        self.set_param("inputData.NUMBER_SUBSTRUCTURE", 2)
        self.set_param("inputData.WAVELENGTH", 1.54179)
        self.set_param("inputData.FPRIME", -0.79)
        self.set_param("inputData.FDPRIME", 7.36)

        self.run_and_wait()

        # Validate outputs
        self.validate_pdb("n_REFMAC5.pdb")
        self.validate_mtz("FPHOUT_HL.mtz")

        # Check quality metrics
        log = self.read_output_file("log.txt")
        foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
        assert max(foms) > 0.7


@pytest.mark.usefixtures("file_based_db")
class TestCrank2API(APITestBase):
    """API tests for Crank2 automated phasing pipeline."""

    task_name = "crank2"
    timeout = 600

    def test_gamma_dmfull(self, gamma_mtz, gamma_asu_xml):
        """Test Crank2 phasing with density modification."""
        self.create_project("test_crank2")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGFanom", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.SEQIN", gamma_asu_xml)

        self.set_param("inputData.NUMBER_SUBSTRUCTURE", 2)
        self.set_param("inputData.ATOM_TYPE", "Xe")
        self.set_param("inputData.FPRIME", -0.79)
        self.set_param("inputData.FDPRIME", 7.36)
        self.set_param("inputData.WAVELENGTH", 1.54179)
        self.set_param("inputData.END_PIPELINE", "dmfull")

        self.run_and_wait()

        # Validate outputs
        self.validate_mtz("FPHOUT_DIFFANOM.mtz")
        self.validate_mtz("FPHOUT_HL.mtz")
        self.validate_mtz("FPHOUT.mtz")
        self.validate_mtz("FREEROUT.mtz")

        # Check FOM from log
        log = self.read_output_file("log.txt")
        foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
        assert max(foms) > 0.7


@pytest.mark.usefixtures("file_based_db")
class TestPhaserEPAPI(APITestBase):
    """API tests for Phaser experimental phasing."""

    task_name = "phaser_EP"
    timeout = 300

    def test_gamma(self, gamma_mtz, gamma_heavy_atoms_pdb, gamma_asu_xml):
        """Test Phaser EP with known heavy atom positions."""
        self.create_project("test_phaser_ep")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )
        self.upload_file("inputData.XYZIN_PARTIAL", gamma_heavy_atoms_pdb)
        self.upload_file("inputData.XYZIN_HA", gamma_heavy_atoms_pdb)
        self.upload_file("inputData.ASUFILE", gamma_asu_xml)

        self.set_param("inputData.COMP_BY", "ASU")
        self.set_param("inputData.WAVELENGTH", 1.542)
        self.set_param("inputData.LLGC_CYCLES", 20)
        self.set_param("inputData.ELEMENTS", "Xe")
        self.set_param("controlParameters.RUNPARROT", False)

        self.run_and_wait()

        # Validate outputs
        self.validate_mtz("ABCDOUT_1.mtz")
        self.validate_mtz("ABCDOUT_2.mtz")
        self.validate_mtz("MAPOUT_1.mtz")
        self.validate_mtz("MAPOUT_2.mtz")
        self.validate_pdb("PHASER.1.pdb")
        self.validate_pdb("PHASER.1.hand.pdb")

        # Check heavy atom occupancies
        import gemmi
        model = gemmi.read_structure(str(self.get_job_directory() / "PHASER.1.pdb"))[0]
        occs = [cra.atom.occ for cra in model.all()]
        assert sum(occ > 0.11 for occ in occs) == 3
