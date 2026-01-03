"""
API-driven tests for utility tasks.

These tests verify the API workflow for:
- LidiaAcedrgNew: Ligand dictionary generation
- validate_protein: Structure validation
- AUSPEX: Data quality analysis
- csymmatch: Coordinate matching
- coordinate_selector: Atom selection
- coot_find_waters: Water finding
- coot_rsr_morph: Real-space refinement morphing

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import os
import pytest

from .base import APITestBase, download, URLs


@pytest.mark.usefixtures("file_based_db")
class TestAcedrgAPI(APITestBase):
    """API tests for LidiaAcedrgNew ligand dictionary generation."""

    task_name = "LidiaAcedrgNew"
    timeout = 120

    def test_from_smiles(self):
        """Test acedrg from SMILES string."""
        self.create_project("test_acedrg_smiles")
        self.create_job()

        # These parameters are in inputData for LidiaAcedrgNew
        self.set_param("inputData.MOLSMILESORSKETCH", "SMILES")
        self.set_param("inputData.TLC", "HCA")
        # SMILES for hydroxycinnamalaldehyde
        self.set_param("inputData.SMILESIN", "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4")

        self.run_and_wait()

        # Validate outputs - check CIF file
        cif_path = self.assert_file_exists("HCA.cif")

        # Verify it's a valid chemical compound CIF
        import gemmi
        doc = gemmi.cif.read(str(cif_path))
        gemmi.make_chemcomp_from_block(doc[-1])

    def test_from_cif_rcsb(self):
        """Test acedrg from RCSB ligand CIF."""
        with download(URLs.rcsb_ligand_cif("A1LU6")) as cif_path:
            self.create_project("test_acedrg_cif")
            self.create_job()

            self.set_param("inputData.MOLSMILESORSKETCH", "DICT")
            self.set_param("inputData.TLC", "A1LU6")
            self.upload_file("inputData.DICTIN2", cif_path)

            self.run_and_wait()

            # Validate CIF output
            out_cif = self.assert_file_exists("A1LU6.cif")
            import gemmi
            doc = gemmi.cif.read(str(out_cif))
            gemmi.make_chemcomp_from_block(doc[-1])

    def test_from_monomer_library(self):
        """Test acedrg from CCP4 monomer library."""
        clibd_mon = os.environ.get("CLIBD_MON")
        if not clibd_mon:
            pytest.skip("CLIBD_MON not set")

        cif_path = os.path.join(clibd_mon, "a", "A1LU6.cif")
        if not os.path.exists(cif_path):
            pytest.skip(f"Monomer library CIF not found: {cif_path}")

        self.create_project("test_acedrg_monlib")
        self.create_job()

        self.set_param("inputData.MOLSMILESORSKETCH", "DICT")
        self.set_param("inputData.TLC", "A1LU6")
        self.upload_file("inputData.DICTIN2", cif_path)
        self.set_param("controlParameters.USE_COORD", True)

        self.run_and_wait()

        self.assert_file_exists("A1LU6.cif")


@pytest.mark.usefixtures("file_based_db")
class TestValidateProteinAPI(APITestBase):
    """API tests for validate_protein structure validation."""

    task_name = "validate_protein"
    timeout = 180

    def test_7beq_basic(self, cif7beq):
        """Test basic validation without MolProbity."""
        self.create_project("test_validate_basic")
        self.create_job()

        self.upload_file("inputData.XYZIN_1", cif7beq)

        self.run_and_wait()

        # Check program.xml has validation results
        xml = self.read_program_xml()
        assert xml.find(".//Chain_count") is not None
        assert xml.find(".//Iris") is not None
        assert xml.find(".//B_factors") is not None
        assert xml.find(".//Ramachandran") is not None


@pytest.mark.usefixtures("file_based_db")
class TestAuspexAPI(APITestBase):
    """API tests for AUSPEX data quality analysis."""

    task_name = "AUSPEX"
    timeout = 120

    def test_gamma(self, gamma_mtz):
        """Test AUSPEX with gamma data."""
        self.create_project("test_auspex")
        self.create_job()

        self.upload_file_with_columns(
            "inputData.F_SIGF", gamma_mtz,
            column_labels="/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
        )

        self.run_and_wait()

        # Check output plots exist
        self.assert_file_exists("F_plot.png")
        self.assert_file_exists("FSigF_plot.png")
        self.assert_file_exists("I_plot.png")
        self.assert_file_exists("ISigI_plot.png")
        self.assert_file_exists("SigF_plot.png")
        self.assert_file_exists("SigI_plot.png")


@pytest.mark.usefixtures("file_based_db")
class TestCsymmatchAPI(APITestBase):
    """API tests for csymmatch coordinate matching."""

    task_name = "csymmatch"
    timeout = 60

    def test_8xfm(self, cif8xfm):
        """Test csymmatch with same structure (identity)."""
        self.create_project("test_csymmatch")
        self.create_job()

        # Same structure as query and target
        self.upload_file("inputData.XYZIN_QUERY", cif8xfm)
        self.upload_file("inputData.XYZIN_TARGET", cif8xfm)

        self.run_and_wait()

        # Validate outputs
        self.assert_file_exists("XYZOUT.cif")

        # Should have identity transformation (origin shift of 0,0,0)
        xml = self.read_program_xml()
        change_of_origin = xml.find(".//ChangeOfOrigin")
        if change_of_origin is not None:
            # Output format varies - check for zeros in any format
            text = change_of_origin.text.strip()
            assert "0" in text, f"Expected zero origin shift, got: {text}"


@pytest.mark.usefixtures("file_based_db")
class TestCoordinateSelectorAPI(APITestBase):
    """API tests for coordinate_selector atom selection."""

    task_name = "coordinate_selector"
    timeout = 60

    def test_8xfm(self, cif8xfm):
        """Test selecting atoms from 8xfm."""
        self.create_project("test_coord_selector")
        self.create_job()

        # Upload with selection
        self.upload_file("inputData.XYZIN", cif8xfm)
        # Select chain A using mmdb CID syntax (bare identifier = chain)
        self.set_param("inputData.XYZIN.selection.text", "A")

        self.run_and_wait()

        # Validate output
        self.assert_file_exists("XYZOUT.cif")


@pytest.mark.usefixtures("file_based_db")
class TestCootFindWatersAPI(APITestBase):
    """API tests for coot_find_waters water finding."""

    task_name = "coot_find_waters"
    timeout = 180

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test finding waters in 8xfm."""
        self.create_project("test_find_waters")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.upload_file_with_columns(
            "inputData.FPHIIN", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )

        self.run_and_wait()

        # Validate output has waters
        self.assert_file_exists("XYZOUT.cif")

        # Check for water residues
        import gemmi
        st = gemmi.read_structure(str(self.get_job_directory() / "XYZOUT.cif"))
        has_water = False
        for model in st:
            for chain in model:
                for res in chain:
                    if res.name == "HOH":
                        has_water = True
                        break
        assert has_water, "No water molecules found in output"


@pytest.mark.usefixtures("file_based_db")
class TestCootRsrMorphAPI(APITestBase):
    """API tests for coot_rsr_morph real-space refinement morphing."""

    task_name = "coot_rsr_morph"
    timeout = 180

    def test_8xfm(self, cif8xfm, mtz8xfm):
        """Test RSR morphing with 8xfm."""
        self.create_project("test_rsr_morph")
        self.create_job()

        self.upload_file("inputData.XYZIN", cif8xfm)
        self.upload_file_with_columns(
            "inputData.FPHIIN", mtz8xfm,
            column_labels="/*/*/[FWT,PHWT]"
        )

        self.run_and_wait()

        # Validate output
        self.assert_file_exists("XYZOUT.cif")


@pytest.mark.usefixtures("file_based_db")
class TestMrparseAPI(APITestBase):
    """API tests for mrparse homolog search."""

    task_name = "mrparse"
    timeout = 300  # Network dependent

    def test_gamma(self, demo_data_dir):
        """Test mrparse with gamma sequence."""
        pir_path = demo_data_dir / "gamma" / "gamma.pir"
        if not pir_path.exists():
            pytest.skip(f"PIR file not found: {pir_path}")

        self.create_project("test_mrparse")
        self.create_job()

        self.upload_file("inputData.SEQIN", str(pir_path))
        # These parameters are in 'options' container for mrparse
        self.set_param("options.DATABASE", "PDB")
        self.set_param("options.USEAPI", True)

        self.run_and_wait()

        # Check for homolog results
        job_dir = self.get_job_directory()
        homologs_path = job_dir / "mrparse_0" / "homologs.json"
        if homologs_path.exists():
            import json
            with homologs_path.open() as f:
                homologs = json.load(f)
                # Should find itself with perfect identity
                identities = [h.get("seq_ident", 0) for h in homologs]
                assert any(ident == 1.0 for ident in identities)


@pytest.mark.usefixtures("file_based_db")
class TestProvideAsuContentsAPI(APITestBase):
    """API tests for ProvideAsuContents ASU definition."""

    task_name = "ProvideAsuContents"
    timeout = 60

    def test_beta_blip(self, demo_data_dir):
        """Test defining ASU contents."""
        self.create_project("test_asu_contents")
        self.create_job()

        # Set up ASU content - this mimics the complex parameter structure
        beta_blip_dir = demo_data_dir / "beta_blip"
        if not beta_blip_dir.exists():
            pytest.skip("beta_blip demo data not found")

        # Read actual sequence from demo data
        beta_seq_file = beta_blip_dir / "beta.seq"
        if not beta_seq_file.exists():
            pytest.skip("beta.seq not found in demo data")

        # Parse the FASTA-like sequence file
        lines = beta_seq_file.read_text().strip().split('\n')
        sequence = ''.join(line for line in lines if not line.startswith('>'))

        # Add ASU content items with real sequence
        self.set_param("inputData.ASU_CONTENT", [{
            "sequence": sequence,
            "nCopies": 1,
            "name": "beta",
            "description": "Beta-lactamase",
            "polymerType": "PROTEIN"
        }])

        self.run_and_wait()

        # Check output XML
        self.assert_file_exists("ASUCONTENTFILE.asu.xml")
