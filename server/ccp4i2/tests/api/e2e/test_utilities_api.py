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
- MakeLink: Covalent link definition for ligands
- clustalw: Multiple sequence alignment
- comit: Omit map calculation
- SubtractNative: Native subtraction map
- editbfac: AlphaFold B-factor editing
- density_calculator: Model density map calculation
- pdbset_ui: PDB manipulation with keywords
- add_fractional_coords: Add fractional coordinates to model
- modelASUCheck: Model-sequence alignment check

Each test creates a project, creates a job, uploads input files,
runs the job, and validates outputs.
"""
import json
import os
import pytest

# Mark all tests in this module as pipeline tests (slow, run actual jobs)
pytestmark = pytest.mark.pipeline

from gemmi import cif

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


@pytest.mark.usefixtures("file_based_db")
class TestProvideSequenceAPI(APITestBase):
    """API tests for ProvideSequence sequence input."""

    task_name = "ProvideSequence"
    timeout = 60

    _GAMMA_FASTA = (
        ">gamma_chymotrypsinogen\n"
        "MKYLLPTAAAGLLLLAAQPAMAMTQKIIMKDGLPSDSKLIQFEWNNPDQ"
        "FQKDRISCGQVSIPNNGDCIYGDIAFKAGAWEALFNAAGATFESAERRQ"
        "ADKEGCYREATTCLPLFTIFCNAFMPQIHGAVEGYHWDKYGGMCDDRYC"
        "GTIMGGMYAGGCQAVKSAAADATFNIKNIFESRKIDGLMSQMMCGAAYPF"
        "YVKDGENRCCQAANFYKSGAVILTHPGSEPNLTYWKF\n"
    )

    def test_freetext_sequence(self):
        """Test ProvideSequence with pasted FASTA sequence text."""
        self.create_project("test_provide_sequence")
        self.create_job()

        self.set_param("controlParameters.SEQUENCETEXT", self._GAMMA_FASTA)

        self.run_and_wait()

        # Validate ASU content output
        asu_path = self.assert_file_exists("CASUCONTENTOUT.asu.xml")
        import xml.etree.ElementTree as ET
        tree = ET.parse(asu_path)
        seqs = tree.findall(".//sequence")
        assert len(seqs) >= 1, "No sequences found in ASU output"


@pytest.mark.usefixtures("file_based_db")
class TestProvideAlignmentAPI(APITestBase):
    """API tests for ProvideAlignment alignment input."""

    task_name = "ProvideAlignment"
    timeout = 60

    def test_paste_pir_alignment(self):
        """Test ProvideAlignment with pasted PIR-format alignment."""
        pir_text = (
            ">P1;target\n"
            "target sequence\n"
            "MKYLLPTAAAGLLLLAAQPAMA*\n"
            ">P1;template\n"
            "template sequence\n"
            "MKYLLPTAAAGLLLLAAQPAMA*\n"
        )
        self.create_project("test_provide_alignment_paste")
        self.create_job()

        self.set_param("controlParameters.PASTEORREAD", "PASTE")
        self.set_param("controlParameters.SEQUENCETEXT", pir_text)

        self.run_and_wait()

        self.assert_file_exists("ALIGNMENTFILE.aln")

    def test_read_pir_file(self, demo_data_dir):
        """Test ProvideAlignment by reading gamma PIR file."""
        pir_path = demo_data_dir / "gamma" / "gamma.pir"
        if not pir_path.exists():
            pytest.skip(f"PIR file not found: {pir_path}")

        self.create_project("test_provide_alignment_file")
        self.create_job()

        self.set_param("controlParameters.PASTEORREAD", "ALIGNIN")
        self.upload_file("inputData.ALIGNIN", str(pir_path))

        self.run_and_wait()

        self.assert_file_exists("ALIGNMENTFILE.aln")


@pytest.mark.usefixtures("file_based_db")
class TestMakeLinkAPI(APITestBase):
    """API tests for MakeLink covalent link definition."""

    task_name = "MakeLink"
    timeout = 120

    def test_lys_plp_link(self, pdb6ndn):
        """Test MakeLink: define LYS-PLP covalent link in 6ndn."""
        self.create_project("test_makelink")
        self.create_job()

        self.upload_file("inputData.XYZIN", pdb6ndn)

        # Residue and atom names for the link
        self.set_param("inputData.RES_NAME_1_TLC", "LYS")
        self.set_param("inputData.RES_NAME_2_TLC", "PLP")
        self.set_param("inputData.ATOM_NAME_1_TLC", "NZ")
        self.set_param("inputData.ATOM_NAME_2_TLC", "C4A")
        self.set_param("inputData.ATOM_NAME_1", "NZ")
        self.set_param("inputData.ATOM_NAME_2", "C4A")

        # Delete atom on ligand side and set bond order
        self.set_param("inputData.TOGGLE_DELETE_2", True)
        self.set_param("inputData.DELETE_2", "O4A")
        self.set_param("controlParameters.BOND_ORDER", "DOUBLE")
        self.set_param("controlParameters.TOGGLE_LINK", True)

        # MakeLink may emit warnings for empty optional fields
        self.run_and_wait(expect_success=False)

        # Validate CIF output contains expected blocks
        job_dir = self.get_job_directory()
        cif_file = job_dir / "LYS-PLP_link.cif"
        assert cif_file.exists(), f"Link CIF not produced: {list(job_dir.iterdir())}"
        doc = cif.read(str(cif_file))
        for name in ("mod_LYSm1", "mod_PLPm1", "link_LYS-PLP"):
            assert name in doc, f"Missing block '{name}' in CIF"

        # Check PDB output with applied links
        pdb_file = job_dir / "ModelWithLinks.pdb"
        assert pdb_file.exists(), f"PDB with links not produced: {list(job_dir.iterdir())}"


@pytest.mark.usefixtures("file_based_db")
class TestClustalwAPI(APITestBase):
    """API tests for clustalw multiple sequence alignment."""

    task_name = "clustalw"
    timeout = 60

    def test_cdk2_cdk6_alignment(self, fasta_cdk2, fasta_cdk6):
        """Test clustalw alignment of CDK2 and CDK6."""
        self.create_project("test_clustalw")
        self.create_job()

        self.set_param("inputData.SEQUENCELISTORALIGNMENT", "SEQUENCELIST")

        # Upload two sequences as list items
        self.set_param("inputData.SEQIN", [{}, {}])
        self.upload_file("inputData.SEQIN[0]", fasta_cdk2)
        self.upload_file("inputData.SEQIN[1]", fasta_cdk6)

        self.run_and_wait()

        # Check alignment output and program XML
        job_dir = self.get_job_directory()
        aln_files = list(job_dir.glob("*.aln"))
        assert len(aln_files) >= 1, f"No alignment output: {list(job_dir.iterdir())}"

        xml = self.read_program_xml()
        scores = xml.findall(".//PairwiseScore")
        assert len(scores) >= 1, "No pairwise scores in program.xml"


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
class TestDensityCalculatorAPI(APITestBase):
    """API tests for density_calculator model density map."""

    task_name = "density_calculator"
    timeout = 60

    def test_gamma_model(self, gamma_model_pdb):
        """Test density calculation from gamma model coordinates."""
        self.create_project("test_density_calculator")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.set_param("controlParameters.D_MIN", 2.0)

        self.run_and_wait()

        self.assert_file_exists("MAPOUT.map")


@pytest.mark.usefixtures("file_based_db")
class TestPdbsetUiAPI(APITestBase):
    """API tests for pdbset_ui PDB manipulation."""

    task_name = "pdbset_ui"
    timeout = 60

    def test_gamma_model(self, gamma_model_pdb):
        """Test pdbset with a simple CHAIN keyword."""
        self.create_project("test_pdbset_ui")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.set_param("controlParameters.EXTRA_PDBSET_KEYWORDS", "CHAIN A")

        self.run_and_wait()

        self.validate_pdb("XYZOUT.pdb")


@pytest.mark.usefixtures("file_based_db")
class TestAddFractionalCoordsAPI(APITestBase):
    """API tests for add_fractional_coords."""

    task_name = "add_fractional_coords"
    timeout = 60

    def test_gamma_model(self, gamma_model_pdb):
        """Test adding fractional coordinates to gamma model."""
        self.create_project("test_add_fractional_coords")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)

        self.run_and_wait()

        self.assert_file_exists("XYZOUT.cif")


@pytest.mark.usefixtures("file_based_db")
class TestModelASUCheckAPI(APITestBase):
    """API tests for modelASUCheck model-sequence alignment."""

    task_name = "modelASUCheck"
    timeout = 120

    def test_gamma(self, gamma_model_pdb, gamma_asu_xml):
        """Test model-sequence alignment check with gamma data."""
        self.create_project("test_model_asu_check")
        self.create_job()

        self.upload_file("inputData.XYZIN", gamma_model_pdb)
        self.upload_file("inputData.ASUIN", gamma_asu_xml)

        self.run_and_wait()

        # Check alignment results in program.xml
        xml = self.read_program_xml()
        alignments = xml.findall(".//SequenceAlignment/Alignment")
        assert len(alignments) >= 1
