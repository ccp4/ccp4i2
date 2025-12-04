import xml.etree.ElementTree as ET
import gemmi
import pytest
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)

# NOTE: All phaser tests run FIRST (order="first") to avoid RDKit pickle contamination
# RDKit (imported by acedrg tests) modifies pickle module's dispatch table,
# causing phaser's pickle.dump() to fail. Running phaser tests before acedrg
# ensures pickle module is clean when phaser needs it.

@pytest.mark.order("first")
def test_gamma_basic():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--RUNREFMAC", "False"]
    args += ["--RUNSHEETBEND", "False"]
    with i2run(args) as job:
        for name in ["PHASER.1"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        assert max(llgs) > 1000

@pytest.mark.order("first")
def test_gamma_sheetbend():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--RUNREFMAC", "False"]
    args += ["--RUNSHEETBEND", "True"]
    with i2run(args) as job:
        for name in ["PHASER.1", "XYZOUT_SHEETBEND"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        assert max(llgs) > 1000

@pytest.mark.order("first")
def test_gamma():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["PHASER.1", "XYZOUT_REFMAC", "XYZOUT_SHEETBEND"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        assert min(rworks) < 0.27
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        assert max(llgs) > 1000


@pytest.mark.order("first")
def test_no_solution():
    """Test that phaser correctly reports no solution when data doesn't match model."""
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("rnase", "rnase_model.pdb")]
    args += ["--RESOLUTION_HIGH", "5.0"]
    # Use allow_errors=True because this test expects phaser to fail (no solution found)
    with i2run(args, allow_errors=True) as job:
        assert not (job / "PHASER.1.pdb").exists()
