import gemmi
import pytest
from .utils import demoData, i2run


# NOTE: All phaser tests run FIRST (order=1) to avoid RDKit pickle contamination
# RDKit (imported by acedrg tests) modifies pickle module's dispatch table,
# causing phaser's pickle.dump() to fail. Running phaser tests before acedrg
# ensures pickle module is clean when phaser needs it.

@pytest.mark.order("first")
def test_phaser_ep():
    args = ["phaser_EP"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_PARTIAL", demoData("gamma", "heavy_atoms.pdb")]
    args += ["--XYZIN_HA", demoData("gamma", "heavy_atoms.pdb")]
    args += ["--ASUFILE", demoData("gamma", "gamma.asu.xml")]
    args += ["--COMP_BY", "ASU"]
    args += ["--WAVELENGTH", "1.542"]
    args += ["--LLGC_CYCLES", "20"]
    args += ["--ELEMENTS", "Xe"]
    args += ["--RUNPARROT", "False"]
    with i2run(args) as job:
        for name in ["ABCDOUT_1", "ABCDOUT_2", "MAPOUT_1", "MAPOUT_2"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        for name in ["PHASER.1", "PHASER.1.hand"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        model = gemmi.read_structure(str(job / "PHASER.1.pdb"))[0]
        occs = [cra.atom.occ for cra in model.all()]
        assert sum(occ > 0.11 for occ in occs) == 3
