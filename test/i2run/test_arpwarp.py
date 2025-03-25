import gemmi
from .utils import demoData, i2run


def test_arpwarp():
    args = ["arp_warp_classic"]
    args += ["--AWA_FOBS", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--AWA_FREE", demoData("gamma", "freeR.mtz")]
    args += ["--AWA_PHINI", demoData("gamma", "initial_phases.mtz")]
    args += ["--AWA_PHREF", demoData("gamma", "initial_phases.mtz")]
    args += ["--AWA_SEQIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--AWA_BIG_CYCLES", "1"]
    args += ["--AWA_SMALL_CYCLES", "1"]
    with i2run(args) as directory:
        for name in ["XYZDUM", "XYZOUT"]:
            gemmi.read_pdb(str(directory / f"{name}.pdb"))
        for name in ["DIFFPHIOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(directory / f"{name}.mtz"))
