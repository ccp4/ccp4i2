import gemmi
from .utils import demoData, i2run


def test_mrbump():
    args = ["mrbump_basic"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--ASUIN", f"seqFile={demoData('gamma', 'gamma.pir')}"]
    args += ["--MRMAX", "5"]
    args += ["--REDUNDANCYLEVEL", "95"]
    args += ["--PJOBS", "2"]
    args += ["--NCYC", "10"]
    with i2run(args) as directory:
        gemmi.read_pdb(str(directory / "output_mrbump_1.pdb"))
        gemmi.read_mtz_file(str(directory / "output_mrbump_1.mtz"))
