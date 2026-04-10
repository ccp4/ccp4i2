from sys import platform
import re
from gemmi import read_mtz_file, read_pdb
from pytest import mark
from .utils import demoData, i2run


pytestmark = mark.skipif(platform == "win32", reason="Not supported on Windows")


# TODO: Test long ligand names (e.g. 8xfm)


def test_arpwarp():
    args = ["arp_warp_classic"]
    args += ["--AWA_FOBS", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--AWA_FREE", demoData("gamma", "freeR.mtz")]
    args += ["--AWA_PHINI", demoData("gamma", "initial_phases.mtz")]
    args += ["--AWA_PHREF", demoData("gamma", "initial_phases.mtz")]
    args += ["--AWA_SEQIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--AWA_BIG_CYCLES", "1"]
    args += ["--AWA_SMALL_CYCLES", "1"]
    with i2run(args) as job:
        for name in ["XYZDUM", "XYZOUT"]:
            read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFFPHIOUT", "FPHIOUT"]:
            read_mtz_file(str(job / f"{name}.mtz"))
        log = (job / "log.txt").read_text()
        pattern = re.compile(r"R = [\.\d]+ \(Rfree = ([\.\d]+)\)")
        rfrees = [float(x) for x in pattern.findall(log)]
        assert rfrees[-1] < 0.4
