import gemmi
from .utils import demoData, i2run


def test_parrot():
    args = ["parrot"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ABCD", demoData("gamma", "initial_phases.mtz")]
    args += ["--ASUIN", f"seqFile={demoData('gamma', 'gamma.pir')}"]
    with i2run(args) as directory:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(directory / f"{name}.mtz"))
