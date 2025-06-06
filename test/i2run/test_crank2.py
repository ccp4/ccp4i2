import re
import gemmi
from .utils import demoData, i2run


def test_crank2():
    args = ["crank2"]
    args += [
        "--F_SIGFanom",
        f"fullPath={demoData('gamma', 'merged_intensities_Xe.mtz')}",
        "columnLabels=/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
    ]
    args += ["--SEQIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--NUMBER_SUBSTRUCTURE", "2"]
    args += ["--ATOM_TYPE", "Xe"]
    args += ["--FPRIME", "-0.79"]
    args += ["--FDPRIME", "7.36"]
    args += ["--WAVELENGTH", "1.54179"]
    args += ["--END_PIPELINE", "dmfull"]
    with i2run(args) as job:
        for name in ["FPHOUT_DIFFANOM", "FPHOUT_HL", "FPHOUT", "FREEROUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        log = (job / "log.txt").read_text()
        foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
        assert max(foms) > 0.7
