import re
import gemmi
from .programs import requires_program
from .utils import demoData, i2run


@requires_program("shelxc", "shelxd", "shelxe")
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
        # Each step's figures in the job's own program.xml (they were only in
        # the steps' sub-job parameters, whose numbering depends on the route)
        import xml.etree.ElementTree as ET
        steps = ET.parse(job / "program.xml").find("CrankSteps")
        assert steps is not None
        assert float(steps.findtext("Step[@name='substrdet']/CFOM")) > 30
