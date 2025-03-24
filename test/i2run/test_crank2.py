from .utils import demoData, i2run
import gemmi


def test_auspex():
    seqFile = demoData("gamma", "gamma.pir")
    mtzFile = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["crank2"]
    args += [
        "--F_SIGFanom",
        f"fullPath={mtzFile}",
        "columnLabels=/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
    ]
    args += ["--SEQIN", f"seqFile={seqFile}"]
    args += ["--NUMBER_SUBSTRUCTURE", "2"]
    args += ["--ATOM_TYPE", "Xe"]
    args += ["--FPRIME", "-0.79"]
    args += ["--FDPRIME", "7.36"]
    args += ["--WAVELENGTH", "1.54179"]
    args += ["--END_PIPELINE", "dmfull"]
    with i2run(args) as directory:
        for name in ["FPHOUT_DIFFANOM", "FPHOUT_HL", "FPHOUT", "FREEROUT"]:
            gemmi.read_mtz_file(str(directory / f"{name}.mtz"))
