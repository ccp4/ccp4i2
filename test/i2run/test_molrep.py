import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_molrep():
    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["XYZOUT_MOLREP", "XYZOUT_SHEETBEND", "XYZOUT"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.findall(".//Cycle/r_factor")]
        rfrees = [float(e.text) for e in xml.findall(".//Cycle/r_free")]
        assert rworks[-1] < rworks[0]
        assert rfrees[-1] < rfrees[0]
        assert rworks[-1] < 0.26
        assert rfrees[-1] < 0.28
