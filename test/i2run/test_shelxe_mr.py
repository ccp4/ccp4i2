import xml.etree.ElementTree as ET
from gemmi import read_mtz_file, read_pdb
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_gamma():
    args = ["shelxeMR"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--NTCYCLES", "2"]
    args += ["--NMCYCLES", "2"]
    with i2run(args) as job:
        read_pdb(str(job / "shelxrun.pdb"))
        read_mtz_file(str(job / "FPHIOUT.mtz"))
        xml = ET.parse(job / "program.xml")
        assert float(xml.find(".//BestCC").text) > 43
