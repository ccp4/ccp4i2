import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_substitute_ligand():
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--SMILESIN", '"CC(C)OC1=C(C=CC(=C1)OC)C2=NC(C(N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)C"']
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        doc = gemmi.cif.read(str(job / "DICTOUT.cif"))
        gemmi.make_chemcomp_from_block(doc[-1])
        for name in (
            "DIFFPHIOUT",
            "F_SIGF_OUT_asFMEAN",
            "F_SIGF_OUT",
            "FREERFLAG_OUT",
        ):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "selected_atoms.pdb"), format=gemmi.CoorFormat.Mmcif)
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23
        assert rfrees[-1] < 0.25
