from .utils import demoData, i2run


def test_substitute_ligand():
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--SMILESIN", "CC(C)OC1=C(C=CC(=C1)OC)C2=NC(C(N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)C"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        assert False, "Check output files"
