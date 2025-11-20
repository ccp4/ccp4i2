from gemmi import cif, read_pdb
from .urls import pdbe_mmcif, pdbe_pdb, rcsb_ligand_cif
from .utils import download, hasLongLigandName, i2run


def test_6ndn():
    with download(pdbe_pdb("6ndn")) as xyzPath:
        args = ["MakeLink"]
        args += ["--RES_NAME_1_TLC", "LYS"]
        args += ["--RES_NAME_2_TLC", "PLP"]
        args += ["--ATOM_NAME_1_TLC", "NZ"]
        args += ["--ATOM_NAME_2_TLC", "C4A"]
        args += ["--ATOM_NAME_1", "NZ"]
        args += ["--ATOM_NAME_2", "C4A"]
        args += ["--TOGGLE_DELETE_2", "True"]
        args += ["--DELETE_2", "O4A"]
        args += ["--BOND_ORDER", "DOUBLE"]
        args += ["--TOGGLE_LINK", "True"]
        args += ["--XYZIN", xyzPath]
        with i2run(args) as job:
            doc = cif.read(str(job / "LYS-PLP_link.cif"))
            for name in ("mod_LYSm1", "mod_PLPm1", "link_LYS-PLP"):
                assert name in doc
            read_pdb(str(job / "ModelWithLinks.pdb"))

def test_9efj():
    with download(rcsb_ligand_cif("A1BI3")) as dictPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "DICT"]
        args += ["--TLC", "A1BI3"]
        args += ["--DICTIN2", dictPath]
        with i2run(args) as acedrgJob:
            with download(pdbe_mmcif("9efj")) as xyzPath:
                args = ["MakeLink"]
                args += ["--MON_1_TYPE", "TLC"]
                args += ["--MON_2_TYPE", "CIF"]
                args += ["--DICT_2", str(acedrgJob / "A1BI3.cif")]
                args += ["--RES_NAME_1_TLC", "HIS"]
                args += ["--RES_NAME_2_CIF", "A1BI3"]
                args += ["--ATOM_NAME_1_TLC", "NE2"]
                args += ["--ATOM_NAME_2_CIF", "S10"]
                args += ["--ATOM_NAME_1", "NE2"]
                args += ["--ATOM_NAME_2", "S10"]
                args += ["--TOGGLE_DELETE_2", "True"]
                args += ["--DELETE_2", "O1"]
                args += ["--BOND_ORDER", "SINGLE"]
                args += ["--TOGGLE_LINK", "True"]
                args += ["--XYZIN", xyzPath]
                with i2run(args) as job:
                    doc = cif.read(str(job / "HIS-A1BI3_link.cif"))
                    for name in ("mod_HISm1", "mod_A1BI3m1", "link_HIS-A1BI3"):
                        assert name in doc
                    assert hasLongLigandName(job / "ModelWithLinks.cif")
