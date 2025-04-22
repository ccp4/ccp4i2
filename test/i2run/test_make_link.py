from gemmi import cif, read_pdb
from .urls import pdbe_pdb
from .utils import download, i2run


def test_6ndn():
    with download(pdbe_pdb("6ndn")) as pdb:
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
        args += ["--XYZIN", str(pdb)]
        with i2run(args) as job:
            doc = cif.read(str(job / "LYS-PLP_link.cif"))
            for name in ("mod_LYSm1", "mod_PLPm1", "link_LYS-PLP"):
                assert name in doc
            read_pdb(str(job / "ModelWithLinks.pdb"))
