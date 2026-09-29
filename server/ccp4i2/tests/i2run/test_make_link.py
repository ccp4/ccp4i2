from gemmi import cif, read_pdb
from .urls import pdbe_pdb, rcsb_mmcif
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
        args += ["--XYZIN", pdb]
        # allow_errors=True: Allow diagnostic warnings for empty optional strings
        # (CHARGE_1_LIST, DELETE_1_LIST, etc.) - the script handles empty values correctly
        with i2run(args, allow_errors=True) as job:
            doc = cif.read(str(job / "LYS-PLP_link.cif"))
            for name in ("mod_LYSm1", "mod_PLPm1", "link_LYS-PLP"):
                assert name in doc
            # The point of the task: the model comes back carrying the link.
            # Reading it is not enough -- an empty model reads fine too.
            structure = read_pdb(str(job / "ModelWithLinks.pdb"))
            assert structure[0].count_atom_sites() > 0
            links = [c for c in structure.connections if c.link_id == "LYS-PLP"]
            assert links, "no LYS-PLP connection in the output model"


def test_6ndn_mmcif_gets_a_struct_conn():
    """An mmCIF model comes back with a _struct_conn row, not just a LINKR."""
    with download(rcsb_mmcif("6ndn")) as model:
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
        args += ["--XYZIN", model]
        with i2run(args, allow_errors=True) as job:
            block = cif.read(str(job / "ModelWithLinks.cif")).sole_block()
            rows = list(block.find("_struct_conn.", ["conn_type_id", "ccp4_link_id"]))
            assert any(row[0] == "covale" and row[1] == "LYS-PLP" for row in rows), \
                "no covale _struct_conn row carrying ccp4_link_id LYS-PLP"
