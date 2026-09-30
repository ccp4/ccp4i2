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
            # AceDRG also regularises the two monomers joined. That pair is
            # what the subjob's Moorhen view shows, so it has to be a real
            # output and not a path to a file AceDRG never wrote.
            dimers = sorted((job / "job_1").glob("*_for_link.pdb"))
            assert dimers, "AceDRG wrote no linked pair"
            dimer = read_pdb(str(dimers[0]))
            assert dimer[0].count_atom_sites() > 0
            assert sorted((job / "job_1").glob("*_for_link.cif")), \
                "linked-pair dictionary not brought into the job directory"


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


def test_residue_codes_typed_in_mixed_case():
    """"Lys" and "plp" link as LYS and PLP, as the atom dropdown implied.

    AceDRG found LYS.cif for "Lys", then looked inside it for a comp "Lys"
    and stopped, failing every MakeLink job typed that way.
    """
    with download(pdbe_pdb("6ndn")) as pdb:
        args = ["MakeLink"]
        args += ["--RES_NAME_1_TLC", "Lys"]
        args += ["--RES_NAME_2_TLC", "plp"]
        args += ["--ATOM_NAME_1_TLC", "NZ"]
        args += ["--ATOM_NAME_2_TLC", "C4A"]
        args += ["--ATOM_NAME_1", "NZ"]
        args += ["--ATOM_NAME_2", "C4A"]
        args += ["--TOGGLE_DELETE_2", "True"]
        args += ["--DELETE_2", "O4A"]
        args += ["--BOND_ORDER", "DOUBLE"]
        args += ["--TOGGLE_LINK", "True"]
        args += ["--XYZIN", pdb]
        with i2run(args, allow_errors=True) as job:
            doc = cif.read(str(job / "LYS-PLP_link.cif"))
            assert "link_LYS-PLP" in doc
            structure = read_pdb(str(job / "ModelWithLinks.pdb"))
            links = [c for c in structure.connections if c.link_id == "LYS-PLP"]
            assert links, "no LYS-PLP connection in the output model"


def test_several_edits_declared_as_lists():
    """The edit lists drive AceDRG: a deletion and a bond order on GLU, a
    charge on LYS -- three kinds of edit, two monomers, one instruction."""
    args = ["MakeLink"]
    args += ["--RES_NAME_1_TLC", "LYS"]
    args += ["--RES_NAME_2_TLC", "GLU"]
    args += ["--ATOM_NAME_1_TLC", "NZ"]
    args += ["--ATOM_NAME_2_TLC", "CD"]
    args += ["--ATOM_NAME_1", "NZ"]
    args += ["--ATOM_NAME_2", "CD"]
    args += ["--DELETE_ATOMS_2", "OE2"]
    args += ["--BOND_ORDERS_2", "ATOM_1=CD", "ATOM_2=OE1", "ORDER=DOUBLE"]
    args += ["--CHARGES_1", "ATOM=NZ", "CHARGE=0"]
    args += ["--TOGGLE_LINK", "False"]
    with i2run(args, allow_errors=True) as job:
        instruction = (job / "link_instruction.txt").read_text()
        assert "DELETE ATOM OE2 2" in instruction
        assert "CHANGE CHARGE 1 NZ 0" in instruction
        assert "CHANGE BOND CD OE1 DOUBLE 2" in instruction
        doc = cif.read(str(job / "LYS-GLU_link.cif"))
        assert "link_LYS-GLU" in doc
        # OE2 is gone from the modified GLU.
        mod = doc["mod_GLUm1"]
        deleted = [row for row in mod.find("_chem_mod_atom.", ["function", "atom_id"])
                   if row[0] == "delete"]
        assert [row[1] for row in deleted] == ["OE2"]
