"""The MakeLink report reads what the link changes from the dictionary.

The dictionary keeps each monomer as it was and lists the changes in a
_chem_mod block; the shape here is AceDRG's (ccp4-20260702) for a
cysteine added across the vinyl of N-methylacrylamide, cut down.
"""
import pytest

pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.pipelines.MakeLink.script.link_dictionary import describe_bond, link_modifications

LINK_CIF = """\
data_mod_list
loop_
_chem_mod.id
_chem_mod.name
_chem_mod.comp_id
_chem_mod.group_id
CYSm1 CYSTEINE CYS peptide
DRGm1 . DRG non-polymer

data_comp_DRG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
DRG C1 C
DRG C2 C
DRG C3 C
DRG H1 H
DRG H2 H
DRG H3 H
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
_chem_comp_bond.value_dist
DRG C1 C2 DOUBLE n 1.306
DRG C2 C3 SINGLE n 1.480
DRG C1 H1 SINGLE n 1.085
DRG C2 H3 SINGLE n 1.085

data_mod_CYSm1
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_charge
CYSm1 delete HG . H HSH1 0
CYSm1 change SG . S S2 0
loop_
_chem_mod_bond.mod_id
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist
CYSm1 delete SG HG single .
CYSm1 change CB SG single 1.811
loop_
_chem_mod_angle.mod_id
_chem_mod_angle.function
_chem_mod_angle.atom_id_1
_chem_mod_angle.atom_id_2
_chem_mod_angle.atom_id_3
CYSm1 delete CB SG HG
CYSm1 change CA CB SG

data_mod_DRGm1
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_charge
DRGm1 change C1 . C CH2 0
DRGm1 change C2 . C CH2 0
DRGm1 add . H32 H H 0
loop_
_chem_mod_bond.mod_id
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist
DRGm1 change C1 C2 single 1.520
DRGm1 change C1 H1 single 0.983
DRGm1 add C2 H32 single 0.981
loop_
_chem_mod_plane_atom.mod_id
_chem_mod_plane_atom.function
_chem_mod_plane_atom.plane_id
_chem_mod_plane_atom.atom_id
DRGm1 delete plan-1 C1
DRGm1 delete plan-1 C2
DRGm1 delete plan-2 C1

data_link_CYS-DRG
loop_
_chem_link_bond.link_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_id_2
CYS-DRG SG C1
"""


@pytest.fixture
def link_cif(tmp_path):
    path = tmp_path / "CYS-DRG_link.cif"
    path.write_text(LINK_CIF)
    return path


def test_the_vinyl_is_single_in_the_modification_though_double_in_the_monomer(link_cif):
    cys, drg = link_modifications(link_cif)
    assert (drg.comp_id, drg.mod_id) == ("DRG", "DRGm1")
    [vinyl] = drg.bonds   # the C1-H1 length change is a hydrogen's, not shown
    assert (vinyl.atom_1, vinyl.atom_2, vinyl.type, vinyl.dist) == ("C1", "C2", "single", "1.520")
    assert (vinyl.old_type, vinyl.old_dist) == ("DOUBLE", "1.306")
    assert describe_bond(vinyl) == "DOUBLE 1.306 A -> single 1.520 A"


def test_an_added_hydrogen_is_named_with_the_atom_it_hangs_off(link_cif):
    _cys, drg = link_modifications(link_cif)
    [added] = drg.atoms
    assert (added.function, added.name, added.element, added.bonded_to) == ("add", "H32", "H", "C2")
    assert drg.retyped == ["C1", "C2"]
    assert drg.revised == {"planes": 2}


def test_a_library_monomer_has_no_before_values(link_cif):
    cys, _drg = link_modifications(link_cif)
    assert cys.comp_id == "CYS"
    [gone] = cys.atoms
    assert (gone.function, gone.name, gone.bonded_to) == ("delete", "HG", "SG")
    [cb_sg] = cys.bonds
    assert describe_bond(cb_sg) == "single 1.811 A"
    assert cys.revised == {"angles": 2}


def test_no_file_means_no_modifications(tmp_path):
    assert link_modifications(tmp_path / "missing.cif") == []
