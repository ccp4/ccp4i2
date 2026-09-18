"""The Coot harvest must tell a restraint dictionary CIF from a
coordinate model CIF, so a ligand-builder dictionary is filed as DICTOUT
(coot1) / caught regardless of filename (coot_rebuild), not mis-harvested
as a coordinate model. Shared classifier in cootbridge.harvest."""

import pytest

pytest.importorskip("gemmi", reason="CIF classifier uses gemmi")


def _classifier():
    from ccp4i2.cootbridge.harvest import cif_is_restraint_dictionary
    return cif_is_restraint_dictionary


DICT_CIF = """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
ANI ANI
data_comp_ANI
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
ANI N1 N
ANI C1 C
"""

MODEL_CIF = """data_model
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N 1.0 2.0 3.0
ATOM 2 C 1.5 2.5 3.5
"""


def test_restraint_dictionary_detected(tmp_path):
    path = tmp_path / "aniline.cif"
    path.write_text(DICT_CIF)
    assert _classifier()(path) is True


def test_coordinate_model_not_flagged_as_dictionary(tmp_path):
    path = tmp_path / "model.cif"
    path.write_text(MODEL_CIF)
    assert _classifier()(path) is False


def test_unreadable_cif_defaults_to_model(tmp_path):
    """On any read failure, treat as a model (prior behaviour) rather
    than silently diverting a file to DICTOUT."""
    path = tmp_path / "broken.cif"
    path.write_text("this is not valid cif {{{")
    assert _classifier()(path) is False
