"""The chain sequences a model digest reports (digest._extract_chain_sequences).

1H1S has a phosphorylated Thr160 (TPO). gemmi.one_letter_code on a list of
residue names writes X for every non-standard residue, so the AU-contents
table filled from the model held an X and failed validation as an invalid
amino acid (2026-10-31). A modified residue is reported as its parent.
"""
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi (ccp4-python)")

from ccp4i2.lib.utils.files.digest import _extract_chain_sequences, residue_one_letter


@pytest.mark.parametrize("name, code", [
    ("THR", "T"), ("TPO", "T"), ("SEP", "S"), ("MSE", "M"), ("PTR", "Y"),
    ("DA", "A"), ("DC", "C"), ("U", "U"),
    ("UNK", "X"), ("HOH", "X"), ("ATP", "X"), ("NOTARESIDUE", "X"),
])
def test_residue_one_letter(name, code):
    assert residue_one_letter(name) == code


def _structure(residue_names):
    pdb = []
    for i, name in enumerate(residue_names, start=1):
        pdb.append(f"ATOM  {i:5d}  CA  {name:>3s} A{i:4d}    {i * 3.8:8.3f}   0.000   0.000  1.00 20.00           C")
    pdb.append("END")
    structure = gemmi.read_pdb_string("\n".join(pdb))
    structure.setup_entities()
    return structure


def test_modified_residue_keeps_its_parent_letter():
    sequences = _extract_chain_sequences(_structure(["MET", "SER", "TPO", "VAL", "MSE", "ALA"]))
    assert sequences == {"A": "MSTVMA"}
