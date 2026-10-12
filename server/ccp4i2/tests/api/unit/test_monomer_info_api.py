"""API tests for GET /api/ccp4i2/monomer-info/<code>/.

The endpoint must answer on a server with no CCP4 (the slim server), and
must tell "no such monomer" apart from "the library could not be read":
the Make Link interface reports them differently.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.lib.utils.files import monomer_library as ml
from ccp4i2.lib.utils.files.monomer_library import MonomerLibraryError

# A minimal glycine restraint block: enough for gemmi to build a ChemComp.
GLY_CIF = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
GLY GLY GLYCINE peptide 5 5 .

data_comp_GLY
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
GLY N N NT3 1
GLY CA C CH2 0
GLY C C C 0
GLY O O O 0
GLY OXT O OC -1
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
GLY N CA single n
GLY CA C single n
GLY C O double n
GLY C OXT single n
"""


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


@pytest.fixture
def no_ccp4(monkeypatch):
    monkeypatch.delenv("CCP4", raising=False)
    monkeypatch.delenv("CLIBD_MON", raising=False)


def test_local_library(client, monkeypatch, tmp_path):
    (tmp_path / "g").mkdir()
    (tmp_path / "g" / "GLY.cif").write_text(GLY_CIF)
    monkeypatch.setenv("CLIBD_MON", str(tmp_path))

    resp = client.get("/api/ccp4i2/monomer-info/gly/")

    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["code"] == "GLY"
    assert data["atoms"] == ["N", "CA", "C", "O", "OXT"]
    assert {"atom1": "C", "atom2": "O", "type": "double"} in data["bonds"]


def test_no_ccp4_uses_the_published_library(client, no_ccp4, monkeypatch, tmp_path):
    """The slim-server case: no local library, the entry is fetched."""
    fetched = tmp_path / "GLY.cif"
    fetched.write_text(GLY_CIF)
    monkeypatch.setattr(ml, "find_monomer_cif", lambda code, opener=None: fetched)

    resp = client.get("/api/ccp4i2/monomer-info/GLY/")

    assert resp.status_code == 200
    assert resp.json()["data"]["atoms"] == ["N", "CA", "C", "O", "OXT"]


@pytest.mark.parametrize("status, reason", [
    (404, "not_found"),
    (502, "unavailable"),
])
def test_failures_carry_a_reason(client, no_ccp4, monkeypatch, status, reason):
    def fail(code, opener=None):
        raise MonomerLibraryError(status, "nope")

    monkeypatch.setattr(ml, "find_monomer_cif", fail)

    resp = client.get("/api/ccp4i2/monomer-info/CYS/")

    assert resp.status_code == status
    body = resp.json()
    assert body["success"] is False
    assert body["reason"] == reason
