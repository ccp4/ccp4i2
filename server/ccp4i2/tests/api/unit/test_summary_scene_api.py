"""
E2E test for the fragment-campaign summary-scene endpoint.

    GET /api/ccp4i2/projectgroups/{id}/summary_scene/

This drives the whole path the real feature uses: DB rows (ProjectGroup ->
memberships -> Projects -> Jobs -> Files) plus the actual coordinate /
dictionary files on disk, read back by gemmi in
``lib/campaign_scene.build_summary_scene``.

The fixture builds a synthetic campaign:

  * a parent reference (apo, no ligand)         -> parent ribbon
  * member "frag-drg": refined, ligand DRG      -> hit
  * member "frag-lig": refined, ligand LIG      -> hit
  * member "frag-apo": refined, no ligand       -> skipped

so we can assert the parent ribbon, two fragment-stick elements (each scoped
to its own dictionary), and the skip accounting. The shapes are synthetic;
expect to revisit thresholds/labels once this runs against real XChem data,
but the contract (parent + per-hit element + scoped dict + stats) should hold.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets — do NOT add @pytest.mark.django_db here.
"""

import uuid
from pathlib import Path

import gemmi
from rest_framework.test import APIClient

from ccp4i2.db import models


# --------------------------------------------------------------------------
# Minimal on-disk coordinate / dictionary builders (gemmi only)
# --------------------------------------------------------------------------

def _write_pdb(path: Path, ligand_code=None):
    """One-residue alanine 'protein', optionally plus a HET ligand."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")

    res = gemmi.Residue()
    res.name = "ALA"
    res.seqid = gemmi.SeqId("1")
    for nm, el in [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]:
        a = gemmi.Atom()
        a.name = nm
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(1, 2, 3)
        a.occ = 1.0
        a.b_iso = 20.0
        res.add_atom(a)
    chain.add_residue(res)

    if ligand_code:
        lig = gemmi.Residue()
        lig.name = ligand_code
        lig.seqid = gemmi.SeqId("2")
        lig.het_flag = "H"
        a = gemmi.Atom()
        a.name = "C1"
        a.element = gemmi.Element("C")
        a.pos = gemmi.Position(5, 5, 5)
        a.occ = 1.0
        a.b_iso = 30.0
        lig.add_atom(a)
        chain.add_residue(lig)

    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    path.write_text(st.make_pdb_string())


def _write_dict(path: Path, *codes):
    """Minimal refmac-style restraint CIF declaring one or more monomers."""
    list_rows = "\n".join(f"{c} {c} 'monomer' non-polymer 3 3 ." for c in codes)
    comp_blocks = "\n".join(
        f"""data_comp_{c}
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
{c} C1 C
{c} C2 C
{c} O1 O"""
        for c in codes
    )
    path.write_text(
        f"""data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
{list_rows}
{comp_blocks}
"""
    )


# --------------------------------------------------------------------------
# DB + on-disk fixture builders
# --------------------------------------------------------------------------

def _make_refined_project(root: Path, name, pdb_type, dict_type,
                          ligand_code=None, dict_codes=None):
    """A Project with one finished refmac job carrying XYZOUT (+ optional dict)."""
    project = models.Project.objects.create(
        name=name, directory=str(root / name)
    )
    job = models.Job.objects.create(
        uuid=uuid.uuid4(),
        project=project,
        number="1",
        title=f"{name} refinement",
        task_name="refmac",
        status=models.Job.Status.FINISHED,
    )
    job.directory.mkdir(parents=True, exist_ok=True)

    _write_pdb(job.directory / "XYZOUT.pdb", ligand_code=ligand_code)
    models.File.objects.create(
        uuid=uuid.uuid4(),
        name="XYZOUT.pdb",
        directory=models.File.Directory.JOB_DIR,
        type=pdb_type,
        job=job,
        job_param_name="XYZOUT",
    )

    if dict_codes:
        _write_dict(job.directory / "LIBOUT.cif", *dict_codes)
        models.File.objects.create(
            uuid=uuid.uuid4(),
            name="LIBOUT.cif",
            directory=models.File.Directory.JOB_DIR,
            type=dict_type,
            job=job,
            job_param_name="LIBOUT",
        )
    return project


def _build_campaign(root: Path):
    """Parent + 2 hit members + 1 apo member, wired into a fragment_set group."""
    pdb_type, _ = models.FileType.objects.get_or_create(name="chemical/x-pdb")
    dict_type, _ = models.FileType.objects.get_or_create(
        name="application/refmac-dictionary"
    )

    parent = _make_refined_project(root, "parent_ref", pdb_type, dict_type)
    frag_drg = _make_refined_project(
        root, "frag_drg", pdb_type, dict_type, ligand_code="DRG", dict_codes=["DRG"]
    )
    frag_lig = _make_refined_project(
        root, "frag_lig", pdb_type, dict_type, ligand_code="LIG", dict_codes=["LIG"]
    )
    # Apo: refined but no ligand bound; dict present but irrelevant.
    frag_apo = _make_refined_project(
        root, "frag_apo", pdb_type, dict_type, ligand_code=None, dict_codes=["DRG"]
    )

    group = models.ProjectGroup.objects.create(
        name="test_fragment_campaign",
        type=models.ProjectGroup.GroupType.FRAGMENT_SET,
    )
    models.ProjectGroupMembership.objects.create(
        group=group, project=parent,
        type=models.ProjectGroupMembership.MembershipType.PARENT,
    )
    for member in (frag_drg, frag_lig, frag_apo):
        models.ProjectGroupMembership.objects.create(
            group=group, project=member,
            type=models.ProjectGroupMembership.MembershipType.MEMBER,
        )
    return group


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

def test_summary_scene_endpoint(bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    data = response.json()

    # ---- stats ----------------------------------------------------------
    stats = data["stats"]
    assert stats["members_total"] == 3
    assert stats["hits"] == 2               # DRG + LIG, not the apo member
    assert stats["parent_present"] is True
    skipped_reasons = [s["reason"] for s in stats["skipped"]]
    assert len(skipped_reasons) == 1
    assert "no ligand" in skipped_reasons[0].lower()

    # ---- scene shape ----------------------------------------------------
    scene = data["scene"]
    assert scene["version"] == 1
    files = scene["files"]
    elements = scene["elements"]

    # Parent reference: a coordinates file rendered as a ribbon (CRs).
    ref_files = [f for f in files if f["name"] == "reference"]
    assert len(ref_files) == 1
    assert ref_files[0]["fileId"] is not None
    assert ref_files[0]["projectId"]
    ribbon_elements = [
        e for e in elements
        if any(r["style"] == "CRs" for r in e["representations"])
    ]
    assert len(ribbon_elements) == 1
    assert ribbon_elements[0]["file"] == "reference"

    # Two fragment elements: ligand sticks (CBs), each scoped to a dictionary.
    stick_elements = [
        e for e in elements
        if any(r["style"] == "CBs" for r in e["representations"])
    ]
    assert len(stick_elements) == 2

    selections = " ".join(
        r["selection"]
        for e in stick_elements
        for r in e["representations"]
    )
    assert "DRG" in selections
    assert "LIG" in selections

    # Challenge B: every hit element references its own dictionary, and that
    # dictionary is declared in files[] as kind: dictionary.
    dict_names = {f["name"] for f in files if f.get("kind") == "dictionary"}
    assert len(dict_names) == 2
    for element in stick_elements:
        assert element["dictionaries"], f"{element['file']} has no scoped dict"
        for dname in element["dictionaries"]:
            assert dname in dict_names


def test_summary_scene_empty_campaign(bypass_api_permissions, test_project_path):
    """A campaign with no members yields an empty-but-valid scene, not a 500."""
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = models.ProjectGroup.objects.create(
        name="empty_campaign",
        type=models.ProjectGroup.GroupType.FRAGMENT_SET,
    )

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    data = response.json()
    assert data["stats"]["members_total"] == 0
    assert data["stats"]["hits"] == 0
    assert data["stats"]["parent_present"] is False
    assert data["scene"]["files"] == []
    assert data["scene"]["elements"] == []
