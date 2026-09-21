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

import math
import uuid
from pathlib import Path

import gemmi
import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models


# --------------------------------------------------------------------------
# Minimal on-disk coordinate / dictionary builders (gemmi only)
# --------------------------------------------------------------------------

N_RESIDUES = 20
LIG_SHIFT = (1.0, -0.5, 2.0)


def _write_pdb(path: Path, ligand_code=None, shift=(0.0, 0.0, 0.0)):
    """A short poly-alanine helix, optionally plus a HET ligand.

    Twenty residues with distinct CA positions, so the summary scene's CA
    fit onto the reference is well conditioned (it needs at least twelve).
    ``shift`` moves the whole model, standing in for a member deposited in
    a different frame.
    """
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(60, 60, 60, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")

    sx, sy, sz = shift
    for i in range(N_RESIDUES):
        res = gemmi.Residue()
        res.name = "ALA"
        res.seqid = gemmi.SeqId(str(i + 1))
        theta = math.radians(100.0 * i)
        cx = 10.0 + 2.3 * math.cos(theta) + sx
        cy = 10.0 + 2.3 * math.sin(theta) + sy
        cz = 5.0 + 1.5 * i + sz
        for nm, el, (dx, dy, dz) in [
            ("N", "N", (-1.2, 0.3, -0.5)),
            ("CA", "C", (0.0, 0.0, 0.0)),
            ("C", "C", (1.0, 0.8, 0.4)),
            ("O", "O", (1.4, 1.9, 0.2)),
        ]:
            a = gemmi.Atom()
            a.name = nm
            a.element = gemmi.Element(el)
            a.pos = gemmi.Position(cx + dx, cy + dy, cz + dz)
            a.occ = 1.0
            a.b_iso = 20.0
            res.add_atom(a)
        chain.add_residue(res)

    if ligand_code:
        lig = gemmi.Residue()
        lig.name = ligand_code
        lig.seqid = gemmi.SeqId("101")
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
                          ligand_code=None, dict_codes=None,
                          shift=(0.0, 0.0, 0.0)):
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

    _write_pdb(job.directory / "XYZOUT.pdb", ligand_code=ligand_code, shift=shift)
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
    # Deposited in its own frame: what the summary's fit exists to undo.
    frag_lig = _make_refined_project(
        root, "frag_lig", pdb_type, dict_type, ligand_code="LIG", dict_codes=["LIG"],
        shift=LIG_SHIFT,
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

    # ---- superposition ---------------------------------------------------
    # Every hit is fitted onto the reference on all shared CAs and carries
    # the transform, with the evidence it rests on. frag_lig was written in
    # a shifted frame, so its fit must undo exactly that shift; frag_drg is
    # already in frame, so its fit is the identity.
    superpose = {sp["move"]: sp for sp in scene["superpose"]}
    assert set(superpose) == {e["file"] for e in stick_elements}
    for entry in superpose.values():
        assert entry["method"] == "matrix"
        assert len(entry["mat"]) == 9 and len(entry["vec"]) == 3
        assert entry["fitted"]["onto"] == "reference"
        assert entry["fitted"]["atoms"] == N_RESIDUES
        assert entry["fitted"]["rmsd"] == pytest.approx(0.0, abs=2e-3)
        assert "radius" not in entry["fitted"]   # a global fit
        assert entry["mat"] == pytest.approx([1, 0, 0, 0, 1, 0, 0, 0, 1], abs=1e-3)
    assert superpose["frag_drg"]["vec"] == pytest.approx([0, 0, 0], abs=2e-3)
    assert superpose["frag_lig"]["vec"] == pytest.approx(
        [-c for c in LIG_SHIFT], abs=2e-3
    )
    fits = {s["project"]: s for s in stats["superpose"]}
    assert set(fits) == {"frag_drg", "frag_lig"}
    assert all(f["ok"] and f["atoms"] == N_RESIDUES and f["reason"] is None
               for f in fits.values())


def test_summary_scene_promotes_a_hit_when_the_parent_has_no_model(
    bypass_api_permissions, test_project_path
):
    """With no parent coordinates, a hit carries the reference ribbon.

    A campaign whose parent project was never populated used to render as
    fragments floating in space with no reference and no explanation. The
    ribbon now comes from a member, drawn from the same files[] entry that
    already carries its ligand -- so no second copy is downloaded -- and
    stats say whose frame the scene is in.
    """
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)
    # Strip the parent, leaving the members untouched: the state a campaign
    # is in before its reference model has been imported.
    group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.PARENT
    ).delete()

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    data = response.json()

    assert data["stats"]["parent_present"] is False
    assert data["stats"]["reference"]["is_parent"] is False
    exemplar = data["stats"]["reference"]["project"]
    assert exemplar in {"frag_drg", "frag_lig"}

    scene = data["scene"]
    # Exactly one ribbon, and it is a member's.
    ribbons = [
        e for e in scene["elements"]
        if any(r["style"] == "CRs" for r in e["representations"])
    ]
    assert len(ribbons) == 1
    assert ribbons[0]["file"] != "reference"

    # It is the SAME element as that member's ligand: one file, two
    # representations, not a duplicate download.
    styles = [r["style"] for r in ribbons[0]["representations"]]
    assert styles == ["CRs", "CBs"], styles
    assert len([f for f in scene["files"] if f["name"] == ribbons[0]["file"]]) == 1

    # Both hits still draw their ligands.
    assert data["stats"]["hits"] == 2

    # The promoted hit IS the frame, so it is not fitted; the other hit is
    # fitted onto it, and the provenance names the member, not "reference".
    superpose = scene["superpose"]
    assert len(superpose) == 1
    assert superpose[0]["move"] != ribbons[0]["file"]
    assert superpose[0]["fitted"]["onto"] == ribbons[0]["file"]
    assert superpose[0]["fitted"]["atoms"] == N_RESIDUES
    assert [s["project"] for s in data["stats"]["superpose"]] != [exemplar]
    assert len(data["stats"]["superpose"]) == 1


def test_summary_scene_parent_reference_wins_over_a_hit(
    bypass_api_permissions, test_project_path
):
    """With a parent present, the ribbon is the parent's and no hit gains one."""
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    data = response.json()

    assert data["stats"]["reference"] == {
        "project": "parent_ref", "is_parent": True
    }
    ribbons = [
        e for e in data["scene"]["elements"]
        if any(r["style"] == "CRs" for r in e["representations"])
    ]
    assert len(ribbons) == 1
    assert ribbons[0]["file"] == "reference"


def test_summary_scene_carries_the_first_site_camera(
    bypass_api_permissions, test_project_path
):
    """The scene's view comes from the campaign's first binding site.

    A regression test with a history: sites moved out of a ProjectGroup JSON
    field into the CampaignSite table in migration 0024, and the builder went
    on reading the old attribute through a defaulted getattr. Nothing raised;
    scenes simply came back with no camera, so the viewer opened on whatever
    it happened to be looking at. Assert the camera is there, and that it is
    the first site's by display order rather than by insertion.
    """
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)

    # Added second but ordered first: pins ordering, not insertion order.
    models.CampaignSite.objects.create(
        group=group, name="second pocket",
        origin_x=1.0, origin_y=2.0, origin_z=3.0, order=1,
    )
    models.CampaignSite.objects.create(
        group=group, name="main pocket",
        origin_x=12.5, origin_y=-3.25, origin_z=28.75,
        quat=[0.0, 0.0, 0.0, 1.0], zoom=0.35, order=0,
    )

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    view = response.json()["scene"]["view"]

    assert view["origin"] == [12.5, -3.25, 28.75]
    assert view["quat"] == [0.0, 0.0, 0.0, 1.0]
    assert view["zoom"] == 0.35


def test_summary_scene_without_sites_has_no_view(
    bypass_api_permissions, test_project_path
):
    """No sites means no camera key at all -- not a camera at the origin.

    The distinction matters: an absent `view` leaves the viewer's own framing
    alone, where `origin: [0, 0, 0]` would actively point it at the corner of
    the cell.
    """
    test_project_path.mkdir(parents=True, exist_ok=True)
    group = _build_campaign(test_project_path)

    client = APIClient()
    response = client.get(
        f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/"
    )
    assert response.status_code == 200, response.content
    assert "view" not in response.json()["scene"]


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
