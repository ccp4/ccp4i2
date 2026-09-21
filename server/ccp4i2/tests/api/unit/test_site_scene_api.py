"""
E2E tests for the per-site campaign scene endpoint.

    GET /api/ccp4i2/projectgroups/{id}/sites/{site_id}/scene/

Same synthetic campaign as ``test_summary_scene_api`` (parent + two hit
members + one apo member), with ``SiteEvaluation`` rows added per test,
because membership of a site scene is the verdict and nothing else: what
``detect_ligands`` finds decides what a hit *draws*, never whether it is in.

The site origin sits on the members' ligand (5, 5, 5). Around it the helix
supplies 9 CAs within 15 A and 13 within 20 A, so the local fit has to grow
once to reach its floor of 12 -- which pins the radius growth on the real
path, not just in the pure-function tests.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.lib import campaign_scene

from .test_summary_scene_api import LIG_SHIFT, _build_campaign

SITE_ORIGIN = (5.0, 5.0, 5.0)
Verdict = models.SiteEvaluation.Verdict


def _site(group, name="pocket", origin=SITE_ORIGIN, **extra):
    x, y, z = origin
    return models.CampaignSite.objects.create(
        group=group, name=name, origin_x=x, origin_y=y, origin_z=z, **extra
    )


def _project(name):
    return models.Project.objects.get(name=name)


def _evaluate(site, project_name, verdict):
    models.SiteEvaluation.objects.create(
        project=_project(project_name), site=site, verdict=verdict
    )


def _get(group, site, query=""):
    return APIClient().get(
        f"/api/ccp4i2/projectgroups/{group.id}/sites/{site.id}/scene/{query}"
    )


def _sticks(scene):
    """Elements that draw ligand sticks (not the reference's pocket)."""
    return [
        e for e in scene["elements"]
        if e["file"] != "reference"
        and any(r["style"] == "CBs" for r in e["representations"])
    ]


def _ribbons(scene):
    return [
        e for e in scene["elements"]
        if any(r["style"] == "CRs" for r in e["representations"])
    ]


@pytest.fixture
def campaign(bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    return _build_campaign(test_project_path)


def test_two_hits_give_two_ligand_elements_and_one_ribbon(campaign):
    site = _site(campaign, quat=[0.0, 0.0, 0.0, 1.0], zoom=0.35)
    _evaluate(site, "frag_drg", Verdict.HIT)
    _evaluate(site, "frag_lig", Verdict.HIT)
    _evaluate(site, "frag_apo", Verdict.EMPTY)

    response = _get(campaign, site)
    assert response.status_code == 200, response.content
    data = response.json()
    scene, stats = data["scene"], data["stats"]

    assert stats["site"] == {"id": site.id, "name": "pocket"}
    assert stats["hits_claimed"] == 2
    assert stats["hits_drawn"] == 2
    assert stats["empty_verdicts"] == 1
    assert stats["unclear_verdicts"] == 0
    assert stats["skipped"] == []
    assert stats["reference"] == {"project": "parent_ref", "is_parent": True}

    # One ribbon, the parent's; two stick elements, each with its own
    # dictionary and its own colour.
    ribbons = _ribbons(scene)
    assert [e["file"] for e in ribbons] == ["reference"]
    sticks = _sticks(scene)
    assert {e["file"] for e in sticks} == {"frag_drg", "frag_lig"}
    colours = [e["representations"][0]["colour"] for e in sticks]
    assert len(set(colours)) == 2
    assert all(c in campaign_scene.HIT_COLOURS for c in colours)
    assert all(e["dictionaries"] for e in sticks)

    # The pocket: the reference's residues around the origin as sticks, as
    # an explicit residue list. Residue 3's CA is 6 A from the origin;
    # residue 20 is at the far end of the helix.
    ref_reps = ribbons[0]["representations"]
    assert [r["style"] for r in ref_reps] == ["CRs", "CBs"]
    pocket = ref_reps[1]["selection"]
    assert "//A/3" in pocket.split("||")
    assert "//A/20" not in pocket.split("||")
    assert ref_reps[1]["colour"] == campaign_scene.POCKET_STICK_COLOUR
    assert stats["pocket_residues"] == len(pocket.split("||"))

    # The camera is the site's, and the slab hangs off the pocket.
    view = scene["view"]
    assert view["origin"] == list(SITE_ORIGIN)
    assert view["quat"] == [0.0, 0.0, 0.0, 1.0]
    assert view["zoom"] == 0.35
    assert view["slab"] == {
        "file": "reference",
        "selection": pocket,
        "pad": campaign_scene.ENVIRONMENT_RADIUS,
    }

    # Each hit is fitted locally: the sphere had to grow once (9 CAs at
    # 15 A, 13 at 20 A), and frag_lig's shift is undone exactly.
    superpose = {s["move"]: s for s in scene["superpose"]}
    assert set(superpose) == {"frag_drg", "frag_lig"}
    for entry in superpose.values():
        assert entry["method"] == "matrix"
        assert entry["fitted"]["onto"] == "reference"
        assert entry["fitted"]["radius"] == 20.0
        assert entry["fitted"]["atoms"] == 13
        assert entry["fitted"]["rmsd"] == pytest.approx(0.0, abs=2e-3)
    assert superpose["frag_lig"]["vec"] == pytest.approx(
        [-c for c in LIG_SHIFT], abs=2e-3
    )
    fits = {s["project"]: s for s in stats["superpose"]}
    assert all(f["ok"] and f["radius"] == 20.0 and "fallback" not in f
               for f in fits.values())

    # The diagnostic is measured in the frame the scene draws. The fixture
    # shifts frag_lig's protein but leaves its ligand at (5, 5, 5), so the
    # fit that brings the protein home carries the ligand away from the
    # site by exactly the shift -- which is what a nearest-fragment number
    # is for: a fragment that is not where the verdict says, made visible.
    drawn = {d["project"]: d for d in stats["drawn"]}
    assert drawn["frag_drg"]["verdict"] == "hit"
    assert drawn["frag_drg"]["nearest"] == pytest.approx(0.0, abs=2e-3)
    assert drawn["frag_lig"]["nearest"] == pytest.approx(
        sum(c * c for c in LIG_SHIFT) ** 0.5, abs=2e-3
    )


def test_empty_and_unevaluated_members_are_absent(campaign):
    """``empty`` means somebody looked; no row means nobody has. Neither is
    drawn, and neither is the same as being a hit somewhere else."""
    site = _site(campaign)
    _evaluate(site, "frag_drg", Verdict.HIT)
    _evaluate(site, "frag_apo", Verdict.EMPTY)
    # frag_lig has a ligand but no verdict here: not drawn.

    data = _get(campaign, site).json()
    assert [e["file"] for e in _sticks(data["scene"])] == ["frag_drg"]
    assert {f["name"] for f in data["scene"]["files"]} == {
        "reference", "frag_drg", "frag_drg_dict",
    }
    assert data["stats"]["hits_drawn"] == 1
    assert data["stats"]["empty_verdicts"] == 1


def test_include_unclear_adds_them_in_a_muted_colour(campaign):
    site = _site(campaign)
    _evaluate(site, "frag_drg", Verdict.HIT)
    _evaluate(site, "frag_lig", Verdict.UNCLEAR)

    without = _get(campaign, site).json()
    assert [e["file"] for e in _sticks(without["scene"])] == ["frag_drg"]
    assert without["stats"]["unclear_verdicts"] == 1
    assert without["stats"]["unclear_drawn"] == 0

    with_unclear = _get(campaign, site, "?include=unclear").json()
    sticks = {e["file"]: e for e in _sticks(with_unclear["scene"])}
    assert set(sticks) == {"frag_drg", "frag_lig"}
    assert sticks["frag_lig"]["representations"][0]["colour"] == campaign_scene.UNCLEAR_COLOUR
    assert sticks["frag_drg"]["representations"][0]["colour"] != campaign_scene.UNCLEAR_COLOUR
    assert with_unclear["stats"]["hits_drawn"] == 1
    assert with_unclear["stats"]["unclear_drawn"] == 1
    drawn = {d["project"]: d["verdict"] for d in with_unclear["stats"]["drawn"]}
    assert drawn == {"frag_drg": "hit", "frag_lig": "unclear"}


def test_site_with_no_verdicts_is_exemplar_only(campaign):
    """No fall-back to ligand detection: a site nobody has evaluated shows
    the reference and its pocket, no sticks, and says so."""
    site = _site(campaign)

    data = _get(campaign, site).json()
    scene, stats = data["scene"], data["stats"]
    assert stats["hits_claimed"] == 0
    assert stats["hits_drawn"] == 0
    assert [f["name"] for f in scene["files"]] == ["reference"]
    assert [e["file"] for e in scene["elements"]] == ["reference"]
    assert [r["style"] for r in scene["elements"][0]["representations"]] == ["CRs", "CBs"]
    assert "superpose" not in scene
    assert stats["superpose"] == []
    assert scene["view"]["origin"] == list(SITE_ORIGIN)


def test_hit_verdict_on_a_dataset_with_no_ligand_is_skipped_not_drawn(campaign):
    """A verdict ahead of the modelling: nothing to draw, and the reason and
    the absent nearest distance say so rather than the dataset vanishing."""
    site = _site(campaign)
    _evaluate(site, "frag_apo", Verdict.HIT)

    stats = _get(campaign, site).json()["stats"]
    assert stats["hits_claimed"] == 1
    assert stats["hits_drawn"] == 0
    assert stats["skipped"] == [
        {"project": "frag_apo", "reason": "no ligand in refined model", "nearest": None}
    ]


def test_site_of_another_campaign_is_404(campaign):
    other = models.ProjectGroup.objects.create(
        name="other_campaign", type=models.ProjectGroup.GroupType.FRAGMENT_SET
    )
    foreign_site = _site(other)

    response = _get(campaign, foreign_site)
    assert response.status_code == 404, response.content


def test_origin_in_solvent_gives_no_pocket_and_a_global_fallback(campaign):
    """An empty pocket must produce NO pocket representation and no slab --
    an empty ``selection:`` would draw the whole molecule as sticks. And
    with no CAs within the radius cap, the local fit cannot be made; the
    dataset is still fitted globally, and stats name the fallback."""
    site = _site(campaign, origin=(80.0, 80.0, 80.0))
    _evaluate(site, "frag_lig", Verdict.HIT)

    data = _get(campaign, site).json()
    scene, stats = data["scene"], data["stats"]

    reference = _ribbons(scene)[0]
    assert [r["style"] for r in reference["representations"]] == ["CRs"]
    assert stats["pocket_residues"] == 0
    assert "slab" not in scene["view"]
    assert scene["view"]["origin"] == [80.0, 80.0, 80.0]

    assert len(scene["superpose"]) == 1
    entry = scene["superpose"][0]
    assert "radius" not in entry["fitted"]       # a global fit
    assert entry["vec"] == pytest.approx([-c for c in LIG_SHIFT], abs=2e-3)
    fit = stats["superpose"][0]
    assert fit["ok"] is True
    assert fit["fallback"] == "global"
    assert fit["radius"] is None
    assert "within 30 A" in fit["reason"]

    # Fitted or not, the fragment is nowhere near this site, and the
    # diagnostic says how far.
    assert stats["drawn"][0]["nearest"] > 100


def test_superpose_none_draws_every_frame_as_deposited(campaign):
    site = _site(campaign)
    _evaluate(site, "frag_drg", Verdict.HIT)
    _evaluate(site, "frag_lig", Verdict.HIT)

    data = _get(campaign, site, "?superpose=none").json()
    assert "superpose" not in data["scene"]
    assert data["stats"]["superpose"] == []
    assert len(_sticks(data["scene"])) == 2
    # Unfitted, each fragment is measured where it was deposited: both sit
    # on the origin, because the fixture shifts only frag_lig's protein.
    # (Compare the fitted case above, where the fit carries it off.)
    drawn = {d["project"]: d["nearest"] for d in data["stats"]["drawn"]}
    assert drawn["frag_drg"] == pytest.approx(0.0, abs=2e-3)
    assert drawn["frag_lig"] == pytest.approx(0.0, abs=2e-3)


def test_view_is_the_addressed_site_not_the_campaigns_first(campaign):
    _site(campaign, name="first", origin=(1.0, 2.0, 3.0), order=0)
    second = _site(campaign, name="second", origin=(5.0, 5.0, 5.0), order=1)

    view = _get(campaign, second).json()["scene"]["view"]
    assert view["origin"] == [5.0, 5.0, 5.0]


def test_without_a_parent_the_lowest_id_hit_is_the_exemplar(campaign):
    """No parent model: a hit carries the ribbon AND the pocket, from the
    same files[] entry as its ligand, and the frame is recorded as its."""
    campaign.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.PARENT
    ).delete()
    site = _site(campaign)
    # Recorded in the opposite order to the project ids, to pin the rule.
    _evaluate(site, "frag_lig", Verdict.HIT)
    _evaluate(site, "frag_drg", Verdict.HIT)

    data = _get(campaign, site).json()
    scene, stats = data["scene"], data["stats"]

    assert stats["reference"] == {"project": "frag_drg", "is_parent": False}
    assert stats["parent_present"] is False
    ribbons = _ribbons(scene)
    assert [e["file"] for e in ribbons] == ["frag_drg"]
    styles = [r["style"] for r in ribbons[0]["representations"]]
    assert styles == ["CRs", "CBs", "CBs"]   # ribbon, pocket, its own ligand
    # Its own fragment is not in the pocket list: it is drawn as the hit.
    assert "//A/101" not in ribbons[0]["representations"][1]["selection"]
    assert len([f for f in scene["files"] if f["name"] == "frag_drg"]) == 1

    assert [s["move"] for s in scene["superpose"]] == ["frag_lig"]
    assert scene["superpose"][0]["fitted"]["onto"] == "frag_drg"
    assert scene["view"]["slab"]["file"] == "frag_drg"
    assert stats["hits_drawn"] == 2
