"""
E2E tests for the several-sites campaign scene endpoint.

    GET /api/ccp4i2/projectgroups/{id}/sites/scene/?ids=<id>,<id>,...

For a site and its subsites, or neighbouring sites worth seeing together.
Same synthetic campaign as ``test_site_scene_api``: a poly-Ala helix along z
from about 5 A to 33 A, the parent plus two hit members and an apo member.
Site A sits on the members' ligand at (5, 5, 5); site B is 20 A up the
helix at (5, 5, 25), so the pair spans 10 A either side of their centroid
and the one fit they share has to start wider than a single site's.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""

import pytest
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.lib import campaign_scene

from .test_site_scene_api import (
    SITE_ORIGIN, Verdict, _evaluate, _site, _sticks,
)
from .test_summary_scene_api import _build_campaign

# A second site further up the helix, stored the way the app stores it.
SITE_B_POSITION = (5.0, 5.0, 25.0)
SITE_B_ORIGIN = tuple(-c for c in SITE_B_POSITION)


@pytest.fixture
def campaign(bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    return _build_campaign(test_project_path)


@pytest.fixture
def sites(campaign):
    a = _site(campaign, name="pocket_A", quat=[0.0, 0.0, 0.0, 1.0], zoom=0.35)
    b = _site(campaign, name="subsite_B", origin=SITE_B_ORIGIN, order=1)
    return a, b


def _get_many(group, ids, query=""):
    joined = ",".join(str(i) for i in ids)
    return APIClient().get(
        f"/api/ccp4i2/projectgroups/{group.id}/sites/scene/?ids={joined}{query}"
    )


def test_hits_at_either_site_are_drawn_once_each(campaign, sites):
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    _evaluate(b, "frag_lig", Verdict.HIT)
    _evaluate(a, "frag_apo", Verdict.EMPTY)
    _evaluate(b, "frag_apo", Verdict.EMPTY)

    response = _get_many(campaign, [a.id, b.id])
    assert response.status_code == 200, response.content
    data = response.json()
    scene, stats = data["scene"], data["stats"]

    assert stats["sites"] == [
        {"id": a.id, "name": "pocket_A"}, {"id": b.id, "name": "subsite_B"},
    ]
    assert "site" not in stats
    assert {e["file"] for e in _sticks(scene)} == {"frag_drg", "frag_lig"}
    assert stats["hits_claimed"] == 2
    assert stats["hits_drawn"] == 2
    # One dataset, empty at both: counted once, not twice.
    assert stats["empty_verdicts"] == 1

    drawn = {d["project"]: d for d in stats["drawn"]}
    assert drawn["frag_drg"]["sites"] == ["pocket_A"]
    assert drawn["frag_lig"]["sites"] == ["subsite_B"]
    assert scene["scene"].endswith("sites pocket_A + subsite_B")


def test_a_dataset_hit_at_both_is_drawn_once_naming_both(campaign, sites):
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    _evaluate(b, "frag_drg", Verdict.HIT)

    data = _get_many(campaign, [a.id, b.id]).json()
    assert [e["file"] for e in _sticks(data["scene"])] == ["frag_drg"]
    assert data["stats"]["hits_claimed"] == 1
    assert data["stats"]["drawn"][0]["sites"] == ["pocket_A", "subsite_B"]


def test_the_strongest_verdict_across_the_sites_decides(campaign, sites):
    """Empty in one subsite and a hit in the next is a hit of the view;
    unclear in one and empty in the other is unclear."""
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.EMPTY)
    _evaluate(b, "frag_drg", Verdict.HIT)
    _evaluate(a, "frag_lig", Verdict.UNCLEAR)
    _evaluate(b, "frag_lig", Verdict.EMPTY)

    plain = _get_many(campaign, [a.id, b.id]).json()
    assert [e["file"] for e in _sticks(plain["scene"])] == ["frag_drg"]
    assert plain["stats"]["drawn"][0]["sites"] == ["subsite_B"]
    assert plain["stats"]["unclear_verdicts"] == 1
    assert plain["stats"]["empty_verdicts"] == 0

    with_unclear = _get_many(campaign, [a.id, b.id], "&include=unclear").json()
    sticks = {e["file"]: e for e in _sticks(with_unclear["scene"])}
    assert set(sticks) == {"frag_drg", "frag_lig"}
    assert sticks["frag_lig"]["representations"][0]["colour"] == campaign_scene.UNCLEAR_COLOUR
    drawn = {d["project"]: d for d in with_unclear["stats"]["drawn"]}
    assert drawn["frag_lig"]["sites"] == ["pocket_A"]


def test_one_fit_on_the_region_and_one_pocket_for_both(campaign, sites):
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    _evaluate(b, "frag_lig", Verdict.HIT)

    data = _get_many(campaign, [a.id, b.id]).json()
    scene, stats = data["scene"], data["stats"]

    # The sites are 10 A either side of their centroid, so the shared sphere
    # starts at 10 + ENVIRONMENT_RADIUS = 18 A rather than a single site's
    # 15 A, and every hit is fitted on it.
    expected_radius = 10.0 + campaign_scene.ENVIRONMENT_RADIUS
    fits = {s["project"]: s for s in stats["superpose"]}
    assert set(fits) == {"frag_drg", "frag_lig"}
    for fit in fits.values():
        assert fit["ok"] and "fallback" not in fit
        assert fit["radius"] == pytest.approx(expected_radius)
    for entry in scene["superpose"]:
        assert entry["fitted"]["radius"] == pytest.approx(expected_radius)

    # The pocket is every residue around either site, once, in helix order.
    single_a = campaign_scene.build_site_scene(campaign, a)["scene"]
    single_b = campaign_scene.build_site_scene(campaign, b)["scene"]

    def pocket(scene):
        reps = [e for e in scene["elements"] if e["file"] == "reference"][0]["representations"]
        return reps[1]["selection"].split("||")

    union = pocket(scene)
    assert set(union) == set(pocket(single_a)) | set(pocket(single_b))
    assert len(union) == len(set(union))
    assert union == sorted(union, key=lambda cid: int(cid.split("/")[-1]))
    assert stats["pocket_residues"] == len(union)
    assert scene["view"]["slab"]["selection"] == "||".join(union)


def test_view_looks_at_the_centroid_the_way_the_first_site_was_saved(campaign, sites):
    a, b = sites
    data = _get_many(campaign, [a.id, b.id]).json()
    view = data["scene"]["view"]
    # Moorhen's convention, averaged as stored: the negated centroid.
    assert view["origin"] == pytest.approx(
        [(p + q) / 2 for p, q in zip(SITE_ORIGIN, SITE_B_ORIGIN)]
    )
    assert view["quat"] == [0.0, 0.0, 0.0, 1.0]
    assert view["zoom"] == 0.35

    # The order of the ids is the order of the view: B first takes B's
    # (unsaved) orientation and leads the title.
    swapped = _get_many(campaign, [b.id, a.id]).json()
    assert "quat" not in swapped["scene"]["view"]
    assert swapped["scene"]["scene"].endswith("sites subsite_B + pocket_A")


def test_one_id_is_exactly_the_single_site_scene(campaign, sites):
    a, _ = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    _evaluate(a, "frag_lig", Verdict.HIT)

    many = _get_many(campaign, [a.id]).json()
    one = APIClient().get(
        f"/api/ccp4i2/projectgroups/{campaign.id}/sites/{a.id}/scene/"
    ).json()
    assert many == one
    assert one["stats"]["site"] == {"id": a.id, "name": "pocket_A"}
    assert "sites" not in one["stats"]


def test_repeated_ids_count_once(campaign, sites):
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    data = _get_many(campaign, [a.id, b.id, a.id]).json()
    assert [s["name"] for s in data["stats"]["sites"]] == ["pocket_A", "subsite_B"]


@pytest.mark.parametrize("query", ["", "?ids=", "?ids=1,x", "?ids=,,"])
def test_missing_or_malformed_ids_are_400(campaign, query):
    response = APIClient().get(
        f"/api/ccp4i2/projectgroups/{campaign.id}/sites/scene/{query}"
    )
    assert response.status_code == 400


def test_a_site_of_another_campaign_is_404(campaign, sites, test_project_path):
    a, _ = sites
    other = models.ProjectGroup.objects.create(name="other_campaign")
    foreign = _site(other, name="elsewhere")
    response = _get_many(campaign, [a.id, foreign.id])
    assert response.status_code == 404
    assert str(foreign.id) in response.json()["error"]


def test_hits_measure_nearest_from_the_site_that_claims_them(campaign, sites):
    """The diagnostic is distance to the claiming site, not to the centroid:
    frag_drg's fragment sits on site A, 10 A from the middle of the pair."""
    a, b = sites
    _evaluate(a, "frag_drg", Verdict.HIT)
    data = _get_many(campaign, [a.id, b.id]).json()
    drawn = {d["project"]: d for d in data["stats"]["drawn"]}
    assert drawn["frag_drg"]["nearest"] == pytest.approx(0.0, abs=2e-3)
