"""
The scene endpoints remember a built scene while nothing it depends on
has changed (``lib/campaign_scene_cache``).

The test settings run with the dummy cache, so these tests switch a real
in-memory backend on for themselves and clear it between tests. They count
builds by wrapping the builder, so that "served from the cache" is a fact
about work not done, not about a response that happens to look the same.

Uses the api/ conftest (transactional DB, AllowAny) -- do NOT add
@pytest.mark.django_db here.
"""
import pytest
from django.core.cache import cache
from django.test import override_settings
from rest_framework.test import APIClient

from ccp4i2.db import models
from ccp4i2.lib import campaign_scene
from ccp4i2.lib import campaign_scene_cache

from .test_summary_scene_api import _build_campaign

LOCMEM = {
    "default": {
        "BACKEND": "django.core.cache.backends.locmem.LocMemCache",
        "LOCATION": "campaign-scene-cache-tests",
    }
}


@pytest.fixture
def campaign(test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    return _build_campaign(test_project_path)


@pytest.fixture
def build_counter(monkeypatch):
    """Count calls to the summary builder without changing what it returns."""
    calls = []
    real = campaign_scene.build_summary_scene

    def counting(group):
        calls.append(group.id)
        return real(group)

    monkeypatch.setattr(campaign_scene, "build_summary_scene", counting)
    return calls


@pytest.fixture
def real_cache():
    with override_settings(CACHES=LOCMEM):
        cache.clear()
        yield cache
        cache.clear()


def _get_summary(client, group, query=""):
    response = client.get(f"/api/ccp4i2/projectgroups/{group.id}/summary_scene/{query}")
    assert response.status_code == 200, response.content
    return response.json()


def test_second_open_is_served_without_rebuilding(
    bypass_api_permissions, campaign, build_counter, real_cache
):
    client = APIClient()
    first = _get_summary(client, campaign)
    second = _get_summary(client, campaign)

    assert build_counter == [campaign.id]
    assert second == first


def test_a_new_verdict_builds_the_scene_afresh(
    bypass_api_permissions, campaign, build_counter, real_cache
):
    client = APIClient()
    _get_summary(client, campaign)

    site = models.CampaignSite.objects.create(
        group=campaign, name="main pocket", origin_x=-5, origin_y=-5, origin_z=-5
    )
    member = campaign.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.MEMBER
    ).first().project
    models.SiteEvaluation.objects.create(
        project=member, site=site, verdict=models.SiteEvaluation.Verdict.HIT
    )

    after = _get_summary(client, campaign)
    assert build_counter == [campaign.id, campaign.id]
    assert [m["text"] for m in after["scene"]["annotations"]] == ["main pocket: 1 hit"]

    # And that state is remembered in its turn.
    _get_summary(client, campaign)
    assert len(build_counter) == 2


def test_a_new_member_builds_the_scene_afresh(
    bypass_api_permissions, campaign, build_counter, real_cache
):
    client = APIClient()
    before = _get_summary(client, campaign)

    newcomer = models.Project.objects.create(name="frag_new", directory="/nowhere")
    models.ProjectGroupMembership.objects.create(
        group=campaign, project=newcomer,
        type=models.ProjectGroupMembership.MembershipType.MEMBER,
    )

    after = _get_summary(client, campaign)
    assert len(build_counter) == 2
    assert after["stats"]["members_total"] == before["stats"]["members_total"] + 1


def test_refresh_rebuilds_on_request(
    bypass_api_permissions, campaign, build_counter, real_cache
):
    client = APIClient()
    _get_summary(client, campaign)
    _get_summary(client, campaign, "?refresh=1")
    assert len(build_counter) == 2


def test_site_scenes_are_keyed_on_their_options(
    bypass_api_permissions, campaign, real_cache, monkeypatch
):
    calls = []
    real = campaign_scene.build_site_scene

    def counting(group, site, **kwargs):
        calls.append(kwargs)
        return real(group, site, **kwargs)

    monkeypatch.setattr(campaign_scene, "build_site_scene", counting)
    site = models.CampaignSite.objects.create(
        group=campaign, name="main pocket", origin_x=-5, origin_y=-5, origin_z=-5
    )
    client = APIClient()
    base = f"/api/ccp4i2/projectgroups/{campaign.id}/sites/{site.id}/scene/"
    assert client.get(base).status_code == 200
    assert client.get(base).status_code == 200
    assert client.get(base + "?include=unclear").status_code == 200
    assert client.get(base + "?superpose=none").status_code == 200

    assert calls == [
        {"include_unclear": False, "superpose": True},
        {"include_unclear": True, "superpose": True},
        {"include_unclear": False, "superpose": False},
    ]


def test_fingerprint_is_stable_while_nothing_changes(campaign):
    assert (
        campaign_scene_cache.scene_fingerprint(campaign)
        == campaign_scene_cache.scene_fingerprint(campaign)
    )


def test_a_broken_cache_still_serves_the_scene(
    bypass_api_permissions, campaign, build_counter, real_cache, monkeypatch
):
    class Broken:
        def get(self, *args, **kwargs):
            raise RuntimeError("cache is down")

        def set(self, *args, **kwargs):
            raise RuntimeError("cache is down")

    # Replace the module's cache object rather than patching methods on
    # Django's proxy: the proxy forwards setattr to whichever backend is
    # current, and monkeypatch's undo can run after override_settings has
    # already switched the backend back.
    monkeypatch.setattr(campaign_scene_cache, "cache", Broken())
    client = APIClient()
    data = _get_summary(client, campaign)
    assert data["stats"]["hits"] == 2
    assert build_counter == [campaign.id]
