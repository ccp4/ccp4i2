"""API tests for program-location config + Tips-of-the-Day endpoints."""

import json

import pytest
from rest_framework.test import APIClient


@pytest.fixture
def client(bypass_api_permissions):
    return APIClient()


def test_discover_programs_explicit_names(client):
    resp = client.get("/api/ccp4i2/config/discover-programs/?names=ls,zzz_nope")
    assert resp.status_code == 200
    programs = {p["name"]: p for p in resp.json()["data"]["programs"]}
    assert programs["ls"]["path"] is not None
    assert programs["ls"]["source"] == "path"
    assert programs["zzz_nope"]["path"] is None
    assert programs["zzz_nope"]["source"] == "missing"


def test_discover_programs_default_is_registry_derived(client):
    """With no ?names, the list comes from real plugin TASKCOMMANDs, each
    annotated with the tasks that use it — not a hardcoded guess."""
    resp = client.get("/api/ccp4i2/config/discover-programs/")
    assert resp.status_code == 200
    programs = resp.json()["data"]["programs"]
    assert len(programs) > 20
    by_name = {p["name"]: p for p in programs}
    # A real executable name (ctruncate) is present and task-annotated...
    assert "ctruncate" in by_name
    assert "ctruncate" in by_name["ctruncate"]["tasks"]
    # ...and core interpreters / abs-path commands are excluded from the default.
    assert "ccp4-python" not in by_name
    assert all("/" not in name for name in by_name)


def test_program_preferences_desktop_roundtrip(client, monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")

    # set
    resp = client.patch(
        "/api/ccp4i2/config/program-preferences/set/",
        data=json.dumps({"SHELXDIR": "/opt/shelx", "bogus": "x"}),
        content_type="application/json",
    )
    assert resp.status_code == 200
    assert resp.json()["data"]["preferences"] == {"SHELXDIR": "/opt/shelx"}

    # get reflects it, editable True on desktop
    resp = client.get("/api/ccp4i2/config/program-preferences/")
    body = resp.json()["data"]
    assert body["editable"] is True
    assert body["preferences"]["SHELXDIR"] == "/opt/shelx"


def test_program_preferences_cloud_readonly(client, monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)  # cloud

    resp = client.patch(
        "/api/ccp4i2/config/program-preferences/set/",
        data=json.dumps({"SHELXDIR": "/opt/shelx"}),
        content_type="application/json",
    )
    assert resp.status_code == 409
    assert resp.json()["success"] is False

    resp = client.get("/api/ccp4i2/config/program-preferences/")
    assert resp.json()["data"]["editable"] is False


def test_tip_of_the_day(client):
    resp = client.get("/api/ccp4i2/tips/")
    assert resp.status_code == 200
    data = resp.json()["data"]
    assert data["count"] > 0
    assert isinstance(data["id"], int)
    assert "<" in data["html"]  # some HTML body
    # image srcs rewritten to the tips-image route
    if "img" in data["html"]:
        assert "/api/proxy/ccp4i2/tips/image/" in data["html"]


def test_tip_specific_id(client):
    resp = client.get("/api/ccp4i2/tips/?id=1")
    assert resp.status_code == 200
    assert resp.json()["data"]["id"] == 1


def test_tip_image_traversal_blocked(client):
    resp = client.get("/api/ccp4i2/tips/image/..%2fviews.py")
    assert resp.status_code == 404
