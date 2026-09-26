"""GET /version/ declares the deployment's run targets (service contract, ccp4i2-api 0.6.0)."""
import pytest
from django.conf import settings
from django.test import Client

from ccp4i2.lib import dispatch

VERSION_URL = "/api/ccp4i2/version/"


@pytest.fixture
def client():
    return Client()


def _version(client):
    r = client.get(VERSION_URL)
    assert r.status_code == 200, r.content
    return r.json()


def test_a_desktop_lists_only_local(client, monkeypatch):
    monkeypatch.delattr(settings, "CCP4I2_RUN_TARGETS", raising=False)
    monkeypatch.delattr(settings, "CCP4I2_JOB_TARGET", raising=False)
    monkeypatch.delenv("CCP4I2_JOB_TARGET", raising=False)
    j = _version(client)
    assert j["job_target"] == "local"
    assert j["run_targets"] == [{"name": "local", "path": dispatch.DEFAULT_RUN_TARGETS["local"],
                                 "runs_jobs": True, "runs_programs": False}]


class FakeBatch:
    def submit(self, *a): return "h"
    def poll(self, h): return "queued"
    def cancel(self, h): return None
    def logs(self, h): return None


def test_a_deployment_with_a_program_target_says_so(client, monkeypatch):
    monkeypatch.setattr(settings, "CCP4I2_RUN_TARGETS",
                        {**dispatch.DEFAULT_RUN_TARGETS, "batch": f"{__name__}.FakeBatch",
                         "broken": "no.such.Thing"}, raising=False)
    monkeypatch.setattr(settings, "CCP4I2_JOB_TARGET", "local", raising=False)
    by_name = {t["name"]: t for t in _version(client)["run_targets"]}
    assert by_name["batch"]["runs_programs"] is True and by_name["batch"]["runs_jobs"] is False
    assert "error" in by_name["broken"] and by_name["broken"]["runs_programs"] is False
