"""The sweep_staged_uploads management command."""

import pytest
from django.core.management import call_command
from django.utils import timezone

from ccp4i2.db import models


@pytest.mark.django_db
def test_command_reaps_expired(tmp_path, monkeypatch, settings, capsys):
    d = tmp_path / "staging"
    d.mkdir()
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(d))
    settings.CCP4I2_IMPORT_STAGING_TTL_HOURS = 24

    old = models.StagedUpload.objects.create(
        owner="7", filename="m.mrc", size_bytes=8,
        state=models.StagedUpload.State.STAGING)
    old.created_at = timezone.now() - timezone.timedelta(hours=48)
    old.save()
    fresh = models.StagedUpload.objects.create(
        owner="7", filename="m.mrc", size_bytes=8,
        state=models.StagedUpload.State.STAGING)

    call_command("sweep_staged_uploads")
    out = capsys.readouterr().out
    assert "Swept 1" in out
    assert models.StagedUpload.objects.filter(uuid=fresh.uuid).exists()
    assert not models.StagedUpload.objects.filter(uuid=old.uuid).exists()


@pytest.mark.django_db
def test_command_noop_when_disabled(monkeypatch, capsys):
    monkeypatch.delenv("CCP4I2_IMPORT_STAGING_DIR", raising=False)
    call_command("sweep_staged_uploads")
    assert "not enabled" in capsys.readouterr().out
