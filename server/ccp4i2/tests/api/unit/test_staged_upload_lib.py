"""The staged-upload transport (issue: staged import for web).

Covers the lifecycle and every refusal, at the library level (the view is a thin
map to these). Owner isolation is checked by passing distinct owner keys, exactly
as the view derives them from request.user.pk.
"""

import hashlib

import pytest
from django.utils import timezone

from ccp4i2.db import models
from ccp4i2.lib.utils.files import staged_upload as su


@pytest.fixture
def staging(tmp_path, monkeypatch, settings):
    d = tmp_path / "staging"
    d.mkdir()
    monkeypatch.setenv("CCP4I2_IMPORT_STAGING_DIR", str(d))
    settings.CCP4I2_IMPORT_STAGING_CHUNK_BYTES = 8
    settings.CCP4I2_IMPORT_STAGING_MAX_BYTES = 1000
    settings.CCP4I2_IMPORT_STAGING_THRESHOLD_BYTES = 4
    settings.CCP4I2_IMPORT_STAGING_TTL_HOURS = 24
    settings.CCP4I2_IMPORT_STAGING_MAX_INFLIGHT = 3
    return d


def _stage(owner, data, sha256="", filename="map.mrc"):
    row = su.begin(owner, filename, len(data), sha256)
    # 8-byte chunks
    for i, off in enumerate(range(0, len(data), 8)):
        su.write_chunk(row, i, data[off:off + 8])
    return su.finish(su.get_owned(row.uuid, owner))


def test_happy_path_with_hash(staging):
    data = b"cryo-em-map-bytes-x" * 3
    h = hashlib.sha256(data).hexdigest()
    row = _stage("7", data, sha256=h)
    assert row.state == models.StagedUpload.State.READY
    _, path = su.resolve_for_import(row.uuid, "7")
    assert path.read_bytes() == data


def test_repeated_and_out_of_order_chunks_are_idempotent(staging):
    data = b"abcdefghijABCDEFGHIJ"   # 20 bytes -> 3 chunks (8,8,4)
    row = su.begin("7", "m.mrc", len(data), "")
    su.write_chunk(row, 2, data[16:20])       # out of order
    su.write_chunk(row, 0, data[0:8])
    su.write_chunk(row, 0, data[0:8])         # repeat -> overwrite
    su.write_chunk(row, 1, data[8:16])
    assert su.present_indexes(row) == [0, 1, 2]
    assert su.finish(row).state == models.StagedUpload.State.READY
    _, path = su.resolve_for_import(row.uuid, "7")
    assert path.read_bytes() == data


def test_missing_chunk_is_409(staging):
    row = su.begin("7", "m.mrc", 20, "")
    su.write_chunk(row, 0, b"abcdefgh")
    su.write_chunk(row, 2, b"xyz")            # index 1 missing
    with pytest.raises(su.ChunksMissing):
        su.finish(row)


def test_bad_hash_is_422(staging):
    data = b"twelve bytes"
    row = su.begin("7", "m.mrc", len(data), sha256="0" * 64)
    su.write_chunk(row, 0, data[:8])
    su.write_chunk(row, 1, data[8:])
    with pytest.raises(su.HashMismatch):
        su.finish(row)


def test_size_mismatch_is_422(staging):
    row = su.begin("7", "m.mrc", 99, "")      # declares 99, sends 8
    su.write_chunk(row, 0, b"eightby.")
    with pytest.raises(su.HashMismatch):
        su.finish(row)


def test_expiry_is_410(staging):
    row = su.begin("7", "m.mrc", 8, "")
    row.created_at = timezone.now() - timezone.timedelta(hours=25)
    row.save()
    with pytest.raises(su.Expired):
        su.write_chunk(row, 0, b"eightby.")


def test_over_max_is_413_at_begin(staging):
    with pytest.raises(su.TooLarge):
        su.begin("7", "big.mrc", 10_000, "")   # > MAX_BYTES=1000


def test_too_many_in_flight_is_429(staging):
    for _ in range(3):                          # MAX_INFLIGHT=3
        su.begin("7", "m.mrc", 8, "")
    with pytest.raises(su.TooManyInFlight):
        su.begin("7", "m.mrc", 8, "")
    # a different owner is unaffected
    assert su.begin("8", "m.mrc", 8, "")


def test_owner_isolation(staging):
    row = _stage("7", b"owned by 7 only!!", "")
    # user 8 can neither see nor import user 7's handle
    with pytest.raises(su.NotFound):
        su.get_owned(row.uuid, "8")
    with pytest.raises(su.NotFound):
        su.resolve_for_import(row.uuid, "8")


def test_consume_deletes_and_blocks_reuse(staging):
    row = _stage("7", b"import me once!!!", "")
    _, path = su.resolve_for_import(row.uuid, "7")
    assert path.is_file()
    su.consume(row)
    assert row.state == models.StagedUpload.State.CONSUMED
    assert not path.exists()
    with pytest.raises(su.NotFound):        # a consumed handle can't import again
        su.resolve_for_import(row.uuid, "7")


def test_sweep_reaps_expired_and_consumed(staging):
    old = su.begin("7", "m.mrc", 8, "")
    old.created_at = timezone.now() - timezone.timedelta(hours=48)
    old.save()
    fresh = su.begin("7", "m.mrc", 8, "")
    removed = su.sweep()
    assert removed == 1
    assert models.StagedUpload.objects.filter(uuid=fresh.uuid).exists()
    assert not models.StagedUpload.objects.filter(uuid=old.uuid).exists()
