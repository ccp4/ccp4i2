"""Staged import transport for web deployments.

A browser delivers a large file into ``CCP4I2_IMPORT_STAGING_DIR`` in chunks that
stay under every body cap in the deployment (Next middleware, Django, the
ingress), then imports it *by handle* -- an owner-bound, server-generated uuid --
past those caps. This module is the server side of that transport and the gate
that turns the a63 "containment" check into real authorisation.

Security invariants (see docs/staged-import.md):

* The server names every path. The directory is the row's ``uuid`` and the
  filename is sanitised here; no client string is ever joined into a path.
* A staged file is owner-bound and importable only by that owner.
* Disk is bounded: per-upload size, per-owner in-flight count, and expiry.
* Failure is loud: each refusal raises a typed error the view maps to a status
  code, never a silent fall-through.

The functions here do not touch the request; the view passes an ``owner`` key
(``str(request.user.pk)``) and the raw values.
"""

import hashlib
import logging
import shutil
from datetime import timedelta
from pathlib import Path

from django.conf import settings
from django.utils import timezone
from django.utils.text import get_valid_filename

from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


# --- typed errors -> HTTP status in the view -------------------------------

class StagedUploadError(Exception):
    """Base: carries an HTTP status and a message."""
    status = 400

    def __init__(self, message):
        super().__init__(message)
        self.message = message


class TooLarge(StagedUploadError):          # 413
    status = 413


class TooManyInFlight(StagedUploadError):   # 429
    status = 429


class NotFound(StagedUploadError):          # 404 (also foreign owner)
    status = 404


class Expired(StagedUploadError):           # 410
    status = 410


class ChunksMissing(StagedUploadError):     # 409
    status = 409


class HashMismatch(StagedUploadError):      # 422
    status = 422


# --- configuration accessors -----------------------------------------------

def staging_dir():
    """The staging root as a resolved ``Path``, or ``None`` when unset."""
    import os
    raw = os.environ.get("CCP4I2_IMPORT_STAGING_DIR")
    if not raw:
        return None
    try:
        return Path(raw).resolve()
    except (OSError, ValueError):
        return None


def staging_enabled():
    return staging_dir() is not None


def owner_key(request):
    """A stable per-user key for owner-binding a staged upload.

    In a real deployment ``request.user`` is authenticated (a pk is set) by the
    local-session (desktop) or JWT (cloud) middleware. Tests without auth may
    leave it ``None`` or anonymous; those collapse to the string ``"None"``,
    which is still internally consistent (the same request stages and imports).
    """
    user = getattr(request, "user", None)
    return str(getattr(user, "pk", None))


# Read config from settings with a fallback default, so the module works under
# any settings module (the test settings do not define these).
_DEFAULTS = {
    "CCP4I2_IMPORT_STAGING_CHUNK_BYTES": 16 * 1024 * 1024,
    "CCP4I2_IMPORT_STAGING_MAX_BYTES": 2 * 1024 * 1024 * 1024,
    "CCP4I2_IMPORT_STAGING_TTL_HOURS": 24,
    "CCP4I2_IMPORT_STAGING_THRESHOLD_BYTES": 32 * 1024 * 1024,
    "CCP4I2_IMPORT_STAGING_MAX_INFLIGHT": 8,
}


def _conf(name):
    return getattr(settings, name, _DEFAULTS[name])


def chunk_bytes():
    return _conf("CCP4I2_IMPORT_STAGING_CHUNK_BYTES")


def ttl_hours():
    return _conf("CCP4I2_IMPORT_STAGING_TTL_HOURS")


def capability():
    """The ``import_staging`` capability advertised to the client, or ``None``."""
    if not staging_enabled():
        return None
    return {
        "chunk_bytes": _conf("CCP4I2_IMPORT_STAGING_CHUNK_BYTES"),
        "max_bytes": _conf("CCP4I2_IMPORT_STAGING_MAX_BYTES"),
        "threshold_bytes": _conf("CCP4I2_IMPORT_STAGING_THRESHOLD_BYTES"),
    }


# --- path helpers (server names everything) --------------------------------

def _row_dir(row):
    return staging_dir() / str(row.uuid)


def _part_path(row, index: int):
    return _row_dir(row) / f"part.{int(index)}"


def final_path(row):
    """The assembled file's path. ``row.filename`` is already sanitised."""
    return _row_dir(row) / row.filename


# --- lifecycle --------------------------------------------------------------

def begin(owner: str, filename: str, size_bytes: int, sha256: str = ""):
    """Create a staging row + its directory. Raises on the size / in-flight bounds."""
    if staging_dir() is None:
        raise NotFound("Staged upload is not enabled on this deployment")
    if size_bytes is None or int(size_bytes) < 0:
        raise StagedUploadError("size_bytes is required and must be >= 0")
    max_bytes = _conf("CCP4I2_IMPORT_STAGING_MAX_BYTES")
    if int(size_bytes) > max_bytes:
        raise TooLarge(
            f"File is larger than the {max_bytes} byte staging limit")
    inflight = models.StagedUpload.objects.filter(
        owner=owner, state=models.StagedUpload.State.STAGING).count()
    if inflight >= _conf("CCP4I2_IMPORT_STAGING_MAX_INFLIGHT"):
        raise TooManyInFlight("Too many uploads already in progress")

    # Sanitise the filename to a plain basename; never trust it for the path.
    safe = get_valid_filename(Path(str(filename or "upload")).name) or "upload"
    row = models.StagedUpload.objects.create(
        owner=owner, filename=safe, size_bytes=int(size_bytes),
        sha256=(sha256 or "").lower().strip(),
    )
    _row_dir(row).mkdir(parents=True, exist_ok=True)
    logger.info("staged upload begin %s owner=%s size=%s", row.uuid, owner, size_bytes)
    return row


def get_owned(uuid, owner: str):
    """The row with this uuid owned by ``owner``, or raise NotFound.

    A foreign or unknown id is a 404 either way -- we never reveal that an id
    exists but belongs to someone else.
    """
    row = models.StagedUpload.objects.filter(uuid=uuid, owner=owner).first()
    if row is None:
        raise NotFound("No such staged upload")
    return row


def _check_not_expired(row):
    horizon = timedelta(hours=_conf("CCP4I2_IMPORT_STAGING_TTL_HOURS"))
    if timezone.now() - row.created_at > horizon:
        raise Expired("This staged upload has expired")


def write_chunk(row, index: int, data: bytes):
    """Write one chunk. Idempotent: a repeated index overwrites."""
    if row.state != models.StagedUpload.State.STAGING:
        raise NotFound("This upload is no longer accepting chunks")
    _check_not_expired(row)
    if len(data) > _conf("CCP4I2_IMPORT_STAGING_CHUNK_BYTES"):
        raise TooLarge("Chunk is larger than the negotiated chunk size")
    path = _part_path(row, index)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)


def present_indexes(row):
    """Sorted chunk indexes already on disk (for resume)."""
    d = _row_dir(row)
    if not d.is_dir():
        return []
    out = []
    for p in d.iterdir():
        if p.name.startswith("part."):
            try:
                out.append(int(p.name.split(".", 1)[1]))
            except ValueError:
                pass
    return sorted(out)


def finish(row):
    """Assemble the parts into the final file, verify, and mark ``ready``.

    Concatenates ``part.0..part.N`` in order; the client sends contiguous indexes
    from 0. Verifies the total byte count (and the sha256 if one was given),
    deletes the parts, and moves the row to ``ready``.
    """
    if row.state == models.StagedUpload.State.READY:
        return row  # idempotent
    if row.state != models.StagedUpload.State.STAGING:
        raise NotFound("This upload cannot be finished")
    _check_not_expired(row)

    indexes = present_indexes(row)
    if not indexes or indexes != list(range(len(indexes))):
        raise ChunksMissing(
            f"Missing chunks: have {indexes}, need a contiguous run from 0")

    out = final_path(row)
    hasher = hashlib.sha256()
    total = 0
    with open(out, "wb") as dst:
        for i in indexes:
            with open(_part_path(row, i), "rb") as src:
                for block in iter(lambda: src.read(1024 * 1024), b""):
                    dst.write(block)
                    hasher.update(block)
                    total += len(block)

    if total != row.size_bytes:
        out.unlink(missing_ok=True)
        raise HashMismatch(
            f"Assembled size {total} != declared {row.size_bytes}")
    if row.sha256 and hasher.hexdigest() != row.sha256:
        out.unlink(missing_ok=True)
        raise HashMismatch("Assembled sha256 does not match the declared hash")

    for i in indexes:
        _part_path(row, i).unlink(missing_ok=True)
    row.state = models.StagedUpload.State.READY
    row.save(update_fields=["state"])
    logger.info("staged upload ready %s (%s bytes)", row.uuid, total)
    return row


def resolve_for_import(uuid, owner: str):
    """The path to import from a ``ready`` handle, owner-checked and contained.

    Returns the assembled file ``Path``, or raises. This is the cloud analogue of
    ``resolve_importable_path`` -- but the caller cannot name a path, only a
    handle it owns, which closes the "any file in the staging dir" hole.
    """
    row = get_owned(uuid, owner)
    if row.state != models.StagedUpload.State.READY:
        raise NotFound("Staged upload is not ready to import")
    _check_not_expired(row)
    path = final_path(row)
    root = staging_dir()
    # Defence in depth: the assembled path must still be inside the staging root.
    if not (path.resolve().is_relative_to(root) and path.is_file()):
        raise NotFound("Staged file is missing")
    return row, path


def consume(row, delete=True):
    """Mark a row imported so its handle can't be reused.

    ``delete`` removes the staging directory now -- correct when the bytes have
    already been copied into the project (the synchronous upload_file_param
    path). Pass ``delete=False`` when a *detached* importer still needs to read
    the file (project-zip import): the row is marked consumed immediately, and
    the sweeper reaps the directory once it is past the TTL.
    """
    row.state = models.StagedUpload.State.CONSUMED
    row.save(update_fields=["state"])
    if delete:
        shutil.rmtree(_row_dir(row), ignore_errors=True)
    logger.info("staged upload consumed %s (delete=%s)", row.uuid, delete)


def sweep():
    """Reap expired rows (any state) and consumed rows whose file is gone.

    Returns the number of rows removed. Idempotent; safe to run on a schedule.
    """
    if staging_dir() is None:
        return 0
    horizon = timezone.now() - timedelta(
        hours=_conf("CCP4I2_IMPORT_STAGING_TTL_HOURS"))
    removed = 0
    for row in models.StagedUpload.objects.all():
        gone = (row.state == models.StagedUpload.State.CONSUMED
                and not _row_dir(row).exists())
        if row.created_at < horizon or gone:
            shutil.rmtree(_row_dir(row), ignore_errors=True)
            row.delete()
            removed += 1
    return removed
