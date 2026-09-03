"""The checksum of what was uploaded, as distinct from what was imported.

`FileImport.checksum` records the file as it ended up in the project. For
anything derived --- an MTZ split down to the columns one parameter asked for
--- that is the checksum of the derivative, so it cannot answer "have I seen
this file before?". `source_checksum` records the bytes as uploaded, taken
before any splitting or conversion.

A hash not taken at upload time cannot be recovered afterwards, which is why
it is recorded even though nothing yet acts on it automatically.

Pure Python -- no CCP4 binaries and no database needed.
"""
import hashlib

from ccp4i2.lib.utils.files.upload_param import (
    _checksum_of, _existing_imports_of_same_source,
)


def test_checksum_matches_md5_of_the_bytes(tmp_path):
    """Same algorithm as CDataFile.checksum(), so the two are comparable."""
    target = tmp_path / "gamma.mtz"
    target.write_bytes(b"MTZ\x00not really" * 1000)
    assert _checksum_of(target) == hashlib.md5(target.read_bytes()).hexdigest()


def test_checksum_is_chunked_not_slurped(tmp_path):
    """Reflection files run to hundreds of megabytes; correctness across the
    chunk boundary is the whole point of reading them piecewise."""
    target = tmp_path / "big.mtz"
    payload = bytes(range(256)) * 4096  # 1 MiB, spans the 64 KiB chunks
    target.write_bytes(payload)
    assert _checksum_of(target) == hashlib.md5(payload).hexdigest()


def test_an_unreadable_file_yields_no_checksum_rather_than_raising(tmp_path):
    """A hash we could not take must degrade to "this file takes no part in
    duplicate detection", never to a failed upload."""
    assert _checksum_of(tmp_path / "was-never-there.pdb") == ""


def test_a_blank_checksum_matches_nothing(tmp_path):
    """Files whose source could not be hashed --- and every row written before
    source_checksum existed --- must not all collide with each other on "".

    Short-circuits before any query, so no database is involved.
    """
    class _NoJob:
        def __getattr__(self, name):
            raise AssertionError("must not reach the database for a blank checksum")

    assert _existing_imports_of_same_source(_NoJob(), "") == []
