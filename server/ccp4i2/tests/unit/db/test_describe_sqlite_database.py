"""Tell a legacy (Qt) CCP4i2 database from one written by this application.

The legacy importer reads the Qt schema. Pointed at a Django database — one
from a previous alpha of this application — every query fails on a missing
table, and the user gets a stack of OperationalErrors instead of being told
they picked the wrong file. Those projects also carry their own snapshots on
disk, so the folder route handles them properly and this one never should.

The two schemas share no table name, so one lookup settles it.
"""

import sqlite3

from ccp4i2.db.import_sqlite import describe_sqlite_database


def _make(path, tables):
    connection = sqlite3.connect(path)
    for table in tables:
        connection.execute(f"CREATE TABLE [{table}] (id INTEGER PRIMARY KEY)")
    connection.commit()
    connection.close()
    return path


def test_recognises_a_legacy_database(tmp_path):
    db = _make(tmp_path / "legacy.sqlite", ["Projects", "Jobs", "Files"])
    described = describe_sqlite_database(db)
    assert described["kind"] == "legacy"
    assert described["projects"] == 0


def test_recognises_this_applications_database(tmp_path):
    db = _make(tmp_path / "django.sqlite3", ["ccp4i2_project", "ccp4i2_job"])
    described = describe_sqlite_database(db)
    assert described["kind"] == "django"
    # The message must point somewhere useful, not merely say no.
    assert "folder" in described["message"].lower()


def test_a_file_that_is_not_a_database(tmp_path):
    """sqlite3.connect opens lazily, so this only fails on the first read —
    inside the function, where it must not escape."""
    plain = tmp_path / "notes.txt"
    plain.write_text("this is not a database")
    assert describe_sqlite_database(plain)["kind"] == "unknown"


def test_a_missing_file(tmp_path):
    described = describe_sqlite_database(tmp_path / "absent.sqlite")
    assert described["kind"] == "unknown"
    assert "No such file" in described["message"]


def test_an_empty_database(tmp_path):
    assert describe_sqlite_database(_make(tmp_path / "e.sqlite", []))["kind"] == "empty"


def test_an_unrelated_database(tmp_path):
    db = _make(tmp_path / "other.sqlite", ["widgets"])
    assert describe_sqlite_database(db)["kind"] == "unknown"
