"""
The SQLite connection must wait for a contended write lock, not fail.

Several processes write the one database: two uvicorn workers and a
ccp4-python subprocess per running job. With Django's default deferred
transaction, a write that follows a read while another process holds the
write lock fails at once -- SQLite does not consult the busy handler in that
case -- and a New Project import lost jobs and create_task calls to
"database is locked" half a second apart. transaction_mode=IMMEDIATE takes
the write lock at BEGIN, so the transaction queues instead.

Pure Python -- no CCP4, no Django connection needed.
"""
import sqlite3
import threading
import time

import pytest

from ccp4i2.config.sqlite import SQLITE_INIT_COMMAND, sqlite_database


def test_the_settings_say_so():
    db = sqlite_database("/tmp/x.sqlite3")
    assert db["ENGINE"] == "django.db.backends.sqlite3"
    opts = db["OPTIONS"]
    assert opts["transaction_mode"] == "IMMEDIATE"
    assert opts["timeout"] >= 30
    assert "journal_mode=WAL" in opts["init_command"]
    assert "busy_timeout=30000" in opts["init_command"]
    assert opts["init_command"] == SQLITE_INIT_COMMAND


def _contended_write(path, begin):
    """Another process holds the write lock for a second; we read, then write.

    Returns (succeeded, seconds waited).
    """
    def connect():
        c = sqlite3.connect(path, isolation_level=None, timeout=0, check_same_thread=False)
        for stmt in SQLITE_INIT_COMMAND.split(";"):
            if stmt.strip():
                c.execute(stmt)
        return c

    holder = connect()
    holder.execute("BEGIN IMMEDIATE")
    holder.execute("UPDATE job SET status = 3 WHERE id = 1")
    threading.Timer(1.0, lambda: holder.execute("COMMIT")).start()

    other = connect()
    started = time.monotonic()
    try:
        other.execute(begin)
        other.execute("SELECT status FROM job WHERE id = 1").fetchone()
        other.execute("UPDATE job SET status = 2 WHERE id = 1")
        other.execute("COMMIT")
        return True, time.monotonic() - started
    except sqlite3.OperationalError as e:
        assert "locked" in str(e)
        return False, time.monotonic() - started
    finally:
        time.sleep(1.2)
        holder.close()
        other.close()


@pytest.fixture
def db_path(tmp_path):
    path = str(tmp_path / "db.sqlite3")
    c = sqlite3.connect(path, isolation_level=None)
    c.execute("CREATE TABLE job (id INTEGER, status INTEGER)")
    c.execute("INSERT INTO job VALUES (1, 1)")
    c.close()
    return path


def test_a_deferred_transaction_fails_at_once(db_path):
    """The failure this guards against: no wait, whatever the busy_timeout."""
    ok, waited = _contended_write(db_path, "BEGIN")
    assert not ok
    assert waited < 0.5


def test_an_immediate_transaction_waits_its_turn(db_path):
    ok, waited = _contended_write(db_path, "BEGIN IMMEDIATE")
    assert ok
    assert 0.8 < waited < 10
