"""The SQLite connection settings CCP4i2 runs on.

Kept apart from settings.py so that the choice is testable without importing
the settings module, and so that the reason for it is written once.
"""

SQLITE_INIT_COMMAND = (
    "PRAGMA journal_mode=WAL;"
    "PRAGMA synchronous=NORMAL;"
    "PRAGMA busy_timeout=30000;"
)


def sqlite_database(name):
    """The SQLite connection settings, used by both SQLite branches of settings.py.

    Several processes write this database at once: two uvicorn workers, and
    a ccp4-python subprocess per running job that records status and gleans
    files. Django's default transaction is deferred -- it takes a read lock
    on the first SELECT and upgrades to a write lock on the first write. If
    another process already holds the write lock at that moment, SQLite does
    not wait and does not call the busy handler; it fails at once with
    "database is locked". That is what killed jobs and create_task calls
    during a New Project import, half a second apart.

    - transaction_mode IMMEDIATE takes the write lock at BEGIN, so a
      contended transaction waits for its turn instead of failing.
    - busy_timeout / timeout give it 30 s to wait.
    - WAL lets readers and the one writer proceed together, and persists in
      the file once set.
    """
    return {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": name,
        "OPTIONS": {
            "timeout": 30,
            "transaction_mode": "IMMEDIATE",
            "init_command": SQLITE_INIT_COMMAND,
        },
    }
