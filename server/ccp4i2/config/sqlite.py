"""The SQLite connection settings CCP4i2 runs on.

Kept apart from settings.py so that the choice is testable without importing
the settings module, and so that the reason for it is written once.
"""

import django

SQLITE_INIT_COMMAND = "PRAGMA busy_timeout=30000;"

# Django 5.1 taught the sqlite3 backend the "transaction_mode" and
# "init_command" OPTIONS. Older backends pass every OPTIONS key straight to
# sqlite3.connect(), which rejects both with a TypeError on the first connect,
# so an install that resolved Django 4.2 (Materia pins django>=4.2,<5.0 in its
# server images) was dead on arrival with an error naming a sqlite keyword.
BACKEND_KNOWS_TRANSACTION_MODE = django.VERSION >= (5, 1)


def sqlite_options():
    """The OPTIONS dict for this Django, as its sqlite3 backend will accept.

    On 5.1+ the write lock is taken at BEGIN (IMMEDIATE) and busy_timeout is
    set by the init command. Before 5.1 neither key exists; the connection
    degrades to Django's deferred transaction, and the 30 s "timeout" -- which
    every version passes to sqlite3.connect() -- still gives a contended
    write that SQLite *does* wait on its chance. The contention fix is kept
    wherever it is supported; nowhere is the install broken by it.
    """
    options = {"timeout": 30}
    if BACKEND_KNOWS_TRANSACTION_MODE:
        options["transaction_mode"] = "IMMEDIATE"
        options["init_command"] = SQLITE_INIT_COMMAND
    return options


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

    Both OPTIONS keys are Django 5.1+; see sqlite_options() for what an older
    Django gets instead.

    The journal mode is deliberately left alone. WAL would let readers run
    beside the writer, but it needs shared memory between every process
    that opens the file and does not work over NFS or SMB -- and a user's
    home, where this database lives, is often a network share. The
    rollback journal works everywhere, and the two settings above are the
    ones that fix the failure.
    """
    return {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": name,
        "OPTIONS": sqlite_options(),
    }
