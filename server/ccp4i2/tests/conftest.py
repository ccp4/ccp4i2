"""Refuse to run the suite against anything but the isolated test settings.

On 2026-09-04 the API unit suite was run from a shell that had exported
``DJANGO_SETTINGS_MODULE=ccp4i2.config.settings`` for i2run. pytest-django
ranks that environment variable above the ``DJANGO_SETTINGS_MODULE`` line in
pytest.ini, and the ``os.environ[...]`` override at the top of the API conftest
comes too late, Django being configured by then. The suite therefore ran on the
production settings, and its transactional teardown flushed the developer's
live ``~/.ccp4i2-django/db.sqlite3``: every project and job gone, while the
desktop app was running on it.

``--ds`` in pytest.ini's addopts now outranks the variable, and this hook is
the second lock on the same door: whatever the settings module turns out to
be, no test runs while the database it points at is a file the user relies on.
"""
from pathlib import Path

import pytest

TEST_SETTINGS = "ccp4i2.config.test_settings"


def _real_homes():
    home = Path.home().resolve()
    return {home / ".ccp4i2x", home / ".ccp4i2-django", home / ".ccp4i2"}


def _is_under(path: Path, roots) -> bool:
    return any(path == root or root in path.parents for root in roots)


def pytest_sessionstart(session):
    from django.conf import settings

    module = getattr(settings, "SETTINGS_MODULE", None)
    if module != TEST_SETTINGS:
        raise pytest.UsageError(
            f"Refusing to run tests with DJANGO_SETTINGS_MODULE={module!r}. "
            f"The suite must run on {TEST_SETTINGS}; a run on the production "
            "settings flushes the live user database. Unset the variable in "
            "this shell (it is exported for i2run) and run again."
        )

    name = str(settings.DATABASES["default"]["NAME"])
    if name and name != ":memory:" and not name.startswith("file:"):
        if _is_under(Path(name).expanduser().resolve(), _real_homes()):
            raise pytest.UsageError(
                f"Refusing to run tests against {name}: that is the user's "
                "live CCP4i2 database, and the suite flushes whatever it runs on."
            )
