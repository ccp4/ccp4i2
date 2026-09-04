"""A test run must not write into the developer's real CCP4i2 home.

Isolating CCP4I2_PROJECTS_DIR is not enough, and the gap was not theoretical:
project_snapshot resolves its ``project_directories.json`` registry through
``preferences.ccp4i2_home()`` rather than through Django settings, so running
the API suite added an entry naming a ``~/.cache/ccp4i2-tests`` directory to
the real ``~/.ccp4i2x/project_directories.json``.

Anything else that reaches for the home directly -- credentials, preferences,
the server log -- had the same hole.
"""
from pathlib import Path

from django.conf import settings

from ccp4i2.config import preferences


def test_the_home_is_not_the_real_user_home():
    home = preferences.ccp4i2_home()
    for real in (".ccp4i2x", ".ccp4i2-django"):
        assert home != (Path.home().resolve() / real), (
            f"tests are resolving the CCP4i2 home to {home}, the real one"
        )


def test_the_home_is_the_isolated_test_directory():
    assert preferences.ccp4i2_home() == Path(settings.USER_DIR).resolve()


def test_the_registry_would_be_written_inside_it():
    """The specific file that leaked, addressed by the specific path used."""
    from ccp4i2.db.project_snapshot import registry_path

    assert registry_path().parent == preferences.ccp4i2_home()
    assert Path.home().resolve() != registry_path().parent


def test_the_settings_module_is_the_test_one():
    """The lock that would have held on 2026-09-04.

    An exported DJANGO_SETTINGS_MODULE=ccp4i2.config.settings (the i2run
    setting) outranks the pytest.ini line, so the API suite ran against the
    production settings and its transactional teardown flushed the live user
    database. ``--ds`` in addopts outranks the variable; this checks it did.
    """
    assert settings.SETTINGS_MODULE == "ccp4i2.config.test_settings"


def test_the_database_is_not_the_users():
    name = str(settings.DATABASES["default"]["NAME"])
    if name == ":memory:" or name.startswith("file:"):
        return
    resolved = Path(name).expanduser().resolve()
    for real in (".ccp4i2x", ".ccp4i2-django", ".ccp4i2"):
        root = Path.home().resolve() / real
        assert resolved != root and root not in resolved.parents, (
            f"tests would run against the live database {resolved}"
        )
