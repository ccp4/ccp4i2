"""
Shared CCP4i2 user preferences — the single source of truth for desktop config
that must agree across the GUI, the ``i2``/``i2run`` CLI, and the Django control
plane.

Location
--------
``~/.ccp4i2x/preferences.json`` (override the home with the ``CCP4I2_HOME``
environment variable). This sits next to the SQLite database, the project store,
credentials and the project registry, so ONE directory holds all per-user state.
``~/.ccp4i2x`` and not ``~/.ccp4i2``: the latter collides with the legacy Qt home
``~/.CCP4I2`` on case-insensitive filesystems (see ``ccp4i2_home``). An install
predating the rename keeps its ``~/.ccp4i2-django`` home, adopted in place. JSON
(vs the classic XML prefs) so the Electron/JS side reads it as easily as Python.

Resolution precedence (for every setting)
-----------------------------------------
``environment variable  >  preferences.json  >  built-in default``

This keeps cloud deployments (which drive everything from env vars) completely
unaffected — the file is purely the *desktop* persistence layer — while giving the
GUI and CLI one place to agree.

Schema (``preferences.json``)
-----------------------------
All keys are optional; an absent key falls through to the env var then the default.
::

    {
      "version":     1,                                   # schema version
      "ccp4Dir":     "/path/to/ccp4",                     # CCP4 install (for ccp4-python + binaries)
      "projectsDir": "/path/to/projects",                 # project store      (env CCP4I2_PROJECTS_DIR)
      "database":    "sqlite:////abs/path/db.sqlite3",    # DB URL, sqlite or postgres (env DATABASE_URL)
      "ccp4i2Root":  "/path/to/ccp4i2/server/ccp4i2",     # package root       (env CCP4I2_ROOT)
      "userPreferences": { }                              # open bag for future app/scientific prefs
    }

Bootstrap note: ``ccp4Dir``/``projectsDir``/``database`` must live in a *file*, not
the database — you need them to find the database in the first place.

This module is pure stdlib (it is imported while Django settings are still being
assembled) — it must not import Django or any ccp4i2 package.
"""

import json
import os
import tempfile
from pathlib import Path

PREFERENCES_VERSION = 1

# Keys that hold filesystem paths (kept here so callers/tests can reason about them).
PATH_KEYS = ("ccp4Dir", "projectsDir", "ccp4i2Root")


# The user home, and the one it superseded. ``~/.ccp4i2x`` rather than
# ``~/.ccp4i2``: the latter differs from the legacy Qt-era home ``~/.CCP4I2``
# only by case, and on a case-insensitive filesystem (the macOS default) they
# are the SAME directory — the new app would write its db and preferences into
# the tree it is migrating FROM. ``x`` rather than ``-django`` because the home
# should not name an implementation detail the user never chose.
HOME_DIR_NAME = ".ccp4i2x"
LEGACY_HOME_DIR_NAME = ".ccp4i2-django"


def ccp4i2_home() -> Path:
    """The CCP4i2 user home: everything per-user lives under one root.

    Resolution, in order:

    1. ``CCP4I2_HOME`` — wins outright (cloud, tests, deliberate relocation).
    2. An existing ``~/.ccp4i2x``.
    3. An existing ``~/.ccp4i2-django`` — **adopted in place.** Installs from
       3.1.0a26 and earlier keep working exactly where they are; nothing is
       moved, so no project directory is ever relocated behind the user's back.
    4. ``~/.ccp4i2x`` — a fresh install.

    Kept in sync with the Electron side (client/main/ccp4i2-preferences.ts);
    the two MUST agree or the desktop app and the server disagree about where
    the database and projects live, which is exactly the drift this replaced.
    """
    override = os.environ.get("CCP4I2_HOME")
    if override:
        return Path(override).expanduser().resolve()

    home = Path.home()
    current = home / HOME_DIR_NAME
    if current.is_dir():
        return current.expanduser().resolve()

    legacy = home / LEGACY_HOME_DIR_NAME
    if legacy.is_dir():
        return legacy.expanduser().resolve()

    return current.expanduser().resolve()


def default_projects_dir() -> Path:
    """Where new projects go unless the user says otherwise: ``<home>/projects``.

    An adopted pre-a27 home keeps its existing ``CCP4X_PROJECTS`` directory —
    renaming it would strand every absolute path recorded in those projects.
    """
    home = ccp4i2_home()
    legacy = home / "CCP4X_PROJECTS"
    if legacy.is_dir() and not (home / "projects").is_dir():
        return legacy
    return home / "projects"


def preferences_path() -> Path:
    """Full path to ``preferences.json`` inside the CCP4i2 user home."""
    return ccp4i2_home() / "preferences.json"


def is_desktop() -> bool:
    """True on a desktop (Electron) launch, False in a cloud/server deployment.

    ``preferences.json`` and legacy-preference migration are a *desktop* concept:
    in cloud, settings arrive as env vars (see the module docstring) and the
    per-user file/home is meaningless (ephemeral, per-replica). The desktop
    launcher sets ``CCP4I2_LOCAL_SESSION_TOKEN`` (the same signal
    ``settings.py`` uses to pick the local-session auth middleware), so use it
    here to gate file writes and legacy migration.
    """
    return bool(os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN"))


def load_preferences() -> dict:
    """Load ``preferences.json``; return ``{}`` if absent or unreadable.

    Never raises — a missing or corrupt file must not break settings load or the
    CLI; callers fall through to env vars / defaults.
    """
    path = preferences_path()
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (FileNotFoundError, NotADirectoryError, json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def save_preferences(prefs: dict) -> Path:
    """Atomically write ``preferences.json`` (temp file + ``os.replace``).

    Stamps the current schema version. Last writer wins — adequate for a small,
    rarely-contended config file shared between the GUI and the CLI.
    """
    home = ccp4i2_home()
    home.mkdir(parents=True, exist_ok=True)
    payload = {"version": PREFERENCES_VERSION, **prefs}

    fd, tmp = tempfile.mkstemp(dir=str(home), prefix=".preferences-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp, preferences_path())
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return preferences_path()


def resolve(key: str, env: str = None, default=None, prefs: dict = None):
    """Resolve one setting as ``env var > preferences.json > default``.

    Args:
        key:     key in ``preferences.json``.
        env:     environment variable that overrides it (optional).
        default: value if neither env nor file provides one.
        prefs:   pre-loaded preferences dict (avoids re-reading the file when
                 resolving several settings in a row, e.g. from settings.py).
    """
    if env:
        env_val = os.environ.get(env)
        if env_val not in (None, ""):
            return env_val
    if prefs is None:
        prefs = load_preferences()
    file_val = prefs.get(key)
    if file_val not in (None, ""):
        return file_val
    return default


def user_preference(name: str, default=None):
    """Resolve a *functional* preference (from the ``userPreferences`` bag).

    These are the preferences consumed by wrappers/pipelines via the
    ``PREFERENCES()`` accessor — e.g. ``SHELXDIR``, ``BUSTERDIR``,
    ``COOT_EXECUTABLE``, ``PDB_REDO_TOKEN_ID``, ``RETAIN_DIAGNOSTIC_FILES``.

    Precedence is the standard ``env var > preferences.json > default``; the
    environment variable name is the preference name itself, so in cloud these
    arrive as plain env vars / Key Vault secrets and on desktop from the file's
    ``userPreferences`` object.
    """
    bag = load_preferences().get("userPreferences", {})
    bag = bag if isinstance(bag, dict) else {}
    return resolve(name, env=name, default=default, prefs=bag)
