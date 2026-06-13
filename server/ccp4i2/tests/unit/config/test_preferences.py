"""
Tests for the shared CCP4i2 preferences mechanism
(``ccp4i2/config/preferences.py``) and its integration into Django settings.

Covers:
  * the pure-stdlib preferences module (load/save/resolve, home override,
    robustness to missing/corrupt files);
  * end-to-end settings precedence (env var > preferences.json > default),
    verified by importing settings.py in a fresh subprocess so each case picks
    up its own environment and preferences file.

All of this runs on a slim, CCP4-free interpreter.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ccp4i2.config import preferences

SERVER_ROOT = str(Path(__file__).resolve().parents[4])  # .../server


# ---------------------------------------------------------------------------
# preferences module (pure stdlib)
# ---------------------------------------------------------------------------

def test_home_override(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    assert preferences.ccp4i2_home() == tmp_path.resolve()
    assert preferences.preferences_path() == tmp_path.resolve() / "preferences.json"


def test_load_missing_returns_empty(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    assert preferences.load_preferences() == {}


def test_corrupt_file_returns_empty(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    preferences.preferences_path().write_text("{ this is not json", encoding="utf-8")
    assert preferences.load_preferences() == {}


def test_save_load_roundtrip(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    preferences.save_preferences({"ccp4Dir": "/opt/ccp4", "projectsDir": "/data/proj"})
    loaded = preferences.load_preferences()
    assert loaded["ccp4Dir"] == "/opt/ccp4"
    assert loaded["projectsDir"] == "/data/proj"
    # version is stamped, and the file is valid JSON on disk
    assert loaded["version"] == preferences.PREFERENCES_VERSION
    json.loads(preferences.preferences_path().read_text(encoding="utf-8"))


def test_save_is_atomic_no_temp_left(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    preferences.save_preferences({"ccp4Dir": "/opt/ccp4"})
    leftovers = list(tmp_path.glob(".preferences-*.tmp"))
    assert leftovers == []


def test_resolve_precedence(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    prefs = {"projectsDir": "/from/file"}

    # default only
    monkeypatch.delenv("CCP4I2_PROJECTS_DIR", raising=False)
    assert preferences.resolve("missing", default="/the/default", prefs=prefs) == "/the/default"
    # file beats default
    assert preferences.resolve("projectsDir", default="/the/default", prefs=prefs) == "/from/file"
    # env beats file
    monkeypatch.setenv("CCP4I2_PROJECTS_DIR", "/from/env")
    assert (
        preferences.resolve("projectsDir", env="CCP4I2_PROJECTS_DIR",
                             default="/the/default", prefs=prefs)
        == "/from/env"
    )


def test_resolve_ignores_empty_env(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.setenv("CCP4I2_PROJECTS_DIR", "")  # empty string must not win
    prefs = {"projectsDir": "/from/file"}
    assert (
        preferences.resolve("projectsDir", env="CCP4I2_PROJECTS_DIR", prefs=prefs)
        == "/from/file"
    )


# ---------------------------------------------------------------------------
# settings.py precedence (end-to-end, fresh subprocess per case)
# ---------------------------------------------------------------------------

_PROBE = (
    "import json; from ccp4i2.config import settings as s; "
    "print(json.dumps({"
    "'projects': str(s.CCP4I2_PROJECTS_DIR), "
    "'root': str(s.CCP4I2_ROOT), "
    "'engine': s.DATABASES['default']['ENGINE'], "
    "'name': str(s.DATABASES['default']['NAME']), "
    "'host': str(s.DATABASES['default'].get('HOST', '')), "
    "}))"
)


def _probe_settings(home: Path, prefs: dict, extra_env: dict = None) -> dict:
    if prefs is not None:
        home.mkdir(parents=True, exist_ok=True)
        (home / "preferences.json").write_text(json.dumps(prefs), encoding="utf-8")
    env = dict(os.environ)
    env["CCP4I2_HOME"] = str(home)
    env["PYTHONPATH"] = SERVER_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    # isolate from any ambient overrides
    for var in ("DATABASE_URL", "CCP4I2_DB_FILE", "CCP4I2_PROJECTS_DIR", "CCP4I2_ROOT"):
        env.pop(var, None)
    env.update(extra_env or {})
    result = subprocess.run(
        [sys.executable, "-c", _PROBE], env=env, capture_output=True, text=True, timeout=120
    )
    assert result.returncode == 0, f"probe failed:\n{result.stdout}\n{result.stderr}"
    return json.loads(result.stdout.strip().splitlines()[-1])


def test_settings_reads_preferences(tmp_path):
    home = tmp_path / "home"
    out = _probe_settings(home, {
        "projectsDir": str(tmp_path / "myprojects"),
        "ccp4i2Root": str(tmp_path / "myroot"),
        "database": "postgres://u:p@dbhost:5432/mydb",
    })
    assert out["projects"] == str(tmp_path / "myprojects")
    assert out["root"] == str(tmp_path / "myroot")
    assert out["engine"] == "django.db.backends.postgresql"
    assert out["name"] == "mydb"
    assert out["host"] == "dbhost"


def test_settings_default_when_no_preferences(tmp_path):
    home = tmp_path / "home"
    out = _probe_settings(home, {})  # empty preferences.json
    assert out["projects"] == str(home.resolve() / "CCP4X_PROJECTS")
    assert out["engine"] == "django.db.backends.sqlite3"


def test_settings_env_overrides_preferences(tmp_path):
    home = tmp_path / "home"
    out = _probe_settings(
        home,
        {"projectsDir": str(tmp_path / "fromfile")},
        extra_env={"CCP4I2_PROJECTS_DIR": str(tmp_path / "fromenv")},
    )
    assert out["projects"] == str(tmp_path / "fromenv")


# ---------------------------------------------------------------------------
# functional preferences (the userPreferences bag) and the PREFERENCES() accessor
# ---------------------------------------------------------------------------

def test_user_preference_unset_returns_default(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("SHELXDIR", raising=False)
    assert preferences.user_preference("SHELXDIR") is None
    assert preferences.user_preference("SHELXDIR", default="/fallback") == "/fallback"


def test_user_preference_from_file(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("SHELXDIR", raising=False)
    preferences.save_preferences({"userPreferences": {"SHELXDIR": "/opt/shelx"}})
    assert preferences.user_preference("SHELXDIR") == "/opt/shelx"


def test_user_preference_env_overrides_file(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    preferences.save_preferences({"userPreferences": {"SHELXDIR": "/from/file"}})
    monkeypatch.setenv("SHELXDIR", "/from/env")
    assert preferences.user_preference("SHELXDIR") == "/from/env"


def test_user_preference_preserves_json_bool(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("RETAIN_DIAGNOSTIC_FILES", raising=False)
    preferences.save_preferences({"userPreferences": {"RETAIN_DIAGNOSTIC_FILES": False}})
    # a real JSON bool comes back as a bool (not the truthy string "False")
    assert preferences.user_preference("RETAIN_DIAGNOSTIC_FILES") is False


def test_PREFERENCES_accessor(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    monkeypatch.delenv("SHELXDIR", raising=False)
    monkeypatch.delenv("BUSTERDIR", raising=False)
    preferences.save_preferences({"userPreferences": {"SHELXDIR": "/opt/shelx"}})

    from ccp4i2.core import CCP4Modules
    prefs = CCP4Modules.PREFERENCES()
    assert prefs.SHELXDIR == "/opt/shelx"      # from file
    assert prefs.BUSTERDIR is None             # unset -> None, never raises
