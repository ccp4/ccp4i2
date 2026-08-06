"""
Tests for program (binary) discovery and legacy-preference migration
(``ccp4i2/config/program_discovery.py`` + ``preferences_migration.py``).

Slim, CCP4-free: resolution is pure path logic; migration parses an XML file we
construct in a temp dir. Each test isolates CCP4I2_HOME via monkeypatch.
"""

import json
import os
import stat
from pathlib import Path

import pytest

from ccp4i2.config import preferences, program_discovery
from ccp4i2.config import preferences_migration


def _write_prefs(home: Path, user_prefs: dict):
    (home / "preferences.json").write_text(
        json.dumps({"version": 1, "userPreferences": user_prefs})
    )


def _make_exe(path: Path):
    path.write_text("#!/bin/sh\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_resolve_via_suite_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    sdir = tmp_path / "shelxbin"
    sdir.mkdir()
    _make_exe(sdir / "shelxd")
    _write_prefs(tmp_path, {"SHELXDIR": str(sdir)})

    result = program_discovery.discover_program("shelxd")
    assert result["path"] == str(sdir / "shelxd")
    assert result["source"] == program_discovery.SOURCE_SUITE_DIR


def test_resolve_via_explicit_executable(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    coot = tmp_path / "mycoot"
    _make_exe(coot)
    _write_prefs(tmp_path, {"COOT_EXECUTABLE": str(coot)})

    result = program_discovery.discover_program("coot")
    assert result["path"] == str(coot)
    assert result["source"] == program_discovery.SOURCE_EXECUTABLE_PREF


def test_resolve_via_exe_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    d = tmp_path / "extra"
    d.mkdir()
    _make_exe(d / "someprog")
    _write_prefs(tmp_path, {"exePaths": [str(d)]})

    result = program_discovery.discover_program("someprog")
    assert result["path"] == str(d / "someprog")
    assert result["source"] == program_discovery.SOURCE_EXE_PATHS


def test_resolve_via_path(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    # ls is on PATH on any POSIX CI host.
    result = program_discovery.discover_program("ls")
    assert result["path"] is not None
    assert result["source"] == program_discovery.SOURCE_PATH


def test_resolve_missing(monkeypatch, tmp_path):
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    result = program_discovery.discover_program("zzz_definitely_not_a_program")
    assert result["path"] is None
    assert result["source"] == program_discovery.SOURCE_MISSING


def test_precedence_explicit_beats_path(monkeypatch, tmp_path):
    # An explicit COOT_EXECUTABLE must win over a PATH coot.
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    explicit = tmp_path / "explicit_coot"
    _make_exe(explicit)
    _write_prefs(tmp_path, {"COOT_EXECUTABLE": str(explicit)})
    # Even if a 'coot' were on PATH, the explicit pref is checked first.
    result = program_discovery.discover_program("coot")
    assert result["source"] == program_discovery.SOURCE_EXECUTABLE_PREF


def test_env_var_overrides_file(monkeypatch, tmp_path):
    # user_preference precedence is env > file; a SHELXDIR env var must win.
    monkeypatch.setenv("CCP4I2_HOME", str(tmp_path))
    filedir = tmp_path / "filedir"
    filedir.mkdir()
    _make_exe(filedir / "shelxd")
    envdir = tmp_path / "envdir"
    envdir.mkdir()
    _make_exe(envdir / "shelxd")
    _write_prefs(tmp_path, {"SHELXDIR": str(filedir)})
    monkeypatch.setenv("SHELXDIR", str(envdir))

    result = program_discovery.discover_program("shelxd")
    assert result["path"] == str(envdir / "shelxd")


# --- migration ---------------------------------------------------------------


LEGACY_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header><pluginName>guipreferences</pluginName></ccp4i2_header>
  <ccp4i2_body>
    <SHELXDIR>/opt/shelx</SHELXDIR>
    <COOT_EXECUTABLE>/opt/coot/bin/coot</COOT_EXECUTABLE>
    <EXEPATHLIST>
      <CExePath>
        <exeName>moorhen</exeName>
        <exePath>
          <baseName>Moorhen.app</baseName>
          <relPath>/opt/moorhen</relPath>
        </exePath>
      </CExePath>
    </EXEPATHLIST>
  </ccp4i2_body>
</ccp4:ccp4i2>
"""


def test_parse_legacy_program_prefs(tmp_path):
    legacy = tmp_path / "guipreferences.params.xml"
    legacy.write_text(LEGACY_XML)
    parsed = preferences_migration.parse_legacy_program_prefs(legacy)
    assert parsed["SHELXDIR"] == "/opt/shelx"
    assert parsed["COOT_EXECUTABLE"] == "/opt/coot/bin/coot"
    assert parsed["exePaths"] == ["/opt/moorhen"]


def test_migration_desktop_imports_and_is_idempotent(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    legacy = tmp_path / "guipreferences.params.xml"
    legacy.write_text(LEGACY_XML)
    monkeypatch.setenv("CCP4I2_HOME", str(home))
    monkeypatch.setenv("CCP4I2_LEGACY_PREFS", str(legacy))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")  # desktop

    imported = preferences_migration.migrate_legacy_program_prefs()
    assert imported.get("SHELXDIR") == "/opt/shelx"
    bag = preferences.load_preferences()["userPreferences"]
    assert bag["SHELXDIR"] == "/opt/shelx"
    assert bag["exePaths"] == ["/opt/moorhen"]

    # Second run is a no-op (marker set).
    assert preferences_migration.migrate_legacy_program_prefs() == {}


def test_migration_cloud_noop(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    legacy = tmp_path / "guipreferences.params.xml"
    legacy.write_text(LEGACY_XML)
    monkeypatch.setenv("CCP4I2_HOME", str(home))
    monkeypatch.setenv("CCP4I2_LEGACY_PREFS", str(legacy))
    monkeypatch.delenv("CCP4I2_LOCAL_SESSION_TOKEN", raising=False)  # cloud

    assert preferences_migration.migrate_legacy_program_prefs() == {}
    assert not (home / "preferences.json").exists()


def test_migration_non_destructive(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    legacy = tmp_path / "guipreferences.params.xml"
    legacy.write_text(LEGACY_XML)
    monkeypatch.setenv("CCP4I2_HOME", str(home))
    monkeypatch.setenv("CCP4I2_LEGACY_PREFS", str(legacy))
    monkeypatch.setenv("CCP4I2_LOCAL_SESSION_TOKEN", "desktop")
    # A Django-set SHELXDIR must not be clobbered by the legacy value.
    _write_prefs(home, {"SHELXDIR": "/my/django/shelx"})

    preferences_migration.migrate_legacy_program_prefs()
    bag = preferences.load_preferences()["userPreferences"]
    assert bag["SHELXDIR"] == "/my/django/shelx"  # existing value wins
