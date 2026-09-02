"""One user home, resolved the same way by the server and the desktop app.

Testers reported state scattered across ~/.ccp4i2/, ~/.ccp4i2-django/,
~/.config/ccp4i2-django/ and ~/.config/Electron/. Part of that was genuine
sprawl; part was a bug — the Electron app defaulted projects to
~/.ccp4i2/CCP4X_PROJECTS while Django settings said
~/.ccp4i2-django/CCP4X_PROJECTS, so the app and the server it launched
disagreed about where projects lived.

These tests pin the resolution order. The Electron side must match; it is
checked by reading the TypeScript, since there is no JS runtime here.
"""

import os
import re
from pathlib import Path

import pytest

from ccp4i2.config import preferences


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.delenv("CCP4I2_HOME", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    return tmp_path


def test_fresh_install_uses_the_new_home(home):
    assert preferences.ccp4i2_home() == (home / ".ccp4i2x").resolve()


def test_existing_install_is_adopted_in_place(home):
    """An a26-or-earlier install keeps working where it is. Nothing is moved:
    relocating a user's project directories behind their back is a far worse
    failure than an un-consolidated layout."""
    (home / ".ccp4i2-django").mkdir()
    assert preferences.ccp4i2_home() == (home / ".ccp4i2-django").resolve()


def test_new_home_wins_when_both_exist(home):
    (home / ".ccp4i2-django").mkdir()
    (home / ".ccp4i2x").mkdir()
    assert preferences.ccp4i2_home() == (home / ".ccp4i2x").resolve()


def test_env_override_beats_everything(home, monkeypatch):
    (home / ".ccp4i2x").mkdir()
    monkeypatch.setenv("CCP4I2_HOME", str(home / "elsewhere"))
    assert preferences.ccp4i2_home() == (home / "elsewhere").resolve()


def test_home_never_collides_with_the_legacy_qt_home():
    """~/.ccp4i2 would differ from the Qt-era ~/.CCP4I2 only by case, and on a
    case-insensitive filesystem they are the same directory — the new app would
    write its database into the tree it migrates from."""
    assert preferences.HOME_DIR_NAME.lower() != ".ccp4i2"


def test_default_projects_dir_is_under_the_home(home):
    assert preferences.default_projects_dir() == home / ".ccp4i2x" / "projects"


def test_adopted_home_keeps_its_existing_project_store(home):
    """Renaming CCP4X_PROJECTS would strand every absolute path recorded in the
    projects inside it."""
    (home / ".ccp4i2-django" / "CCP4X_PROJECTS").mkdir(parents=True)
    assert preferences.default_projects_dir().name == "CCP4X_PROJECTS"


def test_electron_side_resolves_the_same_way():
    """The two implementations must agree, or the desktop app and the server it
    launches use different databases. Checked by reading the source: a mismatch
    here is the bug this consolidation existed to fix."""
    ts = (
        Path(preferences.__file__).parents[3]
        / "client" / "main" / "ccp4i2-preferences.ts"
    )
    if not ts.is_file():          # server-only checkout
        pytest.skip("client/ not present")
    text = ts.read_text(encoding="utf-8")

    names = dict(re.findall(r'export const (\w+) = "([^"]+)";', text))
    assert names.get("HOME_DIR_NAME") == preferences.HOME_DIR_NAME
    assert names.get("LEGACY_HOME_DIR_NAME") == preferences.LEGACY_HOME_DIR_NAME
    # and it must consult the legacy home rather than jumping straight to the new one
    assert "LEGACY_HOME_DIR_NAME" in text.split("export function ccp4i2Home")[1]
