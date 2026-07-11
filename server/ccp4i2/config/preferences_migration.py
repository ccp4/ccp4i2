"""
One-time import of program-location preferences from the legacy Qt CCP4i2.

Established users have their program locations in the classic GUI's
``~/.CCP4I2/configs/guipreferences.params.xml`` (COOT_EXECUTABLE, SHELXDIR,
EXEPATHLIST …). The Django app keeps its own ``preferences.json`` and does not
share the legacy home (deliberately — see ``preferences.ccp4i2_home``). Without
migration those users would have to re-enter every program path.

This module reads the legacy XML and seeds the *program-location* keys of
``preferences.json`` (only the ones relevant to binary discovery — not the whole
legacy pref set). It is:

- **non-destructive**: never overwrites a key already present in
  ``preferences.json`` (the Django value wins);
- **idempotent**: a ``userPreferences._legacyProgramImport`` marker records that
  the import ran, so it does not re-import on every launch;
- **best-effort**: any parse error is swallowed — migration must never block
  startup.

Legacy serialisation handled:
- scalar ``<KEY>value</KEY>`` for COOT_EXECUTABLE / CCP4MG_EXECUTABLE /
  SHELXDIR / DIALSDIR / BUSTERDIR;
- ``<EXEPATHLIST><CExePath><exeName>…</exeName><exePath><baseName>…</baseName>
  <relPath>…</relPath></exePath></CExePath>…</EXEPATHLIST>`` — each entry's
  directory (relPath, or relPath/baseName) is collected into ``exePaths``.
"""

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional

from ccp4i2.config.preferences import (
    is_desktop,
    load_preferences,
    save_preferences,
)

# Scalar legacy keys copied verbatim into the userPreferences bag.
_SCALAR_KEYS = (
    "COOT_EXECUTABLE",
    "CCP4MG_EXECUTABLE",
    "SHELXDIR",
    "DIALSDIR",
    "BUSTERDIR",
)

_IMPORT_MARKER = "_legacyProgramImport"


def legacy_preferences_path() -> Path:
    """Classic Qt CCP4i2 GUI preferences file (``~/.CCP4I2/configs/…``)."""
    override = os.environ.get("CCP4I2_LEGACY_PREFS")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".CCP4I2" / "configs" / "guipreferences.params.xml"


def _text(elem: Optional[ET.Element]) -> Optional[str]:
    if elem is None or elem.text is None:
        return None
    value = elem.text.strip()
    return value or None


def parse_legacy_program_prefs(path: Path) -> Dict[str, object]:
    """Extract the program-location keys from a legacy prefs XML file.

    Returns a dict with any of the scalar keys plus an ``exePaths`` list
    (directories derived from EXEPATHLIST). Missing/empty values are omitted.
    Returns ``{}`` on any parse failure.
    """
    try:
        tree = ET.parse(str(path))
    except (ET.ParseError, OSError):
        return {}
    root = tree.getroot()
    body = root.find("ccp4i2_body")
    if body is None:
        # Some legacy files namespace or omit the wrapper; fall back to root.
        body = root

    out: Dict[str, object] = {}
    for key in _SCALAR_KEYS:
        val = _text(body.find(key))
        if val:
            out[key] = val

    exe_paths: List[str] = []
    exepathlist = body.find("EXEPATHLIST")
    if exepathlist is not None:
        for cexe in exepathlist.findall("CExePath"):
            exepath = cexe.find("exePath")
            if exepath is None:
                continue
            rel = _text(exepath.find("relPath"))
            base = _text(exepath.find("baseName"))
            if not rel:
                continue
            # relPath is the directory; baseName (when present) is the app/exe
            # inside it. The directory is what feeds exePaths.
            directory = rel
            if directory and directory not in exe_paths:
                exe_paths.append(directory)
    if exe_paths:
        out["exePaths"] = exe_paths

    return out


def migrate_legacy_program_prefs(force: bool = False) -> Dict[str, object]:
    """Seed ``preferences.json`` program-location keys from the legacy GUI.

    Non-destructive (existing keys win), idempotent (records a marker). Returns
    the dict of keys actually imported (empty if nothing to do).

    Args:
        force: re-run even if the import marker is present.
    """
    # Desktop-only: in cloud there is no legacy home, preferences arrive as env
    # vars, and preferences.json is per-replica/ephemeral. Never migrate there.
    if not is_desktop() and not force:
        return {}

    prefs = load_preferences()
    bag = prefs.get("userPreferences")
    bag = bag if isinstance(bag, dict) else {}

    if bag.get(_IMPORT_MARKER) and not force:
        return {}

    legacy = parse_legacy_program_prefs(legacy_preferences_path())

    imported: Dict[str, object] = {}
    for key, value in legacy.items():
        if key == "exePaths":
            existing = bag.get("exePaths")
            existing = list(existing) if isinstance(existing, list) else []
            merged = existing + [p for p in value if p not in existing]
            if merged != existing:
                bag["exePaths"] = merged
                imported["exePaths"] = value
        elif key not in bag:  # scalar: Django value wins if already set
            bag[key] = value
            imported[key] = value

    # Always stamp the marker so we don't re-scan the legacy file each launch.
    bag[_IMPORT_MARKER] = True
    prefs["userPreferences"] = bag
    save_preferences(prefs)
    return imported
