"""
Program (binary) discovery — resolve a task's executable against user
preferences before falling back to ``PATH``.

Django's task execution was a bare ``shutil.which(TASKCOMMAND)``: if a program
(coot, ccp4mg, shelx*, dials, buster, …) was not on ``PATH``, the job failed
with no supported way to point CCP4i2 at it. The legacy Qt GUI let users set
program locations (``COOT_EXECUTABLE``, ``SHELXDIR``, ``EXEPATHLIST`` …); this
module is the Django equivalent.

Resolution order for a program ``name``:

1. an explicit ``{PROG}_EXECUTABLE`` preference (exact executable path) for the
   programs that have one (coot, ccp4mg);
2. a ``{SUITE}DIR`` preference joined with ``name`` (shelx* -> SHELXDIR,
   dials -> DIALSDIR, buster/refine -> BUSTERDIR);
3. each directory in the general ``exePaths`` list (the django name for the
   legacy ``EXEPATHLIST``);
4. ``shutil.which(name)`` — the ``PATH`` default (unchanged behaviour);
5. ``None`` — caller raises an actionable error.

Preferences are read via ``config.preferences.user_preference`` (env var >
``preferences.json`` ``userPreferences`` bag > default), so cloud deployments
can supply them as plain env vars. Pure stdlib — no Django, no CCP4 import — so
it is safe to call from the CCP4-free server (the probe endpoint) and from the
ccp4-python job environment (``CCP4PluginScript``) alike.
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional

from ccp4i2.config.preferences import user_preference

# Programs with an explicit single-executable preference.
_EXECUTABLE_PREF: Dict[str, str] = {
    "coot": "COOT_EXECUTABLE",
    "ccp4mg": "CCP4MG_EXECUTABLE",
}

# Programs whose executable lives in a preference-specified suite directory.
# Maps program name -> the {SUITE}DIR preference key that holds its directory.
_SUITE_DIR_PREF: Dict[str, str] = {
    "shelxc": "SHELXDIR",
    "shelxd": "SHELXDIR",
    "shelxe": "SHELXDIR",
    "shelxl": "SHELXDIR",
    "shelx": "SHELXDIR",
    "dials": "DIALSDIR",
    "dials.import": "DIALSDIR",
    "dials.integrate": "DIALSDIR",
    "buster": "BUSTERDIR",
    "refine": "BUSTERDIR",
}

# The general extra-search-path list preference (legacy EXEPATHLIST).
_EXE_PATHS_KEY = "exePaths"

# Discovery source labels (also surfaced by the probe endpoint).
SOURCE_EXECUTABLE_PREF = "executable_pref"
SOURCE_SUITE_DIR = "suite_dir"
SOURCE_EXE_PATHS = "exe_paths"
SOURCE_PATH = "path"
SOURCE_MISSING = "missing"

# Programs the discovery UI probes by default. Extendable via preferences later.
KNOWN_PROGRAMS: List[str] = [
    "coot", "ccp4mg", "shelxc", "shelxd", "shelxe", "dials", "buster",
    "refmac5", "refmacat", "servalcat", "aimless", "pointless", "ctruncate",
    "freerflag", "phaser", "molrep", "modelcraft", "nautilus", "buccaneer",
]


def _exe_paths() -> List[str]:
    """The user's extra executable search directories (``exePaths`` pref)."""
    value = user_preference(_EXE_PATHS_KEY, default=None)
    if not value:
        return []
    if isinstance(value, str):
        # tolerate an os.pathsep-joined string as well as a JSON list
        return [p for p in value.split(os.pathsep) if p]
    if isinstance(value, (list, tuple)):
        return [str(p) for p in value if p]
    return []


def _is_executable_file(path: Path) -> bool:
    try:
        return path.is_file() and os.access(str(path), os.X_OK)
    except OSError:
        return False


def _candidate_names(name: str) -> List[str]:
    """Executable filenames to try for ``name`` (adds ``.exe`` on Windows)."""
    if os.name == "nt" and not name.lower().endswith(".exe"):
        return [name, f"{name}.exe"]
    return [name]


def resolve_program(name: str) -> Optional[str]:
    """Resolve ``name`` to an absolute executable path, or ``None``.

    Convenience wrapper around :func:`discover_program` returning just the path.
    """
    return discover_program(name)["path"]


def discover_program(name: str) -> Dict[str, Optional[str]]:
    """Resolve ``name`` and report *how* it was found.

    Returns ``{"name", "path", "source"}`` where ``source`` is one of the
    ``SOURCE_*`` constants (``missing`` when unresolved). Read-only: never runs
    the program.
    """
    # 1. explicit {PROG}_EXECUTABLE
    pref_key = _EXECUTABLE_PREF.get(name)
    if pref_key:
        explicit = user_preference(pref_key, default=None)
        if explicit:
            p = Path(explicit)
            if _is_executable_file(p):
                return {"name": name, "path": str(p), "source": SOURCE_EXECUTABLE_PREF}

    # 2. {SUITE}DIR joined with the program name
    dir_key = _SUITE_DIR_PREF.get(name)
    if dir_key:
        suite_dir = user_preference(dir_key, default=None)
        if suite_dir:
            for cand in _candidate_names(name):
                p = Path(suite_dir) / cand
                if _is_executable_file(p):
                    return {"name": name, "path": str(p), "source": SOURCE_SUITE_DIR}

    # 3. general exePaths list
    for d in _exe_paths():
        for cand in _candidate_names(name):
            p = Path(d) / cand
            if _is_executable_file(p):
                return {"name": name, "path": str(p), "source": SOURCE_EXE_PATHS}

    # 4. PATH (unchanged default)
    which = shutil.which(name)
    if which:
        return {"name": name, "path": which, "source": SOURCE_PATH}

    # 5. not found
    return {"name": name, "path": None, "source": SOURCE_MISSING}


def discover_programs(names: Optional[List[str]] = None) -> List[Dict[str, Optional[str]]]:
    """Discover a list of programs (default :data:`KNOWN_PROGRAMS`)."""
    return [discover_program(n) for n in (names if names is not None else KNOWN_PROGRAMS)]
