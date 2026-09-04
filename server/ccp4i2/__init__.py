"""
CCP4i2 - Graphical interface and scripting environment for CCP4

This package provides the core functionality for CCP4i2, including:
- Task wrappers for crystallographic programs
- Pipelines that chain multiple programs
- Data management utilities
- Report generation
"""

from datetime import datetime
from pathlib import Path

MAJOR = 3
MINOR = 1
PATCH = 0
# PEP 440 pre-release segment: "a1" (alpha), "b2" (beta), "rc1", or "" for a
# final release. The alpha line is deliberately a pre-release so plain
# `pip install ccp4i2` never selects it, and the desktop app pins the EXACT
# version (see client/main/ccp4i2-server-version.ts) — an alpha app and its
# backend are strictly bound and cannot mix with the un-suffixed 3.0.x wheels
# published before this discipline.
PRERELEASE = "a45"

__version__ = f"{MAJOR}.{MINOR}.{PATCH}{PRERELEASE}"
__version_date__ = datetime(2026, 9, 4)
__version_info__ = (MAJOR, MINOR, PATCH)

I2_TOP = Path(__file__).resolve().parent
