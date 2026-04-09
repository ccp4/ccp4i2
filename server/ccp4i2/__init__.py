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
MINOR = 0
PATCH = 0

__version__ = f"{MAJOR}.{MINOR}.{PATCH}"
__version_date__ = datetime(2025, 12, 10)
__version_info__ = (MAJOR, MINOR, PATCH)

I2_TOP = Path(__file__).resolve().parent
