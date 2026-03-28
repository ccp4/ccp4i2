# Copyright (C) 2025-2026 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
