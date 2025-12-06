"""
Pytest configuration for ccp4i2 tests.

This file is automatically loaded by pytest and contains test fixtures
and configuration that applies to all tests.
"""

import sys
import os
from pathlib import Path

# Ensure CCP4I2_ROOT is set to project root if not already defined
if not os.environ.get("CCP4I2_ROOT"):
    os.environ["CCP4I2_ROOT"] = str(Path(__file__).resolve().parent.parent)

# Add ccp4i2 to Python path if needed
ccp4i2_root = os.environ.get("CCP4I2_ROOT")
if ccp4i2_root:
    ccp4i2_path = Path(ccp4i2_root)
    # Check if this is actually the ccp4i2 directory (has wrappers/)
    if (ccp4i2_path / "wrappers").exists():
        # This is the actual ccp4i2 directory
        if str(ccp4i2_path) not in sys.path:
            sys.path.append(str(ccp4i2_path))
        print(f"Added ccp4i2 to Python path: {ccp4i2_path}")
