"""
Pytest configuration for cdata-codegen tests.

This file is automatically loaded by pytest and contains test fixtures
and configuration that applies to all tests.
"""

import sys
import os
from pathlib import Path

# Add ccp4i2 to Python path if CCP4I2_ROOT points to it
# IMPORTANT: We add it at the END of sys.path so cdata-codegen modules take precedence
ccp4i2_root = os.environ.get("CCP4I2_ROOT")
if ccp4i2_root:
    ccp4i2_path = Path(ccp4i2_root)
    # Check if this is actually the ccp4i2 directory (has wrappers/)
    if (ccp4i2_path / "wrappers").exists():
        # This is the actual ccp4i2 directory
        if str(ccp4i2_path) not in sys.path:
            sys.path.append(str(ccp4i2_path))  # Use append, not insert(0)
        print(f"Added ccp4i2 to Python path (at end): {ccp4i2_path}")
    else:
        # CCP4I2_ROOT is set to cdata-codegen directory, try to find actual ccp4i2
        parent_dir = ccp4i2_path.parent
        actual_ccp4i2 = parent_dir / "ccp4i2"
        if actual_ccp4i2.exists() and (actual_ccp4i2 / "wrappers").exists():
            if str(actual_ccp4i2) not in sys.path:
                sys.path.append(str(actual_ccp4i2))  # Use append, not insert(0)
            print(f"Found and added ccp4i2 to Python path (at end): {actual_ccp4i2}")
