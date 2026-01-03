"""
CCP4 DLL path initialization for Windows.
Auto-imported to set up DLL directories before other modules load.
"""

import os

if os.name == "nt" and hasattr(os, "add_dll_directory"):
    ccp4_dir = os.environ.get("CCP4")
    if ccp4_dir:
        try:
            os.add_dll_directory(os.path.join(ccp4_dir, "bin"))
        except (OSError, AttributeError):
            pass
