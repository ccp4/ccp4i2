# Coot 1.x startup stub (Python 3).
#
# Passed to Coot as `--script` by the coot1 task wrapper. It carries no
# job data: everything comes from the environment handshake (see
# handshake.py / api_client.py). Runs in Coot's __main__.
#
# This is a real, static file on purpose - the wrapper used to generate
# it as an interpolated string, which was fragile and unlintable.

import os
import sys
import traceback

try:
    try:
        from ccp4i2.cootbridge import coot1_gui
    except ImportError:
        # Non-CCP4 Coot builds may not have ccp4i2 importable; fall back
        # to the package location the handshake published.
        bridge_dir = os.environ.get("CCP4I2_COOTBRIDGE_DIR")
        if bridge_dir:
            server_dir = os.path.dirname(os.path.dirname(bridge_dir))
            if server_dir not in sys.path:
                sys.path.insert(0, server_dir)
        from ccp4i2.cootbridge import coot1_gui
    coot1_gui.start()
except Exception:
    print("[ccp4i2-cootbridge] startup failed:")
    traceback.print_exc()
