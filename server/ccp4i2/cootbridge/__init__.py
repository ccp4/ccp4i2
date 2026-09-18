# The CCP4i2 <-> Coot bridge.
#
# Layout (see docs in each module):
#
#   api_client.py    - shared data-interaction layer. Python 2.7 AND 3.x,
#                      stdlib only, no GUI imports, no other ccp4i2 imports.
#                      Coot 0.9 stubs load this file directly by path
#                      (imp.load_source), so it must stay standalone.
#   coot1_gui.py     - Coot 1.x adapter: initial data load, the "CCP4i2"
#                      menu (save to job), and the project-hierarchy
#                      browser widget. Python 3 / GTK4 (PyGObject).
#   coot09_loader.py - Coot 0.9 adapter: initial data load through the
#                      flat scripting namespace. Python 2.7 compatible,
#                      no GUI (the GTK2 browser is a later tier).
#
# The module loaded into Coot receives ONLY connection details and the
# job identity, via environment variables set by the coot1 task wrapper
# (wrappers/coot1/script/coot1.py). Everything else - what to load, where
# files live, where to save - is fetched from the CCP4i2 REST API (or,
# offline, from the job directory the identity points at).
#
# This package intentionally has an empty-ish __init__: it must be
# importable without Django, CCP4, or GTK present.

import os as _os

_HERE = _os.path.dirname(_os.path.abspath(__file__))

#: The static scripts a wrapper passes to Coot as ``--script`` (real
#: files, not generated source).
COOT1_STARTUP_STUB = _os.path.join(_HERE, "_coot1_stub.py")
COOT09_STARTUP_STUB = _os.path.join(_HERE, "_coot09_stub.py")
