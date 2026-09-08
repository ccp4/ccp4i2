# -*- coding: utf-8 -*-
# Coot 0.9 startup stub (Python 2.7 / runs in Coot's __main__).
#
# Passed to Coot as `--script` by the coot_rebuild task wrapper. Carries
# no job data - the environment handshake does (see handshake.py). A
# real static file: the wrapper used to generate it as an interpolated
# string, which was fragile and unlintable.
#
# Coot 0.9's embedded Python 2.7 cannot import the ccp4i2 package, so the
# bridge modules are loaded by file path from CCP4I2_COOTBRIDGE_DIR.
# Coot's ~2000 scripting functions are globals of this __main__ frame, so
# they are collected here (where an imported module could not see them)
# and handed to the loader.

from __future__ import print_function

import imp
import os
import traceback

_bridge_dir = os.environ.get("CCP4I2_COOTBRIDGE_DIR")

ccp4i2_bridge_controller = None
try:
    _bridge = imp.load_source(
        "ccp4i2_api_client", os.path.join(_bridge_dir, "api_client.py"))
    _loader = imp.load_source(
        "ccp4i2_coot09", os.path.join(_bridge_dir, "coot09_loader.py"))
    _gui = imp.load_source(
        "ccp4i2_coot09_gui", os.path.join(_bridge_dir, "coot09_gui.py"))
    _fns = {}
    for _name in ("read_pdb", "make_and_draw_map", "read_cif_dictionary",
                  "set_map_colour", "set_molecule_name",
                  "coot_menubar_menu", "add_simple_coot_menu_menuitem",
                  "molecule_chooser_gui", "save_coordinates",
                  "molecule_name", "graphics_n_molecules",
                  "is_valid_model_molecule"):
        try:
            _fns[_name] = eval(_name)
        except NameError:
            pass
    ccp4i2_bridge_controller = _gui.install(_fns, _bridge, _loader)
except Exception:
    print("[ccp4i2-cootbridge] startup failed:")
    traceback.print_exc()


# Legacy COOTSCRIPTOUT scripts (from refmac/edstats/...) call methods on a
# ccp4i2Interface object that no longer exists; give them a no-op so they
# degrade politely instead of crashing startup.
class _CCP4i2InterfaceShim(object):
    def __getattr__(self, _name):
        def _noop(*_args, **_kwargs):
            print("[ccp4i2-cootbridge] legacy interface call ignored:", _name)
        return _noop


ccp4i2Interface = _CCP4i2InterfaceShim()

# Optional per-job extras, flagged by the wrapper through the environment.
if os.environ.get("CCP4I2_COOT_KEYBINDINGS"):
    try:
        file_to_preferences("template_key_bindings.py")  # noqa: F821
    except Exception:
        traceback.print_exc()

_script_file = os.environ.get("CCP4I2_COOTSCRIPTFILE")
if _script_file and os.path.isfile(_script_file):
    # A follow-on Coot script (COOTSCRIPTOUT). Executed in this frame so it
    # sees Coot's builtins and the ccp4i2Interface shim; __future__ imports
    # are stripped (only legal at a module top, not in exec'd text) and any
    # failure is contained so it cannot abort startup.
    try:
        _handle = open(_script_file)
        try:
            _source = "".join(
                _line for _line in _handle
                if "from __future__ import" not in _line)
        finally:
            _handle.close()
        exec(_source)
    except Exception:
        traceback.print_exc()
