# -*- coding: utf-8 -*-
"""Coot 0.9 adapter for the CCP4i2 bridge (Python 2.7 compatible).

The 0.9 half of the two-renderer design: the shared data layer
(api_client) does everything toolkit-neutral; this module translates
the load plan into Coot 0.9's flat scripting namespace. No GUI is
built here - the GTK2 browser is a later tier; saving in 0.9 sessions
continues through the classic generated-script menu items.

Coot 0.9's ~2000 scripting calls are globals of __main__ and are NOT
visible to imported modules, so the generated stub script must hand
them in. Expected stub shape (written by a coot_rebuild-successor
task):

    import imp
    bridge = imp.load_source("ccp4i2_api_client", r"<...>/api_client.py")
    loader = imp.load_source("ccp4i2_coot09", r"<...>/coot09_loader.py")
    loader.install_and_load({
        "read_pdb": read_pdb,
        "make_and_draw_map": make_and_draw_map,
        "read_cif_dictionary": read_cif_dictionary,
        "set_map_colour": set_map_colour,
        "set_molecule_name": set_molecule_name,
    }, bridge)

Only the functions actually present need be passed; missing ones
degrade gracefully.
"""

from __future__ import absolute_import, print_function

import traceback


def _log(message):
    print("[ccp4i2-cootbridge-0.9] {0}".format(message))


def install_and_load(coot_functions, bridge_module):
    """Load the launching job's data through Coot 0.9 builtins.

    coot_functions: dict of name -> callable from Coot's flat namespace.
    bridge_module:  the api_client module (loaded via imp.load_source).
    """
    try:
        config = bridge_module.BridgeConfig()
        client = bridge_module.CootBridgeClient(config)
    except Exception:
        _log("failed to configure bridge:\n" + traceback.format_exc())
        return
    try:
        plan = bridge_module.load_plan(config, client)
    except Exception:
        _log("could not build load plan:\n" + traceback.format_exc())
        return
    for item in plan:
        try:
            load_item(item, coot_functions)
        except Exception:
            _log("load failed for {0}:\n{1}".format(
                item.get("path"), traceback.format_exc()))


def load_item(item, coot_functions):
    """Dispatch one {kind, path, label} through 0.9 flat-namespace calls."""
    kind = item["kind"]
    path = str(item["path"])
    label = item.get("label")
    read_pdb = coot_functions.get("read_pdb")
    make_map = coot_functions.get("make_and_draw_map")
    if kind == "coordinates" and read_pdb is not None:
        imol = read_pdb(path)
        set_name = coot_functions.get("set_molecule_name")
        if label and set_name is not None and imol >= 0:
            set_name(imol, label)
    elif kind == "map_2fofc" and make_map is not None:
        make_map(path, "F", "PHI", "PHI", 0, 0)
    elif kind == "map_fofc" and make_map is not None:
        make_map(path, "F", "PHI", "PHI", 0, 1)
    elif kind == "map_anom" and make_map is not None:
        imap = make_map(path, "F", "PHI", "PHI", 0, 1)
        set_colour = coot_functions.get("set_map_colour")
        if set_colour is not None and imap >= 0:
            set_colour(imap, 0.75, 0.9, 0.75)
    elif kind == "dictionary":
        read_dict = coot_functions.get("read_cif_dictionary")
        if read_dict is not None:
            read_dict(path)
    else:
        _log("no handler for kind {0!r}, skipped".format(kind))
