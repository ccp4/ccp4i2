import os
import gtk

from ccp4i2.wrappers.SyncToDjango.script import testTree
for cootFunction in [coot_menubar_menu, add_simple_coot_menu_menuitem, read_pdb, set_molecule_name, make_and_draw_map, auto_read_make_and_draw_maps, read_cif_dictionary, save_coordinates, save_state_file, molecule_chooser_gui, residue_centre, set_rotation_centre]:
    setattr(testTree, cootFunction.__name__, cootFunction)

btve = testTree.DjangoInteractionWidget()
