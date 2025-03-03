#
#  Copyright (C) 2017 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on ideas and code by Nat Echols and Martin Noble.
#

"""Create xia2_xds.def.xml from PHIL parameters"""

import sys
import os

# Nasty trick required to import Xia2DialsTaskCreator when running with
# ccp4-python
this_dir = os.path.dirname(os.path.realpath(__file__))
ccp4i2_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_dir)))
sys.path.append(ccp4i2_dir)
from wrappers.xia2_dials.script.create_def_xml import Xia2DialsTaskCreator


class Xia2XDSTaskCreator(Xia2DialsTaskCreator):
    def __init__(self, debug=False):

        Xia2DialsTaskCreator.__init__(self, debug)
        self.fmt_dic["PLUGINNAME"] = "xia2_xds"

        self._elts_to_remove = ["dials", "strategy", "xia2__settings__input__image"]

    def __call__(self):

        # Modify the pipeline parameter to only allow XDS versions
        for cont in self.phil_tree.iter():
            if cont.get("id") == "xia2__settings__pipeline":
                qual = cont.find("qualifiers")
                qual.find("toolTip").text = (
                    "Select XDS processing pipeline. "
                    "3d: XDS, XSCALE, LABELIT   "
                    "3di: as 3d, but use 3 wedges for indexing  "
                    "3dii: XDS, XSCALE, using all images for autoindexing   "
                    "3dd: as 3d, but use DIALS for indexing"
                )
                qual.find("enumerators").text = "3d,3dd,3di,3dii"
                qual.find("default").text = "3dii"
                break

        super(Xia2XDSTaskCreator, self).__call__()


if __name__ == "__main__":

    x2xtc = Xia2XDSTaskCreator()
    x2xtc()
