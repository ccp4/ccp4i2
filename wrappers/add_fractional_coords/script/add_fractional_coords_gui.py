"""
    add_fractional_coords_gui.py: CCP4 GUI Project

    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""

from ccp4i2.baselayer import QtCore
from qtgui.CCP4TaskWidget import CTaskWidget


class add_fractional_coords_gui(CTaskWidget):
    TASKMODULE = "model_data_utility"
    TASKTITLE = "Add Fractional Coordinates"
    TASKNAME = "add_fractional_coords"
    TASKVERSION = 0.1
    DESCRIPTION = "Calculate fractional coordinates and output them in mmCIF format"
    MGDISPLAYFILES = ["XYZIN"]
    WHATNEXT = []

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction="inputData")

        self.createLine(["subtitle", "Input model"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "XYZIN"])
        self.getWidget("XYZIN").showAtomSelection()
        self.closeSubFrame()
