"""
    density_calculator_gui.py: CCP4 GUI Project

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

from baselayer import QtCore, QtWidgets
from qtgui.CCP4TaskWidget import CTaskWidget


class density_calculator_gui(CTaskWidget):
    TASKMODULE = "expt_data_utility"
    TASKTITLE = "Density Calculator"
    TASKNAME = "density_calculator"
    TASKVERSION = 0.1
    DESCRIPTION = "Calculate a map for a structure"
    MGDISPLAYFILES = ["XYZIN"]
    WHATNEXT = ["coot_rebuild"]

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction="inputData")

        self.createLine(["subtitle", "Structure"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "XYZIN"])
        self.getWidget("XYZIN").showAtomSelection()
        self.closeSubFrame()

        self.createLine(["subtitle", "Options"])
        self.openSubFrame(frame=True)
        self.createLine(["label", "Scattering form factor", "widget", "FORM_FACTOR"])
        self.createLine(["label", "Resolution limit / &#8491;", "widget", "D_MIN"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Advanced Options"])
        self.openSubFrame(frame=True)
        rate_line = self.createLine(["label", "Oversampling rate", "widget", "RATE"])
        self.grid_spacing_label = QtWidgets.QLabel("")
        rate_line.addWidget(self.grid_spacing_label)
        self.createLine(["label", "Blurring mode", "widget", "BLUR_MODE"])
        self.createLine(
            ["label", "Blurring B-factor / &#8491;<sup>2</sup>", "widget", "BLUR"],
            toggle=["BLUR_MODE", "open", "custom"],
        )
        self.createLine(
            [
                "widget",
                "UNBLUR",
                "label",
                "Unblur when calculating the reciprocal-space map coefficients",
            ],
            toggle=["BLUR_MODE", "close", "none"],
        )
        self.createLine(["label", "Density cutoff", "widget", "CUTOFF"])
        self.createLine(
            [
                "widget",
                "MOTT_BETHE",
                "label",
                "Approximate electron scattering factors using the Mott-Bethe formula",
            ],
            toggle=["FORM_FACTOR", "open", "xray"],
        )
        self.closeSubFrame()

        self.container.controlParameters.D_MIN.dataChanged.connect(self.update)
        self.container.controlParameters.RATE.dataChanged.connect(self.update)
        self.update()

    @QtCore.Slot()
    def update(self):
        try:
            d_min = self.container.controlParameters.D_MIN
            rate = self.container.controlParameters.RATE
            spacing = d_min / (2 * rate)
            text = f"<span>(grid spacing: {spacing:.2f} &#8491;)</span>"
            self.grid_spacing_label.setText(text)
        except AttributeError:
            pass
