"""
    modelcraft_gui.py: CCP4 GUI Project

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

from PySide2 import QtCore
from qtgui.CCP4TaskWidget import CTaskWidget


class modelcraft_gui(CTaskWidget):
    TASKMODULE = "model_building"
    TASKTITLE = "Autobuild with ModelCraft, Buccaneer and Nautilus"
    SHORTTASKTITLE = "ModelCraft"
    TASKNAME = "modelcraft"
    TASKVERSION = 0.1
    DESCRIPTION = "Automated model building of protein, nucleic acid and water"
    MGDISPLAYFILES = ["XYZIN"]
    WHATNEXT = ["coot_rebuild"]

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    @QtCore.Slot()
    def update_requirements(self):
        use_model = bool(self.container.controlParameters.USE_MODEL_PHASES)
        self.container.inputData.PHASES.setQualifier("allowUndefined", use_model)
        self.container.inputData.XYZIN.setQualifier("allowUndefined", not use_model)
        phases_widget = self.getWidget("PHASES")
        if hasattr(phases_widget, "jobCombo"):
            phases_combo_text = "..is not used" if use_model else "..must be selected"
            phases_widget.jobCombo.setItemText(0, phases_combo_text)
        phases_widget.validate()
        xyzin_widget = self.getWidget("XYZIN")
        xyzin_widget.setUndefinedAllowedBehaviour(not use_model)
        xyzin_widget.validate()

    @QtCore.Slot()
    def basic_changed(self):
        basic = bool(self.container.controlParameters.BASIC)
        self.container.controlParameters.CYCLES.set(5 if basic else 25)
        self.container.controlParameters.STOP_CYCLES.set(2 if basic else 4)

    def drawContents(self):
        self.openFolder(folderFunction="inputData")

        self.createLine(["subtitle", "Reflection data"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "F_SIGF"])
        self.createLine(["widget", "FREERFLAG"])
        self.createLine(
            [
                "widget",
                "USE_MODEL_PHASES",
                "label",
                (
                    "Get initial phases from refining the starting model"
                    " (uncheck to specify starting phases,"
                    " e.g. from experimental phasing)"
                ),
            ]
        )
        self.createLine(
            ["widget", "PHASES"], toggle=["USE_MODEL_PHASES", "open", [False]]
        )
        self.createLine(
            [
                "widget",
                "UNBIASED",
                "label",
                "Phases are unbiased "
                "and should be used as refinement restraints when the model is poor",
            ],
            toggle=["USE_MODEL_PHASES", "open", [False]],
        )
        self.closeSubFrame()

        self.createLine(["subtitle", "Asymmetric unit contents"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "ASUIN"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Starting model"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "XYZIN"])
        self.getWidget("XYZIN").showAtomSelection()
        self.closeSubFrame()

        self.createLine(["subtitle", "Options"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "BASIC", "label", "Run a quicker basic pipeline"])
        self.createLine(["label", "Run for", "widget", "CYCLES", "label", "cycles"])
        self.createLine(
            [
                "widget",
                "AUTO_STOP",
                "label",
                "Stop automatically if R-free does not improve in",
                "widget",
                "STOP_CYCLES",
                "label",
                "cycles",
            ]
        )
        self.createLine(
            [
                "widget",
                "SELENOMET",
                "label",
                "Build selenomethionine (MSE) instead of methionine (MET)",
            ]
        )
        self.createLine(["widget", "TWINNED", "label", "Use twinned refinement"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Optional pipeline steps"])
        self.openSubFrame(frame=True)
        self.createLine(["advice", "Uncheck to disable the following optional steps:"])
        self.createLine(
            [
                "widget",
                "SHEETBEND",
                "label",
                "Preliminary low-resolution refinement with Sheetbend",
            ]
        )
        self.createLine(
            ["widget", "PRUNING", "label", "Residue and chain pruning"],
            toggle=["BASIC", "open", [False]],
        )
        self.createLine(
            ["widget", "PARROT", "label", "Classical density modification with Parrot"]
        )
        self.createLine(
            [
                "widget",
                "DUMMY_ATOMS",
                "label",
                "Phase improvement through addition and refinement of dummy atoms",
            ],
            toggle=["BASIC", "open", [False]],
        )
        self.createLine(
            ["widget", "WATERS", "label", "Addition of waters"],
            toggle=["BASIC", "open", [False]],
        )
        self.createLine(
            ["widget", "SIDE_CHAIN_FIXING", "label", "Final side-chain fixing"],
            toggle=["BASIC", "open", [False]],
        )
        self.closeSubFrame()

        USE_MODEL_PHASES = self.container.controlParameters.USE_MODEL_PHASES
        USE_MODEL_PHASES.dataChanged.connect(self.update_requirements)
        self.update_requirements()

        BASIC = self.container.controlParameters.BASIC
        BASIC.dataChanged.connect(self.basic_changed)
        self.basic_changed()
