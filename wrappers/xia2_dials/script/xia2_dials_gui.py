#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

from baselayer import QtGui, QtWidgets, QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Container


class xia2_dials_gui(CTaskWidget):

    # Subclass CTaskWidget to give specific task window
    TASKTITLE = "Automated integration of images with DIALS using xia2"
    DESCRIPTION = "Select a directory containing images and integrate them"
    TASKNAME = "xia2_dials"
    TASKMODULE = "data_processing"
    TASKVERSION = 0.0
    RANK = 1
    SHORTTASKTITLE = "xia2/dials"
    TASKVERSION = 0.1

    def drawContents(self):

        # Input data
        self.openFolder(folderFunction="inputData")
        self.createLine(
            [
                "tip",
                "Provide either the parent directory of the "
                "dataset(s), or pick one image from each",
                "subtitle",
                "Locate datasets",
            ]
        )
        self.openSubFrame(frame=True)

        self.createLine(
            ["widget", "-title", "Choose one image from each dataset...", "IMAGE_FILE"]
        )
        self.createLine(
            ["advice", "...Or let xia2 find datasets under a parent directory"]
        )
        self.createLine(["widget", "IMAGE_DIRECTORY"])
        self.closeSubFrame()

        # Basic parameters
        self.createLine(
            ["tip", "Subset of xia2 control parameters", "subtitle", "Basic parameters"]
        )
        self.nestedAutoGenerate(
            container=self.container.controlParameters, expertLevel=["0"]
        )

        # Advanced parameters in a separate folder
        folder = self.openFolder(
            folderFunction="controlParameters",
            title="Advanced parameters",
            drawFolder=self.drawAdvanced,
        )

        self.connectDataChanged("IMAGE_FILE", self.handleImageFile)
        self.connectDataChanged("IMAGE_DIRECTORY", self.handleImageDirectory)
        self.connectDataChanged("dials__index__method", self.handleIndexMethod)
        self.handleImageFile()
        self.handleImageDirectory()
        self.handleIndexMethod()

        return

    @QtCore.Slot()
    def handleImageFile(self):
        image_supplied = [
            e for e in self.container.inputData.IMAGE_FILE if str(e.imageFile).strip()
        ]
        image_supplied = len(image_supplied) > 0
        self.container.inputData.IMAGE_DIRECTORY.setQualifiers(
            {"allowUndefined": image_supplied}
        )
        self.getWidget("IMAGE_DIRECTORY").validate()

    @QtCore.Slot()
    def handleImageDirectory(self):
        if self.container.inputData.IMAGE_DIRECTORY.isSet():
            self.container.inputData.IMAGE_FILE.setQualifiers({"listMinLength": 0})
        else:
            self.container.inputData.IMAGE_FILE.setQualifiers({"listMinLength": 1})
        self.getWidget("IMAGE_FILE").validate()

    @QtCore.Slot()
    def handleIndexMethod(self):
        try:
            unit_cell_required = str(self.container.controlParameters.dials.dials__index.dials__index__method) == "real_space_grid_search"
        except:
            #Default option of index method fft3d. We can get here in xia2_xds
            unit_cell_required = False
        self.container.controlParameters.xia2.xia2__settings.xia2__settings__unit_cell.setQualifiers({'allowUndefined':(not unit_cell_required)})
        self.getWidget("xia2__settings__unit_cell").validate()


    def drawAdvanced(self):
        self.nestedAutoGenerate(
            container=self.container.controlParameters, expertLevel=["0", "1"]
        )

    def nestedAutoGenerate(self, container, expertLevel=["0"]):
        """Autogenerate a GUI from nested containers using only elements at
        the specified expertLevels"""

        dataOrder = container.dataOrder()
        contents = [getattr(container, name) for name in dataOrder]

        # Filter by expertLevel
        expert = [
            c.qualifiers("guiDefinition").get("expertLevel", "0") for c in contents
        ]
        dataOrder = [d for d, e in zip(dataOrder, expert) if e in expertLevel]
        contents = [c for c, e in zip(contents, expert) if e in expertLevel]

        # If this container has any content other than nested containers,
        # autogenerate that content
        data = [
            name
            for name, model in zip(dataOrder, contents)
            if not isinstance(model, CCP4Container.CContainer)
        ]
        subcontainers = [
            model for model in contents if isinstance(model, CCP4Container.CContainer)
        ]
        if data:
            # only open a sub-frame if there is a label available
            make_subframe = container.qualifiers("guiLabel") is not NotImplemented
            if make_subframe:
                self.openSubFrame(frame=True)
            self.autoGenerate(container, selection={"includeParameters": data})
            for c in subcontainers:
                self.nestedAutoGenerate(c, expertLevel=expertLevel)
            if make_subframe:
                self.closeSubFrame()
        else:
            # just drill down without opening a sub frame
            for c in subcontainers:
                self.nestedAutoGenerate(c, expertLevel=expertLevel)
        return
