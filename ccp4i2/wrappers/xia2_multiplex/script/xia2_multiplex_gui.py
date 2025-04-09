#
#  Copyright (C) 2022 UKRI/STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#

from PySide2 import QtCore, QtWidgets

from .... import qtgui
from ....core import CCP4Container
from ....qtgui.CCP4TaskWidget import CTaskWidget
from .find_integrated import find_integrated


class FindIntegratedWorker(QtCore.QObject):
    finished = QtCore.Signal()

    def __init__(self, path):
        self._path = path
        super(FindIntegratedWorker, self).__init__()

    def run(self):
        self.results = find_integrated(self._path)
        self.finished.emit()


class xia2_multiplex_gui(CTaskWidget):

    # Subclass CTaskWidget to give specific task window
    TASKTITLE = "Automated combination of datasets using xia2.multiplex"
    DESCRIPTION = "Select previous xia2 runs"
    TASKNAME = "xia2_multiplex"
    TASKMODULE = "test"
    # TASKMODULE = "data_reduction"
    TASKVERSION = 0.0
    RANK = 1
    SHORTTASKTITLE = "multiplex"
    TASKVERSION = 0.1

    def drawContents(self):

        # Input data
        self.openFolder(folderFunction="inputData")

        self.createLine(
            [
                "tip",
                "Select a directory in which to search for DIALS integrated files",
                "subtitle",
                "Search for DIALS integrated files",
            ]
        )

        self.createLine(
            [
                "widget",
                "-title",
                "Find DIALS integrated files",
                "SEARCH_ROOT_DIR",
            ]
        )

        self.createLine(
            [
                "widget",
                "-title",
                "DIALS integrated.refl files (each must have an associated .expt file)",
                "DIALS_INTEGRATED",
            ]
        )

        self.closeSubFrame()

        # Basic parameters
        self.createLine(
            ["tip", "Subset of xia2 control parameters", "subtitle", "Basic parameters"]
        )
        self.simpleAutoGenerate(container=self.container.controlParameters)

        # Advanced parameters in a separate folder
        self.openFolder(
            folderFunction="controlParameters",
            title="Advanced parameters",
            drawFolder=self.drawAdvanced,
        )

        self.connectDataChanged("SEARCH_ROOT_DIR", self.handleSearchRootDir)

        return

    @QtCore.Slot()
    def handleSearchRootDir(self):

        self.getWidget("SEARCH_ROOT_DIR").setEnabled(False)
        self.getWidget("DIALS_INTEGRATED").setEnabled(False)
        # Set up a thread that is owned by the application instance, to avoid
        # crashes due to thread object going out of scope
        self.find_integrated_thread = QtCore.QThread(
            parent=QtWidgets.QApplication.instance()
        )
        self.find_integrated_worker = FindIntegratedWorker(
            str(self.container.inputData.SEARCH_ROOT_DIR)
        )
        self.find_integrated_worker.moveToThread(self.find_integrated_thread)
        self.find_integrated_thread.started.connect(self.find_integrated_worker.run)
        self.find_integrated_worker.finished.connect(self.handleSearchRootDirFinished)
        self.find_integrated_worker.finished.connect(self.find_integrated_thread.quit)
        self.find_integrated_thread.finished.connect(
            self.find_integrated_thread.deleteLater
        )
        self.find_integrated_thread.start()

    @QtCore.Slot()
    def handleSearchRootDirFinished(self):

        dials_integrated = self.getWidget("DIALS_INTEGRATED")
        self.getWidget("SEARCH_ROOT_DIR").setEnabled(True)
        dials_integrated.setEnabled(True)
        results = self.find_integrated_worker.results
        self.container.inputData.DIALS_INTEGRATED.unSet()
        self.container.inputData.DIALS_INTEGRATED.set([r + ".refl" for r in results])

        dials_integrated.updateViewFromModel()

        # Attempt to set input as valid to get rid of red highlights
        # (does not currently work)
        dials_integrated.validate()
        for i in range(dials_integrated.getNofLines()):
            dials_integrated.setValidity(i)

        # Update the widget to redraw
        dials_integrated.updateViewFromModel()
        dials_integrated.validate()
        # I do not know why this is necessary, but it seems to be.
        childer = dials_integrated.parent().findChildren(QtWidgets.QWidget)
        for child in childer:
            if type(child) is qtgui.CCP4Widgets.CListViewListWidget:
                if hasattr(child.parent(), "isValid"):
                    if child.parent().isValid:
                        child.setStyleSheet("QFrame { background-color:white}")
                    else:
                        child.setStyleSheet(
                            'QFrame [isValid="false"] { background: qlineargradient( x1:0.15 y1:0, x2:1 y2:0, stop:0 transparent, stop:1 #f55); border:0px; }'
                        )
        self.validate()

    def drawAdvanced(self):
        self.nestedAutoGenerate(
            container=self.container.controlParameters, expertLevel=["0", "1"]
        )

    def simpleAutoGenerate(self, container):
        """Autogenerate a specific selection of basic parameters"""

        self.autoGenerate(
            container.symmetry,
            selection={"includeParameters": ["symmetry__space_group"]},
        )
        self.autoGenerate(
            container.resolution,
            selection={"includeParameters": ["resolution__d_max", "resolution__d_min"]},
        )
        self.autoGenerate(
            container.filtering,
            selection={"includeParameters": ["filtering__method", "resolution__d_min"]},
        )
        self.createLine(["advice", "Clustering"])
        self.autoGenerate(
            container,
            selection={"includeParameters": ["max_clusters", "cluster_method"]},
        )

        # TODO Add here:
        # - filtering mode either image group or dataset
        # - tidy up intermediate files checkbox

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
