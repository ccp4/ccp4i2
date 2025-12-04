#
#  Copyright (C) 2024 UKRI/STFC Rutherford Appleton Laboratory, UK.
#
#  Author: Martin Maly, David Waterman
#

from baselayer import QtCore, QtWidgets
from qtgui.CCP4TaskWidget import CTaskWidget
from core import CCP4Container
import qtgui
from .find_expt_refl import find_expt_refl
from dxtbx.serialize import load
from cctbx import uctbx
from scitbx.array_family import flex
import os


class FindIntegratedWorker(QtCore.QObject):
    finished = QtCore.Signal()

    def __init__(self, path, preference):
        self._path = path
        self._preference = preference
        super(FindIntegratedWorker, self).__init__()

    def run(self):
        self.results = find_expt_refl(self._path, self._preference)
        self.finished.emit()


class xia2_ssx_reduce_gui(CTaskWidget):

    # Subclass CTaskWidget to give specific task window
    TASKTITLE = "Reduction of serial datasets using xia2.ssx_reduce"
    DESCRIPTION = "Select integrated data from xia2.ssx or dials.stills_process"
    TASKNAME = "xia2_ssx_reduce"
    TASKMODULE = "data_reduction"
    TASKVERSION = 0.0
    RANK = 1
    SHORTTASKTITLE = "xia2.ssx_reduce"
    TASKVERSION = 0.1

    # def __init__(self,parent):
    #     CTaskWidget.__init__(self,parent)

    def isValid(self):
        if self.container.inputData.DIALS_INTEGRATED.isSet():
            if len(self.container.inputData.DIALS_INTEGRATED) >= 1:
                firstReflPath = self.container.inputData.DIALS_INTEGRATED[0].__str__()
                if os.path.isfile(firstReflPath):
                    if "DataFiles_scaled" in os.path.split(firstReflPath)[-1]:
                        self.container.inputData.SEARCH_PREFERENCE.set('scaled')
                        self.updateScaleMerge1()
        #     if self.getWidget('followFrom') is None: return
        #     followJobId = self.getWidget('followFrom').currentJobId()
        #     if followJobId is not None:
        #         container_follow = CCP4Modules.PROJECTSMANAGER().getJobParams(jobId=followJobId)
        return CTaskWidget.isValid(self)


    def drawContents(self):

        # Input data
        self.openFolder(folderFunction="inputData")

        self.createLine(
            [
                "tip",
                "Select a directory in which to search for xia2.ssx integrated files",
                "subtitle",
                "Search for xia2.ssx integrated files",
            ]
        )
        self.createLine(
            [
                'label', 'Search for',
                'widget',
                '-guiMode',
                'radio',
                'SEARCH_PREFERENCE',
                'label', 'files',
            ]
        )
        self.connectDataChanged('SEARCH_PREFERENCE', self.updateScaleMerge1)
        self.createLine(
            [
                "widget",
                "-title",
                "Search for xia2.ssx integrated files",
                "SEARCH_ROOT_DIR",
            ]
        )
        self.createLine(
            [
                "widget",
                "-title",
                "integrated.refl or scaled.refl files (each must have an associated .expt file)",
                "DIALS_INTEGRATED",
            ]
        )
        self.closeSubFrame()

        # Basic parameters
        self.openSubFrame(frame=[True])
        self.createLine(
            ["tip", "Subset of xia2.ssx_reduce control parameters", "subtitle", "Basic parameters"]
        )
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["d_min"]},
        )
        self.autoGenerate(
            self.container.inputData,
            selection={"includeParameters": ["reference"]},
        )
        self.closeSubFrame()

        # Advanced parameters in a separate folder
        self.openFolder(folderFunction="controlParameters", title="Advanced parameters")
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["partiality_threshold"]},
        )
        self.createLine(
            ['label', 'Keep anomalous pairs separate during scaling', 'widget', 'scaling__anomalous'],
            toggle=['workflow__steps', 'open', ['scale+merge']],
        )
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["dials_cosym_phil_d_min"]},
        )

        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Workflow options'])
        self.createLine(['label', 'Lattice tolerance for triggering mis-indexing assessment', 'widget', 'symmetry__lattice_symmetry_max_delta'])
        self.createLine(['label', 'Workflow: ', 'widget', '-guiMode', 'radio', 'workflow__steps' ] )
        self.connectDataChanged('workflow__steps', self.updateScaleMerge2)
        self.createLine(
            ["advice", "Reduction of batch size and number of processors prevents system memory issues."]
        )
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["reduction_batch_size"]},
        )
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["multiprocessing__nproc"]},
        )
        self.closeSubFrame()

        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Unit cell filtering'])
        self.createLine(['advice', 'A filter is applied based on the median cell and the absolute tolerances.'])
        self.createLine(['label', 'Calculated median cell', 'widget', 'MEDIAN_CELL'])
        self.getWidget('MEDIAN_CELL').setEnabled(False)
        self.getWidget('MEDIAN_CELL').setFixedWidth(300)
        self.createLine(['label', 'Absolute angle tolerance (degrees)', 'widget', 'clustering__absolute_angle_tolerance'])
        self.createLine(['label', 'Absolute length tolerance (A)', 'widget', 'clustering__absolute_length_tolerance'])
        self.createLine(['label', 'Instead of using the median cell, use these central cell values for the cell filtering<br>(Separate unit cell parameters using comma)', 'widget', 'clustering__central_unit_cell'])
        self.getWidget('clustering__central_unit_cell').setFixedWidth(300)
        self.createLine(['label', 'Instead of filtering based on cell values and tolerances, use a clustering approach to select a cell cluster.<br>Clustering threshold (Andrews–Bernstein distance)', 'widget', 'clustering__threshold'])
        self.closeSubFrame()

        self.openSubFrame(frame=[True])
        self.createLine(['label', 'Space group for scaling and merging', 'widget', 'symmetry__space_group'])
        self.closeSubFrame()

        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Reference model options (if generating reference from PDB model)'])
        self.createLine(['label', 'Average solvent density (k-sol)', 'widget', 'reference_model__k_sol'])
        self.createLine(['label', 'Average solvent B-factor (B-sol)', 'widget', 'reference_model__b_sol'])
        self.closeSubFrame()

        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Splitting mixed-condition data'])
        self.createLine(['label', 'Dose series - number of repeated measurements at each point', 'widget', 'dose_series_repeat'])
        self.autoGenerate(
            self.container.inputData,
            selection={"includeParameters": ["grouping"]},
        )
        self.closeSubFrame()

        self.connectDataChanged("SEARCH_ROOT_DIR", self.handleSearchRootDir)

        return

    def updateScaleMerge1(self):
        if self.container.inputData.SEARCH_PREFERENCE == 'scaled':
            self.container.controlParameters.workflow.workflow__steps.set('merge')
            self.getWidget('reference').setEnabled(False)
            self.container.inputData.reference.unSet()
        elif self.container.inputData.SEARCH_PREFERENCE == 'integrated':
            self.container.controlParameters.workflow.workflow__steps.set('scale+merge')
            self.getWidget('reference').setEnabled(True)

    def updateScaleMerge2(self):
        if self.container.controlParameters.workflow.workflow__steps == 'merge':
            self.getWidget('reference').setEnabled(False)
            self.container.inputData.reference.unSet()
        elif self.container.controlParameters.workflow.workflow__steps == 'scale+merge':
            self.getWidget('reference').setEnabled(True)

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
            str(self.container.inputData.SEARCH_ROOT_DIR),
            str(self.container.inputData.SEARCH_PREFERENCE)
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

        # Update median cell
        exptPaths = []
        for i, path in enumerate(results):
            exptPath = results[i] + ".expt"
            if os.path.isfile(exptPath):
                exptPaths.append(exptPath)
        try:
            expts = [load.experiment_list(path, check_format=False) for path in exptPaths]
            uc_params = [flex.double() for _ in range(6)]
            for experiments in expts:
                for c in experiments.crystals():
                    for i, p in enumerate(c.get_unit_cell().parameters()):
                        uc_params[i].append(p)
            medianCell = uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])
            medianCellStr = str(medianCell).replace("(", "").replace(")", "")
            self.container.controlParameters.MEDIAN_CELL.set(medianCellStr)
        except:
            print("WARNING: Median unit cell could not be calculated.")
