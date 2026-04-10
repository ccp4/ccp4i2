"Liz Potterton August 2012 - gui for mtz split"

from dataclasses import dataclass
from functools import partial
from PySide2 import QtCore, QtWidgets
from core.CCP4ErrorHandling import CErrorReport
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets


@dataclass
class DataType:
    description: str
    groupType: str
    contentFlag: int
    columnTypes: str
    labels: list


DATA_TYPES = [
    DataType("Anomalous intensities", "Obs", 1, "KMKM", ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]),
    DataType("Anomalous SFs", "Obs", 2, "GLGL", ["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"]),
    DataType("Mean intensities", "Obs", 3, "JQ", ["I", "SIGI"]),
    DataType("Mean SFs", "Obs", 4, "FQ", ["F", "SIGF"]),
    DataType("HL Phases", "Phs", 1, "AAAA", ["HLA", "HLB", "HLC", "HLD"]),
    DataType("Phi-FOM Phases", "Phs", 2, "PW", ["PHI", "FOM"]),
    DataType("Map coefficients", "MapCoeffs", 1, "FP", ["F", "PHI"]),
]


class CTaskSplitMtz(CCP4TaskWidget.CTaskWidget):
    TASKNAME = "splitMtz"
    TASKVERSION = 0.0
    TASKTITLE = "Import and Split MTZ into experimental data objects"
    TASKMODULE = "data_entry"
    SHORTTASKTITLE = "Import and Split MTZ"
    DESCRIPTION = "Select groups of columns from the MTZ file (csplitmtz)"
    ERROR_CODES = {200: {"description": "There are no selected column groups"}}

    def drawContents(self):
        self.openFolder(folderFunction="inputData", followFrom=False)
        self.createLine(["subtitle", "Select experimental data file"])
        self.createLine(["widget", "HKLIN"])
        self.createLine(["subtitle", "Select column groups"])
        self.createLine(["widget", "COLUMNGROUPLIST"])
        self.container.inputData.HKLIN.dataChanged.connect(self.handleSelectHklin)

        line = self.createLine(["subtitle", "Select individual column data for"])
        self.butGroup = QtWidgets.QButtonGroup(self)
        for idx, dataType in enumerate(DATA_TYPES):
            if idx % 4 == 0:
                line = self.createLine([])
            but = QtWidgets.QRadioButton(dataType.description, self)
            self.butGroup.addButton(but, idx)
            line.addWidget(but)
            but.clicked.connect(partial(self.drawSelectColumns, dataType))

        line = self.createLine([])
        self.selectFrame = QtWidgets.QFrame(self)
        self.selectFrame.setLayout(QtWidgets.QGridLayout())
        line.addWidget(self.selectFrame)

    @QtCore.Slot(DataType)
    def drawSelectColumns(self, dataType):
        for iR in (0, 1):
            for iC in (0, 1):
                layoutItem = self.selectFrame.layout().itemAtPosition(iR, iC)
                if layoutItem is not None:
                    layoutItem.widget().hide()
                    layoutItem.widget().deleteLater()

        HKLIN = self.container.inputData.HKLIN
        self.comboList = []
        for nframes, (cType, text) in enumerate(zip(dataType.columnTypes, dataType.labels)):
            frame = QtWidgets.QFrame(self)
            frame.setLayout(QtWidgets.QHBoxLayout())
            frame.layout().setContentsMargins(0, 0, 0, 0)
            frame.layout().setSpacing(0)
            label = QtWidgets.QLabel(text, frame)
            frame.layout().addWidget(label)
            self.comboList.append(CCP4Widgets.CComboBox(frame))
            self.comboList[-1].setMinimumContentsLength(25)
            self.comboList[-1].setEditable(False)
            self.comboList[-1].addItem("Select column..")
            for col in HKLIN.fileContent.getListOfColumns([cType]):
                self.comboList[-1].addItem(f"{col.dataset}/{col.columnLabel}")
            frame.layout().addWidget(self.comboList[-1])
            self.selectFrame.layout().addWidget(frame, nframes // 2, nframes % 2)

    @QtCore.Slot()
    def handleSelectHklin(self):
        inputData = self.container.inputData
        inputData.HKLIN.loadFile()
        inputData.COLUMNGROUPLIST.set(inputData.HKLIN.fileContent.getColumnGroups())

        # Redo the individual column input
        butId = self.butGroup.checkedId()
        if butId >= 0:
            self.drawSelectColumns(DATA_TYPES[butId])

    def fix(self):
        self.getWidget("COLUMNGROUPLIST").updateModelFromView()
        nGroupsSelected = self.getWidget("COLUMNGROUPLIST").numberSelected()

        userGroup = self.container.inputData.USERCOLUMNGROUP
        userGroup.unSet()
        butId = self.butGroup.checkedId()
        if butId >= 0:
            selectedCols = []
            for c in self.comboList:
                text = str(c.currentText())
                if c.currentIndex() != 0 and text not in selectedCols:
                    selectedCols.append(text)
            if len(selectedCols) == len(self.comboList):
                userGroup.columnGroupType.set(DATA_TYPES[butId].groupType)
                userGroup.contentFlag.set(DATA_TYPES[butId].contentFlag)
                for idx, col in enumerate(selectedCols):
                    dataset, label = col.split("/")
                    columnType = DATA_TYPES[butId].columnTypes[idx]
                    userGroup.columnList.addItem()
                    userGroup.columnList[-1].columnLabel.set(label)
                    userGroup.columnList[-1].columnType.set(columnType)
                    userGroup.columnList[-1].dataset.set(dataset)
                userGroup.dataset.set(userGroup.columnList[0].dataset)
        elif nGroupsSelected == 0:
            return CErrorReport(self.__class__, 200, stack=False)
        return CErrorReport()
