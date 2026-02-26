from qtgui import CCP4TaskWidget


class CTaskAreaimolUI(CCP4TaskWidget.CTaskWidget):
    TASKNAME = "areaimol"
    TASKTITLE = "Solvent accessible surface - AREAIMOL"
    TASKMODULE = "model_data_utility"
    SHORTTASKTITLE = "AREAIMOL"
    PROGRAMHELP = "areaimol"
    DESCRIPTION = "Solvent accessible surface calculation with AREAIMOL"

    def diffmode_changed(self):
        compare = self.container.controlParameters.DIFFMODE == "COMPARE"
        self.container.inputData.XYZIN2.setQualifier("allowUndefined", not compare)
        if not compare:
            self.container.inputData.XYZIN2.unSet()
        xyzin2_widget = self.getWidget("XYZIN2")
        xyzin2_widget.setUndefinedAllowedBehaviour(not compare)
        xyzin2_widget.validate()

    def drawContents(self):
        self.setProgramHelpFile("areaimol")
        self.openFolder(folderFunction="inputData", title="Input Data", followFrom=False)
        self.createLine(["subtitle", "Input"])
        self.openSubFrame(frame=[True])
        self.createLine(["subtitle", "Input Structure"])
        self.createLine(["tip", "This is the structure which the area will be calculated.", "widget", "XYZIN"])
        self.createLine(["tip", "This keyword controls the program function, the data required, and how it is processed and analysed", "widget", "DIFFMODE"])
        self.createLine(["label", "The symmetry of the molecule, e.g. P212121", "widget", "SYMMETRY"], toggle=["DIFFMODE", "open", ["IMOL"]])
        self.createLine(["subtitle", "Second structure to compare with"], toggle=["DIFFMODE", "open", ["COMPARE"]])
        self.createLine(["tip", "A second set of input coordinates, and is only used in DIFFMODE COMPARE.", "widget", "XYZIN2"], toggle=["DIFFMODE", "open", ["COMPARE"]])
        self.createLine(["label", "Output mode", "widget", "OUTPUT_MODE"], toggle=["DIFFMODE", "open", ["OFF"]])
        self.createLine(["label", "Output mode", "widget", "OUTPUT_MODE_COMPARE"], toggle=["DIFFMODE", "open", ["COMPARE", "IMOL"]])
        self.createLine(["label", "Residue names to exclude (space separated)", "widget", "EXCLUDE"])
        self.closeSubFrame()
        self.connectDataChanged("DIFFMODE", self.diffmode_changed)
