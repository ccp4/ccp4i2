from qtgui import CCP4TaskWidget


class NucleoFindGUI(CCP4TaskWidget.CTaskWidget):
    TASKNAME = "nucleofind"
    TASKMODULE = "model_building"
    TASKTITLE = "Predict regions of nucleic acids with NucleoFind"
    DESCRIPTION = "Use NucleoFind to predict the regions of nucleic acid phosphates, sugars and base in an electron density map or Coulomb potential map."
    SHORTTASKTITLE = "NucleoFind"

    def drawContents(self):
        self.setProgramHelpFile("nucleofind")
        self.openFolder(folderFunction="inputData", title="Input Data")

        self.createLine(["subtitle", "Input map"])
        self.openSubFrame(frame=[True])
        self.createLine(["widget", "FPHIIN"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Options", ""])
        self.openSubFrame(frame=[True])
        self.createLine(["label", "CPU threads", "widget", "THREADS", "label", "(higher values require more memory)"])
        self.createLine(["widget", "GPU", "label", "Use GPU acceleration if it is available"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Advanced Options", ""])
        self.openSubFrame(frame=[True])
        self.createLine(["label", "Truncate reflections up to", "widget", "RESOLUTION", "label", "Å resolution"])
        self.createLine(["widget", "FULL_CELL", "label", "Predict over entire unit cell"])
        self.createLine(["label", "Overlap predicted boxes by ", "widget", "OVERLAP", "label", "grid points (lowering this may increase prediction accuracy but will increase runtime)"])
        self.createLine(["label", "Output map values", "widget", "OUTPUT_TYPE"])
        self.closeSubFrame()
