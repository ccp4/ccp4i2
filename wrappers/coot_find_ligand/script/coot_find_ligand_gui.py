from qtgui.CCP4TaskWidget import CTaskWidget


class Ccoot_find_ligand(CTaskWidget):
    TASKNAME = "coot_find_ligand"
    TASKVERSION = 0.0
    TASKMODULE = "wrappers"
    TASKTITLE = "Find a ligand with Coot API"
    DESCRIPTION = "Find and flexibly fit a ligand (non-interactive Coot)"

    def drawContents(self):
        self.setProgramHelpFile("coot_find_ligand")
        self.openFolder(folderFunction="inputData")

        self.createLine(["subtitle", "Input data"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "XYZIN"])
        self.createLine(["widget", "FPHI"])
        self.createLine(["widget", "DICT"])
        self.createLine(["label", "Component ID", "widget", "COMPID"])
        self.closeSubFrame()

        self.createLine(["subtitle", "Options"])
        self.openSubFrame(frame=True)
        self.createLine(["label", "Map RMS threshold", "widget", "THRESHOLD"])
        self.createLine(["widget", "FLEXIBLE", "label", "Use flexible fitting"])
        self.createLine(
            ["label", "Number of conformers", "widget", "CONFORMERS"],
            toggle=["FLEXIBLE", "open", [True]],
        )
        self.closeSubFrame()
