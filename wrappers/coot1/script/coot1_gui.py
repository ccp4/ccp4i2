from qtgui.CCP4TaskWidget import CTaskWidget


class Ccoot1(CTaskWidget):
    TASKNAME = "coot1"

    def drawContents(self):
        self.openFolder(folderFunction="inputData")
        self.createLine(["subtitle", "Coordinates"])
        self.createLine(["widget", "XYZIN_LIST"])
        self.createLine(["subtitle", "Electron density maps"])
        self.createLine(["widget", "FPHIIN_LIST"])
        self.createLine(["subtitle", "Difference density maps"])
        self.createLine(["widget", "DELFPHIIN_LIST"])
        self.createLine(["subtitle", "Anomalous difference density maps"])
        self.createLine(["widget", "DELFPHIINANOM_LIST"])
        self.createLine(["subtitle", "Additional data"])
        self.createLine(["widget", "DICT"])
