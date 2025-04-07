from ....qtgui.CCP4TaskWidget import CTaskWidget


class add_fractional_coords_gui(CTaskWidget):
    TASKMODULE = "model_data_utility"
    TASKTITLE = "Add Fractional Coordinates"
    TASKNAME = "add_fractional_coords"
    TASKVERSION = 0.1
    DESCRIPTION = "Calculate fractional coordinates and output them in mmCIF format"
    MGDISPLAYFILES = ["XYZIN"]
    WHATNEXT = []

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction="inputData")

        self.createLine(["subtitle", "Input model"])
        self.openSubFrame(frame=True)
        self.createLine(["widget", "XYZIN"])
        self.getWidget("XYZIN").showAtomSelection()
        self.closeSubFrame()
