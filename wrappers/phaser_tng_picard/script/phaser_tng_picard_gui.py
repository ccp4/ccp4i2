from qtgui.CCP4TaskWidget import CTaskWidget


class phaser_tng_picard_gui(CTaskWidget):
    TASKNAME = "phaser_tng_picard"
    TASKTITLE = "Molecular replacement with Phaser TNG Picard"
    SHORTTASKTITLE = "PhaserTNG"
    DESCRIPTION = "Automated MR with Phaser (The Next Generation) Picard pipeline"

    def drawContents(self):
        self.openFolder(folderFunction="inputData")
        self.createLine(["subtitle", "Observations"])
        self.createLine(["widget", "OBSIN"])
        self.createLine(["subtitle", "Free-R flag"])
        self.createLine(["widget", "FREERFLAG"])
        self.createLine(["subtitle", "Asymmetric unit contents"])
        self.createLine(["widget", "ASUIN"])
        self.createLine(["subtitle", "Input models"])
        self.createLine(["widget", "XYZIN_LIST"])
