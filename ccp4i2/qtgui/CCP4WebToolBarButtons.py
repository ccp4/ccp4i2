from __future__ import print_function
import sys, os
import functools

from PySide2 import QtCore, QtGui, QtWidgets, QtWebEngine, QtWebEngineWidgets, QtWebChannel

class CCP4WebToolBarButtonBridge(QtCore.QObject):
    buttonClicked = QtCore.Signal(str)
    buttonList = QtCore.Signal(list)
    @QtCore.Slot(str)
    def clicked(self,buttonName):
        self.buttonClicked.emit(buttonName)
    @QtCore.Slot('QJsonArray')
    def namesRequest(self,buttonNames):
        ls = []
        for i in range(buttonNames.count()):
            ls.append(buttonNames.at(i).toString())
        self.buttonList.emit(ls)

class CCP4WebToolBarButtons(QtWebEngineWidgets.QWebEngineView):
    buttonClicked = QtCore.Signal(str)
    buttonList = QtCore.Signal(list)
    editVisible = QtCore.Signal()

    def contextMenuEvent(self,e):
        self.menu = QtWidgets.QMenu()
        act = QtWidgets.QAction("Customize",self)
        self.menu.addAction(act)
        self.menu.popup(self.mapToGlobal(e.pos()))
        act.triggered.connect(self.editVisible.emit)

    def __init__(self,parent=None,fileName=None):
        super(CCP4WebToolBarButtons, self).__init__()

        self.channel = QtWebChannel.QWebChannel()
        self.handler = CCP4WebToolBarButtonBridge()
        self.channel.registerObject('handler', self.handler)
        self.page().setWebChannel(self.channel)

        print(os.path.abspath(fileName))
        local_url = QtCore.QUrl.fromLocalFile(os.path.abspath(fileName))

        self.load(local_url)
        self.handler.buttonClicked.connect(self.buttonClicked.emit)
        self.handler.buttonList.connect(self.buttonList.emit)

    @QtCore.Slot(str,bool)
    def setButtonEnabled(self,button,state):
        if state:
            self.page().runJavaScript("setButtonEnabled('"+button+"',true)")
        else:
            self.page().runJavaScript("setButtonEnabled('"+button+"',false)")

    @QtCore.Slot(str,bool)
    def setButtonVisible(self,button,state):
        if state:
            self.page().runJavaScript("setButtonVisible('"+button+"',true)")
        else:
            self.page().runJavaScript("setButtonVisible('"+button+"',false)")


if __name__== "__main__":

    app = QtWidgets.QApplication(sys.argv)

    webview = CCP4WebToolBarButtons(fileName=sys.argv[1])
    win = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(webview)
    win.setLayout(layout)

    bigButtonLayout = QtWidgets.QHBoxLayout()

    cbTaskMenu = QtWidgets.QCheckBox("Task menu")
    cbExportProject = QtWidgets.QCheckBox("Export project")
    cbRun = QtWidgets.QCheckBox("Run")
    cbRunOnServer = QtWidgets.QCheckBox("Run on server")
    cbClone = QtWidgets.QCheckBox("Clone")
    cbHelp = QtWidgets.QCheckBox("Help")
    cbBibliography = QtWidgets.QCheckBox("Bibliography")
    cbExportMTZ = QtWidgets.QCheckBox("Export MTZ")
    cbShowLogFile = QtWidgets.QCheckBox("Show log file")
    cbShowI2run = QtGui.QCheckBox("Show equivalent i2run command")
    cbViewCoot = QtWidgets.QCheckBox("View Coot")
    cbViewCCP4MG = QtWidgets.QCheckBox("View in CCP4MG")

    buttonLayout = QtWidgets.QVBoxLayout()
    buttonEnableLabel = QtWidgets.QLabel("Enable/disable buttons")
    buttonLayout.addWidget(buttonEnableLabel)
    buttonLayout.addWidget(cbTaskMenu)
    buttonLayout.addWidget(cbExportProject)
    buttonLayout.addWidget(cbRun)
    buttonLayout.addWidget(cbRunOnServer)
    buttonLayout.addWidget(cbClone)
    buttonLayout.addWidget(cbHelp)
    buttonLayout.addWidget(cbBibliography)
    buttonLayout.addWidget(cbExportMTZ)
    buttonLayout.addWidget(cbShowLogFile)
    buttonLayout.addWidget(cbShowI2run)
    buttonLayout.addWidget(cbViewCoot)
    buttonLayout.addWidget(cbViewCCP4MG)

    cbTaskMenu.setCheckState(QtCore.Qt.Checked)
    cbExportProject.setCheckState(QtCore.Qt.Checked)
    cbRun.setCheckState(QtCore.Qt.Checked)
    cbRunOnServer.setCheckState(QtCore.Qt.Checked)
    cbClone.setCheckState(QtCore.Qt.Checked)
    cbHelp.setCheckState(QtCore.Qt.Checked)
    cbBibliography.setCheckState(QtCore.Qt.Checked)
    cbExportMTZ.setCheckState(QtCore.Qt.Checked)
    cbShowLogFile.setCheckState(QtCore.Qt.Checked)
    cbViewCoot.setCheckState(QtCore.Qt.Checked)
    cbViewCCP4MG.setCheckState(QtCore.Qt.Checked)

    cbTaskMenu.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"task_menu"))
    cbExportProject.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"export_project"))
    cbRun.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"run"))
    cbRunOnServer.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"run_remote"))
    cbClone.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"clone"))
    cbHelp.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"task_help"))
    cbBibliography.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"references"))
    cbExportMTZ.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"export_mtz"))
    cbShowLogFile.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"show_log"))
    cbShowI2run.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"show_i2run"))
    cbViewCoot.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"view_coot"))
    cbViewCCP4MG.clicked[bool].connect(functools.partial(webview.setButtonEnabled,"view_ccp4mg"))

    bigButtonLayout.addLayout(buttonLayout)

    cbTaskMenu = QtWidgets.QCheckBox("Task menu")
    cbExportProject = QtWidgets.QCheckBox("Export project")
    cbRun = QtWidgets.QCheckBox("Run")
    cbRunOnServer = QtWidgets.QCheckBox("Run on server")
    cbClone = QtWidgets.QCheckBox("Clone")
    cbHelp = QtWidgets.QCheckBox("Help")
    cbBibliography = QtWidgets.QCheckBox("Bibliography")
    cbExportMTZ = QtWidgets.QCheckBox("Export MTZ")
    cbShowLogFile = QtWidgets.QCheckBox("Show log file")
    cbShowI2run = QtGui.QCheckBox("Show equivalent i2run command")
    cbViewCoot = QtWidgets.QCheckBox("View Coot")
    cbViewCCP4MG = QtWidgets.QCheckBox("View in CCP4MG")

    buttonLayout = QtWidgets.QVBoxLayout()
    buttonEnableLabel = QtWidgets.QLabel("Hide/show buttons")
    buttonLayout.addWidget(buttonEnableLabel)
    buttonLayout.addWidget(cbTaskMenu)
    buttonLayout.addWidget(cbExportProject)
    buttonLayout.addWidget(cbRun)
    buttonLayout.addWidget(cbRunOnServer)
    buttonLayout.addWidget(cbClone)
    buttonLayout.addWidget(cbHelp)
    buttonLayout.addWidget(cbBibliography)
    buttonLayout.addWidget(cbExportMTZ)
    buttonLayout.addWidget(cbShowLogFile)
    buttonLayout.addWidget(cbShowI2run)
    buttonLayout.addWidget(cbViewCoot)
    buttonLayout.addWidget(cbViewCCP4MG)

    cbTaskMenu.setCheckState(QtCore.Qt.Checked)
    cbExportProject.setCheckState(QtCore.Qt.Checked)
    cbRun.setCheckState(QtCore.Qt.Checked)
    cbRunOnServer.setCheckState(QtCore.Qt.Checked)
    cbClone.setCheckState(QtCore.Qt.Checked)
    cbHelp.setCheckState(QtCore.Qt.Checked)
    cbBibliography.setCheckState(QtCore.Qt.Checked)
    cbExportMTZ.setCheckState(QtCore.Qt.Checked)
    cbShowLogFile.setCheckState(QtCore.Qt.Checked)
    cbShowI2run.setCheckState(QtCore.Qt.Checked)
    cbViewCoot.setCheckState(QtCore.Qt.Unchecked)
    cbViewCCP4MG.setCheckState(QtCore.Qt.Unchecked)

    cbTaskMenu.clicked[bool].connect(functools.partial(webview.setButtonVisible,"task_menu"))
    cbExportProject.clicked[bool].connect(functools.partial(webview.setButtonVisible,"export_project"))
    cbRun.clicked[bool].connect(functools.partial(webview.setButtonVisible,"run"))
    cbRunOnServer.clicked[bool].connect(functools.partial(webview.setButtonVisible,"run_remote"))
    cbClone.clicked[bool].connect(functools.partial(webview.setButtonVisible,"clone"))
    cbHelp.clicked[bool].connect(functools.partial(webview.setButtonVisible,"task_help"))
    cbBibliography.clicked[bool].connect(functools.partial(webview.setButtonVisible,"references"))
    cbExportMTZ.clicked[bool].connect(functools.partial(webview.setButtonVisible,"export_mtz"))
    cbShowLogFile.clicked[bool].connect(functools.partial(webview.setButtonVisible,"show_log"))
    cbShowI2run.clicked[bool].connect(functools.partial(webview.setButtonVisible,"show_i2run"))
    cbViewCoot.clicked[bool].connect(functools.partial(webview.setButtonVisible,"view_coot"))
    cbViewCCP4MG.clicked[bool].connect(functools.partial(webview.setButtonVisible,"view_ccp4mg"))
    
    bigButtonLayout.addLayout(buttonLayout)

    layout.addLayout(bigButtonLayout)

    win.show()
    win.raise_()

    def logger(name):
        print(name)
    webview.buttonClicked.connect(logger)

    def setCootMGInivisible():
        webview.setButtonVisible("view_coot",False)
        webview.setButtonVisible("view_ccp4mg",False)

    webview.loadFinished.connect(setCootMGInivisible)

    sys.exit(app.exec_())
