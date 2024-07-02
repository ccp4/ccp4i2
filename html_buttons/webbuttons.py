from __future__ import print_function
import sys, os
import functools


from PySide6 import QtWidgets, QtCore, QtWebEngineWidgets, QtWebChannel

class CCP4WebToolBarButtonBridge(QtCore.QObject):

    buttonClicked = QtCore.Signal(str)

    @QtCore.Slot(str)
    def clicked(self, buttonName):
      self.buttonClicked.emit(buttonName)

class CCP4WebToolBarButtons(QtWebEngineWidgets.QWebEngineView):

    buttonClicked = QtCore.Signal(str)

    def __init__(self):
        super(CCP4WebToolBarButtons, self).__init__()

        self.channel = QtWebChannel.QWebChannel()
        self.handler = CCP4WebToolBarButtonBridge()
        self.channel.registerObject('handler', self.handler)
        self.page().setWebChannel(self.channel)

        file_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "qtgui", "html_i2_buttons.html"))
        local_url = QtCore.QUrl.fromLocalFile(file_path)

        self.load(local_url)
        self.handler.buttonClicked.connect(self.buttonClicked.emit)

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

    webview = CCP4WebToolBarButtons()

    win = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(webview)
    win.setLayout(layout)

    bigButtonLayout = QtWidgets.QHBoxLayout()
    buttonLayout = QtWidgets.QVBoxLayout()

    cbTaskMenu = QtWidgets.QCheckBox("Task menu")
    cbExportProject = QtWidgets.QCheckBox("Export project")
    cbRun = QtWidgets.QCheckBox("Run")
    cbRunOnServer = QtWidgets.QCheckBox("Run on server")
    cbClone = QtWidgets.QCheckBox("Clone")
    cbHelp = QtWidgets.QCheckBox("Help")
    cbBibliography = QtWidgets.QCheckBox("Bibliography")
    cbExportMTZ = QtWidgets.QCheckBox("Export MTZ")
    cbShowLogFile = QtWidgets.QCheckBox("Show log file")
    cbShowI2run = QtWidgets.QCheckBox("Show equivalent i2run command")
    cbViewCoot = QtWidgets.QCheckBox("View Coot")
    cbViewCCP4MG = QtWidgets.QCheckBox("View in CCP4MG")

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
    cbShowI2run.setCheckState(QtCore.Qt.Checked)
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
    cbShowI2run = QtWidgets.QCheckBox("Show equivalent i2run command")
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
