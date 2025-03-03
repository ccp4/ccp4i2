import sys
import functools

from PyQt4 import sip
from PyQt4 import QtCore, QtGui

import icon_rc
from jiffyQt import TaskList, TreeModel
import SimpleTaskManager

class TaskListProxyModel(QtGui.QSortFilterProxyModel):
        def __init__(self,parent=None):
                QtGui.QSortFilterProxyModel.__init__(self,parent)
                self.filterString = ""
                self.filterMode = ()

        def setFilterMode(self,string):
             self.filterMode = string
             self.invalidateFilter()

        def setFilterString(self,string):
             self.filterString = string
             self.invalidateFilter()

        def hasAcceptedChildrenOrIsAccepted(self, item):

            theData = item.data(QtCore.Qt.DisplayRole)
            theDataKeywords = item.keywords()

            filterModeTasks = []
            for m in self.filterMode:
                filterModeTasks.extend(SimpleTaskManager.MODULE_DEFAULTS[m])

            acceptedFilterMode = True

            if len(filterModeTasks) > 0 and not str(item.name()) in filterModeTasks:
                acceptedFilterMode = False

            if acceptedFilterMode and (len(str(self.filterString).strip()) == 0 or str(self.filterString).lower() in theData.lower()):
                return True

            if acceptedFilterMode and (theDataKeywords is not None and any(str(self.filterString).lower() in s for s in theDataKeywords)):
                return True

            accepted = False
            for ic in range(item.childCount()):
                if self.hasAcceptedChildrenOrIsAccepted(item.child(ic)):
                    accepted = True
            return accepted

        def filterAcceptsRow(self, sourceRow, sourceParent):
                source_index = self.sourceModel().index(sourceRow, 0, sourceParent);

                ip = sourceParent.internalPointer()
                theData = source_index.internalPointer().data(QtCore.Qt.DisplayRole)

                return self.hasAcceptedChildrenOrIsAccepted(source_index.internalPointer())

class Window(QtGui.QWidget):
    def createTaskDummy(self,taskName,taskDesc):
        alertLabel = QtGui.QLabel("In a real program this would create the following task:<br/>"+taskName+"<br/>"+taskDesc)
        alertDialog = QtGui.QDialog()
        alertLayout = QtGui.QVBoxLayout()
        alertLayout.addWidget(alertLabel)
        alertDialog.setLayout(alertLayout)
        dbb = QtGui.QDialogButtonBox()
        okButton = dbb.addButton(QtGui.QDialogButtonBox.Ok)
        okButton.clicked.connect(alertDialog.accept)
        alertLayout.addWidget(dbb)
        alertDialog.exec_()

    def __init__(self,parent=None):
        QtGui.QWidget.__init__(self,parent)
        self.setLayout(QtGui.QVBoxLayout())
        self.modeLayout = QtGui.QHBoxLayout()
        modeLabel = QtGui.QLabel("Mode")

        modeLabel.setStyleSheet("QLabel { font-size: 24pt;}")
        self.modeLayout.addWidget(modeLabel)

        modeWidget = QtGui.QWidget()
        modeWidget.setLayout(self.modeLayout)
        self.layout().addWidget(modeWidget)

        refinementButton = QtGui.QRadioButton("Refinement")
        coreButton = QtGui.QRadioButton("Core")
        anyButton = QtGui.QRadioButton("Any")
        mrButton = QtGui.QRadioButton("Molecular\nreplacement")
        epButton = QtGui.QRadioButton("Experimental\nphasing")

        self.modeLayout.addWidget(mrButton)
        self.modeLayout.addWidget(epButton)
        self.modeLayout.addWidget(refinementButton)
        self.modeLayout.addWidget(coreButton)
        self.modeLayout.addWidget(anyButton)

        self.modeStack = QtGui.QStackedWidget()
        self.layout().addWidget(self.modeStack)

        self.modeStack.setMaximumHeight(60)

        self.refinementWidget = QtGui.QWidget()
        self.coreWidget = QtGui.QWidget()
        self.anyWidget = QtGui.QWidget()
        self.mrWidget = QtGui.QWidget()
        self.epWidget = QtGui.QWidget()

        self.modeStack.addWidget(self.refinementWidget)
        self.modeStack.addWidget(self.coreWidget)
        self.modeStack.addWidget(self.anyWidget)
        self.modeStack.addWidget(self.mrWidget)
        self.modeStack.addWidget(self.epWidget)

        refinementButton.clicked.connect(functools.partial(self.modeStack.setCurrentWidget,self.refinementWidget))
        coreButton.clicked.connect(functools.partial(self.modeStack.setCurrentWidget,self.coreWidget))
        anyButton.clicked.connect(functools.partial(self.modeStack.setCurrentWidget,self.anyWidget))
        mrButton.clicked.connect(functools.partial(self.modeStack.setCurrentWidget,self.mrWidget))
        epButton.clicked.connect(functools.partial(self.modeStack.setCurrentWidget,self.epWidget))

        self.modeStack.setCurrentWidget(self.anyWidget)
        anyButton.setChecked(True)
        
        self.refinementWidget.setLayout(QtGui.QHBoxLayout())
        self.coreWidget.setLayout(QtGui.QHBoxLayout())
        self.anyWidget.setLayout(QtGui.QHBoxLayout())
        self.mrWidget.setLayout(QtGui.QHBoxLayout())
        self.epWidget.setLayout(QtGui.QHBoxLayout())

        self.refinementWidget.layout().addWidget(QtGui.QLabel("Pipelines"))
        self.mrWidget.layout().addWidget(QtGui.QLabel("Pipelines"))
        self.epWidget.layout().addWidget(QtGui.QLabel("Pipelines"))

        crankButton = QtGui.QPushButton("Crank2")
        phaserEPButton = QtGui.QPushButton("Phaser")

        self.epWidget.layout().addWidget(crankButton)
        self.epWidget.layout().addWidget(phaserEPButton)

        mrBumpButton = QtGui.QPushButton("MrBUMP")
        self.mrWidget.layout().addWidget(mrBumpButton)

        refmacButton = QtGui.QPushButton("Refmac5")
        self.refinementWidget.layout().addWidget(refmacButton)

        crankButton.clicked.connect(functools.partial(self.createTaskDummy,"crank2","<b>"+SimpleTaskManager.TASK_TITLES["crank2"]+"</b><br/>"+SimpleTaskManager.TASK_DESC["crank2"]))
        phaserEPButton.clicked.connect(functools.partial(self.createTaskDummy,"phaser_EP_AUTO","<b>"+SimpleTaskManager.TASK_TITLES["phaser_EP_AUTO"]+"</b><br/>"+SimpleTaskManager.TASK_DESC["phaser_EP_AUTO"]))
        mrBumpButton.clicked.connect(functools.partial(self.createTaskDummy,"mrbump_basic","<b>"+SimpleTaskManager.TASK_TITLES["mrbump_basic"]+"</b><br/>"+SimpleTaskManager.TASK_DESC["mrbump_basic"]))
        refmacButton.clicked.connect(functools.partial(self.createTaskDummy,"prosmart_refmac","<b>"+SimpleTaskManager.TASK_TITLES["prosmart_refmac"]+"</b><br/>"+SimpleTaskManager.TASK_DESC["prosmart_refmac"]))

        searchLayout =  QtGui.QHBoxLayout()
        searchLayout.addWidget(QtGui.QLabel("Filter"))

        searchBox = QtGui.QLineEdit()
        searchLayout.addWidget(searchBox)
        self.layout().addLayout(searchLayout)

        taskList = TaskList()
        self.layout().addWidget(taskList)

        self.oldModel = taskList.model()

        pm2 = TaskListProxyModel()
        pm2.setSourceModel(self.oldModel)
        taskList.setModel(pm2)
        def clicky(idx):
            taskName = taskList.model().data(idx,QtCore.Qt.UserRole + 2)
            self.createTaskDummy(taskName,taskList.model().data(idx))
        taskList.doubleClicked.connect(clicky)

        searchBox.textChanged.connect(pm2.setFilterString)

        refinementButton.clicked.connect(functools.partial(pm2.setFilterMode,("refinement",)))
        coreButton.clicked.connect(functools.partial(pm2.setFilterMode,("core",)))
        anyButton.clicked.connect(functools.partial(pm2.setFilterMode,()))
        mrButton.clicked.connect(functools.partial(pm2.setFilterMode,("molecular_replacement",)))
        epButton.clicked.connect(functools.partial(pm2.setFilterMode,("expt_phasing",)))

sip.setdestroyonexit(False)
app = QtGui.QApplication(sys.argv)

win = Window()
win.show()
win.raise_()

sys.exit(app.exec_())
