from __future__ import print_function

import sys
import os
import functools
import datetime
import time
import collections
import queue
import re
from lxml import etree

import sqlite3
from PyQt4 import QtGui, QtCore, QtWebKit
from lxml import etree

import SimpleTaskManager
import CLiteOneOffs

def pretty(d, indent=0):
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))

class TreeItem(QtGui.QStandardItem):
    def __init__(self,data,parentItem=None):
        QtGui.QStandardItem. __init__(self)
        self.m_parentItem = parentItem
        self.m_itemData = data
        self.m_childItems = []
        self._icon = None
        self._itemType = 0
        self._name = None
        self._number = None
        self._date = None
        self._keywords = None

    def setName(self,name):
        self._name = name

    def name(self):
        return self._name

    def keywords(self):
        return self._keywords

    def setKeywords(self,keywords):
        self._keywords = keywords

    def setDate(self,date):
        self._date = date

    def setNumber(self,number):
        self._number = number

    def parentItem(self):
        return self.m_parentItem

    def appendChild(self,child):
        self.m_childItems.append(child)

    def child(self,row):
        return self.m_childItems[row]

    def childCount(self):
        return len(self.m_childItems)

    def row(self):
        if self.m_parentItem is not None:
            return self.m_parentItem.m_childItems.index(self);
        return 0

    def columnCount(self):
        return 1#len(self.m_itemData)

    def data(self,column):
        return self.m_itemData#[column]

    def setIcon(self,icon):
        self._icon = icon

    def icon(self):
        return self._icon

class JobTreeItem(TreeItem):
    def __init__(self,data,parentItem=None):
        TreeItem. __init__(self,data,parentItem)

    def columnCount(self):
        return 2#len(self.m_itemData)

    def data(self,column):
        if column == 0:
            return "{} {}".format(self.m_itemData[0], self.m_itemData[1])
        elif column == 1:
            return self.m_itemData[2]
        return self.m_itemData[column]

class TreeModel(QtCore.QAbstractItemModel):

    def __init__(self,parent=None):
        QtCore.QAbstractItemModel.__init__(self,parent)

        self.rootItem = TreeItem("Root")

    def headerData(self, section, orientation, role = QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if section == 0:
                return "Title"
            if section == 1:
                return "Date"

        return QtCore.QVariant()

    def flags(self,index):
        if not index.isValid():
            return 0

        return QtCore.QAbstractItemModel.flags(self, index)

    def data(self, index, role):
        #FIXME - Maybe I should have columns instead of all these roles? Maybe I do already?
        if not index.isValid():
            return QtCore.QVariant()

        item = index.internalPointer()

        if role == QtCore.Qt.UserRole + 5:
            return item._keywords

        if role == QtCore.Qt.UserRole + 4:
            return item._date

        if role == QtCore.Qt.UserRole + 3:
            return item._number

        if role == QtCore.Qt.UserRole + 2:
            return item._name

        if role == QtCore.Qt.UserRole + 1:
            return item._itemType

        if role == QtCore.Qt.DecorationRole:
            return item.icon()

        if role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant();

        return item.data(index.column())

    def columnCount(self,parent = QtCore.QModelIndex()):
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def rowCount(self,parent = QtCore.QModelIndex()):
        if (parent.column() > 0):
            return 0

        if not parent.isValid():
            parentItem = self.rootItem;
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount();

    def parent(self,index):
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parentItem()

        if parentItem is self.rootItem:
            return QtCore.QModelIndex()

        if parentItem is not None:
            return self.createIndex(parentItem.row(), 0, parentItem)

    def index(self,row,column,parent=QtCore.QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem is not None:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

class JobListTreeModel(TreeModel):
    def __init__(self,data,parent=None):
        TreeModel.__init__(self,parent)

        self.rootItem = JobTreeItem("Root")
        self.setupModelData(data,self.rootItem)

    def setupModelData(self,lines,parent):
        lines = [x for x in lines if len(x)>0]

        jobTreeItems = []
        for l in lines:
            jobTreeItems.append(JobTreeItem(l, None))
            jobTreeItems[-1].setDate(l[2])
            jobTreeItems[-1].setNumber(l[0])
            jobTreeItems[-1].setName(l[1])

        for job in jobTreeItems:
            childJobs = [candidate for candidate in jobTreeItems if candidate.m_itemData[4] == job.m_itemData[3]]
            for childJob in childJobs:
                childJob.m_parentItem = job
                job.appendChild(childJob)

        unrootedJobs = [job for job in jobTreeItems if job.m_parentItem is None]
        for job in unrootedJobs:
            job.m_parentItem = parent
            parent.appendChild(job)

        return

class ScrollEater(QtCore.QObject):

    def __init__(self,parent=None):
        QtCore.QObject.__init__(self,parent)
        self.receiver = None

    def setScrollReceiver(self,receiver):
        self.receiver = receiver

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.Wheel:
            #Eat wheel events
            if self.receiver is not None:
                QtGui.QApplication.sendEvent(self.receiver,event)
            return True
        else:
            #standard event processing
            return QtCore.QObject.eventFilter(self, obj, event)

class FileSelector(QtGui.QWidget):
    def __init__(self,parent=None,text=""):
        QtGui.QWidget. __init__(self,parent)
        self.combo = QtGui.QComboBox()
        self.combo.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Maximum)
        layout = QtGui.QHBoxLayout()
        layout.setContentsMargins(3,3,3,3)
        self.setLayout(layout)
        if text:
            print(text)
            text_label = QtGui.QLabel(text)
            text_label.setToolTip(text)
            text_label.setFixedWidth(150)
            layout.addWidget(text_label)
        layout.addWidget(self.combo)
        #self.combo.setStretch(2)
        fileButton = QtGui.QToolButton()
        layout.addWidget(fileButton)
        file_icon = QtGui.QIcon(":/qticons/"+"file_manager.png")
        fileButton.setIcon(file_icon)
        dbButton = QtGui.QToolButton()
        layout.addWidget(dbButton)
        db_icon = QtGui.QIcon(":/qticons/"+"database.png")
        dbButton.setIcon(db_icon)

        def getFile():
            fileName = QtGui.QFileDialog.getOpenFileName()
        fileButton.clicked.connect(getFile)

class HTMLDelegate(QtGui.QStyledItemDelegate):

    def paint(self, painter, option, index):
        options = QtGui.QStyleOptionViewItemV4(option)
        self.initStyleOption(options,index)

        style = QtGui.QApplication.style() if options.widget is None else options.widget.style()

        doc = QtGui.QTextDocument()
        doc.setHtml(options.text)

        options.text = ""
        style.drawControl(QtGui.QStyle.CE_ItemViewItem, options, painter);

        ctx = QtGui.QAbstractTextDocumentLayout.PaintContext()

        # Highlighting text if item is selected
        if options.state & QtGui.QStyle.State_Selected:
            ctx.palette.setColor(QtGui.QPalette.Text, options.palette.color(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText));

        textRect = style.subElementRect(QtGui.QStyle.SE_ItemViewItemText, options)
        painter.save()
        painter.translate(textRect.topLeft())
        painter.setClipRect(textRect.translated(-textRect.topLeft()))
        doc.documentLayout().draw(painter, ctx)

        painter.restore()

    def sizeHint(self, option, index):
        options = QtGui.QStyleOptionViewItemV4(option)
        self.initStyleOption(options,index)

        doc = QtGui.QTextDocument()
        doc.setHtml(options.text)
        doc.setTextWidth(options.rect.width())
        itemType = index.data(QtCore.Qt.UserRole + 1)
        if itemType == 1:
            return QtCore.QSize(doc.idealWidth(), 48)
        return QtCore.QSize(doc.idealWidth(), 32)

class TaskList(QtGui.QTreeView):
    def __init__(self,parent=None):
        QtGui.QTreeView.__init__(self,parent)
        self.setIconSize(QtCore.QSize(32,32))
        self.taskListModel = TreeModel()
        self.setModel(self.taskListModel)
        self.setHeaderHidden(True)
        rootItem = self.taskListModel.rootItem
        icon = QtGui.QIcon(":/qticons/"+"ccp4.png")
        for mod in SimpleTaskManager.MODULE_ORDER:
            moduleItem  = TreeItem("<b>"+SimpleTaskManager.MODULE_TITLES[mod]+"</b>",rootItem)
            modIcon =  QtGui.QIcon(":/qticons/"+SimpleTaskManager.MODULE_ICONS[mod])
            moduleItem.setIcon(modIcon)
            rootItem.appendChild(moduleItem)
            for task in SimpleTaskManager.MODULE_DEFAULTS[mod]:
                taskIcon =  QtGui.QIcon(":/qticons/"+SimpleTaskManager.TASK_ICONS[task])
                if task in SimpleTaskManager.TASK_TITLES:
                    if task in SimpleTaskManager.TASK_DESC:
                        taskItem  = TreeItem("<b>"+SimpleTaskManager.TASK_TITLES[task]+"</b><br/><i>"+SimpleTaskManager.TASK_DESC[task]+"</i>",moduleItem)
                    else:
                        taskItem  = TreeItem("<b>"+SimpleTaskManager.TASK_TITLES[task]+"</b><br/>",moduleItem)
                else:
                    print("No title for",task)
                    taskItem  = TreeItem("<b>"+task+"</b><br/>",moduleItem)
                if task in SimpleTaskManager.TASK_KEYWORDS:
                    taskItem.setKeywords(SimpleTaskManager.TASK_KEYWORDS[task])
                taskItem.setIcon(taskIcon)
                taskItem._itemType = 1
                taskItem.setName(task)
                moduleItem.appendChild(taskItem)

        delegate = HTMLDelegate()
        self.setItemDelegate(delegate)

class WindowContainer(object):
    def __init__(self):
        self.window_list = []

class SimpleDBBrowser(QtGui.QMainWindow):
    def setDBFile(self,fn):
        self.dbfile = fn

    def __init__(self,parent=None,container=None):
        QtGui.QMainWindow.__init__(self,parent)

        self.projModel = QtGui.QStringListModel()
        self.jobModel = JobListTreeModel([])#QtGui.QStringListModel()

        widget = QtGui.QWidget()

        splitter = QtGui.QSplitter()
        self.setCentralWidget(splitter)

        self.projList = QtGui.QTreeView()
        self.projList.setModel(self.projModel)
        self.projList.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)

        self.jobList = QtGui.QTreeView()
        self.jobList.setModel(self.jobModel)
        self.jobList.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)

        self.projList.setHeaderHidden(True)
        #self.jobList.setHeaderHidden(True)

        webSettings = QtWebKit.QWebSettings.globalSettings()
        webSettings.setAttribute(QtWebKit.QWebSettings.DeveloperExtrasEnabled,True)
        webSettings.setAttribute(QtWebKit.QWebSettings.LocalStorageEnabled, True);
        webSettings.setFontSize(QtWebKit.QWebSettings.DefaultFontSize,12)
        webSettings.setFontSize(QtWebKit.QWebSettings.DefaultFixedFontSize,12)
        #webSettings.setAttribute(QtWebKit.QWebSettings.SessionStorageEnabled, True);
        self.jobInfoWidget = QtWebKit.QWebView()
        self.jobInfoWidget.setHtml("<h1>A report will appear here</h1>")

        splitter.addWidget(self.projList)
        self.projList.hide()
        splitter.addWidget(self.jobList)

        self.stack = QtGui.QStackedWidget()
        self.stack.addWidget(self.jobInfoWidget)
        splitter.addWidget(self.stack)

        self.jobInfo = {}

        self.projList.selectionModel().selectionChanged.connect(self.selectedProject)
        self.jobList.selectionModel().selectionChanged.connect(self.selectedJob)

        self.fileMenu = self.menuBar().addMenu("File")
        self.projectMenu = self.menuBar().addMenu("Projects")

        self.selectedProjectName = None

        self.window_container = container
        if container:
            container.window_list.append(self)

        closeAction = self.fileMenu.addAction("Close")
        closeAction.setShortcut(QtGui.QKeySequence.Close)
        closeAction.triggered.connect(self.close)
        self.fileMenu.addSeparator()

        self.taskList = TaskList()
        self.stack.addWidget(self.taskList)
        taskListAction = self.fileMenu.addAction("Task List")

        taskListAction.triggered.connect(self.jobList.selectionModel().clearSelection)
        taskListAction.triggered.connect(functools.partial(self.stack.setCurrentWidget,self.taskList))

        self.taskList.doubleClicked.connect(self.newJob)

        self.taskWidget = QtGui.QLabel()
        self.stack.addWidget(self.taskWidget)

        self.stack.setCurrentWidget(self.taskList)

        toolBar = self.addToolBar("toolbar")
        toolBar.addAction(taskListAction)
        taskListAction.setIcon(QtGui.QIcon(":/qticons/taskmenu.svg"))
        toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)

        self.files = []
        self.responseQueue = queue.Queue()

    def newJob(self, idx):
        taskName = self.taskList.model().data(idx,QtCore.Qt.UserRole + 2)
        print(taskName)
        projectId = self.projIds[self.selectedProjectName]

        from core.CCP4Modules import PROJECTSMANAGER
        pm = PROJECTSMANAGER()
        print(pm.db())

        from dbapi.CCP4DbUtils import COpenJob
        cOpenJob = COpenJob(projectId=projectId)
        cOpenJob.createJob(taskName=taskName)
        #THis should trigger a "jobCreated signal" from the CDbApi
        self.newJobSMcN(idx)

    def jobCreated(self, createdJobDict):
        if createdJobDict['projectId'] != self.projIds[self.selectedProjectName]: return
        print("job has been created", createdJobDict)
        response = self.dbQuery('getProjectsAndJobs')
        if response != 'FailedDbInteraction':
            self.setInfo(response['projSums'], response['projJobs'], response['projDirs'], response['projIds'])

        thisJobInfo = self.jobInfo[self.selectedProjectName]
        textJobList = []
        for job in thisJobInfo:
            t = datetime.datetime.fromtimestamp(job[2])
            if job[7] is not None:
                textJobList.append((job[1],job[7],t.ctime(),job[0],job[11]))
            else:
                if job[9] in SimpleTaskManager.TASK_TITLES:
                    textJobList.append((job[1],SimpleTaskManager.TASK_TITLES[job[9]],t.ctime(),job[0],job[11]))
                else:
                    textJobList.append((job[1],job[9],t.ctime(),job[0],job[11]))

        textJobList.sort(key=lambda s: list(map(int, s[0].split('.'))))

        #This does nice nesting, but gets wrong report because indexing is now wrong! ....
        mod = JobListTreeModel(textJobList)
        self.jobList.setModel(mod)
        self.jobModel = mod
        self.jobList.selectionModel().selectionChanged.connect(self.selectedJob)

    def newJobSMcN(self,idx):
        #FIXME - Most inputs not implemented!
        #FIXME - JobListTreeModel needs more info, icluding job name.
        if self.taskList.model().data(idx,QtCore.Qt.UserRole + 1) == 1:
            if len(self.projList.selectionModel().selectedIndexes()) > 0:
                projId = str(self.projModel.data(self.projList.selectionModel().selectedIndexes()[0],QtCore.Qt.DisplayRole))
                projId = self.projIds[projId]
                print(projId)
                if len(self.files) == 0:
                    t1 = time.time()
                    self.files = self.dbQuery("getProjectFiles?projId={}".format(projId))
                    t2 = time.time()
                    print("Time to get",len(self.files),"files:",t2-t1)

            taskName = self.taskList.model().data(idx,QtCore.Qt.UserRole + 2)
            self.stack.setCurrentWidget(self.taskWidget)
            self.taskWidget.setText(taskName)
            """
            from ccp4i2.qtgui.CCP4ProjectViewer import CTaskInputFrame
            win = QtGui.QWidget()
            win2 = CTaskInputFrame(win)
            win2.createTaskWidget(taskName)
            """

            """
            #Let's not even bother with this, let's just ude def.xml"
            from ccp4i2.wrappers.ccp4mg_general.script import ccp4mg_general_gui
            win = ccp4mg_general_gui.Cccp4mg_general()
            win.draw()
            """

            print(taskName)
            xmlfn = SimpleTaskManager.get_def_xml_for_task(taskName)
            if xmlfn is None:
                return
            print(xmlfn)
            xmlf = open(xmlfn)
            xml = xmlf.read()
            xmlf.close()

            scrollWin = QtGui.QScrollArea()
            scrollEater = ScrollEater(scrollWin)
            scrollEater.setScrollReceiver(scrollWin)
            win = QtGui.QWidget()
            layout = QtGui.QVBoxLayout()
            layout.setContentsMargins(3,3,3,3)
            layout.setSpacing(2)
            win.setLayout(layout)

            headWidget = QtGui.QFrame()
            headLayout = QtGui.QVBoxLayout()
            headWidget.setLayout(headLayout)
            headLayout.setContentsMargins(2,2,2,2)
            headLayout.setSpacing(0)

            qticonsDir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(SimpleTaskManager.__file__))),"qticons")
            print('qticonsDir',qticonsDir)
            headWidget.setObjectName('jobHeaderFrame')
            headWidget.setStyleSheet('''QFrame#jobHeaderFrame { background-image: url("'''+qticonsDir+'''/backgroundJobTitle.png");; padding-left: 10px; padding-bottom: 5px; padding-top: 4px; padding-right: 10px; border: 1px solid; border-radius: 0px; border-color: #ec9; margin-bottom: 8px;}''')

            titleLayout = QtGui.QHBoxLayout()
            titleLayout.setContentsMargins(2,2,2,2)
            titleLabel = QtGui.QLabel("Job title")
            titleLayout.addWidget(titleLabel)
            titleLineEdit = QtGui.QLineEdit()
            titleLineEdit.setText(SimpleTaskManager.TASK_TITLES[taskName])
            titleLayout.addWidget(titleLineEdit)
            headLayout.addLayout(titleLayout)

            jobDataLayout = QtGui.QHBoxLayout()
            jobDataLayout.setContentsMargins(2,2,2,2)
            jobDataLabel = QtGui.QLabel("Use data from jobxx:")
            jobDataLayout.addWidget(jobDataLabel)
            jobDataCombo = QtGui.QComboBox()
            jobDataLayout.addWidget(jobDataCombo)
            jobDataLayout.addStretch(20)
            headLayout.addLayout(jobDataLayout)

            #FIXME - There's probably a better way to do this, by setting up a model for jobDataCombo and also having joblist model return just job number as an option.
            jobList = []
            for i in range (self.jobModel.rowCount()):
                jobList.append(str(self.jobModel.data(self.jobModel.index(i,0),QtCore.Qt.DisplayRole)).split()[0]+" " +self.jobModel.data(self.jobModel.index(i,0),QtCore.Qt.UserRole + 2))
            jobDataCombo.addItems(jobList)

            layout.addWidget(headWidget)

            parser = ET.XMLParser()
            tree = ET.fromstring(xml, parser)
            #print xml
            #print [x[3] for x in self.files if x[3] is not None]

            pdb_labels_pair = [(x[3],x[2]) for x in self.files if x[3] is not None and x[5] == 2]
            pdb_model = QtGui.QStandardItemModel()
            for i,x in zip(list(range(len(pdb_labels_pair))),pdb_labels_pair):
                item = QtGui.QStandardItem(x[0])
                item.setData(x[1],QtCore.Qt.UserRole + 1)
                pdb_model.setItem(i,0,item)

            map_labels_pair = [(x[3],x[2]) for x in self.files if x[3] is not None and x[5] in [10,11,12,13]]
            map_model = QtGui.QStandardItemModel()
            for i,x in zip(list(range(len(map_labels_pair))),map_labels_pair):
                item = QtGui.QStandardItem(x[0])
                item.setData(x[1],QtCore.Qt.UserRole + 1)
                map_model.setItem(i,0,item)

            def printUUID(fs,i):
                idx = fs.combo.model().index(i,0)
                print(fs.combo.model().data(idx), fs.combo.model().data(idx,QtCore.Qt.UserRole + 1))

            inputContainers =  tree.findall("ccp4i2_body/container[@id='inputData']")
            for container in inputContainers:
                contents = container.findall("content")
                for content in contents:
                    if len(content.findall("className")) == 0:
                        continue
                    className = content.findall("className")[0].text
                    #print className
                    text = content.attrib["id"]
                    if className == "CDictDataFile" or className == "CPdbDataFile" or className == "CMapCoeffsDataFile" or className == "CSeqDataFile":
                        if len(content.findall("qualifiers/label")) > 0:
                            text = content.findall("qualifiers/label")[0].text
                        fs = FileSelector(text=text)
                        fs.combo.installEventFilter(scrollEater)
                        if className == "CPdbDataFile":
                            fs.combo.setModel(pdb_model)
                        elif className == "CMapCoeffsDataFile":
                            fs.combo.setModel(map_model)
                        fs.combo.currentIndexChanged[int].connect(functools.partial(printUUID,fs))
                        layout.addWidget(fs)
                    elif className == "CList":

                        subClassName = content.findall("subItem/className")[0].text
                        frame = QtGui.QFrame()
                        frame.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
                        frameLayout = QtGui.QVBoxLayout()
                        frameLayout.setContentsMargins(2,2,2,2)
                        frameLayout.setSpacing(2)
                        layout.addWidget(frame)
                        frame.setLayout(frameLayout)
                        listWidget = QtGui.QTreeView()
                        listWidget.setHeaderHidden(True)
                        listWidget.hide()
                        hideShow = QtGui.QPushButton("Show list")
                        plusMinus = QtGui.QWidget()
                        plusMinus.hide()
                        plusMinusLayout = QtGui.QHBoxLayout()
                        plusMinusLayout.setContentsMargins(2,2,2,2)
                        plusMinusLayout.setSpacing(2)
                        plusMinus.setLayout(plusMinusLayout)
                        plusButton = QtGui.QPushButton()
                        minusButton = QtGui.QPushButton()
                        plusMinusLayout.addWidget(plusButton)
                        plusMinusLayout.addWidget(minusButton)
                        plusMinusLayout.addStretch(20)
                        plus_icon = QtGui.QIcon(":/qticons/"+"list_add_grey.png")
                        minus_icon = QtGui.QIcon(":/qticons/"+"list_delete_grey.png")
                        plusButton.setIcon(plus_icon)
                        minusButton.setIcon(minus_icon)
                        def toggleList(listWidget,plusMinus,hideShow):
                            if listWidget.isVisible():
                                listWidget.hide()
                                plusMinus.hide()
                                hideShow.setText("Show list")
                            else:
                                listWidget.show()
                                plusMinus.show()
                                hideShow.setText("Hide list")
                        hideShow.clicked.connect(functools.partial(toggleList,listWidget,plusMinus,hideShow))
                        hideShowLayout = QtGui.QHBoxLayout()
                        hideShowLayout.setContentsMargins(0,0,0,0)
                        hideShowLayout.setSpacing(2)
                        hideShowLayout.addWidget(hideShow)
                        hideShowLayout.addStretch(20)
                        frameLayout.addLayout(hideShowLayout)
                        frameLayout.addWidget(listWidget)
                        frameLayout.addWidget(plusMinus)
                        listModel = QtGui.QStandardItemModel()
                        item = QtGui.QStandardItem("--")
                        item.setData("--",QtCore.Qt.UserRole + 1)
                        listModel.setItem(0,0,item)
                        listWidget.setModel(listModel)
                        def addToModel(listModel,listWidget,minusButton):
                            row = listModel.rowCount()
                            item = QtGui.QStandardItem("--")
                            item.setData("--",QtCore.Qt.UserRole + 1)
                            listModel.setItem(row,0,item)
                            idx = listModel.index(row,0)
                            listWidget.selectionModel().select(idx,QtGui.QItemSelectionModel.ClearAndSelect)
                            minusButton.setEnabled(True)
                        def removeFromModel(listModel,listWidget,minusButton):
                            if listModel.rowCount() < 2:
                                return
                            if len(listWidget.selectionModel().selectedIndexes()) < 1:
                                return
                            index = listWidget.selectionModel().selectedIndexes()[0]
                            print(index.row())
                            listModel.removeRows(index.row(),1)
                            print("Now select",max(index.row()-1,0))
                            listWidget.selectionModel().select(index.sibling(max(index.row()-1,0),0),QtGui.QItemSelectionModel.ClearAndSelect)
                            if listModel.rowCount() < 2:
                                minusButton.setEnabled(False)
                        plusButton.clicked.connect(functools.partial(addToModel,listModel,listWidget,minusButton))
                        minusButton.clicked.connect(functools.partial(removeFromModel,listModel,listWidget,minusButton))
                        minusButton.setEnabled(False)
                        idx = listModel.index(0,0)
                        def updataListModel(listModel,listWidget,fs,text):
                            idx = fs.combo.model().index(fs.combo.currentIndex(),0)
                            index = listWidget.selectionModel().selectedIndexes()[0]
                            listModel.setData(index,text)
                            listModel.setData(index,fs.combo.model().data(idx,QtCore.Qt.UserRole + 1),QtCore.Qt.UserRole + 1)

                        listWidget.selectionModel().select(idx,QtGui.QItemSelectionModel.ClearAndSelect)
                        if subClassName == "CDictDataFile" or subClassName == "CPdbDataFile" or subClassName == "CMapCoeffsDataFile" or subClassName == "CSeqDataFile":
                            if len(content.findall("qualifiers/label")) > 0:
                                text = content.findall("qualifiers/label")[0].text
                            if len(content.findall("subItem/qualifiers/label")) > 0:
                                text = content.findall("subItem/qualifiers/label")[0].text
                            fs = FileSelector(text=text)
                            fs.combo.installEventFilter(scrollEater)
                            if subClassName == "CPdbDataFile":
                                fs.combo.setModel(pdb_model)
                            elif subClassName == "CMapCoeffsDataFile":
                                fs.combo.setModel(map_model)
                            if fs.combo.count() > 0:
                                updataListModel(listModel,listWidget,fs,fs.combo.currentText())
                            fs.combo.currentIndexChanged[str].connect(functools.partial(updataListModel,listModel,listWidget,fs))
                            fs.combo.currentIndexChanged[int].connect(functools.partial(printUUID,fs))
                            frameLayout.addWidget(fs)

                        def listChanged(fs,listModel,listWidget,idx):
                            print("List changed")
                            if len(idx.indexes())<1:
                                return
                            print(listModel.data(idx.indexes()[0],QtCore.Qt.DisplayRole))
                            uuid = listModel.data(idx.indexes()[0],QtCore.Qt.UserRole + 1)
                            fs.combo.blockSignals(True)
                            items = fs.combo.model().findItems(listModel.data(idx.indexes()[0],QtCore.Qt.DisplayRole))
                            for it in items:
                                if uuid == it.data(QtCore.Qt.UserRole + 1):
                                    fs.combo.setCurrentIndex(it.index().row())
                                    break
                            fs.combo.blockSignals(False)
                        listWidget.selectionModel().selectionChanged.connect(functools.partial(listChanged,fs,listModel,listWidget))

            layout.addStretch(20)
            scrollWin.setWidget(win)
            self.stack.addWidget(scrollWin)
            self.stack.setCurrentWidget(scrollWin)
            scrollWin.setWidgetResizable(True)

    def setInfo(self,projInfo,jinfo,projDirs,projIds):
        self.jobInfo = jinfo
        self.projModel.setStringList(projInfo)
        self.projDirs = projDirs
        self.projIds = projIds
        self.projectMenu.clear()

        if len(projInfo) < 16:
            for proj in projInfo:
                act = self.projectMenu.addAction(proj)
        else:
            sorted_list = sorted(projInfo, key=lambda s: str(s).lower())
            adMenu = self.projectMenu.addMenu("a-d")
            ehMenu = self.projectMenu.addMenu("e-h")
            ilMenu = self.projectMenu.addMenu("i-l")
            mpMenu = self.projectMenu.addMenu("m-p")
            qtMenu = self.projectMenu.addMenu("q-t")
            uzMenu = self.projectMenu.addMenu("u-z")
            numMenu = self.projectMenu.addMenu("other")
            for proj in sorted_list:
                if str(proj).lower().startswith(("a","b","c","d")):
                    act = adMenu.addAction(proj)
                elif str(proj).lower().startswith(("e","f","g","h")):
                    act = ehMenu.addAction(proj)
                elif str(proj).lower().startswith(("i","j","k","l")):
                    act = ilMenu.addAction(proj)
                elif str(proj).lower().startswith(("m","n","o","p")):
                    act = mpMenu.addAction(proj)
                elif str(proj).lower().startswith(("q","r","s","t")):
                    act = qtMenu.addAction(proj)
                elif str(proj).lower().startswith(("u","v","w","x","y","z")):
                    act = uzMenu.addAction(proj)
                else:
                    act = numMenu.addAction(proj)

        self.projectMenu.triggered.connect(self.selectedProject2)

    def getFileListForJob(self,jobId):
        self.files = self.dbQuery('getFileListForJob?jobId={}'.format(jobId))

    def dbQuery(self, queryString):
        response = CLiteOneOffs.cLiteDbThread().performAction(action=queryString, responseQueue=self.responseQueue)
        if response != 'FailedDbInteraction':
            return response
        else:
            raise Exception("Failed Db Interaction","getFileListForJob")

    def selectedJob(self, idx):
        def getIndexOfTuple(l, index, value):
            for pos,t in enumerate(l):
                if t[index] == value:
                    return pos

            # Matches behavior of list.index
            raise ValueError("list.index(x): x not in list")

        projId = str(self.projModel.data(self.projList.selectionModel().selectedIndexes()[0],QtCore.Qt.DisplayRole))
        thisJobInfo = self.jobInfo[projId]
        if len(idx.indexes()) == 0: return

        searchText = self.jobList.model().data(idx.indexes()[0],QtCore.Qt.DisplayRole).split()[0]
        searchRow  = getIndexOfTuple(thisJobInfo,1,searchText)

        jobId = thisJobInfo[searchRow][0]

        files = self.getFileListForJob(jobId)
        t = ""
        for f in self.files:
            t += f[1]
            if f[2] is not None:
                t += f[2]
            t += '\n'

        dir_split = thisJobInfo[searchRow][1].split(".")
        job_dir = []
        for i in dir_split: job_dir.append( "job_"+i)
        job_dir = os.path.join(*job_dir)
        job_dir = os.path.join(self.projDirs[projId],"CCP4_JOBS",job_dir)

        self.stack.setCurrentWidget(self.jobInfoWidget)
        report_file = os.path.join(job_dir,"report.html")

        #Look for report, or try to make it if necessary
        if not os.path.exists(report_file):
            concurrencyMethod = "Threading"
            if concurrencyMethod ==  "Subprocess":
                import subprocess
                thisDir = os.path.dirname(os.path.dirname(os.path.abspath(SimpleTaskManager.__file__)))
                cQuarantinedActionsPath = os.path.join(thisDir, "qs", "CQuarantinedActions.py")
                print(thisDir, cQuarantinedActionsPath)
                a=subprocess.call(['ccp4-python', cQuarantinedActionsPath, 'remakeReport', jobId, 'Finished', CLiteOneOffs.cLiteDbThread().dbfile])
            elif concurrencyMethod == "Threading":
                import threading
                import CQuarantinedActions
                bubblyQueue = queue.Queue()
                self.reportThread=threading.Thread(target=CQuarantinedActions.remakeReport, args=(bubblyQueue, jobId, 'Finished', CLiteOneOffs.cLiteDbThread().dbfile, ))
                self.reportThread.start()
                print(bubblyQueue.get(True, 10))
            elif concurrencyMethod == "Multiprocess":
                import multiprocessing
                import CQuarantinedActions
                bubblyQueue = multiprocessing.Queue()
                a=multiprocessing.Process(target=CQuarantinedActions.remakeReport, args=(bubblyQueue, jobId, 'Finished', CLiteOneOffs.cLiteDbThread().dbfile, ))
                a.start()
                print(bubblyQueue.get(True, 10))
        #Handle existing report file
        if os.path.exists(report_file):
            html = patchedReport(report_file, jobId=jobId, files=self.files, jobs=self.jobInfo[projId])
            urlString = "http://127.0.0.1:{}".format(CLiteOneOffs.cLiteHTTPThread().port)
            baseURL = QtCore.QUrl(urlString)
            self.jobInfoWidget.setHtml(html, baseURL)
        else:
            self.jobInfoWidget.setHtml("<span>No report available</span>")

    def setSelectedRow(self,row):
        self.projList.selectionModel().select(self.projModel.index(0,0).sibling(row,0),QtGui.QItemSelectionModel.ClearAndSelect)

    def selectedProject2(self,projact):
        idx = self.projModel.match(self.projModel.index(0,0),QtCore.Qt.DisplayRole, projact.text())[0]
        self.projList.selectionModel().select(idx,QtGui.QItemSelectionModel.Deselect)
        self.projList.selectionModel().select(idx,QtGui.QItemSelectionModel.ClearAndSelect)

    def selectedProject(self,idx):
        if len(idx.indexes()) == 0: return
        projId = str(self.projModel.data(idx.indexes()[0],QtCore.Qt.DisplayRole))

        tlw = QtGui.QApplication.topLevelWidgets()
        for widget in tlw:
            if widget.objectName() == projId:
                widget.show()
                widget.raise_()
                widget.activateWindow()
                return
            simpleWidgets = widget.findChildren(SimpleDBBrowser)
            for child in simpleWidgets:
                if child.objectName() == projId:
                    child.show()
                    child.raise_()
                    child.activateWindow()
                    return

        if self.selectedProjectName is None:
            self.selectedProjectName = projId
            self.setObjectName(projId)

        elif projId != self.selectedProjectName:
            newWindow = SimpleDBBrowser(container=self.window_container)
            newWindow.show()
            newWindow.raise_()
            projInfo = list(self.projModel.stringList())
            jinfo = self.jobInfo
            projDirs = self.projDirs
            projIds = self.projIds
            newWindow.setInfo(projInfo,jinfo,projDirs,projIds)
            newWindow.setSelectedRow(idx.indexes()[0].row())
            newWindow.setObjectName(projId)
            if hasattr(self, "dbfile"): newWindow.dbfile = self.dbfile
            return

        self.setWindowTitle(projId)
        thisJobInfo = self.jobInfo[projId]
        textJobList = []
        for job in thisJobInfo:
            t = datetime.datetime.fromtimestamp(job[2])
            if job[7] is not None:
                textJobList.append((job[1],job[7],t.ctime(),job[0],job[11]))
            else:
                if job[9] in SimpleTaskManager.TASK_TITLES:
                    textJobList.append((job[1],SimpleTaskManager.TASK_TITLES[job[9]],t.ctime(),job[0],job[11]))
                else:
                    textJobList.append((job[1],job[9],t.ctime(),job[0],job[11]))

        textJobList.sort(key=lambda s: list(map(int, s[0].split('.'))))

        #This does nice nesting, but gets wrong report because indexing is now wrong! ....
        mod = JobListTreeModel(textJobList)
        self.jobList.setModel(mod)
        self.jobModel = mod
        self.jobList.selectionModel().selectionChanged.connect(self.selectedJob)
        """
        #This shows report properly ...
        self.jobModel.setStringList(textJobList)
        """

def patchedReport(reportPath, jobId=None, files=[], jobs=[]):
    try:
        #For some reason, report files seem to have a hardwired root address...THis Must go !
        reportString = open(reportPath).read()
        #
        listOfPortedAddresses = re.findall(r"http:\/\/127.0.0.1:[0-9]*", reportString)
        setOfPortedAddresses = set(listOfPortedAddresses)
        for portedAddress in setOfPortedAddresses:
            reportString = reportString.replace(portedAddress, "")
        #
        reportString = reportString.replace("http://127.0.0.1:None/report_files/","/report_files/")
        reportString = reportString.replace("None/report_files/","/report_files/")
        reportXML = ET.fromstring(reportString)
    except:
        print('Failed to parse report', reportPath)
        #Attempt to remake report
        remade = remakeReport(jobId)
        try:
            reportXML = ET.parse(reportPath)
        except:
            print('Failed to parse report', reportPath)
            return None

    #print 'reportXML parsed',jobId
    # Patch img elements
    for img in reportXML.findall('.//img'):
        try: img.set('src', "/jobId/{}/{}".format(jobId, img.get("src")))
        except: print('Failed patching img', ET.tostring(img))

    # Patch drawn spans
    for drawnSpan in reportXML.findall(".//span[@data-url]"):
        try:
            patchedElement ="/jobId/{}/{}".format(jobId, drawnSpan.get('data-url'))
            #patchedElement ="/ManageCCP4i2Archive/?File=True?jobId="+jobId+"?filePath="+drawnSpan.get('data-url')
            drawnSpan.set('data-url',patchedElement)
        except: print('Failed patching drawnSpan', ET.tostring(drawnSpan))
    #print 'data-is-urls patched'

    # Patch drawn divs
    for drawnDiv in reportXML.findall(".//div[@data-is-urls='True']"):
        try:
            patchedElement ="/jobId/{}/{}".format(jobId, drawnDiv.get('data-data'))
            #patchedElement ="/ManageCCP4i2Archive/?File=True?jobId="+jobId+"?filePath="+drawnDiv.get('data-data')
            drawnDiv.set('data-data',patchedElement)
        except: print('Failed patching drawnDiv', ET.tostring(drawnDiv))
    #print 'data-is-urls patched'

    #Patch downloadObjects
    for qt_object in reportXML.findall(".//object[@class='qt_object']"):
        try:
            newDiv = ET.Element('div')
            newDiv.set('class','qt_object')

            mimeTypeName = qt_object.findall('param[@name="mimeTypeName"]')[0].get("value")
            dbFileId = qt_object.findall('param[@name="dbFileId"]')[0].get("value")
            try: annotation = qt_object.findall('param[@name="annotation"]')[0].get("value")
            except:
                try: annotation = qt_object.findall('param[@name="baseName"]')[0].get("value")
                except: annotation = "Can't find annotation or baseName"
            imgNode = ET.SubElement(newDiv,"img",src="/icon/"+mimeTypeName, height="20", width="20", style="margin:5px;display:inline-block;border:1px solid grey;")

            fileSpec = "?File=True?fileId="+dbFileId
            clickDiv = ET.SubElement(newDiv, "div", style="width:500px;height:25px;float:right;font-size:125%;border:0px solid red;margin:2px;", id="Download_"+dbFileId, onclick = 'window.downloadFile("'+fileSpec+'");')
            clickDiv.set("class","downloadable")

            matchedFiles = [file for file in files if file[0] == dbFileId]
            jobNumber = "?"
            if len (matchedFiles) == 1:
                fileCreationJobId = matchedFiles[0][6]
                matchedJobs = [job for job in jobs if job[0] == fileCreationJobId]
                if len(matchedJobs) == 1:
                    jobNumber = matchedJobs[0][1]

            clickDiv.text = "Job: "+str(jobNumber) + " - "+annotation

            qt_object.getparent().replace(qt_object, newDiv)
        except Exception as err: print('Failed to patch qt_object',etree.tostring(qt_object), err)
    #print 'qt_objects patched'

    #Patch downloadObjects
    for qt_launch in reportXML.findall(".//object[@type='x-ccp4-widget/CDownloadButton']"):
        try:
            downloadsJobId = qt_launch.findall('param[@name="jobId"]')[0].get("value")
            dataName = qt_launch.findall('param[@name="dataName"]')[0].get("value")
            fileSpec = "/jobId/{}/tables_as_csv_files/{}.csv".format(downloadsJobId, dataName);
            newButton = ET.Element('Button', onclick = 'window.downloadFile("'+fileSpec+'");')
            newButton.text = 'Download'

            qt_launch.getparent().replace(qt_launch, newButton)
        except: print('Failed to patch qt_launch',etree.tostring(qt_launch))
    #print 'CDownloadButtons patched'

    #Patch CLauncherButton objects
    for qt_launch in reportXML.findall(".//object[@type='x-ccp4-widget/CLauncherButton']"):
        try:
            executable = qt_launch.findall('param[@name="exe"]')[0].get("value")
            if executable == "CCP4mg":
                sceneJobId = qt_launch.findall('param[@name="jobId"]')[0].get("value")
                sceneName = qt_launch.findall('param[@name="sceneFile"]')[0].get("value")
                #Remove .scene.xml
                pathElements = sceneName.split(".")
                noExtension = pathElements[-3]
                sceneNumber = noExtension.split("_")[-1]
                urlOfScene = "/ManageCCP4i2Archive/SceneForJob/"+sceneJobId+"/"+sceneNumber
                newButton = ET.Element('Button', onclick = 'window.open("'+urlOfScene+'");')
                newButton.text = '3D view'

                qt_launch.getparent().replace(qt_launch, newButton)
            else:
                newButton = ET.Element('Button', disabled='disabled')
                newButton.text = 'Not active'

                qt_launch.getparent().replace(qt_launch, newButton)
        except: print('Failed to patch qt_launch',etree.tostring(qt_launch))
    #print 'Launch MGs patched'

    return ET.tostring(reportXML)

class QSSupervisorApp(QtGui.QApplication):
    def __init__(self, *argv, **kw):
        t1 = time.time()
        import icon_rc

        modArgs = list(argv[0])
        if len(modArgs) < 2:
            from core import CCP4Config
            from core.CCP4Utils import getDotDirectory
            fileName = CCP4Config.DBFILE()
            if fileName is None:
                fileName = os.path.join(getDotDirectory(), 'db', 'database.sqlite')
            print(fileName)
            modArgs.append(fileName)

        from utils.startup import startProjectsManager
        pm = startProjectsManager(dbFileName=modArgs[1])

        super(QSSupervisorApp, self).__init__([modArgs], **kw)


        #Prepare to close HTTP and DBTHreads on exit
        self.aboutToQuit.connect(self.aboutToQuitCallback)

        window_container = WindowContainer()
        win = SimpleDBBrowser(container=window_container)
        win.show()
        win.raise_()

        self.flush()

        pm.db().jobCreated.connect(win.jobCreated)

        self.responseQueue = queue.Queue()
        cLiteDbThread = CLiteOneOffs.cLiteDbThread(parent=self, dbFile=modArgs[1])

        import CLiteHTTPThread
        thisDir = os.path.split(os.path.abspath(__file__))[0]
        parentDir = os.path.split(thisDir)[0]
        reportFilesPath = os.path.join(parentDir,"docs")

        cLiteHTTPThread = CLiteHTTPThread.CLiteHTTPThread(parent=self, dbQueue=cLiteDbThread.queue, parentDir=reportFilesPath, hasStartedQueue = self.responseQueue)
        cLiteHTTPThread.start()
        #Block with timeout until thread has started
        print(self.responseQueue.get(True, 10))

        response = cLiteDbThread.performAction('getProjectsAndJobs', win.responseQueue)
        if response != 'FailedDbInteraction':
            win.setInfo(response['projSums'], response['projJobs'], response['projDirs'], response['projIds'])

        print("Startup in ", time.time() - t1)

        modelData = """1 Hello
        1.1 is it me your
        3.2 Stranger to ...
        4 Captain Caveman
        1.2 looking for?
        2 Poor old Jonny Ray,
        3 Fergus sings the blues
        3.1.1.3 foo
        3.1.1.2 bar
        3.1.1.1 sna
        3.1 in bars of twelve
        2.1 blah, blah
        3.1.1 or less
        3.1.2 or more
        """

        """
        mod = TreeModel(modelData.split("\n"))
        tv = QtGui.QTreeView()
        tv.setModel(mod)
        tv.show()
        tv.raise_()
        """

    def aboutToQuitCallback(self):
        CLiteOneOffs.cLiteDbThread().shutdown()
        CLiteOneOffs.cLiteHTTPThread().shutdown()
        print("About to quit")

if __name__ == "__main__":
    #Here we use this file location to identify the appropriate i2 root
    parentDir = os.path.split(os.path.abspath(__file__))[0]
    grandParentDir = os.path.split(parentDir)[0]
    sys.path.insert(1, grandParentDir)
    from utils import startup

    # Proper rendering of svg files requires this environment variable to be set
    os.environ['QT_PLUGIN_PATH'] = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(QtCore.__file__))))),'qt4','plugins')
    print(os.environ['QT_PLUGIN_PATH'])

    qsSupervisorApp = QSSupervisorApp(sys.argv)
    from utils import QApp
    QApp.MYAPPLICATION = qsSupervisorApp
    sys.exit(qsSupervisorApp.exec_())
