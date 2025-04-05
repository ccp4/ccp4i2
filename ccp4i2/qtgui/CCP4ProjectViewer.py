"""
Copyright (C) 2010 University of York
Liz Potterton June 2011 - Project document in a QtMainWindow
"""

##@package CCP4ProjectViewer (QtGui)   Project document in a QtMainWindow

from pathlib import Path
import functools
import glob
import json
import math
import os
import shutil
import sys
import tempfile
import time

from lxml import etree
from PySide2 import QtCore, QtGui, QtWebChannel, QtWebEngineWidgets, QtWidgets
import clipper
import gemmi

from . import CCP4ProjectWidget
from . import CCP4TaskWidget
from . import CCP4WebBrowser
from . import CCP4WebToolBarButtons
from ..core import CCP4I2Runner
from ..core import CCP4Utils
from ..core.CCP4Config import CONFIG
from ..core.CCP4ErrorHandling import CErrorReport, CException, Severity
from ..core.CCP4TaskManager import TASKMANAGER
from ..core.CCP4WarningMessage import warningMessage
from ..dbapi import CCP4DbApi
from ..dbapi.CCP4DbUtils import COpenJob
from ..utils.QApp import QTAPPLICATION
from .CCP4MessageBox import CMessageBox


KJS_DEVINTER_ON = False
DEFAULT_WINDOW_SIZE = (1400, 800)
ALWAYS_SHOW_SERVER_BUTTON = False


def PROJECTVIEWER(projectId=None, open=False):
    for pv in CProjectViewer.Instances:
        if pv.taskFrame.openJob.projectId == projectId: return pv
    if open:
        pv = CProjectViewer(projectId=projectId)
        return pv
    return None


staticIconMap = {}
def getMenuIcon(parent, name, size=None):
    # Execute a reload if HD icon required (96x96) and has a non-HD pixmap loaded
    reloadHD = False
    if staticIconMap.get(name):
        if staticIconMap.get(name).availableSizes()[0].width() < 96 and size == "HD":
            reloadHD= True
    # Return & cache icons for the GUI (used for the task icons)
    if staticIconMap.get(name) and not reloadHD:
        return staticIconMap[name]
    else:
        pixFile = TASKMANAGER().searchIconFile(name)
        if size == "HD":
            pixFile = os.path.splitext(pixFile)[0] + "_96" + os.path.splitext(pixFile)[1]
            if not os.path.isfile(pixFile): # fallback in case no _96 file.
                pixFile = TASKMANAGER().searchIconFile(name)
        tmpIcn = QtGui.QIcon(pixFile)
        staticIconMap[name] = tmpIcn
        return tmpIcn

def currentlyOpenJobs():
    jobIdList = []
    for pv in CProjectViewer.Instances:
        jobIdList.append(pv.taskFrame.openJob.jobId)
        for tw in pv.findChildren(CTaskMainWindow):
            jobIdList.append(tw.objectName()[3:])
    return jobIdList

class CCP4WebReportBridge(QtCore.QObject):
    reloadDataResult = QtCore.Signal(int)
    jsSignal = QtCore.Signal(dict)
    @QtCore.Slot(str)
    def clicked(self,jsargs):
        launchParams = json.loads(jsargs)
        self.jsSignal.emit(launchParams)
    
    @QtCore.Slot(int)
    def reloadReportDataResult(self,total):
        self.reloadDataResult.emit(total)

class HTMLDelegate(QtWidgets.QStyledItemDelegate):

    helpClicked = QtCore.Signal(QtCore.QModelIndex)

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)
        qticonsDir = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons')
        self.pipelinePix = QtGui.QPixmap(os.path.join(qticonsDir,"gears.png") )
        self.helpPix = QtGui.QIcon(os.path.join(qticonsDir,"help-no-highlight.png")).pixmap(QtCore.QSize(28,28))
        self.helpPixHL = QtGui.QIcon(os.path.join(qticonsDir,"help-highlight.png")).pixmap(QtCore.QSize(28,28))

        self._compact = False
        self.requiredSize = None
        self.lastTextPos = QtCore.QPoint(0,0)
        self.showHelpButton = False
        self.showHelpClose = False
        self.showHelpIndex = None

    @QtCore.Slot()
    def clearShowHelp(self):
        self.showHelpButton = False
        self.showHelpClose = False

    @QtCore.Slot('QSize')
    def setRequiredSize(self,size):
        self.requiredSize = size

    def setCompact(self,comapct):
        self._compact = comapct

    def paint(self, painter, option, index):

        taskType = index.model().mapToSource(index).internalPointer().taskType()
        helpAvailable = index.model().mapToSource(index).internalPointer().helpAvailable()
        icon_size = index.model().mapToSource(index).internalPointer().iconSize()

        options = QtWidgets.QStyleOptionViewItem(option)
        widget = options.widget
        self.initStyleOption(options,index)
        
        style = QtWidgets.QApplication.style() if options.widget is None else options.widget.style()

        options.decorationSize = icon_size

        textWidth,textHeight = self.determineWidthHeight(options)

        doc = QtGui.QTextDocument()

        textOption = doc.defaultTextOption()
        doc.setDefaultTextOption(textOption)

        if options.font.pixelSize()> -1:
            doc.setHtml("<div style='vertical-align: top;font-size:"+str(options.font.pixelSize())+"px;'>"+options.text+"</div>")
            theIconSize = max(options.font.pixelSize(),icon_size.width())
            options.decorationSize = QtCore.QSize(theIconSize,theIconSize)
        elif options.font.pointSize()> -1:
            doc.setHtml("<div style='vertical-align: top;font-size:"+str(options.font.pointSize())+"pt;'>"+options.text+"</div>")

        textOption.setWrapMode(QtGui.QTextOption.WordWrap)
        doc.setTextWidth(options.rect.width());
        options.text = ""

        style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, options, painter);

        ctx = QtGui.QAbstractTextDocumentLayout.PaintContext()

        # Highlighting text if item is selected
        if options.state & QtWidgets.QStyle.State_Selected:
            ctx.palette.setColor(QtGui.QPalette.Text, options.palette.color(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText));

        textRect = style.subElementRect(QtWidgets.QStyle.SE_ItemViewItemText, options)

        if index.column() != 0:
            textRect.adjust(5, 0, 0, 0)
 
        textRect.setTop(textRect.top())

        painter.save()
        painter.translate(textRect.topLeft())

        if taskType == "taskpipe" and not self._compact:
            painter.drawPixmap(options.rect.width()-self.pipelinePix.width()-options.rect.left(),0,self.pipelinePix)

        painter.setClipRect(textRect.translated(-textRect.topLeft()))
        doc.documentLayout().draw(painter, ctx)

        if self.showHelpButton and helpAvailable:
            if index.parent().row() != -1:
                if self.showHelpIndex is not None and  self.showHelpIndex.row() == index.row() and self.showHelpIndex.column() == index.column():
                    if self.showHelpIndex.parent().row() == index.parent().row() and self.showHelpIndex.parent().column() == index.parent().column():
                        if self.showHelpClose:
                            painter.drawPixmap(options.rect.width()-self.helpPixHL.width()-options.rect.left()-50,int((options.rect.height()-self.helpPixHL.height())/2),self.helpPixHL)
                        else:
                            painter.drawPixmap(options.rect.width()-self.helpPix.width()-options.rect.left()-50,int((options.rect.height()-self.helpPix.height())/2),self.helpPix)

        painter.restore()


    def determineWidthHeight(self, options):
        text = options.text
        totalHeight = 0.0

        for t in text.split("<br/>"):
            doc = QtGui.QTextDocument()
            doc.setTextWidth(10000);

            if options.font.pixelSize()> -1:
                doc.setHtml("<div style='vertical-align: middle;font-size:"+str(options.font.pixelSize())+"px;'>"+t+"</div>")
            elif options.font.pointSize()> -1:
                doc.setHtml("<div style='vertical-align: middle;font-size:"+str(options.font.pointSize())+"pt;'>"+t+"</div>")
            nrows = 1
            if self.requiredSize is not None:
                nrows = int(math.ceil(1.0*doc.idealWidth()/(self.requiredSize.width()-64)))
            totalHeight += 1.0*nrows*doc.documentLayout().documentSize().height()
        return doc.idealWidth(),1.0*totalHeight

    def sizeHint(self, option, index):
        #Need to determine if we are going to wrap in here and hence multiply doc.documentLayout().documentSize().height() by number of lines
        options = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(options,index)

        icon_size = index.model().mapToSource(index).internalPointer().iconSize()
        theIconSize = max(options.font.pixelSize(),icon_size.width())

        idealWidth,totalHeight = self.determineWidthHeight(options)

        return QtCore.QSize(idealWidth, max(totalHeight,theIconSize))

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
        self._description = ""
        self._task_type = ""
        self._help_available = False
        self._icon_size = QtCore.QSize(16,16)

    def setIconSize(self,icon_size):
        self._icon_size = icon_size

    def iconSize(self):
        return self._icon_size

    def setHelpAvailable(self,help_available):
        self._help_available = help_available

    def setTaskType(self,task_type):
        self._task_type = task_type

    def taskType(self):
        return self._task_type

    def helpAvailable(self):
        return self._help_available

    def setLongDescription(self,desc):
        self._description = desc

    def longDescription(self):
        return self._description

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

#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def flags(self,index):
        if not index.isValid():
            return 0

        return QtCore.QAbstractItemModel.flags(self, index)

    def data(self, index, role):
        #FIXME - Maybe I should have columns instead of all these roles? Maybe I do already?
        if not index.isValid():
#FIXME PYQT - or maybe None? This used to return QVariant.
            return None

        item = index.internalPointer()

        if role == QtCore.Qt.UserRole + 6:
            print("return description")
            return item._description

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
#FIXME PYQT - or maybe None? This used to return QVariant.
            return None

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


class CTaskListProxyModel(QtCore.QSortFilterProxyModel):
        def __init__(self,parent=None):
            QtCore.QSortFilterProxyModel.__init__(self,parent)
            self.filterString = ""
            self.filterMode = ()

        def setFilterMode(self,string):
            self.filterMode = string
            self.invalidateFilter()

        def setFilterString(self,string):
            self.filterString = string
            self.invalidateFilter()

        def hasAcceptedChildrenOrIsAccepted(self, item):

            theData = item.data(QtCore.Qt.DisplayRole) + item.longDescription()
            theDataKeywords = item.keywords()

            acceptedFilterMode = True

            if len(self.filterMode) > 0:
                taskmodule = TASKMANAGER().getTaskAttribute(item.name(), 'TASKMODULE')
                if type(taskmodule) == list:
                    if len(set(taskmodule).intersection(set(self.filterMode)))==0:
                        acceptedFilterMode = False
                if type(taskmodule) == str:
                    if len(set([taskmodule]).intersection(set(self.filterMode)))==0:
                        acceptedFilterMode = False
                if taskmodule is None:
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


#------------------------------------------------------------------------------------------------------
class CProjectViewer(CCP4WebBrowser.CMainWindow):
#------------------------------------------------------------------------------------------------------

    jobSearch = QtCore.Signal()

    ERROR_CODES = {101 : {'description' : 'Error creating task widget for'},
                   102 : {'description' : 'Wrong task name in imported params file'},
                   103 : {'description' : 'Error attempting to auto-generate task widget for'},
                   104 : {'description' : 'Error writing job parameters file'},
                   120 : {'description' : 'Unknown error attempting to kill remote job'}}
    INPUT_TAB = 0
    OUTPUT_TAB = 1
    COMMENT_TAB = 2
    MARGIN=2
    Instances = []

    def isNotProjectViewer(self):
            return False

    @QtCore.Slot()
    def setToolButtonStyle(self):
        from ..core.CCP4Preferences import PREFERENCES
        if PREFERENCES().TOOLBARBUTTONSSTYLE == 1:
            self.webviewToolBar.page().runJavaScript("setButtonStyle('textOnly')")
            self.dockWidget.setFixedHeight(57)
        elif PREFERENCES().TOOLBARBUTTONSSTYLE == 0:
            self.webviewToolBar.page().runJavaScript("setButtonStyle('iconOnly')")
            self.dockWidget.setFixedHeight(57)
        else:
            self.webviewToolBar.page().runJavaScript("setButtonStyle('textUnderIcon')")
            self.dockWidget.setFixedHeight(60+PREFERENCES().GUI_FONT_SIZE)

    def __init__(self, parent=None, projectId=None, jobId=None, graphical=True):
        CCP4WebBrowser.CMainWindow.__init__(self, parent)
        CProjectViewer.Instances.append(self)
        self.setObjectName('projectViewer')
        #self.destroyed.connect(CProjectViewer.updateInstances)
        self._dictionaryWidget = None
        self.mainwins = []
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        if jobId is None:
            jobIdList = PROJECTSMANAGER().db().getProjectJobList(projectId=projectId, topLevelOnly=True, maxLength=1)
            if len(jobIdList) > 0:
                jobId = jobIdList[0]
        openJob = COpenJob(projectId=projectId, jobId=jobId)
        PROJECTSMANAGER().db().updateProject(projectId=projectId, key='lastaccess')
        self.setWindowTitle(self.version + 'Project Viewer: ' + openJob.projectName)
        self.layout().setContentsMargins(CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN)
        self.layout().setSpacing(CProjectViewer.MARGIN)
        # left side project widget and buttons
        leftFrame = QtWidgets.QFrame(self)
        leftFrame.setLayout(QtWidgets.QVBoxLayout())
        leftFrame.layout().setContentsMargins(CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN)
        leftFrame.layout().setSpacing(CProjectViewer.MARGIN)
        self._projectWidget = CCP4ProjectWidget.CProjectWidget(self, projectId)
        self._projectWidget.setMinimumSize(QtCore.QSize(300, 300))
        leftFrame.layout().addWidget(self._projectWidget)
        # Right side is task frame or task chooser
        self.rightStack = QtWidgets.QStackedWidget(self)
        self.taskFrame = CTaskFrame(self, projectId=openJob.projectId)
        self.taskFrame.buttons.frame.hide()
        self.rightStack.addWidget(self.taskFrame)
        self.taskChooser = CChooseTaskFrame(self)

        self._projectWidget.model().sourceModel().rowsInserted.connect(self.handleNumberOfJobsChanged)
        self._projectWidget.model().sourceModel().rowsRemoved.connect(self.handleNumberOfJobsChanged)

        rootIdx = QtCore.QModelIndex()
        modelRows = self._projectWidget.model().sourceModel().rowCount(rootIdx)
        jobList = PROJECTSMANAGER().db().projectJobSearch(self._projectWidget.projectId)
        if modelRows == 0:
            self.taskChooser.cancelButton.setEnabled(False)

        self.rightStack.addWidget(self.taskChooser)
        self.rightStack.setCurrentIndex(0)
        # left/right splitter
        self.splitterSizes = [320, CCP4TaskWidget.WIDTH + CCP4TaskWidget.PADDING_ALLOWANCE]
        mainWidget = QtWidgets.QSplitter(self)
        mainWidget.setOrientation(QtCore.Qt.Horizontal)
        mainWidget.addWidget(leftFrame)
        mainWidget.addWidget(self.rightStack)
        self.rightStack.layout().setContentsMargins(CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN, CProjectViewer.MARGIN)
        self.rightStack.layout().setSpacing(CProjectViewer.MARGIN)
        # And set the splitter sizes after show ...
        self.setCentralWidget(mainWidget)
        self.taskFrame.outputFrame.webView.csvDownloadRequest.connect(self.handleCsvDownloadRequest)
        from ..core.CCP4Preferences import PREFERENCES
        @QtCore.Slot(int)
        def zoomFactorChanged(val):
            PREFERENCES().REPORT_ZOOM_FACTOR.set(self.taskFrame.outputFrame.webView.zoomFactor())
        self.taskFrame.outputFrame.webView.zoomFactorChanged.connect(zoomFactorChanged)
        self.taskFrame.outputFrame.webView.setZoomFactor(PREFERENCES().REPORT_ZOOM_FACTOR)
        # Signals to open a new task
        self.taskFrame.buttons.button('clone').clicked.connect(self.cloneTask)
        self.taskFrame.buttons.button('help').menu().aboutToShow.connect(self.showHelpMenu)
        self.taskFrame.buttons.button('help').menu().triggered.connect(self.handleHelpMenu)
        self.taskFrame.outputFrame.nextTask.connect(self.handleNextTask)
        self.taskFrame.outputFrame.interrupt.connect(functools.partial(self.stopJob, False, False))
        self.taskFrame.launchJobRequestSignal.connect(self.launchJobRequest)
        self.taskChooser.taskClicked.connect(self.handleChooseTask)
        self.taskChooser.taskTree.taskDropped.connect(self.handleChooseTaskAndFollowFrom)
        # Show chooser
        self.taskChooser.closed.connect(functools.partial(self.showTaskChooser, False))
        # Handle job/project deleted
        PROJECTSMANAGER().db().jobDeleted.connect(self.handleJobDeleted)
        PROJECTSMANAGER().db().jobStatusUpdated.connect(self.taskFrame.titleBar.setStatusBar)
        PROJECTSMANAGER().db().jobStatusUpdated.connect(self.updateActionEnabled)
        PROJECTSMANAGER().db().jobFinished.connect(self.taskFrame.titleBar.setStatusBar)
        PROJECTSMANAGER().db().projectDeleted.connect(self.handleProjectDeleted)
        self.taskFrame.labelEdited.connect(self.redoReport)
        @QtCore.Slot()
        def remakeReport80(openJob):
           print("Attempting to remake report")
           self.taskFrame.outputFrame.showOutput(redo=True, doReload=True, openJob=openJob)
        self.taskFrame.report71.connect(remakeReport80)

#FIXME (SJM 3/10/2019) - I'm hacking at the moment by still creating the CTOolbar on Mac and hiding it. This is because the toolbar and signals are intimately entwined at the moment. I have to disentagle.
        self.addToolBar(CCP4WebBrowser.CToolBar(self, 'project', 'project'))
        toolbar = self.toolBar('project') 
        if toolbar is not None:
            from ..qtcore.CCP4JobController import JOBCONTROLLER
            if ALWAYS_SHOW_SERVER_BUTTON or JOBCONTROLLER().serversEnabled():
                toolbar.extend(['task_menu', 'export_project', 'sep', 'run', 'run_remote', 'clone',
                            'task_help', 'references', 'export_mtz', 'show_log', 'view_coot', 'view_ccp4mg', 'show_i2run', 'sep', 'new_project2'])
            else:
                toolbar.extend(['task_menu', 'export_project', 'sep', 'run', 'clone',
                            'task_help', 'references', 'export_mtz', 'show_log', 'view_coot', 'view_ccp4mg', 'show_i2run', 'sep', 'new_project2'])
            toolbar.setMovable(False)
            toolbar.show()
        if jobId is not None:
            try:
                if PREFERENCES().RESTORE_TO_TASKLIST:
                    self.showTaskChooser(True)
                else:
                    # Note that opening a job like this is slow (avoid on start).
                    openJob = self.taskFrame.openTask(jobId=jobId)
                    self.setSelectedJob(jobId=openJob.jobId)
            except:
                print('ERROR CProjectView.init opening jobId',jobId)
        else:
            #Initialise task chooser open if not jobs in project
            self.showTaskChooser(True)
        if graphical:
            self.show()
        mainWidget.setSizes(self.splitterSizes)
        self.taskFrame.tab.currentChanged[int].connect(self.handleTaskFrameChanged)
        #self.lastTaskFrameMode = CProjectViewer.OUTPUT_TAB
        #self.handleTaskFrameChanged(CProjectViewer.INPUT_TAB)
        self.shortcuts = {}
        self.setShortcuts()
        self.updateActionEnabled()

        label = QtWidgets.QLabel("Am I in the dark?")
        text_hsv_value = label.palette().color(QtGui.QPalette.WindowText).value()
        bg_hsv_value = label.palette().color(QtGui.QPalette.Background).value()
        isDarkMode = text_hsv_value > bg_hsv_value
        if sys.platform == "darwin" and not isDarkMode:
            self.toolBar('project').hide()
            self.webviewToolBar = CCP4WebToolBarButtons.CCP4WebToolBarButtons(fileName=os.path.join(os.path.dirname(__file__),"html_i2_buttons.html"))
            dockWidget = QtWidgets.QDockWidget()
            dockWidget.setTitleBarWidget(QtWidgets.QWidget());
            dockWidget.setWidget(self.webviewToolBar)
            self.addDockWidget(QtCore.Qt.TopDockWidgetArea,dockWidget)
            dockWidget.setFixedHeight(70)    # Don't want any scroolbars or resize decorations.
            dockWidget.setMinimumWidth(1000) # Don't want the toolbar scrunching.
            self.dockWidget = dockWidget
            @QtCore.Slot(str)
            def handleToolBarClick(buttonName):
                if self.findChild(QtWidgets.QAction, buttonName):
                    self.findChild(QtWidgets.QAction, buttonName).trigger()
            self.webviewToolBar.buttonClicked.connect(handleToolBarClick)
            self.webviewToolBar.loadFinished.connect(functools.partial(self.updateActionEnabled,None))
            PREFERENCES().GUI_FONT_SIZE.dataChanged.connect(self.setToolButtonFontSize)

            @QtCore.Slot()
            def editVisible():
                pynames = self.jsnames

                prefWidget = QtWidgets.QDialog()
                prefWidget.setWindowTitle("Edit visible tool buttons")
                listWidget = QtWidgets.QListWidget()
                layout = QtWidgets.QVBoxLayout()
                prefWidget.setLayout(layout)
                layout.addWidget(listWidget)

                for v in pynames:
                    name = v
                    if str(name) in CCP4WebBrowser.CToolBar.toolBarPreferencesMapping:
                        mapping = CCP4WebBrowser.CToolBar.toolBarPreferencesMapping[str(name)]
                        if self.findChild(QtWidgets.QAction, name):
                            child = self.findChild(QtWidgets.QAction, name)
                            item = QtWidgets.QListWidgetItem()
                            item.setText(child.text())
                            item.setIcon(child.icon())
                            listWidget.addItem(item)
                            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
                            item.setCheckState(QtCore.Qt.Unchecked)
                            item.setData(QtCore.Qt.UserRole,name)
                            #MN Trying to make code more compact/readable
                            if mapping == "SHOW_RUN_REMOTE_TOOLBUTTON":
                                from ..qtcore.CCP4JobController import JOBCONTROLLER
                                if ALWAYS_SHOW_SERVER_BUTTON or JOBCONTROLLER().serversEnabled():
                                   val = PREFERENCES().SHOW_RUN_REMOTE_TOOLBUTTON
                                else:
                                   val = False
                            else:
                                val = getattr(PREFERENCES(), mapping)
                            #replaces:
                            if val:
                                item.setCheckState(QtCore.Qt.Checked)

                @QtCore.Slot('QlistWidgetItem')
                def setItemVisibilities(item):
                    checkState = item.checkState()
                    name = item.data(QtCore.Qt.UserRole)
                    if checkState:
                        self.webviewToolBar.setButtonVisible(str(name),True)
                        val = True
                    else:
                        self.webviewToolBar.setButtonVisible(str(name),False)
                        val = False
                    mapping = CCP4WebBrowser.CToolBar.toolBarPreferencesMapping[str(name)]
                    
                    #MN Trying to make code more compact/readable
                    getattr(PREFERENCES(), mapping).set(val)
                    #Replacing
                listWidget.itemChanged.connect(setItemVisibilities)
                prefWidget.exec_()

            self.jsnames = []

            @QtCore.Slot()
            def updateActionVisible():
                
                self.menuBar().updateMenu("File")
                PREFERENCES().TOOLBARBUTTONSSTYLE.dataChanged.connect(self.setToolButtonStyle)
                self.setToolButtonStyle()

                pynames = self.jsnames
                for v in pynames:
                    name = v
                    if str(name) in CCP4WebBrowser.CToolBar.toolBarPreferencesMapping:
                        mapping = CCP4WebBrowser.CToolBar.toolBarPreferencesMapping[str(name)]
                        #MN Trying to increase code compactness/readability
                        if mapping == "SHOW_RUN_REMOTE_TOOLBUTTON":
                            from ..qtcore.CCP4JobController import JOBCONTROLLER
                            if ALWAYS_SHOW_SERVER_BUTTON or JOBCONTROLLER().serversEnabled():
                               val = PREFERENCES().SHOW_RUN_REMOTE_TOOLBUTTON
                            else:
                               val = False
                        else:
                            val = getattr(PREFERENCES(), mapping)
                        self.webviewToolBar.setButtonVisible(str(name),val)

                    else:
                        if CONFIG().developer: print(name,"not in mapping")
            
            @QtCore.Slot(list)
            def handleButtonListUpdated(buttonList):
                self.jsnames = buttonList
                updateActionVisible()

            @QtCore.Slot()
            def getNamesAndUpdateActionVisible():
                self.webviewToolBar.buttonList.connect(handleButtonListUpdated)
                self.webviewToolBar.page().runJavaScript("requestButtonNames()")
            
            self.webviewToolBar.loadFinished.connect(getNamesAndUpdateActionVisible)
            self.webviewToolBar.editVisible.connect(editVisible)

        else:
            self.toolBar('project').show()

        # Handle projectWidget actions
        self._projectWidget.purge.connect(self.purgeJobFiles)
        self._projectWidget.cloneTask.connect(self.cloneTask)
        self._projectWidget.export.connect(self.exportJobFile)
        self._projectWidget.deleteJob.connect(self.deleteJob)
        self._projectWidget.interruptJob.connect(functools.partial(self.stopJob, False, False))
        self._projectWidget.killJob.connect(functools.partial(self.stopJob, True, False))
        self._projectWidget.killDeleteJob.connect(functools.partial(self.stopJob, True, True))
        self._projectWidget.markFailedJob.connect(functools.partial(self.markFailedJob, False))
        self._projectWidget.markFinishedJob.connect(functools.partial(self.markFailedJob, True))
        self._projectWidget.currentJobChanged.connect(self.handleJobListSelectionChange)
        self._projectWidget.nextJob.connect(self.handleNext)
        self._projectWidget.openFrame.connect(self.openTaskMainWindow)
        self._projectWidget.showWhatNextPage.connect(self.showWhatNextPage)
        self._projectWidget.showBibliography.connect(self.showBibliography)
        self._projectWidget.openMergeMtz.connect(functools.partial(self.launchJobRequest, 'mergeMtz'))
        pInfo = PROJECTSMANAGER().db().getProjectInfo(projectId, ['I1ProjectName', 'I1ProjectDirectory'])
        if pInfo['i1projectname'] is not None:
            self.addI1Widget(pInfo['i1projectname'], pInfo['i1projectdirectory'])
        self.checkForRemotePasswords()
 
    def resetZoom(self):
        if hasattr(self,"taskFrame") and hasattr(self.taskFrame,"outputFrame") and hasattr(self.taskFrame.outputFrame,"webView"):
            self.taskFrame.outputFrame.webView.setZoomFactor(1.0)
            from ..core.CCP4Preferences import PREFERENCES
            PREFERENCES().REPORT_ZOOM_FACTOR.set(self.taskFrame.outputFrame.webView.zoomFactor())

    def zoomIn(self):
        if hasattr(self,"taskFrame") and hasattr(self.taskFrame,"outputFrame") and hasattr(self.taskFrame.outputFrame,"webView"):
            self.taskFrame.outputFrame.webView.setZoomFactor(self.taskFrame.outputFrame.webView.zoomFactor()*1.2)
            from ..core.CCP4Preferences import PREFERENCES
            PREFERENCES().REPORT_ZOOM_FACTOR.set(self.taskFrame.outputFrame.webView.zoomFactor())

    def zoomOut(self):
        if hasattr(self,"taskFrame") and hasattr(self.taskFrame,"outputFrame") and hasattr(self.taskFrame.outputFrame,"webView"):
            self.taskFrame.outputFrame.webView.setZoomFactor(self.taskFrame.outputFrame.webView.zoomFactor()/1.2)
            from ..core.CCP4Preferences import PREFERENCES
            PREFERENCES().REPORT_ZOOM_FACTOR.set(self.taskFrame.outputFrame.webView.zoomFactor())
 
    @QtCore.Slot('QModelIndex',int,int)
    def handleNumberOfJobsChanged(self,parentIdx,startRow,endRow):
        #Disable "Cancel" button if there are no jobs in project
        rootIdx = QtCore.QModelIndex()
        modelRows = self._projectWidget.model().sourceModel().rowCount(rootIdx)
        if modelRows == 0:
            self.taskChooser.cancelButton.setEnabled(False)
        else:
            self.taskChooser.cancelButton.setEnabled(True)
            #self.updateActionEnabled()

    def close(self):
        if sys.platform == "darwin":
            from ..core.CCP4Preferences import PREFERENCES
            PREFERENCES().TOOLBARBUTTONSSTYLE.dataChanged.disconnect(self.setToolButtonStyle)
        # Update the last access info for all open projects
        self._projectWidget.model().sourceModel().rowsInserted.disconnect(self.handleNumberOfJobsChanged)
        self._projectWidget.model().sourceModel().rowsRemoved.disconnect(self.handleNumberOfJobsChanged)

        openProjectIdList = []
        for w in self.Instances:
            if w != self:
                openProjectIdList.append(w.taskFrame.openJob.projectId)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().updateProjectLastAcess(self.taskFrame.openJob.projectId, openProjectIdList)
        CCP4WebBrowser.CMainWindow.close(self)

    def projectId(self):
        return self.taskFrame.openJob.projectId

    def addI1Widget(self, projectName, projectDir):
        '''Open an i1 project viewer'''
        try:
            from . import CCP4I1Projects
            frame = QtWidgets.QFrame(self)
            frame.setLayout(QtWidgets.QVBoxLayout())
            frame.layout().setContentsMargins(0, 0, 0, 0)
            frame.layout().setSpacing(0)
            self._projectWidget.tab.addTab(frame, 'CCP4i project')
            self.i1Widget = CCP4I1Projects.CI1ProjectWidget(self)
            frame.layout().addWidget(self.i1Widget)
            #self._projectWidget.tab.addTab(self.i1Widget,'CCP4i project')
            infoList = [projectName, projectName]
            pObj = CCP4I1Projects.CI1TreeItemProject(self.i1Widget.model().rootItem, infoList=infoList, directory=projectDir)
            err=pObj.loadDatabase()
            if err.maxSeverity() > Severity.WARNING:
                warningMessage(err, 'CCP4i project', 'Error loading old CCP4 project data', parent=self)
            else:
                self.i1Widget.model().beginResetModel()
                self.i1Widget.model().rootItem = pObj
                self.i1Widget.model().endResetModel()
        except CException as e:
            print(e)
            return
        try:
            line = QtWidgets.QHBoxLayout()
            line.setContentsMargins(3, 3, 3, 3)
            line.setSpacing(2)
            frame.layout().addLayout(line)
            line.addWidget(QtWidgets.QLabel('CCP4i project: ' + projectName))
            self.i1Reload = QtWidgets.QPushButton(self)
            self.i1Reload.setText('Reload')
            line.addWidget(self.i1Reload)
            self.i1Reload.clicked.connect(self.reloadI1Project)
        except:
            print(e)
        self.i1Widget.showLogFile.connect(functools.partial(self.viewI1File, 'log'))
        self.i1Widget.projectView.fileClicked.connect(functools.partial(self.viewI1File, 'clicked'))
        self.i1Widget.viewInQtrview.connect(functools.partial(self.viewI1File, 'qtrview'))
        self.i1Widget.viewInCoot.connect(functools.partial(self.viewI1Job, 'coot'))
        self.i1Widget.viewInMg.connect(functools.partial(self.viewI1Job, 'ccp4mg'))
        self.i1Widget.viewDefFile.connect(functools.partial(self.viewI1Job, 'def'))
        self.i1Widget.reloadProject.connect(self.reloadI1Project)
        self.i1Widget.viewFile.connect(functools.partial(self.viewI1File, 'view'))
        self.i1Widget.copyFile.connect(self.copyI1File)

    @QtCore.Slot(str,str,str)
    def viewI1Job(self,mode, projectId=None, jobId=None):
        '''View some file from a job from an i1 project'''
        jItem = self.i1Widget.model().getJob1(jobId)
        if jItem is None:
            print('ERROR in viewI1Job - could not find job:', jobId)
            return
        if mode == 'def':
            defFile = os.path.join(self.i1Widget.model().rootItem.directory, 'CCP4_DATABASE', str(jobId) + '_' + jItem.taskName + '.def')
            from .CCP4WebBrowser import WEBBROWSER
            WEBBROWSER().openFile(defFile, format="text/plain")
        if mode in ['coot', 'ccp4mg']:
            fileList = []
            for fItem in jItem.childOutFiles:
                if fItem.broken>0:
                    path = fItem.filePath()
                    if os.path.splitext(path)[1] in ['.mtz', '.pdb']:
                        fileList.append(fItem.filePath())
            if mode == 'coot':
                mode = 'coot0'
            from ..qtcore.CCP4Launcher import LAUNCHER
            LAUNCHER().openInViewer(viewer=mode, fileName=fileList)

    @QtCore.Slot(str,str,str)
    def viewI1File(self, mode, fileName, fileType=None):
        '''View a file from i1project viewer'''
        if fileName is None:
            return
        from ..qtcore.CCP4Launcher import LAUNCHER
        if mode == 'qtrview':
            if os.path.splitext(fileName)[1] == '.html':
                fileName = os.path.splitext(fileName)[0]
            LAUNCHER().launch('logview', [fileName])
        elif mode in ['view', 'clicked']:
            if fileType is None:
                ext = os.path.splitext(fileName)[1]
                if ext == '.mtz':
                    fileType = "application/CCP4-mtz"
                elif ext == '.pdb':
                    fileType = "chemical/x-pdb"
                else:
                    fileType = "text/plain"
            from .CCP4WebBrowser import WEBBROWSER
            if fileType == "application/CCP4-mtz":
                LAUNCHER().launch('hklview', [fileName])
            elif fileType == "chemical/x-pdb":
                WEBBROWSER().openFile(fileName)
            else:
                WEBBROWSER().openFile(fileName)

    @QtCore.Slot()
    def reloadI1Project(self):
        self.i1Widget.model().beginResetModel()
        err = self.i1Widget.model().rootItem.loadDatabase()
        if err.maxSeverity() > Severity.WARNING:
            warningMessage(err, 'CCP4i project', 'Error loading old CCP4 project data', parent=self)
        self.i1Widget.model().endResetModel()

    @QtCore.Slot(str)
    def copyI1File(self, fileName):
        if os.path.exists(fileName):
            urlList = [QtCore.QUrl()]
            urlList[0].setPath(fileName )
            urlList[0].setScheme("file")
            mimeData = QtCore.QMimeData()
            mimeData.setUrls(urlList)
            QTAPPLICATION().clipboard().setMimeData(mimeData)

    @QtCore.Slot(dict)
    def updateActionEnabled(self, jobStatus=None, dummy = None):
        if jobStatus is None:
            jobStatus = self.taskFrame.openJob.status
        elif isinstance(jobStatus, dict):
            if not jobStatus['jobId'] == self.taskFrame.openJob.jobId:
                return
            jobStatus = jobStatus['status']
        if jobStatus is None:
            #No job set - disable all buttons
            for item in ['run', 'clone', 'task_help','show_i2run', 'show_log', 'references', 'export_mtz']:
                if hasattr(self,"webviewToolBar"):
                    self.webviewToolBar.setButtonEnabled(item,False)
                self.findChild(QtWidgets.QAction, item).setEnabled(False)
                for item in ['run_remote']:
                    if hasattr(self,"webviewToolBar"):
                        self.webviewToolBar.setButtonEnabled(item,False)
                    obj = self.findChild(QtWidgets.QAction, item)
                    if obj is not None:
                        obj.setEnabled(False)
        else:
            if self.taskFrame.buttons.mode in ['task', 'input']:
                self.findChild(QtWidgets.QAction, 'run').setEnabled(jobStatus in ['Pending', 'Interrupted'])
                self.findChild(QtWidgets.QAction, 'view_coot').setEnabled(jobStatus in ['Finished', 'Interrupted'])
                self.findChild(QtWidgets.QAction, 'view_ccp4mg').setEnabled(jobStatus in ['Finished', 'Interrupted'])
                if hasattr(self,"webviewToolBar"):
                    self.webviewToolBar.setButtonEnabled('run',jobStatus in ['Pending', 'Interrupted'])
                    self.webviewToolBar.setButtonEnabled('run_remote',jobStatus in ['Pending', 'Interrupted'])
                    self.webviewToolBar.setButtonEnabled('view_coot',jobStatus in ['Finished', 'Interrupted'])
                    self.webviewToolBar.setButtonEnabled('view_ccp4mg',jobStatus in ['Finished', 'Interrupted'])
                obj =  self.findChild(QtWidgets.QAction, 'run_remote')
                if obj is not None:
                    obj.setEnabled(jobStatus in ['Pending', 'Interrupted'])
            else:
                self.findChild(QtWidgets.QAction, 'run').setEnabled(False)
                obj = self.findChild(QtWidgets.QAction,'run_remote')
                if obj is not None:
                    obj.setEnabled(False)
            #self.findChild(QtWidgets.QAction,'view').setEnabled(jobStatus in ['Finished','Interrupted','To delete'])
            self.findChild(QtWidgets.QAction, 'clone').setEnabled(True)
            if hasattr(self,"webviewToolBar"):
                self.webviewToolBar.setButtonEnabled('clone',True)
            refFileList = TASKMANAGER().searchReferenceFile(name=self.taskFrame.openJob.taskName, drillDown=True)
            if len(refFileList) == 0:
                # Is ther acustom biblio for this job
                refFileList = glob.glob(os.path.join(self.taskFrame.openJob.jobDir, '*.medline.txt'))
            self.findChild(QtWidgets.QAction, 'references').setEnabled((len(refFileList) > 0))
            enabled = False
            exportMenu = []
            if jobStatus in ['Finished', 'Interrupted']:
                exportMenu = TASKMANAGER().exportJobFiles(taskName=self.taskFrame.openJob.taskName, jobId=self.taskFrame.openJob.jobId)
                for item in exportMenu:
                    if item[0] == 'complete_mtz':
                        enabled = True
                if not enabled:
                    params = TASKMANAGER().getTaskAttribute(taskName=self.taskFrame.openJob.taskName, attribute='EXPORTMTZPARAMS', default=None)
                    if params is not None:
                        enabled = True
            self.findChild(QtWidgets.QAction,'export_mtz').setEnabled(enabled)
            self.findChild(QtWidgets.QAction,'show_log').setEnabled(jobStatus in ['Running','Finished','Failed','Interrupted','To delete','Unsatisfactory'])
            helpFile = TASKMANAGER().searchHelpFile(name=self.taskFrame.openJob.taskName)
            self.findChild(QtWidgets.QAction,'task_help').setEnabled((helpFile is not None))
            self.findChild(QtWidgets.QAction,'show_i2run').setEnabled(True)
            if hasattr(self,"webviewToolBar"):
                self.webviewToolBar.setButtonEnabled('export_mtz',enabled)
                self.webviewToolBar.setButtonEnabled('show_i2run',jobStatus in ['Pending','Running','Finished','Failed','Interrupted','To delete','Unsatisfactory'])
                self.webviewToolBar.setButtonEnabled('show_log',jobStatus in ['Running','Finished','Failed','Interrupted','To delete','Unsatisfactory'])
                self.webviewToolBar.setButtonEnabled('show_i2run',True)
                self.webviewToolBar.setButtonEnabled('task_help',helpFile is not None)
                self.webviewToolBar.setButtonEnabled('references',(len(refFileList) > 0))
                self.setToolButtonFontSize()

    @QtCore.Slot(str)
    def setToolButtonFontSize(self):
        settings = self.webviewToolBar.settings()
        """
        css = bytes(".button {font-size: "+str(PREFERENCES().GUI_FONT_SIZE)+"px;}","utf-8")
        cssUrl = QtCore.QUrl("data:text/css;charset=utf-8;base64,"+base64.b64encode(css).decode())
        settings.setUserStyleSheetUrl(cssUrl)
        """
        from ..core.CCP4Preferences import PREFERENCES
        settings.setFontSize(QtWebEngineWidgets.QWebEngineSettings.DefaultFontSize,int(PREFERENCES().GUI_FONT_SIZE))
        settings.setFontSize(QtWebEngineWidgets.QWebEngineSettings.DefaultFixedFontSize,int(PREFERENCES().GUI_FONT_SIZE))
        self.setToolButtonStyle()

    @QtCore.Slot(str)
    def runTask(self, mode='Now'):
        self.taskFrame.inputFrame.runTask(mode)

    @QtCore.Slot(str)
    def handleViewTask(self, mode):
        self.taskFrame.handleViewTask(mode)

    def initialiseActionDefinitions(self):
        label = QtWidgets.QLabel("Am I in the dark?")
        text_hsv_value = label.palette().color(QtGui.QPalette.WindowText).value()
        bg_hsv_value = label.palette().color(QtGui.QPalette.Background).value()
        isDarkMode = text_hsv_value > bg_hsv_value

        CCP4WebBrowser.CMainWindow.initialiseActionDefinitions(self)
        self.actionDefinitions['task_menu'] = dict(text = "Task menu", tip="Show/hide the task menu",
                                                   slot = functools.partial(self.showTaskChooser, True), icon='taskmenu')
        self.actionDefinitions['clone'] = dict(text="Clone job", tip="Make another job with same parameters",
                                               slot=self.cloneTask, icon='clone')
        if isDarkMode:
            self.actionDefinitions['run'] = dict(text="Run", tip="Run this task", slot=functools.partial(self.runTask, 'local'), icon = 'running_dark')
            self.actionDefinitions['run_remote'] = dict(text="Run on server", tip="Run this job on a server",
                                                    slot=functools.partial(self.runTask, 'run_remote'), icon='running_dark')
        else:
            self.actionDefinitions['run'] = dict(text="Run", tip="Run this task", slot=functools.partial(self.runTask, 'local'), icon = 'running')
            self.actionDefinitions['run_remote'] = dict(text="Run on server", tip="Run this job on a server",
                                                    slot=functools.partial(self.runTask, 'run_remote'), icon='running')

        self.actionDefinitions['job_search'] = dict(text="Job search", tip="Search job list",
                                                    slot=self.jobSearch.emit, icon='search')
        self.actionDefinitions['view_coot'] = dict(text="View in Coot", tip="Show map & models related to the job in Coot",
                                                   slot=functools.partial(self.handleViewTask, 'view_coot'), icon='coot_rebuild')
        self.actionDefinitions['view_ccp4mg'] = dict(text="View in CCP4mg", tip="Show map & models related to the job in CCP4mg",
                                                     slot=functools.partial(self.handleViewTask,'view_ccp4mg'), icon='ccp4mg_edit_model')
        self.actionDefinitions['export_mtz'] = dict(text="Export MTZ", tip="Export traditional MTZ file with all the data for the job",
                                                    slot=self.exportMtz, icon='MiniMtzDataFile')
        self.actionDefinitions['task_help'] = dict(text="Help",tip="Show documentation for this task",
                                                   slot=functools.partial(self.handleHelpMenu, 'task_help'), icon='help')
        self.actionDefinitions['references'] = dict(text="Bibliography", tip="Show bibliographic references for the task",
                                                    slot=functools.partial(self.handleHelpMenu, 'references'), icon='biblio')
        self.actionDefinitions['show_log'] = dict(text="Show log file", tip="Show complete log file for selected job",
                                                    slot=self.showLog, icon='showlogfile')
        self.actionDefinitions['new_project2'] = dict(text="New project", tip="Create a new project",
                                                    slot = functools.partial(self.handleProjectMenu,'new_project'), icon='new_project')
        self.actionDefinitions['show_i2run'] = dict(text="Show i2run command", tip="Show i2run command that would achieve same outcome",
                                                    slot=self.showI2run, icon='showlogfile')

    def setShortcuts(self):
        for name, key in [['taskmenu', QtCore.Qt.CTRL + QtCore.Qt.Key_M], ['nextproject', QtCore.Qt.CTRL + QtCore.Qt.Key_N]]:
            self.shortcuts[name] = QtWidgets.QShortcut(QtGui.QKeySequence(key),self)
            self.shortcuts[name].activated.connect(functools.partial(self.handleShortcut, name))

    @QtCore.Slot(str)
    def handleShortcut(self, mode):
        if mode == 'taskmenu':
            self.showTaskChooser(True)
        elif mode == 'nextproject':
            idx = self.Instances.index(self) + 1
            if idx >= len(self.Instances):
                idx = 0
            self.Instances[idx].show()
            self.Instances[idx].raise_()

    def projectWidget(self):
        return self._projectWidget

    @QtCore.Slot()
    def showHelpMenu(self):
        menu = self.taskFrame.buttons.button('help').menu()
        menu.clear()
        action = menu.addAction('Task description')
        action.setObjectName('task_help')
        action = menu.addAction('Bibliographic references')
        action.setObjectName('references')
        refFileList = TASKMANAGER().searchReferenceFile(name=self.taskFrame.openJob.taskName, drillDown=True)
        menu.findChild(QtWidgets.QAction, 'references').setEnabled((len(refFileList) > 0))
        helpFile = TASKMANAGER().searchHelpFile(name=self.taskFrame.openJob.taskName)
        menu.findChild(QtWidgets.QAction, 'task_help').setEnabled((helpFile is not None))
        # Add task-specific help
        if hasattr(self,"webviewToolBar"):
            self.webviewToolBar.setButtonEnabled('references',(len(refFileList) > 0))
            self.webviewToolBar.setButtonEnabled('task_help',(helpFile is not None))
        programHelpList = TASKMANAGER().getTaskAttribute(self.taskFrame.openJob.taskName, 'PROGRAMHELP', default=self.taskFrame.openJob.taskName)
        if programHelpList is not None:
            if not isinstance(programHelpList, list):
                programHelpList = [programHelpList]
            for programHelp in programHelpList:
                if programHelp.count('$CCP4I2'):
                    helpPath = CCP4Utils.getCCP4I2Dir() + programHelp[7:]
                elif programHelp.count('$CCP4'):
                    helpPath = CCP4Utils.getCCP4Dir() + programHelp[5:]
                else:
                    helpPath = os.path.join(CCP4Utils.getCCP4Dir(), 'html', programHelp + '.html')
                if os.path.exists(helpPath):
                    action = menu.addAction(os.path.splitext(os.path.split(helpPath)[-1])[0])
                    action.setData(helpPath)

    @QtCore.Slot(str)
    @QtCore.Slot('QObject')
    def handleHelpMenu(self, mode):
        if not isinstance(mode, str):
            mode = str(mode.objectName())
        from .CCP4WebBrowser import WEBBROWSER
        if mode == 'references':
            self.showBibliography(taskNameList=[self.taskFrame.openJob.taskName])
        elif mode == 'task_help':
            fileName = TASKMANAGER().searchHelpFile(name=self.taskFrame.openJob.taskName)
            WEBBROWSER().openFile(fileName)
        else:
            fileName = action.data().__str__()  # KJS : Obvious bug here. Needs sorting.
            WEBBROWSER().openFile(fileName)

    def reportAvailable(self):
        status = self.taskFrame.openJob.status in CCP4DbApi.FINISHED_JOB_STATUS or self.taskFrame.openJob.status in ['Running', 'Running remotely', 'Failed' ]
        return status

    def testReportAvailable(self):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        testReportList = glob.glob(os.path.join(PROJECTSMANAGER().db().getProjectInfo(projectId=self.taskFrame.openJob.projectId, mode='projectdirectory'), 'CCP4_test*.log'))
        return (len(testReportList) > 0)

    def redoReport(self):
        try:
            self.taskFrame.outputFrame.showOutput(redo=True, doReload=True)
        except CException as e:
            warningMessage(e, 'Creating report','Error creating report', parent=self)

    def openTask(self, taskName=None, jobId=None, followJobId=None):
        try:
            openJob = self.taskFrame.openTask(jobId=jobId, taskName=taskName, followJobId=followJobId)
            self.setSelectedJob(jobId=openJob.jobId)
        except:
            print('ERROR CProjectView.init opening jobId,taskName', jobId, taskName)

    def setSelectedJob(self,jobId):
        self._projectWidget.selectJob(jobId)

    def Exit(self):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().updateProject(projectId=self.taskFrame.openJob.projectId, key='lastaccess')
        self.taskFrame.saveStatus()

    @QtCore.Slot('QObject')
    @staticmethod
    def updateInstances(qobj):
        l = []
        for w in CProjectViewer.Instances:
            if CCP4Utils.isAlive(w):
                l.append(w)
        CProjectViewer.Instances = l

    def getProject(self):
        #The current project and job data is held in a COpenJob member of CTaskFrame
        return self.taskFrame.openJob.projectId

    def getOpenJobNumber(self):
        ''' Return the current job number '''
        return self.taskFrame.openJob.jobnumber 

    @QtCore.Slot(bool)
    def showTaskChooser(self, state=None):
        '''Toggle display of task menu'''
        if state is None:
            state = (self.rightStack.currentIndex() == 0)
        if state:
            self.rightStack.setCurrentIndex(1)
            self.taskChooser.taskTree.setFocus(QtCore.Qt.OtherFocusReason)
            if len(self.taskChooser.taskTree.selectedIndexes()) == 0:
                #Select the top one ...
                cols = self.taskChooser.taskTree.model().columnCount()-1
                sel = QtCore.QItemSelection(self.taskChooser.taskTree.model().index(0,0),self.taskChooser.taskTree.model().index(0,cols))
                self.taskChooser.taskTree.selectionModel().select(sel,QtCore.QItemSelectionModel.Select)
        else:
            self.rightStack.setCurrentIndex(0)

        if state:#   KJS - Removed this code until the bug is fixed.
            for item in ['run_remote','run', 'clone','view_coot','view_ccp4mg','export_mtz','show_log','show_i2run']:
                obj = self.findChild(QtWidgets.QAction, item)
                if obj is not None:
                    obj.setEnabled(False)
                if hasattr(self,"webviewToolBar"):
                    self.webviewToolBar.setButtonEnabled(item,False)
        else:
            rootIdx = QtCore.QModelIndex()
            modelRows = self._projectWidget.model().sourceModel().rowCount(rootIdx)
            if modelRows > 0:
                self.updateActionEnabled()

    @QtCore.Slot(str)
    def handleChooseTask(self, taskName):
        self.handleChooseTaskAndFollowFrom(taskName)

    @QtCore.Slot(str,str)
    def handleChooseTaskAndFollowFrom(self, taskName, data=None):
        '''Open a selected task'''
        self.taskChooser.taskTree.setEnabled(False)
        openJob = self.taskFrame.openTask(taskName=taskName)
        self.setSelectedJob(jobId=openJob.jobId)
        self.showTaskChooser(False)
        self.taskChooser.taskTree.setEnabled(True)

    @QtCore.Slot(str,str,str,bool)
    def handleJobListSelectionChange(self, fileId, jobId, pipelineJobId, detatch):
        '''Keep the task display in step with the selection on the project job list'''
        if detatch:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['status', 'taskname'])
            if jobInfo['status'] in CCP4DbApi.FINISHED_JOB_STATUS:
                self.openTaskMainWindow('output',jobId)
            elif jobInfo['status'] in ['Running', 'Running remotely'] and TASKMANAGER().hasReportDefinition(name=jobInfo['taskname'], jobStatus=jobInfo['status']):
                self.openTaskMainWindow('output', jobId)
            else:
                self.openTaskMainWindow('input', pipelineJobId)
        else:
            if self.taskFrame.openJob.jobId is None or jobId != str(self.taskFrame.openJob.jobId):
                # Open task interface for the top level pipeline
                try:
                    openJob = self.taskFrame.openTask(jobId=jobId)
                except CException as e:
                    warningMessage(e, 'Opening new job','Failed opening job', parent=self)
                self.showTaskChooser(False)
                self.deleteSpawnedTaskWindows(jobId)
            elif  jobId == str(self.taskFrame.openJob.jobId) and self.rightStack.currentIndex() == 1:
                self.showTaskChooser(False)

    @QtCore.Slot(int)
    def handleTaskFrameChanged(self, mode):
        '''If necessary create the task input frame when the task frame tab changed to input
          This is part of mechanism to ensure slow task frame drawing is only done when necessary'''
        if mode == CProjectViewer.INPUT_TAB and self.taskFrame.inputFrame.taskWidget is None:
            oJ = self.taskFrame.openJob
            defFile = TASKMANAGER().lookupDefFile(oJ.taskname, oJ.taskversion)
            from ..core import CCP4Container
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            container = CCP4Container.CContainer(parent=self.taskFrame, definitionFile=defFile, guiAdmin=True)
            paramsFile = PROJECTSMANAGER().makeFileName(jobId = oJ.jobId, mode='JOB_INPUT')
            if not os.path.exists(paramsFile): paramsFile = PROJECTSMANAGER().makeFileName(jobId = oJ.jobId,mode='PARAMS')
            container.loadDataFromXml(paramsFile)
            taskWidget = self.taskFrame.inputFrame.createTaskWidget(oJ.taskname, projectId=oJ.projectid, jobId=oJ.jobId,
                                                                    container=container, taskEditable=False)
            taskWidget.launchJobRequestSignal.connect(self.taskFrame.launchJobRequest)
        if mode == CProjectViewer.INPUT_TAB:
            self.splitterSizes = self.centralWidget().sizes()
            w = CCP4TaskWidget.WIDTH + CCP4TaskWidget.PADDING_ALLOWANCE
            if self.splitterSizes[1] < w:
                totWidth = self.splitterSizes[0] + self.splitterSizes[1]
                self.centralWidget().setSizes([totWidth - w, w])
                self.centralWidget().setStretchFactor(2, 0)
        self.lastTaskFrameMode = mode

    def showLog(self):
        jobId = str(self.taskFrame.openJob.jobId)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        jobDirectory = PROJECTSMANAGER().makeFileName(jobId=jobId,mode='ROOT')
        logfiles = []
        for root, subFolders, files in os.walk(jobDirectory):
            for fn in files:
                if fn == "log.txt" and os.path.exists(os.path.join(root,fn)):
                    fileName = os.path.join(root,fn)
                    if os.path.exists(fileName):
                        logfiles.append((fileName, os.path.getmtime(fileName)))
        logfiles.sort(key=lambda tup: tup[1]) 
        text = ""
        for f in logfiles:
            fh = open(f[0])
            text += "================================================================================\n"
            text += f[0] + "\n"
            text += "================================================================================\n"
            text += fh.read()
            fh.close()
        logfile = tempfile.NamedTemporaryFile(suffix='.txt',delete=False)
        logFileName = logfile.name
        logfile.write(bytes(text,"UTF-8"))
        logfile.close()
        from ..qtcore.CCP4Launcher import LAUNCHER
        LAUNCHER().launch('logview',[logFileName])

    def burrowOn(self, etreeElement, root):
        result = ""
        if len(etreeElement) == 0:
            subKeyword = "/".join(root.getpath(etreeElement).split("/")[2:])
            #print(etree.tostring(etreeElement))
            if etreeElement.text is not None and len(etreeElement.text) != 0:
                textToUse = etreeElement.text
                #Here apply some patching to make syntax a bit more compact:
                #print("Before",subKeyword)
                subKeyword=subKeyword.replace('/CPdbEnsembleItem','')
                subKeyword=subKeyword.replace('selection/item','selection')
                #print("After",subKeyword)
                if " " in textToUse:
                    keywordPair = '\t\t"{}={}" \\\n'.format(subKeyword, textToUse)
                    print("Warning: spaces in subKeyworded argument", keywordPair)
                    result += keywordPair
                else:
                    result += '\t\t"{}={}" \\\n'.format(subKeyword, textToUse)
        else:
            for subElement in etreeElement:
                result += self.burrowOn(subElement, root)
        return result

    def i2runForElement(self, element, objectName=None):
        result = ""
        from ..core import CCP4Data
        elementEtree = element.getEtree()
        if len(elementEtree) == 0:
            if elementEtree.text is not None and len(elementEtree.text)>0:
                if isinstance(element,(CCP4Data.CString,)) and hasattr(elementEtree, "text") and elementEtree.text is not None and " " in elementEtree.text:
                    result += '\t--{} "{}" \\\n'.format(objectName, elementEtree.text)
                else:
                    result += '\t--{} {} \\\n'.format(objectName, elementEtree.text)
        else:
            interimResult = ""
            for subElement in elementEtree:
                interimResult += self.burrowOn(subElement, etree.ElementTree(elementEtree))
            print("interimResult [{}]".format(interimResult))
            if len(interimResult) > 0:
                result += '\t--{} \\\n'.format(objectName, elementEtree.text)
                result += interimResult
        return result
                    
    def listContents(self, container, baseXml):
        from ..core import CCP4Container
        from ..core import CCP4Data
        #print("Contents of container", container)
        result = ""
        
        for child in container.children():
            #print(child.objectName())
            #print("Child", child.objectName(), child, child.isSet(), (not hasattr(child,'isDefault')) or ( child.isDefault()))
            if isinstance(child, (CCP4Container.CContainer,)):
                #print ("isContainer")
                if child.objectName() not in ["outputData", "guiAdmin", "guiControls", "patchSelection", "guiParameters", "temporary"]:
                    #print ("willList")
                    result += self.listContents(child, baseXml)
            else:
                #print("Not container")
                if (hasattr(child, "isSet") and child.isSet()) or (isinstance(child,(CCP4Data.CList,)) and len(child)>0):
                    if (not hasattr(child,'isDefault')) or (not child.isDefault()):
                        #print("not default")
                        nameToUse = child.objectName()
                        if len(baseXml.findall('.//content[@id="{}"]'.format(nameToUse))) > 1:
                            nameToUse='.'.join(child.objectPath().split('.')[1:])
                        #print('nameToUse', nameToUse)
                        if isinstance(child,(CCP4Data.CList,)):
                            for listElement in child:
                                result += self.i2runForElement(listElement, objectName=nameToUse)
                        else:
                            result += self.i2runForElement(child, objectName=nameToUse)
        return result
    
    def showI2run(self):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        taskName = PROJECTSMANAGER().db().getJobInfo(jobId=self.taskFrame.openJob.jobId, mode=['taskname'])
        projectName = PROJECTSMANAGER().db().getJobInfo(jobId=self.taskFrame.openJob.jobId, mode=['projectname'])['projectname']
        i2Runner = CCP4I2Runner.CI2Runner(['i2run', taskName])
        
        #Select the inut tab to ensure that container is loaded
        if self.taskFrame.inputFrame.taskWidget is None:
            self.handleTaskFrameChanged(mode=CProjectViewer.INPUT_TAB)
        #Following call might be useful, but might mess up input
        #self.taskFrame.inputFrame.taskWidget.container.removeUnsetListItems()
        
        #print(i2Runner.defXml,'i2Runner.defXml')
        result = "i2run {} --projectName {} \\\n".format(taskName, projectName) + self.listContents(self.taskFrame.inputFrame.taskWidget.container, i2Runner.defXml)

        d = QtWidgets.QDialog()
        okButton = QtWidgets.QPushButton("OK")

        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(okButton)

        vbox = QtWidgets.QVBoxLayout()
        textEdit = QtWidgets.QTextEdit(d)
        textEdit.setText(result)
        vbox.addWidget(textEdit)
        vbox.addLayout(hbox)
        
        d.setLayout(vbox)
        okButton.clicked.connect(d.close)
        d.setWindowTitle("Equivalent i2run command")
        d.setModal(True)
        d.setStyleSheet("QTextEdit {min-width: 600px; height: 600px;}");
        d.exec_()
        
    def exportMtz(self):
        '''Export a monster MTZ of the job output'''
        self.exportJobFile(['complete_mtz', 'MTZ file', 'application/CCP4-mtz'], str(self.taskFrame.openJob.jobId))
        str(self.taskFrame.openJob.jobId)

    @QtCore.Slot(object,str)
    def exportJobFile(self,exportInfo,jobId):
        '''Export a job (by export project mechanism) or export specific output file'''
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        from ..qtcore.CCP4CustomMimeTypes import MIMETYPESHANDLER
        if exportInfo == 'job':
            if getattr(PROJECTSMANAGER(), 'exportThread', None) is not None:
                QtWidgets.QMessageBox.warning(self, 'Export job', 'Please wait - another job currently being exported')
                return
            from . import CCP4FileBrowser
            self.browser = CCP4FileBrowser.CFileDialog(self, title='Save all job files',
                                                       filters=[MIMETYPESHANDLER().getMimeTypeInfo('application/CCP4-compressed-db', 'filter')],
                                                       defaultSuffix=MIMETYPESHANDLER().getMimeTypeInfo('application/CCP4-compressed-db', 'fileExtensions')[0],
                                                       fileMode=QtWidgets.QFileDialog.AnyFile)
            self.browser.selectFile.connect(functools.partial(self.doExportJob, jobId))
            self.browser.show()
        else:
            taskName = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode=['taskname'])
            fromFile = TASKMANAGER().exportJobFiles(jobId=jobId, taskName=taskName, mode=exportInfo[0])
            if fromFile is None or len(fromFile) == 0 and exportInfo[0] == 'complete_mtz':
                fromFile = PROJECTSMANAGER().exportJobMtzFile(jobId=jobId)
                if fromFile is None:
                    QtWidgets.QMessageBox.warning(self, self.windowTitle(), 'No export file for job')
                    return
            if fromFile is not None:
                from . import CCP4FileBrowser
                self.browser = CCP4FileBrowser.CFileDialog(self, title='Save ' + exportInfo[1],
                                                           filters=[MIMETYPESHANDLER().getMimeTypeInfo(exportInfo[2], 'filter')],
                                                           defaultSuffix=MIMETYPESHANDLER().getMimeTypeInfo(exportInfo[2], 'fileExtensions')[0],
                                                           fileMode=QtWidgets.QFileDialog.AnyFile)
                self.browser.selectFile.connect(functools.partial(self.doExportFile, jobId, fromFile))
                self.browser.rejected.connect(functools.partial(PROJECTSMANAGER().cleanupExportFiles, jobId))
                self.browser.show()
    
    @QtCore.Slot(str,str,str,bool)
    def doExportFile(self, jobId, fromFile, toFile, move=False):
        '''Export a file by copying to selected path'''
        if move:
            shutil.move(fromFile, toFile)
        else:
            shutil.copyfile(fromFile, toFile)

    @QtCore.Slot(str,str)
    def doExportJob(self,jobId,toFile):
        '''Export a job by project export'''
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        projectId = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode=['projectid'])
        PROJECTSMANAGER().compressProject(projectId, jobList=[jobId], excludeI2files=False, fileName=toFile)

    @QtCore.Slot(str,str,str)
    def purgeJobFiles(self, jobId, context='temporary'):
        from ..core import CCP4ProjectsManager
        purger = CCP4ProjectsManager.CPurgeProject(jobId=jobId, projectId=self.taskFrame.openJob.projectId)
        purger.purgeJob(jobId=jobId, context=context)

    @QtCore.Slot(str)
    def cloneTask(self, oldJobId=None,suggestedParams=None):
        '''Clone a job'''
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("cloneTask",self.taskFrame.openJob.jobId,oldJobId)
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        if not oldJobId:
            try:
                if PROJECTSMANAGER().db().getJobInfo(jobId=self.taskFrame.openJob.jobId, mode='status') == 'Pending':
                    self.taskFrame.inputFrame.makeJobInputFile(self.taskFrame.openJob.jobId)
            except:
                pass
            oldJobId = self.taskFrame.openJob.jobId
        clonedJobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=oldJobId, mode=['status', 'taskname'])
        if oldJobId is None or clonedJobInfo['taskname'] is None: return
        if clonedJobInfo['status'] == 'Pending':
            #Do we have the Pending job open in a taskFrame - if so then save the params
            for taskWidget in self.findChildren(CCP4TaskWidget.CTaskWidget):
                if taskWidget._jobId == oldJobId:
                    taskWidget.saveToXml()
        try:
            openJob =  self.taskFrame.openTask(taskName=clonedJobInfo['taskname'], cloneJobId=oldJobId, suggestedParams=suggestedParams)
        except CException as e:
            warningMessage(e, 'Opening new job', 'Failed opening ' + str(clonedJobInfo['taskname']) + ' job', parent=self)
        if openJob is not None:
            self.setSelectedJob(jobId=openJob.jobId)
            self.showTaskChooser(False)

    @QtCore.Slot(bool,bool,str)
    def stopJob(self, kill=False, delete=False, jobId=None):
        '''Stop a running job'''
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        taskName = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode='taskname')
        if kill and taskName == 'coot_rebuild':
            QtWidgets.QMessageBox.warning(self, self.windowTitle(), 'Please close Coot from within Coot')
            return
        if kill:
            from ..qtcore.CCP4JobController import JOBCONTROLLER
            err = JOBCONTROLLER().killJobProcess(jobId=jobId)
            if delete:
                try:
                    self.deleteJob(jobIdList=[jobId])
                except CException as e:
                    err.extend(e)
                    warningMessage(err, 'Deleting job', 'Failed deleting job', parent=self)
            else:
                PROJECTSMANAGER().db().updateJobStatus(jobId, CCP4DbApi.JOB_STATUS_FAILED)
        else:
            jobDir = PROJECTSMANAGER().jobDirectory(jobId=jobId, create=False)
            if jobDir is not None:
                CCP4Utils.saveFile( os.path.join(jobDir, 'INTERRUPT'), 'Signal to interrupt job')

    @QtCore.Slot(bool,str)
    def markFailedJob(self, succeeded=False, jobId=None):
        '''Set job status to finished or failed'''
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        if succeeded:
            PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status = CCP4DbApi.JOB_STATUS_FINISHED)
        else:
            PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status = CCP4DbApi.JOB_STATUS_FAILED)

    @QtCore.Slot(dict)
    def handleJobDeleted(self,args):
        '''Update gui on jobDeleted signal from database'''
        if args['jobId'] == self.taskFrame.openJob.jobId:
            newJobId = self.suggestOpenJob()
            if newJobId is not None:
                openJob = self.taskFrame.openTask(jobId=self.suggestOpenJob())
                self.setSelectedJob(jobId=openJob.jobId)
            else:
                pass
        self.deleteSpawnedTaskWindows(args['jobId'])
        
    def deleteSpawnedTaskWindows(self, jobId):
        '''Delete the popout task windows'''
        children = self.findChildren(CTaskMainWindow, 'job' + str(jobId))
        for child in children:
            child.close()
            child.deleteLater()

    def suggestOpenJob(self, lastOpenJobId=None):
        # Suggest an open job when the current one has been deleted. This only get called if the last job in list 
        # (usually the first job that was run) is delete because the Qt widget just moves the selection to the next job down list
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        jobIdList = PROJECTSMANAGER().db().getProjectJobListInfo(projectId=self.taskFrame.openJob.projectId, mode='JobId', topLevelOnly=True, order='DESC')
        if len(jobIdList) > 0:
            return jobIdList[0]['jobid']
        else:
            return None

    @QtCore.Slot(dict)
    def handleProjectDeleted(self,args):
        '''Close window if project has been deleted'''
        if args['projectId'] == self.taskFrame.openJob.projectId:
            self._projectWidget.setForceUpdate(False)
            self.close()

    def showHelp(self, mode='ccp4i2',newTab=True):
        '''Show help pages'''
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().showHelp(mode=mode,newTab=newTab)

    @QtCore.Slot(str,str,str)
    def handleNext(self, nextTask, jobId, patchParamsFile):
        '''Handle user selection of next task (in Job list) by creating new job'''
        openJob = self.taskFrame.openTask(taskName=nextTask, followJobId=jobId, patchParamsFile=patchParamsFile)
        self.setSelectedJob(jobId=openJob.jobId)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().updateJob(jobId=openJob.jobId, key='preceedingjobid', value=jobId)

    @QtCore.Slot(str,str)
    def handleNextTask(self,taskName,patchParamsFile):
        '''Handle Next task selection from buttons under report in CTaskFrame'''
        if taskName == 'clone_rerun':
            self.cloneTask()
        elif taskName == 'clone_suggested':
                from ..core.CCP4ProjectsManager import PROJECTSMANAGER
                outputXmlFile = PROJECTSMANAGER().makeFileName(jobId=self.taskFrame.openJob.jobId,mode='PROGRAMXML')
                parser = etree.XMLParser()
                f = open(outputXmlFile)
                s = f.read()
                f.close()
                tree = etree.fromstring(s, parser)
                suggested = tree.xpath("Verdict/suggestedParameters")
                self.cloneTask(suggestedParams=suggested[0])
        else:
            openJob = self.taskFrame.openTask(taskName=taskName, followJobId=self.taskFrame.openJob.jobId, patchParamsFile=patchParamsFile)
            self.setSelectedJob(jobId=openJob.jobId)

    @QtCore.Slot(str,str,str)
    def showWhatNextPage(self, jobId, taskName=None, projectId=None):
        '''Show a What next page'''
        nextPage = TASKMANAGER().getWhatNextPage(taskName=taskName, jobId=jobId)
        if nextPage is None:
            taskTitle= TASKMANAGER().getTitle(taskName)
            QtWidgets.QMessageBox.warning(self, self.windowTitle(), 'Sorry - no What Next? page for ' + str(taskTitle))
        else:
            from .CCP4WebBrowser import WEBBROWSER
            view = WEBBROWSER().loadWebPage(helpFileName=nextPage)
            if view is not None:
                if projectId is None:
                    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
                    projectId = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode='projectid')
                view.report.setJobInfo(jobInfo= {'jobId' : jobId, 'projectId' : projectId})

    @QtCore.Slot(str,list)
    def showBibliography(self, jobId=None, taskNameList=[]):
        '''Show bibliography for task'''
        from . import CCP4BibliographyViewer
        if jobId is None:
            jobId = self.taskFrame.openJob.jobId
        viewer = CCP4BibliographyViewer.CBibliographyViewer(self)
        viewer.setReferences(jobId=jobId, taskNameList=taskNameList)
        viewer.show()

    def openApplication(self,application):
        '''Pass on request for program logs'''
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().openApplication(application)
    
    def launchJobRequest(self,taskName,args):
        try:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            jobId, pName, jNumber = PROJECTSMANAGER().newJob(taskName=taskName, projectId=self.taskFrame.openJob.projectId)
        except:
            QtWidgets.QMessageBox.warning(self, 'Error in creating new job','Unknown error creating new job' + str(taskName) + '\n')
            return
        # Create an input params file for job
        from ..core import CCP4Container
        container = CCP4Container.CContainer(definitionFile=TASKMANAGER().lookupDefFile(taskName), guiAdmin=True)
        # Load known data to the params file
        if taskName == 'mergeMtz' and args.get('fileName', None) is not None:
            container.inputData.MINIMTZINLIST[0].set({'fileName' : args.get('fileName', None)})
        container.saveDataToXml(fileName=PROJECTSMANAGER().makeFileName(jobId=jobId, mode='JOB_INPUT'))
        # If there is a jobId for a 'parent' job that launched this then it needs a call to CTaskWidget.handleLaunchedJobFinish
        # when the 'child' job finishes
        if args.get('launchJobId', None) is not None:
            # use non-existance of _launchedPopoutJobs attribute as flag to
            if getattr(self, '_launchedPopoutJobs', None) is None:
                self._launchedPopoutJobs = {jobId : args['launchJobId']}
                PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)
        # Open a a popout window
        win = self.openTaskMainWindow('input', jobId)
        # Call 'parent' job to inform of progress and pass over reference to child widget
        childTaskWidget = None
        parentTaskWidget = None
        for taskWidget in self.findChildren(CCP4TaskWidget.CTaskWidget):
            if taskWidget.jobId() == jobId:
                childTaskWidget = taskWidget
            elif taskWidget.jobId() == args['launchJobId']:
                parentTaskWidget = taskWidget
        if parentTaskWidget is not None:
            parentTaskWidget.handleLaunchedJob(jobId=jobId, status=1, taskWidget=childTaskWidget)
    
    def handleJobFinished(self, args):
        '''Update a popout window when job status changed'''
        launchJobId = self._launchedPopoutJobs.get(args['jobId'], None)
        if launchJobId is not None:
            for taskWidget in self.findChildren(CCP4TaskWidget.CTaskWidget):
                if taskWidget.jobId() == launchJobId:
                    taskWidget.handleLaunchedJob(jobId=args['jobId'], status=args.get('status',None))
                    del self._launchedPopoutJobs[args['jobId']]
    
    @QtCore.Slot(str,str)
    def openTaskMainWindow(self, mode, jobId=None):
        '''Open a popout window containing task input or report'''
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        openJob = COpenJob(jobId=jobId)
        if mode == 'input':
            paramsFile = PROJECTSMANAGER().makeFileName(jobId = jobId, mode='JOB_INPUT')
            if not os.path.exists(paramsFile):
                return None
            defFile = TASKMANAGER().lookupDefFile(openJob.taskname, openJob.taskversion)
            if defFile is None:
                return None
            from ..core import CCP4Container
            container = CCP4Container.CContainer(parent=self, definitionFile=defFile, guiAdmin=True)
            container.loadDataFromXml(paramsFile)
            taskEditable = (openJob.status in ['Unknown', 'Pending'])
            try:
                widget = CTaskInputFrame(self)
                widget.createTaskWidget(openJob.taskname, projectId=openJob.projectid, jobId=jobId, container=container, taskEditable=taskEditable)
            except CException as e:
                warningMessage(e, parent=self)
                return None
            except Exception as e:
                QtWidgets.QMessageBox.warning(self, 'Error in creating task input','Unknown error creating task input widget for job number '+str(openJob.jobnumber)+str(e))
        elif mode == 'output':
            if 1:
                widget = CReportView(self)
                widget.showOutput(openJob=openJob)
        elif mode == 'status':
            widget = CJobStatusWidget(self)
            widget.setJob(openJob)
        win = CTaskMainWindow(None,openJob.projectname,openJob.jobId)
        win.buttons.setMode(mode)
        if win.buttons.button('clone') is not None:
            win.buttons.button('clone').clicked.connect(functools.partial(self.cloneTask, jobId))
            win.buttons.button('view').menu().triggered.connect(functools.partial(self.handleViewTask0, jobId))
        if mode == 'input':
            win.buttons.button('run').clicked.connect(win.runTask)

        @QtCore.Slot(dict)
        def handleJobStatusUpdated(status):
            if "projectId" in status and "jobId" in status and "parentJobId" in status and status["projectId"] == openJob.projectid and (status["jobId"] == openJob.jobId or status["parentJobId"] == openJob.jobId):
                if status["status"] == CCP4DbApi.JOB_STATUS_RUNNING:
                    win.buttons.button('run').setEnabled(False)
                    win.centralWidget().setEnabled(False)

        PROJECTSMANAGER().db().jobStarted.connect(handleJobStatusUpdated)

        win.titleBar.setOpenJob(openJob)
        win.buttons.setEnabled(openJob.status)
        win.centralWidget().layout().insertWidget(1, widget)
        win.show()
        win.raise_()
        self.mainwins.append(win)
        return win

    @QtCore.Slot(str,'QAction')
    def handleViewTask0(self,jobId, action=None):
        '''Show job files in ccp4mg or coot'''
        text = str(action.text())
        from ..qtcore.CCP4Launcher import LAUNCHER
        if text.count('4mg'):
            LAUNCHER().openInViewer(viewer='ccp4mg', jobId=jobId, projectId=self.taskFrame.openJob.projectId, guiParent=self)
        elif text.count('oot'):
            LAUNCHER().openInViewer(viewer='coot_job', jobId=jobId, projectId=self.taskFrame.openJob.projectId, guiParent=self)

    def openSendReport(self):
        '''Open window to send developer error report'''
        from . import CCP4ErrorReportViewer
        widget = CCP4ErrorReportViewer.CSendJobError(self, projectId=self.taskFrame.openJob.projectId, projectName=self.taskFrame.openJob.projectName)
        widget.show()

    @QtCore.Slot(list,bool)
    def deleteJob(self, jobIdList, deleteImportFiles=True):
        '''Delete one or more jobs'''
        followOnJobs = []
        allJobTree = []
        if not isinstance(jobIdList, list):
            jobIdList = [jobIdList]
        try:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            jobTreeList = PROJECTSMANAGER().db().getMultiFollowOnJobs(jobIdList=jobIdList, traceImportFiles=deleteImportFiles)
        except CException as e:
            warningMessage(e, parent=self, windowTitle=self.windowTitle(), message='Error attempting to find jobs dependent on output files')
            return
        xtrJobIdList = []
        for  jobTree in jobTreeList:
            for followOnJob in jobTree[2]:
                if followOnJob[0] not in jobIdList:
                    xtrJobIdList.append(followOnJob[0])
        #Are there open pending jobs that might use files?
        jobsToDeleteWithSelectedFiles = []
        if self.taskFrame.openJob.status == 'Pending':
            selectedFileDict = self.taskFrame.inputFrame.taskWidget.container.inputData.inputFilesFileIds()
            if len(selectedFileDict) > 0:
                jobsWithSelectedFiles = PROJECTSMANAGER().db().getJobSourceOfFiles(fileIdList=list(selectedFileDict.keys()))
                jobsToDeleteWithSelectedFilesSet = set(jobIdList+xtrJobIdList) & set(jobsWithSelectedFiles)
                for i in range(len(jobsToDeleteWithSelectedFilesSet)):
                    jobsToDeleteWithSelectedFiles.append(jobsToDeleteWithSelectedFilesSet.pop())
        if len(xtrJobIdList) == 0 and len(jobsToDeleteWithSelectedFiles) == 0:
            # Delete in reverse order
            for iJob in range(len(jobTreeList) - 1, -1, -1):
                delJobId, importFiles, followOnJobs = jobTreeList[iJob]
                try:
                    PROJECTSMANAGER().deleteJob(jobId=delJobId, importFiles=importFiles, projectId=self.taskFrame.openJob.projectId,
                                                deleteImportFiles=deleteImportFiles)
                except CException as e:
                    warningMessage(e, parent=self, windowTitle=self.windowTitle(), message='Error attempting to delete job')
                    return
                except Exception as e:
                    QtWidgets.QMessageBox.warning(self, self.windowTitle(), 'Unrecognised error deleting job\n'+str(e))
                    return
                self.handleJobsDeleted()
        else:
            #QtWidgets.QMessageBox.warning(self,self.windowTitle(),"Deleting multiple jobs with 'knock on' jobs temporarily disabled\n")
            #return
            deleteJobGui = CDeleteJobGui(self, self.taskFrame.openJob.projectId, jobIdList=jobIdList, jobTreeList=jobTreeList,
                                         jobsToDeleteWithSelectedFiles=jobsToDeleteWithSelectedFiles, deleteImportFiles=deleteImportFiles, ifXtrJobs=len(xtrJobIdList) > 0)
            deleteJobGui.jobsDeleted.connect(self.handleJobsDeleted)
            deleteJobGui.show()

    @QtCore.Slot()
    def handleJobsDeleted(self):
        '''When job deleted reset the CDataFile 'job' combo to remove deleted job'''
        if self.taskFrame.inputFrame.taskWidget is not None:
            self.taskFrame.inputFrame.taskWidget.resetJobCombos()
            self.taskFrame.inputFrame.taskWidget.validate()
        rootIdx = QtCore.QModelIndex()
        modelRows = self._projectWidget.model().sourceModel().rowCount(rootIdx)
        if modelRows == 0:
            self.showTaskChooser(True)

    def handleProjectMenuExport(self):
        if CCP4WebBrowser.CMainWindow.projectManagerDialog is None:
            from . import CCP4ProjectManagerGui
            CCP4WebBrowser.CMainWindow.projectManagerDialog = CCP4ProjectManagerGui.CProjectManagerDialog()
            CCP4WebBrowser.CMainWindow.projectManagerDialog.hide()    # KJS : Problem here. Non-existent functions by the looks of it.
        CCP4WebBrowser.CMainWindow.projectManagerDialog.handleExport3(self.taskFrame.openJob.projectId)

    def openManageImportFiles(self):
        '''Open window showing imported files for this project'''
        from . import CCP4ImpExpWidgets
        widget = CCP4ImpExpWidgets.CManageImportFiles(self, projectId=self.taskFrame.openJob.projectId)
        widget.show()

    @QtCore.Slot(str,str)
    def handleCsvDownloadRequest(self, jobId, dataName):
        '''Export a comma-separated-variables file containing data from a report table'''
        if jobId is None:
            jobId = self.taskFrame.openJob.jobId
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        fileName = os.path.join(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='TABLES_DIR'), dataName + '.csv')
        if not os.path.exists(fileName):
            fileName = os.path.join(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='TABLES_DIR'), dataName + '.csv').replace(" ","_")
            if not os.path.exists(fileName):
                return
        from . import CCP4FileBrowser
        self.fileBrowser = CCP4FileBrowser.CFileDialog(self, title='Save csv table file for ' + dataName,
                                                       filters=['csv table file (*.csv)'], defaultSuffix='csv',
                                                       fileMode=QtWidgets.QFileDialog.AnyFile)
        self.fileBrowser.selectFile.connect(functools.partial(self.downloadCsvFile, jobId,dataName))
        self.fileBrowser.show()
    
    @QtCore.Slot(str,str,str)
    def downloadCsvFile(self, jobId, dataName, targetFile):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        sourceFile = os.path.join(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='TABLES_DIR'), dataName + '.csv')
        if not os.path.exists(sourceFile):
            sourceFile = os.path.join(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='TABLES_DIR'), dataName + '.csv').replace(" ","_")
        try:
            shutil.copyfile(sourceFile, targetFile)
        except:
            QtWidgets.QMessageBox.warning(self, self.windowTitle(), 'Failed to copy csv table file to ' + str(targetFile))

    def grabWidget(self):
        '''Developer tool to grab image of task input or report widget'''
        if not hasattr(self, 'grabStatus'):
            self.grabStatus = {'taskName' : None, 'saveDir' : None, 'mode' : None}
        if self.taskFrame.tab.currentIndex() == CTaskFrame.INPUT_TAB:
            wList = self.findChildren(CCP4TaskWidget.CTaskWidget)
            if len(wList) == 0:
                QtWidgets.QMessageBox.warning(self, 'Grab widget', 'Failed to find child task widget')
                return
            taskWidget = wList[0]
        elif self.taskFrame.tab.currentIndex() == CTaskFrame.OUTPUT_TAB:
            wList = self.findChildren(CReportView)
            if len(wList) == 0:
                QtWidgets.QMessageBox.warning(self, 'Grab widget','Failed to find child report view')
                return
            reportWidget = wList[0]
        else:
            return
        taskName = self.taskFrame.openJob.taskName
        if self.grabStatus['taskName'] is None or taskName != self.grabStatus['taskName']:
            rv = QtWidgets.QFileDialog.getExistingDirectory(caption='Select directory for saving screen grabs for task ' + str(taskName))
            if rv is None or len(rv) == 0:
                return
            self.grabStatus['saveDir'] = str(rv)
            self.grabStatus['taskName'] = taskName
        if self.taskFrame.tab.currentIndex() == CTaskFrame.INPUT_TAB:
            self.grabStatus['mode'] = 'task'
        else:
            self.grabStatus['mode'] = 'report'
        fileList = glob.glob(os.path.join(self.grabStatus['saveDir'], self.grabStatus['taskName'] + '_'+self.grabStatus['mode'] + '*.png'))
        if len(fileList) == 0:
            fileName = os.path.join(self.grabStatus['saveDir'], self.grabStatus['taskName'] + '_' + self.grabStatus['mode'] + '_1.png')
        else:
            fileList.sort()
            num = int(fileList[-1][-6:-4].strip('_'))
            fileName = os.path.join(self.grabStatus['saveDir'], self.grabStatus['taskName'] + '_' + self.grabStatus['mode'] + '_' + str(num+1) + '.png')
        pix0 = QtGui.QPixmap()
        if self.grabStatus['mode'] == 'task':
            pix1 = pix0.grabWidget(taskWidget.widget)
        else:
            pix1 = pix0.grabWidget(reportWidget)
        f = QtCore.QFile(fileName)
        f.open(QtCore.QIODevice.WriteOnly)
        pix1.save(f, "PNG")
        f.close()
        print('Saved to ', fileName)
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().openFile(fileName)

    def queryDeleteJob(self, jobId, jobInfo):
        self.queryDeleteWindow = QtWidgets.QDialog(self)
        self.queryDeleteWindow.setWindowTitle(self.windowTitle())
        self.queryDeleteWindow.setModal(True)
        self.queryDeleteWindow.setLayout(QtWidgets.QVBoxLayout())
        self.queryDeleteWindow.layout().addWidget(QtWidgets.QLabel('''By default interactive jobs (such as ''' + TASKMANAGER().getShortTitle(jobInfo['taskname']) + \
        ''') with no output files are deleted.\nYou can change this in the Preferences window.\nJob number ''' + str(jobInfo.get('jobnumber',' ')) + \
        '''will be deleted.''', self))
        self.queryDelete = QtWidgets.QCheckBox('Delete interactive jobs with no output files', self.queryDeleteWindow)
        from ..core.CCP4Preferences import PREFERENCES
        self.queryDelete.setChecked(bool(PREFERENCES().DELETE_INTERACTIVE_JOBS))
        self.queryShow = QtWidgets.QCheckBox('Do not show this warning again', self.queryDeleteWindow)
        self.queryShow.setChecked(True)
        self.queryDeleteWindow.layout().addWidget(self.queryDelete)
        self.queryDeleteWindow.layout().addWidget(self.queryShow)
        butBox = QtWidgets.QDialogButtonBox(self.queryDeleteWindow)
        butBox.addButton(QtWidgets.QDialogButtonBox.Ok)
        self.queryDeleteWindow.layout().addWidget(butBox)
        butBox.buttons()[0].clicked.connect(functools.partial(self.handleQueryDeleteJob, jobId))
        self.queryDeleteWindow.show()
        self.queryDeleteWindow.raise_()

    @QtCore.Slot(str)
    def handleQueryDeleteJob(self, jobId):
        ifShow = not(self.queryShow.isChecked())
        ifDelete = self.queryDelete.isChecked()
        self.queryDeleteWindow.hide()
        self.queryDeleteWindow.deleteLater()
        from ..core.CCP4Preferences import PREFERENCES
        PREFERENCES().DELETE_INTERACTIVE_JOBS.set(ifDelete)
        PREFERENCES().SHOW_DELETE_INTERACTIVE_JOBS.set(ifShow)
        PREFERENCES().save()
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        if ifDelete:
            try:
                PROJECTSMANAGER().db().deleteJob(jobId=jobId, deleteChildren=True)
            except:
                print('Projectviewer handleQueryDeleteJob failed to delete job')
        else:
            PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FINISHED)
            # Function to return list of names of exportable MTZ(s)

    def checkForRemotePasswords(self):
        '''Request remote machine password in order to check job status
          This is usually only relevant when i2 restarted after closing down
          with remote jobs running'''
        # i2 does not save user password on remote machines - check if we need a password in order to check the remote machine
        from ..qtcore.CCP4JobController import JOBCONTROLLER
        requirePass = JOBCONTROLLER().checkForRemotePasswords(self.taskFrame.openJob.projectId)
        if len(requirePass) == 0:
            return
        from . import CCP4JobControlGui
        self.passwordEntry = CCP4JobControlGui.CPasswordEntry(parent=self, label='Please enter password for ' + requirePass[0]['username'] + ' on ' + requirePass[0]['machine'] + '\nTo enable recovering remote jobs')
        self.passwordEntry.passwordEntered.connect(functools.partial(self.handlePasswordEntry,requirePass[0]['jobId']))
        self.passwordEntry.show()
    
    @QtCore.Slot(str,str)
    def handlePasswordEntry(self,jobId,password):
        from ..qtcore.CCP4JobController import JOBCONTROLLER
        JOBCONTROLLER().patchRemotePasswords(jobId,password)

    def widgetIsSaveable(self):
        return False

    def widgetIsPrintable(self):
        return False

    def widgetIsSearchable(self):
        return False

    def widgetIsRunable(self):
        return False

    def handleSave(self):
        pass

    def handlePrint(self):
        pass

    def handleRun(self):
        pass

    def openFind(self):
        pass

    def isFindFrameOpen(self):
        return False

    def deleteTab(self):
        pass

    def historyBack(self):
        pass

    def historyForward(self):
        pass

    def reloadPage(self):
        pass


class CTaskButtonsLayout(QtWidgets.QHBoxLayout):
    def maximumSize(self):
        return QtCore.QSize(CCP4TaskWidget.WIDTH, 50)

    def sizeHint(self):
        return self.maximumSize()


class CTaskButtons(QtWidgets.QButtonGroup):
    '''Buttons now only used in the popout task window'''
    BUTTONS = ['task_menu', 'help', 'view', 'clone', 'run']
    BUTTONTEXT = ['Task menu', 'Help', 'View', 'Clone job', 'Run']
    TOOLTIPS = ['Open another task from the task chooser', 'Open task and program documentation', 'View the output maps and model in molecular graphics',
                'Create another job with the same parameters set', 'Check the input data then run the job']
    MOREINFO = 'More info..'
    DEFAULTMODE = 0
    RUNONLYMODE = 1

    #buttons
    def __init__(self, parent, parentLayout=None, mode=DEFAULTMODE):
        QtWidgets.QButtonGroup.__init__(self,parent)
        self.mode = 'task'  # task/input/output/status
        self.frame = QtWidgets.QFrame(parent)
        self.frame.setObjectName('buttons')
        layout = CTaskButtonsLayout()
        self.frame.setLayout(layout)
        layout.setContentsMargins(CProjectViewer.MARGIN,CProjectViewer.MARGIN,CProjectViewer.MARGIN,CProjectViewer.MARGIN)
        layout.setSpacing(0)
        layout.addStretch(1)
        if mode == CTaskButtons.RUNONLYMODE:
            butRange = [4]
        else:
            butRange = list(range(len(CTaskButtons.BUTTONTEXT)))
        for ii in butRange:
            self.addButton(QtWidgets.QPushButton(CTaskButtons.BUTTONTEXT[ii],parent),ii)
            but = self.button(CTaskButtons.BUTTONS[ii])
            but.setToolTip(CTaskButtons.TOOLTIPS[ii])
            but.setAutoDefault(False)
            layout.addWidget(but)
            layout.addStretch(1)
        #parentLayout.addLayout(layout)
        parentLayout.addWidget(self.frame)
        self.button('run').setDefault(True)
        if mode==CTaskButtons.DEFAULTMODE:
            helpMenu = QtWidgets.QMenu(parent)
            self.button('help').setMenu(helpMenu)
            viewMenu = QtWidgets.QMenu(parent)
            self.button('view').setMenu(viewMenu)
            #nextMenu = QtWidgets.QMenu(parent)
            #self.button('next').setMenu(nextMenu)
            action = viewMenu.addAction('In CCP4mg')
            action = viewMenu.addAction('In Coot')
        if ALWAYS_SHOW_SERVER_BUTTON:
            runMenu =  QtWidgets.QMenu(parent)
            self.button('run').setMenu(runMenu)
            action = runMenu.addAction('Now')
            action = runMenu.addAction('Remotely - test output')
            action = runMenu.addAction('Remotely - import result')

    def setMode(self,mode):
        self.mode = mode

    def button(self,name):
        ii = CTaskButtons.BUTTONS.index(name)
        return QtWidgets.QButtonGroup.button(self, ii)

    def setRunMode(self, editor=None, status='Pending'):
        if editor:
            self.button('run').setText('Save')
        else:
            if status == 'Interrupted':
                self.button('run').setText('Restart')
            else:
                self.button('run').setText('Run')
        
    def setEnabled(self, jobStatus):
        if jobStatus is None:
            #No job set - disable all buttons
            for item in ['run', 'view', 'clone']:
                b = self.button(item)
                if b is not None:
                    b.setEnabled(False)
        else:
            if self.mode in ['task', 'input']:
                self.button('run').setEnabled(jobStatus in ['Pending', 'Interrupted'])
            else:
                self.button('run').setEnabled(False)
            if self.button('view') is not None:
                self.button('view').setEnabled(jobStatus in ['Finished', 'Interrupted', 'To delete'])
                self.button('clone').setEnabled(True)


class CJobStatusWidget(QtWidgets.QFrame):

    MARGIN = 1

    def __init__(self, parent=None):
        QtWidgets.QFrame.__init__(self, parent=None)
        # Save the parent viewer - the Qt widget is reparented to a StackedWidget when it is put in a tab
        self.openJob = COpenJob()
        self.commentId = None
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(CJobStatusWidget.MARGIN, CJobStatusWidget.MARGIN, CJobStatusWidget.MARGIN, CJobStatusWidget.MARGIN)
        self.layout().setSpacing(CJobStatusWidget.MARGIN)
        splitter = QtWidgets.QSplitter(self)
        splitter.setOrientation(QtCore.Qt.Vertical)
        splitter.setSizes([2, 1])
        self.layout().addWidget(splitter)
        # Annotation frame
        frame =  QtWidgets.QFrame(self)
        layout = QtWidgets.QVBoxLayout()
        frame.setLayout(layout)
        line =  QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Rate this job:', self))
        self.thunk = QtWidgets.QComboBox(self)
        self.thunk.setEditable(False)
        for item in CCP4DbApi.JOB_EVALUATION_TEXT:
            self.thunk.addItem(CCP4ProjectWidget.jobIcon(item), item)
        line.addWidget(self.thunk)
        line.addWidget(QtWidgets.QLabel('or change job title..', self))
        line.addStretch(1)
        layout.addLayout(line)
        self.title = QtWidgets.QLineEdit(self)
        self.title.setMaxLength(255)
        layout.addWidget(self.title)
        line =  QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Add comment..', self))
        line.addStretch(1)
        layout.addLayout(line)
        self.annotation = CCommentTextEdit(self)
        layout.addWidget(self.annotation)
        line =  QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Previous comments..', self))
        line.addStretch(1)
        layout.addLayout(line)
        self.history = QtWidgets.QTextEdit(self)
        self.history.setReadOnly(True)
        layout.addWidget(self.history)
        splitter.addWidget(frame) 
        self.thunk.currentIndexChanged[int].connect(self.saveThunk)
        self.title.editingFinished.connect(self.saveTitle)
        self.annotation.newLine.connect(self.saveAnnotation)

    def setJob(self,openJob):
        self.thunk.blockSignals(True)
        self.title.blockSignals(True)
        self.history.blockSignals(True)
        self.annotation.blockSignals(True)
        self.saveAnnotation()
        self.openJob= openJob
        self.commentId = None
        if self.openJob.jobId is None:
            self.thunk.setCurrentIndex(0)
            self.annotation.setPlainText('')
            self.title.setText('')
            self.history.setReadOnly(False)
            self.history.clear()
            self.history.setReadOnly(True)
        else:
            self.thunk.setCurrentIndex(CCP4DbApi.JOB_EVALUATION_TEXT.index(self.openJob.evaluation))
            if self.openJob.jobtitle is not None:
                self.title.setText(self.openJob.jobtitle)
            else:
                self.title.setText('')
            self.annotation.setPlainText('')
            self.history.setReadOnly(False)
            self.history.clear()
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            db = PROJECTSMANAGER().db()
            commentList = db.getComments(jobId=self.openJob.jobId)
            text = '<html><body>'
            for cid,user,time,comment in commentList:
                if db.timeIsToday(time):
                    self.commentId = cid
                    self.annotation.setPlainText(comment)
                else:
                    strTime = db.timeString(time)
                    text = text + '<h5>' + user + ' on ' + strTime+'</h5>\n<p>' + comment + '</p>\n'
            text = text + '</body></html>'
            self.history.setHtml(text)
            self.history.setReadOnly(True)
        self.thunk.blockSignals(False)
        self.title.blockSignals(False)
        self.history.blockSignals(False)
        self.annotation.blockSignals(False)
        self.openJob.jobUpdated.connect(self.handleJobUpdated)

    @QtCore.Slot(int)
    def saveThunk(self, indx):
        if self.openJob.jobId is None:
            return
        try:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            PROJECTSMANAGER().db().updateJob(jobId=self.openJob.jobId, key='evaluation', value=indx)
        except:
            pass

    @QtCore.Slot()
    def saveTitle(self):
        if self.openJob.jobId is None:
            return
        try:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            text = str(self.title.text())
            PROJECTSMANAGER().db().updateJob(jobId=self.openJob.jobId, key='jobTitle', value=text)
            projectId=self.openJob.projectId
            projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=projectId)
            dbxml = os.path.join(projectDir,"DATABASE.db.xml")
            print("Saving",dbxml)
            PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml)
        except:
            pass

    @QtCore.Slot()
    def saveAnnotation(self):
        if self.openJob.jobId is None: return
        try:
            text = str(self.annotation.toPlainText())
            if len(text) > 0:
                from ..core.CCP4ProjectsManager import PROJECTSMANAGER
                if self.commentId is None:
                    self.commentId = PROJECTSMANAGER().db().createComment(jobId=self.openJob.jobId,comment=text)
                else:
                    PROJECTSMANAGER().db().updateComment(commentId=self.commentId, comment=text)
                projectId=self.openJob.projectId
                projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=projectId)
                dbxml = os.path.join(projectDir,"DATABASE.db.xml")
                print("Saving",dbxml)
                PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml)
        except:
            pass

    @QtCore.Slot(str,str)
    def handleJobUpdated(self,key,value):
        if key == 'evaluation':
            self.thunk.blockSignals(True)
            self.thunk.setCurrentIndex(CCP4DbApi.JOB_EVALUATION_TEXT.index(self.openJob.evaluation))
            self.thunk.blockSignals(False)


class CCommentTextEdit(QtWidgets.QTextEdit):

    newLine = QtCore.Signal()

    def keyReleaseEvent(self,event):
        if event.key() == QtCore.Qt.Key_Return:
            self.newLine.emit()
            event.accept()
        else:
            event.ignore()

class keyEnterReceiver(QtCore.QObject):
    def eventFilter(self, obj, event):
        if (event.type()==QtCore.QEvent.KeyPress):
            key = event
            if (key.key()==QtCore.Qt.Key_Enter) or (key.key()==QtCore.Qt.Key_Return):
                if hasattr(obj,"enterSignal"):
                    obj.enterSignal.emit()
            else:
                return QtCore.QObject.eventFilter(self, obj, event);
            return True
        else:
            return QtCore.QObject.eventFilter(self, obj, event);

        return False

class CChooseTaskFrame(QtWidgets.QFrame):

    taskClicked = QtCore.Signal(str)
    closed = QtCore.Signal()

    def closeEvent(self):
        self.closed.emit()

    def __init__(self, parent):
        QtWidgets.QFrame.__init__(self, parent)
        self.setLayout(QtWidgets.QVBoxLayout())

        self.taskTree = CTaskTree(self)
        self.taskTree.setHeaderHidden(True)

        modeLayout = QtWidgets.QHBoxLayout()
        modeLayout.addWidget(QtWidgets.QLabel("Only show tasks related to:"))

        self.refinementButton = QtWidgets.QRadioButton("Refinement")
        #self.coreButton = QtWidgets.QRadioButton("Core")
        self.anyButton = QtWidgets.QRadioButton("Any")
        self.mrButton = QtWidgets.QRadioButton("Molecular\nreplacement")
        self.epButton = QtWidgets.QRadioButton("Experimental\nphasing")

        modeLayout.addWidget(self.mrButton)
        modeLayout.addWidget(self.epButton)
        modeLayout.addWidget(self.refinementButton)
        #modeLayout.addWidget(self.coreButton)
        modeLayout.addWidget(self.anyButton)

        self.anyButton.setChecked(True)

        from ..core.CCP4Preferences import PREFERENCES
        if PREFERENCES().SHOW_TASK_MODE_BUTTONS:
            self.layout().addLayout(modeLayout)

        self.taskSearchBox = QtWidgets.QLineEdit()
        self.taskSearchBox.setPlaceholderText("Only show tasks containing text typed here")
        searchLayout = QtWidgets.QHBoxLayout()
        searchLayout.addWidget(QtWidgets.QLabel("Filter:"))
        searchLayout.addWidget(self.taskSearchBox)
        self.layout().addLayout(searchLayout)
        self.layout().addWidget(self.taskTree)
        self.taskTreeLoaded = False
        self.loadTaskTree()

        buttonGroup = QtWidgets.QButtonGroup(self)
        buttonGroup.addButton(QtWidgets.QPushButton('New job', self), 0)
        self.cancelButton = QtWidgets.QPushButton('Cancel', self)
        buttonGroup.addButton(self.cancelButton, 1)
        layout0 = QtWidgets.QHBoxLayout()
        layout0.addStretch(1)
        for but in buttonGroup.buttons():
            layout0.addWidget(but)
            layout0.addStretch(1)
        self.layout().addLayout(layout0)
        buttonGroup.button(0).clicked.connect(self.handleChooseTask0)
        self.taskTree.doubleClicked.connect(self.handleChooseTask)
        buttonGroup.button(1).clicked.connect(functools.partial(self.window().showTaskChooser, False))
        from ..core.CCP4WorkflowManager import WORKFLOWMANAGER
        WORKFLOWMANAGER().listChanged.connect(self.invalidateTaskTree)
        from ..core.CCP4CustomTaskManager import CUSTOMTASKMANAGER
        CUSTOMTASKMANAGER().listChanged.connect(self.invalidateTaskTree)
        PREFERENCES().preferencesSaved.connect(self.invalidateTaskTree)

        @QtCore.Slot('QModelIndex','QModelIndex')
        def toggleNewButton(current,previous):
            taskName = str(self.taskTree.model().data(current,QtCore.Qt.UserRole + 2))
            if len(taskName) == 0:
                buttonGroup.button(0).setEnabled(False)
            else:
                buttonGroup.button(0).setEnabled(True)
        buttonGroup.button(0).setEnabled(False)
        self.taskTree.taskChanged.connect(toggleNewButton)
        enterReceiver = keyEnterReceiver(self.taskTree)
        self.taskTree.installEventFilter(enterReceiver)
        self.taskTree.enterSignal.connect(self.handleChooseTask0)

    def showEvent(self, theEvent):
        if not self.taskTreeLoaded:
            self.loadTaskTree()
        super(CChooseTaskFrame, self).showEvent(theEvent)

    @QtCore.Slot(str)
    def handleChooseTask(self, idx):
        try:
            taskName = self.taskTree.model().data(idx,QtCore.Qt.UserRole + 2)
            if taskName is None:
                return
            taskName = str(self.taskTree.model().data(idx,QtCore.Qt.UserRole + 2))
        except:
            return
        if len(taskName) == 0:
            return
        self.taskClicked.emit(taskName)

    @QtCore.Slot()
    def handleChooseTask0(self):
        idxs = self.taskTree.selectedIndexes()
        if len(idxs)>0:
            self.handleChooseTask(idxs[0])
    
    @QtCore.Slot()
    def invalidateTaskTree(self):
        self.taskTreeLoaded = False
        if CCP4Utils.isAlive(self) and self.isVisible():
            self.loadTaskTree()
    
    def loadTaskTree(self):
        # KJS: 4k (w,h) approx. double that of 1080p.
        from ..core.CCP4Preferences import PREFERENCES
        compact = PREFERENCES().COMPACT_TASK_MENU
        if PREFERENCES().HD_ICONS:
            setHD = "HD"
            icsize = (72, 96)
        else:
            setHD = None
            icsize = (24, 48)
        taskListModel = TreeModel()
        self.taskTree.setModel(taskListModel)
        rootItem = taskListModel.rootItem

        no_win = ["arcimboldo", "morda_i2", "arp_warp_classic", "xia2_xds"]
        for module, title, taskList in TASKMANAGER().taskTree():

            if module == "developer_tools" and not KJS_DEVINTER_ON:
                continue

            moduleItem  = TreeItem("<b>"+title+"</b>",rootItem)
            modIcon = getMenuIcon(None, module, setHD) 
            moduleItem.setIcon(modIcon)
            rootItem.appendChild(moduleItem)
            moduleItem.setIconSize(QtCore.QSize(icsize[0], icsize[0]))

            for taskName in taskList:
                try:
                    if sys.platform == "win32" and (taskName in no_win):
                        continue
                    taskIcon =  getMenuIcon(None, taskName, setHD)
                    desc = TASKMANAGER().getTaskAttribute(taskName, 'DESCRIPTION')
                    rank = TASKMANAGER().getTaskAttribute(taskName, 'RANK')
                    if compact:
                        name = TASKMANAGER().getTitle(taskName)
                    else:
                        name = "<b>"+TASKMANAGER().getTitle(taskName)+"</b>"
                    if desc is not None and not compact:
                        name += "<br/>" + "<i>"+desc+"</i>"
                    taskItem  = TreeItem(name,moduleItem)
                    if desc is not None:
                        taskItem.setLongDescription(desc)
                    taskItem.setIcon(taskIcon)
                    taskItem._itemType = 1
                    taskItem.setName(taskName)
                    if rank == 1:
                        taskItem.setTaskType('taskpipe')
                    else:
                        taskItem.setTaskType('tasktool')
                    helpFile = TASKMANAGER().searchHelpFile(name=taskName)
                    if helpFile is None:
                        taskItem.setHelpAvailable(False)
                    else:
                        taskItem.setHelpAvailable(True)
                    if compact:
                        taskItem.setIconSize(QtCore.QSize(icsize[0], icsize[0]))
                    else:
                        taskItem.setIconSize(QtCore.QSize(icsize[1], icsize[1]))
                    moduleItem.appendChild(taskItem)
                except:
                    print(taskName,"failed")
            
        self.taskTreeLoaded = True

        delegate = HTMLDelegate()
        if compact:
            delegate.setCompact(True)
        self.taskTree.setItemDelegate(delegate)
        self.taskTree.setMouseTracking(True) 
        self.taskTree.wheel_moved.connect(delegate.clearShowHelp)

        oldModel = self.taskTree.model()
        taskListProxyModel = CTaskListProxyModel()
        taskListProxyModel.setSourceModel(oldModel)
        self.taskTree.setModel(taskListProxyModel)

        self.taskSearchBox.textChanged.connect(taskListProxyModel.setFilterString)
        self.refinementButton.clicked.connect(functools.partial(taskListProxyModel.setFilterMode,("refinement",)))
        #self.coreButton.clicked.connect(functools.partial(taskListProxyModel.setFilterMode,("core",)))
        self.anyButton.clicked.connect(functools.partial(taskListProxyModel.setFilterMode,()))
        self.mrButton.clicked.connect(functools.partial(taskListProxyModel.setFilterMode,("molecular_replacement",)))
        self.epButton.clicked.connect(functools.partial(taskListProxyModel.setFilterMode,("expt_phasing",)))
        self.taskTree.sizeChangedSignal.connect(delegate.setRequiredSize)
        self.taskTree.sizeChangedSignal.connect(self.taskTree.doItemsLayout)

        @QtCore.Slot('QModelIndex')
        def handleHelpClicked(idx):
            try:
                taskName = str(self.taskTree.model().data(idx,QtCore.Qt.UserRole + 2))
                fileName = TASKMANAGER().searchHelpFile(name=taskName)
                from .CCP4WebBrowser import WEBBROWSER
                WEBBROWSER().openFile(fileName)
            except:
                pass
        delegate.helpClicked.connect(handleHelpClicked)

        @QtCore.Slot()
        def expandIfFiltered():
            self.taskTree.collapseAll()
            if len(str(self.taskSearchBox.text()))>0:
                self.taskTree.expandAll()
            elif not self.anyButton.isChecked():
                self.taskTree.expandAll()

        self.taskSearchBox.textChanged.connect(expandIfFiltered)
        self.refinementButton.clicked.connect(expandIfFiltered)
        #self.coreButton.clicked.connect(expandIfFiltered)
        self.anyButton.clicked.connect(expandIfFiltered)
        self.mrButton.clicked.connect(expandIfFiltered)
        self.epButton.clicked.connect(expandIfFiltered)


class CTaskTree(QtWidgets.QTreeView):

    taskChanged = QtCore.Signal('QModelIndex','QModelIndex')
    taskDropped = QtCore.Signal(str,str)
    sizeChangedSignal = QtCore.Signal('QSize')
    enterSignal = QtCore.Signal()
    wheel_moved = QtCore.Signal()

    def __init__(self, parent):
        QtWidgets.QTreeView.__init__(self, parent)
        self.setObjectName('taskTree')
        self.setDragEnabled(True)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DropOnly)
        self.setAcceptDrops(True)

    def currentChanged(self,current,previous):
        self.taskChanged.emit(current,previous)
        QtWidgets.QTreeView.currentChanged(self,current,previous)
        if (self.viewport().size().height()- self.visualRect(current).y()) < 30  :
            #It is extremely bad that Qt ever lets this happen.
            #Keyboard scrolling with items of different height trigger this problem which took me *ages* to fix.
            self.scrollTo(current,QtWidgets.QAbstractItemView.PositionAtBottom)

    def mimeTypes(self):
        typesList = []
        from ..qtcore.CCP4CustomMimeTypes import MIMETYPESHANDLER
        for item in list(MIMETYPESHANDLER().mimeTypes.keys()):
            typesList.append(item)
        typesList.append('jobid')
        return typesList

    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat('FollowFromJob'):
            event.accept()
        else:
            event.ignore()
    
    def dragMoveEvent(self,event):
        dropItem = self.itemAt(event.pos().x(), event.pos().y())
        if event.mimeData().hasFormat('jobid'):
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self,event):
        targetItem = self.itemAt(event.pos().x(), event.pos().y())
        if event.mimeData().hasFormat('FollowFromJob'):
            taskName = str(targetItem.data(0, 101))
            data = str(event.mimeData().data('FollowFromJob').data())
            self.taskDropped.emit(taskName,data)
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def wheelEvent(self,event):
        self.wheel_moved.emit()
        return QtWidgets.QTreeView.wheelEvent(self,event)

    def resizeEvent(self,event):
        self.sizeChangedSignal.emit(self.size())
        QtWidgets.QTreeView.resizeEvent(self,event)


def FILEWATCHER():
    if CFileSystemWatcher.insts is None:
        CFileSystemWatcher.insts = CFileSystemWatcher()
    return CFileSystemWatcher.insts


class CFileSystemWatcher(QtCore.QFileSystemWatcher):

    jobFileChanged = QtCore.Signal(tuple)

    # Sub-class so we keep track of jobId as well as path
    insts = None

    def __init__(self, parent=None):
        if parent is None:
            parent = QTAPPLICATION()
        QtCore.QFileSystemWatcher.__init__(self, parent)
        # Fix for fail on NFS (and perhaps elsewhere) to convert to using internal poller mechanism
        # https://bugreports.qt.io/browse/QTBUG-8351
        from ..core.CCP4Preferences import PREFERENCES
        if PREFERENCES().FILESYSTEMWATCHERPOLLER:
            self.setObjectName("_qt_autotest_force_engine_poller")
        self.jobPaths = {}
        self.jobPendingPaths = {}
        self.fileChanged.connect(self.handleFileChanged)
        self.directoryChanged.connect(self.handleDirectoryChanged)
        self.fileSizes = {}
        self._diagnostic = False
        self.jobsByUpdateInterval = {}

    @QtCore.Slot(str)
    def handleFileChanged(self, path):
        # Small changes introduced to handle a not-uncommon case.  Specifically, if a watched file is overwritten by doing
        # a pseudo-atomic overwrite (write to temporary file, and then rename over the top of existing file), then a
        # sad things happens: the rename causes the original file to be unlinked andhence removed from the
        # observed path list.
        if self._diagnostic:
            print('CFileSystemWatcher.handleFileChanged', path)
            sys.stdout.flush()
        path = str(path)
        for jobId,pathList in list(self.jobPaths.items()):
            if path in pathList:
                # Further processing appropriate only if the path represents a valid (existing) file object
                if os.path.isfile(path):
                    # Ignore change in files that makes them > 100 characters (arbitrary cutoff)smaller (since probably implies they
                    # are in the course of being overwritten)
                    newFileSize = os.stat(path).st_size
                    if path not in self.fileSizes or (path in self.fileSizes and (self.fileSizes[path] < (newFileSize+100))):
                        self.jobFileChanged.emit((jobId,path))
                    self.fileSizes[path] = newFileSize
        if self._diagnostic:
            print('CFileSystemWatcher Currently watching files:',[str(file) for file in self.files()])
            print('CFileSystemWatcher Currently watching directories:',[str(directory) for directory in self.directories()])
        return

    def clear(self):
        self.removePaths(self.files())
        self.jobPaths = {} 
        self.jobPendingPaths = {}
    
    def addJobPath(self,jobId,path,updateInterval=None):
        # Beware the file we want to watch might not exist
        # so put a watch on the (hopefully existing) parent directory
        if self._diagnostic:
            print('FILEWATCHER.addJobPath',jobId,path,updateInterval)
            sys.stdout.flush()
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print('FILEWATCHER.addJobPath',jobId,path,updateInterval)
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        if updateInterval is not None:
            if jobId not in self.jobsByUpdateInterval:
                self.jobsByUpdateInterval[jobId] = { 'path' : path , 'size' : 0 }
            return
        if jobId not in self.jobPaths:
            self.jobPaths[jobId] = []
        if not path in self.jobPaths[jobId]:
            self.jobPaths[jobId].append(path)
        directoryOfPath, fileName = os.path.split(path)
        currentlyWatchedDirectories = [str(directory) for directory in self.directories()]
        if not directoryOfPath in currentlyWatchedDirectories and os.path.exists(directoryOfPath):
            self.addPath(directoryOfPath)
        if not directoryOfPath in currentlyWatchedDirectories and os.path.exists(os.path.dirname(directoryOfPath)):
            self.addPath(os.path.dirname(directoryOfPath))
            if jobId not in self.jobPendingPaths:
                self.jobPendingPaths[jobId] = []
            self.jobPendingPaths[jobId].append(directoryOfPath)
        if os.path.isfile(path):
            self.addPath(path)
        if self._diagnostic:
            print('CFileSystemWatcher Currently watching files:',[str(file) for file in self.files()])
            print('CFileSystemWatcher Currently watching directories:',[str(directory) for directory in self.directories()])

    @QtCore.Slot(str)
    def handleDirectoryChanged(self, path):
        if self._diagnostic:
            print('CFileSystemWatcher.handleDirectoryChanged', path)
            sys.stdout.flush()
        #Here in case of a directory event: go through the
        #files we wish to be watching and add them to watch list *iff* they exist and are
        #not currently being watched
        for jobId, pathsOfJob in list(self.jobPaths.items()):
            if jobId in self.jobPendingPaths:
                added = []
                for pendingDir in self.jobPendingPaths[jobId]:
                    if not pendingDir in self.directories() and os.path.exists(pendingDir):
                        #print("----------------------------------------")
                        #print("addPath 4 (directory)",pendingDir)
                        self.addPath(pendingDir)
                        #print("----------------------------------------")
                        added.append(pendingDir)
                    self.jobPendingPaths[jobId] = [x for x in self.jobPendingPaths[jobId] if not x in added]

        currentlyWatchedPaths = [str(file) for file in self.files()]
        for jobId, pathsOfJob in list(self.jobPaths.items()):
            for pathOfJob in pathsOfJob:
                directoryOfPath, fileName = os.path.split(pathOfJob)
                if directoryOfPath == str(path) and pathOfJob not in currentlyWatchedPaths and os.path.isfile(pathOfJob):
                    self.addPath(pathOfJob)
                    self.handleFileChanged(pathOfJob)
        if self._diagnostic:
            print('CFileSystemWatcher Currently watching files:',[str(file) for file in self.files()])
            print('CFileSystemWatcher Currently watching directories:',[str(directory) for directory in self.directories()])

    def removeJob(self,jobId):
        if  jobId in self.jobPaths:
            self.removePaths(self.jobPaths[jobId])
            del self.jobPaths[jobId]
        if self._diagnostic:
            print('FILEWATCHER.removeJob',jobId,jobId in self.jobsByUpdateInterval)
        if jobId in self.jobsByUpdateInterval:
            del self.jobsByUpdateInterval[jobId]

    def removePath(self,path):
        for jobId,pathList in list(self.jobPaths.items()):
            if path in pathList:
                self.jobPaths[jobId].remove(path)
                QtCore.QFileSystemWatcher.removePath(self,path)
                return

    def triggerJobsByUpdateInterval(self):
        for jobId,value in list(self.jobsByUpdateInterval.items()):
            if os.path.isfile(value['path']):
                newSize = os.stat(value['path']).st_size
                if self._diagnostic:
                    print('FILEWATCHER',jobId,value,newSize)
                if newSize > value['size']:
                    value['size'] = newSize
                    self.jobFileChanged.emit((jobId,value['path']))


class CReportView(QtWidgets.QStackedWidget):

    nextTask = QtCore.Signal(str,str)
    interrupt = QtCore.Signal(str)
    reportAvailable = QtCore.Signal(str,bool)
    labelEdited =  QtCore.Signal()

    def handleFind(self, direction=1):
        widget=self.webView
        if widget is None or not widget.isSearchable():
            return
        text = self.findWidget.text()
        rv = widget.findText(subString=text, direction=direction, caseSensitive=self.findWidget.caseSensitive())

    def handleFindNext(self):
        self.handleFind(1)

    def handleFindPrevious(self):
        self.handleFind(-1)

    def __init__(self,parent=None):
        QtWidgets.QStackedWidget.__init__(self,parent=parent)
        #self.setLayout(QtWidgets.QVBoxLayout())
        #self.layout().setSpacing(0)
        #self.layout().setContentsMargins(1,1,1,1)
        self.setObjectName('reportView')
        self.openJob = COpenJob()
        self.openChildJob = None
        self.generator = None
        self._reportFile = None
        self.findShortcut = QtWidgets.QShortcut(QtGui.QKeySequence.Find, self)
        self.findWidget = CCP4WebBrowser.CFindFrame(self)

        @QtCore.Slot()
        def showFindWidget():
           self.findWidget.show()
           self.findWidget.findTextWidget.setFocus(QtCore.Qt.OtherFocusReason)

        self.findShortcut.activated.connect(showFindWidget)
        self.findWidget.findNext.connect(self.handleFindNext)
        self.findWidget.findPrevious.connect(self.handleFindPrevious)

        self.webFrame = QtWidgets.QFrame(self)
        self.webFrame.setLayout(QtWidgets.QVBoxLayout())
        #self.webFrame.layout().setSpacing(0)
        #self.webFrame.layout().setContentsMargins(1,1,1,1)
        self.webFrame.setObjectName('reportFrame')
        self.linksStack = QtWidgets.QStackedWidget(self)
        self.linksFrame = []
        self._channelIsSet = False
        labelList = [ [], [['Results','results'],['Input data','inputData'],['Output data','outputData'],['Show details','Show details'],['Job details','jobDetails']],
                      [['Error output','errors'],['Terminal output','terminal'],['Diagnostic','diagnostic']]
                    ]
        for ii in range(3):
            self.linksFrame.append(QtWidgets.QFrame(self))
            self.linksFrame[ii].setObjectName('reportLinkFrame'+str(ii))
            self.linksFrame[ii].setFrameShape(QtWidgets.QFrame.StyledPanel)
            self.linksFrame[ii].setMaximumHeight(25)
            self.linksFrame[ii].setLayout(QtWidgets.QHBoxLayout())
            self.linksFrame[ii].layout().setSpacing(0)
            self.linksFrame[ii].layout().setContentsMargins(1,1,1,1)
            for label,link in labelList[ii]:
                but = QtWidgets.QPushButton(label,self)
                but.setFlat(True)
                self.linksFrame[ii].layout().addWidget(but)
                but.clicked.connect(functools.partial(self.handleLinkButton,link))
            self.linksStack.addWidget(self.linksFrame[ii])
        self.webFrame.layout().addWidget(self.linksStack)
        self.backToReportButton = QtWidgets.QPushButton("Back to main report")
        from . import CCP4WebView
        self.webView = CCP4WebView.CWebView(self,blockLoading=True)
        self.webView.loadFinished.connect(self.handleLoadFinished)
        self.webView.page().profile().downloadRequested.connect(self.handleDownload)
        self.webFrame.layout().addWidget(self.backToReportButton)
        self.webFrame.layout().addWidget(self.webView)
        self.webFrame.layout().addWidget(self.findWidget)
        self.findWidget.hide()
        self.webFrame.layout().setStretch(2,10.0)
        self.addWidget(self.webFrame)
        self.backToReportButton.hide()
        @QtCore.Slot()
        def hideBackButton():
            self.backToReportButton.hide()
            self.webFrame.layout().setSpacing(0)
            self.webFrame.layout().setContentsMargins(1,1,1,1)
        @QtCore.Slot()
        def showBackButton():
            self.backToReportButton.show()
            self.webFrame.layout().setSpacing(2)
            self.webFrame.layout().setContentsMargins(2,2,2,2)

        self.webView.searchFound.connect(self.findWidget.searchCallback)

        self.webView.nonReportNavigated.connect(showBackButton)
        self.webView.reportNavigated.connect(hideBackButton)

        self.webView.loadFinished.connect(self.handleFindNext)

        self.backToReportButton.clicked.connect(self.webView.back)
        self.backToReportButton.clicked.connect(hideBackButton)
        self.nextFrame = QtWidgets.QFrame(self)
        self.nextFrame.setObjectName('highlight')
        self.nextFrame.setLayout(QtWidgets.QGridLayout())
        self.nextFrame.layout().setSpacing(0)
        self.nextFrame.layout().setContentsMargins(0,0,0,0)
        self.webFrame.layout().addWidget(self.nextFrame)
        self.textView = None
        self.webView.page().NavigationRequest.connect(self.handleNavigationRequest)
        self.webView.page().CustomMimeTypeRequested.connect(self.handleCustomMimeTypeRequested)
        FILEWATCHER().jobFileChanged.connect(self.handleFileChanged)
        from ..qtcore.CCP4JobController import JOBCONTROLLER
        JOBCONTROLLER().remoteJobUpdated.connect(self.handleFileChanged)
        self.openJob.jobStatusUpdated.connect(self.handleJobStatusUpdated)
        self.openJob.workflowJobStatusUpdated.connect(self.handleWorkflowJobStatusUpdated)
        def actOnJSSignal(args):
            print(args)
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            from ..qtcore.CCP4Launcher import LAUNCHER
            if 'exe' in args and args['exe'] == 'loggraph':
                reportFile = PROJECTSMANAGER().makeFileName(jobId=self.openJob.jobId,mode='REPORT')
                launchArgs = []
                if 'ccp4_data_id' in args:
                    launchArgs.extend(["--select",str(args['ccp4_data_id'])])
                if 'ccp4_data_current_index' in args:
                    launchArgs.extend(["-G",str(args['ccp4_data_current_index'])])
                launchArgs.append(reportFile)
                LAUNCHER().launch(viewer='loggraph',argList=launchArgs)
            elif 'exe' in args and args['exe'] == 'CCP4mg':
                if args.get('sceneFile',None) is not None:
                    sceneFile = os.path.join(PROJECTSMANAGER().jobDirectory(args['jobId']),os.path.split(args['sceneFile'])[-1])
                    LAUNCHER().launch(viewer='CCP4mg',argList=['-nore',sceneFile],projectId=args.get('projectId',None),guiParent=self.parent())
                elif args.get('pictureDefinitionFile',None) is not None:
                    LAUNCHER().launch(viewer='CCP4mg',argList=['-nore','-pict',args['pictureDefinitionFile']],projectId=args.get('projectId',None),guiParent=self.parent())
                elif args.get('jobId',None) is not None:
                    LAUNCHER().openInViewer(viewer='CCP4mg',jobId=args.get('jobId'),projectId=args.get('projectId',None),guiParent=self.parent())
            elif args.get('exe') == 'Coot':
                LAUNCHER().openInViewer(viewer='coot_job',jobId=args.get('jobId',None),projectId=args.get('projectId',None),guiParent=self.parent())
            elif 'action' in args and args['action'] == 'download':
                self.webView.emitDownloadRequest(args["jobId"],args["dataName"])
            elif 'action' in args and args['action'].startswith("view_"):
                fileInfo = PROJECTSMANAGER().db().getFileInfo(args.get('dbFileId',None),mode=['projectid','relpath','filename','annotation','filesubtype','filecontent','filetype'])
                projectDir = PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectdirectory')
                filePath = os.path.join(projectDir,fileInfo['relpath'],fileInfo['filename'])
                if args['action'] == "view_text":
                    if args["filetypeclass"] == "AsuDataFile":
                        tree = CCP4Utils.openFileToEtree(filePath)
                        seqs = tree.xpath("//CAsuContentSeq")
                        html = """<body>"""
                        for seq in seqs:
                            desc = seq.xpath("description")[0].text
                            nCopies = seq.xpath("nCopies")[0].text
                            polymerType = seq.xpath("polymerType")[0].text
                            sequence = seq.xpath("sequence")[0].text
                            if nCopies == 1:
                                html += "<h4>"+nCopies+" copy of '"+desc+"' ("+polymerType+")</h4>\n"
                            else:
                                html += "<h4>"+nCopies+" copies of '"+desc+"' ("+polymerType+")</h4>\n"
                            n = 80
                            seqs = [sequence[i:i+n] for i in range(0, len(sequence), n)]
                            html += '<pre>'
                            for seq in seqs:
                                html += seq +"\n"
                            html += "</pre>\n"
                        win = QtGui.QDialog()
                        widget = QtGui.QTextEdit()
                        widget.setReadOnly(True)
                        widget.setHtml(html)
                        layout = QtGui.QVBoxLayout()
                        layout.addWidget(widget)
                        dbb = QtGui.QDialogButtonBox()
                        layout.addWidget(dbb)
                        closeButton = dbb.addButton(QtGui.QDialogButtonBox.Close)
                        closeButton.clicked.connect(win.accept)
                        win.setLayout(layout)
                        win.exec_()
                    else:
                        from .CCP4WebBrowser import WEBBROWSER
                        WEBBROWSER().openFile(fileName=filePath,format='text/plain')
                else:
                    LAUNCHER().openInViewer(viewer=args['action'][5:],fileName=filePath,projectId=fileInfo['projectid'],guiParent=self.parent())
            elif 'action' in args and args['action'] == "quick_view":
                fileInfo = PROJECTSMANAGER().db().getFileInfo(args.get('dbFileId',None),mode=['projectid','relpath','filename','annotation','filesubtype','filecontent','filetype'])
                projectDir = PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectdirectory')
                filePath = os.path.join(projectDir,fileInfo['relpath'],fileInfo['filename'])
                gemmi_struct = gemmi.read_structure(filePath)
                html = """<body>"""
                html += "<p>Spacegroup: " + gemmi_struct.spacegroup_hm + "</p>"
                cell = gemmi_struct.cell
                html += "<p>Cell: " + str(cell.a) + " " + str(cell.b) + " " + str(cell.c) + " " + str(cell.alpha) + " " + str(cell.beta) + " " + str(cell.gamma) + "</p>"
                html += "<p></p>"
                html += "<p></p>"
                html += "<p>Structure contains " + str(len(gemmi_struct)) + " models</p>"
                html += "<p>containing chains:</p>"
                for gemmi_model in gemmi_struct:
                    for chain in gemmi_model:
                         html += "<p>"+chain.name + " (" + str(len(chain)) + " residues)</p>"
                html += """</body>"""
                win = QtWidgets.QDialog()
                widget = QtWidgets.QTextEdit()
                widget.setReadOnly(True)
                widget.setHtml(html)
                layout = QtWidgets.QVBoxLayout()
                layout.addWidget(widget)
                dbb = QtWidgets.QDialogButtonBox()
                layout.addWidget(dbb)
                closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Close)
                closeButton.clicked.connect(win.accept)
                win.setLayout(layout)
                win.exec_()

            elif 'action' in args and args['action'] == "export_as_map":
                def coefficientsToMap(coefficientsPath, mapPath=None, overSample=1.0):
                    mtz_file = clipper.CCP4MTZfile()
                    hkl_info = clipper.HKL_info()
                    mtz_file.open_read (str(coefficientsPath))
                    mtz_file.import_hkl_info ( hkl_info )
                    sg, cell = hkl_info.spacegroup(), hkl_info.cell()
                    fphidata = clipper.HKL_data_F_phi_float(hkl_info)
                    mtz_file.import_hkl_data( fphidata, str("/*/*/[F,PHI]") );
                    mtz_file.close_read()
                    #Clipper will sample the output map according to Fourier theory and the nominal resolution
                    #for visualisation, it is generally nicer to make things a bit more finely sampled
                    fudgedResolution = hkl_info.resolution()
                    fudgedResolution.init(hkl_info.resolution().limit()/overSample)
                    mygrid=clipper.Grid_sampling ( hkl_info.spacegroup(), hkl_info.cell(), fudgedResolution )
                    mymap = clipper.Xmap_float(hkl_info.spacegroup(), hkl_info.cell(), mygrid )
                    mymap.fft_from(fphidata)

                    mapout = clipper.CCP4MAPfile()
                    if mapPath is None:
                        coefficientsRoot, extension = os.path.splitext(os.path.abspath(coefficientsPath))
                        mapPath = coefficientsRoot+".map"

                    mapout.open_write( mapPath )
                    mapout.export_xmap_float( mymap )
                    mapout.close_write()
                    return mapPath
            
                print("Export as map!")
                fileInfo = PROJECTSMANAGER().db().getFileInfo(args.get('dbFileId',None),mode=['filetype','filename','projectid','relpath'])
                projectDir = PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectdirectory')
                fileType = fileInfo["filetype"]
                fileName = fileInfo["filename"]
                filePath = os.path.join(projectDir,fileInfo['relpath'],fileInfo['filename'])
                mapPath = coefficientsToMap(filePath, overSample=1.0)
                print("mapPath",mapPath)
                fileBrowser = QtWidgets.QFileDialog()
                fileBrowser.setWindowTitle("Save map")
                fileBrowser.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
                fileBrowser.setFileMode(QtWidgets.QFileDialog.AnyFile)
                fileBrowser.setDefaultSuffix(".map")
                fileBrowser.setNameFilters(["Map files (*.map *.ccp4)","Any files (*)"])
                fileBrowser.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
                if fileBrowser.exec_():
                   if len(fileBrowser.selectedFiles()) > 0:
                       path = fileBrowser.selectedFiles()[0]
                       if path:
                           shutil.copyfile(mapPath,path)

            elif 'action' in args and args['action'] == "export":
                fileInfo = PROJECTSMANAGER().db().getFileInfo(args.get('dbFileId',None),mode=['filetype','filename','projectid','relpath'])
                projectDir = PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectdirectory')
                fileType = fileInfo["filetype"]
                fileName = fileInfo["filename"]
                filePath = os.path.join(projectDir,fileInfo['relpath'],fileInfo['filename'])
                from ..qtcore.CCP4CustomMimeTypes import MIMETYPESHANDLER
                filters = MIMETYPESHANDLER().getMimeTypeInfo(fileType,'filter')
                suffices =  MIMETYPESHANDLER().getMimeTypeInfo(fileType,'fileExtensions')
                if len(suffices)>0:
                    defaultSuffix =  suffices[0]
                else:
                    defaultSuffix =  ""
                from . import CCP4FileBrowser
                fileBrowser = CCP4FileBrowser.CFileDialog(parent=self,
                                      title='Export '+fileName,
                                      filters = [filters],
                                      defaultSuffix = defaultSuffix,
                                      fileMode = QtWidgets.QFileDialog.AnyFile)
                fileBrowser.selectFile.connect(functools.partial(self.exportData,filePath,args.get('dbFileId',None)))
                fileBrowser.show()

            elif 'action' in args and args['action'] == "editLabel":
                from . import CCP4Widgets
                d = CCP4Widgets.CEditFileLabel(parent=self,fileId=args.get('dbFileId',None))
                d.accepted.connect(self.labelEdited.emit)

            elif 'action' in args and args['action'] == "WebGL":
                from ..report.CCP4ReportParser import WEBGLSOURCES, MTZToB64Map
                htmlBase = args['htmlBase']
                print("Now WebGL it")
                onLoadText = """<script>
     var queue = [];
     var loaded = false;

     function enqueue(callback,fogdiv,clipdiv)
     {
         if(!loaded) queue.push([callback,fogdiv,clipdiv]);
         else callback(fogdiv,clipdiv);
     }

     var getElementsByRegex = function(node,pattern,recurse){
        var arrElements = [];   // to accumulate matching elements
        var re = new RegExp(pattern);   // the regex to match with

        function findRecursively(aNode) { // recursive function to traverse DOM
            if (!aNode) 
                return;
            if (aNode.tagName !== undefined && aNode.tagName.search(re) != -1)
                arrElements.push(aNode);  // FOUND ONE!
            for (var idx in aNode.childNodes) // search children...
                findRecursively(aNode.childNodes[idx]);
        };

        findRecursively(node); // initiate recursive matching
        return arrElements; // return matching elements
    };
    function loadData(gl,pdbatoms,enerLib,bruteForceHB,wizard){
        var ribbons = [];
        //FIXME - This is hackery ... need more wizards to be written.
        if(wizard.startsWith("Site and ribbons by ")){
           ribbons = wizardSiteAndRibbonsByChain(pdbatoms,enerLib,bruteForceHB);
        } else if(wizard.startsWith("colour chains")) {
           ribbons = wizardRibbonsByChain(pdbatoms,enerLib,bruteForceHB);
        } else {
           ribbons = wizardBonds(pdbatoms,enerLib,bruteForceHB);
        }
        //ribbons = wizardWorms(pdbatoms,enerLib,bruteForceHB);

        var appendStart = new Date().getTime();
        for(var i=0;i<ribbons.length;i++){
            // FIXME - Should be done in wizardBonds call????
            ribbons[i]["symmetry"] = pdbatoms["symmetry"];
            gl.appendOtherData(ribbons[i],true);
        }
        console.log("Time to do appendData(actual): "+(new Date().getTime()-appendStart));
        gl.buildBuffers();
        console.log("Time to done buildBuffers after append : "+(new Date().getTime()-appendStart));
        gl.drawScene();

        var hier = pdbatoms["atoms"];
        gl.setOrigin(hier[0].centre());
        console.log(hier[0].centre());

    }

    function _base64ToArrayBuffer(base64) {
        var binary_string =  window.atob(base64);
        var len = binary_string.length;
        var bytes = new Uint8Array( len );
        for (var i = 0; i < len; i++)        {
            bytes[i] = binary_string.charCodeAt(i);
        }
        return bytes.buffer;
    }

    function doParseAndLoad(gl,enerLib,xmlDoc){
        var extent = 0.0;
            var atoms_ids = {};
            if(xmlDoc.getElementsByTagName("ccp4i2_body").length>0){
                var ccp4i2_body = xmlDoc.getElementsByTagName("ccp4i2_body")[0];
                if(ccp4i2_body.getElementsByTagName("data").length>0){
                    var data = ccp4i2_body.getElementsByTagName("data")[0];
                    var molData = data.getElementsByTagName("MolData");
                    for (var i=0; i < molData.length; i++) {
                        var customResCIFFiles = molData[i].getElementsByTagName("customResCIFFiles");
                        for(var k=0;k<customResCIFFiles.length;k++){
                            var cifmonomers = customResCIFFiles[k].getElementsByTagName("cifmonomer");
                            for(var l=0;l<cifmonomers.length;l++){
                                console.log("Have a cifmonomer");
                                var name = cifmonomers[l].getElementsByTagName("name");
                                var ciffiledatanodes = cifmonomers[l].getElementsByTagName("filedata");
                                if(name.length>0&&ciffiledatanodes.length>0){
                                    console.log("Have some text");
                                    var nameText = getNodeText(name[0]);
                                    var cifText = getNodeText(ciffiledatanodes[0]);
                                    console.log(nameText);
                                    console.log(cifText);
                                    enerLib.addCIFAtomTypes(nameText,cifText);
                                    enerLib.addCIFBondTypes(nameText,cifText);
                                }
                            }
                        }
                        var modid = molData[i].getAttribute("id");
                        var fileDatas = molData[i].getElementsByTagName("filedata");
                        var fileData;
                        for(var ifd=0;ifd<fileDatas.length;ifd++){
                            console.log(ifd);
                            if(fileDatas[ifd].parentNode===molData[i]){
                                fileData = fileDatas[ifd];
                                break;
                            }
                        }
                        var pdbatoms = parsePDB(fileData.innerHTML.split(String.fromCharCode(10)));
                        var hier = pdbatoms["atoms"];
                        atoms_ids[modid] = pdbatoms;
                        var molDisp = molData[i].getElementsByTagName("MolDisp");
                        for (var j=0; j < molDisp.length; j++) {
                            var selection = molDisp[j].getElementsByTagName("select")[0].innerHTML;
                            var style = molDisp[j].getElementsByTagName("style")[0].innerHTML;
                            if(true){
                                var model = hier[0]; //FIXME - Loop over models?
                                var selAtoms = model.getAtoms(selection);
                                //FIXME - Now what do I do?
                                if(style==="BALLSTICK"||style==="BONDS"||style==="THINBONDS"||style==="FATBONDS"||style==="SPHERES"||style==="CYLINDERS"){
                                    //FIXME - pick up the colour scheme ....
                                    var colourScheme = new ColourScheme(pdbatoms);
                                    var atomColours = colourScheme.colourByAtomType();
                                    var contactsAndSingletons = model.getBondsContactsAndSingletons(selAtoms);
                                    var contacts = contactsAndSingletons["contacts"];
                                    var singletons = contactsAndSingletons["singletons"];
                                    if(style==="BONDS"||style==="FATBONDS"){
                                        var singletonPrimitiveInfo = singletonsToLinesInfo(singletons,4,atomColours);
                                        gl.appendOtherData(singletonPrimitiveInfo,true);
                                        var ligandAtoms =  model.getAtoms(selection+" and ligands");
                                        var multipleBonds = getMultipleBonds(ligandAtoms,enerLib,4,atomColours,"LINES");
                                        for(var imbo=0;imbo<multipleBonds.length;imbo++){
                                            gl.appendOtherData(multipleBonds[imbo]);
                                        }
                                        if(ligandAtoms.length<selAtoms.length){
                                            var linePrimitiveInfo = contactsToLinesInfo(contacts,4,atomColours);
                                            gl.appendOtherData(linePrimitiveInfo,true);
                                        }
                                    } else if(style==="THINBONDS"){
                                        var singletonPrimitiveInfo = singletonsToLinesInfo(singletons,2,atomColours);
                                        gl.appendOtherData(singletonPrimitiveInfo,true);
                                        var ligandAtoms =  model.getAtoms(selection+" and ligands");
                                        var multipleBonds = getMultipleBonds(ligandAtoms,enerLib,2,atomColours,"LINES");
                                        for(var imbo=0;imbo<multipleBonds.length;imbo++){
                                            gl.appendOtherData(multipleBonds[imbo]);
                                        }
                                        if(ligandAtoms.length<selAtoms.length){
                                            var linePrimitiveInfo = contactsToLinesInfo(contacts,2,atomColours);
                                            gl.appendOtherData(linePrimitiveInfo,true);
                                        }
                                    } else if(style==="BALLSTICK"){
                                        var spheres = atomsToSpheresInfo(selAtoms,0.4,atomColours);
                                        gl.appendOtherData(spheres,true);
                                        var ligandAtoms =  model.getAtoms(selection+" and ligands");
                                        var multipleBonds = getMultipleBonds(ligandAtoms,enerLib,0.2,atomColours,"CYLINDERS");
                                        for(var imbo=0;imbo<multipleBonds.length;imbo++){
                                            gl.appendOtherData(multipleBonds[imbo]);
                                        }
                                        if(ligandAtoms.length<selAtoms.length){
                                            var cylinderPrimitiveInfo = contactsToCylindersInfo(contacts,0.2,atomColours);
                                            gl.appendOtherData(cylinderPrimitiveInfo,true);
                                        }
                                    } else if(style==="CYLINDERS"){
                                        var spheres = atomsToSpheresInfo(selAtoms,0.2,atomColours);
                                        gl.appendOtherData(spheres,true);
                                        var ligandAtoms =  model.getAtoms(selection+" and ligands");
                                        var multipleBonds = getMultipleBonds(ligandAtoms,enerLib,0.2,atomColours,"CYLINDERS");
                                        for(var imbo=0;imbo<multipleBonds.length;imbo++){
                                            gl.appendOtherData(multipleBonds[imbo]);
                                        }
                                        if(ligandAtoms.length<selAtoms.length){
                                            var cylinderPrimitiveInfo = contactsToCylindersInfo(contacts,0.2,atomColours);
                                            gl.appendOtherData(cylinderPrimitiveInfo,true);
                                        }
                                    } else if(style==="SPHERES"){
                                        var spheres = atomsToSpheresInfo(selAtoms,1.0,atomColours);
                                        gl.appendOtherData(spheres,true);
                                    }
                                } else if(style==="SPLINE"||style==="WORM"){
                                    var flagBulge = true;
                                    CalcSecStructure(hier,flagBulge);
                                    //FIXME - pick up the colour scheme ....
                                    var colourScheme = new ColourScheme(pdbatoms);
                                    //var atomColours = colourScheme.colourByAtomType();
                                    var atomColours = colourScheme.colourBySecondaryStructure({"strand":[0.0,0.0,1.0,1.0],"helix":[1.0,0.0,0.0,1.0]});
                                    //FIXME - ignores selection. Can we select in spline yet?
                                    var coloured_splines_info;
                                    if(pdbatoms["models"].length>0){
                                        for(var modidx=0;modidx<pdbatoms["models"].length;modidx++){
                                            if(style==="SPLINE"){
                                                coloured_splines_info = GetSplinesColoured(pdbatoms,atomColours,false,pdbatoms["models"][modidx]);
                                            } else {
                                                coloured_splines_info = GetWormColoured(pdbatoms,atomColours,false,pdbatoms["models"][modidx]);
                                            }
                                            for(var i=0;i<coloured_splines_info.length;i++){
                                                gl.appendOtherData(coloured_splines_info[i],true);
                                            }
                                        }
                                    } else {
                                        if(style==="SPLINE"){
                                            coloured_splines_info = GetSplinesColoured(pdbatoms,atomColours);
                                        } else {
                                            coloured_splines_info = GetWormColoured(pdbatoms,atomColours);
                                        }
                                        for(var i=0;i<coloured_splines_info.length;i++){
                                            gl.appendOtherData(coloured_splines_info[i],true);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    var mapData = data.getElementsByTagName("MapData");
                    for (var i=0; i < mapData.length; i++) {
                        var filedata = mapData[i].getElementsByTagName("filedata")[0];
                        var theData = getNodeText(filedata);
 
                        var isDiffNode = mapData[i].getElementsByTagName("isDifferenceMap")[0];
                        var isDiff = false;
                        if(typeof(isDiffNode)==='undefined'){
                        } else {
                            var isDiffText = getNodeText(isDiffNode);
                            if(isDiffText==="True"||isDiffText==="true"||isDiffText==="1"){
                                isDiff = true;
                            }
                        }
                        var mapDisp = mapData[i].getElementsByTagName("MapDisp");
                        var sigma;
                        var units = mapData[i].getElementsByTagName("contourUnits")[0].innerHTML;
                        var contourMult = 1.0;

                        var sigmaNode = mapData[i].getElementsByTagName("sigma")[0];
                        if(typeof(sigmaNode)==='undefined'){
                        } else {
                            var sigmaText = getNodeText(sigmaNode);
                            sigma = parseFloat(sigmaText);
                        }

                        var unitsNode = mapData[i].getElementsByTagName("contourUnits")[0];
                        if(typeof(unitsNode)==='undefined'){
                        } else {
                            var unitsText = getNodeText(unitsNode);
                            if(unitsText==="sigma"){
                                contourMult = sigma;
                            }
                        }

                        var contents = _base64ToArrayBuffer(theData);
                        var map = readMapFromArrayBuffer(contents);
                        var mapGrid = mapToMapGrid(map);

                        for (var j=0; j < mapDisp.length; j++) {
                            if(isDiff){
                                var contourLevel = contourMult * parseFloat(mapDisp[j].getElementsByTagName("contourLevel")[0].innerHTML);
                                var mapTriangleData = {"mapColours": [[0.0,1.0,0.0,1.0]], "mapContourLevels":[contourLevel], "mapGrids":[mapGrid],"col_tri":[[]], "norm_tri":[[]], "vert_tri":[[]], "idx_tri":[[]] , "prim_types":[[]] };
                                gl.appendOtherDataIfReady(mapTriangleData);
                                var contourLevel = -contourMult * parseFloat(mapDisp[j].getElementsByTagName("contourLevel")[0].innerHTML);
                                var mapTriangleData = {"mapColours": [[1.0,0.0,0.0,1.0]], "mapContourLevels":[contourLevel], "mapGrids":[mapGrid],"col_tri":[[]], "norm_tri":[[]], "vert_tri":[[]], "idx_tri":[[]] , "prim_types":[[]] };
                                gl.appendOtherDataIfReady(mapTriangleData);
                            } else {
                                var contourLevel = contourMult * parseFloat(mapDisp[j].getElementsByTagName("contourLevel")[0].innerHTML);
                                var mapTriangleData = {"mapColours": [[0.0,0.0,1.0,1.0]], "mapContourLevels":[contourLevel], "mapGrids":[mapGrid],"col_tri":[[]], "norm_tri":[[]], "vert_tri":[[]], "idx_tri":[[]] , "prim_types":[[]] };
                                gl.appendOtherDataIfReady(mapTriangleData);
                            }
                        }
                    }
                }
                if(ccp4i2_body.getElementsByTagName("wizard").length>0){
                    var wizard = ccp4i2_body.getElementsByTagName("wizard")[0];
                    var template = wizard.getElementsByTagName("template")[0].innerHTML.split(":")[1];
                    var parameters = wizard.getElementsByTagName("parameters")[0];
                    var molIds = getElementsByRegex(parameters,'^MolData.*');
                    for (var i=0; i < molIds.length; i++) {
                        loadData(gl,atoms_ids[molIds[i].innerHTML],enerLib,true,template);
                    }
                }
                gl.buildBuffers();
                gl.drawScene();
                gl.setOrigin(hier[0].centre());

                var atoms = hier[0].getAllAtoms();
                var minX = 1e+8;
                var minY = 1e+8;
                var minZ = 1e+8;
                var maxX = 1e-8;
                var maxY = 1e-8;
                var maxZ = 1e-8;
                for(var iat=0;iat<atoms.length;iat++){
                    var x = atoms[iat].x();
                    var y = atoms[iat].y();
                    var z = atoms[iat].z();
                    if(x>maxX) maxX = x;
                    if(y>maxY) maxY = y;
                    if(z>maxZ) maxZ = z;
                    if(x<minX) minX = x;
                    if(y<minY) minY = y;
                    if(z<minZ) minZ = z;
                }
                extent = 0.0;
                if((maxX-minX)>extent) extent = maxX-minX;
                if((maxY-minY)>extent) extent = maxY-minY;
                if((maxZ-minZ)>extent) extent = maxZ-minZ;
                gl.setFog([-0.25*extent,0.75*extent]);
                gl.set_clip_range(-0.25*extent,extent*.25,true);
                
                if(ccp4i2_body.getElementsByTagName("View").length>0){
                    var View = ccp4i2_body.getElementsByTagName("View")[0];
                }
            } else {
                console.log("no ccp4i2_body");
            }
            //addSymmetry(pdbatoms);
            return extent;
            }

    window.onload = function()
    {
        loaded = true;
        for(var i = 0; i < queue.length; i++)
        {
            queue[i][0](queue[i][1],queue[i][2]);
        }
    }
    
</script>        """
                html = ""
                html += """<!DOCTYPE html>
<html>
<head>
<title>
CCP4I2 3D View
</title>\n"""
                html += onLoadText
                html += '<link rel="stylesheet" type="text/css" href="'+htmlBase+'/jquery-ui-1.10.4.custom.min.css">\n'
                for webgl_source in WEBGLSOURCES:
                    html += '<script src="'+htmlBase+'/MGWebGL/'+webgl_source+'">type="text/javascript"></script>\n'
                html += '<script src="'+htmlBase+'/jquery.min.js">type="text/javascript"></script>\n'
                html += '<script src="'+htmlBase+'/jquery-ui-1.10.4.custom.min.js">type="text/javascript"></script>\n'
                html += '<body style="font-family: Helvetica, Arial, sans-serif;">\n'
                html += "<table>\n"
                html += "<tr>\n"
                html += "<td colspan=\"2\">\n"
                html += '<div id="webgldiv" style="border: none; width: 600px; height: 600px; float: left;"></div>\n'
                html += "</td>\n"
                html += "</tr>\n"
                html += "<tr>\n"
                html += "<td>Fog</td>\n"
                html += '<td style="border: none; width: 500px;">\n'
                html += '<div id="fog_slider" style="border: none; width: 500px; float: left;"></div>\n'
                html += "</td>\n"
                html += "</tr>\n"
                html += "<tr>\n"
                html += "<td>Clip</td>\n"
                html += '<td style="border: none; width: 500px;">\n'
                html += '<div id="clip_slider" style="border: none; width: 500px; float: left;"></div>\n'
                html += "</td>\n"
                html += "</tr>\n"
                html += "</table>\n"

                tree = CCP4Utils.openFileToEtree(args["picdef"])
                MolDatas = tree.xpath("ccp4i2_body/scene/data/MolData")

                for m in MolDatas:
                    fn = m.xpath("filename")
                    with open(fn[0].text) as f:
                        fnt = f.read()
                    fd = etree.SubElement(m,"filedata")
                    fd.text = fnt.replace("<","&lt;").replace(">","&gt;")
                    m.remove(fn[0])

                MapDatas = tree.xpath("ccp4i2_body/scene/data/MapData")

                for m in MapDatas:
                    fn = m.xpath("filename")
                    fd = etree.SubElement(m,"filedata")
                    b64map, mean,std_dev = MTZToB64Map(fn[0].text)
                    fd.text = b64map
                    sigma = etree.SubElement(m,"sigma")
                    sigma.text = str(std_dev)
                    mean = etree.SubElement(m,"mean")
                    mean.text = str(mean)
                    m.remove(fn[0])

                sceneText = etree.tostring(tree,pretty_print=True).decode("utf-8")

                fog_slider_range_uuid_str = "fog_slider"
                clip_slider_range_uuid_str = "clip_slider"

                scriptText = """<script>
    function webglfun() {

        var gl = new MGWebGL("webgldiv",true,false); 
        var enerLib = new EnerLib();

        gl.setBackground([1.0,1.0,1.0,1.0]);
        var contents = """+'`'+sceneText+'`'+""";
        if (window.DOMParser) {
            var parser=new DOMParser();
            xmlDoc=parser.parseFromString(contents,"text/xml");
        } else {
            xmlDoc=new ActiveXObject("Microsoft.XMLDOM");
            xmlDoc.async=false;
            xmlDoc.loadXML(contents);
        }
        var extent = doParseAndLoad(gl,enerLib,xmlDoc);
        gl.setShowAxes(true);
  $( function() {
    $( "#"""+fog_slider_range_uuid_str+"""\" ).slider({
      range: true,
      min: -extent,
      max: extent,
      values: [ -0.25*extent,0.75*extent ],
      slide: function( event, ui ) {
        gl.setFog(ui.values);
      }
    });
    $( "#amount" ).val( "$" + $( "#"""+fog_slider_range_uuid_str+"""\" ).slider( "values", 0 ) +
      " - $" + $( "#"""+fog_slider_range_uuid_str+"""\" ).slider( "values", 1 ) );
  } );
  
  $( function() {
    $( "#"""+clip_slider_range_uuid_str+"""\" ).slider({
      range: true,
      min: -extent,
      max: extent,
      values: [ -0.25*extent,extent*.25 ],
      slide: function( event, ui ) {
        gl.set_clip_range(ui.values[0],ui.values[1],true);
      }
    });
    $( "#amount" ).val( "$" + $( "#"""+clip_slider_range_uuid_str+"""\" ).slider( "values", 0 ) +
      " - $" + $( "#"""+clip_slider_range_uuid_str+"""\" ).slider( "values", 1 ) );
  } );
        };
    enqueue(webglfun,"#"""+fog_slider_range_uuid_str+"\",\"#"+clip_slider_range_uuid_str+"""\");
    </script>
                """
                html += scriptText
                
                html += "</body>\n"

                fout = tempfile.NamedTemporaryFile(suffix=".html",delete=False)
                tfilename = fout.name
                print("Writing",tfilename); sys.stdout.flush()
                fout.write(html.encode())
                fout.close()

                self.web_view = QtWebEngineWidgets.QWebEngineView()
                self.web_view.load(QtCore.QUrl("file://"+tfilename))
                self.web_view.resize(630,660)
                self.web_view.show()
                self.web_view.raise_()
 
        if not self._channelIsSet:
            self.channel = QtWebChannel.QWebChannel()
            self.handler = CCP4WebReportBridge()
            self.channel.registerObject('handler', self.handler)
            self.webView.page().setWebChannel(self.channel)
            self.handler.jsSignal.connect(actOnJSSignal)
            self._channelIsSet = True

        self.generator = None
        self._diagnostic = True

    def exportData(self,myFileName,fileId,exportFileName):
      if os.path.splitext(exportFileName) !=  os.path.splitext(myFileName):
        exportFileName = os.path.splitext(exportFileName)[0] +  os.path.splitext(myFileName)[1]
      try:
        shutil.copyfile(myFileName,exportFileName)
      except:
        from . import CCP4Widgets
        e = CException(CCP4Widgets.CDataFileView,100,'From: '+str(myFileName)+' to: '+str(exportFileName),name=self.modelObjectPath())
        warningMessage(e, 'Copying file',parent=self)
      else:
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().createExportFile(fileId=fileId,exportFilename=exportFileName)
        fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['jobid','projectname'])
        from ..dbapi import CCP4DbUtils
        CCP4DbUtils.makeJobBackup(jobId=fileInfo['jobid'],projectName=fileInfo['projectname'])

    def setReportLinks(self,labelList=[]):
        # Set the self.linksFrame[1] which is the option for finished reports
        for widget in self.linksFrame[1].findChildren(QtWidgets.QPushButton):
            self.linksFrame[1].layout().removeWidget(widget)
            widget.deleteLater()
        for label,brief in labelList:
            if brief is not None:
                but = QtWidgets.QPushButton(brief,self)
                but.setToolTip(label)
            else:
                but = QtWidgets.QPushButton(label,self)
            but.setFlat(True)
            self.linksFrame[1].layout().addWidget(but)
            but.clicked.connect(functools.partial(self.handleLinkButton,label))

    def clear(self):
        self.webView.clear()
    
    def setLinkButtons(self):
        try:
            linksObj= self.webView.report.getDataObject(id='links')
        except:
            pass

    def setNextButtons(self,taskName=None,jobId=None,status=None,interruptLabel=None):
        #print 'CReportView.setNextButtons',taskName,jobId,status,interruptLabel
        buts = self.nextFrame.findChildren(QtWidgets.QPushButton)
        for but in buts:
            but.hide()
            but.deleteLater()
        if status is not None:
            if status == 'Finished' or status == 'Unsatisfactory':
                nextList = TASKMANAGER().whatNext(taskName,jobId)
                ii = 0
                for link,label,defFile in nextList:
                    but = QtWidgets.QPushButton(label,self)
                    self.nextFrame.layout().addWidget(but,ii/3,ii%3)
                    but.clicked.connect(functools.partial(self.nextTask.emit,link,defFile))
                    ii = ii + 1

#I'm protecting this in try/except because I am trying to be very careful not to break other the whole program.
                try:
                    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
                    outputXmlFile = PROJECTSMANAGER().makeFileName(jobId=self.openJob.jobId,mode='PROGRAMXML')
                    parser = etree.XMLParser()
                    f = open(outputXmlFile,'rb')
                    s = f.read()
                    f.close()
                    tree = etree.fromstring(s, parser)
                    suggested = tree.xpath("Verdict/suggestedParameters")
                    if len(suggested)>0 and len(suggested[0])>0:
                        but = QtWidgets.QPushButton("Re-run with suggested parameters",self)
                        self.nextFrame.layout().addWidget(but,ii/3+1,0)
                        but.clicked.connect(functools.partial(self.nextTask.emit,'clone_suggested',None))
                except:
                    pass
                    """
                    print("Problem adding the clone with suggested button")
                    exc_type, exc_value, exc_tb = sys.exc_info()[:3]
                    sys.stderr.write(str(exc_type) + '\n')
                    sys.stderr.write(str(exc_value) + '\n')
                    traceback.print_tb(exc_tb)
                    """

            elif status == 'Failed':
                but = QtWidgets.QPushButton('Clone job to rerun',self)
                self.nextFrame.layout().addWidget(but,0,0)
                but.clicked.connect(functools.partial(self.nextTask.emit,'clone_rerun',None))
            elif status in [ 'Running','Running remotely'] and taskName is not None:
                interruptLabel = TASKMANAGER().getTaskAttribute(taskName=taskName,attribute='INTERRUPTLABEL',default=None,script=True)
                #print 'setNextButtons interruptLabel',interruptLabel
                if interruptLabel is not None:
                    but = QtWidgets.QPushButton(interruptLabel,self)
                    self.nextFrame.layout().addWidget(but,0,0)
                    but.clicked.connect(functools.partial(self.handleInterruptButton,jobId))

    @QtCore.Slot(str)
    def handleInterruptButton(self,jobId):
        self.setNextButtons(jobId=jobId)
        self.interrupt.emit(jobId)

    @QtCore.Slot(str)
    def handleLinkButton(self,link):
        #print 'handleLinkButton',link
        url = self.webView.url()
        url.setFragment(link)
        # This should not be necessary..
        #self.webView.setTarget(link)
        self.webView.load(url)
        
    @QtCore.Slot('QUrl')
    def handleNavigationRequest(self,url):
        # Intercept navigation requests and open links in another window 
        print('CReportView.handleNavigationRequest',str(url.toLocalFile()),str(url.fragment()))
        sys.stdout.flush()
        path = str(url.path())
        print('CReportView.handleNavigationRequest 2'); sys.stdout.flush()
        if str(url.toLocalFile()) == str(CCP4Utils.getCCP4I2Dir())+"/docs/":
            # probably an internal link so show in the report window
            newUrl = QtCore.QUrl.fromLocalFile(self._reportFile)
            newUrl.setFragment(url.fragment())
            self.webView.load(newUrl)
            return
        print('CReportView.handleNavigationRequest 3'); sys.stdout.flush()
        print(len(path), path[-11:])
        print(path[-11:]== 'report.html')
        print(path.endswith('report.html'))
        if len(path)>11 and path[-11:] == 'report.html':
            if str(url.fragment()) != "":
                self.webView.load(url)
                return
            jobId = self.webView.report.getLinkId(path)
            print(jobId)
            if jobId is not None:
                openJob = COpenJob(jobId=jobId)
                self.openJob.jobStatusUpdated.connect(self.handleJobStatusUpdated)
                self.openJob.workflowJobStatusUpdated.connect(self.handleWorkflowJobStatusUpdated)
                if not os.path.exists(path):
                    # Requesting report file that is not yet created
                    if self.generator is None:
                        from ..report import CCP4ReportGenerator
                        self.generator = CCP4ReportGenerator.CReportGenerator(jobId=openJob.jobId,jobStatus=openJob.status,jobNumber=openJob.jobNumber)
                        #self.generator.FinishedPictures.connect(self.handleMgFinished)
                    try:
                        reportFile, newPageOrNewData = self.generator.makeReportFile()
                    except CException as e:
                        if  e.maxSeverity()>Severity.WARNING:
                            warningMessage(e, windowTitle=self.parent().windowTitle(),message='Failed creating job report',parent=self)
                    except Exception as e:
                        QtWidgets.QMessageBox.warning(self,self.parent().windowTitle(),'Unknown error creating report file for job number '+str(openJob.jobNumber))
                    if os.path.exists(reportFile):
                        err = self.generator.mergeIntoParent(parentFile=self._reportFile)
                        if err.maxSeverity() <= Severity.WARNING:
                            self.webView.reload()
                return
        print('CReportView.handleNavigationRequest 4',url); sys.stdout.flush()
        # See if the web browser can do something
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().loadPage(url=url,newTab=True)
        WEBBROWSER().raise_()
        print('CReportView.handleNavigationRequest end')
        sys.stdout.flush()

    @QtCore.Slot('QUrl')
    def handleCustomMimeTypeRequested(self,url):
        #print 'CReportView.handleCustomMimeTypeRequest',url
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().CustomMimeTypeRequested(url)

    def showOutput(self,openJob=None,reportFile=None,reportErr=True,redo=False,doReload=False):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        from ..report import CCP4ReportGenerator
        #print 'showOutput',openJob,reportFile,redo
        #traceback.print_stack(limit=8)
        # If not report file then try making it from an outputXmlFile
        # Otherwise use a log or diagnostic file
        if openJob is None: openJob = self.openJob
        linkList = None
        newPageOrNewData = "NEWPAGE"
        # Do we have a log file ot fall back to?
        logFile = PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='LOG')
        if not os.path.exists(logFile):
            logFile = logFile + '.html'
            if not os.path.exists(logFile):
                logFile = None
        #print 'showOutput logFile',logFile,'reportFile',reportFile
        if reportFile is None:
            hasReportDef = TASKMANAGER().hasReportDefinition(openJob.taskname, openJob.status)
            #print 'showOutput hasReportDefinition',openJob.taskname, openJob.status,hasReportDef
            if openJob.isWorkflow:
                childOpenJob = openJob.childOpenJob(-1)
                #print 'showOutput workflow child',childOpenJob
                if childOpenJob is not None: openJob = childOpenJob
            #print 'showOutput',openJob.status,CCP4DbApi.FINISHED_JOB_STATUS,openJob.status in CCP4DbApi.FINISHED_JOB_STATUS
            if openJob.status == 'Failed':
                self.generator = CCP4ReportGenerator.CReportGenerator(jobId=openJob.jobId,jobStatus=openJob.status,jobNumber=openJob.jobNumber)
                reportFile = self.generator.makeFailedReportFile(redo=redo)
            elif  (openJob.status in CCP4DbApi.FINISHED_JOB_STATUS or (openJob.status in ['Running','Running remotely'] ) ) and hasReportDef:
                if self.generator is None or openJob.jobId != self.generator.jobId:
                    self.generator = CCP4ReportGenerator.CReportGenerator(jobId=openJob.jobId,jobStatus=openJob.status,jobNumber=openJob.jobNumber)
                    #self.generator.FinishedPictures.connect(self.handleMgFinished)
                else:
                    self.generator.setJobStatus(openJob.status)
                # Comment out to ensure errors are trapped
                if CONFIG().developer:
                    reportFile, newPageOrNewData = self.generator.makeReportFile(redo=redo,doReload=doReload,useGeneric=(logFile is None))
                else:
                    try:
                        #print 'showOutput makeReportFile'
                        if openJob.status == 'Failed':
                            reportFile = self.generator.makeFailedReportFile(redo=redo)
                        else:
                            reportFile, newPageOrNewData = self.generator.makeReportFile(redo=redo,doReload=doReload,useGeneric=(logFile is None))
                    except CException as e:
                        # Dont report lack for report definition file
                        if reportErr and e.maxSeverity()>Severity.WARNING and e.code != 3:
                            warningMessage(e, windowTitle=self.parent().windowTitle(),message='Failed creating job report',parent=self)
                        reportFile = None
                    except Exception as e:
                        if reportErr:
                            QtWidgets.QMessageBox.warning(self,self.parent().windowTitle(),'Unknown error creating report file for job number '+str(openJob.jobNumber))
                        reportFile = None
            if reportFile is None:
                reportFile = PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='REPORT')
                if not os.path.exists(reportFile):
                    if logFile is not None:
                        reportFile = logFile
                    else:
                        reportFile = PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='DIAGNOSTIC')
                        if not os.path.exists(reportFile):
                            reportFile = None
        if reportFile is None:
            reportFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'docs','report_files','blank.html')
        #print 'CProjectViewer.showOutput reportFile',openJob.jobId,openJob.jobnumber,openJob.taskname,reportFile
        if newPageOrNewData == 'NEWDATA' and reportFile is not None:
            #Here if we believe there is an existing report in the webview and we simply want new data
            #to be loaded into it
#We cannot get nElementsUpdated synchronously with QtWebEngine
#So I guess we're going to have to have use QWebChannel
            loadedUrl=self.webView.url()
            reportUrl =  QtCore.QUrl.fromLocalFile(reportFile)

            if not self._channelIsSet:
                self.channel = QtWebChannel.QWebChannel()
                self.handler = CCP4WebReportBridge()
                self.channel.registerObject('handler', self.handler)
                self.webView.page().setWebChannel(self.channel)
                self._channelIsSet = True
            
            if reportUrl == loadedUrl:
                def reloadReportDataResultHandler(total):
                     if total>0:
                         self.showOutput0(openJob,reportFile,linkList,True)
                self.handler.reloadDataResult.connect(reloadReportDataResultHandler)
                self.webView.page().runJavaScript('''
        try {
            reloadReportData();
        } catch(e) {
            window.reloadDataBridge.reloadReportDataResult(0);
        }
                ''');
                return reportFile
            else:
                return self.showOutput0(openJob,reportFile,linkList,False)
        else:
            return self.showOutput0(openJob,reportFile,linkList,False)

    @QtCore.Slot(object,str,bool)
    def showOutput0(self,openJob,reportFile,linkList,reloadedData):
        if not reloadedData:
            if reportFile is None or os.path.splitext(reportFile)[1] in ['.html','.htm']:
                self.setCurrentIndex(0)
                if reportFile is None:
                    self.webView.load( QtCore.QUrl())
                    self.linksStack.setCurrentIndex(0)
                else:
                    if os.path.basename(reportFile) == "report.html":
                        #FIXME - This causes server problems if we view sub jobs.
                        url = QtCore.QUrl("http://127.0.0.1:43434/database/?getProjectJobFile?projectId="+openJob.projectid+"?fileName="+os.path.basename(reportFile)+"?jobNumber="+openJob.jobnumber)
                        #print(url)
                    else:
                        url =  QtCore.QUrl.fromLocalFile(reportFile)
                    loadedUrl=self.webView.url()
#This is a bit of a hack to try to make sure that arcimboldo.html/PAIREF_project.html is updated if open.
                    if (os.path.dirname(url.path()) == os.path.dirname(loadedUrl.path()) and os.path.basename(url.path()) == "report_tmp.html" and os.path.basename(loadedUrl.path()) == "arcimboldo.html") or (os.path.dirname(url.path()) == os.path.dirname(os.path.dirname(loadedUrl.path())) and os.path.basename(url.path()) == "report_tmp.html" and os.path.basename(loadedUrl.path()) == "PAIREF_project.html"):
                       self.webView.load(loadedUrl)
                    else:
                       self.webView.load(url)
                    if openJob.status == 'Failed':
                        self.linksStack.setCurrentIndex(2)
                    else:
                        if linkList is not None:
                            self.setReportLinks(linkList)
                        self.linksStack.setCurrentIndex(1)
            else:
                # Xml file - treat as text
                self.linksStack.setCurrentIndex(0)
                if self.textView is None:
                    from . import CCP4TextViewer
                    self.textView = CCP4TextViewer.CTextViewer(self)
                    self.addWidget(self.textView)
                self.setCurrentIndex(1)
                try:
                    self.textView.open( reportFile)
                    self.textView.setFont(style='fixed_width')
                except:
                    pass
        #print 'CReportView.showOutput',reportFile
        self._reportFile = reportFile
        if openJob.jobId != self.openJob.jobId:
            FILEWATCHER().removeJob(self.openJob.jobId)
            self.addWatchFile(openJob)
            self.openJob = openJob
            self.openJob.jobStatusUpdated.connect(self.handleJobStatusUpdated)
            self.openJob.workflowJobStatusUpdated.connect(self.handleWorkflowJobStatusUpdated)
        self.reportAvailable.emit(self.openJob.jobId, (reportFile is not None))
        return reportFile

#Adapted from https://stackoverflow.com/questions/55963931/how-to-download-csv-file-with-qwebengineview-and-qurl
    @QtCore.Slot("QWebEngineDownloadItem*")
    def handleDownload(self,download):
        old_path = os.path.basename(download.url().path())
        suffix = QtCore.QFileInfo(old_path).suffix()

        fileBrowser = QtWidgets.QFileDialog()
        fileBrowser.setWindowTitle("Save file")
        fileBrowser.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        fileBrowser.setFileMode(QtWidgets.QFileDialog.AnyFile)
        fileBrowser.setDefaultSuffix(suffix)
        fileBrowser.selectFile(old_path)
        fileBrowser.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        if fileBrowser.exec_():
           if len(fileBrowser.selectedFiles()) > 0:
               path = fileBrowser.selectedFiles()[0]
               if path:
                   download.setPath(path)
                   download.accept()

    @QtCore.Slot(bool)
    def handleLoadFinished(self,ok):
        if not ok:
            return
        if hasattr(self.webView.report,"topFolds"):
            self.setReportLinks(self.webView.report.topFolds)

    @QtCore.Slot(str)
    def handleMgFinished(self,jobId):
        if jobId == self.openJob.jobId:
            self.webView.attemptImageLoad()

    @QtCore.Slot(tuple)
    def handleFileChanged(self,args):
        if self._diagnostic:
            print('CReportView.handleFileChanged', args)
            sys.stdout.flush()
        jobId,outputXmlFile = args
        if jobId == self.openJob.jobId:
            self.showOutput(redo=True)
        elif self.openJob.childjobs is not None and jobId in self.openJob.childjobs:
            self.showOutput(openJob= self.openChildJob,redo=True)

    @QtCore.Slot(str)
    def handleJobStatusUpdated(self,status):
        if status in ['Running']:
            self.addWatchFile(self.openJob)
        elif status in ['Finished','Failed']:
            FILEWATCHER().removeJob(self.openJob.jobId)
            #self.showOutput(redo=True)

    @QtCore.Slot(tuple)
    def handleWorkflowJobStatusUpdated(self,argList):
        if not argList[0] == self.openJob.jobId:
            return
        # NEED TO REMOVE FILEWATCH on subjobs too!
        FILEWATCHER().removeJob(self.openJob.jobId)
        if len(self.openJob.childjobs) > 1:
            FILEWATCHER().removeJob(self.openJob.childjobs[-2])
        self.openChildJob = COpenJob(jobId=argList[1])
        self.addWatchFile(self.openChildJob)
        try:
            self.showOutput(openJob=self.openChildJob,redo=True)
        except:
            pass

    def addWatchFile(self,openJob):
        if not openJob.status in ['Running', 'Running remotely']:
            return
        # If there is specification for a 'running' report then watch the task program output
        runningXrt = TASKMANAGER().searchXrtFile(openJob.taskName,jobStatus='Running')
        if runningXrt is None:
            runningXrt = TASKMANAGER().getReportClass(openJob.taskName,jobStatus='Running')
        if runningXrt is not None:
            updateInterval = TASKMANAGER().getReportAttribute(openJob.taskName,'UPDATE_INTERVAL')
            watchFile = TASKMANAGER().getReportAttribute(openJob.taskName,'WATCHED_FILE')
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            if watchFile is None:
                outputXmlFile = PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='PROGRAMXML')
            else:
                outputXmlFile = os.path.normpath(os.path.join(PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='ROOT'),watchFile))
            FILEWATCHER().addJobPath(jobId=openJob.jobId,path=outputXmlFile,updateInterval=updateInterval)


class CTaskInputFrame(QtWidgets.QFrame):

    lastOpenFolder = {}
    ERROR_CODES = {101 : { 'description' : 'Error attempting to run internal plugin' }}

    def __init__(self,parent):
        QtWidgets.QFrame.__init__(self,parent)
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(0,CProjectViewer.MARGIN,0,CProjectViewer.MARGIN)
        self.layout().setSpacing(CProjectViewer.MARGIN)
        self.taskWidget = None
        from ..core.CCP4Preferences import PREFERENCES
        from ..qtcore.CCP4JobController import JOBCONTROLLER
        PREFERENCES().TASK_WINDOW_LAYOUT.dataChanged.connect(self.redraw)
        JOBCONTROLLER().serverJobFailed.connect(self.handleServerJobFail)

    def createTaskWidget(self,taskName,projectId=None,jobId=None,container=None,taskEditable=True,followJobId=None,excludeInputData=False):
        # Create task widget
        taskWidgetClass = TASKMANAGER().getTaskWidgetClass(taskName)
        if taskWidgetClass is not None:
            try:
                taskWidget = taskWidgetClass(self)
            except CException as e:
                warningMessage(e, 'Error opening task window: '+str(taskName),parent=self)
                raise e
            except Exception as e:
                mess = QtWidgets.QMessageBox.warning(self,'Error opening task window: '+str(taskName),str(e))
                raise CException(self.__class__,101,taskName)
            taskWidget.folderAttributes.setAttribute(attribute='editable',folderFunction='all',value=taskEditable)
            taskWidget.setContainer(container)
        else:
            # Try for an auto-generated gui?
            from . import CCP4ContainerView
            taskWidget = CCP4ContainerView.CContainerView(self,container=container)
            taskWidget.folderAttributes.setAttribute(attribute='editable',folderFunction='all',value=taskEditable)
            #  raise CException(self.__class__,103,taskName)    
        taskWidget.setProjectId(projectId)
        taskWidget.setJobId(jobId)
        taskWidget.setDefaultParameters()
        taskWidget.excludeInputData = excludeInputData
        e = taskWidget.draw()
        if len(e) > 0:
            print('drawTaskWidget errors', e.report())
        if followJobId is not None:
            taskWidget.setFollowJobId(followJobId)
        taskWidget.setDefaultFiles()
        if self.taskWidget is not None:
            CTaskInputFrame.lastOpenFolder[self.taskWidget.jobId()] = self.taskWidget.visibleFolder()
        self.closeTaskWidget()
        # Stick new widget in window
        self.taskWidget = taskWidget
        if CTaskInputFrame.lastOpenFolder.get(jobId,None) is not None:
            self.taskWidget.setVisibleFolder(CTaskInputFrame.lastOpenFolder[jobId])
        self.layout().insertWidget(1,self.taskWidget)
        self.taskWidget.show()
        invalidList = self.taskWidget.isValid()
        return self.taskWidget

    def closeTaskWidget(self):
        if self.taskWidget is not None:
            # Close any exisiting task widget. Probably should be saving a def file
            self.taskWidget.close()
            self.taskWidget.deleteLater()
            self.taskWidget = None

    @QtCore.Slot(str)
    def runTask(self,runMode):
        # Update control file path in db and set status to queue the job
        #print 'CTaskInputFrame.runTask',runMode
        if isinstance(runMode,QtWidgets.QAction): runMode = runMode.text().__str__()
        if self.taskWidget is None: return
        self.taskWidget.updateModelFromView()
        # The method to fix any issues in user input
        rv = self.taskWidget.fix()
        if len(rv) > 0:
            warningMessage(rv, parent=self, windowTitle='Running task', message='Error in input data')
            return
        # Check validity
        invalidList = self.taskWidget.isValid()
        taskErrors = self.taskWidget.taskValidity()
        if len(invalidList)>0 or len(taskErrors)>0:
            text = ''
            if len(invalidList)>0:
                text = text + 'Can not run task - missing or invalid data\nInvalid items are highlighted on the interface\n'
                for item in invalidList:
                    try:
                        itemName = item.objectPath()
                    except:
                        itemName = 'Unknown'
                    if isinstance(item,str):
                        text = text + item + '\n'
                    else:
                        for xtra in ['label','toolTip']:
                            label = item.qualifiers(xtra)
                            if label is not None and label is not NotImplemented:
                                itemName = itemName +': '+label
                        text = text + itemName + '\n'
                    print('Invalid data value',itemName,item)
            if len(taskErrors)>0:
                text = text = taskErrors.report(user=True,ifStack=False)
            if not hasattr(self,'messageBox'):
                self.messageBox = QtWidgets.QMessageBox(self)
                self.messageBox.setWindowTitle(str(self.taskWidget.title()))
                but = self.messageBox.addButton('Close',QtWidgets.QMessageBox.RejectRole)
                but.clicked.connect(self.messageBox.close)
                self.messageBox.setDefaultButton(but)
                self.messageDetailsButton = self.messageBox.addButton('Details',QtWidgets.QMessageBox.AcceptRole)
                self.messageDetailsButton.clicked.connect(functools.partial(self.showInvalidDetails,invalidList))
            else:
                self.messageDetailsButton.show()
            self.messageBox.setText(text)
            self.messageBox.exec_()
            return
        jobId = self.taskWidget.jobId()
        projectId=self.taskWidget.projectId()
        taskName = self.taskWidget.taskName()
        container = self.taskWidget.getContainer()
        # Check if the task is blocked from local running
        print('blockLocal',taskName, TASKMANAGER().getTaskAttribute(taskName,'blockLocal'))
        if runMode.count('run_remote')==0 and TASKMANAGER().getTaskAttribute(taskName,'blockLocal'):
            title = TASKMANAGER().getTitle(taskName)
            QtWidgets.QMessageBox.warning(self,"Running "+title, "The task "+title+"\nand other large jobs should not be run on this computer.\nPlease use 'Run on server' to run on a bigger computer." )
            return
        # Remove unset items from lists
        container.removeUnsetListItems()
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        ifImportFile, errors = PROJECTSMANAGER().importFiles(jobId=jobId, container=container)
        rv = self.makeJobInputFile(jobId)
        if not rv:
            QtWidgets.QMessageBox.warning(self, str(self.taskWidget.title()), 'Job failed writing input parameters file')
            return
        # This is used in the 'workflow' automation 
        subJobWidgets = getattr( self.taskWidget, 'subJobTaskWidgets', {})
        for jobName, taskWidget in list(subJobWidgets.items()):
            taskWidget.saveToXml()
        # Preceeding job concept not now used to track job - this is pretty much redundant
        preceedingjobId = self.taskWidget.getContainer().guiAdmin.followFrom
        PROJECTSMANAGER().db().updateJob(jobId=jobId, key='preceedingjobid',value=preceedingjobId.pyType())
        #Record input files in database
        PROJECTSMANAGER().db().gleanJobFiles(jobId=jobId, container=container,projectId=projectId, roleList=[CCP4DbApi.FILE_ROLE_IN])
        if self.taskWidget.isEditor():
            PROJECTSMANAGER().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FINISHED)
        else:
            # Delete a pre-existing report - assume we are restarting job
            try:
                os.remove(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='REPORT'))
            except:
                pass
            if runMode.count('run_remote') > 0:
                self.runRemotely(jobId,projectId)
                return
            elif TASKMANAGER().isInternalPlugin(taskName):
                try:
                    PROJECTSMANAGER().runInternalTask(jobId=self.taskWidget.jobId(), projectId=self.taskWidget.projectId(),
                                                        taskName=self.taskWidget.taskName())
                except Exception as e:
                    err = CException(self.__class__,101,taskName,str(e))
                    warningMessage(err, 'Running internaltask','Failed running task',parent=self)
            else:
                PROJECTSMANAGER().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_QUEUED)
        self.redrawTaskWidget()

    def redrawTaskWidget(self,editable=False):
        # Redraw task gui as non-editable
        container = self.taskWidget.container
        taskName = self.taskWidget.taskName()
        jobId = self.taskWidget.jobId()
        projectId=self.taskWidget.projectId()
        from ..core.CCP4Preferences import PREFERENCES
        if isinstance(self.taskWidget,CCP4TaskWidget.CTaskWidget) and PREFERENCES().TASK_WINDOW_LAYOUT == 'FOLDER':
            scrollDisp = self.taskWidget.getScrollDisplacement()
            folderOpenStatus = self.taskWidget.widget.getFolderOpenStatus()
            taskWidget = self.createTaskWidget(taskName,container=container,taskEditable=editable,jobId=jobId,projectId=projectId)
            taskWidget.widget.setFolderOpenStatus(folderOpenStatus)
            taskWidget.setScrollDisplacement(scrollDisp)
        else:
            taskWidget = self.createTaskWidget(taskName,container=container,taskEditable=editable,jobId=jobId,projectId=projectId)

    @QtCore.Slot(list)
    def showInvalidDetails(self,invalidList):
        #print 'showInvalidDetails',invalidList
        text = 'Can not run task - missing or invalid data\nInvalid items are highlighted on the interface\n'
        for modelObj in invalidList:
            text = text + str(modelObj)+'\n'
            if modelObj is not None:
                text = text + modelObj.validity(modelObj.get()).report(ifStack=False) + '\n'
        self.messageDetailsButton.hide()
        self.messageBox.setText(text)
        self.messageBox.show()
          
    def makeJobInputFile(self, jobId=None):
        if self.taskWidget is None:
            return False
        #Beware -- this is trying to save status of previous task widget which may have been deleted 
        try:
          from ..core.CCP4ProjectsManager import PROJECTSMANAGER
          status =  PROJECTSMANAGER().db().getJobInfo(jobId=self.taskWidget._jobId,mode=['status'])
        except:
            print('makeJobInputFile NOT saving input_params.xml file - db query fail')
            return False
        if status not in ['Pending','Interrupted']:
            # Ensure do not overwrite a job that is already started
            print('makeJobInputFile NOT saving input_params.xml file')
            return False
        self.taskWidget.saveToXml()
        return True

    def makeJobBall(self,jobId,projectId,mechanism='ssh_shared'):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        jobNumber = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber'])
        projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)
        #dbxml = JOBCONTROLLER().getServerParam(jobId,'dbXml')
        jobDir = PROJECTSMANAGER().db().jobDirectory(jobId=jobId)
        if mechanism not in ['ssh_shared','test','qsub_local','qsub_shared']:
            dbxml = os.path.join( projectInfo['projectdirectory'],'CCP4_TMP','DATABASE'+str(int(time.time()))+'.db.xml')
        else:
            dbxml = os.path.join( jobDir, 'DATABASE.db.xml' )
            from ..qtcore.CCP4JobController import JOBCONTROLLER
            JOBCONTROLLER().setServerParam(jobId,'dbXml',dbxml)
        inputFilesList,inputFileIdList,fromJobList,errReport =  PROJECTSMANAGER().getJobInputFiles(projectDir=projectInfo['projectdirectory'],jobIdList=[jobId],jobNumberList=[jobNumber])
        print('runRemotely inputFilesList',inputFilesList,'fromJobList',fromJobList)
        print('runRemotely  errReport',errReport.report())
        fromJobIdList = []
        fromJobNumberList = []
        for item in fromJobList:
            fromJobIdList.append(item['jobid'])
            fromJobNumberList.append(item['jobnumber'])
        jobNumberList,errReport = PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml,jobList=[jobId],inputFileList=inputFileIdList,inputFileFromJobList=fromJobIdList)
        if errReport.maxSeverity()>Severity.WARNING:
            warningMessage(errReport, "title",'Error creating XML database file',parent=self)
            return False
        if mechanism in ['ssh_shared','qsub_local','qsub_shared']:
            self.runRemotely1(jobId,projectId)
        else:
            from ..qtcore import CCP4Export
            tarball = os.path.join( projectInfo['projectdirectory'],'CCP4_TMP','job_'+jobNumber+'_setup.ccp4db.zip')
            if os.path.exists(tarball):
                os.remove(tarball)
            self.exportThread = CCP4Export.ExportProjectThread(self,projectDir=projectInfo['projectdirectory'],dbxml=dbxml,target=tarball,jobList=[jobNumber],inputFilesList=inputFilesList,directoriesList=[],extraJobList=fromJobNumberList)
            self.exportThread.jobsWithInputFiles = [jobNumber]
            self.exportThread.finished.connect(functools.partial(self.runRemotely1,jobId,projectId))
            print('makeJobBall starting',jobId,projectId,tarball)
            self.exportThread.start()

    def runRemotely(self,jobId,projectId,message=None):
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        from .CCP4JobControlGui import JOBCONTROLLERGUI
        dialog = JOBCONTROLLERGUI()
        jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId,['jobnumber','taskname'])
        dialog.setWindowTitle('Run '+jobInfo['jobnumber']+ ' ' + TASKMANAGER().getTitle(jobInfo['taskname']))
        dialog.setInfo(message)
        rv = dialog.exec_()
        print('runRemotely',rv,dialog.valid(),dialog.get('mechanism'))
        if rv == QtWidgets.QDialog.Accepted and dialog.valid():
            from ..qtcore.CCP4JobController import JOBCONTROLLER
            JOBCONTROLLER().createServerParams(jobId,dialog.getParams())
            print('project viewer runRemotely runningReport',self.taskWidget.taskName(),TASKMANAGER().getReportClass(self.taskWidget.taskName(),jobStatus='Running'))
            JOBCONTROLLER().setServerParam(jobId,'runningReport', (TASKMANAGER().getReportClass(self.taskWidget.taskName(),jobStatus='Running') is not None) )
            self.makeJobBall(jobId,projectId,dialog.get('mechanism'))

    @QtCore.Slot(str,str)
    def runRemotely1(self,jobId,projectId):
        print('runRemotely1',jobId,projectId)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_QUEUED)
        self.redrawTaskWidget()

    @QtCore.Slot(tuple)
    def handleServerJobFail(self,args):
        jobId,projectId,exception = args
        print('handleServerJobFail',jobId,projectId,exception.report(user=True,ifStack=False))
        if projectId != self.taskWidget.projectId(): return
        if jobId == self.taskWidget.jobId():
            self.redrawTaskWidget(editable=True)
        message = 'Failed to start remote job\n'+exception.report(mode=2,user=True,ifStack=False)
        if exception[0]['code'] == 331:
            message = message + '\nPlease check that work directory exists on remote machine'
        self.runRemotely(jobId,self.taskWidget.projectId(),message)

    def clear(self):
        if self.taskWidget is None: return
        self.taskWidget.close()
        self.taskWidget.deleteLater()
        self.taskWidget = None
        
    @QtCore.Slot()
    def redraw(self):
        if self.taskWidget is None: return
        taskName = self.taskWidget.taskName()
        jobId = self.taskWidget.jobId()
        projectId = self.taskWidget.projectId()
        container= self.taskWidget.getContainer()
        taskEditable = self.taskWidget.folderAttributes.attribute('editable')
        # ? followJobId copied in container ?
        taskWidget = self.createTaskWidget(taskName,projectId=projectId,jobId=jobId,container=container,taskEditable=taskEditable)


class CTaskTitleBarLayout(QtWidgets.QHBoxLayout):

    def minimumSize(self):
        return QtCore.QSize(CCP4TaskWidget.WIDTH,25)

    def sizeHint(self):
        return self.minimumSize()

class CTaskTitleBar(QtWidgets.QFrame):

    MARGIN = 2

    def __init__(self,parent):
        QtWidgets.QFrame.__init__(self,parent)
        self.jobId= None
        self.setLayout(CTaskTitleBarLayout())
        self.layout().setContentsMargins(0,0,0,0)
        self.layout().setSpacing(0)
        self.title = QtWidgets.QLabel(self)
        self.title.setObjectName('jobTitle')
        self.layout().addWidget(self.title)
        self.status =  QtWidgets.QLabel(self)
        self.status.setObjectName('jobStatus')
        self.layout().addWidget(self.status)
        self.layout().setStretchFactor(self.title,3.0)
        self.layout().setStretchFactor(self.status,3.0)
        movie = QtGui.QMovie(os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','running_1.gif')))
        self.icon = QtWidgets.QLabel(self)
        movie.setScaledSize(QtCore.QSize(24,24))
        self.icon.setMovie(movie)
        movie.start()
        self.icon.setObjectName('jobStatus')
        self.layout().addWidget(self.icon)
        self.llfont = 10  # Lower & upper font size limits for bar
        self.ulfont = 20
        self.WinScale = 0.5

    def setOpenJob(self, openJob):
        # Sets up the blue bar with Job status
        self.ftsize = self.llfont
        from ..core.CCP4Preferences import PREFERENCES
        if PREFERENCES().GUI_FONT_SIZE >= self.llfont:
            self.ftsize = min(PREFERENCES().GUI_FONT_SIZE, self.ulfont)
        if sys.platform == "win32":
            self.ftsize = round(self.ftsize*self.WinScale) # Win Fonts not the same size.
        fnt = QtGui.QFont("Sans Serif", self.ftsize, QtGui.QFont.Bold, True)
        self.title.setFont(fnt);
        self.title.setText(openJob.title)
        self.jobId = openJob.jobId
        self.setStatusBar({'jobId' : openJob.jobId, 'status' : openJob.info['status'] } )

    @QtCore.Slot(dict)
    def setStatusBar(self, info={}):
        if info.get('jobId') != self.jobId:
            return
        status = info.get('status',0)
        if isinstance(status,int): 
            status = CCP4DbApi.JOB_STATUS_TEXT[status]
        fnt = QtGui.QFont("Sans Serif", self.ftsize, QtGui.QFont.Bold, True)
        self.status.setFont(fnt);
        self.status.setText('The job is '+status)
        if status in  ['Running','Running remotely']:
            self.icon.show()
        else:
            self.icon.hide()


class CTaskFrame(QtWidgets.QFrame):

    launchJobRequestSignal = QtCore.Signal(str,dict)
    labelEdited =  QtCore.Signal()
    report71 = QtCore.Signal(object)

    INPUT_TAB = 0
    OUTPUT_TAB = 1
    COMMENT_TAB =2
    MARGIN = 0
    ERROR_CODES = {100 : {'description' : 'Unknown error drawing task widget'}}

    def __init__(self,parent,projectId=None):
        QtWidgets.QFrame.__init__(self,parent)
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(CTaskFrame.MARGIN,CTaskFrame.MARGIN,CTaskFrame.MARGIN,CTaskFrame.MARGIN)
        self.layout().setSpacing(CTaskFrame.MARGIN)
        self.openJob = COpenJob(projectId=projectId)
        self.titleBar = CTaskTitleBar(self)
        self.layout().addWidget(self.titleBar)
        self.tab = QtWidgets.QTabWidget(self)
        self.tab.currentChanged.connect(self.setZoomActionsEnabled)
        self.inputFrame = CTaskInputFrame(self)
        self.tab.addTab(self.inputFrame,'Input')
        self.outputFrame = CReportView(self)
        self.tab.addTab(self.outputFrame,'Results')
        self.statusFrame = CJobStatusWidget(self)
        self.tab.addTab(self.statusFrame,'Comments')
        self.setTaskTab('input')
        self.layout().addWidget(self.tab)
        bottomLayout = QtWidgets.QHBoxLayout()
        bottomLayout.setContentsMargins(0,0,0,0)
        bottomLayout.setSpacing(0)
        self.layout().addLayout(bottomLayout)
        self.buttons = CTaskButtons(self,parentLayout=bottomLayout)
        bottomLayout.addStretch()
        self.buttons.button('run').clicked.connect(functools.partial(self.inputFrame.runTask,'Now'))
        self.buttons.button('run').clicked.connect(functools.partial(self.buttons.button('run').setDefault,False))
        if ALWAYS_SHOW_SERVER_BUTTON:
            self.buttons.button('run').menu().triggered.connect(self.inputFrame.runTask)
        self.buttons.button('view').menu().triggered.connect(self.handleViewTask)
        self.buttons.button('task_menu').clicked.connect(functools.partial(self.window().showTaskChooser,True))
        self.outputFrame.reportAvailable.connect(self.handleReportAvailable)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)
        PROJECTSMANAGER().db().jobStarted.connect(self.handleJobStarted)
        self.outputFrame.labelEdited.connect(self.labelEdited.emit)

    def saveStatus(self,openJob=None):
        if hasattr(openJob,"jobId"):
            self.inputFrame.makeJobInputFile(openJob.jobId)
        else:
            self.inputFrame.makeJobInputFile()
        from ..dbapi import CCP4DbUtils
        if openJob is None:
            openJob = self.openJob
        if openJob.jobId is not None:
            try:
                CCP4DbUtils.makeJobBackup(jobId=openJob.jobId,projectName=openJob.projectName)
            except Exception as e:
                pass

    @QtCore.Slot()
    def setZoomActionsEnabled(self):
        ar = self.parent().findChild(QtWidgets.QAction, 'resetZoom')
        ap = self.parent().findChild(QtWidgets.QAction, 'zoomIn')
        am = self.parent().findChild(QtWidgets.QAction, 'zoomOut')
        if hasattr(self.parent(),"parent") and hasattr(self.parent().parent(),"parent") and hasattr(self.parent().parent().parent(),"findChild"):
            ar = self.parent().parent().parent().findChild(QtWidgets.QAction, 'resetZoom')
            ap = self.parent().parent().parent().findChild(QtWidgets.QAction, 'zoomIn')
            am = self.parent().parent().parent().findChild(QtWidgets.QAction, 'zoomOut')
        if not ar:
            return
        if self.tab.currentIndex() == 1:
            ar.setEnabled(True)
            ap.setEnabled(True)
            am.setEnabled(True)
        else:
            ar.setEnabled(False)
            ap.setEnabled(False)
            am.setEnabled(False)

    def setTaskTab(self,mode=None):
        if mode is None: return
        mode = mode.lower()
        if mode == 'status':
            self.tab.setCurrentIndex(self.COMMENT_TAB)
        elif mode == 'input':
            self.tab.setCurrentIndex(self.INPUT_TAB)
        elif mode == 'output':
            self.tab.setCurrentIndex(self.OUTPUT_TAB)

    @QtCore.Slot(dict)
    def handleJobStarted(self,args):
        if args.get('jobId','') != self.openJob.jobId and args.get('parentJobId','') != self.openJob.jobId :
            return
        self.outputFrame.setNextButtons(args.get('taskName',None), args['jobId'],status='Running')

    @QtCore.Slot(dict)
    def handleJobFinished(self,args):
        jobId = args.get('jobId','')
        status = args.get('status',None)
        if jobId == self.openJob.jobId:
            if isinstance(status,int):
                status = CCP4DbApi.JOB_STATUS_TEXT[status]
            if status != self.openJob.status:
                self.openJob.status = status
                self.buttons.setRunMode(status=self.openJob.status)
                self.buttons.setEnabled(self.openJob.status)
                self.window().updateActionEnabled()
                self.outputFrame.showOutput(self.openJob,redo=True)
                self.outputFrame.setNextButtons(self.openJob.taskname,self.openJob.jobId,status)
                self.outputFrame.setLinkButtons()
        elif args.get('parentJobId','') == self.openJob.jobId:
            # Have finished a sub-job that might have been interruptable so clear the 'next' buttons
            self.outputFrame.setNextButtons(jobId=jobId)


    def updateTaskFrame(self, openJob=None):
        reportFile = None
        if openJob is not None:
            self.openJob = openJob
        self.statusFrame.setJob(self.openJob)
        if self.openJob.jobId is None:
            self.buttons.button('next').menu().clear()
            self.buttons.setEnabled(self.openJob.status)
            self.inputFrame.clear()
            self.outputFrame.clear()
        else:
            try:
                reportFile = self.outputFrame.showOutput(self.openJob, reportErr=False)
            except CException as e:
                print(e)
                warningMessage(e, self.window().windowTitle(), 'Error creating job report', parent=self)
            except Exception as e:
                print(e)
                CMessageBox(self,message='Error creating report for job number '+str(self.openJob.jobnumber),exception=e,openJob=self.openJob)
            self.buttons.setEnabled(self.openJob.status)
            if self.openJob.status == 'Running':
                from ..core.CCP4ProjectsManager import PROJECTSMANAGER
                runningSubJob = PROJECTSMANAGER().db().getRunningSubJob(jobId=self.openJob.jobId)
                if runningSubJob is not None:
                    self.outputFrame.setNextButtons(runningSubJob['taskName'],runningSubJob['jobId'],'Running')
                else:
                    self.outputFrame.setNextButtons(self.openJob.taskname,self.openJob.jobId,self.openJob.status)
            else:
                self.outputFrame.setNextButtons(self.openJob.taskname,self.openJob.jobId,self.openJob.status)
            self.outputFrame.setLinkButtons()
            if self.inputFrame.taskWidget is not None:
                self.buttons.setRunMode(editor=self.inputFrame.taskWidget.isEditor(),status=self.openJob.status)
        if (self.openJob.status in CCP4DbApi.FINISHED_JOB_STATUS or self.openJob.status in ['Running','Failed']) and reportFile is not None:
            self.setTaskTab('output')
            self.tab.setTabEnabled(self.OUTPUT_TAB, True)
        else:
            self.setTaskTab('input')
            if reportFile is None:
                self.tab.setTabEnabled(self.OUTPUT_TAB, False)
        self.window().updateActionEnabled(self.openJob.status)

    def handleViewTask(self,mode):
        if not isinstance(mode,str):
            mode = str(mode.text())
        from ..qtcore.CCP4Launcher import LAUNCHER
        if mode.count('4mg'):
            LAUNCHER().openInViewer(viewer='ccp4mg',jobId=self.openJob.jobId,projectId=self.openJob.projectId,guiParent=self)
        elif mode.count('oot'):
            LAUNCHER().openInViewer(viewer='coot_job',jobId=self.openJob.jobId,projectId=self.openJob.projectId,guiParent=self)

    @QtCore.Slot(str,bool)
    def handleReportAvailable(self,jobId,status):
        if jobId == self.openJob.jobId or (self.openJob.childjobs is not None and jobId in self.openJob.childjobs):
            print("Report available?")
            pass
        else:
            print("************************************************************")
            print("Report not available?")
            print(jobId)
            print(self.openJob.jobId)
            print(self.openJob.childjobs)
            print("************************************************************")
            return
        if status:
            self.tab.setTabEnabled(self.OUTPUT_TAB,True)
            self.setTaskTab('output')
        else:
            self.tab.setTabEnabled(self.OUTPUT_TAB,False)

    def openTask(self,taskName=None,jobId=None,cloneJobId=None,followJobId=None,patchParamsFile=None,suggestedParams=None):
        from ..core import CCP4Container, CCP4File
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        #print 'CTaskFrame.openTask',taskName,jobId,'followJobId',followJobId,'cloneJobId',cloneJobId
        # If there is jobid try to get the paramsFile and ensure consistent taskName
        # ??? Should we be concerned about version number ???
        time1 = time.time()
        # Check we have a clone params file before creating new job
        cloneParamsFile = None
        if cloneJobId is not None:
            cloneParamsFile = PROJECTSMANAGER().makeFileName(jobId=cloneJobId, mode='JOB_INPUT')
            if not os.path.exists(cloneParamsFile):
                cloneParamsFile = PROJECTSMANAGER().makeFileName(jobId=cloneJobId, mode='PARAMS')
            if not os.path.exists(cloneParamsFile):
                QtWidgets.QMessageBox.warning(self,self.windowTitle(),'No parameter file found for task to clone')
                return None
        #print 'openTask',cloneJobId,cloneParamsFile
        paramsFile = None
        # If we are opening a pre-existing job then just ensure can read the params.def.xml file
        # or create a new job in the database and with a job directory
        if jobId is not None:
            paramsFile = PROJECTSMANAGER().makeFileName(jobId = jobId,mode='JOB_INPUT')
            # We could be loading a sub-job without a input params file
            if not os.path.exists(paramsFile): paramsFile = PROJECTSMANAGER().makeFileName(jobId = jobId,mode='PARAMS')
            #print 'CProjectViewer.openTask paramsFile',paramsFile,os.path.exists(paramsFile)
            if os.path.exists(paramsFile):
                header = CCP4File.xmlFileHeader(paramsFile)
                if taskName is None:
                    taskName = str(header.pluginName)
                elif taskName != str(header.pluginName):
                    err = CException(self.__class__,102,'Suggested: '+str(taskName)+' File: '+str(paramsFile)+' Contains: '+str(header.pluginName))
                    warningMessage(err, 'Error loading params file','Error loading params file for job id: '+jobId,parent=self)
                    return None
            else:
                paramsFile = None
            ifNewJob= False
        else:
            jobId,pName,jNumber = PROJECTSMANAGER().newJob(taskName=taskName,projectId=self.openJob.projectId)      
            ifNewJob= True
        # Create an COpenJob instance to hold the meta-data for this job
        openJob=COpenJob(jobId=jobId,projectId=self.openJob.projectId)
        taskEditable =  ( openJob.status in ['Unknown','Pending'] )
        from ..core.CCP4Preferences import PREFERENCES
        if openJob.status == "Finished" and PREFERENCES().AUTO_UPDATE_REPORT80:
            stamp80 = os.path.join(PROJECTSMANAGER().makeFileName(jobId = openJob.jobId,mode='ROOT'),"REP8STMP")
            reportFile = PROJECTSMANAGER().makeFileName(jobId = openJob.jobId,mode='REPORT')
            redoReportFor80 = False
            if not os.path.exists(stamp80):
                #This report has never been created in 8.0
                redoReportFor80 = True
                print("########################################")
                print(stamp80,"does not exist. I am remaking report file.")
                print("########################################")
            elif os.path.exists(stamp80) and os.path.exists(reportFile):
                reportTime = os.path.getmtime(reportFile)
                stampTime = os.path.getmtime(stamp80)
                if reportTime > stampTime:
                    redoReportFor80 = True
                    print("########################################")
                    print(stamp80,"is older than report file. I am remaking the report file.")
                    print("########################################")
            if redoReportFor80:
                self.report71.emit(openJob)
                Path(stamp80).touch()
                print("Created stamp",stamp80)
        # For a cloned job copy the params.def.xml
        if cloneJobId is not None:
            print("%%%%%%%%%%%%%%%%%%%% COPYING")
            if paramsFile is None:
                paramsFile = PROJECTSMANAGER().makeFileName(jobId = openJob.jobId,mode='JOB_INPUT')
            try:
                CCP4File.cloneI2XmlFile(cloneParamsFile,paramsFile,{ 'jobId' : str(openJob.jobnumber) }, taskFrame=self, taskName=taskName, suggestedParams=suggestedParams )
                openJob.clonedFromJobId = cloneJobId
                #print 'CTaskFrame.openTask cloneFromJobId',cloneJobId, openJob.clonedFromJobId
            except:
                print('ERROR cloning params file',cloneParamsFile)
                paramsFile = None
            else:
                splitCloneParamsFile = os.path.split(cloneParamsFile)      
                subJobParamsFiles = glob.glob(os.path.join(splitCloneParamsFile[0],'job_*_'+splitCloneParamsFile[1]))
                for subFile in subJobParamsFiles:
                    CCP4File.cloneI2XmlFile(subFile,os.path.join(os.path.split(paramsFile)[0],os.path.split(subFile)[1]),{ 'jobId' : str(openJob.jobnumber) }, taskFrame=self, taskName=taskName,suggestedParams=suggestedParams )
        # Find the task def file and create a CContainer with data contents based on the def file
        defFile = TASKMANAGER().lookupDefFile(openJob.taskname,openJob.taskversion)
        if defFile is None:
            print('Failed to find def file in openTask')
            print('CTaskFrame.openTask defFile',openJob.taskname,type(openJob.taskname),defFile,openJob.taskversion)
            return self.openJob
        # Set up data container
        container = CCP4Container.CContainer(parent=self,definitionFile=defFile,guiAdmin=True)

        if patchParamsFile is not None:
            # expect patch params file to be passed thru from a whatNext command
            if patchParamsFile.startswith('$CCP4I2'):
                patchParamsFile = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),patchParamsFile[8:]))
            elif patchParamsFile.startswith('$PROJECT'):
                patchParamsFile = os.path.normpath(os.path.join(openJob.projectDir,patchParamsFile[9:]))
            try:
                container.loadDataFromXml(patchParamsFile)
            except:
                print('ERROR loading patch params file',patchParamsFile)
        # If it is existing or cloned job then paramsFile exists and is loaded to container
        if paramsFile is not None:
            print('CProjectViewer.openTask loading',paramsFile)
            try:
                container.loadDataFromXml(paramsFile)
            except:
                print('ERROR loading params file',paramsFile)
        if cloneParamsFile is not None:
            print('CProjectViewer.openTask loading cloned file',cloneParamsFile)
            try:
                #MAJOR CHECKME SJM why on earth would I load clonedParamsFilethis if I already used paramsFile above?
                if cloneParamsFile is not None and paramsFile is None:
                    container.loadDataFromXml(cloneParamsFile)
                if cloneJobId is not None and container.guiAdmin.get('jobTitle') is not None:
                    jobTitle = PROJECTSMANAGER().db().getJobInfo(jobId=cloneJobId,mode='jobtitle')
                    if jobTitle is not None:
                        container.guiAdmin.jobTitle.set(jobTitle)
            except:
                print('ERROR loading cloned params file',cloneParamsFile)
            try:
                if container.guiAdmin.jobTitle.isSet():
                    PROJECTSMANAGER().db().updateJob(jobId=jobId,key='jobTitle',value=container.guiAdmin.jobTitle.__str__())
            except:
                print('ERROR saving cloned jobTitle')
        # Automatically set output filenames in the container (to overwrite those from cloned job) and save to params.xml
        PROJECTSMANAGER().setOutputFileNames(container=container,projectId=openJob.projectId,
                                             jobNumber=openJob.jobnumber,force=ifNewJob or (openJob.clonedFromJobId is not None ))
        self.saveStatus(openJob=openJob)
        # Make a guess here if there is no followJobId and it is a new job
        # and its not a clone with params already set!
        if openJob.clonedFromJobId is not None:
            followJobId = None
        elif ifNewJob and followJobId is None:
            followJobId = PROJECTSMANAGER().db().getProjectFollowFromJobId(projectId=openJob.projectId)
        # Draw the task input widget 
        try:
            t2 = time.time()
            if taskEditable or self.tab.currentIndex() == self.INPUT_TAB:
                taskWidget = self.inputFrame.createTaskWidget(openJob.taskname,projectId=openJob.projectid,jobId=openJob.jobId,
                                                              container=container,taskEditable=taskEditable,followJobId=followJobId)
                if hasattr(taskWidget,"launchJobRequestSignal"):
                    taskWidget.launchJobRequestSignal.connect(self.launchJobRequest)
            else:
                self.inputFrame.closeTaskWidget()
        except CException as e:
            warningMessage(e, self.windowTitle(),'Error drawing task widget',parent=self)
        except Exception as e:
            err = CErrorReport(self.__class__,999,details=str(e),exc_info=sys.exc_info())
            warningMessage(err, self.windowTitle(),'Unknown error drawing task widget\n\n'+str(e),parent=self)
        else:
            t3 = time.time()
            self.openJob = openJob
            # Update the task title bar and the 'Next' buttons
            self.titleBar.setOpenJob(self.openJob)
            self.updateTaskFrame()
            print('opentask times total',time.time()-time1)
            print('opentask times drawing',t3-t2)
            return self.openJob
        # Only get here if failed to create taskWidget properly so delete it
        try:
            taskWidget.deleteLater()
        except:
            pass
        return self.openJob

    def projectId(self):
        return self.openJob.projectId

    @QtCore.Slot(str,dict)
    def launchJobRequest(self,taskName,args):
        self.launchJobRequestSignal.emit(taskName,args)


class CTaskMainWindow(CCP4WebBrowser.CMainWindow):

    windowAboutToClose = QtCore.Signal()

    '''Popout window for task input and report'''

    def openApplication(self,application):
        '''Pass on request for program logs'''
        from .CCP4WebBrowser import WEBBROWSER
        WEBBROWSER().openApplication(application)

    def handleProjectMenuExport(self):
        pass

    def openSendReport(self):
        '''Open window to send developer error report'''
        from . import CCP4ErrorReportViewer
        widget = CCP4ErrorReportViewer.CSendJobError(self, projectId=self.taskFrame.openJob.projectId, projectName=self.taskFrame.openJob.projectName)
        widget.show()

    def showHelp(self,mode='ccp4i2',newTab=True):
        pass

    def openManageImportFiles(self):
        pass

    def isProjectViewer(self):
        return False

    def widgetIsSearchable(self):
        return False

    def widgetIsRunable(self):
        return False

    def handleSave(self):
        pass

    def widgetIsSaveable(self):
        return False

    def handlePrint(self):
        pass

    def widgetIsPrintable(self):
        return False

    def handleRun(self):
        pass

    def openFind(self):
        pass

    def isFindFrameOpen(self):
        return False

    def deleteTab(self):
        pass

    def historyBack(self):
        pass

    def historyForward(self):
        pass

    def reloadPage(self):
        pass

    def resetZoom(self):
        if hasattr(self,"widget") and hasattr(self.widget(),"webView") and hasattr(self.widget().webView,"setZoomFactor"):
            self.widget().webView.setZoomFactor(1.0)

    def zoomIn(self):
        if hasattr(self,"widget") and hasattr(self.widget(),"webView") and hasattr(self.widget().webView,"setZoomFactor"):
            self.widget().webView.setZoomFactor(self.widget().webView.zoomFactor()*1.2)

    def zoomOut(self):
        if hasattr(self,"widget") and hasattr(self.widget(),"webView") and hasattr(self.widget().webView,"setZoomFactor"):
            self.widget().webView.setZoomFactor(self.widget().webView.zoomFactor()/1.2)

    def __init__(self,parent,projectName,jobId,version=''):
        #QtWidgets.QMainWindow.__init__(self,parent)
        CCP4WebBrowser.CMainWindow.__init__(self, parent)
        self.version=version
        self.setWindowTitle(self.version+'job from project: '+projectName)
        self.setWindowIcon(CCP4WebBrowser.mainWindowIcon())
        frame = QtWidgets.QFrame(self)
        frame.setLayout( QtWidgets.QVBoxLayout())
        self.setObjectName('job'+str(jobId))
        self.titleBar = CTaskTitleBar(self)
        frame.layout().addWidget(self.titleBar)
        self.buttons = CTaskButtons(self,frame.layout(),mode=CTaskButtons.RUNONLYMODE)
        self.setCentralWidget(frame)
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        PROJECTSMANAGER().db().jobStatusUpdated.connect(self.titleBar.setStatusBar)
        PROJECTSMANAGER().db().jobFinished.connect(self.titleBar.setStatusBar)

    def widget(self):
        return self.centralWidget().layout().itemAt(1).widget()

    def closeEvent(self,event):
        self.deleteLater()
        self.windowAboutToClose.emit()
        event.accept()

    def getTaskName(self):
        return None

    @QtCore.Slot()
    def runTask(self):
        self.widget().runTask('Local')


class CDeleteJobGui(QtWidgets.QDialog):

    jobsDeleted = QtCore.Signal()

    '''Show knock-on effect of deleting a job and get user confirmation to delete all'''

    def __init__(self,parent=None,projectId=None,jobIdList=None,jobTreeList=[],deleteImportFiles=False,jobsToDeleteWithSelectedFiles=[],label=None,ifXtrJobs=False):
        QtWidgets.QDialog.__init__(self,parent)
        self.setWindowTitle('Delete jobs')
        self.projectId = projectId
        self.jobIdList = jobIdList
        self.jobTreeList = jobTreeList
        self.importFileList = []
        self.deleteImportFiles = deleteImportFiles
        self.jobsToDeleteWithSelectedFiles = jobsToDeleteWithSelectedFiles
        self.setModal(True)
        self.setLayout(QtWidgets.QVBoxLayout())
        if len(jobsToDeleteWithSelectedFiles)>0:
            line = QtWidgets.QHBoxLayout()
            lab = QtWidgets.QLabel( self)
            lab.setPixmap( QtGui.QPixmap(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','list_delete.png')  ).scaled(16,16) )
            line.addWidget(lab)
            lab = QtWidgets.QLabel('marked jobs have files used in the current Job Input',self)
            lab.setObjectName('emphasise')
            line.addWidget(lab)
            line.addStretch(5)
            self.layout().addLayout(line)
        if ifXtrJobs:
            lab = QtWidgets.QLabel('Unselected jobs highlighted in pink',self)
            lab.setObjectName('emphasise')
            self.layout().addWidget(lab)
        nImportedFiles = 0
        for jobTree in self.jobTreeList:
            nImportedFiles += len(jobTree[1])
        if nImportedFiles >0:
            self.deleteImportWidget = QtWidgets.QCheckBox('Delete files imported by subsequent jobs - this may imply deleting additional jobs',self)
            if deleteImportFiles:
                self.deleteImportWidget.setCheckState(QtCore.Qt.Checked)
            self.layout().addWidget(self.deleteImportWidget)
            self.deleteImportWidget.stateChanged.connect(self.handleDeleteImportChanged)
        '''
        if len(importFileList)>0:
          for importId,fileName in importFileList:
            self.layout().addWidget(QtWidgets.QLabel(fileName,self))
            self.importFileList.append([importId,fileName])
        '''
        self.tree = CJobTree(self)
        self.layout().addWidget(self.tree)
        self.tree.clear()
        for jobTree in self.jobTreeList:
            self.tree.load(jobTree,jobsToDeleteWithSelectedFiles=self.jobsToDeleteWithSelectedFiles)
        buttonBox = QtWidgets.QDialogButtonBox(self)
        but = buttonBox.addButton('Delete all these jobs',QtWidgets.QDialogButtonBox.ApplyRole)
        but.setAutoDefault(0)
        but.clicked.connect(self.deleteJobs)
        but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
        but.setAutoDefault(0)
        but.clicked.connect(self.close)
        self.layout().addWidget(buttonBox)

    @QtCore.Slot(bool)
    def handleDeleteImportChanged(self,deleteImportFiles):
        self.deleteImportFiles = deleteImportFiles
#FIXME There is no self.jobId, what is this meant to do and what is self.jobId meant to be?
#Possibly getMultiFollowOnJobs(jobIdList=self.jobIdList,...) ?
#And why is this self.followOnJobs? Why member of class?
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        self.followOnJobs = PROJECTSMANAGER().db().getFollowOnJobs(jobId=self.jobId,traceImportFiles=(deleteImportFiles>0))
        self.tree.load(self.followOnJobs,jobsToDeleteWithSelectedFiles=self.jobsToDeleteWithSelectedFiles)
        
    @QtCore.Slot()
    def deleteJobs(self):
        '''
        # Beware importFilePath() returns a name for a new import so this will not work
        for importId,fileName in self.importFileList:
          filePath = PROJECTSMANAGER().importFilePath(projectId=self.projectId,baseName=fileName)
          try:
            os.remove(filePath)
            PROJECTSMANAGER().db().deleteImportFile(importId=importId)
          except:
            print 'ERROR deleting file',filePath
        '''
        for jobTree in self.jobTreeList:
            self.deleteJobs0(jobTree,deleteImportFiles=self.deleteImportFiles)
        self.jobsDeleted.emit()
        self.close()

    def deleteJobs0(self,jobTree=None,deleteImportFiles=False):
        jobId,importFiles,descendents = jobTree
        if jobId is not None:
            from ..core.CCP4ProjectsManager import PROJECTSMANAGER
            PROJECTSMANAGER().deleteJob(jobId=jobId,importFiles=importFiles,projectId=self.projectId,deleteImportFiles=deleteImportFiles)
        for childJobTree in descendents:
            self.deleteJobs0(childJobTree,deleteImportFiles=False)


class CJobTree(QtWidgets.QTreeWidget):
    ''' Sub-widget of CDeleteJobGui window'''
    def __init__(self,parent):
        QtWidgets.QTreeWidget.__init__(self,parent)
        self.setColumnCount(3)
        self.setHeaderLabels(['Job number', 'Task name', 'Status'])
        self.setColumnWidth(0, 100)
        self.setColumnWidth(1, 300)
        self.icon = QtGui.QIcon( QtGui.QPixmap(os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', 'list_delete.png')).scaled(16,16))

    def load(self, jobTree, treeParentId=None, jobsToDeleteWithSelectedFiles=[]):
        jobId, importFiles, childJobTree = jobTree
        from ..core.CCP4ProjectsManager import PROJECTSMANAGER
        jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber', 'taskname', 'status'])
        taskTitle = TASKMANAGER().getTitle(jobInfo['taskname'])
        # Only add job once - beware a jobs could be child of multiple preceeding jobs so appear in jobTree more than once
        if len(self.findItems(jobInfo['jobnumber'],QtCore.Qt.MatchExactly)) == 0:
            item = QtWidgets.QTreeWidgetItem([jobInfo['jobnumber'], taskTitle,str(jobInfo['status'])])
            if treeParentId is None:
                self.addTopLevelItem(item)
            else:
                #treeParentId.addChild(item)
                self.addTopLevelItem(item)
                for col in 0,1,2:
                    item.setBackground(col,QtGui.QBrush(QtGui.QColor('pink')))
            if jobId in jobsToDeleteWithSelectedFiles:
                item.setIcon(0,self.icon)
            for childJob in childJobTree:
                self.load(childJob,treeParentId=item,jobsToDeleteWithSelectedFiles=jobsToDeleteWithSelectedFiles)
