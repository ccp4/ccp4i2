from __future__ import print_function
from PySide2 import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from core.CCP4Modules import PROJECTSMANAGER
from core.CCP4ErrorHandling import *
from . import CCP4i2DjangoSession
import functools
import sys, os

#-------------------------------------------------------------------
class SyncToDjango(CTaskWidget):
#-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'SyncToDjango'
    TASKVERSION = 0.1
    TASKMODULE=[ 'developer_tools' ]
    TASKTITLE='Exchange projects with Django Back end'
    EDITOR = True  
    DESCRIPTION='''Synchronization utility for CCP4i2 projects in locations which have a django-based project archive set up'''
    RANK = 2
    MAINTAINER = 'martin.noble@ncl.ac.uk'
    DIRECTION="NeitherUpNorDown"
    
    def drawContents(self):
        from lxml import etree
        self.openFolder(folderFunction='inputData',followFrom=False)
        from core import CCP4Modules
        container = CCP4Modules.SERVERSETUP()
        container.load()
        serversEtree = container.getEtree(excludeUnset=True)
        allServerNames = serversEtree.xpath("//CHostname/text()")
        archiveServerNames = [name for name in allServerNames if "ManageCCP4i2Archive" in name]
        
        self.createLine(['label','Select server','widget','serverName'])
        self.container.inputData.serverName.setQualifiers({'enumerators':archiveServerNames,'menuText':archiveServerNames})
        self.container.inputData.serverName.set(archiveServerNames[0])
        try:
            self.getWidget('serverName').populateComboBox(self.container.inputData.serverName)
            self.getWidget('serverName').updateViewFromModel()
        except:
            pass

        subFrame = self.openSubFrame(frame=True)
        hbox = QtWidgets.QHBoxLayout()
        label = QtWidgets.QLabel("Username")
        label.setFixedWidth(200)
        hbox.addWidget(label)
        self.usernameField = QtWidgets.QLineEdit()
        hbox.addWidget(self.usernameField)
        subFrame.layout().addLayout(hbox)
        
        hbox = QtWidgets.QHBoxLayout()
        label = QtWidgets.QLabel("Password")
        label.setFixedWidth(200)
        hbox.addWidget(label)
        self.passwordField = QtWidgets.QLineEdit()
        self.passwordField.setEchoMode(QtWidgets.QLineEdit.Password)
        hbox.addWidget(self.passwordField)
        subFrame.layout().addLayout(hbox)

        hbox = QtWidgets.QHBoxLayout()
        listRemote = QtWidgets.QPushButton("List remote")
        listRemote.clicked.connect(self.onListRemoteClicked)
        hbox.addWidget(listRemote)
        hbox.addStretch(1)
        listLocal = QtWidgets.QPushButton("List local")
        listLocal.clicked.connect(self.onListLocalClicked)
        hbox.addWidget(listLocal)
        subFrame.layout().addLayout(hbox)
        self.projectsTree = QtWidgets.QTreeView()
        self.projectsTree.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection);
        subFrame.layout().addWidget(self.projectsTree)
        
        hbox = QtWidgets.QHBoxLayout()
        filterLabel = QtWidgets.QLabel("Filter by text:")
        filterLabel.setFixedWidth(200)
        hbox.addWidget(filterLabel)
        self.filterText = QtWidgets.QLineEdit()
        self.filterText.textChanged.connect(self.onFilterTextChanged)
        hbox.addWidget(self.filterText)
        subFrame.layout().addLayout(hbox)
        
        syncSelectedButton = QtWidgets.QPushButton("Sync selected")
        syncSelectedButton.clicked.connect(self.syncSelected)
        subFrame.layout().addWidget(syncSelectedButton)
        
        openSelectedButton = QtWidgets.QPushButton("Open selected")
        openSelectedButton.clicked.connect(self.openSelected)
        subFrame.layout().addWidget(openSelectedButton)

        cootSelectedButton = QtWidgets.QPushButton("Coot selected")
        cootSelectedButton.clicked.connect(self.cootSelected)
        subFrame.layout().addWidget(cootSelectedButton)

    def ccp4i2DjangoSession(self):
        if not hasattr(self,"_ccp4i2DjangoSession"):
            serverName = self.container.inputData.serverName.__str__()
            transportType = "http://"
            if "443" in serverName: transportType = "https://"
            self._ccp4i2DjangoSession = CCP4i2DjangoSession.CCP4i2DjangoSession(transportType+serverName,self.usernameField.text(), self.passwordField.text(), projectsManager=PROJECTSMANAGER())
        return self._ccp4i2DjangoSession

    @QtCore.Slot(str)
    def onListRemoteClicked(self, item):
        serverName = self.container.inputData.serverName.__str__()
        transportType = "http://"
        if "443" in serverName: transportType = "https://"
        response = self.ccp4i2DjangoSession().getURLWithValues(transportType+serverName+"?listProjects",{})
        responseText = response.read()
        import json
        self.directoryList = json.loads(responseText)
        projectHierarchy = ProjectHierarchy(directoryList = self.directoryList)
        self.projectsTree.setModel(projectHierarchy)
        self.DIRECTION = "downloaded from"

    @QtCore.Slot(str)
    def onListLocalClicked(self, item):
        self.directoryList = PROJECTSMANAGER().db().listProjects()
        projectHierarchy = ProjectHierarchy(directoryList = self.directoryList)
        self.projectsTree.setModel(projectHierarchy)
        self.DIRECTION = "uploaded to"
    
    @QtCore.Slot(str)
    def onFilterTextChanged(self, item):
        if hasattr(self,"directoryList"):
            itemAsString = str(item)
            if len(itemAsString) == 0:
                projectHierarchy = ProjectHierarchy(directoryList = self.directoryList)
                self.projectsTree.setModel(projectHierarchy)
            else:
                catenatedSelfAndDescendents = {}
                for project in self.directoryList:
                    catenatedSelfAndDescendents[project[0]] = self.catenatedSelfAndDescendents(project)
                filteredList = [project for project in self.directoryList if itemAsString in catenatedSelfAndDescendents[project[0]]]
                projectHierarchy = ProjectHierarchy(directoryList = filteredList)
                self.projectsTree.setModel(projectHierarchy)

    def catenatedSelfAndDescendents(self, project):
        result = project[1]
        descendents = [subProject for subProject in self.directoryList if subProject[3] == project[0]]
        for descendent in descendents:
            result += ("_"+self.catenatedSelfAndDescendents(descendent))
        return result
    
    @QtCore.Slot(str)
    def openSelected(self, item):
        from qtgui.CCP4ProjectManagerGui import openProject
        for anIndex in self.projectsTree.selectionModel().selection().indexes():
            openProject(anIndex.internalPointer().ref.projectTuple[0])

    @QtCore.Slot(str)
    def cootSelected(self, item):
        for anIndex in self.projectsTree.selectionModel().selection().indexes():
            projectId = anIndex.internalPointer().ref.projectTuple[0]
            from core.CCP4Modules import PREFERENCES
            archiveName = self.ccp4i2DjangoSession().lastLigandPipelineExportForProjectId(projectId)
            import subprocess
            currentDirectory = os.path.realpath(PROJECTSMANAGER().db().jobDirectory(self.jobId()))
            import zipfile
            zip_ref = zipfile.ZipFile(archiveName, 'r')
            zip_ref.extractall(currentDirectory)
            filenames = [os.path.basename(filePath) for filePath in zip_ref.namelist()]
            dirnames = [os.path.dirname(filePath) for filePath in zip_ref.namelist()]
            zip_ref.close()
            exePath = 'coot'
            if sys.platform == 'win32':
                exePath = str(PREFERENCES().COOT_EXECUTABLE)
            setupCootPath = os.path.join(os.path.dirname(__file__),"SetupCoot.py")
            cootArgs = [exePath,'--python','--no-state-script','--script',setupCootPath,'XYZOUT.pdb','hklout.mtz']
            if "DICTOUT.cif" in filenames: cootArgs += ["--dict","DICTOUT.cif"]
            subprocess.call(cootArgs, cwd=os.path.join(currentDirectory, dirnames[0]))

    @QtCore.Slot(str)
    def syncSelected(self, item):
        print("Sync selected", self.DIRECTION)
        if self.DIRECTION == "NeitherUpNorDown":
            pass
        else:
            messageText = "Confirm that following projects should be synced with server "+self.container.inputData.serverName.__str__()+":\n"
            projectsToSync = []
            for anIndex in self.projectsTree.selectionModel().selection().indexes():
                #projectName = self.projectsTree.model().data(anIndex, QtCore.Qt.DisplayRole)
                messageText += ("\t"+anIndex.internalPointer().ref.projectTuple[1] + "\n")
                projectsToSync.append(anIndex.internalPointer().ref.projectTuple)
            button = QtWidgets.QMessageBox.question(self, "Sync selected", messageText,
                                                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                                                
            if button == QtWidgets.QMessageBox.Yes:
                print("Okay...I'll do it")
                for project in projectsToSync:
                    self.ccp4i2DjangoSession().syncProject(project[0])
                return


from .TreeModel import TreeNode, TreeModel

class NamedElement(object): # your internal structure
    def __init__(self, name, subelements, projectTuple):
        self.name = name
        self.subelements = subelements
        self.projectTuple = projectTuple

class NamedNode(TreeNode):
    def __init__(self, ref, parent, row):
        self.ref = ref
        TreeNode.__init__(self, parent, row)

    def _getChildren(self):
        return [NamedNode(elem, self, index)
            for index, elem in enumerate(self.ref.subelements)]

class ProjectHierarchy(TreeModel):
    def __init__(self, directoryList):
        self.directoryList = directoryList
        self.elements = []
        for project in self.directoryList:
            self.elements.append(NamedElement(project[1], None, project))
        for element in self.elements:
            subelements = [subelement for subelement in self.elements if subelement.projectTuple[3] == element.projectTuple[0]]
            element.subelements = subelements
        self.rootElements = [element for element in self.elements if element.projectTuple[3] is None]
        TreeModel.__init__(self)

    def _getRootNodes(self):
        return [NamedNode(elem, None, index)
            for index, elem in enumerate(self.rootElements)]

    def columnCount(self, parent):
        return 1

    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == QtCore.Qt.DisplayRole and index.column() == 0:
            return node.ref.name
        return None

    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole \
            and section == 0:
            return 'Name'
        return None
    
#--------------------------------------------------------------------------------------------------
    def flags(self, index=QtCore.QModelIndex()):
#--------------------------------------------------------------------------------------------------
        #print "Evaluating flags for", str(index.column()), str(index.row())
        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled




