"""
Copyright (C) 2010 University of York
Liz Potterton Feb 2010 - Create CCP4ProjectManager mostly as placeholder
Liz Potterton May 2010 - Rename CCP4ProjectManagerGui 
"""

##@package CCP4ProjectManagerGui (QtGui) Project handling tools for web browser (currently just a Project menu)

import copy
import datetime
import functools
import glob
import os
import shutil
import sys
import time

from PySide2 import QtCore, QtGui, QtWidgets

from . import CCP4Widgets
from .. import I2_TOP
from ..core import CCP4Utils
from ..core.CCP4Config import CONFIG
from ..core.CCP4ErrorHandling import CException, Severity
from ..core.CCP4Modules import MIMETYPESHANDLER
from ..core.CCP4Modules import PROJECTSMANAGER
from ..core.CCP4Modules import WEBBROWSER
from ..core.CCP4WarningMessage import warningMessage


DIAGNOSTIC = True

def openProject(projectId=None,projectName=None):
    print('openProject',projectId,projectName, end=' ')
    if projectId is None: projectId = PROJECTSMANAGER().db().getProjectId(projectName=projectName)
    if projectId is None: return None
    from . import CCP4ProjectViewer
    for window in CCP4ProjectViewer.CProjectViewer.Instances:
      if  window.projectId() == projectId:
        window.show()
        window.raise_()
        return window

    p = CCP4ProjectViewer.CProjectViewer(projectId=projectId)
    p.show()
    p.raise_()
    return p


def formatDate(fTime):
  from . import CCP4ProjectWidget
  if fTime is None: return ''
  try:
    date = time.strftime(CCP4ProjectWidget.CTreeItemJob.DATE_FORMAT,time.localtime(int(fTime)))
    if date == CCP4ProjectWidget.CTreeItemJob.TODAY:
      date = time.strftime(CCP4ProjectWidget.CTreeItemJob.TIME_FORMAT,time.localtime(int(fTime)))
    elif date.split(' ')[-1] == CCP4ProjectWidget.CTreeItemJob.THISYEAR:
      date=date.rsplit(' ',1)[0]
    else:
      date=date.split(' ',1)[1]
    return  date
  except:
    return ''

class ProjectsViewDelegate(QtWidgets.QStyledItemDelegate):
    def paint(self, painter, option, index):
        options = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(options,index)

        style = QtWidgets.QApplication.style() if options.widget is None else options.widget.style()

        if index.column() != 0:
            options.icon = QtGui.QIcon()

        if index.column() == 2 or index.column() == 3:
            try:
                dt = datetime.datetime.utcfromtimestamp(float(options.text))
            except ValueError as err:
                dt = datetime.datetime(2011, 11, 4, 0, 0)
            options.text = dt.strftime("%d %b %y")
            style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, options, painter)
        else:
            style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, options, painter)


class CNewProjectGui(QtWidgets.QDialog):

  projectCreated = QtCore.Signal(str)

  insts = None

  def __init__(self,parent=None):
    QtWidgets.QDialog.__init__(self,parent)
    # This makes this window stay on top even when we want the file browser to be on top
    #self.setModal(True)
    self.setWindowTitle('Create a New Project')
    self.setLayout(QtWidgets.QVBoxLayout())

    from ..core import CCP4Data
    from ..core import CCP4DataManager
    from ..core import CCP4File
    self.name = CCP4Data.CString(parent=self)
    self.directory = CCP4File.CDataFile(parent=self,qualifiers={'isDirectory' : True, 'mustExist': False, 'fromPreviousJob' : False  })


    MARGIN = 1

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('<b>Name of project/folder</b><br/><i>(alphanumeric, _ or -)</i>',self))
    self.nameWidget = CCP4DataManager.DATAMANAGER().widget(model=self.name,parentWidget=self)
    self.nameWidget.setToolTip('Project name - may be used by programs that only allow\none word of alphanumeric characters, dash or underscore.')
    line.addWidget(self.nameWidget)
    self.nameWidget.widget.setValidator(QtGui.QRegExpValidator(QtCore.QRegExp("[A-Za-z0-9_-]{0,255}"), self ));
    self.layout().addLayout(line)

    for textLine in [   "By default all projects go in the 'CCP4I2_PROJECTS' directory in your",
                     "home area - click 'Select directory' to choose an alternative.",
                     "Hint to organise your projects: in the 'Manage projects' window",
                     "you can use a project as a folder and drag other projects into it" ]:

      line = QtWidgets.QHBoxLayout()
      line.addWidget(QtWidgets.QLabel(textLine,self))
      self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    self.descriptionWidget = CProjectDescription(self,projectId=None,projectName='New project')
    line.addWidget(self.descriptionWidget)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    buttonBox = QtWidgets.QDialogButtonBox(self)

    self.createButton = buttonBox.addButton('Create project',QtWidgets.QDialogButtonBox.AcceptRole)
    self.createButton.setFocusPolicy(QtCore.Qt.NoFocus)
    self.createButton.setDefault(True)
    self.createButton.clicked.connect(functools.partial(self.addProject,None))
    selectDirButton = buttonBox.addButton('Select directory',QtWidgets.QDialogButtonBox.ActionRole)
    selectDirButton.setFocusPolicy(QtCore.Qt.NoFocus)
    selectDirButton.clicked.connect(functools.partial(self.createButton.setEnabled,False))
    selectDirButton.clicked.connect(self.selectDir)
    cancelButton = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    cancelButton.setFocusPolicy(QtCore.Qt.NoFocus)
    cancelButton.clicked.connect(functools.partial(self.createButton.setEnabled,True))
    cancelButton.clicked.connect(self.close)
    helpButton = buttonBox.addButton(QtWidgets.QDialogButtonBox.Help)
    helpButton.setFocusPolicy(QtCore.Qt.NoFocus)
    helpButton.clicked.connect(self.help)

    line.addStretch(0.5)
    line.addWidget(buttonBox)
    line.addStretch(0.5)
    self.layout().addLayout(line)

    self.directory.dataChanged.connect(self.setProjectName)
    self.rejected.connect(self.close)


  @QtCore.Slot()
  def setProjectName(self):
    if not self.directory.isSet() or self.name.isSet(): return
    tail = os.path.split(self.directory.__str__())[1]
    #print 'setProjectName',tail
    self.name.set(tail)
    self.nameWidget.updateViewFromModel()

  @QtCore.Slot()
  def help(self):
    WEBBROWSER().loadWebPage(helpFileName='general/tutorial' , target='projects')

  def clear(self):
    self.name.unSet()
    self.directory.unSet()
    #self.directoryWidget.updateViewFromModel()
    self.nameWidget.updateViewFromModel()

  @QtCore.Slot()
  def selectDir(self):
    self.nameWidget.updateModelFromView()
    name = self.name.get()
    if name is None:
      self.createButton.setEnabled(True)
      QtWidgets.QMessageBox.warning(self,'No project name','You must provide a project name')
      return
    elif PROJECTSMANAGER().projectNameStatus(name) != 0:
      self.createButton.setEnabled(True)
      QtWidgets.QMessageBox.warning(self,'Project exists','A project of that name already exists')
      return

    rv = QtWidgets.QFileDialog.getExistingDirectory(caption='Select directory for saving project '+str(name))
    #print 'selectDir getExistingDirectory',rv,len(rv)
    if rv is not None and len(rv)>0:
        self.addProject(str(rv))
    self.createButton.setEnabled(True)

  @QtCore.Slot(str)
  def addProject(self,directory = None):
    #print 'CNewProjectGui.addProject',directory,'*'
    self.nameWidget.updateModelFromView()
    name = self.name.get()
    #print 'CNewProjectGui.addProject',name 
    if name is None:
      QtWidgets.QMessageBox.warning(self,'No project name','You must provide a project name')
      return
    elif PROJECTSMANAGER().projectNameStatus(name) != 0:
      QtWidgets.QMessageBox.warning(self,'Project exists','A project of that name already exists')
      return

    name0 = CCP4Utils.safeOneWord(name)
    if name0 != name:
      print('Fixing name to ',name0)

    if directory is None:
      directory = CCP4Utils.getProjectDirectory()
      if directory is None:
        QtWidgets.QMessageBox.warning(self,'CCP4I2_PROJECTS directory','Failed creating a CCP4I2_PROJECTS directory in your home directory')
        return
      directory = os.path.join(directory,name0)
    else:
      directory = os.path.normpath(directory)
      if os.path.isfile(directory):
        QtWidgets.QMessageBox.warning(self,'Create project','Selected directory is a file')
        return
      elif os.path.samefile(directory,CCP4Utils.getHOME()):
        QtWidgets.QMessageBox.warning(self,'Create project','Selected directory is users home area')
        return
      elif os.path.samefile(directory,CCP4Utils.getProjectDirectory()):
        QtWidgets.QMessageBox.warning(self,'Create project','Selected directory is the CCP4 master project directory')
        return



    altName = PROJECTSMANAGER().aliasForDirectory(directory)
    if altName is not None:
         QtWidgets.QMessageBox.warning(self,'Directory already a project','This directory is already used by the project: '+altName)
         return

    parentProject,relpath,parentProjectId =  PROJECTSMANAGER().interpretDirectory(directory)
    if parentProject is not None:
      rv = QtWidgets.QMessageBox.question(self,'Make sub-project',
          'This is sub-directory of another project \nMake this project a sub-project of '+parentProject,
          QtWidgets.QMessageBox.Ok|QtWidgets.QMessageBox.Cancel)
      if rv == QtWidgets.QMessageBox.Cancel: return

    if not os.path.exists(directory):
      parentDir = os.path.split(directory)[0]
      if not os.path.exists(parentDir):
        QtWidgets.QMessageBox.warning(self,'Directory does not exist','The parent directory for the project directory does not exist: ' +parentDir )
        return

      try:
        os.mkdir(directory)
      except:
         QtWidgets.QMessageBox.warning(self,'Error making directory','Error making directory: '+directory)
         return

    try:
      projectId = PROJECTSMANAGER().createProject(projectName=name0,projectPath=directory)
    except CException as e:
      warningMessage(e, parent=self)
      return
    except:
      QtWidgets.QMessageBox.warning(self,'Error creating project','Unknown error creating project')
      return

    self.descriptionWidget.save(projectId)

    self.projectCreated.emit(projectId)
    PROJECTSMANAGER().backupDBXML()
    pView = openProject(projectId)

    #print 'CNewProjectGui.addProject to close'
    self.close()

  def close(self):
    for child in self.findChildren(QtWidgets.QDialog):
      child.close()
    QtWidgets.QDialog.close(self)

class CExportOptionsDialog(QtWidgets.QDialog):

  exportSignal = QtCore.Signal(int)

  def __init__(self,parent=None,projectId=None,projectName=None):
    QtWidgets.QDialog.__init__(self,parent)
    #self.setModal(True)
    self.setWindowTitle('Options for saving project: '+str(projectName))
    self.setLayout(QtWidgets.QVBoxLayout())

    self.modeGroup = QtWidgets.QButtonGroup(self)
    self.modeGroup.setExclusive(True)

    line = QtWidgets.QHBoxLayout()
    but = QtWidgets.QRadioButton('Export entire project',self)
    self.modeGroup.addButton(but,1)
    but.setChecked(True)
    line.addWidget(but)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    but = QtWidgets.QRadioButton('Export every job after previous import/export',self)
    self.modeGroup.addButton(but,2)
    line.addWidget(but)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    self.impExpWidget = CImportExportListWidget(self,projectId=projectId)
    line.addSpacing(20)
    line.addWidget(self.impExpWidget)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    but = QtWidgets.QRadioButton('Select jobs to export',self)
    self.modeGroup.addButton(but,3)
    line.addWidget(but)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    self.selectionWidget = CCP4Widgets.CJobSelectionLineEdit(self,projectId)
    line.addSpacing(20)
    line.addWidget(self.selectionWidget)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton('Export',QtWidgets.QDialogButtonBox.AcceptRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.handleExport)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.close)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Help)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.help)
    line.addStretch(0.5)
    line.addWidget(buttonBox)
    line.addStretch(0.5)
    self.layout().addLayout(line)

  @QtCore.Slot()
  def handleExport(self):
    self.exportSignal.emit(self.modeGroup.checkedId())

  @QtCore.Slot()
  def help(self):
    pass

class CExportJobSelection(QtWidgets.QLineEdit):
    def __init__(self,parent,projectId=None):
      QtWidgets.QLineEdit.__init__(self,parent)
      self.projectId = projectId
      self.setToolTip("Enter list of jobs e.g. '27-29,31'")

    def getSelection(self):
      errList = []
      seleList = []
      text = str(self.text())
      splitList = text.split(',')
      #print 'CExportJobSelection.getSelection',text,type(text)
      for item in splitList:
        #print 'CExportJobSelection.getSelection split',item
        rSplit = item.split('-')
        if len(rSplit)==1:
          try:
            jobId = PROJECTSMANAGER().db().getJobId(projectId=self.projectId,jobNumber=item.strip())
          except:
            errList.append(item)
          else:
            seleList.append(jobId)
        elif len(rSplit)==2:
          try:
            jobList = PROJECTSMANAGER().db().getJobsInRange(projectId=self.projectId,jobNumberRange=[rSplit[0].strip(),rSplit[1].strip()])
          except:
            errList.append(item)
          else:
            if len(jobList)==0:
              errList.append(item)
            else:
              seleList.extend(jobList)
        else:
           errList.append(item)
      #print 'CExportJobSelection.getSelection', seleList,errList
      return seleList,errList

class CImportExportListWidget(QtWidgets.QComboBox):

    def __init__(self,parent,projectId=None):
      QtWidgets.QComboBox.__init__(self,parent)
      self.projectId = projectId

      self.setEditable(False)
      self.load()

    def load(self):
      exportList = PROJECTSMANAGER().db().getProjectExportInfo(projectId=self.projectId)
      importList = PROJECTSMANAGER().db().getProjectImportInfo(projectId=self.projectId)
      #print 'CImportExportListWidget.load',len(exportList),len(importList)
      def formatImport(data):
        imDate = time.strftime( '%a %d %b %H:%M' , time.localtime(data[1]))
        exDate = time.strftime( '%a %d %b %H:%M' , time.localtime(data[2]))
        text = 'Import '+imDate+' from '+importList[iIm][3]+ ' on '+exDate
        return text,data[0]
      def formatExport(data):
        exDate = time.strftime( '%a %d %b %H:%M' , time.localtime(data[1]))
        text = 'Export '+exDate
        return text,data[0]
      iEx = 0
      iIm = 0
      while ( iEx<len(exportList) or iIm<len(importList) ):
        #print 'CImportExportListWidget.load',iEx,iIm
        if iEx>=len(exportList):
          text,qVar = formatImport(importList[iIm])
          iIm += 1
        elif iIm>=len(importList):
          text,qVar = formatExport(exportList[iEx])
          iEx += 1
        elif importList[iIm][1]>exportList[iEx][1]:
          text,qVar = formatImport(importList[iIm])
          iIm += 1
        else:
          text,qVar = formatExport(exportList[iEx])
          iEx += 1
        self.addItem(text,qVar)

    def getCurrentItem(self):
      qVar = self.itemData(self.currentIndex())
      #print 'CImportExportListWidget.getCurrentItem',qVar
      return qVar.__str__()

    def getTime(self):
      tid = self.getCurrentItem()
      try:
        info = PROJECTSMANAGER().db().getProjectExportInfo(projectExportId=tid)
      except:
        info = {}
      if info.get('projectexporttime',None) is not None:
        return info['projectexporttime']
      else:
        info = PROJECTSMANAGER().db().getProjectImportInfo(projectImportId=tid)
        return info['projectimporttime']


class CImportNewProjectDialog(QtWidgets.QDialog):

  importSignal = QtCore.Signal()

  def __init__(self,parent=None):
    MARGIN = 2
    QtWidgets.QDialog.__init__(self,parent)
    #self.setModal(True)
    self.setWindowTitle('Import a New Project?')
    self.setLayout(QtWidgets.QVBoxLayout())

    from ..core import CCP4DataManager
    from ..core import CCP4File

    self.projectInfoDisplay = CProjectInfoDisplay(self)
    self.layout().addWidget(self.projectInfoDisplay)

    line = QtWidgets.QHBoxLayout()
    lab = QtWidgets.QLabel('The project id in the file does not match any existing project.')
    lab.setObjectName('emphasise')
    line.addWidget(lab)
    self.layout().addLayout(line)
    self.modeButs = QtWidgets.QButtonGroup(self)
    self.modeButs.setExclusive(True)

    lab =  QtWidgets.QLabel('Create a new project directory',self)
    but =  QtWidgets.QRadioButton('Create a new project directory',self)
    but.setChecked(True)
    but.setVisible(False)
    self.modeButs.addButton(but,1)
    self.layout().addWidget(lab)
    self.layout().addWidget(but)

    line = QtWidgets.QFrame()
    line.setLayout(QtWidgets.QHBoxLayout())
    line.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
    line.layout().setSpacing(MARGIN)
    line.layout().addSpacing(30)
    self.directory0 = CCP4File.CDataFile(parent=self,qualifiers={ 'isDirectory' : True,'mustExist': False,'fromPreviousJob' : False  })
    self.directoryWidget0 = CCP4DataManager.DATAMANAGER().widget(model=self.directory0,parentWidget=self,qualifiers={'jobCombo': False, 'projectBrowser' : False} )
    self.directoryWidget0.fileLineEdit.setCharWidth(60,mode='minimum')
    line.layout().addWidget(self.directoryWidget0)
    self.layout().addWidget(line)

    line = QtWidgets.QHBoxLayout()
    line.addSpacing(30)
    line.addWidget(QtWidgets.QLabel('Project name',self))
    self.projectNameWidget = QtWidgets.QLineEdit(self)
    line.addWidget(self.projectNameWidget)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton('Import',QtWidgets.QDialogButtonBox.AcceptRole)
    but.clicked.connect(self.handleAccept)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.setDefault(True)
    but.clicked.connect(self.close)
    line.addStretch(0.1)
    line.addWidget(buttonBox)
    line.addStretch(0.1)
    self.layout().addLayout(line)

  def reset(self,importProjectInfo):
    self.directory0.unSet()
    self.projectInfoDisplay.load(importProjectInfo)
    self.projectNameWidget.setText(importProjectInfo['projectName'])
    directory = CCP4Utils.getProjectDirectory()
    if directory is not None:
        projName = CCP4Utils.safeOneWord(importProjectInfo['projectName'])
        projDir = os.path.join(directory,projName)
        #FIXME - This only puts a new name into the widget. The user is still free to override this choice with the name
        #of an existing file. We do not trap such behaviour. We should and forbid it somehow.
        if os.path.exists(projDir):
            projDir += "_"
            for i in range(10000):
                if not os.path.exists(projDir+str(i)):
                    projDir = projDir+str(i)
                    self.directory0.set(projDir)
                    break
        else:
            self.directory0.set(projDir)

  @QtCore.Slot()
  def handleAccept(self):
    if self.modeButs.checkedId() == 1:
      if not self.directory0.isSet():
        QtWidgets.QMessageBox.warning(self,'Import a New Project?','Please enter a new project directory')
        return
      projectName = self.projectName()
      if len(projectName)==0:
        QtWidgets.QMessageBox.warning(self,'Import a New Project?','Please enter a name for the new project')
        return
      try:
        projectId = PROJECTSMANAGER().db().getProjectId(projectName = projectName)
      except:
        pass
      else:
        QtWidgets.QMessageBox.warning(self,'Import a New Project?','There is already a project in the database called '+projectName)
        return
    else:
      pid = self.dbTreeWidget.selectedProjectId()
      if pid is None:
        QtWidgets.QMessageBox.warning(self,'Import a New Project?','Please select a project to merge with')
        return
    self.importSignal.emit()

  def mode(self):
    return self.modeButs.checkedId()

  def newDirectory(self):
    return  self.directory0.__str__()

  def projectName(self):
    name = str(self.projectNameWidget.text()).strip()
    if len(name)==0: name = None
    return name

  def mergeProject(self):
    return self.dbTreeWidget.selectedProjectId()

class CImportExistingProjectDialog(QtWidgets.QDialog):

  importSignal = QtCore.Signal()

  def __init__(self,parent=None,projectName=None):
    MARGIN = 2
    QtWidgets.QDialog.__init__(self,parent)
    #self.setModal(True)
    self.setWindowTitle('Import to an Existing Project?')
    self.setLayout(QtWidgets.QVBoxLayout())

    self.projectInfoDisplay = CProjectInfoDisplay(self)
    self.layout().addWidget(self.projectInfoDisplay)

    line = QtWidgets.QHBoxLayout()
    lab = QtWidgets.QLabel('The project id in the file matches project '+projectName)
    lab.setObjectName('emphasise')
    line.addWidget(lab)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget( QtWidgets.QLabel('Import jobs to that existing project'))
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton('Import',QtWidgets.QDialogButtonBox.AcceptRole)
    but.clicked.connect(self.importSignal.emit)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.setDefault(True)
    but.clicked.connect(self.close)
    line.addStretch(0.1)
    line.addWidget(buttonBox)
    line.addStretch(0.1)
    self.layout().addLayout(line)

  def reset(self,importProjectInfo):
    self.projectInfoDisplay.load(importProjectInfo)


class CProjectInfoDisplay(QtWidgets.QFrame):
  def __init__(self,parent,projectInfo={}):
    QtWidgets.QFrame.__init__(self,parent)
    self.setObjectName('highlight')
    self.setLayout(QtWidgets.QVBoxLayout())
    self.widgets = {}
    for item in ['projectName','hostName','userId','creationTime']:
       self.widgets[item] = QtWidgets.QLabel(self)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('The compressed file contains data for project: ',self))
    line.addWidget(self.widgets['projectName'])
    line.addStretch(0.5)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Created on',self))
    line.addWidget(self.widgets['hostName'])
    line.addWidget(QtWidgets.QLabel('by:',self))
    line.addWidget(self.widgets['userId'])
    line.addWidget(QtWidgets.QLabel('on:',self))
    line.addWidget(self.widgets['creationTime'])
    line.addStretch(0.5)
    self.layout().addLayout(line)

    self.load(projectInfo)

  def load(self,projectInfo):
    #print 'CProjectInfoDisplay.load',projectInfo
    for item in ['projectName','hostName','userId']:
       self.widgets[item].setText(str(projectInfo.get(item,'')))
    t = projectInfo.get('creationTime',None)
    if t is None:
      text = ''
    else:
      text = time.strftime(PROJECTSMANAGER().db().TIMEFORMAT,time.localtime(float(t)))
    self.widgets['creationTime'].setText(text)

class ProjectTreeItem(QtGui.QStandardItem):
    #FIXME - Reparenting is required.
    def __init__(self,data,parentItem=None):
        QtGui.QStandardItem. __init__(self)
        self.m_parentItem = parentItem
        self.m_itemData = data
        self.m_childItems = []
        self._icon = None
        self._itemType = 0
        self._pid = -1
        self._name = ""
        self._directory = ""
        self._createdDate = None
        self._accessDate = None
        self._tags = ""
        self._annotation = ""

    def setName(self,name):
        self._name = name

    def name(self):
        return self._name

    def pid(self):
        return self._pid

    def tags(self):
        return self._tags

    def setTags(self,tags):
        self._tags = tags

    def setAnnotation(self,annotation):
        self._annotation = annotation

    def setAccessDate(self,date):
        self._accessDate = date

    def setCreatedDate(self,date):
        self._createdDate = date

    def setId(self,pid):
        self._pid = pid

    def setDirectory(self,directory):
        self._directory = directory

    def parentItem(self):
        return self.m_parentItem

    def setParentItem(self,parent):
        self.m_parentItem = parent

    def removeChild(self,child):
        self.m_childItems.remove(child)

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
        return 6

    def data(self,column):
        if column == 0:
            return self.m_itemData
        elif column == 1:
            return self._directory
        elif column == 2:
            return self._createdDate
        elif column == 3:
            return self._accessDate
        elif column == 4:
            return self._tags
        elif column == 5:
            return self._annotation
        elif column == 6:
            return self._pid

    def setIcon(self,icon):
        self._icon = icon

    def icon(self):
        return self._icon

class CProjectsViewModel(QtCore.QAbstractItemModel):
    def __init__(self,parent=None):
        QtCore.QAbstractItemModel.__init__(self,parent)
        self.rootItem = ProjectTreeItem("Root")

    def headerData(self, section, orientation, role = QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if section == 0:
                return "Name"
            elif section == 1:
                return "Directory"
            elif section == 2:
                return "Created"
            elif section == 3:
                return "Last active"
            elif section == 4:
                return "Tags"
            elif section == 5:
                return "Annotation"

#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def flags(self,index):

        defaultFlags = QtCore.QAbstractItemModel.flags(self,index);

        if index.isValid():
            return QtCore.Qt.ItemIsDragEnabled | QtCore.Qt.ItemIsDropEnabled | defaultFlags;
        else:
            return QtCore.Qt.ItemIsDropEnabled | defaultFlags;

    def data(self, index, role):
        #FIXME - Maybe I should have columns instead of all these roles? Maybe I do already?
        if not index.isValid():
#FIXME PYQT - or maybe None? This used to return QVariant.
            return None

        item = index.internalPointer()

        if role == QtCore.Qt.UserRole + 8:
            return item._pid

        if role == QtCore.Qt.UserRole + 7:
            return item._annotation

        if role == QtCore.Qt.UserRole + 6:
            return item._tags

        if role == QtCore.Qt.UserRole + 5:
            return item._accessDate

        if role == QtCore.Qt.UserRole + 4:
            return item._createdDate

        if role == QtCore.Qt.UserRole + 3:
            return item._directory

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

class CProjectsViewProxyModel(QtCore.QSortFilterProxyModel):
    def __init__(self,parent=None):
        QtCore.QSortFilterProxyModel.__init__(self,parent)
        self.filterString = ""
        self.filterStartDate = None
        self.filterEndDate = None

    def setFilterString(self,string):
        self.filterString = string
        self.invalidateFilter()

    def setFilterDates(self,startDate,endDate):
        self.filterStartDate = startDate
        self.filterEndDate = endDate
        self.invalidateFilter()

    def checkDates(self,item):
        startDate = item.data(2)
        endDate = item.data(3)
        accepted = True
        if self.filterStartDate is not None and self.filterEndDate is not None:
            t = datetime.datetime.strptime(self.filterStartDate,"%Y-%m-%d %H:%M:%S")
            stampStart = time.mktime(t.timetuple())
            if startDate < stampStart:
                accepted = False
            t = datetime.datetime.strptime(self.filterEndDate,"%Y-%m-%d %H:%M:%S")
            stampEnd = time.mktime(t.timetuple())
            if endDate > stampEnd and stampEnd > stampStart: # Second clause because sometimes end date is undefined.
                accepted = False
        return accepted

    def hasAcceptedChildrenOrIsAccepted(self, item):

            theData = item.data(0)
            theAnnotation = item.data(5)
            theTags = item.data(4)

            theData += " " + theAnnotation + " " + theTags

            if (len(str(self.filterString).strip()) == 0 or str(self.filterString).lower() in theData.lower()):
                return self.checkDates(item)

            accepted = False
            for ic in range(item.childCount()):
                if self.hasAcceptedChildrenOrIsAccepted(item.child(ic)):
                    accepted = True
            if accepted:
                return self.checkDates(item)
            return accepted

    def filterAcceptsRow(self, sourceRow, sourceParent):
        source_index = self.sourceModel().index(sourceRow, 0, sourceParent);
        ip = sourceParent.internalPointer()

        return self.hasAcceptedChildrenOrIsAccepted(source_index.internalPointer())

class CProjectsTreeView(QtWidgets.QTreeView):

    selectionIndexesChanged = QtCore.Signal(int)
    DROP_ITEM_TOL = 6

    def __init__(self,parent=None):
        QtWidgets.QTreeView.__init__(self,parent)
        self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.setUniformRowHeights(True)
        self.dropSite = None

    def selectionChanged(self,sel,desel):
        self.selectionIndexesChanged.emit(len(self.selectionModel().selectedIndexes()))
        QtWidgets.QTreeView.selectionChanged(self,sel,desel)

    def sizeHint(self):
        return QtCore.QSize(600,200)

    def mimeData(self,indexes):
        encodedData = QtCore.QByteArray()
        for widgetIndex in indexes:
            pid = self.model().sourceModel().data(widgetIndex,QtCore.Qt.UserRole+8)
            encodedData.append(bytes(str(pid)+",","utf-8"))
        # With mime type as text the data can be dropped on desktop
        # but the type of the data is lost
        mimeData = QtCore.QMimeData()
        mimeData.setData('project',encodedData)
        return mimeData

    def dragEnterEvent(self,event):
        if event.mimeData().hasFormat('project'):
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self,event):
        self.dropSite = event.answerRect()
        self.repaint()
        QtWidgets.QTreeView.dragMoveEvent(self,event)
        idx = self.indexAt(event.pos())
        if not idx.isValid():
            targetProjectId = None
        else:
            idx_s = idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0)
            item = idx_s.internalPointer()
            targetProjectId = self.model().sourceModel().data(idx_s,QtCore.Qt.UserRole+8)
        movedProject = event.mimeData().data('project').data()
        if movedProject == targetProjectId:
            event.setDropAction(QtCore.Qt.IgnoreAction)
            event.ignore()
            return


        event.setDropAction(QtCore.Qt.MoveAction)
        event.accept()

    def paintEvent(self,event):
        QtWidgets.QTreeView.paintEvent(self,event)

        if self.dropSite is not None:
            painter = QtGui.QPainter(self.viewport())
            painter.save()
            point = QtCore.QPoint(self.dropSite.x(),self.dropSite.y())
            idx = self.indexAt(point)
            if idx.isValid():
                arect = self.visualRect(idx);
                b = arect.y()
                if abs(b - self.dropSite.y()) < self.DROP_ITEM_TOL:
                    painter.drawLine ( 1, b, self.viewport().width()-3, b )
                    painter.drawLine ( 1, b-4, 1, b+4 )
                    painter.drawLine ( self.viewport().width()-3, b-4, self.viewport().width()-3, b+4 )
                else:
                    painter.drawRect ( 1, b, self.viewport().width()-3, self.rowHeight(idx) )
            painter.restore()

        event.accept()  

    def dropEvent(self,event):
        idx = self.indexAt(event.pos())

        if not idx.isValid():
            targetProjectId = None
        else:
            if self.dropSite is not None:
                arect = self.visualRect(idx);
                b = arect.y()
                if abs(b - self.dropSite.y()) < self.DROP_ITEM_TOL:
                    targetProjectId = None
                else:
                    idx_s = idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0)
                    item = idx_s.internalPointer()
                    targetProjectId = self.model().sourceModel().data(idx_s,QtCore.Qt.UserRole+8)
            else:
                #Just in case dropSite is None.
                idx_s = idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0)
                item = idx_s.internalPointer()
                targetProjectId = self.model().sourceModel().data(idx_s,QtCore.Qt.UserRole+8)

        self.dropSite = None

        movedProjects = event.mimeData().data('project').data().decode()
        movedProjects_a = movedProjects.strip(",").split(",")

        for movedProject in movedProjects_a:
            if movedProject == targetProjectId:
                continue
            try:
                PROJECTSMANAGER().db().updateProject(movedProject,key='parentProjectId',value=targetProjectId)
            except Exception as e:
                print('Failed to move project',e)

        event.setDropAction(QtCore.Qt.MoveAction)
        event.accept()
        self.repaint()

    def startDrag(self,dropActions):

        if len(self.selectedIndexes()) == 0:
            return

        idxs = []
        for idx in  self.selectedIndexes():
            if idx.column() == 0:
                idx_s = idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0)
                item = idx_s.internalPointer()
                pid = self.model().sourceModel().data(idx_s,QtCore.Qt.UserRole+8)
                idxs.append(idx_s)

        mimeData = self.mimeData(idxs)
        drag = QtGui.QDrag(self)
        drag.setMimeData(mimeData)
        pixmap = item.icon().pixmap(18,18)
        drag.setHotSpot(QtCore.QPoint(9,9))
        drag.setPixmap(pixmap)

        if drag.exec_(QtCore.Qt.MoveAction) == QtCore.Qt.MoveAction:
            pass#print 'Dragging...'

    def selectedProjectId(self):
        if len(self.selectionModel().selectedIndexes()):
            idx = self.selectionModel().selectedIndexes()[0]
            oldName = self.model().sourceModel().data(idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0),QtCore.Qt.DisplayRole)
            pid = self.model().sourceModel().data(idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0),QtCore.Qt.UserRole+8)
            return pid

        return None

class CProjectsTreeWidget(QtWidgets.QTreeWidget):


  rightMousePress = QtCore.Signal('QMouseEvent')

  PROJECTICON = None

  def __init__(self,parent):
    QtWidgets.QTreeWidget.__init__(self,parent)
    self.projInfo = {}
    self.setColumnCount(2)
    self.expandAll()
    self.setHeaderLabels(['Name','Directory','Created','Last active','Tags'])
    self.setItemsExpandable(True)
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    self.setDragEnabled(True)
    self.setAcceptDrops(True)
    self.header().resizeSection(0,300)
    #self.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
    PROJECTSMANAGER().db().projectTagsChanged.connect(self.handleProjectTagsChanged)

  def projectIcon(self):
    if CProjectsTreeWidget.PROJECTICON is None:
      fileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','project.png')
      CProjectsTreeWidget.PROJECTICON = QtGui.QIcon(QtGui.QPixmap(fileName))
    return CProjectsTreeWidget.PROJECTICON

  def populate(self,args={}):
    self.clear()
    proj_dir_list0=PROJECTSMANAGER().db().getProjectDirectoryList()
    #print 'CProjectsView.populate',proj_dir_list
    proj_dir_list = []
    self.projInfo = {}
    for pid,pName,pDir,parent,created,lastAccess in proj_dir_list0:
      if parent is None: proj_dir_list.append(pid)
      self.projInfo[pid] = { 'name' : pName,
                        'dir' : pDir,
                        'children' : [],
                        'created' : created,
                        'lastAccess' : lastAccess,
                        'tags' : '',
                        'annotation' : '' }

    for pid,pName,pDir,parent,created,lastAccess in proj_dir_list0:
      if parent is not None:
        self.projInfo[parent]['children'].append(pid)

    pTagList = PROJECTSMANAGER().db().getProjectTagList()
    for pid,tagText in pTagList:
      self.projInfo[pid]['tags'] += ','+tagText

    commentList = PROJECTSMANAGER().db().getProjectCommentsList()
    for pid,comment in commentList:
      self.projInfo[pid]['annotation'] = comment
    #print 'CProjectsTreeWidget',self.projInfo
    self.drawTree(projectList=proj_dir_list,parent=None)

  @QtCore.Slot(str)
  def handleProjectTagsChanged(self,projectId):
    #print 'CProjectsTreeWidget.handleProjectTagsChanged',projectId,self.projInfo[projectId]['name']
    tagList = PROJECTSMANAGER().db().getProjectTags(projectId,tagText=True)
    text = ''
    for tid,tagText in tagList:
      text += ','+tagText
    treeWidgetItemList = self.findItems(self.projInfo[projectId]['name'],QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive,0)
    #print 'handleProjectTagsChanged',treeWidgetItemList,text
    if len(treeWidgetItemList)==1:
      treeWidgetItemList[0].setData(4,QtCore.Qt.DisplayRole,text[1:])
      self.update(self.indexFromItem(treeWidgetItemList[0],4))

  def makeTreeWidgetItem(self,pid):

    t1 = formatDate(self.projInfo[pid]['created'])
    t2 = formatDate(self.projInfo[pid]['lastAccess'])
    item = QtWidgets.QTreeWidgetItem([self.projInfo[pid]['name'],self.projInfo[pid]['dir'],t1,t2,self.projInfo[pid]['tags'][1:]],1001)
    item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsDragEnabled)
    item.setData(0,QtCore.Qt.UserRole,pid)
    item.setData(0,QtCore.Qt.DecorationRole,self.projectIcon())
    item.setData(0,QtCore.Qt.ToolTipRole,self.projInfo[pid]['annotation'])
    return item

  def drawTree(self,projectList=None,parent=None):
    for pid in projectList:
      item = self.makeTreeWidgetItem(pid)
      if parent is None:
        self.addTopLevelItem(item)
      else:
        parent.addChild(item)
      if len(self.projInfo[pid]['children'])>0:
        self.drawTree(self.projInfo[pid]['children'],item)

  def setHighlights(self,treeItem=None,highlightList=[],colour=['red','pink']):
    #print 'setHighlights',highlightList
    nHits = 0
    if treeItem is None:
      for n in range(self.topLevelItemCount()):
        item = self.topLevelItem(n)
        childHits = self.setHighlights(item,highlightList,colour)
        if highlightList.count(item.data(0,QtCore.Qt.UserRole).__str__()):
          item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor(colour[0])))
          nHits += 1
        elif childHits>0:
         item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor(colour[1])))
        else:
          item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor('black')))
        #if childHits>0: item.setExpanded(True)
        nHits+= childHits
      self.update()
    else:
      for n in range(treeItem.childCount()):
        item = treeItem.child(n)
        childHits = self.setHighlights(item,highlightList,colour)
        if highlightList.count(item.data(0,QtCore.Qt.UserRole).__str__()):
          item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor(colour[0])))
          nHits += 1
        elif childHits>0:
         item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor(colour[1])))
        else:
          item.setData(0,QtCore.Qt.ForegroundRole,QtGui.QBrush(QtGui.QColor('black')))
        #if childHits>0: item.setExpanded(True)
        nHits+= childHits
    #print 'setHighlights nHits',nHits
    return nHits 


  def sizeHint(self):
    return QtCore.QSize(600,200)


  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
      self.rightMousePress.emit(event)      
    QtWidgets.QTreeWidget.mousePressEvent(self,event)

  def selectedItem(self):
    seleList = self.selectedItems()
    #print 'contentsTree.selectedItem',seleList
    if len(seleList) == 0:
      return None
    elif len(seleList) == 1:
      return seleList[0]
    else:
      #print 'Error in CContentsTree - more than one selected item'
      return seleList[0]

  def mimeData(self,widgetItemList):
    from ..core import CCP4Data
    #print 'CProjectsTreeWidget.mimeData',widgetItemList
    encodedData = QtCore.QByteArray()
    for widgetItem in widgetItemList:
      pid = CCP4Data.varToUUID(widgetItem.data(0,QtCore.Qt.UserRole))
      #print 'CProjectsTreeWidget.mimeData',pid
      encodedData.append(str(pid))
    # With mime type as text the data can be dropped on desktop
    # but the type of the data is lost
    mimeData = QtCore.QMimeData()
    mimeData.setData('project',encodedData)
    #mimeData.setText('project '+str(pid))
    return mimeData

  def dragEnterEvent(self,event):
    if event.mimeData().hasFormat('project'):
      event.accept()
    else:
      event.ignore()

  def dragMoveEvent(self,event):
    dropItem = self.itemAt(event.pos().x(),event.pos().y())
    #print 'dragMoveEvent',event.mimeData().hasFormat('project')
    #if event.mimeData().hasFormat('project') and dropItem is not None:
    if event.mimeData().hasFormat('project'):
      from ..dbapi import CCP4DbApi
      movedProject = CCP4DbApi.UUIDTYPE(event.mimeData().data('project').data())
      targetProjectId = self.item2ProjectId(dropItem)
      if movedProject == targetProjectId:
        event.ignore()
      else:
        event.setDropAction(QtCore.Qt.MoveAction)
        event.accept()
    else:
      event.ignore()

  def dropEvent(self,event):
    targetItem = self.itemAt(event.pos().x(),event.pos().y())
    #print 'CProjectsTreeWidget.dropEvent',event,targetItem  
    if event.mimeData().hasFormat('project'):
      from ..dbapi import CCP4DbApi
      movedProject = CCP4DbApi.UUIDTYPE(event.mimeData().data('project').data())
      targetProjectId = self.item2ProjectId(targetItem)
      if movedProject == targetProjectId:
        print('ERROR trying to drop project in itself')
        return
      #print 'dropEvent targetProject',targetItem,targetProjectId
      try:
        PROJECTSMANAGER().db().updateProject(movedProject,key='parentProjectId',value=targetProjectId)
      except CException as e:
        print('Failed to move project',e.report())
        event.setDropAction(QtCore.Qt.IgnoreAction)
        event.ignore()
        return
      except Exception as e:
        print('Failed to move project',e)
        event.setDropAction(QtCore.Qt.IgnoreAction)
        event.ignore()
        return
      #self.drawTree([movedProject],targetItem)
      #if targetProjectId is not None:
      #  self.projInfo[targetProjectId]['children'].append(movedProject)
      self.populate()
      event.setDropAction(QtCore.Qt.MoveAction)
      event.accept()
    else:
      event.ignore()

  def startDrag(self,dropActions):
    item = self.currentItem()
    projectId = self.item2ProjectId(item)
    mimeData = self.mimeData([item])
    drag = QtGui.QDrag(self)
    drag.setMimeData(mimeData)
    pixmap = item.icon(0).pixmap(18,18)
    drag.setHotSpot(QtCore.QPoint(9,9))
    drag.setPixmap(pixmap)
    if drag.exec_(QtCore.Qt.MoveAction) == QtCore.Qt.MoveAction:
      #self.takeItem(self.row(item))
      try:
        if item.parent() is None:
          indx = self.indexOfTopLevelItem(item)
          self.takeTopLevelItem(indx)
        else:
          parent = item.parent()
          parent.removeChild(item)
          parentId = self.item2ProjectId(parent)
          self.projInfo[parentId]['children'].remove(projectId)
      except:
          pass
          #print 'WARNING from CProjectsTreeWidget.startDrag'

  def selectedProjectId(self):
    selList = self.selectedItems()
    #print 'CProjectsTreeWidget.selectedProjectId',selList
    if len(selList)==0: return None
    return self.item2ProjectId(selList[0])

  def item2ProjectId(self,item):
    if item is None: return None
    from ..core import CCP4Data
    projectId = CCP4Data.varToUUID(item.data(0,QtCore.Qt.UserRole))
    return  projectId


class CProjectsStatusBar(QtWidgets.QStatusBar):
  def __init__(self,parent):
    QtWidgets.QStatusBar.__init__(self,parent)
    self.setObjectName('statusWidget')
    self.setSizeGripEnabled(False)
    self.progressWidget =QtWidgets.QProgressBar(self)
    self.progressWidget.hide()

  def showMessage(self,text=None,timeout=0):
    #print 'CProjectsStatusBar.showMessage',text
    QtWidgets.QStatusBar.showMessage(self,text,timeout*1000)
    self.show()

  @QtCore.Slot()
  def clear(self):
    self.clearMessage()
    self.hideProgress()

  def showProgress(self,max=None,min=None):
    self.addPermanentWidget(self.progressWidget)
    self.progressWidget.show()
    self.progressWidget.setMaximum(max)

  @QtCore.Slot(int)
  def setProgress(self,value):
    self.progressWidget.setValue(value)

  def hideProgress(self):
    self.removeWidget(self.progressWidget)


def repopulateTreeViewNew(proxyModel,projectsView):
       projectsViewModel = CProjectsViewModel()

       rootItem = projectsViewModel.rootItem

       iconFileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','project.png')
       projectIcon = QtGui.QIcon(iconFileName)

       oldRow = -1
       oldName = ''

       projectsView.selectionIndexesChanged.emit(0)

       if proxyModel.sourceModel() is not None:
           idxs = projectsView.selectionModel().selectedIndexes()
           for idx in idxs:
               if idx.column() == 0:
                   oldName = projectsView.model().sourceModel().data(idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0),QtCore.Qt.DisplayRole)
                   oldRow = idx.row()

       proxyModel.setSourceModel(projectsViewModel)
       projectsView.setModel(proxyModel)

       proj_dir_list0=PROJECTSMANAGER().db().getProjectDirectoryList()
       pTagList = PROJECTSMANAGER().db().getProjectTagList()
       commentList = PROJECTSMANAGER().db().getProjectCommentsList()

       projectItems = {}
       for pid,pName,pDir,parent,created,lastAccess in proj_dir_list0:
           projectItem  = ProjectTreeItem(pName,rootItem)
           projectItem.setId(pid)
           projectItem.setDirectory(pDir)
           projectItem.setIcon(projectIcon)
           projectItem.setCreatedDate(created)
           projectItem.setAccessDate(lastAccess)
           rootItem.appendChild(projectItem)
           projectItems[pid] = projectItem

       for pidt,tagText in pTagList:
           try:
               if len(projectItems[pidt].tags()) == 0:
                   projectItems[pidt].setTags(tagText)
               else:
                   projectItems[pidt].setTags(projectItems[pidt].tags()+","+tagText)
           except:
               print("Failed to set tags. Possibly orphaned tag with no existing project.")

       for pidt,tagText in commentList:
           try:
               projectItems[pidt].setAnnotation(tagText)
           except:
               print("Failed to set comments. Possibly orphaned comment with no existing project.")

       for pid,pName,pDir,parent,created,lastAccess in proj_dir_list0:
           if parent is not None:
               item = projectItems[pid]
               parentItem = projectItems[parent]
               rootItem.removeChild(item)
               parentItem.appendChild(item)
               item.setParentItem(parentItem)

       proxyModel.beginResetModel()

       if oldRow >-1 and len(oldName) > 0:
           idx = projectsView.model().index(oldRow,0)
           if idx.row() > -1 and idx.column() > -1:
               newName = projectsView.model().sourceModel().data(idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0),QtCore.Qt.DisplayRole)
               if oldName == newName:
                   projectsView.selectionModel().select(QtCore.QItemSelection(idx,idx.sibling(idx.row(),5)),QtCore.QItemSelectionModel.Select)
       proxyModel.endResetModel()

       projectsView.header().resizeSection(0,250)
       projectsView.setSortingEnabled(True)
       projectsView.sortByColumn(0,QtCore.Qt.AscendingOrder)


class CProjectManagerDialog(QtWidgets.QDialog):

  doneSavingJobDataSignal = QtCore.Signal()
  startSavingJobDataSignal = QtCore.Signal()

  def __init__(self,parent=None):
    QtWidgets.QDialog.__init__(self,parent)
    self.newProjectGui = None
    self.importProjectManager = None

    layout = QtWidgets.QGridLayout()
    self.setLayout(layout)
    self.setWindowTitle("Manage projects")

    @QtCore.Slot('QModelIndex')
    def myDoubleClickHandler(idx):
        projName = self.projectsView.model().sourceModel().data(idx.model().mapToSource(idx).sibling(idx.model().mapToSource(idx).row(),0),QtCore.Qt.DisplayRole)
        openProject(projectName=projName)

    self.projectsView = CProjectsTreeView()
    self.projectsView.doubleClicked.connect(myDoubleClickHandler)
    if sys.platform == "darwin": 
        self.projectsView.setIconSize(QtCore.QSize(16,16))
    self.projectsViewProxyModel = CProjectsViewProxyModel()
    repopulateTreeViewNew(self.projectsViewProxyModel,self.projectsView)

    pvd = ProjectsViewDelegate()
    self.projectsView.setItemDelegate(pvd)

    self.buttonFrame = QtWidgets.QDialogButtonBox(QtCore.Qt.Vertical)

    openButton = self.buttonFrame.addButton('Open',QtWidgets.QDialogButtonBox.ActionRole)
    openButton.setDefault(True)
    openButton.setEnabled(False)
    openButton.clicked.connect(self.openProject)
    addButton = self.buttonFrame.addButton('Add project or folder',QtWidgets.QDialogButtonBox.ActionRole)
    addButton.clicked.connect(self.addProject)
    renameButton = self.buttonFrame.addButton('Rename project',QtWidgets.QDialogButtonBox.ActionRole)
    renameButton.setEnabled(False)
    renameButton.clicked.connect(self.handleRenameProject)
    editDescriptionButton = self.buttonFrame.addButton('Edit description',QtWidgets.QDialogButtonBox.ActionRole)
    editDescriptionButton.setEnabled(False)
    editDescriptionButton.clicked.connect(self.handleEditDescription)
    moveButton = self.buttonFrame.addButton('Move directory',QtWidgets.QDialogButtonBox.ActionRole)
    moveButton.setEnabled(False)
    moveButton.clicked.connect(self.handleMoveProject)
    deleteButton = self.buttonFrame.addButton('Delete',QtWidgets.QDialogButtonBox.ActionRole)
    deleteButton.setEnabled(False)
    deleteButton.clicked.connect(self.handleDeleteProject)
    cleanupButton = self.buttonFrame.addButton('Cleanup files',QtWidgets.QDialogButtonBox.ActionRole)
    cleanupButton.clicked.connect(self.handleCleanup)
    exportButton = self.buttonFrame.addButton('Export',QtWidgets.QDialogButtonBox.ActionRole)
    exportButton.setEnabled(False)
    exportButton.clicked.connect(self.handleExport1)
    importButton = self.buttonFrame.addButton('Import',QtWidgets.QDialogButtonBox.ActionRole)
    importButton.clicked.connect(self.handleImportProject)
    exportAlltButton = self.buttonFrame.addButton('Export all projects',QtWidgets.QDialogButtonBox.ActionRole)
    exportAlltButton.clicked.connect(self.handleExportAll)
    openButton.setAutoDefault(False)
    renameButton.setAutoDefault(False)
    editDescriptionButton.setAutoDefault(False)
    moveButton.setAutoDefault(False)
    deleteButton.setAutoDefault(False)
    cleanupButton.setAutoDefault(False)
    exportButton.setAutoDefault(False)

    if CONFIG().developer:
        rerunButton = self.buttonFrame.addButton('Rerun test project',QtWidgets.QDialogButtonBox.ActionRole)
        rerunButton.setEnabled(False)
        rerunButton.clicked.connect(self.handleRerun)
        rerunButton.setAutoDefault(False)

    @QtCore.Slot(int)
    def checkButtonStates(lenSel):
        #FIXME - Want to know if I have a single row selected. lenSel / rowCount is what is wanted.
        if lenSel == 6:
            enab = True
        else:
            enab = False
        openButton.setEnabled(enab)
        renameButton.setEnabled(enab)
        editDescriptionButton.setEnabled(enab)
        moveButton.setEnabled(enab)
        deleteButton.setEnabled(enab)
        exportButton.setEnabled(enab)
        if CONFIG().developer:
            rerunButton.setEnabled(enab)
        openButton.setDefault(True)

    self.projectsView.selectionIndexesChanged.connect(checkButtonStates)

    self.layout().addWidget(self.buttonFrame,1,1)

    self.statusWidget = CProjectsStatusBar(self)
    self.layout().addWidget(self.statusWidget,4,0,1,2)

    self.openButton = QtWidgets.QPushButton(self,text='Open')
    self.openButton.setMaximumWidth(100)
    self.openButton.clicked.connect(self.openProject)
    self.layout().addWidget(self.openButton,5,0)

    self.layout().addWidget(self.projectsView,1,0)

    @QtCore.Slot()
    def expandIfFiltered():
        self.projectsView.collapseAll()
        if len(str(self.searchFilterEdit.text()))>0:
            self.projectsView.expandAll()

    searchFilterWidget = QtWidgets.QWidget()
    self.searchFilterEdit = QtWidgets.QLineEdit()
    self.searchFilterEdit.setPlaceholderText("Only show projects with name, annotation or tag containing text typed here")
    self.searchFilterEdit.textChanged.connect(self.projectsViewProxyModel.setFilterString)
    self.searchFilterEdit.textChanged.connect(expandIfFiltered)
    searchFilterWidget.setLayout(QtWidgets.QHBoxLayout())
    searchFilterLabel = QtWidgets.QLabel("Filter:")
    searchFilterWidget.layout().addWidget(searchFilterLabel)
    searchFilterWidget.layout().addWidget(self.searchFilterEdit)
    advButton = QtWidgets.QPushButton()
    if sys.platform == "darwin":
        advButton.setIconSize(QtCore.QSize(12,12))
    upArrow = QtGui.QIcon(str(I2_TOP / 'qticons' / "up.png"))
    downArrow = QtGui.QIcon(str(I2_TOP / 'qticons' / "down.png"))
    advButton.setIcon(downArrow)
    searchFilterWidget.layout().addWidget(advButton)
    searchFilterWidget.layout().setContentsMargins(0,0,0,0)

    advSearchLayout = QtWidgets.QGridLayout()
    dateLabel = QtWidgets.QLabel("Created/active between:")
    dateLabelAnd = QtWidgets.QLabel("and:")
    qdt = QtCore.QDateTime.fromString("1970-01-01 00:00:00","yyyy-MM-dd HH:mm:ss")
    dateStart = QtWidgets.QDateTimeEdit(qdt)
    dateEnd = QtWidgets.QDateTimeEdit(QtCore.QDateTime.currentDateTime())
    advSearchLayout.addWidget(dateLabel,1,0)
    advSearchLayout.addWidget(dateLabelAnd,2,0)
    advSearchLayout.addWidget(dateStart,1,1)
    advSearchLayout.addWidget(dateEnd,2,1)
    dateStartReset = QtWidgets.QPushButton("Reset")
    dateEndReset = QtWidgets.QPushButton("Now")
    dateStart.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Maximum)
    dateEnd.setSizePolicy(QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Maximum)
    advSearchLayout.addWidget(dateStartReset,1,2)
    advSearchLayout.addWidget(dateEndReset,2,2)
    @QtCore.Slot()
    def resetEnd():
        dateEnd.setDateTime(QtCore.QDateTime.currentDateTime())
    dateStartReset.clicked.connect(functools.partial(dateStart.setDateTime,qdt))
    dateEndReset.clicked.connect(resetEnd)

    advWidget = QtWidgets.QWidget()
    advWidget.setLayout(advSearchLayout)
    advWidget.hide()

    @QtCore.Slot()
    def dateChanged():
       startDate = str(dateStart.dateTime().toString("yyyy-MM-dd HH:mm:ss"))
       endDate = str(dateEnd.dateTime().toString("yyyy-MM-dd HH:mm:ss"))
       self.projectsViewProxyModel.setFilterDates(startDate,endDate)
    dateStart.dateTimeChanged.connect(dateChanged)
    dateEnd.dateTimeChanged.connect(dateChanged)

    self.layout().addWidget(searchFilterWidget,2,0)
    self.layout().addWidget(advWidget,3,0)

    @QtCore.Slot()
    def toggleAdvanced():
        if dateLabel.isVisible():
            advWidget.hide()
            advButton.setIcon(downArrow)
        else:
            advWidget.show()
            advButton.setIcon(upArrow)
    advButton.clicked.connect(toggleAdvanced)

    self.exportThread= None

    @QtCore.Slot()
    def handleRepopulateTreeViewNew():
        repopulateTreeViewNew(self.projectsViewProxyModel,self.projectsView)

    PROJECTSMANAGER().db().projectsListChanged.connect(handleRepopulateTreeViewNew)
    PROJECTSMANAGER().db().projectUpdated.connect(handleRepopulateTreeViewNew)
    PROJECTSMANAGER().db().projectCommentEdited.connect(handleRepopulateTreeViewNew)
    PROJECTSMANAGER().db().projectCommentDeleted.connect(handleRepopulateTreeViewNew)

    self.projectsView.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
    self.projectsView.setDragEnabled(True)
    self.projectsView.setAcceptDrops(True)
    self.projectsView.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)

  def setMode(self,mode = 'all'):
    if mode == 'all':
      self.projectsView.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
      self.projectsView.setDragEnabled(True)
      self.projectsView.setAcceptDrops(True)
      self.projectsView.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
      self.buttonFrame.show()
      self.statusWidget.show()
      self.openButton.hide()
    else:    
      self.projectsView.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
      self.projectsView.setDragEnabled(False)
      self.projectsView.setAcceptDrops(False)
      self.buttonFrame.hide()
      self.statusWidget.hide()
      self.openButton.show()

  @QtCore.Slot()
  def addProject(self):
    if self.newProjectGui is None:
      self.newProjectGui = CNewProjectGui(parent=self)
    #self.newProjectGui.start()
    self.newProjectGui.show()

  @QtCore.Slot()
  def handleEditDescription(self):
    projectId=self.projectsView.selectedProjectId()
    if projectId is None: return
    for w in self.findChildren(CProjectDescriptionDialog):
      if w.widget.projectId == projectId:
        w.show()
        w.raise_()
        return
    pName = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectname')
    descriptionDialog = CProjectDescriptionDialog(self,projectId,pName)
    descriptionDialog.show()
    descriptionDialog.raise_()


  def closeEvent(self,event):
    # On closing project manager delete project description children
    #print 'CProjectManagerDialog.closeEvent'
    for w in self.findChildren(CProjectDescriptionDialog):
      w.deleteLater()
    return QtWidgets.QDialog.closeEvent(self,event)


  @QtCore.Slot()
  def handleRenameProject(self):
    projectId=self.projectsView.selectedProjectId()
    if projectId is None: return
    pName = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectname')

    self.renameDialog = QtWidgets.QDialog(self)
    self.renameDialog.setModal(True)
    self.renameDialog.setWindowTitle('Rename project?')
    self.renameDialog.setLayout(QtWidgets.QVBoxLayout())
    self.renameDialog.layout().addWidget(QtWidgets.QLabel('Rename project: '+pName,self))
    self.renameWidget = QtWidgets.QLineEdit(self)
    self.renameWidget.setText(pName)
    self.renameDialog.layout().addWidget(self.renameWidget)
    butBox = QtWidgets.QDialogButtonBox(self)
    but = butBox.addButton('Rename',QtWidgets.QDialogButtonBox.ApplyRole)
    but.setDefault(False)
    but.clicked.connect(functools.partial(self.renameProject,projectId))
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.renameDialog.close)
    self.renameDialog.layout().addWidget(butBox)
    self.renameDialog.show()
    self.renameDialog.raise_()

  @QtCore.Slot(str)
  def renameProject(self,projectId):
    self.renameDialog.hide()   
    #print 'renameProject',projectId,self.renameWidget.text()
    name0 = self.renameWidget.text().__str__()
    if len(name0)==0: return
    name = CCP4Utils.safeOneWord(name0)
    if name0 != name:
      print('Fixing name to ',name)
    try:
      PROJECTSMANAGER().db().getProjectId(name)
    except:
      # Failed finding project of this name - good!
      pass
    else:
      QtWidgets.QMessageBox.warning(self,'Rename project','A project called '+name+' already exists')
      return

    try:
      PROJECTSMANAGER().db().updateProject(projectId,'projectname',name)
    except CException as e:
      warningMessage(e, parent=self,windowTitle=self.windowTitle(),message='Error changing project name')


  @QtCore.Slot()
  def openProject(self):
    projectId=self.projectsView.selectedProjectId()
    if projectId is not None: openProject(projectId)


  def handleDoubleClick(self,item,col=None):
    #print 'handleDoubleCLick',item,col
    nameVar = item.data(0,QtCore.Qt.DisplayRole)
    name = nameVar.__str__()
    openProject(projectName=name)


  @QtCore.Slot()
  def handleExport1(self):
    if getattr(self,'exportThread',None) is not None:
        QtWidgets.QMessageBox.warning(self,'Export project','Please wait - another project export currently in progress')
        return
    projectId=self.projectsView.selectedProjectId()
    if projectId is None:
      #print 'handleExport calling statusBar'
      self.statusWidget.showMessage("Select a project before selecting 'Export'",5)
      return
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)

    if not hasattr(self,'exportOptions'):
      self.exportOptionsWidget = CExportOptionsDialog(self,projectId,projectInfo['projectname'])
      self.exportOptionsWidget.exportSignal.connect(lambda mode: self.handleExport2(projectId,mode))
    self.exportOptionsWidget.show()

  def handleExport3(self,projectId):
    if getattr(self,'exportThread',None) is not None:
        QtWidgets.QMessageBox.warning(self,'Export project','Please wait - another project export currently in progress')
        return
    if projectId is None:
      #print 'handleExport calling statusBar'
      self.statusWidget.showMessage("Select a project before selecting 'Export'",5)
      return
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)

    if not hasattr(self,'exportOptions'):
      self.exportOptionsWidget = CExportOptionsDialog(self,projectId,projectInfo['projectname'])
      self.exportOptionsWidget.exportSignal.connect(lambda mode: self.handleExport2(projectId,mode))
    self.exportOptionsWidget.show()

    progressBar = QtWidgets.QProgressBar()
    progressBar.setWindowTitle("Exporting "+str(projectInfo['projectname']))
    progWin = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout()
    progWin.setLayout(layout)
    layout.addWidget(progressBar)
    label = QtWidgets.QLabel("Exporting "+str(projectInfo['projectname']))
    layout.addWidget(label)

    @QtCore.Slot(tuple)
    def setLabelProgress(ret):
        job,done = ret
        label.setText("Saving job "+str(job))

    @QtCore.Slot()
    def showProgressBar():
        progressBar.setMaximum(0)
        progressBar.setMinimum(0)
        progWin.show()
        if hasattr(self,"exportThread"):
           self.exportThread.savingJobData.connect(setLabelProgress)
    @QtCore.Slot()
    def hideProgressBar():
        progWin.close()
    self.startSavingJobDataSignal.connect(showProgressBar)
    self.doneSavingJobDataSignal.connect(hideProgressBar)

  @QtCore.Slot(str,int)
  def handleExport2(self,projectId,mode):
    #print('handleExport2',projectId,mode)
    if mode == 0: return
    if mode == 1:
      after = None
      jobList = None
    if mode == 2:
      after = self.exportOptionsWidget.impExpWidget.getTime()
      jobList = None
    elif mode == 3:
      after = None
      jobList,errList = self.exportOptionsWidget.selectionWidget.getSelection()
      if len(errList)>0:
        errText = 'Unable to interpret selection: '
        for err in errList: errText = errText +"'"+err+"' "
        mess = QtWidgets.QMessageBox.warning(self,'Export project',errText)
        return

    excludeI2files = False

    #print 'handleExport2',after,jobList
    self.exportOptionsWidget.close()
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)

    safeName = CCP4Utils.safeOneWord(projectInfo['projectname'])

    from . import CCP4FileBrowser
    from ..qtcore import CCP4Export
    self.browser = CCP4FileBrowser.CFileDialog(self,
           title='Save project data to compressed file',
           filters= ['CCP4 Project Database (*.'+CCP4Export.COMPRESSED_SUFFIX+')'],
           defaultSuffix=CCP4Export.COMPRESSED_SUFFIX,
           defaultFileName=safeName,
           fileMode=QtWidgets.QFileDialog.AnyFile  )
    self.browser.selectFile.connect(functools.partial(self.compressProject,projectId,after,jobList,excludeI2files))
    self.browser.show()

  def export(self,projectId,fileName):
    jobNumberList,errReport = PROJECTSMANAGER().db().exportProjectXml(projectId,fileName)
    if len(errReport)>0:
       warningMessage(errReport, parent=self,windowTitle='Errors exporting project',message='Some errors in exporting project to file:\n'+fileName)


  @QtCore.Slot(str,float,list,list,str)
  def compressProject(self,projectId,after=None,jobList=None,excludeI2files=False,fileName=None):
    #print 'CProjectManagerDialog.compressProject',after,jobList,excludeI2files,fileName
    self.browser.hide()
    self.browser.deleteLater()

    projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)
    # If there is a limited set of jobs then find the input jobs that are not output by jobs on that list
    inputFilesList,inputFileIdList,fromJobList,errReport =  PROJECTSMANAGER().getJobInputFiles(projectDir=projectInfo['projectdirectory'],jobIdList=jobList,useDb=True,excludeI2files=excludeI2files)
    #print 'CProjectManagerDialog.compressProject inputFilesList,fromJobIdList',inputFilesList,fromJobIdList
    fromJobIdList = []
    fromJobNumberList = []
    for item in fromJobList:
      fromJobIdList.append(item[0])
      fromJobNumberList.append(item[1])

    dbxml = os.path.join( projectInfo['projectdirectory'],'CCP4_TMP','DATABASE'+str(int(time.time()))+'.db.xml')
    self.statusWidget.showMessage('Creating XML database:'+dbxml)
    #print 'Creating temporary database xml file in:',dbxml
    # exportProjectXml returns list of TOP-LEVEL jobNumbers for the export
    jobNumberList,errReport = PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml,recordExport=True,status='exportable',after=after,jobList=jobList,inputFileList=inputFileIdList,inputFileFromJobList=fromJobIdList)
    #print 'CProjectManagerDialog.compressProject jobNumberList',jobNumberList
    if errReport.maxSeverity()>Severity.WARNING:
        self.statusWidget.clearMessage()
        warningMessage(errReport, 'Export project','Error creating XML database file',parent=self)
        return

    directoriesList = ['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES']
    from ..qtcore import CCP4Export
    self.exportThread = CCP4Export.ExportProjectThread(self,projectDir=projectInfo['projectdirectory'],dbxml=dbxml,target=fileName,jobList=jobNumberList,inputFilesList=inputFilesList,directoriesList=directoriesList,)
    self.exportThread.savingJobData.connect(self.updateSavingJobData)
    self.exportThread.startSavingJobData.connect(self.progressSavingJobData)
    self.exportThread.finished.connect(functools.partial(self.doneSavingJobData,projectInfo['projectname'],fileName))
    self.startSavingJobDataSignal.emit()
    self.exportThread.start()

  @QtCore.Slot(int)
  def progressSavingJobData(self,numberOfJobs):
    self.statusWidget.showMessage('Exporting project: Saving '+str(numberOfJobs)+' jobs to compressed file')
    self.statusWidget.showProgress(numberOfJobs)

  @QtCore.Slot(tuple)
  def updateSavingJobData(self,ret):
    jobNumber,jobsDone = ret    
    #print 'updateSavingJobData',jobNumber,jobsDone
    if jobNumber in ['IMPORT','DATABASE','PROJECT']:
      label = ['imported files','database file','project data'][ ['IMPORT','DATABASE','PROJECT'].index(jobNumber)]
      self.statusWidget.clear()
      if not jobsDone:
        self.statusWidget.showMessage('Exporting project: Saving '+label)
      else:
        self.statusWidget.showMessage('Exporting project: Finished saving '+label)
    else:
      self.statusWidget.setProgress(jobsDone)

  @QtCore.Slot(str,str)
  def doneSavingJobData(self,projectName,filename):
    report = self.exportThread.errorReport
    if report.maxSeverity()>Severity.WARNING:
      warningMessage(report, 'Saving job data','Error saving data files for export',parent=self)
    else:
      self.statusWidget.hideProgress()
      self.statusWidget.showMessage('Project '+str(projectName)+' saved to: '+str(filename))
    self.exportThread.deleteLater()
    self.exportThread = None
    self.doneSavingJobDataSignal.emit()

  @QtCore.Slot()
  def handleExportAll(self):
      projects = PROJECTSMANAGER().db().getProjectDirectoryList()

      cwd = os.getcwd()

      rv = QtWidgets.QFileDialog.getExistingDirectory(caption='Select directory for saving all projects')
      if not rv:
          return

      origLenProjects = len(projects)

      progressWindow = QtWidgets.QWidget()
      progressLabel = QtWidgets.QLabel("Exporting "+str(len(projects))+" projects")
      progressBar = QtWidgets.QProgressBar()
      progressBar.setValue(0)
      progressLayout = QtWidgets.QVBoxLayout()
      progressLayout.addWidget(progressLabel)
      progressLayout.addWidget(progressBar)
      progressDBB = QtWidgets.QDialogButtonBox()
      cancelButton = progressDBB.addButton(QtWidgets.QDialogButtonBox.Cancel)
      doneButton = progressDBB.addButton(QtWidgets.QDialogButtonBox.Ok)
      doneButton.setEnabled(False)
      progressLayout.addWidget(progressDBB)
      progressWindow.setLayout(progressLayout)
      progressWindow.show()
      progressWindow.raise_()


      @QtCore.Slot(bool)
      def doneClicked(dum=False):
          progressWindow.close()
          return
      doneButton.clicked.connect(doneClicked)

      self.cancelExportAll = False

      @QtCore.Slot(bool)
      def cancelClicked(dum=False):
          self.cancelExportAll = True
          cancelButton.setEnabled(False)
          progressLabel.setText("Cancelling...")
      cancelButton.clicked.connect(cancelClicked)

      iproject = 0

      @QtCore.Slot(int)
      def exportNextProject(iproject):
          try:
              if self.cancelExportAll:
                  progressWindow.close()
                  return
              name = projects.pop()[1]
              pc = int(100. * iproject / origLenProjects)
              progressBar.setValue(pc)
              progressLabel.setText("Exporting project "+name+" "+str(iproject+1)+" of "+str(origLenProjects))
              from ..utils import ExportAllProjects
              compClass = ExportAllProjects.CompressClass(projectName=name,fileName=os.path.join(rv,name+".ccp4_project.zip"))
              compClass.doneSignal.connect(functools.partial(exportNextProject,iproject+1))
              compClass.run()
          except:
              progressLabel.setText("Finished exporting "+str(origLenProjects)+" projects")
              progressBar.setValue(100)
              #progressWindow.close()
              doneButton.setEnabled(True)
              return

      exportNextProject(iproject)

  @QtCore.Slot()
  def handleImportProject(self):
    from ..dbapi import CCP4DbApi
    CCP4DbApi.CDbXml.updateInstances()
    if len(CCP4DbApi.CDbXml.Instances)>0:
      QtWidgets.QMessageBox.warning(self,'Import project','Another import is already in progress - please wait')
      return

    from . import CCP4FileBrowser
    self.browser = CCP4FileBrowser.CFileDialog(self,
           title='Import project compressed file',
           filters = [ MIMETYPESHANDLER().getMimeTypeInfo("application/CCP4-compressed-db",'filter') ],
           defaultSuffix=  MIMETYPESHANDLER().getMimeTypeInfo("application/CCP4-compressed-db",'fileExtensions')[0],
           fileMode=QtWidgets.QFileDialog.ExistingFile  )
    self.browser.selectFile.connect(self.handleImportProject1)
    self.browser.show()

  @QtCore.Slot(str)
  def handleImportProject1(self,compressedFile):

    if hasattr(self,'browser'):
      try:
        self.browser.hide()
        self.browser.deleteLater()
        del self.browser
      except:
        pass

    self.statusWidget.showMessage('Checking contents of compressed file')

    try:
      xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedFile)
    except CException as e:
      warningMessage(e, 'Import project','Failed extracting database XML file from compressed file',parent=self)
      self.statusWidget.clear()
    except Exception as e:
      QtWidgets.QMessageBox.warning(self,'Import project','Error extracting database xml file from '+str(compressedFile))
      self.statusWidget.clear()
      return

    try:
      from ..dbapi import CCP4DbApi
      self.dbImport = CCP4DbApi.CDbXml(db=PROJECTSMANAGER().db(),xmlFile=xmlFile)
      #self.dbImport.setDiagnostic(True)
      importProjectInfo = self.dbImport.loadProjectInfo()  
    except:
      QtWidgets.QMessageBox.warning(self,'Import project','Error attempting to read database file in\n'+str(compressedFile))
      return
    try:
      projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=self.dbImport.projectId)
    except:
      projectInfo = None

    #print 'handleImportProject1 importProjectInfo',importProjectInfo

    if projectInfo is None:
      # There is no matching projectId in db- query user to create new project
      if not hasattr(self,'importNewDialog'): self.importNewDialog = CImportNewProjectDialog(self)
      self.importNewDialog.reset(importProjectInfo=self.dbImport.headerInfo())
      self.importNewDialog.importSignal.connect(functools.partial(self.importProject0,compressedFile))
      self.importNewDialog.show()

    else:

      self.importExistingDialog=CImportExistingProjectDialog(self,projectName=projectInfo['projectname'])
      self.importExistingDialog.reset(importProjectInfo=self.dbImport.headerInfo())
      self.importExistingDialog.importSignal.connect(functools.partial(self.importProject1,compressedFile,projectInfo['projectdirectory'],self.dbImport.projectId))
      self.importExistingDialog.show()     

  @QtCore.Slot(str)
  def importProject0(self,compressedFile):
    # Close importNewDialog
    if not hasattr(self,'importNewDialog'): return
    self.importNewDialog.close()

    if self.importNewDialog.mode() == 1:
      # Create a new project
      dirName = self.importNewDialog.newDirectory()
      projectName = self.importNewDialog.projectName()
      self.importNewDialog.deleteLater()
      del self.importNewDialog
      self.importProject(compressedFile,dirName,projectName=projectName)
    else:
      # Merge into existing project despite incompatible projectId
      projectId = self.importNewDialog.mergeProject()
      self.importNewDialog.deleteLater()
      del self.importNewDialog
      projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)
      self.dbImport.projectId = projectId
      self.dbImport.projectDirectory = projectInfo['projectdirectory']
      #self.dbImport.projectName = self.importNewDialog.projectName()
      self.importProject(compressedFile=compressedFile,dirName=self.dbImport.projectDirectory,existingProject=projectId,forceProjectId=True)

  @QtCore.Slot(str,str,str)
  def importProject1(self,compressedFile,dirName,existingProject=None):
    self.importExistingDialog.close()
    self.importExistingDialog.deleteLater()
    del self.importExistingDialog

    # Append to existing project
    self.importProject(compressedFile=compressedFile,dirName=dirName,existingProject=existingProject)

  def importProject(self,compressedFile,dirName,existingProject=None,forceProjectId=False,projectName=None):
    if DIAGNOSTIC:
        print('CProjectManagerDialog.importProject',compressedFile,dirName,existingProject,projectName)
    # Load the database.xml into temporary tables in db
    self.dbImport.projectDirectory = dirName
    if projectName is not None: self.dbImport.projectName = projectName
    if existingProject is None:
      ret = self.dbImport.createProject()
      if DIAGNOSTIC:
          print('dbImport.createProject()',ret.report())
      if ret.maxSeverity()>Severity.WARNING:
        print(ret.report())
        warningMessage(ret, parent=self,windowTitle='Error creating project in database',ifStack=False)
        return


    self.dbImport.createTempTables()
    if forceProjectId:
      self.dbImport.loadTempTable(resetProjectId=existingProject)
    else:
      self.dbImport.loadTempTable()
    # If loading jobs to an existing project flag up jobs in temp tables that
    # are already in db
    if existingProject is not None:
      self.dbImport.setExclInTempTables()
    # Flag imported files to be imported (there is no checking yet that they exist)
    self.dbImport.setExclImportedFiles()

    if DIAGNOSTIC: print('CProjectManagerDialog.importProject setting Temp Tables',self.dbImport.errReport.report())

    if self.dbImport.errReport.maxSeverity()>Severity.WARNING:
      if DIAGNOSTIC: print('Error report from the import process..')
      if DIAGNOSTIC: print(self.dbImport.errReport.report())
      warningMessage(self.dbImport.errReport, parent=self,windowTitle='Error loading data from project export file',ifStack=False)

    # Make project directory if necessary
    if not os.path.exists(dirName):
      try:
        os.mkdir(dirName)
      except:
        QtWidgets.QMessageBox.warning(self,'Import project','Failed to create directory:'+dirName)
        return

    self.statusWidget.showMessage('Unpacking project files to '+dirName)
    from ..qtcore import CCP4Export
    # Unpack project files from the tar file (possibly in separate thread) 
    # Pass import thread dbImport to enable query database and flagging loaded jobs/files
    if DIAGNOSTIC: print('CProjectManagerDialog.importProject creating import thread')
    self.importThread = CCP4Export.ImportProjectThread(self,projectDir=dirName,compressedFile=compressedFile,
                                                       dbImport=self.dbImport,diagnostic=DIAGNOSTIC)
    #self.importThread.finished.connect(self.doneImportProgress)
    errReport = self.importThread.run()
    if DIAGNOSTIC: print('CProjectManagerDialog.importProject import thread running',errReport.report())

    self.dbImport.cleanupTempTables()
    #self.dbImport.listTempJobs('TempJobs after cleanup')
    #self.dbImport.listTempFiles('TempFiles after cleanup')
    #self.dbImport.listTempFileUses('TempFileUses after cleanup')
    stats = self.dbImport.importStats()
    if DIAGNOSTIC:
      for key,value in list(stats.items()):
        if key == 'failedFiles':
          if len(value)>0: print('Failed to import files..(probably already present)')
          for item in value:
              print('Job_'+str(item[4]),item[2])
        else:
          print('CProjectManagerDialog.importProject stats', key,value)
    if errReport.maxSeverity()>Severity.WARNING:
      self.dbImport.removeTempTables()
      text = 'ERRORS UNPACKING DATA FILES\n'
      for err in errReport: text = text + err['details'] + '\n'
      QtWidgets.QMessageBox.warning(self,'Import failed',text)
      self.statusWidget.showMessage('Error unpacking project files to '+dirName)
      return

    self.statusWidget.showMessage('Loading project data to database')
    self.dbImport.importTempTables()
    #self.dbImport.importProjectCommentsTempTables()
    self.dbImport.removeTempTables()
    self.statusWidget.showMessage('Finished importing project')

    self.dbImport.db.projectReset.emit({'projectId':self.dbImport.projectId})

    if stats['jobsTotal']>stats['jobsImported'] or stats['filesTotal']>stats['filesImported']:
      text = 'Some of the jobs/files are already in the database in project '+str(self.dbImport.projectName)+'\n' + \
      'Imported '+str(stats['jobsImported'])+' new jobs from '+str(stats['jobsTotal'])+' in file \n' + \
     'and '+str(stats['filesImported'])+' new data files from '+str(stats['filesTotal'])+' in file\n'
    else:
      text = 'Successfully imported '+str(stats['jobsImported'])+' jobs and '+str(stats['filesImported'])+' data files'

    if stats.get('incrJobNumber',0) > 0:
        text = text +'\nImporting jobs '+str(stats['importMin'])+' to '+str(stats['importMax'])+' have been renumbered\n'+str(int(stats['importMin'])+int(stats['incrJobNumber']))+' to '+str(int(stats['importMax'])+int(stats['incrJobNumber'])) +' to avoid clash with existing jobs'
    if len(text)>0:  QtWidgets.QMessageBox.information(self,'Import complete',text)
    PROJECTSMANAGER().backupDBXML()
    projectId = copy.deepcopy(self.dbImport.projectId)
    openProject(projectId=projectId)

  def createImportProgress(self):
    pass

  def handleImportProgress(self,ret):
    #print 'handleImportProgress',ret
    jobNo,done = ret

  @QtCore.Slot()
  def doneImportProgress(self):
    pass

  @QtCore.Slot()
  def handleDeleteProject(self):
    projectId=self.projectsView.selectedProjectId()
    if projectId is None:
        self.statusWidget.showMessage("Select a project before selecting 'Delete'",5)
        return
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)

    self.deleteDialog = QtWidgets.QDialog(self)
    #self.deleteDialog.setModal(True)
    self.deleteDialog.setWindowTitle('Delete project?')
    self.deleteDialog.setLayout(QtWidgets.QVBoxLayout())
    self.deleteDialog.layout().addWidget(QtWidgets.QLabel('Delete project from database: '+projectInfo['projectname'],self))
    self.deleteDirectory = QtWidgets.QCheckBox('Delete directory: '+projectInfo['projectdirectory'],self)
    self.deleteDialog.layout().addWidget(self.deleteDirectory)
    butBox = QtWidgets.QDialogButtonBox(self)
    but = butBox.addButton('Delete project',QtWidgets.QDialogButtonBox.ApplyRole)
    but.setDefault(False)
    but.clicked.connect(functools.partial(self.deleteProject,projectId))
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.deleteDialog.close)
    self.deleteDialog.layout().addWidget(butBox)
    self.deleteDialog.show()
    self.deleteDialog.raise_()

  @QtCore.Slot(str)
  def deleteProject(self,projectId):
    deleteDirectory = self.deleteDirectory.isChecked()
    self.deleteDialog.close()
    if projectId is None: return
    e = PROJECTSMANAGER().deleteProject(projectId=projectId,deleteDirectory=deleteDirectory)
    if e.maxSeverity()>Severity.WARNING:
       warningMessage(e, 'Deleting project',parent=self)
    PROJECTSMANAGER().backupDBXML()

  @QtCore.Slot()
  def handleMoveProject(self):
    projectId=self.projectsView.selectedProjectId()
    #print 'handleMoveProject projectId',projectId,type(projectId)
    if projectId is None: return
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)

    self.moveProjectId = projectId
    self.moveDialog = QtWidgets.QDialog(self)
    self.moveDialog.setModal(True)
    self.moveDialog.setWindowTitle('Move project directory?')
    self.moveDialog.setLayout(QtWidgets.QVBoxLayout())
    self.moveDialog.layout().addWidget(QtWidgets.QLabel('Move the project directory for project: '+projectInfo['projectname'],self))
    self.moveDialog.layout().addWidget(QtWidgets.QLabel('Current directory: '+projectInfo['projectdirectory'],self))
    # Mode button group
    self.moveMode = -1
    self.moveModeGroup = QtWidgets.QButtonGroup(self)
    self.moveModeGroup.setExclusive(True)
    but =  QtWidgets.QRadioButton('Move the directory and register with database',self)
    self.moveModeGroup.addButton(but,1)
    self.moveDialog.layout().addWidget(but)
    but =  QtWidgets.QRadioButton('Directory is already moved - just register with database',self)
    self.moveModeGroup.addButton(but,0)
    self.moveDialog.layout().addWidget(but)
    self.moveModeGroup.buttonReleased[int].connect(self.handleMoveModeChange)

    from ..core import CCP4DataManager
    from ..core import CCP4File
    self.moveStack = QtWidgets.QStackedLayout()
    self.moveDialog.layout().addLayout(self.moveStack)

    self.moveStack.addWidget(QtWidgets.QLabel('Choose mode above',self))


    self.moveDirectory0 =  CCP4File.CDataFile(parent=self,qualifiers= {  'isDirectory' : True, 'mustExist' : True } )
    self.moveDirectory0Widget = CCP4DataManager.DATAMANAGER().widget(model=self.moveDirectory0,parentWidget=self,qualifiers={'jobCombo': False, 'projectBrowser' : False } )
    self.moveDirectory0Widget.updateViewFromModel()
    self.moveDirectory0Widget.fileLineEdit.setCharWidth(60,mode='minimum')
    self.moveStack.addWidget(self.moveDirectory0Widget)

    self.moveDirectory1 =  CCP4File.CDataFile(parent=self,qualifiers= { 'isDirectory' : True,'mustExist' : False } )
    self.moveDirectory1Widget = CCP4DataManager.DATAMANAGER().widget(model=self.moveDirectory1,parentWidget=self,qualifiers={'jobCombo': False, 'projectBrowser' : False} )
    self.moveDirectory1Widget.updateViewFromModel()
    self.moveDirectory1Widget.fileLineEdit.setCharWidth(60,mode='minimum')
    self.moveStack.addWidget(self.moveDirectory1Widget)

    butBox = QtWidgets.QDialogButtonBox(self)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Apply)
    but.setDefault(False)
    but.clicked.connect(self.moveProject)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.moveDialog.close)
    self.moveDialog.layout().addWidget(butBox)

    self.moveStack.setCurrentIndex(self.moveMode -1 )
    self.moveDialog.show()
    self.moveDialog.raise_()

  @QtCore.Slot(int)
  def handleMoveModeChange(self,mode):
    mode = int(mode)
    if mode != self.moveMode:
      self.moveMode = mode
      self.moveDirectory0Widget.closeBrowser()
      self.moveDirectory1Widget.closeBrowser()
      self.moveStack.setCurrentIndex(self.moveMode + 1 )
      #print 'handleMoveModeChange',self.moveMode,self.moveStack.currentIndex()

  @QtCore.Slot()
  def moveProject(self):
    errMess = None
    if self.moveMode < 0:
        errMess = 'Choose if moving the directory or just updating database'
    if self.moveMode == 0:
      if not self.moveDirectory0.isSet():
        errMess = 'Choose directory'
      else:
        newDir = str(self.moveDirectory0)
        if not os.path.exists(newDir) or not os.path.isdir(newDir):
          errMess = 'Directory does not exist or is not directory'
    else:
      if not self.moveDirectory1.isSet():
        errMess = 'Choose new directory path'
      else:
        newDir = str(self.moveDirectory1)
        if os.path.exists(newDir):
          errMess = 'New directory should not exist already'
        #else:
        #  if not os.path.exists(os.path.split(newDir)[0]):
        #    errMess = 'Parent directory for new directory must exist'
    if errMess is not None:
       QtWidgets.QMessageBox.warning(self,'Can not move directory',errMess)
       return

    #print 'moveProject',newDir,self.moveMode
    projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=self.moveProjectId)
    if self.moveMode == 1:
      try:
        subDirList = glob.glob(os.path.join(projectInfo['projectdirectory'], 'CCP4_*'))
        for subDir in subDirList:
          shutil.move(subDir, newDir)
        if len(glob.glob(os.path.join(projectInfo['projectdirectory'], '*'))) > 0:
          errMess = 'Project directory contains files/directories not created by CCP4i2. These have not been moved'
          QtWidgets.QMessageBox.warning(self,'Warning unexpected files or directories in project directory',errMess)

      except Exception as e:
         QtWidgets.QMessageBox.warning(self,'Error moving directory',str(e))
         return
    try:
      PROJECTSMANAGER().db().updateProject(self.moveProjectId,key='projectdirectory',value=newDir)
    except CException as e:
      warningMessage(e, windowTitle='Error updating database',parent=self)
      return
    except Exception as e:
      QtWidgets.QMessageBox.warning(self,'Error updating database',str(e))
      return

    self.moveDialog.close()

  @QtCore.Slot()
  def handleCleanup(self):
    projectId=self.projectsView.selectedProjectId()
    #print 'handleCleanup projectId',projectId
    if projectId is None:
      ret = QtWidgets.QMessageBox.question(self,self.windowTitle(),'Remove temporary files from all projects?',QtWidgets.QMessageBox.Ok|QtWidgets.QMessageBox.Cancel)
      if ret != QtWidgets.QMessageBox.Ok: return
      PROJECTSMANAGER().cleanupAllProjects(context='temporary')
    else:
      msgBox = QtWidgets.QMessageBox()
      msgBox.setWindowTitle('Cleanup files')
      msgBox.setText('Do you want to delete only temporary files or temporary files and intermediate data files from sub-jobs')
      b = msgBox.addButton(QtWidgets.QMessageBox.Cancel)
      b = msgBox.addButton('Temporary files only',QtWidgets.QMessageBox.ApplyRole)
      b.clicked.connect(functools.partial(self.handleCleanup2,projectId,'temporary'))
      b = msgBox.addButton('Temporary and intermediate files',QtWidgets.QMessageBox.ApplyRole)
      b.clicked.connect(functools.partial(self.handleCleanup2,projectId,'intermediate'))
      msgBox.exec_()

  @QtCore.Slot(str,str)
  def handleCleanup2(self,projectId,context):
      #print 'handleCleanup2',projectId,context
      from ..core import CCP4ProjectsManager
      cleanup = CCP4ProjectsManager.CPurgeProject(projectId = projectId)
      cleanup.purgeProject(context=context)

  @QtCore.Slot()
  def handleRerun(self):
    projectId=self.projectsView.selectedProjectId()
    if projectId is None: return
    self.rerunProjectId = projectId
    filter_list = []

    from . import CCP4FileBrowser
    self.rerunDialog = CCP4FileBrowser.CFileDialog(self,fileMode=QtWidgets.QFileDialog.Directory,
      title='Directory for rerun output')
    self.rerunDialog.selectFile.connect(self.rerunProject)
    self.rerunDialog.show()

  @QtCore.Slot(str)
  def rerunProject(self,outputDirectory):
    self.rerunDialog.close()
    if not os.path.exists(outputDirectory) or not os.path.isdir(outputDirectory):
      QtWidgets.QMessageBox.warning(self,'No suitable directory for project rerun','You must provide a directory in which the rerun project directory will be created')
      return
    from ..core import CCP4ProjectBasedTesting
    self.projectBasedTest = CCP4ProjectBasedTesting.CProjectBasedTesting(sourceProjectList=[self.rerunProjectId],outputDirectory=outputDirectory,useCurrentDb=True)
    #self.projectBasedTest.start()
    self.projectBasedTest.runTests()
    logFile = self.projectBasedTest.logFiles[-1]
    print('Re-running project log file:',logFile)
    widget = WEBBROWSER().openFile(logFile)
    self.projectBasedTest.reportUpdated.connect(WEBBROWSER().reloadFile)

  def recoverProject(self,projectDir):
    #print 'recoverProject',projectDir
    if not os.path.exists(projectDir) or not os.path.isdir(projectDir):
      QtWidgets.QMessageBox.warning(self,'No project directory selected','You must select a project directory')
      return
    if not os.path.exists(os.path.join(projectDir,'CCP4_JOBS')):
      QtWidgets.QMessageBox.warning(self,'Not a project directory?','This does not appear to be a CCP4i2 project directory which should contain a CCP4_JOBS sub-directory')
      return

    from ..dbapi import CCP4DbUtils
    self.makeDbXml = CCP4DbUtils.CMakeProjectDbXml(self,projectDir=projectDir,projectName=os.path.split(projectDir)[-1])
    self.makeDbXml.jobLoaded.connect(self.statusWidget.setProgress)

    self.statusWidget.showMessage('Making database file from project directory')
    self.statusWidget.showProgress(self.makeDbXml.numberOfJobs())

    errReport = self.makeDbXml.loadProject()

    if errReport.maxSeverity()>Severity.WARNING:
      warningMessage(errReport, parent=self,windowTitle='Error recovering project directory',
                               message='Some errors are reported',ifStack=False,minSeverity=Severity.ERROR)
      return

    xmlFile = self.makeDbXml.saveXmlFile()
    if len(errReport)>0:
       warningMessage(errReport, parent=self,windowTitle='Project database information recovered',
         message='The project database information has been recovered and saved to '+xmlFile+ \
        '\nThere are some warnings which can probably be ignored.\nTo seee click Show Details' ,ifStack=False)

    self.statusWidget.showMessage('Loading database file to database')
    self.importXmlDatabase(xmlFile,projectDir=projectDir)
    self.statusWidget.clear()
    self.statusWidget.removeProgress()

  def importXmlDatabase(self,xmlFile,projectDir=None):
    from ..dbapi import CCP4DbApi
    self.dbImport = CCP4DbApi.CDbXml(db=PROJECTSMANAGER().db(),xmlFile=xmlFile)
    projectInfo = self.dbImport.loadProjectInfo()
    # Check that we have a projectId and it is not already in DB
    #print 'importXmlDatabase projectInfo',projectInfo
    if self.dbImport.projectId is None:
      print('Imported project XML file does not appear to have a projectId')
      return
    try:
      dbProjectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=self.dbImport.projectId,mode=['projectname','projectdirectory'])
    except:
      pass
    else:
      QtWidgets.QMessageBox.warning(self,'Error importing database XML','There is aready a project with the same projectID and project name '+str(dbProjectInfo.get('projectname','unknown'))+' and it has project directory '+str(dbProjectInfo.get('projectdirectory','unknown')))
      return
    if projectDir is not None: self.dbImport.projectDirectory = projectDir


    ret = self.dbImport.createProject()
    # Expect an error from createProject if the project is already in db
    #print 'CProjectManagerDialog.importProject from createProject',ret.report()
    self.dbImport.createTempTables()
    self.dbImport.loadTempTable()

    if self.dbImport.errReport.maxSeverity()>Severity.WARNING:
      print('Error report from the import process..')
      warningMessage(self.dbImport.errReport, windowTitle='Project manager- importing XML',
                                message='Failed importing database from XML file',ifStack=False,parent=self)
    else:
      self.dbImport.setAllToImport()
      #self.dbImport.listTempJobs('Temp jobs')
      #self.dbImport.listTempFiles('Temp files')
      self.dbImport.importTempTables()
      self.dbImport.removeTempTables()

      openProject(projectId=self.dbImport.projectId)

class CTagLineEdit( QtWidgets.QLineEdit ):

    enterKeyPress = QtCore.Signal('QKeyEvent')

    def keyPressEvent(self,event):
      if event.key() in [QtCore.Qt.Key_Return,QtCore.Qt.Key_Enter]:
        self.enterKeyPress.emit(event)
        event.accept()
      else:
        QtWidgets.QLineEdit.keyPressEvent(self,event)

class CProjectDescription(QtWidgets.QFrame):
  TAGSPERROW = 3
  TAGLABEL = 'Choose tag..'

  def __init__(self,parent,projectId,projectName=''):
    QtWidgets.QFrame.__init__(self,parent)
    self.projectId = projectId
    self.dbTags = PROJECTSMANAGER().db().getTagList()
    #print 'CProjectDescription self.dbTags',self.dbTags

    self.setLayout(QtWidgets.QGridLayout())
    MARGIN = 1
    self.layout().setSpacing(MARGIN)
    self.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
    self.layout().addWidget(QtWidgets.QLabel('Description of project',self),0,0,1,1)
    self.annotationWidget = QtWidgets.QTextEdit(self)
    self.annotationWidget.setToolTip('Enter description of project')
    self.layout().addWidget(self.annotationWidget,1,0,4,self.TAGSPERROW)

    #self.layout().addWidget(QtWidgets.QLabel('Tag the project..',self),5,0,1,self.TAGSPERROW)
    from . import CCP4GuiUtils
    icon = CCP4GuiUtils.createIcon('list_add_grey')
    self.moreTagsBut = QtWidgets.QToolButton(self)
    self.moreTagsBut.setIcon(icon)
    #self.moreTagsBut = QtWidgets.QPushButton('Add row of tags',self)
    self.layout().addWidget(self.moreTagsBut,5,0)
    self.moreTagsBut.clicked.connect(self.drawTags)

    # Edit frame
    self.newTagFrame = QtWidgets.QFrame(self)
    self.newTagFrame.setLayout(QtWidgets.QHBoxLayout())
    self.newTagFrame.layout().setSpacing(MARGIN)
    self.newTagFrame.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
    self.editTagLabel = QtWidgets.QLabel('New tag',self)
    self.newTagFrame.layout().addWidget(self.editTagLabel)
    self.newTagWidget = CTagLineEdit(self)
    self.newTagWidget.setToolTip('Enter new tag name')
    self.newTagFrame.layout().addWidget(self.newTagWidget)

    # Save/delete 
    newTagBut = QtWidgets.QPushButton('Save',self)
    self.newTagFrame.layout().addWidget(newTagBut)


    # Tool button
    icon =CCP4GuiUtils.createIcon(name='gears2')    
    but = QtWidgets.QToolButton(self)
    but.setIcon(icon)
    but.setToolTip('Tag handling tools')
    but.setPopupMode(QtWidgets.QToolButton.InstantPopup)
    but.setMenu(QtWidgets.QMenu(self))
    but.menu().addAction('Edit tag',self.editTag)
    but.menu().addAction('Delete tag',self.deleteTag)
    but.menu().addAction('Delete unused tags',self.deleteUnusedTags)
    self.newTagFrame.layout().addWidget(but)

    self.newTagWidget.enterKeyPress.connect(self.saveNewTag)
    newTagBut.clicked.connect(self.saveNewTag)
    self.layout().addWidget(self.newTagFrame,6,0,1,self.TAGSPERROW)

    PROJECTSMANAGER().db().tagCreated.connect(functools.partial(self.handleTagEditSignal,'tagCreated'))
    PROJECTSMANAGER().db().tagUpdated.connect(functools.partial(self.handleTagEditSignal,'tagUpdated'))
    PROJECTSMANAGER().db().tagDeleted.connect(functools.partial(self.handleTagEditSignal,'tagDeleted'))
    PROJECTSMANAGER().db().unusedTagsDeleted.connect(functools.partial(self.handleTagEditSignal,'unusedTagsDeleted'))

    self.editTagCombo = None
    self.deleteTagCombo = None
    self.blockTagWidgetChanged = False
    self.nTagRows = 0
    self.tagWidgets = []
    self.commentId = None
    self.drawTags()
    self.load()

  @QtCore.Slot()
  def saveNewTag(self):
    text = str(self.newTagWidget.text()).strip()
    #print 'saveNewTag',text
    self.newTagWidget.clear()
    if len(text)==0:
      return
    try:
      tagId = PROJECTSMANAGER().db().newTag(self.projectId,text=text)
    except CException as e:
      warningMessage(e, 'Save project tag','Failed saving new project tag',self)
      return
    except Exception as e:
      print('Failed creating tag',e)
      return

  @QtCore.Slot(str,dict)
  def handleTagEditSignal(self,mode,args={}):
    #print 'handleTagEditSignal',mode,args
    self.blockTagWidgetChanged = True
    if mode == 'tagCreated':
      tagId = args['tagId']
      text = args['text']
      # A tagsCreated signal so can just add new tag to combo boxes
      self.dbTags.append([tagId,None,text])
      # Load new tag to tagWidgets combos and set the first unset combo
      # to the new tag
      for tagWidget in self.tagWidgets:
        tagWidget.addItem(text,tagId)
      if self.editTagCombo is not None: self.editTagCombo.addItem(text,tagId)
      if self.deleteTagCombo is not None: self.deleteTagCombo.addItem(text,tagId)
      loaded = False
      for tagWidget in self.tagWidgets:
        if tagWidget.currentIndex() == 0:
          tagWidget.setCurrentIndex(tagWidget.count()-1)
          loaded = True
          break
      if not loaded:
        # No free tag combos - create another row and set first in the row
        self.drawTags()
        self.tagWidgets[-self.TAGSPERROW].setCurrentIndex(tagWidget.count()-1)
    elif mode in [ 'tagDeleted', 'tagUpdated','unusedTagsDeleted']:
      self.dbTags = PROJECTSMANAGER().db().getTagList()
      comboList = []
      comboList.extend(self.tagWidgets)
      if self.editTagCombo is not None:comboList.append(self.editTagCombo)
      if self.deleteTagCombo is not None:comboList.append(self.deleteTagCombo)
      for tagWidget in self.tagWidgets:
        currentTagId = tagWidget.itemData(tagWidget.currentIndex(),QtCore.Qt.UserRole).__str__()
        tagWidget.clear()
        tagWidget.addItem('Choose tag..',0)
        for tid,parentTid,text in self.dbTags:
          tagWidget.addItem(text,tid)
        tagWidget.setCurrentIndex(max(0,tagWidget.findData(currentTagId,QtCore.Qt.UserRole)))
    self.blockTagWidgetChanged = False

  @QtCore.Slot()
  def drawTags(self):
    self.nTagRows += 1
    self.layout().addItem(self.layout().takeAt(self.layout().indexOf(self.newTagFrame)),6+self.nTagRows,0,1,self.TAGSPERROW)
    self.layout().addItem(self.layout().takeAt(self.layout().indexOf(self.moreTagsBut)),5+self.nTagRows,0)
    for iC in range(self.TAGSPERROW):
      self.tagWidgets.append(QtWidgets.QComboBox(self))
      self.tagWidgets[-1].setEditable(False)
      self.tagWidgets[-1].addItem(CProjectDescription.TAGLABEL,0)
      self.tagWidgets[-1].currentIndexChanged[int].connect(functools.partial(self.handleTagWidgetChanged,len(self.tagWidgets)-1))
      #for tid,parentTid,text,nUses in self.dbTags:
      #  print 'drawTags',tid,parentTid,text,nUses
      for tid,parentTid,text in self.dbTags:
        #print 'drawTags',tid,parentTid,text
        self.tagWidgets[-1].addItem(text,tid)
      self.layout().addWidget(self.tagWidgets[-1],4+self.nTagRows,iC)

  @QtCore.Slot(int,int)
  def handleTagWidgetChanged(self,widgetIndex,comboIndex):
    #print 'handleTagWidgetChanged',self.blockTagWidgetChanged,widgetIndex,comboIndex
    if self.blockTagWidgetChanged: return
    projectTagList = []
    for tW in self.tagWidgets:
      #if tW.currentIndex() != 0: tagList.append((str(tW.currentText()),str(tW.data(tW.currentIndex()))  ))
      if tW.currentIndex() != 0:
        tid = str(tW.itemData(tW.currentIndex()))
        if not tid in projectTagList: projectTagList.append(tid )
    #print 'CProjectDescription.handleTagWidgetChanged',projectTagList
    PROJECTSMANAGER().db().resetProjectTags(self.projectId,projectTagList)

  def load(self):
    if self.projectId is not None:
      commentInfo = PROJECTSMANAGER().db().getCommentInfo(projectId=self.projectId)
    else:
      commentInfo = {}
    #print 'CProjectDescription.load',commentInfo
    if len(commentInfo)==0:
      self.annotationWidget.setPlainText('')
      self.commentId = None
    else:
      self.annotationWidget.setPlainText(commentInfo['comment'])
      self.commentId = commentInfo['commentid']
    if self.projectId is not None:
      projectTagList = PROJECTSMANAGER().db().getProjectTags(projectId=self.projectId)
    else:
      projectTagList = []
    nR = (1 + ((len(projectTagList)-1)//self.TAGSPERROW)) - self.nTagRows
    for iR in range(nR): self.drawTags()
    for iT in range(len(projectTagList)):
      idx = self.tagWidgets[iT].findData(projectTagList[iT])
      if idx>=0:
        self.tagWidgets[iT].setCurrentIndex(idx)
      else:
        self.tagWidgets[iT].setCurrentIndex(0)
    for iT in range(len(projectTagList),len(self.tagWidgets)):
      self.tagWidgets[iT].setCurrentIndex(0)

  @QtCore.Slot(str)
  def save(self,projectId=None):
    if projectId is not None: self.projectId = projectId
    annotation = str(self.annotationWidget.toPlainText())
    #print 'CProjectDescription.save',annotation
    if self.commentId is None:
      PROJECTSMANAGER().db().createComment(projectId=self.projectId,comment=annotation)
    else:
      PROJECTSMANAGER().db().updateComment(self.commentId,annotation)
    projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=self.projectId)
    dbxml = os.path.join(projectDir,"DATABASE.db.xml")
    print("Saving",dbxml)
    PROJECTSMANAGER().db().exportProjectXml(self.projectId,fileName=dbxml)

  def deleteTag(self):
    if self.deleteTagCombo is None:
      win = QtWidgets.QDialog(self)
      win.setLayout(QtWidgets.QVBoxLayout())
      # Select tag to delete
      self.deleteTagCombo =  QtWidgets.QComboBox(self)
      win.layout().addWidget(self.deleteTagCombo)
      self.deleteTagCombo.currentIndexChanged[int].connect(self.selectDeleteTag)
      self.deleteTagWarning = QtWidgets.QLabel(self)
      win.layout().addWidget(self.deleteTagWarning)
      but = QtWidgets.QPushButton('Delete',self)
      win.layout().addWidget(but)
      but.clicked.connect(self.applyDeleteTag)

    self.deleteTagCombo.clear()
    self.deleteTagCombo.addItem('Delete tag..',0)
    for tid,parentTid,text in self.dbTags:
        self.deleteTagCombo.addItem(text,tid)
    self.deleteTagCombo.window().show()
    self.deleteTagCombo.window().raise_()



  @QtCore.Slot(int)
  def selectDeleteTag(self,indx):
    if indx == 0:
      self.deleteTagWarning.setText('')
      return
    tagId = self.deleteTagCombo.itemData(self.deleteTagCombo.currentIndex(),QtCore.Qt.UserRole).__str__()
    taggedProjects = PROJECTSMANAGER().db().getProjectsWithTag(tagId)
    if len(taggedProjects)>0:
      text = 'Tag used: '+ taggedProjects[0][1]
      for tP in taggedProjects[1:]: text+= ' '+tP[1]
      self.deleteTagWarning.setText(text)
    else:
      self.deleteTagWarning.setText('Tag unused')

  @QtCore.Slot()
  def applyDeleteTag(self):
    if self.deleteTagCombo.currentIndex() == 0: return
    tagId = self.deleteTagCombo.itemData(self.deleteTagCombo.currentIndex(),QtCore.Qt.UserRole).__str__()
    PROJECTSMANAGER().db().deleteTag(tagId)
    self.deleteTagCombo.window().hide()

  def deleteUnusedTags(self):
    PROJECTSMANAGER().db().deleteUnusedTags()

  def editTag(self):
    if self.editTagCombo is None:
      win = QtWidgets.QDialog(self)
      win.setLayout(QtWidgets.QVBoxLayout())
      # Select tag to edit
      self.editTagCombo =  QtWidgets.QComboBox(self)
      win.layout().addWidget(self.editTagCombo)
      self.editTagCombo.currentIndexChanged[int].connect(self.selectEditTag)
      line = QtWidgets.QHBoxLayout()
      line.addWidget(QtWidgets.QLabel('Change to',self))
      self.editTagLineEdit = QtWidgets.QLineEdit(self)
      line.addWidget(self.editTagLineEdit)
      but = QtWidgets.QPushButton('Save',self)
      line.addWidget(but)
      but.clicked.connect(self.saveEditTag)
      win.layout().addLayout(line)

    self.editTagCombo.clear()
    self.editTagCombo.addItem('Edit tag..',0)
    for tid,parentTid,text in self.dbTags:
        self.editTagCombo.addItem(text,tid)
    self.editTagLineEdit.clear()
    self.editTagCombo.window().show()
    self.editTagCombo.window().raise_()

  @QtCore.Slot(int)
  def selectEditTag(self,indx):
    if indx == 0:
      self.editTagLineEdit.setText('')
    else:
      self.editTagLineEdit.setText(self.editTagCombo.currentText())

  @QtCore.Slot()
  def saveEditTag(self):
    #print 'saveEditTag',self.editTagCombo.currentIndex()
    if self.editTagCombo.currentIndex() == 0: return
    tagId = self.editTagCombo.itemData(self.editTagCombo.currentIndex(),QtCore.Qt.UserRole).__str__()
    new = self.editTagLineEdit.text().__str__()
    #print 'saveEditTag',new
    if len(new)==0 or new == self.editTagCombo.currentText(): return
    try:
      PROJECTSMANAGER().db().updateTag(tagId,'text',new)
    except Exception as e:
      print('ERROR editing tag', e)
    self.editTagCombo.window().hide()


class CProjectDescriptionDialog(QtWidgets.QDialog):

  def __init__(self,parent,projectId,projectName):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('Edit annotation and tags for project:'+projectName)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.widget = CProjectDescription(self,projectId)
    self.layout().addWidget(self.widget)

    buttonBox = QtWidgets.QDialogButtonBox(self)
    self.layout().addWidget(buttonBox)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Close)
    but.clicked.connect(self.save)

  def save(self):
    ret = self.widget.save()
    self.close()
    self.deleteLater()
