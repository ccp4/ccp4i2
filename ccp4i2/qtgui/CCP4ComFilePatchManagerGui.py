from __future__ import print_function

"""
     CCP4ComFilePatchManagerGui.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
     Liz Potterton July 2013 - create and manage com file patches
"""

import os
from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4Data,CCP4Container
from qtgui import CCP4CustomisationGui
from core.CCP4ErrorHandling import *
from core.CCP4Modules import COMFILEPATCHMANAGER,WEBBROWSER,PROJECTSMANAGER

def openGui():
  if CComFilePatchManagerGui.insts is None:
    CComFilePatchManagerGui.insts = CComFilePatchManagerGui()
  CComFilePatchManagerGui.insts.show()
  CComFilePatchManagerGui.insts.raise_()


class CComFilePatchManagerGui(CCP4CustomisationGui.CCustomisationGui):

  insts = None


  def __init__(self,parent=None):
    CCP4CustomisationGui.CCustomisationGui.__init__(self,parent=parent,mode='comfilepatch',title='Task Parameter and Patches Manager')

  def manager(self):
    return COMFILEPATCHMANAGER()

  def handleNew(self):
   createWidget = CCreatePatchDialog(self)
   createWidget.show()
   createWidget.raise_()

  def handleEdit(self,selected=None):
    if selected is None:
      selected = self.customListView.selectedItem()
      if selected is None: return
    createWidget = CCreatePatchDialog(self,new=False)
    createWidget.loadPatch(selected)
    createWidget.show()
    createWidget.raise_()
    
class CCreatePatchDialog(QtWidgets.QDialog):

  created = QtCore.Signal(str)

  def __init__(self,parent=None,new=True):
    from qtgui import CCP4Widgets
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('Set custom task parameters')
    self.selectedJob = None
    self.setLayout(QtWidgets.QVBoxLayout())

    if not new:
      line = QtWidgets.QHBoxLayout()
      self.taskLabel = QtWidgets.QLabel('Apply patch to tasks: ',self)
      line.addWidget(self.taskLabel)
      self.layout().addLayout(line)
                   
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Edit command file from project:'))
    if new:
      self.projectCombo = QtWidgets.QComboBox(self)
      self.projectCombo.currentIndexChanged[int].connect(self.handleProjectChanged)
      line.addWidget(self.projectCombo)
    else:
      self.projectLabel = QtWidgets.QLabel(self)
      line.addWidget(self.projectLabel)
    line.addWidget(QtWidgets.QLabel('job:'))
    if new:
      self.jobsCombo = CCP4Widgets.CJobSelectionCombo(self)
      self.jobsCombo.currentIndexChanged[int].connect(self.handleJobSelection)
      line.addWidget(self.jobsCombo)
    else:
      self.jobsLabel =  QtWidgets.QLabel(self)
      line.addWidget(self.jobsLabel)
    line.addStretch(1)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Name for task parameters'))
    self.nameLineEdit= QtWidgets.QLineEdit(self)
    self.nameLineEdit.setToolTip('Enter a name for the custom task parameters')
    self.nameLineEdit.setMinimumWidth(300)
    line.addWidget(self.nameLineEdit)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Title to be used in interface'))
    self.titleLineEdit= QtWidgets.QLineEdit(self)
    self.titleLineEdit.setToolTip('Enter a short descriptive description')
    self.titleLineEdit.setMinimumWidth(300)
    line.addWidget(self.titleLineEdit)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    self.useControlParamsWidget = QtWidgets.QCheckBox('Use the control parameters set for this job',self)
    self.useControlParamsWidget.setChecked(True)
    line.addWidget(self.useControlParamsWidget)
    self.layout().addLayout(line)
    
    line = QtWidgets.QHBoxLayout()
    box = QtWidgets.QVBoxLayout()
    box.addWidget(QtWidgets.QLabel('Original command file'))
    self.originalTextEdit = QtWidgets.QTextEdit(self)
    self.originalTextEdit.setReadOnly(True)
    self.originalTextEdit.setObjectName('readOnly')
    box.addWidget(self.originalTextEdit)
    line.addLayout(box)
    box = QtWidgets.QVBoxLayout()
    box.addWidget(QtWidgets.QLabel('Edit the command file'))
    self.textEdit = QtWidgets.QTextEdit(self)
    box.addWidget(self.textEdit)
    line.addLayout(box)
    self.layout().addLayout(line)

    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton('Help',QtWidgets.QDialogButtonBox.HelpRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.help)
    but = buttonBox.addButton('Cancel',QtWidgets.QDialogButtonBox.RejectRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.cancel)
    but = buttonBox.addButton('Save task parameters',QtWidgets.QDialogButtonBox.AcceptRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.accept)
    self.layout().addWidget(buttonBox)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    
    self.reset()

  @QtCore.Slot()
  def handleProjectChanged(self):
    self.selectedProject = self.projectCombo.itemData(self.projectCombo.currentIndex()).__str__()
    self.jobsCombo.setProjectId(self.selectedProject)

  @QtCore.Slot(int)
  def handleJobSelection(self,indx):
    jobId = self.jobsCombo.getSelection()
    #print 'CCreatePatchDialog.handleJobSelection',jobList,errList
    self.selectedJob = jobId
    self.originalText = self.parent().manager().getComFileText(jobId=self.selectedJob)
    self.originalTextEdit.setReadOnly(False)
    self.originalTextEdit.setPlainText(self.originalText)
    self.originalTextEdit.setReadOnly(True)
    self.textEdit.setPlainText(self.originalText)
    #print 'CCreatePatchDialog.handleJobSelection',self.selectedJob
    
  def reset(self):
    self.textEdit.clear()
    self.originalTextEdit.clear()
    # get list of projectId,projectName,projectDir,parentId
    if hasattr(self,'projectCombo'):
      projectList =  PROJECTSMANAGER().db().listProjects(order='name')
      for project in projectList:
        item = self.projectCombo.addItem(project[1],project[0])
      self.jobsCombo.setProjectId(projectList[0][0])
    #self.titleLineEdit.setReadOnly(False)
    self.nameLineEdit.clear()
    self.titleLineEdit.clear()
    self.loadedFrom= None

    
  @QtCore.Slot()
  def help(self):
    WEBBROWSER().loadWebPage(helpFileName='customisation')

  @QtCore.Slot()
  def cancel(self):
    self.hide()

  def accept(self):

    if self.selectedJob is None and self.jobsLabel.text().__str__()!='Unknown':
      QtWidgets.QMessageBox.warning(self,'Task parameters','No job has been selected')
      return
    
    name = self.nameLineEdit.text().__str__()
    if len(name)<1:
      QtWidgets.QMessageBox.warning(self,'Task parameters','Please provide a unique name for task parameters')
      return

    import re
    name0 = re.sub('[^a-zA-Z0-9_-]','_',name)
    if name0 != name:
      self.nameLineEdit.setText(name0)
      QtWidgets.QMessageBox.warning(self,'Create task parameters','The task parameters name will be used as a directory name\n' +
                                'it has been changed to limited chacter set a-z,A-Z,0-9,_,-\n' +
                                "Please click 'Save task parameters' again" )
      return

    diry = COMFILEPATCHMANAGER().getDirectory(name=name)
    if os.path.exists(diry):
      import functools
      msgBox = QtWidgets.QMessageBox()
      msgBox.setWindowTitle('Create task parameter')
      msgBox.setText('There is already a task parameter directory called\n'+os.path.split(diry)[1])
      b = msgBox.addButton(QtWidgets.QMessageBox.Cancel)
      b = msgBox.addButton('Overwrite',QtWidgets.QMessageBox.ApplyRole)
      b.clicked.connect(functools.partial(self.handleOverwrite,name))
      msgBox.exec_()
      return
    
    self.createPatch(name)
    self.hide()

  @QtCore.Slot(str)
  def handleOverwrite(self,name):
    self.createPatch(name,True)
    self.hide()

  def createPatch(self,name,overwrite=False):
    title = None
    jobInfoList = []
    taskNameList = []
    #print 'createPatch selectedJob',self.selectedJob,self.loadedFrom
    if self.loadedFrom is not None:
      container,errMess = self.patchContainer(self.loadedFrom)
      #print 'createPatch loadedFrom',self.loadedFrom, container,errMess
      if container is not None:
        for item in container.taskNameList: taskNameList.append(item.__str__())
    elif  self.selectedJob is not None:
      jobId = self.selectedJob
      while jobId is not None:
        try:
          jobInfoList.append( PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['parentjobid','taskname']) )
          jobId = jobInfoList[-1]['parentjobid']
        except:
          jobId = None
      for jobInfo in jobInfoList: taskNameList.append(jobInfo['taskname'])
    #print 'createPatch taskNameList', jobInfoList,taskNameList
    text2 = self.textEdit.toPlainText().__str__()
    title = self.titleLineEdit.text().__str__()
    useControlParams = self.useControlParamsWidget.isChecked()
    print('createPatch title',title)
    err = COMFILEPATCHMANAGER().createPatch(name,title,taskNameList,self.selectedProject,self.selectedJob,self.originalText,text2,useControlParams,overwrite=True)
    if err.maxSeverity()>SEVERITY_WARNING:
      err.warningMessage('Create command file patch','Error saving patch',parent=self)
      return
    self.created.emit(name)

  def patchContainer(self,name):
    from core import CCP4ComFilePatchManager
    container = CCP4ComFilePatchManager.CPatchDefinition(parent=self,name=name)
    fileName=COMFILEPATCHMANAGER().getCustomFile(name=name)
    #print 'patchContainer fileName',fileName
    if fileName is None:
      return None,'Task parameter file does not exist: '+COMFILEPATCHMANAGER().getCustomFile(name=name,mustExist=False)
    try:      
      container.loadDataFromXml(fileName=fileName,loadHeader=True)
    except:
      return None,'Probably error in file: '+fileName
    return container,None

  def loadPatch(self,name):
    container,errMess = self.patchContainer(name)
    if errMess is not None:
      QtWidgets.QMessageBox.warning(self,'Error loading task parameter information',errMess)
      return
    if hasattr(self,'projectLabel'):
      self.selectedProject = container.projectId.__str__()
      try:
        projectName = PROJECTSMANAGER().db().getProjectInfo(projectId=self.selectedProject,mode='projectname')
      except:
        projectName = None
      if projectName is not None:
        self.projectLabel.setText(projectName)
      else:
        self.projectLabel.setText('Unknown')
      self.selectedJob = container.jobId.__str__()
      try:
        jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=self.selectedJob,mode=['jobnumber','taskname'])
      except:
        jobInfo = None
      if jobInfo is not None:
        self.jobsLabel.setText(jobInfo['jobnumber']+' '+jobInfo['taskname'])
      else:
        self.jobsLabel.setText('Unknown')
    if getattr(self,'taskLabel') is not None:
      label = 'Apply patch to tasks: '
      from core import CCP4TaskManager
      for taskName in container.taskNameList:
        title = CCP4TaskManager.TASKMANAGER().getTitle(taskName=taskName)
        label = label + title + ' '
      self.taskLabel.setText(label)
    self.originalText = str(container.text1)
    self.originalTextEdit.setReadOnly(False)
    self.originalTextEdit.setPlainText(self.originalText)
    self.originalTextEdit.setReadOnly(True)
    self.nameLineEdit.setText(name)
    self.textEdit.setPlainText(str(container.text2))
    title = container.header.pluginTitle.__str__()
    self.titleLineEdit.setText(title)
    self.useControlParamsWidget.setChecked(container.controlParameters.isSet())
    self.loadedFrom= name
    #self.titleLineEdit.setReadOnly(True)
    

