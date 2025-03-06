from __future__ import print_function

"""
     tasks/workflow/CTaskWorkflow.py: CCP4 GUI Project
     Copyright (C) 2014 STFC

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
     Liz Potterton July 2014
"""

import functools,os
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from core import CCP4WorkflowManager,CCP4Container
from qtgui import CCP4ProjectViewer
from core.CCP4Modules import WORKFLOWMANAGER,PROJECTSMANAGER
from core.CCP4TaskManager import TASKMANAGER

class CTaskWorkflow(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'workflow'
  TASKVERSION = 0.0
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)


  def setDefaultParameters(self):
    self.workflowName = self.container.header.pluginName.__str__()
    #print 'CTaskWorkflow.setDefaultParameters',self.workflowName
    self.workflowDef = CCP4WorkflowManager.CWorkflowDefinition(self,name=self.workflowName)
    fileName = WORKFLOWMANAGER().getCustomFile(self.workflowName)
    self.workflowDef.loadDataFromXml(fileName,function='WORKFLOW')
      

  def drawContents(self):
    folder = self.openFolder(folderFunction='inputData',title='Input Data')    
    self.autoGenerate(self.container.inputData)
    folder = self.openFolder(folderFunction='controlParameters',title='Control parameters for sub-tasks')

    self.createLine(['advice','Edit parameters for the sub-jobs in the workflow..'])
    line = self.createLine()
    self.listWidget = QtWidgets.QListWidget(self)
    self.listWidget.itemClicked.connect(self.openSubJobTaskWidget)
                            
    line.addWidget(self.listWidget)
    
    for jobName in self.workflowDef.jobDef.dataOrder()[1:]:
      taskName = self.workflowDef.jobDef.get(jobName).taskName.__str__()
      item = QtWidgets.QListWidgetItem(jobName.split('_')[1]+ ' ' + TASKMANAGER().getTitle(taskName))
      item.setData(QtCore.Qt.UserRole,jobName)
      self.listWidget.addItem( item )

  def paramsFilePath(self,jobName):
    splitPath =  os.path.split(PROJECTSMANAGER().makeFileName(jobId=self.jobId(),mode='JOB_INPUT'))
    return os.path.join(splitPath[0],jobName+'_'+splitPath[1])


  @QtCore.Slot('QListWidgetItem')
  def openSubJobTaskWidget(self,listWidgetItem):
    # This follows the method of CProjectViewer.openTaskMainWindow() to create a task input in a separate window
    jobName = listWidgetItem.data(QtCore.Qt.UserRole).__str__()
    if jobName in self.subJobTaskWidgets:
      self.subJobTaskWidgets[jobName].window().show()
      self.subJobTaskWidgets[jobName].window().raise_()
      return

    taskName = self.workflowDef.jobDef.get(jobName).taskName.__str__()
    
    container = CCP4Container.CContainer(parent=self,name=taskName,
                          definitionFile=os.path.join(WORKFLOWMANAGER().getDirectory(self.workflowName),jobName+'.def.xml'))
    container.loadDataFromXml(os.path.join(WORKFLOWMANAGER().getDirectory(self.workflowName),jobName+'.params.xml'))
    if os.path.exists(self.paramsFilePath(jobName)): container.loadDataFromXml(self.paramsFilePath(jobName))
    
    #print 'CTaskWorkflow.openControlParameters container',container.dataOrder()
    taskInp = CCP4ProjectViewer.CTaskInputFrame(self)
    taskInp.createTaskWidget(taskName=taskName,container=container,excludeInputData=True)
    self.subJobTaskWidgets[jobName] = taskInp.taskWidget
    widget = self.subJobTaskWidgets[jobName].widget
    if isinstance(widget,QtWidgets.QTabWidget):
      if str(widget.tabText(0)) == 'Input Data': widget.removeTab(0)
    else:
      if widget.layout().itemAt(0).widget().title() == 'Input Data': widget.layout().itemAt(0).widget().hide()
    self.subJobTaskWidgets[jobName].setParamsFileName(self.paramsFilePath(jobName))
    projectName=PROJECTSMANAGER().db().getProjectInfo(projectId=self.projectId(),mode='projectname')
    win = CCP4ProjectViewer.CTaskMainWindow(self,projectName=projectName,jobId=self.jobId())
    win.jobTitle.setText( 'Job '+PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId(),mode='jobnumber')+'.'+jobName.split('_')[1]+
                          ' '+self.workflowName+ ' : ' + TASKMANAGER().getTitle(taskName) )
    for but in ['run','view','clone']:
      win.buttons.button(but).hide()
    win.centralWidget().layout().insertWidget(1,self.subJobTaskWidgets[jobName].parent())
    win.windowAboutToClose.connect(functools.partial(self.handleClosingSubTaskWindow,jobName))
    win.show()
    win.raise_()

