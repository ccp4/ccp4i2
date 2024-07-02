from __future__ import print_function

"""
     CCP4WorkflowManagerGui.py: CCP4 GUI Project
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
     Liz Potterton July 2013 - create and manage workflows
"""

import os
import re
import functools

from PySide6 import QtGui, QtWidgets,QtCore
from core import CCP4Container
from qtgui import CCP4CustomisationGui
from core.CCP4ErrorHandling import *
from core.CCP4Modules import WORKFLOWMANAGER,WEBBROWSER,PROJECTSMANAGER

def openWorkflowManagerGui():
    if CWorkflowManagerGui.insts is None:
        CWorkflowManagerGui.insts = CWorkflowManagerGui()
    CWorkflowManagerGui.insts.show()
    CWorkflowManagerGui.insts.raise_()


class CWorkflowManagerGui(CCP4CustomisationGui.CCustomisationGui):

    insts = None
    ERROR_CODES = {201 : {'description' : 'Unknown error creating workflow'},
                   202 : {'description' : 'Unknown error attempting to edit workflow'}}

    def __init__(self, parent=None):
        CCP4CustomisationGui.CCustomisationGui.__init__(self, parent=parent, mode='workflow', title='Workflow Manager')

    def manager(self):
        return WORKFLOWMANAGER()

    def handleNew(self):
        if self.createWidget is None:
            self.createWidget = CCreateWorkflowDialog(self)
            self.createWidget.workflowCreated.connect(self.handleEdit)
        else:
            self.createWidget.reset()
        self.createWidget.show()
        self.createWidget.raise_()

    @QtCore.Slot('QListWidgetItem')
    def handleEdit(self,selected=None):
        if selected is None:
            selected = self.customListView.selectedItem()
        if selected is None:
            return
        try:
            editor = CWorkflowEditDialog(self, name=selected)
        except CException as e:
            e.warningMessage('Create workflow', 'Error opening workflow editor.\nPossibly a file is corrupted.\n' + self.manager().getDirectory(selected), parent=self)
            return
        except Exception as e:
            e = CException(self.__class__, 202, exc_info=sys.exc_info())
            e.warningMessage('Create workflow', 'Error opening workflow editor.\nPossibly a file is corrupted.\n' + self.manager().getDirectory(selected), parent=self)
            return
        editor.show()


class CCreateWorkflowDialog(QtWidgets.QDialog):

    workflowCreated = QtCore.Signal(str)

    ERROR_CODES = {201 : {'description' : 'Unkown error saving workflow'}}

    def __init__(self, parent=None):
        from qtgui import CCP4Widgets
        QtWidgets.QDialog.__init__(self, parent)
        self.setWindowTitle('Create a workflow')
        self.setLayout(QtWidgets.QVBoxLayout())
        self.projectCombo = QtWidgets.QComboBox(self)
        self.projectCombo.setEditable(False)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Select from project'))
        line.addWidget(self.projectCombo)
        self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Select jobs'))
        self.jobsLineEdit = CCP4Widgets.CJobSelectionLineEdit(self)
        line.addWidget(self.jobsLineEdit)
        self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Name'))
        self.nameLineEdit= QtWidgets.QLineEdit(self)
        self.nameLineEdit.setToolTip('Enter a one-word name for the workflow')
        line.addWidget(self.nameLineEdit)
        self.layout().addLayout(line)
        buttonBox = QtWidgets.QDialogButtonBox(self)
        but = buttonBox.addButton('Help', QtWidgets.QDialogButtonBox.HelpRole)
        but.setFocusPolicy(QtCore.Qt.NoFocus)
        but.clicked.connect(self.help)
        but = buttonBox.addButton('Cancel', QtWidgets.QDialogButtonBox.RejectRole)
        but.setFocusPolicy(QtCore.Qt.NoFocus)
        but.clicked.connect(self.cancel)
        but = buttonBox.addButton('Create workflow', QtWidgets.QDialogButtonBox.AcceptRole)
        but.setFocusPolicy(QtCore.Qt.NoFocus)
        but.clicked.connect(self.accept)
        self.layout().addWidget(buttonBox)
        self.reset()

    def reset(self):
        # get list of projectId, projectName, projectDir, parentId
        projectList =  PROJECTSMANAGER().db().listProjects(order='name')
        for project in projectList:
            item = self.projectCombo.addItem(project[1], project[0])
        self.jobsLineEdit.clear()
        self.nameLineEdit.clear()

    @QtCore.Slot()
    def help(self):
        WEBBROWSER().loadWebPage(helpFileName='customisation')

    @QtCore.Slot()
    def cancel(self):
        self.hide()

    @QtCore.Slot()
    def accept(self):
        projectId = self.projectCombo.itemData(self.projectCombo.currentIndex()).__str__()
        self.jobsLineEdit.setProjectId(projectId)
        jobList,errList = self.jobsLineEdit.getSelection()
        if len(errList) > 0:
            text = 'Error interpreting job selection:'
            for err in errList:
                text = text + '\n' + err
            QtWidgets.QMessageBox.warning(self, 'Create workflow', text)
            return
        name = self.nameLineEdit.text().__str__()
        if len(name) < 1:
            QtWidgets.QMessageBox.warning(self, 'Create workflow','Please provide a unique one-word name for workflow')
            return
        one_word_name = re.sub('[^a-zA-Z0-9_-]', '_', name)
        if one_word_name != name:
            self.nameLineEdit.setText(one_word_name)
            QtWidgets.QMessageBox.warning(self,'Create workflow','The workflow name will be used as a directory name\n' +
                                      'it has been changed to limited chacter set a-z,A-Z,0-9,_,-\n' +
                                      "Please click 'Create workflow' again" )
            return
        diry = WORKFLOWMANAGER().getDirectory(name=name)
        if os.path.exists(diry):
            msgBox = QtWidgets.QMessageBox()
            msgBox.setWindowTitle('Create workflow')
            msgBox.setText('There is already a workflow directory called\n' + os.path.split(diry)[1])
            b = msgBox.addButton(QtWidgets.QMessageBox.Cancel)
            b = msgBox.addButton('Overwrite', QtWidgets.QMessageBox.ApplyRole)
            b.clicked.connect(functools.partial(self.handleOverwrite, (projectId, jobList, name)))
            msgBox.exec_()
            return
        self.createWorkflow(projectId, jobList, name)
        self.hide()
    
    @QtCore.Slot(tuple)
    def handleOverwrite(self, args):
        projectId, jobList, name = args
        self.createWorkflow(projectId, jobList, name, True)
        self.hide()

    def createWorkflow(self,projectId,jobList,name,overwrite=False,title=None):
        try:
            err = WORKFLOWMANAGER().createWorkflow(projectId=projectId, jobList=jobList, name=name, overwrite=overwrite, title=title)
        except CException as err:
            err.warningMessage('Create workflow', 'Error saving workflow', parent=self)
            return
        except Exception as e:
            err = CException(self.__class__,201,exc_info=sys.exc_info())
            err.warningMessage('Create workflow', 'Error saving workflow', parent=self)
            return
        if err.maxSeverity()>SEVERITY_WARNING:
            err.warningMessage('Create workflow', 'Error saving workflow', parent=self)
            return
        self.workflowCreated.emit(name)


class CWorkflowEditDialog(QtWidgets.QDialog):

    def __init__(self,parent,name):
        QtWidgets.QDialog.__init__(self,parent)
        self.setWindowTitle('Edit workflow interface:'+str(name))
        self.setLayout(QtWidgets.QVBoxLayout())
        MARGIN = 5
        self.layout().setSpacing(MARGIN)
        self.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
        label = QtWidgets.QLabel('Short title for workflow task',self)
        label.setObjectName('italic')
        self.layout().addWidget(label)
        self.setWindowTitle('Edit workflow: '+name)
        self.name = name
        directory = WORKFLOWMANAGER().getDirectory(self.name)
        from core import CCP4WorkflowManager
        self.workflowDef = CCP4WorkflowManager.CWorkflowDefinition(self,name=self.name)
        fileName = WORKFLOWMANAGER().getCustomFile(self.name)
        #print 'CWorkflowEditDialog workflowDef',self.name,self.workflowDef.header.pluginName,fileName
        self.workflowDef.loadDataFromXml(fileName,function='WORKFLOW')
        self.jobDefs = CCP4WorkflowManager.CWorkflowContainerList()
        self.jobDefs.loadContentsForWorkflow(directory)
        #print 'CWorkflowEditDialog.init jobDefs',self.jobDefs['0'].dataOrder()
        self.titleWidget = QtWidgets.QLineEdit(self)
        self.layout().addWidget(self.titleWidget)
        label = QtWidgets.QLabel('Label for input files in workflow gui',self)
        label.setObjectName('italic')
        self.layout().addWidget(label)
        self.inputWidgets = []
        for i in range(len(self.jobDefs['job_0'].inputData.dataOrder())):
            self.inputWidgets.append(CWorkflowEditInputWidget(self))
            self.layout().addWidget(self.inputWidgets[-1])
        label = QtWidgets.QLabel('Select overall workflow output files and file names to appear in gui',self)
        label.setObjectName('italic')
        self.layout().addWidget(label)
        self.outputFrame = CWorkflowOutputChooser(self)
        self.layout().addWidget(self.outputFrame)
        line = QtWidgets.QHBoxLayout()
        line.addStretch(2)
        butBox = QtWidgets.QDialogButtonBox(self)
        but = butBox.addButton('Save workflow', QtWidgets.QDialogButtonBox.ActionRole)
        but.setDefault(False)
        but.clicked.connect(self.apply)
        but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
        but.clicked.connect(self.close)
        line.addWidget(butBox)
        line.addStretch(2)
        self.layout().addLayout(line)
        self.outputFrame.load(self.workflowDef, self.jobDefs)
        self.updateViewFromModel()

    @QtCore.Slot()
    def close(self):
        self.hide()
        self.deleteLater()

    @QtCore.Slot()
    def apply(self):
        self.updateModelFromView()
        self.save()
        self.close()

    def updateViewFromModel(self):
        self.titleWidget.setText(self.jobDefs['job_0'].header.pluginTitle.__str__())
        inputFileList = self.jobDefs.inputFileList(workflowDef=self.workflowDef)
        ii = -1
        for inputFile in inputFileList:
            ii += 1
            try:
                self.inputWidgets[ii].desc.setText(inputFile['className'] + ' ' + inputFile['name'] + ' ' + inputFile.get('taskName', ''))
                self.inputWidgets[ii].label.setText(inputFile['label'])
                self.inputWidgets[ii].setObjectName(inputFile['name'])
            except:
                print('ERROR in CWorkflowEditDialog.updateViewFromModel for inputFile',inputFile)
        for cb in self.outputFrame.findChildren(QtWidgets.QCheckBox):
            cb.setChecked(False)
        for outFile in self.workflowDef.jobDef.job_0.output:
            cbName = str(outFile.fromJob) + '_' + str(outFile.fromKey)
            cb = self.outputFrame.findChild(QtWidgets.QCheckBox, cbName)
            if cb is not None:
                cb.setChecked(True)
            le = self.outputFrame.findChild(QtWidgets.QLineEdit, cbName)
            if le is not None:
                #le.setText(self.jobDefs.__getattr__(str(outFile.fromJob)).outputData.__getattr__(str(outFile.fromKey)).annotation.__str__())
                le.setText(outFile.annotation.__str__())

    def updateModelFromView(self):
        from core import CCP4WorkflowManager
        self.jobDefs['job_0'].header.pluginTitle.set(str(self.titleWidget.text()))
        for inputWidget in self.inputWidgets:
            label = str(inputWidget.label.text()).strip()
            if len(label) == 0:
                label = None
            self.jobDefs['job_0'].inputData.__getattr__(str(inputWidget.objectName())).setQualifier('guiLabel', label)
        nOutputFiles = 0
        workflowOutput = CCP4WorkflowManager.CWorkflowDataFlowList(parent=self)
        workflowOutput.setObjectName('output')
        jobZeroOutputParams = CCP4Container.CContainer(parent=self, name='outputData')
        #print 'CWorkflowEditDialog.updateModelFromView jobDefs',self.jobDefs.dataOrder()
        for jobKey in self.workflowDef.jobDef.dataOrder():
            job = self.workflowDef.jobDef.__getattr__(jobKey)
            #print 'CWorkflowEditDialog.updateModelFromView jobDefs',self.jobDefs[str(jobKey)].dataOrder()
            for fileObj in job.allOutputFiles:
                obj = self.jobDefs[jobKey].outputData.__getattr__(fileObj.key.__str__())
                cb = self.outputFrame.findChild(QtWidgets.QCheckBox, jobKey + '_' + fileObj.key.__str__() )
                le = self.outputFrame.findChild(QtWidgets.QLineEdit, jobKey + '_' + fileObj.key.__str__() )
                if cb.isChecked():
                    nOutputFiles += 1
                    outKey = 'OUTPUT' + str(nOutputFiles)
                    annotation = le.text().__str__()
                    if len(annotation) == 0:
                        annotation = None
                    workflowOutput.append({'toKey': outKey, 'fromJob' : jobKey, 'fromKey' : fileObj.key.__str__(), 'annotation' : annotation } )
                    qualifiers = obj.qualifiers(default=False, custom=True)
                    for item in ['contentFlag', 'subType', 'sameCrystalAs']:
                        if item in qualifiers:
                            del qualifiers[item]
                    newObj = obj.__class__(parent=jobZeroOutputParams,name=outKey,qualifiers=qualifiers)
                    jobZeroOutputParams.addObject(newObj)
        self.jobDefs['job_0'].__dict__['_value']['outputData'] =jobZeroOutputParams
        self.workflowDef.jobDef['job_0'].__dict__['_value']['output'] = workflowOutput

    def save(self):
        workflowDir = WORKFLOWMANAGER().getDirectory(self.name)
        self.workflowDef.header.pluginTitle = self.titleWidget.text().__str__()
        self.workflowDef.saveDataToXml(WORKFLOWMANAGER().getCustomFile(self.name,mustExist=False),function='WORKFLOW')
        self.jobDefs['job_0'].saveContentsToXml(fileName=os.path.join(workflowDir,'job_0.def.xml'))
        for jobKey in self.jobDefs.dataOrder()[1:]:
            self.jobDefs[jobKey].saveDataToXml(fileName=os.path.join(workflowDir,jobKey+'.params.xml'))
        WORKFLOWMANAGER().listChanged.emit()

class CWorkflowEditInputWidget(QtWidgets.QFrame):

    def __init__(self,parent=None):
        QtWidgets.QFrame.__init__(self,parent=None)
        self.setLayout(QtWidgets.QHBoxLayout())
        MARGIN = 1
        self.layout().setSpacing(MARGIN)
        self.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
        self.desc = QtWidgets.QLabel(self)
        self.desc.setMinimumWidth(200)
        self.layout().addWidget(self.desc)
        self.label = QtWidgets.QLineEdit(self)
        self.label.setMinimumWidth(200)
        self.layout().addWidget(self.label)

class CWorkflowOutputChooser(QtWidgets.QScrollArea):

    def __init__(self,parent=None):
        QtWidgets.QScrollArea.__init__(self,parent)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        
        frame = QtWidgets.QFrame(self)
        frame.setLayout(QtWidgets.QGridLayout())
        MARGIN = 1
        frame.layout().setSpacing(MARGIN)
        frame.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)    
        self.setWidget(frame)
        self.setWidgetResizable(1)

    def load(self,workflowDef,jobDefs):
        layout = self.widget().layout()
        #print 'CWorkflowOutputChooser.load',jobDefs.dataOrder()
        indx = -1
        for jobKey in workflowDef.jobDef.dataOrder()[1:]:
            indx += 1
            job = workflowDef.jobDef.__getattr__(jobKey)
            jobOutput = jobDefs.__getattr__(jobKey).outputData
            layout.addWidget(QtWidgets.QLabel(job.taskName.__str__(),self),indx,0,1,2)
            for fileObj in job.allOutputFiles:
                indx += 1
                key = fileObj.key.__str__()
                line = QtWidgets.QHBoxLayout()
                cb = QtWidgets.QCheckBox(key,self)
                cb.setObjectName(jobKey+'_'+key)
                #print 'CWorkflowOutputChooser.load',fileObj.key,fileObj.ifOverallOutput
                #if fileObj.ifOverallOutput: cb.setChecked(True)
                layout.addWidget(cb,indx,0)
                le = QtWidgets.QLineEdit(self)
                le.setObjectName(jobKey+'_'+key)
                layout.addWidget(le,indx,1)
        self.widget().show()

