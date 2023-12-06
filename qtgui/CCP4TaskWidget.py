from __future__ import print_function

"""
     CCP4TaskWidget.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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
     GNU Lesser General Public License for more details."""

"""
     Liz Potterton Jan 2010 - Create CCP4TaskWidget prototype
"""

##@package CCP4TaskWidget (QtGui)  Widget to view CCP4 tasks
import os
import re
import traceback
import functools
from PySide2 import QtGui, QtWidgets, QtCore
from qtgui import CCP4Widgets
from core.CCP4Modules import *
from core.CCP4ErrorHandling import *
from core.CCP4DataManager import DATAMANAGER
from core.CCP4Config import DEVELOPER

MARGIN = 1
WIDTH = 600
PADDING_ALLOWANCE = 0
LINEHEIGHT = 20
USER_NOVICE = 0
USER_REGULAR = 1
USER_EXPERT = 2

def whatNext(self, jobId=None, childJobs=[], childTaskName=None):
    return []


class CFolderAttributes:

    def __init__(self):
        self.atts = {'all' : {'editable' : True}}
        #print 'CFolderAttributes.__init__'

    def setAttribute(self, attribute='editable', folderFunction='all', value=None):
        if not attribute in ['editable'] or folderFunction is None or value is None:
            print('Invalid input to CFolderAttributes.setAttribute', attribute, folderFunction, value)
            return
        if folderFunction not in self.atts:
            self.atts[folderFunction] = {}
        self.atts[folderFunction][attribute] = value
        #print 'CFolderAttributes.setAttribute',self.atts

    def attribute(self, attribute='editable', folderFunction='all'):
        if folderFunction not in self.atts:
            folderFunction = 'all'
        if attribute in self.atts[folderFunction]:
            #print 'CFolderAttributes.attribute',attribute,folderFunction,self.atts[folderFunction][attribute]
            return  self.atts[folderFunction][attribute]
        else:
            return NotImplemented

    def allAttributes(self,folderFunction='all'):
        if folderFunction not in self.atts:
            return  self.atts['all']
        else:
            atts = {}
            atts.update(self.atts['all'])
            atts.update(self.atts[folderFunction])
            return atts

class CTaskWidgetDrawMethods:

    ERROR_CODES = {100 : {'description' : 'No definition found for data name'},
                   101 : {'description' : 'Parameter name not recognised in setMenuText'},
                   102 : {'description' : 'Wrong number of menu text items in setMenuText'},
                   103 : {'description' : 'Reference to unknown enumerator value(s) in setMenuText'},
                   104 : {'description' : 'Parameter name not recognised in setTooltip'},
                   105 : {'description' : 'Parameter name not recognised in setLabel'}}
    STACK_LAYOUT = True

    def __init__(self):
        self.toggleList = []
        self.subFrame = None
        self.stack = None
        self.tabWidget = None
        self.currentFolder = None
        self.currentFolderLayout = None
        self.jobTitleWidget = None
        self.ignoreInput = False

    def openTabFrame(self, title=None, toolTip=None):
        if self.tabWidget is None:
            self.tabWidget = QtWidgets.QTabWidget(self)
            self.currentFolderLayout.addWidget( self.tabWidget)
        else:
            self.tabWidget.widget(self.tabWidget.count()-1).layout().addStretch(2)
        frame = QtWidgets.QFrame()
        frame.setLayout(QtWidgets.QVBoxLayout())
        frame.layout().setSpacing(MARGIN)
        frame.layout().setContentsMargins(MARGIN, MARGIN, MARGIN, MARGIN)
        self.tabWidget.addTab(frame, title)
        if toolTip is not None:
            self.tabWidget.setTabToolTip(self.tabWidget.count() - 1, toolTip)
        else:
            #Worth setting title astool tip so it is readable if title display is trucated
            self.tabWidget.setTabToolTip(self.tabWidget.count() - 1, title)

    def closeTabFrame(self):
        self.tabWidget.widget(self.tabWidget.count() - 1).layout().addStretch(2)
        self.tabWidget = None

    def openSubFrame(self, toggle=[], toggleFunction=[], frame=False, title=None, tip=''):
        if self.ignoreInput:
            return
        self.closeSubFrame()
        if self.currentFolder is None:
            print('ERROR in openSubFrame: there is not current open folder')
            return [1, 'ERROR in openSubFrame: there is not current open folder']
        if title is not None:
            self.subtitle = self.createLine(['subtitle', title, tip])
        self.subFrame = CSubFrame(self)
        self.subFrame.setObjectName('tasksubframe')
        self.currentFolderLayout.addWidget(self.subFrame)
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(MARGIN)
        layout.setContentsMargins(MARGIN, MARGIN, MARGIN, MARGIN)
        self.subFrame.setLayout(layout)
        self.noChildFrame = False
        #print 'openSubFrame args',self.subFrame,arg
        if len(toggle) == 2:
            self.setToggle(target=self.subFrame, parameter=toggle[0], state=toggle[1])
            if title is not None:
                self.setToggle(target=self.subtitle, parameter=toggle[0], state=toggle[1])
        #print 'setToggle',target
        elif len(toggle) == 3:
            self.setToggle(target=self.subFrame, parameter=toggle[0], state=toggle[1], values=toggle[2])
            if title is not None:
                self.setToggle(target=self.subtitle, parameter=toggle[0], state=toggle[1], values=toggle[2])
        elif len(toggleFunction) == 2:
            self.setToggleFunction(target=self.subFrame, toggleFunction=toggleFunction[0], parameterList=toggleFunction[1])
            if title is not None:
                self.setToggleFunction(target=self.subtitle, toggleFunction=toggleFunction[0], parameterList=toggleFunction[1]) 
        if frame:
            self.subFrame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            self.noChildFrame = True
        return self.subFrame

    def closeSubFrame(self):
        if self.ignoreInput:
            return
        if self.subFrame is not None and self.noChildFrame:
            children = self.subFrame.findChildren(CCP4Widgets.CComplexLineWidget)
            for child in children:
                child.setFrameShape(QtWidgets.QFrame.NoFrame)
            self.subFrame.layout().update()
        self.subFrame = None

    def openStack(self, controlVar=None):
        if self.ignoreInput:
            return
        if self.stack is not None:
            self.closeStack()
        if self.STACK_LAYOUT:
            self.stack = CStackedLayout(self, controlVar)
        else:
            self.stack = CStackedWidget(self, controlVar)
        return self.stack

    def closeStack(self):
        if self.ignoreInput:
            return
        if self.stack is not None:
            #print 'closeStack',self.stack,self.currentFolderLayout
            if self.STACK_LAYOUT:
                if self.subFrame is not None:
                    self.subFrame.layout().addLayout(self.stack)
                else:
                    self.currentFolderLayout.addLayout(self.stack)
            else:
                if self.subFrame is not None:
                    self.subFrame.layout().addWidget(self.stack)
                else:
                    self.currentFolderLayout.addWidget(self.stack)
                #self.stack.update()
            self.stack = None

    def createTitleLine(self, name='TITLE'):
        if self.ignoreInput:
            return
        self.createLine(['label', 'Comment', 'widget', name])

    def createLine(self, definition=[], appendLine=None, toggle=[], toggleFunction=[]):
        if self.ignoreInput:
            return
        container = self.parentTaskWidget()._container
        line = CTaskLine(self)
        if len(definition) > 0:
            #if self.parentTaskWidget().layoutMode == 'TAB':
            #  rv = line.draw(definition=definition,dataContainer=container,setupFolder=self.currentFolder.setupUpdateFolder,attributes=self.currentFolder._attributes)
            #else:
            rv = line.draw(definition=definition, dataContainer=container, attributes=self.currentFolder._attributes)
            self.myException.extend(rv, stack=False)
        else:
            rv = []
        if len(rv) > 0:
            line.deleteLater()
        else:
            if appendLine is not None:
                appendLine.addWidget(line)
            elif self.stack is not None:
                self.stack.addWidget(line)
            elif self.subFrame is not None:
                self.subFrame.layout().addWidget(line)
            elif self.tabWidget is not None:
                self.tabWidget.widget(self.tabWidget.count() - 1).layout().addWidget(line)
            else:
                self.currentFolderLayout.addWidget(line)
        #print 'CTaskWidgetDrawMethods.createLine toggle',self.subFrame, definition
        if len(toggle) > 0:
            if len(toggle) == 2:
                self.setToggle(line, toggle[0], toggle[1], [True])
            else:
                self.setToggle(line, toggle[0], toggle[1], toggle[2])
        if len(toggleFunction) > 0:
            self.setToggleFunction(line, toggleFunction[0], toggleFunction[1])
        return line

    def applyToggles(self):
        if self.ignoreInput:
            return
        for toggle in self.toggleList:
            if not toggle.connected:
                toggle.handleSignal()
                toggle.makeConnection()
        childStacks = self.findChildren(CStackedLayout)
        for stack in childStacks: stack.update()

    def setToggle(self, target=None, parameter=None, state='open', values =[]):
        if self.ignoreInput:
            return
        #print 'setToggle',target
        if target is None or parameter is None or len(values) == 0:
            return
        self.toggleList.append(CToggle(self, target, parameter, state, values))
    
    def setToggleFunction(self, target=None, toggleFunction=None, parameterList=[]):
        if self.ignoreInput:
            return
        if target is None or toggleFunction is None or len(parameterList) == 0:
            return
        self.toggleList.append(CToggle(self, target, parameterList=parameterList, toggleFunction=toggleFunction))

    def createRadioGroup(self, label=None, itemList=[], objectName=None):
        if self.ignoreInput:
            return
        if objectName is not None:
            varObj = self.parentTaskWidget()._container.getObject(objectName)
            itemList = varObj.qualifiers('enumerators')
        else:
            varObj = None
        line = QtWidgets.QFrame(self)
        line.setLayout(QtWidgets.QHBoxLayout())
        if label is not None:
            line.layout().addWidget(QtWidgets.QLabel(label, self))
        group = QtWidgets.QButtonGroup(self)
        idx = -1
        for item in itemList:
            idx = idx + 1
            button = QtWidgets.QRadioButton(item, self)
            group.addButton(button, idx)
            line.layout().addWidget(button)
        self.currentFolderLayout.addWidget(line)
        return group

    def createJobTitle(self, followFrom=True):
        from core import CCP4TaskManager
        self.jobHeaderFrame = self.parentTaskWidget().openSubFrame(frame=[False])
        self.jobHeaderFrame.setObjectName('jobHeaderFrame') # so that it can be styled
        if self.ignoreInput:
            return
        #print 'createJobTitle editable',self.parentTaskWidget().folderAttributes.attribute('editable'),self.parentTaskWidget().jobId()
        line = self.parentTaskWidget().createLine(['label', 'Job title', 'widget', 'jobTitle'])
        try:
            self.jobTitleWidget = line.layout().itemAt(1).widget().widget
        except:
            return
        t = self.getWidget('jobTitle').model
        if t is None: return
        #print self.getWidget('jobTitle'), t
        if not self.parentTaskWidget().folderAttributes.attribute('editable'):
            #Its an old job and job title may have been changed in database since saving params.xml
            title = PROJECTSMANAGER().db().getJobInfo(jobId=self.parentTaskWidget().jobId(), mode='jobtitle')
            #print 'createJobTitle title from db',title
            if title is not None:
                t.set(title) 
        if not t.isSet() or len(t) == 0:
            t.set(CCP4TaskManager.TASKMANAGER().getShortTitle(self.parentTaskWidget().taskName()))
#FIXME - This is major puzzle - I do not know how a CBoldLabel can emit editingFinished
        if hasattr(self.jobTitleWidget,"editingFinished"):
            self.jobTitleWidget.editingFinished.connect(self.saveTitle)
        self.jobTitleWidget.setToolTip('Short description of job to appear in Job list')
        if followFrom and self.parentTaskWidget().folderAttributes.attribute('editable'):
            self.parentTaskWidget().createLine(['widget','followFrom'])
        self.parentTaskWidget().closeSubFrame()

    def createPatchFrame(self):
        if self.ignoreInput:
            return
        taskWidget = self.parentTaskWidget()
        # Force 'fix' of the patchSelection object before drawing
        try:
            patchSelectionObj = taskWidget._container.guiAdmin.patchSelection
        except:
            return
        patchSelectionObj.set(patchSelectionObj.fix({'taskName' : taskWidget.taskName(), 'patch' : patchSelectionObj.get('patch') }))
        if len(patchSelectionObj.getPatchList()) > 0:
            taskWidget.createLine(['widget', 'patchSelection'])
            patchSelectionObj.dataChanged.connect(self.handlePatchChange)
        self.handlePatchChange()

    @QtCore.Slot()
    def handlePatchChange(self):
        taskWidget = self.parentTaskWidget()
        #print 'CTaskWidget.handlePatchChange patch',taskWidget._container.guiAdmin.patchSelection.patch
        try:
            if not taskWidget._container.guiAdmin.patchSelection.patch.isSet():
                return
        except:
            return
        controlParameters = COMFILEPATCHMANAGER().getComFileControlParameters(name=taskWidget._container.guiAdmin.patchSelection.patch.__str__())
        #print 'CTaskWidget.handlePatchChange',controlParameters
        taskWidget._container.controlParameters.update(controlParameters)
        #print 'CTaskWidget after update',taskWidget._container.controlParameters
        for key in controlParameters.dataOrder():
            widget = self.findChild(CCP4Widgets.CViewWidget,key)
            if widget is not None:
                widget.updateViewFromModel()

    @QtCore.Slot()
    def saveTitle(self):
        from core import CCP4Modules
        jobId = self.parentTaskWidget().jobId()
        if jobId is None:
            return
        text = str(self.jobTitleWidget.text())
        print('CTaskWidgetDrawMethods.saveTitle', jobId,text)
        try:
            CCP4Modules.PROJECTSMANAGER().db().updateJob(jobId=jobId, key='jobTitle', value=text)
        except:
            pass

class CTaskWidget(QtWidgets.QFrame):

    doFix = QtCore.Signal()
    followFromJobUpdated = QtCore.Signal(str,str)
    launchJobRequestSignal = QtCore.Signal(str,dict)

    EDITOR = False
    AUTOPOPULATEINPUT = True
    USEPLUGIN = None
    ERROR_CODES = {101 : {'description' : 'Error updating view from model'},
                   102 : {'description' : 'Error updating model from view'},
                   103 : {'description' : 'Error in cootFix'}}

    def __init__(self, parent=None, title=None, projectId=None, jobId=None, layoutMode=None):
        QtWidgets.QFrame.__init__(self, parent)
        from core import CCP4Container
        from core import CCP4Modules
        self.setTitle(title)
        self.helpFile = ''
        self.programHelpFile = ''
        # Dict of sub-job widgets (eg workflows)
        self.subJobTaskWidgets = {}
        self.paramsFileName = None
        # For special cases (eg workflows) do not draw input data frame
        self.excludeInputData = False
        self.paramsList = []
        self.widgetLookup = {}
        self.layoutMode = layoutMode
        if layoutMode is not None and layoutMode in ['TAB', 'FOLDER']:
            self.layoutMode = layoutMode
        else:
            self.layoutMode = str(CCP4Modules.PREFERENCES().TASK_WINDOW_LAYOUT)
        # Project is projectName as this is used by the CDataFileView most frequently
        self._projectId = projectId
        self._jobId = jobId
        self._jobNumber = None
        self._container = CCP4Container.CContainer()
        self.folderAttributes = CFolderAttributes()
        self.contextMenuWidget = QtWidgets.QMenu(self)
        self.widgetWithContextMenu = []
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        self.layout().setSpacing(0)
        if self.layoutMode == 'TAB':
            self.widget = CTabTaskWidget(parent=self)
            self.layout().addWidget(self.widget)
        else:
            self.scrollArea= QtWidgets.QScrollArea(self)
            self.layout().addWidget(self.scrollArea)
            self.widget = CFolderTaskWidget(parent=self)
            self.scrollArea.setWidget(self.widget)
            self.scrollArea.setWidgetResizable(1)
            self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            #print 'CTaskWidget self',self,'scrollArea',self.scrollArea,'widget',self.widget,'widget.parent',self.widget.parent(),self.widget.parent().parent()
        self.messageScrollArea= QtWidgets.QScrollArea(self)
        self.messageScrollArea.setMinimumHeight(40)
        self.messageScrollArea.setMaximumHeight(40)
        #self.messageScrollArea.setMinimumWidth(WIDTH)
        #self.messageScrollArea.setMaximumWidth(WIDTH)
        self.messageScrollArea.setObjectName('messageScrollArea')
        self.message = QtWidgets.QLabel(self)
        self.message.setObjectName('errorMessage')
        self.message.setWordWrap(True)
        self.messageScrollArea.setWidget(self.message)
        self.messageScrollArea.setWidgetResizable(1)
        self.layout().addWidget(self.messageScrollArea)
        self.messageStack = []

    def setMessage(self, text='', parameter=None):
        if text is None:
            text = ''
        self.message.setText(text)
        self.messageStack.append([parameter,text])
        #print 'setMessage',text,parameter,self.messageStack
    
    def unsetMessage(self, parameter=None):
        if len(self.messageStack) == 0:
            return
        if self.messageStack[-1][0] == parameter:
            del self.messageStack[-1]
        else:
            delFrom = None
            for ii in range(len(self.messageStack)-2,-1,-1):
                if self.messageStack[ii][0] == parameter:
                    delFrom = ii
            if delFrom is not None:
                #print 'unsetMessage deleting multi levels of stack',self.messageStack,delFrom
                del self.messageStack[delFrom:-1]
        #print 'unsetMessage',parameter,self.messageStack
        if len(self.messageStack) > 0:
            self.message.setText(self.messageStack[-1][1])

    def updateMessage(self, text='', parameter=None):
        for ii in range(len(self.messageStack)):
            if self.messageStack[ii][0] == parameter:
                if text is None:
                    text = ''
                self.messageStack[ii][1] = text
        if len(self.messageStack) > 0:
            self.message.setText(self.messageStack[-1][1])
        #print 'updateMessage',text,parameter,self.messageStack

    def setDefaultParameters(self):
        # Reimplement in tasks to set initial parameters
        pass

    def setDefaultFiles(self):
        # Reimplement in tasks to set initial files
        # This is called after setFollowJobId() so can override
        # and files set by the system
        pass

    def setParamsFileName(self, fileName):
        self.paramsFileName = fileName

    def setSelection(self, selection):
        self.selection = selection

    def setFollowJobId(self, jobId, force=True):
        self._container.guiAdmin.followFrom = jobId
        widget = self.getWidget('followFrom')
        if widget is not None:
            widget.updateInputFiles(jobId, self._projectId, force)
            widget.updateViewFromModel()

    def loadControlParameters(self,jobId):
        print('CTaskWidget.loadControlParameters', jobId)
        from core import CCP4Container
        defFile = PROJECTSMANAGER().makeFileName(jobId, 'PARAMS')
        if not os.path.exists(defFile):
            defFile1 = PROJECTSMANAGER().makeFileName(jobId, 'JOB_INPUT')
            if not os.path.exists(defFile1):
                print('Can not copy parameters from other job - parameters file ', defFile, 'does not exist')
                return
            else:
                defFile = defFile1
        tmpContainer = CCP4Container.CContainer()
        tmpContainer.loadDataFromXml(defFile, check=False, loadHeader=False)
        self._container.controlParameters.copyData(tmpContainer.controlParameters)
        for key in self._container.controlParameters.dataOrder():
            widget = self.getWidget(key)
            if widget is None:
                print('loadControlParameters no widget for', key)
            else:
                widget.updateViewFromModel()
                print('loadControlParameters updated', key)

    def isEditable(self):
        return self.folderAttributes.attribute('editable')

    def draw(self):
        rv = self.widget.draw()
        self.widget.show()
        return rv

    def visibleFolder(self):
        if self.layoutMode == 'TAB':
            return str(self.widget.tabText(self.widget.currentIndex()))
        else:
            return None

    def setVisibleFolder(self,title=None,index=None):
        if self.layoutMode != 'TAB': return
        if index is None:
            for ii in range(self.widget.count()):
                if title == str(self.widget.tabText(ii)):
                    index = ii
                    break
        if index is not None:
            self.widget.setCurrentIndex(index)

    def openFolder(self,folderFunction='general', title=None, toggle=[], toggleFunction=[], level=USER_NOVICE,
                   autoSelection={}, drawFolder=None, **kw):
        #print 'CtaskWidget.openFolder',folderFunction,self.excludeInputData
        if self.excludeInputData and folderFunction == 'inputData':
            self.widget.ignoreInput = True
            return
        else:
            self.widget.ignoreInput = False
        if title == None:
            if folderFunction == 'protocol':
                title = self.title()
            elif folderFunction == 'inputData':
                title = 'Input data'
            elif folderFunction == 'outputData':
                title = 'Output data'
        attribs = self.folderAttributes.allAttributes(folderFunction)
        attribs['level'] = level
        rv = self.widget.openFolder(folderFunction=folderFunction, title=title, toggle=toggle,
                                    toggleFunction=toggleFunction, attributes=attribs, drawFolder=drawFolder, keywords=kw)
        if len(autoSelection) > 0:
            self.autoGenerate(container=self.container.controlParameters, selection=autoSelection)
        return rv

    def closeFolder(self):
        self.widget.closeFolder()

    def createLine(self, definition=[], appendLine=None, toggle=[], toggleFunction=[]):
        return self.widget.createLine(definition=definition, appendLine=appendLine, toggle=toggle, toggleFunction=toggleFunction)

    def openStack(self, controlVar=None):
        return self.widget.openStack(controlVar=controlVar)

    def closeStack(self):
        self.widget.closeStack()

    def openTabFrame(self, title, toolTip=None):
        self.widget.openTabFrame(title=title, toolTip=toolTip)

    def closeTabFrame(self):
        self.widget.closeTabFrame()

    def createRadioGroup(self, label=None, itemList=[]):
        return self.widget.createRadioGroup(label=label, itemList=itemList)

    def createTitleLine(self, name='TITLE'):
        self.widget.createTitleLine(name=name)

    def openSubFrame(self, toggle=[], toggleFunction=[], frame=False, title=None, tip=''):
        #print 'CTaskWidget.openSubFrame',args
        return self.widget.openSubFrame(toggle=toggle, toggleFunction=toggleFunction, frame=frame, title=title, tip=tip)

    def closeSubFrame(self):
        self.widget.closeSubFrame()

    def setContainer(self, container=None):
        self._container = container

    def getContainer(self):
        return self._container

    container = property(getContainer, setContainer)

    def setProjectId(self, projectId=None):
        self._projectId = projectId

    def projectId(self):
        return self._projectId

    def setJobId(self, jobId):
        self._jobId = jobId
        self._jobNumber= PROJECTSMANAGER().db().getJobInfo(jobId = jobId, mode='jobnumber')

    def jobId(self):
        return self._jobId

    def jobNumber(self):
        return self._jobNumber

    def taskName(self):
        return self.TASKNAME

    def setTitle(self, title):
        self._title = title

    def title(self):
        if self._title is not None:
            return self._title
        elif hasattr(self, 'TASKTITLE'):
            return self.TASKTITLE
        else:
            return self.__class__.__name__

    def isEditor(self):
        return self.EDITOR

    def autoPopulateInput(self):
        return self.AUTOPOPULATEINPUT
      
    def getParams(self, paramValues={}):
        from core import CCP4Data
        for key, value in list(paramValues.items()):
            widget = self.getWidget(key)
            if widget is not None:
                if isinstance(value, CCP4Data.CData):
                    ret_value = widget.model.getDataObjects()
                else:
                    ret_value = widget.getValue()
                    paramValues[key] = ret_value

    def setDefaultParams(self):
        if self._container is not None:
            return 0
        self.setParams(self._container.getDataObjects())

    def setParams(self,paramValues={}):
        for key,value in list(paramValues.items()):
            widget = self.getWidget(key)
            if widget is not None:
                if isinstance(widget, CCP4Widgets.CViewWidget):
                    widget.setModel(value)
                else:
                    widget.setValue(value)

    def setProgramHelpFile(self, helpFile):
        self.programHelpFile = helpFile

    def setHelpFile(self, helpFile):
        self.helpFile = helpFile

    @QtCore.Slot(str,str,str,int,int)
    def populateContextMenu(self, name, helpFile, helpTarget, globalX, globalY):
        #print 'CTaskWidget.populateContextMenu',name,helpTarget,globalX,globalY
        from qtgui import CCP4GuiUtils
        self.widgetWithContextMenu = [name, helpFile, helpTarget]
        self.contextMenuWidget.clear()
        CCP4GuiUtils.populateMenu(self, self.contextMenuWidget, ['help'], default_icon='')
        self.contextMenuWidget.popup(QtCore.QPoint(globalX, globalY))

    def getActionDef(self,name=''):
        if name == 'help':
            return dict(text = 'Help', tip = 'Help', slot = self.help)

    def help(self):
        #print 'CTaskWidget.help',self.widgetWithContextMenu
        try:
            helpPath = os.path.join(os.environ['CCP4'], 'docs', 'sphinx', 'build', 'html')
        except:
            QtWidgets.QMessageBox.warning(None, self.windowTitle(), 'Can not access help files - CCP4 not setup')
        helpPath = os.path.join(helpPath, self.widgetWithContextMenu[1] + '.html')
        if not os.path.exists(helpPath):
            QtWidgets.QMessageBox.warning(None, self.windowTitle(), 'Help file not found: ' + helpPath)
            return
        WEBBROWSER().loadWebPage(helpPath)

    '''
    def updateModelFromView(self,textOnly=False):
        # Not convinced this is needed and sledge-hammer updateModelFromView() of some
        # widgets (eg CImportUnmerged) very bad as trashes gui-generalted data (such as cell params)
        return

    def updateModelFromView(self,textOnly=False):
        rv = CErrorReport()
        from qtgui.CCP4Widgets import CViewWidget
        from qtgui.CCP4ContainerView import CContainerView
        toggledFrame = {}
        for t in self.widget.toggleList:
          toggledFrame[t.target] = t.getTargetVisibility()
        widgetList = self.findChildren(CViewWidget)
        #print 'CTaskWidget.updateModelFromView',textOnly,widgetList
        for widget in widgetList:
          #print widget.model.objectName(),
          w = widget.parent()
          isVis = True
          while isVis and not isinstance(w,(CTaskWidget,CContainerView,QtWidgets.QDialog,QtWidgets.QMainWindow)):
            #print w,toggledFrame.get(w,True) ,'*',
            if not toggledFrame.get(w,True): isVis = False
            w = w.parent()
          #print 'isVis',isVis,widget.getValue()
          if isVis:
            try:
              if textOnly:
                widget.updateModelFromText()
              else:
                widget.updateModelFromView()
            except:
              rv.append(self.__class__,102,str(widget.objectName()),stack=False)
        return rv
    '''

    def updateModelFromView(self, textOnly=False):
        # Called when user clicks run - try to fix any text widget with focus that has not updated
        # the model since has not had the 'finished' signal from user hitting return or moving focus
        rv = CErrorReport()
        fW = QTAPPLICATION().focusWidget()
        if fW is not None:
            if isinstance(fW,CCP4Widgets.CLineEdit):
                fW.editingFinished.emit()
            elif isinstance(fW,CCP4Widgets.CTextEdit):
                fW.textChanged.emit()
        return rv
        '''
        w = fW.parent()
        while w is not None:
            if isinstance(w,(CCP4Widgets.CStringView,CCP4Widgets.CFloatView,CCP4Widgets.CIntView)):
                try:
                    print('CTaskWidget.updateModelFromView updating model from view for',w.model.objectName(),'widget value:',w.getValue(),'old model value:',w.model)
                    w.updateModelFromView()
                except:
                    pass
                return rv
            else:
                w = w.parent()
        '''

    def updateViewFromModel(self):
        rv = CErrorReport()
        for param in self.paramsList:
            try:
                widgetList = self.findWidget(param)
                #print 'CTaskWidget.updateViewFromModel',widget,widget.model.get0()
                for widget in widgetList:
                    widget.blockSignals(True)
                    widget.updateViewFromModel()
                    widget.blockSignals(False)
            except:
                rv.append(self.__class__, 101, str(widget.objectName()))
        #print 'CTaskWidget.updateViewFromModel errors:',rv.report()
        return rv

    def setup(self):
        return 0

    def createWidget(self, name=None, widgetQualifiers={}, helpTarget=None):
        model = self._container.find(name)
        if model is None:
            raise CException(self.__class__, 100, name)
        wQualifiers = {}
        wQualifiers.update(widgetQualifiers)
        wQualifiers.update(model.qualifiers())
        #print 'widgetQualifiers',name,wQualifiers
        widget = DATAMANAGER().widget(model=model, parentWidget=self, qualifiers=wQualifiers, name=name)
        #print 'CTaskLine widget',name,widget
        if widget is not None:
            widget.setObjectName(name)
            #modelTip = model.qualifiers('toolTip')
            #if modelTip is not NotImplemented:
            #  widget.setToolTip(modelTip)
            if model.qualifiers('toolTip') is not NotImplemented:
                widget.setToolTip(model.qualifiers('toolTip'))
            '''
            if widget.STRETCH > 0:
                self.layout().setStretch(self.layout().count(),widget.stretchFactor())
                doneStretch = True
            elif widgetQualifiers.get('charWidth',1) < 0:
                doneStretch = True
            '''
            if helpTarget is not None:   # KJS: Change to functools version 
                widget.contextMenuRequest.connect(functools.partial(self.populateContextMenu, name, self.programHelpFile, helpTarget))
        return widget

    def getWidget(self, name=None):
        return self.widgetLookup.get(name, [None])[0]

    def setMenuText(self, parameter=None, menuText=None):
        '''
        Set data object qualifier menuText from Python
        parameter is the name of the data object
        menuText can be a list which should be same length as the enumerators qualifier
        or it can be a dict with the enumerator values as keys
        '''
        dataObj = self._container.find(parameter)
        if dataObj is None:
            raise CException(self.__class__, 101, str(parameter))
        enumerators = dataObj.qualifiers('enumerators')
        #print 'setMenuText enumerators',parameter,dataObj,enumerators,len(enumerators),len(menuText)
        if isinstance(menuText,list):
            if len(menuText) != len(enumerators):
                raise CException(self.__class__, 102, str(parameter), stack=False)
            dataObj.setQualifier('menuText', menuText)
        elif isinstance(menuText, dict):
            enums = []
            menu = []
            unrecog = []
            for item in enumerators:
                if item in menuText:
                    enums.append(item)
                    menu.append(menuText[item])
                else:
                    unrecog.append(item)
            #print 'setMenuText new menu',enums,menu
            dataObj.setQualifiers({'enumerators' : enums, 'menuText': menu})
            if len(unrecog) > 0:
                return CException(self.__class__, 103, str(parameter)+' ' + str(unrecog))

    def setToolTip(self, parameter=None, toolTip=None):
        dataObj = self._container.find(parameter)
        if dataObj is None:
            raise CException(self.__class__, 104, str(parameter))
        dataObj.setQualifier('toolTip', toolTip)

    def setLabel(self, parameter=None, label=None):
        dataObj = self._container.find(parameter)
        if dataObj is None:
            raise CException(self.__class__, 105, str(parameter))
        dataObj.setQualifier('label', label)

    @QtCore.Slot()
    def validate(self):
        if self.isEditable():
            return self.widget.validate()
        else:
            return 0

    def resetJobCombos(self):
        return self.widget.resetJobCombos()

    def fix(self):
        # Dummy method to be reimplemented in sub-class
        # Called after user clicks run button and before validate()
        # Probably should be used mostly to ensure that if there are
        # possible alternate input files then only one has a set value
        # (this is to prevent unused files being recorded in database)
        self.doFix.emit()
        rv = CErrorReport()
        for key, win in list(self.subJobTaskWidgets.items()):
            rv.extend(win.fix())
        return rv

    def cootFix(self):
        from core import CCP4Modules
        from qtgui import CCP4FileBrowser
        path = str(CCP4Modules.PREFERENCES().COOT_EXECUTABLE)
        #print 'CTaskWidget.cootFix',path
        try:
            if  path is not None and os.path.exists(path):
                return CErrorReport()
        except:
            pass
        '''
        path = CCP4Utils.findCootPath()
        if path is not None:
          CCP4Modules.PREFERENCES().COOT_EXECUTABLE.set(fullPath=path)
          return CErrorReport()
        '''
        self.cootFixDialog = CCP4FileBrowser.CFileDialog(self, 'Find Coot Path', filters=[' (*)'], projectCombo=False)
        label = QtWidgets.QLabel("""Sorry - failed to find Coot. Please enter the Coot executable and then 'Run' again.\nBeware you are probably using a CCP4 nightly build that does not include Coot.\nThe Coot executable can also be set in Preferences.""",self)
        label.setStyleSheet("QLabel { font-weight: bold;  border: 2px solid} ")
        self.cootFixDialog.addWidget(label)
        self.cootFixDialog.selectFile.connect(self.handleCootFix)
        self.cootFixDialog.show()
        self.cootFixDialog.raise_()
        return CErrorReport(self.__class__, 103)

    @QtCore.Slot(str)
    def handleCootFix(self, filePath):
        #self.cootFixWidget.updateModelFromView()
        self.cootFixDialog.hide()
        self.cootFixDialog.deleteLater()
        if filePath is not None and os.path.exists(filePath):
            from core import CCP4Modules
            CCP4Modules.PREFERENCES().COOT_EXECUTABLE.set(filePath)
            CCP4Modules.PREFERENCES().save()
        else:
            pass

    def isValid(self):
        invalidList = []
        for key, widgetList in list(self.widgetLookup.items()):
            for widget in widgetList:
                if widget.isValid is None:
                    widget.validate()
                if widget.isValid is not None and not widget.isValid:
                    if getattr(widget, 'model', None) is not None:
                        if widget.model.objectName() != 'fileContent':
                            #print 'CTaskWidget.isValid',widget.model.objectName(),widget.isValid
                            invalidList.append(widget.model)
                    else:
                        invalidList.append(str(widget))
        # Apply to any sub-job widgets
        for key, win in list(self.subJobTaskWidgets.items()):
            invalidList.extend(win.isValid())
        #print 'CTaskWidget.isValid',invalidList
        return invalidList

    def taskValidity(self):
        return CErrorReport()

    def getScrollDisplacement(self):
        scrollArea = getattr(self, 'scrollArea', None)
        if scrollArea is None:
            return 0
        else:
            return self.scrollArea.verticalScrollBar().value()

    def setScrollDisplacement(self, value):
        scrollArea = getattr(self, 'scrollArea', None)
        if scrollArea is not None:
            self.scrollArea.verticalScrollBar().setValue(value)

    def autoGenerate(self, container=None, selection={}, subFrame=False):
        from core import CCP4File
        from core import CCP4Data
        from core import CCP4Container
        label = container.qualifiers('guiLabel')
        defn = container.qualifiers('guiDefinition')
        if defn is None or defn is NotImplemented:
            defn = {}
        if subFrame:
            frame = self.widget.openSubFrame(frame=True)
            #print 'subFrame',frame,container.__dict__['_qualifiers'],label,defn, defn.get('toggleParameter',None)
            if defn.get('toggleParameter', None) is not None:
                self.widget.setToggle(target=frame, parameter=defn['toggleParameter'],
                                      state=defn.get('toggleState', 'open'), values=defn['toggleValues'])
        if label is not NotImplemented:
            self.widget.createLine(definition=['advice', label])
        #print 'autoGenerate container',repr(container),defn
        includeParameters = selection.get('includeParameters', [])
        #print 'CTaskWidget.autoGenerate includeParameters',includeParameters
        includeWildcards = []
        for item in includeParameters:
            if item.count('*'):
                includeWildcards.append(item)
        excludeParameters = selection.get('excludeParameters', [])
        keyValues = selection.get('keyValues', {})
        #--------------------------------------------------------------------
        def matchesKeyValues(definition):
            if len(keyValues) == 0:
                return True
            for key,value in list(keyValues.items()):
                if key not in definition or definition[key] != value:
                    return False
            return True
        def matchesWildcard(name):
            for wildcard in includeWildcards:
                if re.match(wildcard, name):
                    return True
            return False
        #--------------------------------------------------------------------
        #print 'CTaskWidget.autoGenerate includeWildcards',includeWildcards
        for name in container.dataOrder():
            # If includeParameters is defined then name should be in that list or match a wildcard in that list
            # if includeParameters is not defined then name just needs to not be in excludeParameters list
            #print 'CTaskWidget.autoGenerate',name,len(includeParameters),name in includeParameters,matchesWildcard(name)
            if (len(includeParameters) > 0 and (name in includeParameters or matchesWildcard(name))) or \
                        (len(includeParameters) == 0 and name not in selection.get('excludeParameters', [])):
                model=getattr(container,name)
                if isinstance(model, CCP4Container.CContainer):
                    self.autoGenerate(model, selection=selection, subFrame=True)
                else:
                    label = model.qualifiers('guiLabel')
                    defn = model.qualifiers('guiDefinition')
                    if defn is None or defn is NotImplemented:
                        defn = {}
                    if label is None or label is NotImplemented:
                        label = name
                        #print 'autoGenerate',name,repr(model),label,defn
                    if matchesKeyValues(defn):
                        if isinstance(model, CCP4File.CDataFile):
                            line = self.widget.createLine(definition=['widget', name ])
                            widget = line.findChildren(QtWidgets.QWidget)[0]
                        elif isinstance(model, CCP4Data.CBoolean):
                            line = self.widget.createLine(definition=['widget', name, 'label' , label])
                            widget = line.findChildren(QtWidgets.QWidget)[0]
                        else:
                            line = self.widget.createLine(definition=['label', label, 'widget', name])
                            widget = line.findChildren(QtWidgets.QWidget)[1]
                            #widget.setMaximumWidth(WIDTH-116)
                        if defn.get('toggleParameter',None) is not None:
                            self.widget.setToggle(target=line, parameter=defn['toggleParameter'],
                                                  state=defn.get('toggleState', 'open'), values=defn['toggleValues'])
                        toolTip = model.qualifiers('toolTip')
                        if toolTip is not None and toolTip is not NotImplemented:
                            #print 'toolTip',widgetIndex,line.findChildren(QtWidgets.QWidget)
                            widget.setToolTip(toolTip)
        if subFrame:
            self.closeSubFrame()

    @QtCore.Slot(str)
    def handleClosingSubTaskWindow(self, jobName):
        if jobName in self.subJobTaskWidgets:
            fileNameSplit = os.path.split(PROJECTSMANAGER().makeFileName(jobId=self._jobId, mode='JOB_INPUT'))
            self.subJobTaskWidgets[jobName].saveToXml(fileName=os.path.join(fileNameSplit[0], jobName + '_' + fileNameSplit[1]))
            del self.subJobTaskWidgets[jobName]

    def slugify(self, value):
        """
        Normalizes string, converts to lowercase, removes non-alpha characters,
        and converts spaces to hyphens.
        """
        import unicodedata
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
        if sys.version_info >= (3,0):
            try:
                value = value.decode()
            except:
                pass
            value = str(re.sub('[^\w\s-]', '', value).strip().lower())
            value = str(re.sub('[-\s]+', '-', value))
        else:
            value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
            value = unicode(re.sub('[-\s]+', '-', value))
        
        return value
    
    def patchOutputFilePaths(self, jobInfo, fileName, projectDirectory=None):
        from core import CCP4File
        from core import CCP4Data
        from core import CCP4ModelData
        jobDirectory = os.path.split(os.path.normpath(fileName))[0]
        dataList = self.container.outputData.dataOrder()
        for objectName in dataList:
            dobj = self.container.outputData.find(objectName)
            partPath = os.path.join(jobInfo['jobnumber']+"_"+jobInfo['projectname']+"_"+objectName+"_"+jobInfo['taskname'])
            if sys.version_info >= (3,0):
                partPath = str(self.slugify(str(partPath)))
            else:
                partPath = str(self.slugify(unicode(partPath)))
            if isinstance(dobj,CCP4File.CDataFile):
                if isinstance(dobj,CCP4ModelData.CPdbDataFile) and dobj.contentFlag.isSet() and int(dobj.contentFlag) == 2:
                    fullPath = os.path.join(jobDirectory, partPath+"." + dobj.fileExtensions()[1])
                else:
                    fullPath = os.path.join(jobDirectory, partPath+"." + dobj.fileExtensions()[0])
                #MN 2020-07-09
                #The digestion of a fullpath into a relpath in setFullPath is done only if checkDb is True
                #I don't know why that was not the case, so I am switching the default behaviour to be
                #manually setting elements of the files location, failing over to checkDb False
                if projectDirectory is not None:
                    fileDirectory, fileName = os.path.split(fullPath)
                    dobj.relPath.set(os.path.relpath(fileDirectory, projectDirectory))
                    dobj.baseName.set(fileName)
                    dobj.project.set(jobInfo['projectid'])
                else:
                    try:
                        dobj.setFullPath(fullPath, checkDb=True)
                    except:
                        dobj.setFullPath(fullPath, checkDb=False)
            elif isinstance(dobj,CCP4Data.COutputFileList):
                #Empty the current output file list array
                while len(dobj) > 0: dobj.remove(dobj[-1])
                for iItem in range(dobj.qualifiers()['listMaxLength']):
                    dobj.append(dobj.makeItem())
                    fullPath = os.path.join(jobDirectory, "{}_{}.{}".format(partPath, str(iItem),dobj[-1].fileExtensions()[0]))
                    #MN 2020-07-09
                    #The digestion of a fullpath into a relpath in setFullPath is done only if checkDb is True
                    #I don't know why that was not the case, so I am switching the default behaviour to be
                    #manually setting elements of the files location, failing over to checkDb False
                    if projectDirectory is not None:
                        fileDirectory, fileName = os.path.split(fullPath)
                        dobj[-1].relPath.set(os.path.relpath(fileDirectory, projectDirectory))
                        dobj[-1].baseName.set(fileName)
                        dobj[-1].project.set(jobInfo['projectid'])
                    else:
                        try:
                            dobj[-1].setFullPath(fullPath, checkDb=True)
                        except:
                            dobj[-1].setFullPath(fullPath, checkDb=False)

    def saveToXml(self, fileName=None, jobInfo={}):
        # Assume the updateModelFromView() has been called earlier by project viewer before the
        # validation of input
        #self.updateModelFromView()
        from core import CCP4File
        from core import CCP4Modules
        from core import CCP4TaskManager
        if fileName is None:
            if self.paramsFileName is not None:
                fileName= self.paramsFileName
            else:
                fileName = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self._jobId, mode='JOB_INPUT')
        jobInfo = {}
        if self._jobId is not None:
            try:
                jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self._jobId, mode=['taskname', 'jobnumber', 'projectname', 'status', 'projectid'])
            except:
                pass
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        cHeader = self.container.getHeader()
        if cHeader is not None:
            f.header.set(cHeader)
        f.header.setCurrent()
        f.header.function.set('PARAMS')
        f.header.jobId.set(self._jobId)
        f.header.projectId.set(self._projectId)
        try:
            f.header.pluginName = CCP4TaskManager.TASKMANAGER().getTaskData(self.taskName())['taskName']
        except:
            print('Error getting plugin name for taskName',self.taskName())
        if jobInfo.get('jobnumber', None) is not None:
            f.header.jobNumber.set(jobInfo['jobnumber'])
        if jobInfo.get('projectname', None) is not None:
            f.header.projectName.set(jobInfo['projectname'])
                
        #MN set output file paths here, where project, and jobNumber are all known
        #MNDoing a single lookup of projectDirectory can speed up the patching of file names later
        projectDirectory=CCP4Modules.PROJECTSMANAGER().db().getProjectDirectory(jobInfo['projectid'])
        self.patchOutputFilePaths(jobInfo, fileName, projectDirectory=projectDirectory)
        
        bodyEtree = self.container.getEtree()
        f.saveFile(bodyEtree=bodyEtree)
        #print 'CTaskWidget.saveTofile DONE',fileName,jobInfo

    def handleLaunchedJob(self, jobId=None, status=None, taskWidget=None):
        print('Dummy method for reimplementation in task that has launched a popout task', jobId, status, taskWidget)
        pass

    def connectDataChanged(self, name, handle):
        obj = self.container.find(name)
        if obj is None:
            return
        obj.dataChanged.connect(handle)

    def findWidget(self, name):
        return self.widgetLookup.get(name, None)


#---------------------------------------------------------------------
class CFolderTaskWidget(QtWidgets.QFrame, CTaskWidgetDrawMethods):
#---------------------------------------------------------------------

    MARGIN = 4
    ERROR_CODES = {101 : {'description' : 'Error updating GUI widget with model value'},
                   105 : {'description' : 'Internal error handing file - no task container'}}

    def __init__(self, parent=None):
        QtWidgets.QFrame.__init__(self,parent)
        CTaskWidgetDrawMethods.__init__(self)
        layout = QtWidgets.QVBoxLayout()
        #layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        layout.setContentsMargins(CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN)
        layout.setSpacing(CFolderTaskWidget.MARGIN)
        self.setLayout(layout)

    def parentTaskWidget(self):
        return self.parent().parent().parent()

    def draw(self):
        self.myException = CErrorReport()
        self.parentTaskWidget().drawContents()
        self.parentTaskWidget().updateViewFromModel()
        self.closeFolder()
        self.layout().addStretch(5)
        '''
        for ii in range(0,self.layout().count()):
          w = self.layout().itemAt(ii).widget()
          print 'CTaskWidget.draw',ii,w
          if w is not None:
            w.contents.show()
            w.show()
        '''
        #print 'CTaskWidget.finishDraw',self.isVisible()
        self.parentTaskWidget().setDefaultParams()
        self.applyToggles()
        e = CErrorReport()
        #print 'CTaskWidget.draw', len(e)
        e.extend(self.myException, stack=False)
        del self.myException
        return e

    def openFolder(self, folderFunction='general', title='CCP4 Task Folder', toggle=[], toggleFunction=[],
                   attributes={}, keywords={}, drawFolder=None,**kw):
        if self.subFrame is not None:
            self.closeSubFrame()
        self.closeFolder()
        self.currentFolder = CTaskFolder(self, folderFunction=folderFunction, title=title, toggle=toggle,
                                         toggleFunction=toggleFunction, attributes=attributes)
        self.currentFolderLayout = QtWidgets.QVBoxLayout()
        self.currentFolderLayout.setSpacing(CFolderTaskWidget.MARGIN)
        self.currentFolderLayout.setContentsMargins(CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN, CFolderTaskWidget.MARGIN)
        self.currentFolder.folderToggled.connect(functools.partial(self.folderToggled, self.currentFolder))
        if folderFunction == 'inputData':
            self.createJobTitle(followFrom=keywords.get('followFrom', True))
        return self.currentFolder

    def closeFolder(self):
        if self.subFrame is not None:
            self.closeSubFrame()
        if self.tabWidget is not None:
            self.closeTabFrame()
        if self.currentFolder is not None:
            if self.currentFolder.folderFunction == 'inputData':
                self.createPatchFrame()
            self.currentFolderLayout.addStretch(1)
            self.currentFolder.setContentsLayout(self.currentFolderLayout)
            self.layout().addWidget(self.currentFolder)
        self.currentFolder = None
        self.currentFolderLayout = None

    @QtCore.Slot('CTaskFolder')
    def folderToggled(self, folder):
        self.layout().update()

    def getWidget(self, name):
        # findChild only finds stuff in top tab??
        print("CCP4TaskManager.CFolderTaskWidgetgetWidget")
        return self.widgetLookup.get(name,[None])[0]

    def getFolderOpenStatus(self):
        status = []
        for i in range(self.layout().count()):
            w = self.layout().itemAt(i).widget()
            #print 'getFolderOpenStatus',i,w
            if w is not None:
                status.append(w.isOpen())
        return status

    def setFolderOpenStatus(self, openStatus):
        for i in range(min(len(openStatus), self.layout().count())):
            w = self.layout().itemAt(i).widget()
            if openStatus[i]:
                w.openFolder()
            else:
                w.closeFolder()

#---------------------------------------------------------------------
class CToggle(QtCore.QObject):
#---------------------------------------------------------------------
    ERROR_CODES = {101 : {'description' : 'Can not set up toggle widget undefined for'},
                   102 : {'description' : 'Can not set up toggle model undefined for'},
                   103 : {'description' : 'Can not change visibility, widget undefined for'},
                   104 : {'description' : 'Can not change visibility, model undefined for'}}

    def __init__(self, parent, target=None, parameter=None, state='open', values=[], parameterList=[], toggleFunction=None):
        QtCore.QObject.__init__(self, parent)
        #print 'CToggle.__init__',parameter,parent,target
        self.target = target
        self.parameter = parameter
        self.state = state
        self.values = values
        self.parameterList = []
        self.parameterList.extend(parameterList)
        self.toggleFunction = toggleFunction
        self.connected = False

    def makeConnection(self):
        container = self.parent().parentTaskWidget().container
        if container is None:
            return
        if self.parameter is not None:
            obj = container.find(self.parameter)
            if obj is None:
                return
            obj.dataChanged.connect(self.handleSignal)
        else:
            for param in self.parameterList:
                obj = container.find(param)
                if obj is not None:
                    obj.dataChanged.connect(self.handleSignal)
        self.connected = True

    @QtCore.Slot()
    def handleSignal(self):
        if self.parameter is not None:
            self.setTargetVisibility()
        else:
            self.applyVisibilityFunction()

    def setTargetVisibility(self):
        '''
        widget = self.parent().findChild(QtWidgets.QWidget,self.parameter)
        #print 'setTargetVisibility',self.parent(),self.parameter,widget
        if widget is None:
          raise CException(self.__class__,103,self.parameter)    
        if widget.model is None:
          raise CException(self.__class__,104,self.parameter)
        '''
        obj = self.parent().parentTaskWidget().container.find(self.parameter)
        if obj is None:
            return
        value = obj.get()
        #print 'CToggle.setTargetVisibility',self.parameter,value,self.values,self.values.count(value),self.state
        #if isinstance(self.target,CTaskLine) or isinstance(self.target,CSubFrame) or isinstance(self.target,CTaskFolder):
        if 1:
            if self.values.count(value):
                if self.state == 'open':
                    self.target.show()
                else:
                    self.target.hide()
            else:
                if self.state == 'open':
                    self.target.hide()
                else:
                    self.target.show()

    def getTargetVisibility(self):
        try:
            obj = self.parent().parentTaskWidget().container.find(self.parameter)
        except:
            return True
        if obj is None:
            return True
        value = obj.get()
        if 1:
            if self.values.count(value):
                return self.state == 'open'
            else:
                return self.state != 'open'

    def applyVisibilityFunction(self):
        vis = self.toggleFunction()
        #print 'applyVisibilityFunction',vis
        if vis:
            self.target.show()
        else:
            self.target.hide()


#---------------------------------------------------------------------
class CTaskFolder(CCP4Widgets.CFolder):
#---------------------------------------------------------------------

    def __init__(self,parent=None, folderFunction='general', title='Folder', copen=1, toggle=[], toggleFunction=[], attributes={}):
        self.folderFunction = folderFunction
        if ['protocol'].count(folderFunction):
            titleBar = 0
        else:
            titleBar = 1
        #print 'CTaskFolder.__init__',title
        self._attributes = attributes
        CCP4Widgets.CFolder.__init__(self, parent, title, copen, titleBar, toggle=toggle)
        from qtgui import CCP4StyleSheet
        self.setTitleColour(CCP4StyleSheet.LOWLIGHTCOLOUR)
        if self.titleBar is not None:
            self.titleBar.setMaximumHeight(30)

    def setupUpdateStatus(self, *args):
        pass


#---------------------------------------------------------------------
class CSubFrame(QtWidgets.QFrame):
#---------------------------------------------------------------------

    def __init__(self,parent):
        QtWidgets.QFrame.__init__(self,parent)

#---------------------------------------------------------------------
class CTaskLine(QtWidgets.QFrame):
#---------------------------------------------------------------------
    ERROR_CODES = {101 : {'description' : 'No data container for task line'},
                   102 : {'description' : 'No model found for task line item'},
                   103 : {'description' : 'Error creating widget for task line item'}}
    MARGIN = 0

    def __init__(self, parent=None):
        QtWidgets.QFrame.__init__(self, parent)
        self.setObjectName('taskLine')
        layout = QtWidgets.QHBoxLayout()
        layout.setSpacing(CTaskLine.MARGIN)
        layout.setContentsMargins(CTaskLine.MARGIN, CTaskLine.MARGIN, CTaskLine.MARGIN, CTaskLine.MARGIN)
        self.setLayout(layout)
        self.tip = None
        self.helpTarget = ''

    def parentTaskWidget(self):
        return self.parent().parentTaskWidget()

    def addWidget(self, widget):
        lastItem = self.layout().itemAt(self.layout().count()-1)
        #print 'CTaskLine.addWidget lastItem',lastItem
        if lastItem is not None and isinstance(lastItem, QtWidgets.QSpacerItem):
            self.layout().insertWidget(self.layout().count()-1,widget)
        else:
            self.layout().addWidget(widget)

    def draw(self, definition=[], dataContainer=None, setupFolder=None, attributes={}):
        #print 'CTaskLine.draw',definition
        doneStretch = False
        if dataContainer is None:
            raise CException(self.__class__, 101)
        myException = CErrorReport()
        widgetQualifiers = {}
        widgetQualifiers.update(attributes)
        #print 'CTaskLine.draw widgetQualifiers',widgetQualifiers
        pDef = -1
        definition.append('')
        ifFullLine = False      # KJS : This loop needs looked at. Investigate later.
        while (pDef < len(definition) - 2):
            widget = None
            pDef = pDef + 1
            if definition[pDef] in ['message', 'tip']:
                pDef = pDef + 1
                self.tip = definition[pDef]
            if definition[pDef] == 'help':
                pDef = pDef + 1
                self.helpTarget = definition[pDef]
            if definition[pDef] == 'label':
                pDef = pDef + 1
#FIXME - SJM 30/05/2017 - I have no idea why the colour is hardwired. It means that setEnabled(False) shows no "greying out" effect.
                #lab = QtWidgets.QLabel('<FONT color=black>' + definition[pDef] + '</FONT>')
                lab = QtWidgets.QLabel('<span>'+definition[pDef]+'</span>')
                if self.tip is not None: lab.setToolTip(self.tip)
                self.layout().addWidget(lab)
            if definition[pDef] == 'subtitle':
                pDef = pDef + 1
                subtitle = QtWidgets.QLabel(definition[pDef])
                pDef = pDef + 1
                subtitle.setToolTip ( '<FONT color=black>' + definition[pDef] + '</FONT>' )
                subtitle.setObjectName ( 'subtitle' )
                self.layout().addWidget ( subtitle )
                if self.tip is not None: subtitle.setToolTip(self.tip)
            if definition[pDef] in [ 'advice','warning']:
                pDef = pDef + 1
                #label = QtWidgets.QLabel('<FONT color=black>' + definition[pDef] + '</FONT>')
                label = QtWidgets.QLabel('<span>'+definition[pDef]+'</span>')
                label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
                self.layout().setContentsMargins(CFolderTaskWidget.MARGIN,CFolderTaskWidget.MARGIN,CFolderTaskWidget.MARGIN,CFolderTaskWidget.MARGIN)
                self.layout().setSpacing(CFolderTaskWidget.MARGIN)
                if definition[pDef-1] == 'advice':
                    label.setObjectName('italic')
                    self.layout().addWidget(label)
                else:
                    label.setObjectName('warning')
                    label.setWordWrap(True)
                    self.layout().addWidget(label)
                if self.tip is not None: lab.setToolTip(self.tip)
            if definition[pDef] == 'spacing':
                pDef = pDef + 1
                self.layout().addSpacing( definition[pDef])
                ifFullLine = True
            if definition[pDef] == 'stretch':
                self.layout().addStretch( 1)
                ifFullLine = True
            if definition[pDef] == 'launchButton':
                pDef = pDef + 1
                taskName = definition[pDef]
                widget = QtWidgets.QPushButton(TASKMANAGER().getTaskAttribute(taskName,'TASKTITLE'))
                widget.setObjectName(taskName)
                if self.tip is not None:
                    widget.setToolTip(self.tip)
                widget.clicked.connect(functools.partial(self.parentTaskWidget().launchJobRequestSignal.emit,taskName,{'launchJobId':self.parentTaskWidget().jobId()}))
                self.layout().addWidget(widget)
            if definition[pDef] == 'widget':
                pDef = pDef + 1
                while definition[pDef][0] == '-':
                    widgetQualifiers[definition[pDef][1:]] = definition[pDef+1]
                    pDef = pDef + 2
                model = dataContainer.find(definition[pDef])
                #print 'CTaskLine.draw',model,definition[pDef]
                #if len(widgetQualifiers)>0: print 'CTaskLine.draw widgetQualifiers',definition[pDef],widgetQualifiers
                if model is None:
                    myException.append(self.__class__, 102, definition[pDef], stack=False)
                else:
                    if not (self.parent().parent().paramsList.count(definition[pDef])):
                        self.parent().parent().paramsList.append(definition[pDef])
                    widget = None
                    qualifiers = {}
                    qualifiers.update(model.qualifiers())
                    qualifiers.update(widgetQualifiers)
                    #print 'widgetQualifiers',definition[pDef],model.qualifiers(),widgetQualifiers
                    if DEVELOPER():
                        widget = DATAMANAGER().widget(model=model, parentWidget=self.parent(), qualifiers=widgetQualifiers, name=definition[pDef])
                    else:
                        try:
                            widget = DATAMANAGER().widget(model=model, parentWidget=self.parent(), qualifiers=widgetQualifiers, name=definition[pDef])
                        except CException as e:
                            e.appendDetails(definition[pDef])
                            myException.extend(e)
                        except:
                            myException.append(self.__class__, 103, definition[pDef])
                    #print 'CTaskLine widget',definition[pDef],widget,widget.STRETCH
                    #widget = dataContainer.widget(name=par,parentWidget=self,widgetQualifiers=widgetQualifiers)
                    if widget is not None:
                        widget.setObjectName(definition[pDef])
                        #modelTip = model.qualifiers('toolTip')
                        #if modelTip is not NotImplemented:
                        #  widget.setToolTip(modelTip)
                        if self.tip is not None:
                            widget.setToolTip(self.tip + ' (' + str(definition[pDef]) + ')')
                        else:
                            tip = model.qualifiers('toolTip')
                            if tip is not NotImplemented and tip is not None:
                                widget.setToolTip(model.qualifiers('toolTip') + ' (' + str(definition[pDef]) + ')')
                            else:
                                widget.setToolTip(str(definition[pDef]))
                        self.layout().addWidget(widget)
                        if widget.STRETCH > 0:
                            #print 'CTaskLine.draw stretching',model.objectName()
                            self.layout().setStretch(self.layout().count() - 1, widget.STRETCH)
                            doneStretch = True
                        elif widgetQualifiers.get('charWidth', 1) < 0:
                            doneStretch = True
                        tW = self.parentTaskWidget()    # KJS: Change to functools version
                        if hasattr(widget,"contextMenuRequest"):
                            widget.contextMenuRequest.connect(functools.partial(tW.populateContextMenu, definition[pDef], tW.programHelpFile, self.helpTarget))
                        if definition[pDef] not in tW.widgetLookup:
                            tW.widgetLookup[definition[pDef]] = []
                        tW.widgetLookup[definition[pDef]].append(widget)
                        if setupFolder is not None:
                            setupFolder(model)
                        if isinstance(widget,CCP4Widgets.CComplexLineWidget):
                            ifFullLine = True
            if definition[pDef] == 'format':
                pass
            '''
            if ['toggle','toggle_display'].count(definition[pDef]):
                if pDef+3<len(definition) and ['open','close'].count(definition[pDef+2]):
                  self.parent().setToggle(self,definition[pDef+1],definition[pDef+2],definition[pDef+3])
                  pDef = pDef+3
            '''
        if not ifFullLine and not doneStretch:
            self.layout().addStretch(2)
        #if doneStretch:
        #  self.layout().addSpacing(1)
        #self.layout().addStretch(0.1)
        #else:
        #  self.layout().addStretch(2)
        return myException


#---------------------------------------------------------------------
class CTabTaskWidget(QtWidgets.QTabWidget, CTaskWidgetDrawMethods):
#---------------------------------------------------------------------

    MARGIN = 0
    ERROR_CODES = {101 : {'description' : 'Error updating GUI widget with model value'}}

    def __init__(self, parent=None):
        QtWidgets.QTabWidget.__init__(self, parent)
        self.currentChanged[int].connect(self.handleTabChanged)
        #print 'CTabTaskWidget.tabBar',QtWidgets.QTabWidget.tabBar(self)
        QtWidgets.QTabWidget.tabBar(self).setUsesScrollButtons(True)
        CTaskWidgetDrawMethods.__init__(self)

    def handleTabChanged(self, tabIndex):
        #print 'CTabTaskWidget.handleTabChanged',tabIndex,self.widget(tabIndex)._drawFolder
        if self.widget(tabIndex)._drawFolder is None:
            return
        self.parentTaskWidget().updateModelFromView()
        self.currentFolder = self.widget(tabIndex)
        self.currentFolderLayout = self.widget(tabIndex).frame.layout()
        self.myException = CErrorReport()
        try:
            self.widget(tabIndex)._drawFolder()
            print('DONE drawing folder', self.tabText(tabIndex))
        except Exception as e:
            print('FAILED drawing folder', self.tabText(tabIndex))
            print(e)
            traceback.print_exc()
        if self.currentFolderLayout is not None:
            self.currentFolderLayout.addStretch(10)
        self.widget(tabIndex)._drawFolder = None
        if len(self.myException) > 0:
            self.myException.report()
        self.currentFolder = None
        self.currentFolderLayout = None
        self.parentTaskWidget().updateViewFromModel()
        self.applyToggles()

    def parentTaskWidget(self):
        return self.parent()

    def draw(self):
        self.myException = CErrorReport()
        self.parentTaskWidget().drawContents()
        self.parentTaskWidget().updateViewFromModel()
        self.closeFolder()
        self.parentTaskWidget().setDefaultParams()
        self.applyToggles()
        e = CErrorReport()
        #print 'CTaskWidget.draw',len(e)
        e.extend(self.myException, stack=False)
        del self.myException
        return e

    def openFolder(self, folderFunction='general', title='CCP4 Task Folder', toggle=[], toggleFunction=[],
                   attributes={}, keywords={}, drawFolder=None, **kw):
        if self.subFrame is not None:
            self.closeSubFrame()
        self.closeFolder()
        self.currentFolder = CTabFrame(self, folderFunction=folderFunction, title=title, attributes=attributes, drawFolder=drawFolder)
        if len(toggle) > 0:
            if len(toggle) == 1:
                self.setToggle(self.currentFolder, toggle[0], 'open', [True])
            elif len(toggle) == 2:
                self.setToggle(self.currentFolder, toggle[0], toggle[1], [True])
            else:
                self.setToggle(self.currentFolder, toggle[0], toggle[1], toggle[2])
        if len(toggleFunction) > 0:
            self.setToggle(self.currentFolder, toggleFunction=toggle[0], parameterList=toggle[1])
        self.currentFolderLayout = QtWidgets.QVBoxLayout()
        self.currentFolderLayout.setSpacing(MARGIN)
        self.currentFolderLayout.setContentsMargins(MARGIN, MARGIN, MARGIN, MARGIN)
        self.currentFolder.frame.setLayout(self.currentFolderLayout)
        if folderFunction == 'inputData':
            self.createJobTitle(followFrom=keywords.get('followFrom', True))
        return self.currentFolder

    def closeFolder(self):
        if self.subFrame is not None:
            self.closeSubFrame()
        if self.tabWidget is not None:
            self.closeTabFrame()
        if self.currentFolder is not None:
            if self.currentFolder.folderFunction == 'inputData':
                self.createPatchFrame()
            self.currentFolderLayout.addStretch(1)
            self.currentFolder.frame.show()
            #self.currentFolder.setLayout(self.currentFolderLayout)
            self.addTab(self.currentFolder,self.currentFolder.title)
            self.setTabToolTip(self.count()-1,self.currentFolder.title)
            if ['inputData'].count(self.currentFolder.folderFunction) != 0:
                self.currentFolder.updateStatus()
        self.currentFolder = None
        self.currentFolderLayout = None

    def getWidget(self,name):
        # findChild only finds stuff in top tab??
        return self.parent().getWidget(name)

    def validate(self):
        # This is probably broken!!
        totInvalid = 0
        for i in range(self.count()):
            frame = self.widget(i)
            nInvalid = 0
            widgetList = frame.findChildren(CCP4Widgets.CViewWidget)
            #print 'CTabFolder.validate widgetList',widgetList
            for widget in widgetList:
                if widget.validate() is not None:
                    nInvalid + widget.validate()
            #print 'CTabFolder.validate nInvalid',nInvalid
            if nInvalid > 0:
                print('INVALID FRAME')
                frame.setTabColour('red')
                totInvalid = totInvalid + nInvalid
            else:
                frame.setTabColour('black')
        return totInvalid

    def resetJobCombos(self):
        #print 'CTabTaskWidget.resetJobCombos'
        for i in range(self.count()):
            frame = self.widget(i)
            widgetList = frame.findChildren(CCP4Widgets.CViewWidget)
            for widget in widgetList:
                if isinstance(widget, CCP4Widgets.CDataFileView):
                    widget.loadJobCombo()
                    widget.updateJobCombo()
                    widget.updateModelFromView()
      
#---------------------------------------------------------------------
class CTabFrame(QtWidgets.QFrame):
#---------------------------------------------------------------------

    def __init__(self, parent=None, folderFunction='general', title='Folder', attributes={}, drawFolder=None):
        QtWidgets.QFrame.__init__(self, parent)
        self.setObjectName(title)
        self._drawFolder = drawFolder
        self.title = title
        self.folderFunction = folderFunction
        self._attributes = attributes
        self.dataObjects = {}
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(0, 0, 0, 0)
        self.layout().setSpacing(0)
        self.scrollArea= QtWidgets.QScrollArea(self)
        self.layout().addWidget(self.scrollArea)
        self.frame = QtWidgets.QFrame()
        self.scrollArea.setWidget(self.frame)
        self.scrollArea.setWidgetResizable(1)
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)

    def hide(self):
        tab = self.parent().parent()
        indx = tab.indexOf(self)
        if tab.currentIndex() == indx:
            tab.setCurrentIndex(0)
        tab.setTabEnabled(indx, False)

    def show(self):
        tab = self.parent().parent()
        indx = tab.indexOf(self)
        tab.setTabEnabled(indx, True)

    def setTabColour(self,colour):
        tab = self.parent().parent()
        indx = tab.indexOf(self)
        tab.tabBar().setTabTextColor(indx, QtGui.QColor(colour))

    def updateStatus(self, dataObjectName=None):
        #print 'CTabFrame.updateStatus',dataObjectName
        allSet = True
        for key, obj in list(self.dataObjects.items()):
            #print 'CTabFrame.updateStatus testing',key,obj.validity().report()
            if obj.validity(obj.get()).maxSeverity() > SEVERITY_WARNING:
                allSet = False
                break
        if allSet:
            self.setTabColour('black')
        else:
            self.setTabColour('red')

#---------------------------------------------------------------------
class CStackedWidget(QtWidgets.QStackedWidget):
#---------------------------------------------------------------------
    # This did not work - nothing displayed despite all diagnostic output correct

    def __init__(self, parent=None, controlVar=None):
        self.controlVar = controlVar
        QtWidgets.QStackedWidget.__init__(self, parent)
        controlWidget = self.parent().getWidget(self.controlVar)
        #print 'CStackedWidget.__init__',controlWidget.model
        if controlWidget is not None:
            controlWidget.model.dataChanged.connect(self.update)

    def parentTaskWidget(self):
        parent = self.parent()
        while 1:
            if isinstance(parent, CTaskWidget):
                return parent
            if isinstance(parent, QtWidgets.QMainWindow):
                return None
            parent=parent.parent()

    @QtCore.Slot()
    def update(self):
        # I'd not expected to have to go up so far - has Qt reparented?
        controlWidget = self.parentTaskWidget().getWidget(self.controlVar)
        #print 'CStackedWidget.update',self.controlVar,controlWidget
        if controlWidget is None:
            return
        value = controlWidget.model.get()
        valueList = controlWidget.model.qualifiers('enumerators')
        #print 'CStackedWidget.update',value,valueList
        if value in valueList:
            indx = valueList.index(value)
            #print 'CStackedWidget.update',value,indx
            self.setCurrentIndex(indx)


#---------------------------------------------------------------------
class CStackedLayout(QtWidgets.QStackedLayout):
#---------------------------------------------------------------------

    def __init__(self, parent=None, controlVar=None):
        self.controlVar = controlVar
        QtWidgets.QStackedLayout.__init__(self)
        controlWidget = parent.getWidget(self.controlVar)
        #print 'CStackedLayout.__init__',controlWidget.model
        if controlWidget is not None:
            controlWidget.model.dataChanged.connect(self.update)

    def parentTaskWidget(self):
        parent = self.parentWidget()
        while 1:
            if isinstance(parent, CTaskWidget):
                return parent
            if isinstance(parent, QtWidgets.QMainWindow):
                return None
            parent = parent.parent()

    @QtCore.Slot()
    def update(self):
        # I'd not expected to have to go up so far - has Qt reparented?
        #print 'CStackedLayout.update',self.parentWidget()
        controlWidget = self.parentTaskWidget().getWidget(self.controlVar)
        #print 'CStackedLayout.update',self.controlVar,controlWidget
        if controlWidget is None:
            return
        value = controlWidget.model.get()
        valueList = controlWidget.model.qualifiers('enumerators')
        #print 'CStackedWidget.update',value,valueList
        if value in valueList:
            indx = valueList.index(value)
            #print 'CStackedLayout.update',value,indx
            self.setCurrentIndex(indx)


#===========================================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testTaskWidget)
    return suite

def runAllTests():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

#-------------------------------------------------------------------
class CSummatTask(CFolderTaskWidget):
#-------------------------------------------------------------------
# Subclass CTaskWidget to give specific task window

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent=None)

    def drawContents(self):
        #self.setProgramHelpFile('fft')
        #self.openFolder(folderFunction='protocol')
        self.openFolder(title='Test folder')
        #self.createTitleLine()
        self.createLine(['widget', 'nCycles', 'label', 'cycles with cutoff', 'widget', 'cutoff'])
        self.createLine(['label', 'Range of gubbins', 'widget', 'range' ] )


#class testTaskWidget(unittest.TestCase):
class testTaskWidget():

    def __init__(self):
        self.setUp()

    def setUp(self):
        from core.CCP4Modules import QTAPPLICATION
        self.app = QTAPPLICATION()
        self.window = QtWidgets.QMainWindow()

    def test1(self):
        from core.CCP4Container import CContainer
        summat = CContainer(parent=self.app, definitionFile='/Users/lizp/Desktop/dev/ccp4i2/sandpit/summat.contents.xml')
        self.assertEqual(len(summat.CONTENTS), 4, 'Container - Wrong content length')
        self.assertEqual(summat.gubbins.startValue, 12.0, 'Container - sub-container CFloat wrong initial value')

    def test2(self):
        from core.CCP4Container import CContainer
        import sys
        self.container = CContainer(parent=self.app, definitionFile='/Users/lizp/Desktop/dev/ccp4i2/sandpit/summat.contents.xml')
        self.task = CSummatTask(self.window)
        self.task.setContainer(self.container)
        self.task.draw()
        self.window.setCentralWidget(self.task)
        self.window.show()
        sys.exit(self.app.exec_())

    def tearDown(self):
        self.app.quit()
