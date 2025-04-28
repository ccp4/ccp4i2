"""
Copyright (C) 2013 STFC
Liz Potterton July 2013 - create and manage custom tasks
"""

import copy
import functools
import os

from PySide2 import QtCore, QtGui, QtWidgets

from . import CCP4CustomisationGui, CCP4Widgets
from ..core import CCP4CustomTaskManager
from ..core import CCP4Data
from ..core.CCP4ErrorHandling import Severity
from ..core.CCP4Modules import CUSTOMTASKMANAGER, WEBBROWSER
from ..core.CCP4WarningMessage import warningMessage


def openGui():
  if CCustomTaskManagerGui.insts is None:
    CCustomTaskManagerGui.insts = CCustomTaskManagerGui()
  CCustomTaskManagerGui.insts.show()
  CCustomTaskManagerGui.insts.raise_()


class CCustomTaskManagerGui(CCP4CustomisationGui.CCustomisationGui):

  insts = None


  def __init__(self,parent=None):
    CCP4CustomisationGui.CCustomisationGui.__init__(self,parent=parent,mode='customtask',title='Custom Task Manager')
    self.createWidget= None
    

  def manager(self):
    return CUSTOMTASKMANAGER()

  def handleNew(self):
    if self.createWidget is None:
      self.createWidget = CCreateCustomTaskDialog(self)
      self.customTaskCreated.connect(self.customListView.populate)
    else:
      self.createWidget.updateViewFromModel()
      
    self.createWidget.show()
    self.createWidget.raise_()


  def handleEdit(self,selected=None):
    if selected is None: selected = self.customListView.selectedItem()
    if selected is None: return
    editor = CCreateCustomTaskDialog(self,name=selected)
    editor.show()

class CCustomComFileView(CCP4Widgets.CComplexLineWidget):

  MODEL_CLASS =CCP4CustomTaskManager.CCustomComFile
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    qualis['vboxLayout'] = True
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.layout().takeAt(0)
    self.iconButton.deleteLater()
    self.widgets['text'] = CCP4Widgets.CTextEdit(self)
    self.widgets['text'].textChanged.connect(self.updateText)
    self.widgets['name'] = CCP4Widgets.CLineEdit(self)
    self.widgets['name'].textChanged.connect(self.updateName)
    self.layout().addWidget(QtWidgets.QLabel('Command file content',self))
    self.layout().addWidget(self.widgets['text'])
    self.layout().addWidget(QtWidgets.QLabel('Command file name',self))
    self.layout().addWidget(self.widgets['name'])
    self.setModel(model)

  @QtCore.Slot()
  def updateText(self):
    self.model.text.set(self.widgets['text'].toPlainText())

  @QtCore.Slot(str)
  def updateName(self,txt):
    self.model.name.set(str(txt))

  def updateModelFromView(self):
    CCP4Widgets.CComplexLineWidget.updateModelFromView(self)
    #print 'CCustomComFileView.updateModelFromView',repr(self.model),self.model
    
    
class CCustomComFileListView(CCP4Widgets.CListView):
  MODEL_CLASS =CCP4CustomTaskManager.CCustomComFileList
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis =  { 'mode' : 'table',
                'tableItems' : ['comFileName','comFileText'],
                'columnHeaders':['Name','Text'],
                'title' : 'Command file(s)'
               }   
    qualis.update(qualifiers)
    CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)


    
class CCustomTaskParamView(CCP4Widgets.CComplexLineWidget):

  showingMergeTo = QtCore.Signal()

  MODEL_CLASS =CCP4CustomTaskManager.CCustomTaskParam
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    qualis['vboxLayout'] = True
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.layout().takeAt(0)
    self.iconButton.deleteLater()
    
    self.widgets['name'] = CCP4Widgets.CStringView(self)
    self.widgets['dataType'] = CCP4Widgets.CStringView(self,qualifiers={ 'modelClass' : CCP4Data.CI2DataType })
    self.widgets['dataType'].widget.addItem('More..')
    self.widgets['dataType'].widget.setCurrentIndex(0)
    self.widgets['dataType'].widget.setEditable(True)
    self.widgets['dataType'].widget.currentIndexChanged[int].connect(self.handleDataTypeChanged)
    self.widgets['label'] =  CCP4Widgets.CStringView(self)
    self.widgets['obligatory'] = CCP4Widgets.CBooleanView(self)
    self.widgets['saveDataToDb'] = CCP4Widgets.CBooleanView(self)
    self.widgets['function'] = CCP4Widgets.CStringView(self,qualifiers= { 'modelClass' : CCP4CustomTaskManager.CCustomTaskFileFunction} )
    self.widgets['function'].widget.setCurrentIndex(0)
    self.mergeToWidget = QtWidgets.QComboBox( self )
    self.mergeToWidget.currentIndexChanged[int].connect(self.updateMergeToModel)
    self.widgets['splitColumns'] = CCP4Widgets.CStringView(self)
    self.widgets['requiredContentType'] = CMiniMtzRequiredContentWidget(parent=self)
    self.widgets['outputFilePath'] = CCP4Widgets.CStringView(self)
    self.widgets['outputFilePath'].setToolTip("'glob' style relative file path ")
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['name'] )
    line.addWidget(QtWidgets.QLabel('interface label',self))
    line.addWidget(self.widgets['label'])
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Data type',self))
    line.addWidget(self.widgets['dataType'] )
    line.addStretch(1)
    line.addWidget(QtWidgets.QLabel('file function',self))
    line.addWidget(self.widgets['function'])
    line.addStretch(1)
    self.layout().addLayout(line)
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['obligatory'])
    line.addWidget(QtWidgets.QLabel('Parameter is obligatory',self))
    line.addStretch(1)
    line.addWidget(self.widgets['saveDataToDb'])
    line.addWidget(QtWidgets.QLabel('Record in database',self))
    line.addStretch(1)
    self.layout().addLayout(line)

    self.fileFrame  = QtWidgets.QFrame(self)
    self.fileFrame.setLayout(QtWidgets.QHBoxLayout())
    self.fileFrame.layout().setContentsMargins(0,0,0,0)
    self.fileFrame.layout().setSpacing(0)
    self.fileFrame.layout().addWidget(QtWidgets.QLabel("Input file name/ output file search path"))
    self.fileFrame.layout().addWidget(self.widgets['outputFilePath'])
    self.fileFrame.layout().setStretchFactor(self.widgets['outputFilePath'],5)
    self.layout().addWidget(self.fileFrame)
    
    self.mergeToFrame = QtWidgets.QFrame(self)
    self.mergeToFrame.setLayout(QtWidgets.QVBoxLayout())
    self.mergeToFrame.layout().setContentsMargins(0,0,0,0)
    self.mergeToFrame.layout().setSpacing(0)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Merge to or split from monster MTZ:'))
    line.addWidget(self.mergeToWidget)
    line.addStretch(1)
    self.mergeToFrame.layout().addLayout(line)
    self.columnLabelFrame = QtWidgets.QFrame(self)
    self.columnLabelFrame.setLayout(QtWidgets.QHBoxLayout())
    self.columnLabelFrame.layout().setContentsMargins(0,0,0,0)
    self.columnLabelFrame.layout().setSpacing(0)
    self.columnLabelFrame.layout().addWidget(QtWidgets.QLabel('Comma-separated list of column labels'))
    self.columnLabelFrame.layout().addWidget(self.widgets['splitColumns'])    
    self.mergeToFrame.layout().addWidget(self.columnLabelFrame)
    self.requiredContentFrame= QtWidgets.QHBoxLayout()
    self.requiredContentFrame.addWidget(self.widgets['requiredContentType'] )
    self.mergeToFrame.layout().addLayout(self.requiredContentFrame)
    self.layout().addWidget(self.mergeToFrame)
    self.mergeToFrame.hide()
    self.handleFunctionChange()

    self.setModel(model)


  def setModel(self,model=None):
    if model is not None:
      mergeTo = model.mergeTo.__str__()
    else:
      mergeTo = ''
    CCP4Widgets.CComplexLineWidget.setModel(self,model=model)
    self.handleDataTypeChanged(self.widgets['dataType'].widget.currentIndex())
    self.showingMergeTo.emit()
    indx = self.mergeToWidget.findText(mergeTo)
    if indx<0: indx = 0
    self.mergeToWidget.setCurrentIndex(indx)
    if model is not None:
        self.model.dataType.dataChanged.connect(self.setRequiredContentTypeMode)
    self.handleFunctionChange()

  @QtCore.Slot()
  def setRequiredContentTypeMode(self):
    #print 'CCustomTaskParamView.setRequiredContentTypeMode',self.model.dataType.__str__()
    self.widgets['requiredContentType'].setMode(self.model.dataType.__str__())

  def handleFunctionChange(self,indx=None):
    function = self.widgets['function'].widget.currentText().__str__()
    #print 'CCustomTaskParamView.handleFunctionChange',function,type(function)
    if function == 'input':
      self.widgets['requiredContentType'].show()
      self.columnLabelFrame.hide()
    elif function == 'output':
      self.widgets['requiredContentType'].hide()
      self.columnLabelFrame.show()
    else:
      self.widgets['requiredContentType'].hide()
      self.columnLabelFrame.hide()

      
  @QtCore.Slot('QModelIndex')
  def handleDataTypeChanged(self,indx):
    dType = self.widgets['dataType'].widget.itemData(indx).__str__()
    #print 'CCustomTaskParamView.handleDataTypeChanged',dType
    if dType.count('DataFile')>0:
      if not self.mergeToFrame.isVisible():
        self.showingMergeTo.emit()
        self.mergeToFrame.show()
      self.handleFunctionChange()
      self.widgets['requiredContentType'].setMode( dType )
      self.fileFrame.show()
    else:
      if self.mergeToFrame.isVisible():
        self.mergeToFrame.hide()
      self.widgets['requiredContentType'].setMode(None)
      self.fileFrame.hide()
      

  @QtCore.Slot('QModelIndex')
  def updateMergeToModel(self,indx=None):
    #print 'CCustomTaskParamView.updateMergeToModel dataType',self.model.dataType.__str__() , self.mergeToWidget.currentIndex()
    if self.model.dataType.__str__() in ['CObsDataFile', 'CPhsDataFile', 'CMapCoeffsDataFile','CFreeRDataFile']:
      text = self.mergeToWidget.itemText(self.mergeToWidget.currentIndex())
      if len(text)>0:
        self.model.mergeTo.set(text)
      else:
        self.model.mergeTo.unSet()
    else:
      self.model.mergeTo.unSet()


class CMiniMtzRequiredContentWidget(CCP4Widgets.CViewWidget):

  def __init__(self,parent=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    CCP4Widgets.CViewWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setLayout(QtWidgets.QStackedLayout())
    self.obsWidgets = []
    self.phsWidgets = []

    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QHBoxLayout())
    frame.layout().setContentsMargins(0,0,0,0)
    frame.layout().setSpacing(0)
    frame.layout().addWidget(QtWidgets.QLabel(' ',self))
    self.layout().addWidget(frame)
    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QHBoxLayout())
    frame.layout().setContentsMargins(0,0,0,0)
    frame.layout().setSpacing(0)
    frame.layout().addWidget(QtWidgets.QLabel('Observed data allowed content type:',self))
    for item in ['I+/-','F+/-','Imean','Fmean']:
      self.obsWidgets.append(QtWidgets.QCheckBox(item,self))
      frame.layout().addWidget(self.obsWidgets[-1])
    self.layout().addWidget(frame)
    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QHBoxLayout())
    frame.layout().setContentsMargins(0,0,0,0)
    frame.layout().setSpacing(0)
    frame.layout().addWidget(QtWidgets.QLabel('Phase restraints allowed content type:',self))
    for item in ['Hendrickson-Lattmann coefficients','Phi+FOM']:
      self.phsWidgets.append(QtWidgets.QCheckBox(item,self))
      frame.layout().addWidget(self.phsWidgets[-1])
    self.layout().addWidget(frame)
    self.layout().setCurrentIndex(0)
    self.mode = None



  def setMode(self,mode):
    #print 'CMiniMtzRequiredContentWidget.setMode',mode,self.mode
    if mode == self.mode: return
    self.mode = copy.deepcopy(mode)
    if self.mode == 'CObsDataFile':
      self.layout().setCurrentIndex(1)
    elif self.mode == 'CPhsDataFile':
      self.layout().setCurrentIndex(2)
    else:
      self.layout().setCurrentIndex(0)
      

  def getValue(self):
    v = []
    if self.mode == 'CObsDataFile':
      for w in self.obsWidgets: v.append(w.isChecked())
    elif self.mode == 'CPhsDataFile':
      for w in self.phsWidgets: v.append(w.isChecked())
    #print 'CMiniMtzRequiredContentWidget.getValue',v
    return v

  def setValue(self,value=[]):
    if value is None or len(value) == 0:
      #self.setMode(None)
      pass
    elif len(value) == 2:
      #self.setMode('CPhsDataFile')
      for ii in range(2):
        self.phsWidgets[ii].setChecked(bool(value[ii]))
    elif len(value) == 4:
      for ii in range(4):
        self.obsWidgets[ii].setChecked(bool(value[ii]))
      #self.setMode('CObsDataFile')
    

      
class CCustomTaskParamListView(CCP4Widgets.CListView):
  MODEL_CLASS =CCP4CustomTaskManager.CCustomTaskParamList
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    qualis =  { 'mode' : 'table',
                'tableItems' : ['name','label','dataType','function','obligatory','saveDataToDb','mergeTo','outputFilePath','requiredContentType','splitColumns'],
                'columnHeaders':['Name','Label','Data type','File function','Oblig','Record','Merge to MTZ','File path','Representation','Output columns'],
                'title' : 'Parameters in command line or files'
               }   
    qualis.update(qualifiers)
    CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)
    #print 'CCustomTaskParamListView.__init__ editor',self.editor
    self.editor.showingMergeTo.connect(self.populateEditorMergeTo)
    self.listWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Interactive)
    self.listWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)



  @QtCore.Slot()
  def populateEditorMergeTo(self):
    #print 'CCustomTaskParamListView.populateEditorMergeTo'
    self.editor.mergeToWidget.clear()
    self.editor.mergeToWidget.addItem('')
    for obj in self.model:
      if obj != self.editor.model:
        if obj.dataType == 'CMtzDataFile':
          self.editor.mergeToWidget.addItem(obj.name.__str__())


class CCreateCustomTaskDialog(QtWidgets.QDialog):

  customTaskCreated = QtCore.Signal(str)

  def __init__(self,parent=None,name=None):
    QtWidgets.QDialog.__init__(self,parent)

    self.setWindowTitle('Create a custom task')
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setContentsMargins(2,2,2,2)
    self.layout().setSpacing(2)
    
    self.widgets = {}
    self.model = CCP4CustomTaskManager.CCustomTaskDefinition(self,name=name)
    if name is not None:
      self.model.loadDataFromXml(fileName=os.path.join(CUSTOMTASKMANAGER().getDirectory(name),'task.xml'),loadHeader=True,check=False)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Custom task name',self))
    self.widgets['name'] = CCP4Widgets.CStringView(parent=self,model=self.model.name)
    self.widgets['name'].setToolTip('Unique single-word identifier for task')
    line.addWidget(self.widgets['name'])
    line.addWidget(QtWidgets.QLabel('title',self))
    self.widgets['title'] = CCP4Widgets.CStringView(parent=self,model=self.model.title)
    self.widgets['title'].setToolTip('Task name that appears in the user interface')
    line.addWidget(self.widgets['title'])
    line.setStretchFactor(self.widgets['title'],5)
    line.addStretch(1)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Parameter names in command line/file begin with tag:',self))
    self.widgets['tagCharacter'] = CCP4Widgets.CStringView(parent=self,model=self.model.tagCharacter)
    self.widgets['tagCharacter'].setToolTip('Begin parameter names with this tag')
    self.widgets['tagCharacter'].setMaximumWidth(50)
    line.addWidget(self.widgets['tagCharacter'])
    line.addStretch(1)
    self.layout().addLayout(line)

    line =  QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Command line',self))
    self.widgets['comLine'] = CCP4Widgets.CStringView(parent=self,model=self.model.comLine)
    self.widgets['comLine'].setToolTip("Command line with variables such as file names entered as '#HKLIN'")
    line.addWidget(self.widgets['comLine'])
    self.layout().addLayout(line)
    
    self.widgets['comFileList'] = CCustomComFileListView(parent=self,model=self.model.comFileList)
    self.widgets['comFileList'].setToolTip("Command file with variables such as file names entered as '#HKLIN'")
    self.layout().addWidget(self.widgets['comFileList'])

    line = QtWidgets.QHBoxLayout()
    line.addStretch(1)
    update = QtWidgets.QPushButton('Update parameters',self)
    update.clicked.connect(self.updateParametersList)
    line.addWidget(update)
    line.addStretch(1)
    self.layout().addLayout(line)

    self.widgets['paramsList'] = CCustomTaskParamListView(parent=self,model=self.model.paramList)
    self.widgets['paramsList'].setMinimumHeight(400)
    self.layout().addWidget(self.widgets['paramsList'])

    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton('Help',QtWidgets.QDialogButtonBox.HelpRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.help)
    but = buttonBox.addButton('Cancel',QtWidgets.QDialogButtonBox.RejectRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.cancel)
    but = buttonBox.addButton('Save custom task',QtWidgets.QDialogButtonBox.AcceptRole)
    but.setFocusPolicy(QtCore.Qt.NoFocus)
    but.clicked.connect(self.accept)
    self.layout().addWidget(buttonBox)

    self.updateViewFromModel()
    
    #print 'CCreateCustomTaskDialog.init',len(self.model.comFileList)

  def updateViewFromModel(self):
    for key in list(self.widgets.keys()):
      self.widgets[key].updateViewFromModel()
    print('CCreateCustomTaskDialog.updateViewFromModel pluginTitle',self.model.header.pluginTitle)
    self.widgets['name'].setValue(self.model.header.pluginName)
    self.widgets['title'].setValue(self.model.header.pluginTitle)

  def updateModelFromView(self):
    for key in list(self.widgets.keys()):
      self.widgets[key].updateModelFromView()
    print('CCreateCustomTaskDialog.updateModelFromView title',self.widgets['title'].getValue())
    self.model.header.pluginName.set(self.widgets['name'].getValue())
    self.model.header.pluginTitle.set(self.widgets['title'].getValue())

  @QtCore.Slot()
  def help(self):
    WEBBROWSER().loadWebPage(helpFileName='customisation')

  @QtCore.Slot()
  def cancel(self):
    self.hide()

  @QtCore.Slot()
  def accept(self):
    if not self.model.name.isSet():
       QtWidgets.QMessageBox.warning(self,'Create Custom Task','Custom task name must be entered')
       return
    if self.model.name.validity(self.model.name.__str__()).maxSeverity() > Severity.WARNING:
      QtWidgets.QMessageBox.warning(self,'Create Custom Task','Custom task name must be single word')
      return
    if not self.model.comLine.isSet():
       QtWidgets.QMessageBox.warning(self,'Create Custom Task','Custom command line must be entered')
       return
    tag = self.model.tagCharacter.__str__().strip()
    if len(tag)<0 or len(tag)>2:
      QtWidgets.QMessageBox.warning(self,'Create Custom Task','Tag character must be one or two characters')
      return

    paramNameList = self.paramNameList(tag)
    for param in paramNameList:
      obj,indx = self.model.paramList.getItemByName(param)
      if obj is None:
        QtWidgets.QMessageBox.warning(self,'Create Custom Task',"There are undefined parameters in the command line/file.\nPlease click 'Update parameters' and check the definitions in the parameter list")
        return

    mergedMtzs = CUSTOMTASKMANAGER().getMergedMtzs(self.model.paramList)
    #print 'CCreateCustomTaskDialog.accept mergedMtzs',mergedMtzs
    for key,value in list(mergedMtzs.items()):
      if len(value)<=1:
         QtWidgets.QMessageBox.warning(self,'Create Custom Task',"The 'monster' MTZ: "+key+" is created from only one 'mini' MTZ: "+value[0]+'. ' + \
                                   "A 'monster' MTZ is expected to combine two or more 'mini' MTZs.")
         return

    defFile = CUSTOMTASKMANAGER().getCustomFile(name=self.model.name.__str__(),mustExist=True)
    if defFile is not None:
       msgBox = QtWidgets.QMessageBox()
       msgBox.setWindowTitle('Create custom task')
       msgBox.setText('There is already a custom task directory called '+self.model.name.__str__())
       b = msgBox.addButton(QtWidgets.QMessageBox.Cancel)
       b = msgBox.addButton('Overwrite',QtWidgets.QMessageBox.ApplyRole)
       b.clicked.connect(functools.partial(self.createCustomTask,True))
       msgBox.exec_()
    else:        
      self.createCustomTask()
      self.hide()

  @QtCore.Slot(bool)
  def createCustomTask(self,overwrite=False):
    self.hide()
    #try:
    err = CUSTOMTASKMANAGER().createCustomTask(name=self.model.name.__str__(),title=self.model.title.__str__(),container=self.model,overwrite=overwrite)

    if err.maxSeverity()>Severity.WARNING:
      warningMessage(err, 'Create custom task','Error saving custom task',parent=self)
      return
    
    self.customTaskCreated.emit(self.model.name.__str__())

  @QtCore.Slot()
  def updateParametersList(self):

    try:
      paramListCurrentRow = self.widgets['paramsList'].editItemIndex
      paramListCurrentRowName = self.model.paramList[paramListCurrentRow].name.__str__()
    except:
      paramListCurrentRowName = None
      
    self.updateModelFromView()
    tag = self.model.tagCharacter.__str__().strip()
    if len(tag)<0 or len(tag)>2:
      QtGui.MessageBox.warning(self,'Create Custom Task','Tag character must be one or two characters')
      return
    else:
      # Reset just in case the strip was necessary
      self.model.tagCharacter.set(tag)

    paramNameList = self.paramNameList(tag)

    for param in paramNameList:
      obj,indx = self.model.paramList.getItemByName(param)
      if obj is None:
        self.model.paramList.append( { 'name' : param } )

    self.widgets['paramsList'].updateViewFromModel()
    if paramListCurrentRowName is not None:
      obj,indx = self.model.paramList.getItemByName(paramListCurrentRowName)
    else:
      indx = None
    if indx is None: indx = 0
    self.widgets['paramsList'].handleRowChange(row=indx,force=True)
    

  def paramNameList(self,tag='#'):
    paramNameList = CUSTOMTASKMANAGER().extractParams(tag,self.model.comLine.__str__())
    for ii in range(len(self.model.comFileList)):
      paramNameList.extend(CUSTOMTASKMANAGER().extractParams(tag,self.model.comFileList[ii].text.__str__()))
    return paramNameList
