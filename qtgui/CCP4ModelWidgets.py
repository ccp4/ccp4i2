from __future__ import print_function

"""
     qtgui/CCP4ModelWidgets.py: CCP4 Gui Project
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
     GNU Lesser General Public License for more details.
"""

##@package CCP4ModelWidgets (QtGui) Collection of widgets for model data types

import os
import sys
import re
import functools
from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4ModelData
from core import CCP4Modules,CCP4Utils
from core.CCP4ErrorHandling import *
from qtgui import CCP4Widgets
from qtgui import CCP4SequenceList
from qtgui import CCP4RefmacMultiAtomSelection

class CResidueRangeView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4ModelData.CResidueRange
    ERROR_CODES = { }

    def __init__(self,parent=None,model=None,qualifiers={}):
      qualis = qualifiers
      CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
      if model is not None: self.setModel(model)

      self.layout().addWidget(QtWidgets.QLabel('Chain'))
      if self.editable:
        self.widgets['chainId'] = CCP4Widgets.CComboBox(self,charWidth=2,dragType=self.dragType())
      else:
        self.widgets['chainId'] = CCP4Widgets.CLabel(self,charWidth=4,uneditable=True,dragType=self.dragType())
      self.layout().addWidget(self.widgets['chainId'])

      for item,label in [['firstRes','residue'],['lastRes','to']]:
        self.layout().addWidget(QtWidgets.QLabel(label))
        if self.editable:
          self.widgets[item] = CCP4Widgets.CLineEdit(self,charWidth=10,dragType=self.dragType())
        else:
          self.widgets[item] = CCP4Widgets.CLabel(self,charWidth=10,uneditable=True,dragType=self.dragType())
        self.layout().addWidget(self.widgets[item])

      if self.editable:
        for item in list(self.widgets.keys()):
          self.widgets[item].editingFinished.connect(self.updateModelFromView)
          self.widgets[item].acceptDropDataSignal.connect(self.acceptDropData)
        pdbFileObj = self.model.parent().getDataByKey('pdbFileKey')
        if pdbFileObj is not None:
          self.loadChains()
          pdbFileObj.dataChanged.connect(self.loadChains)
      self.layout().addStretch(5)

    @QtCore.Slot()
    def loadChains(self):
      self.widgets['chainId'].clear()
      for chn in self.model.parent().getChains():
        self.widgets['chainId'].addItem(chn)
          

    def getMenuDef(self):
      return ['clear','copy','paste','help']

class CSequenceEdit:
  def __init__(self,parent=None,model=None):
    pass
      

class CSequenceView(CCP4Widgets.CComplexLineWidget):

  ERROR_CODES = { }
  MODEL_CLASS = CCP4ModelData.CSequence
  
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'gridLayout':True}
    qualis.update(qualifiers)
    
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
    if model is not None: self.setModel(model)
    #if  qualis.get('stackButton',None) is not None:
    #  self.layout().addWidget(qualis['stackButton'],0,1)
    self.layout().addWidget(QtWidgets.QLabel('Description:',self),1,0)
    if self.editable:
        self.widgets['identifier'] = CCP4Widgets.CLineEdit(self,qualifiers=qualis)
        self.widgets['reference'] = CCP4Widgets.CLineEdit(self,qualifiers=qualis)
        self.widgets['referenceDb'] = CCP4Widgets.CComboBox(self,qualifiers= CCP4ModelData.CSequence.CONTENTS['referenceDb']['qualifiers'] ) 
    else:
        self.widgets['identifier'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
        self.widgets['reference'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
        self.widgets['referenceDb'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
    self.widgets['identifier'].setToolTip('Name of the sequence')
    self.widgets['reference'].setToolTip('Database identifier')
    self.widgets['referenceDb'].setToolTip('Sequence database')    
    self.layout().addWidget(self.widgets['identifier'],1,1,1,2)
    self.layout().addWidget(QtWidgets.QLabel('Database reference:',self),0,0)
    self.layout().addWidget(self.widgets['referenceDb'],0,1)
    self.layout().addWidget(self.widgets['reference'],0,2)
    #if self.editable:
    #  self.widgets['moleculeType']=CCP4Widgets.CComboBox(self,qualifiers= CCP4ModelData.CSequence.CONTENTS['moleculeType']['qualifiers'])
    #else:
    #  self.widgets['moleculeType']=CCP4Widgets.CLabel(self,qualifiers=qualis)
    #self.layout().addWidget(self.widgets['moleculeType'],0,5)

    self.widgets['sequence'] = CCP4Widgets.CTextEdit(self,qualifiers=qualis)
    self.layout().addWidget(self.widgets['sequence'],2,0,1,3)

    for item in CCP4ModelData.CSequence.CONTENTS_ORDER:
      widget = self.widgets.get(item,None)
      if widget is not None:
        #print 'CSequenceView.__init__',item,repr(self.model.getDataObjects(item)),type(self.model.getDataObjects(item))
        if model is not None:
          tT = self.model.getDataObjects(item).qualifiers('toolTip')
          if tT is NotImplemented: tT = ''
          widget.setToolTip(tT)
        if self.editable:
          widget.editSignal.connect(self.updateModelFromView)
          widget.acceptDropDataSignal.connect(self.acceptDropData)


  def getMenuDef(self):
    menu = CCP4Widgets.CComplexLineWidget.getMenuDef(self)
    menu.insert(0,'content')
    return menu

  def updateViewFromModel(self):
    #print 'CSequenceView.updateViewFromModel model',self.model
    #print 'CSequenceView.updateViewFromModel',self.model.get()
    if len(self.widgets) == 0: return
    if self.model is None:
      for item in ['identifier','reference','referenceDb','sequence']:
        self.widgets[item].setValue('')
    else:
      for item in ['identifier','reference','referenceDb','sequence']:
        #print 'CSequenceView.updateViewFromModel',item,self.model.__dict__['_value'][item].__str__()
        self.widgets[item].setValue(self.model.__dict__['_value'][item].__str__())
    self.validate()

class CSeqDataFileView(CCP4Widgets.CDataFileView):
  MODEL_CLASS = CCP4ModelData.CSeqDataFile
  
  loadMultiple = QtCore.Signal(list)

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'vboxLayout':True}
    qualis.update(qualifiers)
    self.enableEdit = qualifiers.get('enableEdit',True)
    CCP4Widgets.CDataFileView.__init__(self,parent=parent,model=model,qualifiers=qualis)
    
    qualis['iconButton'] = False  
    if model is not None:
      self.widgets['sequence'] = CSequenceView(parent=parent,model=model.getFileContent(),qualifiers=qualis)
    else:
      self.widgets['sequence'] = CSequenceView(parent=parent,qualifiers=qualis)   
    self.layout().addWidget(self.widgets['sequence'])
    self.widgets['sequence'].hide()
    self.setModel(model)
    parentTask = self.parentTaskWidget()
    if parentTask is not None:
      parentTask.doFix.connect(self.saveSequence)
                             

  def getMenuDef(self):
    if self.editable:
      menu = ['clear',['View','view_text'],'sep','copy','paste','help']
      if self.enableEdit: menu.insert(2,'edit')
    else:
      menu = [['View','view_text'],'sep','copy','editLabel','export','help']
    return menu

  def getActionDef(self,name):
    if name == 'edit':
      return dict (
        text = self.tr("View/edit sequence"),
        tip = self.tr('You can cut-n-paste sequence into the widget'),
        slot = self.showEditor
      )
    else:
      return CCP4Widgets.CDataFileView.getActionDef(self,name)


  def showEditor(self,visible=None):
    if visible is None: visible = not self.widgets['sequence'].isVisible()
    #print 'CSeqDataFile.showEditor',visible,self.model.exists(),self.widgets['sequence'].model
    if visible and self.model.exists():
        self.model.loadFile()
        #print 'CSeqDataFile.showEditor after loadFile',repr(self.widgets['sequence'].model),self.widgets['sequence'].model
        if self.widgets['sequence'].model is None:
          self.widgets['sequence'].setModel(model.fileContent)
        self.widgets['sequence'].updateViewFromModel()
    self.widgets['sequence'].setVisible(visible)

  def setModel(self,model):
    CCP4Widgets.CDataFileView.setModel(self,model)
    if self.widgets.get('sequence',None) is None: return
    if model is not None:
      if hasattr(model,'fileContent'):
        self.widgets['sequence'].setModel(model.fileContent)
      model.dataChanged.connect(self.loadSequence)
      self.loadSequence()
    else:
      self.widgets['sequence'].setModel(None)

  @QtCore.Slot()
  def loadSequence(self):
    if not self.editable: return
    if self.model.isSet():
      try:
        self.model.loadFile()
      except:
        self.model.fileContent.unSet()
    else:
      self.model.fileContent.unSet()
    self.updateViewFromModel()
    #self.widgets['sequence'].updateViewFromModel()
    
  @QtCore.Slot()
  def saveSequence(self):
    widgetValue = self.widgets['sequence'].widgets['sequence'].getValue()
    #print 'CSeqDataFileView.saveSequence',widgetValue
    if widgetValue is None or len(widgetValue.strip())==0: return
    jobId = self.parentTaskWidget().jobId()
    self.widgets['sequence'].updateModelFromView()
    if not self.model.baseName.isSet():
      self.model.saveSequence(jobId=jobId)
      self.validate()
    else:
      if self.model.fileContent.sequence.isSet():
        # User has a filename set and a sequence set - test if they match
        fileContent = CCP4ModelData.CSequence()
        #try:
        if 1:
          fileContent.loadFile(str(self.model),format=self.model.__dict__['format'])
          err = self.model.fileContent.assertSame(fileContent)
          #print 'CSeqDataFileView.saveSequence',self.model.fileContent,fileContent,err.report()
          if err.maxSeverity()>SEVERITY_WARNING:
            mess = QtWidgets.QMessageBox(self)
            mess.setWindowTitle('Sequence file')
            mess.setText('Save sequence or reference data to a new file')
            mess.addButton('Save to file', QtWidgets.QMessageBox.AcceptRole)
            mess.addButton('Keep existing file',QtWidgets.QMessageBox.ActionRole)
            mess.addButton('Cancel', QtWidgets.QMessageBox.RejectRole)
            mess.show()
            ret = mess.exec_()
            if ret == QtWidgets.QMessageBox.RejectRole:
              return
            elif ret == QtWidgets.QMessageBox.ActionRole:
              self.model.loadFile()
            else:

              self.model.saveSequence(jobId=jobId)
            self.validate()
              
        
  def validate(self,isValid=None,reportMessage=True):
    dataFileIsValid = CCP4Widgets.CDataFileView.validate(self,isValid=isValid,reportMessage=reportMessage)
    if self.editable:
      try:
        sequenceIsValid =  self.widgets['sequence'].validate(isValid=isValid,reportMessage=reportMessage)
        if self.model.qualifiers('allowUndefined') and sequenceIsValid.maxSeverity()<SEVERITY_WARNING:
          #Allow the sequence to be undefined if the sequence file is allowed to be undefined
          isValid = dataFileIsValid
        else:
          isValid = dataFileIsValid and sequenceIsValid
      except:
        isValid = dataFileIsValid
    else:
      isValid = dataFileIsValid
    if isValid != self.isValid:
      self.isValid = isValid
      self.setProperty("isValid",isValid)
      self.updateValidityIndicator()
    return self.isValid

  @QtCore.Slot(str,str,dict)
  def handleBrowserOpenFile(self,filename,downloadInfo,**kw):
    kw['validate'] = False
    kw['updateView'] = False
    # Imperative to clear all model related data before trying to load new file
    self.model.unSet()
    self.model.blockSignals(True)
    CCP4Widgets.CDataFileView.handleBrowserOpenFile(self,filename,downloadInfo,**kw)
    self.model.blockSignals(False)
    #print 'CSeqDataFileView.handleBrowserOpenFile identifiers',self.model,self.model.isSet(),self.model.__dict__.get('identifiers','No identifiers'),self.model.fileContent.__dict__.get('loadWarning','No loadWarning'),self.model.__dict__['format']
    if not self.model.isSet():
      self.validate()
      return
    if self.model.__dict__['format'] == 'unknown':
      self.showEditor(True)
      QtWidgets.QMessageBox.warning(self,'Reading sequence file','The format of the file has not been recognised. Please edit the sequence to remove any non-sequence information.' )
    elif self.model.fileContent.__dict__.get('loadWarning',None) is not None:
       #print 'CSeqDataFileView.handleBrowserOpenFile loadWarning',self.model.fileContent.__dict__.get('loadWarning')
       if len(self.model.fileContent.__dict__['loadWarning'])>0 and self.model.fileContent.__dict__['loadWarning'][0]['code'] == 108:
         win = self.makeImportReport('Importing '+str(filename),self.acceptFixedPir)
         label = QtWidgets.QLabel(self)
         label.setText( 'The input file is not correct PIR format but has been fixed and is shown below.\nIf this is still incorrect please correct the file and select it again.\n'  + CCP4ModelData.PIR_DESCRIPTION + '\nYour corrected file:')
         #label.setReadOnly(True)
         win.layout().insertWidget(0,label)
         label = QtWidgets.QTextEdit(self)
         label.setText(CCP4Utils.readFile(self.model.fileContent.__dict__['loadWarning'][0]['details']))
         label.setReadOnly(True)
         win.layout().insertWidget(1,label)
         win.show()
         return
       elif self.model.fileContent.__dict__['loadWarning'].maxSeverity()>SEVERITY_WARNING:
         self.model.fileContent.__dict__['loadWarning'].warningMessage(windowTitle='Importing sequence file',parent=self,message='Failed importing sequence file')
         self.model.unSet()
         self.updateViewFromModel()
         return
    elif len(self.model.__dict__.get('identifiers',[]))>1:
      #print 'handleBrowserOpenFile drawing idChooser',getattr(self,'idChooser',None)
      self.idChooser = QtWidgets.QListWidget(self)
      if isinstance(self.model.parent(), CCP4ModelData.CAsuContentSeqList):
          self.idChooser.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
          title = 'Select one or more sequences from the file'
      else:
          self.idChooser.setSelectionMode(QtWidgets.QListWidget.SingleSelection)
          title = 'Select one of the sequences in the file'
      win = self.makeImportReport('Importing '+str(filename),self.handleIdChooser)
      win.layout().insertWidget(0,self.idChooser)
      win.layout().insertWidget(0,QtWidgets.QLabel(title, self))
      row = 0
      for item in self.model.__dict__['identifiers']:
          self.idChooser.addItem(item)
          if isinstance(self.model.parent(), CCP4ModelData.CAsuContentSeqList):
              self.idChooser.item(row).setSelected(True)
          row += 1
      win.show()
      win.raise_()
      return

    if self.model.__dict__['format'] is None or self.model.__dict__['format'] != 'internal':
      #self.model.importFile(jobId=self.parentTaskWidget().jobId(),jobNumber=self.parentTaskWidget().jobNumber())
      self.doImportFile()
    else:
      self.model.dataChanged.emit()
      
    # No problems just display it!
    self.updateViewFromModel()
    self.validate()

  def makeImportReport(self,title,callBack=None):
    win = QtWidgets.QDialog(self)
    win.setWindowTitle(title)
    win.setLayout(QtWidgets.QVBoxLayout())
    line = QtWidgets.QHBoxLayout()
    box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
    win.layout().addLayout(line)
    line.addStretch(1)
    line.addWidget(box)
    line.addStretch(1)
    if callBack is not None: box.clicked.connect(callBack)
    return win

  def doImportFile(self,label=None):
      if self.model.fileContent.identifier.isSet() and len(self.model.fileContent.identifier)>0:
        anno = self.model.fileContent.identifier
      elif label is not None:
        anno = label+' imported from '+ str(self.model.baseName)
      else:
        anno = self.model.qualifiers('guiLabel')+' imported from '+ str(self.model.baseName)

      #print 'doImportFile anno',anno,'label',label

      #self.model.blockSignals(True)
      validatedFile=self.model.fileContent.__dict__.get('validatedFile',None)
      self.model.importFile(jobId=self.parentTaskWidget().jobId(),annotation=anno,validatedFile=validatedFile,jobNumber=self.parentTaskWidget().jobNumber())
      self.model.annotation = anno
      #self.model.blockSignals(False)
      self.updateViewFromModel()
      if CCP4Modules.PREFERENCES().AUTO_INFO_ON_FILE_IMPORT and not self.model.dbFileId.isSet():        
          self.openInfo(label=self.model.qualifiers('guiLabel').lower(),sourceFileAnnotation=self.model.__dict__.get('sourceFileAnnotation',''))

  def acceptFixedPir(self,but):
      #print 'CSeqDataFileView.acceptFixedPir',but.text().__str__()
      if but.text().__str__() == 'Cancel':
          self.model.unSet()
          self.updateViewFromModel()
      else:
          self.doImportFile()
      win = but.window()
      win.close()
      self.updateViewFromModel()

  @QtCore.Slot('QPushButton')
  def handleIdChooser(self,but):
    #print 'CSeqDataFileView.handleIdChooser', but.text().__str__()
    if self.idChooser.selectionMode() == QtWidgets.QListWidget.SingleSelection:
        recordList = [self.idChooser.currentRow()]
    else:
        recordList = []
        for item in self.idChooser.selectedItems():
            recordList.append(self.idChooser.row(item))

    win = but.window()
    win.setModal(False)
    win.close()
    
#FIXME - Surely this should be if but == QtWidgets.QDialogButtonBox.Cancel ?
    if but.text().__str__() == 'Cancel':
      self.model.unSet()
      self.updateViewFromModel()
      return
    else:
      #print 'CSeqDataFileView.handleIdChooser',record; import sys; sys.stdout.flush()
      if len(recordList) == 1:
          self.model.fileContent.loadExternalFile(self.model.__str__(),self.model.__dict__['format'],record=recordList[0])
          self.doImportFile()
          print('handleIdChooser', self.model.fileContent.identifier,self.model.fileContent.sequence)
          self.updateViewFromModel()
      else:
          self.loadMultiple.emit(recordList)

  def filterText(self):
    # make the filters text for QFileDialog to include the alignment extensions
    textList = []
    for desc,extList in [ [CCP4ModelData.CSeqDataFile.QUALIFIERS['mimeTypeDescription'] , CCP4ModelData.EXTLIST ],
                          ]:
      text = desc + ' ('
      for ext in list(extList.keys()):
        text = text + '*.'+ext+' '
      textList.append( text[0:-1]+')' )    
    return textList

'''
This reimplemetation of CDataFileView provides tools for user to merge Refmac dictionary files
It is no longer presented on the gui but the issues are handled in the appropriate scripts
eg.CPluginScript.mergeDictToProjectLib() is used.

class CDictDataFileView(CCP4Widgets.CDataFileView):
  MODEL_CLASS = CCP4ModelData.CDictDataFile
  PROJECT_FILE_LABEL = 'Ideal ligand geometry file for project'

  def __init__(self,parent=None,model=None,qualifiers={},**kw):
    qualis = {}
    qualis.update(qualifiers)
    qualis.update(kw)
    CCP4Widgets.CDataFileView.__init__(self,parent=parent,model=model,qualifiers=qualis)
    

  def getActionDef(self,name):
    if name == 'manage_dict':
      return dict (
        text = self.tr("View/edit dictionary"),
        tip = self.tr('View/edit the currently selected dictionary'),
        slot = self.manageProjectDictionary
      )
    else:
      return CCP4Widgets.CDataFileView.getActionDef(self,name)

  def manageProjectDictionary(self):
    if not hasattr(self,'dictManager'):     
      #self.dictManager = CDictDataDialog(model=self.model.fileContent)
      self.dictManager = CDictDataDialog(parent=self.parentTaskWidget(),projectId=self.parentTaskWidget().projectId())
    self.dictManager.show()
    self.dictManager.raise_()
    
  @QtCore.Slot()
  def handleMerge(self):
    from qtgui import CCP4FileBrowser
    self.mergeFileBrowser = CCP4FileBrowser.CFileDialog(parent=self,title='Select dictionary file to merge into project dictionary',
          defaultSuffix=CCP4Modules.MIMETYPESHANDLER().getMimeTypeInfo(name='application/refmac-dictionary',info='fileExtensions'),
                                 fileMode=QtWidgets.QFileDialog.ExistingFiles,saveButtonText='Merge these files')
    self.mergeFileBrowser.show()
    self.mergeFileBrowser.selectFiles.connect(self.mergeFiles)

  def mergeFiles(self,selectedFiles=[]):   
    # 'CDictDataFileView.mergeFiles',selectedFiles
    if hasattr(self,'mergeFileBrowser'):
      self.mergeFileBrowser.close()
      self.mergeFileBrowser.deleteLater()
    
    workDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(
        jobId=self.parentTaskWidget().jobId(),projectId=self.parentTaskWidget().projectId())
    newFile,err = self.model.mergeInDictFiles(dictFileList=selectedFiles,parentWorkDirectory=workDirectory)
    
    for fileName in selectedFiles:
      err = self.model.fileContent.mergeFile(fileName=fileName,overwrite=True)
      #print 'CDictDataFileView.mergeFiles',err.report(),err.maxSeverity()
      if err.maxSeverity()>SEVERITY_WARNING:
        err.warningMessage(windowTitle='Error merging Refmac dictionaries',parent=self,message='Error attempting to merge Refmac dictionaries')
      return

  """
  def loadJobCombo(self):
    CCP4Widgets.CDataFileView.loadJobCombo(self)
    if self.jobCombo.findText(CDictDataFileView.PROJECT_FILE_LABEL)<0:
      # Call the defaultProjectDict even though result is not used it will ensure that
      # a file exists (by creting emptly file if necessary)
      if  self.model is None: return
      dictFile = self.model.defaultProjectDict(projectId=self.parentTaskWidget().projectId())
      self.jobCombo.insertItem(1,CDictDataFileView.PROJECT_FILE_LABEL)
      self.jobCombo.setCurrentIndex(1)
  """

  def handleFollowFrom(self,contextJobId,projectId):
    self.jobCombo.setCurrentIndex(1)

  def handleJobComboChange(self,indx=None):
    self.connectUpdateViewFromModel(False)
    if indx is None: indx = 0
    if str(self.jobCombo.itemText(indx)) == CDictDataFileView.PROJECT_FILE_LABEL:
      self.model.set(self.model.defaultProjectDict(projectId=self.parentTaskWidget().projectId()))
      self.model.annotation = CDictDataFileView.PROJECT_FILE_LABEL
    else:
      CCP4Widgets.CDataFileView.handleJobComboChange(self,indx0=indx)
    
'''
      
class CMonomerView(CCP4Widgets.CComplexLineWidget):

  MODEL_CLASS = CCP4ModelData.CMonomer

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'gridLayout':True}
    qualis.update(qualifiers)
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
    if model is not None: self.setModel(model)
    self.layout().addWidget(QtWidgets.QLabel('Identifier:',self),0,1)
    if self.editable:
      self.widgets['identifier'] = CCP4Widgets.CStringView(self,qualifiers=qualis)
    else:
      self.widgets['identifier'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
    self.layout().addWidget(self.widgets['identifier'],0,2)
    self.layout().addWidget(QtWidgets.QLabel('Formula:',self),0,3)
    if self.editable:
      self.widgets['formula'] = CCP4Widgets.CStringView(self,qualifiers=qualis)
    else:
      self.widgets['formula'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
    self.layout().addWidget(self.widgets['formula'],0,4)

    self.layout().addWidget(QtWidgets.QLabel('Dictionary name:',self),1,1)
    if self.editable:
      self.widgets['dictionaryName'] = CCP4Widgets.CStringView(self,qualifiers=qualis)
    else:
      self.widgets['dictionaryName'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
    self.layout().addWidget(self.widgets['dictionaryName'],1,2)
    
    self.layout().addWidget(QtWidgets.QLabel('Smiles string:',self),1,3)
    if self.editable:
      self.widgets['smiles'] = CCP4Widgets.CStringView(self,qualifiers=qualis)
    else:
      self.widgets['smiles'] = CCP4Widgets.CLabel(self,qualifiers=qualis)
    self.layout().addWidget(self.widgets['smiles'],1,4)


    '''
    toolTip = self.model.qualifiers('toolTip')
    if toolTip is NotImplemented:
      toolTip = ''
    else:
      toolTip = toolTip
    
      
    for item in self.model.CONTENTS_ORDER:
      widget = self.widgets[item]
      tT = self.model.getDataObjects(item).qualifiers('toolTip')
      if tT is NotImplemented: tT = ''
      widget.setToolTip(toolTip+tT)
      if self.editable:
        widget.editSignal.connect(self.updateModelFromView)
        widget.acceptDropDataSignal.connect(self.acceptDropData)
    '''

class CPdbEnsembleItemView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4ModelData.CPdbEnsembleItem

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'vboxLayout':True,'iconButton':False}
    qualis.update(qualifiers)
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
    if model is not None: self.setModel(model)
    self.widgets = {}
    self.widgets['structure'] = CPdbDataFileView(self,qualifiers={'iconButton':True,'iconName':'PdbDataFile','ifAtomSelection' : True})
    self.layout().addWidget(self.widgets['structure'])
    self.widgets['identity_to_target'] = CCP4Widgets.CFloatView(self)
    self.widgets['rms_to_target'] = CCP4Widgets.CFloatView(self)
    self.widgets['number'] = CCP4Widgets.CIntView(self, qualifiers = { 'editable' : self.editable, 'guiMode' : 'combo', 'enumerators' : [i for i in range(21)], 'menuText' : [str(i) for i in range(21)]})
    self.widgets['number'].setMaximumWidth(60)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Copies:',self))
    line.addWidget(self.widgets['number'])
    line.addWidget(QtWidgets.QLabel('Sequence identity (fractional):',self))
    line.addWidget(self.widgets['identity_to_target'])
    line.addWidget(QtWidgets.QLabel('OR RMS difference:',self))
    line.addWidget(self.widgets['rms_to_target'])
    self.layout().addLayout(line)

  def setModel(self, model):
    ensemble = None
    if model is not None and hasattr(model,'parent') and model.parent() is not None:
      ensemble = model.parent().parent()
    cEnsembleLabelView = self.parent().findChildren(CEnsembleLabelView)[0]
    if self.model is not None:
      for item in ['structure','identity_to_target','rms_to_target']:
        try:
            self.model.get(item).dataChanged.disconnect(self.validate)
        except:
            print("disconnect fail")
      if ensemble is not None:
        for item in ['number']:
          try:
              ensemble.get(item).dataChanged.disconnect(self.numberChanged)
          except:
              print("disconnect fail")
  
    if model is None or isinstance(model,self.MODEL_CLASS):
      from qtgui.CCP4Widgets import CViewWidget
      CViewWidget.setModel(self,model)
      
      if model is not None:
        model.dataChanged.connect(self.validate)
        toolTip = model.qualifiers('toolTip')
        if toolTip is not NotImplemented and toolTip is not None  and self.iconButton is not None:
          self.iconButton.setToolTip(toolTip+'\n'+self.iconButton.toolTip())
        for key,w in list(self.widgets.items()):
          if key in ['structure','identity_to_target','rms_to_target']:
            if isinstance(w,CViewWidget): w.setModel(model.get(key))
          elif key in ['number']:
            if ensemble is not None:
              if isinstance(w,CViewWidget): w.setModel(ensemble.get(key))
      else:
        for key,w in list(self.widgets.items()):
           if isinstance(w,CViewWidget): w.setModel(None)
  
    # Set allowUndefined True for the first CPdbEnsembleItem in the first CEnsemble where the CEnsembleList
    # is allowed zero length
    if model is not None and ensemble is not None:
      if isinstance(ensemble,CCP4ModelData.CEnsemble) and ensemble.qualifiers('allowUndefined'):
        if model.parent().index(model) == 0:
          model.setQualifier('allowUndefined',True)
          model.structure.setQualifier('allowUndefined',True)

    if model is not None and self.editable:
      for item in ['structure','identity_to_target','rms_to_target']:
        model.get(item).dataChanged.connect(self.validate)
      for item in ['number']:
        if ensemble is not None:
          ensemble.get(item).dataChanged.connect(self.numberChanged)

    if self.widgets['structure'].widgets.get('selection',None) is None:
      self.widgets['structure'].showAtomSelection()

  def updateViewFromModel(self):
    for key,widget in list(self.widgets.items()):
      widget.updateViewFromModel()
    #self.parent().findChildren(CEnsembleLabelView)[0].updateViewFromModel()
      
  @QtCore.Slot(dict)
  def numberChanged(self, **kw):
    cEnsembleLabelView = self.parent().findChildren(CEnsembleLabelView)[0]
    self.parent().parent().updateViewFromModel(**kw)
    cEnsembleLabelView.validate()
  
  '''
  def getMenuDef(self):
    return self.widgets['structure'].getMenuDef()

  def getActionDef(self,name):
    return self.widgets['structure'].getActionDef(name)

  def showAtomSelection(self):
    self.widgets['structure'].showAtomSelection()

  def viewContents(self):
    self.widgets['structure'].viewContents()
  '''

class CEnsembleView(CCP4Widgets.CComplexLineWidget):
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'vboxLayout':True}
    qualis.update(qualifiers)
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
    if model is not None: self.setModel(model)
    self.widgets = {}
    line = QtWidgets.QHBoxLayout()
    iconWidgetItem = self.layout().takeAt(0)
    line.addWidget(iconWidgetItem.widget())
    line.addWidget(QtWidgets.QLabel('Label for ensemble:',self))
    self.widgets['label'] = CCP4Widgets.CStringView(self)
    line.addWidget(self.widgets['label'])
    self.layout().addLayout(line)
    self.widgets['pdbItemList'] = CCP4Widgets.CListView(self,model=self.model.pdbItemList,qualifiers={
        'editable' : self.editable,
        'mode' : 'table',
        'tableItems': ['structure','identity_to_target','rms_to_target'],
        'columnHeaders':['Filename','Identity','RMS']  })
    self.layout().addWidget(self.widgets['pdbItemList'])

class CEnsembleLabelView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4ModelData.CEnsemble
  def __init__(self,parent=None,model=None,qualifiers={}):
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualifiers=qualifiers)
    self.widgets['use'] = CCP4Widgets.CBooleanView(self, qualifiers = { 'editable' : self.editable } )
    self.widgets['number'] = CCP4Widgets.CIntView(self, qualifiers = { 'editable' : self.editable, 'guiMode' : 'combo', 'enumerators' : [i for i in range(21)], 'menuText' : [str(i) for i in range(21)]})
    self.widgets['number'].setMaximumWidth(60)
    self.widgets['label'] = CCP4Widgets.CStringView(self, qualifiers = { 'editable' : self.editable } )
    self.layout().addWidget( self.widgets['use'] )
    self.layout().addWidget(QtWidgets.QLabel('Find',self))
    self.layout().addWidget( self.widgets['number'] )
    self.layout().addWidget(QtWidgets.QLabel('of search ensemble called',self))
    self.layout().addWidget( self.widgets['label'] )
    self.setModel(model)

  def validate(self,isValid=None,reportMessage=True):
    # Just validate the number and label
    if isValid is None:
      if self.model is None:
        isValid = False
      else:
        v = self.model.number.validity(self.model.number.get())
        v.extend(self.model.label.validity(self.model.label.get()))
        isValid = (v.maxSeverity()<=SEVERITY_WARNING)
    
    CCP4Widgets.CComplexLineWidget.validate(self,isValid=isValid,reportMessage=reportMessage)

class CEnsembleListView(CCP4Widgets.CTreeView):
  MODEL_CLASS = CCP4ModelData.CEnsembleList
  def __init__(self,parent=None,model=None,qualifiers={}):
    displayRole = QtCore.Qt.DisplayRole
    qualis = { 'editors' : [ { 'modelClass' : CCP4ModelData.CEnsemble, 'name' : 'ensemble', 'label':'ensemble' } ,
                             { 'modelClass' : CCP4ModelData.CPdbEnsembleItem , 'name' : 'pdbEnsembleItem','label':'structure in ensemble' } ],

               'columnHeaders':[  {displayRole:'Ensemble/Filename','width':240},
                                  {displayRole:'Selection','width':240},
                                  {displayRole:'Identity','width':50},
                                  {displayRole:'RMS','width':50}  ]
               }
    """
        'hierarchy' : [ {'name':'self','label':'ensemble','editorClass':CEnsembleLabelView, 'grey':True,  'list':
                                   [{'name':'pdbItemList','label':'structure in ensemble','editorClass':CPdbEnsembleItemView }]
                             } ],
    """
    #print 'INTO CEnsembleListView',model
    qualis.update(qualifiers)
    super(CEnsembleListView,self).__init__(parent,model=model,qualifiers=qualis)
    self.listWidget.setMinimumHeight(150)

class CPdbDataFileView(CCP4Widgets.CDataFileView):
  MODEL_CLASS = CCP4ModelData.CPdbDataFile
  def __init__(self,parent=None,model=None,qualifiers={}):
    self.ifAtomSelection = qualifiers.get('ifAtomSelection',False)
    qualis = { 'vboxLayout' : True }
    qualis.update(qualifiers)
    #print 'CPdbDataFileView parent',parent
    #CCP4Widgets.CDataFileView.__init__(self,parent=parent,model=model,qualifiers=qualis)
    self.nChainsCol = 0
    super(CPdbDataFileView,self).__init__(parent=parent,model=model,qualifiers=qualis)
    #if self.ifAtomSelection: self.showAtomSelection()
  
  def setUndefinedAllowedBehaviour(self,allowUndefined):
      """
      This allows us to avoid the problem with the selected atoms widget being
      because no data is set even when we have "allowUndefined" true.
      """
      if not allowUndefined:
          if hasattr(self,"jobCombo"):
              self.jobCombo.setItemText(0,"..must be selected")
              if hasattr(self,"widgets"):
                  if "selection" in self.widgets:
                      self.widgets["selection"].allowUndefined = False
      else:
          if hasattr(self,"jobCombo"):
              self.jobCombo.setItemText(0,"..is not used")
              if hasattr(self,"widgets"):
                  if "selection" in self.widgets:
                      self.widgets["selection"].allowUndefined = True

  def getMenuDef(self):
    if self.editable:
      #menu = ['clear','view','annotate','sep','copy','paste','help']
      menu = ['clear',['View','quick_view','view_text','view_CCP4mg','view_Coot','view_PdbView'],'sep','copy','paste','help']
      if self.ifAtomSelection: menu.insert(menu.index('sep'),'select')
    else:
      if self.role is not None and self.role == 'output':
        menu = [['View','quick_view','view_text','view_CCP4mg','view_Coot','view_PdbView'],'sep','copy','editLabel','export','help']
      else:
        menu = [['View','quick_view','view_text','view_CCP4mg','view_Coot','view_PdbView'],'sep','copy','export','help']
    if self._stacked: menu.insert(0,'handleStack')
    if self.ifInfo: menu.insert(menu.index('sep'),'annotate')
    return menu
  
  def getActionDef(self,name):
    if name == 'select':
      def iC():
        try:
            return self.widgets['selection'].isVisible()
        except:
            return False
      return dict (
        text = self.tr("Select atoms"),
        tip = self.tr('Select limited set of atoms'),
        slot = self.showAtomSelection,
        checkable = True,
        checked = iC
      )
    elif name == 'quick_view':
      def e():  return (self.model is not None and self.model.exists())
      return dict (
        text = self.tr("Quick view"),
        tip = self.tr('Quick view of model data'),
        slot = self.viewContents,
        enabled = e
      )
      
    else:
      return CCP4Widgets.CDataFileView.getActionDef(self,name)

  def openViewer(self,mode):
    if mode == 'view_text':
      CCP4Modules.WEBBROWSER().openFile(fileName=self.model.__str__(),toFront=True)
    else:
      CCP4Widgets.CDataFileView.openViewer(self,mode)

  def setModel(self,model=None):
    #print 'CPdbDataFileView.setModel',repr(model),self.widgets.has_key('selection')
    #if model is not None: print model.objectName()
    if self.model is not None and 'selection' in self.widgets:
        try:
            self.model.dataChanged.disconnect(self.widgets['selection'].applySelection)
        except:
            print("CPdbDataFileView.setModel self.model.dataChanged.disconnect failed")
    
    super(CPdbDataFileView,self).setModel(model=model)
    
    if self.ifAtomSelection and model is not None:
        if 'selection' in self.widgets:
            self.widgets['selection'].setModel(model.selection)
            self.widgets['selection'].updateViewFromModel()
    elif model is not None and model.selection.isSet():
        if 'selection' not in self.widgets: self.showAtomSelection()
        self.widgets['selection'].setModel(model.selection)
        self.widgets['selection'].updateViewFromModel()
            
  def getValue(self):
    val = CCP4Widgets.CDataFileView.getValue(self)
    if getattr(self,'selectionLine',None) is not None:
      val['selection'] = self.widgets['selection'].getValue()
    #print 'CPdbDataFileView.getValue',val
    return val

  def setValue(self,value):
    #print 'CPdbDataFileView.setValue',self.model.objectName(),value
    CCP4Widgets.CDataFileView.setValue(self,value)
    if value.get('selection',None) is not None:
      self.showAtomSelection()
      self.widgets['selection'].setValue(value.get('selection',{}))
    elif  getattr(self,'selectionLine',None) is not None:
      self.widgets['selection'].clear()
    else:
      return
    self.widgets['selection'].applySelection()
      

  def showAtomSelection(self):
    #print 'CPdbDataFileView.showAtomSelection'
    hideSimpleSelections = True
    if self.widgets.get('selection',None) is None:
      self.selectionLine = QtWidgets.QFrame(self)
      self.selectionLine.setLayout(QtWidgets.QHBoxLayout())
      self.selectionLine.layout().setSpacing(0)
      self.selectionLine.layout().setContentsMargins(0,0,0,0)
      self.selectionLine.layout().addWidget(QtWidgets.QLabel('Atom selection'))
      selection = None
      if hasattr(self,'model') and self.model is not None:
          selection = self.model.selection
      self.widgets['selection'] = CAtomSelectionView(parent=self,model=selection, qualifiers= { 'editable' :self.editable} )
      
      if hasattr(self,'model') and self.model is not None:
          self.model.dataChanged.connect(self.widgets['selection'].applySelection)
      self.selectionLine.layout().addWidget(self.widgets['selection'])
      self.layout().addWidget(self.selectionLine)

      self.simpleSelectionsCheck = QtWidgets.QCheckBox('Simple selections ...')
      self.layout().addWidget(self.simpleSelectionsCheck)

      expanded_cb_path = os.path.abspath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','toc-minus.png'))
      unexpanded_cb_path = os.path.abspath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','toc-plus.png'))

      expandStyle = """
      QCheckBox::indicator:unchecked {
         image: url("""+unexpanded_cb_path+""");
      }

      QCheckBox::indicator:checked {
          image: url("""+expanded_cb_path+""");
      }
      """
      self.simpleSelectionsCheck.setStyleSheet(expandStyle)

      self.selectionTypes = QtWidgets.QFrame(self)
      self.layout().addWidget(self.selectionTypes)
      self.chainsWidget = QtWidgets.QFrame(self)
      self.layout().addWidget(self.chainsWidget)
      self.excludeEdit = QtWidgets.QLineEdit(self)
      if self.editable:
          self.selectionTypes.setLayout(QtWidgets.QGridLayout())
          self.chainsWidget.setLayout(QtWidgets.QGridLayout())
          self.selectionTypes.layout().addWidget(QtWidgets.QLabel('Residue types'),0,0)
          peptideWidget = QtWidgets.QCheckBox('Peptide')
          nucleicWidget = QtWidgets.QCheckBox('Nucleic acid')
          ligandWidget = QtWidgets.QCheckBox('Ligands')
          waterWidget = QtWidgets.QCheckBox('Water')
          soluteWidget = QtWidgets.QCheckBox('Solute')
          saccharideWidget = QtWidgets.QCheckBox('Saccharide')
          nucleotideWidget = QtWidgets.QCheckBox('Nucleotide')
          metalWidget = QtWidgets.QCheckBox('Metal')
          @QtCore.Slot()
          def interpretSimpleSelections():
              sel = ''
              if peptideWidget.isChecked():
                  sel += 'peptide'
              if nucleicWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or nucleic'
                  else:
                      sel = 'nucleic'
              if ligandWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or ligands'
                  else:
                      sel = 'ligands'
              if waterWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or water'
                  else:
                      sel = 'water'
              if soluteWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or solute'
                  else:
                      sel = 'solute'
              if saccharideWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or saccharide'
                  else:
                      sel = 'saccharide'
              if nucleotideWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or nucleotide'
                  else:
                      sel = 'nucleotide'
              if metalWidget.isChecked():
                  if len(sel)>0:
                      sel += ' or metal'
                  else:
                      sel = 'metal'

              chains = self.chainsWidget.findChildren(QtWidgets.QCheckBox)
              chainsel = ''
              for ch in chains:
                  #print ch.isChecked(), ch.text()
                  if ch.isChecked():
                      if len(chainsel)>0:
                          chainsel += " or " + str(ch.text())+"/"
                      else:
                          chainsel = str(ch.text())+"/"
              #print chainsel; sys.stdout.flush()

              if len(chainsel)>0 and len(sel)>0:
                 sel = "{"+sel+"} and {"+chainsel+"}"
              elif len(sel)==0 and len(chainsel)>0:
                 sel = chainsel

              if len(str(self.excludeEdit.text()).strip())>0:
                  theElements = str(self.excludeEdit.text()).strip().split(",")
                  theNewElements = []
                  for ele in theElements:
                      if not ele in theNewElements:
                          theNewElements.append(ele)
                      if ele.upper() != ele:
                          if not ele.upper() in theNewElements:
                              theNewElements.append(ele.upper())
                  theNewElements = ",".join(theNewElements)

                  if len(sel)>0:
                      sel += ' and {not /*/*/*/*['+theNewElements+']:*}'
                  else:
                      sel = '{not /*/*/*/*['+theNewElements+']:*}'

              val = {'selection':{'text':sel}}
              CCP4Widgets.CDataFileView.setValue(self,val)
              self.widgets['selection'].setValue(val)
              self.widgets['selection'].applySelection()
              # Surely should not have to do this? But applySelection0 seems not to work...
              self.widgets['selection'].connectUpdateViewFromModel(False)
              self.widgets['selection'].model.text.set(sel)
              self.widgets['selection'].connectUpdateViewFromModel(True)
    
          peptideWidget.stateChanged.connect(interpretSimpleSelections)
          nucleicWidget.stateChanged.connect(interpretSimpleSelections)
          ligandWidget.stateChanged.connect(interpretSimpleSelections)
          waterWidget.stateChanged.connect(interpretSimpleSelections)
          soluteWidget.stateChanged.connect(interpretSimpleSelections)
          saccharideWidget.stateChanged.connect(interpretSimpleSelections)
          nucleotideWidget.stateChanged.connect(interpretSimpleSelections)
          metalWidget.stateChanged.connect(interpretSimpleSelections)
          self.selectionTypes.layout().addWidget(peptideWidget,1,0)
          self.selectionTypes.layout().addWidget(nucleicWidget,2,0)
          self.selectionTypes.layout().addWidget(ligandWidget,3,0)
          self.selectionTypes.layout().addWidget(waterWidget,4,0)
          self.selectionTypes.layout().addWidget(soluteWidget,1,1)
          self.selectionTypes.layout().addWidget(saccharideWidget,2,1)
          self.selectionTypes.layout().addWidget(metalWidget,3,1)
          self.selectionTypes.layout().addWidget(nucleotideWidget,4,1)

          def fillChainsWidget():
              #print "fillChainsWidget"; sys.stdout.flush()
              self.nChainsCol = 0
              if hasattr(self,"model") and hasattr(self.model,"fileContent") and hasattr(self.model.fileContent,"composition") and hasattr(self.model.fileContent.composition,"chains"):
                  if len(self.model.fileContent.composition.chains)<70: # and len(self.model.fileContent.composition.chains)>1:
                      self.chainsWidget.layout().addWidget(QtWidgets.QLabel('Chains'),0,0)
                      for i in range(len(self.model.fileContent.composition.chains)):
                          col = i // 8
                          row = i % 8 + 1
                          chainWidget = QtWidgets.QCheckBox(self.model.fileContent.composition.chains[i])
                          self.chainsWidget.layout().addWidget(chainWidget,row,col)
                          chainWidget.stateChanged.connect(interpretSimpleSelections)
                          self.nChainsCol = col

          def fillExcludeWidget():
              self.excludeLabel = QtWidgets.QLabel('Exclude (comma separated) atoms of type:')
              self.excludeEdit.setToolTip("Comma separated list of atoms to not include, e.g. H,C,N,Ag")
              #This line adds a blank middle column which is somewhat more visually appealing (to SJM) than the long text entry box.
              for i in range(self.nChainsCol):
                  self.chainsWidget.layout().setColumnStretch(i,0)
              self.chainsWidget.layout().setColumnStretch(self.nChainsCol,10)
              self.chainsWidget.layout().addWidget(self.excludeLabel,0,self.nChainsCol+1)
              self.chainsWidget.layout().addWidget(self.excludeEdit,1,self.nChainsCol+1)
              self.excludeEdit.textEdited.connect(interpretSimpleSelections)

          fillChainsWidget()
          fillExcludeWidget()

          def clearLayout(layout):
              #print "clearLayout"; sys.stdout.flush()
              #print self.model.fileContent.molHnd
              while layout.count():
                  child = layout.takeAt(0)
                  if child.widget():
                      child.widget().deleteLater()

          @QtCore.Slot()
          def dataIsChanged():
              #print "dataIsChanged"; sys.stdout.flush()
              if not self.model.fileContent.molHnd is self._old_molHnd:
                  clearLayout(self.chainsWidget.layout())
                  fillChainsWidget()
                  self.excludeEdit = QtWidgets.QLineEdit(self)
                  fillExcludeWidget()
                  self._old_molHnd = self.model.fileContent.molHnd

          if hasattr(self,"model") and hasattr(self.model,"fileContent") and hasattr(self.model.fileContent,"composition") and hasattr(self.model.fileContent.composition,"chains"):
              self._old_molHnd = self.model.fileContent.molHnd
          else:
              self._old_molHnd = None

          self.model.dataChanged.connect(dataIsChanged)
 
    elif not self.widgets['selection'].isVisible():# and self.simpleSelectionsCheck.isChecked():
      self.selectionLine.show()
      if hasattr(self,"simpleSelectionsCheck"):
          self.simpleSelectionsCheck.show()
          if self.simpleSelectionsCheck.isChecked():
              hideSimpleSelections = False
              self.selectionTypes.show()
              self.chainsWidget.show()
    else:
      self.selectionLine.hide()
      self.selectionTypes.hide()
      self.chainsWidget.hide()
      if hasattr(self,"simpleSelectionsCheck"):
          self.simpleSelectionsCheck.hide()

    if not self.editable:
      if hasattr(self,"excludeEdit"):
          self.excludeEdit.hide()
      if hasattr(self,"excludeLabel"):
          self.excludeLabel.hide()
      self.selectionTypes.hide()
      self.chainsWidget.hide()
      self.simpleSelectionsCheck.hide()

    self.simpleSelectionsCheck.stateChanged.connect(self.selectionTypes.setVisible)
    self.simpleSelectionsCheck.stateChanged.connect(self.chainsWidget.setVisible)
    if hideSimpleSelections:
        self.selectionTypes.hide()
        self.chainsWidget.hide()

  def viewContents(self):
    if not hasattr(self,'contentsViewer'):
      self.contentsViewer = CPdbContentsViewer(self)
      try:
        fileContent = self.model.fileContent
      except:
        pass
      self.model.fileContent.dataChanged.connect(self.contentsViewer.load)
      self.contentsViewer.load()
    self.contentsViewer.show()
    
  def setFocus1(self,reason,index=None):
    #print 'CPdbDataFileView.setFocus1',index
    if index is None:
      CCP4Widgets.CDataFileView.setFocus(self,reason)
    elif index==1 and self.ifAtomSelection:
      self.showAtomSelection()
      self.selectionLine.setFocus(reason)
    else:
      self.jobCombo.setFocus(reason)

  def handleBrowserOpenFile(self,filename,downloadInfo,**kw):
    kw['validate'] = False
    kw['updateView'] = False
    CCP4Widgets.CDataFileView.handleBrowserOpenFile(self,filename,downloadInfo,**kw)
    err = self.model.importFile(jobId=self.parentTaskWidget().jobId(),jobNumber=self.parentTaskWidget().jobNumber())
    if err is not None and len(err)>0:
      if err.maxSeverity()>SEVERITY_WARNING:
        message = 'File failed validation test'
        err.warningMessage(windowTitle='Importing coordinate file',parent=self,message=message)
        self.model.unSet()
    self.updateViewFromModel()
    self.validate()
    


class CPdbContentsViewer(QtWidgets.QDialog):
  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.textEdit = QtWidgets.QTextEdit(self)
    self.textEdit.setReadOnly(True)
    self.layout().addWidget(self.textEdit)
    butBox = QtWidgets.QDialogButtonBox(self)
    self.layout().addWidget(butBox)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Close)
    but.clicked.connect(self.close)

  @QtCore.Slot()
  def load(self):
    #print 'CPdbContentsViewer.load', self.parent().model.fileContent.composition.chains
    self.textEdit.setReadOnly(False)
    self.textEdit.clear()
    if self.parent().model.annotation.isSet():
      self.setWindowTitle(self.parent().model.annotation.__str__())
    else:
      self.setWindowTitle(self.parent().model.__str__())
    text = '\nSpace group: ' + str(self.parent().model.fileContent.mmdbManager.GetSpaceGroup())
    try:
      cell = self.parent().model.fileContent.mmdbManager.GetCell()
      text = text + '\nCell:'
      for ii in range(1,7):
        text = text + '  ' + str(cell[ii])
    except:
      pass
     
    if len( self.parent().model.fileContent.composition.chains)>0:
      text = text + '\n\nChains in model:\n'
      comp = self.parent().model.fileContent.composition
      for i in range(len(comp.chains)):
        text = text + str(comp.chains[i])+ '  ' + comp.chainInfo[i][1].split('/')[3] + ' - ' +  comp.chainInfo[i][2].split('/')[3] + \
           '  ('+str(comp.chainInfo[i][0]) +' residues)\n'
    
    if len( self.parent().model.fileContent.composition.monomers)>0:
      text = text + '\n\nMonomers in model:\n'+self.parent().model.fileContent.composition.monomers[0]
      for mon in self.parent().model.fileContent.composition.monomers[1:]:
        text = text + ', ' + str(mon)
    self.textEdit.document().setPlainText(text)
    self.textEdit.setReadOnly(True)
  
class CPdbDataFileListView(CCP4Widgets.CListView):
  MODEL_CLASS = CCP4ModelData.CPdbDataFileList
  def __init__(self,parent=None,model=None,qualifiers={}):

    #print 'CPdbDataFileListView'
    qualis = { 'mode' : 'table',
               'tableItems' : ['fullPath','selection'] ,
               'columnHeaders':['Coordinate file','Selection'],
               }
    qualis.update(qualifiers)
    CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)
    # Tis broke! self.editor.showAtomSelection()

class CAtomSelectionView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4ModelData.CAtomSelection
  ERROR_CODES = { 101 : { 'description' : 'No CPdbDataFile when applying atom selection' },
                  102 :  { 'description' : 'Error applying atom selection' }
                  }

  def __init__(self,parent=None,model=None,qualifiers={}):
    super(CAtomSelectionView,self).__init__(parent=parent,qualifiers=qualifiers)

    self.allowUndefined = False
    
    if self.editable:
      self.widgets['text'] = CCP4Widgets.CLineEdit(self)
      self.widgets['text'].setToolTip("Enter selection command e.g. 'A/ or B/11-21'")
      self.setFocusProxy(self.widgets['text'])
    else:
      self.widgets['text'] = CCP4Widgets.CLabel(self)
    self.layout().addWidget(self.widgets['text'])
    self.label = QtWidgets.QLabel(' ( 0 atoms)',self)
    self.layout().addWidget(self.label)
    self.help = QtWidgets.QPushButton('Help',self)
    self.help.setToolTip('Details of atom selection syntax')
    self.layout().addWidget(self.help)
    if hasattr(self.widgets['text'],"textEdited"):
        self.widgets['text'].textEdited.connect(self.applySelection0)
    self.help.clicked.connect(self.showHelp)
    if self.parent().model is not None and hasattr(self.parent().model,"dataLoaded"):
        self.parent().model.dataLoaded.connect(self.applySelection)
    if model is not None:
      self.setModel(model)
      self.updateViewFromModel()

  def validate(self,isValid=None,reportMessage=True):
      self.applySelection()

  @QtCore.Slot()
  def showHelp(self):
    CCP4Modules.WEBBROWSER().loadWebPage(helpFileName='general/atom_selection.html')


  def setValue(self,value={}):
    CCP4Widgets.CComplexLineWidget.setValue(self,value=value)
    try:
      self.applySelection()
    except:
      self.label.setText('No atoms')


  def updateViewFromModel(self):
    #print 'CAtomSelectionView.updateViewFromModel'
    #import traceback
    #traceback.print_stack(limit=10)
    self.widgets['text'].blockSignals(True)
    self.widgets['text'].setText(self.model.__str__())
    self.widgets['text'].blockSignals(False)
    try:
      self.applySelection()
    except:
      self.label.setText('No atoms')
      
  def clear(self):
    CCP4Widgets.CComplexLineWidget.clear(self)
    self.widgets['text'].clear()
    self.label.setText('No atoms')

  @QtCore.Slot(str)
  def  applySelection0(self,text=None):    
    try:
      self.applySelection(text=text)
    except:
      pass
    else:
      self.connectUpdateViewFromModel(False)
      self.model.text.set(text)
      self.connectUpdateViewFromModel(True)
    
  @QtCore.Slot(str)
  def applySelection(self,text=None):
    #print 'CSelectionLine.applySelection',text
    #import traceback
    #traceback.print_stack(limit=5)
    if text is None: text = self.widgets['text'].text()
    if isinstance(self.model.parent(),CCP4ModelData.CPdbDataFile):
      pdbDataObj = self.model.parent()
    else:
      pdbDataObj = self.model.getDataByKey('pdbFileKey')
    #print 'CSelectionLine.applySelection',text,pdbDataObj
    if pdbDataObj is None:
      try:
        pdbDataObj = self.model.parent()
      except:
        self.label.setText('No atoms')
        raise CException(self.__class__,101,name=self.modelObjectPath())
    if pdbDataObj.isSet():
      try:
        # Call to loadFile removed from here - CPdbDataFile.updateData() unsets fileContent
        # so this accessing fileContent should cause a loadFile() if necessary
        nSelAtoms,selHnd = pdbDataObj.fileContent.interpretSelection(str(text))
        if selHnd is not None: pdbDataObj.fileContent.molHnd.DeleteSelection(selHnd)
      except Exception as e:
        #print 'CSelectionLine.applySelection loadFile fail\n',e
        #self.label.setText('No atoms')
        #raise
        #raise CException(self.__class__,102,name=self.modelObjectPath())
        nSelAtoms = 0
    else:
      nSelAtoms = 0

    self.label.setText(' ( '+str(nSelAtoms)+' atoms)')
    
    if (nSelAtoms == 0 and not self.allowUndefined) or (nSelAtoms == 0 and pdbDataObj.isSet()):
      self.setProperty("isValid",False)
      self.setProperty("hasWarning",True)
      self.isValid = False
      self.validityMessage = 'No atoms selected'
    else:
      self.setProperty("isValid",True)
      self.setProperty("hasWarning",False)
      self.isValid = True
      self.validityMessage = None
    self.updateValidityIndicator()



class CDictDataDialog(QtWidgets.QDialog):
  def __init__(self,parent=None,model=None,projectId=None):
    if parent is None:
      parent = CCP4Modules.QTAPPLICATION()
    QtWidgets.QDialog. __init__(self,parent)
    self.setLayout(QtWidgets.QHBoxLayout())
    self.setModal(False)
    if projectId is not None:
      pName = ' for '+CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectname')
    else:
      pName=''
    self.setWindowTitle('Manage project dictionary'+pName)
    
      
    self.frame = CDictDataView(self,model=model,projectId=projectId)
    self.layout().addWidget(self.frame)

  def closeEvent(self,event):
    #print 'CDictDataDialog.closeEvent',self.frame.showWidget
    if self.frame.showWidget is not  None and CCP4Utils.isAlive(self.frame.showWidget):
      CCP4Modules.WEBBROWSER().tab().deleteTabWidget(widget=self.frame.showWidget)
    event.accept()

class CDictDataView(CCP4Widgets.CViewWidget):

  def __init__(self,parent=None,model=None,projectId=None):
    CCP4Widgets.CViewWidget.__init__(self,parent=parent)
    self.showWidget = None
    if model is None:
      self.dictDataFile = CCP4ModelData.CDictDataFile()
      #self.dictDataFile.set(self.dictDataFile.defaultProjectDict(projectId=projectId))
      self.dictDataFile.loadFile()
      self.model = self.dictDataFile.fileContent
      from qtgui import CCP4ProjectViewer
      CCP4ProjectViewer.FILEWATCHER().addJobPath(os.path.split( os.path.split(str(self.dictDataFile))[0])[1],self.dictDataFile.__str__())
      CCP4ProjectViewer.FILEWATCHER().fileChanged.connect(self.handleFileChanged)
    else:
      self.model = model
      self.dictDataFile = model.parent()
    #print 'CDictDataView',self.dictDataFile,self.model,self.model.monomerList
    QtWidgets.QFrame.__init__(self,parent)
    self.setLayout(QtWidgets.QHBoxLayout())
    self.layout().setContentsMargins(1,1,1,1)
    self.layout().setSpacing(1)
    self.monomerListWidget = CDictDataList(self)   
    self.layout().addWidget(self.monomerListWidget )
    butLayout = QtWidgets.QVBoxLayout()
    for label,connect in [['Show',self.handleShow],['Import geometry file',self.handleMerge],['Delete',self.handleDelete]]:
      but = QtWidgets.QPushButton(label,self)
      but.clicked.connect(connect)
      butLayout.addWidget(but)
    self.layout().addLayout(butLayout)
    self.updateViewFromModel()
    self.model.dataChanged.connect(self.updateViewFromModel)
    self.monomerListWidget.itemDoubleClicked.connect(self.handlDoubleClick)

  @QtCore.Slot(str)
  def handleFileChanged(self,path):
    #print 'CDictDataView.handleFileChanged',path
    if str(path) == str(self.dictDataFile):
      #print 'CDictDataView.handleFileChanged updating'
      self.dictDataFile.fileContent.loadFile(str(self.dictDataFile))
      self.updateViewFromModel()


  @QtCore.Slot('QListWidgetItem')
  def handlDoubleClick(self,item):
    #print 'CDictDataView.handlDoubleClick',item
    idd = item.data(0).toStrint().__str__()
    self.handleShow(idd)

    
  @QtCore.Slot()
  def updateViewFromModel(self):
    self.monomerListWidget.load(self.model.monomerList)
    if self.showWidget is not  None and CCP4Utils.isAlive(self.showWidget):
      self.showWidget.reload()
      
  @QtCore.Slot()
  def updateModelFromView(self):
    # All updates should be handled by the edit methods
    pass

  @QtCore.Slot()
  def handleShow(self,idd=None):
    if idd is None:
      idd = self.monomerListWidget.currentSelection()
    #print 'CDictDataView.handleShow',idd,self.model,self.model.parent()
    if self.showWidget is not None: print('alive?',CCP4Utils.isAlive(self.showWidget))
    if self.showWidget is None or (not CCP4Utils.isAlive(self.showWidget)):
      self.showWidget = CCP4Modules.WEBBROWSER().openFile(self.model.parent().__str__())
    else:
      self.showWidget.browserWindow().show()
      self.showWidget.browserWindow().raise_()
    if idd is not None:
      self.showWidget.findText(subString='data_comp_'+idd)


  def handleMerge(self):
    from qtgui import CCP4FileBrowser
    self.mergeFileBrowser = CCP4FileBrowser.CFileDialog(parent=self,title='Select geometry file to merge into project geometry file',
          defaultSuffix=CCP4Modules.MIMETYPESHANDLER().getMimeTypeInfo(name='application/refmac-dictionary',info='fileExtensions'),
          filters = ['Geometry file for refinement(*.cif)'], fileMode=QtWidgets.QFileDialog.ExistingFiles,saveButtonText='Merge these files')
    self.mergeFileBrowser.show()
    self.mergeFileBrowser.selectFiles.connect(self.mergeFiles)

  @QtCore.Slot(list)
  def mergeFiles(self,selectedFiles):
    if hasattr(self,'mergeFileBrowser'):
      self.mergeFileBrowser.close()
      self.mergeFileBrowser.deleteLater()
    #print 'CDictDataFileView.mergeFiles',selectedFiles
    
    for fileName in selectedFiles:
      err = self.model.mergeFile(fileName=fileName,overwrite=True)
      #print 'CDictDataFileView.mergeFiles',err.report(),err.maxSeverity()
      if err.maxSeverity()>SEVERITY_WARNING:
        err.warningMessage(windowTitle='Error merging geometry files',parent=self,message='Error attempting to merge geometry files')
      return

  @QtCore.Slot()
  def handleDelete(self):
    idd = self.monomerListWidget.currentSelection()
    #print 'CDictDataView.handleDelete',idd
    if idd is None: return
    err = self.model.delete(idd)
    if err.maxSeverity()>SEVERITY_WARNING:
      err.warningMessage(windowTitle='Error deleting item in geometry file',parent=self,message='Error attempting to edit geometry file')
    return


class CDictDataList(QtWidgets.QTreeWidget):
  def __init__(self,parent):
    QtWidgets.QTreeWidget.__init__(self,parent)
    self.setMinimumWidth(400)
    self.setColumnCount(3)
    self.setHeaderLabels(['Id','Code','Name'])
    self.setColumnWidth(0,80)
    self.setColumnWidth(1,80)
    self.setColumnWidth(2,200)
    self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
                          
  @QtCore.Slot(list)
  def load(self,monomerList):
    # monomerList is a Clist of CChemComp
    #print 'CDictDataList.load',monomerList
    self.clear()
    for monomer in monomerList:
      #print 'CDictDataList.load',monomer.get('id')
      qList = []
      for item in ['id', 'three_letter_code', 'name']:
          qList.append(str(monomer.get(item)))
      item = QtWidgets.QTreeWidgetItem(qList)
      treeId=self.addTopLevelItem(item)

  def handleMonomerDeleted(self,id):
    #print 'CDictDataList.handleMonomerDeleted',id
    if self.model().rowCount()==0: return
    modInxList = self.model().match(id)
    if len(modInxList)==0:
      #print 'CDictDataList.handleMonomerDeleted no match to id',id
      pass
    else:
      self.model().removeRow(modInxList[0].row())

  def handleMonomerAdded(self,id):
    pass
      
  def currentSelection(self):
    indices = self.selectionModel().selectedRows()
    if len(indices) == 0:
        return None
    idd = str(indices[0].data())
    return idd

class CSeqEdit(CCP4Widgets.CTextEdit):

  def __init__(self,parent):
    CCP4Widgets.CTextEdit.__init__(self,parent)
    self.setStyleSheet("CSeqEdit { font-family: Courier ;}")
    
  def setValue(self,text):
    self.blockSignals(True)
    if text is None:
      CCP4Widgets.CTextEdit.setValue(self,'')
      return
    text = re.sub(' ','',text)
    text = re.sub('\n','',text)
    nGap = 10
    nLine = 60
    seqList = []
    while len(text)>0:
      line = text[0:nLine]
      text = text[nLine:]
      while len(line)>0:
        seqList.append(line[0:nGap]+' ')
        line = line[nGap:]
      seqList[-1] = seqList[-1]+'\n'
    newText = ''
    newText = newText.join(seqList)
    CCP4Widgets.CTextEdit.setValue(self,newText)
    self.blockSignals(False)

class PositiveIntegerDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        validator = QtWidgets.QRegExpValidator(QtCore.QRegExp("^[1-9]+\d*"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class SpaceDelimitedPositiveIntegerDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("(^[1-9]+\d*)([ ][1-9]+\d*)+"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class AlphanumericDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None,maxchars=0):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)
        self.maxchars = maxchars

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        regexp = "[A-Za-z0-9]"
        if self.maxchars > 0:
           regexp = "[A-Za-z0-9]{1,"+str(self.maxchars)+"}"
        validator = QtGui.QRegExpValidator(QtCore.QRegExp(regexp), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class SpaceDelimitedAlphanumericDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None,maxchars=0):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)
        self.maxchars = maxchars

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        regexp = "[A-Za-z0-9]([ ][A-Za-z0-9])*"
        if self.maxchars > 0:
           regexp = "[A-Za-z0-9]{1,"+str(self.maxchars)+"}([ ][A-Za-z0-9]{1,"+str(self.maxchars)+"})*"
        validator = QtGui.QRegExpValidator(QtCore.QRegExp(regexp), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class ResnumDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("-?[0-9]+[A-Za-z]?"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class SingleCharDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("[A-Za-z0-9]?"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class CAtomRefmacSelectionListView(CCP4Widgets.CViewWidget):

    MODEL_CLASS = CCP4ModelData.CAtomRefmacSelectionList

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CViewWidget.__init__(self, parent, qualis)
        self.STRETCH = 2
        self.setLayout(QtWidgets.QHBoxLayout())
        line = QtWidgets.QVBoxLayout()
        self.widgets['atomSelTable'] = CCP4Widgets.CMultiAtomSelection(self)
        self.widgets['atomSelTable'].setToolTip("Manually specify rigid group selections, by chain identifier and residue range.\nMultiple ranges can be part of the same group, by assigning them the same integer Group ID.\nIf groups are not manually specified then each chain (or segment) in the model will be treated as a seperate rigid group.")

        groupidDelegate = PositiveIntegerDelegate(self.widgets['atomSelTable'].table)
        chainDelegate = AlphanumericDelegate(self.widgets['atomSelTable'].table,3)
        res1Delegate = ResnumDelegate(self.widgets['atomSelTable'].table)
        res2Delegate = ResnumDelegate(self.widgets['atomSelTable'].table)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(0,groupidDelegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(1,chainDelegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(2,res1Delegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(3,res2Delegate)

        line.addWidget(self.widgets['atomSelTable'])
        self.layout().addLayout(line)
        self.model = model
        #privateModel = CCP4RefmacMultiAtomSelection.AtomSelectionListModel()
        privateModel = CCP4RefmacMultiAtomSelection.RigidBodySelectionListModel()
        self.widgets['atomSelTable'].setModel(privateModel)
        privateModel.setSelectionData(model)
        def myDataChanged(tl,br):
            #print(privateModel.selectionData())
            self.model.blockSignals(True)
            self.model.set(privateModel.selectionData())
            self.model.blockSignals(False)
            self.model.emitDataChanged()

        privateModel.dataChanged.connect(myDataChanged)

class COccRefmacSelectionListView(CCP4Widgets.CViewWidget):

    MODEL_CLASS = CCP4ModelData.COccRefmacSelectionList

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CViewWidget.__init__(self, parent, qualis)
        self.STRETCH = 2
        self.setLayout(QtWidgets.QHBoxLayout())
        line = QtWidgets.QVBoxLayout()
        self.widgets['atomSelTable'] = CCP4Widgets.CMultiAtomSelection(self)
        self.widgets['atomSelTable'].setToolTip("Manually specify partial occupancy group selections, by chain identifier, residue range and optionally atom name and/or alt code.\nMultiple chains can be specified by providing a space-seperated list.\nMultiple selections can be part of the same group, by assigning them the same integer Group ID.")
        
        groupidDelegate = PositiveIntegerDelegate(self.widgets['atomSelTable'].table)
        chainsDelegate = SpaceDelimitedAlphanumericDelegate(self.widgets['atomSelTable'].table,3)
        res1Delegate = ResnumDelegate(self.widgets['atomSelTable'].table)
        res2Delegate = ResnumDelegate(self.widgets['atomSelTable'].table)
        atomDelegate = AlphanumericDelegate(self.widgets['atomSelTable'].table,4)
        altDelegate = SingleCharDelegate(self.widgets['atomSelTable'].table)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(0,groupidDelegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(1,chainsDelegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(2,res1Delegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(3,res2Delegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(4,atomDelegate)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(5,altDelegate)

        line.addWidget(self.widgets['atomSelTable'])
        self.layout().addLayout(line)
        self.model = model
        privateModel = CCP4RefmacMultiAtomSelection.OccupancySelectionListModel()
        self.widgets['atomSelTable'].setModel(privateModel)
        privateModel.setSelectionData(model)
        def myDataChanged(tl,br):
            #print(privateModel.selectionData())
            self.model.blockSignals(True)
            self.model.set(privateModel.selectionData())
            self.model.blockSignals(False)
            self.model.emitDataChanged()

        privateModel.dataChanged.connect(myDataChanged)

class COccRelationRefmacListView(CCP4Widgets.CViewWidget):

    MODEL_CLASS = CCP4ModelData.COccRelationRefmacList

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CViewWidget.__init__(self, parent, qualis)
        self.STRETCH = 2
        self.setLayout(QtWidgets.QHBoxLayout())
        line = QtWidgets.QVBoxLayout()
        self.widgets['atomSelTable'] = CCP4Widgets.CMultiAtomSelection(self)
        self.widgets['atomSelTable'].setToolTip("Specify which partial occupancy groups are mutually exclusive / correspond to conformers that exist in different parts of the crystal.")
        
        groupidsDelegate = SpaceDelimitedPositiveIntegerDelegate(self.widgets['atomSelTable'].table)
        self.widgets['atomSelTable'].table.setItemDelegateForColumn(0,groupidsDelegate)

        line.addWidget(self.widgets['atomSelTable'])
        self.layout().addLayout(line)
        self.model = model
        privateModel = CCP4RefmacMultiAtomSelection.OccupancyRelationListModel()
        self.widgets['atomSelTable'].setModel(privateModel)
        privateModel.setSelectionData(model)
        def myDataChanged(tl,br):
            #print(privateModel.selectionData())
            self.model.blockSignals(True)
            self.model.set(privateModel.selectionData())
            self.model.blockSignals(False)
            self.model.emitDataChanged()

        privateModel.dataChanged.connect(myDataChanged)

class CAsuContentSeqListView(CCP4Widgets.CViewWidget):

    MODEL_CLASS = CCP4ModelData.CAsuContentSeqList

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CViewWidget.__init__(self, parent, qualis)
        self.STRETCH = 2
        self.setLayout(QtWidgets.QHBoxLayout())
        line = QtWidgets.QVBoxLayout()
        self.widgets['seqTable'] = CCP4Widgets.CSequenceTable(self)
        self.widgets['seqTable'].setToolTip("Sequence list view. Click '+' to add a sequence.\n\nOr drag and drop sequence/coordinate files into here.")
        line.addWidget(self.widgets['seqTable'])
        self.layout().addLayout(line)
        self.model = model
        self.redraw()
        self.widgets['seqTable'].loadFromSeqFileSignal.connect(self.actOnLoadSeqFile)
        self.widgets['seqTable'].loadFromCoorFileSignal.connect(self.actOnLoadCoorFile)
        self.widgets['seqTable'].urlsDroppedSignal.connect(self.actOnUrlsDropped)
        self.widgets['seqTable'].textDroppedSignal.connect(self.actOnTextDropped)

    def actOnTextDropped(self,text):
#
        try:
            fileInfo = CCP4Modules.PROJECTSMANAGER().db().getFileInfo(text,mode=['projectid','relpath','filename','annotation','filesubtype','filecontent','filetype'])
            projectDir = CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectdirectory')
            filePath = os.path.join(projectDir,fileInfo['relpath'],fileInfo['filename'])
            url = QtCore.QUrl.fromLocalFile(filePath)
            self.actOnUrlsDropped([url])
        except:
            print(text,"is not a database id for a CCP4i2 input/output file")
            raise

    def setValue(self, value=None):
        self.redraw()

    def redraw(self):
        privateModel = CCP4SequenceList.SequenceModel()
        self.widgets['seqTable'].setModel(privateModel)
        for mod in self.model:
            privateModel.addItem((str(mod.name),int(mod.nCopies),str(mod.description),str(mod.polymerType),str(mod.sequence)))
        self.validate()
        @QtCore.Slot('QModelIndex','QModelIndex')
        def privateModelDataChanged(tl,br):
            newValues = []
            for i in range(privateModel.rowCount()):
                name, copies, desc, seqType, seq = privateModel.getRow(i)
                newSeq = CCP4ModelData.CAsuContentSeq()
                newSeq.sequence.set(str(seq))
                newSeq.name.set(CCP4Utils.safeOneWord(str(name)))
                newSeq.description.set(str(desc))
                newSeq.nCopies.set(int(copies))
                newSeq.polymerType.set(str(seqType))
                newValues.append(newSeq)
            self.model.blockSignals(True)
            self.model.set(newValues)
            self.model.blockSignals(False)
            self.model.dataChanged.emit()
            self.validate()
        privateModel.dataChanged.connect(privateModelDataChanged)
        if not self.editable:
            privateModel.setEditable(False)

    @QtCore.Slot(list)
    def actOnUrlsDropped(self,urls):
        for url in urls:
            if url.isLocalFile():
                fname = url.toLocalFile()
                if fname.lower().endswith(".pdb") or fname.lower().endswith(".cif") or fname.lower().endswith(".mmcif") or fname.lower().endswith(".ent"):
                    pdbFile = CCP4ModelData.CPdbData(parent=self)
                    pdbFile.loadFile(fname)
                    self.loadPdbFile(pdbFile,os.path.abspath(fname),os.path.basename(fname))
                elif "." in fname.lower() and (fname[fname.rfind(".")+1:]).lower() in CCP4ModelData.EXTLIST:
                    seqFile = CCP4ModelData.CSeqDataFile(parent=self.model)
                    seqFile.setFullPath(fname)
                    seqFile.loadFile()
                    if seqFile.__dict__['format'] == 'unknown':
                        QtWidgets.QMessageBox.warning(self,'Reading sequence file','The format of the file has not been recognised.' )
                    if len(seqFile.__dict__.get('identifiers',[])) == 0:
                        return
                    elif len(seqFile.__dict__.get('identifiers',[])) == 1:
                        entry = CCP4ModelData.CAsuContentSeq()
                        entry.nCopies.set(1)
                        entry.sequence.set(seqFile.fileContent.sequence)
#FIXME - Do not understand why my data (rcsb 1ibm) is invalid! It is valid in import sequence task!
                        print(seqFile.fileContent.identifier)
                        print(type(seqFile.fileContent.identifier))
                        entry.name.set(CCP4Utils.safeOneWord(str(seqFile.fileContent.identifier)))
                        entry.description.set(seqFile.fileContent.description)
                        entry.autoSetPolymerType()
                        self.model.blockSignals(True)
                        self.model.append(entry)
                        self.model.blockSignals(False)
                        self.redraw()
                    else:
                        win = QtWidgets.QDialog()
                        title = 'Select one or more sequences from the file'
                        label = QtWidgets.QLabel(title, self)
                        layout = QtWidgets.QVBoxLayout()
                        win.setLayout(layout)
                        idChooser = QtWidgets.QListWidget(self)
                        layout.addWidget(label)
                        layout.addWidget(idChooser)
                        idChooser.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
                        row = 0
                        for item in seqFile.__dict__['identifiers']:
                            idChooser.addItem(item)
                            idChooser.item(row).setSelected(True)
                            row += 1
                        dbb = QtWidgets.QDialogButtonBox()
                        closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Ok)
                        discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Cancel)
                        closeButton.clicked.connect(win.accept)
                        discardButton.clicked.connect(win.reject)
                        layout.addWidget(dbb)
                        ret = win.exec_()
                        if ret == QtWidgets.QDialog.Accepted:
                            recordList = []
#FIXME Need to call autoSetPolymerType() for each entry
                            for item in idChooser.selectedItems():
                                recordList.append(idChooser.row(item))
                            rv = self.model.extendSeqList(seqFile, recordList)
                            self.redraw()
                else:
                    QtWidgets.QMessageBox.warning(self,'Unhandled file suffix','File name suffix unrecognised as coordinate or sequence data.')
            else:
                QtWidgets.QMessageBox.warning(self,'Unhandled URL type','Only local files can be dropped here at present.')
                

    @QtCore.Slot()
    def actOnLoadCoorFile(self):
#TODO - Make OK button greyed if nothing selected ...
        win = QtWidgets.QDialog(self)
        win.setModal(True)
        layout = QtWidgets.QVBoxLayout(win)
        win.setLayout(layout)
        dbb = QtWidgets.QDialogButtonBox()
        closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Ok)
        discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Cancel)
        closeButton.clicked.connect(win.accept)
        discardButton.clicked.connect(win.reject)
        pdbFile = CCP4ModelData.CPdbDataFile(parent=self)
        pdbWidget = CPdbDataFileView(self,model=pdbFile,qualifiers={'iconButton':True,'iconName':'PdbDataFile','ifAtomSelection' : False})
        layout.addWidget(pdbWidget)
        layout.addWidget(dbb)
        win.setMinimumWidth(400)
        win.show()
        @QtCore.Slot()
        def onCoorFileAccept():
            fullPath = pdbFile.fullPath
            label = pdbFile.guiLabel()
            self.loadPdbFile(pdbFile.fileContent,fullPath,label)
        win.accepted.connect(onCoorFileAccept)

    @QtCore.Slot(str,str,str)
    def loadPdbFile(self, pdbFile,fullPath,label):
#TODO - Deal With k being " ".
        if pdbFile is None or pdbFile.sequences is None:
            return
        fname = os.path.basename(str(fullPath))
        if fname.lower().endswith(".pdb"):
            fname = fname[:-4]
        elif fname.lower().endswith(".cif"):
            fname = fname[:-4]
        elif fname.lower().endswith(".mmcif"):
            fname = fname[:-6]
        elif fname.lower().endswith(".ent"):
            fname = fname[:-4]
        if len(pdbFile.sequences) == 0:
            QtWidgets.QMessageBox.warning(self,'Reading sequence from coordinate file','No sequences were read from the file\n')
            return
        elif len(pdbFile.sequences)==1:
            k = (list(pdbFile.sequences.keys()))[0]
            v = (list(pdbFile.sequences.values()))[0]
            entry = CCP4ModelData.CAsuContentSeq()
            entry.nCopies.set(1)
            entry.sequence.set(v)
            entry.name.set(CCP4Utils.safeOneWord(fname+"_"+k))
            entry.description.set(label+"_"+k)
            entry.autoSetPolymerType()
            self.model.blockSignals(True)
            self.model.append(entry)
            self.model.blockSignals(False)
            self.redraw()
            self.model.dataChanged.emit()
        else:
            idChooser = QtWidgets.QDialog(self)
            layout = QtWidgets.QVBoxLayout()
            idChooser.setLayout(layout)
            layout.addWidget(QtWidgets.QLabel("Please select one or more sequences to import"))
            dbb = QtWidgets.QDialogButtonBox()
            closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Ok)
            discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Cancel)
            closeButton.clicked.connect(idChooser.accept)
            discardButton.clicked.connect(idChooser.reject)
            listWidget = QtWidgets.QListWidget()
            listWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            items = {}
            for k,v in pdbFile.sequences.items():
                items[k] = listWidget.addItem("Chain "+k)
                #items[k].setSelected(True)
            layout.addWidget(listWidget)
            layout.addWidget(dbb)
            ret = idChooser.exec_()
            if ret == QtWidgets.QDialog.Rejected:
                return
            else:
                irow=0
                for k,v in pdbFile.sequences.items():
                    if listWidget.item(irow).isSelected():
                        entry = CCP4ModelData.CAsuContentSeq()
                        entry.nCopies.set(1)
                        entry.sequence.set(pdbFile.sequences[k])
                        entry.name.set(CCP4Utils.safeOneWord(fname+"_"+k))
                        entry.description.set(label+"_"+k)
                        entry.autoSetPolymerType()
                        self.model.blockSignals(True)
                        self.model.append(entry)
                        self.model.blockSignals(False)
                    irow += 1
                self.redraw()
                self.model.dataChanged.emit()

    @QtCore.Slot()
    def actOnLoadSeqFile(self):
#TODO - DEAL With multiple sequences in seq file in a more elegant way.
        win = QtWidgets.QDialog(self)
        layout = QtWidgets.QVBoxLayout(win)
        win.setLayout(layout)
        dbb = QtWidgets.QDialogButtonBox()
        closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Ok)
        discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Cancel)
        closeButton.clicked.connect(win.accept)
        discardButton.clicked.connect(win.reject)
        seqFile = CCP4ModelData.CSeqDataFile(parent=self.model)
        seqWidget = CSeqDataFileView(self,model=seqFile,qualifiers={'iconButton':True,'iconName':'SeqDataFile'})
        @QtCore.Slot(list)
        def loadMultipleFromFile(recordList):
            win.close()
            rv = self.model.extendSeqList(seqFile, recordList)
            self.redraw()
        seqWidget.loadMultiple.connect(loadMultipleFromFile)
        layout.addWidget(seqWidget)
        layout.addWidget(dbb)
        win.setMinimumWidth(400)
        win.setModal(True)
        win.show()
        @QtCore.Slot()
        def onSeqFileAccept():
            entry = CCP4ModelData.CAsuContentSeq()
            entry.nCopies.set(1)
            entry.sequence.set(seqFile.fileContent.sequence)
#FIXME - Do not understand why my data (rcsb 1ibm) is invalid! It is valid in import sequence task!
            entry.name.set(CCP4Utils.safeOneWord(str(seqFile.fileContent.identifier)))
            entry.description.set(seqFile.fileContent.description)
            entry.autoSetPolymerType()
            self.model.blockSignals(True)
            self.model.append(entry)
            self.model.blockSignals(False)
            self.redraw()
            self.model.dataChanged.emit()
            
        win.accepted.connect(onSeqFileAccept)

    def validate(self,isValid=None,reportMessage=True):
        if hasattr(self,"widgets") and "seqTable" in self.widgets and hasattr(self.widgets["seqTable"],"table") and hasattr(self.widgets["seqTable"].table,"model") and hasattr(self.widgets["seqTable"].table.model(),"rowCount"):
            if self.widgets['seqTable'].table.model().rowCount() == 0:
                return CCP4Widgets.CViewWidget.validate(self,isValid=False,reportMessage="No sequences in AU contents")
            else:
                return CCP4Widgets.CViewWidget.validate(self,isValid=True,reportMessage="Sequences OK")
        return CCP4Widgets.CViewWidget.validate(self,isValid=True,reportMessage="Sequences OK")

class CAsuContentSeqView(CCP4Widgets.CViewWidget):

    MODEL_CLASS = CCP4ModelData.CAsuContentSeq
    ERROR_CODES = {101 : { 'description' : 'There is no valid sequence name'},
                   102 : { 'description' : 'There is no valid sequence'}}
 
    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CViewWidget.__init__(self, parent, qualis)
        print(self.widgets)

class CAsuDataFileView(CCP4Widgets.CDataFileView):
    # Ok, this is the class that gets drawn for the ASU datatype. It also includes the modal interface code....
    MODEL_CLASS = CCP4ModelData.CAsuDataFile
    NOTICE_OLD = ("CCP4i2 now uses a 'AU content' file that contains all chain sequences in "
              "a crystal instead of separate sequence files. Where necessary you will be able "
              "to select chain(s) from the file. In future projects please use the "
              "' Define AU contents/' (in the ' Import ../' module). "
              "Please enter all known chain sequences for your structure.")
    NOTICE = ("It is recommended that you enter an abbreviated\ndescription of your asymmetric unit below in the\n"
              "text box (this will appear in drop down menus).\nA generic description will be "
              "generated if left blank.")
    ANNO = 'ASU Data File'

    def __init__(self, parent=None, model=None, qualifiers={}):
        self.selectionMode = qualifiers.get('selectionMode', 0)
        qualis = { 'vboxLayout' : True, 'browseFiles' : False }
        qualis.update(qualifiers)
        self.selTextList = []
        super(CAsuDataFileView, self).__init__(parent=parent, model=None, qualifiers=qualis)
        self.infoline = None
        if self.editable:
            #FIXME - Need to provide an API to hide this.
            self.icona = QtGui.QIcon(CCP4Modules.PIXMAPMANAGER().getPixmap("list_add_grey"))
            self.iconm = QtGui.QIcon(CCP4Modules.PIXMAPMANAGER().getPixmap("list_delete_grey"))
            self.buttonAS = QtWidgets.QToolButton(self)
            self.buttonAS.setIcon(self.icona)
            self.buttonAS.setMaximumHeight(16)
            self.buttonAS.setMaximumWidth(16)
            #button.setText(item[1])
            self.buttonAS.setToolTip("Unlock button for Quick ASU Create (can also be done through task menu & Create ASU Contents)")
            # The original part
            self.defineWidget = QtWidgets.QWidget()
            defineLayout = QtWidgets.QHBoxLayout()
            self.defineWidget.setLayout(defineLayout)
            defineLayout.setContentsMargins(0, 0, 0, 0)
            self.layout().addWidget(self.defineWidget)
            self.defineButton = QtWidgets.QPushButton("New AU contents")
            self.defineButton.setToolTip("Open the dialog to specify AU contents information")
            defineLayout.addWidget(self.buttonAS)
            self.labhlp = QtWidgets.QLabel("If a suitable ASU is not available above, you can press the cross & then button to quickly create one.")
            defineLayout.addWidget(self.labhlp)
            defineLayout.addWidget(self.defineButton)
            self.defineButton.hide()
            defineLayout.addStretch(2)
            # KJS : The orig passed a blank string (using the partial) in order to to pop a modal gui up in an
            # overridden function which is used for browsing files. Unclear why.
            self.defineButton.clicked.connect(functools.partial(self.handleBrowserOpenFile,"",""))
            self.buttonAS.clicked.connect(self.HideB)
        if self.selectionMode > 0:
            if self.selectionMode == 1:
                self.selectionLabel = QtWidgets.QLabel('Select one sequence', self)
            else:
                self.selectionLabel = QtWidgets.QLabel('Select one or more sequences', self)
            self.selectionLabel.setObjectName('italic')
            self.layout().addWidget(self.selectionLabel)
            self.selectionFrame = QtWidgets.QFrame(self)
            self.selectionFrame.setLayout(QtWidgets.QGridLayout())
            self.layout().addWidget(self.selectionFrame)
            self.selectionGroup = QtWidgets.QButtonGroup(self)
            if self.selectionMode == 1:
                self.selectionGroup.setExclusive(True)
            else:
                self.selectionGroup.setExclusive(False)
            self.selectionButtons = []
            self.selectionGroup.buttonClicked[int].connect(self.updateModelSelection)
            if self.selectionMode == 2:
                self.selectionFrame.hide()
                self.selectionLabel.hide()
        self.setModel(model)

    @QtCore.Slot()
    def HideB(self):
        if self.defineButton.isHidden():
            self.buttonAS.setIcon(self.iconm)
            self.defineButton.show()
            self.labhlp.hide()
        else:
            self.buttonAS.setIcon(self.icona)
            self.defineButton.hide()
            self.labhlp.show()

    def hideDefineWidget(self):
        if hasattr(self,"defineWidget"):
            self.defineWidget.hide()

    def setModel(self,model):
        if self.model is not None: 
            self.model.dataChanged.disconnect(self.drawSelection)
        CCP4Widgets.CDataFileView.setModel(self,model)
        if self.model is not None:
            self.model.dataChanged.connect(self.drawSelection)
            self.model.selection.dataChanged.connect(self.drawSelection)
      
    def filterText(self):
        # make the filters text for QFileDialog
        textList = CCP4Widgets.CDataFileView.filterText(self)  
        extList = CCP4ModelData.CSeqDataFile().qualifiers('fileExtensions')
        line = 'Sequence files (*.' + extList[0]
        for ext in extList:
            line = line + ' *.' + ext
        line += ')'
        textList.insert(0, line)
        return textList

    def dropTypes(self):
        # Enable drag-n-drop of old sequence data
        return ['AsuDataFile', 'SeqDataFile']

    @QtCore.Slot('QMimeData')
    def acceptDropData(self,textData):
        #print 'CAsuDataFileView.acceptDropData',textData
        if textData.count('application/CCP4-seq'):
          self.handleBrowserOpenFile(textData.split(' ')[0])
        else:
          from lxml import etree
          tree = etree.fromstring(textData)
          if tree.tag == 'SeqDataFile':
            fileName = os.path.join(tree.xpath('./relPath')[0].text,tree.xpath('./baseName')[0].text)
            self.handleBrowserOpenFile(fileName,None)
          else:
            CCP4Widgets.CDataFileView.acceptDropData(self,textData)

    def handleBrowserOpenFile(self, fileName, downloadInfo, autoInfoOnFileImport=True,
                              validate=True,updateView=True):
        if fileName.endswith('asu.xml'):
            CCP4Widgets.CDataFileView.handleBrowserOpenFile(self, fileName, downloadInfo,
                                                            autoInfoOnFileImport=autoInfoOnFileImport,
                                                            validate=validate, updateView=updateView)
        else:
            #projectId = self.parentTaskWidget().projectId()
            self.importDialog = QtWidgets.QDialog(self)
            self.importDialog.setObjectName('createAsu')
            self.importDialog.setLayout(QtWidgets.QVBoxLayout())
            self.importDialog.setWindowTitle('Define AU content - import sequences:'+os.path.split(fileName)[1])
            self.importDialog.setModal(False)
            self.importAsuContent = CCP4ModelData.CAsuContentSeqList(parent=self)
            self.importAsuContentView = CAsuContentSeqListView(parent=self.importDialog,model=self.importAsuContent)
            self.importAsuContentView.setVisible(True)
            self.importAsuContent.set([])
            self.importAsuContentView.setMinimumWidth(550)
            self.importAsuContentView.updateViewFromModel()
            self.importDialog.layout().addWidget(self.importAsuContentView)
            #self.importDialog.layout().addWidget(QtWidgets.QLabel(self.NOTICE,self.importDialog))
            self.importDialog.show()

            #Need a CMtzDataFileView for Matthews calc...
            from core import CCP4XtalData
            from qtgui import CCP4XtalWidgets
            scl = CCP4Widgets.CLabel(self)
            scl.setText("<b>Solvent content analysis</b>")
            self.hklFile = CCP4XtalData.CMtzDataFile(parent=self.importAsuContentView.model)
            hklinWidget = CCP4XtalWidgets.CMtzDataFileView(parent=self,model=self.hklFile,qualifiers={'iconButton':True,'iconName':'MTZDataFile'})
            self.importDialog.layout().addWidget(scl)
            self.importDialog.layout().addWidget(hklinWidget)
            resultWidget = QtWidgets.QTextEdit(self)
            resultWidget.setMinimumHeight(100)
            resultWidget.setReadOnly(True)
            self.importDialog.layout().addWidget(resultWidget)
            resultWidget.setHtml("<p>Solvent content analysis will appear here when there is a valid sequence list and reflection file above.</p></body></html>")
            def runMatthews():
                 if (self.importAsuContent.validity(self.importAsuContent.get()).maxSeverity() <= SEVERITY_WARNING) and self.hklFile.isSet() and len(self.importAsuContent) > 0:
                     totWeight = 0.0
                     text = "<table><tr><th style=\"text-align:left\">Name</th><th>Number of copies</th><th>Molecular weight</th></tr>"
                     for seqObj in self.importAsuContent:
                         totWeight = totWeight + seqObj.molecularWeight(seqObj.polymerType)
                         text = text + '<tr><td>{0}</td><td>{1} </td><td>{2:.2f}</td></tr>'.format(seqObj.name,seqObj.nCopies,float(seqObj.molecularWeight(seqObj.polymerType))) + '\n'
                     text += "</table><br/>"

                     if totWeight < 1e-6:
                         resultWidget.setHtml('<b style="color: red;">Matthew\'s analysis suggests the current composition and cell volume are incompatible.</b>.')
                         return
                     text +=  '<br>Total sequence weight: {0:.2f}</br>'.format(float(totWeight))
                     rv = self.hklFile.fileContent.matthewsCoeff(molWt=totWeight)
                     vol = rv.get('cell_volume','Unkown')
                     nmol=[]
                     solv = []
                     matt=[]
                     prob=[]
                     if vol == 'Unkown':
                         text = text + '<p>Cell volume = Unknown</p>'
                     else:
                         text = text + '<p>Cell volume = {0:.1f}<p>\n'.format(float(vol)) 
                     headText = ""
                     if len(rv.get('results',[])) == 0:
                         headText = '<b style="color: red;">Matthews analysis suggests the current composition and cell volume are incompatible.</b>'
                         self.matthewsInvalid = True
                     else:
                         headText = headText +    '<table><tr><th>  Nmol  </th><th>  %solvent  </th><th>  Matthews  </th><th>  prob(Matthews)  </th></tr>\n'
                     for result in rv.get('results',[]):
                         headText = headText + '<tr><td>  {0}  </td><td>  {1:.2f}  </td><td>  {2:.2f}  </td><td>  {3:.2f}  </td></tr>'.format(result.get('nmol_in_asu'),result.get('percent_solvent'),result.get('matth_coef'),result.get('prob_matth')) + '\n'
                         nmol.append(result.get('nmol_in_asu'))
                         solv.append(result.get('percent_solvent'))
                         matt.append(result.get('matth_coef'))
                         prob.append(result.get('prob_matth'))
                     headText = headText +    '</table><br/>'
                     text = headText + text
                     resultWidget.setHtml(text)
                 else:
                     resultWidget.setHtml("<p>Solvent content analysis will appear here when there is a valid sequence list and reflection file above.</p></body></html>")

            self.hklFile.dataChanged.connect(runMatthews)
            self.importAsuContent.dataChanged.connect(runMatthews)

            qled = QtWidgets.QLineEdit(self)
            qled.textChanged[str].connect(self.infochan)
            #self.importDialog.layout().addWidget(qled)
            butBox = QtWidgets.QDialogButtonBox(self)
            self.importDialog.layout().addWidget(butBox)
            but = butBox.addButton(QtWidgets.QDialogButtonBox.Save)
            but.clicked.connect(self.doImport)
            but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
            but.clicked.connect(self.importDialog.hide)
            but = butBox.addButton(QtWidgets.QDialogButtonBox.Help)
            but.clicked.connect(self.showHelp)

    @QtCore.Slot(str)
    def infochan(self, text):
      self.infoline = text

    @QtCore.Slot()
    def doImport(self):
        #print 'CAsuDataFileView.doImport validate',self.importAsuContent,self.importAsuContent.validity(self.importAsuContent.get())
        #print 'CAsuDataFileView.doImport jobCombo',self.jobCombo.count()
        if self.importAsuContent.validity(self.importAsuContent.get()).maxSeverity() <= SEVERITY_WARNING:
            self.importDialog.hide()
            from dbapi import CCP4DbUtils
            openJob = CCP4DbUtils.COpenJob(projectId=self.parentTaskWidget().projectId())
            openJob.createJob(taskName='ProvideAsuContents')
            print(openJob.container)
            jobId = CCP4Modules.PROJECTSMANAGER().db().getJobId(projectId=self.parentTaskWidget().projectId(),jobNumber=openJob.jobNumber)
            openJob.container.inputData.ASU_CONTENT.set(self.importAsuContent)
            if self.infoline:
                ntxt = str(self.infoline)
                CAsuDataFileView.ANNO = ntxt
                openJob.container.guiAdmin.jobTitle.set(ntxt)
            # Do a cleanup of the sequences before saving
            for seqObj in openJob.container.inputData.ASU_CONTENT:
                clean,err = seqObj.cleanupSequence(str(seqObj.sequence))
                seqObj.sequence.set(clean)
                if hasattr(seqObj,"autoSetPolymerType"):
                    seqObj.autoSetPolymerType()
            if self.hklFile.isSet():
                openJob.container.inputData.HKLIN.set(self.hklFile)
            err = openJob.runJob()
            if err.maxSeverity() > SEVERITY_WARNING:
                err.warningMessage(windowTitle='Importing sequence file to AU contents',parent=self,message='Failed creating AU content file')
                return
            self.resetJobCombo = openJob.jobNumber
            print("End of doImport", openJob.container, openJob.container)
            if self.infoline:
                print(ntxt, openJob.container.guiAdmin)
                CCP4Modules.PROJECTSMANAGER().db().updateJob(jobId=jobId, key='jobTitle', value=ntxt)

    def handleJobFinished(self, args):
        # After the CDataFileView.handleJobFinished() has updated the jobCombo make sure the new asu file is selected in it
        CCP4Widgets.CDataFileView.handleJobFinished(self,args)
        #print 'CAsuDataFileView.handleJobFinished',getattr(self,'resetJobCombo',None)
        if getattr(self,'resetJobCombo',None) is not None:
            idx = self.jobCombo.findText(self.resetJobCombo,QtCore.Qt.MatchStartsWith)
            #print 'CAsuDataFileView Attempting to reset resetJobCombo',idx
            if idx >= 0:
                self.jobCombo.setCurrentIndex(idx)
                del self.resetJobCombo
    
    def getMenuDef(self):
        menu = ['clear','view','exportSeq','editLabel','sep','copy','paste','help']
        # For now try selection always open
        if self.selectionMode>0: menu.insert(menu.index('sep'),'select')
        return menu
  
    def getActionDef(self,name):
        if name == 'select':
            def iC():
                try:
                    return self.selectionList.isVisible()
                except:
                    return False
            return dict(text=self.tr("Select chains"), tip=self.tr('Select limited set of sequences'),
                        slot=self.showSelection, checkable=True, checked=iC)
        elif name == 'exportSeq':
            return dict(text=self.tr("Export sequence file"), tip=self.tr('Export sequence file'),
                        slot=self.exportSeq)
        elif name == 'editLabel':
            def e(): 
                  return hasattr(self, "jobLabel")
            return dict(text=self.tr("Edit label"), tip=self.tr('Edit file label'),
                        slot=self.handleEditLabel, enabled=e,)
        return CCP4Widgets.CDataFileView.getActionDef(self, name)

    def openViewer(self, mode):
        if not self.model.exists():
            return
        from qtgui import CCP4TextViewer
        text = ''
        for seqObj in self.model.fileContent.seqList:
            text += str(seqObj.nCopies) + ' copies of ' + str(seqObj.name) + '\n'
            text += 'Description: ' + str(seqObj.description) + '\n'
            if seqObj.source.isSet():
                text += 'Loaded from file: ' + str(seqObj.source)  + '\n'
            text += '\n' + seqObj.formattedSequence() + '\n\n'
        self.viewDialog = QtWidgets.QDialog(self.parentTaskWidget())
        self.viewDialog.setLayout(QtWidgets.QVBoxLayout())
        self.viewDialog.setWindowTitle(str(self.model.annotation))
        self.viewDialog.setMinimumWidth(600)
        self.viewDialog.setMinimumHeight(300)
        viewer = CCP4TextViewer.CTextViewer(parent=self.viewDialog)
        viewer.loadText(text,fixedWidthFont=True)
        self.viewDialog.layout().addWidget(viewer)
        self.viewDialog.show()
        self.viewDialog.raise_()
    
    def showSelection(self, visible=None):
        if visible is None:
            visible = (not self.selectionFrame.isVisible())
        self.selectionFrame.setVisible(visible)
        self.selectionLabel.setVisible(visible)
    
    def getValue(self):
        rv = CCP4Widgets.CDataFileView.getValue(self)
        if self.selectionMode>0:
            selDict = {}
            for idx in range(len(self.model.fileContent.seqList)):
                selDict[str(self.model.fileContent.seqList[idx].name)] = self.selectionGroup.button(idx).isChecked()
            rv['selection'] = selDict
            #print 'CAsuDataFileView.getValue',rv['selection']
        return rv

    def setValue(self,value):
        #print 'CAsuDataFileView.setValue',value['selection']
        CCP4Widgets.CDataFileView.setValue(self, value)
        if self.selectionMode == 0:
            return
        self.selectionGroup.blockSignals(True)
        for idx in range(len(self.model.fileContent.seqList)):
            self.selectionGroup.button(idx).setChecked(value['selection'][str(self.model.fileContent.seqList[idx].name)])
        self.selectionGroup.blockSignals(False)

    def handleJobComboChange(self, indx0=None):
        CCP4Widgets.CDataFileView.handleJobComboChange(self, indx0=indx0)
        self.updateSelectionViewFromModel()
        self.validate()
    
    @QtCore.Slot()
    def updateModelFromView(self):
        CCP4Widgets.CDataFileView.updateModelFromView(self)
        if self.model is None:
            return
        if self.selectionMode == 0:
            return
        self.updateModelSelection()

    @QtCore.Slot()
    def updateModelSelection(self):
        selDict = {}
        if len(self.selectionButtons) == 0:
            return
        self.connectUpdateViewFromModel(False)
        for idx in range(len(self.model.fileContent.seqList)):
            selDict[str(self.model.fileContent.seqList[idx].name)] = self.selectionGroup.button(idx).isChecked()
        self.model.selection.set(selDict)
        self.connectUpdateViewFromModel(True)

    def updateViewFromModel(self):
        #print 'CAsuFileView.updateViewFromModel'
        CCP4Widgets.CDataFileView.updateViewFromModel(self)
        self.updateSelectionViewFromModel()

    def updateSelectionViewFromModel(self):
        self.drawSelection()
        if self.model is None:
            return
        #print 'CAsuFileView.updateViewFromModel nResidues',self.model.numberOfResidues(),self.model.molecularWeight()
        if self.selectionMode == 1 and len(self.model.fileContent.seqList) > 1:
            self.showSelection(True)
        elif self.selectionMode == 2:
            for key, value in list(self.model.selection.items()):
                #print 'CAsuFileView.updateViewFromModel selection',key,value
                if not value:
                    self.showSelection(True)
                    break

    @QtCore.Slot()
    def drawSelection(self):
        if  self.selectionMode == 0:
            return
        # Is the selection text up-to-date?
        selTextList = []
        if self.model is not None and self.model.isSet():
            for asuObject in self.model.fileContent.seqList:
                selTextList.append(str(asuObject.name) + '    ' +  str(asuObject.description))
        #print 'drawSelection',selTextList,selTextList == self.selTextList
        if selTextList == self.selTextList:
            return
        self.selectionGroup.blockSignals(True)
        # Delete current selection buttons
        for but in self.selectionButtons:
            self.selectionGroup.removeButton(but)
            but.deleteLater()
        self.selectionButtons = []
        # Update the selTextList and redraw the buttons
        self.selTextList = selTextList
        idx = 0
        for text in self.selTextList:
            try:
                if self.selectionMode == 1:
                    self.selectionButtons.append(QtWidgets.QRadioButton(text, self))
                else:
                    self.selectionButtons.append(QtWidgets.QCheckBox(text, self))
                self.selectionGroup.addButton(self.selectionButtons[-1], idx)
                self.selectionFrame.layout().addWidget(self.selectionButtons[-1], idx, 0)
                self.selectionGroup.button(idx).setChecked(self.model.selection[str(self.model.fileContent.seqList[idx].name)])
                idx += 1
            except Exception as e:
                print('Error drawing asu content selection', e)
        self.selectionGroup.blockSignals(False)

    def showHelp(self):
        CCP4Modules.WEBBROWSER().loadWebPage(helpFileName='model_data')

    def exportSeq(self):
        from qtgui import CCP4FileBrowser
        filters = ["FASTA format (*.fasta *.fa *.ffn *.ffa *.fna *.frn)","EMBL format (*.embl)"]
        self.exportBrowser = CCP4FileBrowser.CFileDialog(self, title='Save selected sequences to file',
                                                         filters=filters, fileMode=QtWidgets.QFileDialog.AnyFile)
        self.exportBrowser.selectFile.connect(self.doExportSeq)
        self.exportBrowser.show()

    def getBioPyVersion(self):
        from Bio import __version__ as bioversion
        try:
            version=float(bioversion)
            return version
        except:
            return None

    @QtCore.Slot(str)
    def doExportSeq(self, fileName):
        theFilter = "fasta"
        if "." in fileName:
            theFilter = fileName[fileName.rfind(".") + 1:]
            if theFilter.lower().startswith("em"):
                theFilter = "embl"
            else:
                theFilter = "fasta"
        elif not CCP4Modules.PREFERENCES().NATIVEFILEBROWSER:
            if str(self.exportBrowser.widget.fileDialog.selectedFilter()).startswith("EMBL"):
                theFilter = "embl"
                fileName += ".embl"
            else:
                theFilter = "fasta"
                fileName += ".fasta"
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        bioversion=getBioPyVersion()
        if bioversion is not None:
            if bioversion < 1.79:
                from Bio.Alphabet import generic_protein
                from Bio.Alphabet import generic_nucleotide
        seqs = []
        for seqObj in self.model.fileContent.seqList:
            if bioversion is not None:
                if bioversion < 1.79:
                    seqs.append(SeqRecord(Seq(str(seqObj.sequence), generic_protein), name=str(seqObj.description), 
                                  id=str(seqObj.name), description=str(seqObj.description)))
                else:
                    seqs.append(SeqRecord(Seq(str(seqObj.sequence), id="prot", annotations={"molecule_type": "protein"}), name=str(seqObj.description), 
                                  id=str(seqObj.name), description=str(seqObj.description)))
            else:
                seqs.append(SeqRecord(Seq(str(seqObj.sequence), id="prot", annotations={"molecule_type": "protein"}), name=str(seqObj.description), 
                                  id=str(seqObj.name), description=str(seqObj.description)))
        with open(fileName,'w+') as outputFileHandle:
            SeqIO.write(seqs,outputFileHandle, theFilter)
    
