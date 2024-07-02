from __future__ import print_function

"""
     CCP4ReportExternalManagerGui.py: CCP4 GUI Project
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
     Liz Potterton July 2013 - report jobs run external to ccp4i2
"""

import os, functools
from PySide6 import QtGui, QtWidgets,QtCore
from core import CCP4Data,CCP4Container
from qtgui import CCP4CustomisationGui,CCP4Widgets
from core import CCP4ImportedJobManager
from core.CCP4ErrorHandling import *
from core.CCP4Modules import IMPORTEDJOBMANAGER,WEBBROWSER,PROJECTSMANAGER

def openGui():
  if CImportedJobManagerGui.insts is None:
    CImportedJobManagerGui.insts = CImportedJobManagerGui()
  CImportedJobManagerGui.insts.show()
  CImportedJobManagerGui.insts.raise_()


class CImportedJobManagerGui(CCP4CustomisationGui.CCustomisationGui):

  insts = None
  def __init__(self,parent=None):
    CCP4CustomisationGui.CCustomisationGui.__init__(self,parent=parent,mode='importedjobs',title='Import Job Manager')
    self.createWidget= None
    

  def manager(self):
    return IMPORTEDJOBMANAGER()

  def handleNew(self):
    if self.createWidget is None:
      self.createWidget = CCreateImportedJobDialog(self)
      self.createWidget.importedJobCreated.connect(self.customListView.populate)
    else:
      self.createWidget.updateViewFromModel()
    self.createWidget.show()
    self.createWidget.raise_()


  def handleEdit(self,selected=None):
    if selected is None: selected = self.customListView.selectedItem()
    if selected is None: return
    editor = CCreateImportedJobDialog(self,name=selected)
    editor.show()

class CImportedJobDataView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4ImportedJobManager.CImportedJobData
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'vboxLayout' : True }
    qualis.update(qualifiers)
    #print 'CListView qualis',qualis
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setModel(model)
    line =  QtWidgets.QHBoxLayout()
    line.addWidget(self.layout().takeAt(0).widget())
    self.widgets['dataType'] = CCP4Widgets.CStringView(parent=self,model=self.model.dataType)
    self.model.dataType.dataChanged.connect(self.handleDataTypeChanged)
    line.addWidget(self.widgets['dataType'] )
    self.widgets['label'] = CCP4Widgets.CStringView(parent=self,model=self.model.label)
    line.addWidget(self.widgets['label'] )
    line.setStretchFactor(self.widgets['label'],5)
    self.layout().addLayout(line)
    
    self.fileNameLayout =  QtWidgets.QHBoxLayout()
    qualis = { 'autoInfoOnFileImport' : False, 'jobCombo' : False, 'browseDb' : False , 'toolTip' :'Select imported job directory' }
    self.widgets['fileName'] = CCP4Widgets.CDataFileView(parent=self,model=self.model.fileName,qualifiers=qualis)
    self.fileNameLayout.addWidget(self.widgets['fileName'])
    self.layout().addLayout(self.fileNameLayout)
    
    self.layout().addStretch(5)
    self.handleDataTypeChanged()
    self.model.dataChanged.connect(self.validate)


  @QtCore.Slot()
  def handleDataTypeChanged(self):
    from core import CCP4DataManager
    dataType = self.model.dataType.__str__()
    cls = CCP4DataManager.DATAMANAGER().getClass(dataType)
    if cls is None: return
    self.model.resetFileNameClass(cls)
    self.model.label.set(self.model.fileName.qualifiers('guiLabel'))
    
    widget = CCP4DataManager.DATAMANAGER().widget(model=self.model.fileName,parentWidget=self,qualifiers={'vboxLayout' : False }) 
    print('CImportedJobDefinitionView.handleDataTypeChanged',dataType,widget)
    if widget is not None:
      try:
        self.widgets['fileName'].deleteLater()
        self.fileNameLayout.takeAt(0)
      except:
        pass
      self.widgets['fileName'] = widget
      self.fileNameLayout.addWidget(widget)
    self.validate()  
    
class CImportedJobDataListView(CCP4Widgets.CListView):
  MODEL_CLASS = CCP4ImportedJobManager.CImportedJobDataList
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = { 'mode' : 'table',
               'tableItems' : ['name','dataType','label'] ,
               'columnHeaders':['Name','Data type','Label on interface'],
               }
    qualis.update(qualifiers)
    CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)
    
class CCreateImportedJobDialog(QtWidgets.QDialog):

  def __init__(self,parent=None,name=None):
    QtWidgets.QDialog.__init__(self,parent)

    self.setWindowTitle('Report an external job')
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setContentsMargins(2,2,2,2)
    self.layout().setSpacing(2)
    
    self.widgets = {}
    self.model = CCP4ImportedJobManager.CImportedJobDefinition(self,name=name)
    if name is not None:  
      self.model.loadDataFromXml(fileName=os.path.join(IMPORTEDJOBMANAGER().getDirectory(name),'task.xml'),loadHeader=True)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Name of imported job',self))
    self.widgets['name'] = CCP4Widgets.CStringView(parent=self,model=self.model.name)
    self.widgets['name'].setToolTip('Unique single-word identifier for task')
    line.addWidget(self.widgets['name'])
    line.addWidget(QtWidgets.QLabel('title',self))
    self.widgets['title'] = CCP4Widgets.CStringView(parent=self,model=self.model.title)
    self.widgets['title'].setToolTip('Task name that appears in the user interface')
    line.addWidget(self.widgets['title'])
    line.setStretchFactor(self.widgets['title'],5)
    self.widgets['title'].setMinimumWidth(300)
    line.addStretch(1)
    self.layout().addLayout(line)
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('Directory to copy to CCP4i2 project directory (optional)'))
    self.layout().addLayout(line)   
    line = QtWidgets.QHBoxLayout()
    self.widgets['commandFile'] = CCP4Widgets.CDataFileView(parent=self,model=self.model.commandFile)
    line.addWidget(self.widgets['commandFile'])
    self.layout().addLayout(line)

    self.inputFrame = QtWidgets.QFrame()
    self.inputFrame.setFrameShape(QtWidgets.QFrame.Box)
    self.inputFrame.setLayout(QtWidgets.QVBoxLayout())
    self.inputFrame.layout().setSpacing(1)
    self.inputFrame.layout().setContentsMargins(1,1,1,1)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('Input files to job'))
    self.inputFrame.layout().addLayout(line)
    print('CCreateImportedJobDialog.__init__',self.model,self.model.dataOrder())
    self.widgets['inputFileDefinitionList'] = CImportedJobDataListView(parent=self,model=self.model.inputFileDefinitionList)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['inputFileDefinitionList'])
    self.inputFrame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addStretch(4)
    but = QtWidgets.QPushButton('Add input file',self)
    but.clicked.connect(functools.partial(self.addDataObject,'input'))
    line.addWidget(but)
    line.addStretch(1)
    self.inputFrame.layout().addLayout(line)
    self.layout().addWidget(self.inputFrame)

    self.outputFrame = QtWidgets.QFrame()
    self.outputFrame.setFrameShape(QtWidgets.QFrame.Box)
    self.outputFrame.setLayout(QtWidgets.QVBoxLayout())
    self.outputFrame.layout().setSpacing(1)
    self.outputFrame.layout().setContentsMargins(1,1,1,1)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('Output files'))
    self.outputFrame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    self.widgets['outputFileDefinitionList'] = CImportedJobDataListView(parent=self,model=self.model.outputFileDefinitionList)
    line.addWidget(self.widgets['outputFileDefinitionList'])
    self.outputFrame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addStretch(4)
    but = QtWidgets.QPushButton('Add output file',self)
    but.clicked.connect(functools.partial(self.addDataObject,'output'))
    line.addWidget(but)
    line.addStretch(1)
    self.outputFrame.layout().addLayout(line)
    self.layout().addWidget(self.outputFrame)
    
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

  @QtCore.Slot(str)
  def addDataObject(self,mode):
    if mode == 'input':
      defList = self.model.inputFileDefinitionList
      name = 'inputFileDefinitionList'
      frame = self.inputFrame
    else:
      defList = self.model.outputFileDefinitionList
      name = 'outputFileDefinitionList'
      frame = self.outputFrame
    print('CCreateImportedJobDialog.addDataObject',mode,frame.layout().count(),defList.__len__())
    defList.addItem()
    indx = defList.__len__() - 1
    line = QtWidgets.QHBoxLayout()
    self.widgets[name+str(indx)] = CImportedJobDefinitionView(parent=self,model=defList[indx])
    line.addWidget(self.widgets[name+str(indx)])
    frame.layout().insertWidget(indx,self.widgets[name+str(indx)])

  def updateViewFromModel(self):
    for key in list(self.widgets.keys()):
      self.widgets[key].updateViewFromModel()
    
  @QtCore.Slot()
  def help(self):
    WEBBROWSER().loadWebPage(helpFileName='customisation')

  @QtCore.Slot()
  def cancel(self):
    self.hide()

  @QtCore.Slot()
  def accept(self):
    if not self.model.isValid(): return
    name = self.model.name.__str__()
    if os.path.exists(os.path.join(IMPORTEDJOBMANAGER().getDirectory(name))):
      QtWidgets.QMessageBox.warning(self,self.windowTitle(),'There is already an imported job called: '+name+'.  Please enter alternative name')
      return
    self.model.saveDataToXml(fileName=os.path.join(IMPORTEDJOBMANAGER().getDirectory(name),'task.xml'))
    self.hide()
