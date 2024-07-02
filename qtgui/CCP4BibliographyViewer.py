"""
     CCP4BibliographyViewer.py: CCP4 GUI Project
     Copyright (C) 2014STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.sstac
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Sept 2014 - View and export bibliography
"""

import functools
from PySide6 import QtGui, QtWidgets,QtCore
from core.CCP4ErrorHandling import *
from core import CCP4Annotation, CCP4Modules
from qtgui import CCP4Widgets

class CBibReferenceView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4Annotation.CBibReference

  
  def __init__(self,parent=None,model=None,qualifiers={}):
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers={'vboxLayout' : True })
    self.setMinimumWidth(500)
        
    # Beware changeof widget type may need change in CCP4StyleSheet
    #self.widgets['selection'] = QtWidgets.QCheckBox(self)
    self.widgets['selection'] = QtWidgets.QLabel(self)
    self.widgets['selection'].setObjectName('bibreference_selection')
    for item in ['authorList','source']:
      self.widgets[item] = QtWidgets.QLabel(self)
      self.widgets[item].setObjectName('bibreference_'+item)
      self.widgets[item].setWordWrap(True)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['selection'])
    #line.addWidget(self.widgets['title'])
    #line.addStretch(1)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['authorList'])
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['source'])
    self.layout().addLayout(line)

    self.setModel(model)

  def setModel(self,model):
    if model is None:
      for item in ['selection','authorList','source']:
        self.widgets[item].setText('')
    else:
      self.widgets['selection'].setText(str(model.title))
      self.widgets['source'].setText(str(model.source))
      if len(model.authorList)>0:
        authorText = str(model.authorList[0])
        for author in model.authorList[1:]: authorText = authorText + ',' + str(author)
        self.widgets['authorList'].setText(authorText)
    
class CBibReferenceGroupView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4Annotation.CBibReferenceGroup

  def __init__(self,parent=None,model=None,qualifiers={}):
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers={'vboxLayout' : True, 'iconName' : 'book' })
    iconWidgetItem = self.layout().takeAt(0)
    if iconWidgetItem is not None:
        line = QtWidgets.QHBoxLayout()        
        line.addWidget(iconWidgetItem.widget())
        self.layout().addLayout(line)

    self.widgets['title']  = QtWidgets.QLabel(self)
    line.addWidget(self.widgets['title'])
    self.referencesLayout = QtWidgets.QVBoxLayout()
    self.layout().addLayout(self.referencesLayout)

    self.setModel(model)
    self.fileBrowser = None

  def setModel(self,model):
    #CCP4Widgets.CComplexLineWidget.setModel(self,model)
    self.model = model
    if model is None:
      self.widgets['title'].setText('')
    else:
      self.widgets['title'].setText(str(model.title))
    if model is not None and len(model.references)>0:                                            
      for refObj in model.references:
        #print 'CBibReferenceGroupView.setModel',refObj                  
        widget = CBibReferenceView(parent=self,model=refObj)
        self.referencesLayout.addWidget(widget)
       
  def getMenuDef(self):
    #return ['select_all','deselect_all','export_medline','export_bibtxt']
    return ['export_medline','export_bibtex']
      
  def getActionDef(self,name):
    if name == 'select_all':
      return dict (
        text = self.tr("Select all"),
        tip = self.tr('Select all of the listed references'),
        slot = self.selectAll,
      )
    elif name == 'deselect_all':
      return dict (
        text = self.tr("Deselect all"),
        tip = self.tr('Deselect all of the listed references'),
        slot = functools.partial(self.selectAll,False),
      )
    elif name == 'export_medline':
      return dict (
        text = self.tr("Export Medline format"),
        tip = self.tr('Export all selected references in Medline format'),
        slot = functools.partial(self.handleExport,'medline')
      )
    elif name == 'export_bibtex':
      return dict (
        text = self.tr("Export Bibtex format"),
        tip = self.tr('Export all selected references in Bibtex format'),
        slot = functools.partial(self.handleExport,'bibtex')
      )

  def selectAll(self,mode=True):
    for i in range(self.referencesLayout.count()):
      self.referencesLayout.itemAt(i).widget().widgets['selection'].setChecked(mode)

  def handleExport(self,format=None):
    #for i in range(self.referencesLayout.count()):
    #  if self.referencesLayout.itemAt(i).widget().widgets['selection'].isChecked():
    if self.fileBrowser is None:
      from qtgui import CCP4FileBrowser
      self.fileBrowser = CCP4FileBrowser.CFileDialog(parent=self,
                                      title='Export bibliography',
                                      filters = ['.'+format+'.txt'],
                                      defsaultSuffix = format+'.txt',
                                      fileMode = QtWidgets.QFileDialog.AnyFile)
      self.fileBrowser.selectFile.connect(functools.partial(self.export,format))
    self.fileBrowser.show()

  @QtCore.Slot(str,str)
  def export(self,format,fileName):
    err = self.model.export(cformat=format,fileName=fileName)
    if err.maxSeverity()>SEVERITY_WARNING:
      err.warningMessage(parent=self,windowTitle='Error attempting to export bibliography')
    

class CBibliographyViewer(QtWidgets.QMainWindow):
  def __init__(self,parent=None):
    self.referenceGroups = []
    self.referenceGroupViews = []
    QtWidgets.QMainWindow.__init__(self,parent=parent)
    #self.scroll = QtWidgets.QScrollArea(self)
    self.frame = QtWidgets.QFrame(self)
    self.setCentralWidget(self.frame)
    #self.scroll.setWidget(self.frame)
    self.frame.setLayout(QtWidgets.QVBoxLayout())
    

  def setReferences(self,jobId=None,taskNameList=[]):
    import os,glob
    pm = CCP4Modules.PROJECTSMANAGER()
    if jobId is not None:
      taskBiblio = glob.glob(os.path.join(pm.jobDirectory(jobId=jobId),'*.medline.txt'))
      if len(taskBiblio)>0:
        taskName = pm.db().getJobInfo(jobId,'taskname')
        self.referenceGroups.append(CCP4Annotation.CBibReferenceGroup(self))
        self.referenceGroups[-1].loadFromMedline(fileNameList=taskBiblio,taskName=taskName)
        self.referenceGroupViews.append(CBibReferenceGroupView(self))
        self.referenceGroupViews[-1].setModel(self.referenceGroups[-1])
        self.frame.layout().addWidget(self.referenceGroupViews[-1])
        self.setWindowTitle('Bibliography for '+CCP4Modules.TASKMANAGER().getTitle(taskName))
        return
  
    if len(taskNameList) == 0: taskNameList.append(pm.db().getJobInfo(jobId=jobId,mode='taskname') )

    for taskName in taskNameList:
      self.referenceGroups.append(CCP4Annotation.CBibReferenceGroup(self))
      self.referenceGroups[-1].loadFromMedline(taskName=taskName)
      self.referenceGroupViews.append(CBibReferenceGroupView(self))
      #print 'setReferences referenceGroups',self.referenceGroups[-1]
      self.referenceGroupViews[-1].setModel(self.referenceGroups[-1])
      self.frame.layout().addWidget(self.referenceGroupViews[-1])
    self.setWindowTitle('Bibliography for '+CCP4Modules.TASKMANAGER().getTitle(taskNameList[0]))
    self.frame.show()
