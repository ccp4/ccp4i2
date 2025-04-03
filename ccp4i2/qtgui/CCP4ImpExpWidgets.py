"""
Copyright (C) 2012 STFC
"""

# May 2012 Tools to handle import/export of files

import os

from PySide2 import QtCore, QtWidgets

from . import CCP4Widgets
from ..core import CCP4Annotation, CCP4Data, CCP4File
from ..core.CCP4DataManager import DATAMANAGER


class CExportedFileCombo(CCP4Widgets.CComplexLineWidget):

  MODEL_CLASS = CCP4File.CExportedFile

  def __init__(self,parent=None,model=None,qualifiers={}):
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualifiers=qualifiers)
    self.widget = CCP4Widgets.CComboBox(self)
    self.widget.setEditable(False)
    self.layout().addWidget(self.widget)
    self.fileTypeClass = qualifiers.get('fileTypeClass',None)
    self.before = qualifiers.get('before',None)
    self.projectId = None
    self.setProjectId(qualifiers.get('projectId',None))
    #print 'CExportedFileCombo',qualifiers    
    self.setModel(model)
    self.updateViewFromModel()
    self.widget.currentIndexChanged[int].connect(self.updateModelFromView0)

  def setProjectId(self,projectId=None):
    self.projectId = projectId
    self.widget.clear()
    if self.projectId is not None: self.load()

  def setValue(self,value):
    #print 'CExportedFileCombo.setValue',value
    expId = value.get('exportId',None)
    if expId is None or expId < 0:
      self.widget.setCurrentIndex(0)
    else:
      indx = self.widget.findData(expId)
      if indx>=0: self.widget.setCurrentIndex(indx)

  def load(self):
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    exportedFiles = PROJECTSMANAGER().db().getExportFilesByFileType(fileTypeClass=self.fileTypeClass,before=self.before,projectId=self.projectId)
    self.widget.addItem('',-1)
    for exFile in exportedFiles:
      from ..core.CCP4TaskManager import TASKMANAGER
      title = TASKMANAGER().getTitle(exFile['taskname'])
      if len(exFile['exportfilename'])>22:
        self.widget.addItem('..'+exFile['exportfilename'][-20:]+' '+exFile['jobnumber']+' '+title,exFile['exportid'])
      else:
        self.widget.addItem(exFile['exportfilename']+' '+exFile['jobnumber']+' '+title,exFile['exportid'])

  def getValue(self):
    #print 'CExportedFileCombo.getValue',self.widget.currentIndex(),item,type(item)
    indx = CCP4Data.varToUUID( self.widget.itemData(self.widget.currentIndex()) )
    return { 'exportId' : indx }

  @QtCore.Slot(int)
  def updateModelFromView0(self,indx):
    self.updateModelFromView()


class CExportedFileListView(CCP4Widgets.CListView):

  MODEL_CLASS = CCP4File.CExportedFileList
  def __init__(self,parent,model=None,qualifiers={},**kw):
    qualis = { 
               'editorClassName' : 'CExportedFileCombo',
               'dragType' : 'file_info' 
               }
    qualis.update(qualifiers)
    qualis.update(kw)
    CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis,editorQualifiers={'before':qualis.get('before',None)})


class CImportInfo(QtWidgets.QFrame):

  def __init__(self,parent,model=None,importId=None,importFileName=None,label=None,sourceFileAnnotation=''):
    QtWidgets.QFrame.__init__(self,parent)
    self.setMinimumWidth(500)
    self.model = model
    self.importId = importId
    if importId is not None:
      from ..core.CCP4ProjectsManager import PROJECTSMANAGER
      importFileInfo = PROJECTSMANAGER().db().getImportFileInfo(importId=importId)
      before = importFileInfo.get('creationtime')
      if importFileInfo.get('annotation',None) is not None:
        sourceFileAnnotation = sourceFileAnnotation + importFileInfo['annotation']
      sourceFileName = importFileInfo['sourcefilename']
    else:
      before = None
      if importFileName is not None:
        sourceFileName = importFileName
      else:
        try:
          sourceFileName = model.getSourceFileName()
        except:
          sourceFileName = model.__str__()
    #print 'CImportInfo',importId,before,sourceFileName
    try:
      self.parent().setWindowTitle('Provenance of file: '+os.path.split(sourceFileName)[1])
    except:
      pass
 
    self.setLayout(QtWidgets.QVBoxLayout())
    if label is None:
      label = 'Provide info on the origin of the file: '
    else:
      label = 'Provide info on the origin of '+label+' from the file:'
    self.layout().addWidget(QtWidgets.QLabel(label))
    self.layout().addWidget(QtWidgets.QLabel(sourceFileName))

    self.annotation = CCP4Annotation.CAnnotation(parent=self)
    self.annotation.text.set(sourceFileAnnotation)
    self.annotationView = DATAMANAGER().widget(model=self.annotation,parentWidget=self, qualifiers = { 'multiLine' : True, 'title' : 'Describe source of this file' } )
    #print 'CImportInfo.__init__',self.annotationView,self.annotationView.model,repr(self.annotation),self.annotation
    self.annotationView.updateViewFromModel()
    self.annotationView.setMaximumHeight(150)
    self.layout().addWidget(self.annotationView)

    self.exportedFiles = CCP4File.CExportedFileList(parent=self)
    self.exportedFilesView = DATAMANAGER().widget(model=self.exportedFiles,parentWidget=self , qualifiers = { 'title' :'Select any previously exported files used to derive this file', 'before' : before } )
    self.exportedFilesView.setMaximumHeight(150)
    self.layout().addWidget(self.exportedFilesView)

  def save(self):
    self.annotationView.updateModelFromView()
    #print 'CImportInfo.save',self.annotation.get('text').__str__()
    if self.importId is not None:
      if self.annotation.isSet():
        anno = self.annotation.get('text').__str__()
      else:
        anno = ''
      from ..core.CCP4ProjectsManager import PROJECTSMANAGER
      PROJECTSMANAGER().db().updateImportFile(importId=self.importId,key='annotation',value=anno)
    else:
      anno = self.annotation.get('text').__str__()
      if len(anno)>0:
        self.model.__dict__['sourceFileAnnotation'] = anno


class CImportInfoDialog(QtWidgets.QDialog):
  def __init__(self,parent,model=None,importId=None,importFileName=None,label=None,sourceFileAnnotation=''):
    # Expect input model to be set to a CDataFile if this is called for a freshly imported file
    # for a job that has not yet been run (so no files and importFiles yet saved to db)
    # Expect database importId and importFileName if eiting the import annotation after the
    # job has been run
    #print 'CImportInfoDialog',model,importFileName
    QtWidgets.QDialog.__init__(self,parent)
    self.setModal(True)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.showAutoWidget = QtWidgets.QCheckBox('Show this window automatically when importing a file (can be changed via Preferences)',self)
    self.showAutoWidget.toggled.connect(self.handleShowAuto)
    self.layout().addWidget(self.showAutoWidget)
    #self.layout().addWidget(QtWidgets.QLabel('This window accessible from the file icon menu',self))
    from ..core.CCP4Preferences import PREFERENCES
    if PREFERENCES().AUTO_INFO_ON_FILE_IMPORT:
      self.showAutoWidget.setCheckState(QtCore.Qt.Checked)
    else:
      self.showAutoWidget.setCheckState(QtCore.Qt.Unchecked)
    self.infoWidget = CImportInfo(self,model=model,importId=importId,importFileName=importFileName,label=label,
                                  sourceFileAnnotation=sourceFileAnnotation)
    self.layout().addWidget(self.infoWidget)
    butBox = QtWidgets.QDialogButtonBox(self)
    self.layout().addWidget(butBox)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Ok)
    but.setAutoDefault(0)
    but.clicked.connect(self.saveInfo)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.setAutoDefault(0)
    but.clicked.connect(self.close)

  def setProjectId(self,projectId):
    self.infoWidget.exportedFilesView.editor.setProjectId(projectId)

  @QtCore.Slot(str)
  def saveInfo(self):
    self.infoWidget.save()
    self.close()

  @QtCore.Slot(bool)
  def handleShowAuto(self,auto):
    from ..core.CCP4Preferences import PREFERENCES
    PREFERENCES().AUTO_INFO_ON_FILE_IMPORT = auto
    PREFERENCES().save()

class CManageImportFiles(QtWidgets.QDialog):
  def __init__(self,parent,projectId):
    QtWidgets.QDialog.__init__(self,parent)
    self.setModal(False)
    self.projectId = projectId
    self.setLayout(QtWidgets.QHBoxLayout())
    self.fileList = CImportFileList(self)   
    self.layout().addWidget(self.fileList )
    self.fileList.load(self.projectId)
    butLayout = QtWidgets.QVBoxLayout()
    for label,connect in [['Edit info',self.editInfo],['Delete',self.deleteFile]]:
      but = QtWidgets.QPushButton(label,self)
      but.clicked.connect(connect)
      butLayout.addWidget(but)
    self.layout().addLayout(butLayout)
    self.importInfoGui = None
    self.deleteFilesGui = None
    PROJECTSMANAGER().db().jobDeleted.connect(self.fileList.handleJobDeleted)

  @QtCore.Slot()
  def editInfo(self):
    rv = self.fileList.currentSelection()
    if rv is None: return
    importId,fileId,jobId,fileName = rv
    if self.importInfoGui is not None:
      self.importInfoGui.hide()
      self.importInfoGui.deleteLater()
    self.importInfoGui = CImportInfoDialog(self,importFileName=fileName,importId=importId)
    self.importInfoGui.show()

  @QtCore.Slot()
  def deleteFile(self):
    rv = self.fileList.currentSelection()
    if rv is None: return
    importId,fileId,jobId,fileName = rv
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    jobTree = PROJECTSMANAGER().db().getFollowOnJobs(jobId=jobId)
    delJobId,importFiles,followOnJobs = jobTree

    from . import CCP4ProjectViewer
    self.deleteJobGui = CCP4ProjectViewer.CDeleteJobGui(self,projectId=self.projectId,jobIdList=[jobId],jobTreeList=[jobTree], 
                                           label='Delete jobs that use imported file:'+fileName,deleteImportFiles=True)
    self.deleteJobGui.show()

class CImportFileList(QtWidgets.QTreeWidget):
  def __init__(self,parent):
    QtWidgets.QTreeWidget.__init__(self,parent)
    self.setColumnCount(5)
    self.setHeaderLabels(['Source file','Saved as','Date','Job number','Task name'])
    self.setColumnWidth(0,200)
    self.setColumnWidth(1,200)
    self.setColumnWidth(2,60)
    self.setColumnWidth(3,60)
    self.setColumnWidth(4,200)
    self.setMinimumWidth(460)
    self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    self.projectId = None
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    PROJECTSMANAGER().db().importFileDeleted.connect(self.handleImportFileDeleted)

  def load(self,projectId):
    # This returns list of: ImportId,JobId,FileID,FileTypeId,Filename,Annotation, JobNumber TaskName SourceFileName CreateTime
    #                          0        1       2        3         4           5        6        7         8            9
    self.projectId = projectId
    self.clear()
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    fileInfoList = PROJECTSMANAGER().db().getProjectImportFiles(projectId=projectId,ifJobInfo=True)
    for fileInfo in fileInfoList:
      from ..core.CCP4TaskManager import TASKMANAGER
      taskTitle = TASKMANAGER().getTitle(fileInfo[7])
      qList = []
      for item in [fileInfo[8],fileInfo[4],fileInfo[9],fileInfo[6],taskTitle]: qList.append(str(item))
      item = QtWidgets.QTreeWidgetItem(qList)
      item.setData(0,QtCore.Qt.UserRole,fileInfo[0])
      item.setData(1,QtCore.Qt.UserRole,fileInfo[1])
      item.setData(2,QtCore.Qt.UserRole,fileInfo[2])
      treeId=self.addTopLevelItem(item)

  @QtCore.Slot(dict)
  def handleImportFileDeleted(self,args):
    #print 'CImportFileList.handleImportFileDeleted',args
    if self.model().rowCount()==0: return
    modInxList = self.model().match(self.model().index(0,0),QtCore.Qt.UserRole,args['importId'],1)
    if len(modInxList)==0:
      print('CImportFileList.handleImportFileDeleted no match to importId',args)
    else:
      self.model().removeRow(modInxList[0].row())

  @QtCore.Slot(dict)
  def handleJobDeleted(self,args):
    #print 'CImportFileList.handleJobDeleted',args
    if args['projectId'] == self.projectId: self.load(self.projectId)

  def currentSelection(self):
    indices = self.selectionModel().selectedRows()
    #print 'CImportFileList.currentSelection',indices,len(indices)
    if len(indices) == 0: return None
    #print 'CImportFileList.currentSelection row',indices[0].row()
    importId = CCP4Data.varToUUID(indices[0].data(QtCore.Qt.UserRole))
    jobId = CCP4Data.varToUUID(indices[0].sibling(indices[0].row(),1).data(QtCore.Qt.UserRole))
    fileId = CCP4Data.varToUUID(indices[0].sibling(indices[0].row(),2).data(QtCore.Qt.UserRole))
    fileName = indices[0].data().__str__()
    #print 'CImportFileList.currentSelection', fileId,ok
    return importId,fileId,jobId,fileName
