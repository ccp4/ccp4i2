from __future__ import print_function


"""
     CCP4ProjectWidget.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.highlightL 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton April 2011 - list database
"""

##@package CCP4ProjectWidget View a project
                            
from PySide2 import QtGui, QtWidgets,QtCore,QtSvg
from core.CCP4Modules import WEBBROWSER,PROJECTSMANAGER,MIMETYPESHANDLER,QTAPPLICATION,LAUNCHER,PREFERENCES
from core.CCP4TaskManager import TASKMANAGER
from core.CCP4ErrorHandling import *
from dbapi import CCP4DbApi
from qtgui import CCP4StyleSheet
from core import CCP4Utils, CCP4File
import os,sys, time, datetime
import functools

_PROJECTMODEL = {}
_BROWSERMODE = 0
DRAGICONSIZE=16
_JOBICON = {}

def PROJECTMODEL(projectId):
    if projectId not in _PROJECTMODEL:
        #print QTAPPLICATION()
        _PROJECTMODEL[projectId] = CProjectModel(parent=QTAPPLICATION(),projectId=projectId)
    return _PROJECTMODEL[projectId]

def loadSvg(fileName,size=24):
    svg = QtSvg.QSvgRenderer()
    svg.load(fileName)
    pixmap = QtGui.QPixmap(size,size)
    pixmap.fill(QtGui.QColor(0,0,0,0))
    painter = QtGui.QPainter(pixmap)
    svg.render(painter)
    painter.end()
    return pixmap

def jobIcon(style='job'):
    if style not in _JOBICON:
      fileName = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons',style))
      if os.path.exists(fileName+'.png'):
        _JOBICON[style] = QtGui.QIcon(QtGui.QPixmap(fileName+'.png'))
        #print 'jobIcon',fileName+'.png'
      elif os.path.exists(fileName+'.svg'):
        _JOBICON[style] = QtGui.QIcon(loadSvg(fileName+'.svg'))
        #print 'jobIcon',fileName+'.svg'
      elif os.path.exists(fileName+'.gif'):
        movie = QtGui.QMovie(fileName+'.gif');
        label = QtWidgets.QLabel()
        label.setMovie(movie)
        movie.start()
        _JOBICON[style] = label
      else:
        fileName = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','evaluations',style))
        if os.path.exists(fileName+'.png'):
          _JOBICON[style] = QtGui.QIcon(QtGui.QPixmap(fileName+'.png'))
        elif os.path.exists(fileName+'.svg'):
          _JOBICON[style] = QtGui.QIcon(loadSvg(fileName+'.svg'))
    #print 'CCP4ProjectModel.jobIcon',JOBICON
    return  _JOBICON[style]


class CTreeItem:
  def __init__(self,parent=None):
    self._parent = parent
    self.name = ''  
  def getName(self,jobNumber=False):
    if isinstance(self.name,str):
      return self.name
    else:
      return self.name.__str__()
  def parent(self):
    return self._parent
  def setParent(self,parent):
    self._parent = parent
  def columnCount(self):
    return CProjectModel.NCOLUMNS
  def child(self,row):
    return None
  def childCount(self):
    return 0
  def canFetchMore(self):
    return False

  def isJob(self):
    return False
  def isFile(self):
    return False
  def isProject(self):
    return False
  def isFolder(self):
    return False
  def root(self):
   p = self
   while not isinstance(p.parent(),QtCore.QAbstractItemModel):
     p = p.parent()
   return p

  def model(self):
   p = self
   while not isinstance(p,QtCore.QAbstractItemModel):
     p = p.parent()
   return p
    
     
class CTreeItemFolder(CTreeItem):
  ERROR_CODES = { 101 : { 'description' : 'Failed deleting project from folder' } }
  def __init__(self,parent=None,info={}):
    self._parent = parent
    self.name = info.get('name')
    self.childProjects = []
    self.childFolders = []

  def isFolder(self):
    return True

  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.name
    elif role == QtCore.Qt.DecorationRole:
      if column == 0:
        return jobIcon('fileopen') 
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None
  
  def child(self,row):
    if row<len(self.childFolders):
      return self.childFolders[row]
    else:
      row = row-len(self.childFolders)
      if row<len(self.childProjects):
        return self.childProjects[row]
    return None

  def columnCount(self):
    return 1

  def rowCount(self):
    return self.childCount()

  def childCount(self,mode='all'):
    if mode == 'folders':
      return len(self.childFolders)
    elif mode == 'projects':
      return  len(self.childProjects)
    else:
      return len(self.childProjects) + len(self.childFolders)

  def appendChildProject(self,child):
    self.childProjects.append(child)

  def appendChildFolder(self,child):
    self.childFolders.append(child)

  def removeChildProject(self,projectName):
    for ii in range(len(self.childProjects)):
      if self.childProjects[ii].refProject == projectName:
        del self.childProjects[ii]
        return CErrorReport()
   
    return CErrorReport(self.__class__,101,projectName)

  def removeChildFolder(self,folderName):
    for ii in range(len(self.childFolders)):
      if self.childFolders[ii].name == folderName:
        del self.childFolders[ii]
        return CErrorReport()
   
    return CErrorReport(self.__class__,101,folderName)

  def row(self):
    try:
      return self._parent.childFolders.index(self)
    except:
      print('CTreeItemFolder.row() failed')
      return -1

class CTreeItemProject(CTreeItem):
  # The project tree item currently only supports having child files
  # for use in the database browser
  def __init__(self,parent=None,infoList=[]):
    CTreeItem.__init__(self,parent=parent)
    self.projectId =  infoList[0]
    self.projectName = infoList[1]
    self.childProjects = []
    self.childFiles = []
    self.childJobs = []
    self.childFilesLoaded = False

  def isProject(self):
    return True

  def addDummyFile(self):
    # The files in a prject are loaded when user opens the project folder
    # but the project will only appear as a folder if we give it a dummy childfile
    self.appendChildFile(CTreeItemFile(self,[' ',' ',None,0,' ',' ','']))

  def clearFiles(self):
     self.childFiles = []

  def getProjectId(self):
    return self.projectId.__str__()

  def getProjectName(self):
    return self.projectName.__str__()

  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.projectName
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier

  def canFetchMore(self):
    return (not self.childFilesLoaded)
      
  def appendChildProject(self,child):
    self.childProjects.append(child)
    
  def appendChildFile(self,child):
    self.childFiles.append(child)
    
  def appendChildJob(self,child):
    self.childJobs.append(child)

  def child(self,row):
    if row<len(self.childProjects):
      return self.childProjects[row]
    else:
      row = row-len(self.childProjects)
      if row<len(self.childJobs):
        return self.childJobs[row]
      else:
        row = row-len(self.childJobs)
        if row<len(self.childFiles):
          return self.childFiles[row]     
    return None

  def columnCount(self):
    return 1

  def rowCount(self):
    return self.childCount()

  def childCount(self,mode='all'):
    if mode == 'files':
      return len(self.childFiles)
    elif mode == 'projects':
      return  len(self.childProjects)
    elif mode == 'jobs':
      return  len(self.childJobs)
    else:
      #print 'CTreeItemProject.childCount',len(self.childProjects) + len(self.childFiles)
      return len(self.childProjects) + len(self.childFiles) + len(self.childJobs)

  '''
  def childFileIndex(self,fileNode=None):
    if fileNode is not None and self.childFiles.count(fileNode):
      return self.childFiles.index(fileNode) 
    else:
      return -1
  '''

  def data(self,column,role):
    #print 'CTreeItemJob.data',column,role
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.projectName
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.projectId
      '''
      elif role == QtCore.Qt.BackgroundRole:
        if self.highlight:
          return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.HIGHLIGHTCOLOUR))
        
        elif CProjectModel.CONTRAST_MODE == 1 and self.top:
          return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.LOWLIGHTCOLOUR))
                             
      elif role == QtCore.Qt.FontRole and CProjectModel.CONTRAST_MODE == 2:
        if self.top:
          return CTreeItemJob.boldFont
        else:
          return CTreeItemJob.italicFont
      '''

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def row(self):
    try:
      return self._parent.childProjects.index(self)
    except:
      print('CTreeItemProject.row() failed')
      return -1

  
class CTreeItemFile(CTreeItem):
  def __init__(self,parent=None,infoList=[],displayColumn=0,maxChar=40,displayJobNumber=False):
    CTreeItem.__init__(self,parent=parent)
    #jid,fid,impid,ftype,fname,annotation = infoList
    self.displayColumn = displayColumn
    self.fileId = infoList[1]
    self.fileName = infoList[4]
    self.fileType = infoList[3]
    self.ifImport = infoList[2] is not None
    if len(infoList)>6:
      self.jobNumber = infoList[6]
    else:
      self.jobNumber = ''
    self.identifier = self.fileId
    try:
      # Beware this could be broken if using older version of code
      mimeType = CCP4DbApi.FILETYPES_TEXT[self.fileType]
      if self.ifImport:
        self.icon = MIMETYPESHANDLER().icon(mimeType,modifier='import')
      else:
        self.icon = MIMETYPESHANDLER().icon(mimeType)
    except:
      pass
    self.setName(infoList[5],maxChar=maxChar,displayJobNumber=displayJobNumber)

  def yAdjust(self):
    return 0

  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == self.displayColumn:
        return "<span class=\"fileItem\">"+str(self.name)+"</span>"
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier
    elif role == QtCore.Qt.UserRole + 1:
        return 1
    elif role == QtCore.Qt.UserRole + 2:
        return self.name
    elif role == QtCore.Qt.DecorationRole:
      if column == self.displayColumn:
        return self.icon
    elif role == QtCore.Qt.EditRole:
      if column == 0: return True
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def updateName(self,annotation,mimeType=None,maxChar=40,displayJobNumber=False):
    if displayJobNumber:
      text = self.jobNumber + ' '
    else:
      text = ''
    if annotation is not None and len(annotation)>0:
      text = text + annotation
    else:
      if mimeType is None:
        try:
          mimeType = CCP4DbApi.FILETYPES_TEXT[self.fileType]
        except:
          pass
      if mimeType is not None:
        text = text +  MIMETYPESHANDLER().getMimeTypeInfo(mimeType,'description')
      else:
        text = text + self.fileName        
    self.name = text[0:maxChar]

    try:
        PROJECTSMANAGER().db().updateFile(self.fileId,key='annotation',value=annotation)
    except:
        print("Failed to update fileId")
    

  def setName(self,annotation,mimeType=None,maxChar=40,displayJobNumber=False):
    if displayJobNumber:
      text = self.jobNumber + ' '
    else:
      text = ''
    if annotation is not None and len(annotation)>0:
      text = text + annotation
    else:
      if mimeType is None:
        try:
          mimeType = CCP4DbApi.FILETYPES_TEXT[self.fileType]
        except:
          pass
      if mimeType is not None:
        text = text +  MIMETYPESHANDLER().getMimeTypeInfo(mimeType,'description')
      else:
        text = text + self.fileName        
    self.name = text[0:maxChar]

  
  def row(self):
    ii = self._parent.len(self._parent.childInFiles) - 1
    for item in self._parent.childOutFiles:
      ii += 1
      if item == self: return ii
    return -1

  def getJobId(self):
    return self._parent.jobId

  def getFileId(self):
    return self.fileId

  def getFilename(self):
    return self.fileName
 
  def isJob(self):
    return False

  def isFile(self):
    return True

  def mimeData(self):
    from lxml import etree
    urlList = []
    mimeType = CCP4DbApi.FILETYPES_CLASS[self.fileType]
    root,err = PROJECTSMANAGER().db().getFileEtree(fileId=self.fileId)
    info = PROJECTSMANAGER().db().getFileInfo(fileId=self.fileId,mode='projectid')
    #print 'CTreeItemFile.mimeData projectId',projectId
    relPath = root.find('relPath')
    #print 'CTreeItemFile.mimeData projectId',info,PROJECTSMANAGER().getProjectDirectory(projectId=info['projectid']),relPath.text
    relPath.text = os.path.normpath(os.path.join(PROJECTSMANAGER().getProjectDirectory(projectId=info['projectid']),str(relPath.text)))
    dragText = etree.tostring(root,pretty_print=False)
    urlList = [QtCore.QUrl()]
    urlList[0].setPath( PROJECTSMANAGER().db().getFullPath(fileId=self.fileId ) )
    urlList[0].setScheme("file")
    #print 'CProjectModel.mimeData',mimeType,dragText
    encodedData = QtCore.QByteArray()
    encodedData.append(dragText)
    mimeData = QtCore.QMimeData()
    mimeData.setData(mimeType,encodedData)
    if len(urlList)>0:
      mimeData.setUrls(urlList)
    return mimeData

def formatDate(intTime):
#FIXME PYQT - or maybe None? This used to return QVariant.
  if intTime is None: ""
  try:
    if PREFERENCES().JOB_LIST_DATE_TIME:
      return time.strftime(CTreeItemJob.DATE_FORMAT+" "+CTreeItemJob.TIME_FORMAT,time.localtime(intTime))
    date = time.strftime(CTreeItemJob.DATE_FORMAT,time.localtime(intTime))
    if date == CTreeItemJob.TODAY:
      date = time.strftime(CTreeItemJob.TIME_FORMAT,time.localtime(intTime))
    elif date.split(' ')[-1] == CTreeItemJob.THISYEAR:
      date=date.rsplit(' ',1)[0]
    else:
      date=date.split(' ',1)[1]
    return  date
  except:
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None
           

class JobListHTMLDelegate(QtWidgets.QStyledItemDelegate):

    editingFinishedSignal = QtCore.Signal()

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate. __init__(self,parent)
        self.editorWidget = None
        qticonsDir = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons')
        self.donePix = QtGui.QPixmap(os.path.join(qticonsDir,"job.png") )

    def createEditor(self,parent,option,modelIndex):
        #if self.parent().model().data(modelIndex,QtCore.Qt.DisplayRole):
        window = QtWidgets.QDialog(parent)
        window.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.Dialog)
        layout = QtWidgets.QVBoxLayout()
        layout.setMargin(CProjectWidget.MARGIN)
        layout.setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
        window.setLayout(layout)
        editor = CJobEditLineEdit(window)
        editor.setMinimumWidth(350)
        editor.setObjectName('editor')
        window.layout().addWidget(editor)
        self.editorWidget = window
        self.editorWidget.setFocus()
        editor.losingFocus.connect(self.closeEditor)
        return window

    @QtCore.Slot()
    def closeEditor(self):
        if self.editorWidget is not None:
          self.editingFinishedSignal.emit()
          self.editorWidget.deleteLater()
          self.editorWidget = None
    
    def setEditorData(self,editor,modelIndex):
        data = self.parent().model().sourceModel().nodeFromIndex(modelIndex.model().mapToSource(modelIndex)).getName(jobNumber=False)
        if data is not None:
          editor.findChild(QtWidgets.QLineEdit,'editor').setText(data)
  
    def setModelData(self,editor,model,modelIndex):
        text = str(editor.findChild(QtWidgets.QLineEdit,'editor').text())
        node = self.parent().model().sourceModel().nodeFromIndex(modelIndex.model().mapToSource(modelIndex))
        if node.isJob():
          try:
            PROJECTSMANAGER().db().updateJob(node.getJobId(),key='jobtitle',value=text)
          except:
            pass
        elif node.isFile():
          try:
            PROJECTSMANAGER().db().updateFile(node.getFileId(),key='annotation',value=text)
          except:
            pass
 
    def updateEditorGeometry(self,editor,option,modelIndex):
        #Put the edit line in the right place
        editor.show()
        view = editor.parent().parent()
        # mapToGlobal() does not allow for the header - need to add header height
        # and want to not overlap the icon (and the job number for jobs)
        node = self.parent().model().sourceModel().nodeFromIndex(self.parent().model().mapToSource(modelIndex))
        rightShift = 25 + node.isJob()*25
        pos=view.mapToGlobal(view.visualRect(modelIndex).topLeft())
        editor.move(pos+QtCore.QPoint(rightShift,20))

    def paint(self, painter, option, index):
        options = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(options,index)

        style = QtWidgets.QApplication.style() if options.widget is None else options.widget.style()

        doc = QtGui.QTextDocument()

        textOption = doc.defaultTextOption()
        textOption.setWrapMode(QtGui.QTextOption.NoWrap)
        doc.setDefaultTextOption(textOption)

        if options.font.pixelSize()> -1:
            pixSize = str(options.font.pixelSize())
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
            font-size: """+pixSize+"""px;
        }
        .jobStatusPending {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .jobStatusRunning {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown  {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
            font-size: """+pixSize+"""px;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
            font-size: """+pixSize+"""px;
        }
        .jobKpi {
            font-size: """+pixSize+"""px;
        }
        .subJob {
            font-size: """+pixSize+"""px;
        }
        </style>
        </head><body>
            """
        elif options.font.pointSize()> -1:
            ptSize = str(options.font.pointSize())
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusPending {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusRunning {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
            font-size: """+ptSize+"""pt;
        }
        .jobKpi {
            font-size: """+ptSize+"""pt;
        }
        .subJob {
            font-size: """+ptSize+"""pt;
        }
        </style>
        </head><body>
            """
        else:
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
        }
        .jobStatusPending {
            font-style: italic;
        }
        .jobStatusRunning {
            font-style: italic;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown {
            font-style: italic;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
        }
        .jobKpi {
        }
        </style>
        </head><body>
            """

        html += options.text + "</body>"

        if options.font.pixelSize()> -1:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX",str(options.font.pixelSize()))
        elif options.font.pointSize()> -1:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX",str(options.font.pointSize()))
        else:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX","13")

        if options.font.pixelSize()> -1:
            theIconSize = max(options.font.pixelSize(),14)
            options.decorationSize = QtCore.QSize(theIconSize,theIconSize)

        doc.setHtml(html)

        doc.setTextWidth(options.rect.width());
        options.text = ""
        
        style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, options, painter);

        color1 = QtGui.QColor(220,220,255)
        color2 = QtGui.QColor(240,240,255)

        ctx = QtGui.QAbstractTextDocumentLayout.PaintContext()

        # Highlighting text if item is selected
        if options.state & QtWidgets.QStyle.State_Selected:
            ctx.palette.setColor(QtGui.QPalette.Text, options.palette.color(QtGui.QPalette.Active, QtGui.QPalette.HighlightedText));
            color1 = QtGui.QColor(0,0,255)
            color2 = QtGui.QColor(0,0,255)

        textRect = style.subElementRect(QtWidgets.QStyle.SE_ItemViewItemText, options)

        yAdjust = index.model().mapToSource(index).internalPointer().yAdjust()
        textRect.adjust(0, yAdjust, 0, 0)
 
        textRect.setTop(textRect.top())
        painter.save()
        
        painter.translate(textRect.topLeft())
        painter.setClipRect(textRect.translated(-textRect.topLeft()))

        doc.documentLayout().draw(painter, ctx)

        painter.restore()

    def sizeHint(self, option, index):
        options = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(options,index)

        doc = QtGui.QTextDocument()

        textOption = doc.defaultTextOption()
        textOption.setWrapMode(QtGui.QTextOption.NoWrap)
        doc.setDefaultTextOption(textOption)

        if options.font.pixelSize()> -1:
            pixSize = str(options.font.pixelSize())
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
            font-size: """+pixSize+"""px;
        }
        .jobStatusPending {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .jobStatusRunning {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown {
            font-style: italic;
            font-size: """+pixSize+"""px;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
            font-size: """+pixSize+"""px;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
            font-weight: bold;
            font-size: """+pixSize+"""px;
        }
        .jobKpi {
            font-size: """+pixSize+"""px;
        }
        .subJob {
            font-size: """+pixSize+"""px;
        }
        </style>
        </head><body>
            """
        elif options.font.pointSize()> -1:
            ptSize = str(options.font.pointSize())
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusPending {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusRunning {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown {
            font-style: italic;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
            font-size: """+ptSize+"""pt;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
            font-weight: bold;
            font-size: """+ptSize+"""pt;
        }
        .jobKpi {
            font-size: """+ptSize+"""pt;
        }
        .subJob {
            font-size: """+ptSize+"""pt;
        }
        </style>
        </head><body>
            """
        else:
            html = """<html><head>
        <style>
        .jobHeader {
            font-weight: bold;
        }
        .jobStatusPending {
            font-style: italic;
        }
        .jobStatusRunning {
            font-style: italic;
        }
        .fileItem, .jobStatusFileHolder, .jobStatusUnsatisfactory, .jobStatusToDelete, .jobStatusUnknown {
            font-style: italic;
        }
        .jobStatusFailed {
            font-style: italic;
            color: red;
        }
        .jobStatusSuccess, .jobStatusInterrupted {
            font-weight: bold;
        }
        .jobKpi {
        }
        </style>
        </head><body>
            """

        html += options.text + "</body>"

        if options.font.pixelSize()> -1:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX",str(options.font.pixelSize()))
        elif options.font.pointSize()> -1:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX",str(options.font.pointSize()))
        else:
            html = html.replace("XXXX_IMG_HEIGHT_XXXX","13")

        if options.font.pixelSize()> -1:
            theIconSize = max(options.font.pixelSize(),14)
            options.decorationSize = QtCore.QSize(theIconSize,theIconSize)

        doc.setHtml(html)

        doc.setTextWidth(options.rect.width());

        itemType = index.data(QtCore.Qt.UserRole + 1)

        mult = 1.0
        if options.font.pixelSize()> -1:
            mult = options.font.pixelSize()/12.
        if options.font.pointSize()> -1:
            mult = options.font.pointSize()/12.

        if itemType == 1:
            return QtCore.QSize(doc.idealWidth(), doc.size().height())
        elif itemType == 2:
            return QtCore.QSize(doc.idealWidth(), doc.size().height())

        return QtCore.QSize(doc.idealWidth(), doc.size().height())
 
class CTreeItemJob(CTreeItem):

  statusList = ['Finished','Interrupted','Failed','Running','File holder','To delete','Unsatisfactory','Running remotely']
  statusIconList = ['job','job_interrupted','sad','running','fileopen','list_delete','unsatisfactory','running']
  boldFont = None
  italicFont = None
  TIME_FORMAT = '%H:%M'
  DATE_FORMAT = '%a %d %b %y'
  TODAY = None
  THISYEAR = None
  qticonsDir = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons')
  DONE_PIX = os.path.join(qticonsDir,"green-tick.png")
  PENDING_PIX = os.path.join(qticonsDir,"undone.png")
  BEST_PIX = os.path.join(qticonsDir,"evaluations","Best.png")
  GOOD_PIX = os.path.join(qticonsDir,"evaluations","Good.png")
  REJEJECTED_PIX = os.path.join(qticonsDir,"evaluations","Rejected.png")
  FAILED_PIX = os.path.join(qticonsDir,"red-cross.png")
  BLANK_PIX = os.path.join(qticonsDir,"blank.png")
  RUNNING_PIX = os.path.join(qticonsDir,"running.png")
  RUNNING_DARK_PIX = os.path.join(qticonsDir,"running_dark.png")
  FOLDER_PIX = os.path.join(qticonsDir,"file_manager2.png")
  TO_DELETE_PIX = os.path.join(qticonsDir,"dustbin.png")
  PAUSE_PIX = os.path.join(qticonsDir,"pause.png")
  UNKNOWN_PIX = os.path.join(qticonsDir,"question.png")
  IMG_HEIGHT_STR = "XXXX_IMG_HEIGHT_XXXX"
  
  def __init__(self,parent=None,info={}):
    CTreeItem.__init__(self,parent=parent)
    #import traceback
    #print 'CTreeItemJob',info.get('jobnumber',None),info.get('taskname',None)
    #traceback.print_stack()
    #print '\n\n'
    self.childJobs = []
    self.childInFiles = []
    self.childOutFiles = []
    self.highlight = False
    self.jobId = None
    self.taskName =  info.get('taskname')
    if len(info)==0: return
    self.jobId = info.get('jobid')
    self.jobNumber = info.get('jobnumber')
    self.identifier = self.jobId
    self.performance = None
    if info.get('performance',None) is not None:
      performanceClass = TASKMANAGER().getPerformanceClass(self.taskName)
      if performanceClass is not None:
        self.performance = performanceClass()
        self.performance.set(info['performance'])
    #print 'CTreeItemJob.__init__ jobtitle',info.get('jobtitle',None), TASKMANAGER().getShortTitle( info.get('taskname'),substitute=False)
    self.setName(info['jobtitle'])
    if info['parentjobid'] is not None:
#FIXME PYQT - or maybe None? These used to set QVariant.
      self.evaluation =  ""
      self.followFrom = ""
      self.top=False
    else:
      self.evaluation = jobIcon(info.get('evaluation'))
      self.followFrom = jobIcon('greendot')
      self.top=True
    self.statusStr = info.get('status')
    self.evaluationStr = info.get('evaluation')
    if info.get('status') in CTreeItemJob.statusList:
      self.status = jobIcon(CTreeItemJob.statusIconList[CTreeItemJob.statusList.index(info.get('status'))])
      status = info['status']+ ' '
    else:
      self.status = jobIcon('Unknown')
      status = 'Job pending '
    toolTip =TASKMANAGER().getTitle( info.get('taskname'))
    if info.get('jobtitle',None) is not None and info['jobtitle'] != toolTip:
      toolTip =  "Job "+str(self.jobNumber)+' '+toolTip + '\n' + info['jobtitle'] + '\n'
    else:
      toolTip =  "Job "+str(self.jobNumber)+' '+toolTip + '\n'
    if self.performance is not None: toolTip = toolTip + self.performance.__str__() + '\n'

    self.dateTime =  formatDate(info['finishtime'])
    
    self.toolTip = toolTip + status + '\nRight mouse click for options'
    self.colour = None

    #print 'CTreeItemJob.__init__',self.name.__str__()
    if CTreeItemJob.boldFont is None: 
      CTreeItemJob.boldFont = QtGui.QFont()
      CTreeItemJob.boldFont.setItalic(True)
      CTreeItemJob.boldFont.setBold(True)
      CTreeItemJob.italicFont = QtGui.QFont()
      CTreeItemJob.italicFont.setItalic(True)

  def updateName(self,jobTitle=None):
    if jobTitle is not None and len(jobTitle)>0:
      self.name = self.jobNumber +' '+ jobTitle
    else:
      self.name = self.jobNumber +' '+ TASKMANAGER().getShortTitle( self.taskName )
    try:
      PROJECTSMANAGER().db().updateJob(self.jobId,key='jobtitle',value=jobTitle)
      from dbapi import CCP4DbUtils
      CCP4DbUtils.makeJobBackup(jobId=self.jobId)
    except:
      print('ERROR in editing job name')
      
  def setName(self,jobTitle=None):
    #import traceback
    #traceback.print_stack(limit=5)
    if jobTitle is not None and len(jobTitle)>0:
      self.name = self.jobNumber +' '+ jobTitle
    else:
      self.name = self.jobNumber +' '+ TASKMANAGER().getShortTitle( self.taskName )
      
  def update(self,info={}):
    #print 'CTreeItemJob.update',info
    if info['parentjobid'] is None:
      self.evaluation = jobIcon(info.get('evaluation'))
    if info.get('status') in CTreeItemJob.statusList:
      self.status = jobIcon(CTreeItemJob.statusIconList[CTreeItemJob.statusList.index(info.get('status'))])
    else:
      #print 'CTreeItemJob.update status unknown',info.get('status')
      self.status = jobIcon('Unknown')
    self.setName(info.get('jobtitle',None))

  def set(self,key,value):
    #print 'CTreeItemJob.set',key,value
    key = key.lower()
    if key == 'evaluation':
      value = CCP4DbApi.JOB_EVALUATION_TEXT[value]
      self.evaluation = jobIcon(value)
      self.evaluationStr = value
    elif key == 'status':
      if value in CTreeItemJob.statusList:
        #print 'CTreeItemJob.set status',value,CTreeItemJob.statusIconList[CTreeItemJob.statusList.index(value)]
        self.statusStr = value
        self.status = jobIcon(CTreeItemJob.statusIconList[CTreeItemJob.statusList.index(value)])
      else:
        self.status = jobIcon('Unknown')
    elif key == 'jobtitle':
      self.setName(value)
    elif key == 'highlight':
      self.highlight = value
    elif key == 'performance':
      cls =  TASKMANAGER().getPerformanceClass(self.taskName)
      if cls is not None:
        self.performance = cls()
        self.performance.set(value)
    elif key == 'finishtime':
      self.dateTime =  formatDate(value)

  def getName(self,jobNumber=True):
    if isinstance(self.name,str):
      name = self.name
    else:
      name = self.name.__str__()
    if jobNumber:
      return name
    else:
      return name[len(self.jobNumber)+1:]
    

  def appendChildJob(self,child):
    self.childJobs.append(child)

  def removeChildJob(self,row):
    if row>=0 and row < len(self.childJobs):
      del self.childJobs[row]
      
  def removeChildFile(self,row=None,fileNode=None):
    if row is not None:
      if row < len(self.childInFiles):
        del self.childInFiles[row]
      else:
        row = row - len(self.childInFiles)
        if row < len(self.childOutFiles):
          del self.childOutFiles[row]
        else:
          return False
    elif fileNode is not None:
      if self.childInFiles.count(fileNode):
        ii = self.childInFiles.index(fileNode)
        del self.childInFiles[ii]
      elif self.childOutFiles.count(fileNode):
        ii = self.childOutFiles.index(fileNode)
        del self.childOutFiles[ii]
      else:
        return False
    else:
      return False

  def removeChildOutputFiles(self):
    for indx in range(len(self.childOutFiles)-1,-1,-1):
      del self.childOutFiles[indx]

  def removeChildJobs(self):
    for indx in range(len(self.childJobs)-1,-1,-1):
      del self.childJobs[indx]
       
  def appendChildInputFile(self,child):
    self.childInFiles.append(child)
    
  def appendChildOutputFile(self,child):
    self.childOutFiles.append(child)

  def insertChildJob(self,pos,child):
    self.childJobs.insert(pos,child)

  def child(self,row):
    if row<len(self.childInFiles):
      return self.childInFiles[row]
    else:
      row = row-len(self.childInFiles)
      if row<len(self.childOutFiles):
        return self.childOutFiles[row]
      else:
        row = row-len(self.childOutFiles)
        if row<len(self.childJobs):
          return self.childJobs[row]
        else:
          return None

  def childCount(self,mode='jobs'):
    if mode == 'infiles':
      return len(self.childInFiles)
    elif mode == 'outfiles':
      return  len(self.childInFiles) + len(self.childOutFiles)
    else:
      return len(self.childJobs) + len(self.childInFiles) + len(self.childOutFiles)

  def childFileIndex(self,fileNode=None):
    if fileNode is not None:
      if self.childInFiles.count(fileNode):
        return self.childInFiles.index(fileNode)
      elif self.childOutFiles.count(fileNode):
        return self.childOutFiles.index(fileNode) + len(self.childInFiles)
      else:
        return -1

  def yAdjust(self):
      return 0

  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if "." in str(self.name).split()[0]:
          if self.statusStr == "Failed":
              bigStr = "<span class=\"subJob\">" + str(self.name) + "</span>&nbsp;&nbsp;<span class=\"jobStatusFailed\">" + self.statusStr + "</span>"
          elif self.statusStr == "Pending":
              bigStr = "<span class=\"subJob\">" + str(self.name) + "</span>&nbsp;&nbsp;<span class=\"jobStatusPending\">" + self.statusStr + "</span>"
          else:
              if self.performance is not None:
                  bigStr = "<span class=\"subJob\">" + str(self.name) + "&nbsp;&nbsp;" + str(self.performance) + "</span>"
              else:
                  bigStr = "<span class=\"subJob\">" + str(self.name) + "</span>"
          return bigStr
      if column == 0:
        bigStr = '<table width="100%" cellspacing="0" cellpadding="0">'
        bigStr += "<tr><td class=\"jobHeader\">" + str(self.name) + "</tr></td>"
        bigStr += "<tr>"
        if self.evaluationStr in ["Best","Good","Rejected"] and self.statusStr == "Finished":
            if self.evaluationStr == "Best":
                bigStr += " <td class=\"jobStatusSuccess\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.BEST_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " (" + self.evaluationStr + ") " + str(self.dateTime) + "</span></td>"
            if self.evaluationStr == "Good":
                bigStr += " <td class=\"jobStatusSuccess\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.GOOD_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " (" + self.evaluationStr + ") " + str(self.dateTime) + "</span></td>"
            if self.evaluationStr == "Rejected":
                bigStr += " <td class=\"jobStatusSuccess\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.REJEJECTED_PIX+"\"/>&nbsp;&nbsp;<span class=\"jobStatusSuccess\">" + self.statusStr + " (" + self.evaluationStr + ") " + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "Failed":
            bigStr += " <td class=\"jobStatusFailed\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.FAILED_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "Pending" or self.statusStr == "Queued":
            bigStr += " <td class=\"jobStatusPending\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.PENDING_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + "</span></td>"
        elif self.statusStr == "Finished":
            bigStr += " <td class=\"jobStatusSuccess\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.DONE_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"

        elif self.statusStr == "Unknown":
            bigStr += " <td class=\"jobStatusUnknown\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.UNKNOWN_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "Interrupted":
            bigStr += " <td class=\"jobStatusInterrupted\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.PAUSE_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "File holder":
            bigStr += " <td class=\"jobStatusFileHolder\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.FOLDER_PIX+"\"/>&nbsp;&nbsp;<span>" + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "To delete":
            bigStr += " <td class=\"jobStatusToDelete\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.TO_DELETE_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
        elif self.statusStr == "Unsatisfactory":
                bigStr += " <td class=\"jobStatusUnsatisfactory\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.REJEJECTED_PIX+"\"/>&nbsp;&nbsp;<span class=\"jobStatusSuccess\">" + self.statusStr + " (" + self.evaluationStr + ") " + str(self.dateTime) + "</span></td>"

        else:
            label = QtWidgets.QLabel("Am I in the dark?")
            text_hsv_value = label.palette().color(QtGui.QPalette.WindowText).value()
            bg_hsv_value = label.palette().color(QtGui.QPalette.Background).value()
            isDarkMode = text_hsv_value > bg_hsv_value
            if isDarkMode:
                bigStr += " <td class=\"jobStatusRunning\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.RUNNING_DARK_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
            else:
                bigStr += " <td class=\"jobStatusRunning\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.RUNNING_PIX+"\"/>&nbsp;&nbsp;<span>" + self.statusStr + " " + str(self.dateTime) + "</span></td>"
        bigStr += "</tr>"
        if self.performance is not None and not self.statusStr == "File holder":
            bigStr += "<tr><td class=\"jobKpi\"><img height=\""+self.IMG_HEIGHT_STR+"\" src=\""+self.BLANK_PIX+"\"/>&nbsp;&nbsp;" + str(self.performance) + "</td></tr>"
        bigStr += "</table>"
        return bigStr
      elif column == 2 and self.performance is not None:
        return ""
        return self.performance.data(QtCore.Qt.DisplayRole)
      elif column == 3:
        return ""
        return self.dateTime
    if role == QtCore.Qt.EditRole:
       try:
           justName = self.name
           jns = self.jobNumber
           lenNum = len(jns) + 1
           return justName[lenNum:]
       except:
           print("Setting edit name failed")

    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier
    elif role == QtCore.Qt.UserRole + 1:
        if "." in str(self.name).split()[0]:
           return 1
        elif self.performance is not None:
           return 0
        else:
           return 2
    elif role == QtCore.Qt.UserRole + 2:
      return self.name
    elif role == QtCore.Qt.DecorationRole:
      if column == 0:
        return None
        if self.top:
          if self.evaluationStr in ["Best","Good","Rejected"] and self.statusStr == "Finished":
              return self.evaluation
          return self.status
    elif role == QtCore.Qt.BackgroundRole:
      #print 'CJobTreeItem.data', self.statusStr
      if self.statusStr in  ['Running','Running remotely']:
         return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.RUNCOLOUR))
      if self.highlight:
        return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.LOWLIGHTCOLOUR))
      
      elif CProjectModel.CONTRAST_MODE == 1 and self.top:
        return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.LOWLIGHTCOLOUR))

    elif role == QtCore.Qt.ForegroundRole:
      try:
        colour = self.root().highlightList.get(self.jobNumber,None)
      except:
        #print 'ERROR CTreeItemJob.data root',self.name,self.root()
        colour = None
      if colour is None:
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None
      else:
        return QtGui.QBrush(QtGui.QColor(colour))
      
    elif role == QtCore.Qt.FontRole and CProjectModel.CONTRAST_MODE == 2:
      if self.top:
        return CTreeItemJob.boldFont
      else:
        return CTreeItemJob.italicFont
    elif role == QtCore.Qt.ToolTipRole:
      if column == 0:
        return "Right mouse click for menu of options"
      else:
        return self.toolTip

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def row(self):
    ii = -1
    for item in self._parent.childJobs:
      ii += 1
      if item == self: return ii
    return -1

  def getTopJob(self):
    # Track up to the top (pipeline) job
    #print 'CTreeItemJob.getTopJob',self,self.jobId
    if not isinstance(self._parent,CTreeItemJob):
      return self
    if self._parent.jobId is None:
      # its the root
      return self
    else:
      return self._parent.getTopJob()

  def getJobId(self):
    return self.jobId

  def getTaskName(self):
    return self.taskName

  def getStatus(self):
    if self.jobId is not None:
      return PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId,mode='status')
    else:
      return None
      
  
  def isJob(self):
    return True

  def isTopJob(self):
    return isinstance(self._parent._parent,CProjectModel)

  def isFile(self):
    return False

  def mimeData(self):
    from lxml import etree
    urlList = []
    mimeType = 'FollowFromJob'
    root = etree.Element('jobId')
    root.text = str(self.jobId)
    dragText = etree.tostring(root,pretty_print=False)
    sceneFiles = PROJECTSMANAGER().getSceneFiles(jobId=self.jobId)
    if len(sceneFiles)>0:
      urlList = [QtCore.QUrl()]
      urlList[0].setPath( sceneFiles[0] )
      urlList[0].setScheme("file")
    #print 'CProjectModel.mimeData',mimeType,dragText
    encodedData = QtCore.QByteArray()
    encodedData.append(dragText)
    mimeData = QtCore.QMimeData()
    mimeData.setData(mimeType,encodedData)
    if len(urlList)>0:
      mimeData.setUrls(urlList)
    return mimeData


class CProjectModel(QtCore.QAbstractItemModel):

  redrawSignal = QtCore.Signal()

  CONTRAST_MODE = 2

  NCOLUMNS = 1
  COLUMNS = ['name']
  COLUMNHEADERS = { 'name':'Job/File' }
  # performance needs to be last in this list due to workings of CDbApi.jobInfoDict()
  JOBINFOITEMS = ['jobid','jobnumber','jobtitle','status','evaluation','taskname','parentjobid','finishtime','performance']
  
  def __init__(self, parent=None, projectId=None):
    QtCore.QAbstractItemModel.__init__(self,parent)
    self._projectId = projectId
    self.parents=[]

    self.rootItem = CTreeItemJob(self)
    self.rootItem.taskName = 'root'
    self.rootItem.highlightList = {}
    #print 'CProjectModel.__init__ rootItem',self.rootItem
    self.setupModelData()

    PROJECTSMANAGER().db().jobFinished.connect(self.updateFinishedJob)
    PROJECTSMANAGER().db().jobUpdated.connect(self.updateJob)
    PROJECTSMANAGER().db().jobCreated.connect(self.createJob)
    PROJECTSMANAGER().db().jobStarted.connect(self.createJob)
    PROJECTSMANAGER().db().jobDeleted.connect(self.deleteJob)
    PROJECTSMANAGER().db().fileUpdated.connect(self.updateFile)
    PROJECTSMANAGER().db().projectReset.connect(self.resetAll)
    PROJECTSMANAGER().db().followFromJobChanged.connect(self.updateFollowFrom)
    PROJECTSMANAGER().db().setJobToImportSignal.connect(self.setJobToImport)


  @QtCore.Slot(dict)
  def resetAll(self,args):
    #print 'CProjectModel.resetAll',args
    if not args['projectId'] ==  self._projectId: return
    #print 'CProjectModel.resetAll resetting'
    self.beginResetModel()
    self.rootItem = CTreeItemJob(self)
    self.rootItem.taskName = 'root'
    
    self.setupModelData()
    self.endResetModel()
    #print 'from resetAll'
    
  def setData(self, index, value, role):
    if index.isValid() and role == QtCore.Qt.EditRole:
      prev_value = self.getValue(index,role)
      item = index.internalPointer()
      if hasattr(item,"updateName"):
          item.updateName(str(value))
      else:
          item.setData(str(value))
      return True
    else:
      return False
    
  def removeRows(self, position=0, count=1, parent=QtCore.QModelIndex()):
    node = self.nodeFromIndex(parent)
    self.beginRemoveRows(parent, position, position + count - 1)
    node.childItems.pop(position)
    self.endRemoveRows()
    
  def nodeFromIndex(self, index):
    if index is not None and index.isValid():
      #print 'nodeFromIndex',index,index.internalPointer()
      return index.internalPointer()
    else:
      return self.rootItem
    
  def modelIndexFromJob(self,jobId=None):
    indexList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.UserRole,jobId,1)
    #print 'modelIndexFromJob',jobId,indexList
    if len(indexList)>0:
      return indexList[0]
    else:
      jobPath = PROJECTSMANAGER().db().getJobAncestors(jobId=jobId)
      if len(jobPath)<=1: return None
      modelIndex = QtCore.QModelIndex()
      for jid in jobPath:
        indexList = self.match(self.index(0,0,modelIndex),QtCore.Qt.UserRole,jid,1)
        if len(indexList)==0: return None
        modelIndex = indexList[0]
      return modelIndex
    
  def modelIndexFromFile(self,fileId=None,jobId=None,jobIndex=None):
    if jobIndex is None:
      if jobId is None:
        jobId = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode='jobid')
      jobIndex = self.modelIndexFromJob(jobId)
    #print 'modelIndexFromFile jobIndex',jobId,jobIndex
    if jobIndex is None: return None
    indexList = self.match(self.index(0,0,jobIndex),QtCore.Qt.UserRole,fileId,1)
    #print 'modelIndexFromFile',fileId,indexList
    if len(indexList)>0:
      return indexList[0]
    else:
      return None
    
  def getValue(self, index, role):
    item = index.internalPointer()
    return item.data(index.column(),role)
  
  def columnCount(self, parent):
    return CProjectModel.NCOLUMNS
  
  def data(self, index, role):
    if not index.isValid():
      return None    
    item = index.internalPointer()
    return item.data(index.column(),role)
  

  def flags(self,modelIndex):
    if not modelIndex.isValid(): return QtCore.Qt.NoItemFlags
    ic = modelIndex.column()
    #print 'CProjectModel.flags',ic
    if ic == 0:
      return QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|QtCore.Qt.ItemIsEditable
    else:
      return  QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled
  
    
  def headerData(self,section,orientation,role=QtCore.Qt.DisplayRole):
    if orientation == QtCore.Qt.Horizontal:
      if role==QtCore.Qt.DisplayRole:
        if self.COLUMNHEADERS[self.COLUMNS[section]] is not None:
          return self.COLUMNHEADERS[self.COLUMNS[section]]

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def index(self, row, column, parent):
    if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
      return QtCore.QModelIndex()
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    childItem = parentItem.child(row)
    #print 'CProjectModel.index',row, column, parent,childItem.getName()
    if childItem:
      return self.createIndex(row, column, childItem)
    else:
      return QtCore.QModelIndex()
    
  def parent(self, index):
    if not index.isValid():
      return QtCore.QModelIndex()
    childItem = index.internalPointer()
    parentItem = childItem.parent()
    if parentItem == self.rootItem:
      return QtCore.QModelIndex()
    return self.createIndex(parentItem.row(), 0, parentItem)
  
  def rowCount(self, parent):
    #if parent.column() > 0:
    #  return 0
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    return parentItem.childCount()

  
  def setupModelData(self):
    CTreeItemJob.TODAY = time.strftime(CTreeItemJob.DATE_FORMAT,time.localtime())
    CTreeItemJob.THISYEAR = time.strftime('%y',time.localtime())
    jobInfoList = PROJECTSMANAGER().db().getProjectJobListInfo(mode=self.JOBINFOITEMS,projectId=self._projectId)
    for jobInfo in jobInfoList:
      if jobInfo.get('parentjobid',None) is None:
        item = CTreeItemJob(parent=self.rootItem,info=jobInfo)
        #self.rootItem.appendChildJob(item)
        self.rootItem.insertChildJob(0,item)
      else:
        parent = self.getJobTreeItem(jobInfo['parentjobid'],self.rootItem)
        if parent is not None:
          item = CTreeItemJob(parent=parent,info=jobInfo)
          parent.appendChildJob(item)
          #print('setupModelData',parent.getName(),item.getName())
    fileList = PROJECTSMANAGER().db().getProjectFiles(projectId=self._projectId,topLevelOnly=False)
    for infoList in fileList:
      #jid,fid,impid,ftype,fname,annotation = infoList
      parent = self.getJobTreeItem(infoList[0],self.rootItem)
      if parent is not None:
        item = CTreeItemFile(parent=parent,infoList=infoList)
        if infoList[2] is not None:
          parent.appendChildInputFile(item)
        else:
          parent.appendChildOutputFile(item)
    followFrom = PROJECTSMANAGER().db().getProjectFollowFromJobId(projectId=self._projectId)
    if followFrom is not None: self.updateFollowFrom([self._projectId,None,followFrom])

  def getJobTreeItem(self,jobId,parentTreeItem):
    for item in parentTreeItem.childJobs:
      if item.jobId == jobId:
        return item
      elif len(item.childJobs)>0:
        rv = self.getJobTreeItem(jobId,item)
        if rv is not None: return rv
    return None
        
  def currentProjectId(self):
    return self._projectId

  @QtCore.Slot(dict)
  def createJob(self,args):
    #print 'CProjectModel.createJob',argList
    if args['projectId']  != self._projectId: return
    if self.getJobTreeItem(args['jobId'],self.rootItem) is not None: return    
    jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=args['jobId'],mode=CProjectModel.JOBINFOITEMS)
    #print 'CProjectModel.createJob',jobId,jobInfo
    if jobInfo['parentjobid'] is None:
      self.beginInsertRows(QtCore.QModelIndex(),0,0)
      item = CTreeItemJob(self.rootItem,jobInfo)
      self.rootItem.insertChildJob(0,item)
      self.endInsertRows()
      jobIndex = self.index(0,0,QtCore.QModelIndex())
    else:
      parentIndex = self.modelIndexFromJob(jobInfo['parentjobid'])
      if parentIndex is None:
        # Assume parent job started externally and not finished so not reported
        # by CCP4DbApi.getRecentlyFinishedJobs() - need to do a createJob() on it
        parentIndex = self.createJob({'jobId':jobInfo['parentjobid'],'projectId':args['projectId']})
      parentNode = self.nodeFromIndex(parentIndex)
      parentRows = parentNode.childCount()
      self.beginInsertRows(parentIndex,parentRows,parentRows)
      item = CTreeItemJob(parentNode,jobInfo)
      parentNode.appendChildJob(item)
      self.endInsertRows()
      jobIndex = self.index(parentRows,0,parentIndex)

    self.createImportFiles(args['jobId'],jobIndex,item)
    self.createOutputFiles(args['jobId'],jobIndex,item)

    #print 'from createJob'
    return jobIndex

  def createChildJobs(self,argList):
    #print 'createJob',jobId
    jobId,projectId = argList
    if projectId != self._projectId: return
    childJobs = PROJECTSMANAGER().db().getChildJobs(jobId,descendents=True)
    #print 'createChildJobs',childJobs
    def recurseCreateJob(jobList,model):
      for jid,childJobList in jobList:
        model.createJob({'jobId':jid,'projectId':projectId})
        if len(childJobList)>0: recurseCreateJob(childJobList,model)
    recurseCreateJob(childJobs,self)
    

  def createImportFiles(self,jobId=None,jobIndex=None,jobNode=None):
    # Get imported files
    fileInfoList =  PROJECTSMANAGER().db().getJobImportFiles(jobId=jobId)
    jobIndex = None
    for fileInfo in fileInfoList:
      #jobid,fileid,importid,filetype,filename,annotation
      if jobIndex is None:
        jobIndex = self.modelIndexFromJob(jobId)
        jobNode = self.nodeFromIndex(jobIndex)
      index = self.modelIndexFromFile(fileInfo[1],jobIndex=jobIndex)
      #print 'CProjectModel.createImportFiles',fileInfoList,index
      if index is None:
        item = CTreeItemFile(parent=jobNode,infoList=fileInfo)
        if fileInfo[2] is not None:
          nFiles = jobNode.childCount('infiles')
          self.beginInsertRows(jobIndex,nFiles,nFiles)
          jobNode.appendChildInputFile(item)
        else:
          nFiles = jobNode.childCount('outfiles')
          self.beginInsertRows(jobIndex,nFiles,nFiles)
          jobNode.appendChildOutputFile(item)
        self.endInsertRows()

  def createOutputFiles(self,jobId=None,jobIndex=None,jobNode=None):
      
    # Get list of fileIds
    fileInfoList = PROJECTSMANAGER().db().getJobFiles(jobId=jobId,mode='all')
    #print 'CProjectsModel.createOutputFiles',jobId,fileInfoList
    if len(fileInfoList) ==0 : return
    if jobIndex is None: jobIndex = self.modelIndexFromJob(jobId)
    if jobNode is None: jobNode = self.nodeFromIndex(jobIndex)
    # Move any files we already know about (assume this is pipeline job adopting an output file from child job)
    # Assume it is a child of a child job
    newFileList = []
    for fileInfo in fileInfoList:
      fileIndex = self.modelIndexFromFile(fileInfo[1],jobId=fileInfo[0])
      #print 'CProjectsModel.createOutputFiles',jobId,fileInfo,fileIndex
      if fileIndex is not None:
        #print 'CProjectsModel.createOutputFiles parent',fileIndex.parent().data(QtCore.Qt.UserRole).__str__()
        if fileIndex.parent().data(QtCore.Qt.UserRole).__str__() != fileInfo[0]:
          fileNode = self.nodeFromIndex(fileIndex)
          fromJobNode = self.nodeFromIndex(fileIndex.parent())
          fromRowIndex = fromJobNode.childFileIndex(self.nodeFromIndex(fileIndex))
          toRowIndex = jobNode.childCount(mode='outfiles')      
          #print 'CProjectModel.createOutputFiles moving file from',fromJobNode.jobNumber,'to',jobId
          self.beginMoveRows(fileIndex.parent(),fromRowIndex,fromRowIndex,jobIndex,toRowIndex)
          try:
            fromParentNode.removeChildFile(fileNode=fileNode)
          except:
            #print 'createJobOutputFiles error removing fileId from old job'
            pass
          newFileNode= CTreeItemFile(jobNode,fileInfo)
          jobNode.appendChildOutputFile(newFileNode)
          self.endMoveRows()
      else:
        newFileList.append(fileInfo)

    #print 'createJobOutputFiles',jobId,newFileIdList
    if len(newFileList)==0: return

    nFiles = jobNode.childCount('outfiles')
    self.beginInsertRows(jobIndex,nFiles,nFiles+len(newFileList)-1)
    for fileInfo in newFileList:
      item = CTreeItemFile(parent=jobNode,infoList=fileInfo)
      jobNode.appendChildOutputFile(item)
    self.endInsertRows()

  @QtCore.Slot(dict)
  def deleteJob(self,args):
    if args['projectId'] != self._projectId: return
    #print 'CProjectModel.deleteJob',args
    #print 'deleteJob in _jobIndexList',len(self._jobIndexList),self._jobIndexList
    
    modelIndex = self.modelIndexFromJob(jobId=args['jobId'])
    if modelIndex is None: return
    node=self.nodeFromIndex(modelIndex)
    row=modelIndex.row()
    parent = node.parent()
    self.beginRemoveRows(QtCore.QModelIndex(),row,row)
    parent.removeChildJob(row)
    self.endRemoveRows()
    #print 'from deleteJob'

  @QtCore.Slot(dict)
  def setJobToImport(self,args):
    #print 'setJobToImport',args
    if args['projectId'] != self._projectId: return
    modelIndex = self.modelIndexFromJob(jobId=args['jobId'])
    if modelIndex is None: return
    node=self.nodeFromIndex(modelIndex)
    #print 'CProjectModel.setJobToImport',node,node.childCount('infiles'),node.childCount()
    self.beginRemoveRows(modelIndex,node.childCount('infiles'),node.childCount())
    node.removeChildOutputFiles()
    node.removeChildJobs()
    self.endRemoveRows()
    #node.taskName = 'import_files'
    #print 'from setJobToImport'


  def backupDB(self,args):
    PROJECTSMANAGER().backupDB()

  @QtCore.Slot(dict)
  def updateFinishedJob(self,args):
    #print 'CProjectModel.updateFinishedJob',args
    newArgs = {}
    newArgs.update(args)
    newArgs['value'] = args['status']
    newArgs['key'] = 'status'
    self.updateJob(newArgs)
    
      
  @QtCore.Slot(dict)
  def updateJob(self,args):
    #print 'CProjectModel.updateJob',args
    if args.get('projectId','') != self._projectId: return
    jobId = args.get('jobId')
    key = args.get('key')
    value = args.get('value')
    #print 'CProjectModel.updateJob',jobId,key,value
    index = self.modelIndexFromJob(jobId)
    if index is None:
      index = self.createJob({'jobId':jobId,'projectId':args['projectId']})
    else:
      if key == 'status' and isinstance(value,int): value = CCP4DbApi.JOB_STATUS_TEXT[value]
      node = self.nodeFromIndex(index)
      #print 'CProjectModel.updateJob node',node,key,value
      if node is not None:
        if key == 'jobtitle':
          node.setName( value )
        else:
          node.set(key,value)
          if key == 'status':
            if value in ['Queued']:
              self.createImportFiles(jobId)
            elif value in ['Finished','Failed','Interrupted','To delete']:
              self.createOutputFiles(jobId)
              if value in ['Finished','Interrupted']:
                # What about evaluation?
                perfDict=PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode='performance')
                #print 'CProjectModel.updateJob perfDict',perfDict
                if perfDict is not None: node.set('performance',perfDict)
        # If there is jobTitle info (from CDbApi jobFinished signal) then update title
        if args.get('jobTitle',None) is not None:
          node.setName( args['jobTitle'] )
        if args.get('finishTime',None) is not None: node.set('finishtime',args['finishTime'])
        self.dataChanged.emit(index,index)
        self.redraw()
    #print 'from updateJob'

  @QtCore.Slot(dict)
  def updateFile(self,args):
    #print 'CProjectModel.updateFile',args
    if args['projectId'] != self._projectId or args['key'] != 'annotation' : return
    index = self.modelIndexFromFile(args['fileId'],jobId=args['jobId'])  
    if index is None: return
    node = self.nodeFromIndex(index)   
    node.setName(args['value'])
    self.dataChanged.emit(index,index)
    self.redraw()
    #print 'from updateFile'

    
  @QtCore.Slot(dict)
  def updateFollowFrom(self,args):
    #print 'CProjectModel.updateFollowFrom',args
    projectId,previousFollowFromJobId,jobId = args
    if projectId != self._projectId: return
    if previousFollowFromJobId is not None:
      index = self.modelIndexFromJob(jobId=previousFollowFromJobId)
      #print 'CProjectModel.updateFollowFrom previousFollowFromJobId index',index
      if index is not None:
        node = self.nodeFromIndex(index)
        if node is not None: node.set('followFrom',False)
        #ffindex = index.sibling(index.row(),0)
        self.dataChanged.emit(index,index)
    if jobId is not None:
      index = self.modelIndexFromJob(jobId=jobId)
      if index is not None:
        node = self.nodeFromIndex(index)
        if node is not None: node.set('followFrom',True)
        #ffindex = index.sibling(index.row(),0)
        self.dataChanged.emit(index,index)
    self.redraw()
    #print 'from updateFollowFrom'

  def updateHighlight(self,jobId=None):
    for jid,highlight in [ [ self.highlightJob,False ] , [ jobId , True ] ]:
      if jid is not None:
        index = self.modelIndexFromJob(jobId=jid)
        if index is not None:
          node = self.nodeFromIndex(index)
          node.set('highlight',highlight)
          self.dataChanged.emit(index,index)
    self.highlightJob = jobId
    self.redraw()
    
  def fileLabel(self,modelIndex=None,maxLength=None):
    node = self.nodeFromIndex(modelIndex)
    if node is None: return ''
    #print 'CProjectModel.fileLabel node',node
    return node.getName()

  def redraw(self):
    pass
    #QtCore.QAbstractItemModel.reset(self)
    self.redrawSignal.emit()

  def unsetJobColour(self):
    '''
    for job in self.rootItem.childJobs:
      job.colour = None
    '''
    self.rootItem.highlightList = {}

  def setJobColour(self,jobList,colour,unset=True):
    if unset: self.rootItem.highlightList= {}
    for job in jobList:
      self.rootItem.highlightList[job] = colour
    

class CProjectView(QtWidgets.QTreeView):

  rightMousePress = QtCore.Signal('QMouseEvent')
  jobClicked = QtCore.Signal('QModelIndex')
  fileClicked = QtCore.Signal('QModelIndex')

  """
  def changeEvent(self,e):
      #FIXME - *THIS* is the cause of my performance woes. I have to not react to every change event, just last one. Hmm ...
      print("changeEvent")
      import traceback
      traceback.print_stack()
      return QtWidgets.QTreeView.changeEvent(self,e)

      if e.type() == QtCore.QEvent.StyleChange or QtCore.QEvent.FontChange:
          self.resizeColumnToContents(0)
  """

  def __init__(self,parent=None):
    QtWidgets.QTreeView.__init__(self,parent)
    self.setObjectName('projectWidget')
    self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
    self.setDragEnabled(True)
    #self.setAcceptDrops(True)
    self.setDragDropMode(QtWidgets.QAbstractItemView.DragOnly)
    self.setExpandsOnDoubleClick(False)
    self.setRootIsDecorated(True)
    self.setIconSize(QtCore.QSize(16,16))
    self.setEditTriggers(QtWidgets.QAbstractItemView.EditKeyPressed)
    #self.setItemDelegateForColumn(2,CJobEditDelegate(self))
    #self.setFocusPolicy(QtCore.Qt.NoFocus)
    self.setToolTip('Right mouse click for options to view jobs and files')
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)
    PREFERENCES().TABLES_ALTERNATING_COLOR.dataChanged.connect(self.resetAlternatingRowColors)
    self.forceUpdate = sys.platform.count('inux') or sys.platform.count('arwin')
    #print 'CProjectView.__init__ forceUpdate',self.forceUpdate

    delegate = JobListHTMLDelegate(self)
    @QtCore.Slot()
    def editLabelsDone():
        if hasattr(self,"model") and hasattr(self.model(),"sourceModel") and hasattr(self.model().sourceModel(),"_projectId"):
            projectId = self.model().sourceModel()._projectId
            projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=projectId)
            dbxml = os.path.join(projectDir,"DATABASE.db.xml")
            print("Saving",dbxml)
            PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml)
    delegate.editingFinishedSignal.connect(editLabelsDone)
    self.setItemDelegate(delegate)
    self.setUniformRowHeights(False)

  def update(self):
    if not self.forceUpdate: return
    #print 'CProjectView.update'
    #QtWidgets.QTreeView.update(self)
    # Why am I doing this???  This is broken on Windows
    # Maybe this fixes the failure to update on Linux
    frameRect = self.frameRect()
    self.setDirtyRegion(QtGui.QRegion(frameRect.x(),frameRect.y(),frameRect.width(),frameRect.height()))
    #print 'from CProjectView.update'
  
    
  def resetAlternatingRowColors(self):
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)

  def keyPressEvent(self,event=None):
    '''
    From Qt Docs (../Qt/html/qt.html#KeyboardModifier-enum)
    Note: On Mac OS X, the ControlModifier value corresponds to the Command keys on the Macintosh keyboard,
    and the MetaModifier value corresponds to the Control keys.
    The KeypadModifier value will also be set when an arrow key is pressed as the arrow keys are considered part of the keypad.
    Note: On Windows Keyboards, Qt::MetaModifier and Qt::Key_Meta are mapped to the Windows key.
    Other possible values: QtCore.Qt.NoModifier
    modifier test by: event.modifiers() == QtCore.Qt.ControlModifier
    '''
    #print 'CProjectView.keyPressEvent',event.key()
    if event.key() == 16777265:
      # mapFromGlobal() seems to not allow for the header so need to subtract header height 
      modelIndex = self.indexAt(self.mapFromGlobal(self.cursor().pos()-QtCore.QPoint(0,self.header().height())))
      #print 'CProjectView.keyPressEvent modelIndex',modelIndex,modelIndex.isValid(),modelIndex.row()
      if modelIndex.isValid() and modelIndex.column()==2:
        event.accept()
        self.edit(modelIndex)
        return
    QtWidgets.QTreeView.keyPressEvent(self,event)
    
    
  def mousePressEvent(self,event=None):
    #print 'mousePressEvent'
    if event.button() == QtCore.Qt.RightButton:
      self.rightMousePress.emit(event)
      event.accept()
      
      return
    else:     
      mousePressX = event.x()
      modelIndex = self.indexAt(event.pos())
      r = self.visualRect(modelIndex)
      if r.isValid():
        # Weird: changed how handle selection and now don't seem to need to allow for indent!
        #indent = self.nestedLevel(modelIndex)*self.indentation()
        #print 'mousePressEvent indent',indent
        indent = 0
        #print 'mousePressEvent mousePressX',mousePressX,'left',r.left(),(r.left() + self.iconSize().width() + 4)
        if mousePressX  >r.left()+indent and mousePressX < (r.left() + self.iconSize().width() + 4 + indent):
          self.startDrag(modelIndex=modelIndex)
          event.accept()
          return
      if modelIndex.model() is None: return
      node = self.model().sourceModel().nodeFromIndex(modelIndex.model().mapToSource(modelIndex))      
      if node.isJob():
        self.jobClicked.emit(modelIndex)
      else:
        self.fileClicked.emit(modelIndex)
        event.accept()
        return
    #print 'calling QTreeView.mousePressEvent'
    QtWidgets.QTreeView.mousePressEvent(self,event)

  def nestedLevel(self,modelIndex):
    if not modelIndex.isValid(): return 0
    level = -1
    while level<5:
      level = level + 1
      modelIndex = modelIndex.parent()
      if not modelIndex.isValid(): return level
    return level

  def nodeFromEvent(self,event):
    #print("mapToSource 2")
    modelIndex = self.model().mapToSource(self.indexAt(QtCore.QPoint(event.x(),event.y())))
    col = self.model().sourceModel().COLUMNS[modelIndex.column()]
    return modelIndex,self.model().sourceModel().nodeFromIndex(modelIndex),col

    return [None,None,None,indx]
  
  def selectRow(self,modelIndex=None):
    #print 'CProjectView.selectRow',modelIndex.row()
    return
    sel = QtCore.QItemSelection( modelIndex.sibling(modelIndex.row(),0) , modelIndex.sibling(modelIndex.row(),CProjectModel.NCOLUMNS) )
    self.selectionModel().select(sel,QtCore.QItemSelectionModel.ClearAndSelect)
    #print 'CProjectView.selectRow DONE'

  
  def startDrag(self,dropActions=None,modelIndex=None):
    if modelIndex is None:
      modelIndex = self.currentIndex()
    if modelIndex is None: return
    #print("mapToSource 3")
    node = self.model().sourceModel().nodeFromIndex(self.model().mapToSource(modelIndex))
    if isinstance(node,CTreeItemFile) or node.getStatus() in ['Finished','noStatus']:
      drag = QtGui.QDrag(self)
      drag.setMimeData(node.mimeData())
    
      icon = node.data(CProjectModel.COLUMNS.index('name'),QtCore.Qt.DecorationRole)
      if icon is not None:
        try:
          pixmap = icon.pixmap(18,18)
          drag.setHotSpot(QtCore.QPoint(9,9))
          drag.setPixmap(pixmap)
        except:
          pass
      drag.exec_(QtCore.Qt.CopyAction)
     
  def getSelectedJobs(self):
    selectedRows = self.selectionModel().selectedRows()
    selectedJobs = []
    for row in selectedRows:
      #print("mapToSource 4")
      node = self.model().sourceModel().nodeFromIndex(self.model().mapToSource(row))
      #print 'getSelectedJobs', row,node
      if node is not None:
        selectedJobs.append(node.getJobId())
    #print 'getSelectedJobs',selectedJobs
    return selectedJobs

class CProjectProxyModel(QtCore.QSortFilterProxyModel):
        def __init__(self,parent=None):
            QtCore.QSortFilterProxyModel.__init__(self,parent)
            self.filterString = ""
            self.filterMode = ()
            self.filterStartDate = None
            self.filterEndDate = None

        def setFilterDates(self,startDate,endDate):
            self.filterStartDate = startDate
            self.filterEndDate = endDate
            self.invalidateFilter()

        def setFilterMode(self,string):
            self.filterMode = string
            self.invalidateFilter()

        def setFilterString(self,string):
            self.filterString = string
            self.invalidateFilter()

        def checkDates(self,item):
            #print "checkDates",item,self.filterStartDate,self.filterEndDate,item.data(4,QtCore.Qt.DisplayRole)
            accepted = True
            if hasattr(item,"getTopJob") and self.filterStartDate is not None and self.filterEndDate is not None:
                jobId = item.getTopJob().getJobId()
                info = PROJECTSMANAGER().db().getJobInfo(jobId)
                if "creationtime" in info and info["creationtime"] is not None:
                    t = datetime.datetime.strptime(self.filterStartDate,"%Y-%m-%d %H:%M:%S")
                    stamp = time.mktime(t.timetuple())
                    if info['creationtime']<stamp:
                        accepted = False
                if "finishtime" in info and info["finishtime"] is not None:
                    t = datetime.datetime.strptime(self.filterEndDate,"%Y-%m-%d %H:%M:%S")
                    stamp = time.mktime(t.timetuple())
                    if info['finishtime']>stamp:
                        accepted = False
                elif "creationtime" in info and info["creationtime"] is not None:
                    t = datetime.datetime.strptime(self.filterEndDate,"%Y-%m-%d %H:%M:%S")
                    stamp = time.mktime(t.timetuple())
                    if info['creationtime']>stamp:
                        accepted = False
            return accepted

        def hasAcceptedChildrenOrIsAccepted(self, item):

            #FIXME, this requires more comprehensive data searching.
            #print(item.data(0,QtCore.Qt.DisplayRole))
            #print(item.data(1,QtCore.Qt.DisplayRole))
            #print(item.data(0,QtCore.Qt.UserRole))
            theData = ""
            ur2 = item.data(0,QtCore.Qt.UserRole+2)
            if ur2 is not None and len(ur2)>0:
                theData = ur2
            theDataKeywords = ""#item.keywords()

            acceptedFilterMode = True

            if acceptedFilterMode and (len(str(self.filterString).strip()) == 0 or str(self.filterString).lower() in theData.lower()):
                return self.checkDates(item)

            if acceptedFilterMode and (theDataKeywords is not None and any(str(self.filterString).lower() in s for s in theDataKeywords)):
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

 
class CProjectWidget(QtWidgets.QFrame):

  interruptJob = QtCore.Signal(str)
  killJob = QtCore.Signal(str)
  killDeleteJob = QtCore.Signal(str)
  markFailedJob = QtCore.Signal(str)
  markFinishedJob = QtCore.Signal(str)
  openFrame = QtCore.Signal(str,str)
  purge = QtCore.Signal(str,str)
  export = QtCore.Signal(object,str)
  nextJob = QtCore.Signal(str,str,str)
  deleteJob = QtCore.Signal(list,bool)
  stopJob = QtCore.Signal(str)
  cloneTask = QtCore.Signal(str)
  showBibliography = QtCore.Signal(str)
  showWhatNextPage = QtCore.Signal(str)
  openMergeMtz = QtCore.Signal(dict)
  currentJobChanged = QtCore.Signal(str,str,str,bool)

  MARGIN = 0

  ERROR_CODES = { 100 : { 'description' : 'Error copying file' } }
  
  def  __init__(self,parent=None,projectId=None):
    QtWidgets.QFrame.__init__(self)
    self.projectId = projectId

    layout = QtWidgets.QVBoxLayout()
    layout.setMargin(CProjectWidget.MARGIN)
    layout.setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    #layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
    self.setLayout(layout)

    iconDir = os.path.join(os.environ['CCP4I2_TOP'],'qticons')
    upArrow = QtGui.QIcon(os.path.join(iconDir,"up.png"))
    downArrow = QtGui.QIcon(os.path.join(iconDir,"down.png"))
    self.taskSearchBox = QtWidgets.QLineEdit()

    searchLayout = QtWidgets.QHBoxLayout()
    self.taskSearchBox.setPlaceholderText("Only show jobs containing text typed here")
    searchLayout.addWidget(QtWidgets.QLabel("Filter:"))
    searchLayout.addWidget(self.taskSearchBox)
    advButton = QtWidgets.QPushButton()
    advButton.setIcon(downArrow)
    if sys.platform == "darwin":
        advButton.setIconSize(QtCore.QSize(12,12))
    searchLayout.addWidget(advButton)

    advSearchLayout = QtWidgets.QGridLayout()
    dateLabel = QtWidgets.QLabel("Created/finished between:")
    dateLabelAnd = QtWidgets.QLabel("and:")
    projectInfo = PROJECTSMANAGER().db().getProjectInfo(self.projectId)
    creationDate = projectInfo['projectcreated']
    dt = datetime.datetime.utcfromtimestamp(creationDate)
    qdt = QtCore.QDateTime.fromString(dt.strftime("%Y-%m-%d %H:%M:%S"),"yyyy-MM-dd HH:mm:ss")
    dateStart = QtWidgets.QDateTimeEdit(qdt)
    dateEnd = QtWidgets.QDateTimeEdit(QtCore.QDateTime.currentDateTime())
    advSearchLayout.addWidget(dateLabel,1,0)
    advSearchLayout.addWidget(dateLabelAnd,2,0)
    advSearchLayout.addWidget(dateStart,1,1)
    advSearchLayout.addWidget(dateEnd,2,1)
    dateStartReset = QtWidgets.QPushButton("Project start")
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
    def toggleAdvanced():
        if dateLabel.isVisible():
            advWidget.hide()
            advButton.setIcon(downArrow)
        else:
            advWidget.show()
            advButton.setIcon(upArrow)
    advButton.clicked.connect(toggleAdvanced)
    advButton.setMaximumSize(QtCore.QSize(22,22))

    self.tab = QtWidgets.QTabWidget(self)
    self.layout().addWidget(self.tab)

    model = PROJECTMODEL(projectId=projectId)
    print(model)
    projectProxyModel = CProjectProxyModel()
    projectProxyModel.setSourceModel(model)
    self.projectView=CProjectView(self)
    self.projectView.setModel(projectProxyModel)
    model.redrawSignal.connect(self.projectView.update)
    #self.projectView.model().reset()
    @QtCore.Slot('QModelIndex','QModelIndex')
    def resizeColumns(topLeft,bottomRight):
        self.projectView.resizeColumnToContents(0)
    model.dataChanged.connect(resizeColumns)
    self.projectView.setHeaderHidden(True)

    self.taskSearchBox.textChanged.connect(projectProxyModel.setFilterString)
    @QtCore.Slot()
    def dateChanged():
       startDate = str(dateStart.dateTime().toString("yyyy-MM-dd HH:mm:ss"))
       endDate = str(dateEnd.dateTime().toString("yyyy-MM-dd HH:mm:ss"))
       projectProxyModel.setFilterDates(startDate,endDate)
    dateStart.dateTimeChanged.connect(dateChanged)
    dateEnd.dateTimeChanged.connect(dateChanged)

    #self.projectView.setColumnWidth(model.COLUMNS.index('name'),2000)

    self.projectView.rightMousePress.connect(self.showJobListPopup)

    searchLayout.setMargin(0)
    searchLayout.setContentsMargins(0,0,0,0)
    advSearchLayout.setMargin(0)
    advSearchLayout.setContentsMargins(0,0,0,0)
    advSearchLayout.setSpacing(0)

    self.historyGui = CHistoryGui(self)
    self.historyGui.unSet()
    self.historyGui.clearButton.clicked.connect(self.clearHistory)
    self.historyGui.slider.valueChanged[int].connect(self.handleHistorySlider)
    if parent is not None: parent.jobSearch.connect(self.showJobSearch)
    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QVBoxLayout())
    frame.layout().addLayout(searchLayout)
    frame.layout().addWidget(advWidget)
    frame.layout().setMargin(CProjectWidget.MARGIN)
    frame.layout().setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    frame.layout().addWidget(self.historyGui)
    frame.layout().addWidget(self.projectView)
    self.tab.addTab(frame,'Job list')

    self.dirModel = CProjectDirModel(self,projectId)
    self.proxyModel = CProjectDirProxyModel(self)
    self.proxyModel.setSourceModel(self.dirModel)

    self.dirView = CProjectDirView(self)
    self.dirView.setModel(self.proxyModel)
    self.dirView.setSortingEnabled(True)
    self.dirView.sortByColumn(0,QtCore.Qt.DescendingOrder)
    projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=self.projectId)
    #print 'CProjectDirModel',projectName,projectDir
    modelIndex = self.proxyModel.mapFromSource(self.dirModel.index(os.path.join(projectDir,'CCP4_JOBS')))
    if modelIndex.isValid():
      self.dirView.setRootIndex(modelIndex)
    self.projectView.fileClicked.connect(self.handleFileClicked)
    self.projectView.selectionModel().selectionChanged.connect(self.handleSelectionChanged)
    self.projectView.doubleClicked.connect(self.handleDoubleClick)
    self.dirView.doubleClicked.connect(self.showFileFromDirView)
    self.tab.addTab(self.dirView,'Project directory')    
    if sys.platform == "darwin": self.dirView.setIconSize(QtCore.QSize(14,14))

    self.popupMenu = QtWidgets.QMenu(self)
    self.doubleClicked = None
    self.projectView.resizeColumnToContents(0)
    @QtCore.Slot()
    def resizeColumnsFromPrefs():
        try:
            self.projectView.resizeColumnToContents(0)
        except RuntimeError:
            #This is to avoid some PyQt horribleness, sorry.
            pass
    QTAPPLICATION().prefsChanged.connect(resizeColumnsFromPrefs)

  def model(self):
      return self.projectView.model()

  @QtCore.Slot()
  def clearHistory(self):
    self.model().sourceModel().unsetJobColour()
    self.projectView.update()
    self.historyGui.unSet()

  @QtCore.Slot(int)
  def handleHistorySlider(self,value):
    #print 'handleHistorySlider',value
    self.setJobColour(self.historyLeafList[value-1],self.historySelection)
    self.projectView.update()
    
  @QtCore.Slot(bool)
  def showJobSearch(self,show=True):
    if not show:
     if getattr(self,'jobSeachGui',None) is not None:
        self.jobSeachGui.hide()
     return
    if getattr(self,'jobSeachGui',None) is None:
      projectName = PROJECTSMANAGER().db().getProjectInfo(self.projectId,mode='projectname')
      self.jobSeachGui = CJobSearchDialog(self,projectName=projectName)
    self.jobSeachGui.show()
    self.jobSeachGui.raise_()

  def setForceUpdate(self,value):
    self.projectView.forceUpdate = value
    
  @QtCore.Slot('QMouseEvent')
  def showJobListPopup(self,event):
    from core.CCP4TaskManager import TASKMANAGER
    modelIndex,node,column = self.projectView.nodeFromEvent(event)
    #print 'showJobListPopup',node.getName(),column
    position = QtCore.QPoint(event.globalX(),event.globalY())
    
    if node is None: return
    itemName = node.getName()
    selectedJobs = []
    if column in ['evaluation','name','performance','dateTime']:
      if node.isJob():
        try:
          selectedJobs = self.projectView.getSelectedJobs()
        except:
          selectedJobs = []
        try:
          taskname = node.getTaskName()
          jobStatus = node.getStatus()
          jobId=node.jobId
          importFiles = node.childCount('infiles')
        except:
          return
        if jobId is None: return
          
        #print 'showJobListPopup',taskname,type(taskname),jobStatus
        self.popupMenu.setTitle(itemName)
        self.popupMenu.clear()
        if taskname == 'import_files':
          action = self.popupMenu.addAction('Delete job and imported files')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'Delete+import',jobId,itemName,modelIndex,None,selectedJobs))
        else:
          if jobStatus in ['Running','Running remotely']:
            popupStopSubMenu = self.popupMenu.addMenu('Stop job')
            if jobStatus != 'Running remotely' and TASKMANAGER().getTaskAttribute(taskName=taskname,attribute='INTERRUPTABLE',default=False,script=True) and taskname != 'coot_rebuild': 
                action = popupStopSubMenu.addAction('Interrupt safely')
                action.triggered.connect(functools.partial(self.interruptJob.emit,jobId))
            if taskname != 'coot_rebuild':
                action = popupStopSubMenu.addAction('Kill immediately')
                action.triggered.connect(functools.partial(self.killJob.emit,jobId))
            action = popupStopSubMenu.addAction('Kill and delete')
            action.triggered.connect(functools.partial(self.killDeleteJob.emit,jobId))
            action = popupStopSubMenu.addAction('Mark as failed')
            action.triggered.connect(functools.partial(self.markFailedJob.emit,jobId))
            action = popupStopSubMenu.addAction('Mark as finished')
            action.triggered.connect(functools.partial(self.markFinishedJob.emit,jobId))

          else:
            if len(selectedJobs)>1 and jobId in selectedJobs:
              action = self.popupMenu.addAction('Delete selected jobs')
              action.triggered.connect(functools.partial(self.handleJobListPopup,'Delete+selected',jobId,itemName,modelIndex,None,selectedJobs))
            elif importFiles == 0:
              action = self.popupMenu.addAction('Delete')
              action.triggered.connect(functools.partial(self.handleJobListPopup,'Delete',jobId,itemName,modelIndex,None,selectedJobs))
            else:
              popupDelMenu = self.popupMenu.addMenu('Delete')
              for label,mode in [['Delete job - save imported files','Delete'],['Delete job - delete imported files','Delete+import']]:
                action = popupDelMenu.addAction(label)
                action.triggered.connect(functools.partial(self.handleJobListPopup,mode,jobId,itemName,modelIndex,None,selectedJobs))
          # Clone - disable if not CLONEABLE
          action = self.popupMenu.addAction('Clone')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'Clone',jobId,itemName,modelIndex,None,selectedJobs))
          action.setEnabled(TASKMANAGER().getTaskAttribute(taskName=taskname,attribute='CLONEABLE',default=True) )

          action = self.popupMenu.addAction('Copy parameters')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'Copy parameters',jobId,itemName,modelIndex,None,selectedJobs))
          action.setEnabled(node.isTopJob())

          action = self.popupMenu.addAction('Edit label')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'Edit label',jobId,itemName,modelIndex,position,selectedJobs))
          popupOpenSubMenu = self.popupMenu.addMenu('Open')
          for frame,label in [['input','Input'],['output','Output']]:
            action = popupOpenSubMenu.addAction(label)
            action.setEnabled(jobStatus != CCP4DbApi.JOB_STATUS_FILE_HOLDER)
            action.triggered.connect(functools.partial(self.openFrame.emit,frame,jobId))
          action = popupOpenSubMenu.addAction('Comments')
          action.triggered.connect(functools.partial(self.openFrame.emit,'status',jobId))
          popupViewSubMenu = self.popupMenu.addMenu('View')
          action = popupViewSubMenu.addAction('Job report')
          #print 'showJobListPopup jobStatus',jobStatus
          if jobStatus in ['Finished','Failed','Interrupted']:
            action.triggered.connect(functools.partial(self.handleJobListPopup,'Job report',jobId,itemName,modelIndex,None,selectedJobs))
          else:
            action.setEnabled(False)
          action = popupViewSubMenu.addAction('Command file')
          if os.path.exists(PROJECTSMANAGER().makeFileName(jobId,'COM')):
            action.triggered.connect(functools.partial(self.handleJobListPopup,'Command file',
                                                                    jobId,itemName,modelIndex,None,selectedJobs))
          else:
            action.setEnabled(False)
          action = popupViewSubMenu.addAction('Log file')
          action1 = popupViewSubMenu.addAction('Log graphs')
          logFile = PROJECTSMANAGER().makeFileName(jobId,'LOG')
          #print 'showJobListPopup logFile',logFile
          if os.path.exists(logFile) or os.path.exists(logFile+'.html'):
            action.triggered.connect(functools.partial(self.handleJobListPopup,'Log file',jobId,itemName,modelIndex,None,selectedJobs))
            action1.triggered.connect(functools.partial(self.handleJobListPopup,'Log graphs',jobId,itemName,modelIndex,None,selectedJobs))
          else:
            action.setEnabled(False)
            action1.setEnabled(False)
          action = popupViewSubMenu.addAction('Diagnostic')
          if jobStatus in ['Finished','Failed','Interrupted']:
            #diagFile = PROJECTSMANAGER().makeFileName(jobId,'DIAGNOSTIC')
            #if os.path.exists(diagFile):  
            action.triggered.connect(functools.partial(self.handleJobListPopup,'Diagnostic',jobId,itemName,modelIndex,None,selectedJobs))
          else:
            action.setEnabled(False)
          for label,prog in [['CCP4mg','ccp4mg'],['Coot','coot']]:
            action = popupViewSubMenu.addAction('In '+label)
            if jobStatus in ['Finished','Failed','Interrupted']:
              action.triggered.connect(functools.partial(self.handleJobListPopup,prog,jobId,itemName,modelIndex,None,selectedJobs))
            else:
              action.setEnabled(False)
          action = popupViewSubMenu.addAction('Job directory')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'job_dir',jobId,itemName,modelIndex,None,selectedJobs))
          action = popupViewSubMenu.addAction('Database entry')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'job_db',jobId,itemName,modelIndex,None,selectedJobs))
          action = popupViewSubMenu.addAction('Bibliography')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'job_references',jobId,itemName,modelIndex,None,selectedJobs))
          action.setEnabled(True)

          self.popupHistoryMenu = self.popupMenu.addMenu('Data history')
          self.loadHistoryMenu(jobId)
          
          popupPurgeMenu = self.popupMenu.addMenu('Cleanup files')
          action = popupPurgeMenu.addAction('Delete temporary files')
          action.triggered.connect(functools.partial(self.purge.emit,jobId,'temporary'))
          action = popupPurgeMenu.addAction('Delete temporary and sub-job data files')
          action.triggered.connect(functools.partial(self.purge.emit,jobId,'intermediate'))          
          if jobStatus in ['Finished','Interrupted','Failed'] :
            popupViewSubMenu = self.popupMenu.addMenu('Export')
            action = popupViewSubMenu.addAction('All job files - compressed' )         
            action.triggered.connect(functools.partial(self.export.emit,'job',jobId))
            if jobStatus in ['Finished','Interrupted'] :
              exportMenu = TASKMANAGER().exportJobFiles(taskname,jobId=jobId)
              for item in exportMenu:
                action = popupViewSubMenu.addAction(item[1])         
                action.triggered.connect(functools.partial(self.export.emit,item,jobId))
          
          popupViewSubMenu = self.popupMenu.addMenu('Next task..')
          # What next...
          nextOptionList = TASKMANAGER().whatNext(taskname,jobId=jobId)
          for nextTask,nextTitle,nextDefFile in nextOptionList:
            action =  popupViewSubMenu.addAction(nextTitle)
            action.triggered.connect(functools.partial(self.nextJob.emit,nextTask,jobId,nextDefFile))

          popupViewSubMenu.addSeparator()
          action = popupViewSubMenu.addAction('Help..')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'What next?',jobId,itemName,modelIndex,None,selectedJobs))

          popupMenuEvalution = self.popupMenu.addMenu('Evaluation')
          for item in CCP4DbApi.JOB_EVALUATION_TEXT:
              action = popupMenuEvalution.addAction(jobIcon(item),item)
              action.triggered.connect(functools.partial(self.handleEvaluationPopup,item,jobId))

          popupMenuFollowFrom = self.popupMenu.addMenu('Follow from status')
          for menuItem,icon in [['Follow from',jobIcon('greenarrowsup')],['Clear',jobIcon('greendot')]]:
              action = popupMenuFollowFrom.addAction(icon,menuItem)
              action.triggered.connect(functools.partial(self.handleFollowFromPopup,menuItem,jobId))

          if jobStatus != 'Finished':
              popupMenuEvalution.setEnabled(False)
              popupMenuFollowFrom.setEnabled(False)

      elif node.isFile():
        fileId = node.fileId
        self.popupMenu.setTitle(itemName)
        self.popupMenu.clear()
        popupViewSubMenu = self.popupMenu.addMenu('View')
        #ext = node.getFilename().split('.')[-1]
        if node.fileType == 2:
          subMenuDef =   [['text','As text'],['ccp4mg','In CCP4mg'],['coot','In Coot'],['db','Database entry']]
        elif node.fileType in [4,5,6,10,11,12,13,16]:
          subMenuDef =   [['viewhkl','In ViewHKL'],['db','Database entry'],['coot','In Coot']]
        elif node.fileType in [ 13]:
          subMenuDef =   [['viewhkl','In ViewHKL'],['ccp4mg','In CCP4mg'],['coot','In Coot'],['db','Database entry']]
        elif node.fileType == 1:
          subMenuDef =   [['text','As text'],['db','Database entry']]
        else:
          subMenuDef =   [['text','As text'],['db','Database entry']]
        for alias,label in subMenuDef:
          action = popupViewSubMenu.addAction(label)
          action.triggered.connect(functools.partial(self.handleJobListPopup,'file_'+alias,fileId,itemName,modelIndex,None,selectedJobs))
        action = self.popupMenu.addAction('Edit label')
        action.triggered.connect(functools.partial(self.handleJobListPopup,'file_edit_label',fileId,itemName,modelIndex,position,selectedJobs))
        if  node.fileType in [4,5,6,10,11,12,13,16]:
          popupExportMenu = self.popupMenu.addMenu('Export file')
          action = popupExportMenu.addAction('Only this file')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'file_export',fileId,itemName,modelIndex,None,selectedJobs))
          action = popupExportMenu.addAction('All input/output exptal data for this job')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'file_exptdata_export',fileId,itemName,modelIndex,None,selectedJobs))
          action = popupExportMenu.addAction('Select data to export')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'file_select_exptdata_export',fileId,itemName,modelIndex,None,selectedJobs))
        else:
          action = self.popupMenu.addAction('Export file')
          action.triggered.connect(functools.partial(self.handleJobListPopup,'file_export',fileId,itemName,modelIndex,None,selectedJobs))
      self.popupMenu.popup(QtCore.QPoint(event.globalX(),event.globalY()))

  def loadHistoryMenu(self,jobId):
    self.popupHistoryMenu.clear()
    fileList = PROJECTSMANAGER().db().getJobFiles(jobId,CCP4DbApi.FILE_ROLE_IN,mode='all')
    fileList0 = PROJECTSMANAGER().db().getJobFiles(jobId,CCP4DbApi.FILE_ROLE_OUT,mode='all')
    for f0 in fileList0:
      new=True
      for f in fileList:
        if f0[3] == f[3]:
          new = False
          break
      if new: fileList.append(f0)
    for f in fileList:
        if f[3] is not None:
          action = self.popupHistoryMenu.addAction(self.historyFileTypeText(f[3]))
          action.triggered.connect(functools.partial(self.showHistory,jobId,f[3],role=None))

  @QtCore.Slot(str,str,str,str)
  def showHistory(self,jobId,fileType,role=None,selection=None):
    import copy
    #print 'showHistory',jobId,fileType,role,selection
    jobNum = PROJECTSMANAGER().db().getJobInfo(jobId,mode='jobnumber')
    self.historyLeafList = []
    self.historySelection = selection
    preceedingJobList = []
    tmpLeafList = []
    def sortPreceedingJobs(jTree,history):
        #print 'job ',jTree[0],jTree[1],history
        preceedingJobList.append(jTree[1])
        if len(jTree[2]) > 0:
          for item in jTree[2]:
            sortPreceedingJobs(item,history+[jTree[1]])
    def sortSuceedingJobs(jTree,history):
        #print 'job ',jTree[0],jTree[1],history
        if len(jTree[2]) == 0:
          tmpLeafList.append(history+[jTree[1]])
        else:
          for item in jTree[2]:
            sortSuceedingJobs(item,history+[jTree[1]])
    def compareHistoryLeaf(aList,bList):
      return int(aList[-1]) - int(bList[-1])
    if role == CCP4DbApi.FILE_ROLE_IN:
      jobsTree = PROJECTSMANAGER().db().getPreceedingJobs(jobId,fileType=fileType)
      sortPreceedingJobs(jobsTree,[jobNum])
      #print 'showHistory preceedingJobList',preceedingJobList
      self.setJobColour(preceedingJobList,self.historySelection)
      self.historyGui.set('History of job '+str(jobNum)+' '+self.historyFileTypeText(fileType))
      self.historyGui.unsetSlider()
    elif role == CCP4DbApi.FILE_ROLE_OUT or role is None:
      #if role is None:
      if 1:
        jobsTree = PROJECTSMANAGER().db().getPreceedingJobs(jobId,fileType=fileType)
        sortPreceedingJobs(jobsTree,[jobNum])
        preceedingJobList.reverse()
      jobsTree = PROJECTSMANAGER().db().getSucceedingJobs(jobId,fileType=fileType)
      sortSuceedingJobs(jobsTree,preceedingJobList)
      # historyLeafList is list of list of jobs to each final 'leaf' job
      if sys.version_info > (3,0):
        self.historyLeafList = sorted(tmpLeafList,key=functools.cmp_to_key(compareHistoryLeaf))
      else:
        self.historyLeafList = sorted(tmpLeafList,cmp=compareHistoryLeaf)
      #print 'showHistory leafList',self.historyLeafList
      if len(self.historyLeafList)<=1:
        self.historyGui.set('History of job '+str(jobNum)+' '+self.historyFileTypeText(fileType))
        self.historyGui.unsetSlider()
        self.setJobColour(self.historyLeafList[0],self.historySelection)
      else:
        # Append a list of all successive jobs to historyLeafList
        allBranchesList = copy.deepcopy(self.historyLeafList[0])
        for branch in self.historyLeafList[1:]:
          for jN in branch:
            if not jN in allBranchesList: allBranchesList.append(jN)
        self.historyLeafList.append(sorted(allBranchesList))
        #print 'allBranchesList',self.historyLeafList[-1]
        self.historyGui.set('History of job '+str(jobNum)+' '+self.historyFileTypeText(fileType))
        self.historyGui.setSlider(self.historyLeafList)
        self.setJobColour(self.historyLeafList[0],self.historySelection)
    self.projectView.update()
    return

  def setJobColour(self,jobList,selection=None):
    if selection is None:
      self.model().sourceModel().setJobColour(jobList,'red')
    else:
      l0 = []
      l1 = []
      for j in jobList:
        if j in  selection:
          l0.append(j)
        else:
          l1.append(j)
      self.model().setJobColour(l0,'red')
      self.model().setJobColour(l1,'pink',False)
      
  def historyFileTypeText(self,fileTypeId):
    return CCP4DbApi.FILETYPELIST[fileTypeId][2]
    

  @QtCore.Slot(str,str,bool)
  def handleEvaluationPopup(self,evaluation,jobId,triggerBool=None):
    #print 'handleEvaluationPopup',evaluation,jobId
    try:
      PROJECTSMANAGER().db().updateJob(jobId,key='evaluation',value=evaluation)
      from dbapi import CCP4DbUtils
      CCP4DbUtils.makeJobBackup(jobId=jobId)
    except:
      print('ERROR in handleEvaluationPopup')

  @QtCore.Slot(str,str,bool)
  def handleFollowFromPopup(self,action,jobId,triggerBool=None):
    #print 'handleFollowFromPopup',action,jobId
    if action == 'Follow from':
      PROJECTSMANAGER().db().setProjectFollowFromJobId(self.model().sourceModel()._projectId,jobId)
    elif action == 'Clear':
      PROJECTSMANAGER().db().setProjectFollowFromJobId(self.model().sourceModel()._projectId,jobId,clear=True)
      
  @QtCore.Slot(str,str,str,'QModelIndex','QPos',list,bool)
  def handleJobListPopup(self,action,jobId,itemName,modelIndex=None,position=None,selectedJobs=None,triggerBool=None):
    #print 'CProjectWidget.handleJobListPopup',action,jobId,selectedJobs
    #print 'handleJobListPopup current', self.projectView.selectionModel().currentIndex().row(), self.projectView.selectionModel().selectedRows(),self.projectView.selectionModel().selectedIndexes()
    if action == 'Delete':
      self.deleteJob.emit([jobId],False)
    elif action == 'Delete+import':
      self.deleteJob.emit([jobId],True)
    elif action == 'Delete+selected':
      #print 'handleJobListPopup Delete+selected'
      self.deleteJob.emit(selectedJobs,True)
    elif action == 'Stop job':
      self.stopJob.emit(jobId)
    elif action == 'Clone':
      self.cloneTask.emit(jobId)
    elif action == 'Copy parameters':
      # Set 'taskParameters' data on the application clipboard
      # to enable it to be pasted elsewhere
      from lxml import etree
      root = etree.Element('taskParameters')
      jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId,mode=['taskname','jobnumber','projectname','projectid'])
      for name,value in [[ 'jobId' , jobId],['taskName',jobInfo['taskname']],['jobNumber',jobInfo['jobnumber']],['projectName',jobInfo['projectname']],['projectId',jobInfo['projectid']]]:
        e = etree.SubElement(root,name)
        e.text = value
      dragText = etree.tostring(root,pretty_print=True)
      data = QtCore.QByteArray()
      data.append(dragText)
      mimeData = QtCore.QMimeData()
      mimeData.setData('taskParameters_'+jobInfo['taskname'],data)
      QTAPPLICATION().clipboard().setMimeData(mimeData)
    elif action == 'Edit label':
      #from qtgui import CCP4Widgets
      #d = CCP4Widgets.CEditFileLabel(parent=self,jobId=jobId)
      #d.move(position-QtCore.QPoint(20,20))
      #print 'handleJobListPopup modelIndex',modelIndex.row(),modelIndex.column()
      #self.projectView.edit(modelIndex)
      self.projectView.edit(self.projectView.model().mapFromSource(modelIndex))
    elif action == 'Command file':
      comFile = PROJECTSMANAGER().makeFileName(jobId,'COM')
      if os.path.exists(comFile): WEBBROWSER().openFile(comFile,toFront=True)
    elif action == 'Log file':
      logFile = PROJECTSMANAGER().makeFileName(jobId,'LOG')
      #print 'CProjectViewer.handleJobListPopup logFile',logFile
      if os.path.exists(logFile):
        LAUNCHER().launch(viewer='logview',argList=[logFile])
        #WEBBROWSER().openFile(logFile,format="text/plain",toFront=True)
      elif os.path.exists(logFile+'.html'):
        WEBBROWSER().openFile(logFile+'.html',toFront=True)
    elif action == 'Log graphs':
      logFile = PROJECTSMANAGER().makeFileName(jobId,'LOG')
      if os.path.exists(logFile):
        LAUNCHER().launch(viewer='loggraph',argList=[logFile])
    elif action == 'Job report':
      reportFile = PROJECTSMANAGER().makeFileName(jobId,'REPORT')
      if os.path.exists(reportFile):
        WEBBROWSER().openFile(reportFile,toFront=True)
      else:
        try:
          from report import CCP4ReportGenerator
          jobNumber =  PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode='jobnumber')
          generator = CCP4ReportGenerator.CReportGenerator(jobId=jobId,jobStatus='Finished',jobNumber=jobNumber)
          reportFile, newPageOrNewData = generator.makeReportFile()
        except CException as e:
          print(e.report())
        except Exception as e:
          print(e)
        if os.path.exists(reportFile):
          webView = WEBBROWSER().openFile(reportFile,toFront=True)
          #if webView is not None: generator.FinishedPictures.connect(webView.attemptImageLoad)
    elif action == 'Diagnostic':
      #diagFile = PROJECTSMANAGER().makeFileName(jobId,'DIAGNOSTIC')
      jobInfo =  PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber','status'])
      from report import CCP4ReportGenerator
      generator = CCP4ReportGenerator.CReportGenerator(jobId=jobId,jobStatus=jobInfo['status'],jobNumber=jobInfo['jobnumber'])
      reportFile = generator.makeFailedReportFile(redo=False)
      WEBBROWSER().openFile(reportFile,toFront=True)
    elif action in ['ccp4mg','coot']:
      if action == 'coot': action='coot_job'
      LAUNCHER().openInViewer(viewer=action.lower(),jobId=jobId,projectId=self.projectView.model().sourceModel()._projectId,guiParent=self)
    elif action == 'job_dir':
      self.dirView.focusOn(jobNumber=itemName)
      self.tab.setCurrentIndex(1)
    elif action == 'job_db':
      self.showDatabaseEntry(jobId=jobId)
    elif action == 'job_references':
      self.showBibliography.emit(jobId)
    elif action == 'What next?':
      self.showWhatNextPage.emit(jobId)
    elif action[0:4] == 'file':
        fileId = jobId
        fileName = PROJECTSMANAGER().db().getFullPath(fileId)
        if fileName is not None:
          if action == 'file_text':
            WEBBROWSER().openFile(fileName,toFront=True)
          elif action == 'file_db':
            self.showDatabaseEntry(fileId=fileId)
          elif action in ['file_ccp4mg','file_coot','file_viewhkl']:
            if action == 'file_coot': action='file_coot_job'
            LAUNCHER().openInViewer(viewer=action[5:],fileName=str(fileName),projectId=self.projectView.model().sourceModel()._projectId,guiParent=self)
          elif action == 'file_edit_label':
            #self.projectView.edit(modelIndex)
            self.projectView.edit(self.projectView.model().mapFromSource(modelIndex))
            #fileLabel = self.projectView.model().fileLabel(modelIndex=modelIndex,maxLength=None)
            #from qtgui import CCP4Widgets
            #d = CCP4Widgets.CEditFileLabel(parent=self,fileId=fileId)
            #d.move(position-QtCore.QPoint(20,20))
          
          elif action in ['file_export','file_exptdata_export']:
            fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['filetype','jobid','jobnumber'])
            filters = MIMETYPESHANDLER().getMimeTypeInfo(fileInfo['filetype'],'filter')
            defaultSuffix =  MIMETYPESHANDLER().getMimeTypeInfo(fileInfo['filetype'],'fileExtensions')[0]
            if action == 'file_export':
              title = 'Export '+fileName
            else:
              title = 'Export all experimental data associated with job '+str(fileInfo['jobnumber'])
            #print 'file_export',fileType,filters,defaultSuffix
            from qtgui import CCP4FileBrowser
            fileBrowser = CCP4FileBrowser.CFileDialog(parent=self,
                                      title=title,
                                     filters = [filters],
                                      defaultSuffix = defaultSuffix,
                                      fileMode = QtWidgets.QFileDialog.AnyFile)
            if action == 'file_export':
              fileBrowser.selectFile.connect(functools.partial(self.exportData,fileName,fileId,None))
            else:
              fileBrowser.selectFile.connect(functools.partial(self.exportData,fileName,None,fileInfo['jobid']))
            fileBrowser.show()
          elif action == 'file_select_exptdata_export':
            self.openMergeMtz.emit({'fileName':fileName,'fileId':fileId,'jobId':jobId})
    else:
      pass
    #print 'CProjectWidget.handleJobListPopup DONE'


  def launchViewer(self,filePath):
    from core import CCP4Modules
    format = MIMETYPESHANDLER().formatFromFileExt(fileName=filePath)
    viewerList = MIMETYPESHANDLER().getViewers(format)
    #print 'CProjectWidget.launchViewer',filePath,format,viewerList
    if len(viewerList)<=0:
      CCP4Modules.WEBBROWSER().openFile(filePath,toFront=True)
    elif isinstance(viewerList[0],str):
      CCP4Modules.LAUNCHER().openInViewer(viewer=viewerList[0],fileName=filePath,projectId=self.projectView.model().sourceModel()._projectId,guiParent=self)
    else:
      from qtgui import CCP4WebBrowser
      CCP4WebBrowser.OPENFILE(filePath,toFront=True)

  def showDatabaseEntry(self,jobId=None,fileId=None):
    import copy
    win = QtWidgets.QDialog(self)
    win.setWindowTitle('Database entry')
    win.setLayout(QtWidgets.QVBoxLayout())
    label = QtWidgets.QTextEdit()
    win.layout().addWidget(label)
    
    importInfo = {}
    importOrder = ['sourcefilename','annotation']
    if jobId is not None:
      order = copy.deepcopy(PROJECTSMANAGER().db().JOBITEMS)
      order.sort()
      order.extend(['relpath','projectid','projectname','parentjobnumber','childjobs'])
      info = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=order)
      title = 'job id: '+str(jobId) 
    else:
      order = copy.deepcopy(PROJECTSMANAGER().db().FILEITEMS)
      order.sort()
      order.extend(['jobnumber','taskname','projectname','projectid'])  
      info = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=order)
      title = 'file id: '+str(fileId)

      try:
        importInfo = PROJECTSMANAGER().db().getImportFileInfo(fileId=fileId,mode=importOrder)
      except:
        pass
    
    text = 'Database entry for '+title+'\n'
    for item in order:
      try:
        text = text + item+': '+str(info[item])+'\n'
      except:
        pass

    if len(importInfo)>0:
      for item in importOrder:
        try:
          text = text + item+': '+str(importInfo[item])+'\n'
        except:
          pass

    label.setPlainText(text)
    label.setReadOnly(True)
    label.setMinimumHeight(300)
    label.setMinimumWidth(500)
    win.show()
    win.raise_()
 
  @QtCore.Slot('QModelIndex')
  def showFileFromDirView(self,modelIndex):
    filePath = str(self.dirModel.filePath(self.proxyModel.mapToSource(modelIndex)))
    self.launchViewer(filePath)

  @QtCore.Slot(str,str,str,str)
  def exportData(self,myFileName,fileId,jobId,exportFileName):
    #print 'exportData',myFileName,fileId,jobId,exportFileName
    if fileId is not None or (jobId is not None and os.path.exists(os.path.join(PROJECTSMANAGER().db().jobDirectory(jobId=jobId),'hklout.mtz'))) :
      if os.path.splitext(exportFileName) !=  os.path.splitext(myFileName):
          exportFileName = os.path.splitext(exportFileName)[0] +  os.path.splitext(myFileName)[1]
      if fileId is None:
        myFileName = os.path.join(PROJECTSMANAGER().db().jobDirectory(jobId=jobId),'hklout.mtz')
      import shutil
      try:
        shutil.copyfile(myFileName,exportFileName)
      except:
        e = CException(self.__class__,100,'From: '+str(myFileName)+' to: '+str(exportFileName))
        e.warningMessage('Copying file',parent=self)
      else:
        PROJECTSMANAGER().db().createExportFile(fileId=fileId,exportFilename=exportFileName)
        fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['jobid','projectname'])
        from dbapi import CCP4DbUtils
        CCP4DbUtils.makeJobBackup(jobId=fileInfo['jobid'],projectName=fileInfo['projectname'])
    elif jobId is not None:
      # Use the mergeMtz plugin to merge all input and output data objects
      from core import CCP4TaskManager
      from dbapi import CCP4DbApi
      taskObj = CCP4TaskManager.TASKMANAGER().getPluginScriptClass('mergeMtz')(self)
      #print 'CProjectWidget.exportData taskObj',taskObj
      taskObj.container.outputData.HKLOUT.setFullPath(exportFileName)
      for role in [CCP4DbApi.FILE_ROLE_IN,CCP4DbApi.FILE_ROLE_OUT]:
        fileIdList = PROJECTSMANAGER().db().getJobFiles(jobId=jobId,role=role,fileTypes=CCP4DbApi.MINIMTZFILETYPES)
        #print 'CProjectWidget.exportData getJobFiles',role,fileIdList
        for fileId in fileIdList:
          name = PROJECTSMANAGER().db().getFullPath(fileId=fileId)
          taskObj.container.inputData.MINIMTZINLIST.append({'fileName' : name })
          taskObj.container.inputData.MINIMTZINLIST[-1].columnTag.set(PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode='JobParamName'))
          taskObj.container.inputData.MINIMTZINLIST[-1].setColumnNames(mode='applyTag')
      #print 'CProjectWidget.exportData MINIMTZINLIST',taskObj.container.inputData.MINIMTZINLIST
      taskObj.process()
      #print 'CProjectWidget.exportData',exportFileName,os.path.exists(exportFileName)

  @QtCore.Slot('QModelIndex')
  def handleFileClicked(self,modelIndex):
    #print("mapToSource 5",self.projectView.model(),modelIndex.model())
    node = self.projectView.model().sourceModel().nodeFromIndex(self.projectView.model().mapToSource(modelIndex))
    filePath = PROJECTSMANAGER().db().getFullPath(node.fileId)
    self.launchViewer(filePath)
    
  def handleJobClicked(self,modelIndex):
    node = self.projectView.model().sourceModel().nodeFromIndex(self.projectView.model().mapToSource(modelIndex))
    try:
      jobId = node.jobId
    except:
      return
    print(jobId)
    pipelineJobNode = node.getTopJob()
    self.currentJobChanged.emit(None,jobId,pipelineJobNode.jobId,False)

  @QtCore.Slot('QItemSelection','QItemSelection')
  def handleSelectionChanged(self,selected,deselected):
    indices = selected.indexes()
    #print 'CProjectWidget.handleSelectionChanged',indices
    if len(indices)>0:
      if indices[0].model() is self.projectView.model():
          self.handleJobClicked(indices[0])
      else:
          self.handleJobClicked(self.projectView.model().mapFromSource(indices[0]))
      

  @QtCore.Slot('QModelIndex')
  def handleDoubleClick(self,modelIndex):
    #print 'CProjectWidget.handleDoubleClick',double
    #print("mapToSource 7")
    node = self.projectView.model().sourceModel().nodeFromIndex(self.projectView.model().mapToSource(modelIndex))
    if node.isJob():
      fileId = None
      jobId = node.jobId
      pipelineJobNode = node.getTopJob()
    else:
      fileId = node.fileId
      jobId = node.parent().jobId
      pipelineJobNode = node.parent().getTopJob()
     
    if fileId is None:
      #if not double: self.projectView.model().setHighlightRow(currentCacheIndex)
      # Add state of 'double' which indicates job view should be detatched
      self.currentJobChanged.emit(fileId,jobId,pipelineJobNode.jobId,True)
      #elif double and fileId is not None:
    else:
      filePath = PROJECTSMANAGER().db().getFullPath(fileId)
      self.launchViewer(filePath)
    #print 'done handleClick'
        
  def getHeaderState(self):
    return str(self.projectView.header().saveState().data())

  def setHeaderState(self,state):
    stateByteArray = QtCore.QByteArray()
    stateByteArray.append(state)
    self.projectView.header().restoreState(stateByteArray)

  def selectJob(self,jobId,selectionFlags=QtCore.QItemSelectionModel.ClearAndSelect):
    #print 'CProjectWidget.selectJob',jobId
    modelIndex = self.projectView.model().sourceModel().modelIndexFromJob(jobId)
    if modelIndex is not None:
      sel = QtCore.QItemSelection( modelIndex.sibling(modelIndex.row(),0) , modelIndex.sibling(modelIndex.row(),CProjectModel.NCOLUMNS-1) )
      self.projectView.selectionModel().select(sel,selectionFlags)


    
class CFileIconProvider(QtWidgets.QFileIconProvider):

  def icon(self,fileInfo):
    return QtGui.QIcon()
    #print 'CFileIconProvider.icon',fileInfo
    if isinstance(fileInfo,QtCore.QFileInfo):
      suffix = str(fileInfo.completeSuffix())
      #print 'CFileIconProvider',suffix
      if suffix == 'mtz':
        content = fileInfo.baseName().__str__().split(CCP4File.CDataFile.SEPARATOR)[-1]
        #print 'CFileIconProvider content',content
        mimeType=MIMETYPESHANDLER().formatFromFileExt(ext=suffix,contentLabel=content)
      else:
        mimeType=MIMETYPESHANDLER().formatFromFileExt(ext=suffix)
      #print 'CFileIconProvider.icon',suffix,mimeType
      if mimeType is not None:
        icon = MIMETYPESHANDLER().icon(mimeType)
        if icon is not None:
          #print 'CFileIconProvider providing',icon
          return icon
    return QtGui.QIcon()
    
  
    
class CProjectDirProxyModel(QtCore.QSortFilterProxyModel):
    def lessThan(self,left,right):
        leftData = self.sourceModel().data(left,QtCore.Qt.DisplayRole)
        rightData = self.sourceModel().data(right,QtCore.Qt.DisplayRole)
        if len(leftData)> 4 and len(rightData) > 4 and leftData[:4] == "job_" and rightData[:4] == "job_": 
            try:
                return int(leftData[4:].split()[0]) < int(rightData[4:].split()[0])
            except:
                print("FAIL",leftData,rightData)
                raise
                return QtCore.QSortFilterProxyModel.lessThan(self,left,right)
        return QtCore.QSortFilterProxyModel.lessThan(self,left,right)

class CProjectDirModel(QtWidgets.QFileSystemModel):

  def __init__(self,parent,projectId):
    QtWidgets.QFileSystemModel.__init__(self,parent)
    self.projectId = projectId
    projectDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=self.projectId)
    self.setRootPath(projectDir)
    self.setIconProvider(CFileIconProvider())
    #self.sort(0,QtCore.Qt.DescendingOrder)
    lookup = PROJECTSMANAGER().db().getTaskNameLookup(projectId=self.projectId)
    self.jobNoList = []
    self.taskNameList = []
    for j,t in lookup:
      self.jobNoList.append(j)
      self.taskNameList.append(TASKMANAGER().getShortTitle(t))
    PROJECTSMANAGER().db().jobCreated.connect(self.updateLookup)
    PROJECTSMANAGER().db().jobFinished.connect(self.updateLookup1)
    PROJECTSMANAGER().db().jobDeleted.connect(self.updateLookup2)

  @QtCore.Slot(dict)
  def updateLookup(self,args):
    if  args['projectId'] != self.projectId: return
    info = PROJECTSMANAGER().db().getJobInfo(jobId=args['jobId'],mode=['jobnumber','taskname'])
    #print 'CProjectDirModel.updateLookup',info
    self.jobNoList.append(info['jobnumber'])
    self.taskNameList.append(TASKMANAGER().getShortTitle(info['taskname']))

  @QtCore.Slot(dict)
  def updateLookup1(self,args):
    if args.get('projectId') != self.projectId: return    
    info = PROJECTSMANAGER().db().getJobInfo(jobId=args.get('jobId'),mode=['jobnumber','taskname'])
    #print 'CProjectDirModel.updateLookup1',info
    if not info['jobnumber'] in self.jobNoList:
      self.jobNoList.append(info['jobnumber'])
      self.taskNameList.append(TASKMANAGER().getShortTitle(info['taskname']))

  @QtCore.Slot(dict)
  def updateLookup2(self,args):
    if args.get('projectId') != self.projectId: return
    if args['jobNumber'] in self.jobNoList:
      idx = self.jobNoList.index(args['jobNumber'])
      self.jobNoList.pop(idx)
      self.taskNameList.pop(idx)

    
  def data(self,modelIndex,role):
    ret = QtWidgets.QFileSystemModel.data(self,modelIndex,role)
    if role != QtCore.Qt.DisplayRole: return ret
    retStr = ret.__str__()
    if retStr[0:4] != 'job_': return ret
    
    userData = QtWidgets.QFileSystemModel.data(self,modelIndex,QtCore.Qt.UserRole)
    if not userData:
      try:
        pathSplit = (str(self.filePath(modelIndex))).split('/')
        jobNumber = pathSplit[pathSplit.index('CCP4_JOBS')+1][4:]
        for item in pathSplit[pathSplit.index('CCP4_JOBS')+2:]:
          jobNumber = jobNumber + '.' + item[4:]
        #print 'CProjectDirModel.data path',pathSplit,jobNumber
        if self.jobNoList.count(jobNumber)>0:
          taskName = self.taskNameList[self.jobNoList.index(jobNumber)]
          QtWidgets.QFileSystemModel.setData(self,modelIndex, taskName,QtCore.Qt.UserRole)
        else:
          return ret
      except:
        return ret
    else:
      taskName = userData.__str__()
    #print 'CProjectDirModel.data',retStr,taskName
    return retStr+' '+taskName

    
  def supportedDragActions(self):
    return QtCore.Qt.CopyAction

  def mimeTypes(self):
    typesList = []
    for item in list(MIMETYPESHANDLER().mimeTypes.keys()): typesList.append(item)
    #print 'CProjectDirModel.mimeTypes',typesList
    return typesList

  def projectRootIndex(self):
    import os
    return self.index(os.path.join(str(self.rootPath()),'CCP4_JOBS'),0)

  def jobNumberToModelIndex(self,jobNumber):
    jobNumber = str(jobNumber)
    if jobNumber[0] != 'j': jobNumber = 'job_'+jobNumber
    startIndex = self.index(0,0,self.projectRootIndex())
    modelIndexList = self.match(startIndex,QtCore.Qt.DisplayRole,jobNumber)
    #print('jobNumberToModelIndex',jobNumber,modelIndexList)
    if len(modelIndexList)>0:
      return modelIndexList[0]
    else:
      return QtCore.QModelIndex()

  '''
  def mimeData(self,modelIndexList):
    modelIndex = modelIndexList[0]
    qFileInfo = self.fileInfo(modelIndex)
    path = str(qFileInfo.absoluteFilePath())
    mimeType = MIMETYPESHANDLER().formatFromFileExt(fileName=path)
    if mimeType is None: return None
      
    from core import CCP4File
    fileObj = CCP4File.CDataFile(fullPath=path)
    dragText = fileObj.xmlText(pretty_print=False)
    print 'mimeData dragText',mimeType,dragText
    
    encodedData = QtCore.QByteArray()
    encodedData.append(dragText)
    mime = QtCore.QMimeData()
    # With mime type as text the data can be dropped on desktop
    # but the type of the data is lost
    #mimeData.setData('text/plain',data)
    mime.setData(mimeType,encodedData)
    return mime
  '''

class CProjectDirView(QtWidgets.QTreeView):
  
  def __init__(self,parent):
    QtWidgets.QTreeView.__init__(self,parent)
    self.setDragEnabled(True)
    self.setDragDropMode(QtWidgets.QAbstractItemView.DragOnly)
    self._initLayout= True
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)
    PREFERENCES().TABLES_ALTERNATING_COLOR.dataChanged.connect(self.resetAlternatingRowColors)

  @QtCore.Slot()
  def resetAlternatingRowColors(self):
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)
    

  def reset(self):
    QtWidgets.QTreeView.reset(self)
    if self._initLayout:
      self.setColumnWidth(0,300)
      self._initLayout=False

  def focusOn(self,jobNumber=None):
    modelIndex = self.model().sourceModel().jobNumberToModelIndex(jobNumber)
    if modelIndex.isValid():
      self.setExpanded(self.model().mapFromSource(modelIndex),True)
      self.scrollTo(self.model().mapFromSource(modelIndex),QtWidgets.QAbstractItemView.PositionAtTop)

  def startDrag(self,dropActions=None,modelIndex=None):
    if modelIndex is None:
      modelIndex = self.currentIndex()
    #print 'CProjectDirView.startDrag',dropActions,modelIndex,str(modelIndex.data(QtCore.Qt.DisplayRole))
    if modelIndex is None: return

    fileName = ''
    indx = modelIndex
    while indx.isValid():
      fileName = os.path.normpath(os.path.join(str(indx.data(QtCore.Qt.DisplayRole)),fileName))
      indx = indx.parent()
    if len(fileName)<=0: return

    fileName = fileName[0:-1]
    fileId,jobId = PROJECTSMANAGER().db().matchFileName(fileName=fileName)
    #print 'CProjectDirView.startDrag fileName',fileName,fileId
    if fileId is None: return
 
    from lxml import etree   
    mimeType = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode='fileclass')
    root,err = PROJECTSMANAGER().db().getFileEtree(fileId=fileId)
    dragText = etree.tostring(root)
    #print 'CProjectDirView.startDrag',mimeType,dragText
    
    encodedData = QtCore.QByteArray()
    encodedData.append(dragText)
    mimeData = QtCore.QMimeData()
    # With mime type as text the data can be dropped on desktop
    # but the type of the data is lost
    #mimeData.setData('text/plain',data)
    mimeData.setData(mimeType,encodedData)
    url = QtCore.QUrl()
    url.setPath(fileName)
    url.setScheme("file")
    mimeData.setUrls([url])
    
    drag = QtGui.QDrag(self)
    drag.setMimeData(mimeData)
    
    iconVar = modelIndex.data(QtCore.Qt.DecorationRole)
#FIXME - Almost certainly do not want toPyObject any more.
    #icon = iconVar.toPyObject()
    icon = iconVar
    pixmap = icon.pixmap(18,18)
    drag.setHotSpot(QtCore.QPoint(9,9))
    drag.setPixmap(pixmap)

    #print 'CProjectDirView.startDrag exec_'
    drag.exec_(QtCore.Qt.CopyAction)
      
    
    
class CEvaluationDelegate(QtWidgets.QItemDelegate):

  def __init__(self,parent=None):
    QtWidgets.QItemDelegate. __init__(self,parent)
    

  def createEditor(self,parent,option,modelIndex):
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    if 'evaluation' in self.parent().model()._dataCache[cacheIndex]:
      iconCombo = QtWidgets.QComboBox(parent)
      iconCombo.setEditable(False)
      for item in CCP4DbApi.JOB_EVALUATION_TEXT:
        iconCombo.addItem(jobIcon(item),item)
      iconCombo.show()
      return iconCombo
    else:
      return None
    
  def setEditorData(self,editor,modelIndex):
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    data = self.parent().model()._dataCache[cacheIndex]['evaluation']
    #print 'CEvaluationDelegate.setEditorData',data
    editor.setCurrentIndex(CCP4DbApi.JOB_EVALUATION_TEXT.index(data))

  def setModelData(self,editor,model,modelIndex):
    evaluation = str(editor.currentText())
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    self.parent().model()._dataCache[cacheIndex]['evaluation'] = evaluation
    jobid = self.parent().model()._dataCache[cacheIndex]['jobid']
    try:
      PROJECTSMANAGER().db().updateJob(jobid,key='evaluation',value=evaluation)
    except:
      pass

  def updateEditorGeometry(self,editor,option,modelIndex):
    r = option.rect
    r.setHeight(editor.geometry().height())
    r.setWidth(editor.geometry().width())
    editor.setGeometry(r)

class CJobEditLineEdit(QtWidgets.QLineEdit):

  losingFocus = QtCore.Signal()

  def focusOutEvent(self,event):
    #print 'CJobEditLineEdit.focusOutEvent'
    self.losingFocus.emit()
    QtWidgets.QLineEdit.focusOutEvent(self,event)
    
  
class CJobEditDelegate(QtWidgets.QStyledItemDelegate):

  def __init__(self,parent=None):
    QtWidgets.QStyledItemDelegate. __init__(self,parent)
    #print 'CJobEditDelegate.__init__'
    self.editorWidget = None

  def createEditor(self,parent,option,modelIndex):
    #if self.parent().model().data(modelIndex,QtCore.Qt.DisplayRole):
    #print 'CJobEditDelegate.createEditor'
    window = QtWidgets.QDialog(parent)
    window.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.Dialog)
    layout = QtWidgets.QVBoxLayout()
    layout.setMargin(CProjectWidget.MARGIN)
    layout.setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    window.setLayout(layout)
    editor = CJobEditLineEdit(window)
    editor.setMinimumWidth(350)
    editor.setObjectName('editor')
    window.layout().addWidget(editor)
    self.editorWidget = window
    self.editorWidget.setFocus()
    editor.losingFocus.connect(self.closeEditor)
    return window

  @QtCore.Slot()
  def closeEditor(self):
    #print 'CJobEditDelegate.closeEditor',self.editorWidget
    if self.editorWidget is not None:
      self.editorWidget.deleteLater()
      self.editorWidget = None
    
  def setEditorData(self,editor,modelIndex):
    data = self.parent().model().nodeFromIndex(modelIndex).getName(jobNumber=False)
    if data is not None:
      editor.findChild(QtWidgets.QLineEdit,'editor').setText(data)
  
  def setModelData(self,editor,model,modelIndex):
    #print 'CJobEditDelegate.setModelData'
    text = str(editor.findChild(QtWidgets.QLineEdit,'editor').text())
    #print("mapToSource 8")
    node = self.parent().model().sourceModel().nodeFromIndex(self.parent().model().mapToSource(modelIndex))
    if node.isJob():
      try:
        PROJECTSMANAGER().db().updateJob(node.getJobId(),key='jobtitle',value=text)
      except:
        pass
    elif node.isFile():
      try:
        PROJECTSMANAGER().db().updateFile(node.getFileId(),key='annotation',value=text)
      except:
        pass
 
  def updateEditorGeometry(self,editor,option,modelIndex):
    #Put the edit line in the right place
    editor.show()
    view = editor.parent().parent()
    # mapToGlobal() does not allow for the header - need to add header height
    # and want to not overlap the icon (and the job number for jobs)
    #print("mapToSource 9")
    node = self.parent().model().sourceModel().nodeFromIndex(self.parent().model().mapToSource(modelIndex))
    rightShift = 25 + node.isJob()*25
    pos=view.mapToGlobal(view.visualRect(modelIndex).topLeft())
    editor.move(pos+QtCore.QPoint(rightShift,20))


    
class CFollowFromDelegate(QtWidgets.QItemDelegate):

  def __init__(self,parent=None):
    QtWidgets.QItemDelegate. __init__(self,parent)
    

  def createEditor(self,parent,option,modelIndex):
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    if 'status' in self.parent().model()._dataCache[cacheIndex]:
      iconCombo = QtWidgets.QComboBox(parent)
      iconCombo.setEditable(False)
      for text,icon in [['Unused','greendot'],['Current','greenarrowsup']]:
        iconCombo.addItem(jobIcon(icon),text)
      iconCombo.show()
      return iconCombo
    else:
      return None
    
  def setEditorData(self,editor,modelIndex):
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    if self.parent().model()._followFromCacheId == cacheIndex:
      editor.setCurrentIndex(1)
    else:
      editor.setCurrentIndex(0)

  def setModelData(self,editor,model,modelIndex):
    status = str(editor.currentText())
    cacheIndex = self.parent().model().cacheFromModelIndex(modelIndex)
    if status == 'Current':
      if cacheIndex == self.parent().model()._followFromCacheId:
        return
      else:
        jobid = self.parent().model()._dataCache[cacheIndex]['jobid']
        try:
          PROJECTSMANAGER().db().setProjectFollowFromJobId(self.parent().model()._projectId,jobid)
        except:
          print('Error setting db setProjectFollowFromJobId')

  def updateEditorGeometry(self,editor,option,modelIndex):
    r = option.rect
    r.setHeight(editor.geometry().height())
    r.setWidth(editor.geometry().width())
    editor.setGeometry(r)


class CHistoryGui(QtWidgets.QFrame):
  def __init__(self,parent):
    QtWidgets.QFrame.__init__(self,parent)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setMargin(CProjectWidget.MARGIN)
    self.layout().setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    line = QtWidgets.QHBoxLayout()
    self.label = QtWidgets.QLabel('History label',self)
    self.label.setObjectName('jobListLabelRed')
    line.addWidget(self.label)
    self.clearButton = QtWidgets.QPushButton('Clear',self)
    line.addWidget(self.clearButton)
    self.layout().addLayout(line)
    self.sliderFrame = QtWidgets.QFrame()
    self.sliderFrame.setLayout(QtWidgets.QGridLayout())
    self.slider = QtWidgets.QSlider(self)
    self.slider.setOrientation(QtCore.Qt.Horizontal)
    self.slider.setTickInterval(1)
    self.slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
    self.sliderFrame.layout().addWidget(QtWidgets.QLabel('Show route to final job..',self),0,0,1,2)
    self.sliderFrame.layout().addWidget(self.slider,1,0,1,1)
    self.sliderFrame.layout().setMargin(CProjectWidget.MARGIN)
    self.sliderFrame.layout().setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    self.layout().addWidget(self.sliderFrame)
    

  def setSlider(self,optionList):
    self.sliderFrame.show()
    self.slider.setMinimum(1)
    self.slider.setMaximum(len(optionList))
    self.slider.setValue(1)
    # Unset previous layout and delete previous labels
    sliderLayoutItem = self.sliderFrame.layout().takeAt(self.sliderFrame.layout().indexOf(self.slider))
    item = self.sliderFrame.layout().takeAt(1)
    while ( item is not None):
      #print 'setSlider deleting',item
      item.widget().deleteLater()
      item = self.sliderFrame.layout().takeAt(1)
    # Put labels on slider
    for ix in range(len(optionList)-1):
      #print 'setSlider label',ix,optionList[ix]
      self.sliderFrame.layout().setColumnStretch(ix,1)
      self.sliderFrame.layout().addWidget(QtWidgets.QLabel(str(optionList[ix][-1]),self),2,ix)
    # The last element of optionList is a list of all successive jobs
    self.sliderFrame.layout().addWidget(QtWidgets.QLabel('All',self),2,ix+1)
    self.sliderFrame.layout().setColumnStretch(ix+1,0.5)
    self.sliderFrame.layout().addItem(sliderLayoutItem,1,0,1,len(optionList))
    
  def unsetSlider(self):
    self.sliderFrame.hide()
    
  def set(self,text=None):
    if text is None:
      self.hide()
    else:
      self.label.setText(text)
      self.show()

  def unSet(self):
    self.set()

class CJobSearchWidget(QtWidgets.QFrame):

  '''
  DATATYPES =  { 2 : 'Coordinate model data',
                 1 :'Model sequence',
                 11:'Observed intensities and structure factors',
                 12:'Phases',
                 13:'Map coefficients',
                 10:'FreeR flag' }
  '''
  DATATYPESORDER = [ 2,1,11,12,13,10]
  PARAMNAMETEXT = 'Choose a parameter name'
  MARGIN = 1

  def __init__(self,parent=None,projectName=''):
    from core import CCP4Annotation
    from qtgui import CCP4AnnotationWidgets
    QtWidgets.QFrame.__init__(self,parent=parent)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setMargin(CJobSearchWidget.MARGIN)
    self.layout().setContentsMargins(CJobSearchWidget.MARGIN,CJobSearchWidget.MARGIN,
                                CJobSearchWidget.MARGIN,CJobSearchWidget.MARGIN)

    line = QtWidgets.QHBoxLayout()
    self.layout().addLayout(line)
    line.addWidget(QtWidgets.QLabel('Show jobs',self))
    self.taskNameWidget = QtWidgets.QComboBox(self)
    line.addWidget(self.taskNameWidget)
    
    line.addWidget(QtWidgets.QLabel('run',self))
    self.dataRangeCheck = QtWidgets.QCheckBox(self)
    line.addWidget(self.dataRangeCheck)
    self.dateRangeModel=CCP4Annotation.CDateRange(parent=self)
    self.dateRangeWidget = CCP4AnnotationWidgets.CDateRangeView(parent=self,model=self.dateRangeModel)
    self.dateRangeWidget.layout().takeAt(0).widget().deleteLater()
    self.dateRangeWidget.setStyleSheet("QFrame { border : 0px };")
    line.addWidget(self.dateRangeWidget)

    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Title/annotation',self))
    self.textSearchModeWidget = QtWidgets.QComboBox(self)
    self.textSearchModeWidget.setEditable(False)
    for lab in [ 'is', 'starts with', 'contains' ]: 
      self.textSearchModeWidget.addItem(lab)
    self.textSearchModeWidget.setCurrentIndex(2)
    line.addWidget(self.textSearchModeWidget)
    self.textSearchWidget = QtWidgets.QLineEdit(self)
    line.addWidget(self.textSearchWidget)
    but = QtWidgets.QPushButton('Help',self)
    but.clicked.connect(functools.partial(self.showHelp,'text_search_syntax'))
    line.addWidget(but)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    self.modeButtonGroup = QtWidgets.QButtonGroup(self)
    self.modeButtonGroup.setExclusive(True)
    butId = 0
    for item in ['Data history','Control parameters','Imported file']:
      but = QtWidgets.QRadioButton(item,self)
      butId += 1
      self.modeButtonGroup.addButton(but,butId)
      line.addWidget(but)
    line.addStretch(2)
    self.modeButtonGroup.button(1).setChecked(True)
    self.modeButtonGroup.buttonClicked[int].connect(self.handleModeChange)
    self.layout().addLayout(line)

    self.modeStack = QtWidgets.QStackedLayout()
    self.layout().addLayout(self.modeStack)

    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QVBoxLayout())
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Data history of',self))
    self.datatypeHistoryWidget = QtWidgets.QComboBox(self)
    for dT in self.DATATYPESORDER:
      self.datatypeHistoryWidget.addItem(CCP4DbApi.FILETYPELIST[dT][2],dT)
    self.datatypeHistoryWidget.currentIndexChanged[int].connect(self.loadJobHistoryWidget)
    line.addWidget(self.datatypeHistoryWidget)
    line.addWidget(QtWidgets.QLabel('to/from job',self))
    self.historyJobWidget = QtWidgets.QComboBox(self)
    line.addWidget(self.historyJobWidget)
    line.addStretch(2)
    frame.layout().addLayout(line)
    self.modeStack.addWidget(frame)

    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QGridLayout())
    frame.layout().addWidget(QtWidgets.QLabel('Control parameter(s)',self),0,0)
    helpBut = QtWidgets.QPushButton('Help',self)
    helpBut.clicked.connect(functools.partial(self.showHelp,'control_params'))
    frame.layout().addWidget(helpBut,1,0)
    self.paramNameWidgets = []
    self.paramValueWidgets = []
    for row in [0,1,2]:
      self.paramNameWidgets.append(QtWidgets.QComboBox(self))
      self.paramNameWidgets[-1].addItem(self.PARAMNAMETEXT)
      self.paramNameWidgets[-1].setStyleSheet("QComboBox { combobox-popup: 0; }")
      self.paramNameWidgets[-1].setMaxVisibleItems(20)
      frame.layout().addWidget(self.paramNameWidgets[-1],row,1)
      frame.layout().addWidget(QtWidgets.QLabel('has value',self),row,2)
      self.paramValueWidgets.append(QtWidgets.QLineEdit(self))
      frame.layout().addWidget(self.paramValueWidgets[-1],row,3)
    self.modeStack.addWidget(frame)

    frame = QtWidgets.QFrame(self)
    frame.setLayout(QtWidgets.QVBoxLayout())
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Enter name or use browser to select imported file',self))
    frame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    from core import CCP4File
    from core import CCP4DataManager
    self.importedFileObj = CCP4File.CDataFile(parent=self)
    self.importedFileWidget = CCP4DataManager.DATAMANAGER().widget(model=self.importedFileObj,parentWidget=self)
    line.addWidget(self.importedFileWidget)
    frame.layout().addLayout(line)
    self.modeStack.addWidget(frame)
    
    self.loadTaskName()
    self.loadJobHistoryWidget(0)
    self.handleModeChange(1)
    #self.handleProjectModeChange(1)

  @QtCore.Slot(int)
  def handleModeChange(self,mode):
    #print 'handleModeChange',mode
    if self.modeStack.currentIndex() == 2:
      self.taskNameWidget.currentIndexChanged[int].disconnect(self.loadControlParameters)
    self.modeStack.setCurrentIndex(mode-1)
    if mode == 2:
      self.loadControlParameters()
      self.taskNameWidget.currentIndexChanged[int].connect(self.loadControlParameters)
    
    
  def handleProjectModeChange(self,mode):
    #print 'handleProjectModeChange',mode
    pass
    #self.projectSelectionFrame.setVisible((mode==3))

  @QtCore.Slot(int)
  def loadControlParameters(self,indx=None):
    if self.modeStack.currentIndex() != 1: return
    taskName = str(self.taskNameWidget.itemData(self.taskNameWidget.currentIndex()))
    #print 'loadControlParameters',taskName
    if taskName == 'anyTask':
      QtWidgets.QMessageBox.warning(self,self.parent().windowTitle(),'Select a task to load control parameters')
      return
    container = self.parent().loadDefFile(taskName)
    pList = container.controlParameters.dataOrder()
    pList.sort()
    for ii in range(len(self.paramNameWidgets)):
      self.paramNameWidgets[ii].clear()
      self.paramNameWidgets[ii].addItem(self.PARAMNAMETEXT)
      for item in pList:
        #dataType = container.controlParameters.find(item).PYTHONTYPE
        #print 'loadControlParameters',item,dataType
        subList =  container.controlParameters.__getattr__(item).dataOrder()
        if len(subList)==0:
          self.paramNameWidgets[ii].addItem(item)
        else:
          for subItem in subList:
            self.paramNameWidgets[ii].addItem(item+'.'+subItem)

    
  @QtCore.Slot(int)
  def loadJobHistoryWidget(self,indx):
    try:
      currentVar = self.historyJobWidget.itemData(self.historyJobWidget.currentIndex())
    except:
      currentVar = None
    #print 'loadJobHistoryWidget current',currentVar
    
    jobList = PROJECTSMANAGER().db().getProjectJobsByFileType(self.parent().parent().projectId,self.DATATYPESORDER[indx])
    self.historyJobWidget.clear()
    self.historyJobWidget.addItem('No job selected')
    for jNum,jId,tName,cTime in jobList:
      self.historyJobWidget.addItem(jNum+' '+TASKMANAGER().getShortTitle(tName),jId)
    if currentVar is not None:
      idx = max(0,self.historyJobWidget.findData(currentVar))
    else:
      idx = 0
    self.historyJobWidget.setCurrentIndex(idx)
    
  def loadTaskName(self):
    self.taskNameWidget.clear()
    self.projectInfo = PROJECTSMANAGER().db().getProjectJobSearchInfo(self.parent().parent().projectId)
    self.taskNameWidget.addItem('Any task','anyTask')
    for taskName in self.projectInfo['taskNameList']:
      self.taskNameWidget.addItem(TASKMANAGER().getShortTitle(taskName),taskName)
    self.taskNameWidget.setCurrentIndex(0)

  def clear(self):
    self.taskNameWidget.setCurrentIndex(0)
    self.textSearchWidget.clear()
    self.datatypeHistoryWidget.setCurrentIndex(0)
    for ii in range(len(self.paramNameWidgets)):
      self.paramNameWidgets[ii].setCurrentIndex(0)
      self.paramValueWidgets[ii].clear()

  def get(self):
    ret = {}
    ifText = False
    ret['taskName'] = str(self.taskNameWidget.itemData(self.taskNameWidget.currentIndex()))
    if ret['taskName'] == 'anyTask':
      ret['taskName'] = None
    else:
      ifText = True
      
    if self.dataRangeCheck.isChecked():
      self.dateRangeWidget.updateModelFromView()
      ret['minTime'],ret['maxTime'] = self.dateRangeModel.epochRange()
      ifText = True

    ret['searchText'] = str(self.textSearchWidget.text())
    if len(ret['searchText']) == 0:
      ret['searchText'] = None
    else:
      ifText = True

    ret['textSearchMode'] = self.textSearchModeWidget.currentIndex()
      
    ifHistory = False
    ifControlParams = False
    ifImportedFile = False
    ret['historyJob'] = None
    ret['controlValues'] = []
    ret['importedFile'] = None
    if self.modeStack.currentIndex() == 0:
      ret['historyFileType'] = self.DATATYPESORDER[self.datatypeHistoryWidget.currentIndex()]
      idx = self.historyJobWidget.currentIndex()
      if idx>0:
        ret['historyJob'] = str(self.historyJobWidget.itemData(idx))
        ifHistory = True
    elif self.modeStack.currentIndex() == 1:
      for ii in range(len(self.paramNameWidgets)):
        if str(self.paramNameWidgets[ii].currentText()) != self.PARAMNAMETEXT and \
                         len(str(self.paramValueWidgets[ii].text()).strip()) > 0:
          ifControlParams = True
          ret['controlValues'].append( [ str(self.paramNameWidgets[ii].currentText()),
                                 str(self.paramValueWidgets[ii].text()).strip() ] )
    elif self.modeStack.currentIndex() == 2:
      self.importedFileWidget.updateModelFromView()
      ret['importedFile'] = str(self.importedFileObj)
      ifImportedFile = len(ret['importedFile'] )>0
    
    return ifText,ifHistory,ifControlParams,ifImportedFile, ret

  
  @QtCore.Slot(str)
  def showHelp(self,target):
    WEBBROWSER().loadWebPage(helpFileName='searchTools',target=target)

class CJobSearchDialog(QtWidgets.QDialog):
  OPLIST = [ '>=', '<=', '>' , '<', '==','=' ]
  ERROR_CODES = { 101 : { 'description' : 'Input value not appropriate type' }
                  }
  def __init__(self,parent,projectName=''):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('Search jobs in project: '+projectName)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setMargin(CProjectWidget.MARGIN)
    self.layout().setContentsMargins(CProjectWidget.MARGIN,CProjectWidget.MARGIN,
                                CProjectWidget.MARGIN,CProjectWidget.MARGIN)
    self.widgets = []
    self.widgets.append(CJobSearchWidget(self,projectName))
    self.layout().addWidget(self.widgets[0])
    
    line = QtWidgets.QHBoxLayout()
    line.addStretch(1)
    buttonBox = QtWidgets.QDialogButtonBox(self)
    button = buttonBox.addButton('Search',QtWidgets.QDialogButtonBox.ApplyRole)
    button.clicked.connect(self.handleApply)
    button = buttonBox.addButton(QtWidgets.QDialogButtonBox.Close)
    button.clicked.connect(self.hide)
    button = buttonBox.addButton('Clear',QtWidgets.QDialogButtonBox.ActionRole)
    button.clicked.connect(self.handleClear)
    line.addWidget(buttonBox)
    line.addStretch(1)
    self.layout().addLayout(line)

  @QtCore.Slot(str)
  def handleApply(self):
    
    ifText,ifHistory,ifControlParams,ifImportedFile, searchParams = self.widgets[0].get()
    #print 'handleApply searchParams',ifText,ifHistory,ifImportedFile,searchParams

    if ifText:
      PROJECTSMANAGER().loadDbTaskTitles()
      try:
        jobList = PROJECTSMANAGER().db().projectJobSearch(self.parent().projectId,searchParams=searchParams)
      except CException as e:
        e.warningMessage('Searching database',parent=self)
        return
      except Exception as e:
        print(e)
        return
      numList = []
      for id,num in jobList: numList.append(num)
    else:
      jobList = []
      numList = []
      
    if ifHistory:
      self.parent().showHistory(searchParams['historyJob'],searchParams['historyFileType'],selection=numList)
      return
  
    if ifImportedFile:
      # returned subList: importId, fileId, checksum, annotation, jobId
      importedFileList = PROJECTSMANAGER().db().getImportedFile(sourceFileName=searchParams['importedFile'],projectId=self.parent().projectId)
      importFileJobList = []
      for iF in importedFileList:
        jidnum = [ iF[4],iF[5] ]
        if not importFileJobList.count(jidnum): importFileJobList.append(jidnum)

      numList = []
      if ifText:       
        for jidnum in importFileJobList:
          if jidnum in jobList: numList.append(jidnum[1])
      else:
        for jid,num in  importFileJobList: numList.append(num)
            
    elif ifControlParams:
      jobList = self.searchControlParams(jobList,searchParams)
      numList = []
      for id,num in jobList: numList.append(num)
              
    self.parent().setJobColour(numList)
    self.parent().projectView.update()
    self.parent().historyGui.set('All selected jobs')
    self.parent().historyGui.unsetSlider()
      

  @QtCore.Slot()
  def handleClear(self):
    for w in self.widgets: w.clear()
    self.parent().model().sourceModel().unsetJobColour()
    self.parent().projectView.update()
    self.parent().historyGui.unSet()
    
  def searchControlParams(self,inputJobList,searchParams):
    errReport = CErrorReport()
    container = self.loadDefFile(searchParams['taskName'])
    targetList = []
    for ii in range(len(searchParams['controlValues'])):
      paramName = searchParams['controlValues'][ii][0]
      if paramName.count('.'):
        paramName,paramName0 = paramName.split('.')
        dataType = container.controlParameters.find(paramName).__getattr__(paramName0).PYTHONTYPE
      else:
        paramName0 = None
        dataType = container.controlParameters.find(paramName).PYTHONTYPE
      words = searchParams['controlValues'][ii][1].split()
      #print 'searchControlParams',ii,dataType,words
      for word in words:
        done = False
        for op in self.OPLIST:
          if word.startswith(op):
            try:
              value = dataType(word[len(op):])
              targetList.append([paramName,paramName0,op,value])
            except Exception as e:
              errReport.append(self.__class__,101,str(word)+' Error:'+str(e))
            else:
              done = True
              break
        if not done:
          if word.count(','):
            wList = []
            for w in word.split(','):
              try:
                wList.append(dataType(w))
              except:
                errReport.append(self.__class__,101,str(w)+' in '+str(word)+' Error:'+str(e))
            targetList.append([paramName,paramName0,'in',wList])
          else:
            try:
              value = dataType(word)
            except:
              errReport.append(self.__class__,101,str(word)+' Error:'+str(e))
            else:
              targetList.append([paramName,paramName0,'=',value])
        
    #print 'searchControlParams targetList',targetList
    if len(targetList)==0: return []
    jobList = []
    projDir = os.path.split(PROJECTSMANAGER().db().jobDirectory(inputJobList[0][0]))[0]
    #print 'searchControlParams projDir',projDir
    for jid,num in inputJobList:
      ok = True
      try:
        body = self.loadParamFile(os.path.join(projDir,'job_'+str(num),'input_params.xml'))
      except:
        print('Failed loading data for job',num)
        ok = False
      else:
        for param,param0,op,target in targetList:
          try:
            ele = body.find('controlParameters').find(param)
            if param0 is not None: ele = ele.find(param0)
            value = dataType(ele.text)
          except:
            print('Failed finding',param,'in job jumber',num)
            ok = False
          else:
            if op in ['=','==']:
              if value != target: ok=False
            elif op == '!=':
              if value == target: ok=False
            elif op == '>':
              if value <= target: ok=False
            elif op == '<':
              if value >= target: ok=False
            elif op == '>=':
              if value < target: ok=False
            elif op == '<=':
              if value > target: ok=False
            elif op == 'in':
              if value not in target: ok=False
            else:
              ok = False
      if ok: jobList.append([jid,num])
    #print 'searchControlParams jobList',jobList
    return jobList
          

  def loadDefFile(self,taskName):
    defFile = TASKMANAGER().lookupDefFile(taskName)
    if defFile is None: defFile = TASKMANAGER().searchDefFile(taskName)
    if defFile is None: return None
    from core import CCP4Container
    container = CCP4Container.CContainer()
    container.loadContentsFromXml(defFile)
    return container

  def loadParamFile(self,fileName):
    from core import CCP4File
    xmlFile = CCP4File.CI2XmlDataFile(fullPath=fileName)
    body = xmlFile.getBodyEtree()
    return body
