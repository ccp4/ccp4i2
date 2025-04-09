"""
Copyright (C) 20011 University of York
Liz Potterton June 2011 - Viewer for file system based on QDirModel
"""
##@package CCP4WebView (QtWebKit) Web browser 'plugin' to view web pages

from PySide2 import QtCore, QtWidgets

from . import CCP4AbstractViewer


class CFileSystemView(CCP4AbstractViewer.CAbstractViewer):

  def __init__(self,parent,fileName=None):
    #print 'CFileSystemView.__init__',parent,fileName
    CCP4AbstractViewer.CAbstractViewer.__init__(self,parent)     
    self.setLayout(QtWidgets.QVBoxLayout())

    self.model = QtWidgets.QFileSystemModel(self)
    self.view = CFileSystemWidget(self)
    self.view.setModel(self.model)
    self.view.setColumnWidth(0,300)
    self.view.doubleClicked.connect(self.showFileFromView)
    self.layout().addWidget(self.view)

  def open(self,fileName):
    #print 'CFileSystemView.open',str(self.model.rootPath()),fileName,str(self.model.rootDirectory())
    modelIndex = self.model.setRootPath(fileName)
    self.view.setRootIndex(modelIndex)
    #self.view.model().reset()

  def handleRootChanged(self,newDir):
    print('handleRootChanged',str(newDir))


  @QtCore.Slot('QModelIndex')
  def showFileFromView(self,modelIndex):
    filePath = self.model.filePath(modelIndex)
    #print 'CCP4FileSystemView.showFileFromView',filePath
    from . import CCP4WebBrowser
    CCP4WebBrowser.OPENFILE(filePath)

  def title(self):
    return str(self.model.rootPath())
    
  def isPrintable(self):
    return 0

  def isSaveable(self):
    return 0

  def isScaleable(self):
    return 0


class CFileSystemWidget (QtWidgets.QTreeView):
  def __init__(self,parent=None):
    QtWidgets.QTreeView.__init__(self,parent=parent)
