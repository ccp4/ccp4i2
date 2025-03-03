from __future__ import print_function

"""
     CCP4FileSystemView.py: CCP4 GUI Project
     Copyright (C) 20011 University of York

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
     Liz Potterton June 2011 - Viewer for file system based on QDirModel
"""
##@package CCP4WebView (QtWebKit) Web browser 'plugin' to view web pages
from PySide2 import QtGui, QtWidgets,QtCore
from core.CCP4ErrorHandling import *
from qtgui import CCP4AbstractViewer


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
    from qtgui import CCP4WebBrowser
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
