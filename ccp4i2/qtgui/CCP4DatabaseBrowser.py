"""
     qtgui/CCP4DatabaseBrowser.py: CCP4 Gui Project
     Copyright (C) 2014 STFC

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

'''
Liz Potterton Mar 14 - Browse projects and files
'''

'''
This module provides a dialog box to select a file of a given fileType from any project.
It is normally used from CDataFileView.
It reimplements QTreeView and QAbstractItemModel and uses the same CTreeItem classes (to support
the abstract item model) as CCP4ProjectWidget.CProjectWidget.
It displays only the project heirarchy and a list of the files in the project (ie no jobs)
The files are only loaded into the tree widget when the user opens a project folder
There is no special handling of the current project so theer is a possiblity of selecting a file
in the current project.
'''

from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4Modules
from qtgui import CCP4ProjectWidget



class CDatabaseBrowser(QtWidgets.QDialog):

  fileSelected = QtCore.Signal(str)

  def __init__(self,parent=None,title='Open file from another project',fileType=None,projectId=None):
    QtWidgets.QDialog.__init__(self,parent=parent)
    self.setWindowTitle(title)
    self.viewer = CDatabaseBrowserView(self)
    self.viewer.setMinimumWidth(400)
    model = CDatabaseBrowserModel(self,fileType=fileType)
    model.loadModel()
    self.viewer.setModel(model)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().addWidget(self.viewer)

    self.viewer.clicked.connect(self.handleClicked)
    self.viewer.doubleClicked.connect(self.handleDoubleClicked)
    self.show()

  @QtCore.Slot('QModelIndex')
  def handleClicked(self,modelIndex):
    pass

  @QtCore.Slot('QModelIndex')
  def handleDoubleClicked(self,modelIndex):
    treeNode = modelIndex.internalPointer()
    #print 'handleDoubleClicked',treeNode,isinstance(treeNode,CCP4ProjectWidget.CTreeItemFile)
    if isinstance(treeNode,CCP4ProjectWidget.CTreeItemFile):
      self.fileSelected.emit(treeNode.fileId)
  


class CDatabaseBrowserView(QtWidgets.QTreeView):

    def __init__(self,parent=None):
      QtWidgets.QTreeView.__init__(self,parent=parent)

class CDatabaseBrowserModel(QtCore.QAbstractItemModel):

  def __init__(self,parent=None,fileType=None):
    QtCore.QAbstractItemModel.__init__(self,parent=parent)
    self.rootItem = CCP4ProjectWidget.CTreeItemProject(self,['',''])
    self.fileType = fileType
      
  def columnCount(self,parent):
    return 1

  def headerData(self,section,orientation,role):
    return 'Projects/Files'

  def childCount(self):
    return self.rootItem.childCount('projects')

  def rowCount(self,parent):
    if not parent.isValid():
      return self.childCount()
    else:
      return parent.internalPointer().childCount()

    return self.rootItem.childCount('projects')

  def child(self,row):
    return self.rootItem.child(row)

  def data(self, index, role):
    if not index.isValid():
      return None
    else:
      return index.internalPointer().data(index.column(),role)

  def index(self, row, column, parent):
    #print 'CDatabaseBrowserModel.index',row, column,parent.isValid(),parent.internalPointer()
    if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
      return QtCore.QModelIndex()
    treeNode = parent.internalPointer()
    if treeNode is None:
      childNode = self.rootItem.child(row)
    else:
      childNode = parent.internalPointer().child(row)
    return self.createIndex(row, column, childNode)

  def parent(self,index):
    if not index.isValid():
      return QtCore.QModelIndex()
    childItem = index.internalPointer()
    parentItem = childItem.parent()
    #print 'parent',childItem,parentItem
    if parentItem == self.rootItem:
      return QtCore.QModelIndex()
    else:
      return self.createIndex(parentItem.row(), 0, parentItem)

  def canFetchMore(self,parent):
    if not parent.isValid():
      return False
    else:
      return  parent.internalPointer().canFetchMore()

  def fetchMore(self,parent):
    # Load the files for any opened project in the tree model
    if not self.canFetchMore(parent): return
    parentNode = parent.internalPointer()
    projectId = parentNode.getProjectId()
    fileList = CCP4Modules.PROJECTSMANAGER().db().getProjectFiles(projectId=projectId,
                                                           fileType=self.fileType,topLevelOnly=True)
    #print 'CDatabaseBrowserModel.fetchMore fileList',fileList
    parentNode.clearFiles()
    for fileInfo in fileList:
      treeItem = CCP4ProjectWidget.CTreeItemFile(parent=parentNode,infoList=fileInfo,displayColumn=0,maxChar=200,displayJobNumber=True)
      parentNode.appendChildFile(treeItem)
    parentNode.childFilesLoaded = True
    

  def loadModel(self):
    # Get list of projects in reverse alphabetic order so can run through the list in reverse
    # (so can delete project from list once it is in the tree)
    projectList = CCP4Modules.PROJECTSMANAGER().db().getProjectDirectoryList(order='DESC')
    #print 'CDatabaseBrowserModel.loadModel projectList',projectList
    # This is a list of projects with [projectId,projectName,projectDir,parentProjectId] for each project
    lookup = {}
    previousProjectListLen = len(projectList)
    for indx in range(len(projectList)-1,-1,-1):
      if projectList[indx][3] is None:
        treeItem = CCP4ProjectWidget.CTreeItemProject(self.rootItem,[projectList[indx][0],projectList[indx][1]])
        treeItem.addDummyFile()
        self.rootItem.appendChildProject(treeItem)
        lookup[projectList[indx][0]] = treeItem
        del projectList[indx]

    #print 'loadModel after first loop len(projectList)', len(projectList),lookup

    while len(projectList)>0 and len(projectList)<previousProjectListLen:
      #print 'CDatabaseBrowserModel.loadModel projectList length',len(projectList)
      previousProjectListLen = len(projectList)
      for indx in range(len(projectList)-1,-1,-1):
        parent = lookup.get(projectList[indx][3],None)
        if parent is not None:
          treeItem = CCP4ProjectWidget.CTreeItemProject(parent,[projectList[indx][0],projectList[indx][1]])
          treeItem.addDummyFile()
          #print 'loadModel appending to',parent.getProjectName(),treeItem.getProjectName()
          parent.appendChildProject(treeItem)
          lookup[projectList[indx][0]] = treeItem
          del projectList[indx]
      
     

      
        
        
                             
