"""
    MakeProjectsAndDoLigandPipeline_gui.py: CCP4 GUI Project
    
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

from __future__ import print_function

from ccp4i2.baselayer import QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class MakeProjectsAndDoLigandPipeline_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'MakeProjectsAndDoLigandPipeline' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'developer_tools' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='Make projects and do ligand pipeline'
    TASKTITLE='Make projects and do ligand pipeline'
    DESCRIPTION = '''A task to generate projects and apply ligand pipeline for multiple datasets'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        from ccp4i2.baselayer import QtGui, QtWidgets,QtCore
        self.openFolder(folderFunction='inputData',followFrom=False)

        self.createLine(['subtitle','Root directory'])
        subFrame = self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Directory in which projects can be found','widget','ROOT_DIRECTORY' ] )
        self.container.inputData.ROOT_DIRECTORY.dataChanged.connect(self.rootDirChanged)
        self.closeSubFrame()
        
        self.createLine(['subtitle','Compound IDs - Smiles - and Relative paths'])
        subFrame = self.openSubFrame(frame=True)
        addButton = QtWidgets.QPushButton("Add",self)
        addButton.clicked.connect(self.on_add_button_clicked)
        
        subFrame.layout().addWidget(addButton)
        removeButton = QtWidgets.QPushButton("Remove",self)
        subFrame.layout().addWidget(removeButton)
        removeButton.clicked.connect(self.on_remove_button_clicked)
        
        self.editor = QtWidgets.QTableView(self)
        self.editor.horizontalHeader().setStretchLastSection(True)
        self.editor.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.editor.setModel(StartTableModel(self))
        
        subFrame.layout().addWidget(self.editor)
        self.closeSubFrame()
        
        self.openSubFrame(frame=True)
        self.createLine(['label','For rigid body refinement use','widget','PIPELINE'])
        self.closeSubFrame()
        
        self.createLine(['subtitle','Input coordinates'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input model','widget','XYZIN' ] )
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()

        self.createLine( [ 'subtitle', 'Free-R flags' ])
        self.openSubFrame(frame=True)
        self.createLine(['widget','-browseDb','True','FREERFLAG'])
        self.closeSubFrame()

        self.openFolder(folderFunction='inputData',followFrom=False,title="Lists as lists")
        self.createLine ( [ 'widget','PATH_LIST' ] )
        self.createLine ( [ 'widget','PROJECTNAME_LIST' ] )
        self.createLine ( [ 'widget','SMILES_LIST' ] )
    
    @QtCore.Slot()
    def on_add_button_clicked(self):
        self.editor.model().insertRows(self.editor.model().rowCount(),1, QtCore.QModelIndex())
        self.getWidget('PROJECTNAME_LIST').updateViewFromModel()
        self.getWidget('SMILES_LIST').updateViewFromModel()
        self.getWidget('PATH_LIST').updateViewFromModel()
        
        self.editor.repaint()
    
    @QtCore.Slot()
    def on_remove_button_clicked(self):
        selection = self.editor.selectionModel().selection()
        rowsToDelete = []
        for anIndex in selection.indexes():
            if not anIndex.row() in rowsToDelete: rowsToDelete.append(anIndex.row())
        
        self.container.inputData.PROJECTNAME_LIST.setQualifiers({'listMinLength':0})
        self.container.inputData.PATH_LIST.setQualifiers({'listMinLength':0})
        self.container.inputData.SMILES_LIST.setQualifiers({'listMinLength':0})
        
        invertedRowsToDelete = sorted(rowsToDelete,reverse=True)
        
        #reverse sort rows so that serial deletion of them does not complicate matters
        for row in invertedRowsToDelete:
            self.editor.model().deleteRows(row, 1, QtCore.QModelIndex())
        self.container.inputData.PROJECTNAME_LIST.setQualifiers({'listMinLength':1})
        self.container.inputData.PATH_LIST.setQualifiers({'listMinLength':1})
        self.container.inputData.SMILES_LIST.setQualifiers({'listMinLength':1})
        
        self.getWidget('PROJECTNAME_LIST').updateViewFromModel()
        self.getWidget('SMILES_LIST').updateViewFromModel()
        self.getWidget('PATH_LIST').updateViewFromModel()
        
        self.editor.repaint()

    @QtCore.Slot()
    def rootDirChanged(self):
        import fnmatch
        import os
        matches = []
        rootLength = len(self.container.inputData.ROOT_DIRECTORY.__str__())
        for root, dirnames, filenames in os.walk(self.container.inputData.ROOT_DIRECTORY.__str__()):
            for filename in fnmatch.filter(filenames, '*INTEGRATE.mtz'):
                if "dials-run" in root:
                    matches.append(os.path.join(root, filename)[rootLength:])
            for filename in fnmatch.filter(filenames, '*INTEGRATE.HKL'):
                if ("3d-run" in root ) or ("3dii-run" in root) or ("3daii-run" in root) or ("fast_dp" in root):
                    matches.append(os.path.join(root, filename)[rootLength:])
            for filename in fnmatch.filter(filenames, 'aimless_unmerged.mtz'):
                if "autoPROC" in root:
                    matches.append(os.path.join(root, filename)[rootLength:])
        
        #Clear existing table
        self.container.inputData.PROJECTNAME_LIST.setQualifiers({'listMinLength':0})
        self.container.inputData.PATH_LIST.setQualifiers({'listMinLength':0})
        self.container.inputData.SMILES_LIST.setQualifiers({'listMinLength':0})
        
        rowsToDelete = [i for i in range (len(self.container.inputData.PROJECTNAME_LIST))]
        
        #reverse sort rows so that serial deletion of them does not complicate matters
        invertedRowsToDelete = sorted(rowsToDelete,reverse=True)
        for row in invertedRowsToDelete:
            self.editor.model().deleteRows(row, 1, QtCore.QModelIndex())
        
        self.container.inputData.PROJECTNAME_LIST.setQualifiers({'listMinLength':1})
        self.container.inputData.PATH_LIST.setQualifiers({'listMinLength':1})
        self.container.inputData.SMILES_LIST.setQualifiers({'listMinLength':1})
        
        for match in matches:
            #Insert a row: this will generate default entried in each of the associated lists
            self.editor.model().insertRows(self.editor.model().rowCount(),1, QtCore.QModelIndex())
            newPath = self.container.inputData.PATH_LIST[-1] = match

        self.getWidget('PROJECTNAME_LIST').updateViewFromModel()
        self.getWidget('SMILES_LIST').updateViewFromModel()
        self.getWidget('PATH_LIST').updateViewFromModel()

        self.editor.repaint()

from ccp4i2.baselayer import QtCore, QtGui, QtWidgets
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class StartTableModel(QtCore.QAbstractTableModel):
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
    def __init__(self, forGui):
        super(StartTableModel,self).__init__()
        self.forGui = forGui
        self.headersAsStrings = ["Id","SMILES","RelativePath"]
        self.headersAsQStrings = self.headersAsStrings
        '''
            self.setHorizontalHeaderLabels(headersAsStrings)
        print(headersAsQStrings)
        for iString, headerQString in enumerate(headersAsQStrings):
            self.setHeaderData(iString, QtCore.Qt.Horizontal, headerQString)
        '''

#------------------------------------------------------------------------------------------------------
    def rowCount(self,index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
        return len(self.forGui.container.inputData.PATH_LIST)
  
#------------------------------------------------------------------------------------------------------
    def columnCount(self,index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
        return 3
  
#------------------------------------------------------------------------------------------------------
    def data(self, index, role = QtCore.Qt.DisplayRole):
#------------------------------------------------------------------------------------------------------
        text = "None"
        if index.column() == 0:
            if index.row() < len(self.forGui.container.inputData.PROJECTNAME_LIST):
                text = self.forGui.container.inputData.PROJECTNAME_LIST[index.row()].__str__()
            else:
                text = "Undefined"
        elif index.column() == 1:
            if index.row() < len(self.forGui.container.inputData.SMILES_LIST):
                text = self.forGui.container.inputData.SMILES_LIST[index.row()].__str__()
            else:
                text = "Undefined"
        elif index.column() == 2:
            import os
            if index.row() < len(self.forGui.container.inputData.PATH_LIST):
                text = self.forGui.container.inputData.PATH_LIST[index.row()].__str__()
                if role == QtCore.Qt.DisplayRole:
                    pathElements = text.split(os.path.sep)
                    try:
                        centralElement = pathElements.index("autoPROC")
                        firstElement = max(0,centralElement-1)
                        lastElement = min(len(pathElements), centralElement+2)
                        text = "::".join(pathElements[firstElement:lastElement])
                    except ValueError as err:
                        pass
                    try:
                        centralElement = pathElements.index("xia2")
                        firstElement = max(0,centralElement-1)
                        lastElement = min(len(pathElements), centralElement+2)
                        text = "::".join(pathElements[firstElement:lastElement])
                    except ValueError as err:
                        pass
                    try:
                        centralElement = pathElements.index("fast_dp")
                        firstElement = max(0,centralElement-1)
                        lastElement = min(len(pathElements), centralElement+1)
                        text = "::".join(pathElements[firstElement:lastElement])
                    except ValueError as err:
                        pass
            else:
                text = "Undefined"

        #Following line avoids the installation of a checkbox in every table cell
#FIXME PYQT - or maybe None? This used to return QVariant.
        if role == QtCore.Qt.CheckStateRole: return None
        if role == QtCore.Qt.BackgroundRole: return QtCore.Qt.red;

        return text

#------------------------------------------------------------------------------------------------------
    def insertRows(self, row, rows=1, index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
        self.beginInsertRows(QtCore.QModelIndex(), row, row + rows - 1)
        
        newPath = self.forGui.container.inputData.PATH_LIST.makeItem()
        #print type(newPath)
        self.forGui.container.inputData.PATH_LIST.append(newPath)
        
        newProjectName = self.forGui.container.inputData.PROJECTNAME_LIST.makeItem()
        self.forGui.container.inputData.PROJECTNAME_LIST.append(newProjectName)
        
        newSmiles = self.forGui.container.inputData.SMILES_LIST.makeItem()
        newSmiles.set("c1ccccc1Cl")
        self.forGui.container.inputData.SMILES_LIST.append(newSmiles)
        
        self.endInsertRows()
        return True

#------------------------------------------------------------------------------------------------------
    def deleteRows(self, row, rows=1, index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
        self.beginRemoveRows(QtCore.QModelIndex(), row, row + rows - 1)
        
        inputData = self.forGui.container.inputData
        
        projectNamesToRemove = []
        pathsToRemove = []
        smilesToRemove = []
        #Make a list of items that will be removed
        for i in range(row, row + rows):
            if i < len(inputData.PROJECTNAME_LIST):
                projectNamesToRemove.append(inputData.PROJECTNAME_LIST[i])
            if i < len(inputData.PATH_LIST):
                pathsToRemove.append(inputData.PATH_LIST[i])
            if i < len(inputData.SMILES_LIST):
                smilesToRemove.append(inputData.SMILES_LIST[i])
    
        for i in range(rows):
            try:
                inputData.PROJECTNAME_LIST.remove(projectNamesToRemove[i])
            except Exception as err:
                print('Error deleting projectName',err)
            try:
                inputData.PATH_LIST.remove(pathsToRemove[i])
            except Exception as err:
                print('Error deleting path',err)
            try:
                inputData.SMILES_LIST.remove(smilesToRemove[i])
            except Exception as err:
                print('Error deleting smiles',err)
        
        self.endRemoveRows()
        return True

#------------------------------------------------------------------------------------------------------
    def flags(self, index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
        #print "Evaluating flags for", str(index.column()), str(index.row())
        return QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled

#------------------------------------------------------------------------------------------------------
    def setData(self, index, variant, role = QtCore.Qt.EditRole):
#------------------------------------------------------------------------------------------------------
        if index.column() == 0: text = self.forGui.container.inputData.PROJECTNAME_LIST[index.row()]=str(variant.toPyObject())
        elif index.column() == 1: text = self.forGui.container.inputData.SMILES_LIST[index.row()]=str(variant.toPyObject())
        elif index.column() == 2: self.forGui.container.inputData.PATH_LIST[index.row()]=str(variant.toPyObject())
        self.dataChanged.emit(index, index)
        return True

#------------------------------------------------------------------------------------------------------
    def headerData(self, section, orientation, role):
#------------------------------------------------------------------------------------------------------
        wouldHaveBeen = super(StartTableModel,self).headerData(section, orientation, role)
        resultQStr=wouldHaveBeen
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            resultQStr = self.headersAsQStrings[section]
        return resultQStr

