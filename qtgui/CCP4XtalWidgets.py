from __future__ import print_function

"""
     qtgui/CCP4XtalWidgets.py: CCP4 Gui Project
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

##@package CCP4XtalWidgets (QtGui) Collection of widgets for crystallographic data types

import os
import functools

from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4XtalData
from qtgui import  CCP4Widgets
from core.CCP4Utils import safeFloat
from core.CCP4ErrorHandling import *
from core import CCP4Modules


class CSpaceGroupsAbstractItemModel(QtCore.QAbstractItemModel):

    def __init__(self,parent):
        QtCore.QAbstractItemModel.__init__(self,parent)
        chiralSpaceGroups= CCP4XtalData.SYMMETRYMANAGER().chiralSpaceGroups
        crystalSystems = CCP4XtalData.SYMMETRYMANAGER().crystalSystems
        model = []
        for key in crystalSystems:
            model.append((key,chiralSpaceGroups[key]))
        self.rootItem = CTreeItem(self,'root',model)

    def columnCount(self,parent):
        return 1

    def childCount(self):
        return self.rootItem.childCount()

    def rowCount(self,index):
        if not index.isValid():
            return self.childCount()
        else:
            return index.internalPointer().childCount()

    def child(self,row):
        return self.rootItem.child(row)

    def data(self, index, role):
        if not index.isValid():
            return None
        item = index.internalPointer()
        #return item.data(index.column(),role)
        return item.data(role)

    def index(self, row, column, parent):
        if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
            return QtCore.QModelIndex()
        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()
        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def parent(self,index):
        if not index.isValid():
            return QtCore.QModelIndex()
        childItem = index.internalPointer()
        parentItem = childItem.parent
        if parentItem == self.rootItem:
            return QtCore.QModelIndex()
        return self.createIndex(parentItem.row(), 0, parentItem)

class CSpaceGroupTreeView(QtWidgets.QTreeView):

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            index = self.indexAt(event.pos())
            #print 'CSpaceGroupTreeView.mousePressEvent', index.row(), index.parent().isValid(), self.parent().parent()
            if index.isValid() and not index.parent().isValid():
                #print 'CSpaceGroupTreeView.mousePressEvent expanding'
                self.setExpanded(index,(not self.isExpanded(index)))
                # Put temporary block on cling the popup menu
                self.parent().parent().blockClose = True
                event.accept()
                return

class CSpaceGroupCombo(QtWidgets.QComboBox):

    def __init__(self,parent):
        QtWidgets.QComboBox.__init__(self, parent)
        self.blockClose = False

    def hidePopup(self):
        #print 'CSpaceGroupCombo.hidePopup'
        if self.blockClose:
            self.blockClose = False
            return
        QtWidgets.QComboBox.hidePopup(self)


class CSpaceGroupView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CSpaceGroup

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self, parent=parent, qualifiers=qualis)
        wQualis = {}
        wQualis.update(qualis)
        wQualis['charWidth'] = 12
        self.layout().addWidget(QtWidgets.QLabel(qualifiers.get('label','Space group'),self))
        if self.editable:
            self.widget = CCP4Widgets.CSearchBox()
            self.widget.setHitList(CCP4XtalData.SYMMETRYMANAGER().nonEnantiogenicList())
            self.widget.changed.connect(self.updateModelFromView)
        else:
            self.widget = QtWidgets.QLabel(self)
        self.layout().addWidget(self.widget)
        self.setModel(model)
        if hasattr(self.widget,"dataSelected"):
            self.widget.dataSelected.connect(self.updateModelFromView)



class CCellView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CCell
    ERROR_CODES = {}

    STRETCH = 5

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = {'gridLayout' : True}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent, qualifiers=qualis)
        #print 'CCellView.__init__',self.layout()
        for row,rowItems in [[0, ['a', 'b', 'c']], [1, ['alpha', 'beta', 'gamma']]]:
            col = 0
            for item in rowItems:
                if item in ['alpha', 'beta', 'gamma']:
                    label = QtWidgets.QLabel('&'+item+';')
                else:
                    label = QtWidgets.QLabel(item)
                label.setTextFormat(QtCore.Qt.RichText)
                col = col+1
                self.layout().addWidget(label,row,col)
                if self.editable:
                    if model is not None:
                        wQualis = model.get(item).qualifiers()
                    else:
                        wQualis = {}
                    if model is not None:
                        wQualis['toolTip']=model.getDataObjects(item).qualifiers('toolTip')
                    self.widgets[item] = CCP4Widgets.CFloatView(self,qualifiers=wQualis)
                    self.widgets[item].editSignal.connect(self.updateModelFromView)
                    self.widgets[item].acceptDropDataSignal.connect(self.acceptDropData)
                    #w.setReadOnly(1-self.editable)
                else:
                    self.widgets[item] = CCP4Widgets.CLabel(self,qualis,uneditable=True)
                col = col+1
                self.layout().addWidget(self.widgets[item],row,col)
                #self.layout().setStretchFactor(self.widgets[item],2)
        if self.editable and model is not None:
            self.setModel(model)

    def getMenuDef(self):
        return ['copy','paste','help']

    def updateModelFromText(self):
        self.updateModelFromView()


class CSpaceGroupCellView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CSpaceGroupCell

    def __init__(self, parent=None, model=None, qualifiers={}):
        qualis = { 'gridLayout' : True,
                     'iconName' : 'cell' }
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self, parent, qualifiers=qualis)
        wQualis = { 'charWidth':12 }
        if qualis.get('editable',True):
            self.widgets['spaceGroup'] = CCP4Widgets.CSearchBox()
            self.widgets['spaceGroup'].setHitList(CCP4XtalData.SYMMETRYMANAGER().nonEnantiogenicList())
            self.widgets['spaceGroup'].changed.connect(self.updateModelFromView)
        else:
            self.widgets['spaceGroup'] = CCP4Widgets.CLabel(self)
        #self.widgets['spaceGroup'] = CSpaceGroupView(self,qualifiers=wQualis)
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(QtWidgets.QLabel('SpaceGroup',self))
        layout.addWidget(self.widgets['spaceGroup'])
        self.layout().addLayout(layout,1,1)
        if self.editable:
            self.loadButton = QtWidgets.QPushButton('Load from another file',self)
            self.loadButton.setToolTip('Take cell parameters from a different PDB or MTZ file')
            self.layout().addWidget(self.loadButton,0,1)
            self.loadButton.clicked.connect(self.handleLoadButton)
        self.widgets['cell'] = CCellView(parent=self.parent(),qualifiers=qualis)
        self.widgets['cell'].iconButton.deleteLater()
        self.layout().addWidget(self.widgets['cell'],0,2,2,1)
        self.setModel(model)

    @QtCore.Slot()
    def handleLoadButton(self):
        #print 'CSpaceGroupCellView.handleLoadButton'
        filterText = [ 'Coordinate file (*.pdb *.cif *.ent)' ,'Experimental data file (*.mtz *.cif *.ent)']
        from qtgui import CCP4FileBrowser
        self.fileBrowser = CCP4FileBrowser.CFileDialog(self, title='Select experimental data file or coordinate file',
                                                       filters = filterText, defaultFileName='')
        self.fileBrowser.setStyleSheet("")
        self.fileBrowser.setDownloadMode(modeList=['ebiPdb','ebiSFs'],projectId=self.parentTaskWidget().projectId())
        self.fileBrowser.selectFile.connect(self.loadFromFile)
        self.fileBrowser.show()

    @QtCore.Slot(str)
    def loadFromFile(self, fileName):
        ext = os.path.splitext(fileName)[1]
        self.model.blockSignals(True)
        if ext == '.pdb':
            from core import CCP4ModelData
            obj = CCP4ModelData.CPdbDataFile(fileName)
            #print 'CSpaceGroupCellView.loadFromFile',obj.fileContent.mmdbManager.GetCell(),obj.fileContent.mmdbManager.GetSpaceGroup()
            self.model.unSet()
            try:
                self.model.cell.set(obj.fileContent.mmdbManager.GetCell()[1:7] )
            except:
                pass
            try:
                self.model.spaceGroup.set(obj.fileContent.mmdbManager.GetSpaceGroup() )
            except:
                pass
        elif ext == '.mtz':
            obj = CCP4XtalData.CMtzDataFile(fileName)
            self.model.unSet()
            try:
                self.model.spaceGroup.set(obj.fileContent.spaceGroup)
            except:
                pass
            try:
                self.model.cell.set(obj.fileContent.cell)
            except:
                pass
        self.model.blockSignals(False)
        self.updateViewFromModel()
        self.validate()

    def dropTypes(self):
        return ['MtzDataFile', 'MiniMtzDataFile', 'ObsDataFile', 'PhsDataFile',
                'FreeRDataFile', 'MapCoeffsDataFile', 'PdbDataFile']

    def acceptDropData(self, textData):
        from lxml import etree
        tree = etree.fromstring(textData)
        try:
            path = os.path.join(tree.xpath('//relPath')[0].text.__str__(), tree.xpath('//baseName')[0].text.__str__())
        except:
            pass
        else:
            self.loadFromFile(path)

class CTreeItem:

    def __init__(self, parent=None, data={}, children=[]):
        #if len(data)>0:
        #  print 'CTreeItem.__init__',self,parent,str(data[QtCore.Qt.DisplayRole])
        self.parent = parent
        self.myData = {}
        if isinstance(data,dict):
            self.myData.update(data)
        else:
            self.myData[QtCore.Qt.DisplayRole] = data
        self.children = []
        for c in children:
            if isinstance(c,(list,tuple)):
                self.appendChild(CTreeItem(self, c[0], c[1]))
            else:
                self.appendChild(CTreeItem(self, c))

    def childCount(self):
        return len(self.children)

    def child(self,row):
        if row >= 0 and row < len(self.children):
            return self.children[row]
        else:
            return None

    def appendChild(self,item):
        self.children.append(item)

    def row(self):
        if self.parent is not None:
            return self.parent.children.index(self)
        else:
            return 0

    def data(self,role):
#FIXME PYQT ???? No idea what should be in place of QVariant() here.
        return self.myData.get(role,"")

    def findChild(self,role,value):
        for child in self.children:
            if role in child.myData and child.myData[role] == value:
                return child
        return None

class CProgramColumnGroupModel(QtCore.QAbstractItemModel):

    def __init__(self,parent=None,mtzData=None):
        QtCore.QAbstractItemModel.__init__(self,parent)
        self.clear()

    def clear(self):
        self.rootItem = CTreeItem()

    def setMtzData(self,mtzData=None, programGroupTypes=None, programColumnGroup=None):
        self.clear()
        datasetIndex = -1
        for dataset in mtzData.datasets:
            datasetIndex = datasetIndex + 1
            #print 'CProgramColumnGroupModel.setMtzData',dataset.name, datasetIndex
            columnGroupIndex = -1
            datasetListed = False
            for columnGroup in dataset.columnGroups:
                columnGroupIndex = columnGroupIndex + 1
                if self.isRightType(columnGroup,programGroupTypes,programColumnGroup):
                    if not datasetListed:
                        ds = CTreeItem (self.rootItem, {QtCore.Qt.DisplayRole: str(dataset.name) } )
                        self.rootItem.appendChild( ds )
                        datasetListed = True
                    ds.appendChild(CTreeItem(ds, { QtCore.Qt.DisplayRole : columnGroup.guiLabel(),
                                                QtCore.Qt.UserRole :  str(datasetIndex)+'.'+str(columnGroupIndex) }) )
        #print 'CProgramColumnGroupModel.setMtzData', self.rootItem.childCount(), self.nRows

    def isRightType(self, columnGroup, programGroupTypes, programColumnGroup):
        if len(programGroupTypes) > 0:
            return (str(columnGroup.groupType) in programGroupTypes)
        else:
            if len(columnGroup.columns) != len(programColumnGroup):
                return False
            for i in range(len(columnGroup.columns)):
                if columnGroup.columns[i].columnType not in programColumnGroup[i].columnType:
                    return False
            return True

    def flags(self,modelIndex):
        #if not modelIndex.parent().isValid() and not modelIndex.row() == 0:
        if not modelIndex.parent().isValid():
            return QtCore.Qt.ItemIsEnabled
        else:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def data(self,modelIndex,role):
        if not modelIndex.isValid():
            return None
        item = modelIndex.internalPointer()
        return item.data(role)

    def headerData(self, section, orientation, role):
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def rowCount(self, modelIndex):
        if modelIndex.column() > 0:
            return 0
        if not modelIndex.isValid():
            parentItem = self.rootItem
        else:
            parentItem = modelIndex.internalPointer()
        return parentItem.childCount()

    def columnCount(self,modelIndex):
        return 1

    def index(self, row, column, parent):
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()
        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()
        childItem = parentItem.child(row)
        if childItem is not None:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def parent(self, modelIndex):
        if not modelIndex.isValid():
            return QtCore.QModelIndex()
        childItem = modelIndex.internalPointer()
        parentItem = childItem.parent
        if parentItem is None:
            return None
        elif parentItem == self.rootItem:
            return QtCore.QModelIndex()
        else:
            return self.createIndex(parentItem.row(), 0, parentItem)


class CProgramColumnGroupQtView(QtWidgets.QTreeView):

    def __init__(self, parent=None):
        QtWidgets.QTreeView.__init__(self,parent)
        self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.header().setVisible(False)

class  CProgramColumnGroupCombo(QtWidgets.QComboBox):

    def __init__(self, parent):
        QtWidgets.QComboBox.__init__(self,parent)
        self.setModel(CProgramColumnGroupModel(self))
        self.setView(CProgramColumnGroupQtView(self))
        self.setEditable(False)

    def beep(self):
        print('CProgramColumnGroupCombo.beep', self.currentIndex(), self.view().selectionModel().currentIndex().row())

    def currentColumnGroupIndices(self):
        modelIndex = self.view().selectionModel().currentIndex()
        if modelIndex.isValid():
            data = modelIndex.internalPointer().data(QtCore.Qt.UserRole)
            if data is None:
                return -1, -1
            val = data.__str__().split('.')
            return int(val[0]), int(val[1])
        else:
            return -1,-1

    def setCurrentColumnGroup(self,dataset,columnGroupText):
        #print 'CProgramColumnGroupCombo.setCurrentColumnGroup',dataset,columnGroupText
        ds = self.model().rootItem.findChild(QtCore.Qt.DisplayRole,dataset)
        #print 'CProgramColumnGroupCombo.setCurrentColumnGroup ds',ds
        if ds is None: return
        cg = ds.findChild(QtCore.Qt.DisplayRole,columnGroupText)
        #print 'CProgramColumnGroupCombo.setCurrentColumnGroup cg',cg
        if cg is not None:
            dsi = self.model().index(ds.row(),0, QtCore.QModelIndex())
            cgi = self.model().index(cg.row(),0,dsi)
            self.setRootModelIndex( QtCore.QModelIndex())
            self.setCurrentIndex(cgi.row())


class CProgramColumnGroupView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CProgramColumnGroup
    MARGIN = 0
    STRETCH = 5
    ERROR_CODES = {101 : { 'description' : 'CProgramColumnGroup model not attached to a CMtzDataFile'},
                   102 : { 'description' : 'No CProgramColumnGroup model assigned for displaying MTZ column selection'},
                   103 : { 'description' : 'Error attempting to update display of MTZ columns'},
                   104 : { 'description' : 'Error number of MTZ columns in GUI does not match number in model'},
                   105 : { 'description' : 'MTZ file contents not available'}}

    def __init__(self, parent=None, model=None, qualifiers={},**kw):
        qualis = {}
        qualis.update(model.qualifiers())
        qualis.update(qualifiers)
        qualis.update(kw)
        qualis['gridLayout'] = True
        CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
        iconWidgetItem = self.layout().takeAt(0)
        iconWidgetItem.widget().deleteLater()
        self.widgets = []
        self.setModel(model)

    def setModel(self, model=None):
        if model == self.model:
            return
        if model is None or isinstance(model,self.MODEL_CLASS):
            if self.model is not None: CCP4Widgets.CViewWidget.unsetModel(self)
            self.model = model


    def drawColumns(self,editable=None):
        if self.model is None: raise CException(self.__class__,102,name=self.modelObjectPath())
        n = -1
        toolTip = self.model.qualifiers('toolTip')
        if editable is not None: self.editable = editable
        #print 'CProgramColumnGroupView.drawColumns columnGroupNames', self.model.columnGroupNames()
        for name in self.model.columnGroupNames():
            n = n + 1
            label = QtWidgets.QLabel(name,self)
            if self.editable:
                self.widgets.append( CCP4Widgets.CComboBox(self,qualifiers=  { 'dragType' : self.dragType() }) )
                self.widgets[-1].setMinimumContentsLength(25)
                self.widgets[-1].setEditable(False)
                if toolTip is not NotImplemented: self.widgets[-1].setToolTip(toolTip)
            else:
                self.widgets.append(CCP4Widgets.CLabel(self , { 'charWidth' : 25, 'editable' :False,  'dragType' : self.dragType()  } ))
            row = n/2
            column = (n%2) * 2
            self.layout().addWidget(label,row,column+1)
            self.layout().addWidget(self.widgets[n],row,column+2)
            if self.editable:
                self.widgets[-1].acceptDropDataSignal.connect(self.acceptDropData)
                # Just attempt to update all columns in sync with the first column
                if n == 0:
                    self.widgets[-1].currentIndexChanged[int].connect(self.updateFromFirstColumn)
                else:
                    self.widgets[-1].currentIndexChanged[int].connect(self.updateModelFromView)
        if n%2 == 1:
            pass


    def populateColumns(self, listOfColumnGroups=[]):
        for widget in self.widgets:
            widget.blockSignals(True)
            widget.clear()
        if self.model is None:
            for widget in self.widgets:
                widget.blockSignals(False)
            return
        if len(listOfColumnGroups)==1:
            #load into a label widget
            for n in range(0,self.model.nColumns()):
                label = str(listOfColumnGroups[0].columnList[n].dataset)+'/'+str(listOfColumnGroups[0].columnList[n].columnLabel)
                self.widgets[n].setText(label)
        else:
            for n in range(0,self.model.nColumns()):
                for item in listOfColumnGroups:
                    label = str(item.columnList[n].dataset)+'/'+str(item.columnList[n].columnLabel)
                    self.widgets[n].addItem(label,label)
        for widget in self.widgets:
            widget.blockSignals(False)
        self.updateFromFirstColumn()

    @QtCore.Slot()
    def updateFromFirstColumn(self):
        if self.model is None:
            raise CException(self.__class__, 102, name=self.modelObjectPath())
        data = {}
        nameList = self.model.columnGroupNames()
        #print 'updateFromFirstColumn',self.model,type(self.model),nameList
        if self.editable:
            data[nameList[0]] = str(self.widgets[0].currentText())
        else:
            data[nameList[0]] = str(self.widgets[0].text())
        for name in nameList[1:]:
            data[name] = None
            #self.model.__setattr__(name,columnLabel)
        #print 'CProgramColumnGroupView.updateFromFirstColumn',self.model.objectName(),data
        self.model.set(data,fix=True)
        self.updateViewFromModel()

    def updateViewFromModel(self):
        if self.model is None:
            raise CException(self.__class__,102,name=self.modelObjectPath())
        #print 'CProgramColumnGroupView.updateViewFromModel',self.model.objectName(),self.model
        if len(self.widgets)<self.model.nColumns(): raise CException(self.__class__,103,name=self.model.objectName())
        columnGroupItems = self.model.columnGroupNames()
        data = self.model.get('guiLabels')
        #print 'CProgramColumnGroupView.updateViewFromModel',self.model.objectName(),data
        n = -1
        for name in columnGroupItems:
            n = n + 1
            self.widgets[n].blockSignals(True)
            #print 'CProgramColumnGroupView.updateViewFromModel',name,n, data.get(name,None),'widget:',self.widgets[n]
            if data.get(name,None) is not None:
                if self.editable:
                    indx = self.widgets[n].findText(data[name])
                    if indx>=0: self.widgets[n].setCurrentIndex(indx)
                else:
                    #print 'CProgramColumnGroupView.updateViewFromModel',n,name,data[name]
                    self.widgets[n].setValue(data[name])
            else:
                if self.editable:
                    self.widgets[n].setCurrentIndex(self.widgets[n].count()-1)
                else:
                    self.widgets[n].setValue('')
            self.widgets[n].blockSignals(False)

    def updateModelFromView(self, **kw):
        if self.model is None:
            raise CException(self.__class__, 102, name=self.modelObjectPath())
        if len(self.widgets) < self.model.nColumns():
            raise CException(self.__class__, 104, name=self.modelObjectPath())
        n = -1
        data = {}
        for name in self.model.columnGroupNames():
            n = n + 1
            if self.editable:
                columnLabel = str(self.widgets[n].currentText())
            else:
                columnLabel = str(self.widgets[n].text())
            if columnLabel.count('/'):
                data[name] = columnLabel.split('/')[-1]
            else:
                data[name] = columnLabel
            #self.model.__setattr__(name,columnLabel)
        #print 'CProgramColumnGroupView.updateModelFromView',self.model.objectName(),data
        self.connectUpdateViewFromModel(False)
        self.model.set(data,fix=True)
        self.connectUpdateViewFromModel(True)
        self.validate()

    def getMenuDef(self):
        return ['copy', 'paste', 'help']

    def validate(self,**kw):
        # Dummy routine for widget that should now be redundant
        return True


class CDatasetView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CDataset

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {'vboxLayout' : True, 'iconName' : 'ObsDataFile', 'dragType' : 'ObsDataFile'}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        self.widgets = {}
        if model is None:
            selected = None
            obsDataFile = None
            crystalName = None
            datasetName = None
            formFactors = None
            formFactorSource = None
        else:
            selected = model.selected
            obsDataFile = model.obsDataFile
            crystalName = model.crystalName
            datasetName = model.datasetName
            formFactors = model.formFactors
            formFactorSource = model.formFactorSource
        self.widgets['selected']  = CCP4Widgets.CBooleanView(parent=self,model=selected)
        self.widgets['obsDataFile']  =CMiniMtzDataFileView(parent=self,model=obsDataFile,qualifiers={'jobCombo' : True,'dragType':'ObsDataFile'})
        self.widgets['crystalName'] = CCP4Widgets.CStringView(parent=self,model=crystalName)
        self.widgets['datasetName'] = CCP4Widgets.CStringView(parent=self,model=datasetName)
        self.widgets['formFactors'] = CFormFactorView(parent=self,model=formFactors)
        self.widgets['formFactorSource'] = CCP4Widgets.CStringView(parent=self,model=formFactorSource)
        self.layout().removeWidget(self.iconButton)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.iconButton)
        line.addWidget(self.widgets['obsDataFile'])
        self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.widgets['selected'])
        line.addWidget(QtWidgets.QLabel('Use this dataset.',self))
        line.addWidget(QtWidgets.QLabel('Crystal name',self))
        line.addWidget(self.widgets['crystalName'])
        line.addWidget(QtWidgets.QLabel('dataset name',self))
        line.addWidget(self.widgets['datasetName'])
        line.addStretch(1)
        self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.widgets['formFactors'])
        line.addWidget(QtWidgets.QLabel('set from (not working yet)',self))
        line.addWidget(self.widgets['formFactorSource'])
        line.addStretch(1)
        self.layout().addLayout(line)

class CDatasetListView(CCP4Widgets.CListView):

    MODEL_CLASS = CCP4XtalData.CDatasetList

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis =  { 'mode' : 'table',
                    'tableItems' : ['file','crystalName','datasetName','formFactors'],
                    'columnHeaders':['Observed data','Crystal','Dataset','Form factors'],
                   }
        qualis.update(qualifiers)
        CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)


class CAsuComponentView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CAsuComponent

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {'vboxLayout' : True, 'iconName' : 'SeqDataFile', 'dragType' : 'SeqDataFile'}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        self.widgets = {}
        #self.widgets['seqfile'] = CCP4Widgets.CSeqDataFileView(self,model=model.seqfile,qualifiers={'jobCombo' : False})
        if model is None:
            seqObj = None
            copiesObj = None
        else:
            seqObj = model.seqFile
            copiesObj = model.numberOfCopies
        self.widgets['seqFile']  = CCP4Widgets.CDataFileView(parent=self,model=seqObj,qualifiers={'jobCombo' : True,'dragType':'SeqDataFile'})
        self.widgets['numberOfCopies'] = CCP4Widgets.CIntView(parent=self,model=copiesObj,qualifiers= { 'editable' : self.editable, 'guiMode' : 'combo', 'enumerators' : [i for i in range(21)], 'menuText' : [str(i) for i in range(21)] })
        self.layout().removeWidget(self.iconButton)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.iconButton)
        #MN qualifiers dictionary provided has guiLabel efined as NotImplemented...intercept and provide value
        guiLabel = qualifiers.get('guiLabel','Number of copies in asymmetric unit:')
        if guiLabel is NotImplemented: guiLabel = 'Number of copies in asymmetric unit:'
        line.addWidget(QtWidgets.QLabel(guiLabel,self))
        line.addWidget(self.widgets['numberOfCopies'])
        line.addWidget(QtWidgets.QLabel('of sequence:',self))
        line.addStretch(2)
        self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.widgets['seqFile'])
        self.layout().addLayout(line)
        self.setModel(model)

    def setModel(self,model):
        #print 'CAsuComponentView.setModel',model,self.model
        if self.model is not None:
            for item in ['seqFile','numberOfCopies']:
                self.widgets[item].connectUpdateViewFromModel(False)
        self.model = model
        if model is None:
            for item in ['seqFile','numberOfCopies']:
                self.widgets[item].setModel(None)
                self.widgets[item].setValue(None)
        else:
            for item in ['seqFile','numberOfCopies']:
                #print 'CAsuComponentView.setModel updating',item,type(self.widgets[item])
                self.widgets[item].setModel(self.model.get(item))
                self.widgets[item].updateViewFromModel()
                self.widgets[item].connectUpdateViewFromModel(True)

    def acceptDropData(self,textData):
        from lxml import etree
        tree = etree.fromstring(textData)
        #print 'CAsuComponentView.acceptDropData',textData,tree.tag
        self.widgets['seqFile'].connectUpdateViewFromModel(False)
        if tree.tag in ['SeqDataFile','seqFile']:
            self.model.seqFile.setEtree(tree)
        else:
            self.model.setEtree(tree)
        self.widgets['seqFile'].connectUpdateViewFromModel(True)
        self.updateViewFromModel()
        self.widgets['seqFile'].validate()

    '''
    def createMimeData(self):
        return self.widgets['seqFile'].createMimeData()
    '''

class CFormFactorView(CCP4Widgets.CComplexLineWidget):
    MODEL_CLASS = CCP4XtalData.CFormFactor

    def __init__(self,parent=None,model=None,qualifiers={}):
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualifiers)
        self.layout().takeAt(0)
        self.widgets['Fp'] = CCP4Widgets.CFloatView(self)
        self.widgets['Fpp'] = CCP4Widgets.CFloatView(self)
        if model is not None:
            print('CFormFactorView.__init__',model)
            self.widgets['Fp'].setModel(model.Fp)
            self.widgets['Fpp'].setModel(model.Fpp)
        self.layout().addWidget(QtWidgets.QLabel("Form factor f'",self))
        self.layout().addWidget(self.widgets['Fp'])
        self.layout().addWidget(QtWidgets.QLabel("f''",self))
        self.layout().addWidget(self.widgets['Fpp'])


class CImportUnmergedView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CImportUnmerged

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = { 'vboxLayout' : True, 'iconName' : 'UnmergedDataFile' , 'dragType' :'UnmergedDataFile','dropTypes' : ['UnmergedDataFile','ImportUnmerged'] }
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        self.widgets['crystalName'] = CCP4Widgets.CStringView(parent=self,qualifiers={})
        self.widgets['crystalName'].setToolTip('Short name for the crystal used to identify it throughout project')
        self.widgets['dataset'] = CCP4Widgets.CStringView(parent=self,qualifiers={})
        self.widgets['dataset'].setToolTip('Short name for the dataset used to identify it throughout project')
        self.widgets['file'] =CCP4Widgets.CDataFileView(parent=self,qualifiers={ 'ifInfo': True, 'mustExist' : True })
        self.widgets['cell'] = CCellView(parent=self)
        self.widgets['cell'].setToolTip('Enter the crystal cell that is not provided in the data file')
        self.widgets['wavelength'] = CCP4Widgets.CFloatView(parent=self,qualifiers={'charWidth' : 10 })
        self.widgets['wavelength'].setToolTip('Enter the wavelength that is not provided in the data file')
        self.widgets['excludeSelection'] = CCP4Widgets.CStringView(parent=self,qualifiers={})
        self.widgets['excludeSelection'].setToolTip("Enter lists and/or ranges of batches in form '7,9,12-13,21-23'")
        self.widgets['batchsInFile'] = CCP4Widgets.CLabel(self)
        self.widgets['batchsInFile'].setMinimumWidth(100)
        self.sameAsCombo = QtWidgets.QComboBox(self)
        self.sameAsCombo.currentIndexChanged[int].connect(self.handleSameAsCombo)
        line = QtWidgets.QHBoxLayout()
        iconWidgetItem = self.layout().takeAt(0)
        line.addWidget(iconWidgetItem.widget())
        line.addWidget(self.widgets['file'])
        self.layout().addLayout(line)
        self.cellFrame = QtWidgets.QFrame(self)
        self.cellFrame.setLayout(QtWidgets.QVBoxLayout())
        self.cellFrame.layout().setContentsMargins(0,0,0,0)
        self.cellFrame.layout().setSpacing(0)
        #self.cellFrame.setFrameShape(QtWidgets.QFrame.Box)
        self.cellFrame.layout().addWidget(QtWidgets.QLabel('Cell parameters'))
        self.cellFrame.layout().addWidget(self.widgets['cell'])
        self.layout().addWidget(self.cellFrame )
        self.cellFrame.hide()
        self.wLFrame =  QtWidgets.QFrame(self)
        self.wLFrame.setLayout(QtWidgets.QHBoxLayout())
        self.wLFrame.layout().setContentsMargins(0,0,0,0)
        self.wLFrame.layout().setSpacing(0)
        self.wLFrame.layout().addWidget(QtWidgets.QLabel('Wavelength'))
        self.wLFrame.layout().addWidget(self.widgets['wavelength'])
        self.wLFrame.layout().addStretch(2)
        self.layout().addWidget(self.wLFrame )
        self.wLFrame.hide()
        #line = QtWidgets.QHBoxLayout()
        #line.addWidget(QtWidgets.QLabel(
        #  'Assign crystal and dataset names which will be used to identify data in subsequent tasks',self))
        #self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Crystal name',self))
        line.addWidget(self.widgets['crystalName'])
        line.addWidget(QtWidgets.QLabel('dataset name',self))
        line.addWidget(self.widgets['dataset'])
        self.sameAsFrame = QtWidgets.QFrame(self)
        self.sameAsFrame.setLayout(QtWidgets.QHBoxLayout())
        self.sameAsFrame.layout().setContentsMargins(0,0,0,0)
        self.sameAsFrame.layout().setSpacing(0)
        self.sameAsFrame.layout().addWidget(QtWidgets.QLabel('OR same dataset as',self))
        self.sameAsFrame.layout().addWidget(self.sameAsCombo)
        line.addWidget(self.sameAsFrame)
        line.addStretch(1)
        self.layout().addLayout(line)
        self.batchFrame =  QtWidgets.QFrame(self)
        self.batchFrame.setLayout(QtWidgets.QHBoxLayout())
        self.batchFrame.layout().setContentsMargins(0,0,0,0)
        self.batchFrame.layout().setSpacing(1)
        label = QtWidgets.QLabel('Batches in file:', self)
        italicFont = label.font()
        italicFont.setItalic(True)
        label.setFont(italicFont)
        self.batchFrame.layout().addWidget(label)
        self.batchFrame.layout().addWidget(self.widgets['batchsInFile'])
        self.layout().addWidget(self.batchFrame )
        self.batchFrame.hide()
        self.batchFrame2 =  QtWidgets.QFrame(self)
        self.batchFrame2.setLayout(QtWidgets.QHBoxLayout())
        self.batchFrame2.layout().setContentsMargins(0,0,0,0)
        self.batchFrame2.layout().setSpacing(1)
        label = QtWidgets.QLabel('Exclude batches from calculations and output', self)
        italicFont = label.font()
        italicFont.setItalic(True)
        label.setFont(italicFont)
        self.batchFrame2.layout().addWidget(label)
        self.batchFrame2.layout().addWidget(self.widgets['excludeSelection'])
        self.batchFrame2.layout().setStretchFactor(self.widgets['excludeSelection'],2)
        self.layout().addWidget(self.batchFrame2 )
        self.batchFrame2.hide()
        self.setModel(model)

    def setModel(self,model):
        #print 'CImportUnmergedView.setModel',repr(model),model
        itemList =  ['crystalName','dataset','file','excludeSelection','cell','wavelength']
        #itemList =  ['file','cell']
        if self.model is not None:
            for item in itemList:
                self.widgets[item].connectUpdateViewFromModel(False)
            """
            I am removing all these, since I see no corresponding connect statements.
            self.model.dataChanged.disconnect(self.validateNames)
            self.model.file.dataChanged.disconnect(self.model.file.loadFile)
            self.model.file.dataChanged.disconnect(self.model.loadDatasetName)
            self.model.file.fileContent.knowncell.dataChanged.disconnect(self.toggleCellView)
            self.model.file.fileContent.knownwavelength.dataChanged.disconnect(self.toggleWavelengthView)
            self.model.file.fileContent.batchs.dataChanged.disconnect(self.toggleBatchView)
            """
        if model is None:
            for item in  itemList:
                self.widgets[item].setModel(None)
                self.widgets[item].setValue(None)
        else:
            #data = model.get0()
            for item in itemList:
                self.widgets[item].setModel(model.get(item))
            #self.widgets['cell'].setModel(model.file.fileContent.cell)
            for item in itemList:
                self.widgets[item].updateViewFromModel()
                self.widgets[item].connectUpdateViewFromModel(True)
            self.sameAsFrame.setVisible(model.parent().index(model)!=0)
        self.model = model
        if self.model is not None:
            self.loadSameAsCombo()
            self.model.file.dataChanged.connect(self.handleFileChanged)
            self.toggleCellView()
            self.toggleWavelengthView()
            self.toggleBatchView()

    @QtCore.Slot()
    def handleFileChanged(self):
        self.model.blockSignals(True)
        self.model.file.loadFile()
        self.model.loadDatasetName()
        self.toggleCellView()
        self.toggleWavelengthView()
        self.toggleBatchView()
        self.updateViewFromModel()
        self.model.blockSignals(False)
#SJM 15/2/2022 - I put this in to force revalidation of the whole unmerged data file view widget and to ensure that crystal and data names are updated if there is a table. 
        self.model.dataChanged.emit()

    def validateNames(self):
        self.widgets['crystalName'].validate()
        self.widgets['dataset'].validate()

    def validate(self,isValid=None,reportMessage=True):
        excludeWidgets= []
        #print 'CImportUnmergedView.validate',isValid,self.model.file.fileContent.knowncell
        if self.model.file.fileContent.knowncell:
            excludeWidgets.append('cell')
            self.widgets['cell'].validate(True)
        if self.model.file.fileContent.knownwavelength:
            excludeWidgets.append('wavelength')
            self.widgets['wavelength'].validate(True)
        CCP4Widgets.CViewWidget.validate(self,isValid,excludeWidgets=excludeWidgets,reportMessage=reportMessage)

    def toggleCellView(self):
        if self.model.file.fileContent.knowncell.isSet() and not(self.model.file.fileContent.knowncell):
            self.cellFrame.show()
        else:
            self.cellFrame.hide()

    def toggleWavelengthView(self):
        if self.model.file.fileContent.knownwavelength.isSet() and not(self.model.file.fileContent.knownwavelength):
            self.wLFrame.show()
        else:
            self.wLFrame.hide()

    def hasBatchInfo(self):
        # return True if the file has batch information
        # format = ['unk','mtz' ,'xds', 'sca', 'saint', 'shelx']
        if ((self.model.file.fileContent.format == 'xds') or
            (self.model.file.fileContent.format == 'saint')):
            return True
        if ((self.model.file.fileContent.format == 'sca') or
            (self.model.file.fileContent.format == 'mtz')):
            # merged = ['unk','merged' ,'unmerged']
            if (self.model.file.fileContent.merged == 'unmerged'):
                return True
            return False
        # default no batch info
        return False
        #if (self.model.file.fileContent.batchs.isSet()
        #and (self.model.file.fileContent.batchs != '-1 - 1')):

    def toggleBatchView(self):
        #print "toggleBatchView", self.model.file.fileContent.batchs
        #print "toggleBatchView", self.model.file.fileContent.format
        #print "toggleBatchView", self.model.file.fileContent.merged
        if (self.hasBatchInfo()):
            self.batchFrame.show()
            self.batchFrame2.show()
        else:
            self.batchFrame.hide()
            self.batchFrame2.hide()

    def loadSameAsCombo(self):
        #print 'CImportUnmergedView.loadSameAsCombo',self.model
        if self.model is None: return
        indx = self.model.parent().index(self.model)
        self.sameAsCombo.blockSignals(True)
        self.sameAsCombo.clear()
        self.sameAsCombo.blockSignals(False)
        if indx<1: return
        self.sameAsCombo.addItem(' ')
        for item in self.model.parent()[0:len(self.model.parent())-1]:
            self.sameAsCombo.addItem(item.getTextItem())

    @QtCore.Slot(int)
    def handleSameAsCombo(self,indx):
        #print 'CImportUnmergedView.handleSameAsCombo',indx
        if indx==0: return
        try:
            self.model.dataset.set(self.model.parent()[indx-1].dataset.__str__())
            self.model.crystalName.set(self.model.parent()[indx-1].crystalName.__str__())
        except:
            pass

    def updateViewFromModel(self):
        #print 'CImportUnmergedView.updateViewFromModel'
        CCP4Widgets.CComplexLineWidget.updateViewFromModel(self)
        self.widgets['file'].updateViewFromModel()
        self.widgets['crystalName'].validate()
        self.widgets['dataset'].validate()
        self.widgets['wavelength'].validate()
        self.widgets['excludeSelection'].validate()
        if self.model is None:
            self.widgets['batchsInFile'].setValue('')
        else:
            self.widgets['batchsInFile'].setValue(self.model.file.fileContent.batchs)

    def getMenuDef(self):
        menu = self.widgets['file'].getMenuDef()
        #menu[menu.index('view')] = ['View','view_ViewHKL','view_text']
        return menu

    def openViewer(self,mode):
        if not self.model.file.exists(): return
        ext = self.model.file.getExt()
        if ext is None:
            return
        if ext == '.mtz' or mode == 'view_ViewHKL':
            CCP4Modules.LAUNCHER().openInViewer('viewHKL',fileName=self.model.file.__str__())
        else:
            CCP4Modules.WEBBROWSER().openFile(self.model.file.__str__(),format='text/plain',toFront=True)

    # Just drag the file name
    def dragData(self):
        if self.model.file.isSet():
            tree = self.model.file.getEtree()
            tree.tag = 'UnmergedDataFile'
            from lxml import etree
            text = etree.tostring(tree,pretty_print=False)
            return text
        else:
            return None

    def dragType(self):
        return 'UnmergedDataFile'

    def acceptDropData(self,textData):
        from lxml import etree
        #print 'CImportUnmergedView.acceptDropData',textData
        tree = etree.fromstring(textData)
        self.connectUpdateViewFromModel(False)
        if tree.tag == 'UnmergedDataFile':
            self.model.file.unSet()
            self.model.file.setEtree(tree)
        elif  tree.tag == 'ImportUnmerged':
            self.model.unSet()
            self.model.setEtree(tree)
        self.connectUpdateViewFromModel(True)
        self.updateViewFromModel()
        self.validate()


class CImportUnmergedListView(CCP4Widgets.CListView):

    MODEL_CLASS = CCP4XtalData.CImportUnmergedList

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {'mode' : 'table','tableItems' : ['file','crystalName','dataset','excludeSelection'],
                  'columnHeaders':['Filename','Crystal','Dataset','Exclude batches'],}
        qualis.update(qualifiers)
        CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)


class CXia2ImageSelectionView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CXia2ImageSelection

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = { 'vboxLayout':True, 'iconButton':False}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        self.widgets['imageFile'] = CCP4Widgets.CDataFileView(parent=self,qualifiers={})
        self.widgets['imageFile'].setToolTip('xia2 will automatically find '
          'other images in the dataset, optionally limited to the specified'
          ' range. For HDF5 datasets, select the master file.')
        self.widgets['imageStart'] = CCP4Widgets.CIntView(parent=self,qualifiers={})
        self.widgets['imageEnd'] = CCP4Widgets.CIntView(parent=self,qualifiers={})
        self.imageFileFrame = QtWidgets.QFrame(self)
        self.imageFileFrame.setLayout(QtWidgets.QHBoxLayout())
        self.imageFileFrame.layout().addWidget(QtWidgets.QLabel("Image file",self))
        self.imageFileFrame.layout().addWidget(self.widgets['imageFile'])
        self.imageFileFrame.layout().setContentsMargins(0,0,0,0)
        self.imageFileFrame.layout().setSpacing(0)
        self.layout().addWidget(self.imageFileFrame)
        self.imageRangeFrame = QtWidgets.QFrame(self)
        self.imageRangeFrame.setLayout(QtWidgets.QHBoxLayout())
        self.imageRangeFrame.layout().addWidget(QtWidgets.QLabel("Image range",self))
        self.imageRangeFrame.layout().addWidget(self.widgets['imageStart'])
        self.imageRangeFrame.layout().addWidget(QtWidgets.QLabel('to',self))
        self.imageRangeFrame.layout().addWidget(self.widgets['imageEnd'])
        self.imageRangeFrame.layout().setContentsMargins(0,0,0,0)
        self.imageRangeFrame.layout().setSpacing(0)
        self.layout().addWidget(self.imageRangeFrame)


class CXia2ImageSelectionListView(CCP4Widgets.CListView):

    MODEL_CLASS = CCP4XtalData.CXia2ImageSelectionList

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {'mode' : 'table', 'tableItems' : ['imageFile','imageStart','imageEnd'],
                  'columnHeaders':['Image filename','Start','End'],}
        qualis.update(qualifiers)
        CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)


class CMmcifReflDataFileView(CCP4Widgets.CDataFileView):

    MODEL_CLASS = CCP4XtalData.CMmcifReflDataFile

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CDataFileView.__init__(self,parent=parent,model=model,qualifiers=qualis)

    def getMenuDef(self):
        return ['clear','view','sep','copy','editLabel','export','help']


class CMtzDataFileView(CCP4Widgets.CDataFileView):

    MODEL_CLASS = CCP4XtalData.CMtzDataFile

    def setModel(self,model=None):
        CCP4Widgets.CDataFileView.setModel(self,model=model)
        if self.model is not None:
            otherMtz = self.model.getDataByKey('sameCrystalAs')
            #print 'CMtzDataFileView.setModel connecting otherMtz',model.objectName(),repr(otherMtz)
            if otherMtz is not None:
                otherMtz.dataChanged.connect(self.validate)

    def getMenuDef(self):
        #print 'CMtzDataFileView.getMenuDef'
        if isinstance(self.model,(CCP4XtalData.CObsDataFile,CCP4XtalData.CPhsDataFile,CCP4XtalData.CFreeRDataFile)):
            viewSubMenu = ['View','view_ViewHKL']
        else:
            viewSubMenu = ['View','view_ViewHKL','view_CCP4mg','view_Coot']
        if self.editable:
            menu = ['clear',viewSubMenu,'sep','copy','paste','help']
        else:
            if getattr(self,'role',None) is not None and self.role == 'output':
                menu = [viewSubMenu,'sep','copy','editLabel','export','help']
            else:
                menu = [viewSubMenu,'sep','copy','export','help']
        if getattr(self,'_stacked',False):
            menu.insert(0,'handleStack')
        return menu


class CMiniMtzDataFileView(CMtzDataFileView):

    MODEL_CLASS = CCP4XtalData.CMiniMtzDataFile

    def dropTypes(self):
        return ["MtzDataFile"]

    def jobNumber(self):
        return self.parentTaskWidget().jobNumber()

    def handleBrowserOpenFile(self,filename,downloadInfo={}):
        self.model.blockSignals(True)
        CMtzDataFileView.handleBrowserOpenFile(self,filename,downloadInfo=downloadInfo,autoInfoOnFileImport=False,validate=False)
        self.model.blockSignals(False)
        if self.model.getExt() in ['.cif','.ent']:
            err = self.model.importFromCif(jobId=self.parentTaskWidget().jobId())
            if err.maxSeverity()>SEVERITY_WARNING:
                err.warningMessage('Importing experimental data','Failed loading cif format file',parent=self)
                self.model.unSet()
                return
        errors = self.model.validColumns()
        #print 'CMiniMtzDataFileView.handleBrowserOpenFile',errors.report()
        #print 'CMiniMtzDataFileView.handleBrowserOpenFile sourceFileAnnotation',self.model.__dict__.get('sourceFileAnnotation','')
        if errors.maxSeverity()>SEVERITY_WARNING:
            #print errors.report()
            if errors.count(cls=self.model.__class__,code=206)>0:
                QtWidgets.QMessageBox.warning(self,'Error in selected MTZ file','This file contains unmerged data - use Data Reduction task to import it')
                self.model.unSet()
                return
            elif errors.count(cls=self.model.__class__,code=203)>0:
                QtWidgets.QMessageBox.warning(self,'Error in selected MTZ file','Selected MTZ file does not contain correct type of data')
                self.model.unSet()
                return
            elif errors.count(cls=self.model.__class__,code=204)>0 or errors.count(cls=self.model.__class__,code=205)>0:
                #applyNow = (errors.count(cls=self.model.__class__,code=205)>0)
                #print 'handleBrowserOpenFile applyNow',applyNow,filename
                try:
                    self.dialog=CSelectColumnsWidget(parent=self,model=self.model,applyNow=False,filename=filename)
                except CException:
                    mess = 'This file '+ filename + '\ndoes not contain the appropriate data: '
                    for rC in self.model.requiredContent():
                        mess = mess + self.model.CONTENT_ANNOTATION[rC-1]+','
                    QtWidgets.QMessageBox.warning(self.parentTaskWidget(),'Error in selected MTZ file',mess[0:-1])
                    self.model.unSet()
                    return
                else:
                    self.dialog.applySignal.connect(self.handleDialogApply)
                    self.dialog.cancelSignal.connect(self.handleDialogCancel)
        else:
            # Selected file has correct complement of columns but needs to be imported to same
            # target filename as used by self.model.splitMtz()
            self.model.importFile(jobId=self.parentTaskWidget().jobId(),jobNumber=self.parentTaskWidget().jobNumber())
            if CCP4Modules.PREFERENCES().AUTO_INFO_ON_FILE_IMPORT and not self.model.dbFileId.isSet():
                # If this is a downloaded file then it will have some provenance info in
                # self.model.sourceFileAnnotation put there by CDataFileView.handleBrowserOpenFile()
                self.openInfo(label=self.model.qualifiers('guiLabel').lower(),sourceFileAnnotation=self.model.__dict__.get('sourceFileAnnotation',''))

    @QtCore.Slot()
    def handleDialogApply(self):
        selectedColumns = self.dialog.getSelection()
        #print 'handleDialogApply',selectedColumns
        try:
            self.dialog.close()
        except:
            pass
        jobId = self.parentTaskWidget().jobId()
        projectId = self.parentTaskWidget().projectId()
        sourceFileName = self.model.__str__()
        error = self.model.splitMtz(jobId=jobId,projectId=projectId,contentFlag=selectedColumns[0],i2Labels=selectedColumns[1],columnLabels=selectedColumns[3])
        #print 'CMiniMtzDataFileView.handleDialogApply',error
        self.model.dataChanged.emit()
        if error.maxSeverity()==SEVERITY_WARNING and error[0]['code']==212:
            mess = QtWidgets.QMessageBox.warning(self,self.windowTitle(),'This data is already imported as\n'+error[0]['details'])
            self.loadJobCombo()
            self.updateJobCombo()
            self.validate()
        elif error.maxSeverity()>=SEVERITY_WARNING:
            if error[0]['code']==211:
                mess = QtWidgets.QMessageBox.warning(self,self.windowTitle(),'No column data selected')
            else:
                error.warningMessage(windowTitle='Splitting MTZ: '+sourceFileName,jobId=jobId,parent=self)
            self.model.unSet()
            self.updateJobCombo()
            self.validate()
        else:
            self.loadJobCombo()
            self.updateJobCombo()
            self.validate()
            if CCP4Modules.PREFERENCES().AUTO_INFO_ON_FILE_IMPORT:
                #print 'CMiniMtzDaaFile.handleDialogApply sourceFileReference',self.model.__dict__.get('sourceFileReference','')
                # If this is a downloaded file then it will have some provenance info in
                # self.model.sourceFileAnnotation put there by CDataFileView.handleBrowserOpenFile()
                self.openInfo(label=self.model.qualifiers('guiLabel').lower()+' '+self.model.__dict__.get('sourceFileReference',''),
                              sourceFileAnnotation=self.model.__dict__.get('sourceFileAnnotation',''))
        self.dialog.deleteLater()

    @QtCore.Slot()
    def handleDialogCancel(self):
        self.dialog.close()
        self.dialog.deleteLater()
        self.model.unSet()


class CMergeMiniMtzView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS =CCP4XtalData.CMergeMiniMtz

    def __init__(self,parent=None,model=None,qualifiers={},**kw):
        qualis = {'vboxLayout' : True, 'iconName' : 'MiniMtzDataFile',
                  'dragType' : ['MiniMtzDataFile','ObsDataFile','PhsDataFile','FreeRDataFile','MapCoeffsDataFile'],
                  'iconButton' : True }
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis,**kw)
        self.setFrameShape(QtWidgets.QFrame.Box)
        iconWidgetItem = self.layout().takeAt(0)
        qualis['iconButton'] =  False
        del qualis['vboxLayout']
        self.widgets['fileName'] = CMiniMtzDataFileView(parent=self,qualifiers=qualis)
        self.widgets['columnTag'] = CCP4Widgets.CStringView(parent=self,qualifiers=qualis)
        self.widgets['columnTag'].widget.setMaximumWidth(300)
        self.widgets['columnTag'].widget.setMinimumWidth(300)
        self.widgets['columnNames'] = CCP4Widgets.CStringView(parent=self,qualifiers=qualis)
        self.widgets['fileName'].setToolTip('Select an experimental data object')
        self.widgets['columnNames'].setToolTip('Comma-separated list of output column names - spaces will be ignored')
        line = QtWidgets.QHBoxLayout()
        line.addWidget(iconWidgetItem.widget())
        line.addWidget( self.widgets['fileName'] )
        self.layout().addLayout( line )
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Tag for column names'))
        line.addWidget(self.widgets['columnTag'] )
        line.addWidget(QtWidgets.QLabel('(maximum 20 characters)'))
        line.addStretch(2)
        self.layout().addLayout( line )
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Output column names'))
        line.addWidget(self.widgets['columnNames'] )
        self.layout().addLayout( line )
        self.setModel(model)

    def openViewer(self,mode):
        self.widgets['fileName'].openViewer(mode)

    def setModel(self,model):
        CCP4Widgets.CComplexLineWidget.setModel(self,model)
        if model is None:
            self.widgets['fileName'].setModel(None)
            self.widgets['columnTag'].setModel(None)
            self.widgets['columnNames'].setModel(None)
        else:
            self.widgets['fileName'].setModel(model.fileName)
            self.widgets['columnTag'].setModel(model.columnTag)
            self.widgets['columnNames'].setModel(model.columnNames)
        # Beware the dataChanged signal from CDataFile may come before dbFileId is set
        self.model.fileName.dbFileId.dataChanged.connect(functools.partial(self.updateColumnNames,True))
        self.widgets['columnTag'].widget.editingFinished.connect(functools.partial(self.updateColumnNames,False))
        #self.model.setColumnTag(overwrite=False)
        self.model.setColumnNames('applyTag',overwrite=False)

    def updateViewFromModel(self):
        if self.model is None: return
        for item in ['fileName','columnTag','columnNames']:
            self.widgets[item].updateViewFromModel()

    def updateModelFromView(self):
        if self.model is None: return
        for item in ['fileName','columnTag','columnNames']:
            self.widgets[item].updateModelFromView()

    @QtCore.Slot(bool)
    def updateColumnNames(self,setColumnTag=False):
        if setColumnTag: self.model.setColumnTag(overwrite=True)
        self.model.setColumnNames('applyTag',overwrite=True)

    def getMenuDef(self):
        if self.editable:
            menu = ['clear','view_ViewHKL','sep','copy','paste','help']
        else:
            menu = ['view_ViewHKL','sep','copy','help']
        if self._stacked: menu.insert(0,'handleStack')
        return menu

    @QtCore.Slot('QMimeData')
    def acceptDropData(self,textData):
        #print 'CMergeMiniMtzView.acceptDropData',textData
        from lxml import etree
        tree = etree.fromstring(textData)
        tree.tag = self.dragType()
        #print 'CMergeMiniMtzView.acceptDropData',textData,tree.tag
        self.connectUpdateViewFromModel(False)
        self.model.fileName.unSet()
        self.model.fileName.setEtree(tree)
        #print 'CMergeMiniMtzView.acceptDropData model',self.model.fileName
        self.connectUpdateViewFromModel(True)
        self.updateViewFromModel()
        #print 'CMergeMiniMtzView.acceptDropData validity',self.model.validity(self.model.get()).report()
        self.validate()



class CMergeMiniMtzListView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS =CCP4XtalData.CMergeMiniMtzList

    def __init__(self, parent=None, model=None, qualifiers={},**kw):
        '''
        qualis =  { 'mode' : 'table',
                    'tableItems' : ['fileName','columnNames'],
                    'columnHeaders':['Filename','Column names']
                   }

        qualis.update(qualifiers)
        CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis)
        self.listWidget.setColumnWidth(0,270)
        self.listWidget.setColumnWidth(1,270)
        '''
        qualis = {'vboxLayout' : True, 'iconButton' : False}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis,**kw)
        self.widgets = []
        if self.editable:
            but = QtWidgets.QPushButton('Add input experimental data object',self)
            but.clicked.connect(self.handleAddObject)
            self.layout().addWidget(but)
        self.setModel(model)

    def setModel(self, model=None):
        self.model = model
        if model is None:
            return
        for item in model.__dict__['_value']:
            #print 'CMergeMiniMtzListView.setModel make widget for item',item
            self.widgets.append(CMergeMiniMtzView(self,model=item,qualifiers={'editable':self.editable}))
            self.layout().insertWidget(self.layout().count()-1,self.widgets[-1])

    @QtCore.Slot()
    def handleAddObject(self):
        self.model.addItem()
        self.widgets.append(CMergeMiniMtzView(self,model=self.model[-1],qualifiers={'editable':self.editable}))
        self.layout().insertWidget(self.layout().count()-1,self.widgets[-1])

    def setValue(self,value=[]):
        #print 'CMergeMiniMtzListView.setValue',value
        if len(value)>len(self.widgets):
            for mod in value[len(self.widgets):]:
                self.widgets.append(CMergeMiniMtzView(self,model=mod,qualifiers={'editable':self.editable}))
                self.layout().insertWidget(self.layout().count()-1,self.widgets[-1])
                #print 'CMergeMiniMtzListView.setValue',len(self.widgets),mod
        indx = 0
        for item in value:
            self.widgets[indx].setValue(item)
            indx += 1
        self.validate()

    def getValue(self,index=-1):
        #print 'CMergeMiniMtzListView.getValue NOT IMPLEMENTED'
        value = []
        if index<0:
            first = 0
            last = self.body.layout().count()
        else:
            first = index
            last = index + 1
        return value

    def updateViewFromModel(self):
        for w in self.widgets: w.updateViewFromModel()
        self.validate()

    def updateModelFromView(self):
        for w in self.widgets: w.updateModelFromView()


class CRunBatchRangeView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS =CCP4XtalData.CRunBatchRange

    def __init__(self,parent=None,model=None,qualifiers={},**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(**kw)
        #print('CRunBatchRangeView',qualis)
        CCP4Widgets.CComplexLineWidget. __init__(self,parent=parent,qualifiers=qualis,**kw)
        # Remove the icon which is not doing much for us
        icon = self.layout().takeAt(0).widget()
        icon.deleteLater()
        # Phil Evans November 2022: comment out connects here as they do not work
        # and do not seem to be needed
        for item in ['runNumber','batchRange0','batchRange1','fileNumber']:
            self.widgets[item] = CCP4Widgets.CIntView(parent=self,qualifiers=qualis)
            # self.widgets[item].editingFinished.connect(self.updateModelFromView)
            # self.widgets[item].enterKeyPress.connect(self.updateModelFromView)
        item = 'resolution'
        self.widgets[item] = CCP4Widgets.CFloatView(parent=self,qualifiers=qualis)
        # self.widgets[item].editingFinished.connect(self.updateModelFromView)
        # self.widgets[item].enterKeyPress.connect(self.updateModelFromView)

        self.fileNumberLabel = QtWidgets.QLabel('file number',self)
        self.layout().addWidget(QtWidgets.QLabel('Run number',self))
        self.layout().addWidget( self.widgets['runNumber'] )
        self.layout().addWidget(self.fileNumberLabel)
        self.layout().addWidget( self.widgets['fileNumber'] )
        self.layout().addWidget(QtWidgets.QLabel('for batch range',self))
        self.layout().addWidget(self.widgets['batchRange0'])
        self.layout().addWidget(QtWidgets.QLabel('to',self))
        self.layout().addWidget(self.widgets['batchRange1'])
        self.layout().addWidget(QtWidgets.QLabel('resolution limit',self))
        self.layout().addWidget(self.widgets['resolution'])
        self.layout().addStretch(1)
        self.showFileNumber(qualis.get('showFileNumber',False))
        if model is not None:
            self.setModel(model)

    def showFileNumber(self, mode):
        self.fileNumberLabel.setVisible(mode)
        self.widgets['fileNumber'].setVisible(mode)

    def setModel(self, model):
        # Reimplement to ensure this is validiated when any data item changes
        # (why is this not done in the base class?)
        if self.model is not None:
            for item in ['runNumber','batchRange0','batchRange1','resolution']:
                self.model.__getattr__(item).dataChanged.disconnect(self.validate)
        CCP4Widgets.CComplexLineWidget.setModel(self,model)
        if self.model is not None:
            for item in ['runNumber','batchRange0','batchRange1','resolution']:
                self.model.__getattr__(item).dataChanged.connect(self.validate)


class CRunBatchRangeListView(CCP4Widgets.CListView):

    MODEL_CLASS =CCP4XtalData.CRunBatchRangeList

    def __init__(self, parent=None, model=None, qualifiers={},**kw):
        if qualifiers.get('showFileNumber',False):
            qualis = {'mode' : 'table', 'tableItems' : ['runNumber', 'fileNumber',
                                                        'batchRange0', 'batchRange1','resolution'],
                      'columnHeaders' : ['Run number', 'File number', 'Batch start', 'Batch end','Resolution']}
        else:
            qualis = {'mode' : 'table', 'tableItems' : ['runNumber',
                                                        'batchRange0', 'batchRange1','resolution'],
                      'columnHeaders' : ['Run number', 'Batch start', 'Batch end','Resolution']}
        qualis.update(qualifiers)
        qualis.update(**kw)
        #print('CRunBatchRangeListView',qualis)
        CCP4Widgets.CListView.__init__(self,parent,model=model,qualifiers=qualis,editorQualifiers={ 'showFileNumber' :qualifiers.get('showFileNumber',False) } )


class CReindexOperatorView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4XtalData.CReindexOperator

    def __init__(self,parent=None,model=None,qualifiers={},**kw):
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualifiers)
        #self.layout().takeAt(0).widget().deleteLater()
        self.layout().addWidget(QtWidgets.QLabel('Reindex operator',self))
        for key in ['h','k','l']:
            self.widgets[key] = CCP4Widgets.CLineEdit(self)
            self.widgets[key].setToolTip('Enter operator for '+key)
            self.widgets[key].editingFinished.connect(self.updateModelFromView)
            #self.widgets[key].enterKeyPress.connect(self.updateModelFromView)
            self.layout().addWidget(QtWidgets.QLabel(key+'=',self))
            self.layout().addWidget(self.widgets[key])
            self.layout().setStretchFactor(self.widgets[key],2)
            if model is not None:
                self.setModel(model)

    def help(self):
        from core import CCP4Utils
        page = os.path.join(CCP4Utils.getCCP4Dir(),'html','reindexing.html')
        if not os.path.exists(page):
            page = os.path.join(CCP4Utils.getCCP4Dir(),'docs','reindexing.html')
        #print 'CReindexOperatorView.help',page,os.path.exists(page)
        if os.path.exists(page):
            CCP4Modules.WEBBROWSER().loadWebPage(fileName=page)


class CColumnGroupListView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS =CCP4XtalData.CColumnGroupList

    def __init__(self,parent=None,model=None,qualifiers={},**kw):
        qualis = { 'vboxLayout' : True }
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(self.layout().takeAt(0).widget())
        line.addWidget(QtWidgets.QLabel('Select column groups to be imported as separate data objects',self))
        self.layout().addLayout(line)
        self.grid = QtWidgets.QGridLayout()
        self.layout().addLayout(self.grid)
        col=0
        for item in ['Dataset','Type','Column labels and types']:
            l = QtWidgets.QLabel(item,self)
            l.setObjectName('bold')
            col = col + 1
            self.grid.addWidget(l,0,col)
        line = QtWidgets.QHBoxLayout()
        #for item in ['Select all','Deselect all','Unique observed data']:
        for item in ['Select all importable','Deselect all']:
            button = QtWidgets.QPushButton(item,self)
            line.addWidget(button)
            button.clicked.connect(functools.partial(self.handleButton,item))
        line.addStretch(1)
        self.layout().addLayout(line)
        if model is not None:
            self.setModel(model)

    @QtCore.Slot(str)
    def handleButton(self,mode):
        #print 'CColumnGroupListView.handleButton',mode
        if mode == 'Unique observed data':
            pass
        else:
            flag = (mode == 'Select all importable')
            for row in range(1,self.grid.rowCount()):
                w = self.grid.itemAtPosition(row,0).widget()
                if not isinstance(w,CCP4Widgets.CCheckBoxUneditable):
                    w.setChecked(flag)

    def clear(self):
        for row in range(1,self.grid.rowCount()):
            for col in (0, 1, 2, 3):
                item = self.grid.itemAtPosition(row,col)
                if item is not None:
                    item.widget().hide()
                    item.widget().deleteLater()

    def updateViewFromModel(self):
        #print 'CColumnGroupListView.updateViewFromModel'
        self.clear()
        for row in range(1,len(self.model.__dict__['_value'])+1):
            value = self.model.__dict__['_value'][row-1]
            if self.editable and len(value.columnGroupType)>0:
                cb = QtWidgets.QCheckBox(self)
                cb.setChecked(bool(value.selected))
                cb.clicked.connect(functools.partial(self.handleClick,row))
            else:
                cb = CCP4Widgets.CCheckBoxUneditable(parent=self)
            self.grid.addWidget(cb,row,0)
            self.grid.addWidget(QtWidgets.QLabel(value.dataset.__str__(),self),row,1)
            self.grid.addWidget(QtWidgets.QLabel(value.columnGroupType.__str__(),self),row,2)
            self.grid.addWidget(QtWidgets.QLabel(value.columnListStr(),self),row,3)

    @QtCore.Slot(int,str)
    def handleClick(self,row,status):
        print('handleClick',row,status)
        self.model.blockSignals(True)
        self.model.__dict__['_value'][row-1].selected.set(status)
        self.model.blockSignals(False)

    @QtCore.Slot()
    def updateModelFromView(self):
        for row in range(1,self.grid.rowCount()):
            #print 'updateModelFromView', row,self.grid.itemAtPosition(row,0).widget()
            self.model.__dict__['_value'][row-1].selected.set(self.grid.itemAtPosition(row,0).widget().isChecked())

    def numberSelected(self):
        n = 0
        for row in range(1,len(self.model.__dict__['_value'])+1):
            if self.grid.itemAtPosition(row,0).widget().isChecked():
                n += 1
        return n


class CGenericReflDataFileView(CMtzDataFileView):

    MODEL_CLASS = CCP4XtalData.CGenericReflDataFile

    def getMenuDef(self):
        viewSubMenu = ['View','view_ViewHKL','view_text']
        if self.editable:
            menu = ['clear',viewSubMenu,'sep','copy','paste','help']
        else:
            if self.role is not None and self.role == 'output':
                menu = [viewSubMenu,'sep','copy','editLabel','export','help']
            else:
                menu = [viewSubMenu,'sep','copy','export','help']
        if self._stacked:
            menu.insert(0,'handleStack')
        return menu


class CSelectColumnsWidget(QtWidgets.QDialog):

    applySignal = QtCore.Signal()
    cancelSignal = QtCore.Signal()

    ERROR_CODES = {300 : {'description' : 'There is no data of required type in MTZ file'}}

    def __init__(self,parent=None,model=None,applyNow=False,filename=None):
        from qtgui import CCP4TaskWidget
        from core import CCP4DataManager
        QtWidgets.QDialog.__init__(self,parent)
        self.model = model
        baseName = str(self.model.baseName)
        if len(baseName)>0:
            baseName = os.path.splitext(baseName)[0]
        self.setModal(True)
        self.setWindowModality(QtCore.Qt.WindowModal)
        self.setWindowTitle('Select columns from MTZ file')
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(3,3,3,3)
        self.layout().setSpacing(3)
        self.setMaximumWidth(CCP4TaskWidget.WIDTH)
        self.layout().addWidget(QtWidgets.QLabel('Extracting data from the file: '+filename,self))
        self.layout().addWidget(QtWidgets.QLabel('Select one set of columns for '+self.model.qualifiers('toolTip'),self))
        self.buttonGroup = QtWidgets.QButtonGroup(self)
        self.buttonGroup.setExclusive(True)
        columnGroupModelList = self.model.columnGroup()
        columnGroupsInFile = self.model.fileContent.getColumnGroups()
        #for item in columnGroupModelList:
        #  print 'CSelectColumnsWidget.__init__ columnGroupModelList',item.get()
        colGpType = self.model.__class__.__name__[1:-8]
        orLabel = 'Select '
        requiredContent = self.model.requiredContent()
        #print 'CSelectColumnsWidget.__init__ requiredContent',requiredContent
        if requiredContent is not None and 0 in requiredContent: requiredContent = None
        contentFlag = 0
        checked = True
        idxcolumnGroup = 0
        self.columnGroupWidgets = []
        self.contentFlags = []  # for each widget

        # Phil Evans Feb 2022
        # In previous versions, the widgets and buttons were indexed by
        # columnGroupModel, ie group type, eg F,SIGF etc.
        # This did not allow for the possibily of multiple datasets with
        # same data types.
        # The widgets etc are now just indexed by serial number
        for columnGroupModel in columnGroupModelList:
            contentFlag += 1
            cglidx = []
            for i in range(len(columnGroupsInFile)):
                item = columnGroupsInFile[i]
                if item.columnGroupType == colGpType and item.contentFlag == contentFlag:
                    # list of indices into columnGroupsInFile
                    cglidx.append(i)

            for idx in cglidx:
                # Phil Evans Feb 2022: we need to loop here to allow for
                #  more than one columnGroup in the file belonging to this
                #  columnGroupModel type
                columnGroupList = []
                if (requiredContent is None or contentFlag in requiredContent) and len(cglidx)>0:
                    #print( 'CSelectColumnsWidget.__init__ columnGroupModel',columnGroupModel)
                    self.layout().addWidget(QtWidgets.QLabel( orLabel + columnGroupModel.qualifiers('guiLabel').__str__().lower(),self))
                    columnGroupWidget = CCP4DataManager.DATAMANAGER().getWidgetClass(model=columnGroupModel)(self,model=columnGroupModel)
                    columnGroupWidget.setObjectName('contentFlag_'+str(contentFlag))
                    radio = QtWidgets.QRadioButton(self)
                    self.buttonGroup.addButton(radio)
                    # Buttons are now indexed by serial number rather than contentFlag
                    self.buttonGroup.setId(radio,idxcolumnGroup)
                    radio.setChecked(checked)
                    columnGroupWidget.layout().addWidget(radio,0,0)
                    # columnGroupList here always contains a single item,
                    #  as the code to handle multiple items is broken
                    columnGroupWidget.drawColumns(False)
                    columnGroupList.append(columnGroupsInFile[idx])
                    columnGroupWidget.populateColumns(columnGroupList)
                    self.layout().addWidget(columnGroupWidget)
                    self.columnGroupWidgets.append(columnGroupWidget)
                    # ... but we do need to store the contentFlag
                    self.contentFlags.append(contentFlag)
                    orLabel = '..or '
                    checked= False
                    idxcolumnGroup += 1

        ##
        #print 'CMiniMtzDataFileView.handleBrowserOpenFile',type(self.model),self.model.CONTENTS
        if idxcolumnGroup == 0:
            #There is no content of required type
            raise CException(self.__class__,300)
        if len(self.model.CONTENTS['subType']['qualifiers']['enumerators'])>0:
            line = QtWidgets.QHBoxLayout()
            subTypeWidget = CCP4DataManager.DATAMANAGER().getWidgetClass(model=self.model.subType)(self,model=self.model.subType,qualifiers=self.model.CONTENTS['subType']['qualifiers'])
            subTypeWidget.updateViewFromModel()
            line.addWidget(QtWidgets.QLabel('This is',self))
            line.addWidget(subTypeWidget)
            line.addStretch(1)
            self.layout().addLayout(line)
        line = QtWidgets.QHBoxLayout()
        self.annotationWidget = CCP4DataManager.DATAMANAGER().getWidgetClass(model=self.model.annotation)(self,model=self.model.annotation)
        self.annotationWidget.updateViewFromModel()
        line.addWidget(QtWidgets.QLabel('File label',self))
        line.addWidget(self.annotationWidget)
        self.layout().addLayout(line)
        self.setAnnotation()
        self.buttonGroup.buttonClicked[int].connect(self.setAnnotation)
        for columnGroupWidget in self.findChildren(CProgramColumnGroupView):
            columnGroupWidget.model.dataChanged.connect(self.setAnnotation)
        buttonLine = QtWidgets.QHBoxLayout()
        buttonLine.addStretch(.5)
        buttonBox = QtWidgets.QDialogButtonBox(self)
        button = buttonBox.addButton(QtWidgets.QDialogButtonBox.Apply)
        button.clicked.connect(self.applySignal.emit)
        button = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
        button.clicked.connect(self.cancelSignal.emit)
        button = buttonBox.addButton(QtWidgets.QDialogButtonBox.Help)
        button.clicked.connect(self.handleHelp)
        buttonLine.addWidget(buttonBox)
        buttonLine.addStretch(.5)
        self.layout().addLayout(buttonLine)
        if not applyNow:
            self.show()
        else:
            self.applySignal.emit()
            # Correct column types - wrong column labels

    @QtCore.Slot()
    def setAnnotation(self):
        baseName = str(self.model.baseName)
        if len(baseName)>0: baseName = os.path.splitext(baseName)[0]
        rv = self.getSelection()
        if rv is None: return
        contentFlag, i2Names, dataset, colLabels = self.getSelection()
        label = ''
        for item in colLabels: label = label + item + ','
        #label =self.model.qualifiers('guiLabel')
        jN = self.parent().jobNumber()
        if jN is not None:
            jN = ' imported by job '+jN
        else:
            jN = ''
        #print('setAnnotation',baseName,'*',dataset,'*',label,'*',jN)
        #self.model.annotation.set(baseName + ' ' + dataset + '/'+label[0:-1] + jN)
        self.model.annotation.set(baseName + ': ' + dataset + jN)
        self.annotationWidget.updateViewFromModel()

    def getSelection(self):
        buttonId =  self.buttonGroup.checkedId()  # index into lists
        columnGroupWidget = self.columnGroupWidgets[buttonId]
        contentFlag = self.contentFlags[buttonId]
        #print("getSelection_1", buttonId, contentFlag, columnGroupWidget)
        if columnGroupWidget is None: return None
        colLabels = []
        for w in columnGroupWidget.widgets:
            #print 'CSelectColumnsWidget.getSelection currentText',w.currentIndex(),w.currentText()
            if isinstance(w,CCP4Widgets.CLabel):
                txt = w.text().__str__().split('/')
            else:
                txt = w.currentText().__str__().split('/')
            colLabels.append(txt[-1])
        if len(txt)>1:
            dataset = txt[0]
        else:
            dataset = ''
        #print 'getSelection',columnGroupWidget.model.columnGroupNames(),dataset,colLabels
        return contentFlag, columnGroupWidget.model. columnGroupNames(), dataset, colLabels

    @QtCore.Slot()
    def handleHelp(self):
        CCP4Modules.WEBBROWSER().loadWebPage(helpFileName='data_files')

