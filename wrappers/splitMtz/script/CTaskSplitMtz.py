from __future__ import print_function

"""
     tasks/splitMtz/CTaskSplitMtz.py: CCP4 GUI Project
     Copyright (C) 2012 STFC

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
     Liz Potterton August 2012 - gui for mtz split
"""

import functools
from PySide6 import QtGui, QtWidgets,QtCore
from core.CCP4ErrorHandling import *
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskSplitMtz(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'splitMtz'
  TASKVERSION = 0.0
  TASKMODULE='data_entry'
  TASKTITLE='Import and Split MTZ into experimental data objects'
  SHORTTASKTITLE='Import and Split MTZ'
  DESCRIPTION = 'Select groups of columns from the MTZ file (csplitmtz)'
  WHATNEXT = []

  ERROR_CODES = { 200 : { 'description' : 'There are no selected column groups' } }
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    self.dataTypeList = \
      [ [ 'Anomalous intensities' , 'Obs' , 1 ],
                         [ 'Anomalous SFs' , 'Obs' , 2 ],
                         [ 'Mean intensities' , 'Obs' , 3 ],
                         [ 'Mean SFs' , 'Obs' , 4 ],
                         ['HL Phases' ,'Phs', 1 ],
                         ['Phi-FOM Phases' ,'Phs', 2 ],
                         ['Map coefficients', 'MapCoeffs' , 1]
                         ]


  def drawContents(self):                    
    folder = self.openFolder(folderFunction='inputData',title='Input Data',followFrom=False)
    
    self.createLine(  [ 'subtitle', 'Select experimental data file'] )
    self.createLine ( [ 'widget','HKLIN' ] )
    self._container.inputData.HKLIN.dataChanged.connect(self.handleSelectHklin)
    self.createLine(  [ 'subtitle', 'Select column groups'] )
    self.createLine([ 'widget', 'COLUMNGROUPLIST'])
    
    # Could put this in a  different tab
    #folder = self.openFolder(folderFunction='controlParameters',title='Advanced')

    line = self.createLine(['subtitle','Select individual column data for'])
    self.butGroup = QtWidgets.QButtonGroup(self)
    idx = 0
    for label,dType,subtype in self.dataTypeList:
      if idx%4==0 : line = self.createLine([])
      but = QtWidgets.QRadioButton(label,self)
      self.butGroup.addButton(but,idx)
      line.addWidget(but)
      but.clicked.connect(functools.partial(self.drawSelectColumns,dType,subtype))
      idx+=1
      
    line = self.createLine([])
    self.selectFrame = QtWidgets.QFrame(self)
    self.selectFrame.setLayout(QtWidgets.QGridLayout())
    line.addWidget(self.selectFrame)


  @QtCore.Slot(str,str)
  def drawSelectColumns(self,dataType,contentFlag):
    from core import CCP4XtalData
    for iR in (0,1):
      for iC in (0,1):
        layoutItem = self.selectFrame.layout().itemAtPosition(iR,iC)
        if not layoutItem is None:
          layoutItem.widget().hide()
          layoutItem.widget().deleteLater()

    # Beware for each of these data types there is a list of possible 'correctColumns' patterns
    if dataType == 'Obs':
     self.currentColumnTypes = CCP4XtalData.CObsDataFile.QUALIFIERS['correctColumns'][contentFlag-1]
    elif dataType == 'Phs':
     self.currentColumnTypes = CCP4XtalData.CPhsDataFile.QUALIFIERS['correctColumns'][contentFlag-1]
    elif dataType == 'MapCoeffs':
     self.currentColumnTypes = CCP4XtalData.CMapCoeffsDataFile.QUALIFIERS['correctColumns'][contentFlag-1]

    print('selectColumns columnTypes', self.currentColumnTypes)
    nframes = 0
    self.comboList = []
    for cType in self.currentColumnTypes:
      frame = QtWidgets.QFrame(self)
      frame.setLayout(QtWidgets.QHBoxLayout())
      frame.layout().setContentsMargins(0,0,0,0)
      frame.layout().setSpacing(0)
      label = QtWidgets.QLabel(cType,frame)
      frame.layout().addWidget(label)
      self.comboList.append( CCP4Widgets.CComboBox(frame) )
      self.comboList[-1].setMinimumContentsLength(25)
      self.comboList[-1].setEditable(False)
      self.comboList[-1].addItem('Select column..')
      for col in self.container.inputData.HKLIN.fileContent.getListOfColumns([cType]):
        print('cType',cType,col)
        self.comboList[-1].addItem(str(col.dataset)+'/'+str(col.columnLabel))
      frame.layout().addWidget(self.comboList[-1])
      self.selectFrame.layout().addWidget(frame,nframes/2,nframes%2)
      nframes+=1
      
  @QtCore.Slot()
  def handleSelectHklin(self):
    #print 'CTaskSplitMtz.handleSelectHklin',self._container.inputData.HKLIN
    self._container.inputData.HKLIN.loadFile()
    self._container.inputData.COLUMNGROUPLIST.set(self._container.inputData.HKLIN.fileContent.getColumnGroups())

    # Redo the individual column input
    butId = self.butGroup.checkedId()
    #print 'CTaskSplitMtz.handleSelectHklin butId',butId
    if butId>=0: self.drawSelectColumns(self.dataTypeList[butId][1],self.dataTypeList[butId][2])
    
  def fix(self):
     self.getWidget('COLUMNGROUPLIST').updateModelFromView()
     nGroupsSelected = self.getWidget('COLUMNGROUPLIST').numberSelected()

     userSelected = False
     self.container.inputData.USERCOLUMNGROUP.unSet()
     userGroup = self.container.inputData.USERCOLUMNGROUP
     if self.butGroup.checkedId()>=0:
       selectedCols = []
       for c in self.comboList:
         if c.currentIndex()!=0 and selectedCols.count(str(c.currentText()))==0:
           selectedCols.append(str(c.currentText()))
       if len(selectedCols)==len(self.comboList):
         userGroup.columnGroupType = self.dataTypeList[self.butGroup.checkedId()][1]
         userGroup.contentFlag = self.dataTypeList[self.butGroup.checkedId()][2]
         idx = 0
         for col in selectedCols:
           userGroup.columnList.addItem()
           userGroup.columnList[-1].dataset = col.split('/')[0]
           userGroup.columnList[-1].columnLabel = col.split('/')[1]
           userGroup.columnList[-1].columnType = self.currentColumnTypes[idx]
           idx += 1
         userGroup.dataset.set(userGroup.columnList[0].dataset)
         userSelected = True
           
     if nGroupsSelected == 0 and not userSelected:
       return CErrorReport(self.__class__,200,stack=False)
     else:
       return  CErrorReport()
