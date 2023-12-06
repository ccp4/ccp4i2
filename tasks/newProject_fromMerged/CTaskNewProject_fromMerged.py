from __future__ import print_function

"""
     tasks/demo.py: CCP4 GUI Project
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

"""
     Liz Potterton June 2011 - task to import data files
"""

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from core import CCP4XtalData
from core.CCP4Modules import WEBBROWSER,PROJECTSMANAGER,MIMETYPESHANDLER

def whatNext(jobId):
  container = PROJECTSMANAGER().getJobParams(jobId)
  if container is None: return []
  mode = container.find('MODE').get()
  # allowed modes..unknown,molrep,phasing,solved
  if mode == 'molrep':
    return ['molrep_mr','phaser_pipeline']
  elif mode == 'phasing':
    return []
  elif mode == 'solved':
    return ['refmac_basic']
  else:
    return []
  

#-------------------------------------------------------------------
class CTaskNewProject_fromMerged(CCP4TaskWidget.CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'newProject_fromMerged'
  TASKVERSION = 0.0
  TASKMODULE='test'
  TASKTITLE='Import project starting data from merged data'
  EDITOR = True
  AUTOPOPULATEINPUT = False


  def drawContents(self):

    self.openFolder(title='Experimental data')

    self.createLine( ['advice','Enter the file containing the data (supported formats: MTZ, mmCIF)'] )

    self.createLine( ['widget','-dropTypes',
                     ['MtzDataFile','MmcifReflDataFile'],'HKLIN'])
    self.container.inputData.HKLIN.dataChanged.connect(self.updateFromInputFile)

    # if HKLIN doesn't have both spacegroup and symmetry information, then ask for it
    self.createLine ( [ 'widget', 'SPACEGROUPCELL'],
                      toggle=['NEED_USER_CELLSYMM','open',[True]])

    # if HKLIN doesn't have a freeR column, then advise that we will generate one
    self.createLine( ['advice','Your data file has no freeR flags, so we will generate them.'],
                     toggle=['FREER_IN_HKLIN','open',[False]] )

    self.createLine( ['label','How do you expect to solve the structure?','widget','MODE','stretch'])

    self.openFolder(title='Composition for MR',toggle=['MODE','open',['molrep']])
    self.createLine( ['advice','Enter sequence of structure to be solved'])
    #self.createLine( ['widget','SEQUENCE'])
    self.createLine( ['widget','SEQFILEIN'])
    self.createLine( ['widget','CONTAINS_SEMET',
                      'label','Structure contains selenomethionine'])

    self.createLine( ['advice','Enter suggested homologous structures'])
    '''
    group = self.createRadioGroup(itemList=['No homologs','Use structures in directory','Enter structure files'])
    group.button(0).setChecked(True)
    '''
    self.createLine( ['widget','-guiMode','radio','MR_MODE'])
    stack = self.openStack(controlVar='MR_MODE')
    self.createLine(['label','No homologous structures'])
    self.createLine(['widget','MR_MODEL_DIR'])
    self.createLine(['widget','MR_MODEL_LIST'])
    self.closeStack()
    mr_model_list = self.getWidget('MR_MODEL_LIST')
    if mr_model_list.getNofLines()<1: mr_model_list.setNofLines(1)
    mr_model_list.setMinimumWidth(700)

    print('CTaskNewProject_fromMerged find MR_MODEL_DIR',self.findChildren(QtWidgets.QWidget,'MR_MODEL_DIR'))
    self.openFolder(title='Composition for ab initio',toggle=['MODE','open',['phasing']])

    
    #self.createLine( ['widget','COMPOSITION'])


  @QtCore.Slot()
  def updateFromInputFile(self):
    """
      Update knowledge of HKLIN contents whenever file is selected/changed.
    """

    #print "updateFromInputFile() called"
    importMimeType = MIMETYPESHANDLER().formatFromFileExt(self.container.inputData.HKLIN.fullPath)

    if importMimeType == 'chemical/x-cif':

      cifdata = CCP4XtalData.CMmcifReflData()
      cifdata.loadFile(self.container.inputData.HKLIN.fullPath)

      if cifdata.__dict__['_value']['cell'].get('a') == None or cifdata.__dict__['_value']['spaceGroup'] == '':
         self.container.inputData.NEED_USER_CELLSYMM = True
      else:
         self.container.inputData.NEED_USER_CELLSYMM = False

      if cifdata.__dict__['_value']['haveFreeRColumn']:
         self.container.inputData.FREER_IN_HKLIN = True
      else:
         self.container.inputData.FREER_IN_HKLIN = False

    else:
      self.container.inputData.NEED_USER_CELLSYMM = False
      self.container.inputData.FREER_IN_HKLIN = True

