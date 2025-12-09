"""
     tasks/scaleit/CTaskscaleit.py
     Copyright (C) 2011 STFC
     Author: Phil Evans

"""
from ccp4i2.baselayer import QtCore
from qtgui.CCP4TaskWidget import CTaskWidget

from pipelines.aimless_pipe.script.aimless_pipe_utils import CellCheck, CellFormat, colourText
from ccp4i2.wrappers.scaleit.script.scaleit_utils import DatalistCheck
from pipelines.import_merged.script.dybuttons import MyMessageBox


class CTaskscaleit(CTaskWidget):

  TASKNAME = 'scaleit'
  TASKVERSION = 0.0
  TASKMODULE= ['expt_data_utility' ]
  TASKTITLE='Compare two or more datasets'
  DESCRIPTION="Compare two or more datasets by running SCALEIT"
  
  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.openFolder(folderFunction='inputData')
    self.createLine (['label',
                      'Assign two or more datasets in the list: the first one will be treated as "native"'])

    self.createLine(['widget', '-listVisible', True, 'MERGEDFILES'])
    self.container.inputData.MERGEDFILES.dataChanged.connect\
                           (self.handleSelectFile)

    self.createLine( ['label','Maximum resolution used (\xc5)',
                      'widget', 'RESOLUTION_MAX'])             

    self.createLine (['label',
                      'The program SCALEIT will be run with default options to compare datasets to the first one (treated as "native")'])
    self.createLine (['label',
                      'The other "derivative" datasets will be scaled to the first one with a scale and anistropic B-factor'])
    self.createLine (['label',
                      'Statistics are analysed by resolution, and by comparison of amplitudes F, not intensities'])
    self.createLine (['label',
                      'If there is just one "derivative" (ie two datasets), the differences are also analysed by a Normal Probability plot'])
    self.createLine (['label',
                      'Data are analysed as amplitudes F even if the input contains intensities'])
    self.createLine (['label',
                      'Datasets from outside ccp4i2 should first be imported using the "Import merged reflection data" task'])


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  @QtCore.Slot()
  def handleSelectFile(self):
    # nothing to do at present
    pass

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def isValid(self):
    invalidElements = super().isValid()
    # Check whether invocation is from runTask  (cf prosmart_refmac_gui)
    import traceback
    functionNames = [a[2] for a in traceback.extract_stack()]
    # print("\n*** isValid call ", functionNames)
    if functionNames[-2] != 'runTask':
      #  ignore other entry points eg createTaskWidget
      return invalidElements
      
    datalistcheck = DatalistCheck(self.container.inputData.MERGEDFILES)
    OK = datalistcheck.popup()
    if OK == 'Onefile' or OK == 'No':
      invalidElements.append(self.container.inputData.MERGEDFILES)
      
    #  OK 
    return invalidElements

