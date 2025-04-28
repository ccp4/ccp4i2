"""
Copyright (C) 2021 University of York
"""

from ......qtgui import CCP4TaskWidget


class Cmrbump_model_prep(CCP4TaskWidget.CTaskWidget):
  TASKNAME = 'mrbump_model_prep'
  TASKVERSION = 0.0
  TASKTITLE='MrBUMP model preparation'
  SHORTTASKTITLE = 'MrBUMP model prep'
  DESCRIPTION = 'Model preparation with MrBUMP for data reduction/MR/model build pipeline'
  RANK = 1

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):
      folder = self.openFolder(folderFunction='inputData',title='Input Data')
      self.createLine( [ 'widget','ASUIN' ] )
      self.createLine( [ 'widget','F_SIGF' ] )
      self.createLine( [ 'widget','FREERFLAG' ] )
      self.createLine( [ 'advice', 'Non-redundancy level for homologue search:'] )
      self.createLine( [ 'widget', 'REDUNDANCYLEVEL' ] )
      self.createLine( [ 'advice', 'Maximum no. of search models to create:'] )
      self.createLine( [ 'widget', 'MRMAX' ] )
