"""
     dr_mr_modelbuild_gui.py: CCP4 GUI Project
     Copyright (C) 2021 University of York

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

from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

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
