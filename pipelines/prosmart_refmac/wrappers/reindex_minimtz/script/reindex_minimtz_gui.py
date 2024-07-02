
"""
     tasks/reindex_minimtz/Creindex_minimtz.py: CCP4 GUI Project
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
Maritn Noble messed around with this
"""

from PySide6 import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets


class Creindex_minimtz(CCP4TaskWidget.CTaskWidget):

  TASKTITLE='Reindex any miniMTZ'
  TASKNAME = 'reindex_minimtz'
  TASKVERSION = 0.0
  TASKMODULE = 'test'
  MGDISPLAYFILES = []
  AUTOPOPULATEINPUT = False

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('reindex_minimtz')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data')
    self.createLine( [ 'advice',' '] )
    self.createLine( [ 'advice',' '] )
    self.createLine( [ 'advice', 'Minimal inputs :' ])
    self.createLine( [ 'widget', 'HKLIN' ] )

#-   --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='controlParameters',title='Options')
    self.createLine( [ 'label', 'Reindexing operator', 'stretch', 'widget', 'OPERATION' ] )
    self.createLine( [ 'label', 'New space group', 'stretch', 'widget', 'NEWSPACEGROUP' ] )
