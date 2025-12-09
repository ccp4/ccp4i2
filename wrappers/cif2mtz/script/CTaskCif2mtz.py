"""
     tasks/cif2mtz/CTaskCif2mtz.py
     Copyright (C) 2011 STFC
     Author: Martyn Winn

"""

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.baselayer import QtCore

#-------------------------------------------------------------------
class CTaskCif2mtz(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'cif2mtz'
  TASKVERSION = 0.1
  TASKMODULE='test'
  TASKTITLE='Import merged mmCIF X-ray data'
  DESCRIPTION = 'Import merged mmCIF X-ray data e.g. from the PDB database'

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)


  @QtCore.Slot()
  def updateFromCifFile(self):
    #print 'CTaskCif2mtz.updateFromCifFile',self.container.inputData.HKLIN.fileContent
    self.container.inputData.SPACEGROUPCELL.cell.set(self.container.inputData.HKLIN.fileContent.cell)
    self.container.inputData.SPACEGROUPCELL.spaceGroup.set(self.container.inputData.HKLIN.fileContent.spaceGroup)

  def drawContents(self):

    self.setProgramHelpFile('cif2mtz')

    self.openFolder(folderFunction='inputData',followFrom=False)

    self.createLine ( [ 'widget','HKLIN' ] )
    self.createLine ( [ 'widget', 'SPACEGROUPCELL' ] )

    self.updateFromCifFile()
    self.container.inputData.HKLIN.dataChanged.connect(self.updateFromCifFile)
