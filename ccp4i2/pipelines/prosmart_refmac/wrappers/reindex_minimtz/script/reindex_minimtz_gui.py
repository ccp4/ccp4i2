"""
     tasks/reindex_minimtz/Creindex_minimtz.py: CCP4 GUI Project




Maritn Noble messed around with this
"""

from ......qtgui import CCP4TaskWidget


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
