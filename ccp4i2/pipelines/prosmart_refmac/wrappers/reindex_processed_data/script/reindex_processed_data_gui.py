
"""
     tasks/reindex_processed/Creindex_processed_data.py: CCP4 GUI Project



"""

"""
Martin Noble messed around with this
"""

from ......qtgui import CCP4TaskWidget


class Creindex_processed_data(CCP4TaskWidget.CTaskWidget):

  TASKTITLE='Reindex observed and FreeR data'
  DESCRIPTION= 'Reindex processed data to correspond to a master dataset (Reindex)'
  TASKNAME = 'reindex_processed_data'
  TASKVERSION = 0.0
  TASKMODULE = 'test'
  AUTOPOPULATEINPUT = True

  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('reindex')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data')
    self.createLine( [ 'widget', 'OBSIN' ] )
    self.createLine( [ 'widget', 'FREERIN' ] )

    self.createLine( [ 'widget', 'OPERATION' ] )
    self.createLine( [ 'widget', '-label', 'New space group' ,'NEWSPACEGROUP' ] )
