"""
    mergeMtz_gui.py
    Copyright (C) 2015 SFTC    
    """

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class mergeMtz_gui(CTaskWidget):
#-------------------------------------------------------------------
    
    TASKMODULE = 'export'
    TASKTITLE = 'Merge experimental data objects to MTZ'
    SHORTTASKTITLE = 'Merge to MTZ'
    TASKVERSION = 0.1
    DESCRIPTION = "Export 'old' style MTZ file"
    TASKNAME = 'mergeMtz'
    

    def drawContents(self):        
      self.openFolder(folderFunction='inputData',followFrom=False)
      self.createLine(['subtitle','Select data objects'])
      self.createLine(['widget','MINIMTZINLIST'])
       
      self.createLine(['subtitle','Enter the export file name'])
      self.createLine(['widget','HKLOUT'])

