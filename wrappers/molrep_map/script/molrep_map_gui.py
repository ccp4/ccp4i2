"""
tasks/molrep_map
"""

from PySide2 import QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class Cmolrep_map(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'molrep_map'
    TASKVERSION = 0.1
    TASKMODULE='molecular_replacement'
    TASKTITLE  = 'Molrep into cryo EM map'
    SHORTTASKTITLE  = 'Molrep into cryo EM map'
    DESCRIPTION='''Molrep into cryo EM map'''
      
    def drawContents(self):
        
        self.setProgramHelpFile('molrep_map')
        folder = self.openFolder(folderFunction='inputData',title='Basic data')
        
        self.createLine( [ 'subtitle', 'Model:' ] )
        self.createLine( [ 'widget', 'XYZIN' ] )
        self.createLine( [ 'subtitle', 'Density maps' ] )
        self.createLine( [ 'widget', 'MAPIN' ] )
        self.closeSubFrame()
        
        for scriptName in self.container.controlParameters.contents():
            self.autoGenerate(container=getattr(self.container.controlParameters, scriptName), subFrame=False)

