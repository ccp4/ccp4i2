"""
    nautilus_build_refine_gui.py: CCP4 GUI Project
    
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

from PySide2 import QtCore
from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class nautilus_build_refine_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'nautilus_build_refine'
    TASKVERSION = 0.1
    TASKMODULE = [ 'model_building' ]
    SHORTTASKTITLE='Autobuild RNA - NAUTILUS'
    TASKTITLE='Autobuild RNA/DNA (Nautilus pipeline)'
    DESCRIPTION = 'Iterations of model building (Nautilus) and refinement (Refmac5)'
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile('nautilus')

        self.openFolder(folderFunction='inputData',followFrom=False)
        
        self.createLine( [ 'subtitle', 'Select experimental data', 'Enter your observations here.' ] )
        self.openSubFrame( frame=[True] )
        self.createLine( [ 'widget', 'F_SIGF' ] )
        self.createLine( [ 'widget', 'ABCD' ] )
        self.createLine( [ 'widget', 'FREERFLAG' ] )
        self.closeSubFrame()

        self.createLine( [ 'subtitle', 'Enter the AU content containing the structure sequence(s)' ] )
        self.openSubFrame( frame=[True] )
        self.createLine( [ 'widget','ASUIN' ] )
        self.closeSubFrame()

        self.createLine( [ 'tip', 'Your model should contain confidently-built metals, ligands and/or protein parts', 'widget', 'XYZIN_MODE', 'label', 'Start from a partially built model' ] )
        self.container.controlParameters.XYZIN_MODE.dataChanged.connect(self.updateRequirements )        
        self.openSubFrame( frame=[True], toggle=[ 'XYZIN_MODE', 'open', [True] ] )
        self.createLine( [ 'widget','XYZIN' ] )
        self.closeSubFrame()

        self.closeFolder()

# a control parameter tab with basic options starts here

        folder = self.openFolder( folderFunction='controlParameters', title='Options' )

        self.createLine ( [ 'subtitle', 'Pipeline control', '<b>Warning</b>: default values work best for most situations.' ] )
        self.openSubFrame ( frame = [True] )
        self.createLine( [ 'label', 'Iterate', 'widget', 'ITERATIONS', 'label', 'times over', 'widget', 'NAUTILUS_CYCLES', 'label', 'cycles of build and', 'widget', 'REFMAC_CYCLES', 'label', 'of refinement' ] )
        self.closeSubFrame ( )

        self.createLine( [ 'widget', 'NAUTILUS_ANISOTROPY_CORRECTION','label','Apply anisotropy correction' ] )

        self.closeFolder()

# a control parameter tab with more advanced options starts here

        folder = self.openFolder( folderFunction='controlParameters', title='Advanced Nautilus Options' )

        self.openSubFrame ( frame=[True])
        self.createLine( [ 'tooltip', '(Using data beyond 2A is slow and seldom helps)', 'label','Leave out reflection data beyond a','widget','NAUTILUS_RESOLUTION', 'label', '&#8491; resolution limit'] )
        self.closeSubFrame ( )

        self.closeFolder()



    @QtCore.Slot()
    def updateRequirements( self ):
        if self.container.controlParameters.XYZIN_MODE :
          self.container.inputData.XYZIN.setQualifier ( 'allowUndefined', False )
          self.container.inputData.XYZIN.setQualifier ( 'toolTip', 'An initial atomic model is required' )
        else :
          self.container.inputData.XYZIN.setQualifier ( 'allowUndefined', True )
        self.getWidget ( 'XYZIN' ).validate()
        return

