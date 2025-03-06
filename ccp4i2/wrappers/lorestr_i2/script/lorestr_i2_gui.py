"""
    lorestr_i2_gui.py: CCP4 GUI Project
    
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

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class lorestr_i2_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'lorestr_i2' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'refinement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='LORESTR'
    TASKTITLE='Low Resolution Refinement Pipeline (LORESTR)'
    DESCRIPTION = '''Automated Low Resolution Structure Refinement Pipeline (LORESTR)'''
    MGDISPLAYFILES = ['XYZIN', 'XYZOUT','FPHIOUT','DIFFPHIOUT']
    WHATNEXT = ['coot_rebuild', 'prosmart_refmac', 'edstats']
    AUTOPOPULATEINPUT = True



    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
#        self.openFolder(folderFunction='inputData',followFrom=False)
        
#        self.createLine(['subtitle','Input coordinates'])
        self.setProgramHelpFile('refmac')

        indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

        #-  --------------------          --------------------          --------------------

        folder = self.openFolder(folderFunction='inputData',title='Input Data')
        #self.createLine( [ 'advice',' '] )
        self.createLine( [ 'subtitle', 'Main inputs' ])
        self.openSubFrame(frame=[True])
        self.createLine( [ 'widget', '-browseDb', True, 'XYZIN' ] )
        self.createLine( [ 'widget', '-browseDb', True, 'F_SIGF' ])
#        self.createLine( [ 'label','Use anomalous data for ', 'widget', 'USEANOMALOUSFOR', 'stretch', 'label', 'Wavelength', 'widget', 'WAVELENGTH'],toggleFunction=[self.anomalousAvailable,['F_SIGF']])
    #    if self.isEditable():
    #        if not self.container.controlParameters.WAVELENGTH.isSet(): self.getWavelength()
        self.createLine( [ 'widget', '-browseDb', True, 'FREERFLAG' ] )
        self.closeSubFrame()
        #self.createLine( [ 'advice',' '] )
        self.createLine( [ 'subtitle', 'Optional additional inputs' ])
        self.openSubFrame(frame=[True])
        self.createLine( [ 'widget', '-browseDb', True, 'TLSIN'] )
        self.createLine( [ 'widget', '-browseDb', True, 'DICT'] )
#        self.createLine( [ 'widget', '-browseDb', True, 'REFERENCE_MODEL' ] )
        self.createLine ( [ 'widget', 'REFERENCE_LIST' ] )
        self.closeSubFrame()

        #-   --------------------          --------------------          --------------------

        folder = self.openFolder(folderFunction='controlParameters',title='Options', drawFolder=self.drawControlParameters)

        #-   --------------------          --------------------          --------------------

        folder = self.openFolder(folderFunction='controlParameters',title='Advanced Options', drawFolder=self.drawAdvancedOptions)

#-  --------------------          --------------------          --------------------


    def drawControlParameters( self ):

        indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

        self.createLine( [ 'subtitle', 'Main options'] )
        self.openSubFrame(frame=[True])

        self.createLine( ['label', 'Automatically fetch homologues from following databases: ', 'widget', 'AUTO'])
        self.createLine( ['widget', 'OVB', 'label', 'Refine overall B-factors (for very low resolution)'])
        self.createLine( [ 'widget', 'DNA', 'label', 'Generate external restraints for DNA/RNA chains' ] )
        self.createLine( [ 'label', 'Number of CPUs to use:', 'widget', 'CPU' ] ) # 'stretch'

        self.closeSubFrame()

        self.createLine( [ 'subtitle', 'After Molecular Replacement'] )
        self.openSubFrame(frame=[True])

        self.createLine( [ 'widget', 'MR', 'label', 'Run 100-200 cycles of jelly body refinement first (if model is straight after MR)' ] )

        self.closeSubFrame()


    def drawAdvancedOptions( self ):
        indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'

        self.createLine( [ 'widget', 'SS', 'label', 'Save disk space by removing excessive ProSMART output' ] )
        self.createLine( [ 'label', 'Download and use homologues with resolution better than', 'widget', 'MINRES', 'label', 'Angstrom' ] ) # 'stretch'
        self.createLine( [ 'label', 'Use up to ', 'widget', 'NH', 'label', 'homologues for restraints' ] ) # 'stretch'
        self.createLine( [ 'label', 'Use up to ', 'widget', 'NC', 'label', 'chains to generate restraints' ] ) # 'stretch'

