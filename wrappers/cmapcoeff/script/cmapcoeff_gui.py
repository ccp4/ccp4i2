"""
    cmapcoeff_gui.py
    Copyright (C) 2014-2019 Jon Agirre & University of York
    Author: Jon Agirre

"""

from ccp4i2.baselayer import QtGui, QtWidgets,QtCore
from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.core import CCP4XtalData

class cmapcoeff_gui(CTaskWidget):

    TASKNAME = 'cmapcoeff'
    TASKVERSION = 0.1
    TASKMODULE='expt_data_utility'
    TASKTITLE='Calculate unusual map coefficients'
    SHORTTASKTITLE='Calculate map coefficients'
    DESCRIPTION = 'Compute map coefficients from set of observations and phases (cmapcoeff)'

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def single_visible ( self ) :
        self.updateRequirements ( )
        if self.container.controlParameters.MAPTYPE == "fobsfobs" :
            return True
        else :
            return False

    def drawContents(self):
        self.setProgramHelpFile ( 'cmapcoeff' )

        self.openFolder ( folderFunction='inputData' )

        self.createLine ( [ 'subtitle', 'Type of map to create', 'Select one of the available map types. Only one set of map coefficients will be created.'] )
        self.openSubFrame ( frame=[True] )
        self.createLine ( [ 'widget', '-guiMode', 'radio', 'MAPTYPE' ] )
        self.closeSubFrame ( )
        self.createLine ( [ 'subtitle', 'First dataset', 'Anomalous differences can be computed for just one dataset.' ] )
        self.openSubFrame ( frame=[True] )
        self.createLine ( [ 'widget', 'F_SIGF1' ] )
        self.createLine ( [ 'widget', 'ABCD1' ] )
        self.closeSubFrame ( )

        #self.container.controlParameters.MAPTYPE.dataChanged.connect(self.updateRequirements )

        self.createLine ( [ 'subtitle', 'Second dataset', 'If you do not provide phases, those from the first dataset will be used.' ], toggleFunction=[ self.single_visible, [ 'MAPTYPE' ] ] )
        self.openSubFrame ( frame=[True], toggleFunction=[ self.single_visible, [ 'MAPTYPE' ] ] )
        self.createLine ( [ 'widget', 'F_SIGF2' ] )
        self.createLine ( [ 'widget', 'ABCD2' ] )
        self.closeSubFrame ( )

        self.createLine ( [ 'subtitle', 'Basic operations on map coefficients', 'All these operations are done in reciprocal space.' ] )
        self.openSubFrame ( frame = [True] )
        self.createLine ( [ 'tip', 'Use a negative B-factor to sharpen the map, or a positive one to blur it', 'label', 'Sharpen or blur the output map using an isotropic B-factor: ', 'widget', 'B_VALUE' ] )
        self.createLine ( [ 'widget', 'SCALE', 'label', 'scale the observations' ], toggleFunction=[ self.single_visible, [ 'MAPTYPE' ] ] )

        self.createLine ( [ 'spacing', 30, 'widget', '-guiMode', 'radio', 'F1_TO_F2' ], toggle=[ 'SCALE', 'open', [True] ] )
        self.createLine ( [ 'label', 'Choose a different resolution limit', 'widget', 'RESOLUTION' ] )
        self.closeSubFrame ( )

        self.createLine ( [ 'subtitle', 'Output options', 'Map files are not used within i2. Use options with care, defaults are right for most cases.' ] )
        self.openSubFrame ( frame = [True] )
        self.createLine ( [ 'widget', 'MAP_OUTPUT', 'label', 'produce a map file too' ] )
        self.createLine ( [ 'advice', 'Grid sampling should not be changed unless you know what you are doing.' ], toggle=[ 'MAP_OUTPUT', 'open', [True] ] )
        self.createLine ( [ 'label', 'u ', 'widget', 'INDEX_U', 'stretch', 'label', 'v ', 'widget', 'INDEX_V', 'stretch', 'label', 'w ', 'widget', 'INDEX_W' ], toggle=[ 'MAP_OUTPUT', 'open', [True] ] )
        self.closeSubFrame ( )

        self.closeFolder ( )

        self.updateRequirements ( )
#        self.updateViewFromModel ( )


    @QtCore.Slot()
    def updateRequirements ( self ) :

        if self.container.controlParameters.MAPTYPE == "anom" :
            self.container.inputData.F_SIGF1.setQualifier ( 'requiredContentFlag', [ 2 ] )
            self.container.inputData.F_SIGF2.setQualifier ( 'allowUndefined', True )
            self.getWidget ( 'F_SIGF2' ).validate ( )
            self.getWidget ( 'F_SIGF2' ).loadJobCombo ( )
            self.container.inputData.ABCD2.setQualifier ( 'allowUndefined', True )
            self.getWidget ( 'ABCD2' ).validate ( )
            self.getWidget ( 'ABCD2' ).loadJobCombo ( )
            self.container.inputData.F_SIGF1.setQualifier ( 'toolTip', 'Computing an anomalous difference maps requires the data to contain intensity or amplitude pairs' )
            self.container.inputData.F_SIGF1.loadFile ( )
            self.getWidget ( 'F_SIGF1' ).loadJobCombo ( )
            self.getWidget ( 'F_SIGF1' ).validate ( )
            self.container.controlParameters.SCALE = False
            self.getWidget ( 'SCALE' ).validate ( )

        elif self.container.controlParameters.MAPTYPE == "fobsfobs" :
            self.container.inputData.F_SIGF1.setQualifier ( 'requiredContentFlag', [ 1, 2, 3, 4 ] )
            self.container.inputData.F_SIGF2.setQualifier ( 'allowUndefined', False )
            self.container.inputData.F_SIGF1.loadFile ( )
            self.getWidget ( 'F_SIGF1' ).loadJobCombo ( )
            self.getWidget ( 'F_SIGF1' ).validate ( )
            self.container.controlParameters.SCALE = True
            self.getWidget ( 'SCALE' ).validate ( )

        else : # simple fobs map calculation
            self.container.inputData.F_SIGF1.setQualifier ( 'requiredContentFlag', [ 1, 2, 3, 4 ] )
            self.container.inputData.F_SIGF2.setQualifier ( 'allowUndefined', True  )
            self.container.inputData.ABCD1.setQualifier   ( 'allowUndefined', False )
            self.container.inputData.ABCD2.setQualifier   ( 'allowUndefined', True  )
            self.container.inputData.F_SIGF1.loadFile ( )
            self.getWidget ( 'F_SIGF1' ).loadJobCombo ( )
            self.getWidget ( 'F_SIGF1' ).validate ( )
            self.container.inputData.ABCD1.loadFile ( )
            self.getWidget ( 'ABCD1' ).loadJobCombo ( )
            self.getWidget ( 'ABCD1' ).validate ( )
            self.container.controlParameters.SCALE = False
            self.getWidget ( 'SCALE' ).validate ( )
        return
