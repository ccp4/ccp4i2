"""
    arcimboldo_gui.py
    Copyright (C) 2015 The University of York
    Author: Jon Agirre

"""

import os
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui.CCP4TaskWidget import CTaskWidget

class arcimboldo_gui(CTaskWidget):

    TASKNAME = 'arcimboldo'
    TASKVERSION = 0.1
    TASKMODULE=['alpha_fold', 'molecular_replacement']
    TASKTITLE='Ab initio phasing and chain tracing - ARCIMBOLDO (LITE, BORGES, SHREDDER)'
    SHORTTASKTITLE='Arcimboldo'
    DESCRIPTION = 'Structure solution from ideal molecular fragments using PHASER and SHELXE'


    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def customBor( self ):
        controlParameters = self.container.controlParameters
        return str(controlParameters.RUN_MODE) == 'CUSTOM' and str(controlParameters.ARCIMBOLDO_RUN) != 'multiprocessing'

    def checkBorgesLib(self):
        ccp4_master_home = os.environ.get ( "CCP4_MASTER", "not_set" )
        return os.path.isdir(os.path.join(ccp4_master_home,'BORGES_LIBS'))

    @QtCore.Slot()
    def changeCoiled(self):
        if self.container.controlParameters.COIL_COILED:
            self.container.controlParameters.TNCS.set(False)
            self.container.controlParameters.TNCS_T.set(False)
        else:
            self.container.controlParameters.TNCS.set(False)
            self.container.controlParameters.TNCS_T.set(True)

    @QtCore.Slot()
    def changeRunMode(self):
        if str(self.container.developerOptions.DEVELOPER_MODE) == 'EXISTING':
            self.container.developerOptions.EXISTING_FOLDER.setQualifiers({'allowUndefined':False,'mustExist':True})
        else:
            self.container.developerOptions.EXISTING_FOLDER.setQualifiers({'allowUndefined':True,'mustExist':False})

        self.getWidget('EXISTING_FOLDER').updateViewFromModel()
        self.validate()

    @QtCore.Slot()
    def changeFixedFragment(self):
        if self.container.controlParameters.LITE_PARTIAL and str(self.container.controlParameters.ARCIMBOLDO_OPTIONS) == 'LITE':
            self.container.inputData.LITE_FIXED.setQualifiers({'allowUndefined':False,'mustExist':True})
        else:
            self.container.inputData.LITE_FIXED.setQualifiers({'allowUndefined':True,'mustExist':False})

        self.getWidget('LITE_FIXED').updateViewFromModel()
        self.validate()


    def showSpherical(self):
        if str(self.container.controlParameters.ARCIMBOLDO_OPTIONS) == 'SHREDDER' and \
        str (self.container.controlParameters.SHREDDER_OPTIONS) == 'spherical':
            return True

        return False

    @QtCore.Slot()
    def changeMode(self):
        self.container.inputData.PDB_SHREDDER.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.container.inputData.BORGES_CUSTOM.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.container.inputData.PDB_LITE.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.container.inputData.LITE_CUSTOMS_LIST.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.container.inputData.LITE_HELICES_LIST.setQualifiers({'allowUndefined':True,'mustExist':False})

        if str(self.container.controlParameters.ARCIMBOLDO_OPTIONS) == 'LITE':
            if str(self.container.controlParameters.LITE_MODELS) == 'CUSTOM':
                self.container.inputData.PDB_LITE.setQualifiers({'allowUndefined':False,'mustExist':True})
            elif str(self.container.controlParameters.LITE_MODELS) == 'HELICES':
                self.container.inputData.LITE_HELICES_LIST.setQualifiers({'allowUndefined':False,'mustExist':True})
            elif str(self.container.controlParameters.LITE_MODELS) == 'CUSTOMS':
                self.container.inputData.LITE_CUSTOMS_LIST.setQualifiers({'allowUndefined':False,'mustExist':True})

        elif str(self.container.controlParameters.ARCIMBOLDO_OPTIONS) == 'BORGES':
            if str(self.container.controlParameters.BORGES_LIBRARY) == 'CUSTOM':
                self.container.inputData.BORGES_CUSTOM.setQualifiers({'allowUndefined':False,'mustExist':True})

        else: #SHREDDER
            self.container.inputData.PDB_SHREDDER.setQualifiers({'allowUndefined':False,'mustExist':True})

        self.getWidget('PDB_SHREDDER').updateViewFromModel()
        self.getWidget('BORGES_CUSTOM').updateViewFromModel()
        self.getWidget('PDB_LITE').updateViewFromModel()
        self.getWidget('LITE_CUSTOMS_LIST').updateViewFromModel()
        self.getWidget('LITE_HELICES_LIST').updateViewFromModel()
        self.validate()

    def drawContents(self):
        self.setProgramHelpFile ( 'arcimboldo' )

        self.openFolder ( folderFunction='inputData' )
        self.createLine ( [ 'tip', 'Choose ARCIMBOLDO program', 'label','Run ARCIMBOLDO','widget', 'ARCIMBOLDO_OPTIONS', 'tip', 'Define where to run ARCIMBOLDO', 'label', 'on', 'widget', 'ARCIMBOLDO_RUN'] )

        #Coil coiled mode just for Arcimboldo Lite
        self.createLine ( [ 'tip', 'Run in coil coiled mode', 'widget','COIL_COILED', 'label','Run in coil coiled mode'])
        self.createLine ( [ 'tip', 'Run in predicted model mode', 'widget','SHREDDER_PREDICTED', 'label','Run in predicted model mode'], toggleFunction=[ self.showSpherical, [ 'ARCIMBOLDO_OPTIONS', 'SHREDDER_OPTIONS' ] ] )
        
        self.createLine ( [ 'label','Grid configuration','widget','-guiMode','radio', 'RUN_MODE'], toggle=['ARCIMBOLDO_RUN','open',[ 'local_grid','remote_grid' ] ] )
        self.createLine ( [ 'tip', 'This file will be referred to in the job\'s bor-file in the instruction setup_bor_path = ...', 'label','Local path of the configuration file', 'widget','CONFIG_FILE'], toggleFunction=[self.customBor,['RUN_MODE','ARCIMBOLDO_RUN'] ] )

        self.createLine ( [ 'subtitle', 'Input data', 'Enter reflection data and specify fragment search' ] )
        self.openSubFrame ( frame=[True])
        self.createLine ( [ 'widget', 'F_SIGF' ] )
        self.createLine ( [ 'tip', 'Composition of the asymmetric unit', 'label', 'Asymmetric unit contains', 'widget', 'N_COMPONENTS', 'label', 'components of molecular weight', 'widget',
            'MOLECULAR_WEIGHT', 'label', 'Daltons'] )

        self.closeSubFrame()

        #Models Arcimboldo Lite
        self.createLine ( [ 'subtitle', 'Model', 'Enter the model information' ], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'LITE' ] ]  )
        self.openSubFrame ( frame=[True], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'LITE'] ]  )
        self.createLine ( [ 'tip', 'Define type of search models and its expected r.m.s.d. from target', 'label', 'Use', 'widget', 'LITE_MODELS', 'label', 'assuming rmsd from target', 'widget', 'LITE_RMSD', 'label', 'A'])
        self.createLine ( [ 'tip', 'Define the number of helical fragments to search for and the length of the fragment', 'label', 'Search for', 'widget', 'N_FRAGMENTS', 'label', 'copies of a helix containing', 'widget', 'HELIX_LENGTH', 'label', 'residues' ], toggle=['LITE_MODELS','open',[ 'HELIX'] ] )
        self.createLine ( [ 'tip', 'Define the number of custom fragments to search for','label', 'Search for', 'widget', 'N_FRAGMENTS', 'label', 'copies of custom model' ], toggle=['LITE_MODELS','open',[ 'CUSTOM'] ] )
        self.createLine ( [ 'widget', 'PDB_LITE' ], toggle=['LITE_MODELS','open',[ 'CUSTOM'] ] )
        self.createLine ( [ 'label', 'Search for the following models in specified order:' ], toggle=['LITE_MODELS','open',[ 'CUSTOMS','HELICES'] ] )
        self.createLine ( [ 'tip', 'Define length of an Ideal helical fragment', 'widget', 'LITE_HELICES_LIST' ], toggle=['LITE_MODELS','open',[ 'HELICES'] ] )
        self.createLine ( [ 'widget', 'LITE_CUSTOMS_LIST' ], toggle=['LITE_MODELS','open',[ 'CUSTOMS'] ] )
        self.createLine ( [ 'tip', 'Start from known partial structure', 'widget','LITE_PARTIAL', 'label','Start from known partial structure']   )
        self.createLine ( [ 'label', 'Fixed in', 'widget','LITE_FIXED'], toggle=['LITE_PARTIAL','open',[ True ] ]  ) 
        self.closeSubFrame()

        #Libraries Arcimboldo Borges
        self.createLine ( [ 'subtitle', 'Borges libraries and Model Handling Parameters', 'Enter the model information' ], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'BORGES' ] ]  )
        self.openSubFrame ( frame=[True], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'BORGES'] ]  )
        self.createLine ( [ 'tip','Orientations of secondary structure elements are shown as u (up) and d (down)' ,'label', 'Define library to be used with Arcimboldo Borges', 'widget', 'BORGES_LIBRARY', 'label', '(topology: u - up, d - down)'],toggleFunction=[self.checkBorgesLib,['ARCIMBOLDO_OPTIONS'] ] )
        self.createLine ( [ 'label', 'Library in', 'widget', 'BORGES_CUSTOM' ], toggle=['BORGES_LIBRARY','open',[ 'CUSTOM'] ])
        self.createLine ( [ 'label', 'Fragment refinement against experimental data:' ] )
        self.createLine ( [ 'tip', 'Use GYRE option when running Phaser\'s Rotation Function', 'widget','BORGES_GYRE', 'label','Switch Phaser GYRE option', 'widget', 'BORGES_GYRE_T']  )
        self.createLine ( [ 'tip', 'Use GIMBLE option when running Phaser\'s Rigid Body Refinement', 'widget','BORGES_GIMBLE', 'label','Switch Phaser GIMBLE option', 'widget', 'BORGES_GIMBLE_T']  )
        self.createLine ( [ 'tip', 'Use MULTICOPY mode', 'widget','BORGES_MULTICOPY', 'label','Switch MULTICOPY option', 'widget', 'BORGES_MULTICOPY_T']  )        
        self.closeSubFrame()

        #Shredder models
        self.createLine ( [ 'subtitle', 'Shredder models', 'Enter the model information' ], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'SHREDDER' ] ]  )
        self.openSubFrame ( frame=[True], toggle=['ARCIMBOLDO_OPTIONS','open',[ 'SHREDDER'] ]  )
        self.createLine ( [ 'widget', 'PDB_SHREDDER' ] )
        self.createLine ( [ 'widget','SHREDDER_RMSD', 'label', 'Similarity of PDB to the target structure: rmsd difference is', 'widget', 'SHREDDER_RMSD_T', 'label', 'A'])
        self.createLine ( [ 'tip', 'Convert input model to polyalanine','widget','SHREDDER_CONVERT', 'label','Convert to polyalanine', 'widget', 'SHREDDER_CONVERT_T']  )
        self.createLine ( [ 'tip', 'Make all B-factors equal', 'widget','SHREDDER_MAKE', 'label','Make all B-factors equal', 'widget', 'SHREDDER_MAKE_T']  )
        self.createLine ( [ 'tip', 'Define Shredder Mode', 'label','Shredder mode','widget', 'SHREDDER_OPTIONS'] )
        self.createLine ( [ 'tip', 'Maintain coil in the model', 'widget','SHREDDER_COIL', 'label','Maintain coil in the model', 'widget', 'SHREDDER_COIL_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )
        self.createLine ( [ 'tip', 'Use GYRE option when running Phaser\'s Rotation Function', 'widget','SHREDDER_GYRE', 'label','Perform gyre refinement', 'widget', 'SHREDDER_GYRE_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )
        self.createLine ( [ 'tip', 'Use GIMBLE option when running Phaser\'s Rigid Body Refinement', 'widget','SHREDDER_GIMBLE', 'label','Perform gimble refinement', 'widget', 'SHREDDER_GIMBLE_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )
        self.createLine ( [ 'tip', 'Perform Phaser\'s LLG-guided pruning', 'widget','SHREDDER_LLG', 'label','Perform LLG-guided pruning', 'widget', 'SHREDDER_LLG_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )
        self.createLine ( [ 'tip', 'Combine partial solutions using ALIXE - better phases but runs longer', 'widget','SHREDDER_COMBINE', 'label','Combine phases with alixe', 'widget', 'SHREDDER_COMBINE_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )
        self.createLine ( [ 'tip', 'Use MULTICOPY mode', 'widget','SHREDDER_MULTICOPY', 'label','Switch MULTICOPY option', 'widget', 'SHREDDER_MULTICOPY_T'], toggle=['SHREDDER_OPTIONS','open',[ 'spherical' ] ]  )        
        self.closeSubFrame()
        self.closeFolder ( )
        self.openFolder ( folderFunction='advancedData', title='Advanced data' )
        self.createLine ( [ 'tip','Number of amino acids in each fragment', 'widget','FRAGMENT_SIZE', 'label','Fragment size', 'widget', 'FRAGMENT_SIZE_T'], toggleFunction=[ self.showSpherical, [ 'ARCIMBOLDO_OPTIONS', 'SHREDDER_OPTIONS' ] ]  )
        self.createLine ( [ 'tip','Use TNCS option when running Phaser','widget','TNCS', 'label','Switch Phaser TNCS option', 'widget', 'TNCS_T']  )
        self.createLine ( [ 'tip','Replace default SHELXE options','label','shelxe_line =','widget','SHELXE_LINE'] )
        self.createLine ( [ 'tip', 'A line key = value to be added to bor-file', 'label','Add lines to bor-file'] )
        self.createLine ( [ 'tip','A line key = value to be added to bor-file', 'widget', '-guiMode', 'multiLine', 'KEYWORDS' ] )
        self.closeFolder ( )
        self.openFolder ( folderFunction='developerOptions', title='Developer options' )
        self.createLine ( [ 'tip', 'Run mode', 'label','Select run mode','widget', 'DEVELOPER_MODE'] )
        self.createLine ( [ 'label', 'Existing run', 'widget', 'EXISTING_FOLDER'], toggle=['DEVELOPER_MODE', 'open',[ 'EXISTING' ] ] )
        
        self.closeFolder ( )
        #Options (all)
        #self.createLine ( [ 'subtitle', 'Options', 'Enter the model information' ] )
        #self.openSubFrame ( frame=[True])
        #self.closeSubFrame()

        self.container.controlParameters.ARCIMBOLDO_OPTIONS.dataChanged.connect( self.changeMode )
        self.container.controlParameters.ARCIMBOLDO_OPTIONS.dataChanged.connect( self.changeFixedFragment )
        self.container.controlParameters.BORGES_LIBRARY.dataChanged.connect( self.changeMode )
        self.container.controlParameters.LITE_MODELS.dataChanged.connect( self.changeMode )
        self.changeMode()

        self.container.controlParameters.COIL_COILED.dataChanged.connect( self.changeCoiled )
        self.container.controlParameters.LITE_PARTIAL.dataChanged.connect ( self.changeFixedFragment )

        #Advanced options (all)
        self.container.developerOptions.DEVELOPER_MODE.dataChanged.connect( self.changeRunMode )

        #Developer (all)
