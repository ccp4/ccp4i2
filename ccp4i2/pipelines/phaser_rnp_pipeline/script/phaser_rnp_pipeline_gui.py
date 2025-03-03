"""
    pipelines/phaser_rnp_pipeline/phaser_rnp_pipeline_gui.py
    Copyright (C) 2014 Newcastle University
    Author: Martin Noble
    
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from PySide2 import QtCore
from pipelines.phaser_pipeline.wrappers.phaser_MR_RNP.script import phaser_MR_RNP_gui

#-------------------------------------------------------------------
class phaser_rnp_pipeline_gui(phaser_MR_RNP_gui.phaser_MR_RNP_gui):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_rnp_pipeline'
    TASKVERSION = 0.1
    TASKMODULE='refinement'
    TASKTITLE="Rigid body refinement - PHASER"
    SHORTTASKTITLE="Rigid body PHASER"
    DESCRIPTION = '''Define rigid bodies for refinement (Phaser), fill partial residues (Coot) and refine (Refmac)'''
    
    def __init__(self,parent):
        super(phaser_rnp_pipeline_gui,self).__init__(parent)
    
    def setDefaultParameters(self):
        # Reimplement in tasks to set initial parameters
        # read ccp4/src/phaser/source/phaser/defaults for defaults
        pass
    
    def drawContents(self):
        self.setProgramHelpFile('phaser')
        self.drawFrontPage()
        
        self.drawPhaserKeywordsFolder()
        #Open a folder to deal with UI for extra steps
        self.openFolder(folderFunction='extraSteps', title='Extra steps')
        self.createLine(['advice',''])
        self.createLine(['advice','Optionally run COOT to fill and fit partial residues'])
        self.createLine(['widget','RUNCOOT','label','Run COOT to fit and fill partial residues'])
        self.createLine(['advice',''])
        self.createLine(['advice','Optionally run a simple 10-cycle refmac job'])
        self.createLine(['widget','RUNREFMAC','label','Run REFMAC on top output solution'])

    def drawFrontPage(self):
        self.openFolder(folderFunction='inputData')
        
        self.createLine( ['subtitle', 'Reflections' ])
        self.openSubFrame()
        self.createLine ( ['widget', 'F_SIGF' ] )
        self.createLine(['tip','Provide a FREERFLAG data object: useful for the optional REFMAC step (see "Extra steps" input folder) or if PHASER identifies a change in space group','widget','FREERFLAG'])
        self.createLine (['label','Resolution range','stretch','widget','RESOLUTION_LOW','widget','RESOLUTION_HIGH'])
        self.closeSubFrame()
        
        self.createLine( ['tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit','subtitle', 'Composition' ])
        self.openSubFrame(frame=True)
        self.setMenuText('COMP_BY',{ 'DEFAULT':'Provide no guidance to Phaser',
                         'MW':'Provide the mol. weight of protein and nucleic acids',
                         'ASU':'Provide a list of sequences (plus multiplicities)'
                         })
        self.createLine ( ['tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit','label', 'For estimating asymmetric unit contents:','widget','COMP_BY'])
        self.container.inputData.COMP_BY.dataChanged.connect(self.handleCOMP_BY)
        self.handleCOMP_BY()
        self.createLine ( ['label', '\n '], toggle=['COMP_BY','open',['DEFAULT']])
         
        self.createLine ( [ 'widget', '-title','Contents of crystal', 'ASUFILE' ], toggle=['COMP_BY','open',['ASU']])
        self.createLine ( [ 'label','Molecular weight (Da) of protein in the ASU','stretch','widget', 'ASU_PROTEIN_MW' ], toggle=['COMP_BY','open',['MW']])
        self.createLine ( [ 'label','Molecular weight (Da) of nucleic acid','stretch','widget', 'ASU_NUCLEICACID_MW' ], toggle=['COMP_BY','open',['MW']])
        self.closeSubFrame()
        
        self.createLine( ['tip','Phaser will tweak the position and orientation of the specified coordinate set, after breaking it down into the specified fragments','subtitle', 'Fragments' ])
        #Widgets to get input PDB coordinates and set of selection strings to apply
        self.openSubFrame()
        self.createLine ( [  'widget', '-browseDb','True','XYZIN_PARENT' ] )
        self.createLine ( [  'widget', '-title' , 'Atom selections for regions of protein to use as rigid bodies', 'SELECTIONS' ] )
        self.getWidget('SELECTIONS').setListVisible(visible=True)
        self.getWidget('SELECTIONS').setEditorVisible(visible=True)
        self.closeSubFrame()

