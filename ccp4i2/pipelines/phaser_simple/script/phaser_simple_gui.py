"""
    pipelines/phaser_simple/phaser_simple_gui.py
    Copyright (C) 2016 Newcastle University
    Author: Martin Noble
    
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.baselayer import QtCore
from ccp4i2.pipelines.phaser_pipeline.script import phaser_pipeline_gui

#-------------------------------------------------------------------
class phaser_simple_gui(phaser_pipeline_gui.phaser_pipeline_gui):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_simple'
    TASKVERSION = 0.1
    TASKMODULE='molecular_replacement'
    TASKTITLE='Basic Molecular Replacement - PHASER'
    SHORTTASKTITLE='Basic MR - PHASER'
    DESCRIPTION = '''Simple MR with optional refinement and rebuilding (Phaser)'''

    def __init__(self,*args,**kws):
        super(phaser_simple_gui,self).__init__(*args, **kws)
        
 
    def drawContents(self):
        self.setProgramHelpFile('phaser')
        self.drawPhaserFrontPage()
        self.openFolder(title='Simple options',drawFolder=self.drawSimpleOptionsFolder)
        self.openFolder(folderFunction='extraSteps', title='Extra steps',drawFolder=self.drawExtraStepsFolder)
        self.openFolder(folderFunction='keywords', title='Keywords', drawFolder=self.drawPhaserKeywordsFolder)
    
    #Override the phaser_pipeline front page drawer so as not to preent Ensemble List
    def drawPhaserFrontPage(self):
        self.openFolder(folderFunction='inputData')
        self.refineWarningLineWidget = self.createLine( [ 'label', '<b>N.B. Please be aware that the default options for this task have changed to run <br/>shift field refinement (<em>sheetbend</em>) and refinement (<em>refmac5</em>) after molecular replacement.<br/> You can revert to the old behaviour by turning them off in the "Additional steps" section.</b>'], toggleFunction=[self.showFreeR,['RUNSHEETBEND','RUNREFMAC']] )
        self.drawReflectionPanel()
        self.drawCompositionPanel()
        self.createLine(['subtitle','Search model'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','XYZIN','label','Copies:','widget','NCOPIES'])
        self.container.inputData.XYZIN.dataChanged.connect(self.XYZINchanged)
        id_rms = self.createLine( [ 'label', 'Similarity of ensemble to target:', 'widget', 'ID_RMS', 'label'] )
        self.createLine( [ 'label', '', 'widget', 'SEARCHSEQUENCEIDENTITY' ], toggle = ['ID_RMS', 'open', [ 'ID' ] ], appendLine=id_rms )
        self.createLine( [ 'label', '', 'widget', 'SEARCHRMS' ], toggle = ['ID_RMS', 'open', [ 'RMS' ] ], appendLine=id_rms )
       
         
        # self.createLine(['label','Sequence identity to target (in range 0.0-1.0)','stretch','widget','SEARCHSEQUENCEIDENTITY'])
        self.closeSubFrame()
        self.createLine(['subtitle','Already placed coordinates','widget','INPUT_FIXED'])
        self.openSubFrame(frame=True,toggle=['INPUT_FIXED','open',[True]])
        self.createLine(['widget','XYZIN_FIXED'])
        self.container.inputData.INPUT_FIXED.dataChanged.connect(self.INPUT_FIXEDchanged)
        #This line has to come after the creation of the XYZIN_FIXED widget ...Duh !
        self.INPUT_FIXEDchanged()
        id_rms = self.createLine( [ 'label', 'Similarity of fixed ensemble to target:', 'widget', 'FIXED_ID_RMS', 'label'] )
        self.createLine( [ 'label', '', 'widget', 'FIXEDSEQUENCEIDENTITY' ], toggle = ['FIXED_ID_RMS', 'open', [ 'ID' ] ], appendLine=id_rms )
        self.createLine( [ 'label', '', 'widget', 'FIXEDRMS' ], toggle = ['FIXED_ID_RMS', 'open', [ 'RMS' ] ], appendLine=id_rms )
    
    def drawReflectionPanel(self):
        self.createLine( ['subtitle', 'Reflections' ])
        self.openSubFrame(frame=True)
        self.drawReflections()
        self.closeSubFrame()
    
    def drawReflections(self):
        #A simple reflections panel, so that some choices can be moved to the "options" panel
        self.createLine ( [ 'widget', 'F_SIGF' ] )
        self.createLine (['label','Use Intensity (I) or amplitude (F) ML target','stretch','widget', 'F_OR_I'])
        self.container.inputData.F_SIGF.dataChanged.connect(self.F_SIGFchanged)
        self.F_SIGFchanged()
        self.createLine ( [ 'widget', 'FREERFLAG' ], toggleFunction=[self.showFreeR,['RUNSHEETBEND','RUNREFMAC']])
    
    def drawSimpleOptionsFolder(self):
        self.createLine( ['subtitle', 'Simple options' ])
        self.openSubFrame(frame=True)
        self.createLine (['label','Resolution','stretch','widget','RESOLUTION_LOW','widget','RESOLUTION_HIGH'])
        self.createLine (['label', 'Spacegroups','stretch','widget', 'SGALT_SELECT',] )
        self.drawSGChoice()
        self.closeSubFrame()
    
    def drawExtraStepsFolder( self ):
        self.createLine(['subtitle','Reference structure','If you know a structure in a similar crystal form, you may provide a reference structure,'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','-browseDb', True, 'XYZIN_TARGET'])
        self.closeSubFrame()
        self.createLine(['subtitle','Fit and fill in coot','Optionally run COOT to fill and fit partial residues'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNCOOT','label','Run COOT to fill and fit partial residues'])
        self.closeSubFrame()
        self.createLine(['subtitle','Run sheetbend and/or refmac','-toolTip', 'Run Sheetbend and/or REFMAC on top output solution'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNSHEETBEND','label','Run Sheetbend on top output solution'])
        self.createLine(['widget','RUNREFMAC','label','Run REFMAC on top output solution'])
        self.closeSubFrame()

    @QtCore.Slot()
    def INPUT_FIXEDchanged(self):
        if self.container.inputData.INPUT_FIXED:
            self.container.inputData.XYZIN_FIXED.setQualifiers({'allowUndefined':False,'mustExist':True})
        else:
            self.container.inputData.XYZIN_FIXED.setQualifiers({'allowUndefined':True,'mustExist':False})
        self.getWidget('XYZIN_FIXED').validate()
        
    @QtCore.Slot()
    def XYZINchanged(self):
        # check for REMARK PHASER cards in XYZIN
        cardFound = False
        if self.container.inputData.XYZIN.isSet():
            with open (self.container.inputData.XYZIN.__str__(),'r') as pdbFile:
                for line in pdbFile:
                    if "REMARK PHASER ENSEMBLE MODEL" in line:
                        cardFound = True
                        break
                    elif line[0:4] == "ATOM":
                        break
        if cardFound:
            self.container.inputData.ID_RMS.setQualifiers({'enumerators':['ID','RMS','CARD'], 'menuText':['sequence identity (in range 0.0-1.0)','RMS','read from header of PDB']})
            self.container.inputData.ID_RMS.set('CARD')
        else:
            self.container.inputData.ID_RMS.setQualifiers({'enumerators':['ID','RMS']})
            self.container.inputData.ID_RMS.set('ID')
        self.getWidget('ID_RMS').populateComboBox(self.container.inputData.ID_RMS)
        self.getWidget('ID_RMS').updateViewFromModel()

    def showFreeR(self):
        if self.container.inputData.RUNSHEETBEND or self.container.inputData.RUNREFMAC:
            return True
        return False

