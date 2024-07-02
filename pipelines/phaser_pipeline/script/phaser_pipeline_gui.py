"""
    pipelines/phaser_pipeline/wrappers/phaser_MR/script/phaser_MR_gui.py
    Copyright (C) 2014 Newcastle University
    Author: Martin Noble
    
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script import phaser_MR_AUTO_gui
from PySide6 import QtCore
from core import CCP4ErrorHandling

#-------------------------------------------------------------------
class phaser_pipeline_gui(phaser_MR_AUTO_gui.phaser_MR_AUTO_gui):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'phaser_pipeline'
    TASKVERSION = 0.1
    TASKMODULE='molecular_replacement'
    TASKTITLE='Expert Mode Molecular Replacement - PHASER'
    SHORTTASKTITLE = 'Expert MR - PHASER'
    DESCRIPTION = '''Advanced MR options followed by refinement and rebuilding (Phaser, Refmac5, Coot)'''
    RANK=1
    EXPORTMTZPARAMS = [ ['F_SIGF','I_SIGI'], 'MAPOUT_REFMAC' , ['DIFMAPOUT_REFMAC', 'diff'] ]

    def __init__(self,parent):
        super(phaser_pipeline_gui,self).__init__(parent)
        #Here deal with clones of legacy tasks where the logic of selecting I vs F was different POSSIBILITIES:
        '''if self.container.inputData.F_OR_I.isSet() and self.container.inputData.F_OR_I.__str__() == 'I' and self.container.inputData.I_SIGI.isSet():
            self.container.inputData.inputData.F_SIGF = self.container.inputData.inputData.I_SIGI
            self.container.inputData.inputData.I_SIGI.unSet()
        '''


    def setDefaultParameters(self):
        # Reimplement in tasks to set initial parameters
        # read ccp4/src/phaser/source/phaser/defaults for defaults
        pass
    
    def drawReflectionPanel(self):
        self.createLine( ['subtitle', 'Reflections' ])
        self.openSubFrame(frame=True)
        self.drawReflections()
        self.drawSGChoice()
        self.createLine(['tip','Free R flag will be used in the (optional) REFMAC step','widget','FREERFLAG'])
        self.closeSubFrame()
    
    def drawExtraStepsFolder(self):
        self.createLine(['subtitle','Reference structure','-toolTip','If you know a structure in a similar crystal form, you may provide a reference structure,'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','-browseDb', True, 'XYZIN_TARGET'])
        self.closeSubFrame()
        self.createLine(['subtitle','Fit and fill in coot','-toolTip','Optionally run COOT to fill and fit partial residues'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNCOOT','label','Run COOT to fill and fit partial residues'])
        self.closeSubFrame()
        self.createLine(['subtitle','Run sheetbend','-toolTip', 'Run Sheetbend on top output solution'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNSHEETBEND','label','Run Sheetbend on top output solution'])
        self.closeSubFrame()
        self.createLine(['subtitle','Run refmac','-toolTip', 'Run REFMAC on top output solution'])
        self.openSubFrame(frame=True)
        self.createLine(['widget','RUNREFMAC','label','Run REFMAC on top output solution'])
        self.closeSubFrame()
    
    def drawContents(self):
        self.setProgramHelpFile('phaser')
        self.drawPhaserFrontPage()
        self.openFolder(folderFunction='knownStructure', title='Known structure',drawFolder=self.drawPhaserKnownStructureFolder)
        self.openFolder(folderFunction='keywords', title='Keywords', drawFolder=self.drawPhaserKeywordsFolder)
        self.openFolder(folderFunction='extraSteps', title='Extra steps',drawFolder=self.drawExtraStepsFolder)
    
    def fix(self):
        # Avoid recording unused data in db
        self.container.inputData.I_SIGI.unSet()
        return CCP4ErrorHandling.CErrorReport()
