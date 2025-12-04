#=======================================================================================
#
#    CTaskShelxeMR.py : CTaskShelxeMR(CCP4TaskWidget)
#    
#    Author  : Kyle Stevenson,STFC
#    Created : 16th Sep. 2015, KJS
#
#    Gui Class for refinement of MR solutions using Shelxe
#
#=======================================================================================


from baselayer import QtGui, QtWidgets,QtCore
from qtgui import CCP4TaskWidget
from qtgui import CCP4Widgets

class CTaskAcorn(CCP4TaskWidget.CTaskWidget):
    
    TASKNAME    = 'acorn'
    TASKVERSION = 1.0
    TASKMODULE  ='density_modification'
    
    TASKTITLE       = "ACORN - Phase Refinement with Dynamic Density Modification"
    SHORTTASKTITLE  = "ACORN"
    
    DESCRIPTION     = "Un-biased improvement of initial phases for high resolution data (1.5 Angstoms and better)"
    ERROR_CODES = {  200 : { 'description' : 'Space group of reflection and phases file does not match' } } 
    WHATNEXT = ['coot_rebuild','prosmart_refmac','modelcraft'] 
    MGDISPLAYFILES  = ['FPHIOUT']


    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        import functools
        self.setProgramHelpFile('acorn')  # Comes in two varieties; version 1 needs the reflections & the model; version 2 also needs the phase info.
        
        folder1 = self.openFolder(folderFunction='inputData',title='Input Data')
        self.createLine( [ 'label', 'Run ACORN with ', 'widget', '-guiMode', 'radio', 'ACORN_PHSIN_TYPE' ] )
        self.openSubFrame( frame=[True] )
        self.createLine(['subtitle', 'Reflection Data','Integrated reflection data for experimental phasing.'])
        self.createLine( [ 'widget', 'F_SIGF' ] )
        self.createLine( [ 'widget', 'ABCD' ], toggle=[ 'ACORN_PHSIN_TYPE', 'open', 'phases'] )
        self.closeSubFrame()
    
        self.openSubFrame( frame=[True], toggle=[ 'ACORN_PHSIN_TYPE', 'open', 'model'] )
        self.createLine(['subtitle', 'Model for approximate co-ordinates','Model for approximate co-ordinates (if known)'])
        self.createLine( ['widget', 'XYZIN'] )
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()
                
        folder2 = self.openFolder(folderFunction='acornGenParam',title='Advanced Acorn Parameters')
        self.createLine( ['subtitle', 'General ACORN phase improvement parameters','General Parameter Choice'] )
        self.createLine( [ 'label', 'Number of trials: ', 'widget' ,'ACOPH_TRIALS' ] )
        countObject = getattr(self.container.controlParameters, 'ACOPH_TRIALS')
        countObject.dataChanged.connect( functools.partial(self.handleNObjectsChanged))
        self.createLine( [ 'label', 'Define custom parameters for each trial', 'widget', 'ACOPH_CUSTOM' ] )
        self.openSubFrame( frame=[True], toggle=[ 'ACOPH_CUSTOM', 'open', [True]] )
        self.createLine( [ 'label', 'Define DDM type for each trial', 'widget', 'ACOPH_CUSTDDM' ] )
        customDDM = getattr(self.container.controlParameters, 'ACOPH_CUSTDDM')
        customDDM.dataChanged.connect( self.changeDDM)
        for i in range(10):
            toggleValues = [10-j for j in range(10-i)]
            trial_line = self.createLine( [ 'tip', 'NCDDM keyword', 'label', 'Trial %s ' % str(i+1), 'widget', 'ACOPH_NCDDM_%s' % str(i+1) ],toggle = ['ACOPH_TRIALS','open', toggleValues] )
            self.createLine( [ 'label', ' cycles of ', 'widget', 'ACOPH_DDMK_%s' % str(i+1) ], appendLine=trial_line )
            self.createLine( [ 'label', ' ', 'widget', 'ACOPH_REFINE_%s' % str(i+1) ], appendLine=trial_line )     
        self.handleNObjectsChanged
        self.closeSubFrame()
        
        # for trial in range(1, self.container.controlParameters.ACOPH_TRIALS.__int__() + 1 ):
            # trial_line = self.createLine( [ 'tip', 'NCDDM keyword', 'label', 'Trial %d ' % trial, 'widget', 'ACOPH_NCDDM_%d' % trial ],toggle=[ toggleParam[trial - 1],'close', ['close'] ] )
            # self.createLine( [ 'label', ' cycles of ', 'widget', 'ACOPH_DDMK_%d' % trial ], appendLine=trial_line )
            # self.createLine( [ 'label', ' ', 'widget', 'ACOPH_REFINE_%d' % trial ], appendLine=trial_line )
        self.createLine( [ 'tip', 'PSFINISH keyword', 'label', 'Cease DDM cycling if phase shift between consecutive cycles is less than','widget', 'ACOPH_PSFINISH' ] )
        
        self.closeSubFrame()
    
        #folder3 = self.openFolder(folderFunction='acornReflectParam',title='Selection of Reflection Data')
        self.createLine(['subtitle', 'Selection of Reflection Data','Reflection Data Settings'])
        self.createLine( [ 'tip', 'Use reflections within resolution range if ticked, all data used otherwise', 'widget', 'ACORN_BRESOL', 'label', 'User defined resolution range for reflection data' ] )
        self.openSubFrame( frame=[True], toggle=[ 'ACORN_BRESOL', 'open', [True] ] )
        self.createLine( [ 'tip', 'RESOLUTION keyword, Low Resolution Limit', 'label','Use reflections within resolution limit of ', 'widget', 'ACOREF_RESOLL', 'label', 'Angstrom' ] )
        self.createLine( [ 'tip', 'RESOLUTION keyword, High Resolution Limit', 'label','to a limit of', 'widget', 'ACOREF_RESOLU', 'label', 'Angstrom' ] )
        self.closeSubFrame()
        
        self.createLine( [ 'tip', 'The reflections with FP less than cut*SIGFP will be rejected and treated as extended reflections', 'widget', 'ACORN_BEXCLUDE', 'label', 'Exclude low SIGFP reflections' ] )
        self.openSubFrame( frame=[True], toggle=[ 'ACORN_BEXCLUDE', 'open', [True] ] )
        self.createLine( [ 'tip', 'EXCLUDE keyword', 'label','Reject reflections that are less than ','widget', 'ACOREF_EXCLUDE' ,'label', '* the standard deviation (sigma)' ] )
        self.closeSubFrame()
        self.createLine( [ 'tip', 'The reflections with E-values greater than a certain cut will be rejected', 'widget', 'ACORN_BECUT', 'label', 'Exclude reflections that are high E-value outliers' ] )
        self.openSubFrame( frame=[True], toggle=[ 'ACORN_BECUT', 'open', [True] ] )
        self.createLine( [ 'tip', 'ECUT keyword', 'label','Reject observed reflections with E-values greater than ','widget', 'ACOREF_ECUT' ] )
        self.closeSubFrame()
        
        #folder4 = self.openFolder(folderFunction='acornPhasePar',title='Acorn Phase Parameters')
        self.createLine(['subtitle', 'Advanced Settings','Advanced Settings'])
        self.createLine( [ 'tip', 'CUTDDM keyword', 'label', 'Upper density limit for Dynamic Density Modification (DDM)','widget', 'ACOPH_CUTDDM' ] )
        self.createLine( [ 'tip', 'Custom ', 'widget', 'ACORN_BGRID', 'label', 'Choose a user defined grid size' ] )
        self.openSubFrame( frame=[True], toggle=[ 'ACORN_BGRID', 'open', [True] ] )
        self.createLine( [ 'tip', 'GRID keyword', 'label','Grid Size (Angstroms)','widget', 'ACOGEN_GRID' ] )
        self.closeSubFrame()
        # only applies to ACORN-MR
        # self.createLine( [ 'tip', 'Info', 'widget', 'ACORN_BSEED', 'label', 'Choose a user defined random number seed for ACORN' ] )
        # self.openSubFrame( frame=[True], toggle=[ 'ACORN_BSEED', 'open', [True] ] )
        # self.createLine( [ 'tip', 'SEED keyword', 'label','Random Number Seed','widget', 'ACOGEN_SEED' ] )
        
        
    @QtCore.Slot()
    def handleNObjectsChanged(self):
        # newNObjects = getattr(self.container.controlParameters, 'ACOPH_TRIALS').__int__()
        self.validate()
    
    @QtCore.Slot()
    def changeDDM(self):
        if self.container.controlParameters.ACOPH_CUSTDDM:
            for i in range(10):
                getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ).setQualifiers({'enumerators':['DDM0','DDM1','DDM2'], 'menuText':['DDM0','DDM1','DDM2']})
                getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ).set('DDM0')
                self.getWidget('ACOPH_DDMK_%s' % str(i+1)).populateComboBox(getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ))
                self.getWidget('ACOPH_DDMK_%s' % str(i+1)).updateViewFromModel()
        else:
            for i in range(10):
                getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ).setQualifiers({'enumerators':['DDM'], 'menuText':['DDM']})
                getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ).set('DDM')
                self.getWidget('ACOPH_DDMK_%s' % str(i+1)).populateComboBox(getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1) ))
                self.getWidget('ACOPH_DDMK_%s' % str(i+1)).updateViewFromModel()
        
                
        
    def taskValidity(self):
        from core import CCP4ErrorHandling
        rv = CCP4ErrorHandling.CErrorReport()
        # Check the space group is same in both input Mini-MTZ files 
        if self.container.controlParameters.ACORN_PHSIN_TYPE == "phases":
            if self.container.inputData.F_SIGF.exists():
                fsig_sg = str(self.container.inputData.F_SIGF.fileContent.spaceGroup)
            else:
                fsig_sg = None
            if self.container.inputData.ABCD.exists():
                phifom_sg = str(self.container.inputData.ABCD.fileContent.spaceGroup)
            else:
                phifom_sg = None
            if fsig_sg is not None and phifom_sg is not None and fsig_sg != phifom_sg:
                rv.append(self.__class__,200,details='Reflections file space group:'+ fsig_sg +' Phases file space group:'+ phifom_sg ,stack=False)
        else:
          #unset any input phases due to "use data from"
            self.container.inputData.ABCD.unSet()
        return rv        
