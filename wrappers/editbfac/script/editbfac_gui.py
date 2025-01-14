
from ....qtgui import CCP4TaskWidget


class editbfac_gui(CCP4TaskWidget.CTaskWidget):

# Subclass CTaskWidget to give specific task window
    TASKNAME = 'editbfac'
    TASKVERSION = 0.0
    TASKMODULE = ['alpha_fold', 'model_data_utility' ]
    TASKTITLE = 'Process Predicted Models'
    SHORTTASKTITLE = 'Process Predicted Models'
    DESCRIPTION = 'Process Predicted Models - automatically process predicted models'
    WHATNEXT = ['phaser_simple', 'phaser_pipeline', 'molrep_pipe']
    PROGRAMHELP = 'editbfac'
  
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)


    def drawContents(self):
        self.setProgramHelpFile('editbfac')
        # Remove the 'followFrom' widget cos there can be no preceeding jobs - (this may be wrong thing to do)
        # Main Input Tab
        folder = self.openFolder(folderFunction='inputData', title='Input Data', followFrom=False)
        self.createLine(['subtitle', 'Atomic model'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'XYZIN'])
        self.closeSubFrame()
        
        self.createLine(['subtitle', 'Options'])
        self.openSubFrame(frame=[True])
        self.createLine(['advice', 'Select B-factor treatment option - it is important this is set correctly' ])
        self.setMenuText('BTREATMENT', {'plddt': 'AlphaFold model - convert pLDDT scores to B-factors',
                                        'rmsd': 'RoseTTAFold model - convert rmsd estimates to B-factors',
                                        'b_value': 'Model is already using B-factors, rather than plddt or rmsd values' })
        self.createLine(['label', '         ', 'widget', '-guiMode', 'multiLineRadio', 'BTREATMENT']) 
        self.createLine(['advice', 'Options for residues and regions' ] )
        self.createLine(['widget', 'CONFCUT', 'label', 'Remove low confidence residues (lddt or rmsd) (recommended)'])
        self.createLine(['widget', 'COMPACTREG', 'label', 'Split model into compact regions (recommended)'])
        self.closeSubFrame()
        
        self.createLine( ['subtitle', 'Optional Files'])
        self.openSubFrame(frame=[True])
        self.createLine( ['advice', 'PAE file'])
        self.createLine( ['widget', 'PAEIN'])
        self.createLine( ['advice', 'Distance model file'])
        self.createLine( ['widget', 'XYZDISTMOD'])
        self.closeSubFrame()
        # Define Options Tab
        fladd = self.openFolder(folderFunction='addSettings', title='Additional Settings', followFrom=False)
        self.createLine(['subtitle', 'Options'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'MAXDOM', 'label', 'Maximum domains to obtain.'])
        self.createLine(['widget', 'DOMAINSIZE', 'label', 'Approximate size of domains to be found (Angstroms)'])
        self.createLine(['widget', 'MINDOML', 'label', 'Minimum domain length (residues)'])
        self.createLine(['widget', 'MAXFRACCL', 'label', 'Maximum fraction close'])
        self.createLine(['widget', 'MINSEQRESI', 'label', 'Minimum sequential residues'])
        self.createLine(['widget', 'MINREMSEQL', 'label', 'Minimum remainder sequence length'])
        self.createLine(['widget', 'MINLDDT', 'label', 'Minimum LDDT for removing residues'])
        self.createLine(['widget', 'MAXRMSD', 'label', 'Maximum RMSD for removing residues'])
        self.createLine(['widget', 'SUBMINB', 'label', 'Subtract the lowest B-value from all B-values'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'PAE File Options (if PAE matrix supplied)'])
        self.createLine(['widget', 'PAEPOWER', 'label', 'PAE power'])
        self.createLine(['widget', 'PAECUTOFF', 'label', 'PAE cutoff'])
        self.createLine(['widget', 'PAEGRAPHRES', 'label', 'PAE graph resolution'])
        self.closeSubFrame()
        self.createLine( [ 'subtitle', 'Distance Model Options (if suitable model provided)' ])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'WEIGHTCA', 'label', 'Weight by CA-CA distance (if distance_model supplied)'])
        self.createLine(['widget', 'DISTPOW', 'label', 'Distance power (for weighting by CA-CA distance)'])
        self.closeSubFrame()
