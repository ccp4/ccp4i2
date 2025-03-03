from qtgui import CCP4TaskWidget

class CTaskpairef(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'pairef'
    TASKVERSION = 1.1
    TASKMODULE =['refinement']
    TASKTITLE = 'Pairef'
    SHORTTASKTITLE = "Pairef"
    DESCRIPTION = 'Paired Refinement with Pairef'
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        tab = '&nbsp;&nbsp;&nbsp;&nbsp;'
        self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['label', 'Run Pairef with ', 'widget', '-guiMode', 'radio', 'SH_TYPE' ] )
        self.createLine(['widget', 'USE_PREREF', 'label', 'Use Pre-Refinement of model (before main paired refinement)'])
        self.openSubFrame(frame=[True])
        self.createLine(['advice', "Advice: to specify reflections, import intensities using the 'Import merged' task before a pairef run"])
        self.createLine(['tip', 'Input Reflections', 'widget', 'F_SIGF'])
        self.createLine(['tip', 'FreeR Flags', 'widget', 'FREERFLAG'])
        self.createLine(['tip', 'Model File', 'widget', 'XYZIN'])
        self.createLine(['tip', 'Unmerged Reflections (optional)', 'widget', 'UNMERGED'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['tip', 'Ligand dictionary', 'widget', 'DICT'])
        self.createLine(['tip', 'Command file for REFMAC5', 'widget', 'REFMAC_KEYWORD_FILE'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True], toggle=['SH_TYPE', 'open', 'semi'])
        self.createLine(['advice', 'Define Resolution Shells (semi-automated)'])
        self.createLine(['label', 'Add ', 'widget', 'NSHELL', 'label', 'resolution shells of width ', 'widget', 'WSHELL', 'label', 'Angstroms'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True], toggle=['SH_TYPE', 'open', 'manual'])
        self.createLine(['advice', 'Define Resolution Shells (manually)'])
        self.createLine(['label', 'Explicitly define shells', 'widget', 'MANSHELL', 'label', 'using comma separated format (e.g. 1.6,1.5,1.4,1.3)'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True], toggle=['USE_PREREF', 'open', [True]])
        self.createLine(['advice', 'Pre-Refinement Options'])
        self.createLine(['label', 'Number of pre-refinement cycles is', 'widget', 'NPRECYCLES'])
        self.createLine(['widget', 'USE_SHAKE', 'label', 'Shake model co-ordinates on/off'])
        self.createLine(['label', tab,'label', 'Shake/Randomize the co-ordinate sites using a given mean error of ', 'widget', 'SHAKE'], toggle=['USE_SHAKE', 'open', [True]])
        self.createLine([ 'widget', 'RESETBFAC', 'label', 'Reset B-factors to mean value'])
        self.createLine(['widget', 'COMPLETE', 'label', 'Perform complete cross-validation, using all available free reference sets.'])
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['advice', 'Pairef Options'])
        self.createLine(['label', 'Manually set the initial resolution (use if required)', 'widget', 'INIRES'])
        self.createLine(['label', 'Number of refinement cycles to perform', 'widget', 'NCYCLES'])
        self.createLine(['widget', 'AUTO_WGT', 'label', 'Use automatic weighting'])
        self.createLine(['label', tab, 'label', 'Set weighting factor', 'widget', 'WGT_TRM'], toggle=['AUTO_WGT', 'open', [False]])
        self.createLine(['widget', 'FIXED_TLS', 'label', 'Input fixed TLS parameters'])
        self.createLine(['label', tab, 'widget', '-browseDb', True, 'TLSIN', 'tip', 'Input Fixed TLS Parameters'], toggle=['FIXED_TLS', 'open', [True]])
        self.createLine(['label', tab, 'label', 'Number of TLS Refinement Cycles', 'widget', 'TLSCYC'], toggle=['FIXED_TLS', 'open', [True]])
        self.closeSubFrame()
