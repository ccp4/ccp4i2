#=======================================================================================
#
#    CTaskPhaserSingleMR.py : CTaskPhaserSingleMR(CCP4TaskWidget.CTaskWidget)
#
#    Author  : Kyle Stevenson,STFC
#    Created : 14th April 2018, KJS
#
#    Gui Class for single MR with Phaser
#
#=======================================================================================

from qtgui import CCP4TaskWidget
from PySide2 import QtCore

class CTaskPhaserSingleMR(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'phaser_singleMR'
    TASKVERSION = 0.0
    TASKMODULE = 'molecular_replacement'
    TASKTITLE = 'Single Atom Molecular Replacement'
    SHORTTASKTITLE = 'Single Atom MR'
    DESCRIPTION = 'Perform Single Atom MR using Phaser'
    WHATNEXT = ['coot_rebuild', 'prosmart_refmac', 'modelcraft', 'arp_warp_classic']
    MGDISPLAYFILES = ['FPHIOUT']

    def __init__(self, parent):
        CCP4TaskWidget.CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.setProgramHelpFile('phaserSingleMR')
        folder = self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters')
        self.createLine(['subtitle', 'Select input data', 'Observed structure factors and protein model details are required'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'F_SIGF'])
        self.createLine(['widget', 'FREERFLAG'])
        self.createLine(['tip', 'Do something', 'widget', 'RESOL_ON', 'label', 'Use resolution range : ',
                         'widget', 'RESOL_LO', 'label', ' Angstroms to ','widget', 'RESOL_HI'])
        self.closeSubFrame()
        # Also need the ASU stuff in here.
        self.createLine(['tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit',
                         'subtitle', 'Composition'])
        self.openSubFrame(frame=True)
        self.setMenuText('COMP_BY', {'DEFAULT':'Use protein average solvent content',
                                     'MW':'Provided as molecular weights of protein and nucleic acid',
                                     'ASU':'Provided as full specification by sequence'})
        self.createLine(['tip','Phaser uses this for maximum likelihood and to suggest number of copies in the asymmetric unit',
                         'label', 'Composition of asymmetric unit:', 'stretch', 'widget', 'COMP_BY'])
        self.createLine(['label', '\n '], toggle=['COMP_BY', 'open', ['DEFAULT']])
        self.createLine(['widget','ASUFILE'], toggle=['COMP_BY', 'open', ['ASU']])
        self.createLine(['label','Molecular weight (Da) of protein in the ASU', 'stretch', 'widget', 'ASU_PROTEIN_MW' ],
                        toggle=['COMP_BY','open', ['MW']])
        self.createLine(['label','Molecular weight (Da) of nucleic acid', 'stretch', 'widget', 'ASU_NUCLEICACID_MW' ],
                        toggle=['COMP_BY', 'open', ['MW']])
        self.closeSubFrame()
        self.createLine(['tip', 'Define parameters for the single MR search & Log-likelihood completion', 'subtitle', 'Define Search'])
        self.openSubFrame()
        self.createLine(['tip', 'Type of atoms used for single MR',
                         'label', 'Search for','widget', 'SINGLE_ATOM_NUM', 'label', ' atom(s) of type ',
                         'widget', 'SINGLE_ATOM_TYPE'])
        self.createLine(['tip', 'Log-likelihood Gain Completion On/Off',
                         'label', 'Turn on Log-likelihood Completion', 'widget', 'LLG_COMPL_ON'])
        self.createLine(['tip', 'Log-likelihood Gain Completion', 'label', 'Log-likelihood Gain Map (LLG) Completion with',
                         'widget', 'LLG_COMPL_ATOMTYP', 'label', 'atoms'], toggle=['LLG_COMPL_ON', 'open', [True]])
        self.createLine(['tip', 'Log-likelihood Gain Completion', 'widget', 'LLG_COMPL_SIGCO_ON',
                         'label', 'LLG-Map : sigma cutoff for adding new atom sites','widget', 'LLG_COMPL_SIGCO'],
                        toggle=['LLG_COMPL_ON','open',[True]])
        self.createLine(['tip', 'Log-likelihood Gain Completion', 'widget', 'LLG_COMPL_ATMSEP_ON',
                         'label', 'LLG-Map : atomic separation distance cut-off','widget', 'LLG_COMPL_ATMSEP'],
                        toggle=['LLG_COMPL_ON', 'open', [True]])
        self.createLine(['tip', 'Log-likelihood Gain Completion', 'widget', 'LLG_COMPL_MAXCYC_ON',
                         'label', 'LLG-Map : max number of completion cycles to use', 'widget', 'LLG_COMPL_MAXCYC'],
                        toggle=['LLG_COMPL_ON', 'open', [True]])
        self.closeSubFrame()
        self.drawExpertTab()

    def drawExpertTab(self):
        folderexp = self.openFolder(folderFunction='expertInputData', title='Expert Settings')
        self.createLine(['subtitle', 'Expert Inputs', 'Expert settings for advanced use of Phaser'])
        self.createLine(['tip', 'Manually set options for packing criteria',
                         'widget', 'EXP_PACKCT_ON', 'label', 'Set packing criteria : ',
                         'widget', 'EXP_PACKCT_TYPE', 'label', ' at ', 'widget', 'EXP_PACKCT_AMT'])
        self.createLine(['tip', 'Use Quick Packing (not relevant if all solutions selected)',
                         'label', '\t','widget', 'EXP_QKPACK_ON', 'label', 'Use Quick Packing'],
                        toggle=['EXP_PACKCT_ON','open',[True]])
        self.createLine(['tip', 'Manually set translation peak selection',
                         'widget', 'EXP_TRAN_SRCHPK_ON', 'label', 'Translation search peak criteria : ',
                         'widget', 'EXP_TRAN_SRCHPK_TYPE', 'label', ' at ', 'widget', 'EXP_TRAN_SRCHPK_AMT'])
        self.createLine(['tip', 'Purge translational peaks',
                         'widget', 'EXP_PURGE_TRANPK_ON', 'label','Purge translation peaks with % cutoff set to ',
                         'widget', 'EXP_PURGE_TRANPK_PER','label',' but with maximum number ','widget', 'EXP_PURGE_TRANPK_NUM'])
        self.createLine(['tip', 'Set expected Least Likelihood Gain (LLG) criteria',
                         'widget', 'EXP_EXPLLG_CRT_ON', 'label','Use expected LLG criteria, where the targeted gain value is ',
                         'widget', 'EXP_EXPLLG_CRT'])
        self.createLine(['tip', 'Set resolution range', 'widget', 'EXP_RESRAN_HRREFINE_ON',
                         'label', 'Use resolution range : ','widget', 'EXP_RESRAN_HRREFINE_LO',
                         'label', ' Angstroms to ','widget', 'EXP_RESRAN_HRREFINE_HI'])
        self.createLine(['tip', 'Use translational non-crystallographic symmetry if possible',
                         'widget' , 'EXP_TRAN_NCS_ON', 'label', 'Use translational NCS if present'])
