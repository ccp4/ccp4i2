#=======================================================================================
#
#    CTaskShelxeMR.py : CTaskShelxeMR(CCP4TaskWidget.CTaskWidget)
#    Copyright (C) 2015 STFC
#    Author  : Kyle Stevenson,STFC
#    Created : 14th April 2016, KJS
#
#    Gui Class for refinement of MR solutions using Shelxe
#
#=======================================================================================

from ....qtgui import CCP4TaskWidget


class CTaskShelxeMR(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'shelxeMR'
    TASKVERSION = 0.0
    TASKMODULE ='model_building'
    TASKTITLE = 'Model building from Molecular Replacement solution using Shelxe'
    SHORTTASKTITLE = "SHELXE-MR"
    DESCRIPTION = 'Use Shelxe to attempt to improve (or verify) a solution from Molecular Replacement'
    WHATNEXT = ['coot_rebuild', 'prosmart_refmac', 'modelcraft', 'arp_warp_classic']
    MGDISPLAYFILES = ['FPHIOUT']

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile('shelxeMR')
        folder = self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters')
        self.createLine(['subtitle', 'Select input data', 'Observed structure factors and protein model details are required'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'F_SIGF' ] )
        self.createLine(['widget', 'FREERFLAG'] )
        self.closeSubFrame()
        self.openSubFrame(frame=[True] )
        self.createLine(['subtitle', 'Model used for Molecular Replacement', 'Model used to evaluate phases during Molecular Replacement calculations.'])
        self.createLine(['widget', 'XYZIN'] )
        self.closeSubFrame()
        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Solvent Content', 'The amount of solvent - estimated or measured - in the asymmetric unit is required for the calculation.' ] )
        self.createLine(['label', 'Fraction of solvent in the asymmetric unit:', 'stretch', 'widget', 'FSOLVENT' ])
        self.closeSubFrame()
        self.closeFolder()
        folder2 = self.openFolder(folderFunction='controlParameters', title='Run Options')
        self.openSubFrame(frame=[True])
        self.createLine(['subtitle', 'Shelxe Run Options', 'Recommended to run with defaults.'])
        self.createLine(['label', 'Number of shelxe tracing cycles:', 'stretch', 'widget', 'NTCYCLES'])
        self.createLine(['label', 'Number of density modification cycles (per trace cycle):', 'stretch', 'widget', 'NMCYCLES'])
        self.createLine(['label', 'Time factor for helix and peptide search (Recommended setting: 1-4', 'stretch', 'widget', 'TIMEFAC'])
        self.createLine(['label', ', where 1 is the quickest, up to 4 with increasingly thorough, but slower, searches)'])
        self.createLine(['label', 'Search for Alpha Helices', 'stretch', 'widget', 'SALPHELICE'])
        self.createLine(['label', 'Search for Beta Sheets (parallel)', 'stretch', 'widget', 'SBETA'])
        self.createLine(['label', 'Search for Beta Sheets (anti-parallel)', 'stretch', 'widget', 'SANTIBETA'])
        self.createLine(['label', 'Optimize correlation coefficient for input model:', 'stretch', 'widget', 'PRUNRES'])
        self.createLine(['label', 'Apply n-fold symmetry (NC) to C-alpha traces:', 'stretch', 'widget', 'USENFOLD'])
        self.closeSubFrame()
        self.closeFolder()
