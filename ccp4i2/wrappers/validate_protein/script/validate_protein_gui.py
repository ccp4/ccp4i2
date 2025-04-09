"""
Copyright (C) 2022
Author: William Rochira and Jon Agirre
"""

from PySide2 import QtCore

from ....qtgui.CCP4TaskWidget import CTaskWidget


class validate_protein_gui(CTaskWidget):
    TASKNAME = 'validate_protein'
    TASKVERSION = 0.2
    TASKMODULE='validation'
    TASKTITLE='Multimetric model validation - Iris'
    SHORTTASKTITLE='Multimetric validation'
    DESCRIPTION = 'Calculates per-residue metrics including B-factors, density fit quality, Ramachandran plots, rotamer outliers, and clashes (Iris-Validation & MolProbity)'
    RANK = 1 # This means the task's entry will carry a gears badge, signifying that it is a pipeline

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile('validate_protein')

        self.openFolder(folderFunction='inputData')

        self.createLine([ 'subtitle', 'Model data', 'Specify coordinates and reflections data of the model to be validated' ])
        self.openSubFrame(frame=[True])
        self.createLine([ 'widget', '-browseDb', True, 'XYZIN_1' ])
        self.createLine([ 'widget', '-browseDb', True, 'F_SIGF_1' ])
        self.createLine([ 'label', 'Label dataset 1 as ', 'widget', 'NAME_1'])
        self.closeSubFrame()

        self.createLine( [ 'widget', 'TWO_DATASETS', 'label', 'Compare against a second dataset, e.g. before and after model building/refinement'])
        
        self.openSubFrame ( frame=[True], toggle=[ 'TWO_DATASETS', 'open', [True]] )
        self.createLine([ 'widget', '-browseDb', True, 'XYZIN_2' ])
        self.createLine([ 'widget', '-browseDb', True, 'F_SIGF_2' ])
        self.createLine([ 'label', 'Label dataset 2 as ', 'widget', 'NAME_2'])
        self.closeSubFrame()

        self.container.controlParameters.TWO_DATASETS.dataChanged.connect(self.updateRequirements )

        self.createLine([ 'subtitle', 'Which validation and analysis tools should be run?', '' ])
        self.openSubFrame(frame=[True])
        self.createLine([ 'widget', 'DO_IRIS', 'label', 'Iris: interactive multimetric validation charts' ])
        self.createLine([ 'widget', 'DO_MOLPROBITY', 'label', 'MolProbity: checks including C-betas, side-chain flips, omega angles, and clashes' ])
        self.createLine([ 'widget', 'DO_TORTOIZE', 'label', 'Tortoize: calculate conformation-dependent Z-scores for backbone geometry' ])
        self.createLine([ 'widget', 'DO_BFACT', 'label', 'B-factors: average B-factor graphs and breakdown tables for publications' ])
        self.createLine([ 'widget', 'DO_RAMA', 'label', 'Ramachandran plots' ])
        self.closeSubFrame()

        self.closeFolder()

    @QtCore.Slot()
    def updateRequirements ( self ) :
        if self.container.controlParameters.TWO_DATASETS :
            self.container.inputData.XYZIN_2.setQualifier ( 'allowUndefined', False )
        else :
            self.container.inputData.XYZIN_2.setQualifier ( 'allowUndefined', True )
