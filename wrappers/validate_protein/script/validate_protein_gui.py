"""
    validate_protein_gui.py
    Copyright (C) 2022
    Author: William Rochira and Jon Agirre
"""

from qtgui.CCP4TaskWidget import CTaskWidget


class validate_protein_gui(CTaskWidget):
    TASKNAME = 'validate_protein'
    TASKVERSION = 0.1
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
        self.closeSubFrame()

        self.createLine([ 'subtitle', 'Model data (comparison)', 'Specify coordinates and reflections data from a previous revision of the model to enable comparison in the Iris panel' ])
        self.openSubFrame(frame=[True])
        self.createLine([ 'widget', '-browseDb', True, 'XYZIN_2' ])
        self.createLine([ 'widget', '-browseDb', True, 'F_SIGF_2' ])
        self.closeSubFrame()

        self.createLine([ 'subtitle', 'Which validation checks are needed?', '' ])
        self.openSubFrame(frame=[True])
        self.createLine([ 'widget', 'DO_IRIS', 'label', 'Iris: interactive multimetric validation charts' ])
        self.createLine([ 'widget', 'DO_MOLPROBITY', 'label', 'MolProbity: checks including C-betas, side-chain flips, omega angles, and clashes' ])
        self.createLine([ 'widget', 'DO_BFACT', 'label', 'B-factors: average B-factor graphs and breakdown tables for publications' ])
        self.createLine([ 'widget', 'DO_RAMA', 'label', 'Ramachandran plots' ])
        self.closeSubFrame()

        self.closeFolder()
