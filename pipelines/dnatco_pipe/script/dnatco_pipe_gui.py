"""
    dnatco_pipe_gui.py: CCP4 GUI Project
    Copyright (C) 2025 MRC-LMB
    Author: Martin Maly
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class dnatco_pipe_gui(CTaskWidget):
#-------------------------------------------------------------------

    TASKNAME = 'dnatco_pipe'
    TASKVERSION = 0.1
    TASKMODULE = ['refinement']
    TASKTITLE = 'DNATCO'
    SHORTTASKTITLE = 'DNATCO'
    DESCRIPTION = 'Restrain and validate nucleic acid structures'
    WHATNEXT = ['servalcat_pipe']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):

        # self.setProgramHelpFile('dnatco_pipe')

        self.openFolder(folderFunction="inputData")

        self.openSubFrame(frame=[True], title="Input data")
        self.autoGenerate(
           self.container.inputData,
           selection={"includeParameters": ["XYZIN1"]},
        )
        self.openSubFrame(frame=[True], toggle=['controlParameters.TOGGLE_XYZIN2', 'open', [False]])
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["GENERATE_RESTRAINTS"]},
        )
        self.closeSubFrame()
        self.closeSubFrame()

        self.createLine(['widget', 'TOGGLE_XYZIN2', 'label', 'Compare with another structure model'])
        self.openSubFrame(frame=[True], toggle=['controlParameters.TOGGLE_XYZIN2', 'open', [True]])
        self.autoGenerate(
           self.container.inputData,
           selection={"includeParameters": ["XYZIN2"]},
        )
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["GENERATE_RESTRAINTS"]},
        )
        self.closeSubFrame()

        self.openSubFrame(frame=[True], title="Parameters for restraints generation", toggle=['controlParameters.GENERATE_RESTRAINTS', 'open', [True]])
        self.autoGenerate(
            self.container.controlParameters,
            selection={"includeParameters": ["MAX_RMSD", "RESTRAINTS_SIGMA"]},
        )
        self.closeSubFrame()
        self.closeFolder()

