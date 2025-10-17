"""
    modelASUCheck_gui.py: CCP4 GUI Project

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


class modelASUCheck_gui(CTaskWidget):
    TASKNAME = 'modelASUCheck'
    TASKVERSION = 0.1
    TASKMODULE = ['model_data_utility']
    SHORTTASKTITLE = 'Check model against AU contents'
    TASKTITLE = 'Check model against AU contents'
    DESCRIPTION = 'Align sequences in model to those in AU contents to validate model.'

    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)

        self.createLine(['subtitle', 'Input coordinates'])
        self.openSubFrame(frame=True)
        self.createLine (['tip', 'Input model', 'widget', 'XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()

        self.createLine(['subtitle', 'Input AU Contents'])
        self.openSubFrame(frame=True)
        self.createLine (['tip', 'Input AU Contents', 'widget', 'ASUIN'])
        self.closeSubFrame()
