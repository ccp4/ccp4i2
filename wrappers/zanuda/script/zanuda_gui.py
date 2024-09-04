'''
    zanuda_gui.py: CCP4 GUI Project

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
'''

from qtgui.CCP4TaskWidget import CTaskWidget


class zanuda_gui(CTaskWidget):
    TASKMODULE = 'refinement'
    TASKTITLE = 'Zanuda'
    TASKNAME = 'zanuda'
    TASKVERSION = 0.1
    DESCRIPTION = 'Space group validation'
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['aimless_pipe', 'coot_rebuild']

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData')

        self.createLine(['subtitle', 'Reflection data'])
        self.openSubFrame(frame=True)
        self.createLine(['widget', 'F_SIGF'])
        self.createLine(['widget', 'FREERFLAG'])
        self.closeSubFrame()

        self.createLine(['subtitle', 'Starting model'])
        self.openSubFrame(frame=True)
        self.createLine(['widget', 'XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        self.closeSubFrame()

        self.createLine(['subtitle', 'Options'])
        self.openSubFrame(frame=True)
        self.createLine(
            [
                'widget',
                'AVERAGE',
                'label',
                'Symmetryse input model before further transformations'
            ]
        )
        self.closeSubFrame()

