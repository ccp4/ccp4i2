"""
    morda_i2_gui.py: CCP4 GUI Project
    
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

import os, sys
from qtgui.CCP4TaskWidget import CTaskWidget
from core.CCP4ErrorHandling import CErrorReport
try:
    from morda import version as mrd_version

except:
    mrd_version = None

class morda_i2_gui(CTaskWidget):

    TASKNAME = 'morda_i2'
    TASKVERSION = 0.0
    TASKMODULE = 'test' if sys.platform.startswith('win') else 'molecular_replacement'
    SHORTTASKTITLE = 'MORDA'
    TASKTITLE = 'Automated molecular replacement - MORDA'
    DESCRIPTION = 'Molecular Replacement with Domains and Assemblies'
    MGDISPLAYFILES = ['XYZOUT', 'FPHIOUT', 'DIFFPHIOUT']
    WHATNEXT = ['prosmart_refmac', 'shelxeMR', 'buccaneer_build_refine_mr', 'coot_rebuild']
    AUTOPOPULATEINPUT = True
    PROGRAMHELP = ['morda']
    RANK=1
    not_installed_msg = 'Press Help button in the taskbar above to find installation instructions'

    def drawContents(self):
        if not mrd_version:
            self.openFolder(folderFunction='inputData', title='Input Data')
            self.createLine(['subtitle', self.not_installed_msg])
            self.WHATNEXT = []
            return

        self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['subtitle', 'Reflection data'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'F_SIGF'])
        self.createLine(['widget', 'FREERFLAG'])
        self.createLine(['label', 'Check alternative space groups', 'widget', 'ALTSG'])
        self.closeSubFrame()
        self.createLine(['subtitle', 'Model preparation'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'ASUIN'])
        self.createLine(['label', 'Number of homologous structures to try:', 'widget', 'NSTRUCT'])
        self.closeSubFrame()
        self.createLine(['label', ''])
        self.createLine(['label', 'Number of CPUs to use:', 'widget', 'NCPU'])
        b1 = self.container.inputData.MOCK_DIR.isSet()
        b2 = self.container.inputData.MOCK_PAUSE.isSet()
        b3 = self.container.inputData.CHECK_PAUSE.isSet()
        if b1 or b2 or b3 or 'MRD_DEVEL' in os.environ:
            self.openFolder(folderFunction='controlParameters', title='Developer')
            tooltip = 'There is a bug if you see this section'
            self.createLine(['subtitle', 'Developer options', tooltip])
            self.openSubFrame(frame=[True])
            self.createLine(['label', 'Start from a copy of existing work directory:'])
            self.createLine(['widget', 'MOCK_DIR'])
            lab1 = 'Set for each mock search to take'
            self.createLine(['label', lab1, 'widget', 'MOCK_PAUSE', 'label', 'seconds'])
            lab1 = 'Check for finished searches each'
            self.createLine(['label', lab1, 'widget', 'CHECK_PAUSE', 'label', 'seconds'])
            self.closeSubFrame()

    def taskValidity(self):
        rv = CErrorReport()
        if not mrd_version:
            rv._reports.append({
                'class': type(self),
                'name': '',
                'label': '',
                'code': 0,
                'details': self.not_installed_msg,
                'description': 'not installed\n',
            })

        return rv        

