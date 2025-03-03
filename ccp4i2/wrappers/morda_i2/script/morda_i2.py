"""
    morda_i2.py: CCP4 GUI Project
    
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

import os, sys, shutil, re
from core import CCP4PluginScript
from core import CCP4ErrorHandling
from core import CCP4XtalData

class morda_i2(CCP4PluginScript.CPluginScript):
    TASKNAME = 'morda_i2'
    TASKVERSION= 0.1
    MAINTAINER = 'andrey.lebedev@stfc.ac.uk'
    TASKCOMMAND = sys.executable
    PERFORMANCECLASS = 'CRefinementPerformance'
    INTERRUPTABLE = True
    INTERRUPTLABEL = 'Stop and keep current best solution'
    
    def processInputFiles(self):
        content_flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN
        self.hklin, error = self.makeHklin([['F_SIGF', content_flag], 'FREERFLAG'])
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return self.FAILED

        assert self.container.inputData.ASUIN.isSet()
        self.seqin_file = os.path.join(self.workDirectory,'SEQIN.fasta')
        self.container.inputData.ASUIN.writeFasta(self.seqin_file)
        return self.SUCCEEDED

    def makeCommandAndScript(self, **kw):
        self.appendCommandLine(['-m', 'morda'])
        self.appendCommandLine(['-f', self.hklin])
        if self.container.inputData.ALTSG:
            self.appendCommandLine(['-a'])

        self.appendCommandLine(['-s', self.seqin_file])
        self.appendCommandLine(['-n', str(self.container.inputData.NSTRUCT)])
        self.appendCommandLine(['-c', str(self.container.inputData.NCPU)])
        self.appendCommandLine(['-w', str(self.getWorkDirectory())])
        self.appendCommandLine(['--mp', '--no-stdout'])
        self.appendCommandLine(['-x', self.makeFileName('PROGRAMXML')])
        self.appendCommandLine(['-i', 'INTERRUPT'])
        self.appendCommandLine(['-e'])
        if self.container.inputData.MOCK_DIR.isSet():
            self.appendCommandLine(['-v'])
            self.appendCommandLine(['-z'])
            source = os.path.join(str(self.container.inputData.MOCK_DIR), 'morda')
            target = os.path.join(self.getWorkDirectory(), 'morda')
            shutil.copytree(source, target)

        if self.container.inputData.MOCK_PAUSE.isSet():
            self.appendCommandLine(['--mock-pause', float(self.container.inputData.MOCK_PAUSE)])

        if self.container.inputData.CHECK_PAUSE.isSet():
            self.appendCommandLine(['-p', float(self.container.inputData.CHECK_PAUSE)])

        return self.SUCCEEDED

    def processOutputFiles(self):
        xyzout = os.path.join(self.getWorkDirectory(), 'output', 'morda.pdb')
        hklout = os.path.join(self.getWorkDirectory(), 'output', 'morda.mtz')
        if os.path.isfile(xyzout) and os.path.isfile(hklout):
            grp1 = ['FPHIOUT', 'DIFFPHIOUT']
            grp2 = ['FWT,PHWT', 'DELFWT,PHDELWT']
            error = self.splitHklout(grp1, grp2, infile=hklout)
            if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                return self.FAILED

            try:
                shutil.copy(xyzout, str(self.container.outputData.XYZOUT.fullPath))
                with open(os.path.join(self.getWorkDirectory(), 'output', 'morda.res')) as istream:
                  rc, rf = istream.read().split()

                self.container.outputData.PERFORMANCE.RFactor = float(rc)
                self.container.outputData.PERFORMANCE.RFree = float(rf)

            except:
                return self.FAILED

            return self.SUCCEEDED

        elif os.path.isfile(os.path.join(self.getWorkDirectory(), 'INTERRUPT')):
            return self.INTERRUPTED

        else:
            return self.UNSATISFACTORY

