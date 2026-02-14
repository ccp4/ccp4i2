import os
import shutil
import sys

from ccp4i2.core import CCP4ErrorHandling, CCP4PluginScript, CCP4XtalData


class morda_i2(CCP4PluginScript.CPluginScript):
    TASKNAME = 'morda_i2'
    TASKCOMMAND = sys.executable
    PERFORMANCECLASS = 'CRefinementPerformance'
    
    def processInputFiles(self):
        content_flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN
        self.hklin, error = self.makeHklin([['F_SIGF', content_flag], 'FREERFLAG'])
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return self.FAILED

        assert self.container.inputData.ASUIN.isSet()
        self.seqin_file = os.path.join(self.workDirectory,'SEQIN.fasta')
        self.container.inputData.ASUIN.writeFasta(self.seqin_file)
        return self.SUCCEEDED

    def makeCommandAndScript(self):
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

