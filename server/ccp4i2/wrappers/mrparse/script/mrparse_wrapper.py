import json
import os
import shutil

from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class mrparse(CPluginScript):
    TASKTITLE = 'mrparse'
    TASKNAME = 'mrparse'
    TASKCOMMAND = 'mrparse'
    WHATNEXT = ['phaser_simple', 'phaser_pipeline', 'molrep_pipe']
    PERFORMANCECLASS = 'CExpPhasPerformance'

    def __init__(self, *args, **kwargs):
        self.seqin = None
        self.hklin = None
        CPluginScript.__init__(self, *args, **kwargs)

    def processInputFiles(self):
        self.seqin = self.container.inputData.SEQIN
        if self.container.inputData.F_SIGF.isSet():
            self.hklin, error = self.makeHklin([['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
            if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        pdb_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "homologs.json")
        af_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "af_models.json")
        esm_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "esm_models.json")

        if os.path.exists(pdb_json):
            with open(pdb_json, 'r') as f:
                pdb_data = json.load(f)
                if self.hklin:
                    pdb_data = sorted(pdb_data, key=lambda k: k['ellg'], reverse=True)
                else:
                    pdb_data = sorted(pdb_data, key=lambda k: k['seq_ident'], reverse=True)
                for i in pdb_data:
                    if i['pdb_file'] is not None:
                        xyz_in = os.path.join(self.getWorkDirectory(), "mrparse_0", i['pdb_file'])
                        xyz_out = os.path.join(self.getWorkDirectory(), os.path.basename(i['pdb_file']))
                    if os.path.isfile(xyz_in):
                        shutil.copy(xyz_in, xyz_out)
                    self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                    self.container.outputData.XYZOUT[-1].setFullPath(xyz_out)
                    self.container.outputData.XYZOUT[-1].annotation = "PDB hit: {}".format(i['name'])

        if os.path.exists(af_json):
            with open(af_json, 'r') as f:
                pdb_data = json.load(f)
                pdb_data = sorted(pdb_data, key=lambda k: k['seq_ident'], reverse=True)
                for i in pdb_data:
                    xyz_in = os.path.join(self.getWorkDirectory(), "mrparse_0", i['pdb_file'])
                    xyz_out = os.path.join(self.getWorkDirectory(), os.path.basename(i['pdb_file']))
                    if os.path.isfile(xyz_in):
                        shutil.copy(xyz_in, xyz_out)
                    self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                    self.container.outputData.XYZOUT[-1].setFullPath(xyz_out)
                    self.container.outputData.XYZOUT[-1].annotation = "AFDB hit: {}".format(i['name'])

        if os.path.exists(esm_json):
            with open(esm_json, 'r') as f:
                pdb_data = json.load(f)
                pdb_data = sorted(pdb_data, key=lambda k: k['seq_ident'], reverse=True)
                for i in pdb_data:
                    xyz_in = os.path.join(self.getWorkDirectory(), "mrparse_0", i['pdb_file'])
                    xyz_out = os.path.join(self.getWorkDirectory(), os.path.basename(i['pdb_file']))
                    if os.path.isfile(xyz_in):
                        shutil.copy(xyz_in, xyz_out)
                    self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                    self.container.outputData.XYZOUT[-1].setFullPath(xyz_out)
                    self.container.outputData.XYZOUT[-1].annotation = "ESM hit: {}".format(i['name'])

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        self.appendCommandLine("--seqin")
        self.appendCommandLine(str(self.seqin))
        if self.hklin:
            self.appendCommandLine("--hklin")
            self.appendCommandLine(str(self.hklin))
        if self.container.options.MAXHITS:
            self.appendCommandLine("--max_hits")
            self.appendCommandLine(str(self.container.options.MAXHITS))
        if self.container.options.DATABASE:
            self.appendCommandLine("--database")
            self.appendCommandLine((str(self.container.options.DATABASE)).lower())
        if self.container.options.USEAPI == 'True':
            self.appendCommandLine("--use_api")
        if self.container.options.PDBLOCAL.isSet():
            self.appendCommandLine("--pdb_local")
            self.appendCommandLine(str(self.container.options.PDBLOCAL))
        if self.container.options.PDBSEQDB.isSet():
            self.appendCommandLine("--pdb_seqdb")
            self.appendCommandLine(str(self.container.options.PDBSEQDB))
        if self.container.options.AFDBSEQDB.isSet():
            self.appendCommandLine("--afdb_seqdb")
            self.appendCommandLine(str(self.container.options.AFDBSEQDB))
        if self.container.options.NPROC:
            self.appendCommandLine("--nproc")
            self.appendCommandLine(str(self.container.options.NPROC))
        if self.container.options.DO_CLASSIFY == 'True':
            self.appendCommandLine("--do_classify")
        self.appendCommandLine("--ccp4cloud")
        return CPluginScript.SUCCEEDED
