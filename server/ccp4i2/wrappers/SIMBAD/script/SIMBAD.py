import os
import shutil

from simbad.util import SIMBAD_DIRNAME
from simbad.util.simbad_results import SimbadResults

from ccp4i2.core.CCP4PluginScript import CPluginScript


class SIMBAD(CPluginScript):
    TASKNAME = 'SIMBAD'
    MAINTAINER = 'hlasimpk@liv.ac.uk'
    ERROR_CODES = {1: {'description' : 'Something not very good has happened.'}}
    WHATNEXT = ['prosmart_refmac', 'modelcraft', 'coot_rebuild']
    TASKCOMMAND="simbad"

    def processInputFiles(self):
        #Preprocess reflections to generate an "HKLIN" file
        """
        #makeHklin0 takes as arguments a list of sublists
        #Each sublist comprises 1) A reflection data object identifier (one of those specified in the inputData container
        #                           the task in the corresponding .def.xml
        #                       2) The requested data representation type to be placed into the file that is generated
        #
        #makeHklin0 returns a tuple comprising:
        #                       1) the file path of the file that has been created
        #                       2) a list of strings, each of which contains a comma-separated list of column labels output from
        #                       the input data objects
        #                       3) A CCP4 Error object
        """
        from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
        self.hklin, error = self.makeHklin([['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])

        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        if self.hklin is None:
            return CPluginScript.FAILED
        return self.SUCCEEDED

    def makeCommandAndScript(self):
        params = self.container.inputData
        if params.SIMBAD_SEARCH_LEVEL == 'Lattice':
            self.TASKCOMMAND = 'simbad-lattice'
        elif params.SIMBAD_SEARCH_LEVEL == 'Contaminants':
            self.TASKCOMMAND = 'simbad-contaminant'
        elif params.SIMBAD_SEARCH_LEVEL == 'Lattice + contaminants':
            self.TASKCOMMAND = 'simbad'
        else: assert False
            
        # General flags
        self.appendCommandLine(['-ccp4i2_xml', self.makeFileName('PROGRAMXML')])
        self.appendCommandLine(['-nproc', str(params.SIMBAD_NPROC)])
        
        # Program-specific
        if params.SIMBAD_SEARCH_LEVEL != 'Lattice' and params.SIMBAD_ORGANISM != 'ALL':
            self.appendCommandLine(['-organism', params.SIMBAD_ORGANISM])
        
        # Add the advanced options
        if params.SIMBAD_SEARCH_LEVEL != 'Lattice':
            self.appendCommandLine(['-rot_program', params.SIMBAD_ROT_PROGRAM])
        self.appendCommandLine(['-mr_program', params.SIMBAD_MR_PROGRAM])
        self.appendCommandLine(['-nmol', params.SIMBAD_NMOL])
        if params.SIMBAD_PROCESS_ALL:
            self.appendCommandLine(['--process_all'])
        self.appendCommandLine(['-sga', params.SIMBAD_SGALTERNATIVE])
            
        # Finally add the mtz file
        self.appendCommandLine([self.hklin])
        return self.SUCCEEDED

    def processOutputFiles(self):
        """
        Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        """
        top_files = SimbadResults(os.path.join(self.getWorkDirectory(),SIMBAD_DIRNAME)).top_files()
        if top_files:
            for fo in top_files:
                # Need to copy the files into the actual project directory - cannot be a sub-directory. Not entirely sure why but...
                xyz = os.path.join(self.getWorkDirectory(),os.path.basename(fo.ref_pdb))
                mtz = os.path.join(self.getWorkDirectory(),os.path.basename(fo.ref_mtz))
                # REMOVE CHECK AS FILES SHOULD EXIST
                if os.path.isfile(fo.ref_pdb): shutil.copy2(fo.ref_pdb, xyz)
                if os.path.isfile(fo.ref_mtz): shutil.copy2(fo.ref_mtz, mtz)
                self.container.outputData.XYZOUT.append(xyz)
                self.container.outputData.XYZOUT[-1].annotation = fo.ref_pdb_annotation
                self.container.outputData.HKLOUT.append(mtz)
                self.container.outputData.HKLOUT[-1].annotation = fo.ref_mtz_annotation

        return self.SUCCEEDED
