import os
import shutil
import xml.etree.ElementTree as ET

from ample.constants import AMPLE_PKL
from ample.util import mrbump_util
from ample.util.ample_util import I2DIR

from ....core import CCP4ErrorHandling
from ....core.CCP4PluginScript import CPluginScript


class AMPLE(CPluginScript):
    TASKNAME = 'AMPLE'  # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION = 0.1  # Version of this plugin
    MAINTAINER = 'jens.thomas@liv.ac.uk'
    ERROR_CODES = {
        1: {
            'description': 'Something not very good has happened.'
        },
    }
    WHATNEXT = ['prosmart_refmac', 'modelcraft', 'coot_rebuild']
    #     PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
    #                        ['log_mtzjoin.txt', 0]
    #                        ]
    TASKCOMMAND = "ample"

    # Andre's stuff for a clean shutdown - a file caled INTERUPT will be created.
    #     INTERRUPTABLE = True
    #     INTERRUPTLABEL = 'Stop and keep current best solution'

    def __init__(self, *args, **kws):
        super(AMPLE, self).__init__(*args, **kws)

    def processInputFiles(self):
        from ....core import CCP4XtalData
        # No idea why we need the 'AMPLE_F_SIGF' bit...
        self.hklin, self.columns, error = self.makeHklin0(
            [['AMPLE_F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
        if error.maxSeverity() > CCP4ErrorHandling.Severity.WARNING:
            return CPluginScript.FAILED
        if self.hklin is None: return CPluginScript.FAILED

        self.F, self.SIGF = self.columns.split(',')
        self.fasta = self.container.inputData.AMPLE_SEQIN

        return self.SUCCEEDED

    def makeCommandAndScript(self):
        params = self.container.inputData

        run_type = None
        # Switch for run_type
        ABINITIO = 0
        IMPORT_MODELS = 1
        IMPORT_HOMOLOGS = 2
        NMR_REMODEL = 3
        NMR_IMPORT = 4
        IDEAL_HELICES = 5
        ROSETTA_TM = 6
        ROSETTA = 7
        ROSETTA_CC = 8
        HELICAL_ENSEMBLES = 9

        # Calculate the run_type
        run_type = None
        if params.AMPLE_EXISTING_MODELS == 'True':
            if params.AMPLE_MODEL_TYPE == 'abinitio':
                run_type = IMPORT_MODELS
            elif params.AMPLE_MODEL_TYPE == 'multiple_homologs':
                run_type = IMPORT_HOMOLOGS
            elif params.AMPLE_MODEL_TYPE == 'nmr_ensemble':
                if params.AMPLE_NMR_REMODEL == 'nmr_remodel_true':
                    run_type = NMR_REMODEL
                elif params.AMPLE_NMR_REMODEL == 'nmr_remodel_false':
                    run_type = NMR_IMPORT
                else:
                    assert False, "Unrecognised Parameter: {0}".format(
                        params.AMPLE_NMR_REMODEL)
            else:
                assert False, "Unrecognised Parameter: {0}".format(
                    params.AMPLE_MODEL_TYPE)
        else:
            # No models
            if params.AMPLE_MODEL_GENERATION == 'ideal_helices':
                run_type = IDEAL_HELICES
            elif params.AMPLE_MODEL_GENERATION == 'helical_ensembles':
                run_type = HELICAL_ENSEMBLES
            elif params.AMPLE_MODEL_GENERATION == 'rosetta':
                if params.AMPLE_PROTEIN_CLASS == 'transmembrane':
                    run_type = ROSETTA_TM
                elif params.AMPLE_PROTEIN_CLASS == "coiled-coil":
                    run_type = ROSETTA_CC
                elif params.AMPLE_PROTEIN_CLASS == 'globular':
                    run_type = ROSETTA
                else:
                    assert False, "Unrecognised Parameter: {0}".format(
                        params.AMPLE_PROTEIN_CLASS)
            else:
                assert False, "Unrecognised Parameter: {0}".format(
                    params.AMPLE_MODEL_GENERATION)

        # Sort out the model file
        if params.AMPLE_MODELS_SOURCE == 'directory':
            models_file = params.AMPLE_MODELS_DIR
        elif params.AMPLE_MODELS_SOURCE == 'file':
            models_file = params.AMPLE_MODELS_FILE
        else:
            assert False, "Unrecognised Parameter: {0}".format(
                params.AMPLE_MODELS_FILE)

        # Add modelling parameters shared by all run_types
        #self.appendCommandLine(self.getWorkDirectory())
        self.appendCommandLine(['-fasta', self.fasta])
        #         self.appendCommandLine( params.AMPLE_SEQIN)
        #         self.appendCommandLine('-mtz')
        #         self.appendCommandLine(params.AMPLE_F_SIGF)
        #         self.appendCommandLine('-F')
        #         self.appendCommandLine( self.columnsAsArray[0])
        #         self.appendCommandLine('-SIGF')
        #         self.appendCommandLine( self.columnsAsArray[1])
        self.appendCommandLine(['-mtz', self.hklin])
        self.appendCommandLine(['-F', self.F])
        self.appendCommandLine(['-SIGF', self.SIGF])

        # Model source if using existing models
        if run_type in [IMPORT_MODELS, IMPORT_HOMOLOGS]:
            self.appendCommandLine(['-models', models_file])
        elif run_type in [NMR_REMODEL, NMR_IMPORT]:
            self.appendCommandLine(['-nmr_model_in', models_file])

        # Generating models with rosetta
        if run_type in [ROSETTA, ROSETTA_CC, ROSETTA_TM, NMR_REMODEL]:
            self.appendCommandLine(['-rosetta_dir', params.AMPLE_ROSETTA_DIR])
            self.appendCommandLine(
                ['-frags_3mers', params.AMPLE_ROSETTA_FRAGS3])
            self.appendCommandLine(
                ['-frags_9mers', params.AMPLE_ROSETTA_FRAGS9])
            if params.AMPLE_CONTACT_FILE and params.AMPLE_CONTACT_FILE.exists(
            ):
                self.appendCommandLine(
                    ['-contact_format', params.AMPLE_CONTACT_FORMAT])
                self.appendCommandLine(
                    ['-contact_file', params.AMPLE_CONTACT_FILE])

        # Runtype-specific flags
        if run_type == IDEAL_HELICES:
            self.appendCommandLine(['-ideal_helices', 'True'])
        elif run_type == HELICAL_ENSEMBLES:
            self.appendCommandLine(['-helical_ensembles', 'True'])
        elif run_type == ROSETTA:
            pass  # Nothing to do currently
        elif run_type == ROSETTA_TM:
            self.appendCommandLine(['-transmembrane', 'True'])
        elif run_type == ROSETTA_CC:
            self.appendCommandLine(['-coiled_coil', 'True'])
        elif run_type == IMPORT_HOMOLOGS:
            self.appendCommandLine(['-homologs', 'True'])
        elif run_type == NMR_REMODEL:
            self.appendCommandLine(['-nmr_remodel', 'True'])

        # Stuff that applies to all runtypes
        if params.AMPLE_ENSEMBLING_TM:
            self.appendCommandLine(['-cluster_method', 'spicker_tm'])
        self.appendCommandLine(['-use_shelxe', str(params.AMPLE_USE_SHELXE)])
        self.appendCommandLine(
            ['-shelxe_rebuild',
             str(params.AMPLE_SHELXE_REBUILD)])
        if params.AMPLE_REFINE_REBUILD:
            self.appendCommandLine(['-refine_rebuild_arpwarp', 'True'])

        if len(params.AMPLE_EXTRA_FLAGS):
            # Split the text in the box on whitespace and newline characters
            self.appendCommandLine([
                b for a in params.AMPLE_EXTRA_FLAGS.split(os.linesep)
                for b in a.split() if b
            ])
        # General flags
        self.appendCommandLine(['-nproc', str(params.AMPLE_NPROC)])
        self.appendCommandLine(
            ['-ccp4i2_xml', self.makeFileName('PROGRAMXML')])
        return self.SUCCEEDED

    def flushXML(self):
        tmpFilename = self.makeFileName('PROGRAMXML') + '_tmp'
        with open(tmpFilename, 'wb') as xmlFile:
            ET.indent(self.xmlroot)
            xmlFile.write(ET.tostring(self.xmlroot))
        if os.path.exists(self.makeFileName('PROGRAMXML')):
            os.remove(self.makeFileName('PROGRAMXML'))
        os.rename(tmpFilename, self.makeFileName('PROGRAMXML'))

    def processOutputFiles(self):
        '''
        Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        '''
        #debug_console()
        # programXML file is generated by pyrvapi so we only handle the results specific for I2 here.
        pkl = os.path.join(self.getWorkDirectory(), I2DIR, AMPLE_PKL)
        top_files = mrbump_util.ResultsSummary(results_pkl=pkl).topFiles()
        if top_files:
            for indx, file_info in enumerate(top_files):
                # Need to copy the files into the actual project directory - cannot be a sub-directory. Not entirely sure why but...
                xyz = os.path.join(self.getWorkDirectory(),
                                   os.path.basename(file_info['pdb']))
                mtz = os.path.join(self.getWorkDirectory(),
                                   os.path.basename(file_info['mtz']))
                if os.path.isfile(file_info['pdb']):
                    shutil.copy2(file_info['pdb'], xyz)
                if os.path.isfile(file_info['mtz']):
                    shutil.copy2(file_info['mtz'], mtz)
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(xyz)
                self.container.outputData.XYZOUT[-1].annotation = 'Coordinates for solution {} ({}): {}'.format(
                    indx + 1, file_info['name'], file_info['info'])
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(mtz)
                self.container.outputData.HKLOUT[-1].annotation = 'MTZ for solution {} ({}): {}'.format(
                    indx + 1, file_info['name'], file_info['info'])
            self.splitHkloutList(
                miniMtzsOut=['FPHIOUT'],
                programColumnNames=['FWT,PHWT'],
                outputBaseName=['FPHIOUT'],
                infileList=self.container.outputData.HKLOUT)
            for indx in range(len(self.container.outputData.FPHIOUT)):
                file_info = top_files[indx]
                self.container.outputData.FPHIOUT[indx].annotation = 'Map for solution {} ({}): {}'.format(
                    indx + 1, file_info['name'], file_info['info'])
                self.container.outputData.FPHIOUT[indx].contentFlag = 1
                self.container.outputData.FPHIOUT[indx].subType = 1

        return self.SUCCEEDED
