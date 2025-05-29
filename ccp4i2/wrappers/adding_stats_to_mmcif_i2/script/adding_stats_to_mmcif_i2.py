import logging
import os
import shutil
import sys
import time
import xml.etree.ElementTree as ET

from onedep import __apiUrl__
from onedep.api.Validate import Validate
from adding_stats_to_mmcif.__main__ import run_process

from ....core import CCP4ErrorHandling
from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)
conversions = {
    'H': 'index_h',
    'K': 'index_k',
    'L': 'index_l',
    'FREER': 'pdbx_r_free_flag',
    'Iplus': 'pdbx_I_plus',
    'SIGIplus': 'pdbx_I_plus_sigma',
                'Iminus': 'pdbx_I_minus',
                'SIGIminus': 'pdbx_I_minus_sigma',
                'I': 'pdbx_intensity_meas',
                'SIGI': 'pdbx_intensity_sigma',
                'Fplus': 'pdbx_F_plus',
                'SIGFplus': 'pdbx_F_plus_sigma',
                'Fminus': 'pdbx_F_minus',
                'SIGFminus': 'pdbx_F_minus_sigma',
                'F': 'F_meas_au',
                'SIGF': 'F_meas_sigma_au',
                'FWT': 'pdbx_FWT',
                'PHWT': 'pdbx_PHWT',
                'DELFWT': 'pdbx_DELFWT',
                'DELPHWT': 'pdbx_DELPHWT',
}


class adding_stats_to_mmcif_i2(CPluginScript):
    # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKNAME = 'adding_stats_to_mmcif_i2'
    TASKVERSION = 0.1               # Version of this plugin
    MAINTAINER = 'martin.noble@ncl.ac.uk'
    ERROR_CODES = {201: {'description': 'Failed to analyse output files'},
                   202: {'description': 'Failed applying selection ot PDB file'},
                   203: {'description': 'Input XML file contains neither AIMLESS nor AIMLESS_PIPE nodes'},
                   204: {'description': 'Failed with adding_stats_to_mmcif'}
                   }
    PURGESEARCHLIST = [['hklin.mtz', 0],
                       ['log_mtzjoin.txt', 0]
                       ]
    RUNEXTERNALPROCESS = False

    def __init__(self, *args, **kws):
        super(adding_stats_to_mmcif_i2, self).__init__(*args, **kws)
        self.xmlroot = ET.Element('adding_stats_to_mmcif')
        #print("db=", PROJECTSMANAGER().db())

    def processInputFiles(self):

        self.fastaFilePath = os.path.join(
            self.getWorkDirectory(), "Sequences.fasta")
        self.container.inputData.ASUCONTENT.writeFasta(self.fastaFilePath)

        self.coordinatesToUse = self.container.inputData.XYZIN
        return CPluginScript.SUCCEEDED

    def startProcess(self, *args, **kwargs):
        self.createReflectionsCif()
        #print("Imported adding_stats_to_mmcif")
        if self.container.controlParameters.USEAIMLESSXML:
            aimless_xml_file = str(
                self.container.inputData.AIMLESSXML.fullPath)
        else:
            aimless_xml_file = ''
        #print("aimless_xml_file is", aimless_xml_file)
        if aimless_xml_file is not '':
            aimlessXML = ET.parse(aimless_xml_file).getroot()
            try:
                lastAimlessNode = aimlessXML.xpath('.//AIMLESS_PIPE')[-1]
            except IndexError:
                try:
                    lastAimlessNode = aimlessXML.xpath('.//AIMLESS')[-1]
                except IndexError:
                    self.appendErrorReport(203, aimless_xml_file)
                    return CPluginScript.FAILED
            self.xmlroot.append(lastAimlessNode)
            # Here create a re-rooted aimless XML file incase the aimless_pipe or aimless modules
            # was nested...this should maybe be nahdled in adding_stats_to_mmcif core code.
            if aimlessXML.tag != 'AIMLESS' and aimlessXML.tag != 'AIMLESS':
                aimless_xml_file = os.path.join(
                    self.workDirectory, 'TempAimlessXML.xml')
                CCP4Utils.writeXml(lastAimlessNode, aimless_xml_file)

        output_mmcif = str(self.container.outputData.MMCIFOUT.fullPath)

        processArgs = {"input_mmcif": str(self.container.inputData.XYZIN.fullPath),
                       "output_mmcif": output_mmcif,
                       "fasta_file": self.fastaFilePath,
                       "xml_file": aimless_xml_file}
        worked = run_process(**processArgs)
        #print(f'Worked is [{worked}]')

        if not worked:
            self.appendErrorReport(self.__class__, 204,
                                   "Failed to adding_stats_to_mmcif")
            return self.reportStatus(CPluginScript.FAILED)

        CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))

        if self.container.controlParameters.SENDTOVALIDATIONSERVER:
            self.performOnedepValidation()

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        with open(str(self.container.outputData.COOTSCRIPTOUT.fullPath), "w") as cootScript:
            cootScript.write('''
try:
    xml_string = open('{}').read()
    valid_inf = parse_wwpdb_validation_xml(xml_string)
    if valid_inf:
        entry_validation_info = valid_inf[0]
        subgroups = valid_inf[1]
        ss = sort_subgroups(subgroups)
        validation_to_gui(entry_validation_info, ss, 0)
    else:
        message = 'problem with valid_inf'
        print(message)
        add_status_bar_text(message)
except IOError as e_mess:
    print('load_validate_xml IOERROR: {}'.format(e_mess))
except Exception as err:
                print("Exception {}".format(err))'''.format(str(self.container.outputData.VALIDATIONXML.fullPath), '{}', '{}'))

        #Set annotations
        self.container.outputData.COOTSCRIPTOUT.annotation= "Output coot script to review issues"
        self.container.outputData.COOTSCRIPTOUT.annotation= "Output coot script to review issues"
        self.container.outputData.COOTSCRIPTOUT.annotation= "Output coot script to review issues"

        return CPluginScript.SUCCEEDED

    def createReflectionsCif(self):
        # Create the (potentially merged) file from which to generate the data that will be spat into the PDB
        if self.container.controlParameters.USEANOMALOUS:
            if self.container.controlParameters.USE_TWIN:
                refmacContentFlag = 1
            else:
                refmacContentFlag = 2
        else:
            if self.container.controlParameters.USE_TWIN:
                refmacContentFlag = 3
            else:
                refmacContentFlag = 4

        pathForExperimentalDataToCifify = str(
            self.container.inputData.F_SIGF.fullPath)
        if self.container.inputData.F_SIGF.contentFlag != refmacContentFlag:
            pathForExperimentalDataToCifify = os.path.join(
                self.getWorkDirectory(), 'mergedForMakingCif.mtz')
            self.hklin, self.columns, error = self.makeHklin0(miniMtzsIn=[['F_SIGF', int(self.container.inputData.F_SIGF.contentFlag)], [
                                                              'F_SIGF', refmacContentFlag], ['FREERFLAG', 0], ['FPHIOUT', 1], ['DIFFPHIOUT', 1]], hklin='mergedForMakingCif', ignoreErrorCodes=[])
            if error.maxSeverity() > CCP4ErrorHandling.Severity.WARNING:
                self.reportStatus(CPluginScript.FAILED)

        # Now convert refmac input to mmcif
        self.hklin2cifPlugin = self.makePluginObject('hklin2cif')
        print(f'scaledunmerged path is {self.container.inputData.SCALEDUNMERGED.fullPath}')
        self.hklin2cifPlugin.container.inputData.SCALEDUNMERGED = self.container.inputData.SCALEDUNMERGED.getFullPath()

        hklinPath = os.path.normpath(pathForExperimentalDataToCifify)
        self.hklin2cifPlugin.container.inputData.HKLIN.setFullPath(hklinPath)
        rv = self.hklin2cifPlugin.process()
        if rv != CPluginScript.SUCCEEDED:
            self.reportStatus(rv)
        srcFile = os.path.join(
            self.hklin2cifPlugin.getWorkDirectory(), 'Reflections.cif')
        shutil.copyfile(srcFile, str(
            self.container.outputData.CIFREFLECTIONS.fullPath))
        return

    def display_status(self, sD, exitOnError=True):
        if 'onedep_error_flag' in sD and sD['onedep_error_flag']:
            print("OneDep error: %s\n" % sD['onedep_status_text'])
            if exitOnError:
                self.reportStatus(CPluginScript.FAILED)
        else:
            if 'status' in sD:
                print("OneDep status: %s\n" % sD['status'])

    def performOnedepValidation(self):
        val = Validate(apiUrl=__apiUrl__)
        ret = val.newSession()
        self.display_status(ret)
        ret = val.inputModelXyzFile(
            str(self.container.outputData.MMCIFOUT.fullPath))
        self.display_status(ret)
        ret = val.inputStructureFactorFile(
            str(self.container.outputData.CIFREFLECTIONS.fullPath))
        self.display_status(ret)
        ret = val.run()
        self.display_status(ret)
        #
        #   Poll for service completion -
        #
        it = 0
        sl = 2
        val_status = None
        while True:
            #    Pause -
            it += 1
            pause = it * it * sl
            time.sleep(pause)
            ret = val.getStatus()
            if 'onedep_error_flag' in ret and ret['onedep_error_flag']:
                ret['status'] = 'failed'
            if ret['status'] in ['completed', 'failed']:
                val_status = ret['status']
                print('validation {}'.format(val_status))
                break
            print("[%4d] Pausing for %4d (seconds)\n" % (it, pause))
            sys.stdout.flush()

        output_pdf_file_name = str(self.container.outputData.PDFOUT.fullPath)
        output_xml_file_name = str(
            self.container.outputData.VALIDATIONXML.fullPath)
        output_svg_file_name = os.path.join(
            self.getWorkDirectory(), 'validation.svg')
        file_name_of_logfile = os.path.join(
            self.getWorkDirectory(), 'validation.log')
        if val_status == 'completed':
            logging.info('getting validation report {}'.format(
                output_pdf_file_name))
            ret = val.getReport(output_pdf_file_name)
            self.container.outputData.PDFOUT.annotation.set(
                'Validation report')
            self.display_status(ret)
            logging.debug('getting report status: {}'.format(ret))
            logging.info('getting validation xml {}'.format(
                output_xml_file_name))
            ret = val.getReportData(output_xml_file_name)
            self.display_status(ret)
            logging.debug('getting xml status: {}'.format(ret))

            validationXML = ET.parse(output_xml_file_name).getroot()
            self.xmlroot.append(validationXML)
            CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))

            if output_svg_file_name:
                ret = val.getOutputByType(
                    output_svg_file_name, contentType="validation-report-slider")
                self.display_status(ret, exitOnError=False)
                logging.debug('getting svg status: {}'.format(ret))
        else:
            logging.error('validation run status: {}'.format(val_status))
            ret = val.getReportLog(file_name_of_logfile)
            logging.debug('getting report log status: {}'.format(ret))
            logging.error('log file: "{}"'.format(file_name_of_logfile))
