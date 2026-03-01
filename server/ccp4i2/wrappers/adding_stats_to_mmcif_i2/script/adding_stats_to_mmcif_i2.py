import logging
import os
import shutil
import sys
import time

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript

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
    TASKNAME = 'adding_stats_to_mmcif_i2'
    ERROR_CODES = {203: {'description': 'Input XML file contains neither AIMLESS nor AIMLESS_PIPE nodes'},
                   204: {'description': 'Failed with adding_stats_to_mmcif'}}

    def __init__(self, *args, **kws):
        super(adding_stats_to_mmcif_i2, self).__init__(*args, **kws)
        self.xmlroot = etree.Element('adding_stats_to_mmcif')

    def trace_xyzin_lineage(self, file_uuid):
        """Trace XYZIN file back through refinement -> scaling job chain.

        Given the UUID of a coordinates file, finds the refinement job that
        produced it, then collects all related input/output files and control
        parameters needed for deposition.

        Returns a JSON-serializable dict with keys:
            success, files, paths, params, warnings,
            refinement_job_id, refinement_task_name
        """
        from pathlib import Path

        from lxml import etree as ET

        from ccp4i2.db.models import File, FileUse, Job

        result = {
            "success": False,
            "files": {},
            "paths": {},
            "params": {},
            "warnings": [],
        }

        # 1. Find which job produced this file
        try:
            the_file = File.objects.get(uuid=file_uuid)
        except File.DoesNotExist:
            result["warnings"].append("File not found in database")
            return result

        ref_job = the_file.job
        if ref_job is None:
            result["warnings"].append("File has no producing job")
            return result

        # 2. Validate it's a refinement job
        refinement_tasks = ["prosmart_refmac", "servalcat_pipe", "buster"]
        if ref_job.task_name not in refinement_tasks:
            result["warnings"].append(
                f"Coordinates should come from a refinement job, "
                f"not '{ref_job.task_name}'"
            )
            return result

        if ref_job.status != Job.Status.FINISHED:
            result["warnings"].append(
                "Source refinement job did not finish successfully"
            )
            return result

        result["success"] = True
        result["refinement_job_id"] = ref_job.id
        result["refinement_task_name"] = ref_job.task_name

        # 3. Get refinement job's INPUT files via FileUse
        input_uses = FileUse.objects.filter(
            job=ref_job, role=FileUse.Role.IN
        ).select_related("file")
        for use in input_uses:
            file_uuid_str = str(use.file.uuid)
            param = use.job_param_name
            # Map input reflections to canonical "F_SIGF" regardless of
            # refinement task naming conventions:
            #   prosmart_refmac / buster → "F_SIGF"
            #   servalcat_pipe → "HKLIN"
            if param in ("F_SIGF", "HKLIN"):
                result["files"]["F_SIGF"] = file_uuid_str
            elif param in ("FREERFLAG", "TLSIN"):
                result["files"][param] = file_uuid_str
            elif param in ("DICT", "DICT_LIST"):
                result["files"].setdefault("DICT_LIST", [])
                result["files"]["DICT_LIST"].append(file_uuid_str)

        # 4. Get refinement job's OUTPUT files (may override inputs)
        output_files = File.objects.filter(job=ref_job)
        for f in output_files:
            param = f.job_param_name
            if param == "TLSOUT":
                result["files"]["TLSIN"] = str(f.uuid)
            elif param in ("FPHIOUT", "DIFFPHIOUT"):
                result["files"][param] = str(f.uuid)

        # 5. REFMACINPUTPARAMSXML = input_params.xml from refinement job dir
        input_params_path = Path(ref_job.directory) / "input_params.xml"
        if input_params_path.exists():
            result["paths"]["REFMACINPUTPARAMSXML"] = str(input_params_path)

        # 6. Get USEANOMALOUS and USE_TWIN from refinement job's params
        #    (buster doesn't have these — skip gracefully)
        if ref_job.task_name != "buster":
            try:
                from ccp4i2.lib.utils.plugins.get_plugin import get_job_plugin

                ref_plugin = get_job_plugin(ref_job)
                ctrl = ref_plugin.container.controlParameters
                if hasattr(ctrl, "USEANOMALOUS"):
                    result["params"]["USEANOMALOUS"] = bool(ctrl.USEANOMALOUS)
                if hasattr(ctrl, "USE_TWIN"):
                    result["params"]["USE_TWIN"] = bool(ctrl.USE_TWIN)
            except Exception as e:
                result["warnings"].append(
                    f"Could not read refinement control parameters: {e}"
                )

        # 7. Trace F_SIGF back to its producing (scaling) job
        f_sigf_uuid = result["files"].get("F_SIGF")
        if f_sigf_uuid:
            try:
                f_sigf_file = File.objects.select_related("job").get(
                    uuid=f_sigf_uuid
                )
                scaling_job = f_sigf_file.job
                if scaling_job:
                    result["scaling_job_id"] = scaling_job.id

                    # Look for XMLOUT from scaling job (aimless statistics)
                    xmlout_files = File.objects.filter(
                        job=scaling_job, job_param_name="XMLOUT"
                    )
                    if xmlout_files.exists():
                        xmlout = xmlout_files.first()
                        xmlout_path = Path(scaling_job.directory) / xmlout.name
                        if xmlout_path.exists():
                            tree = ET.parse(str(xmlout_path))
                            root = tree.getroot()
                            has_aimless = (
                                root.findall(".//AIMLESS")
                                or root.findall(".//AIMLESS_PIPE")
                            )
                            if has_aimless:
                                result["files"]["AIMLESSXML"] = str(xmlout.uuid)
                                result["params"]["USEAIMLESSXML"] = True
                            else:
                                result["warnings"].append(
                                    "Scaling job XML does not contain "
                                    "AIMLESS statistics"
                                )
                    else:
                        # Fallback: try program.xml on disk
                        program_xml_path = (
                            Path(scaling_job.directory) / "program.xml"
                        )
                        if program_xml_path.exists():
                            tree = ET.parse(str(program_xml_path))
                            root = tree.getroot()
                            has_aimless = (
                                root.findall(".//AIMLESS")
                                or root.findall(".//AIMLESS_PIPE")
                            )
                            if has_aimless and root.tag != "IMPORT_MERGED":
                                result["paths"]["AIMLESSXML"] = str(
                                    program_xml_path
                                )
                                result["params"]["USEAIMLESSXML"] = True

                    # Look for UNMERGEDOUT — stored as UNMERGEDOUT[0],
                    # UNMERGEDOUT[1], etc. for multi-wavelength datasets
                    scaling_outputs = File.objects.filter(
                        job=scaling_job,
                        job_param_name__startswith="UNMERGEDOUT",
                    )
                    count = scaling_outputs.count()
                    if count >= 1:
                        # Take the first (usually only) unmerged output
                        result["files"]["SCALEDUNMERGED"] = str(
                            scaling_outputs[0].uuid
                        )
                        result["params"]["INCLUDEUNMERGED"] = True
                        if count > 1:
                            result["warnings"].append(
                                f"Multiple unmerged outputs found — "
                                f"using first of {count}"
                            )
                    else:
                        result["warnings"].append(
                            "Scaling job did not output unmerged data"
                        )
            except File.DoesNotExist:
                result["warnings"].append(
                    "Could not trace F_SIGF to its producing job"
                )

        return result

    def processInputFiles(self):

        self.fastaFilePath = os.path.join(
            self.getWorkDirectory(), "Sequences.fasta")
        self.container.inputData.ASUCONTENT.writeFasta(self.fastaFilePath)

        self.coordinatesToUse = self.container.inputData.XYZIN
        return CPluginScript.SUCCEEDED

    def startProcess(self):
        rv = self.createReflectionsCif()
        if rv is not None:
            return rv

        from adding_stats_to_mmcif.__main__ import run_process

        #print("Imported adding_stats_to_mmcif")
        if self.container.controlParameters.USEAIMLESSXML:
            aimless_xml_file = str(
                self.container.inputData.AIMLESSXML.fullPath)
        else:
            aimless_xml_file = ''
        #print("aimless_xml_file is", aimless_xml_file)
        if aimless_xml_file != '':
            with open(aimless_xml_file, "r") as aimlessXMLFile:
                aimlessXML = etree.fromstring(aimlessXMLFile.read())
                try:
                    lastAimlessNode = aimlessXML.xpath('.//AIMLESS_PIPE')[-1]
                except IndexError as err:
                    try:
                        lastAimlessNode = aimlessXML.xpath('.//AIMLESS')[-1]
                    except IndexError as err:
                        self.appendErrorReport(203, aimless_xml_file)
                        return CPluginScript.FAILED
                self.xmlroot.append(lastAimlessNode)
            # Here create a re-rooted aimless XML file incase the aimless_pipe or aimless modules
            # was nested...this should maybe be nahdled in adding_stats_to_mmcif core code.
            if aimlessXML.tag != 'AIMLESS' and aimlessXML.tag != 'AIMLESS':
                aimless_xml_file = os.path.join(
                    self.workDirectory, 'TempAimlessXML.xml')
                with open(aimless_xml_file, 'wb') as tempAimlessXML:
                    tempAimlessXML.write(etree.tostring(lastAimlessNode))

        output_mmcif = str(self.container.outputData.MMCIFOUT.fullPath)

        processArgs = {"input_mmcif": str(self.container.inputData.XYZIN.fullPath),
                       "output_mmcif": output_mmcif,
                       "fasta_file": self.fastaFilePath,
                       "xml_file": aimless_xml_file}
        try:
            worked = run_process(**processArgs)
        except Exception as e:
            print(f'adding_stats_to_mmcif raised exception: {e}')
            import traceback
            traceback.print_exc()
            self.appendErrorReport(204, str(e))
            return CPluginScript.FAILED

        if not worked:
            self.appendErrorReport(204, "Failed in adding_stats_to_mmcif")
            return CPluginScript.FAILED

        with open(self.makeFileName('PROGRAMXML'), 'w') as programXML:
            CCP4Utils.writeXML(programXML, etree.tostring(self.xmlroot))

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
        import subprocess

        # Merge F_SIGF + FREERFLAG + FPHIOUT + DIFFPHIOUT if content flags differ
        if self.container.controlParameters.USEANOMALOUS:
            refmacContentFlag = 1 if self.container.controlParameters.USE_TWIN else 2
        else:
            refmacContentFlag = 3 if self.container.controlParameters.USE_TWIN else 4

        mergedMtzPath = str(self.container.inputData.F_SIGF.fullPath)
        if self.container.inputData.F_SIGF.contentFlag != refmacContentFlag:
            mergedMtzPath = os.path.join(
                self.getWorkDirectory(), 'mergedForMakingCif.mtz')
            self.hklin, self.columns, error = self.makeHklin0(
                miniMtzsIn=[
                    ['F_SIGF', int(self.container.inputData.F_SIGF.contentFlag)],
                    ['F_SIGF', refmacContentFlag],
                    ['FREERFLAG', 1],
                    ['FPHIOUT', 1],
                    ['DIFFPHIOUT', 1],
                ],
                hklin='mergedForMakingCif',
            )
            if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED

        # Rename CCP4i2's prefixed column labels to standard names gemmi expects
        import gemmi as _gemmi
        mtz = _gemmi.read_mtz_file(os.path.normpath(mergedMtzPath))
        COLUMN_RENAME = {
            'FREERFLAG_FREER': 'FREE',
            'F_SIGF_Iplus': 'I(+)',   'F_SIGF_SIGIplus': 'SIGI(+)',
            'F_SIGF_Iminus': 'I(-)',   'F_SIGF_SIGIminus': 'SIGI(-)',
            'F_SIGF_Fplus': 'F(+)',    'F_SIGF_SIGFplus': 'SIGF(+)',
            'F_SIGF_Fminus': 'F(-)',   'F_SIGF_SIGFminus': 'SIGF(-)',
            'F_SIGF_F': 'FP',         'F_SIGF_SIGF': 'SIGFP',
            'F_SIGF_Fmean': 'FP',     'F_SIGF_SIGFmean': 'SIGFP',
            'F_SIGF_I': 'IMEAN',      'F_SIGF_SIGI': 'SIGIMEAN',
            'F_SIGF_Imean': 'IMEAN',  'F_SIGF_SIGImean': 'SIGIMEAN',
            'FPHIOUT_F': 'FWT',       'FPHIOUT_PHI': 'PHWT',
            'DIFFPHIOUT_F': 'DELFWT', 'DIFFPHIOUT_PHI': 'PHDELWT',
        }
        renamed = []
        for col in mtz.columns:
            if col.label in COLUMN_RENAME:
                old = col.label
                col.label = COLUMN_RENAME[old]
                renamed.append(f'{old} -> {col.label}')
        if renamed:
            print(f'Renamed columns: {", ".join(renamed)}')
        renamedMtzPath = os.path.join(self.getWorkDirectory(), 'forDeposition.mtz')
        mtz.write_to_file(renamedMtzPath)

        # Call gemmi mtz2cif directly (no hklin2cif sub-plugin needed)
        outputCifPath = str(self.container.outputData.CIFREFLECTIONS.fullPath)
        cmd = ['gemmi', 'mtz2cif']

        includeUnmerged = (
            self.container.controlParameters.INCLUDEUNMERGED
            and self.container.inputData.SCALEDUNMERGED.isSet()
        )
        if includeUnmerged:
            cmd.append('--depo')

        cmd.append(renamedMtzPath)

        if includeUnmerged:
            cmd.append(str(self.container.inputData.SCALEDUNMERGED.fullPath))

        cmd.append(outputCifPath)

        print(f'Running: {" ".join(cmd)}')
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            cwd=self.getWorkDirectory(),
        )
        if result.stdout:
            print(result.stdout)
        if result.returncode != 0:
            print(f'gemmi mtz2cif failed (exit code {result.returncode}):')
            print(result.stderr)
            return CPluginScript.FAILED
        if not os.path.exists(outputCifPath):
            print('gemmi mtz2cif did not produce output CIF')
            return CPluginScript.FAILED

    def display_status(self, sD, exitOnError=True):
        if 'onedep_error_flag' in sD and sD['onedep_error_flag']:
            print("OneDep error: %s\n" % sD['onedep_status_text'])
            if exitOnError:
                self.reportStatus(CPluginScript.FAILED)
        else:
            if 'status' in sD:
                print("OneDep status: %s\n" % sD['status'])

    def performOnedepValidation(self):
        from onedep import __apiUrl__
        from onedep.api.Validate import Validate
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

            with open(output_xml_file_name, "r") as validationXMLFile:
                validationXML = etree.fromstring(validationXMLFile.read())
                self.xmlroot.append(validationXML)
                with open(self.makeFileName('PROGRAMXML'), 'w') as thisXMLFile:
                    CCP4Utils.writeXML(thisXMLFile, etree.tostring(
                        self.xmlroot, pretty_print=True))

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
