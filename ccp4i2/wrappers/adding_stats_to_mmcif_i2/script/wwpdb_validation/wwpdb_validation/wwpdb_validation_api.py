from onedep import __apiUrl__
from onedep.api.Validate import Validate
import time
import os
import argparse
import logging
import shutil
import tempfile

logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)


def display_status(sD, exitOnError=True):
    if 'onedep_error_flag' in sD and sD['onedep_error_flag']:
        logging.error("OneDep error: %s\n" % sD['onedep_status_text'])
        if exitOnError:
            raise SystemExit()
    else:
        if 'status' in sD:
            logging.info("OneDep status: %s\n" % sD['status'])


def run_validation_api(model_file_path, output_pdf_file_name,
                       output_xml_file_name,
                       structure_factors=None,
                       nmr_chemical_shifts=None,
                       em_volume=None,
                       output_cif_file_name=None,
                       output_svg_file_name=None,
                       output_log_file_name=None,
                       api_input_url=None):
    # Given:
    # modelFilePath contains the path to the model file
    # sfFilePath contains the path to the structure factor file
    if api_input_url:
        api_url = api_input_url
    else:
        api_url = __apiUrl__
    val = Validate(apiUrl=api_url)
    ret = val.newSession()
    display_status(ret)
    ret = val.inputModelXyzFile(model_file_path)
    display_status(ret)
    if structure_factors:
        ret = val.inputStructureFactorFile(structure_factors)
        display_status(ret)
    if nmr_chemical_shifts:
        ret = val.inputNmrChemicalShiftsFile(nmr_chemical_shifts)
        display_status(ret)
    if em_volume:
        ret = val.inputEmVolume(em_volume)
        display_status(ret)
    ret = val.run()
    display_status(ret)
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
        if ret['status'] in ['completed', 'failed']:
            val_status = ret['status']
            logging.info('validation {}'.format(val_status))
            break
        logging.info("[%4d] Pausing for %4d (seconds)\n" % (it, pause))
        #
        #
    lt = time.strftime("%Y%m%d%H%M%S", time.localtime())
    temp_output_dir = None
    if output_log_file_name:
        file_name_of_logfile = output_log_file_name
    else:
        temp_output_dir = tempfile.mkdtemp()
        file_name_of_logfile = os.path.join(temp_output_dir, "xray-report-{}.log".format(lt))
    if val_status == 'completed':
        logging.info('getting validation report {}'.format(output_pdf_file_name))
        ret = val.getReport(output_pdf_file_name)
        display_status(ret)
        logging.debug('getting report status: {}'.format(ret))
        logging.info('getting validation xml {}'.format(output_xml_file_name))
        ret = val.getReportData(output_xml_file_name)
        display_status(ret)
        logging.debug('getting xml status: {}'.format(ret))
        if output_cif_file_name:
            ret = val.getOutputByType(output_cif_file_name, contentType="model")
            display_status(ret)
            logging.debug('getting cif status: {}'.format(ret))
        if output_svg_file_name:
            ret = val.getOutputByType(output_svg_file_name, contentType="validation-report-slider")
            display_status(ret)
            logging.debug('getting svg status: {}'.format(ret))
        ret = val.getReportLog(file_name_of_logfile)
        logging.debug('getting report log status: {}'.format(ret))
        logging.info('log file: "{}"'.format(file_name_of_logfile))
        if os.path.exists(output_pdf_file_name and output_xml_file_name):
            if temp_output_dir:
                shutil.rmtree(temp_output_dir)
            return True
    else:
        logging.error('validation run status: {}'.format(val_status))
        ret = val.getReportLog(file_name_of_logfile)
        logging.debug('getting report log status: {}'.format(ret))
        logging.error('log file: "{}"'.format(file_name_of_logfile))

    return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_mmcif', help='input mmcif model file', type=str, required=True)
    parser.add_argument('--input_sf_file', help='input sf file', type=str)
    parser.add_argument('--input_cs_file', help='input chemical shift file', type=str)
    parser.add_argument('--input_em_volume_file', help='input EM volume file', type=str)
    parser.add_argument('--output_pdf_file_name', help='output pdf file name', type=str, required=True)
    parser.add_argument('--output_xml_file_name', help='output xml file name', type=str, required=True)
    parser.add_argument('--output_cif_file_name', help='output cif file name', type=str)
    parser.add_argument('--output_svg_file_name', help='output svg file name', type=str)
    parser.add_argument('--output_log_file_name', help='output log file name', type=str)
    parser.add_argument('--api_url', help='input api url', type=str)
    parser.add_argument('-d', '--debug', help='debugging', action='store_const', dest='loglevel', const=logging.DEBUG,
                        default=logging.INFO)

    args = parser.parse_args()

    logger.setLevel(args.loglevel)

    worked = run_validation_api(structure_factors=args.input_sf_file, model_file_path=args.input_mmcif,
                                nmr_chemical_shifts=args.input_cs_file,
                                em_volume=args.input_em_volume_file,
                                output_pdf_file_name=args.output_pdf_file_name,
                                output_xml_file_name=args.output_xml_file_name,
                                output_cif_file_name=args.output_cif_file_name,
                                output_svg_file_name=args.output_svg_file_name,
                                output_log_file_name=args.output_log_file_name,
                                api_input_url=args.api_url)
    logging.info('worked: {}'.format(worked))
