import unittest
from .access_test_files import TestFiles
import tempfile
import logging
import os
import shutil

import wrappers.adding_stats_to_mmcif_i2.script.wwpdb_validation.wwpdb_validation.wwpdb_validation_api as wwpdb_validation_api


logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.DEBUG)

class TestValidationReportGeneration(unittest.TestCase):

    def setUp(self):
        self.test_files = TestFiles()

    def test_3zt9(self):
        test_dir = tempfile.mkdtemp()
        logging.debug('test dir: {}'.format(test_dir))
        output_pdf_file = os.path.join(test_dir, 'output.pdf')
        output_xml_file = os.path.join(test_dir, 'output.xml')
        output_log_file = os.path.join(test_dir, 'output.log')
        output_cif_file = os.path.join(test_dir, 'output.cif')
        output_svg_file = os.path.join(test_dir, 'output.svg')
        self.test_files.pdb3zt9()
        worked = wwpdb_validation_api.run_validation_api(
            model_file_path = self.test_files.cif,
            structure_factors= self.test_files.sf,
            output_pdf_file_name=output_pdf_file,
            output_xml_file_name=output_xml_file,
            output_log_file_name=output_log_file,
            )
        self.assertTrue(worked)
        shutil.rmtree(test_dir)


if __name__ == '__main__':
    unittest.main()
