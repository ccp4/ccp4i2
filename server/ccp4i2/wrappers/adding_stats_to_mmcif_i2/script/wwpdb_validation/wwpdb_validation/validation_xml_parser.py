#!/usr/bin/env python
import logging
import pprint
import os
import argparse
from .xml_parsing import parse_xml
from ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.wwpdb_validation.tests.access_test_files import TestFiles

logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)

class validationReport:
    def __init__(self, xml_file):
        """
        :param xml_file: input aimless XML file
        """
        self.xml_file = xml_file
        self.tree = None
        self.root = None
        self.result = dict()
        self.entry = dict()
        self.programs = dict()

    def parse_xml(self):
        """
            checks input file is XML file, parses it.
            prevents parsing twice by checking if self.tree already exists
        :return: True if a parsed aimless XML file, False if not
        """
        if not self.tree:
            self.root = parse_xml(xml_file=self.xml_file)

        if self.root is not None:
            if self.root.tag == 'wwPDB-validation-information':
                logging.debug('is an wwPDB validation file')
                return True
        return False

    def get_entry_info(self):
        nodes = self.root.findall(".//Entry")
        for node in nodes:
            for attrName, attrValue in node.attrib.items():
                self.entry[attrName] = attrValue

    def get_programs(self):
        nodes = self.root.findall(".//program")
        for node in nodes:
            prog_dict = dict()
            for attrName, attrValue in node.attrib.items():
                prog_dict[attrName] = attrValue
            prog_name = prog_dict.get('name', '')
            if prog_name:
                self.programs[prog_name] = prog_dict

    def run_process(self):
        self.parse_xml()
        self.get_entry_info()
        self.get_programs()

if __name__ == '__main__':
    test_data = TestFiles()
    test_data.pdb1cbs()
    logging.info(test_data.xml)

    vr = validationReport(xml_file=test_data.xml)
    vr.run_process()
    print(vr.entry)
    print(vr.programs)

