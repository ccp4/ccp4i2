#!/usr/bin/env python

import logging
import os
import xml.etree.ElementTree as ET


def parse_xml(xml_file):
    """
    checks input file is XML file, parses it.
    :return: root of the xml or None
    """
    root = None
    try:
        if os.path.exists(xml_file):
            tree = ET.parse(xml_file)
            root = tree.getroot()
        else:
            logging.error('xml file: "{}" does not exist'.format(xml_file))
    except Exception as e:
        logging.error(e)

    return root


if __name__ == '__main__':
    FILE_ROOT = os.path.dirname(os.path.realpath(__file__))
    package_path = os.path.dirname(os.path.join(FILE_ROOT, '..', '..', ))
    test_data = os.path.join(package_path, 'test_data')
    example_file = os.path.join(test_data, "good_example.xml")

    data = parse_xml(xml_file=example_file)
    print(data)
