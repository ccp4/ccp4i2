import os
import json
import xml.etree.ElementTree as ET

import sys
import pathlib
import os
import importlib.util
import inspect
import logging
from typing import Dict, Any, Type


CCP4I2_ROOT = str(os.environ["CCP4I2_ROOT"])


def find_defxml_files(root_dir):
    """Recursively find all .def.xml files from root_dir."""
    defxml_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".def.xml"):
                defxml_files.append(os.path.join(dirpath, filename))
    return defxml_files


def parse_defxml_file(file_path):
    """Parse .def.xml file and extract pluginName and pluginVersion."""
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        print(f"Parsing: {file_path}")
        print(f"Root tag: {root.tag}")
        print(f"Children tags: {[child.tag for child in root]}")
        # Try to find header with and without namespace
        header = root.find(".//{http://www.ccp4.ac.uk/ccp4ns}ccp4i2_header")
        if header is None:
            header = root.find(".//ccp4i2_header")
        if header is None:
            print(f"Header not found in {file_path}")
        else:
            print(f"Header found: {header.tag}")
            # Try to find pluginName and pluginVersion with and without namespace
            plugin_name = header.findtext("{http://www.ccp4.ac.uk/ccp4ns}pluginName")
            if plugin_name is None:
                plugin_name = header.findtext("pluginName", default="")
            plugin_version = header.findtext(
                "{http://www.ccp4.ac.uk/ccp4ns}pluginVersion"
            )
            if plugin_version is None:
                plugin_version = header.findtext("pluginVersion", default="")
            print(f"pluginName: {plugin_name}, pluginVersion: {plugin_version}")
            # Make file_path relative to this script's directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            rel_path = os.path.relpath(file_path, script_dir)
            return {
                "file_path": rel_path,
                "pluginName": plugin_name,
                "pluginVersion": plugin_version,
            }
    except Exception as e:
        print(f"Exception parsing {file_path}: {e}")
    return None


def trawl_defxml_files(root_dir=CCP4I2_ROOT):
    """Main function to trawl CCP4I2_ROOT for .def.xml files and extract info."""
    results = []
    def_xml_files = find_defxml_files(root_dir)
    print(def_xml_files)
    for file_path in def_xml_files:
        info = parse_defxml_file(file_path)
        print(file_path, info)
        if info:
            results.append(info)
    return results


if __name__ == "__main__":
    root_directory = CCP4I2_ROOT
    result = trawl_defxml_files(root_directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "defxml_lookup.json")
    with open(output_path, "w") as f:
        f.write(json.dumps(result, indent=2))
