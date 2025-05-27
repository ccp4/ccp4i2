#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on ideas and code by Nat Echols and Martin Noble.
#

from datetime import datetime
from io import StringIO
import getpass
import re
import socket
import sys
import time
import xml.etree.ElementTree as ET


class Phil2Etree(object):
    """
    Convert a phil_scope to an etree using container and content elements, as
    used by ccp4i2. Retain the hierarchical scope structure using containers.
    Attributes of phil scopes and definitions will be mapped to qualifiers as
    follows:

      short_caption --> guiLabel
      help          --> toolTip
      expert_level  --> guiDefinition / expertLevel
      style         --> guiDefinition / style
      caption       --> guiDefinition / caption
      multiple      --> guiDefinition / multiple

    At present, these attributes of phil scopes and definitions are not preserved:

      optional
      call
      sequential_format
      disable_add
      disable_delete
    """

    # define mappings for .type to <className>
    phil_type_as_class_name = {
        "str": "CString",
        "int": "CInt",
        "bool": "CBoolean",
        "ternary": "CString",
        "float": "CFloat",
        "choice": "CString",
        "path": "CDataFile",
    }

    def __init__(self, phil_scope):
        self.phil_scope = phil_scope

    def __call__(self, root_id):

        phil_params = ET.Element("container", id=root_id)
        self.convertScope(self.phil_scope, phil_params)
        return phil_params

    @staticmethod
    def sanitize_text(s):

        return s.replace("<", "&lt;").replace(">", "&gt;")

    def parse_choice_options(self, phil_def):

        s = ",".join([re.sub(r"\*", "", word.value) for word in phil_def.words])
        return self.sanitize_text(s)

    def parse_captions(self, phil_def):
        s = ",".join([re.sub("_", " ", item) for item in phil_def.caption.split()])
        return self.sanitize_text(s)

    def attributesToQualifiers(self, qualifiers, phil_obj, force_label=True):
        # guiLabel
        if phil_obj.short_caption is not None:
            gui_label = ET.SubElement(qualifiers, "guiLabel")
            gui_label.text = self.sanitize_text(phil_obj.short_caption)
        elif force_label:
            gui_label = ET.SubElement(qualifiers, "guiLabel")
            gui_label.text = self.sanitize_text(phil_obj.name)
        # toolTip
        if phil_obj.help is not None:
            toolTip = ET.SubElement(qualifiers, "toolTip")
            toolTip.text = self.sanitize_text(phil_obj.help)
        # guiDefinition
        guiDefinition = ET.SubElement(qualifiers, "guiDefinition")
        if phil_obj.expert_level is not None:
            expert_level = ET.SubElement(guiDefinition, "expertLevel")
            expert_level.text = str(phil_obj.expert_level)
        if phil_obj.style is not None:
            style = ET.SubElement(guiDefinition, "style")
            style.text = self.sanitize_text(str(phil_obj.style))
        if phil_obj.caption is not None:
            caption = ET.SubElement(guiDefinition, "caption")
            caption.text = self.sanitize_text(str(phil_obj.caption))
        if phil_obj.multiple is not None:
            multiple = ET.SubElement(guiDefinition, "multiple")
            multiple.text = self.sanitize_text(str(phil_obj.multiple))
        return

    def definitionToElement(self, keyword, phil_params):
        value = keyword.extract()
        elem = ET.SubElement(phil_params, "content")
        elem.set("id", keyword.full_path().replace(".", "__"))
        elem_class = ET.SubElement(elem, "className")

        # Map phil type to class and qualifiers
        phil_type = keyword.type.phil_type
        if phil_type == "bool" and str(value) not in ["True", "False"]:
            phil_type = "ternary"
        elem_class.text = self.phil_type_as_class_name.get(phil_type, "CString")
        qualifiers = ET.SubElement(elem, "qualifiers")
        self.attributesToQualifiers(qualifiers, keyword)

        # Set defaults for strings and bools
        if (phil_type in ["bool", "str"]) and (value is not None):
            default = ET.SubElement(qualifiers, "default")
            default.text = self.sanitize_text(str(value))

        # Set default for ternary logic, treated like a choice
        elif phil_type == "ternary":
            enum = ET.SubElement(qualifiers, "enumerators")
            enum.text = "True,False," + self.sanitize_text(str(value))
            default = ET.SubElement(qualifiers, "default")
            default.text = self.sanitize_text(str(value))
            only = ET.SubElement(qualifiers, "onlyEnumerators")
            only.text = "True"

        # Set attributes for choices
        elif phil_type == "choice":
            if keyword.type.multi:
                # enumerators do not map to PHIL's multi choice well. In that case
                # just use a string.
                default = ET.SubElement(qualifiers, "default")
                default.text = keyword.as_str().split("=")[1].strip()
            else:
                enum = ET.SubElement(qualifiers, "enumerators")
                enum.text = self.parse_choice_options(keyword)
                default = ET.SubElement(qualifiers, "default")
                if isinstance(value, list):
                    value = " ".join(["*" + self.sanitize_text(str(v)) for v in value])
                default.text = self.sanitize_text(str(value))
                only = ET.SubElement(qualifiers, "onlyEnumerators")
                only.text = "True"

        # Set defaults and limits for numerics
        elif phil_type in ["int", "float"]:
            default = ET.SubElement(qualifiers, "default")
            default.text = self.sanitize_text(str(value))
            if keyword.type.value_min is not None:
                min = ET.SubElement(qualifiers, "min")
                min.text = str(keyword.type.value_min)
            if keyword.type.value_max is not None:
                max = ET.SubElement(qualifiers, "max")
                max.text = str(keyword.type.value_max)
        return elem

    def convertScope(self, scope, container):
        for obj in scope.objects:
            if obj.is_definition:
                keywordElement = self.definitionToElement(obj, container)
                container.append(keywordElement)
            elif obj.is_scope:
                print("Object is scope", obj.full_path())
                sub_container = ET.Element(
                    "container", id=obj.full_path().replace(".", "__")
                )
                qualifiers = ET.SubElement(sub_container, "qualifiers")
                self.attributesToQualifiers(qualifiers, obj, force_label=False)
                container.append(sub_container)
                self.convertScope(obj, sub_container)
        return


class PhilTaskCreator(object):
    def __init__(self, phil_scope, debug=False):

        self.debug = debug

        # Translate the PHIL scope into an etree
        p2e = Phil2Etree(phil_scope)
        self.phil_tree = p2e(root_id="controlParameters")

        # Minimal boilerplate XML string
        self.boilerPlateXML = """<ccp4i2>
    <ccp4i2_header>
        <function>DEF</function>
        <userId>{USERID}</userId>
        <hostName>{HOSTNAME}</hostName>
        <creationTime>{CREATIONTIME}</creationTime>
        <pluginVersion></pluginVersion>
        <ccp4iVersion>{CCP4IVERSION}</ccp4iVersion>
        <pluginName>{PLUGINNAME}</pluginName>
        <OS>{OS}</OS>
        <jobId/>
    </ccp4i2_header>
    <ccp4i2_body id="{PLUGINNAME}">
        <container id="inputData">
        </container>
        <container id="outputData">
        </container>
        <container id="controlParameters">
        </container>
    </ccp4i2_body>
</ccp4i2>
"""

        self.fmt_dic = {}

        self.fmt_dic["USERID"] = getpass.getuser()

        self.fmt_dic["CREATIONTIME"] = datetime.fromtimestamp(time.time()).isoformat()

        self.fmt_dic["HOSTNAME"] = socket.gethostname()

        # FIXME, how to get this properly?
        self.fmt_dic["CCP4IVERSION"] = "0.0.1"

        self.fmt_dic["OS"] = sys.platform

        # overload by derived classes
        self.fmt_dic["PLUGINNAME"] = "PhilGUITask"
        self.inputDataXML = None
        self.outputDataXML = None

        return

    @staticmethod
    def _remove_elements_by_id(etree_xml, id_list):
        to_remove = []
        for cont in etree_xml.iter():
            if cont.get("id") in id_list:
                to_remove.append(cont)
        for cont in to_remove:
            cont.getparent().remove(cont)
        return etree_xml

    @staticmethod
    def _replace_element(etree_xml, element_id, new_element):

        for cont in etree_xml.iter():
            if cont.get("id") == element_id:
                parent = cont.getparent()
                parent.replace(cont, new_element)
        return

    def __call__(self):

        task_xml = ET.fromstring(self.boilerPlateXML.format(**self.fmt_dic))

        # Insert inputData
        if self.inputDataXML is not None:
            self._replace_element(task_xml, "inputData", self.inputDataXML)

        # Replace controlParameters in task_xml with the phil_tree
        self._replace_element(task_xml, "controlParameters", self.phil_tree)

        # Insert outputData
        if self.outputDataXML is not None:
            self._replace_element(task_xml, "outputData", self.outputDataXML)

        # Write out prettified version
        out_file = self.fmt_dic["PLUGINNAME"] + ".def.xml"
        parser = ET.XMLParser(remove_blank_text=True)
        tree = ET.parse(StringIO(ET.tostring(task_xml).decode("utf-8")), parser)
        ET.indent(tree)
        with open(out_file, "w") as f:
            f.write(ET.tostring(tree, xml_declaration=True).decode("utf-8"))
        return
