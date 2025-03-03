from __future__ import print_function

#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on ideas and code by Nat Echols and Martin Noble.
#

"""Create xia2_dials.def.xml from PHIL parameters"""

import sys
import os
from lxml import etree

# Nasty trick required to import PhilTaskCreator when running with ccp4-python
this_dir = os.path.dirname(os.path.realpath(__file__))
ccp4i2_dir = os.path.dirname(os.path.dirname(os.path.dirname(this_dir)))
sys.path.append(ccp4i2_dir)
from utils.phil_handlers import PhilTaskCreator


class Xia2DialsTaskCreator(PhilTaskCreator):
    def __init__(self, debug=False):

        from xia2.Handlers.Phil import master_phil

        PhilTaskCreator.__init__(self, master_phil, debug)
        self.fmt_dic["PLUGINNAME"] = "xia2_dials"

        # Specifically remove:
        # 1) elements related to XDS only
        # 2) the option to set a pipeline other than DIALS
        # 3) settings only relevant to xia2.strategy
        # 4) the image= option as this is handled explicitly in inputData
        self._elts_to_remove = [
            "xds",
            "xia2__settings__xds",
            "xia2__settings__xds_cell_deviation",
            "xia2__settings__pipeline",
            "strategy",
            "xia2__settings__input__image",
        ]

        self.inputDataXML = etree.fromstring(
            """
<container id="inputData">
  <content id="IMAGE_FILE">
    <className>CXia2ImageSelectionList</className>
    <qualifiers>
      <listMinLength>1</listMinLength>
    </qualifiers>
    <subItem>
      <className>CXia2ImageSelection</className>
      <qualifiers>
        <guiLabel>Image file path</guiLabel>
        <toolTip>xia2 will automatically find other images in the dataset, optionally limited to the specified range</toolTip>
      </qualifiers>
    </subItem>
  </content>
  <content id="IMAGE_DIRECTORY">
    <className>CDataFile</className>
    <qualifiers>
      <mustExist>True</mustExist>
      <allowUndefined>True</allowUndefined>
      <isDirectory>True</isDirectory>
    </qualifiers>
  </content>
  <qualifiers/>
</container>
"""
        )

        self.outputDataXML = etree.fromstring(
            """
<container id="outputData">
  <content id="UNMERGEDOUT">
    <className>CList</className>
    <qualifiers/>
    <subItem>
      <className>CUnmergedDataFile</className>
      <qualifiers>
        <contentFlag>
          <min>0</min>
        </contentFlag>
        <baseName>
          <allowedCharacters>*?</allowedCharacters>
        </baseName>
        <relPath>
          <allowedCharacters>*?</allowedCharacters>
        </relPath>
      </qualifiers>
    </subItem>
  </content>
  <content id="HKLOUT">
    <className>CList</className>
    <qualifiers/>
    <subItem>
      <className>CObsDataFile</className>
      <qualifiers>
        <default/>
        <contentFlag>
          <min>0</min>
        </contentFlag>
        <subType>
          <menuText>observed data,derived data,reference data</menuText>
          <onlyEnumerators>True</onlyEnumerators>
          <enumerators>1,2,3</enumerators>
          <default>1</default>
        </subType>
      </qualifiers>
    </subItem>
  </content>
  <content id="FREEROUT">
    <className>CList</className>
    <qualifiers/>
    <subItem>
      <className>CFreeRDataFile</className>
      <qualifiers>
        <contentFlag>
          <min>0</min>
        </contentFlag>
        <subType>
          <onlyEnumerators>True</onlyEnumerators>
        </subType>
      </qualifiers>
    </subItem>
  </content>
  <content id="DIALSJOUT">
    <className>CList</className>
    <subItem>
      <className>CDialsJsonFile</className>
    </subItem>
  </content>
  <content id="DIALSPOUT">
    <className>CList</className>
    <subItem>
      <className>CDialsPickleFile</className>
    </subItem>
  </content>
  <content id="PERFORMANCE">
    <className>CDataReductionPerformance</className>
    <qualifiers/>
  </content>
  <qualifiers/>
</container>
"""
        )

    def __call__(self):

        self._remove_elements_by_id(self.phil_tree, self._elts_to_remove)

        if self.debug:
            print(etree.tostring(self.phil_tree, pretty_print=True))

        super(Xia2DialsTaskCreator, self).__call__()


if __name__ == "__main__":

    x2dtc = Xia2DialsTaskCreator()
    x2dtc()
