#  Create xia2_multiplex.def.xml from PHIL parameters
#  Copyright (C) 2022 UKRI/STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on ideas and code by Nat Echols and Martin Noble.
#

import xml.etree.ElementTree as ET

from xia2.cli.multiplex import phil_scope

from ....utils.phil_handlers import PhilTaskCreator
from ccp4i2.core.CCP4Utils import printXml

class Xia2MultiplexTaskCreator(PhilTaskCreator):
    def __init__(self, debug=False):

        PhilTaskCreator.__init__(self, phil_scope, debug)
        self.fmt_dic["PLUGINNAME"] = "xia2_multiplex"

        self.inputDataXML = ET.fromstring(
            """
<container id="inputData">
  <content id="SEARCH_ROOT_DIR">
    <className>CDataFile</className>
    <qualifiers>
      <mustExist>True</mustExist>
      <allowUndefined>True</allowUndefined>
      <isDirectory>True</isDirectory>
      <guiLabel>Root directory</guiLabel>
      <toolTip>Start search from this directory</toolTip>
    </qualifiers>
  </content>

  <content id="DIALS_INTEGRATED">
    <className>CList</className>
    <qualifiers>
      <guiLabel>DIALS .refl</guiLabel>
      <toolTip>Must have an associated .expt file</toolTip>
      <listMinLength>2</listMinLength>
    </qualifiers>
    <subItem>
      <className>CDataFile</className>
      <qualifiers>
        <saveToDb>False</saveToDb>
        <mustExist>True</mustExist>
        <allowUndefined>True</allowUndefined>
        <isDirectory>False</isDirectory>
        <guiLabel>Bar</guiLabel>
        <toolTip>Foo</toolTip>
      </qualifiers>
    </subItem>
  </content>

  <content id="XIA2_RUN">
    <className>CList</className>
    <qualifiers>
        <guiLabel>Previous xia2 run directories</guiLabel>
        <listMinLength>2</listMinLength>
    </qualifiers>
    <subItem>
      <className>CDataFile</className>
      <qualifiers>
        <mustExist>True</mustExist>
        <allowUndefined>True</allowUndefined>
        <isDirectory>True</isDirectory>
        <guiLabel>Path to a previous xia2 run directory (containing a DataFiles sub-directory)</guiLabel>
        <toolTip>Integration files will be automatically extracted from a previous xia2 run</toolTip>
      </qualifiers>
    </subItem>
  </content>

  <qualifiers/>
</container>
"""
        )

        self.outputDataXML = ET.fromstring(
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
  <content id="PERFORMANCE">
    <className>CDataReductionPerformance</className>
    <qualifiers/>
  </content>
  <qualifiers/>
</container>
"""
        )

    def __call__(self):

        # Remove the option to set unit cell refinement. The geometry from
        # imported MTZs is suspect, so we don't want to refine cells here
        self._remove_elements_by_id(
            self.phil_tree,
            [
                "unit_cell",
            ],
        )

        if self.debug:
            printXml(self.phil_tree)

        super(Xia2MultiplexTaskCreator, self).__call__()


if __name__ == "__main__":

    x2mtc = Xia2MultiplexTaskCreator()
    x2mtc()
