#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

import os

from lxml import etree

from ....core.CCP4PluginScript import CPluginScript
from ...xia2_dials.script import xia2_dials


class Cxia2_xds(xia2_dials.Cxia2_dials):

    TASKTITLE = "Data processing with xia2/xds"
    TASKNAME = "xia2_xds"

    def makeCommandAndScript(self):
        par = self.container.controlParameters
        inp = self.container.inputData

        # Create PHIL file and command line
        self._setCommandLineCore(phil_filename="xia2_xds.phil")

        self.xmlroot = etree.Element("Xia2Xds")

        self.watchFile(
            os.path.normpath(os.path.join(self.getWorkDirectory(), "xia2.txt")),
            self.handleXia2DotTxtChanged,
        )

        return CPluginScript.SUCCEEDED

    def _collect_pickles_and_jsons(self):
        # override to do nothing
        pass

    @staticmethod
    def _get_annotation(prefix, suffix):
        """Form suitable annotation strings"""
        return prefix + " from XDS integration of " + suffix
