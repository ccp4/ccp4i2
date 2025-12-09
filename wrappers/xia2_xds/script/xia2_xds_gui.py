#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

from ccp4i2.baselayer import QtGui, QtWidgets, QtCore

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.core import CCP4Container
from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4Modules

from wrappers.xia2_dials.script import xia2_dials_gui


class xia2_xds_gui(xia2_dials_gui.xia2_dials_gui):

    # Subclass CTaskWidget to give specific task window
    TASKTITLE = "Automated integration of images with XDS using xia2"
    DESCRIPTION = "Select a directory containing images and integrate them"
    TASKNAME = "xia2_xds"
    TASKMODULE = "data_processing"
    TASKVERSION = 0.0
    RANK = 1
    SHORTTASKTITLE = "xia2/xds"
    TASKVERSION = 0.1

    def drawContents(self):

        # Try to find dependencies
        failed = []
        for cmd in ["xds", "xscale"]:
            tst = CCP4Modules.PREFERENCES().EXEPATHLIST.getExecutable(cmd)
            if CCP4Utils.which(tst) is None:
                failed.append(cmd)

        if len(failed) > 0:
            self.openFolder(folderFunction="inputData")
            for cmd in failed:
                self.createLine(["warning", "{0} has not been found.".format(cmd)])

            self.createLine(
                [
                    "warning",
                    (
                        "Please ensure the required software is "
                        "installed and either in the search path "
                        "or added to the Preferences, under Other "
                        "Software"
                    ),
                ]
            )
            return

        super(xia2_xds_gui, self).drawContents()

    @QtCore.Slot()
    def handleIndexMethod(self):
        # This is for DIALS only, so do nothing here
        pass
