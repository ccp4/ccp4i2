#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

import os
from report.CCP4ReportParser import Report
import json

from wrappers.xia2_dials.script import xia2_dials_report


class xia2_xds_report(xia2_dials_report.xia2_dials_report):

    TASKNAME = "xia2_xds"
    RUNNING = True
