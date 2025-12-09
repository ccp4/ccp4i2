#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

import os
from ccp4i2.report.CCP4ReportParser import Report
import json
import re


class xia2_dials_report(Report):

    TASKNAME = "xia2_dials"
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        job_stat = kw.get("jobStatus", None)
        if job_stat is not None and job_stat.lower() == "nooutput":
            return
        elif job_stat is not None and job_stat.lower() == "running":
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)

    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        xia2SummaryFold = parent.addFold(label="xia2 text", initiallyOpen=True)
        xia2TxtNode = self.xmlnode.findall("Xia2Txt")[0]
        if xia2TxtNode is not None:
            xia2SummaryFold.addPre(text=xia2TxtNode.text)

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title
        xia2TxtNode = self.xmlnode.findall("Xia2Txt")[0]
        xia2html = os.path.normpath(os.path.join(self.jobInfo["fileroot"], "xia2.html"))
        if os.path.exists(xia2html):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)

            xia2url = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=xia2.html?jobNumber="
                + jobNumber
            )
            xia2HtmlFold = parent.addFold(label="xia2 report", initiallyOpen=True)
            xia2HtmlFold.append(
                '<span style="font-size:110%">Click on the '
                "following link to display the xia2.html report </span>"
            )
            xia2HtmlFold.append('<a href="{0}">Open Results</a>'.format(xia2url))
        if len(self.xmlnode.findall("Xia2Error"))>0:
            xia2ErrorNode = self.xmlnode.findall("Xia2Error")[0]
            parent.addText(
                style="font-size:125%;color:red;", text="xia2 exited with an error"
            )
            if xia2TxtNode is not None:
                parent.addPre(text=xia2TxtNode.text)
            parent.addPre(text=xia2ErrorNode.text)
        elif xia2TxtNode is not None:
            xia2SummaryFold = parent.addFold(label="xia2 Text", initiallyOpen=False)
            xia2SummaryFold.addPre(text=xia2TxtNode.text)

        # xia2 summary (a table gleaned from xia2.json)
        self.xia2Summaries(parent=parent)

        return

    # Datasets are identified in xia2.json's scalr_statistics dictionary by
    # keys of the form u'["AUTOMATIC", "DEFAULT", "WAVE2"]'. Clean this form up
    # for display in the table
    @staticmethod
    def _clean_scalr_statistics_ids(dataset_id):
        dataset_id = str(dataset_id).split(",")
        dataset_id = [re.sub('[\[\]"]', "", e) for e in dataset_id]
        return "/".join(dataset_id)

    def xia2Summaries(self, parent=None):

        if parent is None:
            parent = self
        runSummaryFold = self.addFold(label="xia2 summary", initiallyOpen=True)
        runTable = runSummaryFold.addTable(title="Summary table", transpose=True)

        names = []
        spaceGroups = []
        cells = []

        from collections import namedtuple

        Stat = namedtuple("Stat", ["key", "title", "fmt", "values"])
        statistics = [
            Stat("High resolution limit", "High res.", "%.2f", []),
            Stat("Low resolution limit", "Low res.", "%.2f", []),
            Stat("Completeness", "Completeness", "%.1f", []),
            Stat("Multiplicity", "Multiplicity", "%.1f", []),
            Stat("I/sigma", "I/sigI", "%.2f", []),
            Stat("Rmeas(I)", "Rmeas (I)", "%.3f", []),
            Stat("Rmeas(I+/-)", "Rmeas (I+/-)", "%.3f", []),
            Stat("Rpim(I)", "Rpim (I)", "%.3f", []),
            Stat("Rpim(I+/-)", "Rpim (I+/-)", "%.3f", []),
            Stat("CC half", "CC½", "%.3f", []),
            Stat("Anomalous completeness", "Anom. Compl.", "%.1f", []),
            Stat("Anomalous multiplicity", "Anom. Mult.", "%.1f", []),
            Stat("Total observations", "#Obs", "%d", []),
            Stat("Total unique", "#Unique", "%d", []),
        ]

        # Get data from xia2.json
        xia2json = os.path.normpath(os.path.join(self.jobInfo["fileroot"], "xia2.json"))
        try:
            with open(xia2json, "r") as f:
                dic = json.loads(f.read())
        except IOError:
            dic = None
        if dic is None:
            parent.addText(
                style="font-size:125%;color:red;", text="File not found: xia2.json"
            )
            return

        # Extract relevant parts for the summary
        crystals = dic["_crystals"]
        if len(crystals) > 1:  # don't cope with > 1 crystal here
            parent.addText(
                style="font-size:125%;color:red;",
                text=("xia2 summary table can only be created " "for a single crystal"),
            )
            return

        k = list(dic["_crystals"].keys())[0]
        dat = dic["_crystals"][k]

        dataset_ids = sorted(dat["_scaler"]["_scalr_statistics"].keys())
        for k in dataset_ids:
            names.append(self._clean_scalr_statistics_ids(k))
            stat = dat["_scaler"]["_scalr_statistics"][k]
            spaceGroups.append(", ".join((dat["_scaler"]["_scalr_likely_spacegroups"])))
            cell = dat["_scaler"]["_scalr_cell"]
            cells.append(
                "{:.2f}, {:.2f}, {:.2f}<br/>{:.2f}, {:.2f}, {:.2f}".format(*cell)
            )

            for row in statistics:
                try:
                    overallValue = row.fmt % stat[row.key][0]
                    outerValue = row.fmt % stat[row.key][2]
                    msg = overallValue + " (" + outerValue + ")"
                except KeyError:
                    msg = "Not known"
                row.values.append(msg)

        runTable.addData(title="Run name", data=names)
        runTable.addData(title="Probable space groups", data=spaceGroups)
        runTable.addData(title="Unit cell", data=cells)
        for row in statistics:
            runTable.addData(title=row.title, data=row.values)
        return
