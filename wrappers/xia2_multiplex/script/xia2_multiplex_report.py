#
#  Copyright (C) 2022 UKRI/STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#

import os
from report.CCP4ReportParser import Report
import json
from collections import OrderedDict
from math import sqrt


class xia2_multiplex_report(Report):

    TASKNAME = "xia2_multiplex"
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
        xia2MultiplexLogFold = parent.addFold(
            label="xia2 multiplex log text", initiallyOpen=True
        )

        xia2MultiplexLogNode = self.xmlnode.findall("Xia2MultiplexLog")[0]
        if xia2MultiplexLogNode is not None:
            xia2MultiplexLogFold.addPre(text=xia2MultiplexLogNode.text)

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title
        xia2MultiplexErrorNode = None
        if len(self.xmlnode.findall("Xia2MultiplexError"))>0:
            xia2MultiplexErrorNode = self.xmlnode.findall("Xia2MultiplexError")[0]
        xia2MultiplexHtml = os.path.normpath(os.path.join(self.jobInfo["fileroot"], "xia2.multiplex.html"))
        xia2MultiplexI2Html = os.path.normpath(os.path.join(self.jobInfo["fileroot"], "xia2.multiplex-i2.html"))
        #Changes to pull jquery and plotly from local server as these cdn servers seem to block us
        new_lines = []
        with open(xia2MultiplexHtml) as orig_html:
            old_lines = orig_html.readlines()
            for l in old_lines:
                if "https://code.jquery.com/jquery-1.12.0.min.js" in l:
                    l = l.replace("https://code.jquery.com/jquery-1.12.0.min.js","/report_files/0.1.0/jquery.min.js")
                if "https://cdn.plot.ly/plotly-latest.min.js" in l:
                    l = l.replace("https://cdn.plot.ly/plotly-latest.min.js","/report_files/0.1.0/plotly.js")
                new_lines.append(l)
        with open(xia2MultiplexI2Html,"w") as new_html:
            new_html.writelines(new_lines)

        if os.path.exists(xia2MultiplexI2Html):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)

            xia2MultiplexUrl = (
                "/database/getProjectJobFile?projectId="
                + projectid
                + "&fileName=xia2.multiplex.html&jobNumber="
                + jobNumber
            )
            xia2MultiplexHtmlFold = parent.addFold(
                label="xia2.multiplex report", initiallyOpen=True
            )
            xia2MultiplexHtmlFold.append(
                '<span style="font-size:110%">Click on the '
                "following link to display the xia2.multiplex.html report </span>"
            )
            xia2MultiplexHtmlFold.append(
                '<a href="{0}">Open Results</a>'.format(xia2MultiplexUrl)
            )
        if xia2MultiplexErrorNode is not None:
            parent.addText(
                style="font-size:125%;color:red;",
                text="xia2.multiplex exited with an error",
            )
            parent.addPre(text=xia2MultiplexErrorNode.text)

        # Merging stats summary (a table gleaned from xia2.multiplex.json)
        self.multiplexSummaries(parent=parent)

        return

    def multiplexSummaries(self, parent=None):

        if parent is None:
            parent = self
        runSummaryFold = self.addFold(
            label="xia2.multiplex summary", initiallyOpen=True
        )
        runTable = runSummaryFold.addTable(title="Summary table", transpose=True)

        names = []
        spaceGroups = []
        cells = []

        statisticDict = OrderedDict(
            [
                ("Low res.", []),
                ("High res.", []),
                ("I/sigI", []),
                ("Completeness", []),
                ("Anom Compl", []),
                ("Rpim (I)", []),
                ("Rpim (I+/-)", []),
                ("Rmeas (I)", []),
                ("Rmeas (I+/-)", []),
                ("Multiplicity", []),
                ("#Obs", []),
                ("#Unique", []),
            ]
        )

        # Get data from xia2.multiplex.json
        xia2json = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "xia2.multiplex.json")
        )
        try:
            with open(xia2json, "r") as f:
                dic = json.loads(f.read())
        except IOError:
            dic = None
        if dic is None:
            parent.addText(
                style="font-size:125%;color:red;",
                text="File not found: xia2.multiplex.json",
            )
            return

        # Extract relevant parts for the summary
        datasets = dic["datasets"]
        dataset_ids = sorted(datasets.keys())

        def _format_stat(dataset, key, fmt, callback=None):
            res_bins = dataset
            overall = dataset["overall"]

            val = overall[key]
            if callback:
                val = callback(val)
            try:
                val = ("{:" + fmt + "}").format(val)
            except (TypeError, ValueError):
                val = "?"
            val_outer = res_bins[key][-1]
            if callback:
                val_outer = callback(val_outer)
            try:
                val_outer = (" ({:" + fmt + "})").format(val_outer)
                val += val_outer
            except (TypeError, ValueError):
                pass
            return val

        for k in dataset_ids:
            names.append(k)

            merging_stats = datasets[k]["merging_stats"]
            merging_stats_anom = datasets[k]["merging_stats_anom"]

            spaceGroups.append(self.xmlnode.findall("Xia2MultiplexSG")[0].text)
            cells.append(self.xmlnode.findall("Xia2MultiplexCell")[0].text)

            low_res = _format_stat(
                merging_stats, "d_star_sq_max", ".2f", lambda x: 1.0 / sqrt(x)
            )
            statisticDict["Low res."].append(low_res)

            hi_res = _format_stat(
                merging_stats, "d_star_sq_min", ".2f", lambda x: 1.0 / sqrt(x)
            )
            statisticDict["High res."].append(hi_res)

            ios = _format_stat(merging_stats, "i_over_sigma_mean", ".2f")
            statisticDict["I/sigI"].append(ios)

            completeness = _format_stat(merging_stats, "completeness", ".1f")
            statisticDict["Completeness"].append(completeness)

            anom_completeness = _format_stat(merging_stats, "anom_completeness", ".1f")
            statisticDict["Anom Compl"].append(anom_completeness)

            r_pim = _format_stat(merging_stats, "r_pim", ".3f")
            statisticDict["Rpim (I)"].append(r_pim)

            r_pim_anom = _format_stat(merging_stats_anom, "r_pim", ".3f")
            statisticDict["Rpim (I+/-)"].append(r_pim_anom)

            r_meas = _format_stat(merging_stats, "r_meas", ".3f")
            statisticDict["Rmeas (I)"].append(r_meas)

            r_meas_anom = _format_stat(merging_stats_anom, "r_meas", ".3f")
            statisticDict["Rmeas (I+/-)"].append(r_meas_anom)

            multiplicity = _format_stat(merging_stats, "multiplicity", ".3f")
            statisticDict["Multiplicity"].append(multiplicity)

            n_obs = _format_stat(merging_stats, "n_obs", "d")
            statisticDict["#Obs"].append(n_obs)

            n_uniq = _format_stat(merging_stats, "n_uniq", "d")
            statisticDict["#Unique"].append(n_uniq)

        runTable.addData(title="Run name", data=names)
        runTable.addData(title="Space group", data=spaceGroups)
        runTable.addData(title="Unit cell", data=cells)
        for k, v in statisticDict.items():
            runTable.addData(title=k, data=v)
        return
