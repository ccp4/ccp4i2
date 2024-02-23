#
#  Copyright (C) 2024 UKRI/STFC Rutherford Appleton Laboratory, UK.
#
#  Author: Martin Maly, David Waterman
#

import os
from report.CCP4ReportParser import Report
import json
from collections import OrderedDict
from math import sqrt
import glob
from pathlib import Path


class xia2_ssx_reduce_report(Report):

    TASKNAME = "xia2_ssx_reduce"
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
        xia2SsxReduceLogFold = parent.addFold(
            label="xia2 ssx_reduce log text", initiallyOpen=True
        )
        xia2SsxReduceLogNode = self.xmlnode.findall("Xia2SsxReduceLog")[0]
        if xia2SsxReduceLogNode is not None:
            xia2SsxReduceLogFold.addPre(text=xia2SsxReduceLogNode.text)

    def defaultReport(self, parent=None):
        def create_I2Html(Html, I2Html):
            #Changes to pull jquery and plotly from local server as these cdn servers seem to block us
            new_lines = []
            with open(Html) as orig_html:
                old_lines = orig_html.readlines()
                for l in old_lines:
                    if "https://code.jquery.com/jquery-1.12.0.min.js" in l:
                        l = l.replace("https://code.jquery.com/jquery-1.12.0.min.js","/report_files/0.1.0/jquery.min.js")
                    if "https://cdn.plot.ly/plotly-latest.min.js" in l:
                        l = l.replace("https://cdn.plot.ly/plotly-latest.min.js","/report_files/0.1.0/plotly.js")
                    new_lines.append(l)
            with open(I2Html,"w") as new_html:
                new_html.writelines(new_lines)

        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title
        xia2SsxReduceErrorNode = None
        if len(self.xmlnode.findall("Xia2SsxReduceError"))>0:
            xia2SsxReduceErrorNode = self.xmlnode.findall("Xia2SsxReduceError")[0]
        projectid = self.jobInfo.get("projectid", None)
        jobNumber = self.jobInfo.get("jobnumber", None)
        xia2SsxReduceHtml = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "LogFiles", "dials.merge.html") # MM
        )
        xia2SsxReduceI2Html = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "dials.merge-i2.html") # MM
        )
        if os.path.exists(xia2SsxReduceHtml):
            create_I2Html(xia2SsxReduceHtml, xia2SsxReduceI2Html)
        if os.path.exists(xia2SsxReduceI2Html):
            xia2SsxReduceUrl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=dials.merge-i2.html?jobNumber="
                + jobNumber
            )
            xia2SsxReduceHtmlFold = parent.addFold(
                label="xia2.ssx_reduce report", initiallyOpen=True
            )
            xia2SsxReduceHtmlFold.append(
                '<span style="font-size:110%">Click on the '
                "following link to display the dials.merge.html report </span>"
            )
            xia2SsxReduceHtmlFold.append(
                '<a href="{0}">Open Results</a>'.format(xia2SsxReduceUrl)
            )

        # Link to dials.cosym.html if exists
        DialsCosymHtml = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "LogFiles", "dials.cosym.html") # MM
        )
        if os.path.exists(DialsCosymHtml):
            DialsCosymI2Html = os.path.normpath(
                os.path.join(self.jobInfo["fileroot"], "dials.cosym-i2.html") # MM
            )
            create_I2Html(DialsCosymHtml, DialsCosymI2Html)
            if os.path.exists(DialsCosymI2Html):
                DialsCosymUrl = (
                    "/database/?getProjectJobFile?projectId="
                    + projectid
                    + "?fileName=dials.cosym-i2.html?jobNumber="
                    + jobNumber
                )
                DialsCosymHtmlFold = parent.addFold(
                    label="dials.cosym report", initiallyOpen=True
                )
                DialsCosymHtmlFold.append(
                    '<span style="font-size:110%">Click on the '
                    "following link to display the dials.cosym.html report </span>"
                )
                DialsCosymHtmlFold.append(
                    '<a href="{0}">Open Results</a>'.format(DialsCosymUrl)
                )

        # When error occured
        if xia2SsxReduceErrorNode is not None:
            parent.addText(
                style="font-size:125%;color:red;",
                text="xia2.ssx_reduce exited with an error",
            )
            parent.addPre(text=xia2SsxReduceErrorNode.text)

        # Merging stats summary (tables gleaned from dials.merge*.json)
        DialsMergeJsonFilesPath = []
        LogFilesPath = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "LogFiles")
        )
        for JsonFilePath in glob.iglob(str(Path(LogFilesPath) / "*.json")):
            DialsMergePrefix = os.path.splitext(os.path.basename(JsonFilePath))[0]
            if "dials.merge" in DialsMergePrefix:
                DialsMergeJsonFilePath = os.path.normpath(
                    os.path.join(self.jobInfo["fileroot"], "LogFiles", DialsMergePrefix + ".json")
                )
                if os.path.isfile(DialsMergeJsonFilePath):
                    DialsMergeJsonFilesPath.append(DialsMergeJsonFilePath)
        for i, DialsMergeJsonPath in enumerate(DialsMergeJsonFilesPath):
            if len(DialsMergeJsonFilesPath) == 1:
                DialsMergeJsonLabel = "xia2.ssx_reduce summary"
            else:
                DialsMergeJsonLabel = "xia2.ssx_reduce summary - data set " + str(i + 1)
            runSummaryFold = self.addFold(
                label=DialsMergeJsonLabel, initiallyOpen=True
            )
            self.ssx_reduceSummaries(JsonFilePath=DialsMergeJsonPath, runSummaryFold=runSummaryFold, parent=parent)

        # xia2.ssx_reduce log
        if parent is None:
            parent = self
        xia2SsxReduceLogFold = parent.addFold(
            label="xia2 ssx_reduce log text", initiallyOpen=True
        )
        xia2SsxReduceLogNode = self.xmlnode.findall("Xia2SsxReduceLog")[0]
        if xia2SsxReduceLogNode is not None:
            xia2SsxReduceLogFold.addPre(text=xia2SsxReduceLogNode.text)

        # dials.merge* logs
        DialsMergeLogNodeMaster = self.xmlnode.findall("DialsMergeLogMaster")[0]
        for i in range(len(DialsMergeLogNodeMaster)):
            if len(DialsMergeLogNodeMaster) == 1:
                DialsMergeLogLabel = "dials.merge log text"
            else:
                DialsMergeLogLabel = "dials.merge log text - data set " + str(i + 1)
            DialsMergeLogFold = parent.addFold(
                label=DialsMergeLogLabel, initiallyOpen=False
            )
            DialsMergeLogNode = self.xmlnode.findall("DialsMergeLogMaster")[0][i]
            if DialsMergeLogNode is not None:
                DialsMergeLogFold.addPre(text=DialsMergeLogNode.text)

        # dials.cosym_reindex log
        DialsCosymLogNodeParent = self.xmlnode.findall("DialsCosymLog")
        print(str(DialsCosymLogNodeParent))
        if DialsCosymLogNodeParent != []:
            DialsCosymLogNode = self.xmlnode.findall("DialsCosymLog")[0]
            if DialsCosymLogNode is not None:
                DialsCosymLogFold = parent.addFold(
                    label="dials.cosym_reindex log text", initiallyOpen=False
                )
                DialsCosymLogFold.addPre(text=DialsCosymLogNode.text)

        return

    def ssx_reduceSummaries(self, JsonFilePath, runSummaryFold, parent=None):

        if parent is None:
            parent = self
        runTable = runSummaryFold.addTable(title="Summary table", transpose=True)

        names = []
        spaceGroups = []
        cells = []

        statisticDict = OrderedDict(
            [
                ("Low res.", []),
                ("High res.", []),
                ("#Obs", []),
                ("#Unique", []),
                ("Multiplicity", []),
                ("Completeness", []),
                ("I/sigI", []),
                ("CC1/2", []),
                ("Rsplit", []),
                # ("Anom Compl", []),
                ("Rpim (I)", []),
                ("Rpim (I+/-)", []),
                ("Rmeas (I)", []),
                ("Rmeas (I+/-)", []),
            ]
        )

        # Get data from dials.merge.json
        xia2json = JsonFilePath
        try:
            with open(xia2json, "r") as f:
                dic = json.loads(f.read())
        except IOError:
            dic = None
        if dic is None:
            parent.addText(
                style="font-size:125%;color:red;",
                text="File not found: dials.merge.json",
            )
            return

        # Extract relevant parts for the summary
        # datasets = dic["datasets"] # MM
        datasets = dic               # MM
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

            spaceGroups.append(self.xmlnode.findall("Xia2SsxReduceSG")[0].text)
            cells.append(self.xmlnode.findall("Xia2SsxReduceCell")[0].text)

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

            cc_half = _format_stat(merging_stats, "cc_one_half", ".3f")
            statisticDict["CC1/2"].append(cc_half)

            completeness = _format_stat(merging_stats, "completeness", ".1f")
            statisticDict["Completeness"].append(completeness)

            # anom_completeness = _format_stat(merging_stats, "anom_completeness", ".1f")
            # statisticDict["Anom Compl"].append(anom_completeness)

            r_pim = _format_stat(merging_stats, "r_pim", ".3f")
            statisticDict["Rpim (I)"].append(r_pim)

            r_pim_anom = _format_stat(merging_stats_anom, "r_pim", ".3f")
            statisticDict["Rpim (I+/-)"].append(r_pim_anom)

            r_meas = _format_stat(merging_stats, "r_meas", ".3f")
            statisticDict["Rmeas (I)"].append(r_meas)

            r_meas_anom = _format_stat(merging_stats_anom, "r_meas", ".3f")
            statisticDict["Rmeas (I+/-)"].append(r_meas_anom)

            r_split = _format_stat(merging_stats, "r_split", ".3f")
            statisticDict["Rsplit"].append(r_split)

            multiplicity = _format_stat(merging_stats, "multiplicity", ".2f")
            statisticDict["Multiplicity"].append(multiplicity)

            n_obs = _format_stat(merging_stats, "n_obs", "d")
            statisticDict["#Obs"].append(n_obs)

            n_uniq = _format_stat(merging_stats, "n_uniq", "d")
            statisticDict["#Unique"].append(n_uniq)

        runTable.addData(title="Wavelength", data=names)
        runTable.addData(title="Space group", data=spaceGroups)
        runTable.addData(title="Unit cell", data=cells)
        for k, v in statisticDict.items():
            runTable.addData(title=k, data=v)
        return
