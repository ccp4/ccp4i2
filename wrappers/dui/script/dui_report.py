import os
import glob
import re
import mmap
import json

from report.CCP4ReportParser import *

class dui_report(Report):

    TASKNAME = 'dui'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self,xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus,**kw)
        self.DUI_Outputlist = []
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        results = self.addResults()
        xia2HtmlFold = parent.addFold(label='Integration Reports', initiallyOpen=True)
        # Make sure we select the right directory for the dui_output folder
        if self.jobInfo['inputfiles']:
            annot = self.jobInfo['inputfiles'][0]['annotation'] # Assumes only one input file.
            reresult = re.search('\(([^\)]*)', annot)
            job_dloc = reresult.group(1)
            useDialsDir = os.path.join(os.path.split(os.path.normpath(self.jobInfo['fileroot']))[0],
                                       job_dloc, "dui_files")
        else:
            useDialsDir = os.path.join(self.jobInfo['fileroot'], "dui_files")
        self.ReadDuiFileListFromJson(useDialsDir)
        projectid = self.jobInfo.get("projectid", None)
        for mfile in self.DUI_Outputlist:
            # mfile[0] is the mtz file path, [1] the html report.
            rfilepath = mfile[1].get('report') # is the dict.
            jobNumber = os.path.basename(os.path.dirname(rfilepath))[4:]
            rfileurl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName="+os.path.basename(rfilepath)+"?jobNumber="
                + jobNumber
            )
            xia2HtmlFold.append('<br></br>')
            xia2HtmlFold.append('<span style="font-size:100%">Report for the integration stage of '
                                    + '<b>' + os.path.split(mfile[0])[1] + '</b>' + '</span>')
            xia2HtmlFold.append('<br></br>')
            xia2HtmlFold.append('<a href="{0}">Open Results</a>'.format(rfileurl))

    def ReadDuiFileListFromJson(self, useDialsDir=None):
        if not useDialsDir:
            return
        jsonin = os.path.join(useDialsDir, 'manifest.json')
        if not os.path.isfile(jsonin):
            return
        with open(jsonin, 'r') as fin:
            jfin = json.load(fin)
            self.DUI_Outputlist = jfin

