import os
import shutil

from ....core.CCP4Modules import PROJECTSMANAGER
from ....report.CCP4ReportParser import Report


class pdb_redo_api_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'pdb_redo_api'
    RUNNING = True
    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )

        if jobStatus in ["Finished"]:
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)
            jobId = self.jobInfo.get("jobid", None)

        #FIXME - Need to copy test-page.html into job directory.

            jobDirectory = PROJECTSMANAGER().jobDirectory(jobId = jobId)
            shutil.copyfile(os.path.join(os.path.dirname(__file__),"test-page.html"),os.path.join(jobDirectory,self.xmlnode.findall('PDB_REDO_RESULTS_DIR')[0].text,"test-page.html"))

            pdbredourl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName="+os.path.join(self.xmlnode.findall('PDB_REDO_RESULTS_DIR')[0].text,"test-page.html")+"?jobNumber="
                + jobNumber
            )

            iframe_style = "display: block;background: #000; margin: 10px; border: none;height: 100vh; width: 95vw;"
            self.append('<iframe style="{1}" src="{0}"></iframe>'.format(pdbredourl,iframe_style))

        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            if len(xmlnode.findall('.//PDB_REDO_JOB_ID'))>0:
                jobNo = xmlnode.findall('.//PDB_REDO_JOB_ID')[0].text
                self.append("<p><b>PDB-REDO job {0} is currently running. Results will appear here when the job finishes.</b></p>".format(jobNo))
            else:
                self.append("<p><b>PDB-REDO job is currently running. Results will appear here when the job finishes.</b></p>")
            return

        htmlText = ""
        projectid = self.jobInfo.get("projectid", None)
        jobNumber = self.jobInfo.get("jobnumber", None)
        
        whatCheckFinalResultsHtml = os.path.normpath(os.path.join(self.xmlnode.findall('PDB_REDO_RESULTS_DIR')[0].text,"wf","index.html"))
        watchCheckFinalUrl = "/database/?getProjectJobFile?projectId=" + projectid + "?fileName="+whatCheckFinalResultsHtml+"?jobNumber=" + jobNumber

        whatCheckOriginalResultsHtml = os.path.normpath(os.path.join(self.xmlnode.findall('PDB_REDO_RESULTS_DIR')[0].text,"wo","index.html"))
        watchCheckOriginalUrl = "/database/?getProjectJobFile?projectId=" + projectid + "?fileName="+whatCheckOriginalResultsHtml+"?jobNumber=" + jobNumber

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        if len(xmlnode.findall('.//PDB_REDO_JOB_ID'))>0:
            jobNo = xmlnode.findall('.//PDB_REDO_JOB_ID')[0].text
            self.append("<p><b>Results for PDB-REDO job {0} <em>(job number on PDB-REDO web site).</em></b></p>".format(jobNo))
            self.append("<p>You can see a report for this job, including plots comparing these results with results for structures with similar resolutions, on the PDB_REDO website for 21 days.</p>".format(jobNo))

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        jobDir = self.jobInfo.get("fileroot", None)

        logFilesFold = self.addFold(label='Log files', brief='Log Files', initiallyOpen=True)

        processFilesFold = logFilesFold.addFold(label='PDB_REDO log', brief='PDB_REDO log', initiallyOpen=False)
        if len(xmlnode.findall('.//PDB_REDO_LOG_FILE'))>0 and len(xmlnode.findall('.//PDB_REDO_LOG_FILE'))>0:
            pdbLogFile = xmlnode.findall('.//PDB_REDO_LOG_FILE')[0].text
            with open(os.path.join(jobDir,pdbLogFile)) as f:
                logText = f.read()
                pdbLogDiv = processFilesFold.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
                pdbLogDiv.addPre(text = logText)

        refmacFilesFold = logFilesFold.addFold(label='Final refmac log', brief='Final refmac log', initiallyOpen=False)
        if len(xmlnode.findall('.//PDB_REDO_FINAL_REFMAC_LOG_FILE'))>0 and len(xmlnode.findall('.//PDB_REDO_FINAL_REFMAC_LOG_FILE'))>0:
            refmacLogFile = xmlnode.findall('.//PDB_REDO_FINAL_REFMAC_LOG_FILE')[0].text
            with open(os.path.join(jobDir,refmacLogFile)) as f:
                logText = f.read()
                refmacLogDiv = refmacFilesFold.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
                refmacLogDiv.addPre(text = logText)
