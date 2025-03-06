from report.CCP4ReportParser import *
import sys
import shutil
from core import CCP4Modules

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

            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
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
        """
        summaryfold = self.addFold(label='Validation metrics from PDB-REDO', brief='Metrics', initiallyOpen=True)

        tableDiv = summaryfold.addDiv(style="width:100%;border-width: 1px; border-color: black; margin:0px; padding:0px;")
        selectString = ".//

        title_data = []
        input_data = []
        pdb_redo_data = []

        if len(xmlnode.findall('.//RCAL'))>0 and len(xmlnode.findall('.//RFIN'))>0:
            input_data.append(xmlnode.findall('.//RCAL')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//RFIN')[0].text)
            title_data.append("R<sub>Work</sub>")

        if len(xmlnode.findall('.//RFCAL'))>0 and len(xmlnode.findall('.//RFFIN'))>0:
            input_data.append(xmlnode.findall('.//RFCAL')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//RFFIN')[0].text)
            title_data.append("R<sub>Free</sub>")

        if len(xmlnode.findall('.//OBRMSZ'))>0 and len(xmlnode.findall('.//FBRMSZ'))>0:
            input_data.append(xmlnode.findall('.//OBRMSZ')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//FBRMSZ')[0].text)
            title_data.append("Bond length RMS Z-score")

        if len(xmlnode.findall('.//OARMSZ'))>0 and len(xmlnode.findall('.//FARMSZ'))>0:
            input_data.append(xmlnode.findall('.//OARMSZ')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//FARMSZ')[0].text)
            title_data.append("Bond angle RMS Z-score")

        crystRefinementTable = tableDiv.addTable(transpose=False,select=selectString)
        crystRefinementTable.addData(title='Crystallographic refinement',data=title_data)
        crystRefinementTable.addData(title='Input',data=input_data)
        crystRefinementTable.addData(title='PDB-REDO',data=pdb_redo_data)

        title_data = []
        input_data = []
        pdb_redo_data = []

        if len(xmlnode.findall('.//TOZRAMA'))>0 and len(xmlnode.findall('.//TFZRAMA'))>0:
            input_data.append(xmlnode.findall('.//TOZRAMA')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFZRAMA')[0].text)
            title_data.append("Ramachandran plot appearance")

        if len(xmlnode.findall('.//TOCHI12'))>0 and len(xmlnode.findall('.//TFCHI12'))>0:
            input_data.append(xmlnode.findall('.//TOCHI12')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFCHI12')[0].text)
            title_data.append("Rotamer normality")

        if len(xmlnode.findall('.//TOZPAK1'))>0 and len(xmlnode.findall('.//TFZPAK1'))>0:
            input_data.append(xmlnode.findall('.//TOZPAK1')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFZPAK1')[0].text)
            title_data.append("Coarse packing")

        if len(xmlnode.findall('.//TOZPAK2'))>0 and len(xmlnode.findall('.//TFZPAK2'))>0:
            input_data.append(xmlnode.findall('.//TOZPAK2')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFZPAK2')[0].text)
            title_data.append("Fine packing")

        if len(xmlnode.findall('.//TOCONFAL'))>0 and len(xmlnode.findall('.//TFCONFAL'))>0:
            TOCONFAL = xmlnode.findall('.//TOCONFAL')[0].text
            TFCONFAL = xmlnode.findall('.//TFCONFAL')[0].text
            if TOCONFAL != "None" and TFCONFAL != "None":
                input_data.append(TOCONFAL)
                pdb_redo_data.append(TFCONFAL)
                title_data.append("Dinucleotide conformation (CONFAL)")

        if len(xmlnode.findall('.//TOBPGRMSZ'))>0 and len(xmlnode.findall('.//TFBPGRMSZ'))>0:
            TOBPGRMSZ = xmlnode.findall('.//TOBPGRMSZ')[0].text
            TFBPGRMSZ = xmlnode.findall('.//TFBPGRMSZ')[0].text
            if TOBPGRMSZ != "None" and TFBPGRMSZ != "None":
                input_data.append(TOBPGRMSZ)
                pdb_redo_data.append(TFBPGRMSZ)
                title_data.append("Base pair conformation")

        if len(xmlnode.findall('.//TOWBMPS'))>0 and len(xmlnode.findall('.//TFWBMPS'))>0:
            input_data.append(xmlnode.findall('.//TOWBMPS')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFWBMPS')[0].text)
            title_data.append("Bump severity")

        if len(xmlnode.findall('.//TOHBSAT'))>0 and len(xmlnode.findall('.//TFHBSAT'))>0:
            input_data.append(xmlnode.findall('.//TOHBSAT')[0].text)
            pdb_redo_data.append(xmlnode.findall('.//TFHBSAT')[0].text)
            title_data.append("Hydrogen bond satisfaction")

        modelQualityTable = tableDiv.addTable(transpose=False,select=selectString)
        title_data.append("WHAT_CHECK")
        input_data.append('<a href="{0}">Report</a>'.format(watchCheckOriginalUrl))
        pdb_redo_data.append('<a href="{0}">Report</a>'.format(watchCheckFinalUrl))
        modelQualityTable.addData(title='Model quality percentile',data=title_data)
        modelQualityTable.addData(title='Input',data=input_data)
        modelQualityTable.addData(title='PDB-REDO',data=pdb_redo_data)

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        changesFold = self.addFold(label='Significant model changes', brief='Model changes', initiallyOpen=True)

        tableChangesDiv = changesFold.addDiv(style="width:100%;border-width: 1px; border-color: black; margin:0px; padding:0px;")
        changesTable = tableChangesDiv.addTable(transpose=False,select=selectString)

        title_data = []
        changes_data = []

        if len(xmlnode.findall('.//NDROTA'))>0:
            changes_data.append(xmlnode.findall('.//NDROTA')[0].text)
            title_data.append("Rotamers changed")

        if len(xmlnode.findall('.//HBFLIP'))>0:
            changes_data.append(xmlnode.findall('.//HBFLIP')[0].text)
            title_data.append("Side chains flipped")

        if len(xmlnode.findall('.//NWATDEL'))>0:
            changes_data.append(xmlnode.findall('.//NWATDEL')[0].text)
            title_data.append("Waters deleted")

        if len(xmlnode.findall('.//NBBFLIP'))>0:
            changes_data.append(xmlnode.findall('.//NBBFLIP')[0].text)
            title_data.append("Peptides flipped")

        if len(xmlnode.findall('.//NCHIRFX'))>0:
            changes_data.append(xmlnode.findall('.//NCHIRFX')[0].text)
            title_data.append("Chiralities fixed")

        if len(xmlnode.findall('.//RSCCB'))>0:
            changes_data.append(xmlnode.findall('.//RSCCB')[0].text)
            title_data.append("Residues fitting density better")

        if len(xmlnode.findall('.//RSCCW'))>0:
            changes_data.append(xmlnode.findall('.//RSCCW')[0].text)
            title_data.append("Residues fitting density worse")

        changesTable.addData(title='Description',data=title_data)
        changesTable.addData(title='Count',data=changes_data)

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        """
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
