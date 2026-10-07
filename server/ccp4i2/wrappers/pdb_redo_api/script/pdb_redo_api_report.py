import os
import shutil

from ccp4i2.core import CCP4Modules
from ccp4i2.report import Report


# PDB-REDO's own measures of the run (its data.json, copied into program.xml
# by the wrapper): (label, before, after). A row shows only when the run
# reports both, so a protein's rotamers or an RNA's base pairs appear only
# where they apply.
REFINEMENT = (
    ("R", "RCAL", "RFIN"),
    ("R-free", "RFCAL", "RFFIN"),
    ("Bond length RMS Z-score", "OBRMSZ", "FBRMSZ"),
    ("Bond angle RMS Z-score", "OARMSZ", "FARMSZ"),
    ("Clashscore", "OCLASH", "FCLASH"),
)
PERCENTILES = (
    ("Ramachandran plot appearance", "TOZRAMA", "TFZRAMA"),
    ("Rotamer normality", "TOCHI12", "TFCHI12"),
    ("Coarse packing", "TOZPAK1", "TFZPAK1"),
    ("Fine packing", "TOZPAK2", "TFZPAK2"),
    ("Bump severity", "TOWBMPS", "TFWBMPS"),
    ("Hydrogen bond satisfaction", "TOHBSAT", "TFHBSAT"),
    ("Clashscore", "TOCLASH", "TFCLASH"),
    ("Dinucleotide conformation (CONFAL)", "TOCONFAL", "TFCONFAL"),
    ("Base pair conformation", "TOBPGRMSZ", "TFBPGRMSZ"),
)
CHANGES = (
    ("Rotamers changed", "NDROTA"),
    ("Side chains flipped", "HBFLIP"),
    ("Peptides flipped", "NBBFLIP"),
    ("Waters deleted", "NWATDEL"),
    ("Chiralities fixed", "NCHIRFX"),
    ("Residues fitting density better", "RSCCB"),
    ("Residues fitting density worse", "RSCCW"),
)


def metric(xmlnode, tag):
    """A PDB-REDO measure as text, or None when the run did not report it."""
    found = xmlnode.find('.//' + tag)
    if found is None or found.text in (None, '', 'None'):
        return None
    return found.text


class pdb_redo_api_report(Report):
    TASKNAME = 'pdb_redo_api'
    RUNNING = True
    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )

        if jobStatus in ["Finished"]:
            jobId = self.jobInfo.get("jobid", None)

        #FIXME - Need to copy test-page.html into job directory.

            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
            resultsDir = self.xmlnode.findall('PDB_REDO_RESULTS_DIR')[0].text
            shutil.copyfile(os.path.join(os.path.dirname(__file__),"test-page.html"),os.path.join(jobDirectory,resultsDir,"test-page.html"))

            testPagePath = os.path.join(jobDirectory, resultsDir, "test-page.html")
            if os.path.isfile(testPagePath):
                self.addFileLink(
                    label='Open PDB-REDO Results',
                    relativePath=os.path.join(resultsDir, "test-page.html"),
                    fileType='html',
                )

        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            if len(xmlnode.findall('.//PDB_REDO_JOB_ID'))>0:
                jobNo = xmlnode.findall('.//PDB_REDO_JOB_ID')[0].text
                self.append("<p><b>PDB-REDO job {0} is currently running. Results will appear here when the job finishes.</b></p>".format(jobNo))
            else:
                self.append("<p><b>PDB-REDO job is currently running. Results will appear here when the job finishes.</b></p>")
            return

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        if len(xmlnode.findall('.//PDB_REDO_JOB_ID'))>0:
            jobNo = xmlnode.findall('.//PDB_REDO_JOB_ID')[0].text
            self.append("<p><b>Results for PDB-REDO job {0} <em>(job number on PDB-REDO web site).</em></b></p>".format(jobNo))
            self.append("<p>You can see a report for this job, including plots comparing these results with results for structures with similar resolutions, on the PDB-REDO website for 21 days.</p>".format(jobNo))

        clearDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        self.addMetrics(xmlnode)

        jobDir = self.jobInfo.get("fileroot", None)

        logFilesFold = self.addFold(label='Log files', brief='Log Files', initiallyOpen=True)

        processFilesFold = logFilesFold.addFold(label='PDB-REDO log', brief='PDB-REDO log', initiallyOpen=False)
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


    def addMetrics(self, xmlnode):
        """Before and after, as PDB-REDO measured them; nothing if it reported none."""
        tables = []
        for heading, rows in (("Refinement", REFINEMENT),
                              ("Model quality percentile", PERCENTILES)):
            shown = [(label, metric(xmlnode, before), metric(xmlnode, after))
                     for label, before, after in rows]
            shown = [row for row in shown if row[1] is not None and row[2] is not None]
            if shown:
                tables.append(((heading, 'Input', 'PDB-REDO'), shown))
        shown = [(label, metric(xmlnode, tag)) for label, tag in CHANGES]
        shown = [row for row in shown if row[1] is not None]
        if shown:
            tables.append((('Model changes', 'Count'), shown))
        if not tables:
            return
        fold = self.addFold(label='What PDB-REDO changed', brief='Metrics', initiallyOpen=True)
        for titles, rows in tables:
            table = fold.addTable()
            for column, title in enumerate(titles):
                table.addData(title=title, data=[row[column] for row in rows])
