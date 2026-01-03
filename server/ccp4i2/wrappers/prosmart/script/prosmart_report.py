import os

from ccp4i2.core import CCP4Modules
from ccp4i2.report import Report


class prosmart_report(Report):
    TASKNAME='prosmart'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus in ["Finished"]:
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)
            jobId = self.jobInfo.get("jobid", None)

            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
            pdbredourl = (
                "/database/?getProjectJobFile?projectId=" + projectid
                + "?fileName=ProSMART_Results_i2.html"+"?jobNumber="
                + jobNumber
            )

            urlprefix= (
                "/database/?getProjectJobFile?projectId=" + projectid
                +"?jobNumber="
                + jobNumber
                + "?fileName="
            )

            t  = ""
            with open(os.path.join(jobDirectory,"ProSMART_Results.html")) as f:
                t = f.read()
            t = (t.replace('window.frames[iframeName].location = url','window.frames[iframeName].location = "'+urlprefix+'"+url')
                .replace("transformations.html","transformations_i2.html")
                .replace("pairwise.html","pairwise_i2.html")
                .replace("scores.html","scores_i2.html")
                .replace("pdb.html","pdb_i2.html")
                .replace("restrain.html","restrain_i2.html")
                .replace("sequence.html","sequence_i2.html"))

            with open(os.path.join(jobDirectory,"ProSMART_Results_i2.html"),"w+") as f:
                f.write(t)

            for fname in ["transformations","pairwise","scores","pdb","restrain","sequence"]:
                t  = ""
                with open(os.path.join(jobDirectory,"Output_Files",".HTML_files",fname+".html")) as f:
                    t = f.read()
                t = t.replace('window.frames[iframeName].location = url','window.frames[iframeName].location = "'+urlprefix+"Output_Files/.HTML_files/"+'"+url')
                with open(os.path.join(jobDirectory,"Output_Files",".HTML_files",fname+"_i2.html"),"w+") as f:
                    f.write(t)

            iframe_style = "display: block;background: #fff; margin: 10px; border: none;height: 100vh; width: 95vw;"
            self.append('<iframe style="{1}" src="{0}"></iframe>'.format(pdbredourl,iframe_style))

        self.addDiv(style="clear:both;")
