from __future__ import print_function
from ccp4i2.report.CCP4ReportParser import *
import sys
import os
from lxml import etree

class AUSPEX_report(Report):

    TASKNAME = 'AUSPEX'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'):
            self.outputXml = False
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        print(self.jobInfo)
        for fname in self.jobInfo['filenames']['IM_OUT']:
            print(fname)
        results = self.addResults()

        def imgToUrl(img):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)

            imgUrl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName="+img+"?jobNumber="
                + jobNumber
                )
            return imgUrl

        inst = "<h2>AUSPEX Plots</h2>"
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "I_plot.png":
                inst += "<h3>Plot of intensities against resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "SigI_plot.png":
                inst += "<h3>Plot of sigma (intensity) against resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "ISigI_plot.png":
                inst += "<h3>Plot of intensity / sigma against resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "F_plot.png":
                inst += "<h3>Plot of amplitudes against resolution.</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:            
            if os.path.basename(img) == "SigF_plot.png":
                inst += "<h3>Plot of sigma (intensity) against resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "FSigF_plot.png":
                inst += "<h3>Plot of amplitude / sigma against resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "score.png":
                inst += "<h3>Plot of icefinder score and intensities vs resolution:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "intensities.png":
                inst += "<h3>Plots for intensities:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        for img in self.jobInfo['filenames']['IM_OUT']:
            if os.path.basename(img) == "amplitudes.png":
                inst += "<h3>Plots for amplitudes:</h3>"
                inst += "<img width=\"650px\" src=\"%s\"></img>"%(imgToUrl(img))
        parent.append(inst)

