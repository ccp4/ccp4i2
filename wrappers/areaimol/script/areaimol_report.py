from report.CCP4ReportParser import *
import sys
import base64

class areaimol_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'areaimol'
    RUNNING = False

    def addProgressGraph(self,parent,xmlnode,internalId="SummaryGraph",tag="SAS"):
        print("##################################################")
        print("##################################################")
        print("addProgressGraph")
        if len(xmlnode.findall(tag))>0:
            print("YES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            progressGraph = parent.addFlotGraph(title="SAS by atom",select=tag,style="height:250px; width:400px;float:left;border:0px;",outputXml=False,internalId=internalId)
            progressGraph.addData(title="Atom",    select="serNo")
            progressGraph.addData(title="Area",    select="area")
            plot = progressGraph.addPlotObject()
            plot.append('title','SAS by atom')
            plot.append('plottype','xy')
            plot.append('yrange', rightaxis='false')
            plot.append( 'xlabel', 'Atom' )
            plot.append( 'xintegral', 'true' )
            plot.append( 'ylabel', 'Area' )
            plotLine = plot.append('plotline',xcol=1,ycol=2,rightaxis='false',colour='blue')
        print("##################################################")
        print("##################################################")

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")

        if jobStatus not in ["Running", "Running remotely"]:
            try:
                summaryText = ""
                if len(self.xmlnode.findall(".//SummaryText"))>0:
                    xmlPath = './/SummaryText'
                    xmlNodes = self.xmlnode.findall(xmlPath)
                    for node in xmlNodes:
                        summaryText += base64.b64decode(node.text).decode()
                if summaryText:
                    fold = self.addFold(label="Summary", initiallyOpen=True)
                    fold.addPre(text=summaryText)
            except:
                pass

            """
            fold = self.addFold(label="SAS by atom", initiallyOpen=True)
            graphDiv = fold.addDiv(style='width:800px; height:270px;overflow:auto;')
            reportNode = self.xmlnode.findall('.//SASValues')[0]
            self.addProgressGraph(graphDiv,reportNode,internalId="SummaryGraph",tag="SAS")
            """

        fold = self.addFold(label="Areaimol log file")
        if len(self.xmlnode.findall(".//LogText"))>0:
            xmlPath = './/LogText'
            xmlNodes = self.xmlnode.findall(xmlPath)
            for node in xmlNodes:
                fold.addPre(text=base64.b64decode(node.text).decode())
