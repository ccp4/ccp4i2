import base64

from ccp4i2.report.CCP4ReportParser import Report


class areaimol_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'areaimol'
    RUNNING = False

    def addGraph(self,parent,xmlnode,internalId="SummaryGraph",tag="SAS"):
        if len(xmlnode.findall(tag))>0:
            progressGraph = parent.addFlotGraph(xmlnode=xmlnode, title="SAS by atom",select=tag,style="height:250px; width:400px;float:left;border:0px;",outputXml=False,internalId=internalId)
            progressGraph.addData(title="Atom",    select="serNo")
            progressGraph.addData(title="Area",    select="area")
            plot = progressGraph.addPlotObject()
            plot.append('title','SAS by atom')
            plot.append('plottype','xy')
            plotLine = plot.append('plotline',xcol=1,ycol=2)
            plotLine.append('colour','blue')
            plotLine.append('symbolsize','0')

    def addAreas(self):
        """The result first: accessible areas and, when two models were
        compared, the residues whose accessible area changes."""
        totals = [n.text for n in self.xmlnode.findall(".//Areas/Total")]
        if not totals:
            return
        fold = self.addFold(label="Accessible surface area", initiallyOpen=True)
        difference = self.xmlnode.findtext(".//Areas/Difference")
        if difference is None:
            fold.append("<p>Total accessible area: %s Å²</p>" % totals[0])
            return
        fold.append(
            "<p>Accessible area of the first model %s Å², of the second %s "
            "Å². Over the atoms the two share, the first has %s "
            "Å² %s accessible area: negative where something present only "
            "in the first model covers the surface.</p>" % (
                totals[0], totals[1] if len(totals) > 1 else "?",
                difference.lstrip("-"),
                "less" if difference.startswith("-") else "more"))
        if self.xmlnode.findall(".//Areas/Residue"):
            fold.append("<p>Residues whose accessible area changes, largest change first:</p>")
            table = fold.addTable(select=".//Areas", transpose=False, id="residue_differences")
            for title, select in (("Residue", "Residue/name"), ("Chain", "Residue/chain"),
                                  ("Number", "Residue/number"),
                                  ("Change in area (Å²)", "Residue/change")):
                table.addData(title=title, select=select)

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")

        if jobStatus not in ["Running", "Running remotely"]:
            self.addAreas()
            try:
                summaryText = ""
                if len(self.xmlnode.findall(".//SummaryText"))>0:
                    xmlPath = './/SummaryText'
                    xmlNodes = self.xmlnode.findall(xmlPath)
                    for node in xmlNodes:
                        summaryText += base64.b64decode(node.text).decode()
                if summaryText:
                    # The numbers are above; this is the program's own account.
                    fold = self.addFold(label="AREAIMOL summary",
                                        initiallyOpen=not self.xmlnode.findall(".//Areas/Total"))
                    fold.addPre(text=summaryText)
            except:
                pass

            breakDiv = self.addDiv(style="clear:both;")
            breakDiv.append('<br/>')
            fold = self.addFold(label="SAS by atom", initiallyOpen=True)
            graphDiv = fold.addDiv(style='width:800px; height:270px;overflow:auto;')
            reportNode = self.xmlnode.findall('.//SASValues')[0]
            self.addGraph(graphDiv,reportNode,internalId="SummaryGraph",tag="SAS")

        if jobStatus in ['Unknown','Interrupted','Failed','Unsatisfactory']:
            breakDiv = self.addDiv(style="clear:both;")
            breakDiv.append('<br/>')
            fold = self.addFold(label="Areaimol log file")
            if len(self.xmlnode.findall(".//LogText"))>0:
                xmlPath = './/LogText'
                xmlNodes = self.xmlnode.findall(xmlPath)
                for node in xmlNodes:
                    fold.addPre(text=base64.b64decode(node.text).decode())
