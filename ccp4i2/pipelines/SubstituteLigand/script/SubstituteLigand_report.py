from ....report.CCP4ReportParser import Report


class SubstituteLigand_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'SubstituteLigand'
    RUNNING = True

    def __init__(self, *args, **kws):
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        self.addDiv(style='clear:both;')
        
        #If there is POINTLESS tags in the XML, then the reflections have been through either aimless_pipe or
        #pointless_reindexToMatch
        reflectionNodes = self.xmlnode.findall('.//POINTLESS')
        if len(reflectionNodes) > 0:
            summaryFold = parent.addFold(label='Key reflection summary', brief='Reflections', initiallyOpen=True)
            pointlessNodes = self.xmlnode.findall('.//POINTLESS')
            if len(pointlessNodes) > 0:
                from ....wrappers.pointless.script.pointless_report import pointless_report
                pointlessreport = pointless_report(pointlessNodes[-1])
                pointlessreport.keyText(summaryFold)
            aimlessNodes = self.xmlnode.findall('.//AIMLESS')
            if len(aimlessNodes) > 0:
                from ....wrappers.aimless.script.aimless_report import aimless_report
                aimlessreport = aimless_report(aimlessNodes[-1],jobNumber='0')
                aimlessreport.keyText(None, parent=summaryFold)
    
        pmaNodes = self.xmlnode.findall('.//PhaserMrResults')
        if len(pmaNodes) > 0:
            from ...phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO_report import phaser_MR_AUTO_report
            pmaNode = pmaNodes[0]
            phaser_MRAReport = phaser_MR_AUTO_report(xmlnode=pmaNode, jobStatus='nooutput')
            if len(self.xmlnode.findall('.//PhaserMrSolutions/Solutions')) > 0:
                compareSolutionsFold = parent.addFold(label='Phaser results',initiallyOpen=True)
                phaser_MRAReport.addResults(parent=compareSolutionsFold)

        # Report here if Dimple's pointless run identified need for a reindexing
        reindexNodes = self.xmlnode.findall(".//REINDEX")
        if len(reindexNodes) > 0:
            newFold = parent.addFold(label="POINTLESS result", initiallyOpen=True)
            newFold.addPre(style="font-size:125%; font-color:red;", text="DIMPLE identified a need to reindex.")
            reindexText = "New reflection and FreeR (if given) have been output with operator {}".format(reindexNodes[0].text)
            newFold.addPre(style="font-size:125%; font-color:red;", text=reindexText)
                
        #phaser_MRAReport.drawContent(jobStatus=self.jobStatus, parent=self)
    
        refmacNodes = self.xmlnode.findall('.//REFMAC')
        if len(refmacNodes) > 0:
            from ....wrappers.refmac_i2.script.refmac_report import refmac_report
            refmacreport = refmac_report(xmlnode=refmacNodes[0], jobStatus='nooutput', jobInfo=self.jobInfo)
            refmacreport.addSummary(parent=self, withTables=False)
            objectMap = {}
            #Use DICTOUT or DICTIN if they have been harvested to define monomer geometry in pictures
            if self.jobInfo['filenames'].get('DICTOUT', None) is not None:
                objectMap['DICT'] = 'DICTOUT'
            elif self.jobInfo['filenames'].get('DICTIN', None) is not None:
                objectMap['DICT'] = 'DICTIN'
#FIXME - XML PICTURE
            """
            if not self.jobStatus.lower().count('running'):
                refmac_report.addRefinementPictures(parent=self, objectNameMap=objectMap, initiallyOpen=True)
            """
