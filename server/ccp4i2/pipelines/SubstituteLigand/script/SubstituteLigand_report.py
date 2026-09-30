from ccp4i2.report import Report


class SubstituteLigand_report(Report):
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
                from ccp4i2.wrappers.pointless.script.pointless_report import pointless_report
                pointlessreport = pointless_report(pointlessNodes[-1])
                pointlessreport.keyText(summaryFold)
            aimlessNodes = self.xmlnode.findall('.//AIMLESS')
            if len(aimlessNodes) > 0:
                from ccp4i2.wrappers.aimless.script.aimless_report import aimless_report
                aimlessreport = aimless_report(aimlessNodes[-1],jobNumber='0')
                aimlessreport.keyText(None, parent=summaryFold)
    
        pmaNodes = self.xmlnode.findall('.//PhaserMrResults')
        if len(pmaNodes) > 0:
            pmaNode = pmaNodes[0]
            if pmaNode.find('Verdict') is not None or pmaNode.find('Modules') is not None:
                # The record of phaser_rnp_pipeline_phil's Phaser task
                from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil_report import (
                    phaser_mr_auto_phil_report,
                )
                phaserFold = parent.addFold(label='Phaser results', initiallyOpen=True)
                phaser_mr_auto_phil_report(xmlnode=pmaNode, jobStatus='nooutput').drawContent(
                    jobStatus=self.jobStatus, parent=phaserFold)
            else:
                # The classic wrapper's, in jobs run before the PHIL pipeline
                from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO_report import (
                    phaser_MR_AUTO_report,
                )
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
    
        # The refinement to report is the pipeline's last, servalcat's: its
        # model and maps are the job's outputs. This read the first REFMAC
        # node, which is Dimple's (or the Phaser pipeline's) intermediate run,
        # so the report showed statistics for a model the job did not output.
        servalcatNodes = self.xmlnode.findall('.//SERVALCAT_FIRST')
        refmacNodes = self.xmlnode.findall('.//REFMAC')
        if len(servalcatNodes) > 0:
            from ccp4i2.wrappers.servalcat.script.servalcat_report import servalcat_report
            servalcatReport = servalcat_report(
                xmlnode=servalcatNodes[-1], jobStatus='nooutput', jobInfo=self.jobInfo)
            cycle_data = servalcatReport.getCycleData(xmlnode=servalcatNodes[-1])
            summaryFold = parent.addFold(
                label='Summary of refinement', initiallyOpen=True, brief='Refinement')
            left, right = summaryFold.addTwoColumnLayout(left_span=5, right_span=7, spacing=2)
            servalcatReport.addTablePerCycle(cycle_data, parent=left, initialFinalOnly=True)
            servalcatReport.addGraphPerCycle(parent=right)
            parent.addDiv(style="clear:both;")
        elif len(refmacNodes) > 0:
            # Jobs from before the pipeline refined with servalcat
            from ccp4i2.wrappers.refmac.script.refmac_report import refmac_report
            refmac_report(xmlnode=refmacNodes[-1], jobStatus='nooutput',
                          jobInfo=self.jobInfo).addSummary(parent=parent, withTables=False)

        self.addLigandFit(parent)

    def addLigandFit(self, parent):
        fitNodes = self.xmlnode.findall('.//LIGAND_FIT')
        if len(fitNodes) == 0:
            return
        fit = fitNodes[-1]
        code = fit.findtext('Code', default='')
        placed = int(fit.findtext('Placed', default='0'))
        residues = [r.text for r in fit.findall('Residue')]
        fold = parent.addFold(label='Ligand fitting', initiallyOpen=True, brief='Ligand')
        if placed == 0:
            fold.addText(text=f'Coot found no site for {code} in the density. '
                              'The output model is the refined starting model.')
            return
        where = ', '.join(residues) if residues else 'the model'
        copies = 'copy' if placed == 1 else 'copies'
        fold.addText(text=f'Coot placed {placed} {copies} of {code}: {where}.')
        fold.addText(text='The placed ligand has not been refined. Check its fit '
                          'to the density, then refine the model with the '
                          'ligand dictionary from this job.')
