import os
import sys
import xml.etree.ElementTree as etree

from ccp4i2.report.svg import inline_svg
from ccp4i2.report.CCP4ReportParser import Report
from ccp4i2.wrappers.aimless.script.aimless_report import aimless_report
from ccp4i2.wrappers.adding_stats_to_mmcif_i2.script.adding_stats_to_mmcif_i2_report import (
    adding_stats_to_mmcif_i2_report,
)
from ccp4i2.wrappers.modelASUCheck.script.modelASUCheck_report import modelASUCheck_report


class prepare_deposition_report(Report):
    """What was done to the model (every sanitising fix, "found X, did Y"),
    the metadata written, the data statistics, the sequence coverage, and the
    wwPDB validation summary."""
    TASKNAME = 'prepare_deposition'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        if jobStatus in ['Running', 'Running remotely']:
            self.addDiv()
            self.addProblems(self)
            self.addSanitiseReport(self)
            self.addMetadataSummary(self)
            self.addSequenceAlignmentReport(self)
            if self.xmlnode.find('.//Validation') is not None:
                self.addText(text="Files sent to the wwPDB validation server; the report appears "
                                  "here when it is ready.", tone="emphasis")
        else:
            self.addDiv()
            self.addDefaultReport(self)

    def addDefaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addProblems(parent)
        self.addSanitiseReport(parent)
        self.addMetadataSummary(parent)
        if self.xmlnode.findall('.//AIMLESS'):
            aimlessNodes = self.xmlnode.findall('.//AIMLESS')
            aimlessreport = aimless_report(aimlessNodes[-1], jobStatus='nooutput',
                                           jobNumber=self.jobNumber)
            aimlessFold = parent.addFold(label="Data reduction statistics", initiallyOpen=True)
            aimlessreport.addAimlessSummary(parent=aimlessFold)
        elif self.xmlnode.findall('.//CifStats'):
            self.addCifStatsTable(parent)
        parent.addDiv()
        self.addSequenceAlignmentReport(parent)
        parent.addDiv()
        self.addValidationReport(parent)

    def addProblems(self, parent):
        problems = self.xmlnode.findall('.//Problems/Problem')
        if not problems:
            return
        fold = parent.addFold(label='Problems', initiallyOpen=True)
        fold.addText(text="The model cannot be deposited as it is:", tone="error")
        for p in problems:
            fold.addText(text=p.text or '', tone="error")

    def addSanitiseReport(self, parent):
        fixes = self.xmlnode.findall('.//Sanitise/Fix')
        if not fixes:
            return
        program = self.xmlnode.findtext('.//Sanitise/Program', 'unknown')
        fold = parent.addFold(label='What was done to the model (%s output)' % program,
                              initiallyOpen=True)
        table = fold.addTable(downloadable=True)
        table.addData(title="Found", data=[f.findtext('Found', '') for f in fixes])
        table.addData(title="Done", data=[f.findtext('Did', '') for f in fixes])

    def addMetadataSummary(self, parent):
        items = self.xmlnode.findall('.//Metadata/Item')
        if not items:
            return
        fold = parent.addFold(label='Metadata written into the model', initiallyOpen=False)
        table = fold.addTable(transpose=True, downloadable=True)
        for item in items:
            table.addData(title=item.get('label', ''), data=[item.text or ''])
        for w in self.xmlnode.findall('.//Metadata/Warning'):
            fold.addText(text=w.text or '', tone="warning")

    def addCifStatsTable(self, parent):
        cifStatsNode = self.xmlnode.findall('.//CifStats')[0]
        stats = cifStatsNode.findall('Stat')
        if not stats:
            return
        statsFold = parent.addFold(label="Data reduction statistics", initiallyOpen=True)
        table = statsFold.addTable(select="CifStats", transpose=True, downloadable=True)
        table.addData(title="", data=["Overall", "Inner shell", "Outer shell"])
        for stat in stats:
            table.addData(title=stat.get('label', ''),
                          data=[stat.get('overall', ''), stat.get('inner', ''), stat.get('outer', '')])

    def addSequenceAlignmentReport(self, parent):
        if not self.xmlnode.findall('.//SequenceAlignment/Alignment') \
                and not self.xmlnode.findall('.//SequenceAlignment/NonAlignedModelChains/Alignment'):
            return
        seqFold = parent.addFold(label='Sequence coverage', brief='Seq. Align.', initiallyOpen=True)
        seqReport = modelASUCheck_report(self.xmlnode, jobStatus='nooutput')
        seqReport.addAlignSummary(self.xmlnode, parent=seqFold)
        seqReport.addAlignReport(self.xmlnode, initiallyOpen=False, parent=seqFold)

    def addValidationReport(self, parent):
        # The wwPDB summary (slider SVG + percentile table) is drawn exactly
        # as the previous task drew it; its report class holds the SVG code.
        helper = adding_stats_to_mmcif_i2_report(self.xmlnode, jobInfo=self.jobInfo,
                                                 jobStatus='nooutput')
        try:
            validationXMLFilename = self.jobInfo['filenames']['VALIDATIONXML']
        except KeyError:
            validationXMLFilename = ''
        sys.stdout.flush()
        validationSvg = os.path.join(self.jobInfo['fileroot'], 'validation.svg')
        summary = []
        if validationXMLFilename and os.path.isfile(validationXMLFilename):
            with open(validationXMLFilename, 'r') as handle:
                validationXML = etree.fromstring(handle.read())
            summary = helper.parseXml(validationXML)
            remake = not os.path.exists(validationSvg)
            if not remake:
                with open(validationSvg, 'r') as handle:
                    remake = any('File download failed' in line for line in handle)
            if remake and summary:
                helper.makeValidationSvg(summary, validationSvg)
        if not summary and not os.path.exists(validationSvg):
            fold = parent.addFold(label="Validation", initiallyOpen=True)
            fold.addText(text="No validation report: the files were not sent to the wwPDB "
                              "validation server, or it returned nothing. OneDep validates "
                              "them again as part of deposition.")
            return
        summaryFold = parent.addFold(label="Validation summary", initiallyOpen=True)
        if os.path.exists(validationSvg):
            summaryFold.append(inline_svg(open(validationSvg).read()))
        else:
            summaryFold.addText(text="No validation summary was downloaded.")
        parent.addDiv()
        statisticsFold = parent.addFold(label="Validation statistics", initiallyOpen=True)
        if summary:
            table = statisticsFold.addTable(transpose=True, downloadable=True)
            table.addData(title="", data=["Value", "Percentile relative to all X-ray structures",
                                          "Percentile relative to X-ray structures of similar resolution"])
            for s in summary:
                if s['value'] is not None:
                    table.addData(title=s['label'], data=[s['value'], s['pc_abs'], s['pc_rel']])
        else:
            statisticsFold.addText(text="No validation statistics were downloaded.")
