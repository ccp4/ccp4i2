from decimal import Decimal
import os
import sys
import xml.etree.ElementTree as ET

from ....report.CCP4ReportParser import Report
from ....wrappers.refmac_i2.script import refmac_report
from ....wrappers.aimless.script.aimless_report import aimless_report


class adding_stats_to_mmcif_i2_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'adding_stats_to_mmcif_i2'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        elif jobStatus in ['Running','Running remotely']:
            clearingDiv = self.addDiv(style="clear:both;")
            self.runningText(self)
        else:
            clearingDiv = self.addDiv(style="clear:both;")
            self.addDefaultReport(self)

    def runningText(self, parent=None) :
        if parent is None: parent = self
        parent.append ( "<p><b>Files sent to the Validation server. Updates will be shown here once validation report ready...</b></p>" )

    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        if len(self.xmlnode.findall('.//AIMLESS')) > 0:
            aimlessNodes = self.xmlnode.findall('.//AIMLESS')
            aimlessreport = aimless_report(aimlessNodes[-1], jobStatus='nooutput',jobNumber=self.jobNumber)
            aimlessFold = parent.addFold(label="Data reduction statistics", initiallyOpen=True)
            aimlessreport.addAimlessSummary(parent=aimlessFold)
        if len(self.xmlnode.findall("REFMAC")) > 0:
            refmacReportNode0 = self.xmlnode.findall("REFMAC")[0]
            refmacReport = refmac_report.refmac_report(xmlnode=refmacReportNode0, jobStatus='nooutput', jobInfo=self.jobInfo)
            refmacReport.addSummary(parent=parent, withTables=False)
        clearingDiv = parent.addDiv(style="clear:both;")

        clearingDiv = parent.addDiv(style="clear:both;")
        scaledUnmergedFold = parent.addFold(label="Scaled unmerged", initiallyOpen=True)
        if len(self.jobInfo['filenames']['SCALEDUNMERGED'].strip()) != 0:
            scaledUnmergedFold.append(f"Unmerged data were incorporated from file {self.jobInfo['filenames']['SCALEDUNMERGED']}")
        else:
            scaledUnmergedFold.append(f"No unmerged data were incorporated")

        try:
            validationXMLFilename = self.jobInfo['filenames']['VALIDATIONXML']
        except KeyError as err:
            validationXMLFilename = ''

        sys.stdout.flush()
        validationSvg = os.path.join(self.jobInfo['fileroot'],'validation.svg')
        remakeValidationSvg = False
        summary = []

        if validationXMLFilename != '' and os.path.isfile(validationXMLFilename): # we did send files to validation server:
            validationXML = ET.parse(validationXMLFilename).getroot()
            summary = self.parseXml(validationXML)
            if os.path.exists(validationSvg):
                with open(validationSvg,'r') as f:
                    for line in f:
                        if 'File download failed' in line:
                            remakeValidationSvg = True
            else:
                remakeValidationSvg = True

            if remakeValidationSvg and summary is not None:
                parent.append ( "<p>Failed to download validation summary. Regenerating svg...</p>" )
                self.makeValidationSvg(summary, validationSvg)

        summaryFold = parent.addFold(label="Validation summary", initiallyOpen=True)
        #htmlCode = '<img src="./validation.svg" alt="Diagram" style="margin-top:20px; margin-right: 30px; float:left;" />'
        if os.path.exists(validationSvg):
            summaryFold.append(open(validationSvg).read())
        else:
            summaryFold.append("<p>No validation SVG downloaded...</p>")

        clearingDiv = parent.addDiv(style="clear:both;")
        statisticsFold = parent.addFold(label="Validation statistics", initiallyOpen=True)
        if len(summary) != 0:
            table = statisticsFold.addTable(style='width:600px;', transpose=True, downloadable=True)
            headers = ["Value", "Percentile relative to all X-ray structures", "Percentile relative to X-ray structures of similar resolution"]
            table.addData(title="", data=headers)
            for s in summary:
                if s['value'] is not None:
                    table.addData(title = s['label'], data = [s['value'], s['pc_abs'], s['pc_rel']])
        else:
            statisticsFold.append ( "<p>No validation XML downloaded...</p>" )

    def pc_to_px(self, percentile):
      return 179.0+3.0*float(percentile)

    def draw_slider(self, idx, label, pc_abs, pc_rel, value):
        svg = '''\
<rect fill="url(#id4) " height="10" width="300" x="180" y="{0}"/>\
<rect height="14" style="stroke:black;fill:black;" width="3.0" x="{1}" y="{2}"/>\
<rect height="14" style="stroke:black;fill:none;" width="3.0" x="{3}" y="{2}"/>\
<text style="font-size:12px;text-anchor:end;" x="175" y="{4}">{5}</text>\
<text style="font-size:12px;text-anchor:start;" x="485" y="{4}"> {6}</text>\
'''.format(str(42+25*idx),str(self.pc_to_px(pc_abs)),str(40+25*idx),str(self.pc_to_px(pc_rel)),str(52+25*idx),label,value)
        return svg

    def add_lines(self, pcs):
        svg = ''
        for idx, pc in enumerate(pcs):
            if idx != len(pcs) - 1:
                svg += '<line style="stroke:black;" x1="{0}" x2="{1}" y1="{2}" y2="{3}"/>'.format(str(self.pc_to_px(pc)+1.5),str(self.pc_to_px(pcs[idx+1])+1.5),str(52+25*idx),str(67+25*idx))
        return svg

    def write_footer(self, num_sliders):
        y_coords = [str(y+25*(num_sliders)) for y in [42,50,58,66,74]]
        svg = '''\
<text style="font-size:14px;text-anchor:end;" x="180" y="28">Metric</text>\
<text style="font-size:14px;text-anchor:start;" x="480" y="28">Value</text>\
<text style="font-size:14px;text-anchor:middle;" x="330" y="28">Percentile Ranks</text>\
<text style="font-size:8px;text-anchor:start;font-style:italic;" x="180" y="{0}">Worse</text>\
<text style="font-size:8px;text-anchor:end;font-style:italic;" x="480" y="{0}">Better</text>\
<rect height="8" style="stroke:black;fill:black;" width="3" x="180" y="{1}"/>\
<text style="font-size:8px;text-anchor:start;" x="186" y="{2}">Percentile relative to all X-ray structures</text>\
<rect height="8" style="stroke:black;fill:none;" width="3" x="180" y="{3}"/>\
<text style="font-size:8px;text-anchor:start;" x="186" y="{4}">Percentile relative to X-ray structures of similar resolution</text></svg>\
'''.format(*y_coords)
        return svg

    def parseXml(self, xmlroot):
        rv = []
        for label, attrib in [('Rfree','DCC_Rfree'),
                              ('Clashscore', 'clashscore'), 
                              ('Ramachandran outliers', 'percent-rama-outliers'),
                              ('Sidechain outliers', 'percent-rota-outliers'),
                              ('RSRZ outliers', 'percent-RSRZ-outliers'),
                              ('RNA backbone', 'RNAsuiteness')
                              ]:
            rv.append({'label':label, 'value':xmlroot[0].get(attrib), 'pc_abs':xmlroot[0].get('absolute-percentile-{}'.format(attrib)), 'pc_rel':xmlroot[0].get('relative-percentile-{}'.format(attrib))})
        return rv

    def makeValidationSvg(self, sliders, filename):
        filtered_sliders = [i for i in sliders if i['value'] is not None]
        no_sliders = len(filtered_sliders)
        # print filtered_sliders
        svg = '''\
<svg xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink" baseProfile="full" height="{0}" version="1.1" width="530">\
<defs><linearGradient id="id4" x1="0" x2="1" y1="0" y2="0"><stop offset="0.0" stop-color="#ff0000"/><stop offset="0.1" stop-color="#ff5a5a"/>\
<stop offset="0.2" stop-color="#ffa0a0"/><stop offset="0.3" stop-color="#ffd2d2"/><stop offset="0.4" stop-color="#ffe0e0"/>\
<stop offset="0.5" stop-color="#eeeeee"/><stop offset="0.6" stop-color="#e0e0ff"/><stop offset="0.7" stop-color="#d2d2ff"/>\
<stop offset="0.8" stop-color="#a0a0ff"/><stop offset="0.9" stop-color="#5a5aff"/><stop offset="1.0" stop-color="#0000ff"/></linearGradient></defs>\
'''.format(str(83+25*no_sliders))
        percentiles = []
        for idx, slider in enumerate(filtered_sliders):
            if slider['label'] in ['Rfree','RNA backbone']:
                svg += self.draw_slider(idx, slider['label'], slider['pc_abs'], slider['pc_rel'], slider['value'])
            elif slider['label'] == 'Clashscore':
                svg += self.draw_slider(idx, slider['label'], slider['pc_abs'], slider['pc_rel'], str(Decimal(slider['value']).quantize(Decimal('1.'))))
            elif slider['label'] in ['Ramachandran outliers', 'Sidechain outliers','RSRZ outliers']:
                svg += self.draw_slider(idx, slider['label'], slider['pc_abs'], slider['pc_rel'], '{}%'.format(str(Decimal(slider['value']).quantize(Decimal('0.1')))) if slider['value'] != '0.00' else '0')
            percentiles.append(slider['pc_abs'])
        svg += self.add_lines(percentiles)
        svg += self.write_footer(no_sliders)
        with open(filename,'w') as f:
            f.write(svg)
