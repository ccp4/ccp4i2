import base64
import os

from ccp4i2.report import Report


class molrep_selfrot_report(Report):
    TASKNAME = 'molrep_selfrot'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent=self

        for structureFactorNode in self.xmlnode.findall('.//StructureFactors'):
            parent.append('<br/>')
            parent.append('<br/>')
            parent.addText(text = 'Conclusion of search for anisotropy:',style='font-size:120%;')
            parent.addText(xmlnode=structureFactorNode, select='INFO',style='font-size:120%;')
            dataTable = parent.addTable(xmlnode = self.xmlnode, select='StructureFactors')
            headings = {'LowResProvided':'Dmax',
                'HighResProvided':'Dmin',
                'Completeness':'Completeness',
                'BOverall':'B-factor',
                'OpticalHighRes':'Optical Dmin',
                'EigenValueRatioH':'Aniso H',
                'EigenValueRatioK':'Aniso K',
                'EigenValueRatioL':'Aniso L'}
            for heading in headings:
                dataTable.addData(select=heading, title=headings[heading])

        for pattersonNode in self.xmlnode.findall ('.//Patterson'):
            #Pull out conclusion of search for translational symmetry
            parent.append('<br/>')
            parent.append('<br/>')
            parent.addText(text = 'Conclusion of search for translational symmetry:',style='font-size:120%;')
            parent.addText(xmlnode=pattersonNode, select='INFO',style='font-size:120%;')

            #Provide list of self Patterson peaks as a table
            parent.append('<br/>')
            pattersonFold=parent.addFold(label='Self Patterson peaks',initiallyOpen=False,style='font-size:110%;')
            peakTable = pattersonFold.addTable(xmlnode = pattersonNode, select='Peak',title='Self-Patterson peaks')
            headings = 'No IX  IY  IZ  Xfrac  Yfrac  Zfrac  Xort   Yort    Zort   Dens Dens_sigma'.split()
            for heading in headings:
                peakTable.addData(select=heading, title=heading)

        for rotationNode in self.xmlnode.findall ('.//SelfRotation'):
            #Provide list of symmetry-expanded rotation funciton peaks as a table
            parent.append('<br/>')
            selfRotationFold=parent.addFold(label='Self Rotation peaks',initiallyOpen=False,style='font-size:110%;')
            peakTable = selfRotationFold.addTable(xmlnode = rotationNode, select='Peak',title='Self-Patterson peaks')
            headings=  'No theta    phi     chi    alpha    beta   gamma      Rf    Rf_sigma'.split()
            for heading in headings:
                peakTable.addData(select=heading, title=heading)

        # Render self-rotation function plot as inline SVG
        self._addRotationFunctionPlot(parent)

    def _addRotationFunctionPlot(self, parent):
        """Convert molrep_rf.ps to SVG and embed in the report with a download link."""
        jobFolder = self.getJobFolder()
        if jobFolder is None:
            return

        ps_path = os.path.join(jobFolder, 'molrep_rf.ps')
        if not os.path.exists(ps_path):
            return

        try:
            from ccp4i2.wrappers.molrep_selfrot.script.ps2svg import ps_to_svg

            with open(ps_path, 'r') as f:
                ps_text = f.read()
            svg_text = ps_to_svg(ps_text)
        except Exception as e:
            parent.append(f'<p>Could not render rotation function plot: {e}</p>')
            return

        # Save SVG to job directory for standalone download
        svg_path = os.path.join(jobFolder, 'molrep_srf.svg')
        with open(svg_path, 'w') as f:
            f.write(svg_text)

        # Build a download link using a data URI
        svg_b64 = base64.b64encode(svg_text.encode('utf-8')).decode('ascii')
        download_link = (
            f'<a href="data:image/svg+xml;base64,{svg_b64}" '
            f'download="self_rotation_function.svg" '
            f'style="display:inline-block; margin:8px 0; padding:4px 12px; '
            f'background:#1976d2; color:white; border-radius:4px; '
            f'text-decoration:none; font-size:13px;">'
            f'Download SVG</a>'
        )

        parent.append('<br/>')
        srfFold = parent.addFold(label='Self Rotation Function Plot', initiallyOpen=True, style='font-size:110%;')
        srfFold.append(download_link)
        srfFold.append(svg_text)
