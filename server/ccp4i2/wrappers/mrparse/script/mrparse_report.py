
import os
import xml.etree.ElementTree as etree
from pathlib import Path

from ccp4i2.report import Report
from ccp4i2.report.embedded_assets import localise_report_assets, vendored_asset


def parse_from_unicode(unicode_str):
    utf8_parser = etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return etree.fromstring(s, parser=utf8_parser)

class mrparse_report(Report):

    TASKNAME = 'mrparse'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'):
            self.outputXml = False
        self.defaultReport()
        return

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.append("<p>Finished running MrParse</p>")
        basepath = self.jobInfo['fileroot']
        mrparse_rep = os.path.join(basepath, "mrparse_0", 'mrparse.html')

        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        if not os.path.exists(mrparse_rep):
            ResultsI2Folder.append('<p>MrParse report not found</p>')
            return

        # MrParse writes its stylesheets and scripts as absolute paths into its
        # own installation. Copy them in beside the report and rewrite the
        # references relative to it, so the report renders when served over
        # HTTP and survives project export, import and relocation.
        #
        # This replaces a rewrite that ran on Windows only (elsewhere the Qt
        # app opened the file straight from disk, where absolute paths worked)
        # and pointed at /database/projectid/N/jobnumber/M/file/... — a route
        # that no longer exists, carrying a project identity that changes when
        # a project is imported somewhere else.
        localised = localise_report_assets(
            Path(mrparse_rep),
            marker='mrparse/html',
            subdirectory='mrparse_html',
            # MrParse loads d3 from d3js.org. Serve our own copy: the report
            # is an archival record that should still draw its feature viewer
            # offline, and on a machine with no route to the internet.
            extra_assets={
                'https://d3js.org/d3.v5.min.js': vendored_asset('d3.v5.min.js'),
            },
        )
        report_name = localised.name if localised else 'mrparse.html'

        ResultsI2Folder.append('<span style="font-size:110%">Click on the '
                               'following link to display the browser report '
                               'for the MrParse job</span>')
        ResultsI2Folder.addFileLink(
            label='Open MrParse Results',
            relativePath=f'mrparse_0/{report_name}',
            fileType='html',
        )

#FIXME - XML PICTURE
        return
        mrparse_xml = os.path.join(basepath, 'params.xml')
        with open(mrparse_xml) as f:
            t = f.read()
            tree = parse_from_unicode(t)
            pdbs = tree.findall(".//ccp4i2_body/outputData/XYZOUT/CPdbDataFile")
            pictureFold = self.addFold(label='Picture', initiallyOpen=False)
            pictureFold.addText(text='View of the models')

            pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='550px', tableWidth='260px', contentWidth='450px')
            pdbidx = 1
            for pdb in pdbs:
                scene = """<?xml version='1.0'?>
<scene>
    <data>"""
                baseName = pdb.findall("baseName")[0].text
                scene += "<MolData id='id{0}'>\n".format(pdbidx)
                scene += '<filename>{0}</filename>\n'.format(os.path.join(basepath,baseName))
                scene += "</MolData>\n"
                scene += """</data>
  <View>
     <scale_auto>1</scale_auto>
     <slab_enabled>0</slab_enabled>
     <centre_MolData>id1</centre_MolData>
     <centre_selection>all</centre_selection>
     <!-- <scale_auto_contacts>1</scale_auto_contacts> -->
     <orientation_auto>
       <selection>all</selection>
       <molData>id1</molData>
     </orientation_auto>
     <slab_enabled>0</slab_enabled>
  </View>
  <wizard><template>Ribbons:colour chains</template>
    <parameters>"""
                scene += "<MolData{0}>id{0}</MolData{0}>\n".format(pdbidx)
                pdbidx += 1
                scene += """</parameters>
  </wizard>
</scene>"""
                pic = pictureGallery.addPicture(label=baseName, scene=scene)

        return
