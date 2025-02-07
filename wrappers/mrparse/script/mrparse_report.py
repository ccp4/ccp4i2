
from report.CCP4ReportParser import *
from core.CCP4Modules import PREFERENCES
import os
import re

import xml.etree.ElementTree as etree

def parse_from_unicode(unicode_str):
    utf8_parser = etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return etree.fromstring(s, parser=utf8_parser)

def change_i2_mrparse_paths(lines):
    """
    This changes the urls for css,js,svg to root of i2 server. These files have to
    be served by i2 on Windows it seems because of stricter browser security.
    """
    svg = ["MrParse-logo-tight.svg"]
    js = ["vue.min.js","lodash.min.js","jquery.js","feature-viewer.bundle.js","mrparse_vue.js","https://d3js.org/d3.v5.min.js"]
    css = ["bootstrap.css","style.css","index.css"]

    lsnew = []
    for l in lines:
        if "const mrparse_html_dir" in l:
            l2 = l[:l.find("=")] + '= "http://127.0.0.1:43434/report_files/0.1.0/mrparse/html";\n'
            lsnew.append(l2)
        elif any(elem in l for elem in svg):
            find = re.compile(r'src="[^"]*"')
            p = "http://127.0.0.1:43434/report_files/0.1.0/mrparse/html/static/" + os.path.basename(find.search(l).group()[4:].strip('"'))
            lsnew.append(find.sub('src="'+p+'"',l))
        elif any(elem in l for elem in js):
            find = re.compile(r'src="[^"]*"')
            p = "http://127.0.0.1:43434/report_files/0.1.0/mrparse/html/" + os.path.basename(find.search(l).group()[4:].strip('"'))
            lsnew.append(find.sub('src="'+p+'"',l))
        elif any(elem in l for elem in css):
            find = re.compile(r'href="[^"]*"')
            p = "http://127.0.0.1:43434/report_files/0.1.0/mrparse/html/" + os.path.basename(find.search(l).group()[5:].strip('"'))
            lsnew.append(find.sub('href="'+p+'"',l))
        else:
            lsnew.append(l)
    return lsnew
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
        #FIXME - consider if PREFERENCES().EXTERNAL_FILES_IN_EXTERNAL_BROWSER before iframing
        if parent is None:
            parent = self
        parent.append("<p>Finished running MrParse</p>")
        basepath = self.jobInfo['fileroot']
        mrparse_rep = os.path.join(basepath, "mrparse_0", 'mrparse.html')

        mrep=open(mrparse_rep, "r")
        mreplines=mrep.readlines()
        mrep.close()
        mrparse_rep_i2_tmp = os.path.join(basepath, "mrparse_0", 'mrparse_i2-tmp.html')
        mrepi2_tmp=open(mrparse_rep_i2_tmp, "w")
        for line in mreplines:
            mrepi2_tmp.write(line.replace("homologs" + os.sep, (os.path.join(basepath, "mrparse_0", 'homologs') + os.sep)).replace("models" + os.sep, (os.path.join(basepath, "mrparse_0", 'models') + os.sep)))
        mrepi2_tmp.close()

        mrparse_rep_i2 = os.path.join(basepath, "mrparse_0", 'mrparse_i2.html')
        with open(mrparse_rep_i2_tmp) as f:
            lines = f.readlines()
            if sys.platform == "win32":
                lsnew = change_i2_mrparse_paths(lines)
            else:
                lsnew = lines
            with open(mrparse_rep_i2 ,"w+") as fout:
                for l in lsnew:
                    fout.write(l)
                fout.write("\n")

        if os.path.exists(mrparse_rep):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)

            mrparseurl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=mrparse_0/mrparse_i2.html?jobNumber="
                + jobNumber
            )
            iframe_style = "display: block;background: #000; margin: 10px; border: none;height: 100vh; width: 95vw;"

        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        if PREFERENCES().EXTERNAL_FILES_IN_IFRAME:
            ResultsI2Folder.append('<br></br>')
            ResultsI2Folder.append('<br></br>')
            ResultsI2Folder.append('<iframe style="{1}" src="{0}"></iframe>'.format(mrparseurl,iframe_style))
        else:
            ResultsI2Folder.append('<span style="font-size:110%">Click on the '
                                    'following link to display the'
                                    'browser report for the MrParse job '
                                    '</span>')
            ResultsI2Folder.append('<br></br>')
            ResultsI2Folder.append('<a href="{0}">Open Results</a>'.format(mrparseurl))
            

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
