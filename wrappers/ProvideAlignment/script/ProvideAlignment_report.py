from report.CCP4ReportParser import *
import sys
import base64

class ProvideAlignment_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ProvideAlignment'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None: parent = self
        
        if self.jobStatus.lower() == 'unsatisfactory':
            parent.addText(style='color:red', text="CCP4i2 was unable to interpret the text you provided as one or more aligned sequences")
        for commentNode in self.xmlnode.findall('.//Commentary'):
            commentFold = parent.addFold(label="Conversion commentary", initiallyOpen=(self.jobStatus.lower() == 'unsatisfactory'))
            commentFold.addText(text='CCP4i2 attempted to interpret the text provided in a range of formats, falling back one to the other.  The following commentary shows the output generated as each format in turn was attempted:')
            commentFold.addPre(text=commentNode.text)
            commentFold.addHelp(ref='$CCP4I2/docs/general/model_data.html',label='Alignment file formats supported by CCP4i2')

        if self.jobStatus.lower() == 'unsatisfactory': return

        parent.addText(text="The alignment seems to be in "+self.xmlnode.findall('.//Format')[-1].text+" format")
        table = parent.addTable(title='Aligned sequences',select='.//Sequence')
        table.addData(title='Index',data=[str(i) for i in range(len(self.xmlnode.findall('.//Sequence')))])
        table.addData(title='Identifier',select='Identifier')
        
        for aliNode in self.xmlnode.findall('.//Alignment'):
            parent.addText(text="Alignment in clustal format:")
            parent.addPre(text=base64.b64decode(aliNode.text).decode())

        for aliNode in self.xmlnode.findall('.//Commentary'):
            commentFold = parent.addFold(label="Conversion commentary", initiallyOpen=False)
            commentFold.addText(select = './/Commentary')
