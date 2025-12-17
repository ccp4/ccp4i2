import xml.etree.ElementTree as etree

from ccp4i2.report import Report


class ProvideSequence_report(Report):
    TASKNAME = 'ProvideSequence'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        
        if self.jobStatus.lower() == 'unsatisfactory':
            parent.addText(style='color:red', text="CCP4i2 was unable to interpret the text you provided as one or more sequence")
        for commentNode in self.xmlnode.findall('./Commentary'):
            commentFold = parent.addFold(label="Conversion commentary", initiallyOpen=(self.jobStatus.lower() == 'unsatisfactory'))
            commentFold.addText(text='CCP4i2 attempted to interpret the text provided in a range of formats, falling back one to the other.  The following commentary shows the output generated as each format in turn was attempted:')
            commentFold.addPre(text=commentNode.text)
            commentFold.addHelp(ref='$CCP4I2/docs/general/model_data.html',label='Sequence file formats supported by CCP4i2')

        print(self.jobStatus.lower())
        if self.jobStatus.lower() == 'unsatisfactory': return

        parent.addText(text="The sequence(s) seemed to be in "+self.xmlnode.findall('.//Format')[-1].text+" format")
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("OK 1")
        table = parent.addTable(title='Aligned sequences',select='./Sequence')
        table.addData(title='Index',data=[str(i) for i in range(len(self.xmlnode.findall('./Sequence')))])
        table.addData(title='Identifier',select='id')
        table.addData(title='Name',select='name')
        table.addData(title='Description',select='description')
        
        print("OK 2")
        parent.addText(text="Sequence(s) in fasta format:")
        print("OK 3")
        print(etree.tostring(self.xmlnode))
        for aliNode in self.xmlnode.findall('./Sequence'):
            print("OK 4?")
            print("Adding ...",aliNode.text)
            parent.addPre(text=aliNode.text)
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
