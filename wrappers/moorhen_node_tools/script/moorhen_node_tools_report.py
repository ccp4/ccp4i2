from report.CCP4ReportParser import *
import sys

class moorhen_node_tools_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'moorhen_node_tools'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        self.addText(text='Happily finished')
        try:
            pdbsWrittenPath = './/number_output_pdbs'
            if len(xmlnode.findall(pdbsWrittenPath)) > 0:
                pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
                self.addText(text='Number of PDBs written: ' + pdbsWrittenString)
            cifsWrittenPath = './/number_output_cifs'
            if len(xmlnode.findall(cifsWrittenPath)) > 0:
                cifsWrittenString = xmlnode.findall(cifsWrittenPath)[0].text
                self.addText(text='Number of CIFs written: ' + cifsWrittenString)
        except:
            pass


