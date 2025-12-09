from ccp4i2.report.CCP4ReportParser import *

class import_files_report(Report):
  TASKNAME = 'import_files'
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    kw['standardise'] = False
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
    self.addText(text='The original job has been deleted but the files imported by that job are still available')
    # Add just the Imported files section from the usual 'standard' sections
    self.children.append(ImportedFiles(jobInfo=self.jobInfo))
        
if __name__ == "__main__":
  import sys
  import_files_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
