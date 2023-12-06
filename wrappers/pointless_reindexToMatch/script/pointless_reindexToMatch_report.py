from __future__ import print_function
import os,sys
try:
  from report.CCP4ReportParser import *
except:
  exec(compile(open(os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc'), "rb").read(), os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc'), 'exec'))
  from report.CCP4ReportParser import *
from lxml import etree
from wrappers.pointless.script.pointless_report import pointless_report
#from pointless_report import pointless_report

class pointless_reindexToMatch_report(pointless_report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'pointless_reindexToMatch'
    RUNNING = False
    
    def __init__(self, *args,**kw):
      kw['jobStatus'] = 'nooutput'
      pointless_report.__init__(self,*args, **kw)

      if 'jobStatus' in args and args.get('jobStatus').lower() == 'nooutput': return

      if self.isFatalError():
        errorDiv = self.addDiv(
          style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;color:red")
        self.Errors(errorDiv, colour=False)

      # Get XML element POINTLESS
      if 'xmlnode' in kw:
        xmlnode = kw['xmlnode']
      else:
        if 'xmlFile' in kw:
          xmlnode = etree.parse(kw['xmlFile'])  # for testing
      if xmlnode is None:
        print("**** no xmlnode or xmlFile ****")
        return
      
      analyse = True
      if xmlnode.find('CopyMessage') is not None:
        # usual copy option, not ANALYSE
        analyse = False

      extratext = None
      if analyse:
        extratext = "Analysing merged file with Pointless"

      projectid = None
      jobNumber = None
      if 'jobInfo' in kw:
          jobInfo = kw.get('jobInfo')
          projectid = jobInfo.get('projectid',None)
          jobNumber = jobInfo.get('jobnumber',None)
      self.justPointless(parent=self, extratext=extratext, projectid=projectid, jobNumber=jobNumber)



###########################################################################
if __name__ == "__main__":
  report = pointless_reindexToMatch_report(xmlFile = sys.argv[1],jobStatus="Finished" )
  tree= report.as_etree()
  #  print etree.tostring(tree,pretty_print=True)
  report.as_html_file(fileName='./test-reindex.html')
