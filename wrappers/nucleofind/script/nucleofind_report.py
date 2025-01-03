
from report.CCP4ReportParser import *
from core import CCP4Utils
import numpy
import os
import shutil

import xml.etree.ElementTree as etree

class nucleofind_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'nucleofind'
  RUNNING = False
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
      clearingDiv = self.addDiv(style="clear:both;")
      self.addDefaultReport(self)
      
  def addDefaultReport(self, parent=None):
      if parent is None: parent=self
      if len(self.xmlnode.findall("LogText")) > 0:
          newFold = parent.addFold(label="Log text", initiallyOpen=True)
          newFold.addPre(text = self.xmlnode.findall("LogText")[0].text)


if __name__ == "__main__":
  import sys
  nucleofind_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
