import os
import sys
#from lxml import etree
from xml.etree import ElementTree as ET

try:
  from report.CCP4ReportParser import *
except:
  exec(compile(open(os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc')).read(), os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc'), 'exec'))
  from report.CCP4ReportParser import *

# - - - - - - - - - - - - - - - - -
class freerflag_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'freerflag'
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
      if self.errorReport().maxSeverity()>SEVERITY_WARNING:
        print('FAILED instantiating FreeRflag report generator')
        self.errorReport().report()
        return

      ##if jobStatus is None or jobStatus.lower() is 'nooutput': return
        
      # 'nooutput' mode would be used by another report class that wanted
      # to use some method(s) from this class for its own report
      if jobStatus is not None and jobStatus.lower() == 'nooutput':
        return

      results = self.addResults()  # Does this do anything?
      # for stand-alone FreeR job (not in pipeline)
      self.justFreeRflag(self)

    # - - - - - - - - - - - - - - - - -
    def justFreeRflag(self, parent=None):
      # for stand-alone FreeR job (not in pipeline)
      text = 'GENERATE FREER SET'
      parent.addText(text=text,
                     style='font-weight:bold; font-size:120%; color:blue;')
      parent.append(' <br/>')
      self.mainReport(parent, self.xmlnode)

    # - - - - - - - - - - - - - - - - -
    def mainReport(self, parent, freerxml):

      elementTag = freerxml.tag
      mode = ''
      if elementTag != 'FREERFLAGINFO':
        pass # do something for pipeline?

      
      if len(freerxml.findall('Mode'))>0:
        mode = freerxml.findall('Mode')[0].text
        parent.append(' <br/>')
        if mode == 'Complete':
          message = 'FreeR flags generated to COMPLETE input set'
        elif mode == 'New':
          fraction = freerxml.findall('Fraction')[0].text
          message = 'New FreeR flags generated with fraction ' + fraction
        else:
          message = 'WARNING: unknown FreeR mode'
        parent.addText(text=message)

        if mode == 'Complete':
          highResD = freerxml.findall('ObservedDataResolution')[0].text
          highResF = freerxml.findall('FreeR_Resolution')[0].text
          cutfreerresolution = (freerxml.findall('CutFreerResolution')[0].text == 'True')
          if cutfreerresolution:
            frcut = freerxml.findall('FreerCutResolution')[0].text
            message = 'Input FreeR data was cut to match the resolution'+\
                      ' of the data, '+frcut+' \xc5'
            parent.append(' <br/>')
            parent.addText(text=message)

      if len(freerxml.findall('GlobalResolutionLimit'))>0:
        glres = freerxml.findall('GlobalResolutionLimit')[0].text
        message = 'For output, the FreeR set was cut to  resolution limit of '+\
                  glres+' \xc5'
        parent.append(' <br/>')
        parent.addText(text=message)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
if __name__ == "__main__":
    report = freerflag_report(xmlFile = sys.argv[1] )
    tree= report.as_etree()
    #  print ET.tostring(tree)
    report.as_html_file(fileName='./test-freer.html')

