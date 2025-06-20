import functools
import os
import xml.etree.ElementTree as ET

from ....report import CCP4RvapiParser
from ....report.CCP4ReportParser import Container, Report


dummy_report =  '<html>\n<head>\n<title>Running Crank2</title>\n</head>\n'
dummy_report += '<body>\n<h3>CRANK2 job running - no report available yet</h3>\n</body>\n</html>\n'


def et(fileName=None):
  try:
    root = ET.parse(fileName).getroot()
    script=root.find('body/script')
    if script is not None:
      script.text=script.text.replace('docURI         = "";', 'docURI         = "{}";'.format(os.path.dirname(fileName)+os.sep))
      #root = ET.Element('iframe')
      #root.set('src', fileName)
  except Exception as e:
    print('Crank2 report failed with error message: {0}'.format(e))
    # show logfile 
    if os.path.isfile(os.path.join(os.path.dirname(fileName),'log.txt')):
      with open(os.path.join(os.path.dirname(fileName),'log.txt')) as f:
        t = f.read()
        root = ET.Element("div")
        clearDiv =  ET.SubElement(root,"div")
        clearDiv.set('style','clear:both;')
        pre = ET.SubElement(root,"pre")
        pre.text = t
    else:
    # the code below opens the presentation in the i2 browser.  could be used as an alternative.
      try:
        root = ET.Element('a')
        root.set('href', fileName)
        root.set('style', "font-size: 130%")
        root.text="Presentation not loaded. You can try to click this link to open it or use the View -> Log file  option."
      except Exception as e:
        print('Crank2 report failed (also returning the error message in report): {0}'.format(e))
        root = ET.fromstring(dummy_report)
  return root


class crank2_report(CCP4RvapiParser.RvapiReport):
#class crank2_report(Report):
  TASKNAME="crank2"
  RUNNING = True
  SEPARATEDATA = True
  WATCHED_FILE = 'i2.xml'
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
    rundir = jobInfo['fileroot']
    if os.path.isfile(os.path.join(rundir,'i2_native_present')):
      crank2_report.WATCHED_FILE = self.WATCHED_FILE ='i2.xml'
      CCP4RvapiParser.RvapiReport.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
    else:
      crank2_report.WATCHED_FILE = self.WATCHED_FILE = 'program.xml'
      Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
      #if os.path.isfile(os.path.join(rundir, "index.html")):
      i=0
      while os.path.isdir(os.path.join(rundir, 'job_'+str(i+1))):  i+=1
      if os.path.isfile( os.path.join(rundir, 'job_'+str(i), "index.html") ):
        html = os.path.join(rundir, 'job_'+str(i), "index.html")
      else:
        html = os.path.join(rundir, "index.html")
      results = Container()
      results.as_etree = functools.partial(et,html)
      self.append( results )
