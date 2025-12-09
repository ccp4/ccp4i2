from __future__ import print_function
from ccp4i2.report.CCP4ReportParser import *
from ccp4i2.report import CCP4RvapiParser
import os,shutil

if sys.version_info >= (3,0):
    from io import StringIO
else:
    from StringIO import StringIO


dummy_report =  '<html>\n<head>\n<title>Running Crank2</title>\n</head>\n'
dummy_report += '<body>\n<h3>CRANK2 job running - no report available yet</h3>\n</body>\n</html>\n'


def et(fileName=None):
  try:
    #root = etree.parse( os.path.join(rundir, "index.html"), parser=parser ).getroot()
    root = etree.parse( fileName ).getroot()
    script=root.find('body/script')
    if script is not None:
      script.text=script.text.replace('docURI         = "";', 'docURI         = "{}";'.format(os.path.dirname(fileName)+os.sep))
      #root = etree.Element('iframe')
      #root.set('src', fileName)
  except Exception as e:
    print('Crank2 report failed with error message: {0}'.format(e))
    # show logfile 
    if os.path.isfile(os.path.join(os.path.dirname(fileName),'log.txt')):
      with open(os.path.join(os.path.dirname(fileName),'log.txt')) as f:
        t = f.read()
        root = etree.Element("div")
        clearDiv =  etree.SubElement(root,"div")
        clearDiv.set('style','clear:both;')
        pre = etree.SubElement(root,"pre")
        pre.text = t
    else:
    # the code below opens the presentation in the i2 browser.  could be used as an alternative.
      try:
        root = etree.Element('a')
        root.set('href', fileName)
        root.set('style', "font-size: 130%")
        root.text="Presentation not loaded. You can try to click this link to open it or use the View -> Log file  option."
      except Exception as e:
        print('Crank2 report failed (also returning the error message in report): {0}'.format(e))
        f = StringIO(dummy_report)
        root = etree.parse( f ).getroot()
        f.close()
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
      import functools
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
