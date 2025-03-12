from __future__ import print_function


import os,sys
from lxml import etree

def getCCP4I2Dir(up=1):
    target = os.path.join(os.path.realpath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    splittarget = abstarget.split()
    if splittarget.count('ccp4i2'):
        splittarget.reverse()
        up = splittarget.index('ccp4i2')
    while up>0:
        abstarget = os.path.dirname(abstarget)
        up = up -1
    return abstarget

if __name__ == '__main__':

    top_path = getCCP4I2Dir()
    print('Running CCP4i2 makeReport from: '+top_path)
    exec(compile(open(os.path.join(top_path,'utils','startup.py')).read(), os.path.join(top_path,'utils','startup.py'), 'exec'))
    setupEnvironment()
    setupPythonpath()
    from report import CCP4ReportGenerator,CCP4ReportParser

    argList = sys.argv[1:]
    print('argList',argList,len(argList))
    if len(argList)<2:
        print('Usage: makeReport xrtFileName xmlDataFileName')
        sys.exit()

    xrt = etree.XML( open( argList[0] ).read() )
    xml = etree.XML( open( argList[1] ).read() )
    jobId = None
    outputFile = './report.html'
    ifPrint = False
    ii = 2
    while ii < len(argList):
        print('arg',ii,argList[ii])
        if argList[ii][0:2] == '-p':
            ifPrint = True
        elif argList[ii][0:2] == '-j':
            ii = ii + 1
            jobId = argList[ii]
        elif argList[ii][0:2] == '-o':
            ii = ii + 1
            outputFile = argList[ii]
            ii = ii + 1
    #print argList,outputFile,jobId,ifPrint
    if jobId is not None:
        exec(compile(open(os.path.join(top_path,'utils','startup.py')).read(), os.path.join(top_path,'utils','startup.py'), 'exec'))
        pm = startProjectsManager()

        rg = CCP4ReportGenerator.CReportGenerator(jobId)
        jobInfo = rg.getJobInfo(jobId)
    else:
        jobInfo = {}

    xreport = xrt.xpath( "/report" )[0]
    report = CCP4ReportParser.Report( xreport, xml,jobInfo=jobInfo )
    text = report.as_html()
    if ifPrint: print(text)

    f = open(outputFile,'w')
    f.write(text)
    f.close()
