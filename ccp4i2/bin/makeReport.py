import sys

from lxml import etree

from .. import I2_TOP
from ..utils.startup import setupEnvironment, startProjectsManager


def main():
    print('Running CCP4i2 makeReport from:', I2_TOP)
    setupEnvironment()
    from ..report import CCP4ReportGenerator, CCP4ReportParser

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
    jobInfo = {}
    if jobId is not None:
        startProjectsManager()
        rg = CCP4ReportGenerator.CReportGenerator(jobId)
        jobInfo = rg.getJobInfo(jobId)

    xreport = xrt.xpath( "/report" )[0]
    report = CCP4ReportParser.Report( xreport, xml,jobInfo=jobInfo )
    text = report.as_html()
    if ifPrint: print(text)

    f = open(outputFile,'w')
    f.write(text)
    f.close()
