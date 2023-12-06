from __future__ import print_function

from report.CCP4ReportParser import *
from core import CCP4ModelData, CCP4Modules
from report import CCP4ReportParser
import sys,os
import xml.etree.ElementTree as etree

class ProvideAsuContents_report(Report):
    TASKNAME = 'ProvideAsuContents'
    RUNNING = False
    USEPROGRAMXML = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)

      #Default file name
      fileName = os.path.join( jobInfo['fileroot'], 'ASUCONTENTFILE.asu.xml')

      #File name from job info
      if hasattr(self,"jobInfo") and self.jobInfo and "outputfiles" in self.jobInfo:
          for outfile in self.jobInfo["outputfiles"]:
              if outfile['jobparamname'] == "ASUCONTENTFILE":
                  fileName = os.path.join( jobInfo['fileroot'], outfile['filename'])

      fileObj = CCP4ModelData.CAsuDataFile()
      fileObj.setFullPath(fileName)
      summaryDiv = self.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      fold = summaryDiv.addFold(label="AU content - sequences", initiallyOpen=True)
      for seqObj in fileObj.fileContent.seqList:
        div = fold.addDiv()
        #print 'ProvideAsuContents report',seqObj.name
        text = '<p>' + str(seqObj.nCopies)+' copies of '+str(seqObj.name) + '<br/>'
        text += 'Description: '+str(seqObj.description) + '<br/>'
        if seqObj.source.isSet(): text += 'Loaded from file: '+str(seqObj.source)  + '<b/>'
        text += '</p>'
        #print 'ProvideAsuContents report', text
        div.append(text)
        '''
        div.addText(text=str(seqObj.nCopies)+' copies of '+str(seqObj.name))
        div.addText(text='Description: '+str(seqObj.description))
        if seqObj.source.isSet(): div.addText(text='Loaded from file: '+str(seqObj.source))
        '''
        div.addPre(text=seqObj.formattedSequence())

#Oh dear, I want USEPROGRAMXML = True, but that breaks reports for old jobs, so I have to reimplement getOutputXml here so that
#I can support old jobs without program.xml and new ones with.
      if hasattr(self,"jobInfo") and self.jobInfo and "jobid" in self.jobInfo:

          outputXmlFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobInfo["jobid"],mode='PROGRAMXML')

          if os.path.exists(outputXmlFile):
              #A new job
              with open(outputXmlFile,'r') as inputFile:
                  text = inputFile.read()
                  if sys.version_info > (3,0):
                      outputXml = etree.fromstring( text.encode('utf-8'))
                  else:
                      outputXml = etree.fromstring( text )
                  matthewsAnalysis = outputXml.findall( ".//matthewsCompositions" )
                  if len(matthewsAnalysis) > 0 and len(matthewsAnalysis[0].getchildren())>0:
                      text = ''
                      fold = summaryDiv.addFold(label="Analysis", initiallyOpen=True)
                      div = fold.addDiv()
                      matthewsComposition = outputXml.findall( ".//matthewsCompositions/composition" )
                      if len(matthewsComposition)> 0:
                          title = "Matthews analysis"
                          fullTable = fold.addTable(select=".//matthewsCompositions/composition",xmlnode=outputXml)
                          for title,select in  [[ "Number of copies" ,"nMolecules" ], [ "Solvent %"  , "solventPercentage" ],[ "Matthews coefficient" , "matthewsCoeff" ],["Matthews probability","matthewsProbability"]]:
                              fullTable.addData(title=title,select=select)

                      entries = outputXml.findall( ".//entries/entry" )
                      if len(entries)>0:
                          title = "Sequence inputs to Matthews analysis"
                          fullTable = fold.addTable(select=".//entries/entry",xmlnode=outputXml)
                          for title,select in  [[ "Name"  , "name" ],[ "Number of copies" ,"copies" ],[ "Molecular weight" , "weight" ]]:
                              fullTable.addData(title=title,select=select)

                      div.append(text)

                      totalWeight = outputXml.findall( ".//totalWeight" )
                      cellVolume = outputXml.findall( ".//cellVolume" )
                      if len(totalWeight) > 0:
                          div.append('<p>Total molecular weight: {0:.1f}</p>'.format(float(totalWeight[0].text)))
                      else:
                          div.append('<p>Total molecular weight unknown</p>')
                      if len(cellVolume) > 0:
                          div.append('<p>Cell volume: {0:.2f}</p>'.format(float(cellVolume[0].text)))
                      else:
                          div.append('<p>Cell volume unknown</p>')
