from __future__ import print_function


#==============================================================================================

import os
import sys
import functools

#from lxml import etree
import xml.etree.ElementTree as etree

from PySide6 import QtCore
from core.CCP4ErrorHandling import *
from dbapi import CCP4DbApi
from core.CCP4Modules import PREFERENCES
from core import CCP4Utils

class CReportGenerator(QtCore.QObject):

  FinishedPictures = QtCore.Signal(str)

  ERROR_CODES = { 1 : { 'description' : 'Report definition file not found',
                        'severity' : SEVERITY_WARNING },
                  2 : { 'description' : 'Task data file not found',
                        'severity' : SEVERITY_WARNING },
                  3 : { 'description' : 'Error loading report definition file' },
                  4 : { 'description' : 'Error loading task data file' },
                  5 : { 'description' : 'Failed to find insert point for sub-job report in parent jobs report file' },
                  6 : { 'description' : 'Error inserting report on sub-job into report file' },
                  7 : { 'description' : 'Error saving report file with inserted report on sub-job' },
                  8 : { 'description' : 'Error reading program data file' },
                  9 : { 'description' : 'Error no program or task data file' },
                 10 : { 'description' : 'No report definition file available' }
               }
  insts = []
  def __init__(self,jobId,jobStatus='Finished',**kw):
    QtCore.QObject.__init__(self)
    from dbapi import CCP4DbApi
    self.jobId = CCP4DbApi.UUIDTYPE(jobId)
    if jobStatus=='To delete':
      self.jobStatus = 'Finished'
    else:
      self.jobStatus = jobStatus
    self.jobNumber = kw.get('jobNumber',None)
    self.reportFile = None

    # Defaults for mg picture generation
    self.centre_xyz = [ 0.0,0.0,0.0]
    self.scale = 100.0
    self.orientation = [ 1.0,0.0,0.0,0.0]
    self.imageWidth = 400
    self.imageHeight = 400
    self.jobInfo = None

  def setJobStatus(self,status):
    # Force resetting data from db
    self.jobStatus = status
    if self.jobStatus in CCP4DbApi.FINISHED_JOB_STATUS:
      self.jobInfo = None
      
  def getReportClass(self,doReload=False):
    from core import CCP4TaskManager,CCP4Modules
    #print 'getReportClass',self.jobId
    taskName = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId,mode='taskname')
    cls = CCP4TaskManager.TASKMANAGER().getReportClass(taskName,jobStatus=self.jobStatus,doReload=doReload)
    xrtFile = CCP4TaskManager.TASKMANAGER().searchXrtFile(taskName,jobStatus=self.jobStatus)
    return taskName,cls,xrtFile

  def getCustomFailedReport(self):
    from core import CCP4Modules, CCP4File, CCP4PluginScript
    reportName = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='DIAGNOSTIC')
    #print 'getCustomFailedReport',reportName
    if not os.path.exists(reportName): return None
    try:
      i2Xml = CCP4File.CI2XmlDataFile(fullPath=reportName)
      errReport = CErrorReport()
      errReport.setEtree(i2Xml.getBodyEtree())
    except:
      print('Error loading diagnostic for creating report',reportName)
      return None
    
    if errReport.count(CCP4PluginScript.CPluginScript,36):
      return CFreeRErrorReport
    else:
      return None

  def htmlBase(self):
    #MN this will avoid starting up a server thread if one is not running
    from qtcore.CCP4HTTPServerThread import CHTTPServerThread
    port=43434
    if CHTTPServerThread.insts is not None:
        port = CHTTPServerThread.insts.port
    #print 'makeReportFile HTTPSERVER',CCP4Modules.HTTPSERVER()
    #print 'makeReportFile port',port
    return 'http://127.0.0.1:'+str(port)+'/report_files'

  def getOutputXml(self):
    from core import CCP4Modules
    from report import CCP4ReportParser
    import os
    outputXml = None
    outputXmlFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='PROGRAMXML')
    #print  'CReportGenerator.makeReportFile outputXmlFile',outputXmlFile,os.path.exists(outputXmlFile)
    #MN: Try to make this all a bit more atomic since programxl file might come and go
    # Liz improved error trapping - Nov 15
    # Jon added support for rvapi's i2.xml - previously, it only generated report if program.xml was present
    # Thanks Pep Trivinyo (Uson group, Barcelona) for pointing this out
    
    #print('getOutputXml',outputXmlFile)
    try:
        with open(outputXmlFile,'r') as inputFile:
            text = inputFile.read()
            if sys.version_info > (3,0):
                outputXml = etree.fromstring( text.encode('utf-8'))
            else:
                outputXml = etree.fromstring( text,CCP4ReportParser.PARSER() )
    except Exception as e:
        #print('Error in getOutputXml',e)
        outputXmlFile =  os.path.join(os.path.split(outputXmlFile)[0],'XMLOUT.xml')
        try:
          with open(outputXmlFile,'r') as inputFile:
              text = inputFile.read()
              outputXml = etree.fromstring( text,CCP4ReportParser.PARSER() )
        except Exception as e :
            outputXmlFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='RVAPIXML')
            try:
             with open(outputXmlFile,'r') as inputFile:
                 text = inputFile.read()
                 outputXml = etree.fromstring( text, CCP4RvapiParser.PARSER() )
            except:
                # no recognised form of XML data has been found
                outputXmlFile=None
                outputXml = None

    return outputXmlFile,outputXml
        
  def makeReportFile(self,redo=False,doReload=False,useGeneric=False,makePictures=True):
      
      #MN changed so as to return a tuple reportFile (a path) and newPageOrNewData: the latter has value "NEWPAGE" if a brand new report page is created and needs to be uploaded, or "NEWDATA" if an existing report.html needs only to be updated with new data for graphs, tables etc.
    from core import CCP4Modules
    import os
    from core import CCP4ErrorHandling
    from report import CCP4ReportParser
    
    newPageOrNewData = "NEWPAGE"
   
    taskName,cls,xrtFile = self.getReportClass(doReload=doReload)
    #print 'makeReportFile',taskName,cls,xrtFile,useGeneric
    if cls is None and xrtFile is None:
      if useGeneric:
        cls = CCP4ReportParser.GenericReport
      else:
        raise CCP4ErrorHandling.CException(self.__class__,10,taskName,stack=False)  

    self.reportFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='REPORT')
    if self.jobStatus in [ 'Running', 'Running remotely'] :
      # Always redo a report for a running job
      fSplit = os.path.splitext(self.reportFile)
      self.reportFile = fSplit[0]+'_tmp'+fSplit[1]
      #print 'makeReportFile',self.jobStatus,self.reportFile,os.path.exists(self.reportFile),redo
    else:
      if os.path.exists(self.reportFile) and not redo: return self.reportFile, newPageOrNewData
  
    if cls.USEPROGRAMXML:
      outputXmlFile,outputXml= self.getOutputXml()
      if outputXml is None:
        
        if  self.jobStatus not in [ 'Running', 'Pending' , 'Running remotely' ]:
          return self.makePleaseWaitReport(jobRunning=False), newPageOrNewData
        else:
          return self.makePleaseWaitReport(), newPageOrNewData
        
        #print 'makeReportFile outputXml',outputXml
    else:
        outputXml = etree.fromstring('<dummy/>')
        
    if self.jobInfo is None:
      self.jobInfo = getReportJobInfo(self.jobId)


    projectId = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId,mode='projectid')
    if cls is not None:
      self.report = cls(xmlnode=outputXml,jobInfo=self.jobInfo, standardise=( self.jobStatus not in [ 'Running','Running remotely']), jobStatus = self.jobStatus , jobNumber=self.jobNumber, projectId =  projectId)
      #print 'CReportGenerator report',report
      self.report.as_html_file(fileName=self.reportFile,htmlBase=self.htmlBase())
      #rtfFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='TABLE_RTF')
      #report.graph_data_as_rtf(fileName=rtfFile)
    elif useGeneric:
      self.report = CCP4ReportParser.GenericReport(xmlnode=outputXml,jobInfo=self.jobInfo, standardise=( self.jobStatus not in [ 'Running','Running remotely']), jobStatus = self.jobStatus , jobNumber=self.jobNumber, taskName = taskName, projectId =  projectId )
      self.report.as_html_file(fileName=self.reportFile,htmlBase=self.htmlBase())
    else:
      try:
        xrtTree = etree.XML( open(xrtFile ).read() )
      except Exception as e:
        print(e)
        raise CCP4ErrorHandling.CException(self.__class__,3,stack=False)  
      # standardise  report (add i/o/ files etc) if self.jobStatus is not Running
      self.report = CCP4ReportParser.Report( xrtTree.xpath( "/report" )[0], outputXml,jobInfo=self.jobInfo, standardise=( self.jobStatus not in [ 'Running','Running remotely']) , jobNumber=self.jobNumber, taskVersion=self.jobInfo.get('taskversion',None) )
      self.report.as_html_file(fileName=self.reportFile,htmlBase=self.htmlBase())
      
    #print 'makeReportFile',self.report.containsPictures(),self.report.pictureQueue
    if self.jobStatus not in [ 'Running','Running remotely'] :
      self.tableDirectory = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='TABLES_DIR')
      if not os.path.exists(self.tableDirectory): os.mkdir(self.tableDirectory)
      self.report.makeCSVFiles(directory=self.tableDirectory)

    jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=self.jobId)
    xmlTableDirectory = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='XML_TABLES_DIR')
    if not os.path.exists(xmlTableDirectory): os.mkdir(xmlTableDirectory)
    
    #Here enquire of the report object whether it thinks a new page or a refrsh of the
    #existing page from new data is called for
    #NB I do this before creating the corresponding XML files, since they
    #may hold the key to whether a page load is needed
    
    newPageOrNewData = self.report.newPageOrNewData()
    
    self.report.makeXMLFiles(directory=jobDirectory)

    if not makePictures: return self.reportFile, newPageOrNewData
    
    """
    if self.report.containsPictures():
      from core import CCP4Config
      try:
        self.mgProcess = self.runMg(pictureQueue=self.report.pictureQueue ,callBack=lambda exitCode,exitStatus: self.handleMgFinished(self.jobId,self.report.pictureQueue[1:],exitCode,exitStatus))
      except:
        print('ERROR making report, running CCP4mg to create pictures')
      else:
        pass
    """
    return self.reportFile, newPageOrNewData

  def makePleaseWaitReport(self,jobRunning=True):
    from report import CCP4ReportParser
    if self.jobInfo is None:
      self.jobInfo = getReportJobInfo(self.jobId)
    doc = CCP4ReportParser.htmlDoc()
    body = doc.getroot().findall('./body')[0]
    h3 = etree.Element('h3')
    if jobRunning:
      h3.text = 'Please wait, no output yet from job '+str(self.jobInfo['jobnumber'])+': '+ self.jobInfo['tasktitle']
    else:
      h3.text = 'No output from job '+str(self.jobInfo['jobnumber'])+': '+ self.jobInfo['tasktitle'] + ' - can not make report'
    body.append(h3)

    # Write to file
    text = etree.tostring(doc.getroot())
    CCP4Utils.saveFile(fileName=self.reportFile,text=text)
    return self.reportFile

  def makeFailedReportFile(self,redo=False):
    import glob
    from report import CCP4ReportParser
    from core import CCP4Modules
    if self.jobStatus != 'Failed':
      self.reportFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='DIAGNOSTIC_REPORT')
      title = 'Diagnostic'
    else:
      title = 'Error Report'
      self.reportFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='REPORT')

    reportClass = self.getCustomFailedReport()
    #print 'makeFailedReportFile getCustomFailedReport',reportClass
    if reportClass is not None:
      if self.jobInfo is None:
        self.jobInfo = getReportJobInfo(self.jobId)
      outputXml = etree.Element('dummy')
      report = reportClass(xmlnode=outputXml,jobInfo=self.jobInfo, standardise=True, jobStatus = self.jobStatus , jobNumber=self.jobNumber )
      #print 'CReportGenerator report',report
      report.as_html_file(fileName=self.reportFile,htmlBase=self.htmlBase())
      return self.reportFile
    
    if os.path.exists(self.reportFile) and not redo: return self.reportFile

    if self.jobInfo is None:
      self.jobInfo = getReportJobInfo(self.jobId)
    self.jobInfo['jobStatus'] = 'Failed'
    # getReportClass will only return the task class if it has class variable FAILED == True
    tName,cls,xrtFile = self.getReportClass()
    report = None
    if cls is not None:
      try:
        report = cls(jobInfo = self.jobInfo,jobStatus='Failed')
      except:
        pass
    if report is None:
      report = CCP4ReportParser.Report(jobInfo = self.jobInfo)
    # If we have some content from the task class then need to insert the title line before that content
    report.insert(0,'<h3 class="error_report">'+title+' for Job '+str(self.jobInfo['jobnumber'])+': '+ self.jobInfo['tasktitle']+'</h3>')
    
    # If there diagnostic output - display it
    reportFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='DIAGNOSTIC')
    if os.path.exists(reportFile):
      t = report.addPre(style='color:red')
      t.text = self.getErrorReport(reportFile)

    logFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='LOG')
    if os.path.exists(logFile):
      fold = CCP4ReportParser.Fold()
      fold.label =  'Log file'
      try:
        if os.path.exists(logFile[0:-7] + 'com.txt' ):
            t = fold.addPre()
            t.text =  CCP4Utils.readFile(logFile[0:-7] + 'com.txt')
        t = fold.addPre()
        t.text =  CCP4Utils.readFile(logFile)
      except Exception as e:
        print('Error attempting to read log file for error report:',logFile)
        print(e)
      report.children.append(fold)
    

    for std,title,target in (('STDERR','Error output from job','errors'),
                      ('STDOUT','Terminal output from job','terminal'),
                      ('DIAGNOSTIC','Diagnostic from script','diagnostic')):
      filePath = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode=std)
      #print 'makeFailedReportFile',fileName,os.path.exists(fileName)
      if os.path.exists(filePath):
        fold = CCP4ReportParser.Fold()       
        fold.label =  os.path.split(filePath)[1]        
        if filePath[-7:] == 'log.txt' and os.path.exists(filePath[0:-7] + 'com.txt' ):
          t = fold.addPre()
          t.text =  CCP4Utils.readFile(filePath[0:-7] + 'com.txt')
        t = fold.addPre()
        t.text =  CCP4Utils.readFile(filePath,limit=PREFERENCES().TEXT_VIEW_LINE_LIMIT*10)
        report.children.append(fold)


    top = os.path.split(filePath)[0]
    topFileList = glob.glob(os.path.join(top,'*.log'))
    topFileList.extend(glob.glob(os.path.join(top,'log*.txt')))
    fullFileList = self.getSubJobLogs(top)
    for n in range(len(topFileList)-1,-1,-1):
      if topFileList[n] != logFile:
        fullFileList.insert(0,[topFileList[n],0])

    for filePath,time in fullFileList:                      
      try:
        fold = CCP4ReportParser.Fold()
        subJobNum,subJobId = self.extractJobId(filePath)
        if subJobId is not None:
          taskName = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=subJobId,mode='taskname')
        else:
          taskName = ''
        fold.label =  filePath.replace(top,'')[1:]  + ' from ' + subJobNum + ' ' +  taskName  
        if filePath[-7:] == 'log.txt' and os.path.exists(filePath[0:-7] + 'com.txt' ):
          t = fold.addPre()
          t.text =  CCP4Utils.readFile(filePath[0:-7] + 'com.txt')
        t = fold.addPre()
        t.text =  CCP4Utils.readFile(filePath,limit=PREFERENCES().TEXT_VIEW_LINE_LIMIT*10)
        #print 'reading',filePath,t.text
      except Exception as e:
        print('Error attempting to read log file for error report:',filePath)
        print(e)
      else:
        report.children.append(fold)

    # Write to file  
    report.as_html_file(fileName=self.reportFile)
    return self.reportFile

  def getSubJobLogs(self,top):
    # Need to get sub-job log files!
    fullFileList = []
    for path,dirList,fileList in os.walk(top):
      for fileName in fileList:
        if fileName[0:3] == 'log' or os.path.splitext(fileName)[1] == '.log':
            try:
              filePath = os.path.join(path,fileName)
              fullFileList.append([filePath,os.path.getmtime(filePath)])
            except:
              pass

    def sortFiles( a , b ):
      if a[1] > b[1] :
        return 1
      else:
        return -1
    if sys.version_info > (3,0):
        fullFileList = sorted(fullFileList, key=functools.cmp_to_key(sortFiles) )
    else:
        fullFileList = sorted(fullFileList, cmp = sortFiles )
    return fullFileList


  def extractJobId(self,path):
    from core import CCP4Modules
    num = ''
    continu = True
    while continu:
      try:
        path,ele = os.path.split(path)
        if ele[0:4]=='job_': num = ele[4:]+'.'+num
      except:
        continu = False
      else:
        if len(ele)==0: continu = False
    if len(num)==0:
      return None
    else:
      projectId = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId,mode='projectid')
      subJobId = CCP4Modules.PROJECTSMANAGER().db().getJobId(projectId=projectId,jobNumber=num[0:-1])
      return num[0:-1],subJobId

      
  def getErrorReport(self,fileName):
    from core import CCP4File
    f = CCP4File.CI2XmlDataFile(fullPath=fileName)
    body = f.getEtreeRoot().find('ccp4i2_body')
    report = CErrorReport()
    report.setEtree(body)
    return report.report(user=True,ifStack=False)


  def mergeIntoParent(self,parentFile=None):
    #parentTree = etree.parse(parentFile)
    parentTree = CCP4Utils.openFileToEtree(parentFile)
    #print 'mergeIntoParent jobId',self.jobId
    # Find the span-link created by CCP4ReportParser.foldLinkLine
    path = 'body/div[@class="sub-job-list"]/span[@class="folder_link"]/a[@id="jobId'+str(self.jobId)+'"]'
    aEle = parentTree.xpath(path)
    #print 'mergeIntoParent aEle',aEle
    if len(aEle) == 0: return CErrorReport(self.__class__,5,parentFile,stack=False)
    aEle = aEle[0]
    label = aEle.text
    #print 'mergeIntoParent',label
    insertEle = aEle.xpath('../..')[0]
    insertIndex = insertEle.index(aEle.xpath('..')[0])

    # Convert body of sub-job report to a hidesection div
    #myTree =  etree.parse(self.reportFile)
    myTree = CCP4Utils.openFileToEtree(self.reportFile)
    body = myTree.xpath('./body')[0]
    body.tag ='div'
    body.set('id','jobId'+str(self.jobId))
    body.set('class','subjob')

    # make a span for the folder title
    # This should do the same as CCP4ReportParser.foldTitleLine()
    span = etree.Element('span')
    span.set('class','folder')
    span.set('onclick','toggleview(this)')
    span.text = label

    # Swap out the span-link and swap in proper folder
    insertEle.remove(aEle.xpath('..')[0])
    insertEle.insert(insertIndex,span)
    insertEle.insert(insertIndex+1,body)

    newFile = os.path.join(os.path.split(parentFile)[0],'tmp_report.html')
    #print 'mergeIntoParent newFile',newFile
    text = etree.tostring(parentTree)
    try:
      CCP4Utils.saveFile(fileName=newFile,text=text)
    except:
      raise CException(self.__class__,7,newFile,stack=False)

    return CErrorReport()

  def runMg(self,pictureQueue=[],callBack=None):
    #print 'CReportGenerator.runMg pictureQueue',pictureQueue
    from core import CCP4Modules,CCP4Config
    import os,re
    mgExe = CCP4Modules.LAUNCHER().getExecutable('CCP4MG')
    if mgExe is None:
      print('ERROR no CCP4mg executable found')
      raise CException(self.__class__,103,stack=False) 

    picDefFile = os.path.normpath(pictureQueue[0])
    options = { "size":str(self.imageWidth)+'x'+str(self.imageHeight) }
    
    argList = []
    for item in ['-nore','-quit']: argList.append(item)
    argList.append(picDefFile)
    argList.append('-R')
    picFile = os.path.normpath(os.path.join(os.path.splitext(os.path.splitext(picDefFile)[0])[0]+'.png'))
    argList.append(picFile)
    argList.append('-RO')
    argList.append(re.sub(': ',':',str(options)))

    from core import CCP4Modules
    process = CCP4Modules.LAUNCHER().launch(viewer='ccp4mg',argList=argList,callBack=callBack)
    return process
              

  @QtCore.Slot(str,list,int,int)
  def handleMgFinished(self,jobId,pictureQueue,exitCode,exitStatus):
    #print 'handleMgFinished',jobId,pictureQueue
    if len(pictureQueue)>0:
      self.mgProcess = self.runMg(pictureQueue,callBack=lambda exitCode2,exitStatus2: self.handleMgFinished(jobId,pictureQueue[1:],exitCode2,exitStatus2))
    else:
      self.mgProcess = None
          
    self.FinishedPictures.emit(jobId)

#======================================================================================================
# NB is function so can work stand-alone
def getReportJobInfo(jobId=None,projectName=None,jobNumber=None):

    """
    print("getReportJobInfo")
    print(jobId,projectName,jobNumber)
    import traceback
    traceback.print_stack()
    """
  
    from core import CCP4Modules,CCP4TaskManager,CCP4Data
    from dbapi import CCP4DbApi
    db = CCP4Modules.PROJECTSMANAGER().db()


    if projectName is None:
        projectNameInfo = db.getJobInfo(jobId=jobId)
        projectId = projectNameInfo['projectid']
        projectName = db.getProjectInfo(projectId,mode='projectname')

    if jobId is None: jobId = db.getJobId(projectName=projectName,jobNumber=jobNumber)

    jobInfo = db.getJobInfo(jobId=jobId,mode=['runtime','status','taskname','taskversion','jobnumber','descendentjobs','projectid','projectname',
                                              'jobtitle','creationtime'])
    # print 'CReportGenerator getJobInfo',jobInfo
    jobInfo['jobid'] = jobId
    jobInfo['tasktitle'] = CCP4TaskManager.TASKMANAGER().getTitle(jobInfo['taskname'])
    jobInfo['fileroot'] = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jobId,mode='ROOT')
    importedfiles = db.getJobImportFiles(jobId=jobId)
    # Returned list JobId,FileID,ImportId,FileTypeId,Filename,Annotation
    jobInfo['importedfiles'] = []
    for imp in importedfiles:
      jobInfo['importedfiles'].append( { 'fileId': imp[1], 'filetypeid':imp[3]  ,'filetype': '', 'filetypeclass' : CCP4DbApi.FILETYPES_CLASS[imp[3]], 'filename': imp[4],'relpath':'' ,'projectname': projectName,'projectid': projectId ,'annotation' : imp[5]} )
    
    #print 'getReportJobInfo importedfiles', jobInfo['importedfiles']

    #jobInfo['childjobinfo'] = {}
    #for jid,childList in jobInfo['descendentjobs']:
      

    for key,role in [['inputfiles',CCP4DbApi.FILE_ROLE_IN],['outputfiles',CCP4DbApi.FILE_ROLE_OUT]]:
      jobInfo[key] = []
      fileIdList = db.getJobFiles(jobId=jobId,role=role)
      for fileId in fileIdList:
        fileInfo = db.getFileInfo(fileId = fileId,mode=['filetypeid','filetype','filename','relpath','projectname','projectid','annotation','jobparamname'])
        fileInfo['filetypeclass'] = CCP4DbApi.FILETYPES_CLASS[fileInfo['filetypeid']]
        fileInfo['fileId'] = fileId
        jobInfo[key].append(fileInfo)
    jobInfo['filenames'] = {}
    """
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print('CReportGenerator.getJobInfo',jobInfo['importedfiles'])
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    """
    try:
      container = db.getParamsContainer(jobId=jobId)
    except:
      pass
    else:
      # Get the data keyed by task parameter name- not necessarilly all files but should all have __str__ methods
      for key in container.inputData.dataOrder():
        if isinstance(container.inputData.find(key),CCP4Data.CList):
          jobInfo['filenames'][key] = []
          for item in container.inputData.find(key):
            jobInfo['filenames'][key].append(item.__str__())
        else:
          jobInfo['filenames'][key] = container.inputData.find(key).__str__()
      for key in container.outputData.dataOrder():
        if isinstance(container.outputData.find(key),CCP4Data.CList):
          jobInfo['filenames'][key] = []
          for item in container.outputData.find(key):
            jobInfo['filenames'][key].append(item.__str__())
        else:
          try:
              jobInfo['filenames'][key] = container.outputData.find(key).__str__()
          except:
              jobInfo['filenames'][key] = ""
    #print 'CReportGenerator filenames',jobInfo['filenames']
    
    return jobInfo


from report.CCP4ReportParser import Report

class CFreeRErrorReport(Report):

  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
    self.addText(text="Your chosen Free-R flag data do not cover all of the experimental data, for example due to insufficient resolution. Choose a different Free-R object, or use the 'Generate a Free-R set' task (X-ray data reduction and analysis) to complete your Free-R flag data.")

    
