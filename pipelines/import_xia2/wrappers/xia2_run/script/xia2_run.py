from __future__ import print_function
"""
     xia2_run.py: CCP4 GUI Project
     Copyright (C) 2013 STFC
"""

import os,shutil,glob
from lxml import etree
from core import CCP4PluginScript
import base64

RUN_TITLES = {
    '2d' : 'Mosflm-Scala integration and processing',
    '3d' : 'XDS-XSCALE integration and processing',
    '3dii' : 'XDS-XSCALE integration with image indexing and processing',
    '2da' : 'Mosflm-Aimless integration and processing',
    '3da' : 'XDS-XSCALE-Aimless integration and processing',
    '3daii' : 'XDS-Aimless integration with image indexing and processing'
    }

                                      
class xia2_run(CCP4PluginScript.CPluginScript):

    TASKMODULE = None      # Where this plugin will appear on the gui
    TASKTITLE = 'Import data processing results from XIA2'
    TASKNAME = 'xia2_run'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin
    ERROR_CODES = { 101 : {'description' : 'XIA2 run directory does not exist' },
                    102 : {'description' : 'The XIA2 job failed with error file' },
                    103 : {'description' : 'No XIA2 run report file found' },
                    104 : {'description' : 'No XIA2 run output experimental data file found' },
                    110 : {'description' : 'Unable to open ispyb.xml' },
                    111 : {'description' : 'Unable to read performance data from ispyb.xml' }
                    }
    WHATNEXT = ['aimless_pipe','molrep_mr','phaser_pipeline']
    CLONEABLE = False
    PERFORMANCECLASS = 'CDataReductionPerformance'


    def runTitle(self,runMode):
      runModeList = self.container.controlParameters.RUN_MODE.qualifiers('enumerators')
      runTitleList = self.container.controlParameters.RUN_MODE.qualifiers('menuText')
      #print 'xia2_run.runTitle',runMode,runModeList,runMode in runModeList
      if runMode in runModeList:
        return 'XIA2: '+runTitleList[runModeList.index(runMode)]
      else:
        return None

    def process(self):
      #print 'xia2_run.process',self.container.inputData.XIA2_DIRECTORY
      runMode = str(self.container.controlParameters.RUN_MODE)
      xdir = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),runMode+'-run')
      if not os.path.exists(xdir):
        xdir = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),runMode)
        if not os.path.exists(xdir):
          self.appendErrorReport(cls=self.__class__,code=101,details=xdir)
          self.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
        
      logfile = os.path.join(xdir,'xia2.html')
      if os.path.exists(logfile):
        shutil.copyfile(logfile,os.path.join(self.workDirectory,'xia2_log.html'))
        shutil.copytree(os.path.join(xdir,'xia2_html'),os.path.join(self.workDirectory,'xia2_html'))
      else:
        self.appendErrorReport(cls=self.__class__,code=103)
        '''
        txtfile =  os.path.join(xdir,'xia2.txt')
        if os.path.exists(txtfile):
          shutil.copyfile(txtfile,os.path.join(self.workDirectory,'report.txt'))
        '''

      xmlfile = os.path.join(xdir,'ispyb.xml')
      #print 'xia2_run.process',xmlfile,os.path.exists(xmlfile)
      if os.path.exists(xmlfile):
        shutil.copyfile(xmlfile,os.path.join(self.workDirectory,'XMLOUT.xml'))

      mtzfiles = glob.glob(os.path.join(xdir,'DATAFILES','*.mtz'))
      self.container.outputData.HKLOUT.setFullPath(os.path.join(self.workDirectory,'hklout.mtz'))
      #print 'xia2_run mtzfiles',mtzfiles
      if len(mtzfiles) > 0:
        shutil.copyfile(mtzfiles[0],self.container.outputData.HKLOUT.__str__())
      else:
        self.appendErrorReport(cls=self.__class__,code=104)

      self.splitHklout0()
      #error = self.splitHklout(['FPHIOUT','DIFFPHIOUT','ABCDOUT'],['FWT,PHWT','DELFWT,PHDELWT','HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB'])


      pluginObj = self.makePluginObject('xia2_integration',pluginTitle=self.runTitle(runMode))
      pluginObj.container.inputData.XIA2_DIRECTORY.setFullPath(xdir)
      pluginObj.process()
      #print 'xia2_run.process',pluginObj.container.outputData.HKLOUT
      if pluginObj.container.outputData.HKLOUT.isSet():
        if len(self.container.outputData.UNMERGEDHKL)==0: self.container.outputData.UNMERGEDHKL.addItem()
        #self.container.outputData.UNMERGEDHKL[0].set(pluginObj.container.outputData.HKLOUT)
        self.container.outputData.UNMERGEDHKL[0].set(os.path.join(self.workDirectory,'UNMERGEDHKL'+ os.path.splitext(pluginObj.container.outputData.HKLOUT.baseName)[1]))
        shutil.copyfile(str( pluginObj.container.outputData.HKLOUT),self.container.outputData.UNMERGEDHKL[0].fullPath.__str__())
      #print 'xia2_run.process UNMERGEDHKL',self.container.outputData.UNMERGEDHKL
        
      for subJob in ['pointless','aimless','ctruncate']:
        pluginObj = self.makePluginObject('xia2_'+subJob)
        pluginObj.container.inputData.XIA2_DIRECTORY.setFullPath(xdir)
        pluginObj.process()

      self.TASKTITLE = self.runTitle(runMode)
      
      self.setPerformance(runMode)
      self.makeXml(runMode)
      
      self.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)
      return CCP4PluginScript.CPluginScript.SUCCEEDED

    def splitHklout0(self):
      import re
      columnList = self.container.outputData.HKLOUT.fileContent.getListOfColumns()
      #print '\n\n** xia2_run columnNames',columnList
      fp = None
      ip = None
      for c in columnList:
        name= c.columnLabel.__str__()
        if fp is None and 'F(+)' in name:
            fp = name
        if ip is None and 'I(+)' in name:
          ip = name
      #print 'xia2_run fp,ip',fp,ip

      self.container.outputData.FREER.set(self.container.outputData.HKLOUT)
      self.container.outputData.FREER.baseName='FreeR.mtz'
      self.container.outputData.FREER.annotation='Free Rfactor'      
      out = [[self.container.outputData.FREER.__str__(),'FreeR_flag','FREER']]

      # output as many data representations as are available
      outputDataObject = None
      relabelColumnFromList = ''
      relabelColumnToList = ''
      contentFlagList = []
      outputFileBaseName = ''
      
      if ip is not None:
        #print '\n\n**ip not none'
        outputDataObject = self.container.outputData.IANOM
        outputDataObject.set(self.container.outputData.HKLOUT)
        outputFileBaseName = 'IANOM-'
        outputDataObject.annotation='Exptl data as anomalous intensities'
        self.container.outputData.IANOM.contentFlag = 1
        self.container.outputData.IANOM.subType = 1
        im = re.sub('\+','-',ip)
        #MN change: because internally, we assume labels of Iplus Iminus etc for this contentFlag and data type, and we have to notify system that this is what we have used
        relabelColumnFromList += (ip+','+'SIG'+ip+','+im+','+'SIG'+im)
        relabelColumnToList += 'Iplus,SIGIplus,Iminus,SIGIminus'
        contentFlagList += [1]
      
      elif fp is not None:
        #print '\n\n**fp not none'
        if outputDataObject is None:
            outputDataObject = self.container.outputData.FANOM
            outputDataObject.set(self.container.outputData.HKLOUT)
            outputFileBaseName = ''
            outputDataObject.annotation = ''
        else:
            relabelColumnFromList += ','
            relabelColumnToList += ','

        outputFileBaseName += 'FANOM-'
        outputDataObject.annotation = 'Exptl data as anomalous SFs'
        self.container.outputData.FANOM.contentFlag = 2
        self.container.outputData.FANOM.subType = 1
        fm = re.sub('\+','-',fp)
        #MN change: because internally, we assume labels of Fplus Fminus etc for this contentFlag and data type, and we have to notify system that this is what representation we have used
        relabelColumnFromList += (fp+','+'SIG'+fp+','+fm+','+'SIG'+fm+',')
        relabelColumnToList += 'Fplus,SIGFplus,Fminus,SIGFminus,'
        contentFlagList += [2]

      else:
        if outputDataObject is None:
          outputDataObject = self.container.outputData.FSIGF
          outputDataObject.set(self.container.outputData.HKLOUT)
          outputFileBaseName = ''
          outputDataObject.annotation = ''
        else:
          relabelColumnFromList += ','
          relabelColumnToList += ','

        outputFileBaseName += 'FSIGF-'
        outputDataObject.annotation += 'Exptl data as SFs'
        relabelColumnFromList += 'F,SIGF'
        relabelColumnToList += 'F,SIGF'
        contentFlagList += [4]

      outputFileBaseName += 'observed_data.mtz'
      outputDataObject.baseName = outputFileBaseName
      
      out.append([outputDataObject.__str__(), relabelColumnFromList, relabelColumnToList ] )
      #MN change this...contentFlag needs to reflect possiblity of containing multiple data representations...should either be a list or a logically ORed set of flags
      outputDataObject.contentFlag = min(contentFlagList)

      err = self.splitMtz(self.container.outputData.HKLOUT.__str__(),out)
      print(err)

      
    def setPerformance(self,runMode):
      ispyb = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),runMode+'-run','ispyb.xml')
      #print 'xia2_run.setPerformance',ispyb
      try:
        from core import CCP4Utils
        xml = CCP4Utils.openFileToEtree(ispyb)
      except:
        self.appendErrorReport(110,ispyb,stack=False)
      else:
        try:
          self.container.outputData.PERFORMANCE.spaceGroup.set(xml.xpath('//spaceGroup')[0].text)
          self.container.outputData.PERFORMANCE.highResLimit.set(float(xml.xpath('//resolutionLimitHigh')[0].text))
          self.container.outputData.PERFORMANCE.rMeas.set(float(xml.xpath('//rMeasAllIPlusIMinus')[0].text))
        except:
          self.appendErrorReport(111,ispyb,stack=False)
         
    def makeXml(self,runMode):
      from lxml import etree
      from core import CCP4Utils
      self.xmlroot = etree.Element('XIA2Import')
      runXML = etree.SubElement(self.xmlroot,'XIA2Run',name=str(runMode))
      ispyb = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),runMode+'-run','ispyb.xml')
      try:
        ispybxml = CCP4Utils.openFileToEtree(ispyb)
      except:
        self.appendErrorReport(110,ispyb,stack=False)
      else:
        runXML.append(ispybxml)
        for programName in ['pointless','aimless','truncate']:
          programEtree = self.harvestLogXML(runMode, programName)
          if programEtree is not None: runXML.append(programEtree)

      CCP4Utils.saveEtreeToFile(self.xmlroot,self.makeFileName('PROGRAMXML'))

    def harvestLogXML(self, runMode, programName):
      logFiles = glob.glob( os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(), runMode+'-run', 'LogFiles','') + "*" + programName + ".log" )
      #print 'harvestLogXML',logFiles
      if len(logFiles)==0: return None        
      root = etree.Element(programName.upper())
      
      try:
        import smartie
      except:
        from core import CCP4Utils
        smartiePath = os.path.join(CCP4Utils.getCCP4I2Dir(),'smartie')
        import sys
        sys.path.append(smartiePath)
        import smartie
                       
      from pimple import MGQTmatplotlib
      logfile = smartie.parselog(logFiles[0])
      for smartieTable in logfile.tables():
        if smartieTable.ngraphs() > 0:
          tableelement = MGQTmatplotlib.CCP4LogToEtree(smartieTable.rawtable() )
          root.append(tableelement)

      summaryCount = logfile.nsummaries()
      for iSummary in range(summaryCount):
        summary = logfile.summary(iSummary)
        summaryTextLines = []
        with open(logFiles[0]) as myLogFile:
            summaryTextLines = myLogFile.readlines()[summary.start():summary.end()-1]
        preElement = etree.SubElement(root,'CCP4Summary')
        preElementText = ''
        for summaryTextLine in summaryTextLines:
          preElementText += (summaryTextLine)
          preElement.text= base64.b64encode(preElementText)
                    
      return root
  
  
