from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
import os, sys
from lxml import etree
from core import CCP4Utils
import shutil
import glob

class AlternativeImportXIA2(CPluginScript):

    TASKNAME = 'AlternativeImportXIA2'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()
        
        from lxml import etree
        self.xmlroot = etree.Element('XIA2Import')

        unmergedOut =  self.container.outputData.UNMERGEDOUT
        obsOut =  self.container.outputData.HKLOUT
        freerOut =  self.container.outputData.FREEROUT

        for runSummary in self.container.controlParameters.runSummaries:
            runName = runSummary.split(':')[0]
            runXML = etree.SubElement(self.xmlroot,'XIA2Run',name=str(runName))
            dirPath = os.path.join(self.container.controlParameters.directoryPath.__str__(), runName)
            destDirPath = self.workDirectory

            # Grab digested ispyb XML
            fileNameIfAny = os.path.join(dirPath, "ispyb.xml")
            if os.path.isfile(fileNameIfAny):
                runXML.append(CCP4Utils.openFileToEtree(fileNameIfAny))
            for programName in ['pointless','aimless','truncate']:
                programEtree = self.harvestLogXML(runName, programName)
                if programEtree is not None: runXML.append(programEtree)
        
            #Grab integrated (unmerged) files
            import sys
            pattern = None
            if runName.startswith('3d'):
                pattern = os.path.join(dirPath,'DataFiles','Integrate','')+'*INTEGRATE.HKL'
            elif runName.startswith('2d') or runName.startswith('dials'):
                pattern = os.path.join(dirPath,'DataFiles','Integrate','')+'*INTEGRATE.mtz'
            possibleFilesToCopy = glob.glob(pattern)
            #At some point these files have been moved out of an Integrate subdirectory
            #and into the "DataFiles" directory...not sure how to know when this may have happened
            #so just try searching in alternative location if none found in Integrate subdirectory
            if len(possibleFilesToCopy) == 0:
                if runName.startswith('3d'):
                    pattern = os.path.join(dirPath,'DataFiles','')+'*INTEGRATE.HKL'
                elif runName.startswith('2d') or runName.startswith('dials'):
                    pattern = os.path.join(dirPath,'DataFiles','')+'*INTEGRATE.mtz'
            possibleFilesToCopy = glob.glob(pattern)
            if len(possibleFilesToCopy) != 0:
                try:
                    srcPath = possibleFilesToCopy[0]
                    srcFilename = os.path.split(srcPath)[1]
                    destPath = os.path.join(destDirPath, runName[0:-4]+'_'+srcFilename)
                    shutil.copyfile(srcPath, destPath)
                    unmergedOut.append(unmergedOut.makeItem())
                    unmergedOut[-1].fullPath = destPath
                    unmergedOut[-1].annotation = runName+' of '+srcFilename[:-13]
                except:
                    print('Unable to import unmerged')
            else:
                print('Unable to find unmerged data to import for run ', runName)
                    
            #Grab merged files
            pattern = os.path.join(dirPath,'DataFiles','')+'*free.mtz'
            possibleFilesToCopy = glob.glob(pattern)
            if len(possibleFilesToCopy) != 0:
                try:
                    srcPath = possibleFilesToCopy[0]
                    srcFilename = os.path.split(srcPath)[1]
                    # Need original file for export
                    allPath = os.path.join(destDirPath, runName[0:-4]+'_'+srcFilename[:-9]+'_all.mtz')
                    shutil.copyfile(srcPath,allPath)
                    obsPath = os.path.join(destDirPath, runName[0:-4]+'_'+srcFilename[:-9]+'_obs.mtz')
                    freerPath = os.path.join(destDirPath, runName[0:-4]+'_'+srcFilename)
                    colin = 'I(+),SIGI(+),I(-),SIGI(-)'
                    colout = 'Iplus,SIGIplus,Iminus,SIGIminus'
                    colfree = 'FreeR_flag'
                    colfreeout = 'FREER'
                    logFile = os.path.join(self.workDirectory,'cmtzsplit.log')
                    status = self.splitMtz(srcPath,[[obsPath,colin,colout],[freerPath,colfree,colfreeout]],logFile)
                    if status == CPluginScript.SUCCEEDED:
                        from core import CCP4XtalData
                        obsOut.append(obsOut.makeItem())
                        obsOut[-1].fullPath = obsPath
                        obsOut[-1].annotation = runName+' of '+srcFilename[:-8]
                        obsOut[-1].contentFlag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
                        freerOut.append(freerOut.makeItem())
                        freerOut[-1].fullPath = freerPath
                        freerOut[-1].annotation = 'FreeR from '+runName+' of '+srcFilename[:-8]
                    else:
                        print('CSplitMTZ Failed')
                except:
                    print('Unable to import merged')
            else:
                print('Unable to find merged data to import for run ', runName)
    
            with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
                CCP4Utils.writeXML(xmlFile,etree.tostring(self.xmlroot, pretty_print=True))

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def harvestLogXML(self, runName, programName):
        pattern = os.path.join(self.container.controlParameters.directoryPath.__str__(), runName, 'LogFiles','') + "*" + programName + ".log"
        candidateFiles = glob.glob(pattern)
        pointlessEtree = None
        if len(candidateFiles) > 0:
            fromPath = candidateFiles[0]
            fileRoot = os.path.split(fromPath)[1]
            try:
                os.mkdir(os.path.join(self.workDirectory, runName))
            except:
                print('Directory already exists:',os.path.join(self.workDirectory, runName))
            toPath = os.path.join(self.workDirectory, runName, fileRoot)
            shutil.copyfile(fromPath, toPath)
            
            smartiePath = os.path.join(CCP4Utils.getCCP4I2Dir(),'smartie')
            sys.path.append(smartiePath)
            import smartie
            
            from lxml import etree
            pointlessEtree = etree.Element(programName.upper())
            
            logfile = smartie.parselog(toPath)
            for smartieTable in logfile.tables():
                if smartieTable.ngraphs() > 0:
                    tableelement = self.xmlForSmartieTable(smartieTable, pointlessEtree)

            summaryCount = logfile.nsummaries()
            for iSummary in range(summaryCount):
                summary = logfile.summary(iSummary)
                summaryTextLines = []
                with open(toPath) as myLogFile:
                    summaryTextLines = myLogFile.readlines()[summary.start():summary.end()-1]
                preElement = etree.SubElement(pointlessEtree,'CCP4Summary')
                preElementText = ''
                for summaryTextLine in summaryTextLines:
                    preElementText += (summaryTextLine)
                preElement.text=etree.CDATA(preElementText)
                    
        return pointlessEtree
                
    def xmlForSmartieTable(self, table, parent):
        from pimple.logtable import CCP4LogToEtree
        tableetree = CCP4LogToEtree(table.rawtable())
        parent.append(tableetree)
        return tableetree

    
# Function called from gui to support exporting MTZ files
def exportJobFile(jobId=None,mode=None,fileInfo={}):
    import os,glob
    from core import CCP4Modules
    #print 'AlternativeImportXIA2.exportJobFile',mode
    if mode == 'complete_mtz':
      if fileInfo.get('fullPath',None) is not None:
        exportFile = fileInfo['fullPath'][0:-7]+'all.mtz'
        print('AlternativeImportXIA2.exportJobFile',exportFile)
        if os.path.exists(exportFile):
          return exportFile
        else:
          return None
      else:
        jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
        allMtzs = glob.glob(os.path.join(jobDir,'*_all.mtz'))
        #print 'AlternativeImportXIA2.exportJobFile',jobDir,allMtzs
        if len(allMtzs)>0:
          return allMtzs[0]
        else:
          return None
    
# Function to return list of names of exportable MTZ(s)
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]


    

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testProvideTLS(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from core.CCP4Modules import QTAPPLICATION
     wrapper = ProvideTLS(parent=QTAPPLICATION(),name='ProvideTLS_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testProvideTLS)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
