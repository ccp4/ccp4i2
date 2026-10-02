import glob
import os
import shutil

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class AlternativeImportXIA2(CPluginScript):

    TASKNAME = 'AlternativeImportXIA2'

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

        for runName, dirPath in self.xia2Runs():
            runXML = etree.SubElement(self.xmlroot,'XIA2Run',name=str(runName))
            destDirPath = self.workDirectory
            filePrefix = runName[:-4] if runName.endswith('-run') else runName

            # Grab digested ispyb XML
            fileNameIfAny = os.path.join(dirPath, "ispyb.xml")
            if os.path.isfile(fileNameIfAny):
                runXML.append(CCP4Utils.openFileToEtree(fileNameIfAny))
            for programName in ['pointless','aimless','truncate']:
                programEtree = self.harvestLogXML(runName, programName, dirPath)
                if programEtree is not None: runXML.append(programEtree)
        
            #Grab integrated (unmerged) files
            import sys
            # By what is there, not by the run's name (only "3d..", "2d.."
            # and "dials.." names were recognised; anything else gave no
            # pattern and a TypeError): in DataFiles/Integrate (older xia2)
            # or DataFiles, MTZ (DIALS, MOSFLM) or XDS's INTEGRATE.HKL.
            possibleFilesToCopy = []
            for where in (os.path.join(dirPath, 'DataFiles', 'Integrate'), os.path.join(dirPath, 'DataFiles')):
                for suffix in ('*INTEGRATE.mtz', '*INTEGRATE.HKL'):
                    possibleFilesToCopy = possibleFilesToCopy or sorted(glob.glob(os.path.join(where, suffix)))
            if len(possibleFilesToCopy) != 0:
                try:
                    srcPath = possibleFilesToCopy[0]
                    srcFilename = os.path.split(srcPath)[1]
                    destPath = os.path.join(destDirPath, filePrefix+'_'+srcFilename)
                    shutil.copyfile(srcPath, destPath)
                    unmergedOut.append(unmergedOut.makeItem())
                    unmergedOut[-1].fullPath = destPath
                    unmergedOut[-1].annotation = 'xia2 run '+runName+': integrated, unmerged ('+srcFilename+')'
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
                    allPath = os.path.join(destDirPath, filePrefix+'_'+srcFilename[:-9]+'_all.mtz')
                    shutil.copyfile(srcPath,allPath)
                    obsPath = os.path.join(destDirPath, filePrefix+'_'+srcFilename[:-9]+'_obs.mtz')
                    freerPath = os.path.join(destDirPath, filePrefix+'_'+srcFilename)
                    colin = 'I(+),SIGI(+),I(-),SIGI(-)'
                    colout = 'Iplus,SIGIplus,Iminus,SIGIminus'
                    colfree = 'FreeR_flag'
                    colfreeout = 'FREER'
                    logFile = os.path.join(self.workDirectory,'cmtzsplit.log')
                    status = self.splitMtz(srcPath,[[obsPath,colin,colout],[freerPath,colfree,colfreeout]],logFile)
                    if status == CPluginScript.SUCCEEDED:
                        from ccp4i2.core import CCP4XtalData
                        obsOut.append(obsOut.makeItem())
                        obsOut[-1].fullPath = obsPath
                        obsOut[-1].annotation = 'xia2 run '+runName+': merged intensities'
                        obsOut[-1].contentFlag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
                        freerOut.append(freerOut.makeItem())
                        freerOut[-1].fullPath = freerPath
                        freerOut[-1].annotation = 'xia2 run '+runName+': free R set'
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

    def xia2Runs(self):
        """(name, directory) of each xia2 run to import.

        The runs listed in runSummaries, under directoryPath (as the Qt
        interface filled them); else, as the new interface leaves them empty,
        found from XIA2_DIRECTORY: the directory itself when it is one xia2
        run (it has DataFiles), else each sub-directory that is."""
        names = [n for n in (str(r).split(':')[0].strip() for r in self.container.controlParameters.runSummaries) if n]
        if names:
            base = str(self.container.controlParameters.directoryPath)
            return [(name, os.path.join(base, name)) for name in names]
        base = str(self.container.inputData.XIA2_DIRECTORY.getFullPath() or
                   self.container.controlParameters.directoryPath)
        if os.path.isdir(os.path.join(base, 'DataFiles')):
            return [(os.path.basename(os.path.normpath(base)), base)]
        return [(name, os.path.join(base, name)) for name in sorted(os.listdir(base))
                if os.path.isdir(os.path.join(base, name, 'DataFiles'))] if os.path.isdir(base) else []

    def harvestLogXML(self, runName, programName, dirPath):
        pattern = os.path.join(dirPath, 'LogFiles', '') + "*" + programName + ".log"
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
            
            from lxml import etree

            from ccp4i2.smartie import smartie
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
        from ccp4i2.pimple.logtable import CCP4LogToEtree
        tableetree = CCP4LogToEtree(table.rawtable())
        parent.append(tableetree)
        return tableetree

    
# Function called from gui to support exporting MTZ files
def exportJobFile(jobId=None,mode=None,fileInfo={}):
    import glob
    import os

    from ccp4i2.core import CCP4Modules

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
