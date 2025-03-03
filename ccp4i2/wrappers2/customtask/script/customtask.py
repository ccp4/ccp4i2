from __future__ import print_function


"""
    customtaskw.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the
     license to address the requirements of UK law.

     You should have received a copy of the modified GNU Lesser General
     Public License along with this library.  If not, copies may be
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""
import os,shutil,glob

from core.CCP4PluginScript import CPluginScript
from core import CCP4CustomTaskManager
from core import CCP4Modules
from core import CCP4Utils


class customtask(CPluginScript):

  TASKMODULE = None
  ERROR_CODES = { 101 : { 'description' : 'ERROR writing command file' },                  
                  102 : { 'description' : 'ERROR merging monster input MTZ' },
                  103 : { 'description' : 'ERROR output file not found' },
                  104 : { 'description' : 'ERROR multiple hits found for output file search' },
                  105 : { 'description' : 'ERROR copying input file' }
                  }

  def process(self):
    manager = CCP4Modules.CUSTOMTASKMANAGER()
    customDef = CCP4CustomTaskManager.CCustomTaskDefinition(self,name=self.TASKNAME)
    customDef.loadDataFromXml(fileName=manager.getCustomFile(name=self.TASKNAME))
    #print 'customtask.process container',customDef

    #convertedMiniMtzs = self.convertInputMiniMtzs()
    #print 'customtask.process convertedMiniMtzs',convertedMiniMtzs
    #Copy any input files to specific filepath
    # Copy any input files requiring specific name
    for paramObj in customDef.paramList:
      if paramObj.function in ['input'] and paramObj.outputFilePath.isSet():
         dobj = self.container.outputData.get(paramObj.name.__str__())
         if dobj is not None and dobj.isSet() and dobj.exists():
           targetPath= os.path.normpath(os.path.join(self.workDirectory, paramObj.outputFilePath.__str__()))
           print('CustomTask: Copying input file from',dobj.fullPath,'to',targetPath)
           try:
             shutil.copyfile(dobj.fullPath,targetPath)
           except:
             self.appendErrorReport(105,dobj.fullPath+' to '+targetPath)
           
    # merge input MTZs
    mergedMtzs = manager.getMergedMtzs(customDef.paramList,function='input')
    mergedMtzFileNames = {}
    for monsterName,miniNameList in list(mergedMtzs.items()):
      miniConvList = []
      for mini in miniNameList:
        miniConvList.append([mini,self.container.inputData.get(mini).qualifiers('requiredContentFlag')])
      print('CustomTask merging miniMTZs to',monsterName,miniConvList)
      fileName,err = self.makeHklin(miniConvList,hklin=monsterName)
      self.extendErrorReport(err)
      if os.path.exists(fileName):
        mergedMtzFileNames[monsterName] = fileName
      else:
        self.appendErrorReport(102,fileName)
    print('CustomTask created merged MTZs',mergedMtzFileNames)

    # Create commandline
    comLine,err = manager.substituteParams(customDef.comLine.__str__(),customDef.tagCharacter.__str__(),self.container,mergedMtzFileNames)
    self.extendErrorReport(err)
    print('CustomTask command line', comLine)
    if len(err)>0: print('CustomTask command line substitution errors',err.report())
    comLineWords = comLine.split()
    self.command = comLineWords[0]
    self.appendCommandLine(comLineWords[1:])


    for comFile in customDef.comFileList:
      if comFile.text.isSet():
        #print 'customtask.process making comFile',comFile.name
        comText,err =  manager.substituteParams(comFile.text.__str__(),customDef.tagCharacter.__str__(),self.container,mergedMtzFileNames)
        if comFile.name.isSet():
          comFileName = os.path.normpath(os.path.join(self.workDirectory,comFile.name.__str__()))
        else:
          comFileName = os.path.normpath(os.path.join(self.workDirectory,'com.txt'))
        print('CustomTask command text',comText)
        print('CustomTask command filename',comFileName)
        try:
          CCP4Utils.saveFile(comFileName,comText)
        except:
          self.appendErrorReport(101,comFileName)

    # Run the task
    status = self.startProcess(self.command,reportStatus=False)
    if not status == CPluginScript.SUCCEEDED:
      # The reportStatus() should be called in CRunPlugin
      #self.reportStatus(status) 
      return status
    
    #shutil.copytree('/Users/lizp/Desktop/test_projects/test7/CCP4_JOBS/job_89/20130823_143704',
    #                os.path.join(self.workDirectory,'20130823_143704'))

    # Find and copy output files
    for paramObj in customDef.paramList:
      if paramObj.function in ['output','log'] and paramObj.outputFilePath.isSet():
        os.chdir(self.workDirectory)
        globFileList = glob.glob(paramObj.outputFilePath.__str__())
        print('CustomTask search for output files',paramObj.name,globFileList)
        if len(globFileList)==0:
          self.appendErrorReport(103,'Parameter:'+paramObj.name.__str__()+' Search pattern:'+paramObj.outputFilePath.__str__())
        else:
          if len(globFileList)>1:
            self.appendErrorReport(104,'Parameter:'+paramObj.name.__str__()+' Search pattern:'+paramObj.outputFilePath.__str__())
          dobj = self.container.outputData.get(paramObj.name.__str__())
          dobj.setOutputPath(jobName=self.TASKNAME,projectId=self.projectId(),relPath=self.getWorkDirectory(ifRelPath=True))
          shutil.copyfile(globFileList[0],dobj.__str__())
          print('CustomTask copied',paramObj.name,dobj.__str__())
          
    # Split output MTZ files
    mergedMtzs = manager.getMergedMtzs(customDef.paramList,function='output')
    for monsterName,miniNameColumnList in list(mergedMtzs.items()):
      monsterFileName = self.container.outputData.get(monsterName).__str__()
      nameList = []
      columnList = []
      for name,col in miniNameColumnList:
        nameList.append(name)
        columnList.append(col)
      print('CustomTask splitting output MTZ',monsterName,nameList,columnList)
      self.splitHklout(miniMtzsOut=nameList,programColumnNames=columnList,infile=monsterFileName,logFile=None)
    
    
    status = CPluginScript.SUCCEEDED
    self.reportStatus(status) 
    return status
