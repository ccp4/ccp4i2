from __future__ import print_function


"""
     workflow.py: CCP4 GUI Project
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
import os,shutil,time

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4WorkflowManager,CCP4Modules
from ccp4i2.core.CCP4ErrorHandling import *
from ccp4i2.core import CCP4Data

from ccp4i2.baselayer import QtCore

class workflow(CPluginScript):

  TASKMODULE = None
  
  ERROR_CODES = { 203 : { 'description' : 'FAILED to load and initialise plugin script' },
                  202 : { 'description' : 'FAILED to read workflow definition file' },
                  201 : { 'description' : 'FAILED to found workflow definition directory' },
                  204 : { 'description' : 'FAILED to copy output file to next program input' },
                  205: { 'severity' : SEVERITY_WARNING, 'description' : 'Can not perform specified file copy' }
                  }
  ASYNCHRONOUS = True
  TIMEOUT_PERIOD = 240
  MAXNJOBS = 4
  def process(self):
    #print 'workflow.process',self.TASKNAME,self.container.inputData.dataOrder()
    self.workflowManager = CCP4WorkflowManager.CWorkflowManager()
    self.workflowDirectory = self.workflowManager.getDirectory(name=self.TASKNAME)
    if not os.path.exists(self.workflowDirectory):
      self.appendErrorReport(201,'Directory: '+self.workflowDirectory)
      self.reportStatus(CPluginScript.FAILED)
      return CPluginScript.FAILED
    #print 'workflow.process',directory
    self.workflowDef = CCP4WorkflowManager.CWorkflowDefinition(name=self.TASKNAME)
    #print 'workflow.process workflowDef',self.workflowDef.dataOrder()
    workflowFile = self.workflowManager.getCustomFile(name=self.TASKNAME)
    err = self.workflowDef.loadDataFromXml(fileName=workflowFile,check=False,function='WORKFLOW')
    if err.maxSeverity()>=SEVERITY_WARNING:
      self.extendErrorReport(err)
      if err.maxSeverity()>SEVERITY_WARNING:
        self.appendErrorReport(202,'File: '+workflowFile)
        self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.FAILED
    #print 'workflow.process workflowDef',self.TASKNAME,self.workflowManager.getCustomFile(name=self.TASKNAME), self.workflowDef.jobDef.dataOrder()
    self.subJobs = {}
    self.maxJobIndex = len(self.workflowDef.jobDef.dataOrder())
    self.currentJobIndex = 0
    self.currentJobKey = None
    self.runNextJob()

  def runNextJob(self):
    # Increment self.currentJobIndex and run that job
    if self.currentJobIndex+1>=self.maxJobIndex: self.finish()
    
    self.currentJobIndex +=1
    self.currentJobKey = self.workflowDef.jobDef.dataOrder()[self.currentJobIndex]
    jobDef = self.workflowDef.jobDef[self.currentJobKey]
    print('workflow.runNextJob jobDef',self.currentJobKey,jobDef.taskName.__str__())
    self.subJobs[self.currentJobKey] = self.makePluginObject(jobDef.taskName.__str__())
    if self.subJobs[self.currentJobKey] is None:
      self.appendErrorReport(203,'For job: '+str(self.currentJobKey))
      self.reportStatus(CPluginScript.FAILED)
      return CPluginScript.FAILED
    
    # Load the control parameters
    err = self.subJobs[self.currentJobKey].container.loadDataFromXml(os.path.join(self.workflowDirectory,self.currentJobKey+'.params.xml'))
    if err.maxSeverity()>SEVERITY_WARNING:
      self.extendErrorReport(err)
      self.reportStatus(CPluginScript.FAILED)
      return CPluginScript.FAILED

    # Is there a paramsFile containing user edit of the control parameters?
    paramsFile = os.path.join(self.workDirectory,self.currentJobKey+'_input_params.xml')
    #print 'workflow.runNextJob paramsFile',paramsFile,os.path.exists(paramsFile)
    if os.path.exists(paramsFile):
      err = self.subJobs[self.currentJobKey].container.loadDataFromXml(paramsFile)
      if err.maxSeverity()>SEVERITY_WARNING: self.extendErrorReport(err)
                                                     
    CCP4Modules.PROJECTSMANAGER().setOutputFileNames(container=self.subJobs[self.currentJobKey].container,projectId=self.projectId(),jobNumber=self.subJobs[self.currentJobKey].getJobNumber())
    #print 'workflow.runNextJob subJobs',self.currentJobKey,jobDef.taskName.__str__(),self.subJobs[self.currentJobKey] 
    inputData = self.subJobs[self.currentJobKey].container.inputData
    for inpDef in jobDef.input:                                             
      print('workflow.runNextJob inpDef',self.currentJobKey,inpDef)
      #try:
      if 1:
        target = inputData.__getattr__(inpDef.toKey.__str__())
        if inpDef.fromJob == 'job_0':
          source = self.container.inputData.__getattr__(inpDef.fromKey.__str__())
        else:
          source = self.subJobs[inpDef.fromJob.__str__()].container.outputData.__getattr__(inpDef.fromKey.__str__())
        #print 'source,target',repr(source),repr(target)
        if isinstance(source,CCP4Data.CList): source = source[0]
        #print 'source,target',repr(source),repr(target)
        target.set(source.get())
      #except:
      #  print 'ERROR in workflow attempting to set',inpDef.fromKey.__str__()

    # Save the parameters to input_params.xml - mostly to help diagnostic
    self.subJobs[self.currentJobKey].container.saveDataToXml(fileName=self.subJobs[self.currentJobKey].makeFileName('JOB_INPUT'))

    #print 'runNextJob',self.currentJobIndex,self.currentJobKey
    self.connectSignal(self.subJobs[self.currentJobKey],'finished',self.handleJobFinished)
    ret = self.subJobs[self.currentJobKey].process()
    #print 'runNextJob done'

  @QtCore.Slot(dict)
  def handleJobFinished(self,statusDict):
    #import time
    #print 'handleJobFinished',self.currentJobIndex,statusDict,time.time()
    finishedJobId = statusDict.get('jobId',None)
    if finishedJobId != self.subJobs[self.currentJobKey].getJobId():
      #print 'workflow.handleJobFinished ignoring finished signal from a CPluginScript that is not the last task'
      import traceback
      traceback.print_stack()
      return
    if statusDict['finishStatus']  == CPluginScript.FAILED: self.finish(CPluginScript.FAILED)
    rv = self.copyOutputFiles()
    if self.currentJobIndex+1>=self.maxJobIndex:
      self.finish(CPluginScript.SUCCEEDED)
    else:
      self.runNextJob()

  def copyOutputFiles(self):
    for outDef in self.workflowDef.jobDef['job_0'].output:
      if outDef.fromJob.__str__() == self.currentJobKey:
        fromObj = self.subJobs[self.currentJobKey].container.outputData.__getattr__(outDef.fromKey.__str__())
        toObj = self.container.outputData.__getattr__(outDef.toKey.__str__())
        #print 'workflow.process fromObj toObj',fromObj,toObj
        if fromObj is not None and toObj is not None:
          #print 'workflow.copyOutputFiles',fromObj.isSet(),fromObj.exists(),os.path.exists(fromObj.__str__()),toObj.isSet()
          if fromObj.isSet() and fromObj.exists() and  toObj.isSet():
            try:
              shutil.copyfile(fromObj.__str__(),toObj.__str__())
              #print 'workflow.copyOutputFiles done copy',toObj.__str__()
            except:
              self.appendErrorReport(204,'From: '+fromObj.__str__()+' to: '+toObj.__str__())
              self.reportStatus(CPluginScript.FAILED)
              return CPluginScript.FAILED
            if fromObj.annotation.isSet():
              toObj.annotation.set(fromObj.annotation.__str__())
          else:
            self.appendErrorReport(205,'From: '+fromObj.__str__()+' to: '+toObj.__str__())
    return CPluginScript.SUCCEEDED

  def finish(self,status):
    #print 'workflow.finish',status
    self.reportStatus(status) 
    

