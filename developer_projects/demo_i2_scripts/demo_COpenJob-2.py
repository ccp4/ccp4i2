from __future__ import print_function

# Run with 
# CCP4I2_TOP=$CCP4/share/ccp4i2 $CCP4/share/ccp4i2/bin/pyi2 pathto/demo_i2_scripts/demo_COpenJob-2.py

import sys,os,shutil
from PyQt4 import QtCore

def bootI2(dbDir):
  # Bootstrap ccp4i2 environment - mimics the gui side but without actual graphics
  # Specify  directory for database and its temp files and make sure it is empty
  if os.path.exists(dbDir):  shutil.rmtree(dbDir)
  os.mkdir(dbDir)
  dbFile = os.path.join(dbDir,'db.sqlite')
  
  ccp4i2_top = os.environ['CCP4I2_TOP']
  sys.path.append(os.path.join(ccp4i2_top,'utils'))
  from startup import setupEnvironment,setupPythonpath,startProjectsManager,startJobController
  setupEnvironment()
  setupPythonpath(top=ccp4i2_top,mode='qtcore')
  from core import CCP4Modules 
  app = CCP4Modules.QTAPPLICATION(graphical=False)
  pm = startProjectsManager(dbFileName=dbFile)
  pm.startCheckForFinishedJobs()
  jc = startJobController()
  jc.setDiagnostic(True)
  if dbFile is not None: jc.setDbFile(dbFile)
  pm.doCheckForFinishedJobs.connect(pm.checkForFinishedJobs)
  from core.CCP4DataManager import DATAMANAGER
  DATAMANAGER().buildClassLookup()
  return app

class Runner:
  def __init__(self,projectName,projectPath,sourceList):
    # Create a project
    from core import CCP4Modules
    if os.path.exists(projectPath): shutil.rmtree(projectPath)
    self.projectId =  CCP4Modules.PROJECTSMANAGER().createProject(projectName='myproject',projectPath=projectPath)
    print('projectId',self.projectId)
    self.sourceList = sourceList
    self.jobIndex = -1

  def run(self):
    from core import CCP4Modules
    CCP4Modules.PROJECTSMANAGER().db().jobFinished.connect(runner.handleJobFinished)
    self.runJob()

  def handleJobFinished(self,args):
    print('handleJobFinished',args)
    if self.jobIndex+1>=len(self.sourceList):
      print('FINISHED ALL JOBS')
      sys.exit()
    else:
      self.runJob()

  def runJob(self):
    from dbapi import CCP4DbUtils
    self.jobIndex = self.jobIndex + 1
    # Create instance of COpenJob which will create/run a job for us
    self.jobHandler = CCP4DbUtils.COpenJob(projectId=self.projectId)
    self.jobHandler.createJob(taskName=self.sourceList[self.jobIndex][0])
    print('jobId',self.jobIndex,self.jobHandler.jobId)

    self.loadData()

    # Run job (as a separate process)
    self.jobHandler.runJob()

  def loadData(self):
    # Set input parameters to job
    sourceDirectory = self.sourceList[self.jobIndex][1]
    inp = self.jobHandler.container.inputData
    inp.F_SIGF.setFullPath(os.path.join(sourceDirectory,'F_SIGF.mtz'))
    inp.FREERFLAG.setFullPath(os.path.join(sourceDirectory,'FREERFLAG.mtz'))
    inp.BUCCANEER_MR_MODE_XYZIN.setFullPath(os.path.join(sourceDirectory,'model.pdb'))
    inp.SEQIN.addItem()
    inp.SEQIN[0].setFullPath(os.path.join(sourceDirectory,'seq.fasta'))
    self.jobHandler.container.controlParameters.ITERATIONS=2


# ==========================================================================================================

# Where we want database and project directory to go
dbDir = '/Users/lizp/Desktop/demo_scripts_db'
projectDirectory = '/Users/lizp/Desktop/demo_scripts_project'
# List of tasks and directories containing input data
# This would assume each directories contains set of input files with the same names
sourceList = [[ 'buccaneer_build_refine_mr','/Users/lizp/Desktop/demo_i2_scripts'] ]

app = bootI2(dbDir)

runner = Runner(projectName='myproject', projectPath=projectDirectory,sourceList=sourceList)
runner.run()

sys.exit(app.exec_())
