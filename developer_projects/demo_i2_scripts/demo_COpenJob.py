from __future__ import print_function

# Run with 
# CCP4I2_TOP=$CCP4/share/ccp4i2 $CCP4/share/ccp4i2/bin/pyi2 pathto/demo_i2_scripts/demo_COpenJob.py

import sys,os,shutil
from PyQt4 import QtCore

def handleJobFinished(args):
  print('handleJobFinished',args)
  sys.exit()

# Hard code ccp4i2 source and a directory containing data
ccp4i2_top = os.environ['CCP4I2_TOP']
sourceDirectory = '/Users/lizp/Desktop/demo_i2_scripts'

# Specify a project directory and database file - and make sure they don't exist
projectDirectory = '/Users/lizp/Desktop/demo_scripts_project'
dbDir = '/Users/lizp/Desktop/demo_scripts_db'
if os.path.exists(projectDirectory): shutil.rmtree(projectDirectory)
if os.path.exists(dbDir):  shutil.rmtree(dbDir)
os.mkdir(dbDir)
dbFile = os.path.join(dbDir,'db.sqlite')

# Bootstrap ccp4i2 environment - mimics the gui side but without actual graphics
# Note this sets the specified database file
sys.path.append(os.path.join(ccp4i2_top,'utils'))
from startup import setupEnvironment,setupPythonpath,startProjectsManager,startJobController
setupEnvironment()
setupPythonpath(top=ccp4i2_top,mode='qtcore')
from core.CCP4Modules import QTAPPLICATION,PROJECTSMANAGER
app = QTAPPLICATION(graphical=False)
pm = startProjectsManager(dbFileName=dbFile)
pm.startCheckForFinishedJobs()
jc = startJobController()
jc.setDiagnostic(True)
if dbFile is not None: jc.setDbFile(dbFile)
#PROJECTSMANAGER().startCheckForFinishedJobs()
pm.doCheckForFinishedJobs.connect(pm.checkForFinishedJobs)
PROJECTSMANAGER().db().jobFinished.connect(handleJobFinished)
from core.CCP4DataManager import DATAMANAGER
DATAMANAGER().buildClassLookup()
from dbapi import CCP4DbUtils
from dbapi import CCP4DbApi

# Create a project
projectId =  PROJECTSMANAGER().createProject(projectName='myproject',projectPath=projectDirectory)
print('projectId',projectId)

# Create instance of COpenJob which will create/run a job for us
jobHandler = CCP4DbUtils.COpenJob(projectId=projectId)
jobHandler.createJob(taskName='buccaneer_build_refine_mr')
print('jobId',jobHandler.jobId)

# Set input parameters to job
inp = jobHandler.container.inputData
inp.F_SIGF.setFullPath(os.path.join(sourceDirectory,'F_SIGF.mtz'))
inp.FREERFLAG.setFullPath(os.path.join(sourceDirectory,'FREERFLAG.mtz'))
inp.BUCCANEER_MR_MODE_XYZIN.setFullPath(os.path.join(sourceDirectory,'model.pdb'))
inp.SEQIN.addItem()
inp.SEQIN[0].setFullPath(os.path.join(sourceDirectory,'seq.fasta'))
jobHandler.container.controlParameters.ITERATIONS=2

# Run job (as a separate process)
jobHandler.runJob()

sys.exit(app.exec_())
