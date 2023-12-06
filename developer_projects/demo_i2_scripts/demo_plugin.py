from __future__ import print_function


import sys,os,shutil

#Run with ..
#> export CCP4I2_TOP=/Users/lizp/Desktop/dev/ccp4i2-devel
#> $CCP4I2_TOP/bin/pyi2 $CCP4I2_TOP/developer_projects/demo_i2_scripts/demo_plugin.py

def handleFinish():
  print('handleFinish')
  sys.exit()

doAsync = False

# Hardcoded ccp4i2 directory and directory containing data
ccp4i2_top = os.environ['CCP4I2_TOP']
sourceDirectory = os.path.join(ccp4i2_top,'developer_projects','demo_i2_scripts')

# A work directory that results are written to
workDirectory = '/Users/lizp/Desktop/demo_plugin'

# Bootstrap i2 environment - NO DATABASE
sys.path.append(os.path.join(ccp4i2_top,'utils'))
from startup import setupEnvironment,setupPythonpath
setupEnvironment()
setupPythonpath(top=ccp4i2_top,mode='qtcore')
if doAsync:
  from core.CCP4Modules import QTAPPLICATION
  app = QTAPPLICATION(graphical=False)

#Need to make the workDirectory
if os.path.exists(workDirectory): shutil.rmtree(workDirectory)
os.mkdir(workDirectory)

# Create a plugin object
import refmac_i2
if doAsync:
  from PyQt4 import QtCore
  wrapper = refmac_i2.refmac_i2(parent=QTAPPLICATION(),name='test_test',workDirectory=workDirectory)
  wrapper.finished.connect(handleFinish)
else:
  wrapper = refmac_i2.refmac_i2(name='test_test',workDirectory=workDirectory)
# Could load default params from xml file
#wrapper.container.loadDataFromXml(sourceDirectory,'buc_input.xml')


#Creating mini-MTZs - if you are working with monster MTZs and need to split to
# input to i2 wrappers.  Create an instance of CMtzDataFile with name of monster mtz
# then use the runMtzSplit() method
'''
from core import CCP4XtalData
m = CCP4XtalData.CMtzDataFile(os.path.join(sourceDirectory,'monster.mtz')
m.runMtzSplit(['FNAT,SIGFNAT','FreeR_flag'],[os.path.join(workDirectory,'F_SIGF.mtz'),os.path.join(workDirectory,'FREERFLAG.mtz')])
'''

# Set the input data
inp = wrapper.container.inputData
inp.F_SIGF.setFullPath(os.path.join(sourceDirectory,'F_SIGF.mtz'))
inp.FREERFLAG.setFullPath(os.path.join(sourceDirectory,'FREERFLAG.mtz'))
inp.XYZIN.setFullPath(os.path.join(sourceDirectory,'model.pdb'))

for paramName in ['XYZOUT','FPHIOUT','DIFFPHIOUT','ABCDOUT']:
  wrapper.container.outputData.get(paramName).setOutputPath(relPath=workDirectory)

# set the plugin running
status = wrapper.process()

print('Finished demo_plugin',status, wrapper.errorReport.report())

if doAsync:
  sys.exit(app.exec_())
else:
  sys.exit()
