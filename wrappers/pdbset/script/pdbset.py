from __future__ import print_function

"""
     pdbset.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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

from core.CCP4PluginScript import CPluginScript
     
class pdbset(CPluginScript):

    TASKMODULE = 'demo'
    TASKTITLE = 'PDBSet'
    TASKNAME = 'pdbset'
    TASKCOMMAND = 'pdbset'
    TASKVERSION= 0.0
    COMLINETEMPLATE = '''1 XYZIN $XYZIN
1 XYZOUT $XYZOUT'''
    COMTEMPLATE = '''1 CELL $CELL.a $CELL.b $CELL.c $CELL.alpha $CELL.beta $CELL.gamma
1 END'''




     
#====================================================================================================

import unittest

class testPdbset(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

   def testPdbset(self):
     from core.CCP4Modules import QTAPPLICATION
     import os

     wrapper = pdbset(parent=QTAPPLICATION(),name='pdbset')
     
     wrapper.container.inputData.XYZIN.set(project='CCP4I2_TOP',relPath='wrappers/pdbset/test_data',baseName='1df7.pdb')
     wrapper.container.inputData.CELL.set(a=100.0,b=120.0,c=30.0,alpha=90.0,beta=90.0,gamma=89.1)
     wrapper.container.outputData.XYZOUT.set(project='CCP4I2_TEST',baseName='1df7_mangled.pdb')
     print('XYZOUT',wrapper.container.outputData.XYZOUT.fullPath)

     # Ensure no output file exists
     if wrapper.container.outputData.XYZOUT.exists():  os.remove(wrapper.container.outputData.XYZOUT.fullPath.get())
     wrapper.process()

     #test if output file created
     self.assertEqual(wrapper.container.outputData.XYZOUT.exists(),True,'Failed to create copied pdb file')                             

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testPdbset)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
