from __future__ import print_function

"""
    demo_copycell.py: CCP4 GUI Project
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

from baselayer import QtCore
from core.CCP4PluginScript import CPluginScript

class demo_copycell(CPluginScript):

    TASKMODULE = 'test'
    TASKTITLE = 'Demo copy cell'
    TASKNAME = 'demo_copycell'
    DESCRIPTION='Demo pipeline - copy cell from MTZ file to PDB file'
    TASKVERSION= 0.0
    ASYNCHRONOUS=True

    def process(self):
      # Check all input files exist
      nonExFiles = self.checkInputData()
      if len(nonExFiles)>0:
        self.reportStatus(CPluginScript.FAILED)
        return
      # Provide default output file names if necessary
      self.checkOutputData()
      # Create instance of mtzdump class (and get it registered with the database)
      self.mtzdump = self.makePluginObject('mtzdump')
      # Run mtzdump to get the cell parameters
      self.mtzdump.container.inputData.HKLIN.set(self.container.inputData.HKLIN)
      self.connectSignal(self.mtzdump,'finished',self.process_1)
      self.mtzdump.process()

    @QtCore.Slot(int)
    def process_1(self,status):
      print('demo_copycell.process_1 status',status)
      if status == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      self.pdbset = self.makePluginObject('pdbset')
      self.pdbset.container.inputData.XYZIN.set(self.container.inputData.XYZIN)
      self.pdbset.container.inputData.CELL.set(self.mtzdump.container.outputData.CELL)
      self.pdbset.container.outputData.XYZOUT.set(self.container.outputData.XYZOUT)
      self.connectSignal(self.pdbset,'finished',self.postProcessWrapper)
      self.pdbset.process()
      
      
#=======================================================================================================
import unittest

class test_demo_copycell(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    from core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

  def test_1(self):
    from core.CCP4Modules import QTAPPLICATION
    import os
    from core.CCP4Utils import getCCP4I2Dir

    # Run the pipeline
    wrapper = demo_copycell(parent=QTAPPLICATION(),name='demo_copycell')
    wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(),'pipelines','demo_copycell','test_data','test_1.params.xml'))
    # Ensure no output file exists already
    xyzout = wrapper.container.outputData.XYZOUT.fullPath.get()
    if xyzout is not None and os.path.exists(xyzout): os.remove(xyzout)
    pid = wrapper.process()

    #test if output file created
    self.assertEqual(os.path.exists( xyzout),1,'Failed to create copied pdb file '+xyzout)                             


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(test_demo_copycell)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
