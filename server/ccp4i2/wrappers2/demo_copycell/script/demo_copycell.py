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

from ccp4i2.baselayer import QtCore
from ccp4i2.core.CCP4PluginScript import CPluginScript

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
