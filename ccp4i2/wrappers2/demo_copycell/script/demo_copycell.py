"""
Copyright (C) 2010 University of York
"""

from PySide2 import QtCore

from ....core.CCP4PluginScript import CPluginScript


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
