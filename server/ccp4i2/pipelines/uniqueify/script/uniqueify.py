import shutil

from ccp4i2.baselayer import QtCore
from ccp4i2.core.CCP4PluginScript import CPluginScript


class uniqueify(CPluginScript):

    TASKMODULE = 'test'      # Where this plugin will appear on the gui
    TASKTITLE = 'Standardise reflection data file' # Short title for gui menu
    TASKNAME = 'uniqueify'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    def process(self):
      # Run mtzheader to get spg, cell, reso
      print("### Running mtzheader ###")

      unsetData = self.checkInputData()
      if len(unsetData)>0:
         self.reportStatus(CPluginScript.FAILED)
         return

      self.mtzheader = self.makePluginObject('mtzheader')
      self.mtzheader.container.inputData.HKLIN = self.container.inputData.HKLIN

      self.connectSignal(self.mtzheader,'finished',self.process_unique)
      self.mtzheader.process()

    @QtCore.Slot(dict)
    def process_unique(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      ### determine workflow
      if self.container.controlParameters.COMPLETE == "Yes":
         self.route = "cad_copy_column"
      elif self.container.controlParameters.COMPLETE == "No":
         self.route = "freerflag"
      elif self.container.controlParameters.COMPLETE == "Auto":
         print("Auto not implemented yet")
         self.reportStatus(CPluginScript.FAILED)
         return

      print("### Running unique ###")

      # Run unique to get full list of reflections
      self.unique = self.makePluginObject('unique')   
      self.unique.container.inputData.SPACEGROUPCELL = self.mtzheader.container.outputData.SPACEGROUPCELL
      self.unique.container.inputData.RESOLUTION = self.mtzheader.container.outputData.MAXRESOLUTION

      if self.route == "freerflag":
         self.connectSignal(self.unique,'finished',self.process_freerflag)
      if self.route == "cad_copy_column":
         self.connectSignal(self.unique,'finished',self.process_cad2)
      self.unique.process()

    @QtCore.Slot(dict)
    def process_freerflag(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running freerflag ###")

      # Run freerflag to attach freeR column
      self.freerflag = self.makePluginObject('freerflag')
      self.freerflag.container.inputData.HKLIN = self.unique.container.outputData.HKLOUT

      self.connectSignal(self.freerflag,'finished',self.process_cad)
      self.freerflag.process()

    @QtCore.Slot(dict)
    def process_cad(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running cad_copy_column ###")

      # Run cad to attach freeR column to original file 
      self.cad_copy_column = self.makePluginObject('cad_copy_column')   
      self.cad_copy_column.container.inputData.HKLIN1 = self.container.inputData.HKLIN
      self.cad_copy_column.container.inputData.HKLIN2 = self.freerflag.container.outputData.HKLOUT

      # for the moment, we assume freerflag has outputted standard label
      self.cad_copy_column.container.inputData.FREERFLAG.FREE = 'FreeR_flag'

      self.connectSignal(self.cad_copy_column,'finished',self.process_freerflag_tidy)
      self.cad_copy_column.process()
      
    @QtCore.Slot(dict)
    def process_freerflag_tidy(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running freerflag ###")

      # Run freerflag to ensure freeR column is complete
      # (experience shows that uniqueify can miss some reflections
      # on resolution limit)
      self.freerflag_tidy = self.makePluginObject('freerflag')
      self.freerflag_tidy.container.inputData.HKLIN = self.cad_copy_column.container.outputData.HKLOUT
      self.freerflag_tidy.container.inputData.FREERFLAG.FREE = 'FreeR_flag'
      self.freerflag_tidy.container.controlParameters.COMPLETE = True

      self.connectSignal(self.freerflag_tidy,'finished',self.process_finish)
      self.freerflag_tidy.process()

    @QtCore.Slot(dict)
    def process_cad2(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running cad_copy_column ###")

      # Run cad to attach unique columns to original file 
      self.cad_copy_column2 = self.makePluginObject('cad_copy_column')   
      self.cad_copy_column2.container.inputData.HKLIN1 = self.container.inputData.HKLIN
      self.cad_copy_column2.container.inputData.HKLIN2 = self.unique.container.outputData.HKLOUT

      # for the moment, we assume unique has outputted standard label
      self.cad_copy_column2.container.inputData.FSIGF.F = 'F'
      self.cad_copy_column2.container.inputData.FSIGF.SIGF = 'SIGF'

      self.connectSignal(self.cad_copy_column2,'finished',self.process_freerflag2)
      self.cad_copy_column2.process()
      
    @QtCore.Slot(dict)
    def process_freerflag2(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running freerflag ###")

      # Run freerflag to ensure freeR column is complete
      self.freerflag2 = self.makePluginObject('freerflag')
      self.freerflag2.container.inputData.HKLIN = self.cad_copy_column2.container.outputData.HKLOUT
      self.freerflag2.container.inputData.FREERFLAG.FREE = 'FreeR_flag'
      self.freerflag2.container.controlParameters.COMPLETE = True

      self.connectSignal(self.freerflag2,'finished',self.process_mtzutils)
      self.freerflag2.process()
      
    @QtCore.Slot(dict)
    def process_mtzutils(self,status):
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.reportStatus(status)
         return

      print("### Running mtzutils ###")

      # Run mtzutils to delete spurious columns
      self.mtzutils = self.makePluginObject('mtzutils')
      self.mtzutils.container.inputData.HKLIN1 = self.freerflag2.container.outputData.HKLOUT
      self.mtzutils.container.inputData.FSIGF.F = 'F'
      self.mtzutils.container.inputData.FSIGF.SIGF = 'SIGF'
      self.mtzutils.container.controlParameters.MODE = 'EXCLUDE'

      self.connectSignal(self.mtzutils,'finished',self.process_finish)
      self.mtzutils.process()
      
    @QtCore.Slot(dict)
    def process_finish(self,status):

      if self.route == "freerflag":
        shutil.copyfile(str(self.freerflag_tidy.container.outputData.HKLOUT),str(self.container.outputData.HKLOUT))
      if self.route == "cad_copy_column":
        shutil.copyfile(str(self.mtzutils.container.outputData.HKLOUT),str(self.container.outputData.HKLOUT))

      # Emit finished signal with whatever status has been emitted 
      # by the previous wrapper
      self.reportStatus(status.get('finshStatus'))    
