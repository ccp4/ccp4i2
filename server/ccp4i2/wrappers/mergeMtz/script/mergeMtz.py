import re

from ccp4i2.core.CCP4PluginScript import CPluginScript


class mergeMtz(CPluginScript):
    TASKTITLE = 'Merge experimental data objects to MTZ'     # A short title for gui menu
    TASKNAME = 'mergeMtz'                                  # Task name - should be same as class name

    def startProcess(self):
      inFiles = []
      for miniMtz in self.container.inputData.MINIMTZINLIST:
        if miniMtz.fileName.isSet() and miniMtz.fileName.exists():
          cls,contentFlag =  miniMtz.fileName.miniMtzType()
          if cls is not None:
            # Create instance of class to use the columnNames() method
            stdColumnNames = cls().columnNames(True,contentFlag)
            userColumnNames = ''
            if miniMtz.columnNames.isSet():
              userColumnNames = re.sub(' ','',miniMtz.columnNames.__str__())
              if not userColumnNames.count(',') == stdColumnNames.count(','): userColumnNames = ''
            if len(userColumnNames)==0:
              inFiles.append([miniMtz.fileName.__str__(),stdColumnNames])
            else:
              inFiles.append([miniMtz.fileName.__str__(),userColumnNames])

      rv = self.joinMtz(self.container.outputData.HKLOUT.fullPath.__str__(),inFiles)
      return rv
