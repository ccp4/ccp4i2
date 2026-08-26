import re

from ccp4i2.core.CCP4PluginScript import CPluginScript


class mergeMtz(CPluginScript):
    TASKNAME = 'mergeMtz'
    ERROR_CODES = {
        201: {'description': 'No input reflection files to merge'},
    }

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

      if not inFiles:
        # joinMtz with nothing to join succeeds and writes no file, so without
        # this the job reports success and produces no HKLOUT. A wrapper that
        # computes the inputs to its own step has to say when that computation
        # comes up empty.
        self.appendErrorReport(
            201,
            'No input reflection files to merge. %d were given; none had a '
            'file name that could be read.'
            % len(self.container.inputData.MINIMTZINLIST))
        return CPluginScript.FAILED

      rv = self.joinMtz(self.container.outputData.HKLOUT.fullPath.__str__(),inFiles)
      return rv
