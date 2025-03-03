
from core.CCP4PluginScript import CPluginScript

  
class mergeMtz(CPluginScript):

    TASKTITLE = 'Merge experimental data objects to MTZ'     # A short title for gui menu
    TASKNAME = 'mergeMtz'                                  # Task name - should be same as class name
    RUNEXTERNALPROCESS = False                                  # There is not external process
    MAINTAINER = 'liz.potterton@york.ac.uk'


    def startProcess(self,command,**kw):
      # Reimplement the method that normally starts an external process to do the processing by
      # calling the CPluginScript.joinMtz() method 
      import re
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
            
          
      #print 'mergeMtz.process inFiles',inFiles
      rv = self.joinMtz(self.container.outputData.HKLOUT.fullPath.__str__(),inFiles)
      return rv

      
     
