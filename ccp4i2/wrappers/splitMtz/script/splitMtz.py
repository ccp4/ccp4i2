import os
import sys

from ....core import CCP4XtalData
from ....core.CCP4PluginScript import CPluginScript


class splitMtz(CPluginScript):

    TASKTITLE = 'Import and Split MTZ to experimental data objects'     # A short title for gui menu
    TASKNAME = 'splitMtz'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    RUNEXTERNALPROCESS = False
    MAINTAINER = 'liz.potterton@york.ac.uk'


    def startProcess(self,command,**kw):
      inp = self.container.inputData
      out = self.container.outputData
      splitCommand = [] 
      for idx in range(len(inp.COLUMNGROUPLIST)):
        if inp.COLUMNGROUPLIST[idx].selected:
          try:
            #print 'columnGroupType',inp.COLUMNGROUPLIST[idx].columnGroupType.__str__()
            columnGroupTypeIdx = ['Obs', 'Phs', 'MapCoeffs', 'FreeR'].index(inp.COLUMNGROUPLIST[idx].columnGroupType.__str__()) 
            cls = (CCP4XtalData.CObsDataFile,CCP4XtalData.CPhsDataFile,CCP4XtalData.CMapCoeffsDataFile, CCP4XtalData.CFreeRDataFile)[ columnGroupTypeIdx ]
          except:
            print('Unable to determine data type',inp.COLUMNGROUPLIST[idx].columnGroupType.__str__())
            cls = None
          if cls is not None:
            # This is a kludge to get a list of objects of mixed CMiniMtzDataFile sub-class
            out.MINIMTZOUTLIST.__dict__['_value'].append(cls(fullPath=os.path.join(self.workDirectory,inp.COLUMNGROUPLIST[idx].dataset.__str__()+'_'+ \
                                               inp.COLUMNGROUPLIST[idx].columnListStr(withTypes=False,splitter='_')+'.mtz') ))
            f = out.MINIMTZOUTLIST[-1].__str__()
            out.MINIMTZOUTLIST[-1].set(os.path.join(os.path.dirname(f),os.path.basename(f).replace("(","").replace(")","")))
            out.MINIMTZOUTLIST[-1].contentFlag.set(inp.COLUMNGROUPLIST[idx].contentFlag.__int__())
            out.MINIMTZOUTLIST[-1].annotation.set(os.path.splitext(str(inp.HKLIN.baseName))[0]+'/'+str(inp.COLUMNGROUPLIST[idx].dataset))
            print('out.MINIMTZOUTLIST',out.MINIMTZOUTLIST[-1].get())
            splitCommand.append( [ out.MINIMTZOUTLIST[-1].__str__(), inp.COLUMNGROUPLIST[idx].columnListStr(withTypes=False,splitter=','),           
                      out.MINIMTZOUTLIST[-1].columnNames(ifString=True) ] )

      if inp.USERCOLUMNGROUP.columnGroupType.isSet():
        try:
          columnGroupTypeIdx = ['Obs', 'Phs', 'MapCoeffs', 'FreeR'].index(inp.USERCOLUMNGROUP.columnGroupType.__str__())
          cls = (CCP4XtalData.CObsDataFile,CCP4XtalData.CPhsDataFile,CCP4XtalData.CMapCoeffsDataFile, CCP4XtalData.CFreeRDataFile)[ columnGroupTypeIdx ]
        except:
          print('Unable to determine data type',inp.USERCOLUMNGROUP.columnGroupType.__str__())
          cls = None
        if cls is not None:
          out.MINIMTZOUTLIST.__dict__['_value'].append(cls(fullPath=os.path.join(self.workDirectory,inp.USERCOLUMNGROUP.dataset.__str__()+'_'+ \
                                               inp.USERCOLUMNGROUP.columnListStr(withTypes=False,splitter='_')+'.mtz') ))
          out.MINIMTZOUTLIST[-1].contentFlag.set(inp.USERCOLUMNGROUP.contentFlag.__int__())
          out.MINIMTZOUTLIST[-1].annotation.set(str(inp.USERCOLUMNGROUP.dataset)+'_'+inp.USERCOLUMNGROUP.columnListStr(withTypes=False,splitter='_'))
          splitCommand.append( [ out.MINIMTZOUTLIST[-1].__str__(), inp.USERCOLUMNGROUP.columnListStr(withTypes=False,splitter=','),         
                      out.MINIMTZOUTLIST[-1].columnNames(ifString=True) ] )
      print('splitMtz.process splitCommand',splitCommand) ; sys.stdout.flush()
      status = self.splitMtz(inp.HKLIN.fullPath.__str__(),splitCommand)
      print('splitMtz.process status',status) ; sys.stdout.flush()
      
      return status
