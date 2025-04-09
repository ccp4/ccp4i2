"""
Copyright (C) 2010 University of York
"""

from ....core.CCP4PluginScript import CPluginScript


class mtzdump(CPluginScript):

    TASKMODULE = 'test'
    MAINTAINER = 'liz.potterton@york.ac.uk'
    TASKTITLE = 'MTZDump'
    TASKNAME = 'mtzdump'
    TASKCOMMAND = 'mtzdump'
    TASKVERSION= 0.0
    COMLINETEMPLATE = '1 HKLIN $HKLIN'
    COMTEMPLATE = '''$HEADER HEADER
1 NREF -1
1 END
'''

    ERROR_CODES = { 101 : { 'description' : 'Log file results not parsable' }
                    }

    def processOutputFiles(self):
      #print 'mtzdump.processOutputFiles'

      # get data fom log file and load into outputData
      data = self.parseMtzdumpLog()
      if len(data['cell'])>0:
        rv = 0
        try:
          self.container.outputData.CELL.set( { 'a' : data['cell'][0],
                                         'b' : data['cell'][1],
                                         'c' : data['cell'][2],
                                         'alpha' : data['cell'][3],
                                         'beta' : data['cell'][4],
                                         'gamma' : data['cell'][5] } )
                                         #'spaceGroup' : data['spaceGroup'] } )
        except:
          rv = 1
      else:
        rv = 1
      if rv>0:
        self.appendErrorReport(101)
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED
      
     
    def parseMtzdumpLog(self):
      # Extract data from log file 
      # Code taken from EDNA example
      logText = self.logFileText()
      pyListLogLines = logText.split("\n")
      cell = []
      listOfColumns = []
      column_name_list = []
      column_type_list = []
      pyStrSpaceGroupName = ''
      iSpaceGroupNumber = -1
      lowerResolutionLimit = None
      upperResolutionLimit = None
      
      for j, pyStrLine in enumerate(pyListLogLines):
            if "* Dataset ID, project/crystal/dataset names, cell dimensions, wavelength:" in pyStrLine:
               try:
                 cell = list(map(float, pyListLogLines[j + 5].split()))
               except:
                 pass
            if " * Space group = " in pyStrLine:
                pyStrSpaceGroupName = pyStrLine.split("'")[1].strip()
                iSpaceGroupNumber = int(pyStrLine.replace("(", " ").replace(")", " ").split()[-1])
            if "*  Resolution Range" in pyStrLine:
                lowerResolutionLimit = float(((pyListLogLines[j + 2].split("(")[1]).split())[0])
                upperResolutionLimit = float(((pyListLogLines[j + 2].split("(")[1]).split())[2])
            if "* Column Labels" in pyStrLine:
                column_name_list = pyListLogLines[j + 2].split()
            if "* Column Types" in pyStrLine:
                column_type_list = pyListLogLines[j + 2].split()

      for j, column_name in enumerate(column_name_list):
            column_type = column_type_list[j]
            listOfColumns.append ( {'name':column_name, 'value': column_type} )
      return { 'cell' : cell,
               'spaceGroup' : pyStrSpaceGroupName,
               'spaceGroupNumber' : iSpaceGroupNumber,
               'lowerResolutionLimit' : lowerResolutionLimit,
               'upperResolutionLimit' : upperResolutionLimit,
               'listOfColumns' : listOfColumns }
