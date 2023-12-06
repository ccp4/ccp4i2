from __future__ import print_function

"""
     mtzdump.scripts.py: CCP4 GUI Project
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

     etc
"""

from core.CCP4PluginScript import CPluginScript
     
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
    

#=====================================================================================================
#=================================test suite=========================================================
#=====================================================================================================

import unittest

# unit testing asynchronous processes potential tricky but QProcess has option to wait for finished
 
class testMtzdump(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    from core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

  def testMtzdump(self):
    from core.CCP4Modules import QTAPPLICATION
    self.wrapper = mtzdump(parent=QTAPPLICATION(),name='test_mtzdump')
    self.wrapper.container.inputData.HKLIN.set({'project':'CCP4I2_TOP','baseName':'gere_nat.mtz','relPath':'test/data'})
    pid = self.wrapper.process()
    print(self.wrapper.container.outputData.CELL)
    if len(self.wrapper.errorReport)>0: self.wrapper.errorReport.report()
    self.assertEqual(self.wrapper.container.outputData.CELL.a,108.742,'Mtzdump output CELL wrong')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testMtzdump)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
