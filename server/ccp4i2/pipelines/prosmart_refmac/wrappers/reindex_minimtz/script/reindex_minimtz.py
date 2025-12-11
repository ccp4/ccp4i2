"""
     refmac.py: CCP4 GUI Project
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

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4ErrorHandling


class reindex_minimtz(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'reindex_minimtz'
    TASKNAME = 'reindex_minimtz'
    TASKCOMMAND = 'reindex'
    TASKVERSION= 0.0
    WHATNEXT = []
    PROGRAMHELP = [ 'reindex','reindexing' ]

    def processInputFiles(self):
      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      from ccp4i2.core import CCP4XtalData
      # Need to set the expected content flag  for phases data
      self.container.outputData.HKLOUT.annotation = 'Reindexed reflections'

      xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
      xmlfile = open( xmlout, "w" )
      xmlfile.write( '<?xml version="1.0" encoding="ASCII" standalone="yes"?>\n<ReindexMiniMTZResult/>')
      xmlfile.close()

      return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):
      import os
      # make refmac command script   
      self.appendCommandLine(['HKLIN',self.container.inputData.HKLIN])
      self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT])

      if self.container.controlParameters.OPERATION.isSet():
          self.appendCommandScript("reindex %s,%s,%s"%(str(self.container.controlParameters.OPERATION.h),str(self.container.controlParameters.OPERATION.k),str(self.container.controlParameters.OPERATION.l)))
      if self.container.controlParameters.NEWSPACEGROUP:
          self.appendCommandScript("symmetry \'%s\'"%str(self.container.controlParameters.NEWSPACEGROUP))
      self.appendCommandScript('END')
      
      return CPluginScript.SUCCEEDED

