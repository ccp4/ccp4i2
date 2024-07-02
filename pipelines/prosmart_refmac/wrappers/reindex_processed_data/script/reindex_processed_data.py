from __future__ import print_function

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
from PySide6 import QtCore
from core.CCP4PluginScript import CPluginScript
import core.CCP4ErrorHandling


class reindex_processed_data(CPluginScript):
    
    TASKMODULE = 'expt_data_util'
    TASKTITLE = 'reindex_processed_data'
    TASKNAME = 'reindex_processed_data'
    TASKCOMMAND = 'reindex'
    TASKVERSION= 0.0
    WHATNEXT = []

    ERROR_CODES = { 201 : { 'description' : 'Expected output file not found' }
                    }
    
    def processInputFiles(self):
        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        from core import CCP4XtalData
        # Need to set the expected content flag  for phases data
        self.container.outputData.HKLOUT.annotation = 'Reindexed reflections'
        
        xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
        xmlfile = open( xmlout, "w" )
        xmlfile.write( '<?xml version="1.0" encoding="ASCII" standalone="yes"?>\n<ReindexMiniMTZResult/>')
        xmlfile.close()
        
        return CPluginScript.SUCCEEDED
    
    def process(self):
        # Check all input files exist
        nonExFiles = self.checkInputData()
        if len(nonExFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
            return
        
        # Provide default output file names if necessary
        self.checkOutputData()
        
        # Create instance of reindex_minimtz class (and get it registered with the database)
        self.reindexObs = self.makePluginObject('reindex_minimtz')
        self.reindexObs.container.controlParameters = self.container.controlParameters
        self.reindexObs.container.inputData.HKLIN = self.container.inputData.OBSIN
        self.connectSignal(self.reindexObs,'finished',self.reindexObsFinished)
        self.reindexObs.process()
    
    
    @QtCore.Slot(dict)
    def reindexObsFinished(self, jobStatus):
        import os, shutil
        if jobStatus['finishStatus'] != CPluginScript.SUCCEEDED:
            self.reportStatus(CPluginScript.FAILED)
        try:
            shutil.copyfile(self.reindexObs.container.outputData.HKLOUT.__str__(),
                            self.container.outputData.OBSOUT.__str__())
        except:
            self.appendErrorReport(201,self.reindexObs.container.outputData.HKLOUT.__str__())
            self.reportStatus(CPluginScript.FAILED)
        
        self.container.outputData.OBSOUT.annotation = "Reindexed observations (%s)"%(str(self.container.controlParameters.NEWSPACEGROUP))
        print(self.container.inputData.FREERIN.__str__())
        if self.container.inputData.FREERIN.isSet() and self.container.inputData.FREERIN.exists():
          self.reindexFreer = self.makePluginObject('reindex_minimtz')
          self.reindexFreer.container.controlParameters = self.container.controlParameters
          self.reindexFreer.container.inputData.HKLIN = self.container.inputData.FREERIN
          self.connectSignal(self.reindexFreer,'finished',self.reindexFreerFinished)
          self.reindexFreer.process()
        else:
          self.postProcessWrapper(jobStatus)
          self.reportStatus(CPluginScript.SUCCEEDED)

    @QtCore.Slot(dict)
    def reindexFreerFinished(self, jobStatus):
        import os, shutil
        if jobStatus['finishStatus'] != CPluginScript.SUCCEEDED:
            self.reportStatus(CPluginScript.FAILED)
        try:
            shutil.copyfile(self.reindexFreer.container.outputData.HKLOUT.__str__(),
                            self.container.outputData.FREEROUT.__str__())
        except:
            self.appendErrorReport(201,self.reindexObs.container.outputData.HKLOUT.__str__())
            self.reportStatus(CPluginScript.FAILED)
        self.container.outputData.FREEROUT.annotation = "Reindexed FREE-R Set (%s)"%(str(self.container.controlParameters.NEWSPACEGROUP))
        
        from core import CCP4XtalData
        xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
        xmlfile = open( xmlout, "w" )
        xmlfile.write( '<?xml version="1.0" encoding="ASCII" standalone="yes"?>\n<ReindexMiniMTZResult/>')
        xmlfile.close()
        
        self.reportStatus(CPluginScript.SUCCEEDED)

