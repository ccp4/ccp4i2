from __future__ import print_function

"""
     buccaneer.py: CCP4 GUI Project
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

import os
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4ErrorHandling


class buccaneer_mr(CPluginScript):

    TASKMODULE          = 'wrappers'      # Where this plugin will appear on the gui  
    TASKTITLE           = 'Buccaneer_mr'  # A short title for gui menu
    TASKNAME            = 'buccaneer_mr'  # Task name - should be same as class name
    TASKCOMMAND         = 'cbuccaneer'    # The command to execute, should be reachable
    TASKVERSION         = 0.0             # Version of this plugin
    MAINTAINER = 'jon.agirre@york.ac.uk'

    def setProgramVersion(self):
        print('buccaneer.getProgramVersion')
        return CPluginScript.setProgramVersion(self,'Buccaneer')


    def processInputFiles(self):
      from ccp4i2.core import CCP4XtalData
      #print 'taskMakeHklin F_SIGF',self.container.inputData.F_SIGF,type(self.container.inputData.F_SIGF),self.container.inputData.F_SIGF.contentFlag
      if self.container.inputData.FWT_PHWT_IN.isSet():
        if self.container.inputData.FREERFLAG.isSet():
          self.hklin,columns,error = self.makeHklInput([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'ABCD', 'FREERFLAG', 'FWT_PHWT_IN' ])
          print('FREERFLAG is set, so joining all data objects')
        else :
          print('FREERFLAG is not set, so joining the rest of the data objects')
          self.hklin,columns,error = self.makeHklInput([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'ABCD', 'FWT_PHWT_IN' ])

        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
          print('ERROR creating input HKLIN with FWT_PHWT_IN')
          print(error.report())
          return CPluginScript.FAILED
      else:
        if self.container.inputData.FREERFLAG.isSet():
          self.hklin,columns,error = self.makeHklInput([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD','FREERFLAG'])
          print('FREERFLAG is set, so joining all data objects')
        else :
          print('FREERFLAG is not set, so joining the rest of the data objects')
          self.hklin,columns,error = self.makeHklInput([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD' ])
        
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
          print('ERROR creating input HKLIN')
          print(error.report())
          return CPluginScript.FAILED

      self.seqin = os.path.join(self.workDirectory,'seqin.fasta')
      self.container.inputData.ASUIN.writeFasta(self.seqin, polymerTypes=["PROTEIN"])
        
      self.refHklin = None
      conPars = self.container.controlParameters
      if conPars.F_SIGF_REF.isSet() and conPars.ABCD_REF.isSet() and conPars.XYZIN_REF.isSet() and \
         conPars.F_SIGF_REF.exists() and conPars.ABCD_REF.exists() and conPars.XYZIN_REF.exists():
        self.refHklin,error = self.makeHklin([['F_SIGF_REF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD_REF'])
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
          self.refHklin = None
      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      self.container.outputData.XYZOUT.annotation = 'Model built with Buccaneer'
      self.container.outputData.XYZOUT.subType=1
#      self.setPerformanceData()

      return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
      from ccp4i2.core import CCP4XtalData
   
      self.appendCommandLine(['-stdin'])

      # INPUT DATA
      self.appendCommandScript("mtzin "+self.hklin)
      self.appendCommandScript("colin-fo F_SIGF_F,F_SIGF_SIGF")
      if self.container.inputData.ABCD.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
        self.appendCommandScript("colin-hl ABCD_HLA,ABCD_HLB,ABCD_HLC,ABCD_HLD")
      else:
        self.appendCommandScript("colin-phifom ABCD_PHI,ABCD_FOM")

      if self.container.inputData.FREERFLAG.isSet( ) :
        self.appendCommandScript("colin-free FREERFLAG_FREER")


      self.appendCommandScript("seqin %s"%(self.seqin))

      if self.container.inputData.XYZIN_MODE:
        if self.container.inputData.XYZIN.isSet():
            self.appendCommandScript("pdbin %s"%(str(self.container.inputData.XYZIN.fullPath)))
            if self.container.controlParameters.KNOWN_STRUCTURE.isSet():
                self.appendCommandScript("known-structure %s"%(self.container.controlParameters.KNOWN_STRUCTURE.__str__()))
      
      if self.container.inputData.XYZIN_SEQ.isSet():
          self.appendCommandScript("pdbin-sequence-prior %s"%(str(self.container.inputData.XYZIN_SEQ.fullPath)))

      if self.container.inputData.FWT_PHWT_IN.isSet():
          self.appendCommandScript( "colin-fc FWT_PHWT_IN_F,FWT_PHWT_IN_PHI" )

      # OUTPUT DATA
      if self.container.outputData.XYZOUT.isSet():
          self.appendCommandScript("pdbout %s"%(str(self.container.outputData.XYZOUT.fullPath)))
      self.appendCommandScript("xmlout %s"%(self.makeFileName('PROGRAMXML')))


      # CONTROL PARAMETERS
      #if self.container.controlParameters.TITLE.isSet():
      #    self.appendCommandScript("title %s"%(str(self.container.controlParameters.TITLE)))
      if self.container.controlParameters.TIDY_ONLY:
          self.appendCommandScript("cycles 1")
          self.appendCommandScript("tidy")
      elif self.container.controlParameters.CYCLES.isSet():
          self.appendCommandScript("cycles %s"%(str(self.container.controlParameters.CYCLES)))
      if self.container.controlParameters.ANISOTROPY_CORRECTION:
          self.appendCommandScript("anisotropy-correction")
      if self.container.controlParameters.FAST:
          self.appendCommandScript("fast")
      if self.container.controlParameters.BUILD_SEMET:
          self.appendCommandScript("build-semet")
      if self.container.controlParameters.FIX_POSITION:
          self.appendCommandScript("fix-position")
#      if self.container.controlParameters.CORRELATION_MODE:
#          self.appendCommandScript("correlation-mode")
#      Kevin said that this option is handled automatically now

      if self.container.controlParameters.RESOLUTION.isSet():
          self.appendCommandScript("resolution %s"%(str(self.container.controlParameters.RESOLUTION)))
      else:
          self.appendCommandScript("resolution 2.0")  # I've added the default value here
      if self.container.controlParameters.SEQUENCE_RELIABILITY.isSet():
          self.appendCommandScript("sequence-reliability %s"%(str(self.container.controlParameters.SEQUENCE_RELIABILITY)))
      if self.container.controlParameters.NEW_RESIDUE_NAME.isSet():
          self.appendCommandScript("new-residue-name %s"%(str(self.container.controlParameters.NEW_RESIDUE_NAME)))
      if self.container.controlParameters.JOBS.isSet():
          self.appendCommandScript("jobs %s"%(str(self.container.controlParameters.JOBS)))
      if self.container.controlParameters.VERBOSE.isSet():
          self.appendCommandScript("verbose %s"%(str(self.container.controlParameters.VERBOSE)))

      if self.container.controlParameters.MODEL_SIGMA.isSet() and not self.container.controlParameters.TIDY_ONLY:
          self.appendCommandScript( "model-filter" )
          self.appendCommandScript( "model-filter-sigma %s"%(str(self.container.controlParameters.MODEL_SIGMA)))
      
      # Now we're going to determine if we need to trigger MR mode in Buccaneer 
    
      # pbond and jonesie changed this to make it work like in i1
      if self.container.controlParameters.PHSIN_TYPE == "mr":
          #self.appendCommandScript("mr-model")
          self.appendCommandScript("pdbin-mr %s"%str(self.container.inputData.MR_MODE_XYZIN.fullPath))

          if self.container.controlParameters.MR_MODE_SIGMA.isSet():
            self.appendCommandScript("mr-model-filter-sigma %s"%(str(self.container.controlParameters.MR_MODE_SIGMA)))

          
          if self.container.controlParameters.MR_MODE == "seed":
              self.appendCommandScript("mr-model-seed")
              self.appendCommandScript("mr-model-filter")
          else:
              if self.container.controlParameters.MR_MODE == "initialmodel":
                  self.appendCommandScript("mr-model-filter")
                
      # REFERENCE DATA
      if self.refHklin is not None:
        self.appendCommandScript("colin-ref-fo F_SIGF_REF_F,F_SIGF_REF_SIGF")
        if self.container.controlParameters.F_SIGF_REF.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
          self.appendCommandScript("colin-ref-hl ABCD_REF_HLA,ABCD_REF_HLB,ABCD_REF_HLC,ABCD_REF_HLD")
        else:
          self.appendCommandScript("colin-ref-phifom PHI,FOM")
        self.appendCommandScript("pdbin-ref %s"%(str(self.container.controlParameters.XYZIN_REF.fullPath)))

      return CPluginScript.SUCCEEDED



