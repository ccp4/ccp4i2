"""
Copyright (C) 2010 University of York
"""

import os
import xml.etree.ElementTree as ET

from ....core import CCP4ErrorHandling
from ....core import CCP4Utils
from ....core import CCP4XtalData
from ....core.CCP4PluginScript import CPluginScript
from ....pimple import MGQTmatplotlib
from ....smartie import smartie


class parrot(CPluginScript):

    TASKMODULE = 'density_modification'                       # Where this plugin will appear on the gui
    TASKTITLE = 'Parrot'                                # A short title for gui menu
    TASKNAME = 'parrot'                                 # Task name - should be same as class name
    TASKCOMMAND = 'cparrot'                             # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    PERFORMANCECLASS = 'CExpPhasPerformance'             # Only FOM is relevent
    MAINTAINER = 'kevin.cowtan@york.ac.uk'


    def processInputFiles(self):

      cols = [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD']
      if self.container.inputData.FREERFLAG.isSet(): cols.append('FREERFLAG')
      if self.container.inputData.F_PHI.isSet(): cols.append('F_PHI')
      self.hklin, __, error = self.makeHklInput(cols)
      if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING: return CPluginScript.FAILED

      self.refHklin = None
      conPars = self.container.controlParameters
      if conPars.F_SIGF_REF.isSet() and conPars.ABCD_REF.isSet() and conPars.XYZIN_REF.isSet() and \
         conPars.F_SIGF_REF.exists() and conPars.ABCD_REF.exists() and conPars.XYZIN_REF.exists():
        self.refHklin, __, error = self.makeHklInput([['F_SIGF_REF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD_REF'],
                                                     extendOutputColnames=False, useInputColnames=False)
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
          self.refHklin = None

      self.seqin = os.path.join(self.workDirectory,'seqin.fasta')
      self.container.inputData.ASUIN.writeFasta(self.seqin)
      
      
      return CPluginScript.SUCCEEDED

    def  processOutputFiles(self):
      # Need to set the expected content flag  for phases data
      self.container.outputData.ABCDOUT.contentFlag.set(CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL)
      print('parrot.processOutputFiles jobNumberString',self.jobNumberString())
      self.container.outputData.ABCDOUT.annotation = self.jobNumberString() + ' Phases from density modification'
      self.container.outputData.FPHIOUT.contentFlag.set(1)
      self.container.outputData.FPHIOUT.subType.set(1)
      self.container.outputData.FPHIOUT.annotation = self.jobNumberString() + ' Map coefficients from density modification'
      
      # extend XML output
      rootNode = ET.parse(self.xmlout).getroot()
      smartieNode = ET.SubElement(rootNode,'SmartieGraphs')
      self.scrapeSmartieGraphs(smartieNode)
      CCP4Utils.writeXml(rootNode, self.xmlout)

      # performance data
      final_fom = float(rootNode.xpath('//ParrotResult/Final/MeanFOM')[0].text)
      self.container.outputData.PERFORMANCE.FOM = final_fom

      # error checking
      error = self.splitHklout(['FPHIOUT','ABCDOUT'],['parrot.F_phi.F,parrot.F_phi.phi','parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D'])
      if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED

    def scrapeSmartieGraphs(self, smartieNode):
        logfile = smartie.parselog(self.makeFileName('LOG'))
        for smartieTable in logfile.tables():
            if smartieTable.ngraphs() > 0:
                tableelement = self.xmlForSmartieTable(smartieTable, smartieNode)
        return

    def xmlForSmartieTable(self, table, parent):
        tableetree = MGQTmatplotlib.CCP4LogToEtree(table.rawtable())
        parent.append(tableetree)
        return tableetree

    def makeCommandAndScript(self):
     
      self.hklout = os.path.join(self.workDirectory,"hklout.mtz")
      self.appendCommandLine(['-stdin'])

      # INPUT DATA
      self.appendCommandScript("mtzin "+self.hklin)
      self.appendCommandScript("colin-fo F_SIGF_F,F_SIGF_SIGF")
      self.container.inputData.ABCD.setContentFlag()
      
      if self.container.inputData.ABCD.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
        self.appendCommandScript("colin-hl ABCD_HLA,ABCD_HLB,ABCD_HLC,ABCD_HLD")
      else:
        self.appendCommandScript("colin-phifom ABCD_PHI,ABCD_FOM")
      if self.container.inputData.FREERFLAG.isSet():
        self.appendCommandScript("colin-free FREERFLAG_FREER")
      if self.container.inputData.F_PHI.isSet():
        self.appendCommandScript("colin-fc F_PHI_F,F_PHI_PHI")

      self.appendCommandScript("seqin %s"%(self.seqin))

      if self.container.inputData.XYZIN_MODE=='ha' and self.container.inputData.XYZIN_HA.isSet():
          self.appendCommandScript("pdbin-ha %s"%(str(self.container.inputData.XYZIN_HA.fullPath)))
      elif self.container.inputData.XYZIN_MODE=='mr' and self.container.inputData.XYZIN_MR.isSet():
          self.appendCommandScript("pdbin-mr %s"%(str(self.container.inputData.XYZIN_MR.fullPath)))
      
      # OUTPUT DATA
      self.xmlout = self.makeFileName('PROGRAMXML')
      self.appendCommandScript("mtzout "+self.hklout)
      self.appendCommandScript("xmlout "+self.xmlout)

      # CONTROL PARAMETERS
      if self.container.controlParameters.CYCLES.isSet():
          self.appendCommandScript("cycles %s"%(str(self.container.controlParameters.CYCLES)))
      if self.container.controlParameters.ANISOTROPY_CORRECTION:
          self.appendCommandScript("anisotropy-correction")
      if self.container.controlParameters.RESOLUTION.isSet():
          self.appendCommandScript("resolution %s"%(str(self.container.controlParameters.RESOLUTION)))
      if self.container.controlParameters.SOLVENT_CONTENT.isSet():
          self.appendCommandScript("solvent-content %s"%(str(self.container.controlParameters.SOLVENT_CONTENT)))
      if self.container.controlParameters.NCS_MASK_FILTER_RADIUS.isSet():
          self.appendCommandScript("ncs-mask-filter-radius %s"%(str(self.container.controlParameters.NCS_MASK_FILTER_RADIUS)))
      if self.container.controlParameters.VERBOSE.isSet():
          self.appendCommandScript("verbose %s"%(str(self.container.controlParameters.VERBOSE)))

      # REFERENCE DATA
      if self.refHklin is not None:
        self.appendCommandScript("colin-ref-fo F,SIGF")
        if self.container.controlParameters.F_SIGF_REF.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
          self.appendCommandScript("colin-ref-hl HLA,HLB,HLC,HLD")
        else:
          self.appendCommandScript("colin-ref-phifom PHI,FOM")      
        self.appendCommandScript("pdbin-ref %s"%(str(self.container.controlParameters.XYZIN_REF.fullPath)))

      return CPluginScript.SUCCEEDED


def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]
