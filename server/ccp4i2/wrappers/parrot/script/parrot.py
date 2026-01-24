import os

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class parrot(CPluginScript):

    TASKMODULE = 'density_modification'                       # Where this plugin will appear on the gui
    TASKTITLE = 'Parrot'                                # A short title for gui menu
    TASKNAME = 'parrot'                                 # Task name - should be same as class name
    TASKCOMMAND = 'cparrot'                             # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    PERFORMANCECLASS = 'CExpPhasPerformance'             # Only FOM is relevent
    MAINTAINER = 'kevin.cowtan@york.ac.uk'


    def processInputFiles(self):

      cols = [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD']
      if self.container.inputData.FREERFLAG.isSet(): cols.append('FREERFLAG')
      if self.container.inputData.F_PHI.isSet(): cols.append('F_PHI')
      self.hklin, __, error = self.makeHklInput(cols)
      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING: return CPluginScript.FAILED

      self.refHklin = None
      conPars = self.container.controlParameters
      if conPars.F_SIGF_REF.isSet() and conPars.ABCD_REF.isSet() and conPars.XYZIN_REF.isSet() and \
         conPars.F_SIGF_REF.exists() and conPars.ABCD_REF.exists() and conPars.XYZIN_REF.exists():
        self.refHklin, __, error = self.makeHklInput([['F_SIGF_REF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD_REF'],
                                                     extendOutputColnames=False, useInputColnames=False)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
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
      from lxml import etree
      rootNode = etree.Element("ParrotResult")
      with open(self.xmlout,'r') as xmlFile:
        rootNode = etree.fromstring(xmlFile.read())
      smartieNode = etree.SubElement(rootNode,'SmartieGraphs')
      self.scrapeSmartieGraphs(smartieNode)
      with open(self.xmlout,'w') as xmlFile:
        CCP4Utils.writeXML(xmlFile,etree.tostring(rootNode,pretty_print=True))

      # performance data
      final_fom = float(rootNode.xpath('//ParrotResult/Final/MeanFOM')[0].text)
      self.container.outputData.PERFORMANCE.FOM = final_fom

      # error checking
      error = self.splitHklout(['FPHIOUT','ABCDOUT'],['parrot.F_phi.F,parrot.F_phi.phi','parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D'])
      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED

    def scrapeSmartieGraphs(self, smartieNode):
        from lxml import etree

        from ccp4i2.smartie import smartie
        
        logfile = smartie.parselog(self.makeFileName('LOG'))
        for smartieTable in logfile.tables():
            if smartieTable.ngraphs() > 0:
                tableelement = self.xmlForSmartieTable(smartieTable, smartieNode)
        return

    def xmlForSmartieTable(self, table, parent):
        from ccp4i2.pimple.logtable import CCP4LogToEtree
        tableetree = CCP4LogToEtree(table.rawtable())
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

     
# ----------------------------------------------------------------------
# Function to return list of names of exportable MTZ(s)
'''
def exportJobFile(jobId=None,mode=None):
    import os
    from ccp4i2.core import CCP4Modules
    from ccp4i2.core import CCP4XtalData

    # Devise name for the merged file and check if it has already been created
    jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
    exportFile = os.path.join(jobDir,'exportMtz.mtz')
    if os.path.exists(exportFile): return exportFile

    # Get the source reflection data either from aimless or an imported file
    # getSourceReflectionFile() returns a dict with elements: fileName, source, jobNumber
    reflnInfo = CCP4Modules.PROJECTSMANAGER().getSourceReflectionFile(jobId = jobId, jobParamName='F_SIGF')
    print 'parrot.exportJobFile getSourceReflectionFile',reflnInfo
    if reflnInfo['fileName'] is None: return None

    # Query database for filenames and job info for the input and ouptput phase objects
    db = CCP4Modules.PROJECTSMANAGER().db()
    abcdInfo = db.getJobFilesInfo(jobId=jobId,jobParamName='ABCD',input=True)
    jN2 = abcdInfo[0]['jobnumber']+'_'+abcdInfo[0]['taskname']
    abcdOutInfo = db.getJobFilesInfo(jobId=jobId,jobParamName='ABCDOUT')
    jN3 = abcdOutInfo[0]['jobnumber']+'_'+abcdOutInfo[0]['taskname']
    print 'parrot.exportJobFile abcdInfo',abcdInfo
    print 'parrot.exportJobFile abcdOutInfo',abcdOutInfo

    print 'parrot.exportJobFile  runCad:',exportFile,abcdInfo[0]['fullPath'],  abcdOutInfo[0]['fullPath']
    # Custom input to runCad to set output column names
    comLines = [ ]

    # Beware the input phases not necessarilly HL
    if abcdInfo[0].get('fileContent',None) is not None:
      fC = abcdInfo[0]['fileContent']
    else:
      p = CCP4XtalData.CPhsDataFile(abcdInfo[0]['fullPath'])
      p.setContentFlag()
      fC = int(p.contentFlag)

    # Use CAD LABOUT line to set column labels in export file
    if fC == 1:
     comLines.append ( 'LABOUT FILENUMBER 2 E1=HLA_'+jN2+'  E2=HLB_'+jN2+' E3=HLC_'+jN2+' E4=HLD_'+jN2 )
    else:
      comLines.append ( 'LABOUT FILENUMBER 2 E1=PHI_'+jN2+'  E2=FOM_'+jN2 )
    # Parrot output is always HLs?
    comLines.append ( 'LABOUT FILENUMBER 3 E1=HLA_'+jN3+'  E2=HLB_'+jN3+' E3=HLC_'+jN3+' E4=HLD_'+jN3 )

    #  Create an CMtzDataFile object and initialise with the refln data file
    m = CCP4XtalData.CMtzDataFile(reflnInfo['fileName'])
    #print m.runCad.__doc__   #Print out docs for the function
    outfile,err = m.runCad(exportFile,[ abcdInfo[0]['fullPath'],  abcdOutInfo[0]['fullPath'] ] ,comLines )
    print 'aimless_pipe.exportJobFile',outfile,err.report()
    return   outfile
'''                                                
 
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]
