from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript


class mrbump_basic(CPluginScript):

    TASKTITLE='MrBUMP Basic'
    TASKNAME = 'mrbump_basic'
    TASKCOMMAND = 'mrbump'
    MAINTAINER = 'ronan.keegan@stfc.ac.uk'

    def validity(self):
        """Override to filter out ENSEMBLES list length error.

        ENSEMBLES is defined in .def.xml with listMinLength=1 but is not used
        by MrBUMP - it's a vestigial element from when phaser was called directly.
        MrBUMP finds its own search models via sequence searches.
        """
        # Get parent validation
        error = super(mrbump_basic, self).validity()

        # Filter out the ENSEMBLES minimum length error (code 101)
        filtered = CCP4ErrorHandling.CErrorReport()
        for err in error.getErrors():
            # Skip error code 101 (min list length) for ENSEMBLES
            if err.get('code') == 101 and 'ENSEMBLES' in err.get('name', ''):
                continue
            filtered.append(
                klass=err.get('class', ''),
                code=err.get('code', 0),
                details=err.get('details', ''),
                name=err.get('name', ''),
                severity=err.get('severity', 0)
            )

        return filtered

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData
        error = None
        self.hklin = None
        dataObjects = []
        #Append Observation with representation dependent on whether we are detwining on Is or not
        dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):

      inp = self.container.inputData
      par = self.container.controlParameters
      mod = self.container.modelParameters
      gui = self.container.guiParameters
      out = self.container.outputData

      # Set the max number of processors
      import multiprocessing
      import os

      from ccp4i2.core import CCP4Utils
      MAXPROC=multiprocessing.cpu_count()  

      keyin = "MAPROGRAM clustalw2\n" 
      keyin += "MRPROGRAM phaser\n" 
      keyin += "DOFASTA False\n" 
      if self.container.controlParameters.LOCALONLY.isSet():
          if self.container.controlParameters.LOCALONLY:
              keyin += "DOPHMMER False\n" 
              keyin += "GESE False\n" 
          else:
              keyin += "DOPHMMER True\n" 
              keyin += "GESE True\n" 
      else:
          keyin += "DOPHMMER True\n" 
          keyin += "GESE True\n" 
      keyin += "PICKLE False\n" 
      keyin += "MDLU False\n" 
      keyin += "MDLD False\n" 
      keyin += "MDLC True\n" 
      keyin += "MDLM False\n" 
      keyin += "MDLP False\n" 
      keyin += "MDLS False\n" 
      keyin += "PHAQ True\n" 
      keyin += "USEPQS False\n" 
      if self.container.controlParameters.NCYC.isSet():
          keyin += "NCYC %s\n" % str(self.container.controlParameters.NCYC) 
      else:
          keyin += "NCYC 100\n" 
      keyin += "TRYALL True\n" 
      keyin += "UPDATE False\n" 
      if self.container.modelParameters.MRMAX:
          keyin += "MRNUM %s\n" % str(self.container.modelParameters.MRMAX)
      if self.container.controlParameters.PJOBS:
          if self.container.controlParameters.PJOBS > MAXPROC:
              keyin += "PJOBS %s\n" % str(MAXPROC) 
          else:
              keyin += "PJOBS %s\n" % str(self.container.controlParameters.PJOBS) 

      if self.container.modelParameters.SEARCH_PDB:
          if self.container.modelParameters.REDUNDANCYLEVEL:
              if str(self.container.modelParameters.REDUNDANCYLEVEL) == '110':
                 keyin += "RLEVEL ALL\n" 
              elif str(self.container.modelParameters.REDUNDANCYLEVEL) == '100':
                 keyin += "RLEVEL 100\n" 
              elif str(self.container.modelParameters.REDUNDANCYLEVEL) == '95':
                 keyin += "RLEVEL 95\n" 
              elif str(self.container.modelParameters.REDUNDANCYLEVEL) == '90':
                 keyin += "RLEVEL 90\n" 
              elif str(self.container.modelParameters.REDUNDANCYLEVEL) == '70':
                 keyin += "RLEVEL 70\n" 
              elif str(self.container.modelParameters.REDUNDANCYLEVEL) == '50':
                 keyin += "RLEVEL 50\n" 
      else:
          keyin += "RLEVEL 0\n" 

      if self.container.modelParameters.SEARCH_AFDB:
          if self.container.modelParameters.AFDBLEVEL is not None:
              if str(self.container.modelParameters.AFDBLEVEL) == '90':
                 keyin += "AFLEVEL 90\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '80':
                 keyin += "AFLEVEL 80\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '70':
                 keyin += "AFLEVEL 70\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '60':
                 keyin += "AFLEVEL 60\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '50':
                 keyin += "AFLEVEL 50\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '40':
                 keyin += "AFLEVEL 40\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '30':
                 keyin += "AFLEVEL 30\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '20':
                 keyin += "AFLEVEL 20\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '10':
                 keyin += "AFLEVEL 10\n" 
              elif str(self.container.modelParameters.AFDBLEVEL) == '0':
                 keyin += "AFLEVEL 0\n" 

      if self.container.modelParameters.PDBLOCAL.isSet():
          keyin += "PDBLOCAL %s\n" % str(self.container.modelParameters.PDBLOCAL)
      if self.container.modelParameters.HHPREDIN.isSet():
          if self.container.controlParameters.LOCALONLY.isSet():
              if self.container.controlParameters.LOCALONLY:
                  keyin += "DOHHPRED False\n" 
              else:
                  keyin += "HHRFILE %s\n" % str(self.container.modelParameters.HHPREDIN)
                  keyin += "DOHHPRED True\n" 
          else:
              keyin += "HHRFILE %s\n" % str(self.container.modelParameters.HHPREDIN)
              keyin += "DOHHPRED True\n" 
      else:
          keyin += "DOHHPRED False\n" 

      if len(self.container.inputData.SELECTEDCHAINS)>0:
          includeLine = "INCLUDE"
          for chain in self.container.inputData.SELECTEDCHAINS:
              includeLine += " %s" % str(chain)
          includeLine += "\n"
          keyin += includeLine

      for iCoordSet, xyzin in enumerate(self.container.modelParameters.XYZIN_LIST):
          if not xyzin.isSelectionSet():
              #self.appendCommandLine(['input.model='+xyzin.fullPath.__str__()])
              keyin += "LOCALFILE %s\n" % xyzin.fullPath.__str__()
          else:
              xyzin.selection.text='{'+xyzin.selection.__str__()+'} and {(ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR)}'
              inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected_'+iCoordSet.__str__()+'.pdb'))
              xyzin.getSelectedAtomsPdbFile(inputCoordPath)
              #self.appendCommandLine(['input.model='+inputCoordPath])
              keyin += "LOCALFILE %s\n" % inputCoordPath
    
      if self.container.controlParameters.BUCC.isSet():
          if self.container.controlParameters.BUCC:
              keyin += "BUCC True\n" 
          else:
              keyin += "BUCC False\n" 
      else:
          keyin += "BUCC False\n" 
      keyin += "END\n" 

      keyfile=os.path.join(self.getWorkDirectory(), "keywords.txt")
      kf=open(keyfile, "w")
      kf.write(keyin)
      kf.close()

      self.appendCommandLine( [ 'HKLIN', self.hklin ] )
      seqFile = os.path.join(self.workDirectory,'SEQIN.fasta')
      inp.ASUIN.writeFasta(fileName=seqFile)
      self.appendCommandLine( [ 'SEQIN', seqFile ] )
      self.appendCommandLine( [ 'KEYIN', str( keyfile ) ] )

      self.appendCommandLine( [ 'XMLOUT', str( self.makeFileName( 'PROGRAMXML' ) ) ] )
      self.appendCommandScript( "" )

      return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):
        import os
        import shutil

        xyzout = os.path.join(self.getWorkDirectory(), "output_mrbump_1.pdb")
        if os.path.exists(xyzout):
            self.container.outputData.XYZOUT=xyzout
        
        hklout = os.path.join(self.getWorkDirectory(), "output_mrbump_1.mtz")
        if os.path.exists(hklout):
            self.container.outputData.HKLOUT=hklout

        import os

        from ccp4i2.core import CCP4File, CCP4XtalData

        # Need to set the expected content flag  for phases data
        self.container.outputData.XYZOUT.annotation = 'Model from MrBump refinement'
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from MrBump refinement'
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from MrBump refinement'

        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        
        outputFiles = ['FPHIOUT','DIFFPHIOUT']
        outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']
        error = self.splitHklout(outputFiles,outputColumns,infile=hklout)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED
