from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4ErrorHandling


class pyphaser_mr(CPluginScript):

    TASKTITLE = 'MR using Phaser (pythonic)'
    TASKNAME = 'pyphaser_mr'

    def process(self):

       unsetData = self.checkInputData()
       if len(unsetData)>0:
         self.reportStatus(CPluginScript.FAILED)
         return

       # No output files, so skip checkOutputData

       status = self.processInputFiles()
       if status == CPluginScript.FAILED:
         self.reportStatus(CPluginScript.FAILED)
         return

       import phaser
       inputData = self.container.inputData

       inp = phaser.InputMR_DAT()
       if self.container.controlParameters.TITLE.isSet():
         inp.setTITL(str(self.container.controlParameters.TITLE))
       inp.setHKLI(str(self.hklin))
       inp.setLABI('F','SIGF')
       inp.setMUTE(True)
       self.initialResults = phaser.runMR_DAT(inp)
       self.xmlfile = None
       self.outputInitialXml()

       if self.container.controlParameters.MODE == "MR_AUTO":
          inp = phaser.InputMR_AUTO()
       elif self.container.controlParameters.MODE == "MR_FRF":
          inp = phaser.InputMR_FRF()
       else:
          print('Chosen mode not supported')
          self.reportStatus(CPluginScript.FAILED)
          return

       inp.setSPAC_HALL(self.initialResults.getSpaceGroupHall())
       inp.setCELL6(self.initialResults.getUnitCell())
       inp.setREFL(self.initialResults.getMiller(),self.initialResults.getFobs(),self.initialResults.getSigFobs())

       controlParameters = self.container.controlParameters
       if controlParameters.COMP_BY == 'DEFAULT':
           #Default is 50% solvent ?
           inp.setCOMP_BY("AVERAGE")
           pass
       elif controlParameters.COMP_BY == 'MW':
           if controlParameters.ASU_PROTEIN_MW.isSet():
               inp.addCOMP_PROT_MW_NUM(float(controlParameters.ASU_PROTEIN_MW), 1.0)
           if controlParameters.ASU_NUCLEICACID_MW.isSet():
               inp.addCOMP_NUCL_MW_NUM(float(controlParameters.ASU_NUCLEICACID_MW), 1.0)
       elif controlParameters.COMP_BY == 'ASU':
           for i in range(len(inputData.ASU_COMPONENTS)):
               inp.addCOMP_PROT_SEQ_NUM(str(inputData.ASU_COMPONENTS[i].seqFile),float(inputData.ASU_COMPONENTS[i].numberOfCopies)) 

       # search models
       for i in range(len(inputData.ENSEMBLES)):
           for j in range(len(inputData.ENSEMBLES[i].pdbItemList)):
              if inputData.ENSEMBLES[i].pdbItemList[j].identity_to_target.isSet():
                inp.addENSE_PDB_ID(str(inputData.ENSEMBLES[i].label),str(inputData.ENSEMBLES[i].pdbItemList[j].structure),float(inputData.ENSEMBLES[i].pdbItemList[j].identity_to_target))
              elif inputData.ENSEMBLES[i].pdbItemList[j].rms_to_target.isSet():
                inp.addENSE_PDB_RMS(str(inputData.ENSEMBLES[i].label),str(inputData.ENSEMBLES[i].pdbItemList[j].structure),float(inputData.ENSEMBLES[i].pdbItemList[j].rms_to_target))

       #define search
       if self.container.controlParameters.SEARCHMODE == 'multidomain':
           for i in range(len(inputData.ENSEMBLES)):
              if bool(inputData.ENSEMBLES[i].use):
                inp.addSEAR_ENSE_NUM(str(inputData.ENSEMBLES[i].label),int(inputData.ENSEMBLES[i].number))

       elif self.container.controlParameters.SEARCHMODE == 'alternatives':
           search_array = []
           for i in range(1,len(inputData.ENSEMBLES)):
              if bool(inputData.ENSEMBLES[i].use):
                search_array.append(str(inputData.ENSEMBLES[i].label))
           inp.addSEAR_ENSE_OR_ENSE_NUMB(search_array,int(inputData.ENSEMBLES[0].number))

       if self.container.controlParameters.SGALT_SELECT.isSet():
           inp.setSGAL_SELE(str(self.container.controlParameters.SGALT_SELECT))
           if self.container.controlParameters.SGALT_SELECT == 'LIST' and self.container.controlParameters.SGALT_TEST.isSet():
              for i in range(len(self.container.controlParameters.SGALT_TEST)):
                inp.addSGAL_TEST(str(self.container.controlParameters.SGALT_TEST[i]))

       if self.container.controlParameters.PERMUTATIONS.isSet():
           if self.container.controlParameters.PERMUTATIONS:
               inp.setPERM(True)
       if self.container.controlParameters.NUM_SOL_OUT.isSet():
           inp.setTOPF(float(self.container.controlParameters.NUM_SOL_OUT))
       if self.container.controlParameters.RESOLUTION_HIGH.isSet():
           inp.setRESO_HIGH(float(self.container.controlParameters.RESOLUTION_HIGH))
       if self.container.controlParameters.NJOBS.isSet():
           inp.setJOBS(int(self.container.controlParameters.NJOBS))
       if self.container.controlParameters.PACK_CUTOFF.isSet():
           #These keywords depend on phaser version. I am assuming 2.5.2
           inp.setPACK_SELE(str(self.container.controlParameters.PACK_SELECT))
           inp.setPACK_CUTO(float(self.container.controlParameters.PACK_CUTOFF))
       if self.container.controlParameters.PEAKS_ROT_CUTOFF.isSet():
           inp.setPEAK_ROTA_CUTO(float(self.container.controlParameters.PEAKS_ROT_CUTOFF))

       # Set root to correct working directory
       import os
       if self.container.controlParameters.ROOT.isSet():
           inp.setROOT(os.path.join(self.getWorkDirectory(),str(self.container.controlParameters.ROOT)))
       else:
           inp.setROOT(os.path.join(self.getWorkDirectory(),'PHASER'))

       #log = open(os.path.join(self.getWorkDirectory(),"phaser.log"), "w")
       if self.container.controlParameters.MODE == "MR_AUTO":
           self.results = phaser.runMR_AUTO(inp)
       elif self.container.controlParameters.MODE == "MR_FRF":
           self.results = phaser.runMR_FRF(inp)
       #close(log)
       print(self.results.summary())

       if self.results.foundSolutions() :
           print("Phaser has found MR solutions")
           print("Top LLG = %f" % self.results.getTopLLG())
           print("Top PDB file = %s" % self.results.getTopPdbFile())
       else:
           print("Phaser has not found any MR solutions")

       self.processOutputFiles()

       # Needed!
       self.reportStatus(CPluginScript.SUCCEEDED)

    def outputInitialXml(self):

       xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
       # open for write, and with large enough buffer
       self.xmlfile = open( xmlout, "w", 4096 )
       self.xmlfile.write( '<?xml version="1.0" encoding="ASCII" standalone="yes"?>\n' )
       self.xmlfile.write( "<PhaserMrResults>\n" )
       self.xmlfile.write( " <PhaserMrStats>\n" )
       self.xmlfile.write( "  <SPG_INITIAL>"+str(self.initialResults.getSpaceGroupName())+"</SPG_INITIAL>\n" )
       self.xmlfile.write( "  <HALL>"+str(self.initialResults.getSpaceGroupHall())+"</HALL>\n" )
       self.xmlfile.write( "  <CELL>\n" )
       self.xmlfile.write( "    <CELL_A>"+str(self.initialResults.getUnitCell()[0])+"</CELL_A>\n" )
       self.xmlfile.write( "    <CELL_B>"+str(self.initialResults.getUnitCell()[1])+"</CELL_B>\n" )
       self.xmlfile.write( "    <CELL_C>"+str(self.initialResults.getUnitCell()[2])+"</CELL_C>\n" )
       self.xmlfile.write( "    <CELL_ALPHA>"+str(self.initialResults.getUnitCell()[3])+"</CELL_ALPHA>\n" )
       self.xmlfile.write( "    <CELL_BETA>"+str(self.initialResults.getUnitCell()[4])+"</CELL_BETA>\n" )
       self.xmlfile.write( "    <CELL_GAMMA>"+str(self.initialResults.getUnitCell()[5])+"</CELL_GAMMA>\n" )
       self.xmlfile.write( "  </CELL>\n" )
       self.xmlfile.write( " </PhaserMrStats>\n" )
       self.xmlfile.flush()

    def processInputFiles(self):
      from ccp4i2.core import CCP4XtalData

      self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
        for report in error._reports:
          if report['code'] == 32:
            report['details'] = 'Observed data has no F/SIGF columns, required by Phaser. Check file import.'
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):

      import os,shutil

      if self.container.controlParameters.NUM_SOL_OUT.isSet():
        num_sol = self.container.controlParameters.NUM_SOL_OUT.isSet()
      else:
        num_sol = 1

      for i in range(1,num_sol+1):
        xyzout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".pdb")
        if os.path.exists(xyzout):
          self.container.outputData.XYZOUT.append(xyzout)
          self.container.outputData.XYZOUT[-1].annotation = 'Positioned coordinates for solution '+str(i)
        hklout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".mtz")
        if os.path.exists(hklout):
          self.container.outputData.HKLOUT.append(hklout)

      self.splitHkloutList(miniMtzsOut=['MAPOUT','DIFMAPOUT'],programColumnNames=['FWT,PHWT','DELFWT,PHDELWT'],
                           outputBaseName=['MAPOUT','DIFMAPOUT'],infileList=self.container.outputData.HKLOUT)
      for indx in range(len(self.container.outputData.MAPOUT)):
        self.container.outputData.MAPOUT[indx].annotation = 'Map for solution '+str(indx+1)
        self.container.outputData.DIFMAPOUT[indx].annotation = 'Difference map for solution '+str(indx+1)
        
      ## doesn't seem to be written now?
      phaser_solfile = os.path.join(self.getWorkDirectory(), "PHASER.sol")

      #print dir(self.results)

      ###dump XML from result object

      if (self.xmlfile is None):
         xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
         self.xmlfile = open( xmlout, "w" )
         self.xmlfile.write( '<?xml version="1.0" encoding="ASCII" standalone="yes"?>\n' )

      self.xmlfile.write( " <PhaserMrSolutions>\n" )

      # write out what we were trying to find
      total_ncomp = 0
      for i in range(len(self.container.inputData.ASU_COMPONENTS)):
          total_ncomp += self.container.inputData.ASU_COMPONENTS[i].numberOfCopies
      self.xmlfile.write( "  <Target>\n" )
      self.xmlfile.write( "    <TotalComps>"+str(total_ncomp)+"</TotalComps>\n" )
      self.xmlfile.write( "    <CompTypes>"+str(len(self.container.inputData.ASU_COMPONENTS))+"</CompTypes>\n" )
      self.xmlfile.write( "  </Target>\n" )

      if not self.results.foundSolutions():
         self.xmlfile.write( "  <solutionsFound>False</solutionsFound>\n" )
      else:
         self.xmlfile.write( "  <solutionsFound>True</solutionsFound>\n" )
         self.xmlfile.write( "  <Solutions>\n" )

         isol = 0
         # available items are in phaser/source/phaser/include/mr_set.h
         # example of usage in phaser/source/phaser/phaser/test_reporter.py
         for solution in self.results.getDotSol():
           self.xmlfile.write( "    <Solution>\n" )
           isol += 1
           self.xmlfile.write( "      <ISOL>"+str(isol)+"</ISOL>\n" )
           self.xmlfile.write( "      <SPG>"+str(solution.getSpaceGroupName())+"</SPG>\n" )
           # full SOLU SET line
           self.xmlfile.write( "      <ANNOTATION>"+str(solution.ANNOTATION)+"</ANNOTATION>\n" )
           # final LLG
           self.xmlfile.write( "      <LLG>"+str(solution.LLG)+"</LLG>\n" )
           # TFZ from actual TF search
           self.xmlfile.write( "      <TFZ>"+str(solution.TFZ)+"</TFZ>\n" )
           # final TFZ after refinement
           self.xmlfile.write( "      <TFZeq>"+str(solution.TFZeq)+"</TFZeq>\n" )
           # number of clashes
           self.xmlfile.write( "      <PAK>"+str(solution.PAK)+"</PAK>\n" )
           ### following maybe not interesting, but included for completeness
           # value of TF from actual search i.e. goes with TFZ
           self.xmlfile.write( "      <TF>"+str(solution.TF)+"</TF>\n" )
           # LLG at start of refinement cycle
           self.xmlfile.write( "      <ORIG_LLG>"+str(solution.ORIG_LLG)+"</ORIG_LLG>\n" )
           # value of RF after refinement
           self.xmlfile.write( "      <R>"+str(solution.R)+"</R>\n" )
           # value of RF before refinement
           self.xmlfile.write( "      <ORIG_R>"+str(solution.ORIG_R)+"</ORIG_R>\n" )
           # number requested?
           self.xmlfile.write( "      <ORIG_NUM>"+str(solution.ORIG_NUM)+"</ORIG_NUM>\n" )
           # number found?
           self.xmlfile.write( "      <NUM>"+str(solution.NUM)+"</NUM>\n" )
           for nd in solution.KNOWN:
              self.xmlfile.write( "      <COMPONENT>\n" )
              self.xmlfile.write( "        <known>"+str(nd.getModlid())+"</known>\n" )
              self.xmlfile.write( "      </COMPONENT>\n" )
           self.xmlfile.write( "    </Solution>\n" )

         self.xmlfile.write( "  </Solutions>\n" )
         self.xmlfile.write( "  <numSolutions>"+str(self.results.numSolutions())+"</numSolutions>\n" )

      self.xmlfile.write( " </PhaserMrSolutions>\n" )
      self.xmlfile.write( "</PhaserMrResults>\n" )
      self.xmlfile.close()

      #print dir(self.results.getTemplatesForSolution(0))
      #print dir(self.results.getDotSol()[0])

      return CPluginScript.SUCCEEDED

