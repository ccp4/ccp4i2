import glob
import os
from time import gmtime, strftime

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pimple.logtable import CCP4LogToEtree
from ccp4i2.smartie import smartie


class phaser_mr(CPluginScript):

    TASKMODULE = 'test' # Where this plugin will appear on gui
    TASKTITLE = 'MR using Phaser' # A short title for gui menu
    TASKNAME = 'phaser_mr'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin
    TASKCOMMAND = 'phaser'   # The command to run the executable
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None 
    COMTEMPLATE = None
    TIMEOUT_PERIOD = 9999999.9
    ASYNCHRONOUS = True

    ERROR_CODES = {
        201 : {'description' : 'Failed killing log poller' },
        202 : {'description' : 'Failed analysing problem' },
        203 : {'description' : 'Failed analysing solutions' },
        204 : {'description' : 'Failed analysing log' },
    }

    def makeCommandAndScript(self):

      # contrary to phaserwiki documentation, this seems to be
      # correct syntax, and output file name is generated using ROOT
      self.appendCommandLine(['--xml'])
      inputData = self.container.inputData
      controlParameters = self.container.controlParameters

      if self.container.controlParameters.TITLE.isSet():
          self.appendCommandScript("TITLE %s"%(str(self.container.controlParameters.TITLE)))
      self.appendCommandScript("MODE %s" % self.container.controlParameters.MODE)
      self.appendCommandScript('HKLIN "%s"' % self.hklin )
      self.appendCommandScript('LABIN F=F SIGF=SIGF') 

      # contents of target asu
      #self.appendCommandScript("COMP BY "+controlParameters.COMP_BY.__str__())
      if controlParameters.COMP_BY == 'DEFAULT':
          #Default is 50% solvent ?
          pass
      elif controlParameters.COMP_BY == 'MW':
          if controlParameters.ASU_PROTEIN_MW.isSet():
            self.appendCommandScript("COMP PROTEIN MW " + controlParameters.ASU_PROTEIN_MW.__str__() + " NUM 1")
          if controlParameters.ASU_NUCLEICACID_MW.isSet():
            self.appendCommandScript("COMP NUCLEIC MW " + controlParameters.ASU_NUCLEICACID_MW.__str__() + " NUM 1")

      elif controlParameters.COMP_BY == 'ASU':
          for i in range(len(inputData.ASU_COMPONENTS)):
              self.appendCommandScript("COMP %s SEQU %s NUM %s" % (inputData.ASU_COMPONENTS[i].moleculeType,
                                                                   inputData.ASU_COMPONENTS[i].seqFile,
                                                                   str(inputData.ASU_COMPONENTS[i].numberOfCopies)))
      # search models
      for i in range(len(inputData.ENSEMBLES)):
         print('['+inputData.ENSEMBLES[i].pdbItemList[0].structure.__str__()+']', len(inputData.ENSEMBLES[i].pdbItemList[0].structure.__str__()))
         if len(inputData.ENSEMBLES[i].pdbItemList)>0 and len(inputData.ENSEMBLES[i].pdbItemList[0].structure.__str__())>0:
             cmd_line = 'ENSE %s' % inputData.ENSEMBLES[i].label
             for j in range(len(inputData.ENSEMBLES[i].pdbItemList)):
                pdbItem = inputData.ENSEMBLES[i].pdbItemList[j]
                if pdbItem.structure.isSelectionSet():
                  selAtomsFile = os.path.join(self.workDirectory,'selectedAtomModel_'+str(i)+'_'+str(j)+'.pdb')
                  print('Creating limited selection of atoms from', pdbItem.structure,'in',selAtomsFile)
                  pdbItem.structure.getSelectedAtomsPdbFile(selAtomsFile)
                  cmd_line += ' PDB "%s"' % selAtomsFile
                else:
                  cmd_line += ' PDB "%s"' % pdbItem.structure
                if pdbItem.identity_to_target.isSet():
                  cmd_line += " IDEN %s" % pdbItem.identity_to_target.__str__()
                else:
                  cmd_line += " RMS %s" % pdbItem.rms_to_target.__str__()
             self.appendCommandScript(cmd_line)

      # fixed models
      if inputData.FIXED_STRUCTURE.isSet():
          #len(inputData.FIXED_STRUCTURE) > 0 and len(inputData.FIXED_STRUCTURE[0].pdbItemList) > 0 and str(inputData.FIXED_STRUCTURE[0].pdbItemList[0].structure) != '':
          for i in range(len(inputData.FIXED_STRUCTURE)):
              label = "ENSE FixedModel_"+str(i)
              cmd_line = label
              for j in range(len(inputData.FIXED_STRUCTURE[i].pdbItemList)):
                  pdbItem = inputData.FIXED_STRUCTURE[i].pdbItemList[j]
                  if pdbItem.structure.isSelectionSet():
                      selAtomsFile = os.path.join(self.workDirectory,'selectedAtomModel_Fixed_'+str(i)+'_'+str(j)+'.pdb')
                      print('Creating limited selection of atoms from', pdbItem.structure,'in',selAtomsFile)
                      pdbItem.structure.getSelectedAtomsPdbFile(selAtomsFile)
                      cmd_line += ' PDB "%s"' % selAtomsFile
                  else:
                      cmd_line += ' PDB "%s"' % pdbItem.structure
                  if pdbItem.identity_to_target.isSet():
                      cmd_line += " IDEN %s" % pdbItem.identity_to_target.__str__()
                  else:
                      cmd_line += " RMS %s" % pdbItem.rms_to_target.__str__()
              self.appendCommandScript(cmd_line)

          # Provide null transformation solution cards for the known structure
          if len(inputData.FIXED_STRUCTURE)>0:
              self.appendCommandScript('SOLUtion SET')
              for i in range(len(inputData.FIXED_STRUCTURE)):
                  label = "FixedModel_"+str(i)
                  cmd_line = "SOLUtion 6DIM ENSEmble %s  EULEr 0 0 0 FRACtional 0 0 0" % label
                  self.appendCommandScript(cmd_line)

      #define search

      if self.container.controlParameters.SEARCHMODE == 'multidomain':
        for i in range(len(inputData.ENSEMBLES)):
          if bool(inputData.ENSEMBLES[i].use):
            self.appendCommandScript("SEARCH ENSE %s NUM %s" % (inputData.ENSEMBLES[i].label,
                                   inputData.ENSEMBLES[i].number.__str__() ))
      elif self.container.controlParameters.SEARCHMODE == 'alternatives':
        cmd_line = "SEARCH ENSE %s" % inputData.ENSEMBLES[0].label
        for i in range(1,len(inputData.ENSEMBLES)):
          if bool(inputData.ENSEMBLES[i].use):
            cmd_line += " OR ENSE %s" % inputData.ENSEMBLES[i].label
        cmd_line += " NUM %s" % inputData.ENSEMBLES[0].number.__str__()
        self.appendCommandScript(cmd_line)

      if self.container.controlParameters.SGALT_SELECT.isSet():
         self.appendCommandScript("SGALTERNATIVE SELECT %s" % self.container.controlParameters.SGALT_SELECT)
         if self.container.controlParameters.SGALT_SELECT == 'LIST' and self.container.controlParameters.SGALT_TEST.isSet():
           for i in range(len(self.container.controlParameters.SGALT_TEST)):
             self.appendCommandScript("SGALTERNATIVE TEST %s" % self.container.controlParameters.SGALT_TEST[i])
      if self.container.controlParameters.PERMUTATIONS.isSet():
         if self.container.controlParameters.PERMUTATIONS:
           self.appendCommandScript("PERMUTATIONS ON")
      if self.container.controlParameters.NUM_SOL_OUT.isSet():
         self.appendCommandScript("TOPFILES %d" % self.container.controlParameters.NUM_SOL_OUT)
      if self.container.controlParameters.RESOLUTION_HIGH.isSet():
         self.appendCommandScript("RESOLUTION HIGH %f" % self.container.controlParameters.RESOLUTION_HIGH)
      if self.container.controlParameters.RESOLUTION_AUTO_HIGH.isSet():
          self.appendCommandScript("RESOLUTION AUTO HIGH %f" % self.container.controlParameters.RESOLUTION_AUTO_HIGH)
      if self.container.controlParameters.NJOBS.isSet():
         self.appendCommandScript("JOBS %d" % self.container.controlParameters.NJOBS)
      if self.container.controlParameters.PACK_CUTOFF.isSet():
         #These keywords depend on phaser version. I am assuming 2.5.2
         self.appendCommandScript("PACK SELECT %s" % self.container.controlParameters.PACK_SELECT)
         self.appendCommandScript("PACK CUTOFF %f" % self.container.controlParameters.PACK_CUTOFF)
      if self.container.controlParameters.PEAKS_ROT_CUTOFF.isSet():
         self.appendCommandScript("PEAKS ROT CUTOFF %f" % self.container.controlParameters.PEAKS_ROT_CUTOFF)

      # Set root to correct working directory
      if self.container.controlParameters.ROOT.isSet():
        self.appendCommandScript('ROOT "%s"'%(os.path.join(self.getWorkDirectory(),str(self.container.controlParameters.ROOT))))
      else:
        self.appendCommandScript('ROOT "%s"'%(os.path.join(self.getWorkDirectory(),'PHASER')))

      return 0

    def processInputFiles(self):
      self.xmlText = None
      self.oldLogLength = 0
      self.watchFile(self.makeFileName('LOG'), self.handleLogChanged)
      
      self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED
      else:      
        return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):

        if self.container.controlParameters.NUM_SOL_OUT.isSet():
            num_sol = self.container.controlParameters.NUM_SOL_OUT.isSet()
        else:
            num_sol = 1

        for i in range(1,num_sol+1):
            xyzout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".pdb")
            if os.path.exists(xyzout):
                self.container.outputData.XYZOUT.append(xyzout)
                self.container.outputData.XYZOUT[-1].annotation.set('Positioned coordinates for solution '+str(i))
                self.container.outputData.XYZOUT[-1].subType.set(1)
            
            hklout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".mtz")
            if os.path.exists(hklout):
                self.container.outputData.HKLOUT.append(hklout)

        # Need to set the expected content flag  for phases data
        
        self.splitHkloutList(miniMtzsOut=['MAPOUT','DIFMAPOUT','PHASEOUT'],programColumnNames=['FWT,PHWT','DELFWT,PHDELWT','PHIC,FOM'],outputBaseName=['MAPOUT','DIFMAPOUT','PHASEOUT'],outputContentFlags=[0,0,CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM],infileList=self.container.outputData.HKLOUT)
    
        for indx in range(len(self.container.outputData.MAPOUT)):
            self.container.outputData.MAPOUT[indx].annotation.set('Map for solution '+str(indx+1))
            self.container.outputData.MAPOUT[indx].subType.set(1)
            self.container.outputData.DIFMAPOUT[indx].annotation.set('Difference map for solution '+str(indx+1))
            self.container.outputData.DIFMAPOUT[indx].subType.set(2)
            self.container.outputData.PHASEOUT[indx].annotation.set('Phases for solution '+str(indx+1))
        
        phaserMRElement = self.generateProgramXML()
        newXml = etree.tostring(phaserMRElement,pretty_print=True)

        with open( self.makeFileName( 'PROGRAMXML' )+'.tmp','w') as aFile:
            CCP4Utils.writeXML(aFile,newXml)
            aFile.flush()
        self.renameFile(self.makeFileName( 'PROGRAMXML' )+'.tmp',self.makeFileName( 'PROGRAMXML' ))

        # Add .txt extension to some files to enable display from the 'Project directory' window
        textFileList = glob.glob(os.path.join(self.getWorkDirectory(),'*.sum'))
        textFileList.extend(glob.glob(os.path.join(self.getWorkDirectory(),'*.sol')))
        for textFile in textFileList:
          self.renameFile(textFile,textFile+'.txt')
        

        return CPluginScript.SUCCEEDED
            
    def handleLogChanged(self, logFilename):
        if (os.stat(logFilename).st_size - self.oldLogLength) > 1000:
            self.oldLogLength = os.stat(logFilename).st_size
            phaserMRElement = self.generateProgramXML()
            newXml = etree.tostring(phaserMRElement,pretty_print=True)
            if self.xmlText is None or len(newXml)>len(self.xmlText):
                with open( self.makeFileName( 'PROGRAMXML' )+'.tmp','w') as aFile:
                    aFile.write( newXml )
                    aFile.flush()
                self.renameFile(self.makeFileName( 'PROGRAMXML' )+'.tmp',self.makeFileName( 'PROGRAMXML' ))
                self.xmlText = newXml

    def generateProgramXML(self):
        phaserMRElement = etree.Element("PhaserMrResult")
        total_ncomp = 0

        formattedTime = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
        try:
            self.analyseProblem(phaserMRElement, total_ncomp)
        except:
            self.appendErrorReport(202,formattedTime)
        #analyse solutions (if any) found
        try:
            self.analyseSolfile(phaserMRElement, total_ncomp)
        except:
            self.appendErrorReport(203,formattedTime)
        #use smartie stuff
        try:
            self.appendSmartieStuff(phaserMRElement)
        except:
            self.appendErrorReport(204,formattedTime)
        return phaserMRElement

    def analyseProblem(self,phaserMRElement,total_ncomp):
        #print '\n\n** in analyseProblem'
        #Analyse the problem as given
        total_ncomp = 0
        for i in range(len(self.container.inputData.ASU_COMPONENTS)):
            total_ncomp += self.container.inputData.ASU_COMPONENTS[i].numberOfCopies
        targetElement = etree.SubElement(phaserMRElement,'Target')
        totalCompsElement = etree.SubElement(targetElement,'TotalComps')
        totalCompsElement.text = str(total_ncomp)
        compTypesElement = etree.SubElement(targetElement,'CompTypes')
        compTypesElement.text =str(len(self.container.inputData.ASU_COMPONENTS))
    
    def analyseSolfile(self, phaserMRElement, total_ncomp):
        #print '\n\n** in analyseSolFile'
        phaser_solfile = os.path.join(self.getWorkDirectory(), "PHASER.sol")
        if os.path.exists(phaser_solfile):
            solfile = open( phaser_solfile )
            nsol = 0
            solutionsElement = etree.SubElement(phaserMRElement,'Solutions')
            for line in solfile.readlines():
                if line.strip().startswith( 'SOLU SPAC' ):
                    spacegroup = line.strip()[10:]
                    spaceGroupElement = etree.SubElement(phaserMRElement,'spaceGroup')
                    spaceGroupElement.text = spacegroup
                elif line.strip().startswith( 'SOLU SET' ):
                    solutionElement = etree.SubElement(solutionsElement,'Solution')
                    nsol = nsol + 1
                    iSolElement = etree.SubElement(solutionElement,'ISOL')
                    iSolElement.text=str(nsol)
                    RFZ = []
                    TFZ = []
                    PAK = []
                    LLG = []
                    # this is the TFZ for the refined component
                    refTFZ = []
                    for word in line.split():
                        if word[0:4] == 'RF++':
                            RFZ.append('++')
                        elif word[0:4] == 'RF*0':
                            RFZ.append('*0')
                        elif word[0:3] == 'RFZ':
                            RFZ.append(word[4:])
                        elif word[0:5] == 'TFZ==':
                            refTFZ.append(word[5:])
                        elif word[0:5] == 'TFZ=*':
                            refTFZ.append('*')
                        elif word[0:4] == 'TF*0':
                            refTFZ.append('*0')
                        elif word[0:3] == 'TFZ':
                            TFZ.append(word[4:])
                        elif word[0:3] == 'PAK':
                            PAK.append(word[4:])
                        elif word[0:3] == 'LLG':
                            LLG.append(word[4:])
                    ncomp = len(RFZ)
                    nComponentsElement = etree.SubElement(solutionElement,'NCOMPONENTS')
                    nComponentsElement.text = str(ncomp)
                    allCompFoundElement = etree.SubElement(solutionElement,'AllCompFound')
                    if ncomp < total_ncomp:
                        allCompFoundElement.text = 'False'
                    else:
                        allCompFoundElement.text = 'True'
                    for icomp in range(ncomp):
                        componentElement=etree.SubElement(solutionElement,'Component')
                        rfzElement = etree.SubElement(componentElement,'RFZ')
                        rfzElement.text = str(RFZ[icomp])
                        tfzElement = etree.SubElement(componentElement,'TFZ')
                        tfzElement.text = str(TFZ[icomp])
                        pakElement = etree.SubElement(componentElement,'PAK')
                        pakElement.text = str(PAK[icomp])
                        llgElement = etree.SubElement(componentElement,'LLG')
                        llgElement.text = str(LLG[icomp])
                    #overall scores printed iff all components found
                    if len(LLG) > ncomp:
                        overallLLGElement = etree.SubElement(solutionElement,'overallLLG')
                        overallLLGElement.text = str(LLG[ncomp])
                    #if len(LLG) > ncomp:
                    #refTFZ does not seem to be output by recent versions of phaser
                    #xmlfile.write( "      <overallTFZ>"+refTFZ[ncomp]+"</overallTFZ>\n" )
                else:
                    pass
                    #print 'Unknown line', line
        return


    def appendSmartieStuff(self, programEtree):
        logfile = smartie.parselog(self.makeFileName( 'LOG' ))
        #Collect a list of graphs of different types
        rotationTables = []
        translationTables = []
        packingTables = []
        refinementTables = []

        rotationTablesDict = {}
        translationTablesDict = {}
        refinementTablesDict = {}
        packingTablesDict = {}

        #print 'Log contained tableCount ',str(len(logfile.tables()))
        for smartieTable in logfile.tables():
            #print '\n\n** found a table'
            if smartieTable.ngraphs() > 0:
                #print '\n\n** found a graph table'
                tableelement = CCP4LogToEtree(smartieTable.rawtable())
                graphTableNodes = tableelement.xpath('.//CCP4Table')
                #print '\n\n** found %d graphTableNodes'%len(graphTableNodes)
                graphTableNode = None
                if len(graphTableNodes)>0: graphTableNode = graphTableNodes[0]
                if graphTableNode is not None:
                    #print '\n\n** made a graph table node'
                    if   'Rotation Function Component' in graphTableNode.get("title"):
                        #print '\n\n** made a rotation graph table node'
                        rotationTables.append(tableelement)
                        hashedNumber = graphTableNode.get("title")[27:].split('(')[0].strip()
                        rotationTablesDict[hashedNumber] = tableelement
                    elif 'Translation Function Component' in graphTableNode.get("title"):
                        #print '\n\n** made a translation graph table node'
                        translationTables.append(tableelement)
                        hashedNumber = graphTableNode.get("title")[30:].split('(')[0].strip()
                        if hashedNumber not in translationTablesDict: translationTablesDict[hashedNumber] = []
                        translationTablesDict[hashedNumber].append(tableelement)
                    elif 'Refinement After Placing Component' in graphTableNode.get("title"):
                        #print '\n\n** made a refinement graph table node'
                        refinementTables.append(tableelement)
                        hashedNumber = graphTableNode.get("title")[34:].split('(')[0].strip()
                        if hashedNumber not in refinementTablesDict: refinementTablesDict[hashedNumber] = []
                        refinementTablesDict[hashedNumber].append(tableelement)
                    elif 'Cell Content Analysis' in graphTableNode.get("title"):
                        #print '\n\n** made a refinement graph table node'
                        cellContentProbabilityNode = etree.SubElement(programEtree,"ContentProbability")
                        cellContentProbabilityNode.append(tableelement)
                    elif 'Intensity distribution for Data' in graphTableNode.get("title"):
                        #print '\n\n** made a refinement graph table node'
                        intensityDistributionNode = etree.SubElement(programEtree,"IntensityDistribution")
                        intensityDistributionNode.append(tableelement)

        summaryCount = logfile.nsummaries()
        searchComponentSummaries = []
        #print '\n\n contains %d summaries'%summaryCount
        for iSummary in range(summaryCount):
            summary = logfile.summary(iSummary)
            summaryTextLines = []
            with open(self.makeFileName('LOG')) as myLogFile:
                summaryTextLines = myLogFile.readlines()[summary.start():summary.end()]
            preElement = etree.Element('CCP4Summary')
            if 'Search Order (next search *):' in "".join(summaryTextLines):
                #print '\n\n** Found a search order summary'
                preElement.text = ''
                interesting = False
                for summaryTextLine in summaryTextLines:
                    if len(summaryTextLine.strip()) == 0:
                        interesting = False
                    if 'Search Order (next search *):' in summaryTextLine:
                        interesting = True
                    if interesting: preElement.text += (summaryTextLine)
                searchComponentSummaries.append(preElement)
            else:
                summaryText = "".join(summaryTextLines)
                #print '\n\n** Found a non-search order summary'
                preElement.text = etree.CDATA(summaryText)
                programEtree.append(preElement)

        #print len(searchComponentSummaries), len(rotationTables), len(translationTables), len(refinementTables)
        if len(searchComponentSummaries) > 0:
            searchesElement = etree.SubElement(programEtree,"Searches")
            iSearch = 0
            for searchComponentSummary in searchComponentSummaries:
                searchElement = etree.SubElement(programEtree,"Search")
                searchElement.append(searchComponentSummary)
                # Identify the sought component in this search
                searchSummaryLines = searchComponentSummary.text.split('\n')[1:]
                hashedSoughtComponent = None
                for searchSummaryLine in searchSummaryLines:
                    if searchSummaryLine.strip().endswith('*'):
                        hashedSoughtComponent = searchSummaryLine.split(':')[0].strip()
                #print '\n\nHashedSoughtComponent', hashedSoughtComponent
                if hashedSoughtComponent is not None:
                    if hashedSoughtComponent in rotationTablesDict: searchElement.append(rotationTablesDict[hashedSoughtComponent])
                    if hashedSoughtComponent in translationTablesDict:
                        for translationTable in translationTablesDict[hashedSoughtComponent]:
                            searchElement.append(translationTable)
                    if hashedSoughtComponent in refinementTablesDict:
                        for refinementTable in refinementTablesDict[hashedSoughtComponent]:
                            searchElement.append(refinementTable)
                                
                #print etree.tostring(searchElement,pretty_print=True)

