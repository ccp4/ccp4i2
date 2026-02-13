import os

from ccp4i2.core.CCP4PluginScript import CPluginScript


class chainsaw(CPluginScript):
    TASKTITLE = 'edit search model'
    TASKNAME = 'chainsaw'
    PERFORMANCECLASS = 'CAtomCountPerformance'
    TASKCOMMAND = 'chainsaw'

    def makeCommandAndScript(self):
      print('chainsaw.makeCommandAndScript')

      if self.container.inputData.XYZIN.isSelectionSet():
        xyzin_file = os.path.join(self.getWorkDirectory(),"XYZIN_sel.pdb")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(xyzin_file)
      else:
        xyzin_file = str(self.container.inputData.XYZIN)
      self.appendCommandLine(['XYZIN',xyzin_file])
      self.appendCommandLine(['ALIGNIN',self.inputAlignmentFileName])
      self.appendCommandLine(['XYZOUT',self.container.outputData.XYZOUT.fullPath])

      ### keywords
      self.appendCommandScript('MODE %s' % self.container.controlParameters.MODE)
      self.appendCommandScript('END')

    def processInputFiles(self):
        self.cryst1card = None
        with open(str(self.container.inputData.XYZIN.fullPath),'r') as inputFile:
            lines = inputFile.readlines()
            for line in lines:
                if line.startswith('CRYST1'):
                    self.cryst1card = line

        # Ensure correct alignment file extension and reorder file if target is not first
        formt,idList = self.container.inputData.ALIGNIN.identifyFile()
        ext = { 'unknown' : 'aln',
                'clustal' : 'aln' ,
                'fasta' : 'fas',
                'pir' : 'pir' ,
                'phylip' : 'phy' }.get(formt,'aln')
        self.inputAlignmentFileName = os.path.normpath(os.path.join(self.getWorkDirectory(),'tempAlignment.'+ext))
        #print 'inputAlignmentFileName',self.inputAlignmentFileName,'TARGETINDEX',self.container.controlParameters.TARGETINDEX
        if self.container.controlParameters.TARGETINDEX.isSet() and self.container.controlParameters.TARGETINDEX != 0:
          formt,idList = self.container.inputData.ALIGNIN.identifyFile()
          self.container.inputData.ALIGNIN.convertFormat(formt,self.inputAlignmentFileName,reorder='reverse')
        else:
          # Forcing the file extension to chainsaw requirement
          import shutil
          shutil.copyfile(self.container.inputData.ALIGNIN.__str__(), self.inputAlignmentFileName)
          return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      if self.cryst1card is not None:
          import os
          tmpFilename = str(self.container.outputData.XYZOUT.fullPath)+'_tmp'
          os.rename(str(self.container.outputData.XYZOUT.fullPath), tmpFilename)
          with open(tmpFilename,'r') as inputFile:
              lines = inputFile.readlines()
              with open(str(self.container.outputData.XYZOUT.fullPath),'w') as outputFile:
                  outputFile.write(self.cryst1card)
                  for line in lines:
                      outputFile.write(line)
      self.container.outputData.XYZOUT.subType=2

      self.container.outputData.XYZOUT.annotation = 'Chainsawed model coordinates'

      # sanity check that chainsaw has produced something
      # -1 means we have not managed to get value out of log file
      self.container.outputData.NO_DELETED = -1
      self.container.outputData.NO_CONSERVED = -1
      self.container.outputData.NO_MUTATED = -1
      logText = self.logFileText()
      pyListLogLines = logText.split("\n")
      for j, pyStrLine in enumerate(pyListLogLines):
         if "conserved" in pyStrLine and "mutated" in pyStrLine:
            self.container.outputData.NO_DELETED = int(pyStrLine.split()[0])
            self.container.outputData.NO_CONSERVED = int(pyStrLine.split()[2])
            self.container.outputData.NO_MUTATED = int(pyStrLine.split()[4])

      if self.container.outputData.NO_CONSERVED > -1 and self.container.outputData.NO_MUTATED > -1:
         self.container.outputData.NUMRES_MODEL = self.container.outputData.NO_CONSERVED + self.container.outputData.NO_MUTATED
         self.container.outputData.SEQID = float(self.container.outputData.NO_CONSERVED) / float(self.container.outputData.NUMRES_MODEL)

      self.container.saveDataToXml(self.makeFileName( 'PROGRAMXML' ))

      # Set a performance parameter mostly used in i2 testing
      self.container.outputData.PERFORMANCE.setFromPdbDataFile(self.container.outputData.XYZOUT)
      #print 'PERFORMANCE',self.container.outputData.PERFORMANCE.nAtoms,self.container.outputData.PERFORMANCE.nResidues

      return CPluginScript.SUCCEEDED
