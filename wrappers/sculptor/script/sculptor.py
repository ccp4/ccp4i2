from __future__ import print_function

"""
     sculptor.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
     Author: Martyn Winn

     Wrapper to phaser.sculptor
"""

from core.CCP4PluginScript import CPluginScript
from xml.etree import ElementTree as ET

class sculptor(CPluginScript):

    TASKMODULE = 'bioinformatics'
    TASKTITLE = 'Truncate search model - SCULPTOR'
    TASKNAME = 'sculptor'  
    TASKVERSION= 0.1
    PERFORMANCECLASS = 'CAtomCountPerformance'

    # used by the base class startProcess()
    TASKCOMMAND = 'phaser.sculptor'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None 
    COMTEMPLATE = None
    
    ERROR_CODES = { 201 : {'description' : 'Unable to convert the provided alignment to clustal (.aln) format' },
                    202 : {'description' : 'Failed reading the alignment file' }, }

    def makeCommandAndScript(self):

      self.appendCommandLine(['--stdin'])


      ### input block
      self.appendCommandScript("input {")
      self.appendCommandScript("model { file_name = %s }" % self.container.inputData.XYZIN)
      if self.container.inputData.ALIGNMENTORSEQUENCEIN.__str__() == 'ALIGNMENT':
          self.appendCommandScript("alignment { file_name = %s \ntarget_index = %d}" % (self.inputAlignmentFileName, int(self.container.controlParameters.TARGETINDEX)+1))
      else:
          self.appendCommandScript("sequence { file_name = %s \nchain_ids = %s}" % (self.container.inputData.SEQUENCEIN, self.container.controlParameters.CHAINIDS))
      self.appendCommandScript("}")

      ### output block
      self.appendCommandScript("output {")
      self.appendCommandScript("job_title = %s" % sculptor.TASKTITLE)
      self.appendCommandScript("folder = %s" % self.workDirectory)
      self.appendCommandScript("root = ''")
      self.appendCommandScript("}")

      ### macromolecule block
      self.appendCommandScript("macromolecule {")
      if self.container.controlParameters.DELETION.isSet():
         self.appendCommandScript("deletion {")
         self.appendCommandScript("use = %s" % self.container.controlParameters.DELETION)
         self.appendCommandScript("}")
      if self.container.controlParameters.POLISHING.isSet():
         self.appendCommandScript("polishing {")
         self.appendCommandScript("use = %s" % self.container.controlParameters.POLISHING)
         self.appendCommandScript("}")
      if self.container.controlParameters.PRUNING.isSet():
         self.appendCommandScript("pruning {")
         self.appendCommandScript("use = %s" % self.container.controlParameters.PRUNING)
         self.appendCommandScript("}")
      if self.container.controlParameters.BFACTOR.isSet():
         self.appendCommandScript("bfactor {")
         self.appendCommandScript("use = %s" % self.container.controlParameters.BFACTOR)
         self.appendCommandScript("}")
      if self.container.controlParameters.RENUMBER.isSet():
         self.appendCommandScript("renumber {")
         self.appendCommandScript("use = %s" % self.container.controlParameters.RENUMBER)
         self.appendCommandScript("}")
      self.appendCommandScript("}")

      return 0

    def processInputFiles(self):
        import shutil,os
        if self.container.inputData.ALIGNMENTORSEQUENCEIN.__str__() == 'ALIGNMENT':
          self.inputAlignmentFileName = os.path.join(self.workDirectory,'alignIn.aln')
          formt,identifiers = self.container.inputData.ALIGNIN.identifyFile()
          print('processInputFiles',formt,identifiers) 
          if formt == 'clustal':
            # clustal file should have .aln extension - unknown format liable to fail but let it try
            shutil.copyfile(self.container.inputData.ALIGNIN.__str__(),self.inputAlignmentFileName)
          elif formt == 'unknown':
            # unknown format liable to fail but let it try
            self.inputAlignmentFileName = os.path.join(self.workDirectory,'alignIn.seq')
            shutil.copyfile(self.container.inputData.ALIGNIN.__str__(),self.inputAlignmentFileName)
          else:
            self.container.inputData.ALIGNIN.convertFormat('clustal',self.inputAlignmentFileName)
        return  CPluginScript.SUCCEEDED
            

    def processOutputFiles(self):
        # Import PDB files that have been output

        import os, glob, shutil
        from core import CCP4Utils
        globPath = os.path.normpath(os.path.join(self.workDirectory,'_*.pdb'))
        outList = glob.glob(globPath)
        xyzoutList = self.container.outputData.XYZOUT
        nGood = 0
        for iFile in range(len(outList)):
          # Beware sculptor seems to create empty pdb files
          txt = CCP4Utils.readFile(outList[iFile])
          print('pdb file length',len(txt))
          if len(txt)<5:
            pass
          else:
            nGood += 1
            fpath,fname = os.path.split(outList[iFile])
            xyzoutList.append(xyzoutList.makeItem())
            outputFilePath = os.path.normpath(os.path.join(self.workDirectory,'XYZOUT_'+str(nGood)+'.pdb'))
            shutil.copyfile(outList[iFile], outputFilePath)
            xyzoutList[-1].setFullPath(outputFilePath)
            if len(outList)>1:
              xyzoutList[-1].annotation = "Edited search model number "+str(nGood)
            else:
              xyzoutList[-1].annotation = "Edited search model"
            xyzoutList[-1].subType = 2

        # Create a trivial xml output file
        from core import CCP4Utils
        #from lxml import etree
        root = ET.Element('sculptor')
        e = ET.Element('number_output_files')
        e.text = str(len(xyzoutList))
        root.append(e)
        with open(self.makeFileName('PROGRAMXML'),"w") as programXMLFile:
            CCP4Utils.writeXML(programXMLFile,ET.tostring(root))

        self.container.outputData.PERFORMANCE.setFromPdbDataFile(self.container.outputData.XYZOUT[0])

        if nGood>0:
            return CPluginScript.SUCCEEDED
        else:
            return CPluginScript.FAILED


