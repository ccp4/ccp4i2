from ccp4i2.core.CCP4PluginScript import CPluginScript


class sculptor(CPluginScript):
    TASKNAME = 'sculptor'  
    PERFORMANCECLASS = 'CAtomCountPerformance'
    TASKCOMMAND = 'phaser.sculptor'
    
    ERROR_CODES = { 201 : {'description' : 'Unable to convert the provided alignment to clustal (.aln) format' },
                    202 : {'description' : 'Failed reading the alignment file' },
                    203 : {'description' : 'The atom selection selected nothing' }, }

    def makeCommandAndScript(self):

      self.appendCommandLine(['--stdin'])


      ### input block
      self.appendCommandScript("input {")
      self.appendCommandScript("model { file_name = %s }" % self.modelFileName)
      if self.container.inputData.ALIGNMENTORSEQUENCEIN.__str__() == 'ALIGNMENT':
          self.appendCommandScript("alignment { file_name = %s \ntarget_index = %d}" % (self.inputAlignmentFileName, int(self.container.controlParameters.TARGETINDEX)+1))
      else:
          self.appendCommandScript("sequence { file_name = %s \nchain_ids = %s}" % (self.container.inputData.SEQUENCEIN, self.container.controlParameters.CHAINIDS))
      self.appendCommandScript("}")

      ### output block
      self.appendCommandScript("output {")
      self.appendCommandScript("job_title = Truncate search model - SCULPTOR")
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
        import os
        import shutil
        # The model as selected: Sculptor was given the whole file whatever
        # the selection said (MDM2 job 23 asked for chain A of a four-copy
        # file and got all four). Its own selection PHIL takes cctbx syntax,
        # so the selected atoms are written out, as chainsaw does.
        self.modelFileName = str(self.container.inputData.XYZIN)
        if self.container.inputData.XYZIN.isSelectionSet():
            selected = os.path.join(self.workDirectory, 'XYZIN_selected.pdb')
            if self.container.inputData.XYZIN.getSelectedAtomsPdbFile(selected) != 0 \
                    or not os.path.isfile(selected):
                self.appendErrorReport(203, str(self.container.inputData.XYZIN.selection.text))
                return CPluginScript.FAILED
            self.modelFileName = selected
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

        import glob
        import os
        import shutil

        from ccp4i2.core import CCP4Utils
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

        # What was made: the number of files, and for each the chains and
        # residues it holds and Sculptor's sequence identity per chain (only
        # in its log until now), so a judgement can tell one trimmed chain
        # from several copies, and Phaser can be given the identity.
        import re
        from lxml import etree

        from ccp4i2.core import CCP4Utils
        root = etree.Element('sculptor')
        e = etree.Element('number_output_files')
        e.text = str(len(xyzoutList))
        root.append(e)
        etree.SubElement(root, 'selection_applied').text = str(
            self.container.inputData.XYZIN.isSelectionSet())
        try:
            import gemmi
            for item in xyzoutList:
                model = gemmi.read_structure(str(item.fullPath))[0]
                out = etree.SubElement(root, 'output', file=os.path.basename(str(item.fullPath)))
                etree.SubElement(out, 'chains').text = str(sum(1 for ch in model if len(ch)))
                etree.SubElement(out, 'residues').text = str(sum(len(ch) for ch in model))
        except Exception as err:  # noqa: BLE001 - counts are extra; the outputs stand
            print('Could not count chains of the output:', err)
        try:
            log = CCP4Utils.readFile(self.makeFileName('LOG'))
            for chain, identity in re.findall(r"chain \(id = '([^']*)'\) -> ([0-9.]+)%", log):
                etree.SubElement(root, 'identity', chain=chain).text = identity
        except Exception as err:  # noqa: BLE001
            print('Could not read the identities from the log:', err)
        CCP4Utils.saveEtreeToFile(root,self.makeFileName('PROGRAMXML'))

        if nGood > 0:
            self.container.outputData.PERFORMANCE.setFromPdbDataFile(self.container.outputData.XYZOUT[0])

        if nGood>0:
            return CPluginScript.SUCCEEDED
        else:
            return CPluginScript.FAILED


