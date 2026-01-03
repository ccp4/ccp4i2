from ccp4i2.core.CCP4PluginScript import CPluginScript


class pointless(CPluginScript):

    TASKNAME = 'pointless'   # Task name - should be same as class name
    ##TASKMODULE = 'expt_data_utility'      # Where this plugin will appear on the gui
    TASKTITLE = 'Analyse unmerged dataset (POINTLESS)' # A short title for gui menu
    TASKVERSION= 0.0               # Version of this plugin

    # used by the base class startProcess()
    TASKCOMMAND = 'pointless'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None 
    COMTEMPLATE = None  
    MAINTAINER = 'pre@mrc-lmb.cam.ac.uk'

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def makeCommandAndScript(self):

      par = self.container.controlParameters

      ###self.appendCommandScript("XMLOUT %s" % str( self.makeFileName( 'PROGRAMXML' )))
      # XMLOUT on command line so that syntax errors go into it
      self.appendCommandLine(['XMLOUT',str( self.makeFileName( 'PROGRAMXML' ) )])
      print(">>XMLOUT<<", str( self.makeFileName( 'PROGRAMXML' ) ))
      
      for i in range(len(self.container.inputData.UNMERGEDFILES)):
        # Note: NAME, CIFBLOCK commands must preceed HKLIN
        # print(">>>HKLIN ", self.container.inputData.UNMERGEDFILES[i], file=logit)

        ndatasets = \
          self.container.inputData.UNMERGEDFILES[i].file.fileContent.numberofdatasets
        if str(ndatasets) == 'None':
            ndatasets = 1
        ndatasets = int(ndatasets)
        
        merged = str(self.container.inputData.UNMERGEDFILES[i].file.fileContent.merged)
        merged = (merged == 'merged')

        if self.container.inputData.UNMERGEDFILES[i].crystalName.isSet() and \
               self.container.inputData.UNMERGEDFILES[i].dataset.isSet():
            if (merged or ndatasets<2):
                self.appendCommandScript("NAME PROJECT %s CRYSTAL %s DATASET %s" % \
                                         (self._dbProjectName,
                                   self.container.inputData.UNMERGEDFILES[i].crystalName,
                                   self.container.inputData.UNMERGEDFILES[i].dataset))
        hklin_command = 'HKLIN'
        # if mmCIF, set blockname
        fformat = self.container.inputData.UNMERGEDFILES[i].file.fileContent.format
        if fformat == 'mmcif':
            if self.container.controlParameters.MMCIF_SELECTED_BLOCK.isSet():
                self.appendCommandScript("CIFBLOCK %s"%(self.container.controlParameters.MMCIF_SELECTED_BLOCK))
            hklin_command = 'CIFIN'    # needed by Pointless to process CIFBLOCK command
        self.appendCommandScript("%s %s" %  (hklin_command, self.container.inputData.UNMERGEDFILES[i].file.fullPath))
      # end loop files

      if par.MODE == 'MATCH':
        if par.REFERENCE_DATASET == 'XYZ':
          self.appendCommandScript("XYZIN %s" % self.container.inputData.XYZIN_REF.fullPath)
        else:
          #print "HKLIN_REF ", self.container.inputData.HKLIN_REF
          self.appendCommandScript("HKLREF %s" % self.container.inputData.HKLIN_REF.fullPath)

      # There are some cases where HKLOUT can be a merged file, 
      # for example when using pointless for reindexing.
      if par.WRITE_HKLOUT:
        # unmerged output
        self.appendCommandScript("HKLOUT %s" % self.container.outputData.MTZUNMERGEDOUT.fullPath)
        if merged and (par.MODE == 'MATCH'):
            #  OUTPUT UNMERGED for Pointless from 1.12.16
            self.appendCommandScript("OUTPUT UNMERGED")
        self.container.outputData.deleteObject('MTZMERGEDOUT')
      else:
        self.container.outputData.deleteObject('MTZUNMERGEDOUT')

      if par.CELL.isSet():
          #print("**par.CELL",par.CELL)
          if par.CELL.alpha.isSet():
              self.appendCommandScript("CELL %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f " %
                                       (par.CELL.a,par.CELL.b,par.CELL.c,
                                        par.CELL.alpha,par.CELL.beta,par.CELL.gamma))
          else:
              self.appendCommandScript("CELL %6.3f %6.3f %6.3f 90.0 90.0 90.0" %  (par.CELL.a,par.CELL.b,par.CELL.c))

      #print 'par.WAVELENGTH.isSet()', par.WAVELENGTH.isSet(), par.WAVELENGTH
      if par.WAVELENGTH.isSet() and abs(par.WAVELENGTH)>1e-4:
          self.appendCommandScript("WAVELENGTH %8.5f " %
                                   (par.WAVELENGTH))

      nfiles = len(self.container.inputData.UNMERGEDFILES)
      for i in range(nfiles):
        #print self.container.inputData.UNMERGEDFILES[i].excludeSelection
        if self.container.inputData.UNMERGEDFILES[i].excludeSelection.isSet():
            batch_list = self.container.inputData.UNMERGEDFILES[i].excludeSelection.split(",")
            batch_list_base = "EXCLUDE "
            if nfiles > 1:
                batch_list_base += "FILE %d " % (i+1)

            batch_list_str = " "

            for batch in batch_list:
                if "-" in batch:
                    self.appendCommandScript(batch_list_base+"BATCH %s" % batch.replace("-"," TO "))
                else:
                    batch_list_str += " "+batch

            if batch_list_str != "":
                self.appendCommandScript(batch_list_base+"BATCH %s" % batch_list_str)

      if par.SET_SETTING != 'DEFAULT':
        self.appendCommandScript("SETTING %s" % par.SET_SETTING)

      if par.POINTLESS_USE_RESOLUTION_RANGE:
          s = self.resolutionRangeCommand()
          if s != "":
              self.appendCommandScript(s)

      if par.ISIGLIMIT.isSet():
          self.appendCommandScript("ISIGLIMIT %f" % par.ISIGLIMIT)

      if par.CCHALFLIMIT.isSet():
          self.appendCommandScript("ISIGLIMIT CCHALF %f" % par.CCHALFLIMIT)

      if par.TOLERANCE.isSet():
          self.appendCommandScript("TOLERANCE %f" % par.TOLERANCE)

      if par.MODE == 'CHOOSE':
        if par.CHOOSE_MODE == 'SOLUTION_NO':
          self.appendCommandScript("choose solution %d" % par.CHOOSE_SOLUTION_NO)
        elif par.CHOOSE_MODE == 'LAUEGROUP':
          self.appendCommandScript("choose lauegroup %s" % par.CHOOSE_LAUEGROUP)
        elif par.CHOOSE_MODE == 'SPACEGROUP':
          if len(par.CHOOSE_SPACEGROUP) > 0:
              self.appendCommandScript("choose spacegroup %s" % par.CHOOSE_SPACEGROUP)
        elif par.CHOOSE_MODE == 'REINDEX_SPACE':
          self.appendCommandScript("spacegroup %s" % par.CHOOSE_SPACEGROUP)
          #if par.REINDEX_OPERATOR.isSet():
          #print "REINDEX_SPACE::REINDEX_OPERATOR isset"
          self.appendCommandScript("reindex %s, %s, %s" % (par.REINDEX_OPERATOR.h,par.REINDEX_OPERATOR.k,par.REINDEX_OPERATOR.l))

      if par.MODE == 'COMBINE':
          self.appendCommandScript("copy")

      if par.RUN_MODE == 'BYFILE':
          self.appendCommandScript("run byfile")

      elif par.RUN_MODE == 'BYRANGE':
          
          if par.RUN_BATCHLIST.isSet():
              nruns = len(par.RUN_BATCHLIST)
              for i in range(nruns):
                  runrange = par.RUN_BATCHLIST[i]
                  # print "runrange", runrange.runNumber, runrange.batchRange0, runrange.batchRange1
                  s = 'RUN ' + "%3d " %  runrange.runNumber +\
                      " BATCH %5d " % runrange.batchRange0 + " TO %5d " % runrange.batchRange1
                  self.appendCommandScript(s)

      if par.REMOVE_LATTICE_CENTERING:
          if par.LATTICE_CENTERING.isSet():
              lattype = str(self.container.controlParameters.LATTICE_CENTERING)
              if lattype != 'P':
                  self.appendCommandScript('lattice '+lattype)

      if par.KEEP_LATTICE_CENTERING:
          self.appendCommandScript('keeplattice')
      else:
          if not par.REMOVE_LATTICE_CENTERING:
              if par.LATTICE_CENTERING_THRESHOLD:
                  self.appendCommandScript('keeplattice '+ str(par.LATTICE_CENTERING_THRESHOLD))
              
      if par.ALLOW_NONCHIRAL:
          self.appendCommandScript('CHIRALITY NONCHIRAL')
          
      return 0

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def resolutionRangeCommand(self):
        
        #print "resolutionRangeCommand",par.RESOLUTION_RANGE
        r1 = self.container.controlParameters.RESOLUTION_RANGE.start
        r2 = self.container.controlParameters.RESOLUTION_RANGE.end
        
        high = 0.0
        low  = 0.0
        if not r1.isSet() and not r2.isSet():
            # nothing set
            s = ""
        elif r1.isSet() and r2.isSet():
            low  = float(r1)
            high = float(r2)
            if low < high:
                low, high = high, low
            s = "RESOLUTION LOW %f HIGH %f" % (low, high)
        else:
            if not r1.isSet():
                high = r2            
            if not r2.isSet():
                high = r1
            s = "RESOLUTION HIGH %f" % high

        return s
