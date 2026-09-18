import sys

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import CellCheck
from ccp4i2.pipelines.import_merged.script.mmcifconvert import ConvertCIF
from ccp4i2.pipelines.import_merged.script.mtzimport import ImportMTZ


class import_merged(CPluginScript):

    TASKNAME = 'import_merged'
    # Note - preserving the HKLOUT by changing severity from the system default of 1 to 5 and
    # beware issues with caseinsensitivity
    PURGESEARCHLIST = [ [ 'HKLIN*.mtz' , 1 ],
                        ['aimless_pipe%*/HKLOUT*.mtz', 1],
                        [ 'hklout.mtz' , 5 ],    
                        [ 'HKLOUT.mtz' , 5 ]
                      ]
    #------------------------------------------------------------------------
    def validity(self):
        error = super(import_merged, self).validity()

        # For an mmCIF file the crystal name, dataset name and cell all come
        # from the chosen reflection block.  Without one the base validity()
        # reports only "CRYSTALNAME is not set" / "DATASETNAME is not set",
        # which does not tell the user where to look - and for mmCIF input
        # those two fields are not even on the page.
        hklin_format = str(self.container.inputData.HKLIN_FORMAT).upper()
        if hklin_format == 'MMCIF' and \
                not self.container.inputData.MMCIF_SELECTED_BLOCK.isSet():
            error.append(
                klass=self.TASKNAME, code=201,
                details='No reflection block has been selected. Choose one '
                        'under "mmCIF Reflection Data" - the crystal name, '
                        'dataset name and cell are taken from it.',
                name=f'{self.TASKNAME}.container.inputData.MMCIF_SELECTED_BLOCK',
                severity=CCP4ErrorHandling.SEVERITY_ERROR)

        # Unmerged data is not importable here -- import_merged is for MERGED
        # data. Block it server-side (the sole validation authority) so RUN is
        # disabled, rather than letting the job run and produce nonsense. This
        # reads the file, an exception to the "validity() does no I/O" guideline;
        # the diagnosis is cached by (path, mtime, size) so repeated validation
        # polls read the file only once, and the read is skipped unless HKLIN is
        # set.
        hklin = self.container.inputData.HKLIN
        if hklin.isSet():
            try:
                from ccp4i2.lib.utils.files.reflection_diagnosis import (
                    diagnose_reflection_file_cached,
                )
                diag = diagnose_reflection_file_cached(str(hklin.fullPath))
                if diag.get('merged') is False:
                    error.append(
                        klass=self.TASKNAME, code=202,
                        details='This looks like UNMERGED data. import_merged is '
                                'for merged reflection data - scale and merge it '
                                'first (e.g. the aimless data-reduction task).',
                        name=f'{self.TASKNAME}.container.inputData.HKLIN',
                        severity=CCP4ErrorHandling.SEVERITY_ERROR)

                # Metadata the file does not carry must be supplied by the user
                # (SHELX: cell + space group; the data-type has a default). Block
                # until they are set, rather than failing at run time in the
                # reader. mmCIF/XDS/scalepack carry their own cell/SG.
                needs = diag.get('needs') or []
                if 'cell' in needs and not self.container.inputData.UNITCELL.isSet():
                    error.append(
                        klass=self.TASKNAME, code=203,
                        details='This format carries no unit cell - enter one '
                                'under "Crystal Information".',
                        name=f'{self.TASKNAME}.container.inputData.UNITCELL',
                        severity=CCP4ErrorHandling.SEVERITY_ERROR)
                if 'spaceGroup' in needs and not self.container.inputData.SPACEGROUP.isSet():
                    error.append(
                        klass=self.TASKNAME, code=204,
                        details='This format carries no space group - enter one '
                                'under "Crystal Information".',
                        name=f'{self.TASKNAME}.container.inputData.SPACEGROUP',
                        severity=CCP4ErrorHandling.SEVERITY_ERROR)
            except Exception:
                pass  # never let the merged probe break validation

        return error

    #------------------------------------------------------------------------
    def process(self):
      self.container.inputData.HKLIN.loadFile()
      # Format by CONTENT, not extension: getFormat() keys off the filename and
      # returns 'unknown' for .sca and mis-classifies a .hkl that is really
      # XDS_ASCII. detect_format peeks the bytes. Returns a plain str, so the
      # historical "CString vs str" ambiguity here also goes away.
      from ccp4i2.lib.utils.files.reflection_diagnosis import detect_format
      self.fformat = detect_format(str(self.container.inputData.HKLIN.fullPath))
      merged = self.container.inputData.HKLIN.getMerged()
      self.isintensity = 0  # unknown I or F
      
      obsout = None

      self.x2mtz = None
      self.mmcifXML = None
      self.resolutioncutoff = False
      #print("IDRR", self.container.inputData.RESOLUTION_RANGE_SET)
      if self.container.inputData.RESOLUTION_RANGE_SET:
          self.resolutioncutoff = True

      # Every format now has a binary-free reader (gemmi ImportMTZ for MTZ,
      # ConvertCIF for mmCIF, read_scalepack/read_shelx for .sca/.hkl, gemmi
      # read_xds_ascii for XDS), so the old convert2mtz (f2mtz/combat) plugin is
      # gone. self.x2mtz stays None; the process1/process2 guards that test it
      # are harmless dead branches kept to minimise churn.
      self.freeRcompleteTried = True
      self.importXML = None
      self.freeout = None
      if self.fformat == 'mtz':
        # Both the resolution-cut and no-cut cases now go through the gemmi
        # ImportMTZ path (importmtz). The legacy no-cut branch shelled out to
        # the `cmtzsplit` binary, which is unavailable on the slim server;
        # ImportMTZ applies a (possibly-null) resolution range with gemmi and
        # produces the same OBSOUT/FREEOUT, so the two engines are collapsed
        # into one.  (columnthings() inside importmtz() auto-picks the best
        # observation group when HKLIN_OBS_COLUMNS is unset.)
        fcontent = self.container.inputData.HKLIN.getFileContent()
        # +1 intensity, -1 amplitude, 0 unknown -- used by the QC/report step.
        self.isintensity = self.isIntensity(
            self.container.inputData.HKLIN_OBS_COLUMNS, fcontent.listOfColumns)
        if len(fcontent.datasets) >= 2:
            self.container.inputData.DATASETNAME = fcontent.datasets[1]
        self.importXML = etree.Element('IMPORT_LOG')
        status = self.importmtz()
        self.makeReportXML(self.importXML)
        self.outputLogXML(self.importXML)
        self.process1(status)
        return self.get_status() if self.get_status() is not None else CPluginScript.SUCCEEDED
      else:
          # not MTZ
          self.importXML = etree.Element('IMPORT_LOG')  # information about the import step
          #  +1 if intensity, -1 if amplitude, 0 if unknown
          self.isintensity = 0
          if str(self.fformat) == 'scalepack':
              self.isintensity = +1  # scalepack files are intensity
          if self.container.inputData.MMCIF_SELECTED_ISINTENSITY:
              self.isintensity = self.container.inputData.MMCIF_SELECTED_ISINTENSITY

          self.makeReportXML(self.importXML)  # add initial stuff for XML into self.importXML
          self.outputLogXML(self.importXML)  # send self.importXML to program.xml

          # mmCIF, direct import (gemmi ConvertCIF)
          if str(self.fformat) == 'mmcif':
              status = self.convertmmcif()
              self.process1(status)
              # Return the status that was set by reportStatus()
              return self.get_status() if self.get_status() is not None else CPluginScript.SUCCEEDED

          # Scalepack .sca: pure-Python/gemmi reader (retires the
          # scalepack2mtz + cmtzsplit binaries; slim-safe).
          if str(self.fformat) == 'scalepack':
              status = self.importscalepack()
              self.process1(status)
              return self.get_status() if self.get_status() is not None else CPluginScript.SUCCEEDED

          # SHELX .hkl: pure-Python reader (retires f2mtz). The file declares
          # neither cell/SG nor whether the data are intensities or amplitudes,
          # so all three come from the user (validity() requires cell + SG).
          if str(self.fformat) == 'shelx':
              status = self.importshelx()
              self.process1(status)
              return self.get_status() if self.get_status() is not None else CPluginScript.SUCCEEDED

          # XDS_ASCII: gemmi reads it natively. Only merged XDS reaches here --
          # unmerged XDS (the usual CORRECT/INTEGRATE output) is rejected by the
          # unmerged block in validity().
          if str(self.fformat) == 'xds':
              status = self.importxds()
              self.process1(status)
              return self.get_status() if self.get_status() is not None else CPluginScript.SUCCEEDED

          # Truly unrecognised content: fail with a clear message.
          print("ERROR: import_merged: unsupported reflection format",
                self.fformat)
          self.appendErrorReport(
              201,
              f'Unsupported reflection format: {self.fformat}',
              severity=CCP4ErrorHandling.SEVERITY_ERROR)
          self.process1(CPluginScript.FAILED)
          return self.get_status() if self.get_status() is not None else CPluginScript.FAILED

    def process1(self,status, completeFreeR=True):
        'if completeFreeR False, always generate new FreeR (for 2nd attempt)'
        #print('process1',type(status),status)
        if status == CPluginScript.FAILED:
            self.reportStatus(status)
            return
      
        # Is FreeR generation switched off?
        if self.container.controlParameters.SKIP_FREER:
            # No freeR generation, leave as is, eg from StarAniso
            self.process2(CPluginScript.SUCCEEDED)

        # HASFREER records whether the *imported* file carried FreeR. It must
        # NOT disable completion when the user supplied a separate FREERFLAG --
        # that external set is exactly what we complete (case 1 below). The old
        # cmtzsplit path never set HASFREER, so completion always ran; the gemmi
        # importmtz path sets it False for a FreeR-less MTZ, so guard on both.
        if not self.container.inputData.HASFREER and \
                not self.container.inputData.FREERFLAG.isSet():
            completeFreeR = False   # no valid FreeR data anywhere

        # Create or complete a freer set
        self.freerflag = self.makePluginObject('freerflag')
        if self.x2mtz is not None:
            self.freerflag.container.inputData.F_SIGF = \
                            self.x2mtz.container.outputData.OBSOUT
        else:
            self.freerflag.container.inputData.F_SIGF = \
                        self.container.outputData.OBSOUT
            
        #print 'import_merged.process1',self.x2mtz.container.outputData.FREEOUT,self.x2mtz.container.outputData.FREEOUT.exists()
        newfreer = 'True'
        freeRsource = None
        if completeFreeR:
            # Cases if completeFreeR == True:
            #  1) inputData.FREERFLAG is set (FreeR from separate object|file), complete this one, or
            #     freeRsource = 'Explicit'
            #  2) FREEOUT.exists from main import, complete this
            #     freeRsource = 'Input'
            #  3) else generate new one
            if self.container.inputData.FREERFLAG.isSet():  # case (1)
                self.freerflag.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG
                freeRsource = 'Explicit'
                #  Check compatible cells etc
                #print "****"
                #print "FR", self.freerflag.container.inputData.FREERFLAG.fileContent
                #print "OBSOUT", self.container.outputData.OBSOUT.fileContent
                #print "HKLOUT", self.container.outputData.HKLOUT.fileContent

                tolerance = None # use default
                cellcheck = \
                          CellCheck(self.freerflag.container.inputData.F_SIGF.fileContent,
                                    self.freerflag.container.inputData.FREERFLAG.fileContent,
                                    tolerance)
                cellsAreTheSame, freerReportXML = cellcheck.checks()
                if (not self.container.controlParameters.OVERRIDE_CELL_DIFFERENCE) and \
                       (not cellsAreTheSame['validity']):
                    # not compatible
                    completeFreeR = False

                if freerReportXML is not None:
                    self.importXML.append(freerReportXML)
            elif self.freeout is not None:
                # FREEOUT from mmcif
                # A freeR set has been imported, so extend/complete it
                self.freerflag.container.inputData.FREERFLAG.set(self.freeout)
                freeRsource = 'Input'

            elif self.x2mtz is not None and self.x2mtz.container.outputData.FREEOUT.exists():
                # A freeR set has been imported (via a converter plugin), so
                # extend/complete it. Guarded: x2mtz is None for the MTZ path now.
                self.freerflag.container.inputData.FREERFLAG = self.x2mtz.container.outputData.FREEOUT
                freeRsource = 'Input'
            else:
                completeFreeR = False

        if completeFreeR:
            self.freerflag.container.controlParameters.GEN_MODE = 'COMPLETE'
            self.freerflag.container.controlParameters.COMPLETE = True
            self.freerflag.container.controlParameters.CUTRESOLUTION = \
                          self.container.controlParameters.CUTRESOLUTION
            self._propagateFreerOverride(self.freerflag)
            newfreer = 'False'
          
        self.freerflag.container.controlParameters.FRAC = \
                         self.container.controlParameters.FREER_FRACTION
        self.addElement(self.importXML, 'newFreeR', newfreer) 
        if freeRsource is not None:
            self.addElement(self.importXML, 'freeRsource', freeRsource) 
        self.outputLogXML(self.importXML)  # send self.importXML to program.xml
        self.freerflag.container.outputData.FREEROUT.setFullPath(str(self.container.outputData.FREEOUT))

        status = self.freerflag.process()
        self.process2(status)
        
    #------------------------------------------------------------------------
    def _propagateFreerOverride(self, plugin):
      """The pipeline's cell-difference override covers both cell gates: the
      CellCheck above and the freerflag wrapper's index-only join of the data
      with the input FreeR set (which would otherwise refuse a set from
      another crystal of the same form)."""
      if self.container.controlParameters.OVERRIDE_CELL_DIFFERENCE:
          plugin.container.controlParameters.OVERRIDE_CELL_DIFFERENCE.set(True)

    #------------------------------------------------------------------------
    def process2(self,status):
      freerOK = True
      doFreeR = True
      if self.container.controlParameters.SKIP_FREER:
        doFreeR = False
      else:
        if status == CPluginScript.FAILED:
          # FreeR run has failed, create error message and continue
          self.addElement(self.importXML, 'FreeRfailed', 'True')
          self.outputLogXML(self.importXML)  # send self.importXML to program.xml
          # try again
          if self.freeRcompleteTried:
              print("trying again")
              self.freeRcompleteTried = False
              self.process1(None, False)  # try to make a new FreeR set
          else:
              # failed a 2nd time, make error and continue
              print("failed again")
              self.addElement(self.importXML, 'FreeRfailed', 'Again')
              self.outputLogXML(self.importXML)  # send self.importXML to program.xml
              freerOK = False
              self.reportStatus(status)
              return

        # If FreeR is OK, then create annotation
        if freerOK:
          #print "FREEROUT content",self.freerflag.container.outputData.FREEROUT.fileContent
          if self.freerflag.container.outputData.FREEROUT.fileContent.spaceGroup.isSet():
              sgname = self.freerflag.container.outputData.FREEROUT.fileContent.spaceGroup.__str__()
          else:
              sgname = 'Unk'

          highresFRformatted = "%7.2f" % float(self.freerflag.container.outputData.FREEROUT.fileContent.resolutionRange.high)
          title ='FreeR - Spg:'+str(sgname).strip()+';Resln:'+highresFRformatted.strip() + "A;"
          try:
              title = title + "Cell:"+self.freerflag.container.outputData.FREEROUT.fileContent.cell.guiLabel()
          except Exception as e:
              print('Error writing cell parameters',e)

          #print "FreeR title:",title
          self.container.outputData.FREEOUT.annotation = title

      # Add x2mtz XML if present
      if self.x2mtz is not None:
          x2mtz_XMLpath = self.x2mtz.makeFileName('PROGRAMXML')
          x2mtz_Etree = etree.parse(x2mtz_XMLpath)
          x2mtz = x2mtz_Etree.getroot()
          self.importXML.append(x2mtz)
      if self.mmcifXML is not None:
          self.importXML.append(self.mmcifXML)

      # add in FreeR XML
      if doFreeR:
          self.importXML.append(self.freerflag.getXML())
      else:
          self.addElement(self.importXML, 'freeRsource', 'None') 

      # Run aimless for a report on data quality
      self.aimlesspipe = self.makePluginObject('aimless_pipe',pluginTitle='DR run for data analysis')
      unmergedList = self.aimlesspipe.container.inputData.UNMERGEDFILES
      #print '\nunmergedList 0',    unmergedList
      if len(unmergedList)==0: unmergedList.addItem()
      # Always do analysis on the file which is saved as pipeline output
      unmergedList[0].file.set(self.container.outputData.OBSOUT.__str__())
      ##  earlier versions in some cases analysed the input file
      #if self.fformat in ['mmcif']:
      #    unmergedList[0].file.set(self.container.outputData.HKLOUT.__str__())
      #      elif  self.fformat in ['mtz']:
      #          unmergedList[0].file.set(self.container.outputData.OBSOUT.__str__())
      #      else:
      #          unmergedList[0].file.set(self.container.inputData.HKLIN.__str__())
      #print 'unmergedList 1',    unmergedList
      xname = self.filteredName(str(self.container.inputData.CRYSTALNAME), 'X')
      #xname = CCP4Utils.safeOneWord(str(self.container.inputData.CRYSTALNAME))
      dname = self.filteredName(str(self.container.inputData.DATASETNAME), 'D')
      #dname = CCP4Utils.safeOneWord(str(self.container.inputData.DATASETNAME))
      unmergedList[0].crystalName.set(xname)
      unmergedList[0].dataset.set(dname)
      #print 'unmergedList 2',    unmergedList
      unmergedList[0].cell.set(self.container.inputData.UNITCELL)
      #print 'unmergedList 3',    unmergedList
      #print 'self.container.inputData',self.container.inputData
      unmergedList[0].wavelength.set(self.container.inputData.WAVELENGTH)

      # parameters for Pointless
      self.aimlesspipe.container.controlParameters.MODE = 'CHOOSE'
      self.aimlesspipe.container.controlParameters.CHOOSE_MODE = 'SPACEGROUP'
      self.aimlesspipe.container.controlParameters.CHOOSE_SPACEGROUP = \
            self.container.inputData.SPACEGROUP
      # parameters for Aimless
      self.aimlesspipe.container.controlParameters.SCALING_PROTOCOL = 'CONSTANT'
      self.aimlesspipe.container.controlParameters.ONLYMERGE = True
      self.aimlesspipe.container.controlParameters.ANALYSIS_MODE = True
      self.aimlesspipe.container.controlParameters.OUTPUT_UNMERGED = False
      self.aimlesspipe.container.controlParameters.SDCORRECTION_OVERRIDE = True
      self.aimlesspipe.container.controlParameters.SDCORRECTION_REFINE = False
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SET = True
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDFAC = 1.0
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDB = 0.0
      self.aimlesspipe.container.controlParameters.SDCORRECTION_SDADD = 0.0
      

#  Probably shouldn't run Phaser but try anyway
#      if self.isintensity < 0:
#         print("* Fs input, don't run Phaser")
#          self.aimlesspipe.container.controlParameters.DOPHASERANALYSIS = False

      tempXML = self.importXML
      self.addElement(tempXML, "DRPIPE_RUNNING", "True") 
      self.outputLogXML(tempXML)  # SEND self.importXML to program.xml
      #  Start data reduction
      print("starting aimless_pipe")
      status = self.aimlesspipe.process()
      self.nearlyDone(status)
      return CPluginScript.SUCCEEDED

    #------------------------------------------------------------------------
    def nearlyDone(self,status):
      print('import_merged.nearlyDone')
      self.container.outputData.OBSOUT.setContentFlag(reset=True)
      try:
          # XML data: We have
          #   a) self.importXML etree element report on the import step
          #   b) aimless pipe report
          # so put these together into an IMPORT_MERGED block
          aimless_pipe_XMLpath = self.aimlesspipe.makeFileName('PROGRAMXML')
          aimless_pipe_Etree = etree.parse(aimless_pipe_XMLpath)
          aimless_pipe = aimless_pipe_Etree.getroot()
          #print 'aimless_pipe', type(aimless_pipe), aimless_pipe
          self.outputLogXML(self.importXML, aimless_pipe)  # and output it
      except:
        pass

      print('import_merged.nearlyDone, finished')
      self.reportStatus(status)

    #------------------------------------------------------------------------
    def isIntensity(self, selectedcolumns, listOfColumns):
        """
        check if selection columns are intensity or amplitudes
        return +1 if intensity, -1 if amplitude, 0 if unknown
        selectedcolumns  string of columns labels
        listOfColumns  list eg
        [{'groupIndex': '1', 'columnLabel': 'F_xe1a', 'columnType': 'F', 'dataset': 'xe1a'},
        """
        #print 'isIntensity',selectedcolumns,type(selectedcolumns)
        if not selectedcolumns.isSet():
          selcol1 = str(listOfColumns[0].get()['columnLabel'])
        else:
          selected = selectedcolumns.split(',')
          selcol1 = selected[0] # first label
        #print 'selected', selectedcolumns, selcol1
        columntype = None
        for col in listOfColumns:
            column = col.get()
            if column['columnLabel'] == selcol1:
                #print 'found', column['columnLabel'], column['columnType']
                columntype = column['columnType']
        isintensity = 0
        if columntype is None:
            print("Unrecognised column " + selcol1)
        elif (columntype == 'J') or (columntype == 'K'):
            isintensity = +1
        elif (columntype == 'F') or (columntype == 'G'):
            isintensity = -1

        return isintensity

    #------------------------------------------------------------------------
    def outputLogXML(self, x1XML, x2XML=None):
      'output x1XML and optionally x2XML to program.xml'
      #print "outputLogXML", x1XML
      rootXML = etree.Element('IMPORT_MERGED') # Global XML for everything
      rootXML.append(x1XML)
      if x2XML is not None:
          rootXML.append(x2XML)
      with open (self.makeFileName('PROGRAMXML'),"w") as outputXML:
          CCP4Utils.writeXML(outputXML,etree.tostring(rootXML,pretty_print=True))

    #------------------------------------------------------------------------
    def makeReportXML(self, containerXML):
        'Make initial report XML'
        #print("mrx 1",  type(self.fformat),  self.fformat)
        filename = str(self.container.inputData.HKLIN)  # Input file
        self.addElement(containerXML, 'filename', filename)
        ffmt = str(self.fformat)
        self.addElement(containerXML, 'fileformat', ffmt)

        #print 'type merged', type(self.container.inputData.HKLIN.getMerged())
        if self.container.inputData.HKLIN.getMerged():
            self.addElement(containerXML, 'merged', 'True')
        else:
            self.addElement(containerXML, 'merged', 'False')
        if self.fformat == 'mtz':
            self.addElement(containerXML, 'columnlabels',
                            self.container.inputData.HKLIN_OBS_COLUMNS.get())
        # = +1 if intensity, -1 if amplitude, 0 if unknown
        if self.isintensity == 0:
            IorFtype = 'Unknown'
        elif self.isintensity > 0:
            IorFtype = 'Intensity'
        elif self.isintensity < 0:
            IorFtype = 'Amplitude'
        self.addElement(containerXML, 'IorFtype', IorFtype)

        if self.container.controlParameters.STARANISO_DATA:
            self.addElement(containerXML, 'StarAniso', 'True')
        
        if self.fformat == 'scalepack':
            # Scalepack
            resorange = self.makeResoRange()
            if resorange is not None:
                # XML version of resolution range, a tuple of (dmax, dmin)
                resoxml = etree.Element('ResolutionRange')
                resoxml.set('id', 'cutresolution')
                if resorange[0] > 0.0:
                    self.addElement(resoxml, 'min', "{:.3f}".format(resorange[0]))
                if resorange[1] > 0.0:
                    self.addElement(resoxml, 'max', "{:.3f}".format(resorange[1]))
                containerXML.append(resoxml)
        # mmCIF things
        if ffmt == 'mmcif':
            if self.container.inputData.MMCIF_SELECTED_BLOCK.isSet():
                self.addElement(containerXML, 'mmcifblock',
                                str(self.container.inputData.MMCIF_SELECTED_BLOCK))
            if self.container.inputData.MMCIF_SELECTED_DETAILS.isSet():
                self.addElement(containerXML, 'mmcifblockdetails',
                                str(self.container.inputData.MMCIF_SELECTED_DETAILS))
            if self.container.inputData.MMCIF_SELECTED_INFO.isSet():
                self.addElement(containerXML, 'mmcifblockinfo',
                                str(self.container.inputData.MMCIF_SELECTED_INFO))
            if self.container.inputData.MMCIF_SELECTED_COLUMNS.isSet():
                self.addElement(containerXML, 'mmcifblockcolumns',
                                str(self.container.inputData.MMCIF_SELECTED_COLUMNS))
            
    #------------------------------------------------------------------------
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = etree.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)

    #------------------------------------------------------------------------
    def filteredName(self,name, default=None):
        """ filtered to remove spaces and other rubbish """
        if name is None: return ""
        if len(name) == 0:
            if default is not None: return default
            return ""
        #  funny things with backslash
        if (name.find('\\')):
            name = name.split('\\')[0]
            # only letters and numbers, no spaces or other stuff
        return CCP4Utils.safeOneWord(name)

    #------------------------------------------------------------------------
    def convertmmcif(self):
        # Convert an mmcif file to a data file (OBSOUT) and
        #  (if present) a FreeR file
        # Uses Gemmi

        filename = str(self.container.inputData.HKLIN)  # Input file
        blockname = None    # selected block
        if self.container.inputData.MMCIF_SELECTED_BLOCK.isSet():
            blockname = str(self.container.inputData.MMCIF_SELECTED_BLOCK)

        outfile = str(self.container.outputData.OBSOUT)

        if self.container.controlParameters.SKIP_FREER:
            freerfile = str(self.container.outputData.FREEOUT)
        else:
            # Temporary place for FreeR in job_1 subdirectory. Bind `freerfile`
            # whether or not job_1 already exists (guarding on `not wdir.exists()`
            # raised NameError on a rerun / pre-created dir).
            wdir = self.workDirectory / 'job_1'
            wdir.mkdir(mode=0o777, exist_ok=True)
            freerfile = str(wdir / 'FREEOUT.mtz')

        self.freeout = freerfile
        reducehkl = True  # for now

        # print("convertmmcif files", outfile, freerfile)
        resorange = self.makeResoRange()
            
        cifcontenttype = str(self.container.inputData.MMCIF_SELECTED_CONTENT)
        convertcif = ConvertCIF(filename, blockname, cifcontenttype,
                                outfile, freerfile, reducehkl, resorange)

        self.mmcifXML = convertcif.getXML()
        status = {'finishStatus':CPluginScript.FAILED}
        if convertcif.getstatus():
            contentFlag = convertcif.contentFlag()
            self.container.outputData.OBSOUT.contentFlag.set(contentFlag)
            status = {'finishStatus':CPluginScript.SUCCEEDED}

        return status

    # ---------------------------------------------------------------------------
    def makeResoRange(self):
        resorange = None
        if self.container.inputData.RESOLUTION_RANGE:
            dmax = 0.0
            dmin = 0.0
            r1 = self.container.inputData.RESOLUTION_RANGE.start
            if r1.isSet():
                dmax = float(r1)
            r2 = self.container.inputData.RESOLUTION_RANGE.end
            if r2.isSet():
                dmin = float(r2)
            if r1.isSet() or r2.isSet():
                resorange = (dmax, dmin)
        return resorange

    #------------------------------------------------------------------------
    def importmtz(self):
        # Import an mtz file to a data file (OBSOUT) and
        #  (if present) a FreeR file, with optional resolution cutoffs
        # Uses Gemmi

        outputData = self.container.outputData
        filename = str(self.container.inputData.HKLIN)  # Input file
        outfile = str(outputData.OBSOUT)

        # Input file column labels for observed data and
        #   FreeR (blank '' if no FreeR)
        obsColLabels, freeRcolumnLabel = self.columnthings(filename)
        obsColLabels = list(obsColLabels.split(','))

        if freeRcolumnLabel == '':
            freerfile = None
            self.container.inputData.HASFREER.set(False)
        else:
            if self.container.controlParameters.SKIP_FREER:
                freerfile = str(self.container.outputData.FREEOUT)
            else:
                # Temporary place for FreeR in job_1 subdirectory. `freerfile`
                # must be bound whether or not job_1 already exists -- binding it
                # only inside `if not wdir.exists()` raised NameError on a rerun
                # or when the dir was pre-created.
                wdir = self.workDirectory / 'job_1'
                wdir.mkdir(mode=0o777, exist_ok=True)
                freerfile = str(wdir / 'FREEOUT.mtz')

        self.freeout = freerfile
        reducehkl = True  # for now
        resorange = self.makeResoRange()

        mtzimport = ImportMTZ(filename, outfile, freerfile,
                              obsColLabels, int(self.contentFlag),
                              freeRcolumnLabel,
                              resorange)

        self.mtzXML = mtzimport.getXML()
        self.importXML.append(self.mtzXML)
        # Honour the import result -- do NOT force SUCCEEDED. A failed ImportMTZ
        # was previously reported as success (the verdict was overwritten
        # unconditionally on the next line), so a broken import looked fine.
        if mtzimport.getstatus():
            return {'finishStatus': CPluginScript.SUCCEEDED}
        return {'finishStatus': CPluginScript.FAILED}

    # -------------------------------------------------------------------------
    def importscalepack(self):
        # Import a merged Scalepack .sca file WITHOUT binaries: the pure-Python
        # read_scalepack reader (parity-pinned to scalepack2mtz) produces a
        # source MTZ, which the common gemmi ImportMTZ path then splits to
        # OBSOUT. Replaces the scalepack2mtz + cmtzsplit chain; slim-safe.
        from ccp4i2.lib.utils.files.reflection_formats import read_scalepack
        from ccp4i2.lib.utils.files.reflection_diagnosis import diagnose_reflection_file

        path = str(self.container.inputData.HKLIN)
        diag = diagnose_reflection_file(path)
        cell = diag.get('cell')
        sgnum = diag.get('spaceGroupNumber')
        anomalous = bool(diag.get('anomalous'))

        # A user-supplied cell / space group overrides the .sca header.
        sgc = self.container.inputData.SPACEGROUPCELL
        if sgc.cell.isSet():
            c = sgc.cell
            cell = [c.a.__float__(), c.b.__float__(), c.c.__float__(),
                    c.alpha.__float__(), c.beta.__float__(), c.gamma.__float__()]
        if sgc.spaceGroup.isSet():
            sgnum = sgc.spaceGroup.number()
        if cell is None or sgnum is None:
            print("ERROR: import_merged.importscalepack: no cell/space group for", path)
            return {'finishStatus': CPluginScript.FAILED}

        srcmtz = read_scalepack(path, cell, sgnum, anomalous=anomalous)
        srcpath = str(self.workDirectory / 'scalepack_source.mtz')
        srcmtz.write_to_file(srcpath)

        if anomalous:
            obsColLabels = ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']
            self.contentFlag = 1   # CObsDataFile I(+/-) anomalous
        else:
            obsColLabels = ['IMEAN', 'SIGIMEAN']
            self.contentFlag = 3   # CObsDataFile Imean
        self.isintensity = +1                              # scalepack is intensity
        self.container.inputData.HASFREER.set(False)       # .sca carries no FreeR
        self.freeout = None

        outfile = str(self.container.outputData.OBSOUT)
        resorange = self.makeResoRange()
        mtzimport = ImportMTZ(srcpath, outfile, None,
                              obsColLabels, int(self.contentFlag),
                              None, resorange)
        self.mtzXML = mtzimport.getXML()
        if self.importXML is not None and self.mtzXML is not None:
            self.importXML.append(self.mtzXML)
        if mtzimport.getstatus():
            return {'finishStatus': CPluginScript.SUCCEEDED}
        return {'finishStatus': CPluginScript.FAILED}

    # -------------------------------------------------------------------------
    def _userCellSpaceGroup(self):
        """(cell6, sgNumber) from the user-supplied SPACEGROUP + UNITCELL (the
        Crystal-Information card), falling back to the SPACEGROUPCELL compound.
        For formats that carry no cell/SG of their own (SHELX)."""
        inp = self.container.inputData
        cell = None
        sgnum = None
        if inp.UNITCELL.isSet():
            c = inp.UNITCELL
            cell = [c.a.__float__(), c.b.__float__(), c.c.__float__(),
                    c.alpha.__float__(), c.beta.__float__(), c.gamma.__float__()]
        if inp.SPACEGROUP.isSet():
            try:
                sgnum = inp.SPACEGROUP.number()
            except Exception:
                sgnum = None
        sgc = inp.SPACEGROUPCELL
        if cell is None and sgc.cell.isSet():
            c = sgc.cell
            cell = [c.a.__float__(), c.b.__float__(), c.c.__float__(),
                    c.alpha.__float__(), c.beta.__float__(), c.gamma.__float__()]
        if sgnum is None and sgc.spaceGroup.isSet():
            sgnum = sgc.spaceGroup.number()
        return cell, sgnum

    # -------------------------------------------------------------------------
    def importshelx(self):
        # Import a SHELX .hkl WITHOUT binaries (retires f2mtz). SHELX carries no
        # cell/SG, and whether the two columns are intensities (HKLF 4) or
        # amplitudes (HKLF 3) is undecidable from the file, so all three are
        # user-supplied: cell + SG from the Crystal-Information card (required by
        # validity()), and SHELX_IS_INTENSITY for the data type.
        from ccp4i2.lib.utils.files.reflection_formats import read_shelx

        path = str(self.container.inputData.HKLIN)
        cell, sgnum = self._userCellSpaceGroup()
        if cell is None or sgnum is None:
            print("ERROR: import_merged.importshelx: SHELX needs a cell and "
                  "space group", path)
            return {'finishStatus': CPluginScript.FAILED}

        intensities = bool(self.container.inputData.SHELX_IS_INTENSITY)
        srcmtz = read_shelx(path, cell, sgnum, intensities=intensities)
        srcpath = str(self.workDirectory / 'shelx_source.mtz')
        srcmtz.write_to_file(srcpath)

        if intensities:
            obsColLabels = ['I', 'SIGI']
            self.contentFlag = 3   # Imean
            self.isintensity = +1
        else:
            obsColLabels = ['F', 'SIGF']
            self.contentFlag = 4   # Fmean
            self.isintensity = -1
        self.container.inputData.HASFREER.set(False)       # .hkl carries no FreeR
        self.freeout = None

        outfile = str(self.container.outputData.OBSOUT)
        resorange = self.makeResoRange()
        mtzimport = ImportMTZ(srcpath, outfile, None,
                              obsColLabels, int(self.contentFlag),
                              None, resorange)
        self.mtzXML = mtzimport.getXML()
        if self.importXML is not None and self.mtzXML is not None:
            self.importXML.append(self.mtzXML)
        if mtzimport.getstatus():
            return {'finishStatus': CPluginScript.SUCCEEDED}
        return {'finishStatus': CPluginScript.FAILED}

    # -------------------------------------------------------------------------
    def importxds(self):
        # Import a MERGED XDS_ASCII file. gemmi reads it natively and converts to
        # an MTZ; XDS carries cell/SG/wavelength itself. Only merged XDS reaches
        # here (unmerged XDS is rejected by the unmerged block in validity()).
        import gemmi

        path = str(self.container.inputData.HKLIN)
        xds = gemmi.read_xds_ascii(path)
        srcmtz = xds.to_mtz()
        srcpath = str(self.workDirectory / 'xds_source.mtz')
        srcmtz.write_to_file(srcpath)

        # XDS is intensities; pick anomalous I(+/-) if present, else mean I.
        labels = [c.label for c in srcmtz.columns]
        if 'I(+)' in labels and 'I(-)' in labels:
            obsColLabels = ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']
            self.contentFlag = 1   # I(+/-) anomalous
        elif 'IMEAN' in labels:
            obsColLabels = ['IMEAN', 'SIGIMEAN']
            self.contentFlag = 3
        else:
            obsColLabels = ['I', 'SIGI']
            self.contentFlag = 3   # Imean
        self.isintensity = +1
        self.container.inputData.HASFREER.set(False)
        self.freeout = None

        outfile = str(self.container.outputData.OBSOUT)
        resorange = self.makeResoRange()
        mtzimport = ImportMTZ(srcpath, outfile, None,
                              obsColLabels, int(self.contentFlag),
                              None, resorange)
        self.mtzXML = mtzimport.getXML()
        if self.importXML is not None and self.mtzXML is not None:
            self.importXML.append(self.mtzXML)
        if mtzimport.getstatus():
            return {'finishStatus': CPluginScript.SUCCEEDED}
        return {'finishStatus': CPluginScript.FAILED}

    # -------------------------------------------------------------------------
    def columnthings(self, filename):
        #  Sort out which columns are wanted, cf x2mtz.py
        print('HKLIN_OBS_COLUMNS', self.container.inputData.HKLIN_OBS_COLUMNS)
        print('HKLIN_OBS_CONTENT_FLAG', self.container.inputData.HKLIN_OBS_CONTENT_FLAG)
        #  HKLIN_OBS_COLUMNS   list of wanted input columns, this should be set
        inputData = self.container.inputData
        outputData = self.container.outputData
        inputData.HKLIN_OBS.fileContent.loadFile(filename)
        columnGroups = \
              inputData.HKLIN_OBS.fileContent.getColumnGroups()
        iBestObs, ifree = self.bestcolumns(columnGroups)
        if inputData.HKLIN_OBS_COLUMNS.isSet():
            obsColLabels = str(self.container.inputData.HKLIN_OBS_COLUMNS)
            self.contentFlag = self.container.inputData.HKLIN_OBS_CONTENT_FLAG
        else:
            # Obs columns not already set (should not happen)
            self.contentFlag = columnGroups[iBestObs].contentFlag
            outputData.OBSOUT.contentFlag.set(self.contentFlag)
            outputData.OBSOUT.annotation.set(\
                outputData.OBSOUT.qualifiers('guiLabel')+' from '+\
                inputData.HKLIN.stripedName())
            obsColLabels = str(columnGroups[iBestObs].columnList[0].columnLabel)
            for col in columnGroups[iBestObs].columnList[1:]:
                obsColLabels += ','+str(col.columnLabel)

        #  +1 if intensity, -1 if amplitude, 0 if unknown
        self.isintensity = +1
        if self.contentFlag%2 == 0:
            # F amplitudes
            self.isintensity = -1
        freeRcolumnLabel = ''
        if ifree >= 0:
            #  Freer column, if present
            outputData.FREEOUT.annotation.set(\
                outputData.FREEOUT.qualifiers(\
                    'guiLabel')+' from '+inputData.HKLIN.stripedName())
            freeRcolumnLabel = str(columnGroups[ifree].columnList[0].columnLabel)
        return obsColLabels, freeRcolumnLabel

    # -----------------------------------------------------------------------
    def bestcolumns(self, columnGroups):
        # Try to figure the 'best' obs data and freer data in the input file
        #for cg in columnGroups:
        #    print('bestcolumns cg',cg.get())
        #   print("* type",cg.columnGroupType)

        ifree = -1
        iBestObs = -1
        for ii in range(len(columnGroups)):
            if columnGroups[ii].columnGroupType == 'FreeR':
                ifree = ii  #  Index for FreeR
            elif columnGroups[ii].columnGroupType == 'Obs':
              # print( "@cg ", columnGroups[ii].contentFlag, columnGroups[iBestObs].contentFlag)
              if iBestObs<0 or \
               columnGroups[ii].contentFlag<columnGroups[iBestObs].contentFlag:
                  # Now the "best" obs column based on content type
                  iBestObs = ii
        return iBestObs, ifree   # indices to columngroups for Obs and Free

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function to return name of exportable MTZ
def exportJobFile(jobId=None,mode=None):
    #  If the input file contained intensities,
    #     then return the output from ctruncate + freer
    #     ie I and F and FreeR
    #  If amplitudes F, then return Fs + FreeR
    #     don't use ctruncate output which has intensities derived from F^2
    #     which would mean truncate applied twice
    import os

    from ccp4i2.core import CCP4Modules

    jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
    exportFile = os.path.join(jobDir,'exportMtz.mtz')
    if os.path.exists(exportFile): return exportFile

    db = CCP4Modules.PROJECTSMANAGER().db()
    #print("DB:", db.getJobFilesInfo(jobId=jobId))
    info = db.getJobFilesInfo(jobId=jobId,jobParamName='OBSOUT')
    obsfileContent = info[0]['fileContent']
    #print("obsfilecontent", obsfileContent)
    # intensity I if fileContent == 1 or 3
    isIntensity = (int(obsfileContent)%2 == 1)
    #print("obsfilecontent", obsfileContent, isIntensity)
    obsfilename =info[0]['fileName']

    truncateOut = None
    if isIntensity:
        # Use truncate output
        childJobs = CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=jobId,details=True)
        #print('import_merged.exportMtz',childJobs)
        for jobNo,subJobId,taskName  in childJobs:
            if taskName == 'aimless_pipe':
                aimlessChildJobs = CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=subJobId,details=True)
                print('import_merged.exportMtz aimlessChildJobs',aimlessChildJobs)
                for jobNo0,subJobId0,taskName0  in aimlessChildJobs:
                    if taskName0 == 'ctruncate':
                        truncateOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=subJobId0,create=False),'HKLOUT.mtz')
                        if not os.path.exists(truncateOut): truncateOut = None
        obsOut = truncateOut
    else:
        # Amplitudes F, use original processed file
        obsOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False), obsfilename)


    # No observed-data MTZ located (e.g. intensity input but ctruncate output
    # missing) -> nothing to reconstruct.
    if not obsOut or not os.path.exists(str(obsOut)):
        return None

    info = db.getJobFilesInfo(jobId=jobId,jobParamName='FREEOUT')
    freerfilename = info[0]['fileName']
    freerflagOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False),freerfilename)
    if not os.path.exists(freerflagOut):
        # No FreeR to add: the observed-data file is the best we can offer.
        return str(obsOut)

    # Combine observed data with the FreeR column into one MTZ. gemmi-based
    # (no cad binary); observed data first so it wins on any column clash.
    from ccp4i2.lib.utils.jobs.export import combine_mtz_files
    try:
        return str(combine_mtz_files([obsOut, freerflagOut], exportFile))
    except Exception as e:
        print('ERROR: import_merged.exportJobFile combine failed:', e)
        return None
 
def exportJobFileMenu(jobId=None):
    print("exportJobFileMenu")
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    print("Result:", [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ])
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]
                                                
