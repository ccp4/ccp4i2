from __future__ import print_function

"""
     aimless_pipe.py: CCP4 GUI Project
     Copyright (C) 2012 STFC
"""

import sys
import os
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
#from lxml import etree as lxml_etree
from xml.etree import ElementTree as ET
from pipelines.aimless_pipe.script.aimless_pipe_utils import *
from core import CCP4Utils
from core.CCP4Data import CString

class aimless_pipe(CPluginScript):

    TASKMODULE = 'data_reduction'      # Where this plugin will appear on the gui
    TASKTITLE = 'Scale and merge data' # A short title for gui menu
    TASKNAME = 'aimless_pipe'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin
    PERFORMANCECLASS = 'CDataReductionPerformance'
    MAINTAINER = 'pre@mrc-lmb.cam.ac.uk'
    ERROR_CODES = {
        201 :{'description':'Pointless failed'},
        202 :{'description':'Aimless failed'},
        203 :{'description':'Ctruncate failed'},
        204 :{'description':'FreeR  failed'},
        205 :{'description':'Data useless'}
        }
    PURGESEARCHLIST =  [[ 'ctruncate%*/hklout.mtz', 0],
                        [ 'aimless%*/*.xmgr', 1]
                        ]
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process(self):
      self.rootXML = ET.Element('AIMLESS_PIPE')
      # Start up
      print("### Starting aimless pipeline ###")

      unsetData = self.checkInputData()
      if len(unsetData)>0:
         self.reportStatus(CPluginScript.FAILED)
         return

      self.fatalError = None

      self.aimless1xml = None
      self.phaser_analysisxml = None
      self.phaser_analysis1xml = None

      # Set parameters to control run
      self.doPhaserAnalysis = self.container.controlParameters.DOPHASERANALYSIS

      self.process_pointless()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process_pointless(self):

      print("### Running pointless ###")
      #print self.container.inputData

      nfiles = len(self.container.inputData.UNMERGEDFILES)
      for i in range(nfiles):
          if not self.container.inputData.UNMERGEDFILES[i].file.fileContent.knowncell:
              self.container.controlParameters.CELL = self.container.inputData.UNMERGEDFILES[i].cell

      #print "SCALING_PROTOCOL",self.container.controlParameters.SCALING_PROTOCOL
      # Run pointless
      self.pointless = self.makePluginObject('pointless')
      self.pointless.container.inputData.copyData(self.container.inputData,['UNMERGEDFILES'])
      if self.container.inputData.HKLIN_REF:
          self.pointless.container.inputData.copyData(self.container.inputData,['HKLIN_REF'])
      if self.container.inputData.HKLIN_REF:
          self.pointless.container.inputData.copyData(self.container.inputData,['XYZIN_REF'])

      self.pointless.container.controlParameters.copyData \
        (self.container.controlParameters,
         ['EXCLUDE_BATCH','EXCLUDED_BATCHES',
          'POINTLESS_USE_RESOLUTION_RANGE', 'RESOLUTION_RANGE',
          'ISIGLIMIT','CCHALFLIMIT','TOLERANCE','MODE','REFERENCE_DATASET','SET_SETTING',
          'CHOOSE_MODE','CHOOSE_SOLUTION_NO','CHOOSE_LAUEGROUP','CHOOSE_SPACEGROUP',
          'REINDEX_OPERATOR','CELL','WAVELENGTH','RUN_MODE','RUN_BATCHLIST',
          'REMOVE_LATTICE_CENTERING','LATTICE_CENTERING','ALLOW_NONCHIRAL',
          'MMCIF_SELECTED_BLOCK'])

      self.connectSignal(self.pointless,'finished',self.process_cycle_aimless)
      self.pointless.process()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    @QtCore.Slot(dict)
    def process_cycle_aimless(self,status):
    # Check Pointless ran OK, if so pick up its XML
    # Loop process_aimless once or twice, depending on AUTOCUTOFF
        print("process_cycle_aimless = ", status.get('finishStatus'))
        if status.get('finishStatus') == CPluginScript.FAILED:
            print("failed in Pointless", status)
            self.fatalError = [201, 'Pointless failed', status]
            #         self.appendErrorReport(201, details='Pointless failed')
            #         self.reportStatus(status)
            if status.get('finishStatus') == CPluginScript.FAILED:
                print("failed")
                self.process_finish(CPluginScript.FAILED)
                self.reportStatus(status)
            return

        try:
            pointlessXMLPath = self.pointless.makeFileName('PROGRAMXML')
            pointlessEtree = ET.parse(pointlessXMLPath)
            self.rootXML.append(pointlessEtree.getroot())
            with open (self.makeFileName('PROGRAMXML'),"w") as outputXML:
                CCP4Utils.writeXML(\
                    outputXML,ET.tostring(self.rootXML,
                                             pretty_print=True,
                                             encoding='unicode'))
        except:
            pass

        # Loop aimless runs
        self.AUTOCUTOFF = self.container.controlParameters.AUTOCUTOFF
        # option XML from 1st Aimless run, set if Aimless was run twice

        if self.AUTOCUTOFF and self.container.controlParameters.RESOLUTION_RANGE:
            # Probably doesn't make sense to do AUTOCUTOFF if resrange is expicit
            pass  ##  place holder if needed

        if not self.doPhaserAnalysis:
            print("Set AUTOCUTOFF off because doPhaserAnalysis is off")
            self.AUTOCUTOFF = False   # can't do AUTOCUTOFF without Phaser
        
        # Becomes True if AUTOCUTOFF set and 2nd Aimless run was done
        #   2nd run is not done if data go to the maximum resolution anyway
        self.cutoffdone = False
        
        aimlessruncount = [2]   # no autocutoff, just do "2nd" Aimless
        if self.AUTOCUTOFF:
            aimlessruncount = [1,2];  # autocutoff, maybe do two runs
            
        # this will be set > 0 if a 2nd Aimless was run
        self.highrescutoff = -1.0
      
        for self.aimlessruncount in aimlessruncount:
            self.process_aimless()

        # Never returns here

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process_aimless(self):
        # aimlessruncount = 1 for initial run to get resolution llimit
        #                 = 2 for final run, even if only one
        print("***### Running aimless ", self.aimlessruncount, "###")

        try:
            # Run aimless to edit model
            self.aimless = self.makePluginObject('aimless')
            # Copy all relevant parameters from GUI to Aimless script
            self.aimless.container.controlParameters.copyData \
             (self.container.controlParameters,
         ['REFERENCE_DATASET', 'REFERENCE_FOR_AIMLESS',
          'EXCLUDE_BATCH','EXCLUDED_BATCHES',
          'ACCEPT_OVERLOADS','ACCEPT_EDGES','ACCEPT_XDS_MISFITS',
          'SCALING_PROTOCOL','BFACTOR_SCALE','SCALES_ROTATION_TYPE',
          'SCALES_ROTATION_SPACING','SCALES_ROTATION_NBINS',
          'SCALES_BROTATION_TYPE','SCALES_BROTATION_SPACING',
          'SCALES_BROTATION_NBINS',
          'SCALES_SECONDARY_NSPHHARMONICS',
          'SCALES_TILETYPE','SCALES_NTILEX','SCALES_NTILEY',
          'OUTLIER_OVERRIDE','OUTLIER_EMAX',
          'OUTLIER_SDMAX','OUTLIER_SDMAX2',
          'OUTLIER_SDMAXALL','OUTLIER_SDMAXALL_ADJUST','OUTLIER_COMBINE',
          'RUN_MODE', 'RUN_BATCHLIST',
          'SDCORRECTION_OVERRIDE','SDCORRECTION_REFINE','SDCORRECTION_FIXSDB',
          'SDCORRECTION_OPTIONS','SDCORRECTION_DAMP','SDCORRECTION_SIMILARITY_SDFAC',
          'SDCORRECTION_SIMILARITY_SDB','SDCORRECTION_SIMILARITY_SDADD',
          'SDCORRECTION_SET',
          'SDCORRECTION_SDFAC','SDCORRECTION_SDB','SDCORRECTION_SDADD',
          'SDCORRECTION_TIESDB_SD',
          'INTENSITIES_OVERRIDE','INTENSITIES_OPTIONS',
          'PARTIALS_TEST','PARTIALS_FRACLOW','PARTIALS_FRACHIGH',
          'PARTIALS_CHECK','PARTIALS_SCALE','PARTIALS_SCALE_MIN',
          'SCALING_DETAILS','CYCLES_FLAG','CYCLES_N',
          'PARALLEL','PARALLEL_MODE','NPROC','FRACPROC',
          'SELECT1','SELECT_IOVSDMIN','SELECT2','SELECT_EMIN',
          'REFINE_REFERENCE',
          'TIE_ROTATION','TIE_ROTATION_SD',
          'TIE_BFACTOR','TIE_BFACTOR_SD',
          'TIE_SURFACE','TIE_SURFACE_SD',
          'TIE_BZERO','TIE_BZERO_SD','ONLYMERGE',
          'CCMINIMUM', 'CCANOMMINIMUM', 'ISIGMINIMUM'
          ])

            self.aimless.container.inputData.copyData \
                (self.container.inputData,
                 ['HKLIN_REF', 'XYZIN_REF'])
            self.aimless.container.controlParameters.OUTPUT_MODE.set([CString('MERGED')])
            if self.container.controlParameters.OUTPUT_UNMERGED:
                if self.container.controlParameters.ORIGINAL_HKL:
                    self.aimless.container.controlParameters.OUTPUT_MODE.set([CString('MERGED'),CString('UNMERGED'),CString('ORIGINAL')])
                else:
                    self.aimless.container.controlParameters.OUTPUT_MODE.set([CString('MERGED'),CString('UNMERGED')])

            #print("**4 ", self.aimless.container.controlParameters.OUTPUT_MODE)

            self.aimless.container.inputData.UNMERGEDFILE.set(self.pointless.container.outputData.MTZUNMERGEDOUT)

            if self.aimlessruncount == 1:
                # Optional 1st run to get resolution limit
                self.aimless.container.controlParameters.copyData (\
                    self.container.controlParameters,
                    ['RESOLUTION_RANGE'])  #  Input range if set
            else:
                # 2nd or only run
                if self.AUTOCUTOFF and self.highrescutoff > 0.0:
                    # set 'RESOLUTION_RANGE' from 1st run
                    self.aimless.container.controlParameters.RESOLUTION_RANGE.end.set(self.highrescutoff)
                else:
                    self.aimless.container.controlParameters.copyData (\
                      self.container.controlParameters,
                      ['RESOLUTION_RANGE'])  #  Input range if set

            print("**starting aimless")
            self.aimless.finished.connect(self.process_post_aimless)
            self.aimless.process()
        except Exception as e:
            # failed in Aimless
            print("Aimless error")
            self.appendErrorReport(CPluginScript,39,str(e))
            self.fatalError = [202, 'Aimless failed', status]
            self.process_finish(CPluginScript.FAILED)
            return

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    @QtCore.Slot(dict)
    def process_post_aimless(self,status):
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"); sys.stdout.flush()
        print('process_post_aimless', status)
        if status.get('finishStatus') == CPluginScript.FAILED:
            print("**aimless failed")
            self.fatalError = [202, 'Aimless failed', status]
            self.process_finish(CPluginScript.FAILED)
            self.reportStatus(status)
            import traceback
            traceback.print_exc()
            return

        print("** process_post_aimless self.aimlessruncount self.aimlessruncount:",
              self.aimlessruncount)

        # aimless has finished, and produced an output file for each dataset
        self.ndatasets = len(self.aimless.container.outputData.MTZMERGEDOUT)

        #  Check for Aimless indicating complete disaster
        docheck = False
        if self.AUTOCUTOFF:
            if self.aimlessruncount == 1:
                docheck = True
        else:
           if self.aimlessruncount == 2:
               docheck = True
        # docheck True if first or only Aimless
        
        if docheck:
            # Test for complete hopelessness
            OK = True
            self.aimless1xml =\
              CCP4Utils.openFileToEtree(self.aimless.makeFileName( 'PROGRAMXML' ) )
            message = self.checkaimlessresult(self.aimless1xml)

            self.disastermessage = None
            if message is not None:
                self.disastermessage = message
                self.fatalError = [205, 'Very poor data', status]
                # Exit pipeline!!!
                self.process_finish(CPluginScript.UNSATISFACTORY)

        if self.aimlessruncount == 2:
            # Final or only run, do phaser analysis
            print("**end of final aimless")
            self.process_phaseranalysis()
        else:
            # Optional 1st Aimless run to get resolution limit
            # Save XML for later use
            print("**1st Aimless")

            try:
                self.aimless1xml = CCP4Utils.openFileToEtree(self.aimless.makeFileName( 'PROGRAMXML' ) )
            except:
                self.aimless1xml = ET.Element('AIMLESS')

            # Run phaser_analysis to get resolution limit estimate
            # Return:
            #   anygood     True if data go to edge, so no 2nd Aimless run needed
            #   maxres      Resolution cutoff for 2nd Aimless run if wanted
            first = True
            anygood, maxres = self.analyseResolution()

            if anygood:
                self.cutoffdone = False
                self.highrescutoff = -1.0  # flag that 2nd Aimless was not done
                self.process_phaseranalysis()  # go directly to phaser_analysis

            # Get high resolution limit for 2nd run,
            # choose the best if > 1 datasets
            self.highrescutoff = maxres
            self.cutoffdone = True
            print("Resolution limit from 1st Aimless run ", self.highrescutoff)
        
        return  # for 2nd Aimless run
        
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process_phaseranalysis(self):
        # If after 2nd or only Aimless run, run phaser_analysis on Aimless
        # If after 2nd or only Aimless run, run phaser_analysis on Aimless
        # output files
        # For autocutoff mode, if the 1st Aimless run was accepted, then
        # there is no need to run phaser_analysis again
        if self.aimlessruncount == 2:
            # After second or only Aimless run
            if self.doPhaserAnalysis:
                self.analyseAgain()

        self.process_cycle_ctruncate()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process_cycle_ctruncate(self):
      print('process_cycle_ctruncate')

      self.ndatasets_processed = 0
      self.ndatasets_failed = 0

      import shutil
      # Copy any unmerged files to main directory
      nunmerged = len(self.aimless.container.outputData.MTZUNMERGEDOUT)
      #print "Nunmerged",nunmerged
      out = self.container.outputData

      for unmergedfile in self.aimless.container.outputData.MTZUNMERGEDOUT:
          # Loop unmerged output files if any, one per dataset
          # filename in subjob
          fname = os.path.split(str(unmergedfile))[1]

          # File content to get datasetName
          from core.CCP4XtalData import CUnmergedDataContent
          unmergedcontent = CUnmergedDataContent()
          unmergedcontent.loadFile(unmergedfile.fullPath)
          #print "CONTENTS:",unmergedcontent.CONTENTS
          dname = str(unmergedcontent.datasetName)
          # Before August 2017: xname is wrongly swapped with pname,
          #   fixed, but don't use it for now
          xname = unmergedcontent.crystalName

          #print "WD:", self.getWorkDirectory()
          # Add to UNMERGEDOUT
          out.UNMERGEDOUT.append(out.UNMERGEDOUT.makeItem())
          outfile = os.path.join(self.workDirectory,fname)
          out.UNMERGEDOUT[-1].setFullPath(outfile)
          out.UNMERGEDOUT[-1].annotation = \
             'Scaled unmerged data: ' + dname
          #print "outfile",out.UNMERGEDOUT[-1]
          #print out.UNMERGEDOUT[-1].fullPath
          shutil.move(str(unmergedfile), str(out.UNMERGEDOUT[-1].fullPath))
          
      self.collectedCtruncateEtree = ET.Element('CTRUNCATES')
      for file in self.aimless.container.outputData.MTZMERGEDOUT:
         self.process_ctruncate(file)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def process_ctruncate(self,infile):

      print("### Running ctruncate on file %s ###" % infile)

      if not os.path.exists(infile.fullPath.get()):
         self.process_post_ctruncate(CPluginScript.FAILED)

      else:
        hkloutList = self.container.outputData.HKLOUT
        hkloutList.append(hkloutList.makeItem())
        filePath = os.path.join(self.workDirectory,'HKLOUT_'+str(self.ndatasets_processed)+'-observed_data.mtz')
        self.ctruncateOutputDataObject = hkloutList[-1]
        hkloutList[-1].setFullPath(filePath)
        from core import CCP4XtalData
        hkloutList[-1].contentFlag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
        try:
          name =  os.path.splitext(os.path.split(self.container.inputData.UNMERGEDFILES[len(hkloutList)-1].file.__str__())[1])[0]
          hkloutList[-1].annotation = 'Observed dataset '+str(len(hkloutList))+' from '+name       
        except:
          hkloutList[-1].annotation = 'Observed dataset '+str(len(hkloutList))

        # Run ctruncate to generate amplitudes
        self.ctruncate = self.makePluginObject('ctruncate')   
        self.ctruncate.container.inputData.HKLIN = infile
        # We assume standard labels from aimless
        self.ctruncate.container.inputData.ISIGI.I = "IMEAN"
        self.ctruncate.container.inputData.ISIGI.SIGI = "SIGIMEAN"
        self.ctruncate.container.inputData.ISIGIanom.Ip = "I(+)"
        self.ctruncate.container.inputData.ISIGIanom.SIGIp = "SIGI(+)"
        self.ctruncate.container.inputData.ISIGIanom.Im = "I(-)"
        self.ctruncate.container.inputData.ISIGIanom.SIGIm = "SIGI(-)"

        self.ctruncate.container.controlParameters.copyData(self.container.controlParameters,['NRES'])
        self.ctruncate.container.controlParameters.OUTPUTMINIMTZ = True
        self.ctruncate.container.controlParameters.OUTPUT_INTENSITIES = True

        self.ctruncate.container.outputData.OBSOUT.setFullPath(filePath)

        self.connectSignal(self.ctruncate,'finished',self.process_post_ctruncate)
        self.ctruncate.process()
        print("DONE Running ctruncate on file %s ###" % infile)
      
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    @QtCore.Slot(dict)
    def process_post_ctruncate(self,status):
      # after each ctruncate run

      imeanoutPath = str(self.ctruncate.container.outputData.OBSOUT1)
      imeanoutList = self.container.outputData.IMEANOUT
      imeanoutList.append(imeanoutList.makeItem())
      msg = 'Mean intensities for twin refinement'
      if self.ndatasets > 1:
        msg += ' (dataset ' + str(len(self.container.outputData.IMEANOUT)) + ')'
      imeanoutList[-1].setFullPath(imeanoutPath)
      imeanoutList[-1].annotation = msg

      self.ndatasets_processed += 1
      if status.get('finishStatus') == CPluginScript.FAILED:
         self.ndatasets_failed += 1
        
      print("Datasets processed so far: ",self.ndatasets_processed)
      print("Datasets failed so far: ",self.ndatasets_failed)

      if self.ndatasets_failed > 0:
          self.fatalError = [203, 'Ctruncate failed', status]
          self.process_finish(CPluginScript.FAILED)
          return
            
      #Catenate Ctruncate output
      latestCTruncateEtree = CCP4Utils.openFileToEtree(self.ctruncate.makeFileName('PROGRAMXML'))
      self.collectedCtruncateEtree.append(latestCTruncateEtree)
      # the last dataset so far
      crystalDatasedId = latestCTruncateEtree.findall('.//CrystalDatasetId')[-1].text 
      self.ctruncateOutputDataObject.annotation = crystalDatasedId
      #print "crystalDatasedId",crystalDatasedId
        
      #
      # this is called by each ctruncate process
      # on the last time, proceed to cad or terminate
      ## todo: cad step not implemented yet
      if self.ndatasets_processed == self.ndatasets:
         
         #for indx in range(len(self.container.outputData.HKLOUT)):
         #    self.container.outputData.HKLOUT[indx].annotation = 'Observed data set '+str(indx+1)
             
         self.process_finish(CPluginScript.SUCCEEDED)

    #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    def process_finish(self,status):
      import os,shutil

      print("process_finish", status)
      xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
      xmlroot = ET.Element("AIMLESS_PIPE")

      # fatalError may contain a list of [errornumber, message, status]
      #  for failures in Pointless or Aimless
      if self.fatalError is not None:
          errorxml = ET.Element('PIPELINE_ERROR')
          errorxml.text = self.fatalError[1]
          xmlroot.append(errorxml)

      # print("HKLIN_FORMAT")
      if self.container.inputData.HKLIN_FORMAT == "MMCIF":
          print("** mmcif")
          importXML = ET.Element('IMPORT_LOG')  # information about the import step
          self.makeMmcifXML(importXML)
          xmlroot.append(importXML)
          
      try:
          pointlessxml = CCP4Utils.openFileToEtree(self.pointless.makeFileName( 'PROGRAMXML' ) )
      except:
          pointlessxml = ET.Element('POINTLESS')
      xmlroot.append (pointlessxml)

      try:
          aimlessxml = CCP4Utils.openFileToEtree(self.aimless.makeFileName( 'PROGRAMXML' ) )
      except:
          aimlessxml = ET.Element('AIMLESS')

      haveaimless = False
      if hasattr(self,"aimless"):  # True if self.aimless exists
          # Add ONLYMERGE flag if no scaling, for report
          if self.aimless.container.controlParameters.SCALING_PROTOCOL == 'ONLYMERGE':
              eom = ET.Element('ONLYMERGE')
              eom.text = "True"
              aimlessxml.append(eom)
          haveaimless = True
      else:
          print('** no aimless object')

      if haveaimless:
          if self.disastermessage is not None:
              self.addElement(aimlessxml, "DisasterMessage", self.disastermessage)
          self.aimlesstwice = False
          if self.highrescutoff > 0.0 and self.aimless1xml is not None:
              # Aimless has been run twice, add in XML for first run
              self.aimless1xml.tag = 'AIMLESS1'
              xmlroot.append (self.aimless1xml)
              self.aimlesstwice = True

          xmlroot.append (aimlessxml)
          # Add in Autocutoff element if present
          if self.AUTOCUTOFF:
              print("** adding autocutoff")
              xmlroot.append(self.autocutoffxml)
            
      # Add in PHASER_ANALYSES
      # 1st run if done
      if self.phaser_analysis1xml is not None:
          if self.highrescutoff > 0.0:
              self.phaser_analysis1xml.tag = 'PHASER_ANALYSES1'
          xmlroot.append (self.phaser_analysis1xml)
      # 2nd or sole run
      if self.phaser_analysisxml is not None:
          xmlroot.append (self.phaser_analysisxml)

      try:
          xmlroot.append (self.collectedCtruncateEtree)
      except:
          pass

      # ANALYSIS_MODE is set for import_merge task, to suppress FreeR generation
      if (self.fatalError is None) and not self.container.controlParameters.ANALYSIS_MODE:
        self.freerflag = None
        freeStatus, freerReportXML = self.runFreerflag()
        
        #Add freerflag status to the xmlout - this can be the usual SUCCEEDED,FAILED or
        #is 100 if failure in the mtzjoin stage (likely due to space group mismatch)
        e = ET.Element('FREERFLAG')
        e1= ET.Element('status')
        e1.text=str(freeStatus)
        e.append(e1)
        if freerReportXML is not None:
            e.append(freerReportXML)
        if self.freerflag is not None:
            freerxml = self.freerflag.getXML()
            if freerxml is not None:
                e.append(freerxml)
        xmlroot.append(e)

      newXml = ET.tostring(xmlroot,encoding='unicode')
      xmlfile = open( xmlout, "w" )
      xmlfile.write(newXml)
      xmlfile.close()
      if self.container.outputData.XMLOUT.isSet(): # i.e. not import_merged
        with open(str(self.container.outputData.XMLOUT.fullPath), 'w') as f:
            f.write(newXml)
        self.container.outputData.XMLOUT.annotation = 'Pipeline XML'
      # Populate the performance indicator
      if haveaimless:
          try:
              eList = xmlroot.findall('AIMLESS/ReflectionFile/SpacegroupName')
              #print('aimless_pipe.process_finish spacegroup',str(eList[0].text).strip())
              spGp = self.container.outputData.PERFORMANCE.spaceGroup
              if len(eList)>0: spGp.set(spGp.fix(str(eList[0].text).strip()))
          except:
              pass

          eList = xmlroot.findall('AIMLESS/Result/Dataset/RmeasOverall/Overall')
          if len(eList)>0: self.container.outputData.PERFORMANCE.rMeas.set(float(str(eList[0].text).strip()))
          eList = xmlroot.findall('AIMLESS/Result/Dataset/ResolutionHigh/Overall')
          if len(eList)>0: self.container.outputData.PERFORMANCE.highResLimit.set(float(str(eList[0].text).strip()))

      # Concatenate log files
      logs = ['job_1/log.txt','job_2/log.txt','job_3/log.txt','job_4/log.txt',
              'job_5/log.txt']
      if haveaimless and self.aimlesstwice:
          logs = logs + ['job_6/log.txt','job_7/log.txt']
      
      files = []
      for log in logs:
          fulllog = os.path.join(self.workDirectory,log)
          if os.path.isfile(fulllog):
              #print("Appending log ",log)
              files.append(fulllog)
      import fileinput
      outfile = os.path.join(self.workDirectory,'log.txt')
      fout = open(outfile, 'w')
      for line in fileinput.input(files):
          fout.write(line)

      # Emit finished signal with whatever status has been emitted 
      # by the last wrapper
      self.cleanup()
      
      if status == CPluginScript.FAILED:
          self.reportStatus(CPluginScript.UNSATISFACTORY)
      else:
          self.reportStatus(status)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def makeMmcifXML(self, containerXML):
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
            
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def cleanup(self):
        # Clean up some files which are not always removed by
        # the PURGE mechanism
        print("*cleanup*")
        try:
            import glob
            wd = self.workDirectory
            filelist = glob.glob(wd+"/*/*.xmgr")
            for fn in filelist:
                os.remove(fn)
                # print(fn,"deleted")
            print("*cleanup done*")
        except:
            print("*cleanup failed*")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def analyseResolution(self):
        # Run phaser_analysis to get resolution limit
        # infothreshold is limit on information content
        # Make self.phaser_analysis1xml XML element
        # Make self.autocutoffxml XML element
        # return anygood, flag for all datasets to maximum resolution
        # and maximum resolution
        print("analyseResolution")
        infothreshold = str(self.container.controlParameters.INFOCONTENTTHRESHOLD)

        self.phaser_analysis1xml = ET.Element('PHASER_ANALYSES')

        # start XML element for recording the results of the 1st Aimless 
        self.autocutoffxml = ET.Element('Autocutoff')
        
        # extract high resolution limits, allowing for multiple datasets
        datasetResultList =  self.aimless1xml.findall('Result/Dataset')

        # We want to put into the autocutoffxml block information from both
        # Aimless and Phaser, so we need to match them here
        # If there are multiple datasets, then the Aimless XML has information
        # for each dataset, and there are multiple output files (identified
        # by Dname) which need to be processed  by phaser_analysis.
        # The entries in Aimless XML and the datasets in the list are
        # probably in the same order, but in case they are not,
        # we match the Dnames here
        
        # Extract information from Aimless xml for each dataset
        dsetinfo = {}
        for dataset in datasetResultList:
            pxdname = dataset.attrib['name']
            dname = pxdname.split('/')[2]
            inputresolution = \
                 dataset.find('ResolutionHigh/Overall').text.strip()
            estimates = dataset.findall('ResolutionLimitEstimate[@type="CChalf"]')
            reslimit = 0.0
            for estimate in estimates:  # Pick out Overall limit
                direction = estimate.find('Direction').text
                if direction =='Overall':
                    reslimit = float(estimate.find('MaximumResolution').text)
                    message = ''
                    if estimate.find('Message') is not None:
                      message =  estimate.find('Message').text

            dsetinfo[dname] = {'pxdname': pxdname, 
                               'InputResolution': inputresolution,
                               'ResolutionEstimateCChalf': str(reslimit),
                               'Message': message}
        
        if len(dsetinfo) != self.ndatasets:
            print("!! Help dsetinfo wrong")
        
        maxres = 10000.0
        numalltomaxres = 0

        # Loop output datasets, run phaser_analysis for each one
        for filename in self.aimless.container.outputData.MTZMERGEDOUT:
            if self.ndatasets > 1:
                # Multiple datasets, extract Dname from filename to match dsetinfo
                #  filenames are "HKLOUT_dname.mtz"
                basename = str(filename).split('/')[-1]
                dname = basename.split('_')[1].split('.')[0]
            else:
                dname = list(dsetinfo.keys())[0]  # Key from sole dict entry
            pxdname = dsetinfo[dname]['pxdname']
            self.phaser_analysis = self.makePluginObject('phaser_analysis')
            # Threshold
            self.phaser_analysis.container.controlParameters.copyData \
              (self.container.controlParameters, ['INFOCONTENTTHRESHOLD'])
            # Run phaser
            self.phaser_analysis.process(filename=filename, pxdname=pxdname)
            allOK, reslimit = self.phaser_analysis.getresolutionlimit()
            if allOK > 0:
                numalltomaxres += 1
                reslimit = float(dsetinfo[dname]['InputResolution'])
            maxres = min(maxres, reslimit)
            #print("***dname,reslimit,maxres",dname,reslimit,maxres)
            # Construct XML element dset for Autocutoff
            # <InputResolution>
            # <ResolutionEstimateCChalf>  limit found by CChalf metric
            # <ResolutionEstimate>  limit found by IC metric
            # <Message>  only if data go to maximum resolution
            # <InformationThreshold> information thresholdfor cutoff
            dset = ET.Element('Dataset')
            dset.set('name', pxdname)
            self.addElement(dset, 'InputResolution',
                            dsetinfo[dname]['InputResolution'])
            self.addElement(dset, 'ResolutionEstimateCChalf', 
                            dsetinfo[dname]['ResolutionEstimateCChalf'])
            self.addElement(dset, 'ResolutionEstimate',
                            "{:6.2f}".format(reslimit))
            message =  dsetinfo[dname]['Message']
            if message != '':
                self.addElement(dset, 'Message', message)
            self.autocutoffxml.append(dset)

            #  phaser_analysis XML
            phaseranalysisxml = self.phaser_analysis.getXML()
            # add InputResolution into Phaser XML since the file resolution
            # from Phaser is unreliable
            self.addResolutiontoPhaser(phaseranalysisxml,
                                       dsetinfo[dname]['InputResolution'])
            phaseranalysisxml.set('name', pxdname)
            self.phaser_analysis1xml.append(phaseranalysisxml)
                        
        # end loop datasets

        anygood = False
        cutoff = 'Yes'
        if numalltomaxres > 0:
            # at least one dataset goes to the maximum resolution
            #  so no need to run second Aimless
            anygood = True
            cutoff = 'No'
        self.addElement(self.autocutoffxml, 'Cutoff', cutoff)

        self.addElement(self.autocutoffxml, 'InformationThreshold', str(infothreshold))            
        # Best resolution
        if maxres < 1000.0:
            self.addElement(self.autocutoffxml,
                            'ResolutionCutoff', "{:6.2f}".format(maxres))
        else:
            maxres = -1.0

        return anygood, maxres

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def analyseAgain(self):
        # Run phaser_analysis after final Aimless, if needed
        # infothreshold is limit on information content for analysis
        # Make self.phaser_analysisxml XML element
        infothreshold = str(self.container.controlParameters.INFOCONTENTTHRESHOLD)

        if self.aimless1xml is None:
            # If Aimless is only run once
            self.aimless1xml = CCP4Utils.openFileToEtree(self.aimless.makeFileName( 'PROGRAMXML' ) )

        datasetResultList =  self.aimless1xml.findall('Result/Dataset')
        dnames = []
        pxdnames = {}
        inputresolutions = {}
        for dataset in datasetResultList:
            pxdname = dataset.attrib['name']
            dname = pxdname.split('/')[2]
            dnames.append(dname)
            pxdnames[dname] = pxdname
            inputresolutions[dname] = dataset.find('ResolutionHigh/Overall').text.strip()

        # This will contain a PHASER_ANALYSIS block for each dataset
        self.phaser_analysisxml = ET.Element('PHASER_ANALYSES')

        # Loop output datasets, run phaser_analysis for each one
        for filename in self.aimless.container.outputData.MTZMERGEDOUT:
            if self.ndatasets > 1:
                # Multiple datasets, extract Dname from filename to match dsetinfo
                #  filenames are "HKLOUT_dname.mtz", dname may contain "_"
                basename = str(filename).split('/')[-1]
                i = basename.find('_')
                if i < 0: i = 0
                j = min(i+1,len(basename)-1)
                dname = basename[j:].split('.')[0]
            else:
                dname = dnames[0]
            pxdname = pxdnames[dname]
        
            self.phaser_analysis = self.makePluginObject('phaser_analysis')
            sys.stdout.flush()
        
            # Threshold
            self.phaser_analysis.container.controlParameters.copyData(self.container.controlParameters, ['INFOCONTENTTHRESHOLD'])
            rstatus = self.phaser_analysis.process(filename=filename, pxdname=pxdname)
            print("after Phaser, status", rstatus, CPluginScript.FAILED)
            if rstatus  == CPluginScript.FAILED:
                message = self.phaser_analysis.resultObject.ErrorMessage()
                phaserelement = ET.Element('PHASER_ANALYSIS')
                phaserfailXML = ET.SubElement(phaserelement,'PhaserFailMessage')
                phaserfailXML.text = message
                self.phaser_analysisxml.append(phaserelement)
                return
            
            allOK, reslimit = self.phaser_analysis.getresolutionlimit()

            #  phaser_analysis XML
            phaseranalysisxml = self.phaser_analysis.getXML()
            # add InputResolution into Phaser XML since the file resolution
            # from Phaser is unreliable
            self.addResolutiontoPhaser(phaseranalysisxml, inputresolutions[dname])
            phaseranalysisxml.set('name', pxdname)
            self.phaser_analysisxml.append(phaseranalysisxml)
        # end loop datasets

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    def addResolutiontoPhaser(self, xmlblock, inputresolution):
        # add InputResolution into Phaser XML since the file resolution
        # from Phaser is unreliable
        res = xmlblock.find('Analysis/Resolution')
        self.addElement(res,'InputResolution', inputresolution),

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    def checkaimlessresult(self, xmlnode):
        # Returns None if OK, else error message
        # Other wise add error message to xmlnode

        message = None
        result = xmlnode.findall("Result")
        datasets = result[0].findall("Dataset")
        if len(datasets) == 0:
            print("no datasets")
            return message

        # Limits for inner shell statistics
        cchalfdisaster = self.container.controlParameters.CCHALFDISASTER
        multiplicitydisaster = \
                   self.container.controlParameters.MULTIPLICITYDISASTER
        
        OK = True
        innercchalf = "0.000"
        multiplicity = "0.0"
        for dataset in datasets:
            # Get statistics for inner shell
            innercchalf = dataset.findall("CChalf/Inner")[0].text

            #  Check that there enough data, get multiplicity
            multiplicity = dataset.findall("Multiplicity/Inner")[0].text

            # Test against disaster limits
            #    Disaster if CC(1/2) very low in inner shell,
            #    and reasonable multiplicity
            #    but allow CChalf == 0.000
            if  (not '0.000' in innercchalf) and \
               (float(innercchalf) < float(cchalfdisaster)):
                # fails CC(1/2) test
                if float(multiplicity) > float(multiplicitydisaster):
                    # ... and good multiplicity
                    OK = False

        if OK:
            return message  # no problem

        #  Bail out if Aimless has reported that the data are really poor
        message = "SERIOUS PROBLEM: the data are VERY poor"
        message += ", as identified"+\
                  " by poor statistics in the inner resolution shell"
        message += " CC(1/2) = "+str(innercchalf)+", less than threshold "+\
                   str(cchalfdisaster)
        message += ", with sufficient multiplicity of "+str(multiplicity)+\
                   ", greater than threshold "+str(multiplicitydisaster)
        
        return message


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def addElement(self, containerXML, elementname, elementtext):
        e2 = ET.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def runFreerflag(self):
      print("runFreerflag")

      complete = False
      if self.container.inputData.FREERFLAG.isSet():
          # Check for cell compatibility with input FreeR set
          # Extending existing FreeR
          tolerance = None # use default
          cellcheck = \
             CellCheck(self.container.outputData.HKLOUT[0].fileContent,
                       self.container.inputData.FREERFLAG.fileContent,
                       tolerance)
          cellsAreTheSame, freerReportXML = cellcheck.checks()
          
          complete = True
          if (not self.container.controlParameters.OVERRIDE_CELL_DIFFERENCE) and \
                 (not cellsAreTheSame['validity']):
              # Not compatible
              status = CPluginScript.FAILED
              print('Cells not compatible', status)
              return status, freerReportXML

      from wrappers.freerflag.script import freerflag
      freerReportXML = None
      self.freerflag = self.makePluginObject('freerflag')
      self.freerflag.container.inputData.F_SIGF = self.container.outputData.HKLOUT[0]
      self.freerflag.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG
      self.freerflag.container.controlParameters.COMPLETE = self.container.controlParameters.COMPLETE
      self.freerflag.container.controlParameters.FRAC = \
                          self.container.controlParameters.FREER_FRACTION
      self.freerflag.container.controlParameters.CUTRESOLUTION = \
                          self.container.controlParameters.CUTRESOLUTION
      filePath = os.path.join(self.workDirectory,'FREERFLAG.mtz')
      self.freerflag.container.outputData.FREEROUT.setFullPath(filePath)
      # Set optional COMPLETE flag: GEN_MODE = 'GEN_NEW' or 'COMPLETE'
      self.freerflag.container.controlParameters.GEN_MODE.set('GEN_NEW')

      if complete:
          self.freerflag.container.controlParameters.GEN_MODE.set('COMPLETE')
      else:
          # Generating new FreeR
          # Fraction for FreeR if specified
          if self.container.controlParameters.FREER_FRACTION.isSet():
              frac = self.container.controlParameters.FREER_FRACTION
              # print("FreeR fraction", type(frac), frac)
              freerReportXML = ET.Element('FreeRfraction')
              freerReportXML.text = frac.__str__()

      # Run the FreeR wrapper
      status = self.freerflag.process()
      if status == CPluginScript.SUCCEEDED:
        self.container.outputData.FREEROUT.setFullPath(filePath)
        # print "FREEROUT",self.container.outputData.FREEROUT.fileContent
        if self.container.outputData.FREEROUT.fileContent.spaceGroup.isSet():
            sgname = self.container.outputData.FREEROUT.fileContent.spaceGroup.__str__()
        else:
            sgname = 'Unk'

        title = "New Freer;"
        if complete:
            title = 'Extended freeR;'

        highresFRformatted = "%7.2f" % float(self.container.outputData.FREEROUT.fileContent.resolutionRange.high)
        title += ' Spg:'+str(sgname).strip()+';Resln:'+highresFRformatted.strip() + "A;"
        try:
          title = title + "Cell:"+self.container.outputData.FREEROUT.fileContent.cell.guiLabel()
        except Exception as e:
          print('Error writing cell parameters',e)

        #print("Title:",title)
        self.container.outputData.FREEROUT.annotation.set(title)
      else:
        nMergeFail = self.freerflag.getErrorReport().count(cls=CPluginScript,code=32)
        # print("nMergeFail", nMergeFail)
        if nMergeFail>0: return 100, freerReportXML
          
      return status, freerReportXML


# ----------------------------------------------------------------------
# Function to return list of names of exportable MTZ(s)
def exportJobFile(jobId=None,mode=None):
    import os
    from core import CCP4Modules
    from core import CCP4XtalData

    jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
    exportFile = os.path.join(jobDir,'exportMtz.mtz')
    if os.path.exists(exportFile): return exportFile

    childJobs = CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=jobId,details=True)
    #print 'aimless.exportMtz',childJobs
    truncateOut = None
    freerflagOut = None
    for jobNo,subJobId,taskName  in childJobs:
      if taskName == 'ctruncate':
        truncateOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=subJobId,create=False),'HKLOUT.mtz')
        if not os.path.exists(truncateOut): truncateOut = None
      elif taskName == 'freerflag':
        freerflagOut = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False),'FREERFLAG.mtz')
        if not os.path.exists(freerflagOut): freerflagOut = None
    if truncateOut is None: return None
    if freerflagOut is None: return truncateOut

    # print('aimless_pipe.exportJobFile  runCad:',exportFile,[ freerflagOut ])
    
    m = CCP4XtalData.CMtzDataFile(truncateOut)
    #print m.runCad.__doc__   #Print out docs for the function
    outfile,err = m.runCad(exportFile,[ freerflagOut ] )
    print('aimless_pipe.exportJobFile',outfile,err.report())
    return   outfile                                                   
 
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ] ]
                                                
# ---------------------------------------------------------------------
'''
import clipper
import math

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # NOT USED # #
class SameCell:
    "Check if two cells are the same"

    def __init__(self):
      self.ccell1 = None   # clipper cells
      self.ccell2 = None
      self.tolerance = None

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def compareCellsinFileContent(self, filecontent1cell, filecontent2cell,
                                  tolerance=None):
        "Compare cell from two fileContent cell objects"
        """
        Return dictionary of:
          'validity'           True if cells are simlar within resolution of tolerance
          'maximumResolution1' maximum allowed resolution in cell1
          'maximumResolution2' maximum allowed resolution in cell2
          'difference'  average cell difference in A        
          'tolerance'   in A
        """
        if tolerance is not None:
            self.tolerance = tolerance

        cell1 = self.extractCell(filecontent1cell.get())
        cell2 = self.extractCell(filecontent2cell.get())

        return self.computeDiffs(cell1, cell2)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def compareCells(self, gcell1, gcell2, tolerance=None):
        "Compare two cell dictionaries (entry for testing)"
        if tolerance is not None:
            self.tolerance = tolerance
        cell1 = self.extractCell(gcell1)
        cell2 = self.extractCell(gcell2)
        return self.computeDiffs(cell1, cell2)
        
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def computeDiffs(self, cell1, cell2):
        "cell1, cell2 are arrays of cell dimensions"
        # store clipper cells
        self.ccell1 = self.makeCell(cell1[0], cell1[1], cell1[2],
                                   cell1[3], cell1[4], cell1[5])
        self.ccell2 = self.makeCell(cell2[0], cell2[1], cell2[2],
                                   cell2[3], cell2[4], cell2[5])
        return self.computeResults()

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def makeCell(self, a, b, c, alpha, beta, gamma):
        "Return clipper::Cell from cell dimensions, angles in degress or radians"
        return clipper.Cell(clipper.Cell_descr(a, b, c, alpha, beta, gamma))

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def computeResults(self):
        "Compute results, return as dictionary"
        #print "computeResults"
        recipdifference = self.reciprocalcellDifference()
        difference = self.cellDifference()
        if self.tolerance is None:
            self.tolerance = 1.0
        #print "tolerance", self.tolerance
        equals = self.ccell1.equals(self.ccell2, self.tolerance)  # Boolean
        #print "Result:",equals, difference, recipdifference
        """
        Return dictionary of:
          'validity'           True if cells are simlar within resolution of tolerance
          'maximumResolution1' maximum allowed resolution in cell1
          'maximumResolution2' maximum allowed resolution in cell2
          'difference'  average cell difference in A
          'tolerance'   in A
        """
        result = {'validity': equals,
                  'maximumResolution1': recipdifference[0],
                  'maximumResolution2': recipdifference[1],
                  'difference': difference, 'tolerance': self.tolerance}
        return result
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def extractCell(self, celldict):
        "extract cell dimensions, return array"
        cell = []
        for item in ['a','b','c','alpha','beta', 'gamma']:
            cell.append(float(celldict[item]))
        return cell

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def reciprocalcellDifference(self):
        "Get resolution at which ccell2 is different from ccell1, two values"
        fracmat1 = self.ccell1.matrix_frac()
        fracmat2 = self.ccell2.matrix_frac()
        s = 0.0
        for j in range(3):
            for i in range(3):
                s += (fracmat1[i][j] - fracmat2[i][j])**2
        s = math.sqrt(s)
        volume1 = self.ccell1.volume()
        volume2 = self.ccell2.volume()
        return s * math.pow(volume1, 0.66666666667), \
               s * math.pow(volume2, 0.66666666667)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def cellDifference(self):
        "Get resolution at which self.ccell2 is different from self.ccell1"
        orthmat1 = self.ccell1.matrix_orth()
        orthmat2 = self.ccell2.matrix_orth()
        #print "orth1\n", orthmat1.format()
        #print "orth2\n", orthmat2.format()
        s = 0.0
        for j in range(3):
            for i in range(3):
                s += (orthmat1[i][j] - orthmat2[i][j])**2
        return math.sqrt(s)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
'''
 
#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testaimless_pipe(unittest.TestCase):

   def setUp(self):
    from core import CCP4Modules
    self.app = CCP4Modules.QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core import CCP4Modules
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from core import CCP4Modules
     import os

     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test1')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = aimless_pipe(parent=CCP4Modules.QTAPPLICATION(),name='test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'pipelines','aimless_pipe','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testaimless_pipe)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
