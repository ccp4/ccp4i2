from __future__ import print_function

"""
    aimless.py: CCP4 GUI Project
    Copyright (C) 2012 STFC
    """

import os
from ccp4i2.core.CCP4PluginScript import CPluginScript

class aimless(CPluginScript):
    
    TASKNAME = 'aimless'   # Task name - should be same as class name
    TASKTITLE = 'Scale and merge dataset (AIMLESS)' # A short title for gui menu
    TASKVERSION= 0.0               # Version of this plugin
    MAINTAINER = 'pre@mrc-lmb.cam.ac.uk'
    
    # used by the base class startProcess()
    TASKCOMMAND = 'aimless'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None
    COMTEMPLATE = None

    ERROR_CODES = { 201 : { 'description' : 'Aimless program failed?  No output reflection files found. See log file.' }
                    }
    # -----------------------------------------------------------------------
    def makeCommandAndScript(self):

        #print("Aimless makeCommandAndScript")

        par = self.container.controlParameters
        
        self.appendCommandLine(['XMLOUT',str( self.makeFileName( 'PROGRAMXML' ) )])
        
        self.appendCommandLine(['HKLIN',self.container.inputData.UNMERGEDFILE.fullPath])
        
        self.appendCommandLine(['SCALES',self.container.outputData.SCALES.fullPath])
        self.appendCommandLine(['ROGUES',self.container.outputData.ROGUES.fullPath])
        self.appendCommandLine(['NORMPLOT',self.container.outputData.NORMPLOT.fullPath])
        self.appendCommandLine(['ANOMPLOT',self.container.outputData.ANOMPLOT.fullPath])
        self.appendCommandLine(['CORRELPLOT',self.container.outputData.CORRELPLOT.fullPath])
        self.appendCommandLine(['ROGUEPLOT',self.container.outputData.ROGUEPLOT.fullPath])

        if par.REFERENCE_FOR_AIMLESS:
            if par.REFERENCE_DATASET == "XYZ":
                # print "XYZIN",self.container.inputData.XYZIN_REF 
                self.appendCommandLine(['XYZIN',self.container.inputData.XYZIN_REF.fullPath])
            elif par.REFERENCE_DATASET == "HKL":
                # print "HKLIN_REF",self.container.inputData.HKLIN_REF
                self.appendCommandLine(['HKLREF',self.container.inputData.HKLIN_REF.fullPath])

        # output options for unmerged & scalepack files etc
        # maybe not useful yet in the ccp4i2 context
        self.outputOptions()

        # scaling protocol and associated parameters
        #        print "calling scalingProtocol"
        self.scalingProtocol()
        self.scalingDetails()
        
        # reject outliers
        self.rejectOutliers()

        # ANALYSIS  thresholds CC(1/2), I/sigI
        self.analysisParameters()

        # SD correction things
        self.SDcorrection()

        # intensities and partials
        self.intensitiesAndPartials()

        # Resolution
        s = self.resolutionRangeCommand()
        if s != "":
            self.appendCommandScript(s)

        # Explicit run definitions, may include run resolutions
        print("Aimless makeCommandAndScript 13")
        self.runExplicit()

        return CPluginScript.SUCCEEDED
    
    # -----------------------------------------------------------------------
    def resolutionRangeCommand(self):
        print("Aimless resolutionRangeCommand")
        s = ""
        par = self.container.controlParameters
        #print "resolutionRangeCommand",par.RESOLUTION_RANGE
        r1 = par.RESOLUTION_RANGE.start
        r2 = par.RESOLUTION_RANGE.end
        high = 0.0
        low  = 0.0
        if not r1.isSet() and  not r2.isSet():
            # nothing set
            s = ""
        elif r1.isSet() and r2.isSet():
            low  = float(r1)
            high = float(r2)
            if r1 == r2:
                s = "RESOLUTION HIGH %f" % high
            else:
                if low < high:
                    low, high = high, low
                s = "RESOLUTION LOW %f HIGH %f" % (low, high)
        else:
            if not r1.isSet():
                high = r2            
            if not r2.isSet():
                high = r1
            s = "RESOLUTION HIGH %f" % high

        #print "s",s
        return s
    # -----------------------------------------------------------------------
    def scalingProtocol(self):
        print("Aimless scalingProtocol")
        # scaling protocol and associated parameters
        par = self.container.controlParameters

        if (bool(par.ONLYMERGE) == True) or (str(par.SCALING_PROTOCOL) == 'ONLYMERGE'):
            #  Onlymerge, no SCALES definition
            self.setOnlymerge()
            return

        if par.SCALING_PROTOCOL == 'DEFAULT':
            return

        s = 'SCALES '
        if par.SCALING_PROTOCOL == 'CONSTANT':
            s += 'CONSTANT '
            if par.BFACTOR_SCALE:
                s += 'BFACTOR ON '
        elif par.SCALING_PROTOCOL == 'BATCH':
            s += 'BATCH '
        else:
            s += 'ROTATION '
            if par.SCALES_ROTATION_TYPE == 'SPACING':
                s += "SPACING %f " % par.SCALES_ROTATION_SPACING
            else:
                s += " %d " % par.SCALES_ROTATION_NBINS

            if par.SCALING_PROTOCOL == 'SECONDARY':
                s += "SECONDARY %d " % par.SCALES_SECONDARY_NSPHHARMONICS
                    
            if par.SCALES_TILETYPE == 'NONE':
                s += 'NOTILE '
            elif par.SCALES_TILETYPE == 'CCD':
                ntx = -1
                nty = -1
                if par.SCALES_NTILEX.isSet():
                    ntx = int(par.SCALES_NTILEX)
                    nty = ntx
                if par.SCALES_NTILEY.isSet():
                    nty = int(par.SCALES_NTILEY)
                if ntx > 0:
                    s += "TILE %d %d " % (ntx,nty)
                s += 'CCD '

            if par.BFACTOR_SCALE:
                s += 'BROTATION '
                if par.SCALES_BROTATION_TYPE == 'SPACING':
                    s += "SPACING %f " % par.SCALES_BROTATION_SPACING
                else:
                    s += " %d " % par.SCALES_BROTATION_NBINS
            else:
                s += 'BFACTOR OFF '

        # print "Adding Scale command: ", s
        self.appendCommandScript(s)
    
    # -----------------------------------------------------------------------
    def setOnlymerge(self):
        print("Aimless setOnlymerge")
        # Turn off outliers and sd optimisation unless set

        self.appendCommandScript('ONLYMERGE')
        par = self.container.controlParameters
        if not par.SDCORRECTION_OVERRIDE:
            # no explicit SDCORRECTION given
            self.appendCommandScript("SDCORRECTION NOREFINE 1.0 0.0 0.0")

        if not par.OUTLIER_OVERRIDE:
            # No explicit REJECT given
            self.appendCommandScript("REJECT NONE")

    # -----------------------------------------------------------------------
    def rejectOutliers(self):
        print("Aimless rejectOutliers")
        # reject outliers
        par = self.container.controlParameters

        if not par.OUTLIER_OVERRIDE:
            return

        if not par.OUTLIER_EMAX.isDefault():
            self.appendCommandScript("REJECT EMAX %f" % par.OUTLIER_EMAX)

        s = "REJECT "
        if par.OUTLIER_COMBINE.isSet():
            s += "COMBINE "
        else:
            s += "SEPARATE "
            
        if (not par.OUTLIER_SDMAX.isDefault()) or (not par.OUTLIER_SDMAX2.isDefault()):
            s += "%f " % par.OUTLIER_SDMAX
            if par.OUTLIER_SDMAX2.isSet():
                s += "%f " % par.OUTLIER_SDMAX2

        if not par.OUTLIER_SDMAXALL.isDefault():
            sda = float(par.OUTLIER_SDMAXALL)
            if par.OUTLIER_SDMAXALL_ADJUST:
                sda = -abs(sda)
            s += "ALL %f " % sda

        # print "Adding Reject command: ", s
        self.appendCommandScript(s)

    # -----------------------------------------------------------------------
    def analysisParameters(self):
        print("Aimless analysisParameters")
        # Analysis by CC(1/2) and Mn(I/sigI)
        par = self.container.controlParameters

        s = "ANALYSIS "

        # Pipeline default may be different from program, so always put this in
        s += "CCMINIMUM %f" % par.CCMINIMUM
        
        if not par.CCANOMMINIMUM.isDefault():
            s += " CCANOMMINIMUM %f " % par.CCANOMMINIMUM

        if not par.ISIGMINIMUM.isDefault():
            s += " ISIGMINIMUM %f " % par.ISIGMINIMUM

        self.appendCommandScript(s)

    # -----------------------------------------------------------------------
    # SD correction things
    def SDcorrection(self):
        print("Aimless SDcorrection")
        par = self.container.controlParameters
        if not par.SDCORRECTION_OVERRIDE:
            return

        s = "SDCORRECTION "
        if par.SDCORRECTION_REFINE:
            s += "REFINE "
        else:
            s += "NOREFINE "

        if par.SDCORRECTION_REFINE:
            if par.SDCORRECTION_OPTIONS.isSet():
                if par.SDCORRECTION_OPTIONS == 'INDIVIDUAL':
                    s += 'INDIVIDUAL '
                elif par.SDCORRECTION_OPTIONS == 'SAME':
                    s += 'SAME '
                elif par.SDCORRECTION_OPTIONS == 'SIMILAR':
                    s += 'SIMILAR '
                    # print "SDcorrection1", s

            if par.SDCORRECTION_OPTIONS == 'SIMILAR':
                sdsim = False
                if par.SDCORRECTION_SIMILARITY_SDFAC.isSet():
                    s += "%7.3f "%par.SDCORRECTION_SIMILARITY_SDFAC
                    sdsim = True
                if par.SDCORRECTION_SIMILARITY_SDB.isSet():
                    s += "%7.3f "%par.SDCORRECTION_SIMILARITY_SDB
                if sdsim:   # Always sdAdd if sdFac parameter given
                    s += "%7.3f "%par.SDCORRECTION_SIMILARITY_SDADD

            if par.SDCORRECTION_FIXSDB:
                s += 'FIXSDB '

        self.appendCommandScript(s)
        #print "Adding SDCORRECTION command: ", s
            
        if par.SDCORRECTION_REFINE and par.SDCORRECTION_TIESDB_SD.isSet():
            s = 'SDCORRECTION TIE SdB 0.0 %7.2f '%par.SDCORRECTION_TIESDB_SD 
            self.appendCommandScript(s)

        if par.SDCORRECTION_REFINE and (not par.SDCORRECTION_DAMP.isDefault()):
            s = "SDCORRECTION DAMP %7.3f "% par.SDCORRECTION_DAMP
            self.appendCommandScript(s)

        if (not par.SDCORRECTION_REFINE) and par.SDCORRECTION_SET:
            s = 'SDCORRECTION '+\
            "%7.4f " % par.SDCORRECTION_SDFAC +\
            "%7.4f " % par.SDCORRECTION_SDB +\
            "%7.4f " % par.SDCORRECTION_SDADD
            self.appendCommandScript(s)

    # ----------------------------------------------------------------------------
    def intensitiesAndPartials(self):
        print("Aimless intensitiesAndPartials")
        par = self.container.controlParameters
        if not par.INTENSITIES_OVERRIDE:
            return

        if par.INTENSITIES_OPTIONS.isSet():
            s = ''
            if par.INTENSITIES_OPTIONS == 'COMBINE':
                s = 'INTENSITIES COMBINE'
            elif par.INTENSITIES_OPTIONS == 'PROFILE':
                s = 'INTENSITIES PROFILE'
            elif par.INTENSITIES_OPTIONS == 'SUMMATION':
                s = 'INTENSITIES SUMMATION'
            if s != '':
                self.appendCommandScript(s)

        if par.PARTIALS_TEST.isSet():
            if par.PARTIALS_TEST:
                s = 'PARTIALS TEST '+\
                    "%7.3f "% (par.PARTIALS_FRACLOW) +\
                    "%7.3f "% (par.PARTIALS_FRACHIGH)
                self.appendCommandScript(s)

        if par.PARTIALS_CHECK:
            s = 'PARTIALS CHECK '
        else:
            s = 'PARTIALS NOCHECK '
        self.appendCommandScript(s)
                
        if par.PARTIALS_SCALE:
            s = 'PARTIALS SCALE ' +\
                "%7.3f "%par.PARTIALS_SCALE_MIN
            self.appendCommandScript(s)

        if par.ACCEPT_OVERLOADS:
            # overload flag
            self.appendCommandScript(['KEEP OVERLOADS'])

        if par.ACCEPT_EDGES:
            # edge flag
            self.appendCommandScript(['KEEP EDGE'])

        if par.ACCEPT_XDS_MISFITS:
            # XDS MISFIT (outlier) flag
            self.appendCommandScript(['KEEP MISFIT'])
            
    # -----------------------------------------------------------------------
    def scalingDetails(self):
        print("Aimless scalingDetails")
        par = self.container.controlParameters

        print("scalingDetails parallel", par.PARALLEL,  par.PARALLEL_MODE)

        if par.PARALLEL:
            print("scalingDetails parallel2", par.PARALLEL,  par.PARALLEL_MODE)

            if par.PARALLEL_MODE == 'AUTO':
                self.appendCommandScript('REFINE PARALLEL AUTO')
                print("scalingDetails parallel auto")
            elif par.PARALLEL_MODE == 'NUMBER':
                self.appendCommandScript('REFINE PARALLEL '+" %4d"%par.NPROC)
                print("scalingDetails parallel number")
            elif par.PARALLEL_MODE == 'FRACTION':
                self.appendCommandScript('REFINE PARALLEL '+" %6.3f"%par.FRACPROC)
                print("scalingDetails parallel fraction")

        if par.REFERENCE_FOR_AIMLESS:
            if par.REFINE_REFERENCE:
                self.appendCommandScript('REFINE REFERENCE')

        if not par.SCALING_DETAILS:
            return

        if par.CYCLES_FLAG:
            self.appendCommandScript('REFINE CYCLES '+"%4d "%par.CYCLES_N)

        if par.SELECT2:
            par.SELECT1 = True   # cannot have 2 without 1

        if par.SELECT1:
            s = 'REFINE SELECT '+"%7.3f " % par.SELECT_IOVSDMIN
            if par.SELECT2:
                s += "%7.3f " % par.SELECT_EMIN
            self.appendCommandScript(s)
            
        s = "TIE "
        if par.TIE_ROTATION:
            s += 'ROTATION '+ "%7.4F " % par.TIE_ROTATION_SD
        if par.TIE_BFACTOR:
            s += 'BFACTOR '+ "%7.4F " % par.TIE_BFACTOR_SD
        if par.TIE_SURFACE:
            s += 'SURFACE '+ "%7.4F " % par.TIE_SURFACE_SD
        if par.TIE_BZERO:
            s += 'ZEROB '+ "%7.4F " % par.TIE_BZERO_SD
        # Should do TIE TILE [5 numbers] also

        if s != "TIE ":
             self.appendCommandScript(s)

    # -----------------------------------------------------------------------
    def runExplicit(self):
        print("Aimless runExplicit")
        par = self.container.controlParameters
        if par.RUN_MODE != 'BYRANGE':
            return

        #print "byrange", par.RUN_BATCHLIST.isSet()
        #print "lenRBL", len(par.RUN_BATCHLIST)
        nruns = len(par.RUN_BATCHLIST)
        if nruns > 0:
            for i in range(nruns):
                runrange = par.RUN_BATCHLIST[i]
                #print( "runrange", runrange)
                #print "runrange", runrange.runNumber, runrange.batchRange0, runrange.batchRange1
                s = 'RUN ' + "%3d " %  runrange.runNumber +\
                    " BATCH %5d " % runrange.batchRange0 + " TO %5d " % runrange.batchRange1
                self.appendCommandScript(s)
                
                if runrange.resolution.isSet():
                    s = 'RESOLUTION RUN ' + "%3d " %  runrange.runNumber +\
                        " HIGH %f" % runrange.resolution
                    self.appendCommandScript(s)
        
    # -----------------------------------------------------------------------
    def outputOptions(self):

        print("outputOptions:", self.container.outputData.HKLOUT_BASENAME)        

        if not self.container.outputData.HKLOUT_BASENAME.isSet():
            self.container.outputData.HKLOUT_BASENAME = os.path.join(self.getWorkDirectory(),"HKLOUT")
        par = self.container.controlParameters

        # par.OUTPUT_MODE is a CList of CStrings, each type only occurs once
        # (enforced when defined)

        #### TESTING
        ####par.OUTPUT_MODE = ['MERGED','UNMERGED']

        print("**** outputOptions OUTPUT_MODE", par.OUTPUT_MODE)
        unmergedout = False
        commandlinehklout = False
        for mode in par.OUTPUT_MODE:
            if mode == 'MERGED':
                self.appendCommandScript("OUTPUT MTZ MERGED")
                self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT_BASENAME + '.mtz'])
                commandlinehklout = True
            elif mode == 'UNMERGED':
                self.appendCommandScript("OUTPUT MTZ UNMERGED")
                if not commandlinehklout:  # not needed if already done
                    self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT_BASENAME + '.mtz'])
                unmergedout = True
            elif mode == 'ORIGINAL':
                if unmergedout:  # only relevant for unmerged output
                    self.appendCommandScript("OUTPUT ORIGINAL")
            elif mode == 'SP_MERGED':
                self.appendCommandScript("OUTPUT SCALEPACK MERGED")
                self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT_BASENAME + '.sca'])
            elif mode == 'SP_UNMERGED':
                self.appendCommandScript("OUTPUT SCALEPACK UNMERGED")
                self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT_BASENAME + '.sca'])
            elif mode == 'NONE':
                self.appendCommandScript("OUTPUT NONE")
        
        return 0
    
    # -----------------------------------------------------------------------
    def processOutputFiles(self):
        print("AIMLESS FINISHED START")
        
        import glob
        nOutFiles = 0
        par = self.container.controlParameters
        out = self.container.outputData

        # out.MTZMERGEDOUT is list of merged MTZ files for processing by ctruncate
        outfilesOK = False   #  will fail if no files unless mode = NONE
        for mode in par.OUTPUT_MODE:
            print("processOutputFiles mode:", mode)
            if mode == 'MERGED':
                # Merged files are HKLOUT[_dname], but not "_unmerged"
                for file in glob.glob(os.path.join(self.getWorkDirectory(),str(out.HKLOUT_BASENAME)+"*.mtz")):
                    if not 'unmerged' in file:
                        nOutFiles += 1
                        outfilesOK = True
                        out.MTZMERGEDOUT.append(file)
                        print("Adding to MTZMERGEDOUT:", file)
                        #print("Type:", type(file))
                        #print("HKLOUT_BASENAME:", str(out.HKLOUT_BASENAME))
                        #print("Work directory", self.getWorkDirectory())
            elif mode == 'UNMERGED':
                for file in glob.glob(os.path.join(self.getWorkDirectory(),str(out.HKLOUT_BASENAME)+"*.mtz")):
                    if 'unmerged' in file:
                        nOutFiles += 1
                        out.MTZUNMERGEDOUT.append(file)
            elif mode == 'SP_MERGED':
                if not file.find('unmerged'):
                    for file in glob.glob(os.path.join(self.getWorkDirectory(),str(out.HKLOUT_BASENAME)+"*.sca")):
                        nOutFiles += 1
                        out.SPMERGEDOUT.append(file)
            elif mode == 'SP_UNMERGED':
                for file in glob.glob(os.path.join(self.getWorkDirectory(),str(out.HKLOUT_BASENAME)+"*.sca")):
                    if 'unmerged' in file:
                        nOutFiles += 1
                        out.SPUNMERGEDOUT.append(file)
            elif mode == 'NONE':
                outfilesOK = True

        print("AIMLESS FINISHED END")
        if outfilesOK:
          return CPluginScript.SUCCEEDED
        else:
          self.appendErrorReport(201)
          return CPluginScript.FAILED
