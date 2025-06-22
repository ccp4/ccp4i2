#=======================================================================================
#
#    acorn.py : acorn(CPluginScript)
#    
#    Author  : Kyle Stevenson,STFC
#    Created : 29th May 2017, KJS
#
#    Complementary Gui Class for ab initio solution of .
#    Handles the command and script processing.
#
#=======================================================================================

import os
import sys
import traceback
import xml.etree.ElementTree as ET

import clipper

from ....core import CCP4XtalData
from ....core.CCP4PluginScript import CPluginScript
from ....core.CCP4Utils import writeXml
from ....smartie import smartie


class acorn(CPluginScript):
    
    TASKMODULE = 'density_modification'                 # Gui-2 Menu Class
    TASKTITLE = 'Acorn'                         # Menu title
    TASKNAME = 'acorn'                          # Task name - should be same as class name
    TASKCOMMAND = 'acorn'                       # The command to run the executable
    TASKVERSION= 1.0                            # Version of this plugin
    COMTEMPLATE = None                          # The program com file template
    COMTEMPLATEFILE = None                      # Name of file containing com file template
    PERFORMANCECLASS = 'CExpPhasPerformance'    # KJS Need to change this
    # ASYNCHRONOUS = False
    ERROR_CODES = {  200 : { 'description' : 'Unexpected reflections in input file, list in stdout.txt'} }
    MAINTAINER = "Kyle.Stevenson@stfc.ac.uk"
    
    def __init__(self,*args, **kwargs):
        CPluginScript.__init__(self,*args, **kwargs)
    
    def processInputFiles(self):
        print("Processing INPUT Files - (Acorn)")
        
        # Need to input the necessary input files for an Acorn Run. This will be conditional.
        self.hklin, __, error = self.makeHklInput([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]],
                                                  extendOutputColnames=False, useInputColnames=False)

        if self.container.controlParameters.ACORN_ECALC :
            try:
                # Using clipper-python to extend data (optionally correct for anisotropy) and calculate Es
                # set options
                if self.container.controlParameters.ACORN_PHSIN_TYPE == "phases":
                    self.hklin_phifom, __, error = self.makeHklInput([['ABCD',CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM]],
                                                                     hklin='phases', extendOutputColnames=False, useInputColnames=False)
                else:
                    self.hklin_phifom = None
                if self.container.controlParameters.ACORN_ANISOTROPY :
                    aniso = True
                else:
                    aniso = False
                # initial parameters
                resfilt = -1.0
                nparm = 12 # as in cecalc
                # open mini-MTZ with f/sigf  
                mtzin_fobs = clipper.CCP4MTZfile()
                mtzin_fobs.open_read(self.hklin)
                # spacegroup and cell from first file
                sg = mtzin_fobs.spacegroup()
                cell = mtzin_fobs.cell()
                # extend to either requested resolution or default (1.0A)
                if self.container.controlParameters.ACORN_EXTEND :
                    try:
                        extres = float(self.container.controlParameters.ACORN_EXTENDRES)
                    except:
                        extres = 1.0
                    reso = clipper.Resolution(extres)
                else:
                    reso = mtzin_fobs.resolution()
            
                hkls = clipper.HKL_info(sg, cell, reso, True)
                # import F/SIGF
                fsig = clipper.HKL_data_F_sigF_float(hkls)
                mtzin_fobs.import_hkl_data(fsig, '*/*/[F, SIGF]')
                # and again if aniso correction required
                if aniso:
                    fsigiso = clipper.HKL_data_F_sigF_float(hkls)
                    mtzin_fobs.import_hkl_data(fsigiso, '*/*/[F, SIGF]')
                # close file to actually read data
                mtzin_fobs.close_read()
                # optionally read in phases
                if self.hklin_phifom is not None:
                    mtzin_phifom = clipper.CCP4MTZfile()
                    mtzin_phifom.open_read(self.hklin_phifom)
                    phifom = clipper.HKL_data_Phi_fom_float(hkls)
                    mtzin_phifom.import_hkl_data(phifom, '*/*/[PHI, FOM]')
                    mtzin_phifom.close_read()
                # optionally correct anisotropy
                if aniso:
                    F = clipper.SFscale_aniso_float.SFscale_aniso_F_float
                    M = clipper.SFscale_aniso_float.SFscale_aniso_NORMAL_float
                    sfscale = clipper.SFscale_aniso_float(3.0, M) # rejection criterion for F/sigF from caniso
                    sfscale(fsigiso, resfilt, 12)
                    print("\n Correcting for anisotropy. Scale factors: \n") 
                    print(str(sfscale.u_aniso_orth(F)))
                # Fs to Es
                esig = clipper.HKL_data_E_sigE_float(hkls)
                if aniso: 
                    esig.compute_from_fsigf(fsigiso)
                else:
                    esig.compute_from_fsigf(fsig)
           
                # now calculate scaling
                initial_params = clipper.DoubleVector(nparm, 1.0)
                basis_f = clipper.BasisFn_spline(esig, nparm, 2.0)
                target_f = clipper.TargetFn_scaleEsq_E_sigE(esig)
                escale = clipper.ResolutionFn(hkls, basis_f, target_f, initial_params)
                # apply scaling
                esig.scaleBySqrtResolution(escale)
                # write file
                mtzout = clipper.CCP4MTZfile()
                output_file = str(os.path.join(self.getWorkDirectory(), 'EXTENDED.mtz'))
                mtzout.open_write(output_file)
                mtzout.export_hkl_info(fsig.hkl_info())
                mtzout.export_hkl_data(fsig, "*/*/[F, SIGF]")
                if aniso:
                    mtzout.export_hkl_data(fsigiso, "*/*/[F_ISO, SIGF_ISO]")
                    mtzout.export_hkl_data(esig, "*/*/[E_ISO, SIGE_ISO]")
                else:
                    mtzout.export_hkl_data (esig, "*/*/[E, SIGE]")
                if self.hklin_phifom is not None:
                    mtzout.export_hkl_data (phifom, "*/*/[PHI, FOM]")
                mtzout.close_write()
                return CPluginScript.SUCCEEDED
            except:
                print("\n \n ***** Error!! please submit bug report with traceback below ***** \n")
                traceback.print_tb(sys.exc_info()[2])
                print("\n \n")
                return CPluginScript.FAILED

        print("Finishing processINPUT (KJS)")

    def processOutputFiles(self):
        print("Process OUPUT (KJS)")
        #if (self.container.controlParameters.ACOMPS_PEAKSEARCH): #  forgot Eleanor told me this isn't used these days (.. leave out)
        #    print "You should not be seeing this (KJS - Acorn)" # leave commented out section in for now.
        #    pido1 = PROCESSMANAGER().startProcess( binfft, arglisto1, logFile=logfile1 )
        #    stat1 = PROCESSMANAGER().getJobData( pid1 )
        #    ex1    = PROCESSMANAGER().getJobData( pid1,'exitCode' )
        
        #    pido2 = PROCESSMANAGER().startProcess( binmapmask, arglisto2, logFile=logfile2 )
        #   stat2 = PROCESSMANAGER().getJobData( pid2 )
        #   ex2   = PROCESSMANAGER().getJobData( pid2,'exitCode' )
        
        #   inTxt3  = "threshold rms %s"%(str(self.container.inputData.ACOMPS_MAXPEAKS)) 
        #   inTxt3 += "numpeaks %s"%(str(self.container.inputData.ACOMP_RMSMULT))
        #   pido3 = PROCESSMANAGER().startProcess( binpeakmax, arglisto3, inputText = inTxt3, logFile=logfile3 )
        #   stat3 = PROCESSMANAGER().getJobData( pid3 )
        #   ex3   = PROCESSMANAGER().getJobData( pid3,'exitCode' )
        
        # Parse the output text file to create an xml file which can then be parsed by the report to make html .......
        rootNode = ET.Element("acorn")
        xmlRI = ET.SubElement(rootNode,"RunInfo")
        
        # Use the ccp4 Smartie Class to parse the ascii log file from Acorn.
        aclfile = self.makeFileName('LOG')
        smfile  = smartie.parselog(aclfile)
        
        tabs = smfile.tables("Correlation Coefficient vs number of cycles")
        vcycnum = tabs[0].col("Cycle_number")
        vccoef = tabs[0].col("CC")
        
        for cn, cc in zip(vcycnum,vccoef):
            xmlcyc = ET.SubElement(xmlRI,"Cycle")
            ET.SubElement(xmlcyc,"NCycle").text          = str(cn)
            ET.SubElement(xmlcyc,"CorrelationCoef").text = str(cc)
        
        writeXml(rootNode, self.xmlout)
        
        #Informative labels
        self.container.outputData.PHSOUT.annotation = 'Phase probabilities for measured reflections only'
        self.container.outputData.FPHIOUT.annotation = 'Map coefficients for measured reflections only'
        self.container.outputData.EPHIOUT.annotation = 'Extended normalised map coefficients'
        
        # generate ouput mini-mtzs using clipper python after removing phases for extended reflections.
        try:
            mtzref = clipper.CCP4MTZfile()
            mtzphifom = clipper.CCP4MTZfile()
            mtzfphi = clipper.CCP4MTZfile()
            mtzephi = clipper.CCP4MTZfile()
            hkls_ref = clipper.HKL_info()
            # open the reference and read in f_sigf
            mtzref.open_read(self.hklin)
            sgref = mtzref.spacegroup()
            cellref = mtzref.cell()
            ref_reso = mtzref.resolution()
            mtzref.import_hkl_info(hkls_ref, False)
            fsig_ref = clipper.HKL_data_F_sigF_float(hkls_ref)
            mtzref.import_hkl_data(fsig_ref, "*/*/[F,SIGF]")
            mtzref.close_read()

            mtzin = clipper.CCP4MTZfile()
            #  read fsig, f_iso, eo_ext and phi/fom from acorn
            acorn_hklout = str(os.path.join(self.getWorkDirectory(), 'hklout.mtz'))
            mtzin.open_read(acorn_hklout)
            # paranoia
            assert mtzin.spacegroup().symbol_hall() == sgref.symbol_hall()
            cell = mtzin.cell()
            # more paranoia
            assert cell.a() == cellref.a() and cell.b() == cellref.b() and cell.c() == cellref.c() and cell.alpha() == cellref.alpha() and cell.beta() == cellref.beta() and cell.gamma() == cellref.gamma()
            hkls = clipper.HKL_info()
            mtzin.import_hkl_info(hkls)
            reso = mtzin.resolution()
            fsig = clipper.HKL_data_F_sigF_float(hkls)
            fsigiso = clipper.HKL_data_F_sigF_float(hkls)
            eoext = clipper.HKL_data_F_sigF_float(hkls)
            phifom = clipper.HKL_data_Phi_fom_float(hkls)
            mtzin.import_hkl_data(fsig, "*/*/[F,SIGF]")
            if self.container.controlParameters.ACORN_ANISOTROPY :
                mtzin.import_hkl_data(fsigiso, "*/*/[F_ISO,SIGF_ISO]")
            else:
                # just use the uncorrected Fs
                mtzin.import_hkl_data(fsigiso, "*/*/[F,SIGF]")
            mtzin.import_hkl_data(eoext, "*/*/[EOEXT,EOEXT]")
            mtzin.import_hkl_data(phifom, "*/*/[PHIOUT,WTOUT]")
            
            if reso.limit() < ref_reso.limit(): 
                # extended data so read in only those reflections in original file
                fsigiso_out = clipper.HKL_data_F_sigF_float(hkls_ref)
                phifom_out = clipper.HKL_data_Phi_fom_float(hkls_ref)
                fphi_out = clipper.HKL_data_F_phi_float(hkls_ref)
                if self.container.controlParameters.ACORN_ANISOTROPY :
                    mtzin.import_hkl_data(fsigiso_out, "*/*/[F_ISO,SIGF_ISO]")
                else:
                    mtzin.import_hkl_data(fsigiso_out, "*/*/[F,SIGF]")
                mtzin.import_hkl_data(phifom_out, "*/*/[PHIOUT,WTOUT]")
            else:
                # user truncated data 
                fsigiso_out = fsigiso
                phifom_out = phifom
                fphi_out = clipper.HKL_data_F_phi_float(hkls)
            mtzin.close_read()
            
            
            if reso.limit() < ref_reso.limit(): 
                # This is so the output Mini-MTZs files have MNFs for reflections missing in input file 
                dummy_flag = clipper.HKL_data_Flag(hkls_ref)
                missing = clipper.HKL_data_Flag(hkls_ref)
                # this sets missing to 0 for missing reflections and 3 for rest. 
                clipper.SetFlagBothIfMissing(missing,fsig_ref,dummy_flag,1)
                # this sets the missing data to NaN. 
                fsigiso_out.mask(missing > 0)
                phifom_out.mask(missing > 0)
            
            # calculate map coefficients
            fphi_out.compute_from_fsigf_phifom(fsigiso_out, phifom_out)
            ephi_out = clipper.HKL_data_F_phi_float(hkls)
            ephi_out.compute_from_fsigf_phifom(eoext, phifom)

            # write out mini-MTZs
            acorn_phases = str(self.container.outputData.PHSOUT.fullPath)
            map_coeffs = str(self.container.outputData.FPHIOUT.fullPath)
            emap_coeffs =  str(self.container.outputData.EPHIOUT.fullPath)

            mtzphifom.open_write(acorn_phases)
            mtzphifom.export_hkl_info(phifom_out.hkl_info())
            mtzphifom.export_hkl_data(phifom_out, "*/*/[PHI,FOM]")
            mtzphifom.close_write()
            mtzfphi.open_write(map_coeffs)
            mtzfphi.export_hkl_info(fphi_out.hkl_info())
            mtzfphi.export_hkl_data(fphi_out, "*/*/[F,PHI]")
            mtzfphi.close_write()
            mtzephi.open_write(emap_coeffs)
            mtzephi.export_hkl_info(ephi_out.hkl_info())
            mtzephi.export_hkl_data(ephi_out, "*/*/[F,PHI]")
            mtzephi.close_write()    

            return CPluginScript.SUCCEEDED        
        
        except:
            print("\n \n ***** Error!! please submit bug report with traceback below ***** \n")
            traceback.print_tb(sys.exc_info()[2])
            print("\n \n")
            return CPluginScript.FAILED

    def makeCommandAndScript(self, container=None):
        print("Constructing Command Script for Acorn (KJS)")

        mainfile = os.path.join(self.getWorkDirectory(), 'EXTENDED.mtz')
        moutfile = os.path.join(self.getWorkDirectory(), 'hklout.mtz')

        #coordfile = self.container.inputData.XYZIN.fullPath.__str__()
        if self.container.inputData.XYZIN.isSelectionSet():
          self.selectedXYZIN = os.path.join(self.workDirectory,'XYZIN_selected_atoms.pdb')
          rv = self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedXYZIN)
        else:
          self.selectedXYZIN = self.container.inputData.XYZIN.fullPath.__str__()
        
        self.appendCommandLine("hklin"); self.appendCommandLine(mainfile) 
        self.appendCommandLine("hklout"); self.appendCommandLine(moutfile)
        if self.container.controlParameters.ACORN_PHSIN_TYPE == "model":
            #self.appendCommandLine("xyzin"); self.appendCommandLine(coordfile)
            self.appendCommandLine("xyzin"); self.appendCommandLine(self.selectedXYZIN)
            # always use all fragment
            self.appendCommandScript("POSI 1")
        
        if self.container.controlParameters.ACORN_ANISOTROPY :
            runacols = "labin  FP=F_ISO SIGFP=SIGF_ISO E=E_ISO"
        else:
            runacols = "labin  FP=F SIGFP=SIGF E=E"
        if self.container.controlParameters.ACORN_PHSIN_TYPE == "phases":
            runacols+= " PHIN=PHI WTIN=FOM"

        self.appendCommandScript(runacols)
        self.appendCommandScript("labout  PHIOUT=PHIOUT WTOUT=WTOUT")
        
        # The general stuff
        if self.container.controlParameters.ACORN_EXTEND:
            self.appendCommandScript("EXTEND") #  nb. toxd example just doesn't work (non matching cell volumes) & never did.
        else:
            self.appendCommandScript("NOEXTEND")
            
        self.appendCommandScript("NTRY %s"%(str(self.container.controlParameters.ACOPH_TRIALS)))
        self.appendCommandScript("PSFINISH %s"%(str(self.container.controlParameters.ACOPH_PSFINISH)))
                  
        # Reflection data
        if self.container.controlParameters.ACORN_BRESOL:
            self.appendCommandScript("RESOLUTION %s %s"%(str(self.container.controlParameters.ACOREF_RESOLL),str(self.container.controlParameters.ACOREF_RESOLU)))
        if self.container.controlParameters.ACORN_BEXCLUDE:
            self.appendCommandScript("EXCLUDE %s"%(str(self.container.controlParameters.ACOREF_EXCLUDE)))
        if self.container.controlParameters.ACORN_BECUT:
            self.appendCommandScript("ECUT %s"%(str(self.container.controlParameters.ACOREF_ECUT)))
        
        # Advanced (Developer?) parameters (Simple ver. for first pass)
        if self.container.controlParameters.ACOPH_PATSUP:
            self.appendCommandScript("SUPP 1") 
        self.appendCommandScript("CUTDDM %s"%(str(self.container.controlParameters.ACOPH_CUTDDM)))
        # if self.container.controlParameters.ACORN_BSEED:
        #     self.appendCommandScript( "SEED %s"%(str(self.container.controlParameters.ACOGEN_SEED)) )
        if self.container.controlParameters.ACORN_BGRID:
            self.appendCommandScript("GRID %s"%(str(self.container.controlParameters.ACOGEN_GRID)))
        if self.container.controlParameters.ACOPH_CUSTOM:
            ncser = "NCSER "
            enhs = "ENHS "
            ddmk = "DDMK "
            ncddm = "NCDDM "
            for i in range(self.container.controlParameters.ACOPH_TRIALS.__int__()):
                ACOPH_REFINE_N = getattr(self.container.controlParameters, 'ACOPH_REFINE_%s' % str(i+1))
                ACOPH_DDMK_N = getattr(self.container.controlParameters, 'ACOPH_DDMK_%s' % str(i+1))
                ACOPH_NCDDM_N = getattr(self.container.controlParameters, 'ACOPH_NCDDM_%s' % str(i+1))
                if ACOPH_REFINE_N == 'NO':
                    ncser += " 0"
                    enhs += " 0"
                elif ACOPH_REFINE_N == 'NCSER':
                    ncser += " 2"
                    enhs += " 0"
                elif ACOPH_REFINE_N == 'ENHS':
                    ncser += " 0"
                    enhs += " 1"
                if self.container.controlParameters.ACOPH_CUSTDDM:
                    if ACOPH_DDMK_N == 'DDM1':
                        ddmk += " 1"
                    elif ACOPH_DDMK_N == 'DDM2':
                        ddmk += " 2"
                    else:
                        ddmk += " 0"
                ncddm += str(ACOPH_NCDDM_N)+" "
            self.appendCommandScript(ncser)
            self.appendCommandScript(enhs)
            if self.container.controlParameters.ACOPH_CUSTDDM:
                self.appendCommandScript(ddmk)
            self.appendCommandScript(ncddm)
        
        # Map Peaks Search ....# post-proc... removed (Eleanor - obsolete)
        # self.appendCommandScript( " %s"%(str(self.container.inputData.ACOPEAK_XXXX)) )
        # self.appendCommandScript( " %s"%(str(self.container.inputData.ACOPEAK_XXXX)) )
        
        self.xmlout = self.makeFileName('PROGRAMXML')

        return CPluginScript.SUCCEEDED
    
