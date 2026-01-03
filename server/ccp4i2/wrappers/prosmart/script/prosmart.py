import glob
import os
import shutil

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript


class prosmart(CPluginScript):
    
    TASKMODULE = 'wrappers' # Where this plugin will appear on gui
    TASKNAME = 'prosmart'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin
    TASKCOMMAND = 'prosmart'   # The command to run the executable
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'

    ERROR_CODES = { 201 : { 'description' : 'No output restraint file from Prosmart' },
                    202 : { 'description' : 'Log file does not report successful job completion' , 'severity' : SEVERITY_WARNING },
                    203 : { 'description' : 'Unable to successfully construct reference model list' },
                    204 : { 'description' : 'Unable to successfully construct target chain list' }
                    }
   

    def processInputFiles(self):
        # Use temp input filename from which prosmart takes output restraints filename
        self.tempFile = os.path.splitext(str(self.container.outputData.RESTRAINTS))[0]+'_TARGET.pdb'
        print('prosmart tempFile',self.tempFile)
        shutil.copyfile(self.container.inputData.TARGET_MODEL.__str__(), self.tempFile)
        return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):

        self.appendCommandLine(['-o',self.workDirectory])
        self.appendCommandLine(['-p1',self.tempFile])
        #self.appendCommandLine(['-p2',self.container.inputData.REFERENCE_MODELS])
        
        p2 = ['-p2']
        for iCoordSet, xyzRef in enumerate(self.container.inputData.REFERENCE_MODELS):
           p2.append(str(xyzRef.fullPath))
        if len(p2) > 1:
           self.appendCommandLine(p2)
        else:
           self.appendErrorReport(203)
           return CPluginScript.FAILED
        
        # always output xml log
        self.appendCommandLine(['-xml',self.makeFileName('PROGRAMXML')])
        
        if self.container.inputData.EXT_FILE.isSet():
            self.appendCommandLine(['-f',self.container.inputData.EXT_FILE])
    
        if self.container.inputData.CHAINLIST_1.isSet():
           c1 = ['-c1']
           for chain in self.container.inputData.CHAINLIST_1.__str__().split():
              c1.append(str(chain))
           if len(c1) > 1:
              self.appendCommandLine(c1)
           else:
              self.appendErrorReport(204)
              return CPluginScript.FAILED

        if self.container.controlParameters.NUCLEIC_ACID:
            self.appendCommandLine(['-dna_rna'])

        ### Program Options
        if self.container.controlParameters.PROGRAM_MODE == 'ALIGN_RESTRAIN':
            self.appendCommandLine(['-a','-r'])
        elif self.container.controlParameters.PROGRAM_MODE == 'ALIGN':
            self.appendCommandLine(['-a'])
        elif self.container.controlParameters.PROGRAM_MODE == 'RESTRAIN':
            self.appendCommandLine(['-r'])
        if not self.container.controlParameters.PROGRAM_MODE == 'RESTRAIN':
            if self.container.controlParameters.ALIGN_MODE == '1':
                self.appendCommandLine(['-a1'])
            elif self.container.controlParameters.ALIGN_MODE == '2':
                self.appendCommandLine(['-a2'])


        ### Fragment Options
        if self.container.controlParameters.LIB_MODE == 'LIB':
            self.appendCommandLine(['-lib'])
        elif self.container.controlParameters.LIB_MODE == 'HELIX':
            self.appendCommandLine(['-helix'])
        elif self.container.controlParameters.LIB_MODE == 'STRAND':
            self.appendCommandLine(['-strand'])
               
        ### Alignment Options
        if self.container.controlParameters.FRAGLEN.isSet():
            self.appendCommandLine(['-len',self.container.controlParameters.FRAGLEN])
        if self.container.controlParameters.ALIGN_THRESHOLD.isSet():
            self.appendCommandLine(['-score',self.container.controlParameters.ALIGN_THRESHOLD])
        # (other options - far less important - may be unnecessary for this gui)
        if self.container.controlParameters.HELIX_CUTOFF.isSet():
            self.appendCommandLine(['-helix_cutoff',self.container.controlParameters.HELIX_CUTOFF])
        if self.container.controlParameters.HELIX_PENALTY.isSet():
            self.appendCommandLine(['-helix_penalty',self.container.controlParameters.HELIX_PENALTY])
        if not self.container.controlParameters.REWARD_SEQ:
            self.appendCommandLine(['-no_reward_seq'])
        if not self.container.controlParameters.ALIGN_REFINE:
            self.appendCommandLine(['-skip_refine'])
    
        ### Superposition Options
        if self.container.controlParameters.SUPERPOSE_THRESHOLD.isSet():
            self.appendCommandLine(['-superpose',self.container.controlParameters.SUPERPOSE_THRESHOLD])

        ### Scoring Options
        if self.container.controlParameters.INCLUDE_MAIN:
            self.appendCommandLine(['-main_dist'])
        if not self.container.controlParameters.PERFORM_FLIPS:
            self.appendCommandLine(['-no_fix_errors'])

        ### Rigid Substructure Identification Options
        if not self.container.controlParameters.CLUSTER_PERFORM:
            self.appendCommandLine(['-cluster_skip'])
        if self.container.controlParameters.CLUSTER_SCORE.isSet():
            self.appendCommandLine(['-cluster_score',self.container.controlParameters.CLUSTER_SCORE])
        if self.container.controlParameters.CLUSTER_ANGLE.isSet():
            self.appendCommandLine(['-cluster_angle',self.container.controlParameters.CLUSTER_ANGLE])
        if self.container.controlParameters.CLUSTER_MIN.isSet():
            self.appendCommandLine(['-cluster_min',self.container.controlParameters.CLUSTER_MIN])
        if self.container.controlParameters.CLUSTER_LINK.isSet():
            self.appendCommandLine(['-cluster_link',self.container.controlParameters.CLUSTER_LINK])
        if self.container.controlParameters.CLUSTER_RIGID.isSet():
            self.appendCommandLine(['-cluster_rigid',self.container.controlParameters.CLUSTER_RIGID])
        #if self.container.controlParameters.PYMOL_CLUSTER_COLOR.isSet():
        #    self.appendCommandLine(['-cluster_color',self.container.controlParameters.CLUSTER_COLOR])
        if self.container.controlParameters.OUTPUT_DM:
            self.appendCommandLine(['-output_dm'])
    
        ### Output Options
        if not self.container.controlParameters.DISPLAY_AS_DEGREES:
            self.appendCommandLine(['-cosine'])
        
        ### Restraint Options
        if self.container.controlParameters.RESTRAIN_SELF:
            self.appendCommandLine(['-self_restrain'])
        if self.container.controlParameters.RESTRAIN_ALL_VS_BEST.isSet():
            if self.container.controlParameters.RESTRAIN_ALL_VS_BEST == 'ALL':
                self.appendCommandLine(['-restrain_all'])
            else:
                self.appendCommandLine(['-restrain_best'])
        if self.container.controlParameters.RESTRAIN_TO_SELF:
            self.appendCommandLine(['-restrain_to_self'])
        if self.container.controlParameters.RESTRAIN_SEQID.isSet():
            self.appendCommandLine(['-restrain_seqid',self.container.controlParameters.RESTRAIN_SEQID])
        if self.container.controlParameters.RESTRAIN_RMAX.isSet():
            self.appendCommandLine(['-rmax',self.container.controlParameters.RESTRAIN_RMAX])
        if self.container.controlParameters.RESTRAIN_RMIN.isSet():
            self.appendCommandLine(['-rmin',self.container.controlParameters.RESTRAIN_RMIN])
        if self.container.controlParameters.RESTRAIN_SIGMA.isSet():
            self.appendCommandLine(['-sigma',self.container.controlParameters.RESTRAIN_SIGMA])
        if self.container.controlParameters.RESTRAIN_MIN_SIGMA.isSet():
            self.appendCommandLine(['-minsigma',self.container.controlParameters.RESTRAIN_MIN_SIGMA])
        if self.container.controlParameters.RESTRAIN_SIGMATYPE.isSet():
            if self.container.controlParameters.RESTRAIN_SIGMATYPE == 'DEFAULT':
                self.appendCommandLine(['-sigmatype','0'])
            elif self.container.controlParameters.RESTRAIN_SIGMATYPE == 'CONSTANT':
                self.appendCommandLine(['-sigmatype','1'])
            else:
                self.appendCommandLine(['-sigmatype','2'])
        if self.container.controlParameters.RESTRAIN_MAIN_CUTOFF.isSet():
            self.appendCommandLine(['-cutoff',self.container.controlParameters.RESTRAIN_MAIN_CUTOFF])
        if self.container.controlParameters.RESTRAIN_SIDE_CUTOFF.isSet():
            self.appendCommandLine(['-side_cutoff',self.container.controlParameters.RESTRAIN_SIDE_CUTOFF])
        if self.container.controlParameters.RESTRAIN_OUTLIER_THRESHOLD.isSet():
            self.appendCommandLine(['-multiplier',self.container.controlParameters.RESTRAIN_OUTLIER_THRESHOLD])
        if self.container.controlParameters.RESTRAIN_SCALE_SIGMAS.isSet():
            self.appendCommandLine(['-weight',self.container.controlParameters.RESTRAIN_SCALE_SIGMAS])

        if self.container.controlParameters.RESTRAIN_BFAC_FILTER.isSet():
           if self.container.controlParameters.RESTRAIN_BFAC_FILTER:
              if self.container.controlParameters.RESTRAIN_BFAC_ALPHA.isSet():
                 self.appendCommandLine(['-bfac_alpha',self.container.controlParameters.RESTRAIN_BFAC_ALPHA])
           else:
              self.appendCommandLine(['-no_bfac_filter'])

        if self.container.controlParameters.RESTRAIN_ALT.isSet():
           if self.container.controlParameters.RESTRAIN_ALT:
              self.appendCommandLine(['-alt'])
        if self.container.controlParameters.RESTRAIN_OCCUP.isSet():
           self.appendCommandLine(['-occup',self.container.controlParameters.RESTRAIN_OCCUP])

        if self.container.controlParameters.RESTRAIN_RM_BONDS:
            self.appendCommandLine(['-rm_bonds'])
        if self.container.controlParameters.RESTRAIN_MAIN_VS_SIDE.isSet():
            if self.container.controlParameters.RESTRAIN_MAIN_VS_SIDE == 'MAIN':
                self.appendCommandLine(['-main'])
            else:
                self.appendCommandLine(['-side'])
        if self.container.controlParameters.RESTRAIN_REFMAC_TYPE.isSet():
            if self.container.controlParameters.RESTRAIN_REFMAC_TYPE == 'TYPE0':
                self.appendCommandLine(['-type','0'])
            elif self.container.controlParameters.RESTRAIN_REFMAC_TYPE == 'TYPE1':
                self.appendCommandLine(['-type','1'])
            else:
                self.appendCommandLine(['-type','2'])
        if self.container.controlParameters.OUTPUT_PDB_CHAIN_RESTRAINTS:
            self.appendCommandLine(['-output_pdb_chain_restraints'])

        ### Generic Bond (H-Bond) Restraint Options
        if self.container.controlParameters.H_BOND:
            self.appendCommandLine(['-h'])
        if self.container.controlParameters.H_HELIX:
            self.appendCommandLine(['-h_helix'])
        if self.container.controlParameters.H_SHEET:
            self.appendCommandLine(['-h_sheet'])
        if self.container.controlParameters.H_310:
            self.appendCommandLine(['-3_10'])
        if self.container.controlParameters.H_ALPHA:
            self.appendCommandLine(['-alpha'])
        if self.container.controlParameters.H_PI:
            self.appendCommandLine(['-pi'])
        if self.container.controlParameters.H_STRICT:
            self.appendCommandLine(['-h_strict'])
        if self.container.controlParameters.H_DIST.isSet():
            self.appendCommandLine(['-bond_dist',self.container.controlParameters.H_DIST])
        if self.container.controlParameters.H_MIN.isSet():
            self.appendCommandLine(['-bond_min',self.container.controlParameters.H_MIN])
        if self.container.controlParameters.H_MAX.isSet():
            self.appendCommandLine(['-bond_max',self.container.controlParameters.H_MAX])
        if self.container.controlParameters.H_MIN_SEP.isSet():
            self.appendCommandLine(['-min_sep',self.container.controlParameters.H_MIN_SEP])
        if self.container.controlParameters.H_MAX_SEP.isSet():
            self.appendCommandLine(['-max_sep',self.container.controlParameters.H_MAX_SEP])
        if self.container.controlParameters.H_BOND_OPT.isSet():
            if self.container.controlParameters.H_BOND_OPT == 'TYPE1':
                self.appendCommandLine(['-bond_opt','1'])
            elif self.container.controlParameters.H_BOND_OPT == 'TYPE2':
                self.appendCommandLine(['-bond_opt','2'])
            else:
                self.appendCommandLine(['-bond_opt','3'])
        if self.container.controlParameters.H_BOND_OVERRIDE.isSet():
            self.appendCommandLine(['-bond_override',self.container.controlParameters.H_BOND_OVERRIDE])

        ### Advanced Options
        if self.container.controlParameters.MERGE_CHAINS:
            self.appendCommandLine(['-merge'])
        if self.container.controlParameters.IS_NMR_MD_ENSEMBLE:
            self.appendCommandLine(['-nmr'])
        if self.container.controlParameters.RENAME_CHAIN.isSet():
            self.appendCommandLine(['-rename_chain',self.container.controlParameters.RENAME_CHAIN])
               
        ### Custom Keywords
        if self.container.controlParameters.KEYWORDS.isSet():
           keys = []
           for key in str(self.container.controlParameters.KEYWORDS).split():
              keys.append(str(key))
           if len(keys) > 0:
              self.appendCommandLine(keys)

    def processOutputFiles(self):
        try:
          #Remove potentially confusing tempFile
          os.remove(self.tempFile)
        except:
          pass
        
        if not self.container.outputData.RESTRAINTS.exists():
          resFileList = glob.glob(os.path.splitext(self.container.outputData.RESTRAINTS.__str__())[0]+'*.txt')
          if len(resFileList)==0:
            self.appendErrorReport(201,str(self.container.outputData.RESTRAINTS))
            return CPluginScript.FAILED
          else:
            shutil.copyfile(resFileList[0],self.container.outputData.RESTRAINTS.__str__())

        self.container.outputData.RESTRAINTS.annotation = 'Restraints for ' + str(self.container.inputData.TARGET_MODEL.annotation)

        # sanity check that prosmart has produced something
        ok = False
        logText = self.logFileText()
        pyListLogLines = logText.split("\n")
        for j, pyStrLine in enumerate(pyListLogLines):
            if "ProSMART completed." in pyStrLine:
                ok = True
        if not ok: self.appendErrorReport(202)
        
        return CPluginScript.SUCCEEDED
