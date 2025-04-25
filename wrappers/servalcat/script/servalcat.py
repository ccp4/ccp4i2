"""
    servalcat.py: CCP4 GUI Project
    Copyright (C) 2024 MRC-LBM, University of Southampton
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,S
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core.CCP4ErrorHandling import *
from core import CCP4Modules, CCP4XtalData, CCP4Utils
from xml.etree import ElementTree as ET
from .json2xml import json2xml
import pathlib
import os
import json
import shutil


class servalcat(CPluginScript):
    
    TASKMODULE = 'wrappers'
    TASKTITLE = 'Refinement (servalcat)'
    TASKNAME = 'servalcat'
    TASKCOMMAND = 'servalcat'
    TASKVERSION= 0.0
    ASYNCHRONOUS = False
    PERFORMANCECLASS = 'CServalcatPerformance'
        
    ERROR_CODES = { 201 : {'description' : 'Refmac returned with non zero status' },
                    202:  {'description': 'New library created but strictly required' },
                    203:  {'description': 'New library created', 'severity':CCP4ErrorHandling.SEVERITY_WARNING},
                    204:  {'description': 'Program completed without generating XMLOUT.' },
                    }
    
    def __init__(self,*args, **kwargs):
        super(servalcat, self).__init__(*args, **kwargs)
        self._readyReadStandardOutputHandler = self.handleReadyReadStandardOutput
        self.xmlroot = ET.Element('SERVALCAT')
        self.xmlLength = 0

    @QtCore.Slot()
    def handleReadyReadStandardOutput(self):
        if not hasattr(self,'logFileHandle'): self.logFileHandle = open(self.makeFileName('LOG'),'w')
        if not hasattr(self,'logFileBuffer'): self.logFileBuffer = ''
        pid = self.getProcessId()
        qprocess = CCP4Modules.PROCESSMANAGER().getJobData(pid,attribute='qprocess')
        availableStdout = qprocess.readAllStandardOutput()
        if sys.version_info > (3,0):
            self.logFileHandle.write(availableStdout.data().decode("utf-8"))
        else:
            self.logFileHandle.write(availableStdout)
        self.logFileHandle.flush()

    def xmlAddRoot(self, xmlText, xmlFilePath=None, xmlRootName=None):
        if xmlRootName:
            xmlText = f"<{xmlRootName}>\n{xmlText}\n</{xmlRootName}>"
        if xmlFilePath:
            with open(xmlFilePath, 'w') as programXmlFile:
                programXmlFile.write(xmlText)
        return xmlText

    def processInputFiles(self):
        error = None
        self.hklin = None
        dataObjects = []
        
        #Apply coordinate selection if set
        self.inputCoordPath = os.path.normpath(self.container.inputData.XYZIN.fullPath.__str__())
        if self.container.inputData.XYZIN.isSelectionSet():
            self.inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected.pdb'))
            self.container.inputData.XYZIN.loadFile()
            if self.container.inputData.XYZIN.isMMCIF():
                self.inputCoordPath = str(pathlib.Path(self.inputCoordPath).with_suffix('.cif'))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.inputCoordPath)

        #Create DICT by merging dictionaries in DICT_LIST
        rv = self.joinDicts(self.container.outputData.DICT, self.container.inputData.DICT_LIST)

        #Append Observation with representation dependent on whether we are detwining on Is or not

        if str(self.container.controlParameters.DATA_METHOD) == 'xtal':
            obsTypeRoot = 'CONTENT_FLAG_F'
            obsPairOrMean = 'MEAN'
            if self.container.inputData.HKLIN.isSet():
                if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
                    obsTypeRoot = 'CONTENT_FLAG_I'
                    obsPairOrMean = 'PAIR'
                elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
                    obsTypeRoot = 'CONTENT_FLAG_I'
                    obsPairOrMean = 'MEAN'
                elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR:
                    obsTypeRoot = 'CONTENT_FLAG_F'
                    obsPairOrMean = 'PAIR'
                elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN:
                    obsTypeRoot = 'CONTENT_FLAG_F'
                    obsPairOrMean = 'MEAN'
            if self.container.controlParameters.F_SIGF_OR_I_SIGI.isSet():
                # overwrite to use F despite available I
                if str(self.container.controlParameters.F_SIGF_OR_I_SIGI) == "F_SIGF":
                    obsTypeRoot = 'CONTENT_FLAG_F'
            if str(self.container.controlParameters.SCATTERING_FACTORS) in ["electron", "neutron"]:
                # overwrite to not use anomalous pairs even if they are present in case of data from electron or neutron diffraction
                obsPairOrMean = 'MEAN'
            obsType = getattr(CCP4XtalData.CObsDataFile, obsTypeRoot+obsPairOrMean)
            dataObjects += [['HKLIN', obsType]]

            #Include FreeRflag if called for
            if self.container.inputData.FREERFLAG.isSet():
                dataObjects += ['FREERFLAG']
            self.hklin,error = self.makeHklin(dataObjects)
            if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED

            return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        if hasattr(self,'logFileHandle'):
            self.logFileHandle.write("JOB TITLE SECTION\n")
            jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId)
            if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                self.logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
            while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
                jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobInfo["parentjobid"])
                if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                    self.logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
            self.logFileHandle.close()
        else:
            self.xmlroot.clear()
        
        # First up check for exit status of the program
        from core.CCP4Modules import PROCESSMANAGER
        exitStatus = 0
        exitCode = 0
        try:
            exitStatus = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitStatus')
        except Exception as e:
            print(e)
            self.appendErrorReport(201, 'Exit status: Unable to recover exitStatus')
            return CPluginScript.FAILED
        if exitStatus != 0:
            self.appendErrorReport(201, 'Exit status: ' + str(exitStatus))
            return CPluginScript.FAILED
        
        # Now the exit codes...I think that non zero means Refmac identified an issue
        try:
            exitCode = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitCode')
        except:
            self.appendErrorReport(201, 'Exit code: Unable to recover exitCode')
            return CPluginScript.FAILED
        if exitCode != 0:
            import os
            try:
                logFileText = open(self.makeFileName('LOG')).read()
            except:
                self.appendErrorReport(201, 'Exit code: '+ str(exitCode))
            return CPluginScript.FAILED

        import os
        from core import CCP4XtalData

        outputCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.mmcif'))
        self.container.outputData.CIFFILE.setFullPath(outputCifPath)
        self.container.outputData.CIFFILE.annotation.set('Model from refinement (mmCIF format)')
        outputPdbPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.pdb'))
        if os.path.isfile(outputPdbPath):
            self.container.outputData.XYZOUT.setFullPath(outputPdbPath)
            self.container.outputData.XYZOUT.annotation.set('Model from refinement (PDB format)')
        self.container.outputData.FPHIOUT.annotation.set('Density map (Fourier coeff.)')
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.annotation.set('Difference density map (Fourier coeff.)')
        self.container.outputData.FPHIOUT.subType = 2
        outputFiles = ['FPHIOUT', 'DIFFPHIOUT']
        outputColumns = ['FWT,PHWT', 'DELFWT,PHDELWT']

        if str(self.container.controlParameters.DATA_METHOD) == "xtal":
            hkloutFilePath = str(os.path.join(self.getWorkDirectory(), "refined.mtz"))
        else:  # spa
            hkloutFilePath = str(os.path.join(self.getWorkDirectory(), "refined_diffmap.mtz"))
        hkloutFile=CCP4XtalData.CMtzDataFile(hkloutFilePath)
        hkloutFile.loadFile()
        columnLabelsInFile = [column.columnLabel.__str__() for column in hkloutFile.fileContent.listOfColumns]
        print('columnLabelsInFile', columnLabelsInFile)
        if str(self.container.controlParameters.DATA_METHOD) == "xtal":
            if 'FAN' in columnLabelsInFile and 'PHAN' in columnLabelsInFile:
                self.container.outputData.ANOMFPHIOUT.annotation.set('Anomalous difference map')
                outputFiles += ['ANOMFPHIOUT']
                outputColumns += ['FAN,PHAN']
        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know contentFlags
        error = self.splitHklout(outputFiles, outputColumns, hkloutFilePath)
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        if str(self.container.controlParameters.DATA_METHOD) == 'spa':
            self.container.outputData.MAP_FO.annotation.set('Density map (in real space)')
            outputMapFoPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined_diffmap_normalized_fo.mrc'))
            self.container.outputData.MAP_FO.setFullPath(outputMapFoPath)
            self.container.outputData.MAP_FOFC.annotation.set('Difference density map (in real space)')
            outputMapFoFcPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined_diffmap_normalized_fofc.mrc'))
            self.container.outputData.MAP_FOFC.setFullPath(outputMapFoFcPath)
            # Write a Coot script with set_contour_level_absolute()
            cootScriptI2FilePath = os.path.join(self.getWorkDirectory(), "refined_coot_i2.py")
            with open(cootScriptI2FilePath, "w") as cootScriptI2File:
                cootScriptI2File.write("set_contour_level_absolute(1, 1.2)\n")
                cootScriptI2File.write("set_contour_level_absolute(2, 3.0)\n")
            self.container.outputData.COOTSCRIPTOUT.annotation.set('Coot script')
            self.container.outputData.COOTSCRIPTOUT.setFullPath(cootScriptI2FilePath)

        # Use output JSON from servalcat as the basis for output XML
        # Get stats from JSON, convert to XML, save and load to self.xmlroot
        rxml = None
        jsonFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.json"))
        if os.path.isfile(jsonFilePath):
            with open(jsonFilePath, "r") as jsonFile:
                jsonText = jsonFile.read()
        else:
            rxml = ET.Element('SERVALCAT')
            self.appendErrorReport(204, self.makeFileName('PROGRAMXML'))
            return CPluginScript.FAILED
        try:
            jsonStats = list(json.loads(jsonText))
            # add d_max_4ssqll and d_min_4ssqll
            for i in range(len(jsonStats)):
                for j in range(len(jsonStats[i]["data"]["binned"])):
                    jsonStats[i]["data"]["binned"][j]['d_min_4ssqll'] = \
                        1 / (list(jsonStats)[i]["data"]["binned"][j]['d_min'] * list(jsonStats)[i]["data"]["binned"][j]['d_min'])
                    jsonStats[i]["data"]["binned"][j]['d_max_4ssqll'] = \
                        1 / (list(jsonStats)[i]["data"]["binned"][j]['d_max'] * list(jsonStats)[i]["data"]["binned"][j]['d_max'])
            xmlText = json2xml(jsonStats, tag_name_subroot="cycle")
            xmlFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.xml"))
            xmlText = self.xmlAddRoot(xmlText, xmlFilePath, xmlRootName="SERVALCAT")
            # self.xmlroot = CCP4Utils.openFileToEtree(xmlFilePath)
            self.xmlroot = ET.fromstring(xmlText)
            rxml = self.xmlroot
            et = ET.ElementTree(self.xmlroot)
            # Write program.xml_tmp and move it to program.xml
            self.flushXml()
        except:
            rxml = ET.Element('SERVALCAT')
            self.appendErrorReport(204, self.makeFileName('PROGRAMXML'))
            return CPluginScript.FAILED

        # Extract performance indicators from JSON
        statsOverall = {}
        try:
            statsOverall["-LL"] = jsonStats[-1]['data']['summary']['-LL']
        except:
            pass
        if str(self.container.controlParameters.DATA_METHOD) == 'spa':
            statsOverall["FSCaverage"] = jsonStats[-1]['data']['summary']['FSCaverage']
            self.container.outputData.PERFORMANCEINDICATOR.FSCaverage.set(statsOverall["FSCaverage"])
        elif self.container.inputData.FREERFLAG.isSet():
            try:
                statsOverall["Rwork"] = jsonStats[-1]['data']['summary']['Rwork']
                statsOverall["Rfree"] = jsonStats[-1]['data']['summary']['Rfree']
                statsOverall["CCFworkavg"] = jsonStats[-1]['data']['summary']['CCFworkavg']
                statsOverall["CCFfreeavg"] = jsonStats[-1]['data']['summary']['CCFfreeavg']
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["Rwork"])
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(statsOverall["Rfree"])
                self.container.outputData.PERFORMANCEINDICATOR.CCFwork_avg.set(statsOverall["CCFworkavg"])
                self.container.outputData.PERFORMANCEINDICATOR.CCFfree_avg.set(statsOverall["CCFfreeavg"])
            except:
                try:
                    statsOverall["R1work"] = jsonStats[-1]['data']['summary']['R1work']
                    statsOverall["R1free"] = jsonStats[-1]['data']['summary']['R1free']
                    statsOverall["CCIworkavg"] = jsonStats[-1]['data']['summary']['CCIworkavg']
                    statsOverall["CCIfreeavg"] = jsonStats[-1]['data']['summary']['CCIfreeavg']
                    self.container.outputData.PERFORMANCEINDICATOR.R1Factor.set(statsOverall["R1work"])
                    self.container.outputData.PERFORMANCEINDICATOR.R1Free.set(statsOverall["R1free"])
                    self.container.outputData.PERFORMANCEINDICATOR.CCIwork_avg.set(statsOverall["CCIworkavg"])
                    self.container.outputData.PERFORMANCEINDICATOR.CCIfree_avg.set(statsOverall["CCIfreeavg"])
                except:
                    pass
        else:
            try:
                statsOverall["R"] = jsonStats[-1]['data']['summary']['R']
                statsOverall["CCFavg"] = jsonStats[-1]['data']['summary']['CCFavg']
                self.container.outputData.PERFORMANCEINDICATOR.R.set(statsOverall["R"])
                self.container.outputData.PERFORMANCEINDICATOR.CCF_avg.set(statsOverall["CCFavg"])
            except:
                try:
                    statsOverall["R1"] = jsonStats[-1]['data']['summary']['R1']
                    statsOverall["CCIavg"] = jsonStats[-1]['data']['summary']['CCIavg']
                    self.container.outputData.PERFORMANCEINDICATOR.R1.set(statsOverall["R1"])
                    self.container.outputData.PERFORMANCEINDICATOR.CCI_avg.set(statsOverall["CCIavg"])
                except:
                    pass

        et = ET.ElementTree(self.xmlroot)
        et.write(self.makeFileName('PROGRAMXML'))
        return CPluginScript.SUCCEEDED


    def makeCommandAndScript(self):
        if str(self.container.controlParameters.DATA_METHOD) == "spa":
            # options only for servalcat refine_spa_norefmac
            self.appendCommandLine(['refine_spa_norefmac'])
            self.appendCommandLine(['--halfmaps',
                                    str(self.container.inputData.MAPIN1.fullPath),
                                    str(self.container.inputData.MAPIN2.fullPath)])
            if self.container.inputData.MAPMASK.isSet():
                self.appendCommandLine(['--mask_for_fofc',
                                        str(self.container.inputData.MAPMASK.fullPath)])
            if self.container.controlParameters.MASK_RADIUS.isSet():
                self.appendCommandLine(['-r', str(float(str(self.container.controlParameters.MASK_RADIUS)))])
            self.appendCommandLine(['-d', str(self.container.controlParameters.RES_MIN)])
            self.appendCommandLine(['--source', "electron"])
            if self.container.controlParameters.PIXEL_SIZE.isSet():
                self.appendCommandLine(['--pixel_size', str(self.container.controlParameters.PIXEL_SIZE)])
            if self.container.controlParameters.POINTGROUP.isSet():
                self.appendCommandLine(['--pg', str(self.container.controlParameters.POINTGROUP)])
            if self.container.controlParameters.TWIST.isSet():
                self.appendCommandLine(['--twist', str(self.container.controlParameters.TWIST)])
            if self.container.controlParameters.RISE.isSet():
                self.appendCommandLine(['--rise', str(self.container.controlParameters.RISE)])
            if self.container.controlParameters.CENTER_X.isSet() and \
                    self.container.controlParameters.CENTER_Y.isSet() and \
                    self.container.controlParameters.CENTER_Z.isSet():
                self.appendCommandLine(['--center',
                                        str(self.container.controlParameters.CENTER_X),
                                        str(self.container.controlParameters.CENTER_Y),
                                        str(self.container.controlParameters.CENTER_Z)])
            if self.container.controlParameters.AXIS1_X.isSet() and \
                    self.container.controlParameters.AXIS1_Y.isSet() and \
                    self.container.controlParameters.AXIS1_Z.isSet():
                self.appendCommandLine(['--axis1',
                                        str(self.container.controlParameters.AXIS1_X),
                                        str(self.container.controlParameters.AXIS1_Y),
                                        str(self.container.controlParameters.AXIS1_Z)])
            if self.container.controlParameters.AXIS2_X.isSet() and \
                    self.container.controlParameters.AXIS2_Y.isSet() and \
                    self.container.controlParameters.AXIS2_Z.isSet():
                self.appendCommandLine(['--axis2',
                                        str(self.container.controlParameters.AXIS2_X),
                                        str(self.container.controlParameters.AXIS2_Y),
                                        str(self.container.controlParameters.AXIS2_Z)])
            if self.container.controlParameters.IGNORE_SYMMETRY:
                self.appendCommandLine(['--ignore_symmetry'])
            if self.container.controlParameters.PIXEL_SIZE.isSet():
                self.appendCommandLine(['--pixel_size', str(float(str(self.container.controlParameters.PIXEL_SIZE)))])
            if self.container.controlParameters.CROSS_VALIDATION:
                self.appendCommandLine(['--cross_validation'])
            if self.container.controlParameters.BLURUSE and self.container.controlParameters.BLUR.isSet():
                self.appendCommandLine(['--blur', str(self.container.controlParameters.BLUR)])
            self.hklout = ""

        elif str(self.container.controlParameters.DATA_METHOD) == "xtal":
            # options only for servalcat refine_xtal_norefmac
            self.appendCommandLine(['refine_xtal_norefmac'])
            self.appendCommandLine(['--hklin', self.hklin])
            if self.container.inputData.FREERFLAG.isSet():
                if self.container.controlParameters.FREERFLAG_NUMBER.isSet():
                    self.appendCommandLine(['--free', str(self.container.controlParameters.FREERFLAG_NUMBER)])
            if self.container.controlParameters.USE_TWIN:  # I,SIGI for optimal results
                self.appendCommandLine(['--twin'])
            self.hklout = os.path.join(self.workDirectory, "refined.mtz")
            # Resolution
            if self.container.controlParameters.RES_CUSTOM:
                if self.container.controlParameters.RES_MIN.isSet():
                    self.appendCommandLine(['--d_min', str(self.container.controlParameters.RES_MIN)])
                if self.container.controlParameters.RES_MAX.isSet():
                    self.appendCommandLine(['--d_max', str(self.container.controlParameters.RES_MAX)])
            self.appendCommandLine(['--source', str(self.container.controlParameters.SCATTERING_FACTORS)])
            if self.container.controlParameters.USE_WORK_IN_EST:
                self.appendCommandLine(['--use_work_in_est'])
            if self.container.controlParameters.NO_SOLVENT:
                self.appendCommandLine(['--no_solvent'])
            if self.container.controlParameters.REFINE_DFRAC:
                self.appendCommandLine(['--refine_dfrac'])

        # options for both servalcat refine_xtal_norefmac and servalcat refine_spa_norefmac
        self.appendCommandLine(['--model', self.inputCoordPath])
        self.appendCommandLine(['-o', 'refined'])
        if self.container.controlParameters.NCYCLES.isSet():
            self.appendCommandLine(["--ncycle", str(self.container.controlParameters.NCYCLES)])
        # Hydrogens
        if self.container.controlParameters.HYDR_USE:
            self.appendCommandLine(['--hydrogen', str(self.container.controlParameters.HYDR_ALL)])
            if self.container.controlParameters.H_REFINE:
                self.appendCommandLine(['--refine_h'])
        else:
            self.appendCommandLine(['--hydrogen', 'no'])
        if self.container.controlParameters.H_OUT:
            self.appendCommandLine(['--hout'])
        if self.container.inputData.DICT_LIST.isSet():
            if len(self.container.inputData.DICT_LIST) > 0:
                self.appendCommandLine(['--ligand'])
                for dictionary in self.container.inputData.DICT_LIST:
                    if os.path.exists(str(dictionary)):
                        self.appendCommandLine([str(dictionary)])
                    else:
                        print(f"WARNING: input restraint dictionary {dictionary} does not exist and so will not be used.")
        # Weight
        if str(self.container.controlParameters.WEIGHT_OPT) == "MANUAL":
            if self.container.controlParameters.WEIGHT.isSet():
                self.appendCommandLine(['--weight', str(self.container.controlParameters.WEIGHT)])
        else:  # WEIGHT_OPT == "AUTO":
            if self.container.controlParameters.WEIGHT_NO_ADJUST:
                self.appendCommandLine(['--no_weight_adjust'])
            elif self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN.isSet() and \
                    self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX.isSet():
                self.appendCommandLine(['--target_bond_rmsz_range',
                                        str(self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN),
                                        str(self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX)])
        # ADP
        if str(self.container.controlParameters.B_REFINEMENT_MODE) == "aniso":
            self.appendCommandLine(['--adp', 'aniso'])
        elif str(self.container.controlParameters.B_REFINEMENT_MODE) == "fix":
            self.appendCommandLine(['--adp', 'fix'])
        else:
            self.appendCommandLine(['--adp', 'iso'])
        # Jelly body
        if self.container.controlParameters.USE_JELLY:
            self.appendCommandLine(['--jellybody', '--jellybody_params'])
            if self.container.controlParameters.JELLY_SIGMA.isSet():
                self.appendCommandLine([str(self.container.controlParameters.JELLY_SIGMA)])
            else:
                self.appendCommandLine(['0.01'])  # default
            if self.container.controlParameters.JELLY_DIST.isSet():
                self.appendCommandLine([str(self.container.controlParameters.JELLY_DIST)])
            else:
                self.appendCommandLine(['4.2'])  # default
            if self.container.controlParameters.JELLY_ONLY:
                self.appendCommandLine(['--jellyonly'])
        if self.container.controlParameters.UNRESTRAINED:
            self.appendCommandLine(['--unrestrained'])
        if self.container.controlParameters.USE_NCS:
            self.appendCommandLine(['--ncsr'])
        if self.container.controlParameters.BFACSETUSE and \
                str(self.container.controlParameters.BFACSET):
            self.appendCommandLine(['--bfactor', str(self.container.controlParameters.BFACSET)])
        if float(self.container.controlParameters.ADPR_WEIGHT) != "1":
            self.appendCommandLine(['--adpr_weight', str(self.container.controlParameters.ADPR_WEIGHT)])
        if self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT.isSet():
            self.appendCommandLine(['--max_dist_for_adp_restraint', str(self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT)])
        if self.container.controlParameters.ADP_RESTRAINT_POWER.isSet():
            self.appendCommandLine(['--adp_restraint_power', str(self.container.controlParameters.ADP_RESTRAINT_POWER)])
        if self.container.controlParameters.ADP_RESTRAINT_NO_LONG_RANGE:
            self.appendCommandLine(['--adp_restraint_no_long_range'])
        if self.container.controlParameters.FIND_LINKS:
            self.appendCommandLine(['--find_links'])
        if self.container.controlParameters.FIX_XYZ:
            self.appendCommandLine(['--fix_xyz'])
        if self.container.controlParameters.KEEP_CHARGES:
            self.appendCommandLine(['--keep_charges'])
        if self.container.controlParameters.RANDOMIZEUSE and \
                str(self.container.controlParameters.RANDOMIZE):
            self.appendCommandLine(['--randomize', str(self.container.controlParameters.RANDOMIZE)])

        keywordFilePath = str(os.path.join(self.getWorkDirectory(), 'keywords.txt'))

        # Occupancy refinement
        if self.container.controlParameters.OCCUPANCY_GROUPS:
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.OCCUPANCY_REFINEMENT:
                    keywordFile.write("OCCUPANCY REFINE\n")
                if self.container.controlParameters.OCCUPANCY_NCYCLE.isSet():
                    keywordFile.write("OCCUPANCY REFINE NCYCLE " + str(self.container.controlParameters.OCCUPANCY_NCYCLE) + "\n")
                occup_groupids = []
                for sel0 in self.container.controlParameters.OCCUPANCY_SELECTION:
                    sel = sel0.get()
                    if sel['groupId'] and sel['chainIds'] and sel['firstRes'] and sel['lastRes']:
                        occupText = "OCCUPANCY GROUP ID %s CHAIN %s RESIDUE FROM %s TO %s" % (str(sel['groupId']),str(sel['chainIds']),str(sel['firstRes']),str(sel['lastRes']))
                        if sel['atoms']:
                            occupText += " ATOM %s" % (str(sel['atoms']))
                        if sel['alt']:
                            occupText += " ALT %s" % (str(sel['alt']))
                        keywordFile.write(occupText + '\n')
                        occup_groupids.append(sel['groupId'])
                if self.container.controlParameters.OCCUPANCY_COMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE:
                        sel = sel0.get()
                        occupText = "OCCUPANCY GROUP ALTS COMPLETE %s" % (str(sel['groupIds']))
                        keywordFile.write(occupText + '\n')  # self.appendCommandScript(occupText + '\n')
                if self.container.controlParameters.OCCUPANCY_INCOMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE:
                        sel = sel0.get()
                        occupText = "OCCUPANCY GROUP ALTS INCOMPLETE %s" % (str(sel['groupIds']))
                        keywordFile.write(occupText + '\n')  # self.appendCommandScript(occupText + '\n')

        if self.container.inputData.METALCOORD_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write("\n@%s"%(str(self.container.inputData.METALCOORD_RESTRAINTS.fullPath)))
        if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.PROSMART_PROTEIN_SGMN.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMN %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_SGMN)))
                if self.container.controlParameters.PROSMART_PROTEIN_SGMX.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_SGMX)))
                if self.container.controlParameters.PROSMART_PROTEIN_ALPHA.isSet():
                    keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_ALPHA)))
                if self.container.controlParameters.PROSMART_PROTEIN_DMAX.isSet():
                    keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_DMAX)))
                keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.fullPath)))
        if self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMN.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMN %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_SGMN)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMX.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_SGMX)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA.isSet():
                    keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_DMAX.isSet():
                    keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_DMAX)))
                keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.fullPath)))
        if self.container.inputData.SERVALCAT_KEYWORD_FILE.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write("\n@%s"%(str(self.container.inputData.SERVALCAT_KEYWORD_FILE.fullPath)))

        if os.path.isfile(keywordFilePath):
            self.appendCommandLine(["--keyword_file", keywordFilePath])
        if self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS.isSet():
            self.appendCommandLine(str(self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS).split())

        jsonFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.json"))
        self.watchFile(jsonFilePath, handler=self.handleJsonChanged, unwatchWhileHandling=True)
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(str)
    def handleJsonChanged(self, jsonFilePath):
        self.xmlroot.clear()
        # Use output JSON from servalcat as the basis for output XML
        if os.path.isfile(jsonFilePath):
            with open(jsonFilePath, "r") as jsonFile:
                jsonText = jsonFile.read()
            try:
                # Get statistics from JSON, convert to XML, save and load to self.xmlroot
                jsonStats = list(json.loads(jsonText))
                # add d_max_4ssqll and d_min_4ssqll
                for i in range(len(jsonStats)):
                    for j in range(len(jsonStats[i]["data"]["binned"])):
                        jsonStats[i]["data"]["binned"][j]['d_min_4ssqll'] = \
                            1 / (list(jsonStats)[i]["data"]["binned"][j]['d_min'] * list(jsonStats)[i]["data"]["binned"][j]['d_min'])
                        jsonStats[i]["data"]["binned"][j]['d_max_4ssqll'] = \
                            1 / (list(jsonStats)[i]["data"]["binned"][j]['d_max'] * list(jsonStats)[i]["data"]["binned"][j]['d_max'])
                xmlText = json2xml(list(jsonStats), tag_name_subroot="cycle")
                xmlFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.xml"))
                xmlText = self.xmlAddRoot(xmlText, xmlFilePath, xmlRootName="SERVALCAT")
                self.xmlroot = ET.fromstring(xmlText)
                et = ET.ElementTree(self.xmlroot)
                # Write program.xml_tmp and move it to program.xml
                self.flushXml()
            except Exception:
                import traceback
                print(traceback.format_exc())

    def flushXml(self):  # assumes self.xmlroot and self.xmlLength are well set
        # Save program.xml if program.xml_tmp has grown
        newXml = ET.tostring(self.xmlroot)
        if len(newXml) > self.xmlLength:
            self.xmlLength = len(newXml)
            with open(self.makeFileName('PROGRAMXML')+'_tmp','w') as programXmlFile:
                CCP4Utils.writeXML(programXmlFile, newXml)
            shutil.move(self.makeFileName('PROGRAMXML')+'_tmp', self.makeFileName('PROGRAMXML'))

    def setProgramVersion(self):
      print('refmac.getProgramVersion')
      return CPluginScript.setProgramVersion(self,'Refmac_5')
