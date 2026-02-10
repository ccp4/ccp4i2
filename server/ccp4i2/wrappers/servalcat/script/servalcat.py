import json
import os
import pathlib
import shutil
import traceback
from xml.etree import ElementTree as ET

from ccp4i2.core import CCP4ErrorHandling, CCP4Modules, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript

from .json2xml import json2xml


class servalcat(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'Refinement (servalcat)'
    TASKNAME = 'servalcat'
    TASKCOMMAND = 'servalcat'
    PERFORMANCECLASS = 'CServalcatPerformance'

    ERROR_CODES = {
        201: {'description': 'Servalcat returned with non zero status'},
        202: {'description': 'New library created but strictly required'},
        203: {'description': 'New library created',
              'severity': CCP4ErrorHandling.SEVERITY_WARNING},
        204: {'description': 'Program completed without generating output statistics'},
        205: {'description': 'Failed to parse output JSON statistics'},
        206: {'description': 'Failed to read output MTZ file'},
        207: {'description': 'Failed to split HKL output'},
    }

    def __init__(self, *args, **kwargs):
        super(servalcat, self).__init__(*args, **kwargs)
        self.xmlroot = ET.Element('SERVALCAT')
        self.xmlLength = 0

    def xmlAddRoot(self, xmlText, xmlFilePath=None, xmlRootName=None):
        if xmlRootName:
            xmlText = f"<{xmlRootName}>\n{xmlText}\n</{xmlRootName}>"
        if xmlFilePath:
            with open(xmlFilePath, 'w') as programXmlFile:
                programXmlFile.write(xmlText)
        return xmlText

    def processInputFiles(self):
        self.hklin = None
        dataObjects = []

        # Apply coordinate selection if set
        self.inputCoordPath = os.path.normpath(
            str(self.container.inputData.XYZIN.fullPath))
        if self.container.inputData.XYZIN.isSelectionSet():
            self.inputCoordPath = os.path.normpath(
                os.path.join(self.getWorkDirectory(), 'selected.pdb'))
            self.container.inputData.XYZIN.loadFile()
            if self.container.inputData.XYZIN.isMMCIF():
                self.inputCoordPath = str(
                    pathlib.Path(self.inputCoordPath).with_suffix('.cif'))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(
                self.inputCoordPath)

        # Create DICT by merging dictionaries in DICT_LIST
        self.joinDicts(self.container.outputData.DICT,
                       self.container.inputData.DICT_LIST)

        # Prepare merged HKL input for X-ray mode
        if str(self.container.controlParameters.DATA_METHOD) == 'xtal':
            if str(self.container.controlParameters.MERGED_OR_UNMERGED) == "merged":
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
                    if str(self.container.controlParameters.F_SIGF_OR_I_SIGI) == "F_SIGF":
                        obsTypeRoot = 'CONTENT_FLAG_F'
                if str(self.container.controlParameters.SCATTERING_FACTORS) in ["electron", "neutron"]:
                    obsPairOrMean = 'MEAN'
                obsType = getattr(CCP4XtalData.CObsDataFile, obsTypeRoot + obsPairOrMean)
                dataObjects += [['HKLIN', obsType]]

                if self.container.inputData.FREERFLAG.isSet():
                    dataObjects += ['FREERFLAG']
                self.hklin, error = self.makeHklin(dataObjects)
                if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                    return CPluginScript.FAILED

            return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        if hasattr(self, 'logFileHandle'):
            self.logFileHandle.write("JOB TITLE SECTION\n")
            try:
                jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(
                    jobId=self.jobId)
                if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                    self.logFileHandle.write(str(jobInfo["jobtitle"]) + "\n")
                while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
                    jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(
                        jobId=jobInfo["parentjobid"])
                    if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                        self.logFileHandle.write(str(jobInfo["jobtitle"]) + "\n")
            except Exception as e:
                print(f"[servalcat] Warning: could not write job title: {e}")
            finally:
                self.logFileHandle.close()
        else:
            self.xmlroot.clear()

        # Check exit status of the program
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER
        try:
            exitStatus = PROCESSMANAGER().getJobData(
                pid=self.getProcessId(), attribute='exitStatus')
        except Exception as e:
            self.appendErrorReport(201,
                f'Unable to recover exitStatus: {e}')
            return CPluginScript.FAILED
        if exitStatus != 0:
            self.appendErrorReport(201, f'Exit status: {exitStatus}')
            return CPluginScript.FAILED

        try:
            exitCode = PROCESSMANAGER().getJobData(
                pid=self.getProcessId(), attribute='exitCode')
        except Exception as e:
            self.appendErrorReport(201,
                f'Unable to recover exitCode: {e}')
            return CPluginScript.FAILED
        if exitCode != 0:
            self.appendErrorReport(201, f'Exit code: {exitCode}')
            return CPluginScript.FAILED

        # Set output file paths and annotations
        outputCifPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), 'refined.mmcif'))
        self.container.outputData.CIFFILE.setFullPath(outputCifPath)
        self.container.outputData.CIFFILE.annotation.set(
            'Model from refinement (mmCIF format)')

        outputPdbPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), 'refined.pdb'))
        if os.path.isfile(outputPdbPath):
            self.container.outputData.XYZOUT.setFullPath(outputPdbPath)
            self.container.outputData.XYZOUT.annotation.set(
                'Model from refinement (PDB format)')

        self.container.outputData.FPHIOUT.annotation.set(
            'Density map (Fourier coeff.)')
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.annotation.set(
            'Difference density map (Fourier coeff.)')
        self.container.outputData.DIFFPHIOUT.subType = 2

        outputFiles = ['FPHIOUT', 'DIFFPHIOUT']
        outputColumns = ['FWT,PHWT', 'DELFWT,PHDELWT']

        # Read output MTZ
        if str(self.container.controlParameters.DATA_METHOD) == "xtal":
            hkloutFilePath = str(os.path.join(
                self.getWorkDirectory(), "refined.mtz"))
        else:  # spa
            hkloutFilePath = str(os.path.join(
                self.getWorkDirectory(), "refined_diffmap.mtz"))

        try:
            hkloutFile = CCP4XtalData.CMtzDataFile(hkloutFilePath)
            hkloutFile.loadFile()
            columnLabelsInFile = [
                str(column.columnLabel) for column in hkloutFile.fileContent.listOfColumns
            ]
        except Exception as e:
            self.appendErrorReport(206,
                f'Failed to read output MTZ: {e}\n{traceback.format_exc()}')
            return CPluginScript.FAILED

        if str(self.container.controlParameters.DATA_METHOD) == "xtal":
            if 'FAN' in columnLabelsInFile and 'PHAN' in columnLabelsInFile:
                self.container.outputData.ANOMFPHIOUT.annotation.set(
                    'Anomalous difference map')
                outputFiles += ['ANOMFPHIOUT']
                outputColumns += ['FAN,PHAN']

        # Split HKL output into individual map coefficient files
        error = self.splitHklout(outputFiles, outputColumns, hkloutFilePath)
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            self.appendErrorReport(207,
                f'Failed to split HKL output')
            return CPluginScript.FAILED

        # Handle SPA-specific outputs
        if str(self.container.controlParameters.DATA_METHOD) == 'spa':
            self.container.outputData.MAP_FO.annotation.set(
                'Density map (in real space)')
            outputMapFoPath = os.path.normpath(os.path.join(
                self.getWorkDirectory(), 'refined_diffmap_normalized_fo.mrc'))
            self.container.outputData.MAP_FO.setFullPath(outputMapFoPath)

            self.container.outputData.MAP_FOFC.annotation.set(
                'Difference density map (in real space)')
            outputMapFoFcPath = os.path.normpath(os.path.join(
                self.getWorkDirectory(), 'refined_diffmap_normalized_fofc.mrc'))
            self.container.outputData.MAP_FOFC.setFullPath(outputMapFoFcPath)

            # Write Coot script with contour levels
            cootScriptI2FilePath = os.path.join(
                self.getWorkDirectory(), "refined_coot_i2.py")
            with open(cootScriptI2FilePath, "w") as cootScriptI2File:
                cootScriptI2File.write("set_contour_level_absolute(1, 1.2)\n")
                cootScriptI2File.write("set_contour_level_absolute(2, 3.0)\n")
            self.container.outputData.COOTSCRIPTOUT.annotation.set('Coot script')
            self.container.outputData.COOTSCRIPTOUT.setFullPath(cootScriptI2FilePath)

        # Parse output JSON statistics and convert to XML
        jsonFilePath = str(os.path.join(
            self.getWorkDirectory(), "refined_stats.json"))
        if not os.path.isfile(jsonFilePath):
            self.appendErrorReport(204,
                'refined_stats.json not found')
            return CPluginScript.FAILED

        try:
            with open(jsonFilePath, "r") as jsonFile:
                jsonText = jsonFile.read()
            jsonStats = json.loads(jsonText)

            # Add d_max_4ssqll and d_min_4ssqll for resolution plotting
            for stat in jsonStats:
                data = stat.get("data", {})
                for p in ("binned", "ml"):
                    if p in data:
                        for entry in data[p]:
                            d_min = entry.get('d_min')
                            d_max = entry.get('d_max')
                            if d_min:
                                entry['d_min_4ssqll'] = 1 / (d_min * d_min)
                            if d_max:
                                entry['d_max_4ssqll'] = 1 / (d_max * d_max)

            xmlText = json2xml(jsonStats, tag_name_subroot="cycle")
            xmlFilePath = str(os.path.join(
                self.getWorkDirectory(), "refined_stats.xml"))
            xmlText = self.xmlAddRoot(xmlText, xmlFilePath,
                                       xmlRootName="SERVALCAT")
            self.xmlroot = ET.fromstring(xmlText)
            self.flushXml()
        except Exception as e:
            self.appendErrorReport(205,
                f'Failed to parse output JSON: {e}\n{traceback.format_exc()}')
            return CPluginScript.FAILED

        # Extract performance indicators from JSON
        self._extractPerformanceIndicators(jsonStats)

        et = ET.ElementTree(self.xmlroot)
        et.write(self.makeFileName('PROGRAMXML'))
        return CPluginScript.SUCCEEDED

    def _extractPerformanceIndicators(self, jsonStats):
        """Extract R-factors and other indicators from final cycle stats."""
        try:
            summary = jsonStats[-1]['data']['summary']
        except (IndexError, KeyError):
            return

        # Log-likelihood
        if '-LL' in summary:
            pass  # Stored in XML already

        if str(self.container.controlParameters.DATA_METHOD) == 'spa':
            if 'FSCaverage' in summary:
                self.container.outputData.PERFORMANCEINDICATOR.FSCaverage.set(
                    summary['FSCaverage'])
        elif self.container.inputData.FREERFLAG.isSet():
            # Try F-based R-factors first, then I-based
            if 'Rwork' in summary:
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(
                    summary['Rwork'])
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(
                    summary.get('Rfree', 0))
                if 'CCFworkavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCFwork_avg.set(
                        summary['CCFworkavg'])
                if 'CCFfreeavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCFfree_avg.set(
                        summary['CCFfreeavg'])
            elif 'R1work' in summary:
                self.container.outputData.PERFORMANCEINDICATOR.R1Factor.set(
                    summary['R1work'])
                self.container.outputData.PERFORMANCEINDICATOR.R1Free.set(
                    summary.get('R1free', 0))
                if 'CCIworkavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCIwork_avg.set(
                        summary['CCIworkavg'])
                if 'CCIfreeavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCIfree_avg.set(
                        summary['CCIfreeavg'])
        else:
            # No FreeR flag - use overall R
            if 'R' in summary:
                self.container.outputData.PERFORMANCEINDICATOR.R.set(
                    summary['R'])
                if 'CCFavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCF_avg.set(
                        summary['CCFavg'])
            elif 'R1' in summary:
                self.container.outputData.PERFORMANCEINDICATOR.R1.set(
                    summary['R1'])
                if 'CCIavg' in summary:
                    self.container.outputData.PERFORMANCEINDICATOR.CCI_avg.set(
                        summary['CCIavg'])

    def makeCommandAndScript(self):
        if str(self.container.controlParameters.DATA_METHOD) == "spa":
            # SPA mode: servalcat refine_spa_norefmac
            self.appendCommandLine(['refine_spa_norefmac'])
            self.appendCommandLine(['--halfmaps',
                                    str(self.container.inputData.MAPIN1.fullPath),
                                    str(self.container.inputData.MAPIN2.fullPath)])
            if self.container.inputData.MAPMASK.isSet():
                self.appendCommandLine(['--mask_for_fofc',
                                        str(self.container.inputData.MAPMASK.fullPath)])
            if self.container.controlParameters.MASK_RADIUS.isSet():
                self.appendCommandLine(['-r', str(float(
                    str(self.container.controlParameters.MASK_RADIUS)))])
            self.appendCommandLine(['-d', str(self.container.controlParameters.RES_MIN)])
            self.appendCommandLine(['--source', "electron"])
            if self.container.controlParameters.PIXEL_SIZE.isSet():
                self.appendCommandLine(['--pixel_size',
                                        str(self.container.controlParameters.PIXEL_SIZE)])
            if self.container.controlParameters.POINTGROUP.isSet():
                self.appendCommandLine(['--pg',
                                        str(self.container.controlParameters.POINTGROUP)])
            if self.container.controlParameters.TWIST.isSet():
                self.appendCommandLine(['--twist',
                                        str(self.container.controlParameters.TWIST)])
            if self.container.controlParameters.RISE.isSet():
                self.appendCommandLine(['--rise',
                                        str(self.container.controlParameters.RISE)])
            if (self.container.controlParameters.CENTER_X.isSet() and
                    self.container.controlParameters.CENTER_Y.isSet() and
                    self.container.controlParameters.CENTER_Z.isSet()):
                self.appendCommandLine(['--center',
                                        str(self.container.controlParameters.CENTER_X),
                                        str(self.container.controlParameters.CENTER_Y),
                                        str(self.container.controlParameters.CENTER_Z)])
            if (self.container.controlParameters.AXIS1_X.isSet() and
                    self.container.controlParameters.AXIS1_Y.isSet() and
                    self.container.controlParameters.AXIS1_Z.isSet()):
                self.appendCommandLine(['--axis1',
                                        str(self.container.controlParameters.AXIS1_X),
                                        str(self.container.controlParameters.AXIS1_Y),
                                        str(self.container.controlParameters.AXIS1_Z)])
            if (self.container.controlParameters.AXIS2_X.isSet() and
                    self.container.controlParameters.AXIS2_Y.isSet() and
                    self.container.controlParameters.AXIS2_Z.isSet()):
                self.appendCommandLine(['--axis2',
                                        str(self.container.controlParameters.AXIS2_X),
                                        str(self.container.controlParameters.AXIS2_Y),
                                        str(self.container.controlParameters.AXIS2_Z)])
            if self.container.controlParameters.IGNORE_SYMMETRY:
                self.appendCommandLine(['--ignore_symmetry'])
            if self.container.controlParameters.CROSS_VALIDATION:
                self.appendCommandLine(['--cross_validation'])
            if (self.container.controlParameters.BLURUSE and
                    self.container.controlParameters.BLUR.isSet()):
                self.appendCommandLine(['--blur',
                                        str(self.container.controlParameters.BLUR)])
            self.hklout = ""

        elif str(self.container.controlParameters.DATA_METHOD) == "xtal":
            # X-ray mode: servalcat refine_xtal_norefmac
            self.appendCommandLine(['refine_xtal_norefmac'])
            self.appendCommandLine(['--hklin'])
            if str(self.container.controlParameters.MERGED_OR_UNMERGED) == "unmerged":
                self.appendCommandLine(
                    [str(self.container.inputData.HKLIN_UNMERGED.fullPath)])
                if self.container.inputData.FREERFLAG.isSet():
                    self.appendCommandLine(
                        ['--hklin_free',
                         str(self.container.inputData.FREERFLAG.fullPath)])
            else:
                self.appendCommandLine([self.hklin])
            if self.container.inputData.FREERFLAG.isSet():
                if self.container.controlParameters.FREERFLAG_NUMBER.isSet():
                    self.appendCommandLine(
                        ['--free',
                         str(self.container.controlParameters.FREERFLAG_NUMBER)])
            if self.container.controlParameters.USE_TWIN:
                self.appendCommandLine(['--twin'])
            self.hklout = os.path.join(self.workDirectory, "refined.mtz")
            # Resolution
            if self.container.controlParameters.RES_CUSTOM:
                if self.container.controlParameters.RES_MIN.isSet():
                    self.appendCommandLine(
                        ['--d_min', str(self.container.controlParameters.RES_MIN)])
                if self.container.controlParameters.RES_MAX.isSet():
                    self.appendCommandLine(
                        ['--d_max', str(self.container.controlParameters.RES_MAX)])
            self.appendCommandLine(
                ['--source',
                 str(self.container.controlParameters.SCATTERING_FACTORS)])
            if self.container.controlParameters.USE_WORK_IN_EST:
                self.appendCommandLine(['--use_in_est', 'work'])
            if self.container.controlParameters.NO_SOLVENT:
                self.appendCommandLine(['--no_solvent'])

        # Common options for both modes
        self.appendCommandLine(['--model', self.inputCoordPath])
        self.appendCommandLine(['-o', 'refined'])
        if self.container.controlParameters.NCYCLES.isSet():
            self.appendCommandLine(
                ["--ncycle", str(self.container.controlParameters.NCYCLES)])

        # Hydrogens
        if self.container.controlParameters.HYDR_USE:
            self.appendCommandLine(
                ['--hydrogen',
                 str(self.container.controlParameters.HYDR_ALL)])
            if self.container.controlParameters.H_REFINE:
                self.appendCommandLine(['--refine_h'])
        else:
            self.appendCommandLine(['--hydrogen', 'no'])
        if self.container.controlParameters.H_OUT:
            self.appendCommandLine(['--hout'])

        # Ligand dictionaries
        if self.container.inputData.DICT_LIST.isSet():
            if len(self.container.inputData.DICT_LIST) > 0:
                self.appendCommandLine(['--ligand'])
                for dictionary in self.container.inputData.DICT_LIST:
                    if os.path.exists(str(dictionary)):
                        self.appendCommandLine([str(dictionary)])
                    else:
                        print(f"[servalcat] WARNING: dictionary {dictionary} "
                              f"does not exist, skipping")

        # Weight
        if str(self.container.controlParameters.WEIGHT_OPT) == "MANUAL":
            if self.container.controlParameters.WEIGHT.isSet():
                self.appendCommandLine(
                    ['--weight',
                     str(self.container.controlParameters.WEIGHT)])
        else:  # AUTO
            if self.container.controlParameters.WEIGHT_NO_ADJUST:
                self.appendCommandLine(['--no_weight_adjust'])
            elif (self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN.isSet() and
                  self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX.isSet()):
                self.appendCommandLine(
                    ['--target_bond_rmsz_range',
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
            self.appendCommandLine([
                str(self.container.controlParameters.JELLY_SIGMA)
                if self.container.controlParameters.JELLY_SIGMA.isSet()
                else '0.01'
            ])
            self.appendCommandLine([
                str(self.container.controlParameters.JELLY_DIST)
                if self.container.controlParameters.JELLY_DIST.isSet()
                else '4.2'
            ])
            if self.container.controlParameters.JELLY_ONLY:
                self.appendCommandLine(['--jellyonly'])

        if self.container.controlParameters.UNRESTRAINED:
            self.appendCommandLine(['--unrestrained'])
        if self.container.controlParameters.USE_NCS:
            self.appendCommandLine(['--ncsr'])
        if (self.container.controlParameters.BFACSETUSE and
                str(self.container.controlParameters.BFACSET)):
            self.appendCommandLine(
                ['--bfactor', str(self.container.controlParameters.BFACSET)])
        if float(self.container.controlParameters.ADPR_WEIGHT) != 1.0:
            self.appendCommandLine(
                ['--adpr_weight',
                 str(self.container.controlParameters.ADPR_WEIGHT)])
        if self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT.isSet():
            self.appendCommandLine(
                ['--max_dist_for_adp_restraint',
                 str(self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT)])
        if self.container.controlParameters.ADP_RESTRAINT_POWER.isSet():
            self.appendCommandLine(
                ['--adp_restraint_power',
                 str(self.container.controlParameters.ADP_RESTRAINT_POWER)])
        if self.container.controlParameters.ADP_RESTRAINT_NO_LONG_RANGE:
            self.appendCommandLine(['--adp_restraint_no_long_range'])
        if self.container.controlParameters.FIND_LINKS:
            self.appendCommandLine(['--find_links'])
        if self.container.controlParameters.FIX_XYZ:
            self.appendCommandLine(['--fix_xyz'])
        if self.container.controlParameters.KEEP_CHARGES:
            self.appendCommandLine(['--keep_charges'])
        if (self.container.controlParameters.RANDOMIZEUSE and
                str(self.container.controlParameters.RANDOMIZE)):
            self.appendCommandLine(
                ['--randomize',
                 str(self.container.controlParameters.RANDOMIZE)])

        # Build keyword file for external restraints and advanced options
        keywordFilePath = str(os.path.join(
            self.getWorkDirectory(), 'keywords.txt'))

        if str(float(self.container.controlParameters.VDWR_WEIGHT)) != "1":
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write(
                    f"VDWR {self.container.controlParameters.VDWR_WEIGHT}\n")

        # Occupancy refinement
        if self.container.controlParameters.OCCUPANCY_GROUPS:
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.OCCUPANCY_REFINEMENT:
                    keywordFile.write("OCCUPANCY REFINE\n")
                if self.container.controlParameters.OCCUPANCY_NCYCLE.isSet():
                    keywordFile.write(
                        f"OCCUPANCY REFINE NCYCLE "
                        f"{self.container.controlParameters.OCCUPANCY_NCYCLE}\n")
                for sel0 in self.container.controlParameters.OCCUPANCY_SELECTION:
                    sel = sel0.get()
                    if (sel['groupId'] and sel['chainIds'] and
                            sel['firstRes'] and sel['lastRes']):
                        occupText = (
                            f"OCCUPANCY GROUP ID {sel['groupId']} "
                            f"CHAIN {sel['chainIds']} "
                            f"RESIDUE FROM {sel['firstRes']} TO {sel['lastRes']}")
                        if sel['atoms']:
                            occupText += f" ATOM {sel['atoms']}"
                        if sel['alt']:
                            occupText += f" ALT {sel['alt']}"
                        keywordFile.write(occupText + '\n')
                if self.container.controlParameters.OCCUPANCY_COMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE:
                        sel = sel0.get()
                        keywordFile.write(
                            f"OCCUPANCY GROUP ALTS COMPLETE {sel['groupIds']}\n")
                if self.container.controlParameters.OCCUPANCY_INCOMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE:
                        sel = sel0.get()
                        keywordFile.write(
                            f"OCCUPANCY GROUP ALTS INCOMPLETE {sel['groupIds']}\n")

        if self.container.inputData.METALCOORD_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write(
                    f"\n@{self.container.inputData.METALCOORD_RESTRAINTS.fullPath}")

        if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.PROSMART_PROTEIN_SGMN.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL WEIGHT SGMN "
                        f"{self.container.controlParameters.PROSMART_PROTEIN_SGMN}")
                if self.container.controlParameters.PROSMART_PROTEIN_SGMX.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL WEIGHT SGMX "
                        f"{self.container.controlParameters.PROSMART_PROTEIN_SGMX}")
                if self.container.controlParameters.PROSMART_PROTEIN_ALPHA.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL ALPHA "
                        f"{self.container.controlParameters.PROSMART_PROTEIN_ALPHA}")
                if self.container.controlParameters.PROSMART_PROTEIN_DMAX.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL DMAX "
                        f"{self.container.controlParameters.PROSMART_PROTEIN_DMAX}")
                keywordFile.write(
                    f"\n@{self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.fullPath}")

        if self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMN.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL WEIGHT SGMN "
                        f"{self.container.controlParameters.PROSMART_NUCLEICACID_SGMN}")
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMX.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL WEIGHT SGMX "
                        f"{self.container.controlParameters.PROSMART_NUCLEICACID_SGMX}")
                if self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL ALPHA "
                        f"{self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA}")
                if self.container.controlParameters.PROSMART_NUCLEICACID_DMAX.isSet():
                    keywordFile.write(
                        f"\nEXTERNAL DMAX "
                        f"{self.container.controlParameters.PROSMART_NUCLEICACID_DMAX}")
                keywordFile.write(
                    f"\n@{self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.fullPath}")

        if self.container.inputData.SERVALCAT_KEYWORD_FILE.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write(
                    f"\n@{self.container.inputData.SERVALCAT_KEYWORD_FILE.fullPath}")

        if os.path.isfile(keywordFilePath):
            self.appendCommandLine(["--keyword_file", keywordFilePath])

        if self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS.isSet():
            self.appendCommandLine(
                str(self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS).split())

        return CPluginScript.SUCCEEDED

    def startProcess(self):
        """Set up real-time progress monitoring, then start the external process."""
        jsonFilePath = os.path.join(self.getWorkDirectory(), "refined_stats.json")
        self.watchFile(jsonFilePath, handler=self.handleJsonChanged,
                       unwatchWhileHandling=True)
        return super().startProcess()

    def handleJsonChanged(self, jsonFilePath):
        """Parse servalcat's output JSON and update program.xml in real-time."""
        self.xmlroot.clear()
        if os.path.isfile(jsonFilePath):
            try:
                with open(jsonFilePath, "r") as f:
                    jsonText = f.read()
                jsonStats = json.loads(jsonText)
                for stat in jsonStats:
                    data = stat.get("data", {})
                    for p in ("binned", "ml"):
                        if p in data:
                            for entry in data[p]:
                                d_min = entry.get('d_min')
                                d_max = entry.get('d_max')
                                if d_min:
                                    entry['d_min_4ssqll'] = 1 / (d_min * d_min)
                                if d_max:
                                    entry['d_max_4ssqll'] = 1 / (d_max * d_max)
                xmlText = json2xml(jsonStats, tag_name_subroot="cycle")
                xmlFilePath = os.path.join(
                    self.getWorkDirectory(), "refined_stats.xml")
                xmlText = self.xmlAddRoot(xmlText, xmlFilePath,
                                          xmlRootName="SERVALCAT")
                self.xmlroot = ET.fromstring(xmlText)
                self.flushXml()
            except Exception:
                traceback.print_exc()

    def flushXml(self):
        """Save program.xml if content has grown (atomic write via tmp+move)."""
        newXml = ET.tostring(self.xmlroot)
        if len(newXml) > self.xmlLength:
            self.xmlLength = len(newXml)
            tmpPath = self.makeFileName('PROGRAMXML') + '_tmp'
            with open(tmpPath, 'w') as f:
                CCP4Utils.writeXML(f, newXml)
            shutil.move(tmpPath, self.makeFileName('PROGRAMXML'))

    def setProgramVersion(self):
        return CPluginScript.setProgramVersion(self, 'Refmac_5')
