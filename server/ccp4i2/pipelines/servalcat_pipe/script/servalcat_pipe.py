import os
import shutil
import sys
import threading
import traceback
from operator import itemgetter

import gemmi
import numpy
from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4ErrorHandling import CErrorReport
from ccp4i2.core.CCP4PluginScript import CPluginScript

from . import monitor_differences


class servalcat_pipe(CPluginScript):

    SHORTTASKTITLE = 'Servalcat'
    TASKTITLE = 'Refinement against diffraction data with optional restraints (ProSMART, MetalCoord)'
    TASKNAME = 'servalcat_pipe'
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'
    WHATNEXT = ['servalcat_pipe', 'coot_rebuild', 'modelcraft']
    ASYNCHRONOUS = True
    PERFORMANCECLASS = 'CServalcatPerformance'
    PURGESEARCHLIST = [
        ['refmac%*/hklout.mtz', 0, "hklout"],
        ['refmac%*/hklout.mtz', 7, "hklout"],
        ['*%*/ANOMFPHIOUT.mtz', 1, "ANOMFPHIOUT"],
        ['*%*/DIFANOMFPHIOUT.mtz', 1, "DIFANOMFPHIOUT"],
    ]

    ERROR_CODES = {
        101: {'description': 'Error copying data file from final job to pipeline directory'},
        102: {'description': 'ProSMART protein restraints failed'},
        103: {'description': 'ProSMART nucleic acid restraints failed'},
        104: {'description': 'MetalCoord restraints failed'},
        105: {'description': 'Servalcat refinement failed'},
        106: {'description': 'Servalcat output file not created'},
        107: {'description': 'Coot find waters failed'},
        108: {'description': 'Post-coot servalcat failed'},
        109: {'description': 'Failed to create sub-plugin'},
        110: {'description': 'Validation failed'},
        111: {'description': 'Failed to convert MetalCoord JSON to restraints'},
        112: {'description': 'Missing required output from sub-plugin'},
        113: {'description': 'ADP analysis failed'},
        114: {'description': 'Coordinate/ADP deviation analysis failed'},
    }

    def __init__(self, *args, **kws):
        super(servalcat_pipe, self).__init__(*args, **kws)
        self.xmlroot = etree.Element("SERVALCAT")

        # Sub-plugin references - for processOutputFiles() to access
        self.prosmartProteinPlugin = None
        self.prosmartNucleicAcidPlugin = None
        self.metalCoordPlugins = []
        self.servalcatPlugin = None
        self.cootPlugin = None
        self.servalcatPostCootPlugin = None
        self.validatePlugin = None

    def validity(self):
        """Validate plugin, adjusting qualifiers for embedded wrappers.

        The metalCoordWrapper embeds metalCoord, which requires XYZIN when run
        standalone. In the pipeline context, XYZIN is programmatically filled
        from the pipeline's own XYZIN, so we set allowUndefined=True.
        """
        if hasattr(self.container, 'metalCoordWrapper'):
            self.container.metalCoordWrapper.inputData.XYZIN.set_qualifier(
                'allowUndefined', True
            )
        return super(servalcat_pipe, self).validity()

    # =========================================================================
    # Main pipeline orchestration
    # =========================================================================

    def startProcess(self):
        """Execute the pipeline synchronously through all phases.

        Phase 1: ProSMART protein restraints (optional)
        Phase 2: ProSMART nucleic acid restraints (optional)
        Phase 3: MetalCoord restraints (optional)
        Phase 4: Servalcat refinement (required)
        Phase 5: Water addition + re-refinement (optional)
        Phase 6: Validation and analysis (optional)
        """
        error = CErrorReport()

        # =================================================================
        # Phase 1: ProSMART protein restraints
        # =================================================================
        phase1_error = self._runProsmartProtein()
        if phase1_error and phase1_error.maxSeverity() >= 4:
            return phase1_error

        # =================================================================
        # Phase 2: ProSMART nucleic acid restraints
        # =================================================================
        phase2_error = self._runProsmartNucleicAcid()
        if phase2_error and phase2_error.maxSeverity() >= 4:
            return phase2_error

        # =================================================================
        # Phase 3: MetalCoord restraints
        # =================================================================
        phase3_error = self._runMetalCoords()
        if phase3_error and phase3_error.maxSeverity() >= 4:
            return phase3_error

        # =================================================================
        # Phase 4: Servalcat refinement
        # =================================================================
        phase4_error = self._runServalcat()
        if phase4_error and phase4_error.maxSeverity() >= 4:
            return phase4_error

        # =================================================================
        # Phase 5: Water addition + re-refinement (optional)
        # =================================================================
        phase5_error = self._runWaterAddition()
        if phase5_error and phase5_error.maxSeverity() >= 4:
            return phase5_error

        # =================================================================
        # Phase 6: Validation and analysis (optional)
        # =================================================================
        phase6_error = self._runValidationAndAnalysis()
        if phase6_error and phase6_error.maxSeverity() >= 4:
            return phase6_error

        return error

    # =========================================================================
    # Phase 1: ProSMART protein restraints
    # =========================================================================

    def _runProsmartProtein(self):
        """Run ProSMART to generate protein external restraints (optional)."""
        error = CErrorReport()

        if not self.container.prosmartProtein.TOGGLE:
            return error

        try:
            self.prosmartProteinPlugin = self.makePluginObject('prosmart')
        except Exception as e:
            self.appendErrorReport(109,
                f'Failed to create prosmart plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 109,
                        f'Failed to create prosmart: {e}', 'prosmart_protein', 4)
            return error

        try:
            plugin = self.prosmartProteinPlugin
            plugin.container.inputData.TARGET_MODEL = self.container.inputData.XYZIN
            plugin.container.inputData.CHAINLIST_1 = self.container.prosmartProtein.CHAINLIST_1
            plugin.container.inputData.REFERENCE_MODELS = self.container.prosmartProtein.REFERENCE_MODELS
            plugin.container.controlParameters.RESTRAIN_ALL_VS_BEST = self.container.prosmartProtein.ALL_BEST
            plugin.container.controlParameters.RESTRAIN_SEQID = self.container.prosmartProtein.SEQID
            plugin.container.controlParameters.RESTRAIN_MAIN_VS_SIDE = self.container.prosmartProtein.SIDE_MAIN
            plugin.container.controlParameters.RESTRAIN_RMIN = self.container.prosmartProtein.RMIN
            plugin.container.controlParameters.RESTRAIN_RMAX = self.container.prosmartProtein.RMAX

            if self.container.prosmartProtein.ADVANCED:
                plugin.container.controlParameters.RESTRAIN_BFAC_FILTER = self.container.prosmartProtein.TOGGLE_BFAC
                plugin.container.controlParameters.RESTRAIN_BFAC_ALPHA = self.container.prosmartProtein.BFAC
                plugin.container.controlParameters.RESTRAIN_ALT = self.container.prosmartProtein.TOGGLE_ALT
                plugin.container.controlParameters.RESTRAIN_OCCUP = self.container.prosmartProtein.OCCUPANCY
                plugin.container.controlParameters.KEYWORDS = self.container.prosmartProtein.KEYWORDS

            print("[servalcat_pipe] Running ProSMART protein restraints...")
            status = plugin.process()

            if status == CPluginScript.FAILED:
                self.appendErrorReport(102, 'ProSMART protein restraints failed')
                error.append(self.__class__.__name__, 102,
                            'ProSMART protein restraints failed', 'prosmart_protein', 4)
                return error

            print("[servalcat_pipe] ProSMART protein restraints completed")

        except Exception as e:
            self.appendErrorReport(102,
                f'Exception in ProSMART protein: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 102,
                        f'Exception in ProSMART protein: {e}', 'prosmart_protein', 4)
            return error

        return error

    # =========================================================================
    # Phase 2: ProSMART nucleic acid restraints
    # =========================================================================

    def _runProsmartNucleicAcid(self):
        """Run ProSMART to generate nucleic acid external restraints (optional)."""
        error = CErrorReport()

        if not self.container.prosmartNucleicAcid.TOGGLE:
            return error

        try:
            self.prosmartNucleicAcidPlugin = self.makePluginObject('prosmart')
        except Exception as e:
            self.appendErrorReport(109,
                f'Failed to create prosmart plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 109,
                        f'Failed to create prosmart: {e}', 'prosmart_na', 4)
            return error

        try:
            plugin = self.prosmartNucleicAcidPlugin
            plugin.container.controlParameters.NUCLEIC_ACID = True
            plugin.container.inputData.TARGET_MODEL = self.container.inputData.XYZIN
            plugin.container.inputData.CHAINLIST_1 = self.container.prosmartNucleicAcid.CHAINLIST_1
            plugin.container.inputData.REFERENCE_MODELS = self.container.prosmartNucleicAcid.REFERENCE_MODELS
            plugin.container.controlParameters.RESTRAIN_ALL_VS_BEST = self.container.prosmartNucleicAcid.ALL_BEST
            plugin.container.controlParameters.RESTRAIN_SEQID = self.container.prosmartNucleicAcid.SEQID
            plugin.container.controlParameters.RESTRAIN_MAIN_VS_SIDE = self.container.prosmartNucleicAcid.SIDE_MAIN
            plugin.container.controlParameters.RESTRAIN_RMIN = self.container.prosmartNucleicAcid.RMIN
            plugin.container.controlParameters.RESTRAIN_RMAX = self.container.prosmartNucleicAcid.RMAX

            if self.container.prosmartNucleicAcid.ADVANCED:
                plugin.container.controlParameters.RESTRAIN_BFAC_FILTER = self.container.prosmartNucleicAcid.TOGGLE_BFAC
                plugin.container.controlParameters.RESTRAIN_BFAC_ALPHA = self.container.prosmartNucleicAcid.BFAC
                plugin.container.controlParameters.RESTRAIN_ALT = self.container.prosmartNucleicAcid.TOGGLE_ALT
                plugin.container.controlParameters.RESTRAIN_OCCUP = self.container.prosmartNucleicAcid.OCCUPANCY
                plugin.container.controlParameters.KEYWORDS = self.container.prosmartNucleicAcid.KEYWORDS

            print("[servalcat_pipe] Running ProSMART nucleic acid restraints...")
            status = plugin.process()

            if status == CPluginScript.FAILED:
                self.appendErrorReport(103, 'ProSMART nucleic acid restraints failed')
                error.append(self.__class__.__name__, 103,
                            'ProSMART nucleic acid restraints failed', 'prosmart_na', 4)
                return error

            print("[servalcat_pipe] ProSMART nucleic acid restraints completed")

        except Exception as e:
            self.appendErrorReport(103,
                f'Exception in ProSMART nucleic acid: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 103,
                        f'Exception in ProSMART nucleic acid: {e}', 'prosmart_na', 4)
            return error

        return error

    # =========================================================================
    # Phase 3: MetalCoord restraints
    # =========================================================================

    def _runMetalCoords(self):
        """Run MetalCoord to generate metal site restraints (optional)."""
        error = CErrorReport()

        if not self.container.metalCoordPipeline.RUN_METALCOORD:
            return error

        if self.container.metalCoordPipeline.GENERATE_OR_USE != "GENERATE":
            return error

        if not shutil.which("metalCoord", mode=os.X_OK):
            print("[servalcat_pipe] WARNING: MetalCoord not installed, skipping")
            return error

        # Run MetalCoord for each selected ligand code
        metalCoordOutputJsonPaths = []
        ligand_codes_selected = self.container.metalCoordPipeline.LIGAND_CODES_SELECTED

        for ligand_code in ligand_codes_selected:
            try:
                metalCoordPlugin = self.makePluginObject('metalCoord')
                self.metalCoordPlugins.append(metalCoordPlugin)
            except Exception as e:
                self.appendErrorReport(109,
                    f'Failed to create metalCoord plugin for {ligand_code}: {e}\n{traceback.format_exc()}')
                error.append(self.__class__.__name__, 109,
                            f'Failed to create metalCoord: {e}', 'metalCoord', 4)
                return error

            try:
                metalCoordPlugin.container.inputData.XYZIN.set(self.container.inputData.XYZIN)
                metalCoordPlugin.container.controlParameters.copyData(
                    self.container.metalCoordWrapper.controlParameters)
                metalCoordPlugin.container.inputData.LIGAND_CODE.set(ligand_code)
                metalCoordPlugin.container.controlParameters.SAVE_PDBMMCIF.set(False)

                if str(self.container.metalCoordPipeline.LINKS) == "KEEP":
                    metalCoordPlugin.container.controlParameters.KEEP_LINKS.set(True)

                print(f"[servalcat_pipe] Running MetalCoord for ligand {ligand_code}...")
                status = metalCoordPlugin.process()

                if status == CPluginScript.FAILED:
                    self.appendErrorReport(104,
                        f'MetalCoord failed for ligand {ligand_code}')
                    # MetalCoord failure for one ligand is non-fatal; continue
                    continue

                # Check for output JSON
                outputJsonFilename = str(ligand_code) + ".json"
                outputJsonPath = os.path.join(
                    metalCoordPlugin.getWorkDirectory(), outputJsonFilename)
                if os.path.isfile(outputJsonPath):
                    metalCoordOutputJsonPaths.append(outputJsonPath)

            except Exception as e:
                self.appendErrorReport(104,
                    f'Exception in MetalCoord for {ligand_code}: {e}\n{traceback.format_exc()}')
                continue

        if not metalCoordOutputJsonPaths:
            print("[servalcat_pipe] WARNING: No MetalCoord output produced")
            return error

        # Convert JSON outputs to restraint keywords
        try:
            outputRestraintsPrefix = "metal_restraints"
            outputRestraintsPath = os.path.join(
                self.getWorkDirectory(), outputRestraintsPrefix + ".txt")
            outputRestraintsMmcifPath = os.path.join(
                self.getWorkDirectory(), outputRestraintsPrefix + ".mmcif")
            outputRestraintsPathPrefix = os.path.join(
                self.getWorkDirectory(), outputRestraintsPrefix)

            if self.container.metalCoordPipeline.LINKS == "UPDATE":
                stPath = str(self.container.inputData.XYZIN.fullPath)
                keep_links = False
            elif self.container.metalCoordPipeline.LINKS == "KEEP":
                stPath = str(self.container.inputData.XYZIN.fullPath)
                keep_links = True
            else:  # NOTTOUCH
                stPath = None
                keep_links = True

            print("[servalcat_pipe] Converting MetalCoord JSON to restraints...")
            from ccp4i2.wrappers.metalCoord.script import json2restraints
            json2restraints.main(
                jsonPaths=metalCoordOutputJsonPaths,
                stPath=stPath,
                outputPrefix=outputRestraintsPathPrefix,
                jsonEquivalentsPath=None,
                keep_links=keep_links)

            if os.path.isfile(outputRestraintsPath):
                self.container.outputData.METALCOORD_RESTRAINTS.setFullPath(
                    outputRestraintsPath)
                self.container.outputData.METALCOORD_RESTRAINTS.annotation = \
                    'Restraints for metal sites'

            if os.path.isfile(outputRestraintsMmcifPath) and stPath:
                self.container.outputData.METALCOORD_XYZ.setFullPath(
                    outputRestraintsMmcifPath)
                self.container.outputData.METALCOORD_XYZ.annotation = \
                    'Structure model with links from MetalCoord (mmCIF format)'

            print("[servalcat_pipe] MetalCoord restraints completed")

        except Exception as e:
            self.appendErrorReport(111,
                f'Failed to convert MetalCoord JSON to restraints: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 111,
                        f'MetalCoord JSON conversion failed: {e}', 'metalCoord', 4)
            return error

        return error

    # =========================================================================
    # Phase 4: Servalcat refinement
    # =========================================================================

    def _runServalcat(self):
        """Run the main servalcat refinement."""
        error = CErrorReport()

        try:
            self.servalcatPlugin = self._createServalcatJob()
        except Exception as e:
            self.appendErrorReport(109,
                f'Failed to create servalcat wrapper: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 109,
                        f'Failed to create servalcat: {e}', 'servalcat', 4)
            return error

        try:
            plugin = self.servalcatPlugin
            plugin.doAsync = False

            # Monitor wrapper's program.xml for real-time cycle-by-cycle progress
            stopMonitor = self._monitorSubPluginXml(plugin, 'SERVALCAT_FIRST')

            print("[servalcat_pipe] Running servalcat refinement...")
            status = plugin.process()

            # Stop monitoring before final sync
            stopMonitor()

            if status == CPluginScript.FAILED:
                self.appendErrorReport(105, 'Servalcat refinement failed')
                error.append(self.__class__.__name__, 105,
                            'Servalcat refinement failed', 'servalcat', 4)
                return error

            if status == CPluginScript.UNSATISFACTORY:
                print("[servalcat_pipe] Servalcat completed with UNSATISFACTORY status")
                # Still continue - unsatisfactory results are still usable

            # Verify critical output exists
            if not os.path.isfile(str(plugin.container.outputData.CIFFILE.fullPath)):
                self.appendErrorReport(106, 'Servalcat did not produce output coordinates')
                error.append(self.__class__.__name__, 106,
                            'Missing CIFFILE output', 'servalcat', 4)
                return error

            # Final sync (gets the complete XML including processOutputFiles data)
            self._appendPluginXml(plugin, tag='SERVALCAT_FIRST')
            print("[servalcat_pipe] Servalcat refinement completed")

        except Exception as e:
            self.appendErrorReport(105,
                f'Exception in servalcat: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 105,
                        f'Exception in servalcat: {e}', 'servalcat', 4)
            return error

        return error

    # =========================================================================
    # Phase 5: Water addition + re-refinement
    # =========================================================================

    def _runWaterAddition(self):
        """Optionally add waters with coot, then re-refine (optional)."""
        error = CErrorReport()

        if not self.container.controlParameters.ADD_WATERS:
            return error

        # --- 5a: Find waters with coot ---
        try:
            self.cootPlugin = self.makePluginObject('coot_find_waters')
        except Exception as e:
            self.appendErrorReport(109,
                f'Failed to create coot_find_waters: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 109,
                        f'Failed to create coot_find_waters: {e}', 'coot', 4)
            return error

        try:
            self.cootPlugin.container.inputData.XYZIN.set(
                self.servalcatPlugin.container.outputData.CIFFILE)
            self.cootPlugin.container.inputData.FPHIIN.set(
                self.servalcatPlugin.container.outputData.FPHIOUT)

            print("[servalcat_pipe] Running coot find waters...")
            status = self.cootPlugin.process()

            if status == CPluginScript.FAILED:
                self.appendErrorReport(107, 'Coot find waters failed')
                error.append(self.__class__.__name__, 107,
                            'Coot find waters failed', 'coot', 4)
                return error

            if not self.cootPlugin.container.outputData.XYZOUT.exists():
                self.appendErrorReport(112, 'Coot did not produce output coordinates')
                error.append(self.__class__.__name__, 112,
                            'Missing coot XYZOUT', 'coot', 4)
                return error

            # Record water count in XML
            self._recordWaterCount()

            print("[servalcat_pipe] Coot find waters completed")

        except Exception as e:
            self.appendErrorReport(107,
                f'Exception in coot find waters: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 107,
                        f'Exception in coot: {e}', 'coot', 4)
            return error

        # --- 5b: Re-refine with servalcat after adding waters ---
        try:
            self.servalcatPostCootPlugin = self._createServalcatJob(
                inputCoordinates=self.cootPlugin.container.outputData.XYZOUT,
                ncyc=int(self.container.controlParameters.NCYCLES_AFTER_ADD_WATERS))
        except Exception as e:
            self.appendErrorReport(109,
                f'Failed to create post-coot servalcat: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 109,
                        f'Failed to create post-coot servalcat: {e}', 'servalcat_post', 4)
            return error

        try:
            self.servalcatPostCootPlugin.doAsync = False

            # Monitor for real-time progress during post-coot refinement
            stopMonitor = self._monitorSubPluginXml(
                self.servalcatPostCootPlugin, 'SERVALCAT_WATERS')

            print("[servalcat_pipe] Running post-coot servalcat refinement...")
            status = self.servalcatPostCootPlugin.process()

            stopMonitor()

            if status == CPluginScript.FAILED:
                self.appendErrorReport(108, 'Post-coot servalcat refinement failed')
                error.append(self.__class__.__name__, 108,
                            'Post-coot servalcat failed', 'servalcat_post', 4)
                return error

            self._appendPluginXml(self.servalcatPostCootPlugin, tag='SERVALCAT_WATERS')
            print("[servalcat_pipe] Post-coot servalcat refinement completed")

        except Exception as e:
            self.appendErrorReport(108,
                f'Exception in post-coot servalcat: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 108,
                        f'Exception in post-coot servalcat: {e}', 'servalcat_post', 4)
            return error

        return error

    # =========================================================================
    # Phase 6: Validation and analysis
    # =========================================================================

    def _runValidationAndAnalysis(self):
        """Run optional validation and analysis steps (non-fatal).

        Note: This runs inside startProcess(), BEFORE processOutputFiles() copies
        files to the pipeline directory. All file references must use the wrapper's
        output paths (which exist on disk), not the pipeline's output paths (which
        don't exist yet).
        """
        error = CErrorReport()

        # Determine which servalcat job produced the final output
        finalServalcatPlugin = self.servalcatPostCootPlugin or self.servalcatPlugin

        # Path to the refined model (in the wrapper's job directory, not pipeline's)
        refinedModelPath = str(finalServalcatPlugin.container.outputData.CIFFILE.fullPath)

        # Validation (non-fatal - errors logged but don't fail pipeline)
        try:
            self._runMultimericValidation(finalServalcatPlugin)
        except Exception as e:
            self.appendErrorReport(110,
                f'Validation failed (non-fatal): {e}\n{traceback.format_exc()}')

        # ADP analysis (non-fatal)
        if self.container.controlParameters.RUN_ADP_ANALYSIS:
            try:
                self._runAdpAnalysis(
                    refinedModelPath,
                    float(self.container.controlParameters.ADP_IQR_FACTOR))
            except Exception as e:
                self.appendErrorReport(113,
                    f'ADP analysis failed (non-fatal): {e}\n{traceback.format_exc()}')

        # Coordinate/ADP deviation monitoring (non-fatal)
        if self.container.monitor.RUN_COORDADPDEV_ANALYSIS:
            try:
                self._runCoordAdpDevAnalysis(
                    str(self.container.inputData.XYZIN.fullPath),
                    refinedModelPath)
            except Exception as e:
                self.appendErrorReport(114,
                    f'Coord/ADP deviation analysis failed (non-fatal): {e}\n{traceback.format_exc()}')

        return error

    # =========================================================================
    # Helper: Create servalcat wrapper job
    # =========================================================================

    def _createServalcatJob(self, withWeight=-1, inputCoordinates=None, ncyc=-1):
        """Create and configure a servalcat wrapper job."""
        result = self.makePluginObject('servalcat')

        # Copy input data from pipeline
        result.container.inputData.copyData(self.container.inputData)

        # Copy control parameters, excluding pipeline-only params
        EXCLUDED_PARAMS = {
            'OPTIMISE_WEIGHT', 'REFMAC_CLEANUP',
            'VALIDATE_IRIS', 'VALIDATE_BAVERAGE',
            'VALIDATE_RAMACHANDRAN', 'VALIDATE_MOLPROBITY',
            'RUN_MOLPROBITY',
        }
        for attr in self.container.controlParameters.dataOrder():
            if attr not in EXCLUDED_PARAMS:
                setattr(result.container.controlParameters, attr,
                        getattr(self.container.controlParameters, attr))

        # MetalCoord restraints
        if self.container.metalCoordPipeline.RUN_METALCOORD:
            if self.container.metalCoordPipeline.GENERATE_OR_USE == "GENERATE":
                result.container.inputData.METALCOORD_RESTRAINTS = \
                    self.container.outputData.METALCOORD_RESTRAINTS
                if (self.container.outputData.METALCOORD_XYZ and
                        self.container.metalCoordPipeline.LINKS != "NOTTOUCH"):
                    if os.path.isfile(str(self.container.outputData.METALCOORD_XYZ.fullPath)):
                        result.container.inputData.XYZIN.set(
                            self.container.outputData.METALCOORD_XYZ)
            else:
                result.container.inputData.METALCOORD_RESTRAINTS = \
                    self.container.metalCoordPipeline.METALCOORD_RESTRAINTS

        # ProSMART protein restraints
        if self.container.prosmartProtein.TOGGLE and self.prosmartProteinPlugin is not None:
            result.container.controlParameters.PROSMART_PROTEIN_SGMN = \
                self.container.prosmartProtein.SGMN
            result.container.controlParameters.PROSMART_PROTEIN_SGMX = \
                self.container.prosmartProtein.SGMX
            result.container.controlParameters.PROSMART_PROTEIN_ALPHA = \
                self.container.prosmartProtein.ALPHA
            result.container.controlParameters.PROSMART_PROTEIN_DMAX = \
                self.container.prosmartProtein.DMAX
            result.container.inputData.PROSMART_PROTEIN_RESTRAINTS = \
                self.prosmartProteinPlugin.container.outputData.RESTRAINTS

        # ProSMART nucleic acid restraints
        if self.container.prosmartNucleicAcid.TOGGLE and self.prosmartNucleicAcidPlugin is not None:
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMN = \
                self.container.prosmartNucleicAcid.SGMN
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMX = \
                self.container.prosmartNucleicAcid.SGMX
            result.container.controlParameters.PROSMART_NUCLEICACID_ALPHA = \
                self.container.prosmartNucleicAcid.ALPHA
            result.container.controlParameters.PROSMART_NUCLEICACID_DMAX = \
                self.container.prosmartNucleicAcid.DMAX
            result.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS = \
                self.prosmartNucleicAcidPlugin.container.outputData.RESTRAINTS

        # Manual weight override
        if withWeight >= 0.:
            result.container.controlParameters.WEIGHT_OPT = 'MANUAL'
            result.container.controlParameters.WEIGHT = withWeight

        # Override coordinates and cycle count (for post-coot refinement)
        if inputCoordinates is not None:
            result.container.inputData.XYZIN.set(inputCoordinates)
            if ncyc > 0:
                result.container.controlParameters.NCYCLES.set(ncyc)

        return result

    # =========================================================================
    # Helper: Record water count in XML
    # =========================================================================

    def _recordWaterCount(self):
        """Extract water count from coot output and record in XML."""
        try:
            cootLogXml = os.path.join(
                os.path.dirname(str(self.cootPlugin.container.outputData.XYZOUT)),
                "program.xml")
            nwaters = "unknown"
            if os.path.isfile(cootLogXml):
                with open(cootLogXml, encoding='utf-8') as f:
                    watersXml = etree.fromstring(f.read())
                    nodes = watersXml.findall(".//WatersFound")
                    if len(nodes) > 0:
                        nwaters = nodes[0].text
            cootNode = etree.SubElement(self.xmlroot, "CootAddWaters")
            cootNode.text = f"Coot added {nwaters} water molecules."
            self._flushXML()
        except Exception as e:
            self.appendErrorReport(107,
                f'Failed to record water count: {e}')

    # =========================================================================
    # Validation methods
    # =========================================================================

    def _runMultimericValidation(self, finalServalcatPlugin):
        """Run geometry validation (Iris, B-factors, Ramachandran, MolProbity).

        Args:
            finalServalcatPlugin: The servalcat wrapper that produced the refined model.
                Uses its output CIFFILE (which exists on disk) rather than the pipeline's
                output CIFFILE (which hasn't been copied yet).
        """
        validate_iris = getattr(self.container.controlParameters, 'VALIDATE_IRIS', False)
        validate_baverage = getattr(self.container.controlParameters, 'VALIDATE_BAVERAGE', False)
        validate_ramachandran = getattr(self.container.controlParameters, 'VALIDATE_RAMACHANDRAN', False)
        validate_molprobity = getattr(self.container.controlParameters, 'VALIDATE_MOLPROBITY', False)

        if not any([validate_iris, validate_baverage, validate_molprobity, validate_ramachandran]):
            return

        self.validatePlugin = self.makePluginObject('validate_protein')
        self.validatePlugin.container.inputData.XYZIN_2.set(
            finalServalcatPlugin.container.outputData.CIFFILE)
        self.validatePlugin.container.inputData.XYZIN_1.set(
            self.container.inputData.XYZIN)
        self.validatePlugin.container.inputData.NAME_2.set("Refined")
        self.validatePlugin.container.inputData.NAME_1.set("Input")

        if (str(self.container.controlParameters.SCATTERING_FACTORS) == "XRAY"
                and str(self.container.controlParameters.MERGED_OR_UNMERGED) == "merged"):
            self.validatePlugin.container.inputData.F_SIGF_2.set(
                self.container.inputData.HKLIN)
            self.validatePlugin.container.inputData.F_SIGF_1.set(
                self.container.inputData.HKLIN)
        else:
            self.validatePlugin.container.inputData.F_SIGF_2.set(None)
            self.validatePlugin.container.inputData.F_SIGF_1.set(None)

        self.validatePlugin.container.controlParameters.TWO_DATASETS.set(True)
        self.validatePlugin.container.controlParameters.DO_IRIS.set(validate_iris)
        self.validatePlugin.container.controlParameters.DO_BFACT.set(validate_baverage)
        self.validatePlugin.container.controlParameters.DO_RAMA.set(validate_ramachandran)
        self.validatePlugin.container.controlParameters.DO_MOLPROBITY.set(validate_molprobity)

        print("[servalcat_pipe] Running validation...")
        self.validatePlugin.process()

        # Extract validation XML into pipeline XML
        try:
            validateXMLPath = self.validatePlugin.makeFileName('PROGRAMXML')
            validateXML = CCP4Utils.openFileToEtree(validateXMLPath)
            xml_validation = etree.SubElement(self.xmlroot, "Validation")

            sections = [
                (True, "//Validate_geometry_CCP4i2/Model_info"),
                (validate_iris, "//Validate_geometry_CCP4i2/Iris"),
                (validate_baverage, "//Validate_geometry_CCP4i2/B_factors"),
                (validate_baverage, "//Validate_geometry_CCP4i2/B_averages"),
                (validate_ramachandran, "//Validate_geometry_CCP4i2/Ramachandran"),
                (validate_molprobity, "//Validate_geometry_CCP4i2/Molprobity"),
            ]
            for enabled, xpath in sections:
                if enabled:
                    nodes = validateXML.xpath(xpath)
                    if len(nodes) > 0:
                        xml_validation.append(nodes[0])

            self._flushXML()
            print("[servalcat_pipe] Validation completed")

        except Exception as e:
            self.appendErrorReport(110,
                f'Failed to extract validation XML: {e}')

    def _runAdpAnalysis(self, modelPath, iqrFactor=2.0):
        """Run ADP analysis on the refined model."""
        import csv

        print("[servalcat_pipe] Running ADP analysis...")

        adp_dict = {}
        adp_per_resi = {}
        adp_dict["All"] = []
        cif_block = gemmi.cif.read(modelPath)[0]
        st = gemmi.make_structure_from_block(cif_block)

        for model in st:
            for chain in model:
                polymer = chain.get_polymer()
                ptype = polymer.check_polymer_type()
                adp_dict[chain.name] = []
                adp_per_resi[chain.name] = {"resi": [], "adp": [], "adp_sidechain": []}
                for residue in chain:
                    adp_this_resi = []
                    adp_this_resi_sidechain = []
                    for atom in residue:
                        if not atom.is_hydrogen() and atom.occ > 0:
                            if atom.aniso.nonzero():
                                adp_atom = gemmi.calculate_b_est(atom)
                            else:
                                adp_atom = atom.b_iso
                            adp_dict["All"].append(adp_atom)
                            adp_dict[chain.name].append(adp_atom)
                            if (residue.entity_type == gemmi.EntityType.Polymer and
                                    ptype in [gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD] and
                                    atom.name not in ["CA", "C", "O", "N", "OXT"]):
                                adp_this_resi_sidechain.append(adp_atom)
                            else:
                                adp_this_resi.append(adp_atom)
                    try:
                        if residue.seqid.num in adp_per_resi[chain.name]["resi"]:
                            continue  # ignoring insertion codes
                        adp_per_resi[chain.name]["resi"].append(residue.seqid.num)
                        if adp_this_resi:
                            adp_per_resi[chain.name]["adp"].append(numpy.mean(adp_this_resi))
                        else:
                            adp_per_resi[chain.name]["adp"].append(None)
                        if adp_this_resi_sidechain:
                            adp_per_resi[chain.name]["adp_sidechain"].append(
                                numpy.mean(adp_this_resi_sidechain))
                        else:
                            adp_per_resi[chain.name]["adp_sidechain"].append(None)
                    except Exception:
                        pass

        if not adp_dict["All"]:
            print("[servalcat_pipe] ADP analysis: no atoms found")
            return

        # Find ADP outliers
        q1 = numpy.quantile(adp_dict["All"], 0.25)
        q3 = numpy.quantile(adp_dict["All"], 0.75)
        iqr = q3 - q1
        adp_limit_low = q1 - iqrFactor * iqr
        adp_limit_high = q3 + iqrFactor * iqr
        adp_low = []
        adp_high = []
        for model in st:
            for cra in model.all():
                if not cra.atom.is_hydrogen and cra.atom.occ > 0:
                    if cra.atom.aniso.nonzero():
                        adp_atom = gemmi.calculate_b_est(cra.atom)
                    else:
                        adp_atom = cra.atom.b_iso
                    if adp_atom < adp_limit_low:
                        adp_low.append({"atom": str(cra), "adp": adp_atom})
                    elif adp_atom > adp_limit_high:
                        adp_high.append({"atom": str(cra), "adp": adp_atom})
        adp_low = sorted(adp_low, key=itemgetter('adp'))
        adp_high = sorted(adp_high, key=itemgetter('adp'), reverse=True)

        # Write CSV and XML
        csvFileName = "adp_analysis.csv"
        csvFilePath = str(os.path.join(self.getWorkDirectory(), csvFileName))
        with open(csvFilePath, "a+", newline='') as csvfile:
            fieldnames = ["chain", "resi", "adp", "adp_sidechain"]
            writer = csv.writer(csvfile, quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow(fieldnames)

        adp_root = etree.Element('ADP_ANALYSIS')
        chains_root = etree.SubElement(adp_root, "chains")
        for ch, values in adp_dict.items():
            chain = etree.SubElement(chains_root, "chain")
            chain_name = etree.SubElement(chain, "name")
            chain_name.text = str(ch)
            chain.set("name", str(ch))
            chain_min = etree.SubElement(chain, "min")
            chain_min.text = "{:.2f}".format(min(values))
            chain_max = etree.SubElement(chain, "max")
            chain_max.text = "{:.2f}".format(max(values))
            chain_med_val = numpy.median(values)
            chain_mad_val = numpy.median(numpy.absolute(values - chain_med_val))
            chain_med = etree.SubElement(chain, "med")
            chain_med.text = "{:.2f}".format(chain_med_val)
            chain_mad = etree.SubElement(chain, "mad")
            chain_mad.text = "{:.2f}".format(chain_mad_val)
            chain_q1 = etree.SubElement(chain, "q1")
            chain_q1.text = "{:.2f}".format(numpy.quantile(values, 0.25))
            chain_q3 = etree.SubElement(chain, "q3")
            chain_q3.text = "{:.2f}".format(numpy.quantile(values, 0.75))
            chain_mean_val = numpy.mean(values)
            chain_std_val = numpy.std(values)
            chain_mean = etree.SubElement(chain, "mean")
            chain_mean.text = "{:.2f}".format(chain_mean_val)
            chain_std = etree.SubElement(chain, "std")
            chain_std.text = "{:.2f}".format(chain_std_val)

            hist, bin_edges = numpy.histogram(values, bins="auto")
            bin_edges = numpy.delete(bin_edges, -1)
            chain_histogram = etree.SubElement(chain, 'histogram')
            for i in range(len(bin_edges)):
                bin_elem = etree.SubElement(chain_histogram, 'bin')
                bin_adp = etree.SubElement(bin_elem, 'adp')
                bin_adp.text = "{:.2f}".format(bin_edges[i])
                bin_count = etree.SubElement(bin_elem, 'count')
                bin_count.text = str(hist[i])

            if ch != "All":
                with open(csvFilePath, "a+", newline='') as csvfile:
                    writer = csv.writer(csvfile, quoting=csv.QUOTE_NONNUMERIC)
                    for i in range(len(adp_per_resi[str(ch)]["resi"])):
                        adp_val = adp_per_resi[str(ch)]["adp"][i]
                        adp_sc_val = adp_per_resi[str(ch)]["adp_sidechain"][i]
                        adp_round = round(adp_val, 2) if adp_val else "-"
                        adp_sc_round = round(adp_sc_val, 2) if adp_sc_val else "-"
                        writer.writerow([ch,
                                        adp_per_resi[str(ch)]["resi"][i],
                                        adp_round,
                                        adp_sc_round])

        if os.path.exists(csvFilePath):
            per_resi_element = etree.SubElement(adp_root, "per_resi")
            per_resi_csv_element = etree.SubElement(per_resi_element, "CSV_FILE")
            per_resi_csv_element.text = str(csvFileName)

        outliers_root = etree.SubElement(adp_root, "outliers")
        etree.SubElement(outliers_root, "adp_limit_low").text = "{:.2f}".format(adp_limit_low)
        etree.SubElement(outliers_root, "adp_limit_high").text = "{:.2f}".format(adp_limit_high)
        etree.SubElement(outliers_root, "iqr_factor").text = "{:.2f}".format(iqrFactor)

        outliers_low = etree.SubElement(outliers_root, "low")
        for outlier in adp_low:
            outlier_elem = etree.SubElement(outliers_low, "data")
            etree.SubElement(outlier_elem, "adp").text = "{:.2f}".format(outlier["adp"])
            etree.SubElement(outlier_elem, "atom").text = str(outlier["atom"])

        outliers_high = etree.SubElement(outliers_root, "high")
        for outlier in adp_high:
            outlier_elem = etree.SubElement(outliers_high, "data")
            etree.SubElement(outlier_elem, "adp").text = "{:.2f}".format(outlier["adp"])
            etree.SubElement(outlier_elem, "atom").text = str(outlier["atom"])

        self.xmlroot.append(adp_root)
        self._flushXML()
        print("[servalcat_pipe] ADP analysis completed")

    def _runCoordAdpDevAnalysis(self, model1Path, model2Path):
        """Monitor changes/shifts of coordinates and ADPs."""
        print("[servalcat_pipe] Running coordinate/ADP deviation analysis...")

        coordDevMinReported = float(self.container.monitor.MIN_COORDDEV)
        ADPAbsDevMinReported = float(self.container.monitor.MIN_ADPDEV)
        csvFileName = "report_coord_adp_dev.csv"
        csvFilePath = str(os.path.join(self.getWorkDirectory(), csvFileName))

        df = monitor_differences.main(
            file1=model1Path, file2=model2Path, output=csvFilePath,
            minCoordDev=coordDevMinReported, minAdpDev=ADPAbsDevMinReported,
            useHydrogens=False)

        if df is None or df.empty or not all(col in df.columns for col in ["CoordDev", "ADPDev"]):
            print(f"[servalcat_pipe] No significant coordinate/ADP deviations found "
                  f"above thresholds ({coordDevMinReported} A, {ADPAbsDevMinReported} A^2)")
            return

        coordDevMean = df["CoordDev"].mean()
        ADPAbsDevMean = df["ADPDev"].mean()

        xmlText = "\n<COORD_ADP_DEV>"
        xmlText += "\n<STATISTICS>"
        xmlText += f"\n<coordDevMean>{round(coordDevMean, 2)}</coordDevMean>"
        xmlText += f"\n<coordDevMinReported>{round(coordDevMinReported, 2)}</coordDevMinReported>"
        xmlText += f"\n<ADPAbsDevMean>{round(ADPAbsDevMean, 2)}</ADPAbsDevMean>"
        xmlText += f"\n<coordADPAbsMinReported>{round(ADPAbsDevMinReported, 2)}</coordADPAbsMinReported>"
        xmlText += "\n</STATISTICS>"
        if os.path.exists(csvFilePath):
            xmlText += f"\n<CSV_FILE>{csvFileName}</CSV_FILE>"
        xmlText += "\n</COORD_ADP_DEV>"

        xmlTree = etree.fromstring(xmlText)
        self.xmlroot.append(xmlTree)
        self._flushXML()
        print("[servalcat_pipe] Coordinate/ADP deviation analysis completed")

    # =========================================================================
    # Output harvesting
    # =========================================================================

    def processOutputFiles(self):
        """Harvest output files from the final servalcat job."""
        error = CErrorReport()

        # Determine which servalcat job produced the final output
        finalServalcatPlugin = self.servalcatPostCootPlugin or self.servalcatPlugin

        if finalServalcatPlugin is None:
            self.appendErrorReport(105, 'No servalcat output available')
            error.append(self.__class__.__name__, 105,
                        'No servalcat output', 'processOutputFiles', 4)
            return error

        try:
            # Copy output files from the final servalcat job
            for attr in self.container.outputData.dataOrder():
                try:
                    wrappersAttr = getattr(finalServalcatPlugin.container.outputData, attr)
                    pipelinesAttr = getattr(self.container.outputData, attr)
                except Exception:
                    continue

                if attr == "PERFORMANCEINDICATOR":
                    setattr(self.container.outputData, attr, wrappersAttr)
                elif os.path.exists(str(wrappersAttr.fullPath)):
                    try:
                        shutil.copyfile(str(wrappersAttr.fullPath),
                                       str(pipelinesAttr.fullPath))
                    except Exception as e:
                        self.appendErrorReport(101,
                            f'{wrappersAttr.fullPath} to {pipelinesAttr.fullPath}: {e}')

            # Apply database annotations
            self._applyAnnotations(finalServalcatPlugin)

            # Optional cleanup of intermediate files
            self._cleanupIntermediateFiles()

            # Flush final XML
            self._flushXML()

        except Exception as e:
            self.appendErrorReport(101,
                f'Exception in processOutputFiles: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 101,
                        f'Exception in processOutputFiles: {e}', 'harvest', 3)

        return error

    def _applyAnnotations(self, servalcatJob):
        """Copy annotations and metadata from servalcat output to pipeline output."""
        try:
            self.container.outputData.XYZOUT.annotation.set(
                servalcatJob.container.outputData.XYZOUT.annotation)
            self.container.outputData.CIFFILE.annotation.set(
                servalcatJob.container.outputData.CIFFILE.annotation)
            self.container.outputData.FPHIOUT.annotation.set(
                servalcatJob.container.outputData.FPHIOUT.annotation)
            self.container.outputData.FPHIOUT.subType = \
                servalcatJob.container.outputData.FPHIOUT.subType
            self.container.outputData.DIFFPHIOUT.annotation.set(
                servalcatJob.container.outputData.DIFFPHIOUT.annotation)
            self.container.outputData.DIFFPHIOUT.subType = \
                servalcatJob.container.outputData.DIFFPHIOUT.subType

            if self.container.outputData.DICT.exists():
                self.container.outputData.DICT.annotation = 'Accumulated monomer dictionary'

            if servalcatJob.container.outputData.COOTSCRIPTOUT.exists():
                if servalcatJob.container.outputData.COOTSCRIPTOUT.annotation.isSet():
                    self.container.outputData.COOTSCRIPTOUT.annotation.set(
                        servalcatJob.container.outputData.COOTSCRIPTOUT.annotation)

            if str(servalcatJob.container.controlParameters.DATA_METHOD) == 'xtal':
                if servalcatJob.container.outputData.ANOMFPHIOUT.exists():
                    if servalcatJob.container.outputData.ANOMFPHIOUT.annotation.isSet():
                        self.container.outputData.ANOMFPHIOUT.annotation.set(
                            servalcatJob.container.outputData.ANOMFPHIOUT.annotation)
                if servalcatJob.container.outputData.DIFANOMFPHIOUT.exists():
                    if servalcatJob.container.outputData.DIFANOMFPHIOUT.annotation.isSet():
                        self.container.outputData.DIFANOMFPHIOUT.annotation.set(
                            servalcatJob.container.outputData.DIFANOMFPHIOUT.annotation)
            elif str(servalcatJob.container.controlParameters.DATA_METHOD) == 'spa':
                self.container.outputData.MAP_FO.annotation.set(
                    servalcatJob.container.outputData.MAP_FO.annotation)
                self.container.outputData.MAP_FO.subType = \
                    servalcatJob.container.outputData.MAP_FO.subType
                self.container.outputData.MAP_FOFC.annotation.set(
                    servalcatJob.container.outputData.MAP_FOFC.annotation)
                self.container.outputData.MAP_FOFC.subType = \
                    servalcatJob.container.outputData.MAP_FOFC.subType
        except Exception as e:
            self.appendErrorReport(101,
                f'Failed to apply annotations: {e}')

    def _cleanupIntermediateFiles(self):
        """Optionally purge intermediate job files."""
        try:
            cleanUpIntermediate = False
            if hasattr(self.container.controlParameters, "REFMAC_CLEANUP"):
                cleanUpIntermediate = self.container.controlParameters.REFMAC_CLEANUP

            if cleanUpIntermediate:
                from ccp4i2.core import CCP4ProjectsManager

                if self.servalcatPlugin is not None:
                    cleanup = CCP4ProjectsManager.CPurgeProject(
                        self.servalcatPlugin._dbProjectId)
                    cleanup.purgeJob(
                        self.servalcatPlugin.jobId,
                        context="extended_intermediate",
                        reportMode="skip")

                if self.servalcatPostCootPlugin is not None:
                    cleanup = CCP4ProjectsManager.CPurgeProject(
                        self.servalcatPostCootPlugin._dbProjectId)
                    cleanup.purgeJob(
                        self.servalcatPostCootPlugin.jobId,
                        context="extended_intermediate",
                        reportMode="skip")
        except Exception as e:
            self.appendErrorReport(101,
                f'Failed to cleanup intermediate files: {e}')

    # =========================================================================
    # XML management helpers
    # =========================================================================

    def _appendPluginXml(self, plugin, tag=None):
        """Safely append sub-plugin XML to our xmlroot, replacing existing tag if present."""
        try:
            pluginRoot = CCP4Utils.openFileToEtree(
                plugin.makeFileName('PROGRAMXML'))
            servalcatXML = pluginRoot.xpath("//SERVALCAT")
            if len(servalcatXML) == 1:
                node = servalcatXML[0]
                if tag:
                    node.tag = tag
                    # Remove existing element with same tag to avoid duplicates
                    existing = self.xmlroot.find(tag)
                    if existing is not None:
                        self.xmlroot.remove(existing)
                self.xmlroot.append(node)
            else:
                self.xmlroot.append(pluginRoot)
            self._flushXML()
        except Exception as e:
            self.appendErrorReport(105,
                f'Failed to append plugin XML: {e}')

    # =========================================================================
    # Real-time progress monitoring
    # =========================================================================

    def _monitorSubPluginXml(self, plugin, tag):
        """Start a background thread monitoring a sub-plugin's program.xml.

        The sub-plugin (servalcat wrapper) updates its program.xml in real-time
        via file watching on refined_stats.json. This method monitors those
        updates and mirrors them into the pipeline's program.xml with the
        appropriate tag renaming (e.g. SERVALCAT -> SERVALCAT_FIRST).

        Returns a stop function to call when monitoring should end.
        """
        stop_event = threading.Event()
        xml_path = plugin.makeFileName('PROGRAMXML')
        last_size = [0]

        def monitor():
            while not stop_event.wait(2.0):
                try:
                    if os.path.isfile(xml_path):
                        current_size = os.path.getsize(xml_path)
                        if current_size > last_size[0]:
                            last_size[0] = current_size
                            self._syncSubPluginXml(xml_path, tag)
                except Exception:
                    pass

        thread = threading.Thread(target=monitor, daemon=True,
                                  name=f"XmlMonitor-{tag}")
        thread.start()

        def stop():
            stop_event.set()
            thread.join(timeout=5)

        return stop

    def _syncSubPluginXml(self, xml_path, tag):
        """Read sub-plugin XML, rename root tag, and update pipeline XML."""
        try:
            pluginRoot = CCP4Utils.openFileToEtree(xml_path)
            servalcatXML = pluginRoot.xpath("//SERVALCAT")
            if len(servalcatXML) == 1:
                node = servalcatXML[0]
                node.tag = tag
                existing = self.xmlroot.find(tag)
                if existing is not None:
                    self.xmlroot.remove(existing)
                self.xmlroot.append(node)
                self._flushXML()
        except Exception:
            pass

    def _flushXML(self):
        """Write current XML to program.xml via atomic tmp+move."""
        try:
            tmpFileName = self.makeFileName('PROGRAMXML') + '_tmp'
            with open(tmpFileName, 'w') as f:
                CCP4Utils.writeXML(f, etree.tostring(self.xmlroot, pretty_print=True))
            shutil.move(tmpFileName, self.makeFileName('PROGRAMXML'))
        except Exception as e:
            self.appendErrorReport(105,
                f'Failed to write program.xml: {e}')


# Function called from gui to support exporting MTZ files
def exportJobFile(jobId=None, mode=None, fileInfo={}):
    from ccp4i2.core import CCP4Modules

    theDb = CCP4Modules.PROJECTSMANAGER().db()
    if mode == 'complete_mtz':
        childJobs = theDb.getChildJobs(jobId=jobId, details=True)
        if childJobs[-1][2] == 'servalcat':
            jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(
                jobId=childJobs[-1][1], create=False)
            if os.path.exists(os.path.join(jobDir, 'refined.mtz')):
                return os.path.join(jobDir, 'refined.mtz')
        elif childJobs[-1][2] == 'validate_protein':
            if childJobs[-2][2] == 'servalcat':
                jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(
                    jobId=childJobs[-2][1], create=False)
                if os.path.exists(os.path.join(jobDir, 'refined.mtz')):
                    return os.path.join(jobDir, 'refined.mtz')


# Function to return list of names of exportable MTZ(s)
def exportJobFileMenu(jobId=None):
    return [['complete_mtz', 'MTZ file', 'application/CCP4-mtz']]
