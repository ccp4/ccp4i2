import glob
import logging
import os
import tempfile
import shutil
import traceback

from lxml import etree

try:  # external Coot API, present only in the execution (worker) env, not the slim API
    import coot_headless_api
except ImportError:
    coot_headless_api = None

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4ErrorHandling import CErrorReport

logger = logging.getLogger(__name__)


class SubstituteLigand(CPluginScript):
    """
    Pipeline for substituting a ligand into a protein structure.

    This pipeline runs synchronously through these phases:
    1. Ligand dictionary generation (LidiaAcedrg) - if ligand provided
    2. Data merging (aimless_pipe) - if unmerged data provided
    3. Rigid body refinement (phaser_rnp_pipeline or i2Dimple)
    4. Servalcat refinement - produces superior maps (normal, difference,
       and anomalous if applicable) for use by coot and as final output
    5. Ligand fitting (coot_headless_api) - if ligand provided
    """

    TASKNAME = 'SubstituteLigand'
    PERFORMANCECLASS = 'CDataReductionRefinementPerformance'
    WHATNEXT = ['coot_rebuild', 'coot1']

    ERROR_CODES = {
        201: {'description': 'Failed in SubstituteLigand'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Failed in LidiaAcedrg ligand generation'},
        204: {'description': 'Failed in aimless pipeline'},
        205: {'description': 'Failed in phaser_rnp_pipeline'},
        206: {'description': 'Failed in i2Dimple'},
        207: {'description': 'Failed in coot ligand fitting'},
        208: {'description': 'Failed in coot postprocessing'},
        209: {'description': 'Failed to create sub-plugin'},
        210: {'description': 'Failed in servalcat refinement'},
        211: {'description': 'Missing required output from sub-plugin'},
        212: {'description': 'Invalid input configuration'},
        213: {'description': 'A ligand in the starting coordinates has no dictionary'},
        214: {'description': 'The fitted ligand\'s code is already used in the starting coordinates'},
        215: {'description': 'The starting coordinates contain no atoms to refine'},
        216: {'description': 'The observed data are not supplied in the form the pipeline was told to expect'},
        217: {'description': 'The ligand chemistry is not supplied in the form the pipeline was told to expect'},
        218: {'description': 'An input has been supplied that the chosen mode will ignore'},
        219: {'description': 'No free-R set supplied with merged data'},
    }

    # The code a fitted ligand is given if the user does not choose one. It
    # used to be hard-coded, which is how a model from an earlier run --- full
    # of DRG --- collided with the DRG this run was about to make.
    DEFAULT_LIGAND_CODE = 'DRG'

    # What each way of providing the ligand's chemistry actually reads: the
    # input field, the menu entry naming it, and a phrase for the thing itself.
    # SKETCH is deliberately absent --- it draws a molecule rather than reading
    # one, so there is no file whose absence means anything. NONE likewise.
    LIGAND_INPUT_FOR_MODE = {
        'MOL': ('MOLIN', 'MDL Mol file', 'an MDL Mol file'),
        'DICT': ('DICTIN', 'REFMAC Dict', 'a REFMAC dictionary'),
        'SMILES': ('SMILESIN', 'SMILES String', 'a SMILES string'),
    }

    def __init__(self, *args, **kws):
        super(SubstituteLigand, self).__init__(*args, **kws)
        self.xmlroot = etree.Element('SubstituteLigand')

        # Intermediate state - set during startProcess()
        self.dictToUse = None
        self.obsToUse = None
        self.freerToUse = None
        self.mapToUse = None
        self.coordinatesForCoot = None
        self.finalCoordinates = None
        self.selAtomsFile = None

        # Sub-plugin references - for processOutputFiles() to access
        self.lidiaAcedrgPlugin = None
        self.aimlessPlugin = None
        self.refinementPlugin = None  # Either rnpPlugin or i2DimplePlugin
        self.servalcatPlugin = None

        # Pipeline state flags
        self._ligandMode = None  # 'DICT', 'SMILES', 'MOL', 'NONE'
        self._dataMode = None    # 'MERGED' or 'UNMERGED'
        self._pipelineMode = None  # 'PHASER' or 'DIMPLE'
        self._hasAnomalous = False
        self._anomalousWavelength = None

        # Clean up unmerged files if not needed
        if self.container.controlParameters.OBSAS.__str__() != 'UNMERGED':
            while len(self.container.inputData.UNMERGEDFILES) > 0:
                self.container.inputData.UNMERGEDFILES.remove(
                    self.container.inputData.UNMERGEDFILES[-1]
                )

    def validity(self):
        """Filter CSMILESString validation errors when mode is not SMILES."""
        from ccp4i2.core import CCP4ErrorHandling
        error = super(SubstituteLigand, self).validity()
        mode = str(self.container.controlParameters.LIGANDAS) if self.container.controlParameters.LIGANDAS.isSet() else ""
        if mode != 'SMILES':
            filtered = CCP4ErrorHandling.CErrorReport()
            for err in error.getErrors():
                err_class = err.get('class', '')
                err_name = err.get('name', '')
                if err_class == 'CSMILESString' and ('SMILES' in err_name or 'SMILESIN' in err_name):
                    continue
                filtered.append(
                    err.get('class', ''),
                    err.get('code', 0),
                    err.get('details', ''),
                    err.get('name', ''),
                    err.get('severity', 0)
                )
            return filtered
        return error

    def processInputFiles(self):
        """Prepare input files and determine pipeline configuration."""
        error = CErrorReport()

        # Validate input data
        invalidFiles = self.checkInputData()
        for invalidFile in invalidFiles:
            if (self.container.controlParameters.OBSAS.__str__() == 'MERGED' and
                invalidFile.__str__() == 'UNMERGEDFILES'):
                invalidFiles.remove(invalidFile)

        if len(invalidFiles) > 0:
            for f in invalidFiles:
                self.appendErrorReport(212, f'Missing required input: {f}')
            error.append(self.__class__.__name__, 212,
                        f'Missing required inputs: {invalidFiles}', 'processInputFiles', 4)
            return error

        # Determine pipeline modes
        self._ligandMode = str(self.container.controlParameters.LIGANDAS)
        self._dataMode = str(self.container.controlParameters.OBSAS)

        inp = self.container.inputData
        if hasattr(inp, 'PIPELINE') and inp.PIPELINE.isSet() and str(inp.PIPELINE) == 'DIMPLE':
            self._pipelineMode = 'DIMPLE'
        else:
            self._pipelineMode = 'PHASER'

        # Extract selected atoms from input structure
        try:
            # Keep the input's format here and let each sub-task convert as
            # its own program requires --- i2Dimple hands dimple the file as
            # it is, phaser gets PDB from its own wrapper. Writing an mmCIF
            # model into a file called selected_atoms.pdb, which is what this
            # did, made every DIMPLE run fail inside dimple's reader.
            selAtomsFilePath = self.container.inputData.XYZIN.getSelectedAtomsFile(
                'selected_atoms', self.getWorkDirectory())
            from ccp4i2.core.CCP4ModelData import CPdbDataFile
            self.selAtomsFile = CPdbDataFile(selAtomsFilePath)
        except Exception as e:
            self.appendErrorReport(201, f'Failed to extract selected atoms: {e}')
            error.append(self.__class__.__name__, 201,
                        f'Failed to extract selected atoms: {e}', 'processInputFiles', 4)
            return error

        return error

    def startProcess(self):
        """
        Execute the pipeline synchronously.

        This runs all sub-plugins in sequence:
        1. LidiaAcedrg (if ligand needed)
        2. aimless_pipe (if unmerged data)
        3. phaser_rnp_pipeline or i2Dimple
        4. servalcat_pipe (always - produces final maps)
        5. coot ligand fitting (if ligand)

        Returns:
            CErrorReport with any errors encountered
        """
        error = CErrorReport()

        # =====================================================================
        # Phase 1: Ligand Dictionary Generation
        # =====================================================================
        if self._ligandMode == 'DICT':
            # User provided dictionary directly
            self.dictToUse = self.container.inputData.DICTIN
            print(f"[SubstituteLigand] Using provided dictionary: {self.dictToUse.fullPath}")

        elif self._ligandMode != 'NONE':
            # Need to generate dictionary via LidiaAcedrg
            ligand_error = self._runLidiaAcedrg()
            if ligand_error and ligand_error.maxSeverity() >= 4:
                return ligand_error

        # =====================================================================
        # Phase 2: Data Merging (if unmerged)
        # =====================================================================
        if self._dataMode == 'UNMERGED':
            merge_error = self._runAimless()
            if merge_error and merge_error.maxSeverity() >= 4:
                return merge_error
        else:
            # Use provided merged data
            self.obsToUse = self.container.inputData.F_SIGF_IN
            self.freerToUse = self.container.inputData.FREERFLAG_IN

        # =====================================================================
        # Phase 3: Rigid Body Refinement
        # =====================================================================
        if self._pipelineMode == 'DIMPLE':
            refine_error = self._runDimple()
        else:
            refine_error = self._runPhaserRnp()

        if refine_error and refine_error.maxSeverity() >= 4:
            return refine_error

        # =====================================================================
        # Phase 4: Servalcat Refinement (always)
        # Produces superior maps; also anomalous map if data supports it
        # =====================================================================
        self._checkAnomalousData()
        servalcat_error = self._runServalcat()
        if servalcat_error and servalcat_error.maxSeverity() >= 4:
            return servalcat_error

        # =====================================================================
        # Phase 5: Ligand Fitting (if ligand mode)
        # =====================================================================
        if self._ligandMode != 'NONE':
            coot_error = self._runCootLigandFitting()
            if coot_error and coot_error.maxSeverity() >= 4:
                return coot_error

        return error

    def _runLidiaAcedrg(self):
        """Run LidiaAcedrg to generate ligand dictionary."""
        error = CErrorReport()

        try:
            self.lidiaAcedrgPlugin = self.makePluginObject('LidiaAcedrgNew')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create LidiaAcedrgNew plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create LidiaAcedrgNew plugin: {e}', 'LidiaAcedrg', 4)
            return error

        try:
            plugin = self.lidiaAcedrgPlugin
            plugin.container.inputData.MOLSMILESORSKETCH = self.container.controlParameters.LIGANDAS

            if self.container.inputData.MOLIN.isSet():
                plugin.container.inputData.MOLIN = self.container.inputData.MOLIN
            if self.container.inputData.SMILESIN.isSet():
                plugin.container.inputData.SMILESIN = self.container.inputData.SMILESIN

            plugin.container.inputData.CONFORMERSFROM = 'RDKIT'
            plugin.container.inputData.TLC = self._generatedLigandCode()

            print(f"[SubstituteLigand] Running LidiaAcedrg...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(203, 'LidiaAcedrg plugin failed')
                error.append(self.__class__.__name__, 203,
                            'LidiaAcedrg plugin failed', 'LidiaAcedrg', 4)
                return error

            # Verify output exists
            if (len(plugin.container.outputData.DICTOUT_LIST) == 0 or
                not os.path.isfile(str(plugin.container.outputData.DICTOUT_LIST[0].fullPath))):
                self.appendErrorReport(211, 'LidiaAcedrg did not produce dictionary output')
                error.append(self.__class__.__name__, 211,
                            'LidiaAcedrg did not produce dictionary output', 'LidiaAcedrg', 4)
                return error

            self.dictToUse = plugin.container.outputData.DICTOUT_LIST[0]
            print(f"[SubstituteLigand] LidiaAcedrg completed: {self.dictToUse.fullPath}")

            # Append XML
            self._appendPluginXml(plugin)

        except Exception as e:
            self.appendErrorReport(203, f'Exception in LidiaAcedrg: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 203,
                        f'Exception in LidiaAcedrg: {e}', 'LidiaAcedrg', 4)
            return error

        return error

    def _runAimless(self):
        """Run aimless_pipe to merge unmerged data."""
        error = CErrorReport()

        try:
            self.aimlessPlugin = self.makePluginObject('aimless_pipe')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create aimless_pipe plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create aimless_pipe plugin: {e}', 'aimless', 4)
            return error

        try:
            plugin = self.aimlessPlugin
            plugin.container.controlParameters.MODE.set('MATCH')
            plugin.container.controlParameters.RESOLUTION_RANGE = self.container.controlParameters.RESOLUTION_RANGE
            plugin.container.controlParameters.SCALING_PROTOCOL.set('DEFAULT')
            plugin.container.controlParameters.ONLYMERGE.set(False)
            plugin.container.controlParameters.REFERENCE_DATASET.set('XYZ')
            plugin.container.controlParameters.AUTOCUTOFF.set(True)
            plugin.container.controlParameters.TOLERANCE.set(10.)

            plugin.container.inputData.copyData(self.container.inputData, ['UNMERGEDFILES'])
            plugin.container.inputData.XYZIN_REF = self.container.inputData.XYZIN

            if self.container.inputData.FREERFLAG_IN.isSet():
                plugin.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG_IN

            print(f"[SubstituteLigand] Running aimless_pipe...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(204, 'Aimless pipeline failed')
                error.append(self.__class__.__name__, 204,
                            'Aimless pipeline failed', 'aimless', 4)
                return error

            # Verify outputs
            aimlessOut = plugin.container.outputData
            if not aimlessOut.FREEROUT.isSet() or not os.path.isfile(str(aimlessOut.FREEROUT.fullPath)):
                self.appendErrorReport(211, 'Aimless did not produce FreeR output')
                error.append(self.__class__.__name__, 211,
                            'Aimless did not produce FreeR output', 'aimless', 4)
                return error

            if len(aimlessOut.HKLOUT) == 0 or not os.path.isfile(str(aimlessOut.HKLOUT[0].fullPath)):
                self.appendErrorReport(211, 'Aimless did not produce merged HKL output')
                error.append(self.__class__.__name__, 211,
                            'Aimless did not produce merged HKL output', 'aimless', 4)
                return error

            self.obsToUse = aimlessOut.HKLOUT[0]
            self.freerToUse = aimlessOut.FREEROUT
            print(f"[SubstituteLigand] aimless_pipe completed")

            # Append XML
            self._appendPluginXml(plugin)

        except Exception as e:
            self.appendErrorReport(204, f'Exception in aimless: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 204,
                        f'Exception in aimless: {e}', 'aimless', 4)
            return error

        return error

    def _runPhaserRnp(self):
        """Run phaser_rnp_pipeline for rigid body refinement."""
        error = CErrorReport()

        try:
            self.refinementPlugin = self.makePluginObject('phaser_rnp_pipeline')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create phaser_rnp_pipeline plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create phaser_rnp_pipeline plugin: {e}', 'phaser_rnp', 4)
            return error

        try:
            plugin = self.refinementPlugin
            plugin.container.inputData.XYZIN_PARENT = self.selAtomsFile
            plugin.container.inputData.F_SIGF = self.obsToUse
            plugin.container.inputData.FREERFLAG = self.freerToUse
            plugin.container.inputData.SELECTIONS.append({
                'text': '/*/*/*/*',
                'pdbFileKey': 'XYZIN_PARENT'
            })

            print(f"[SubstituteLigand] Running phaser_rnp_pipeline...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(205, 'phaser_rnp_pipeline failed')
                error.append(self.__class__.__name__, 205,
                            'phaser_rnp_pipeline failed', 'phaser_rnp', 4)
                return error

            # Store results for harvest
            out = plugin.container.outputData
            self.mapToUse = out.MAPOUT_REFMAC

            if self._ligandMode == 'NONE':
                # No ligand - use refined coordinates directly
                if len(out.XYZOUT) > 0:
                    self.finalCoordinates = out.XYZOUT[0]
                else:
                    self.appendErrorReport(211, 'phaser_rnp did not produce coordinate output')
                    error.append(self.__class__.__name__, 211,
                                'phaser_rnp did not produce coordinate output', 'phaser_rnp', 4)
                    return error
            else:
                # Need coordinates for coot
                self.coordinatesForCoot = out.XYZOUT_REFMAC

            # Update obsToUse/freerToUse if plugin produced them
            if os.path.isfile(str(out.F_SIGF_OUT)):
                self.obsToUse = out.F_SIGF_OUT
            if os.path.isfile(str(out.FREERFLAG_OUT)):
                self.freerToUse = out.FREERFLAG_OUT

            print(f"[SubstituteLigand] phaser_rnp_pipeline completed")

            # Append XML
            self._appendPluginXml(plugin)

        except Exception as e:
            self.appendErrorReport(205, f'Exception in phaser_rnp: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 205,
                        f'Exception in phaser_rnp: {e}', 'phaser_rnp', 4)
            return error

        return error

    def _runDimple(self):
        """Run i2Dimple for rigid body refinement."""
        error = CErrorReport()

        try:
            self.refinementPlugin = self.makePluginObject('i2Dimple')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create i2Dimple plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create i2Dimple plugin: {e}', 'i2Dimple', 4)
            return error

        try:
            plugin = self.refinementPlugin
            plugin.container.inputData.XYZIN = self.selAtomsFile
            plugin.container.inputData.F_SIGF = self.obsToUse
            plugin.container.inputData.FREERFLAG = self.freerToUse

            print(f"[SubstituteLigand] Running i2Dimple...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(206, 'i2Dimple pipeline failed')
                error.append(self.__class__.__name__, 206,
                            'i2Dimple pipeline failed', 'i2Dimple', 4)
                return error

            # Store results for harvest
            out = plugin.container.outputData
            self.mapToUse = out.FPHIOUT

            if self._ligandMode == 'NONE':
                self.finalCoordinates = out.XYZOUT
            else:
                self.coordinatesForCoot = out.XYZOUT

            # Update obsToUse/freerToUse if dimple reindexed
            if os.path.isfile(str(out.F_SIGF_OUT.fullPath)):
                self.obsToUse = out.F_SIGF_OUT
            if os.path.isfile(str(out.FREERFLAG_OUT.fullPath)):
                self.freerToUse = out.FREERFLAG_OUT

            print(f"[SubstituteLigand] i2Dimple completed")

            # Append XML
            self._appendPluginXml(plugin)

        except Exception as e:
            self.appendErrorReport(206, f'Exception in i2Dimple: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 206,
                        f'Exception in i2Dimple: {e}', 'i2Dimple', 4)
            return error

        return error

    def _runServalcat(self):
        """Run servalcat_pipe to refine coordinates and produce maps.

        Always runs after dimple/phaser_rnp to produce superior normal and
        difference maps. Also produces anomalous maps if anomalous data
        is available (determined by prior _checkAnomalousData() call).

        Updates self.mapToUse and self.coordinatesForCoot/self.finalCoordinates
        with servalcat's output for use by subsequent coot ligand fitting.
        """
        error = CErrorReport()

        # Determine coordinates from the rigid body refinement step
        coordinatesForServalcat = self.coordinatesForCoot if self.coordinatesForCoot else self.finalCoordinates

        if coordinatesForServalcat is None:
            self.appendErrorReport(210, 'No coordinates available for servalcat refinement')
            error.append(self.__class__.__name__, 210,
                        'No coordinates for servalcat refinement', 'servalcat', 4)
            return error

        try:
            self.servalcatPlugin = self.makePluginObject('servalcat_pipe')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create servalcat_pipe plugin: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create servalcat_pipe: {e}', 'servalcat', 4)
            return error

        try:
            plugin = self.servalcatPlugin

            # Force synchronous execution
            plugin.doAsync = False

            # Input data
            plugin.container.inputData.XYZIN = coordinatesForServalcat
            plugin.container.inputData.HKLIN = self.obsToUse
            plugin.container.inputData.FREERFLAG = self.freerToUse

            # Ligands that were already in the starting model reach servalcat
            # in these coordinates --- the ligand being substituted does not,
            # because it is fitted afterwards. Servalcat's pre-flight monomer
            # check refuses a residue it has no dictionary for, so without
            # this the pipeline stops at refinement with no way past it: an
            # old ligand that genuinely belongs (a second soak on a crystal
            # that already had one) had nowhere to declare itself.
            for dictionary in self._startingLigandDictionaries():
                plugin.container.inputData.DICT_LIST.append(dictionary)

            # Not the dictionary for the ligand being fitted: that ligand is
            # not in these coordinates --- coot places it after refinement ---
            # and if the starting model happened to hold a residue under the
            # same code, this dictionary would be applied to *that*, which is
            # the mistake the clash check exists to prevent.

            # Refinement configuration
            plugin.container.controlParameters.NCYCLES.set(5)
            plugin.container.controlParameters.DATA_METHOD.set('xtal')
            plugin.container.controlParameters.MERGED_OR_UNMERGED.set('merged')

            # Anomalous map calculation if data supports it
            plugin.container.controlParameters.USEANOMALOUS.set(self._hasAnomalous)
            if self._hasAnomalous:
                plugin.container.controlParameters.USEANOMALOUSFOR.set('OUTPUTMAPS')

            # Disable ProSMART restraints
            plugin.container.prosmartProtein.TOGGLE.set(False)
            plugin.container.prosmartNucleicAcid.TOGGLE.set(False)

            # Disable MetalCoord
            plugin.container.metalCoordPipeline.RUN_METALCOORD.set(False)

            # Disable water finding
            plugin.container.controlParameters.ADD_WATERS.set(False)

            # Disable validation (unnecessary for intermediate step)
            plugin.container.controlParameters.VALIDATE_IRIS.set(False)
            plugin.container.controlParameters.VALIDATE_BAVERAGE.set(False)
            plugin.container.controlParameters.VALIDATE_RAMACHANDRAN.set(False)
            plugin.container.controlParameters.VALIDATE_MOLPROBITY.set(False)
            plugin.container.controlParameters.RUN_ADP_ANALYSIS.set(False)
            plugin.container.monitor.RUN_COORDADPDEV_ANALYSIS.set(False)

            print(f"[SubstituteLigand] Running servalcat refinement (anomalous={self._hasAnomalous})...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(210, 'Servalcat refinement failed')
                error.append(self.__class__.__name__, 210,
                            'Servalcat refinement failed', 'servalcat', 4)
                return error

            # Update map and coordinates for coot
            out = plugin.container.outputData
            self.mapToUse = out.FPHIOUT

            if self._ligandMode == 'NONE':
                self.finalCoordinates = out.XYZOUT
            else:
                self.coordinatesForCoot = out.XYZOUT

            print(f"[SubstituteLigand] Servalcat refinement completed")

            # Append XML
            self._appendPluginXml(plugin)

        except Exception as e:
            self.appendErrorReport(210, f'Exception in servalcat refinement: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 210,
                        f'Exception in servalcat refinement: {e}', 'servalcat', 4)
            return error

        return error

    def runTimeValidity(self):
        """Refuse at the Run dialog what the pipeline would refuse steps in.

        Four things are checked, and they are independent --- fixing one does
        not fix the others.

        *Agreement between the mode menus and the inputs.* Two menus say what
        form the job's inputs take: observed data merged or unmerged, and the
        ligand's chemistry as a Mol file, a dictionary or a SMILES string.
        Neither menu carries the input it names, and neither input is
        individually required, so the pair could disagree freely: 'unmerged'
        with no unmerged file, 'SMILES String' with no string. That got as far
        as aimless being handed an empty file list, or acedrg being asked to
        build a molecule from nothing, and the campaign that found this had
        fifty such jobs fail one at a time. It is now refused here, where the
        menu and the empty field are both on screen.

        *Something to refine.* Coordinates that cannot be read, or that hold no
        atoms once the atom selection has been applied, give molecular
        replacement nothing to place. (Coordinates that are simply not set, or
        name a file that is not there, are already refused by the container's
        own checks, which run first.)

        *Coverage.* Every residue in the starting coordinates that the monomer
        library does not describe needs a dictionary, or servalcat's monomer
        check stops the pipeline after ligand generation, merging and
        molecular replacement have run. Only dictionaries declared as being
        for the starting model count: the one describing the ligand being
        fitted describes a different molecule and covers nothing here.

        *Identity.* If such a residue carries the code the fitted ligand will
        carry, the two would be one name for two molecules. That is not a
        missing dictionary and supplying one does not resolve it.

        Both are checked against the atoms that will actually be refined ---
        the selection, not the file as supplied --- so excluding an unwanted
        ligand with the atom selection is a real answer and clears the error.
        """
        error = super(SubstituteLigand, self).runTimeValidity()

        # Whether to go on to the checks that read files is decided by the
        # container's own verdict alone, and decided *before* the mode checks
        # add to the report. The mode checks cost nothing and are independent
        # of the model, so a disagreeing menu must neither be hidden by a
        # missing model nor hide one: the user should see everything they have
        # to fix in one pass, not one thing per attempt.
        containerBlocked = error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR

        self._checkObservedDataMatchesMode(error)
        self._checkLigandChemistryMatchesMode(error)

        if containerBlocked:
            return error

        xyzin = self.container.inputData.XYZIN
        if not xyzin.isSet():
            return error

        workDirectory = tempfile.mkdtemp(prefix='ccp4i2_substituteligand_check_')
        try:
            model = self._selectedStartingModel(workDirectory)
            if not self._checkThereAreAtomsToRefine(model, error):
                # Coverage and code-clash both ask what is in the model. With
                # nothing in it they would pass, and read as approval.
                return error
            self._checkStartingLigandsAreDescribed(model, error)
            self._checkNoCodeClash(model, error)
        finally:
            shutil.rmtree(workDirectory, ignore_errors=True)
        return error

    def _checkObservedDataMatchesMode(self, error) -> None:
        """The observed data must arrive in the form the menu promises.

        The two branches are not symmetrical. Unmerged data go to aimless,
        which needs at least one file; merged data go straight to refinement,
        which needs the reflections but only wants the free-R set --- the
        pipeline runs without one, it just cannot report R-free, so that is an
        advisory rather than a refusal.
        """
        mode = str(self.container.controlParameters.OBSAS)
        inp = self.container.inputData

        if mode == 'UNMERGED':
            if len(inp.UNMERGEDFILES) == 0:
                error.append(
                    klass=self.TASKNAME,
                    code=216,
                    details=(
                        "'Observed data provided as' is set to unmerged, but no "
                        "unmerged file has been given, so there is nothing for "
                        "aimless to merge.\n\n"
                        "Either add the unmerged reflections, or set that menu to "
                        "merged and supply merged reflections instead."
                    ),
                    name=f'{self.TASKNAME}.container.inputData.UNMERGEDFILES',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR,
                )
            return

        if not inp.F_SIGF_IN.isSet():
            error.append(
                klass=self.TASKNAME,
                code=216,
                details=(
                    "'Observed data provided as' is set to merged, but no merged "
                    "reflection file has been given.\n\n"
                    "Either supply the merged reflections, or set that menu to "
                    "unmerged and supply the unmerged files, which the pipeline "
                    "will merge with aimless for you."
                ),
                name=f'{self.TASKNAME}.container.inputData.F_SIGF_IN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
            return

        if not inp.FREERFLAG_IN.isSet():
            error.append(
                klass=self.TASKNAME,
                code=219,
                details=(
                    "No free-R set has been given with the merged data. Nothing "
                    "generates one on this route --- it is aimless that would, and "
                    "aimless only runs on unmerged data --- so refinement will "
                    "proceed with no cross-validation and the job will report no "
                    "R-free.\n\n"
                    "Supply the free-R set that goes with these reflections where "
                    "there is one; a set made now would not be independent of a "
                    "model already refined against this data."
                ),
                name=f'{self.TASKNAME}.container.inputData.FREERFLAG_IN',
                severity=CCP4ErrorHandling.SEVERITY_WARNING,
            )

    def _checkLigandChemistryMatchesMode(self, error) -> None:
        """The ligand's chemistry must arrive in the form the menu promises.

        Both directions matter. The input the mode reads has to be there, or
        acedrg is asked to build a molecule from nothing; and an input some
        *other* mode would have read is not an error but is worth saying,
        because it is silently ignored and looks for all the world like it was
        used.
        """
        mode = str(self.container.controlParameters.LIGANDAS)
        inp = self.container.inputData
        expected = self.LIGAND_INPUT_FOR_MODE.get(mode)

        if expected is not None:
            field, menuEntry, description = expected
            supplied = getattr(inp, field, None)
            if supplied is None or not supplied.isSet():
                error.append(
                    klass=self.TASKNAME,
                    code=217,
                    details=(
                        f"The ligand's chemistry is set to come from '{menuEntry}', "
                        f"but {description} has not been given.\n\n"
                        "Either supply it, or change how the chemistry is provided "
                        "-- to 'NONE' if no ligand is to be fitted, which runs the "
                        "refinement and leaves the site empty."
                    ),
                    name=f'{self.TASKNAME}.container.inputData.{field}',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR,
                )

        for otherMode, (field, menuEntry, description) in sorted(self.LIGAND_INPUT_FOR_MODE.items()):
            if otherMode == mode:
                continue
            supplied = getattr(inp, field, None)
            if supplied is None or not supplied.isSet():
                continue
            chosen = ("no ligand is to be fitted" if mode == 'NONE'
                      else f"the chemistry is set to come from '{self._ligandMenuEntry(mode)}'")
            error.append(
                klass=self.TASKNAME,
                code=218,
                details=(
                    f"{description[:1].upper()}{description[1:]} has been given, "
                    f"but {chosen}, so it will be ignored.\n\n"
                    f"Set 'Chemistry of the ligand' to '{menuEntry}' to use it, or "
                    "clear the field so the job does not appear to read something "
                    "it does not."
                ),
                name=f'{self.TASKNAME}.container.inputData.{field}',
                severity=CCP4ErrorHandling.SEVERITY_WARNING,
            )

    def _ligandMenuEntry(self, mode) -> str:
        """The menu wording for a ligand mode, for use in messages."""
        if mode in self.LIGAND_INPUT_FOR_MODE:
            return self.LIGAND_INPUT_FOR_MODE[mode][1]
        return 'LIDIA Sketch' if mode == 'SKETCH' else str(mode)

    def _checkThereAreAtomsToRefine(self, model, error) -> bool:
        """Whether the model to be refined holds any atoms.

        Returns False --- and reports --- when it holds none, or cannot be read
        at all. An unreadable file is refused rather than waved through: it is
        exactly the file a half-finished download leaves behind, and the whole
        pipeline is built on placing these coordinates.
        """
        atoms = self._atomCountIn(model)
        if atoms is None:
            error.append(
                klass=self.TASKNAME,
                code=215,
                details=(
                    "The starting coordinates could not be read, so there is no "
                    "model for molecular replacement to place.\n\n"
                    "Check that the file is complete and really is a coordinate "
                    "file."
                ),
                name=f'{self.TASKNAME}.container.inputData.XYZIN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
            return False

        if atoms == 0:
            selected = self.container.inputData.XYZIN.isSelectionSet()
            details = (
                (
                    "The atom selection on the starting coordinates leaves no "
                    "atoms, so there is nothing to refine. The selection is "
                    f"'{self.container.inputData.XYZIN.selection}'.\n\n"
                    "Widen it, or clear it to refine the model as supplied."
                ) if selected else (
                    "The starting coordinates contain no atoms, so there is "
                    "nothing for molecular replacement to place.\n\n"
                    "The file was read, so it is not corrupt --- it simply holds "
                    "no coordinates. An interrupted download and an error page "
                    "saved under a .pdb name both look like this."
                )
            )
            error.append(
                klass=self.TASKNAME,
                code=215,
                details=details,
                name=f'{self.TASKNAME}.container.inputData.XYZIN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )
            return False

        return True

    @staticmethod
    def _atomCountIn(path):
        """How many atoms a coordinate file holds, or None if it cannot be read."""
        try:
            import gemmi
            structure = gemmi.read_structure(str(path))
            return sum(len(residue) for model in structure
                       for chain in model for residue in chain)
        except Exception as err:
            logger.warning(f"[SubstituteLigand] Could not count atoms in {path}: {err}")
            return None

    def _selectedStartingModel(self, workDirectory) -> str:
        """The starting coordinates as they will be refined.

        The pipeline refines ``getSelectedAtomsFile``'s output, so a ligand
        the user has excluded is genuinely not there. Checking the file as
        supplied would refuse a model that has already been fixed --- and the
        atom selection is the remedy these checks recommend.
        """
        xyzin = self.container.inputData.XYZIN
        try:
            return xyzin.getSelectedAtomsFile('runtime_validity_check', workDirectory)
        except Exception as err:
            logger.warning(f"[SubstituteLigand] Could not apply the atom selection "
                           f"for validation, checking the file as supplied: {err}")
            return str(xyzin.fullPath)

    def _checkStartingLigandsAreDescribed(self, model, error) -> None:
        """Every ligand already in the model must have a dictionary."""
        saved = self.errorReport
        self.errorReport = CErrorReport()
        try:
            self.checkMonomeCoverage(model, self._startingLigandDictionaries())
            found = self.errorReport
        finally:
            self.errorReport = saved

        present = self._residueCodesIn(model) or set()
        for entry in found.entries():
            offending = self._codesMentionedIn(entry['details'], present)
            error.append(
                klass=self.TASKNAME,
                code=213,
                details=(
                    f"{entry['details']}\n\n"
                    "Refinement cannot describe a ligand it has no dictionary for. "
                    "Either supply its dictionary under 'Dictionaries for ligands "
                    "already in this model', or exclude it from the starting "
                    "coordinates with the atom selection on that file"
                    f"{self._selectionAdvice(offending)}."
                ),
                name=f'{self.TASKNAME}.container.inputData.XYZIN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )

    def _checkNoCodeClash(self, model, error) -> None:
        """The fitted ligand's code must not already be in use in the model."""
        fitted = self._newLigandCodes()
        if not fitted:
            return
        present = self._residueCodesIn(model)
        if present is None:
            return
        for code in sorted(fitted & present):
            error.append(
                klass=self.TASKNAME,
                code=214,
                details=(
                    f"The starting coordinates already contain a residue called "
                    f"{code}, and {code} is the code the ligand being fitted will "
                    "carry. The two are probably not the same molecule, and one "
                    "name cannot mean both: whichever dictionary was read last "
                    "would describe them both.\n\n"
                    "Give the fitted ligand a different code (LIG is the usual "
                    "second choice where DRG is taken), or remove the existing "
                    f"{code} from the starting coordinates with the atom selection"
                    f"{self._selectionAdvice({code})} -- which is the answer when "
                    "it is the same ligand being repositioned against new data."
                ),
                name=f'{self.TASKNAME}.container.controlParameters.LIGAND_CODE',
                severity=CCP4ErrorHandling.SEVERITY_ERROR,
            )

    @staticmethod
    def _codesMentionedIn(details, present) -> set:
        """Which of the model's residue codes a coverage message names.

        Read back from the message rather than recomputed, but intersected
        with what is actually in the model, so no advice can name a residue
        that is not there.
        """
        import re

        return {code for code in present
                if re.search(rf'\b{re.escape(str(code))}\b', details or '')}

    @staticmethod
    def _selectionAdvice(codes) -> str:
        """The selection that would exclude *codes*, ready to paste.

        "not (DRG) and not (LIG)" and not "not (DRG or LIG)" --- the latter
        parses, selects everything, and would read as advice that did nothing.
        """
        wanted = sorted(str(c) for c in codes if str(c))
        if not wanted:
            return ""
        selection = " and ".join(f"not ({code})" for code in wanted)
        return f" -- '{selection}'"

    def _newLigandCodes(self) -> set:
        """The code(s) the ligand being fitted will carry; empty if none is."""
        mode = str(self.container.controlParameters.LIGANDAS)
        if mode == 'NONE':
            return set()
        if mode == 'DICT':
            dictin = self.container.inputData.DICTIN
            # A supplied dictionary names its own residue --- the pipeline does
            # not get to choose it, so the clash check has to read it.
            return self._dictionaryCodes(str(dictin.fullPath)) if dictin.isSet() else set()
        return {self._generatedLigandCode()}

    def _dictionaryCodes(self, path) -> set:
        """The residue codes a restraint dictionary describes."""
        try:
            import gemmi
            document = gemmi.cif.read(str(path))
            codes = set()
            for block in document:
                for row in block.find('_chem_comp.', ['id']):
                    codes.add(row.str(0).strip().upper())
            return codes
        except Exception as err:
            logger.warning(f"[SubstituteLigand] Could not read residue codes from "
                           f"dictionary {path}: {err}")
            return set()

    def _generatedLigandCode(self) -> str:
        """The code the fitted ligand will be given."""
        chosen = str(self.container.controlParameters.LIGAND_CODE or '').strip().upper()
        return chosen or self.DEFAULT_LIGAND_CODE

    def _residueCodesIn(self, path):
        """Every residue code in a coordinate file, or None if it cannot be read.

        None rather than an empty set on failure: an unreadable file must not
        read as *contains nothing*, which would quietly pass a check whose
        whole purpose is to find something.
        """
        try:
            import gemmi
            structure = gemmi.read_structure(str(path))
            return {residue.name for model in structure
                    for chain in model for residue in chain}
        except Exception as err:
            logger.warning(f"[SubstituteLigand] Could not read residue codes from {path}: {err}")
            return None

    def _startingLigandDictionaries(self):
        """Paths of the dictionaries supplied for ligands already in XYZIN.

        Paths rather than the objects themselves: appending a CDataFile to
        another container's list re-parents it, which would take the input
        out of this pipeline's own container.
        """
        supplied = self.container.inputData.STARTING_DICT_LIST
        return [str(d.fullPath) for d in supplied if d.isSet()]

    def _runCootLigandFitting(self):
        """Fit ligand into density using coot_headless_api."""
        error = CErrorReport()

        xyzin = str(self.coordinatesForCoot.fullPath)
        mtzin = str(self.mapToUse.fullPath)
        dictin = str(self.dictToUse.fullPath)
        xyzout = str(self.container.outputData.XYZOUT.fullPath)

        print(f"[SubstituteLigand] Running coot ligand fitting...")
        print(f"  xyzin = {xyzin}")
        print(f"  mtzin = {mtzin}")
        print(f"  dictin = {dictin}")
        print(f"  xyzout = {xyzout}")

        # Change to job directory so coot's intermediate files are created there
        # (coot_headless_api writes map files to CWD during fit_ligand)
        original_cwd = os.getcwd()
        job_dir = os.path.dirname(xyzout)
        if job_dir:
            os.chdir(job_dir)
            print(f"  Changed to job directory: {job_dir}")

        try:

            # Initialize coot headless API
            mc = coot_headless_api.molecules_container_py(True)
            mc.set_make_backups(False)
            mc.set_use_gemmi(False)

            # Load protein coordinates
            imol_protein = mc.read_pdb(xyzin)
            if imol_protein < 0:
                self.appendErrorReport(207, f'Failed to read protein coordinates: {xyzin}')
                error.append(self.__class__.__name__, 207,
                            f'Failed to read protein coordinates: {xyzin}', 'coot', 4)
                return error

            # Load map from MTZ - try multiple column naming conventions
            imol_map = mc.read_mtz(mtzin, "F", "PHI", "", False, False)
            if imol_map < 0:
                imol_map = mc.read_mtz(mtzin, "2FOFCWT", "PH2FOFCWT", "", False, False)
            if imol_map < 0:
                imol_map = mc.read_mtz(mtzin, "FWT", "PHWT", "", False, False)

            if imol_map < 0:
                self.appendErrorReport(207, f'Failed to read map coefficients: {mtzin}')
                error.append(self.__class__.__name__, 207,
                            f'Failed to read map coefficients: {mtzin}', 'coot', 4)
                return error

            mc.set_imol_refinement_map(imol_map)

            # Import ligand dictionary
            dict_result = mc.import_cif_dictionary(dictin, -999999)
            if dict_result == 0:
                self.appendErrorReport(207, f'Failed to import ligand dictionary: {dictin}')
                error.append(self.__class__.__name__, 207,
                            f'Failed to import ligand dictionary: {dictin}', 'coot', 4)
                return error

            # Get monomer for fitting
            ligandCode = self._generatedLigandCode()
            imol_ligand = mc.get_monomer(ligandCode)
            if imol_ligand < 0:
                self.appendErrorReport(207, f'Failed to get monomer {ligandCode} from dictionary')
                error.append(self.__class__.__name__, 207,
                            f'Failed to get monomer {ligandCode} from dictionary', 'coot', 4)
                return error

            # Fit ligand into density
            ligands_found = mc.fit_ligand(imol_protein, imol_map, imol_ligand, 1.5, True, 20)

            if ligands_found and len(ligands_found) > 0:
                first_lig = ligands_found[0]
                print(f"  Found {len(ligands_found)} ligand positions")

                # Determine how many ligands to merge based on NCS
                n_to_copy = 1
                try:
                    ncs_chains = mc.get_ncs_related_chains(imol_protein)
                    if ncs_chains and len(ncs_chains) > 0:
                        n_to_copy = min(len(ligands_found), len(ncs_chains[0]))
                except Exception:
                    n_to_copy = 1

                # Merge found ligands into protein model
                if n_to_copy > 0:
                    ligand_imols = [lig.imol for lig in ligands_found[:n_to_copy]]
                    ligand_indices = ','.join(str(imol) for imol in ligand_imols)
                    mc.merge_molecules(imol_protein, ligand_indices)
            else:
                print("  No ligands found by fit_ligand")

            # Write output coordinates
            mc.write_coordinates(imol_protein, xyzout)

            # Clean up coot intermediate files
            shutil.rmtree("coot-backup", ignore_errors=True)
            # Remove map files generated by fit_ligand() in CWD
            for pattern in ["molecules-container-fit-ligand-*.map", "wlig.output_map.map"]:
                for mapfile in glob.glob(pattern):
                    try:
                        os.remove(mapfile)
                    except OSError:
                        pass

            # Verify output
            if not os.path.isfile(xyzout):
                self.appendErrorReport(207, 'Coot did not produce output coordinates')
                error.append(self.__class__.__name__, 207,
                            'Coot did not produce output coordinates', 'coot', 4)
                return error

            self.finalCoordinates = self.container.outputData.XYZOUT
            print(f"[SubstituteLigand] Coot ligand fitting completed")

        except Exception as e:
            self.appendErrorReport(207, f'Exception in coot ligand fitting: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 207,
                        f'Exception in coot ligand fitting: {e}', 'coot', 4)
            return error

        finally:
            # Ensure we restore working directory even on error
            if job_dir:
                try:
                    os.chdir(original_cwd)
                except Exception:
                    pass

        return error

    def _checkAnomalousData(self):
        """Check if observation data has anomalous signal and extract wavelength."""
        self._hasAnomalous = False
        self._anomalousWavelength = None

        if self.obsToUse is None:
            return

        try:
            from ccp4i2.core.CCP4XtalData import CObsDataFile

            # Check contentFlag - must be Ipair (1) or Fpair (2)
            contentFlag = None
            if hasattr(self.obsToUse, 'contentFlag') and self.obsToUse.contentFlag:
                cf = self.obsToUse.contentFlag
                while hasattr(cf, 'value'):
                    cf = cf.value
                contentFlag = int(cf) if cf else None

            if contentFlag not in [CObsDataFile.CONTENT_FLAG_IPAIR, CObsDataFile.CONTENT_FLAG_FPAIR]:
                return

            # Load file to access fileContent
            self.obsToUse.loadFile()

            wavelength = None
            if hasattr(self.obsToUse, 'fileContent') and self.obsToUse.fileContent:
                fc = self.obsToUse.fileContent
                if hasattr(fc, 'getListOfWavelengths'):
                    try:
                        wavelengths = fc.getListOfWavelengths()
                        if wavelengths and len(wavelengths) > 0:
                            wavelength = wavelengths[-1]
                    except Exception:
                        pass

                if wavelength is None and hasattr(fc, 'wavelength'):
                    wl = fc.wavelength
                    if hasattr(wl, 'value'):
                        wl = wl.value
                    if wl and float(wl) > 0:
                        wavelength = float(wl)

            # Validate wavelength range
            if wavelength is not None and 0.1 < wavelength < 10.0:
                self._hasAnomalous = True
                self._anomalousWavelength = wavelength
                print(f"[SubstituteLigand] Detected anomalous data with wavelength {wavelength}")

        except Exception as e:
            print(f"[SubstituteLigand] Error checking for anomalous data: {e}")

    def _appendPluginXml(self, plugin):
        """Safely append sub-plugin XML to our xmlroot."""
        try:
            pluginRoot = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
        except Exception as e:
            self.appendErrorReport(201, f'Failed to append plugin XML: {e}')

    def processOutputFiles(self):
        """Harvest output files from sub-plugins."""
        error = CErrorReport()

        try:
            # Harvest dictionary if we generated one
            if self.lidiaAcedrgPlugin is not None and self.dictToUse is not None:
                if os.path.isfile(str(self.dictToUse.fullPath)):
                    self._harvestFile(self.dictToUse, self.container.outputData.DICTOUT)

            # Harvest from aimless if we ran it
            if self.aimlessPlugin is not None:
                aimlessOut = self.aimlessPlugin.container.outputData
                if os.path.isfile(str(aimlessOut.FREEROUT.fullPath)):
                    self._harvestFile(aimlessOut.FREEROUT, self.container.outputData.FREERFLAG_OUT)
                if len(aimlessOut.HKLOUT) > 0 and os.path.isfile(str(aimlessOut.HKLOUT[0].fullPath)):
                    self._harvestFile(aimlessOut.HKLOUT[0], self.container.outputData.F_SIGF_OUT)

            # Harvest F_SIGF_OUT and FREERFLAG_OUT from refinement plugin
            # (these may have been reindexed by dimple/phaser_rnp)
            if self.refinementPlugin is not None:
                out = self.refinementPlugin.container.outputData
                if self._pipelineMode == 'DIMPLE':
                    if os.path.isfile(str(out.F_SIGF_OUT.fullPath)):
                        self._harvestFile(out.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                    if os.path.isfile(str(out.FREERFLAG_OUT.fullPath)):
                        self._harvestFile(out.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
                else:
                    if os.path.isfile(str(out.F_SIGF_OUT)):
                        self._harvestFile(out.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                    if os.path.isfile(str(out.FREERFLAG_OUT)):
                        self._harvestFile(out.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)

            # Harvest maps and coordinates from servalcat
            if self.servalcatPlugin is not None:
                sOut = self.servalcatPlugin.container.outputData
                if os.path.isfile(str(sOut.FPHIOUT.fullPath)):
                    self._harvestFile(sOut.FPHIOUT, self.container.outputData.FPHIOUT)
                if os.path.isfile(str(sOut.DIFFPHIOUT.fullPath)):
                    self._harvestFile(sOut.DIFFPHIOUT, self.container.outputData.DIFFPHIOUT)
                if hasattr(sOut, 'ANOMFPHIOUT') and os.path.isfile(str(sOut.ANOMFPHIOUT.fullPath)):
                    self._harvestFile(sOut.ANOMFPHIOUT, self.container.outputData.ANOMFPHIOUT)
                    print(f"[SubstituteLigand] Harvested anomalous map from servalcat")
                # If no ligand, harvest coordinates from servalcat
                if self._ligandMode == 'NONE' and os.path.isfile(str(sOut.XYZOUT.fullPath)):
                    self._harvestFile(sOut.XYZOUT, self.container.outputData.XYZOUT)

            # Harvest performance indicators
            perf = self.container.outputData.PERFORMANCE
            # Data reduction KPIs from aimless_pipe
            if self.aimlessPlugin is not None:
                drPerf = getattr(self.aimlessPlugin.container.outputData, 'PERFORMANCE', None)
                if drPerf is not None:
                    if hasattr(drPerf, 'spaceGroup') and drPerf.spaceGroup.isSet():
                        perf.spaceGroup.set(drPerf.spaceGroup)
                    if hasattr(drPerf, 'highResLimit') and drPerf.highResLimit.isSet():
                        perf.highResLimit.set(drPerf.highResLimit)
                    if hasattr(drPerf, 'rMeas') and drPerf.rMeas.isSet():
                        perf.rMeas.set(drPerf.rMeas)
            # Refinement KPIs from servalcat
            if self.servalcatPlugin is not None:
                sPerf = getattr(self.servalcatPlugin.container.outputData, 'PERFORMANCEINDICATOR', None)
                if sPerf is not None:
                    if hasattr(sPerf, 'RFactor') and sPerf.RFactor.isSet():
                        perf.RFactor.set(sPerf.RFactor)
                    if hasattr(sPerf, 'RFree') and sPerf.RFree.isSet():
                        perf.RFree.set(sPerf.RFree)

            # Update XML with model composition
            self._updateModelComposition()

            # Flush final XML
            self._flushXML()

        except Exception as e:
            self.appendErrorReport(202, f'Exception in processOutputFiles: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 202,
                        f'Exception in processOutputFiles: {e}', 'harvest', 3)

        return error

    def _harvestFile(self, sourceFile, destFile):
        """Copy file and metadata from source to destination."""
        try:
            shutil.copyfile(str(sourceFile.fullPath), str(destFile.fullPath))
            # Extract primitive values for proper smart assignment
            destFile.annotation = str(sourceFile.annotation)
            destFile.contentFlag = int(sourceFile.contentFlag)
            destFile.subType = int(sourceFile.subType)
        except Exception as e:
            self.appendErrorReport(202, f'Failed to harvest {sourceFile.fullPath} -> {destFile.fullPath}: {e}')

    def _updateModelComposition(self):
        """Update XML with monomer composition from output coordinates."""
        try:
            xyzout_path = str(self.container.outputData.XYZOUT.fullPath)
            if os.path.isfile(xyzout_path):
                from ccp4i2.core.CCP4ModelData import CPdbData
                aCPdbData = CPdbData()
                aCPdbData.loadFile(xyzout_path)

                modelCompositionNode = None
                modelCompositionNodes = self.xmlroot.xpath('//ModelComposition')
                if len(modelCompositionNodes) > 0:
                    modelCompositionNode = modelCompositionNodes[-1]
                else:
                    refmacNodes = self.xmlroot.xpath('//REFMAC')
                    if len(refmacNodes) > 0:
                        modelCompositionNode = etree.SubElement(refmacNodes[-1], "ModelComposition")

                if modelCompositionNode is not None:
                    for monomer in aCPdbData.composition.monomers:
                        etree.SubElement(modelCompositionNode, 'Monomer', id=monomer)
        except Exception as e:
            self.appendErrorReport(208, f'Exception updating model composition: {e}')

    def _flushXML(self):
        """Write current XML to program.xml."""
        try:
            with open(self.makeFileName('PROGRAMXML'), 'w') as programXML:
                CCP4Utils.writeXML(programXML, etree.tostring(self.xmlroot, pretty_print=True))
        except Exception as e:
            self.appendErrorReport(201, f'Failed to write program.xml: {e}')
