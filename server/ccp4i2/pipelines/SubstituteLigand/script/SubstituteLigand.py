import os
import shutil
import traceback

from lxml import etree

import coot_headless_api

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport


class SubstituteLigand(CPluginScript):
    """
    Pipeline for substituting a ligand into a protein structure.

    This pipeline runs synchronously through these phases:
    1. Ligand dictionary generation (LidiaAcedrg) - if ligand provided
    2. Data merging (aimless_pipe) - if unmerged data provided
    3. Rigid body refinement (phaser_rnp_pipeline or i2Dimple)
    4. Ligand fitting (coot_headless_api) - if ligand provided
    5. Anomalous map calculation (prosmart_refmac) - if anomalous data available
    """

    TASKNAME = 'SubstituteLigand'
    WHATNEXT = ['coot_rebuild']
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

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
        210: {'description': 'Failed in anomalous map calculation'},
        211: {'description': 'Missing required output from sub-plugin'},
        212: {'description': 'Invalid input configuration'},
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
        self.anomRefmacPlugin = None

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
            selAtomsFilePath = os.path.normpath(
                os.path.join(self.getWorkDirectory(), 'selected_atoms.pdb')
            )
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(selAtomsFilePath)
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
        4. coot ligand fitting (if ligand)
        5. prosmart_refmac anomalous (if anomalous data)

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
        # Phase 4: Ligand Fitting (if ligand mode)
        # =====================================================================
        if self._ligandMode != 'NONE':
            coot_error = self._runCootLigandFitting()
            if coot_error and coot_error.maxSeverity() >= 4:
                return coot_error

        # =====================================================================
        # Phase 5: Anomalous Map Calculation (if anomalous data)
        # =====================================================================
        self._checkAnomalousData()
        if self._hasAnomalous:
            # Non-fatal - continue even if this fails
            anom_error = self._runAnomalousRefinement()
            if anom_error and anom_error.maxSeverity() >= 4:
                print(f"[SubstituteLigand] Anomalous refinement failed (non-fatal): {anom_error}")
                # Don't return error - anomalous is optional

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
            plugin.container.inputData.TLC = 'DRG'

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

    def _runCootLigandFitting(self):
        """Fit ligand into density using coot_headless_api."""
        error = CErrorReport()

        try:
            xyzin = str(self.coordinatesForCoot.fullPath)
            mtzin = str(self.mapToUse.fullPath)
            dictin = str(self.dictToUse.fullPath)
            xyzout = str(self.container.outputData.XYZOUT.fullPath)

            print(f"[SubstituteLigand] Running coot ligand fitting...")
            print(f"  xyzin = {xyzin}")
            print(f"  mtzin = {mtzin}")
            print(f"  dictin = {dictin}")
            print(f"  xyzout = {xyzout}")

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
            imol_ligand = mc.get_monomer('DRG')
            if imol_ligand < 0:
                self.appendErrorReport(207, 'Failed to get monomer DRG from dictionary')
                error.append(self.__class__.__name__, 207,
                            'Failed to get monomer DRG from dictionary', 'coot', 4)
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

            # Clean up coot backup directory
            shutil.rmtree("coot-backup", ignore_errors=True)

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

    def _runAnomalousRefinement(self):
        """Run prosmart_refmac to calculate anomalous map."""
        error = CErrorReport()

        if not self._hasAnomalous or self._anomalousWavelength is None:
            return error

        # Use coordinates from dimple/phaser_rnp, NOT from coot ligand fitting
        # coordinatesForCoot holds dimple/phaser_rnp output when ligand fitting was expected
        # finalCoordinates holds dimple/phaser_rnp output when LIGANDAS=NONE
        coordinatesForAnomalous = self.coordinatesForCoot if self.coordinatesForCoot else self.finalCoordinates

        if coordinatesForAnomalous is None:
            self.appendErrorReport(210, 'No coordinates available for anomalous refinement')
            error.append(self.__class__.__name__, 210,
                        'No coordinates for anomalous refinement', 'anomalous', 3)  # Warning only
            return error

        try:
            self.anomRefmacPlugin = self.makePluginObject('prosmart_refmac')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create prosmart_refmac plugin: {e}')
            error.append(self.__class__.__name__, 209,
                        f'Failed to create prosmart_refmac: {e}', 'anomalous', 3)
            return error

        try:
            plugin = self.anomRefmacPlugin

            # Force synchronous execution (prosmart_refmac has ASYNCHRONOUS=True)
            plugin.doAsync = False

            plugin.container.inputData.XYZIN = coordinatesForAnomalous
            plugin.container.inputData.F_SIGF = self.obsToUse
            plugin.container.inputData.FREERFLAG = self.freerToUse

            plugin.container.controlParameters.WAVELENGTH.set(self._anomalousWavelength)
            plugin.container.controlParameters.NCYCLES.set(5)
            plugin.container.controlParameters.USEANOMALOUS.set(True)
            plugin.container.controlParameters.USE_NCS.set(True)

            # Disable prosmart restraints
            plugin.container.prosmartProtein.TOGGLE.set(False)
            plugin.container.prosmartNucleicAcid.TOGGLE.set(False)
            plugin.container.platonyzer.TOGGLE.set(False)

            # Disable validation
            plugin.container.controlParameters.VALIDATE_IRIS.set(False)
            plugin.container.controlParameters.VALIDATE_BAVERAGE.set(False)
            plugin.container.controlParameters.VALIDATE_RAMACHANDRAN.set(False)
            plugin.container.controlParameters.VALIDATE_MOLPROBITY.set(False)

            print(f"[SubstituteLigand] Running anomalous refinement with wavelength {self._anomalousWavelength}...")
            status = plugin.process()

            if status != CPluginScript.SUCCEEDED:
                self.appendErrorReport(210, 'Anomalous refinement failed')
                error.append(self.__class__.__name__, 210,
                            'Anomalous refinement failed', 'anomalous', 3)
                return error

            print(f"[SubstituteLigand] Anomalous refinement completed")

            # Append XML
            try:
                pluginRoot = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
                anomNode = etree.SubElement(self.xmlroot, 'AnomalousRefinement')
                refmacNodes = pluginRoot.xpath('//REFMAC')
                if len(refmacNodes) > 0:
                    anomNode.append(refmacNodes[0])
            except Exception as e:
                self.appendErrorReport(210, f'Failed to append anomalous XML: {e}')

        except Exception as e:
            self.appendErrorReport(210, f'Exception in anomalous refinement: {e}\n{traceback.format_exc()}')
            error.append(self.__class__.__name__, 210,
                        f'Exception in anomalous refinement: {e}', 'anomalous', 3)
            return error

        return error

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

            # Harvest from refinement plugin
            if self.refinementPlugin is not None:
                out = self.refinementPlugin.container.outputData

                if self._pipelineMode == 'DIMPLE':
                    # i2Dimple outputs
                    if os.path.isfile(str(out.FPHIOUT.fullPath)):
                        self._harvestFile(out.FPHIOUT, self.container.outputData.FPHIOUT)
                    if os.path.isfile(str(out.DIFFPHIOUT.fullPath)):
                        self._harvestFile(out.DIFFPHIOUT, self.container.outputData.DIFFPHIOUT)
                    if os.path.isfile(str(out.F_SIGF_OUT.fullPath)):
                        self._harvestFile(out.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                    if os.path.isfile(str(out.FREERFLAG_OUT.fullPath)):
                        self._harvestFile(out.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
                    # If no ligand, harvest coordinates from i2Dimple
                    if self._ligandMode == 'NONE' and os.path.isfile(str(out.XYZOUT.fullPath)):
                        self._harvestFile(out.XYZOUT, self.container.outputData.XYZOUT)
                else:
                    # phaser_rnp outputs
                    if hasattr(out, 'MAPOUT_REFMAC') and os.path.isfile(str(out.MAPOUT_REFMAC.fullPath)):
                        self._harvestFile(out.MAPOUT_REFMAC, self.container.outputData.FPHIOUT)
                    if hasattr(out, 'DIFMAPOUT_REFMAC') and os.path.isfile(str(out.DIFMAPOUT_REFMAC.fullPath)):
                        self._harvestFile(out.DIFMAPOUT_REFMAC, self.container.outputData.DIFFPHIOUT)
                    if os.path.isfile(str(out.F_SIGF_OUT)):
                        self._harvestFile(out.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                    if os.path.isfile(str(out.FREERFLAG_OUT)):
                        self._harvestFile(out.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)

                    # If no ligand, harvest coordinates from refinement
                    if self._ligandMode == 'NONE' and len(out.XYZOUT) > 0:
                        if os.path.isfile(str(out.XYZOUT[0].fullPath)):
                            self._harvestFile(out.XYZOUT[0], self.container.outputData.XYZOUT)

            # Harvest anomalous map if available
            if self.anomRefmacPlugin is not None:
                anomOut = self.anomRefmacPlugin.container.outputData.ANOMFPHIOUT
                if anomOut and os.path.isfile(str(anomOut.fullPath)):
                    self._harvestFile(anomOut, self.container.outputData.ANOMFPHIOUT)
                    print(f"[SubstituteLigand] Harvested anomalous map")

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
