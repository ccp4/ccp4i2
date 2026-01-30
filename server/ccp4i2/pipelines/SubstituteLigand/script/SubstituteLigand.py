import os
import shutil

from lxml import etree

import coot_headless_api

from ccp4i2.baselayer import QtCore
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class SubstituteLigand(CPluginScript):
    TASKNAME = 'SubstituteLigand'            # Task name - should be same as class name
    TASKVERSION= 0.0                    # Version of this plugin
    TIMEOUT_PERIOD = 9999999.9
    WHATNEXT = ['coot_rebuild']
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {
        201: {'description': 'Failed in SubstituteLigand'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Exception in LidiaAcedrg ligand generation'},
        204: {'description': 'Exception in aimless pipeline'},
        205: {'description': 'Exception in phaser_rnp_pipeline'},
        206: {'description': 'Exception in i2Dimple'},
        207: {'description': 'Exception in coot ligand fitting'},
        208: {'description': 'Exception in coot postprocessing'},
        209: {'description': 'Failed to create sub-plugin'},
    }

    def __init__(self, *args,**kws):
        super(SubstituteLigand, self).__init__(*args, **kws)
        self.xmlroot = etree.Element('SubstituteLigand')
        self.obsToUse = None
        self.freerToUse = None

        if self.container.controlParameters.OBSAS.__str__() != 'UNMERGED':
            #remove any (potentially invalid) entries from UNMERGED list
            while len(self.container.inputData.UNMERGEDFILES)>0:
                self.container.inputData.UNMERGEDFILES.remove(self.container.inputData.UNMERGEDFILES[-1])

    def validity(self):
        """Filter CSMILESString validation errors when mode is not SMILES."""
        from ccp4i2.core import CCP4ErrorHandling
        error = super(SubstituteLigand, self).validity()
        mode = str(self.container.controlParameters.LIGANDAS) if self.container.controlParameters.LIGANDAS.isSet() else ""
        if mode != 'SMILES':
            # Filter out SMILES validation errors when not using SMILES input
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

    def process(self):
        invalidFiles = self.checkInputData()
        for invalidFile in invalidFiles:
            if self.container.controlParameters.OBSAS.__str__() == 'MERGED' and invalidFile.__str__() == 'UNMERGEDFILES':
                invalidFiles.remove(invalidFile)
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()

        # Chop out the chunk of file we want to use
        selAtomsFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected_atoms.pdb'))
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(selAtomsFilePath)
        from ccp4i2.core.CCP4ModelData import CPdbDataFile
        self.selAtomsFile = CPdbDataFile(selAtomsFilePath)
        
        if self.container.controlParameters.LIGANDAS.__str__() == 'DICT':
            self.dictToUse = self.container.inputData.DICTIN
            self.dictDone()
        elif self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
            self.dictDone()
        elif self.container.controlParameters.LIGANDAS.__str__() != 'NONE':
            try:
                self.lidiaAcedrgPlugin = self.makePluginObject('LidiaAcedrgNew')
            except Exception as e:
                self.appendErrorReport(209, f'Failed to create LidiaAcedrgNew plugin: {e}')
                self.reportStatus(CPluginScript.FAILED)
                return
            self.lidiaAcedrgPlugin.container.inputData.MOLSMILESORSKETCH = self.container.controlParameters.LIGANDAS
            if self.container.inputData.MOLIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.MOLIN=self.container.inputData.MOLIN
            if self.container.inputData.SMILESIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.SMILESIN=self.container.inputData.SMILESIN
            self.lidiaAcedrgPlugin.container.inputData.CONFORMERSFROM = 'RDKIT'
            self.lidiaAcedrgPlugin.container.inputData.TLC='DRG'
            self.connectSignal(self.lidiaAcedrgPlugin,'finished',self.lidiaAcedrg_finished)
            self.lidiaAcedrgPlugin.process()

    @QtCore.Slot(dict)
    def lidiaAcedrg_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.appendErrorReport(203, 'LidiaAcedrg plugin failed')
            self.reportStatus(CPluginScript.FAILED)
            return
        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.lidiaAcedrgPlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()
            self.harvestFile(self.lidiaAcedrgPlugin.container.outputData.DICTOUT_LIST[0], self.container.outputData.DICTOUT)
            self.dictToUse = self.container.outputData.DICTOUT
            self.dictDone()
        except Exception as e:
            self.appendErrorReport(203, 'Exception in lidiaAcedrg_finished: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            
    def dictDone(self):
        if self.container.controlParameters.OBSAS.__str__() == 'UNMERGED':
            self.aimlessPipe()
        else:
            self.obsToUse = self.container.inputData.F_SIGF_IN
            self.freerToUse = self.container.inputData.FREERFLAG_IN
            self.rigidBodyPipeline()

    def aimlessPipe(self):
        try:
            self.aimlessPlugin = self.makePluginObject('aimless_pipe')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create aimless_pipe plugin: {e}')
            self.reportStatus(CPluginScript.FAILED)
            return
        self.aimlessPlugin.container.controlParameters.MODE.set('MATCH')
        self.aimlessPlugin.container.controlParameters.RESOLUTION_RANGE = self.container.controlParameters.RESOLUTION_RANGE
        self.aimlessPlugin.container.controlParameters.SCALING_PROTOCOL.set('DEFAULT')
        self.aimlessPlugin.container.controlParameters.ONLYMERGE.set(False)
        self.aimlessPlugin.container.controlParameters.REFERENCE_DATASET.set('XYZ')
        self.aimlessPlugin.container.controlParameters.AUTOCUTOFF.set(True)
        self.aimlessPlugin.container.inputData.copyData(self.container.inputData,['UNMERGEDFILES'])
        self.aimlessPlugin.container.inputData.XYZIN_REF = self.container.inputData.XYZIN
        self.aimlessPlugin.container.controlParameters.TOLERANCE.set(10.)
        if self.container.inputData.FREERFLAG_IN.isSet():
            self.aimlessPlugin.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG_IN
        self.connectSignal(self.aimlessPlugin,'finished',self.aimlessPlugin_finished)
        self.aimlessPlugin.process()

    @QtCore.Slot(dict)
    def aimlessPlugin_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.appendErrorReport(204, 'Aimless pipeline failed')
            self.reportStatus(CPluginScript.FAILED)
            return

        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.aimlessPlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()

            # Check that aimless produced expected outputs before harvesting
            aimlessOut = self.aimlessPlugin.container.outputData
            if not aimlessOut.FREEROUT.isSet() or not os.path.isfile(str(aimlessOut.FREEROUT.fullPath)):
                self.appendErrorReport(204, 'Aimless did not produce FreeR output')
                self.reportStatus(CPluginScript.FAILED)
                return
            if len(aimlessOut.HKLOUT) == 0 or not os.path.isfile(str(aimlessOut.HKLOUT[0].fullPath)):
                self.appendErrorReport(204, 'Aimless did not produce merged HKL output')
                self.reportStatus(CPluginScript.FAILED)
                return

            self.harvestFile(aimlessOut.FREEROUT, self.container.outputData.FREERFLAG_OUT)
            self.harvestFile(aimlessOut.HKLOUT[0], self.container.outputData.F_SIGF_OUT)
            self.obsToUse = self.container.outputData.F_SIGF_OUT
            self.freerToUse = self.container.outputData.FREERFLAG_OUT
            self.rigidBodyPipeline()
        except Exception as e:
            self.appendErrorReport(204, 'Exception in aimlessPlugin_finished: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)

    def rigidBodyPipeline(self):
        inp = self.container.inputData
        if hasattr(inp,'PIPELINE') and inp.PIPELINE.isSet() and inp.PIPELINE.__str__() == 'DIMPLE':
            self.i2Dimple()
        else:
            self.phaser_rnp_pipeline()

    def phaser_rnp_pipeline(self):
        try:
            self.rnpPlugin = self.makePluginObject('phaser_rnp_pipeline')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create phaser_rnp_pipeline plugin: {e}')
            self.reportStatus(CPluginScript.FAILED)
            return
        self.rnpPlugin.container.inputData.XYZIN_PARENT = self.selAtomsFile
        self.rnpPlugin.container.inputData.F_SIGF = self.obsToUse
        self.rnpPlugin.container.inputData.FREERFLAG = self.freerToUse
        self.rnpPlugin.container.inputData.SELECTIONS.append({'text':'/*/*/*/*','pdbFileKey':'XYZIN_PARENT'})
        self.connectSignal(self.rnpPlugin,'finished',self.rnpPlugin_finished)
        self.rnpPlugin.process()

    @QtCore.Slot(dict)
    def rnpPlugin_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.appendErrorReport(205, 'phaser_rnp_pipeline failed')
            self.reportStatus(CPluginScript.FAILED)
            return
        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.rnpPlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()
            self.harvestFile(self.rnpPlugin.container.outputData.MAPOUT_REFMAC, self.container.outputData.FPHIOUT)
            self.mapToUse = self.container.outputData.FPHIOUT
            self.harvestFile(self.rnpPlugin.container.outputData.DIFMAPOUT_REFMAC, self.container.outputData.DIFFPHIOUT)
            if os.path.isfile(str(self.rnpPlugin.container.outputData.FREERFLAG_OUT)):
                self.harvestFile(self.rnpPlugin.container.outputData.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
                self.freerToUse = self.container.outputData.FREERFLAG_OUT
            if os.path.isfile(str(self.rnpPlugin.container.outputData.F_SIGF_OUT)):
                self.harvestFile(self.rnpPlugin.container.outputData.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
                self.obsToUse = self.container.outputData.F_SIGF_OUT
            if self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
                self.harvestFile(self.rnpPlugin.container.outputData.XYZOUT[0], self.container.outputData.XYZOUT)
                self.reportStatus(CPluginScript.SUCCEEDED)
            else:
                self.coordinatesForCoot = self.rnpPlugin.container.outputData.XYZOUT_REFMAC
                self.cootAddLigand()
        except Exception as e:
            self.appendErrorReport(205, 'Exception in rnpPlugin_finished: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
        
    def i2Dimple(self):
        try:
            self.i2DimplePlugin = self.makePluginObject('i2Dimple')
        except Exception as e:
            self.appendErrorReport(209, f'Failed to create i2Dimple plugin: {e}')
            self.reportStatus(CPluginScript.FAILED)
            return
        self.i2DimplePlugin.container.inputData.XYZIN = self.selAtomsFile
        self.i2DimplePlugin.container.inputData.F_SIGF = self.obsToUse
        self.i2DimplePlugin.container.inputData.FREERFLAG = self.freerToUse
        self.connectSignal(self.i2DimplePlugin,'finished',self.i2Dimple_finished)
        self.i2DimplePlugin.process()

    @QtCore.Slot(dict)
    def i2Dimple_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.appendErrorReport(206, 'i2Dimple pipeline failed')
            self.reportStatus(CPluginScript.FAILED)
            return
        try:
            pluginRoot = CCP4Utils.openFileToEtree(self.i2DimplePlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
            self.flushXML()
            self.harvestFile(self.i2DimplePlugin.container.outputData.FPHIOUT, self.container.outputData.FPHIOUT)
            self.mapToUse = self.container.outputData.FPHIOUT
            self.harvestFile(self.i2DimplePlugin.container.outputData.DIFFPHIOUT, self.container.outputData.DIFFPHIOUT)
            #Adopt the reindexed output of the dimple file if present
            if os.path.isfile(self.i2DimplePlugin.container.outputData.F_SIGF_OUT.fullPath.__str__()):
                self.harvestFile(self.i2DimplePlugin.container.outputData.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
            if os.path.isfile(self.i2DimplePlugin.container.outputData.FREERFLAG_OUT.fullPath.__str__()):
                self.harvestFile(self.i2DimplePlugin.container.outputData.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
            if self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
                self.harvestFile(self.i2DimplePlugin.container.outputData.XYZOUT, self.container.outputData.XYZOUT)
                self.reportStatus(CPluginScript.SUCCEEDED)
            else:
                self.coordinatesForCoot = self.i2DimplePlugin.container.outputData.XYZOUT
                self.cootAddLigand()
        except Exception as e:
            self.appendErrorReport(206, 'Exception in i2Dimple_finished: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
        
    def cootAddLigand(self):
        """Fit ligand into density using coot_headless_api."""
        try:
            # Get file paths
            xyzin = str(self.coordinatesForCoot.fullPath)
            mtzin = str(self.mapToUse.fullPath)
            dictin = str(self.dictToUse.fullPath)
            xyzout = str(self.container.outputData.XYZOUT.fullPath)

            print(f"[COOT DEBUG] cootAddLigand starting...")
            print(f"[COOT DEBUG] xyzin = {xyzin}")
            print(f"[COOT DEBUG] mtzin = {mtzin}")
            print(f"[COOT DEBUG] dictin = {dictin}")
            print(f"[COOT DEBUG] xyzout = {xyzout}")

            # Initialize coot headless API
            print("[COOT DEBUG] Initializing coot_headless_api...")
            mc = coot_headless_api.molecules_container_py(True)
            mc.set_make_backups(False)
            mc.set_use_gemmi(False)
            print("[COOT DEBUG] molecules_container_py initialized")

            # Load protein coordinates
            print(f"[COOT DEBUG] Loading protein coordinates from {xyzin}...")
            imol_protein = mc.read_pdb(xyzin)
            print(f"[COOT DEBUG] imol_protein = {imol_protein}")
            if imol_protein < 0:
                self.appendErrorReport(207, f'Failed to read protein coordinates: {xyzin}')
                self.reportStatus(CPluginScript.FAILED)
                return

            # Load map from MTZ
            # CCP4i2 minimtz format uses F/PHI columns
            print(f"[COOT DEBUG] Loading map from {mtzin}...")

            # Try CCP4i2 minimtz column names (F/PHI)
            imol_map = mc.read_mtz(mtzin, "F", "PHI", "", False, False)
            print(f"[COOT DEBUG] F/PHI attempt: imol_map = {imol_map}")
            if imol_map < 0:
                # Fallback to dimple raw output (2FOFCWT/PH2FOFCWT)
                print("[COOT DEBUG] F/PHI columns failed, trying dimple columns (2FOFCWT/PH2FOFCWT)...")
                imol_map = mc.read_mtz(mtzin, "2FOFCWT", "PH2FOFCWT", "", False, False)
            if imol_map < 0:
                # Fallback to REFMAC columns (FWT/PHWT)
                print("[COOT DEBUG] Dimple columns failed, trying REFMAC columns (FWT/PHWT)...")
                imol_map = mc.read_mtz(mtzin, "FWT", "PHWT", "", False, False)

            print(f"[COOT DEBUG] Final imol_map = {imol_map}")
            if imol_map < 0:
                self.appendErrorReport(207, f'Failed to read map coefficients: {mtzin}')
                self.reportStatus(CPluginScript.FAILED)
                return

            # Set the refinement map
            print("[COOT DEBUG] Setting refinement map...")
            mc.set_imol_refinement_map(imol_map)

            # Import ligand dictionary (use -999999 for IMOL_ENC_ANY)
            print(f"[COOT DEBUG] Importing ligand dictionary from {dictin}...")
            dict_result = mc.import_cif_dictionary(dictin, -999999)
            print(f"[COOT DEBUG] dict_result = {dict_result}")
            if dict_result == 0:
                self.appendErrorReport(207, f'Failed to import ligand dictionary: {dictin}')
                self.reportStatus(CPluginScript.FAILED)
                return

            # Get a monomer for fitting (using 'DRG' as the TLC)
            print("[COOT DEBUG] Getting monomer 'DRG'...")
            imol_ligand = mc.get_monomer('DRG')
            print(f"[COOT DEBUG] imol_ligand = {imol_ligand}")
            if imol_ligand < 0:
                self.appendErrorReport(207, 'Failed to get monomer DRG from dictionary')
                self.reportStatus(CPluginScript.FAILED)
                return

            # Fit ligand into density
            print("[COOT DEBUG] Running fit_ligand...")
            ligands_found = mc.fit_ligand(imol_protein, imol_map, imol_ligand, 1.5, True, 20)
            print(f"[COOT DEBUG] ligands_found = {ligands_found}")

            if ligands_found and len(ligands_found) > 0:
                # fit_ligand returns fit_ligand_info_t objects with imol attribute
                # Log details of the first result for debugging
                first_lig = ligands_found[0]
                print(f"[COOT DEBUG] First ligand info: imol={first_lig.imol}, "
                      f"cluster_idx={first_lig.cluster_idx}, ligand_idx={first_lig.ligand_idx}")

                # Determine how many ligands to merge based on NCS
                n_to_copy = 1  # Default to merging just the best ligand
                try:
                    print("[COOT DEBUG] Checking NCS chains...")
                    ncs_chains = mc.get_ncs_related_chains(imol_protein)
                    print(f"[COOT DEBUG] ncs_chains = {ncs_chains}")
                    if ncs_chains and len(ncs_chains) > 0:
                        # ncs_chains is a nested list like [[A,C], [B,D]]
                        # Use the number of chains in the first group as ligand count
                        n_to_copy = min(len(ligands_found), len(ncs_chains[0]))
                except (AttributeError, Exception) as e:
                    # Method may not be available in older coot versions
                    print(f"[COOT DEBUG] NCS check not available: {e}")
                    n_to_copy = 1

                print(f"[COOT DEBUG] Found {len(ligands_found)} ligand positions, merging {n_to_copy}")

                # Merge found ligands into protein model
                # Extract imol from fit_ligand_info_t objects
                if n_to_copy > 0:
                    ligand_imols = [lig.imol for lig in ligands_found[:n_to_copy]]
                    ligand_indices = ','.join(str(imol) for imol in ligand_imols)
                    print(f"[COOT DEBUG] Merging ligand imols: {ligand_indices}")
                    mc.merge_molecules(imol_protein, ligand_indices)
            else:
                print("[COOT DEBUG] No ligands found by fit_ligand")

            # Write output coordinates
            print(f"[COOT DEBUG] Writing coordinates to {xyzout}...")
            mc.write_coordinates(imol_protein, xyzout)

            # Clean up coot backup directory
            shutil.rmtree("coot-backup", ignore_errors=True)

            # Check output was created
            if not os.path.isfile(xyzout):
                self.appendErrorReport(207, 'Coot did not produce output coordinates')
                self.reportStatus(CPluginScript.FAILED)
                return

            print("[COOT DEBUG] Output file created successfully")

            # Postprocessing: update XML with model composition
            self._updateModelComposition()

            print("[COOT DEBUG] cootAddLigand completed successfully")
            self.finishWithStatus(CPluginScript.SUCCEEDED)

        except Exception as e:
            import traceback
            print(f"[COOT DEBUG] Exception in cootAddLigand: {e}")
            print(f"[COOT DEBUG] Traceback: {traceback.format_exc()}")
            self.appendErrorReport(207, 'Exception in cootAddLigand: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)

    def _updateModelComposition(self):
        """Update XML with monomer composition from output coordinates."""
        try:
            if os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)):
                from ccp4i2.core.CCP4ModelData import CPdbData
                aCPdbData = CPdbData()
                aCPdbData.loadFile(self.container.outputData.XYZOUT.fullPath)
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

            self.appendErrorReport(208, 'Exception in model composition update: ' + str(e))

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
        try:
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            pipelineOutputItem.annotation = pluginOutputItem.annotation
            pipelineOutputItem.contentFlag = pluginOutputItem.contentFlag
            pipelineOutputItem.subType = pluginOutputItem.subType
        except Exception as e:
            self.appendErrorReport(202, 'Failed to harvest file: ' + str(pluginOutputItem.fullPath) + ' -> ' + str(pipelineOutputItem.fullPath) + ': ' + str(e))
            self.finishWithStatus(CPluginScript.FAILED)

    def appendXML(self, changedFile, replacingElementOfType=None):
        newXML = CCP4Utils.openFileToEtree(changedFile)
        oldNodes = self.xmlroot.xpath(replacingElementOfType)
        if len(oldNodes) > 0: oldNodes[0].parent().remove(oldNodes[0])
        self.xmlroot.append(newXML)
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))

    def checkFinishStatus( self, statusDict,failedErrCode,outputFile = None,noFileErrCode= None):
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        try:
            assert outputFile.exists(),'Entity provided is not CDataFile or does not exist'
        except Exception as e:
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile) + ': ' + str(e))
            self.finishWithStatus(CPluginScript.FAILED)

    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)

    def flushXML(self):
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))
