from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline.script import phaser_pipeline


class phaser_rnp_pipeline(phaser_pipeline.phaser_pipeline):

    TASKNAME = 'phaser_rnp_pipeline'

    ERROR_CODES = {
        200: {'description': 'Phaser exited with error status'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'Columns not present'},
        204: {'description': 'Failed in plugin'},
        205: {'description': 'Failed in pointless reindex operation'},
        206: {'description': 'Phaser produced no output coordinates'},
        210: {'description': 'Selection yields no atoms'},
        212: {'description': 'Overlapping atom selections'},
    }
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validity(self):
        """Validate inputs, suppressing the inherited ENSEMBLES list-length
        error (ENSEMBLES is built programmatically in createEnsembleElements).

        Additionally checks that every selection in SELECTIONS resolves to a
        non-empty, non-overlapping set of atoms against XYZIN_PARENT.
        """
        from ccp4i2.core import CCP4ErrorHandling

        error = super().validity()

        # Filter out the ENSEMBLES minimum-length error — it is populated
        # at runtime in createEnsembleElements(), so an empty list is expected.
        filtered = CCP4ErrorHandling.CErrorReport()
        for err in error.getErrors():
            if err.get('code') == 101 and 'ENSEMBLES' in err.get('name', ''):
                continue
            filtered.append(
                err.get('class', ''),
                err.get('code', 0),
                err.get('details', ''),
                err.get('name', ''),
                err.get('severity', 0),
            )

        # If XYZIN_PARENT is not set we cannot validate selections
        xyzin = self.container.inputData.XYZIN_PARENT
        if not xyzin.isSet():
            return filtered

        selections = self.container.inputData.SELECTIONS
        if len(selections) == 0:
            return filtered

        # Load the coordinate file and evaluate each selection
        try:
            content = xyzin.fileContent
            if content is None:
                return filtered
        except Exception:
            return filtered

        atom_sets = []
        for idx, sel in enumerate(selections):
            sel_text = str(sel.text.value) if hasattr(sel.text, 'value') else str(sel.text)
            try:
                n_atoms, selected = content.interpretSelection(sel_text)
            except Exception:
                n_atoms, selected = 0, []

            if n_atoms == 0:
                filtered.append(
                    'CAtomSelection',
                    210,
                    f'Selection {idx + 1} ("{sel_text}") matches no atoms '
                    f'in the parent coordinate file',
                    f'{self.TASKNAME}.container.inputData.SELECTIONS[{idx}]',
                    CCP4ErrorHandling.SEVERITY_ERROR,
                )
            else:
                # Build a hashable set of atom identities for overlap checking
                atom_ids = set()
                for _model, chain, residue, atom in selected:
                    atom_ids.add((
                        chain.name,
                        residue.seqid.num,
                        residue.seqid.icode,
                        residue.name,
                        atom.name,
                        atom.altloc,
                    ))
                atom_sets.append((idx, sel_text, atom_ids))

        # Check pairwise overlap
        for i in range(len(atom_sets)):
            for j in range(i + 1, len(atom_sets)):
                idx_i, text_i, set_i = atom_sets[i]
                idx_j, text_j, set_j = atom_sets[j]
                overlap = set_i & set_j
                if overlap:
                    filtered.append(
                        'CAtomSelection',
                        212,
                        f'Selections {idx_i + 1} ("{text_i}") and '
                        f'{idx_j + 1} ("{text_j}") share '
                        f'{len(overlap)} atom(s)',
                        f'{self.TASKNAME}.container.inputData.SELECTIONS',
                        CCP4ErrorHandling.SEVERITY_ERROR,
                    )

        return filtered

    # ------------------------------------------------------------------
    # Process
    # ------------------------------------------------------------------

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()

        self.xmlroot = etree.Element('PhaserPipeline')

        self.F_SIGF_TOUSE = self.container.inputData.F_SIGF
        self.FREERFLAG_TOUSE = self.container.inputData.FREERFLAG
        # Save flags before sub-plugins replace self.container.inputData
        self.runRefmacFlag = bool(self.container.inputData.RUNREFMAC)
        rv = self.runPointless()
        if rv != CPluginScript.SUCCEEDED:
            return rv
        rv = self.runPhaser(F_SIGF=self.F_SIGF_TOUSE)
        if rv != CPluginScript.SUCCEEDED:
            return rv

        if len(self.container.outputData.XYZOUT) == 0:
            self.appendErrorReport(206, 'No coordinates produced by Phaser')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        XYZIN_TOUSE = self.container.outputData.XYZOUT[0]

        if self.runRefmacFlag:
            self.runRefmac(F_SIGF=self.F_SIGF_TOUSE, FREERFLAG=self.FREERFLAG_TOUSE, XYZIN=XYZIN_TOUSE)

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    # ------------------------------------------------------------------
    # Ensemble construction
    # ------------------------------------------------------------------

    def createEnsembleElements(self):
        from ccp4i2.core.CCP4ModelData import CPdbDataFile, CPdbEnsembleItem
        elements = self.container.inputData.ENSEMBLES
        #Before removing all elements from this list, I have to set its listMinLength to 0
        self.container.inputData.ENSEMBLES.setQualifiers({'listMinLength':0})
        while len(elements) > 0: elements.remove(elements[-1])

        iEntity = 1
        if len(self.container.inputData.SELECTIONS) == 0:
            self.container.inputData.SELECTIONS.append(self.container.inputData.SELECTIONS.makeItem())
            self.container.inputData.SELECTIONS[-1].text.set("/*/*/*.*")
        for selection in self.container.inputData.SELECTIONS:
            self.container.inputData.ENSEMBLES.append(self.container.inputData.ENSEMBLES.makeItem())
            ensemble = self.container.inputData.ENSEMBLES[-1]
            ensemble.number.set(1)
            ensemble.label.set('Fragment_'+str(iEntity))
            elements = ensemble.pdbItemList
            while len(elements) > 1: elements.remove(elements[-1])
            if len(elements) == 0:
                elements.append(elements.makeItem())
            pdbItem = elements[-1]
            pdbItem.structure.set(self.container.inputData.XYZIN_PARENT)
            pdbItem.structure.selection.text.set(str(selection))
            pdbItem.identity_to_target.set(0.9)
            self.container.inputData.USINGSOLELEMENTS.append(self.container.inputData.USINGSOLELEMENTS.makeItem())
            self.container.inputData.USINGSOLELEMENTS[-1].set('Fragment_'+str(iEntity))
            iEntity += 1

    def runPhaser(self, F_SIGF=None):
        import traceback
        try:
            self.createEnsembleElements()
            phaserPlugin = self.makePluginObject('phaser_MR_RNP')
            #This funky arrangement is the way to ensure that the plugin behaves the same
            #when it is a part of the ipeline as it does when it is run alone...something about defaults I guess
            for attrName in phaserPlugin.container.keywords.dataOrder():
                if hasattr(self.container.keywords,attrName):
                    attr = getattr(self.container.keywords,attrName)
                    if hasattr(attr,'isSet') and attr.isSet():
                        setattr(phaserPlugin.container.keywords,attrName,attr)
            phaserPlugin.container.inputData=self.container.inputData
            phaserPlugin.container.inputData.RESOLUTION_LOW.set(25.0)
            phaserPlugin.container.inputData.RESOLUTION_HIGH.set(3.0)
            if F_SIGF is not None: phaserPlugin.container.inputData.F_SIGF=F_SIGF

            phaserPlugin.container.inputData.F_SIGF.loadFile()
            columns = phaserPlugin.container.inputData.F_SIGF.fileContent.getListOfColumns()
            columnStrings = [str(column) for column in columns]
            print(columnStrings)
            if 'F' in columnStrings: phaserPlugin.container.inputData.F_OR_I.set('F')

            rv = phaserPlugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(204,'phaser_MR_RNP')
                self.reportStatus(rv)
                return rv

            pluginOutputs=phaserPlugin.container.outputData
            pipelineOutputs = self.container.outputData
            self.appendXML(phaserPlugin.makeFileName('PROGRAMXML'),'PhaserMrResults')
            self.harvestFile(pluginOutputs.SOLOUT, pipelineOutputs.SOLOUT)
            for outputListType in ['XYZOUT', 'MAPOUT', 'DIFMAPOUT','PHASEOUT']:
                pluginOutputList = getattr(pluginOutputs, outputListType, None)
                pipelineOutputList = getattr(pipelineOutputs, outputListType, None)
                self.harvestList(pluginOutputList, pipelineOutputList)
        except Exception:
            traceback.print_exc()
            self.appendErrorReport(202,'phaser_MR_RNP')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED


    def runPointless(self):
        import traceback
        try:
            pointlessPlugin = self.makePluginObject('pointless_reindexToMatch')
            pointInp = pointlessPlugin.container.inputData
            pointlessPlugin.container.controlParameters.REFERENCE = 'XYZIN_REF'
            pointInp.XYZIN_REF = self.container.inputData.XYZIN_PARENT
            pointInp.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.FREERFLAG.isSet():
                pointInp.FREERFLAG = self.container.inputData.FREERFLAG
            rv = pointlessPlugin.process()
            if rv != CPluginScript.SUCCEEDED:
                # Log the sub-plugin's error report for diagnostics
                if hasattr(pointlessPlugin, 'errorReport'):
                    print(f"pointless_reindexToMatch error report: {pointlessPlugin.errorReport.report()}")
                self.appendErrorReport(205,'pointless_reindexToMatch')
                self.reportStatus(rv)
                return rv

            pluginOutputs = pointlessPlugin.container.outputData
            pipelineOutputs = self.container.outputData

            pluginOutputs.F_SIGF_OUT.loadFile()
            self.container.inputData.F_SIGF.loadFile()
            cellsAreSame = pluginOutputs.F_SIGF_OUT.fileContent.clipperSameCell(self.container.inputData.F_SIGF.fileContent)
            if not cellsAreSame['validity']:
                self.harvestFile(pluginOutputs.F_SIGF_OUT, pipelineOutputs.F_SIGF_OUT)
                self.F_SIGF_TOUSE = pluginOutputs.F_SIGF_OUT
                if self.container.inputData.FREERFLAG.isSet():
                    self.harvestFile(pluginOutputs.FREERFLAG_OUT, pipelineOutputs.FREERFLAG_OUT)
                    self.FREERFLAG_TOUSE = pluginOutputs.FREERFLAG_OUT

            try:
                self.appendXML(pointlessPlugin.makeFileName('PROGRAMXML'),'Pointless')
            except Exception:
                pass
        except Exception:
            traceback.print_exc()
            self.appendErrorReport(205,'pointless_reindexToMatch')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
