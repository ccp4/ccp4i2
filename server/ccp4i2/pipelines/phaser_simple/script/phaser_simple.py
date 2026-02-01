from ccp4i2.pipelines.phaser_pipeline.script import phaser_pipeline
from ccp4i2.core.CCP4PluginScript import CPluginScript

class phaser_simple(phaser_pipeline.phaser_pipeline):
    TASKNAME = 'phaser_simple'

    ERROR_CODES = {
        301: {'description': 'Exception in createEnsembleElements'},
        302: {'description': 'Exception setting up search model ensemble'},
        303: {'description': 'Exception setting up fixed structure ensemble'},
    }

    def process(self):
        self.createEnsembleElements()
        return super(phaser_simple,self).process()
        
    def checkInputData(self):
        invalidFiles = super(phaser_simple,self).checkInputData()
        if (not self.container.inputData.INPUT_FIXED) and ('XYZIN_FIXED' in invalidFiles):
            invalidFiles.remove('XYZIN_FIXED')
        return invalidFiles

    def validity(self):
        """Override to filter out ENSEMBLES list length error.

        ENSEMBLES is intentionally empty at validation time because it's
        populated programmatically by createEnsembleElements() during process().
        """
        from ccp4i2.core import CCP4ErrorHandling

        # Get parent validation
        error = super(phaser_simple, self).validity()

        # Filter out the ENSEMBLES minimum length error (code 101)
        # This error is expected since ENSEMBLES is populated in createEnsembleElements()
        filtered = CCP4ErrorHandling.CErrorReport()
        for err in error.getErrors():
            # Skip error code 101 (min list length) for ENSEMBLES
            if err.get('code') == 101 and 'ENSEMBLES' in err.get('name', ''):
                continue
            filtered.append(
                klass=err.get('class', ''),
                code=err.get('code', 0),
                details=err.get('details', ''),
                name=err.get('name', ''),
                severity=err.get('severity', 0)
            )

        return filtered

    def createEnsembleElements(self):
        try:
            from ccp4i2.core.CCP4ModelData import CPdbDataFile, CAtomSelection, CPdbEnsembleItem
            elements = self.container.inputData.ENSEMBLES
            #Before removing all elements from this list, I have to set its listMinLength to 0
            self.container.inputData.ENSEMBLES.setQualifiers({'listMinLength':0})
            while len(elements) > 0: elements.remove(elements[-1])
            self.container.inputData.ENSEMBLES.append(self.container.inputData.ENSEMBLES.makeItem())
            ensemble = self.container.inputData.ENSEMBLES[-1]
            ensemble.number.set(self.container.inputData.NCOPIES)
            ensemble.label.set('SearchModel')
            elements = ensemble.pdbItemList
            while len(elements) > 1: elements.remove(elements[-1])
            while len(elements) < 1: elements.append(elements.makeItem())
            pdbItem = elements[-1]
            pdbItem.structure.set(self.container.inputData.XYZIN)
            if self.container.inputData.ID_RMS == 'ID':
                pdbItem.identity_to_target.set(self.container.inputData.SEARCHSEQUENCEIDENTITY)
            elif self.container.inputData.ID_RMS == 'RMS':
                pdbItem.rms_to_target.set(self.container.inputData.SEARCHRMS)
            else:
                # Default to identity of 0.9 when ID_RMS is 'CARD' or not set
                # This prevents validity errors requiring either identity or RMS to be set
                pdbItem.identity_to_target.set(0.9)
        except Exception as e:
            self.appendErrorReport(302, 'Exception setting up search model ensemble: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            raise

        try:
            if self.container.inputData.INPUT_FIXED.isSet() and self.container.inputData.INPUT_FIXED and self.container.inputData.XYZIN_FIXED.isSet():
                self.container.inputData.ENSEMBLES.append(self.container.inputData.ENSEMBLES.makeItem())
                ensemble = self.container.inputData.ENSEMBLES[-1]
                ensemble.number.set(0)
                ensemble.label.set('KnownStructure')
                elements = ensemble.pdbItemList
                while len(elements) > 1: elements.remove(elements[-1])
                pdbItem = elements[-1]
                pdbItem.structure.set(self.container.inputData.XYZIN_FIXED)
                if self.container.inputData.FIXED_ID_RMS == 'ID':
                    pdbItem.identity_to_target.set(self.container.inputData.FIXEDSEQUENCEIDENTITY)

                elif self.container.inputData.FIXED_ID_RMS == 'RMS':
                    pdbItem.rms_to_target.set(self.container.inputData.FIXEDRMS)

                self.container.inputData.FIXENSEMBLES.append(self.container.inputData.FIXENSEMBLES.makeItem())
                self.container.inputData.FIXENSEMBLES[-1].set('KnownStructure')
        except Exception as e:
            self.appendErrorReport(303, 'Exception setting up fixed structure ensemble: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            raise

