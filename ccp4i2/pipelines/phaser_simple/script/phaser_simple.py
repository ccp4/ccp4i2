try:
    import ccp4mg
    import mmdb2 as mmdb
except:
    print('FAILED CCP4ModelData imported ccp4mg')
import mmut

from pipelines.phaser_pipeline.script import phaser_pipeline
  
class phaser_simple(phaser_pipeline.phaser_pipeline):
    TASKNAME = 'phaser_simple'                                  # Task name - should be same as class name
    
    def process(self):
        self.createEnsembleElements()
        super(phaser_simple,self).process()
        
    def checkInputData(self):
        invalidFiles = super(phaser_simple,self).checkInputData()
        if (not self.container.inputData.INPUT_FIXED) and ('XYZIN_FIXED' in invalidFiles):
            invalidFiles.remove('XYZIN_FIXED')
        return invalidFiles

    def createEnsembleElements(self):
        from core.CCP4ModelData import CPdbDataFile, CAtomSelection, CPdbEnsembleItem
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
        pdbItem = elements[-1]
        pdbItem.structure.set(self.container.inputData.XYZIN)
        if self.container.inputData.ID_RMS == 'ID':
            pdbItem.identity_to_target.set(self.container.inputData.SEARCHSEQUENCEIDENTITY)
        elif self.container.inputData.ID_RMS == 'RMS':
            pdbItem.rms_to_target.set(self.container.inputData.SEARCHRMS)
        else:
            pdbItem.identity_to_target.set(None)
            pdbItem.rms_to_target.set(None)
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
