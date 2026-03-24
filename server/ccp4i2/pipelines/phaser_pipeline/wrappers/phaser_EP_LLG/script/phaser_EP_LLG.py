from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script import phaser_EP_AUTO

class phaser_EP_LLG(phaser_EP_AUTO.phaser_EP_AUTO):

    TASKNAME = 'phaser_EP_LLG'
    WHATNEXT = ['coot_rebuild',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}
    requiredDefaultList = ['PART_VARI', 'PART_DEVI', 'LLGM']

    def validity(self):
        error = super(phaser_EP_LLG, self).validity()
        xyzin_partial = getattr(self.container.inputData, 'XYZIN_PARTIAL', None)
        if xyzin_partial is not None and xyzin_partial.isSet():
            cf = getattr(xyzin_partial, 'contentFlag', None)
            if cf == 2:  # CONTENT_FLAG_MMCIF — Phaser only works with PDB
                error.append(
                    klass=self.TASKNAME, code=200,
                    details='Phaser apps can only work with PDB format',
                    name=f'{self.TASKNAME}.container.inputData.XYZIN_PARTIAL',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR,
                )
        return error
