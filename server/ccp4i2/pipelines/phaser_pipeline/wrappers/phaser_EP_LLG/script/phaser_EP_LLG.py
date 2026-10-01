from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script import phaser_EP_AUTO

class phaser_EP_LLG(phaser_EP_AUTO.phaser_EP_AUTO):

    TASKNAME = 'phaser_EP_LLG'
    WHATNEXT = ['coot_rebuild',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]

    ERROR_CODES = { 105 : { 'description' : 'Phaser stopped with an error' },
                    201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}

    def processOutputFiles(self):
        status = super(phaser_EP_LLG, self).processOutputFiles()
        out = self.container.outputData
        # Phased by a model, Phaser places no sites: its "sites" file is a
        # header with no atoms, and must not be offered as a structure.
        for item in list(out.XYZOUT):
            path = item.getFullPath() if item.isSet() else None
            if not path or not _has_atoms(path):
                out.XYZOUT.remove(item)
        for item in out.HKLOUT:
            item.annotation.set('Anomalous LLG map, all of Phaser\'s columns')
        # One hand only (the model fixes it): "- original hand" says nothing.
        for name in ('LLGMAPOUT', 'MAPOUT', 'ABCDOUT'):
            items = getattr(out, name)
            if len(items) == 1:
                text = str(items[0].annotation)
                if text.endswith(' - original hand'):
                    items[0].annotation.set(text[:-len(' - original hand')])
        return status

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


def _has_atoms(path):
    try:
        with open(path) as f:
            return any(line.startswith(('ATOM', 'HETATM')) for line in f)
    except OSError:
        return False
