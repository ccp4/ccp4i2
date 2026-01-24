from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script import phaser_EP_AUTO

class phaser_EP_LLG(phaser_EP_AUTO.phaser_EP_AUTO):

    TASKNAME = 'phaser_EP_LLG'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    RUNEXTERNALPROCESS=False
    WHATNEXT = ['coot_rebuild',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}
    requiredDefaultList = ['PART_VARI', 'PART_DEVI', 'LLGM']
