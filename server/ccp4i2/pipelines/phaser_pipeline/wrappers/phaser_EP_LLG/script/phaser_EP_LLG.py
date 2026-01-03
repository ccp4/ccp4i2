
from ccp4i2.core.CCP4PluginScript import CPluginScript
import sys, os
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Modules
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script import phaser_EP_AUTO
from lxml import etree

class phaser_EP_LLG(phaser_EP_AUTO.phaser_EP_AUTO):

    TASKNAME = 'phaser_EP_LLG'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS=False
    WHATNEXT = ['coot_rebuild',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}
    requiredDefaultList = ['PART_VARI', 'PART_DEVI', 'LLGM']

    def __init__(self,*args, **kwargs):
        super(phaser_EP_LLG,self).__init__(*args, **kwargs)
        
