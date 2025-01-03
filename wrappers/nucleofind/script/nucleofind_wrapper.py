from __future__ import print_function
from core.CCP4PluginScript import CPluginScript
from core.CCP4Modules import PROCESSMANAGER
from core.CCP4XtalData import CMapDataFile

from core import CCP4ErrorHandling
import os

class nucleofind(CPluginScript):

    TASKMODULE          = 'model_building' # Where this plugin will appear on the gui
    TASKCOMMAND         = 'nucleofind'  # The command to execute, should be reachable"
    WHATNEXT = [ 'coot_rebuild' ]
    MAINTAINER = 'jordan.dialpuri@york.ac.uk'
    TASKVERSION= 0.1          # Version of this plugin"


    def processInputFiles(self):
      if self.container.inputData.FWT_PHWT.isSet():
          self.hklin, error = self.makeHklin([ 'FWT_PHWT' ])
          if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
              return CPluginScript.FAILED
      return CPluginScript.SUCCEEDED


    def makeCommandAndScript(self):
        print('Making command and script====================')
        self.appendCommandLine("-i "+ self.hklin)
        if self.container.inputData.FWT_PHWT.isSet( ) :
          self.appendCommandLine("-amplitude FWT_PHWT_F")
          self.appendCommandLine("-phase FWT_PHWT_PHI")

        self.appendCommandLine("-o nucleofind")
        self.appendCommandLine("-m core")
        print('Command and script made====================')
        return CPluginScript.SUCCEEDED
      
    def processOutputFiles(self):
        directory = os.path.join(self.getWorkDirectory(), "nucleofind")
        phosphate_map = os.path.join(directory, "nucleofind-phosphate.map")
        sugar_map = os.path.join(directory, "nucleofind-sugar.map")
        base_map = os.path.join(directory, "nucleofind-base.map")
        outputData = self.container.outputData
        outputData.PHOSPHATE_MAP.setFullPath(phosphate_map)
        outputData.PHOSPHATE_MAP.annotation.set("NucleoFind Predicted Phosphate Map")

        outputData.SUGAR_MAP.setFullPath(sugar_map)
        outputData.SUGAR_MAP.annotation.set("NucleoFind Predicted Sugar Map")
        
        
        outputData.BASE_MAP.setFullPath(base_map)
        outputData.BASE_MAP.annotation.set("NucleoFind Predicted Base Map")

        return CPluginScript.SUCCEEDED
