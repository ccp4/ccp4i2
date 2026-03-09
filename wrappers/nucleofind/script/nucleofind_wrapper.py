from __future__ import print_function
from core.CCP4PluginScript import CPluginScript
from core.CCP4Modules import PROCESSMANAGER
from core.CCP4XtalData import CMapDataFile
from lxml import etree

from core import CCP4ErrorHandling, CCP4Utils
import os, shutil

class nucleofind(CPluginScript):

    TASKMODULE          = 'model_building' # Where this plugin will appear on the gui
    TASKCOMMAND         = 'nucleofind'  # The command to execute, should be reachable"
    WHATNEXT = [ 'coot_rebuild' ]
    MAINTAINER = 'jordan.dialpuri@york.ac.uk'
    TASKVERSION= 0.1          # Version of this plugin"


    def processInputFiles(self):
      if self.container.inputData.FWT_PHWT.isSet():
          self.hklin, self.columns, error = self.makeHklin0([ 'FWT_PHWT' ])
          if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
              return CPluginScript.FAILED
      return CPluginScript.SUCCEEDED


    def makeCommandAndScript(self):
        self.appendCommandLine(["-i", self.hklin])
        self.appendCommandLine(["-o", "nucleofind"])
        self.appendCommandLine(["-m", "core"])
        self.appendCommandLine(["-n", self.container.controlParameters.THREADS])
        self.appendCommandLine(["--overlap", self.container.controlParameters.OVERLAP])

        if self.container.controlParameters.SYMMETRY.isSet():
            self.appendCommandLine(["--use-symmetry"])
            
        if self.container.controlParameters.RESOLUTION.isSet():
            self.appendCommandLine(["--resolution", self.container.controlParameters.RESOLUTION])
          
        if self.container.inputData.FWT_PHWT.isSet() :
          # f, p = self.columns
          print(f"{self.columns}")  
          self.appendCommandLine(["--amplitude", "FWT_PHWT_F"])
          self.appendCommandLine(["--phase", "FWT_PHWT_PHI"])

        return CPluginScript.SUCCEEDED
      
    def processOutputFiles(self):
        print("In process output files ")
        wd = self.getWorkDirectory()
        directory = os.path.join(wd, "nucleofind")
        phosphate_map = os.path.join(directory, "nucleofind-phosphate.map")
        sugar_map = os.path.join(directory, "nucleofind-sugar.map")
        base_map = os.path.join(directory, "nucleofind-base.map")
        outputData = self.container.outputData
        
        outputData.PHOSPHATE_MAP.setFullPath(os.path.join(wd, "nucleofind-phosphate.map"))
        shutil.copy(phosphate_map, outputData.PHOSPHATE_MAP.fullPath.__str__())
        outputData.PHOSPHATE_MAP.annotation.set("NucleoFind Predicted Phosphate Map")

        outputData.SUGAR_MAP.setFullPath(os.path.join(wd, "nucleofind-sugar.map"))
        shutil.copy(sugar_map, outputData.SUGAR_MAP.fullPath.__str__())
        outputData.SUGAR_MAP.annotation.set("NucleoFind Predicted Sugar Map")
        
        outputData.BASE_MAP.setFullPath(os.path.join(wd, "nucleofind-base.map"))
        shutil.copy(base_map, outputData.BASE_MAP.fullPath.__str__())
        outputData.BASE_MAP.annotation.set("NucleoFind Predicted Base Map")

        xmlRoot = etree.Element('NucleoFind')
        with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
            xmlString = etree.tostring ( xmlRoot, pretty_print=True )
            CCP4Utils.writeXML(xmlFile,xmlString)
        return CPluginScript.SUCCEEDED
 