"""
    Platonyzer.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

import os
from core.CCP4PluginScript import CPluginScript
import base64

class Platonyzer(CPluginScript):
    TASKNAME = 'Platonyzer'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="platonyzer"
    
    def __init__(self, *args, **kws):
        super(Platonyzer, self).__init__(*args, **kws)

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        if self.container.controlParameters.MODE.__str__() == 'NA_MG':
           self.appendCommandLine(['--create-na-mg-links'])
           if self.container.controlParameters.RM_VDW:
              self.appendCommandLine(['--delete-vdw-rest'])
        self.appendCommandLine([self.container.inputData.XYZIN.__str__()])
        self.appendCommandLine([self.container.outputData.XYZOUT.__str__()])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.container.outputData.RESTRAINTS.setFullPath(os.path.join(self.getWorkDirectory(),"XYZOUT.restraints"))
        if not os.path.exists(self.container.outputData.RESTRAINTS.__str__()):
           return CPluginScript.FAILED
        
        # Split an MTZ file into minimtz data objects
        '''
        outputFilesToMake = ['FPHIOUT','DIFFPHIOUT']
        columnsToTake = ['FWT,PHWT','DELFWT,PHDELWT']
        infile = os.path.join(self.workDirectory,'final.mtz')
        error = self.splitHklout(outputFilesToMake, columnsToTake, infile=infile)
        from core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        '''
        
        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
        '''
        #from lxml import etree
        from xml.etree import ElementTree as ET
        import sys
        import base64
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = ET.Element("i2Dimple")
            logText = ET.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                logText.text = base64.b64encode(logFile.read())
            programXMLFile.write(ET.tostring(xmlStructure))
        '''
        return CPluginScript.SUCCEEDED
