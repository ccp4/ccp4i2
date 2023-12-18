"""
    comit.py: CCP4 GUI Project
    
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
from core import CCP4Utils
import base64
#from lxml import etree
from xml.etree import ElementTree as ET

class comit(CPluginScript):
    TASKNAME = 'comit'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kevin.cowtan@york.ac.uk'
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ], ['log_mtzjoin.txt', 0] ]
    TASKCOMMAND="comit"
    
    def __init__(self, *args, **kws):
        super(comit, self).__init__(*args, **kws)

    def processInputFiles(self):
        from core import CCP4XtalData
        colgrps = [ ['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'F_PHI_IN' ] 
        self.hklin, columns, error = self.makeHklin0( colgrps )
        from core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        inp = self.container.inputData
        out = self.container.outputData
        con = self.container.controlParameters

        import os
        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")

        self.appendCommandLine([ '-stdin' ])
        self.appendCommandScript( 'mtzin ' + self.hklin )
        self.appendCommandScript( 'mtzout ' + self.hklout )
        self.appendCommandScript( 'colout i2' )
        self.appendCommandScript( "colin-fo F_SIGF_F,F_SIGF_SIGF" )
        self.appendCommandScript( "colin-fc F_PHI_IN_F,F_PHI_IN_PHI" )
        self.appendCommandScript( "nomit %s"%(str(self.container.controlParameters.NOMIT)) )
        self.appendCommandScript( "pad-radius %s"%(str(self.container.controlParameters.PAD_RADIUS)) )

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        # Split an MTZ file into minimtz data objects
        infile = os.path.join(self.workDirectory,'hklout.mtz')
        error = self.splitHklout( ['F_PHI_OUT'],['i2.F_phi.F,i2.F_phi.phi'], infile=infile )
        self.container.outputData.F_PHI_OUT.subType = 1
        from core import CCP4ErrorHandling
        if error.maxSeverity ( ) > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
        
        import sys
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = ET.Element("i2comit")
            logText = ET.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                #logText.text = ET.CDATA(logFile.read())
                logText.text = base64.b64encode(logFile.read())
            CCP4Utils.writeXML(programXMLFile,ET.tostring(xmlStructure))
        return CPluginScript.SUCCEEDED
