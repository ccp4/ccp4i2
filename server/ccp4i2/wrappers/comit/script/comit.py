import os

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class comit(CPluginScript):
    TASKNAME = 'comit'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kathryn.cowtan@york.ac.uk'
    TASKCOMMAND="comit"
    
    def __init__(self, *args, **kws):
        super(comit, self).__init__(*args, **kws)

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData
        colgrps = [ ['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'F_PHI_IN' ] 
        self.hklin, columns, error = self.makeHklin0( colgrps )
        from ccp4i2.core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
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
        from ccp4i2.core import CCP4ErrorHandling
        if error.maxSeverity ( ) > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
        import sys

        from lxml import etree
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("i2comit")
            logText = etree.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                logText.text = etree.CDATA(logFile.read())
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))
        return CPluginScript.SUCCEEDED
