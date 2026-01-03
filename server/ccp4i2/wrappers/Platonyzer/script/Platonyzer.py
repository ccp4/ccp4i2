import os

from ccp4i2.core.CCP4PluginScript import CPluginScript


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
        from ccp4i2.core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        '''
        
        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
        '''
        from lxml import etree
        import sys
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("i2Dimple")
            logText = etree.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                logText.text = etree.CDATA(logFile.read())
            programXMLFile.write(etree.tostring(xmlStructure))
        '''
        return CPluginScript.SUCCEEDED
