from ccp4i2.core.CCP4PluginScript import CPluginScript


class ZZPluginNameZZ(CPluginScript):
    TASKNAME = 'ZZPluginNameZZ'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'ZZPluginMaintainerZZ'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="command"
    
    def __init__(self, *args, **kws):
        super(ZZPluginNameZZ, self).__init__(*args, **kws)

    def processInputFiles(self):
        #Preprocess reflections to generate an "HKLIN" file
        '''
        #makeHklin0 takes as arguments a list of sublists
        #Each sublist comprises 1) A reflection data object identifier (one of those specified in the inputData container 
        #                           the task in the corresponding .def.xml
        #                       2) The requested data representation type to be placed into the file that is generated
        #
        #makeHklin0 returns a tuple comprising:
        #                       1) the file path of the file that has been created
        #                       2) a list of strings, each of which contains a comma-separated list of column labels output from
        #                       the input data objects
        #                       3) A CCP4 Error object        
        from ccp4i2.core import CCP4XtalData
        self.hklin, self.columns, error = self.makeHklin0([
            ['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]
        ])
        self.columnsAsArray = self.columns.split(",")
        
        from ccp4i2.core import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        '''
        
        #Preprocess coordinates to extract a subset
        '''
        # The method "getSelectedAtomsPdbFile" applied to a coordinate data object
        # selects those atoms declared in the objects "selectionString" property and writes them into
        # a pruned down file, the name of which is provided in the argument
        self.selectedCoordinatesPath = os.path.join(self.getWorkDirectory(), "selected_xyzin.pdb")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedCoordinatesPath)
        '''
        
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        #self.appendCommandLine(self.getWorkDirectory())
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        #Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        '''
        self.container.outputData.XYZOUT.setFullPath(os.path.join(self.getWorkDirectory(),"final.pdb"))
        '''
        
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
