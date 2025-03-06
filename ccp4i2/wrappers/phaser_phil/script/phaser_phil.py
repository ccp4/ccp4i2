"""
    phaser_phil.py: CCP4 GUI Project

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

class phaser_phil(CPluginScript):
    TASKNAME = 'phaser_phil'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kyle.stevenson@stfc.ac.uk;david.waterman@stfc.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND='ccp4-python'

    def __init__(self, *args, **kws):
        super(phaser_phil, self).__init__(*args, **kws)

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
        from core import CCP4XtalData
        self.hklin, self.columns, error = self.makeHklin0([
            ['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]
        ])
        self.columnsAsArray = self.columns.split(",")

        from core import CCP4ErrorHandling
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

    # Copied here from xia2_dials - if useful should be moved to the base
    # class
    def extract_parameters(self, container):
        """Walk through a container locating parameters that have been set
        and return a list of name, value pairs"""

        result=[]
        dataOrder = container.dataOrder()
        contents = [getattr(container,name) for name in dataOrder]
        for model in contents:
            if isinstance(model, CCP4Container.CContainer):
                result.extend(self.extract_parameters(model))
            elif model.isSet():
                name = model.objectName().replace('__','.')
                # ensure commas are converted to whitespace-separated lists. Only
                # whitespace appears to work correctly with PHIL multiple
                # choice definitions.
                val = str(model.get()).split()
                val = ' '.join([v[:-1] if v.endswith(',') else v for v in val])
                result.append((name, val))
        return result

    def makeCommandAndScript(self,**kw):
        par = self.container.controlParameters
        inp = self.container.inputData

        # PHIL parameters set by the gui
        phil_file = os.path.normpath(os.path.join(
                    self.getWorkDirectory(), 'phaser.phil'))
        with open(phil_file, 'w') as f:
            for (name, val) in self.extract_parameters(par):
                #self.appendCommandLine()
                f.write(name + '={0}\n'.format(val))

        # FIXME it is possibly better, rather than writing out a file of
        # definition=value lines, to fetch this against the Phaser master phil
        # and hence construct a fully formatted phil definition for the job.

        # Get path to phaser's main.py
        from phaser.command_line import main
        self.appendCommandLine([main.__file__])

        self.appendCommandLine([phil_file])

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
        from core import CCP4ErrorHandling
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
