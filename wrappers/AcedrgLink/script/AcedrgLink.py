from __future__ import print_function
"""
    AcedrgLink.py: CCP4 GUI Project
    
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
import base64
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
#from lxml import etree
from xml.etree import ElementTree as ET

class AcedrgLink(CPluginScript):
    TASKNAME = 'AcedrgLink'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="acedrg"
    
    def __init__(self, *args, **kws):
        super(AcedrgLink, self).__init__(*args, **kws)

    def processInputFiles(self):
        print('AceDRG in link mode - processing input')

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

    def makeCommandAndScript(self,**kw):
        self.appendCommandLine(['-L',self.container.inputData.INSTRUCTION_FILE.fullPath])
        self.appendCommandLine(['-o',os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__())])
        if self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS.isSet():
           for kwLine in str(self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS).split('\n'):
              kw = kwLine.lstrip().rstrip()
              #print 'kw','['+str(kw)+']'
              if len(kw)>0:
                 if str(kw)[0] != '#':
                    self.appendCommandLine(kw)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print('AceDRG in link mode - processing output')

        self.container.outputData.CIF_OUT.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_link.cif"))
        if self.container.inputData.ANNOTATION.isSet():
           self.container.outputData.CIF_OUT.annotation.set('Link dictionary: '+self.container.inputData.ANNOTATION.__str__())
        else:
           self.container.outputData.CIF_OUT.annotation.set('Link dictionary: '+self.container.inputData.LINK_ID.__str__())

        #if not os.path.isfile(str(self.container.outputData.CIF_OUT.fullPath)):
        #self.appendErrorReport("AceDRG did not successfully create link CIF dictionary.")
        #return CPluginScript.FAILED

        self.container.outputData.UNL_PDB.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_TMP","UNL_for_link.pdb"))
        self.container.outputData.UNL_CIF.setFullPath(os.path.join(self.getWorkDirectory(),self.container.inputData.LINK_ID.__str__()+"_TMP","UNL_for_link.cif"))


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
        
        import sys
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = ET.Element("acedrg_link")
            logText = ET.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                #logText.text = ET.CDATA(logFile.read())
                logText.text = base64.b64encode(logFile.read())
            CCP4Utils.writeXML(programXMLFile,ET.tostring(xmlStructure))

        return CPluginScript.SUCCEEDED
