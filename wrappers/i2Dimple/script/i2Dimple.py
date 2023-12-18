"""
    i2Dimple.py: CCP4 GUI Project
    
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
from core import CCP4XtalData
from core import CCP4ErrorHandling
import platform
import base64
#from lxml import etree
from xml.etree import ElementTree as ET

class i2Dimple(CPluginScript):
    TASKNAME = 'i2Dimple'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'YOUR_EMAIL_ADDRESS'
    ERROR_CODES = { 201 : {'description' : 'Unable to extract information from dimple.log...run incomplete ?' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                        ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="dimple"
    if platform.system() == 'Windows': TASKCOMMAND = 'dimple.bat'
    PERFORMANCECLASS = 'CRefinementPerformance'


    def __init__(self, *args, **kws):
        super(i2Dimple, self).__init__(*args, **kws)

    def processInputFiles(self):
        inputs = [ ['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN] ]
        if self.container.inputData.FREERFLAG.isSet():
            inputs += [ ['FREERFLAG', 1] ]
        self.hklin,self.columns,error = self.makeHklin0(inputs)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        #makeHklin0 takes as arguments a list of sublists
        #Each sublist comprises 1) An input DATA object identifier as specified ni the inputData container of .de.f.xml
        #                       2) The requested data representation type to be placed into the file that is generated
        #makeHklin returns a tuple comprising:
        #                       1) the file path of the file that has been created
        #                       2) a list of strings each of which contains a comma-separated list of column labels output from
        #                       the input data objects that were input
        #                       3) A CCP4 Error object
        self.columnsAsArray = self.columns.split(",")
        
        import ccp4mg
        import mmdb2 as mmdb
        self.xyzin = os.path.join(self.getWorkDirectory(),"selected_xyzin.pdb")

        #A work around to cast input coordinates into PDB format for dimple if necessary
        testLoader = mmdb.Manager()
        result = testLoader.ReadCIFASCII(str(self.container.inputData.XYZIN.fullPath))
        if result == 0:
            temp_pdb = os.path.join(self.getWorkDirectory(),"selected_xyzin_temp.pdb")
            RC=testLoader.WritePDBASCII(self.xyzin)
            from core.CCP4ModelData import CPdbDataFile
            temporaryObject = CPdbDataFile(fullPath=self.xyzin)
            temporaryObject.loadFile()
            if self.container.inputData.XYZIN.isSelectionSet():
                temporaryObject.selection = self.container.inputData.XYZIN.selection
            temporaryObject.getSelectedAtomsPdbFile(temp_pdb)
            self.xyzin = temp_pdb
        else:
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.xyzin)

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        self.appendCommandLine(["-M",self.container.controlParameters.MR_WHEN_R.__str__()])
        self.appendCommandLine(["--fcolumn", self.columnsAsArray[0],
                                "--sigfcolumn", self.columnsAsArray[1]])
                                
        if self.container.inputData.FREERFLAG.isSet():
            self.appendCommandLine(["--free-r-flags", "-", "--freecolumn", self.columnsAsArray[2]])
        
        self.appendCommandLine([self.xyzin, self.hklin, self.getWorkDirectory()])
        
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        #Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        self.container.outputData.XYZOUT.setFullPath(os.path.join(self.getWorkDirectory(),"final.pdb"))        
        
        # Split an MTZ file into minimtz data objects
        outputFiles = ['FPHIOUT','DIFFPHIOUT']
        outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']
        infile = os.path.join(self.workDirectory,'final.mtz')
        error = self.splitHklout(outputFiles,outputColumns,infile=infile)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
                
        #Create (dummy) PROGRAMXML
        import sys
        with open(self.makeFileName("PROGRAMXML"),"wb") as programXMLFile:
            xmlStructure = ET.Element("i2Dimple")
            
            #Here transform information parsed from dimple.log into program.xml
            #Parse the dimple log (which is compatible with configparser !)
            import configparser
            configParser = configparser.SafeConfigParser()
            configParser.read(os.path.join(self.getWorkDirectory(),"dimple.log"))
            import json
            try:
                #Extract phaser outputs if any
                if "phaser" in configParser.sections():
                    phaserElement = ET.SubElement(xmlStructure,"PHASER")
                    phaserElement.text = configParser.get("phaser","status")
                refmacCycleArrays = {}
                #Extract refmac5 jelly cycles output
                if "refmac5 jelly"  in configParser.sections():
                    for property in ["iter_overall_r", "iter_free_r", "rmsBOND", "rmsANGL","rmsCHIRAL"]:
                        perCycleArray = json.loads(configParser.get("refmac5 jelly",property))
                        refmacCycleArrays[property] = perCycleArray
                #Extract refmac5 jelly cycles output
                if "refmac5 restr"  in configParser.sections():
                    for property in ["iter_overall_r", "iter_free_r", "rmsBOND", "rmsANGL","rmsCHIRAL"]:
                        perCycleArray = json.loads(configParser.get("refmac5 restr",property))
                        if property in refmacCycleArrays: refmacCycleArrays[property] += perCycleArray
                        else: refmacCycleArrays[property] = perCycleArray
                #Identify if pointless has spotted a better reindexing
                if "pointless" in configParser.sections():
                    if "alt_reindex" in [item[0] for item in configParser.items("pointless")]:
                        alt_reindexes = json.loads(configParser.get("pointless","alt_reindex"))
                        bestReindex = ("[h,k,l]",0.)
                        for alt_reindex in alt_reindexes:
                            if alt_reindex["cc"] > bestReindex[1]:
                                bestReindex = (alt_reindex["op"], alt_reindex["cc"])
                        if bestReindex[0] != "[h,k,l]":
                            ET.SubElement(xmlStructure,"REINDEX").text = "{}".format(bestReindex[0])
                            
                            reindexCommandPath = os.path.join(self.getWorkDirectory(),"reindex.txt")
                            with open(reindexCommandPath,"w") as reindexCommandFile:
                                reindexCommandFile.write("reindex {}".format(bestReindex[0]))
                            reindexExecutablePath = os.path.join(os.environ["CBIN"],"reindex")
                            import subprocess

                            initialPath = self.container.inputData.F_SIGF.fullPath.__str__()
                            reindexedPath = self.container.outputData.F_SIGF_OUT.fullPath.__str__()
                            with open(reindexCommandPath,"r") as reindexCommandFile:
                                subprocess.call([reindexExecutablePath, "HKLIN", initialPath,"HKLOUT", reindexedPath],stdin=reindexCommandFile)
                                self.container.outputData.F_SIGF_OUT.annotation = ("Observations reindexed by operator {}".format(bestReindex[0]))

                            if self.container.inputData.FREERFLAG.isSet():
                                initialPath = self.container.inputData.FREERFLAG.fullPath.__str__()
                                reindexedPath = self.container.outputData.FREERFLAG_OUT.fullPath.__str__()
                                import subprocess
                                with open(reindexCommandPath,"r") as reindexCommandFile:
                                    subprocess.call([reindexExecutablePath, "HKLIN", initialPath,"HKLOUT", reindexedPath],stdin=reindexCommandFile)
                                self.container.outputData.FREERFLAG_OUT.annotation = ("FreeR reindexed by operator {}".format(bestReindex[0]))

                if len(refmacCycleArrays) > 0:
                    refmacCyclesNode = ET.SubElement(xmlStructure,"REFMAC")
                    overall_statsNode = ET.SubElement(refmacCyclesNode,"Overall_stats")
                    stats_vs_cycleNode = ET.SubElement(overall_statsNode,"stats_vs_cycle")
                    for iCycle in range(len(refmacCycleArrays["iter_overall_r"])):
                        refmacCycleNode = ET.SubElement(stats_vs_cycleNode,"new_cycle")
                        ET.SubElement(refmacCycleNode,"cycle").text = str(iCycle+1)
                        for property, value in refmacCycleArrays.items():
                            modPropName = property
                            if modPropName == 'iter_overall_r': modPropName = 'r_factor'
                            elif modPropName == 'iter_free_r': modPropName = 'r_free'
                            elif modPropName == 'rmsANGL': modPropName = 'rmsANGLE'
                            propertyNode = ET.SubElement(refmacCycleNode, modPropName)
                            propertyNode.text = str(value[iCycle])
                #Report on blobs found
                if "find-blobs" in configParser.sections():
                    findBlobsNode = ET.SubElement(xmlStructure,"find-blobs")
                    blobPositionsArray = json.loads(configParser.get("find-blobs","blobs"))
                    blobScoresArray = json.loads(configParser.get("find-blobs","scores"))
                    for iBlob, blobScore in enumerate(blobScoresArray):
                        blobPosition = blobPositionsArray[iBlob]
                        blobNode = ET.SubElement(findBlobsNode,"Blob")
                        ET.SubElement(blobNode,"x").text = str(blobPosition[0])
                        ET.SubElement(blobNode,"y").text = str(blobPosition[1])
                        ET.SubElement(blobNode,"z").text = str(blobPosition[2])
                        ET.SubElement(blobNode,"score").text = str(blobScore)
            except Exception as err:
                self.appendErrorReport(201, err.__str__())
                return CPluginScript.FAILED
            
            logText = ET.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"r") as logFile:
                #logText.text = ET.CDATA(logFile.read())
                logText.text = base64.b64encode(logFile.read())
            
            #Extract performanceindictors from XML
            try:
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(xmlStructure.findall(".//Cycle/iter_overall_r")[-1].text)
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(xmlStructure.findall(".//Cycle/iter_free_r")[-1].text)
            except: pass
            
            programXMLFile.write(ET.tostring(xmlStructure))
        return CPluginScript.SUCCEEDED
