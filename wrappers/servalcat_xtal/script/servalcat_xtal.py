from __future__ import print_function

"""
    servalcat_xtal.py: CCP4 GUI Project
    Copyright (C) 2024 MRC-LBM, University of Southampton
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,S
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core.CCP4ErrorHandling import *
from core import CCP4Modules, CCP4XtalData, CCP4Utils
# from lxml import etree
from xml.etree import ElementTree as ET
import pathlib
import os
import json
import shutil

class servalcat_xtal(CPluginScript):
    
    TASKMODULE = 'wrappers'
    TASKTITLE = 'Refinement (servalcat)'
    TASKNAME = 'servalcat_xtal'
    TASKCOMMAND = 'servalcat'  # refine_xtal_norefmac
    # TASKCOMMAND = 'python3'  # -m servalcat refine_xtal_norefmac
    TASKVERSION= 0.0
    WHATNEXT = ['servalcat_xtal','buccaneer_build_refine_mr']
    ASYNCHRONOUS = False
    PERFORMANCECLASS = 'CRefinementPerformance'
        
    ERROR_CODES = { 201 : {'description' : 'Refmac returned with non zero status' },
                    202:  {'description': 'New library created but strictly required' },
                    203:  {'description': 'New library created', 'severity':CCP4ErrorHandling.SEVERITY_WARNING},
                    204:  {'description': 'Program completed without generating XMLOUT.' },
                    }
    
    def __init__(self,*args, **kwargs):
        super(servalcat_xtal, self).__init__(*args, **kwargs)
        self._readyReadStandardOutputHandler = self.handleReadyReadStandardOutput
        self.xmlroot = ET.Element('SERVALCAT')
        self.xmlLength = 0
        #from .refmacLogScraper import logScraper
        #self.logScraper = logScraper(xmlroot=self.xmlroot, flushXML=self.flushXML)

    @QtCore.Slot()
    def handleReadyReadStandardOutput(self):
        if not hasattr(self,'logFileHandle'): self.logFileHandle = open(self.makeFileName('LOG'),'w')
        if not hasattr(self,'logFileBuffer'): self.logFileBuffer = ''
        pid = self.getProcessId()
        qprocess = CCP4Modules.PROCESSMANAGER().getJobData(pid,attribute='qprocess')
        availableStdout = qprocess.readAllStandardOutput()
        if sys.version_info > (3,0):
            self.logFileHandle.write(availableStdout.data().decode("utf-8"))
        else:
            self.logFileHandle.write(availableStdout)
        self.logFileHandle.flush()
        # MM
        #if sys.version_info > (3,0):
        #    self.logScraper.processLogChunk(availableStdout.data().decode("utf-8"))
        #else:
        #    self.logScraper.processLogChunk(str(availableStdout))

    def json2xml(self, json_obj, tag_name=None, tag_name_subroot="subroot"):
    # https://gist.github.com/bewestphal/0d2f884b3327f31b86122b5f6a38b3e2
        def tag_name_clear(tag_name):
            tag_name_clean = tag_name.replace(" ", "_")
            tag_name_clean = tag_name_clean.replace("-", "minus")
            tag_name_clean = tag_name_clean.replace(",", "")
            tag_name_clean = tag_name_clean.replace(".", "")
            tag_name_clean = tag_name_clean.replace("(", "")
            tag_name_clean = tag_name_clean.replace(")", "")
            tag_name_clean = tag_name_clean.replace("/", "")
            tag_name_clean = tag_name_clean.replace("|", "")
            tag_name_clean = tag_name_clean.replace("*", "")
            tag_name_clean = tag_name_clean.replace(":", "_")
            return tag_name_clean
        result_list = list()
        json_obj_type = type(json_obj)
        if tag_name:
            # print(tag_name)
            tag_name_clean = tag_name_clear(tag_name)
        else:
            tag_name_clean = tag_name_subroot
        if json_obj_type is list:
            for sub_elem in json_obj:
                result_list.append("\n<%s>" % (tag_name_clean))
                result_list.append(self.json2xml(sub_elem, tag_name=tag_name))
                # tag_name = re.sub('\s\w+="\w+"', '', tag_name)
                result_list.append("</%s>" % (tag_name_clean))
            return "".join(result_list)
        if json_obj_type is dict:
            for tag_name in json_obj:
                if tag_name: tag_name_clean = tag_name_clear(tag_name)
                sub_obj = json_obj[tag_name]
                if isinstance(sub_obj, list):
                    result_list.append(self.json2xml(sub_obj, tag_name=tag_name))
                elif isinstance(sub_obj, dict):
                    result_list.append("\n<%s>" % (tag_name_clean))
                    result_list.append(self.json2xml(sub_obj, tag_name=tag_name))
                    result_list.append("\n</%s>" % (tag_name_clean))
                else:
                    result_list.append("\n<%s>" % (tag_name_clean))
                    result_list.append(self.json2xml(sub_obj, tag_name=tag_name))
                    # tag_name = re.sub('\s\w+="\w+"', '', tag_name_clean)
                    result_list.append("</%s>" % (tag_name_clean))
            return "".join(result_list)
        return "%s" % json_obj

    def xmlAddRoot(self, xmlText, xmlFilePath=None, xmlRootName=None):
        if xmlRootName:
            xmlText = f"<{xmlRootName}>\n{xmlText}\n</{xmlRootName}>"
        if xmlFilePath:
            with open(xmlFilePath, 'w') as programXmlFile:
                programXmlFile.write(xmlText)
        return xmlText

    def processInputFiles(self):
        error = None
        self.hklin = None
        dataObjects = []
        #print '\n\n\n***contentFlag',self.container.inputData.HKLIN.contentFlag
        #Append Observation with representation dependent on whether we are detwining on Is or not
        
        obsTypeRoot = 'CONTENT_FLAG_F'
        #if self.container.controlParameters.USE_TWIN and self.container.inputData.HKLIN.isSet():
        if self.container.inputData.HKLIN.isSet():
            if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
                obsTypeRoot = 'CONTENT_FLAG_I'
            elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
                obsTypeRoot = 'CONTENT_FLAG_I'
        if self.container.controlParameters.F_SIGF_OR_I_SIGI.isSet():
            if str(self.container.controlParameters.F_SIGF_OR_I_SIGI) == "F_SIGF":
                obsTypeRoot = 'CONTENT_FLAG_F'
        
        obsPairOrMean = 'MEAN'
        if self.container.controlParameters.USEANOMALOUS:
            obsPairOrMean = 'PAIR'
                
        obsType = getattr(CCP4XtalData.CObsDataFile, obsTypeRoot+obsPairOrMean)
        dataObjects += [['HKLIN', obsType]]

        #Include phase estimates if called for
        #if self.container.inputData.ABCD.isSet():
        #    dataObjects += ['ABCD']
        
        #Apply coordinate selection if set
        import os
        self.inputCoordPath = os.path.normpath(self.container.inputData.XYZIN.fullPath.__str__())
        if self.container.inputData.XYZIN.isSelectionSet():
            self.inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected.pdb'))
            self.container.inputData.XYZIN.loadFile()
            if self.container.inputData.XYZIN.isMMCIF():
                self.inputCoordPath = str(pathlib.Path(self.inputCoordPath).with_suffix('.cif'))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.inputCoordPath)

        #Create DICT by merging dictionaries in DICT_LIST
        rv = self.joinDicts(self.container.outputData.DICT, self.container.inputData.DICT_LIST)
        #CPluginScript.joinDicts(self.container.inputData.DICT.fullPath.__str__(), self.container.inputData.DICT_LIST)

        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED
    
    def processOutputFilesToDelete(self):
        try:
            self.processOutputFiles2()
        except Exception as e:
            with open("/tmp/ee.txt", "w") as ee:
                ee.write(e)
            
    def processOutputFiles(self):
        if hasattr(self,'logFileHandle'):
            self.logFileHandle.write("JOB TITLE SECTION\n")
            jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId)
            if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                self.logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
            while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
                jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobInfo["parentjobid"])
                if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                    self.logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
            self.logFileHandle.close()
        else:
            self.xmlroot.clear()
            # self.logScraper.scrapeFile( self.makeFileName('LOG') ) MM
        
        # First up check for exit status of the program
        from core.CCP4Modules import PROCESSMANAGER
        exitStatus = 0
        exitCode = 0
        try:
            exitStatus = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitStatus')
        except Exception as e:
            print(e)
            self.appendErrorReport(201, 'Exit status: Unable to recover exitStatus')
            return CPluginScript.FAILED
        if exitStatus != 0:
            self.appendErrorReport(201, 'Exit status: ' + str(exitStatus))
            return CPluginScript.FAILED
        
        # Now the exit codes...I think that non zero means Refmac identified an issue
        try:
            exitCode = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitCode')
        except:
            self.appendErrorReport(201, 'Exit code: Unable to recover exitCode')
            return CPluginScript.FAILED
        if exitCode != 0:
            import os
            try:
                logFileText = open(self.makeFileName('LOG')).read()
                # if 'Your coordinate file has a ligand which has either minimum or no description in the library' in logFileText and self.container.controlParameters.MAKE_NEW_LIGAND_EXIT.isSet() and self.container.controlParameters.MAKE_NEW_LIGAND_EXIT:
                #     self.appendErrorReport(201,'You did not supply a full ligand geometry file: either make and supply one (Make Ligand task), or set the appropriate flag in the advanced options')
                #     import re
                #     #Example line: * Plotfile: /tmp/martin/refmac5_temp1.64630_new_TM7.ps
                #     plotFiles = re.findall(r'^.*\* Plotfile:.*$',logFileText,re.MULTILINE)
                #     print(plotFiles)
                #     for plotFile in plotFiles:
                #         psfileName = plotFile.split()[-1]
                #         import shutil
                #         shutil.copyfile(psfileName, self.container.outputData.PSOUT.__str__())
                #         self.container.outputData.PSOUT.annotation.set(psfileName+' from REFMAC')
                #     return CPluginScript.UNSATISFACTORY
                # else:
                #     self.appendErrorReport(201,'Exit code: '+str(exitCode))
            except:
                self.appendErrorReport(201, 'Exit code: '+ str(exitCode))
            return CPluginScript.FAILED

        from core import CCP4XtalData
        from core import CCP4File
        import os
        
        # Need to set the expected content flag  for phases data

        self.outputCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.mmcif'))
        self.container.outputData.CIFFILE.setFullPath(self.outputCifPath)
        # TO DO what if PDB not on output
        self.outputPdbPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.pdb'))
        self.container.outputData.XYZOUT.setFullPath(self.outputPdbPath)

        self.container.outputData.CIFFILE.annotation = 'Model from refinement (mmCIF format)'
        self.container.outputData.XYZOUT.annotation = 'Model from refinement (PDB format)'
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from refinement'
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from refinement'
        # self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        # self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        self.container.outputData.TLSOUT.annotation = 'TLS parameters from refinement'
        self.container.outputData.LIBOUT.annotation = 'Generated dictionary from refinement'
        self.container.outputData.ANOMFPHIOUT.annotation = 'Weighted anomalous difference map from refinement'
        self.container.outputData.DIFANOMFPHIOUT.annotation = 'Weighted differences of anomalous difference map'

        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        
        outputFiles = ['FPHIOUT', 'DIFFPHIOUT']
        outputColumns = ['FWT,PHWT', 'DELFWT,PHDELWT']
        """if self.container.controlParameters.PHOUT:
            outputFiles+=['ABCDOUT']
            outputColumns+=['HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB']
        
        if self.container.controlParameters.USEANOMALOUS:
            outputFiles += ['ANOMFPHIOUT']
            outputColumns += ['FAN,PHAN']"""
        
        from core import CCP4XtalData
        import os
        #hkloutFile=CCP4XtalData.CMtzDataFile(os.path.join(self.getWorkDirectory(), "hklout.mtz"))
        hkloutFilePath = str(os.path.join(self.getWorkDirectory(), "refined.mtz"))
        #hkloutFile=CCP4XtalData.CMtzDataFile(os.path.join(self.getWorkDirectory(), "refined.mtz"))
        hkloutFile=CCP4XtalData.CMtzDataFile(hkloutFilePath)
        hkloutFile.loadFile()
        columnLabelsInFile = [column.columnLabel.__str__() for column in hkloutFile.fileContent.listOfColumns]
        print('columnLabelsInFile',columnLabelsInFile)
        
        """if self.container.controlParameters.USEANOMALOUS and 'DELFAN' in columnLabelsInFile and 'PHDELAN' in columnLabelsInFile:
            outputFiles += ['DIFANOMFPHIOUT']
            outputColumns += ['DELFAN,PHDELAN']"""

        error = self.splitHklout(outputFiles,outputColumns,hkloutFilePath)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

        #Use Refmacs XMLOUT as the basis for output XML.  If not existent (probably due to failure), then create a new one
        # from core import CCP4Utils
        # rxml = None
        #try:
        #    rxml = CCP4Utils.openFileToEtree(os.path.normpath(os.path.join(self.getWorkDirectory(),'XMLOUT.xml')))
        #except:
        #    rxml = ET.Element('REFMAC')
        #    self.appendErrorReport(204,self.makeFileName('PROGRAMXML'))
        #    return CPluginScript.FAILED
        # MM
        # rxml = ET.Element('REFMAC')
        # Copy across the segments of the scraped log file into this new xml root
        # for childElement in self.xmlroot: rxml.append(childElement) # MM

        # Use output JSON from servalcat as the basis for output XML
        # Get stats from JSON, convert to XML, save and load to self.xmlroot
        rxml = None
        jsonFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.json"))
        if os.path.isfile(jsonFilePath):
            with open(jsonFilePath, "r") as jsonFile:
                jsonText = jsonFile.read()
        else:
            rxml = ET.Element('SERVALCAT')
            self.appendErrorReport(204, self.makeFileName('PROGRAMXML'))
            return CPluginScript.FAILED
        try:
            jsonStats = list(json.loads(jsonText))
            # add d_max_4ssqll and d_min_4ssqll
            for i in range(len(jsonStats)):
                for j in range(len(jsonStats[i]["data"]["binned"])):
                    jsonStats[i]["data"]["binned"][j]['d_min_4ssqll'] = \
                        1 / (list(jsonStats)[i]["data"]["binned"][j]['d_min'] * list(jsonStats)[i]["data"]["binned"][j]['d_min'])
                    jsonStats[i]["data"]["binned"][j]['d_max_4ssqll'] = \
                        1 / (list(jsonStats)[i]["data"]["binned"][j]['d_max'] * list(jsonStats)[i]["data"]["binned"][j]['d_max'])
            xmlText = self.json2xml(jsonStats, tag_name_subroot="cycle")
            xmlFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.xml"))
            xmlText = self.xmlAddRoot(xmlText, xmlFilePath, xmlRootName="SERVALCAT")
            # self.xmlroot = CCP4Utils.openFileToEtree(xmlFilePath)
            self.xmlroot = ET.fromstring(xmlText)
            rxml = self.xmlroot
            et = ET.ElementTree(self.xmlroot)
            # Write program.xml_tmp and move it to program.xml
            self.flushXml()
        except:
            rxml = ET.Element('SERVALCAT')
            self.appendErrorReport(204, self.makeFileName('PROGRAMXML'))
            return CPluginScript.FAILED

        # Extract performance indicators from JSON
        statsOverall = {}
        try:
            statsOverall["-LL"] = jsonStats[-1]['data']['summary']['-LL']
            statsOverall["Rwork"] = jsonStats[-1]['data']['summary']['Rwork']
            statsOverall["Rfree"] = jsonStats[-1]['data']['summary']['Rfree']
            statsOverall["CCFworkavg"] = jsonStats[-1]['data']['summary']['CCFworkavg']
            statsOverall["CCFfreeavg"] = jsonStats[-1]['data']['summary']['CCFfreeavg']
            self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["Rwork"])
            self.container.outputData.PERFORMANCEINDICATOR.RFree.set(statsOverall["Rfree"])
        except:
            try:
                statsOverall["R2work"] = jsonStats[-1]['data']['summary']['R2work']
                statsOverall["R2free"] = jsonStats[-1]['data']['summary']['R2free']
                statsOverall["CCIworkavg"] = jsonStats[-1]['data']['summary']['CCIworkavg']
                statsOverall["CCIfreeavg"] = jsonStats[-1]['data']['summary']['CCIfreeavg']
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["R2work"]) # TO DO
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(statsOverall["R2free"])   # TO DO
            except:
                pass

        # Perform analysis of output coordinate file composition
        # TO DO what if PDB does not exist
        if False: # MM
            if os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)):
                self.container.outputData.XYZOUT.fileContent.loadFile(self.container.outputData.XYZOUT.fullPath)
                modelCompositionNode = ET.SubElement(rxml,'ModelComposition')
                for chain in self.container.outputData.XYZOUT.fileContent.composition.chains:
                    chainNode = ET.SubElement(modelCompositionNode,'Chain',id=chain)
                for monomer in self.container.outputData.XYZOUT.fileContent.composition.monomers:
                    monomerNode = ET.SubElement(modelCompositionNode,'Monomer',id=monomer)

        #Skim smartie graphs from the log file
        #smartieNode = ET.SubElement(rxml,'SmartieGraphs')
        #self.scrapeSmartieGraphs(smartieNode)
        # martin = ET.SubElement(rxml,'test') # MM
        et = ET.ElementTree(self.xmlroot)
        
        #And write out the XML
        et.write(self.makeFileName('PROGRAMXML'))

        if False: # MM
            with open(self.container.outputData.COOTSCRIPTOUT.fullPath.__str__(),"w") as cootscript:
                #Write a GUI to regions that Refmac has identified as containing duffers
                ##badStretches = self.listOfTransgressingSegments(rxml)
                badStretches = [] # MM
                if len(badStretches) > 0:
                    interestingBitsDef = 'interestingBits = ['
                    for badStretch in badStretches:
                        interestingBitsDef += ('{"chain":"%s","firstResidue":%s,"lastResidue":%s},'%(badStretch['chain'],badStretch['firstResidue'],badStretch['lastResidue']))
                    interestingBitsDef += ']\n'
                    cootscript.write(interestingBitsDef)
                    cootscript.write('ccp4i2Interface.addInterestingBitsMenu(title="Refmac-identified outliers", interestingBits=interestingBits)\n')
        return CPluginScript.SUCCEEDED

    def listOfTransgressingSegments(self, rxml):
        badResidueSet = set()
        badStretches = []
        
        outlierNodes = rxml.findall('.//OutliersByCriteria')
        if len(outlierNodes) == 0: return badStretches
        for outlierTypeNode in outlierNodes[0]:
            key = outlierTypeNode.tag
            for outlier in outlierTypeNode.findall('Outlier'):
                res1Tuple = outlier.get('chainId1'),outlier.get('resId1')
                res2Tuple = outlier.get('chainId2'),outlier.get('resId2')
                badResidueSet.add(res1Tuple)
                badResidueSet.add(res2Tuple)
        badResidueTuple = [tuple for tuple in badResidueSet]
        
        #Sort on the basis of a string formed by adding the chain and residue elements of the tuple
        orderedBadResidues = sorted(badResidueTuple, key=lambda residue: residue[0]+residue[1])
        for orderedBadResidue in orderedBadResidues:
            if len(badStretches) == 0 or badStretches[-1]['chain'] != orderedBadResidue[0] or int(orderedBadResidue[1])-int(badStretches[-1]['lastResidue']) > 1:
                badStretch = {"chain":orderedBadResidue[0], "firstResidue":orderedBadResidue[1], "lastResidue":orderedBadResidue[1]}
                badStretches.append(badStretch)
            else:
                badStretch["lastResidue"] = orderedBadResidue[1]
        return badStretches
            
            
    def scrapeSmartieGraphs(self, smartieNode):
        import sys, os
        from core import CCP4Utils
        smartiePath = os.path.join(CCP4Utils.getCCP4I2Dir(),'smartie')
        sys.path.append(smartiePath)
        import smartie
        
        logfile = smartie.parselog(self.makeFileName('LOG'))
        for smartieTable in logfile.tables():
            if smartieTable.ngraphs() > 0:
                tableelement = self.xmlForSmartieTable(smartieTable, smartieNode)
        
        return
    
    def xmlForSmartieTable(self, table, parent):
        from pimple import MGQTmatplotlib
        tableetree = MGQTmatplotlib.CCP4LogToEtree(table.rawtable())
        parent.append(tableetree)
        return tableetree

    def makeCommandAndScript(self):
        # self.appendCommandLine(['-m', 'servalcat'])
        self.appendCommandLine(['refine_xtal_norefmac'])
        self.appendCommandLine(['-o', 'refined'])
        self.hklout = os.path.join(self.workDirectory, "refined.mtz")
        self.appendCommandLine(['--model', self.inputCoordPath])
        self.appendCommandLine(['--hklin', self.hklin])

        if self.container.inputData.DICT_LIST.isSet():
            self.appendCommandLine(['--ligand'])
            for dictionary in self.container.inputData.DICT_LIST:
                if os.path.exists(str(dictionary)):
                    self.appendCommandLine([str(dictionary)])
                else:
                    print(f"WARNING: input restraint dictionary {dictionary} does not exist and so will not be used.")

        # TO DO - labin, free
        if self.container.controlParameters.F_SIGF_OR_I_SIGI.isSet():
            if str(self.container.controlParameters.F_SIGF_OR_I_SIGI) == "F_SIGF" or not self.container.controlParameters.HKLIN_IS_I_SIGI:
                self.appendCommandLine(['--labin', 'F,SIGF,FREER'])
            else:
                self.appendCommandLine(['--labin', 'I,SIGI,FREER'])
        elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR \
                or self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
            self.appendCommandLine(['--labin', 'I,SIGI,FREER'])
        else:
            self.appendCommandLine(['--labin', 'F,SIGF,FREER'])
        self.appendCommandLine(['--source', str(self.container.controlParameters.SCATTERING_FACTORS)])
        if self.container.controlParameters.NCYCLES.isSet():
            self.appendCommandLine(["--ncycle", str(self.container.controlParameters.NCYCLES)])
        # Hydrogens
        if self.container.controlParameters.HYDR_USE:
            self.appendCommandLine(['--hydrogen', str(self.container.controlParameters.HYDR_ALL)])
            #if self.container.controlParameters.HYDR_ALL:
            #    self.appendCommandLine(['--hydrogen', str(self.container.controlParameters.HYDR_ALL)])
            #else:
            #    self.appendCommandLine(['--hydrogen', 'yes'])
            if self.container.controlParameters.H_REFINE:
                self.appendCommandLine(['--refine_h'])
        else:
            self.appendCommandLine(['--hydrogen', 'no'])
        # Resolution
        if self.container.controlParameters.RES_CUSTOM:
            if self.container.controlParameters.RES_MIN.isSet():
                self.appendCommandLine(['--d_min', str(self.container.controlParameters.RES_MIN)])
            if self.container.controlParameters.RES_MAX.isSet():
                self.appendCommandLine(['--d_max', str(self.container.controlParameters.RES_MAX)])
        # Weight
        if str(self.container.controlParameters.WEIGHT_OPT) == "MANUAL" and \
                str(self.container.controlParameters.WEIGHT):
            self.appendCommandLine(['--weight', str(self.container.controlParameters.WEIGHT)])
        # ADP
        if str(self.container.controlParameters.B_REFINEMENT_MODE) == "aniso":
            self.appendCommandLine(['--adp', 'aniso'])
        elif str(self.container.controlParameters.B_REFINEMENT_MODE) == "fix":
            self.appendCommandLine(['--adp', 'fix'])
        else:
            self.appendCommandLine(['--adp', 'iso'])
        # Jelly body
        if self.container.controlParameters.USE_JELLY:
            self.appendCommandLine(['--jellybody', '--jellybody_params'])
            if self.container.controlParameters.JELLY_SIGMA.isSet():
                self.appendCommandLine([str(self.container.controlParameters.JELLY_SIGMA)])
            else:
                self.appendCommandLine(['0.01'])
            if self.container.controlParameters.JELLY_DIST.isSet():
                self.appendCommandLine([str(self.container.controlParameters.JELLY_DIST)])
            else:
                self.appendCommandLine(['4.2'])
            if self.container.controlParameters.JELLY_ONLY:
                self.appendCommandLine(['--jellyonly'])

        if self.container.controlParameters.USE_NCS:
            self.appendCommandLine(['--ncsr'])
        if self.container.controlParameters.BFACSETUSE and \
                str(self.container.controlParameters.BFACSET):
            self.appendCommandLine(['--bfactor', str(self.container.controlParameters.BFACSET)])
        if float(self.container.controlParameters.ADPR_WEIGHT) != "1":
            self.appendCommandLine(['--adpr_weight', str(self.container.controlParameters.ADPR_WEIGHT)])
        if self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT.isSet():
            self.appendCommandLine(['--max_dist_for_adp_restraint', str(self.container.controlParameters.MAX_DIST_FOR_ADP_RESTRAINT)])
        if self.container.controlParameters.ADP_RESTRAINT_POWER.isSet():
            self.appendCommandLine(['--adp_restraint_power', str(self.container.controlParameters.ADP_RESTRAINT_POWER)])
        # if self.container.controlParameters.ADP_RESTRAINT_EXP_FAC.isSet():
        #     self.appendCommandLine(['--adp_restraint_exp_fac', str(self.container.controlParameters.ADP_RESTRAINT_EXP_FAC)])
        if self.container.controlParameters.ADP_RESTRAINT_NO_LONG_RANGE:
            self.appendCommandLine(['--adp_restraint_no_long_range'])
        # if self.container.controlParametersADP_RESTRAINT_MODE.isSet():
        #     self.appendCommandLine(['--adp_restraint_mode', str(self.container.controlParameters.ADP_RESTRAINT_MODE)])
        if self.container.controlParameters.FIND_LINKS:
            self.appendCommandLine(['--find-links'])
        if self.container.controlParameters.UNRESTRAINED:
            self.appendCommandLine(['--unrestrained'])
        if self.container.controlParameters.FIX_XYZ:
            self.appendCommandLine(['--fix_xyz'])
        if self.container.controlParameters.KEEP_CHARGES:
            self.appendCommandLine(['--keep_charges'])
        if self.container.controlParameters.USE_WORK_IN_EST:
            self.appendCommandLine(['--use_work_in_est'])
        if self.container.controlParameters.RANDOMIZEUSE and \
                str(self.container.controlParameters.RANDOMIZE):
            self.appendCommandLine(['--randomize', str(self.container.controlParameters.RANDOMIZE)])
        if str(self.container.controlParameters.FREERFLAG_NUMBER):
            self.appendCommandLine(['--free', str(self.container.controlParameters.FREERFLAG_NUMBER)])
        if self.container.controlParameters.NO_SOLVENT:
            self.appendCommandLine(['--no_solvent'])

        keywordFilePath = str(os.path.join(self.getWorkDirectory(), 'keywords.txt'))
        if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                # if self.container.controlParameters.PROSMART_PROTEIN_WEIGHT.isSet():
                #     keywordFile.write("\nEXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_WEIGHT)))
                if self.container.controlParameters.PROSMART_PROTEIN_SGMN.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMN %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_SGMN)))
                if self.container.controlParameters.PROSMART_PROTEIN_SGMX.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_SGMX)))
                if self.container.controlParameters.PROSMART_PROTEIN_ALPHA.isSet():
                    keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_ALPHA)))
                if self.container.controlParameters.PROSMART_PROTEIN_DMAX.isSet():
                    keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_DMAX)))
                keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.fullPath)))
        if self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                # if self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT.isSet():
                #     keywordFile.write("\nEXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMN.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMN %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_SGMN)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_SGMX.isSet():
                    keywordFile.write("\nEXTERNAL WEIGHT SGMX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_SGMX)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA.isSet():
                    keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA)))
                if self.container.controlParameters.PROSMART_NUCLEICACID_DMAX.isSet():
                    keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_DMAX)))
                keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.fullPath)))
        if self.container.inputData.SERVALCAT_KEYWORD_FILE.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write("\n@%s"%(str(self.container.inputData.SERVALCAT_KEYWORD_FILE.fullPath)))
        # if self.container.controlParameters.EXTRAREFMACKEYWORDS.isSet():
        #     with open(keywordFilePath, "a+") as keywordFile:
        #         keywordFile.write("\n" + str(self.container.controlParameters.EXTRAREFMACKEYWORDS))

        """if self.container.inputData.SERVALCAT_KEYWORD_FILE.isSet() or self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet() or self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
            keywordFilePath = str(os.path.join(self.getWorkDirectory(), 'keywords.txt'))
            if self.container.inputData.SERVALCAT_KEYWORD_FILE.isSet():
                shutil.copy2(str(self.container.inputData.SERVALCAT_KEYWORD_FILE.fullPath), keywordFilePath)
                if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet() or self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
                    import re
                    # Remove END if was in the input keyword file
                    linesToDelete = []
                    with open(keywordFilePath, "r") as keywordFile:
                        lines = keywordFile.readlines()
                        for i, line in enumerate(lines):
                            if re.search("END", line, re.IGNORECASE):
                                linesToDelete.append(i)
                    if linesToDelete:
                        for i in range(len(linesToDelete)):
                            lines.pop(i)
                        with open(keywordFilePath, "w") as keywordFile:
                            keywordFile.writelines(lines)
            # External restraints from ProSMART
            if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet():
                with open(keywordFilePath, "a+") as keywordFile:
                    if self.container.controlParameters.PROSMART_PROTEIN_WEIGHT.isSet():
                        keywordFile.write("\nEXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_WEIGHT)))
                    if self.container.controlParameters.PROSMART_PROTEIN_ALPHA.isSet():
                        keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_ALPHA)))
                    if self.container.controlParameters.PROSMART_PROTEIN_DMAX.isSet():
                        keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_DMAX)))
                    keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.fullPath)))
            if self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
                with open(keywordFilePath, "a+") as keywordFile:
                    if self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT.isSet():
                        keywordFile.write("\nEXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT)))
                    if self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA.isSet():
                        keywordFile.write("\nEXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA)))
                    if self.container.controlParameters.PROSMART_NUCLEICACID_DMAX.isSet():
                        keywordFile.write("\nEXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_DMAX)))
                    keywordFile.write("\n@%s"%(str(self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.fullPath)))
            if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet() or self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
                with open(keywordFilePath, "a+") as keywordFile:
                    keywordFile.write("\nEND")"""
        if os.path.isfile(keywordFilePath):
            self.appendCommandLine(["--keyword_file", keywordFilePath])

        ##if self.container.controlParameters.EXTRAREFMACKEYWORDS.isSet():
        ##   self.appendCommandLine(["--keywords", self.container.controlParameters.EXTRAREFMACKEYWORDS])
            # for kwLine in str(self.container.controlParameters.EXTRAREFMACKEYWORDS).split('\n'):
            #     #print 'KwLine','['+str(kwLine)+']'
            #     self.appendCommandScript(kwLine.rstrip() + '\n')
        #self.appendCommandScript("END")

        if self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS.isSet():
            self.appendCommandLine(str(self.container.controlParameters.EXTRA_SERVALCAT_OPTIONS).split())

        jsonFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.json"))
        self.watchFile(jsonFilePath, handler=self.handleJsonChanged, unwatchWhileHandling=True)
        return CPluginScript.SUCCEEDED

    @QtCore.Slot(str)
    def handleJsonChanged(self, jsonFilePath):
        self.xmlroot.clear()
        # refmacEtree = CCP4Utils.openFileToEtree(xmlFilename)
        # # MN Here is a for example...with lxml could search for //REFMAC xpath to find REFMAC nodes. Now have a challenge since toplevel nodes are REFMAC nodes
        #if refmacEtree.tag == 'REFMAC':
        #    refmacXML = [refmacEtree]
        #else:
        #    refmacXML = refmacEtree.findall(".//REFMAC")
        #if len(refmacXML) == 1:
        #    refmacXML[0].tag="RefmacInProgress"
        #    self.xmlroot.append(refmacXML[0])
        #else:
        #    rxml = ET.Element('SERVALCAT')
        #    self.appendErrorReport(204, self.makeFileName('PROGRAMXML'))
        #    return CPluginScript.FAILED

        # Use output JSON from servalcat as the basis for output XML
        if os.path.isfile(jsonFilePath):
            with open(jsonFilePath, "r") as jsonFile:
                jsonText = jsonFile.read()
            try:
                # Get statistics from JSON, convert to XML, save and load to self.xmlroot
                jsonStats = list(json.loads(jsonText))
                # add d_max_4ssqll and d_min_4ssqll
                for i in range(len(jsonStats)):
                    for j in range(len(jsonStats[i]["data"]["binned"])):
                        jsonStats[i]["data"]["binned"][j]['d_min_4ssqll'] = \
                            1 / (list(jsonStats)[i]["data"]["binned"][j]['d_min'] * list(jsonStats)[i]["data"]["binned"][j]['d_min'])
                        jsonStats[i]["data"]["binned"][j]['d_max_4ssqll'] = \
                            1 / (list(jsonStats)[i]["data"]["binned"][j]['d_max'] * list(jsonStats)[i]["data"]["binned"][j]['d_max'])
                xmlText = self.json2xml(list(jsonStats), tag_name_subroot="cycle")
                xmlFilePath = str(os.path.join(self.getWorkDirectory(), "refined_stats.xml"))
                xmlText = self.xmlAddRoot(xmlText, xmlFilePath, xmlRootName="SERVALCAT")
                self.xmlroot = ET.fromstring(xmlText)
                # self.xmlroot = CCP4Utils.openFileToEtree(xmlFilePath)
                et = ET.ElementTree(self.xmlroot)
                # Write program.xml_tmp and move it to program.xml
                self.flushXml()
            except Exception:
                import traceback
                print(traceback.format_exc())

    def flushXml(self):  # assumes self.xmlroot and self.xmlLength are well set
        # Save program.xml if program.xml_tmp has grown
        newXml = ET.tostring(self.xmlroot)
        if len(newXml) > self.xmlLength:
            self.xmlLength = len(newXml)
            with open(self.makeFileName('PROGRAMXML')+'_tmp','w') as programXmlFile:
                CCP4Utils.writeXML(programXmlFile, newXml)
            shutil.move(self.makeFileName('PROGRAMXML')+'_tmp', self.makeFileName('PROGRAMXML'))

    def mmaakeCommandAndScript_refmac_ToDelete(self):
        import os
        from core import CCP4Utils
        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")
        # make refmac command script
        self.appendCommandLine(['XYZIN',self.inputCoordPath])
        self.appendCommandLine(['HKLIN',self.hklin])

        if self.container.outputData.DICT.isSet():
            if os.path.exists(str(self.container.outputData.DICT.fullPath)):
                self.appendCommandLine(['LIBIN',self.container.outputData.DICT.fullPath])
            else:
                print("Warning: input restraint dictionary does not exist and so will not be used.")

        if self.container.controlParameters.REFINEMENT_MODE.isSet() and str(self.container.controlParameters.REFINEMENT_MODE) == 'RESTR':
           if self.container.controlParameters.TLSMODE.isSet() and str(self.container.controlParameters.TLSMODE) == 'FILE':
               if self.container.inputData.TLSIN.isSet():
                   self.appendCommandLine(['TLSIN',self.container.inputData.TLSIN.fullPath])
        self.appendCommandLine(['XYZOUT',self.container.outputData.XYZOUT.fullPath])
        self.appendCommandLine(['HKLOUT',self.hklout])
                #self.appendCommandLine(['LIBOUT',self.libout])
        self.appendCommandLine(['LIBOUT',self.container.outputData.LIBOUT.fullPath])
        self.appendCommandLine(['XMLOUT',os.path.normpath(os.path.join(self.getWorkDirectory(),'XMLOUT.xml'))])
        if self.container.controlParameters.SCATTERING_FACTORS.isSet():
           if self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'ELECTRON':
              if self.container.controlParameters.SCATTERING_ELECTRON.isSet():
                 if self.container.controlParameters.SCATTERING_ELECTRON.__str__() == 'GAUSSIAN':
                    self.appendCommandLine(['ATOMSF',os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'atomsf_electron.lib')])
           elif self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'NEUTRON':
              self.appendCommandLine(['ATOMSF',os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'atomsf_neutron.lib')])

        if self.container.controlParameters.TITLE.isSet():
            self.appendCommandScript("TITLE %s"%(str(self.container.controlParameters.TITLE)))

        # Main refinement options, and rigid body
        
        if self.container.controlParameters.REFINEMENT_MODE.isSet():
           if str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID':
              self.appendCommandScript("MODE RIGID")
              if self.container.controlParameters.NCYCRIGID.isSet():
                 self.appendCommandScript("RIGID NCYCLE %s" % (str(self.container.controlParameters.NCYCRIGID)))

              for sel0 in self.container.controlParameters.RIGID_BODY_SELECTION:
                  sel = sel0.get()
                  if sel['groupId'] and sel['chainId'] and sel['firstRes'] and sel['lastRes']:
                      rigidText = "RIGID GROUP %s FROM %s %s TO %s %s" % (str(sel['groupId']),str(sel['firstRes']),str(sel['chainId']),str(sel['lastRes']),str(sel['chainId']))
                      self.appendCommandScript(rigidText + '\n')
           else:
              if self.container.controlParameters.NCYCLES.isSet():
                  self.appendCommandScript("NCYCLES %s"%(str(self.container.controlParameters.NCYCLES)))

              # Occupancy refinement
              if self.container.controlParameters.OCCUPANCY_GROUPS:
                 occup_groupids = []
                 for sel0 in self.container.controlParameters.OCCUPANCY_SELECTION:
                    sel = sel0.get()
                    if sel['groupId'] and sel['chainIds'] and sel['firstRes'] and sel['lastRes']:
                       occupText = "OCCUPANCY GROUP ID %s CHAIN %s RESIDUE FROM %s TO %s" % (str(sel['groupId']),str(sel['chainIds']),str(sel['firstRes']),str(sel['lastRes']))
                       if sel['atoms']:
                          occupText += " ATOM %s" % (str(sel['atoms']))
                       if sel['alt']:
                          occupText += " ALT %s" % (str(sel['alt']))
                       self.appendCommandScript(occupText + '\n')
                       occup_groupids.append(sel['groupId'])
#                    else:
#                       self.appendErrorReport(201,"Error - incorrectly specified occupancy group.")
#                       return CPluginScript.FAILED

                 if self.container.controlParameters.OCCUPANCY_REFINEMENT:
                    self.appendCommandScript("OCCUPANCY REFINE")

                 if self.container.controlParameters.OCCUPANCY_COMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE:
                       sel = sel0.get()
#                       for id in sel['groupIds'].split(' '):
#                          if not id in occup_groupids:
#                             self.appendErrorReport(201,"Error - group ID "+id+" not specified in the list of occupancy groups.")
#                             return CPluginScript.FAILED
                       occupText = "OCCUPANCY GROUP ALTS COMPLETE %s" % (str(sel['groupIds']))
                       self.appendCommandScript(occupText + '\n')

                 if self.container.controlParameters.OCCUPANCY_INCOMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE:
                       sel = sel0.get()
#                       for id in sel['groupIds'].split(' '):
#                          if not id in occup_groupids:
#                             self.appendErrorReport(201,"Error - group ID "+id+" not specified in the list of occupancy groups.")
#                             return CPluginScript.FAILED
                       occupText = "OCCUPANCY GROUP ALTS INCOMPLETE %s" % (str(sel['groupIds']))
                       self.appendCommandScript(occupText + '\n')

        if str(self.container.controlParameters.WEIGHT_OPT) == 'AUTO':
            self.appendCommandScript("WEIGHT AUTO")
        elif self.container.controlParameters.WEIGHT.isSet():
            self.appendCommandScript("WEIGHT MATRIX %s"%(str(self.container.controlParameters.WEIGHT)))
        
        if self.container.controlParameters.USE_TWIN:
            self.appendCommandScript("TWIN")
            
        if self.container.controlParameters.HYDR_USE:
            if str(self.container.controlParameters.HYDR_ALL) == 'ALL':
                self.appendCommandScript("MAKE HYDR ALL")
            else:
                self.appendCommandScript("MAKE HYDR YES")
        else:
            self.appendCommandScript("MAKE HYDR NO")
        
        # Parameters

        if self.container.controlParameters.B_REFINEMENT_MODE.isSet():
            self.appendCommandScript("REFI BREF %s"%(str(self.container.controlParameters.B_REFINEMENT_MODE)))

        if self.container.controlParameters.REFINEMENT_MODE.isSet() and str(self.container.controlParameters.REFINEMENT_MODE) == 'RESTR':
          if self.container.controlParameters.TLSMODE.isSet():
            if str(self.container.controlParameters.TLSMODE) != 'NONE' and self.container.controlParameters.NTLSCYCLES.isSet():
                self.appendCommandScript("REFI TLSC %s"%(str(self.container.controlParameters.NTLSCYCLES)))
                self.appendCommandLine(["TLSOUT",self.container.outputData.TLSOUT.fullPath])
                if self.container.controlParameters.TLSOUT_ADDU:
                    self.appendCommandScript("TLSOUT ADDU")

        if self.container.controlParameters.BFACSETUSE and self.container.controlParameters.BFACSET.isSet():
            self.appendCommandScript("BFAC SET "+str(self.container.controlParameters.BFACSET))

        # Restraints

        if self.container.controlParameters.USE_NCS:
            if str(self.container.controlParameters.NCS_TYPE) == 'L':
                self.appendCommandScript("NCSR LOCAL")
            else:
                self.appendCommandScript("NCSR GLOBAL")
        elif self.container.controlParameters.USE_LOCAL_SYMMETRY: # USE_LOCAL_SYMMETRY is not used by the refinement pipeline, but is used by the bucaneer pipeline
            self.appendCommandScript("NCSR LOCAL")


        if self.container.controlParameters.USE_JELLY:
            if self.container.controlParameters.JELLY_SIGMA.isSet():
                self.appendCommandScript("RIDG DIST SIGM %s"%(str(self.container.controlParameters.JELLY_SIGMA)))
            else:
                self.appendCommandScript("RIDG DIST SIGM 0.01")
            if self.container.controlParameters.JELLY_DIST.isSet():
                self.appendCommandScript("RIDG DIST DMAX %s"%(str(self.container.controlParameters.JELLY_DIST)))
     
        if self.container.controlParameters.MAKE_LINK:
            if self.container.controlParameters.OVERRIDE_LINK:
                self.appendCommandScript("MAKE LINK DEFINE")
            else:
                self.appendCommandScript("MAKE LINK YES")

        # Output options

        if self.container.controlParameters.MAP_SHARP:
            if self.container.controlParameters.MAP_SHARP_CUSTOM and self.container.controlParameters.BSHARP.isSet():
                self.appendCommandScript("MAPC SHAR "+str(self.container.controlParameters.BSHARP))
            else:
                self.appendCommandScript("MAPC SHAR")

        # Advanced options

        if self.container.controlParameters.SCATTERING_FACTORS.isSet():
           if self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'ELECTRON':
              if self.container.controlParameters.SCATTERING_ELECTRON.isSet():
                 if self.container.controlParameters.SCATTERING_ELECTRON.__str__() == 'GAUSSIAN':
                    self.appendCommandScript("SOURCE ELECTRON")
                 elif self.container.controlParameters.SCATTERING_ELECTRON.__str__() == 'MB':
                    self.appendCommandScript("SOURCE ELECTRON MB")
           elif self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'NEUTRON':
              self.appendCommandScript("SOURCE NEUTRON")

        if self.container.controlParameters.SCATTERING_FACTORS.isSet():
          if self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'NEUTRON':
            if self.container.controlParameters.HYDR_USE:
              if self.container.controlParameters.H_REFINE:
                 if self.container.controlParameters.H_REFINE_SELECT.__str__() == 'ALL':
                    self.appendCommandScript("HYDROGEN REFINE ALL")

                 elif self.container.controlParameters.H_REFINE_SELECT.__str__() == 'POLAR':
                    self.appendCommandScript("HYDROGEN REFINE POLAR")
                 elif self.container.controlParameters.H_REFINE_SELECT.__str__() == 'RPOLAR':
                    self.appendCommandScript("HYDROGEN REFINE RPOLAR")
              if self.container.controlParameters.H_TORSION:
                 self.appendCommandScript("RESTRAINT TORSION HYDROGEN INCLUDE ALL")
              if self.container.controlParameters.HD_FRACTION:
                 self.appendCommandScript("REFINEMENT DFRACTION")
                 if self.container.controlParameters.HD_FRACTION_TYPE.__str__() == 'ALL':
                    self.appendCommandScript("HYDROGEN DFRACTION ALL")
                 elif self.container.controlParameters.HD_FRACTION_TYPE.__str__() == 'POLAR':
                    self.appendCommandScript("HYDROGEN DFRACTION POLAR")
                 if str(self.container.controlParameters.HYDR_ALL) == 'YES':
                    if self.container.controlParameters.HD_INIT.__str__() == 'DEUTERIUM':
                       self.appendCommandScript("HYDROGEN DFRACTION INIT")
                    elif self.container.controlParameters.HD_INIT.__str__() == 'MIXTURE':
                       self.appendCommandScript("HYDROGEN DFRACTION INIT REFINEABLE 1 UNREFINEABLE 0")
                 elif str(self.container.controlParameters.HYDR_ALL) == 'ALL':
                    if self.container.controlParameters.HD_INIT_HALL.__str__() == 'DEUTERIUM':
                       self.appendCommandScript("HYDROGEN DFRACTION INIT")
                    elif self.container.controlParameters.HD_INIT_HALL.__str__() == 'MIXTURE':
                       self.appendCommandScript("HYDROGEN DFRACTION INIT REFINEABLE 1 UNREFINEABLE 0")

        if self.container.controlParameters.RES_CUSTOM:
            if self.container.controlParameters.RES_MIN.isSet() and self.container.controlParameters.RES_MAX.isSet():
                self.appendCommandScript("REFI RESO %s %s"%(str(self.container.controlParameters.RES_MIN),str(self.container.controlParameters.RES_MAX)))
            elif self.container.controlParameters.RES_MIN.isSet() and not self.container.controlParameters.RES_MAX.isSet():
                self.appendCommandScript("REFI RESO %s 999"%(str(self.container.controlParameters.RES_MIN)))
            elif not self.container.controlParameters.RES_MIN.isSet() and self.container.controlParameters.RES_MAX.isSet():
                self.appendCommandScript("REFI RESO 0 %s"%(str(self.container.controlParameters.RES_MAX)))
        elif self.container.controlParameters.RESOLUTION.isSet():   # RESOLUTION is not used by the refinement pipeline, but is used by the bucaneer pipeline
            self.appendCommandScript("REFI RESO %s"%(str(self.container.controlParameters.RESOLUTION)))

        if self.container.controlParameters.MAKE_NEW_LIGAND_EXIT:
            self.appendCommandScript("MAKE NEWLIGAND EXIT")
        else:
            self.appendCommandScript("MAKE NEWLIGAND NOEXIT")
        
        #if self.container.controlParameters.SCALETYPE.isSet():
        #    if self.container.controlParameters.SCALETYPE.__str__() == 'REFMACDEFAULT':
        #        pass
        #    if self.container.controlParameters.SCALETYPE.__str__() == 'BABINET':
        #        self.appendCommandScript("SCALE TYPE BULK")
        #        self.appendCommandScript("SOLVENT NO")
        #    if self.container.controlParameters.SCALETYPE.__str__() == 'SIMPLE':
        #        self.appendCommandScript("SCALE TYPE SIMPLE")
        #        self.appendCommandScript("SOLVENT YES")
        #        if self.container.controlParameters.SOLVENT_ADVANCED:
        #            if self.container.controlParameters.SOLVENT_VDW_RADIUS.isSet():
        #               self.appendCommandScript("SOLVENT VDWProb %s"%(str(self.container.controlParameters.SOLVENT_VDW_RADIUS)))
        #            if self.container.controlParameters.SOLVENT_IONIC_RADIUS.isSet():
        #                self.appendCommandScript("SOLVENT IONProb %s"%(str(self.container.controlParameters.SOLVENT_IONIC_RADIUS)))
        #            if self.container.controlParameters.SOLVENT_SHRINK.isSet():
        #                self.appendCommandScript("SOLVENT RSHRink %s"%(str(self.container.controlParameters.SOLVENT_SHRINK)))
        #    if self.container.controlParameters.SCALETYPE.__str__() == 'WILSON':
        #        self.appendCommandScript("SCALE TYPE SIMPLE")
        #        self.appendCommandScript("SOLVENT NO")

        if self.container.controlParameters.SCALE_TYPE.isSet():
            if self.container.controlParameters.SCALE_TYPE.__str__() == 'SIMPLE':
                self.appendCommandScript("SCALE TYPE SIMPLE")
            if self.container.controlParameters.SCALE_TYPE.__str__() == 'BABINET':
                self.appendCommandScript("SCALE TYPE BULK")
        if self.container.controlParameters.SOLVENT_MASK_TYPE.__str__() == 'EXPLICIT':
            self.appendCommandScript("SOLVENT YES")
            if self.container.controlParameters.SOLVENT_ADVANCED:
               if self.container.controlParameters.SOLVENT_VDW_RADIUS.isSet():
                  self.appendCommandScript("SOLVENT VDWProb %s"%(str(self.container.controlParameters.SOLVENT_VDW_RADIUS)))
               if self.container.controlParameters.SOLVENT_IONIC_RADIUS.isSet():
                  self.appendCommandScript("SOLVENT IONProb %s"%(str(self.container.controlParameters.SOLVENT_IONIC_RADIUS)))
               if self.container.controlParameters.SOLVENT_SHRINK.isSet():
                  self.appendCommandScript("SOLVENT RSHRink %s"%(str(self.container.controlParameters.SOLVENT_SHRINK)))
        else:
            self.appendCommandScript("SOLVENT NO")

        if self.container.controlParameters.SCATTERING_FACTORS.__str__() == 'NEUTRON':
           self.appendCommandScript("MAKE HOUT YES")
        elif self.container.controlParameters.OUTPUT_HYDROGENS.isSet():
           if self.container.controlParameters.OUTPUT_HYDROGENS.__str__() == 'YES':
              self.appendCommandScript("MAKE HOUT YES")
           elif self.container.controlParameters.OUTPUT_HYDROGENS.__str__() == 'NO':
              self.appendCommandScript("MAKE HOUT NO")

        if self.container.controlParameters.USEANOMALOUS:
            if self.container.controlParameters.USEANOMALOUSFOR.__str__() == 'OUTPUTMAPS':
                self.appendCommandScript("ANOMALOUS MAPONLY")
            elif self.container.controlParameters.USEANOMALOUSFOR.__str__() == 'SADREFINEMENT':
                self.appendCommandScript("REFINE OREFINE ANOMALOUS")
            if self.container.controlParameters.WAVELENGTH.isSet():
                self.appendCommandScript("ANOMALOUS WAVELENGTH %10.5f"%self.container.controlParameters.WAVELENGTH)
            else:
                self.appendCommandScript("ANOMALOUS WAVELENGTH %10.5f"%self.container.inputData.HKLIN.fileContent.getListOfWavelengths()[-1])
        # Additional input/output
        if self.container.controlParameters.PHOUT:
            self.appendCommandScript("PHOUT")
    
        # Filter outlier monitoring, to avoid bloating the logfile
        self.appendCommandScript("MONI DIST 1000000")
        
        # Maintain the same chain IDs as in the input model.
        #self.appendCommandScript("PDBOUT KEEP TRUE")
        # New Refmac keyword that will keep the user's chain IDs in the report page.
        self.appendCommandScript("PDBOUT KEEP USERS")
        
        # Precedence to platonyzer, so that prosmart-related weights don't interfere with platonyzer restraints.
        if self.container.inputData.PLATONYZER_RESTRAINTS.isSet():
            self.appendCommandScript("@%s"%(str(self.container.inputData.PLATONYZER_RESTRAINTS.fullPath)))

        # Next precedence to external restraints file provided
        if self.container.inputData.EXTERNAL_RESTRAINTS_FILE.isSet():
            self.appendCommandScript("@%s"%(str(self.container.inputData.EXTERNAL_RESTRAINTS_FILE.fullPath)))

        # This general external restraints option is no longer used in the Refinement pipeline, but is maintained for legacy purposes (and if e.g. needed by Lorestr)
        if self.container.inputData.EXTERNALRESTRAINTS.isSet():
            self.appendCommandScript("@%s"%(str(self.container.inputData.EXTERNALRESTRAINTS.fullPath)))
        
        # External restraints options currently available in the Refinement pipeline
        if self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.isSet():
           if self.container.controlParameters.PROSMART_PROTEIN_WEIGHT.isSet():
              self.appendCommandScript("EXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_WEIGHT)))
           if self.container.controlParameters.PROSMART_PROTEIN_ALPHA.isSet():
              self.appendCommandScript("EXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_ALPHA)))
           if self.container.controlParameters.PROSMART_PROTEIN_DMAX.isSet():
              self.appendCommandScript("EXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_PROTEIN_DMAX)))
           self.appendCommandScript("@%s"%(str(self.container.inputData.PROSMART_PROTEIN_RESTRAINTS.fullPath)))
        if self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.isSet():
           if self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT.isSet():
              self.appendCommandScript("EXTERNAL WEIGHT SCALE %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT)))
           if self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA.isSet():
              self.appendCommandScript("EXTERNAL ALPHA %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_ALPHA)))
           if self.container.controlParameters.PROSMART_NUCLEICACID_DMAX.isSet():
              self.appendCommandScript("EXTERNAL DMAX %s"%(str(self.container.controlParameters.PROSMART_NUCLEICACID_DMAX)))
           self.appendCommandScript("@%s"%(str(self.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS.fullPath)))

        labin = "LABIN FP=F SIGFP=SIGF"
        if self.container.controlParameters.USE_TWIN:
            if self.container.inputData.HKLIN.isSet():
                from core import CCP4XtalData
                if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN or self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
                    labin = "LABIN IP=I SIGIP=SIGI"
        else:
            if self.container.controlParameters.USEANOMALOUS:
                labin = "LABIN F+=Fplus SIGF+=SIGFplus F-=Fminus SIGF-=SIGFminus"


#        labin = None
#        if self.container.controlParameters.USE_TWIN and self.container.controlParameters.TWIN_TYPE.isSet() and self.container.controlParameters.TWIN_TYPE=="I":
#            if self.container.controlParameters.USEANOMALOUSFOR.isSet() and self.container.controlParameters.USEANOMALOUSFOR.__str__() != 'NOTHING':
#                labin = "LABIN I+=Iplus SIGI+=SIGIplus I-=Iminus SIGI-=SIGIminus"
#            else:
#                labin = "LABIN IP=I SIGIP=SIGI"
#        else:
#            if self.container.controlParameters.USEANOMALOUSFOR.isSet() and self.container.controlParameters.USEANOMALOUSFOR.__str__() != 'NOTHING':
#                labin = "LABIN F+=Fplus SIGF+=SIGFplus F-=Fminus SIGF-=SIGFminus"
#            else:
#                labin = "LABIN FP=F SIGFP=SIGF"

        if self.container.inputData.ABCD.isSet():
            from core import CCP4XtalData
            if  self.container.inputData.ABCD.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
                labin += " HLA=HLA HLB=HLB HLC=HLC HLD=HLD"
            else:
                labin += " PHIB=PHI FOM=FOM"

        if self.container.inputData.FREERFLAG.isSet(): 
            labin += " FREE=FREER"
            print('FreeR flag set')

        if self.container.inputData.REFMAC_KEYWORD_FILE.isSet():
            self.appendCommandScript("@%s"%(str(self.container.inputData.REFMAC_KEYWORD_FILE.fullPath)))
      
        if self.container.controlParameters.EXTRAREFMACKEYWORDS.isSet():
            for kwLine in str(self.container.controlParameters.EXTRAREFMACKEYWORDS).split('\n'):
                #print 'KwLine','['+str(kwLine)+']'
                self.appendCommandScript(kwLine.rstrip() + '\n')

        self.appendCommandScript(labin)
        self.appendCommandScript('END')
    
        return CPluginScript.SUCCEEDED

    def setProgramVersion(self):
      print('refmac.getProgramVersion')
      return CPluginScript.setProgramVersion(self,'Refmac_5')


