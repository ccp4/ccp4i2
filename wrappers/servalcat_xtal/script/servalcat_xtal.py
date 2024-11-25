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
from .json2xml import json2xml
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

        #Include phase estimates if called for
        #if self.container.inputData.ABCD.isSet():
        #    dataObjects += ['ABCD']
        
        #Apply coordinate selection if set
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

        #print '\n\n\n***contentFlag',self.container.inputData.HKLIN.contentFlag
        #Append Observation with representation dependent on whether we are detwining on Is or not

        if str(self.container.controlParameters.DATA_METHOD) == 'xtal':
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

            #Include FreeRflag if called for
            if self.container.inputData.FREERFLAG.isSet():
                dataObjects += ['FREERFLAG']
            self.hklin,error = self.makeHklin(dataObjects)
            if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED

            return CPluginScript.SUCCEEDED

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

        import os
        from core import CCP4XtalData
        from core import CCP4File
        
        # Need to set the expected content flag  for phases data

        outputCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.mmcif'))
        self.container.outputData.CIFFILE.setFullPath(outputCifPath)
        self.container.outputData.CIFFILE.annotation.set('Model from refinement (mmCIF format)')
        outputPdbPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined.pdb'))
        if os.path.isfile(outputPdbPath):
            self.container.outputData.XYZOUT.setFullPath(outputPdbPath)
            self.container.outputData.XYZOUT.annotation.set('Model from refinement (PDB format)')
        # self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        # self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        # self.container.outputData.TLSOUT.annotation = 'TLS parameters from refinement'
        # self.container.outputData.LIBOUT.annotation = 'Generated dictionary from refinement'
        self.container.outputData.FPHIOUT.annotation.set('Density map (Fourier coeff.)')
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.annotation.set('Difference density map (Fourier coeff.)')
        self.container.outputData.FPHIOUT.subType = 2
        outputFiles = ['FPHIOUT', 'DIFFPHIOUT']
        outputColumns = ['FWT,PHWT', 'DELFWT,PHDELWT']

        if str(self.container.controlParameters.DATA_METHOD) == "xtal":
            self.container.outputData.ANOMFPHIOUT.annotation.set('Weighted anomalous difference map from refinement')
            self.container.outputData.DIFANOMFPHIOUT.annotation.set('Weighted differences of anomalous difference map')
            hkloutFilePath = str(os.path.join(self.getWorkDirectory(), "refined.mtz"))
            #hkloutFile=CCP4XtalData.CMtzDataFile(os.path.join(self.getWorkDirectory(), "hklout.mtz"))
            #hkloutFile=CCP4XtalData.CMtzDataFile(os.path.join(self.getWorkDirectory(), "refined.mtz"))
            # Split out data objects that have been generated. Do this after applying the annotation, and flagging
            # above, since splitHklout needs to know the ABCDOUT contentFlag
            """if self.container.controlParameters.PHOUT:
                outputFiles+=['ABCDOUT']
                outputColumns+=['HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB']
            if self.container.controlParameters.USEANOMALOUS:
                outputFiles += ['ANOMFPHIOUT']
                outputColumns += ['FAN,PHAN']"""
        elif str(self.container.controlParameters.DATA_METHOD) == 'spa':
            hkloutFilePath = str(os.path.join(self.getWorkDirectory(), "refined_diffmap.mtz"))
            self.container.outputData.MAP_FO.annotation.set('Density map (in real space)')
            self.container.outputData.MAP_FOFC.annotation.set('Difference density map (in real space)')
            outputMapFoPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined_diffmap_normalized_fo.mrc'))
            self.container.outputData.MAP_FO.setFullPath(outputMapFoPath)
            outputMapFoFcPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined_diffmap_normalized_fofc.mrc'))
            self.container.outputData.MAP_FOFC.setFullPath(outputMapFoFcPath)
            self.container.outputData.COOTSCRIPTOUT.annotation.set('Coot script')
            outputCootScriptPath = os.path.normpath(os.path.join(self.getWorkDirectory(), 'refined_coot.py'))
            self.container.outputData.COOTSCRIPTOUT.setFullPath(outputCootScriptPath)
        hkloutFile=CCP4XtalData.CMtzDataFile(hkloutFilePath)
        hkloutFile.loadFile()
        columnLabelsInFile = [column.columnLabel.__str__() for column in hkloutFile.fileContent.listOfColumns]
        print('columnLabelsInFile', columnLabelsInFile)
        """if self.container.controlParameters.USEANOMALOUS and 'DELFAN' in columnLabelsInFile and 'PHDELAN' in columnLabelsInFile:
            outputFiles += ['DIFANOMFPHIOUT']
            outputColumns += ['DELFAN,PHDELAN']"""
        error = self.splitHklout(outputFiles, outputColumns, hkloutFilePath)
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
            xmlText = json2xml(jsonStats, tag_name_subroot="cycle")
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
        except:
            pass
        if str(self.container.controlParameters.DATA_METHOD) == 'spa':
            statsOverall["FSCaverage"] = jsonStats[-1]['data']['summary']['FSCaverage']
            self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["FSCaverage"])  # TO DO
            self.container.outputData.PERFORMANCEINDICATOR.RFree.set(0)                             # TO DO
        elif self.container.inputData.FREERFLAG.isSet():
            try:
                statsOverall["Rwork"] = jsonStats[-1]['data']['summary']['Rwork']
                statsOverall["Rfree"] = jsonStats[-1]['data']['summary']['Rfree']
                statsOverall["CCFworkavg"] = jsonStats[-1]['data']['summary']['CCFworkavg']
                statsOverall["CCFfreeavg"] = jsonStats[-1]['data']['summary']['CCFfreeavg']
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["Rwork"])
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(statsOverall["Rfree"])
            except:
                try:
                    statsOverall["R1work"] = jsonStats[-1]['data']['summary']['R1work']
                    statsOverall["R1free"] = jsonStats[-1]['data']['summary']['R1free']
                    statsOverall["CCIworkavg"] = jsonStats[-1]['data']['summary']['CCIworkavg']
                    statsOverall["CCIfreeavg"] = jsonStats[-1]['data']['summary']['CCIfreeavg']
                    self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["R1work"]) # TO DO
                    self.container.outputData.PERFORMANCEINDICATOR.RFree.set(statsOverall["R1free"])   # TO DO
                except:
                    pass
        else:
            try:
                statsOverall["R"] = jsonStats[-1]['data']['summary']['R']
                statsOverall["CCFavg"] = jsonStats[-1]['data']['summary']['CCFavg']
                self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["R"])  # TO DO
                self.container.outputData.PERFORMANCEINDICATOR.RFree.set(0)                    # TO DO
            except:
                try:
                    statsOverall["R1"] = jsonStats[-1]['data']['summary']['R1']
                    statsOverall["CCIavg"] = jsonStats[-1]['data']['summary']['CCIavg']
                    self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(statsOverall["R1"])  # TO DO
                    self.container.outputData.PERFORMANCEINDICATOR.RFree.set(0)                     # TO DO
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
        if str(self.container.controlParameters.DATA_METHOD) == "spa":
            # options only for servalcat refine_spa_norefmac
            self.appendCommandLine(['refine_spa_norefmac'])
            self.appendCommandLine(['--halfmaps',
                                    str(self.container.inputData.MAPIN1.fullPath),
                                    str(self.container.inputData.MAPIN2.fullPath)])
            if self.container.inputData.MAPMASK.isSet():
                self.appendCommandLine(['--mask_for_fofc',
                                        str(self.container.inputData.MAPMASK.fullPath)])
            if self.container.controlParameters.MASK_RADIUS.isSet():
                self.appendCommandLine(['-r', str(float(str(self.container.controlParameters.MASK_RADIUS)))])
            self.appendCommandLine(['-d', str(self.container.controlParameters.RES_MIN)])
            self.appendCommandLine(['--source', "electron"])
            if self.container.controlParameters.PIXEL_SIZE.isSet():
                self.appendCommandLine(['--pixel_size', str(self.container.controlParameters.PIXEL_SIZE)])
            if self.container.controlParameters.POINTGROUP.isSet():
                self.appendCommandLine(['--pg', str(self.container.controlParameters.POINTGROUP)])
            if self.container.controlParameters.TWIST.isSet():
                self.appendCommandLine(['--twist', str(self.container.controlParameters.TWIST)])
            if self.container.controlParameters.RISE.isSet():
                self.appendCommandLine(['--rise', str(self.container.controlParameters.RISE)])
            if self.container.controlParameters.CENTER_X.isSet() and \
                    self.container.controlParameters.CENTER_Y.isSet() and \
                    self.container.controlParameters.CENTER_Z.isSet():
                self.appendCommandLine(['--center',
                                        str(self.container.controlParameters.CENTER_X),
                                        str(self.container.controlParameters.CENTER_Y),
                                        str(self.container.controlParameters.CENTER_Z)])
            if self.container.controlParameters.AXIS1_X.isSet() and \
                    self.container.controlParameters.AXIS1_Y.isSet() and \
                    self.container.controlParameters.AXIS1_Z.isSet():
                self.appendCommandLine(['--axis1',
                                        str(self.container.controlParameters.AXIS1_X),
                                        str(self.container.controlParameters.AXIS1_Y),
                                        str(self.container.controlParameters.AXIS1_Z)])
            if self.container.controlParameters.AXIS2_X.isSet() and \
                    self.container.controlParameters.AXIS2_Y.isSet() and \
                    self.container.controlParameters.AXIS2_Z.isSet():
                self.appendCommandLine(['--axis2',
                                        str(self.container.controlParameters.AXIS2_X),
                                        str(self.container.controlParameters.AXIS2_Y),
                                        str(self.container.controlParameters.AXIS2_Z)])
            if self.container.controlParameters.IGNORE_SYMMETRY:
                self.appendCommandLine(['--ignore_symmetry'])
            if self.container.controlParameters.PIXEL_SIZE.isSet():
                self.appendCommandLine(['--pixel_size', str(float(str(self.container.controlParameters.PIXEL_SIZE)))])
            if self.container.controlParameters.CROSS_VALIDATION:
                self.appendCommandLine(['--cross_validation'])
            if self.container.controlParameters.BLURUSE and self.container.controlParameters.BLUR.isSet():
                self.appendCommandLine(['--blur', str(self.container.controlParameters.BLUR)])
            self.hklout = ""

        elif str(self.container.controlParameters.DATA_METHOD) == "xtal":
            # options only for servalcat refine_xtal_norefmac
            self.appendCommandLine(['refine_xtal_norefmac'])
            self.appendCommandLine(['--hklin', self.hklin])
            if self.container.controlParameters.F_SIGF_OR_I_SIGI.isSet():
                if str(self.container.controlParameters.F_SIGF_OR_I_SIGI) == "F_SIGF" or not self.container.controlParameters.HKLIN_IS_I_SIGI:
                    labin = 'F,SIGF'
                else:
                    labin = 'I,SIGI'
            elif self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR \
                    or self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
                labin = 'I,SIGI'
            else:
                labin = 'F,SIGF'
            if self.container.inputData.FREERFLAG.isSet():
                if self.container.controlParameters.FREERFLAG_NUMBER.isSet():
                    self.appendCommandLine(['--free', str(self.container.controlParameters.FREERFLAG_NUMBER)])
                labin += ",FREER"
            self.appendCommandLine(['--labin', labin])
            if self.container.controlParameters.USE_TWIN:  # I,SIGI for optimal results
                self.appendCommandLine(['--twin'])
            self.hklout = os.path.join(self.workDirectory, "refined.mtz")
            # Resolution
            if self.container.controlParameters.RES_CUSTOM:
                if self.container.controlParameters.RES_MIN.isSet():
                    self.appendCommandLine(['--d_min', str(self.container.controlParameters.RES_MIN)])
                if self.container.controlParameters.RES_MAX.isSet():
                    self.appendCommandLine(['--d_max', str(self.container.controlParameters.RES_MAX)])
            self.appendCommandLine(['--source', str(self.container.controlParameters.SCATTERING_FACTORS)])
            if self.container.controlParameters.USE_WORK_IN_EST:
                self.appendCommandLine(['--use_work_in_est'])
            if self.container.controlParameters.NO_SOLVENT:
                self.appendCommandLine(['--no_solvent'])          

        # options for both servalcat refine_xtal_norefmac and servalcat refine_spa_norefmac
        self.appendCommandLine(['--model', self.inputCoordPath])
        self.appendCommandLine(['-o', 'refined'])
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
        if self.container.controlParameters.H_OUT:
            self.appendCommandLine(['--hout'])
        if self.container.inputData.DICT_LIST.isSet():
            if len(self.container.inputData.DICT_LIST) > 0:
                self.appendCommandLine(['--ligand'])
                for dictionary in self.container.inputData.DICT_LIST:
                    if os.path.exists(str(dictionary)):
                        self.appendCommandLine([str(dictionary)])
                    else:
                        print(f"WARNING: input restraint dictionary {dictionary} does not exist and so will not be used.")
        # Weight
        if str(self.container.controlParameters.WEIGHT_OPT) == "MANUAL":
            if self.container.controlParameters.WEIGHT.isSet():
                self.appendCommandLine(['--weight', str(self.container.controlParameters.WEIGHT)])
        else:  # WEIGHT_OPT == "AUTO":
            if self.container.controlParameters.WEIGHT_NO_ADJUST:
                self.appendCommandLine(['--no_weight_adjust'])
            elif self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN.isSet() and \
                    self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX.isSet():
                self.appendCommandLine(['--target_bond_rmsz_range',
                                        str(self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MIN),
                                        str(self.container.controlParameters.WEIGHT_TARGET_BOND_RMSZ_RANGE_MAX)])
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
                self.appendCommandLine(['0.01'])  # default
            if self.container.controlParameters.JELLY_DIST.isSet():
                self.appendCommandLine([str(self.container.controlParameters.JELLY_DIST)])
            else:
                self.appendCommandLine(['4.2'])  # default
            if self.container.controlParameters.JELLY_ONLY:
                self.appendCommandLine(['--jellyonly'])
        if self.container.controlParameters.UNRESTRAINED:
            self.appendCommandLine(['--unrestrained'])
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
            self.appendCommandLine(['--find_links'])
        if self.container.controlParameters.FIX_XYZ:
            self.appendCommandLine(['--fix_xyz'])
        if self.container.controlParameters.KEEP_CHARGES:
            self.appendCommandLine(['--keep_charges'])
        if self.container.controlParameters.RANDOMIZEUSE and \
                str(self.container.controlParameters.RANDOMIZE):
            self.appendCommandLine(['--randomize', str(self.container.controlParameters.RANDOMIZE)])

        keywordFilePath = str(os.path.join(self.getWorkDirectory(), 'keywords.txt'))

        # Occupancy refinement
        if self.container.controlParameters.OCCUPANCY_GROUPS:
            with open(keywordFilePath, "a+") as keywordFile:
                if self.container.controlParameters.OCCUPANCY_REFINEMENT:
                    keywordFile.write("OCCUPANCY REFINE \n")  # self.appendCommandScript("OCCUPANCY REFINE")
                # if self.container.controlParameters.OCCUPANCY_NCYCLE.isSet():
                #     keywordFile.write("OCCUPANCY REFINE NCYCLE " + str(self.container.controlParameters.OCCUPANCY_NCYCLE) + "\n")
                occup_groupids = []
                for sel0 in self.container.controlParameters.OCCUPANCY_SELECTION:
                    sel = sel0.get()
                    if sel['groupId'] and sel['chainIds'] and sel['firstRes'] and sel['lastRes']:
                        occupText = "OCCUPANCY GROUP ID %s CHAIN %s RESIDUE FROM %s TO %s" % (str(sel['groupId']),str(sel['chainIds']),str(sel['firstRes']),str(sel['lastRes']))
                        if sel['atoms']:
                            occupText += " ATOM %s" % (str(sel['atoms']))
                        if sel['alt']:
                            occupText += " ALT %s" % (str(sel['alt']))
                        keywordFile.write(occupText + '\n')  # self.appendCommandScript(occupText + '\n')
                        occup_groupids.append(sel['groupId'])
#                   else:
#                       self.appendErrorReport(201,"Error - incorrectly specified occupancy group.")
#                       return CPluginScript.FAILED
                if self.container.controlParameters.OCCUPANCY_COMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_COMPLETE_TABLE:
                        sel = sel0.get()
#                       for id in sel['groupIds'].split(' '):
#                          if not id in occup_groupids:
#                             self.appendErrorReport(201,"Error - group ID "+id+" not specified in the list of occupancy groups.")
#                             return CPluginScript.FAILED
                        occupText = "OCCUPANCY GROUP ALTS COMPLETE %s" % (str(sel['groupIds']))
                        keywordFile.write(occupText + '\n')  # self.appendCommandScript(occupText + '\n')
                if self.container.controlParameters.OCCUPANCY_INCOMPLETE:
                    for sel0 in self.container.controlParameters.OCCUPANCY_INCOMPLETE_TABLE:
                        sel = sel0.get()
#                       for id in sel['groupIds'].split(' '):
#                          if not id in occup_groupids:
#                             self.appendErrorReport(201,"Error - group ID "+id+" not specified in the list of occupancy groups.")
#                             return CPluginScript.FAILED
                        occupText = "OCCUPANCY GROUP ALTS INCOMPLETE %s" % (str(sel['groupIds']))
                        keywordFile.write(occupText + '\n')  # self.appendCommandScript(occupText + '\n')

        if self.container.inputData.METALCOORD_RESTRAINTS.isSet():
            with open(keywordFilePath, "a+") as keywordFile:
                keywordFile.write("\n@%s"%(str(self.container.inputData.METALCOORD_RESTRAINTS.fullPath)))
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
                xmlText = json2xml(list(jsonStats), tag_name_subroot="cycle")
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

    def setProgramVersion(self):
      print('refmac.getProgramVersion')
      return CPluginScript.setProgramVersion(self,'Refmac_5')
