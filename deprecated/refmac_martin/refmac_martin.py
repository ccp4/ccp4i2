from __future__ import print_function

"""
    refmac.py: CCP4 GUI Project
    Copyright (C) 2010 University of York
    
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

from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from xml.etree import ElementTree as ET

class refmac_martin(CPluginScript):
    
    TASKMODULE = 'deprecated'
    TASKTITLE = 'Basic refinement with Refmac5'
    TASKNAME = 'refmac_martin'
    TASKCOMMAND = 'refmac5'
    TASKVERSION= 0.0
    WHATNEXT = ['refmac_martin','buccaneer_build_refine_mr']
    ASYNCHRONOUS = False
        
    ERROR_CODES = { 201 : {'description' : 'Refmac returned with non zero status' },
                    202:  {'description': 'New library created but strictly required' },
                    203:  {'description': 'New library created', 'severity':CCP4ErrorHandling.SEVERITY_WARNING},
                    204:  {'description': 'Program completed without generating XMLOUT.' },
                    }

    def processInputFiles(self):
        from core import CCP4XtalData
        error = None
        self.hklin = None
        dataObjects = []
        print('\n\n\n***contentFlag',self.container.inputData.F_SIGF.contentFlag)
        #Append Observation with representation dependent on whether we are detwining on Is or not
        if self.container.controlParameters.TWIN_TYPE.isSet() and self.container.controlParameters.TWIN_TYPE=="ON_I":
            dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN]]
        else:
            dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        #Include phase estimates if called for
        if self.container.inputData.ABCD.isSet():
            dataObjects += ['ABCD']
        
        #Apply coordinate selection if set
        import os
        self.inputCoordPath = os.path.normpath(self.container.inputData.XYZIN.fullPath.__str__())
        if self.container.inputData.XYZIN.isSelectionSet():
            self.inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected.pdb'))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.inputCoordPath)

        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        #First up check for exit status of the program
        from core.CCP4Modules import PROCESSMANAGER
        exitStatus = 0
        exitCode=0
        try:
            exitStatus = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitStatus')
        except:
            self.appendErrorReport(201,'Exit status: Unable to recover exitStatus')
            return CPluginScript.FAILED
        if exitStatus != 0:
            self.appendErrorReport(201,'Exit status: '+str(exitStatus))
            return CPluginScript.FAILED
        
        #Now the exit codes...I think that non zero means Refmac identified an issue
        try:
            exitCode = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitCode')
        except:
            self.appendErrorReport(201,'Exit code: Unable to recover exitCode')
            return CPluginScript.FAILED
        if exitCode != 0:
            self.appendErrorReport(201,'Exit code: '+str(exitCode))
            return CPluginScript.FAILED

        #Now check if LIBOUT has been created
        import os
        if os.path.isfile(str(self.container.outputData.LIBOUT.fullPath)):
            if self.container.controlParameters.MAKE_NEW_LIGAND_EXIT.isSet():
                if self.container.controlParameters.MAKE_NEW_LIGAND_EXIT:
                    self.appendErrorReport(202,'Exiting')
                    return CPluginScript.FAILED
            else:
                self.appendErrorReport(203,'Continuing')

        from core import CCP4XtalData
        from core import CCP4File
        import os
        
        # Need to set the expected content flag  for phases data
        self.container.outputData.XYZOUT.annotation = 'Model from refinement'
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from refinement'
        #self.container.outputData.FPHIOUT.subType = 1
        #self.container.outputData.FPHIOUT.contentFlag = 1
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from refinement'
        #self.container.outputData.DIFFPHIOUT.subType = 2
        #self.container.outputData.DIFFPHIOUT.contentFlag = 1
        self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        self.container.outputData.TLSOUT.annotation = 'TLS parameters from refinement'
        self.container.outputData.LIBOUT.annotation = 'Generated dictionary from refinement'

        # Split out data objects that have been generated. Do this (bizzarely) after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        
        error = self.splitHklout(['FPHIOUT','DIFFPHIOUT','ABCDOUT'],['FWT,PHWT','DELFWT,PHDELWT','HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB'])
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        
        #Use Refmacs XMLOUT as the basis for output XML.  If not existent (probably due to failure), then create a new one
        from core import CCP4Utils
        rxml = None
        try:
            rxml = CCP4Utils.openFileToEtree(str(self.container.outputData.XMLOUT.fullPath))
        except:
            rxml = ET.Element('REFMAC')
            self.appendErrorReport(204,'XML Not generated')
            return CPluginScript.FAILED

        outliersByCriteria={}
        #Do the comprehensive log-greppery to add value to the output XML
        self.skimLog(root=rxml, outliersByCriteria=outliersByCriteria)
        rstats = rxml.findall(".//REFMAC/Overall_stats/stats_vs_cycle")
        if len(rstats)>0:
            for node in rstats[0].findall("new_cycle[last()]/r_factor | new_cycle[last()]/r_free | new_cycle[last()]/rmsBOND |  new_cycle[last()]/rmsANGLE"):
                if node.tag == 'r_factor':
                    self.container.outputData.FINALRFACTOR = float(node.text)
                    if self.container.outputData.FINALRFACTOR>0.0:
                      self.container.outputData.PERFORMANCEINDICATOR.RFactor = float(node.text)
                elif node.tag == 'r_free':
                    self.container.outputData.FINALRFREE = float(node.text)
                    if self.container.outputData.FINALRFREE>0.0:
                      self.container.outputData.PERFORMANCEINDICATOR.RFree = float(node.text)
                elif node.tag == 'rmsBOND':
                    self.container.outputData.FINALRMSBOND = float(node.text)
                elif node.tag == 'rmsANGLE':
                    self.container.outputData.FINALRMSANGLE = float(node.text)
    
        #Perform analysis of output coordinate file composition
        if os.path.isfile(str(self.container.outputData.XYZOUT.fullPath)):
            from core.CCP4ModelData import CPdbData
            aCPdbData = CPdbData()
            aCPdbData.loadFile(self.container.outputData.XYZOUT.fullPath)
            modelCompositionNode = ET.SubElement(rxml,'ModelComposition')
            for chain in aCPdbData.composition.chains:
                chainNode = ET.SubElement(modelCompositionNode,'Chain',id=chain)
            for monomer in aCPdbData.composition.monomers:
                monomerNode = ET.SubElement(modelCompositionNode,'Monomer',id=monomer)

        #Skim smartie graphs from the log file
        smartieNode = ET.SubElement(rxml,'SmartieGraphs')
        self.scrapeSmartieGraphs(smartieNode)
        et = ET.ElementTree(rxml)
        
        #And write out the XML
        et.write(self.container.outputData.XMLOUT.fullPath.__str__())
       
        with open(self.container.outputData.COOTSCRIPTOUT.fullPath.__str__(),"w") as cootscript:
            #Write a GUI to regions that Refmac has identified as containing duffers
            badStretches = self.listOfTransgressingSegments(outliersByCriteria)
            if len(badStretches) > 0:
                interestingBitsDef = 'interestingBits = ['
                for badStretch in badStretches:
                    interestingBitsDef += ('{"chain":"%s","firstResidue":%s,"lastResidue":%s},'%(badStretch['chain'],badStretch['firstResidue'],badStretch['lastResidue']))
                interestingBitsDef += ']\n'
                cootscript.write(interestingBitsDef)
                cootscript.write('ccp4i2Interface.addInterestingBitsMenu(title="Refmac-identified outliers", interestingBits=interestingBits)\n')
        return CPluginScript.SUCCEEDED

    def skimLog(self,root=None, outliersByCriteria={}):
        print('\n\n in skimLog')
        #Pass through logfile and recover some things not currently captured in XMLOUT

        linesRead = open(self.makeFileName('LOG')).readlines()

        lastCycleComplete = False
        parsingAlignment = False
        parsingOutliers = False
        parsingTwin1 = False
        parsingTwin2 = False
        outliersLine = 0
        outliersMode = 'None'
        errorNode = None

        for line in linesRead:
            #recover value used for weighting and store as output data
            if line.startswith(' Weight matrix'):
                terms = line.split()
                self.container.outputData.WEIGHTUSED = terms[2]
            elif 'Things for loggraph' in line:
                lastCycleComplete = True
            elif 'Alignment results' in line:
                parsingAlignment = True
                ncsNode = ET.SubElement(root,'NCS')
            elif 'Bond distance outliers' in line or 'VDW outliers' in line:
                if lastCycleComplete:
                    parsingOutliers = True
                    outliersLine = 0
                    if 'Bond distance outliers' in line:
                        outliersMode = 'Bond'
                    elif 'VDW outliers' in line:
                        outliersMode = 'VDW'
                    outliersByCriteria[outliersMode] = {}
                    outliersByCriteria[outliersMode]['Outliers'] = []
            elif 'Twin operators with Rmerge' in line:
                #print '\nEntering parsingTwin1'
                parsingTwin1 = True
                self.twinningNode = ET.SubElement(root,'Twinning')
                tokens = line.split()
                if len(tokens) > 5:
                    twinningRmergeFilterNode = ET.SubElement(self.twinningNode,'RmergeFilter')
                    twinningRmergeFilterNode.text = tokens[5]
            elif 'Twin domains with fraction' in line:
                parsingTwin2 = True
                tokens = line.split()
                if len(tokens) > 5:
                    twinningFractionFilterNode = ET.SubElement(self.twinningNode,'FractionFilter')
                    twinningFractionFilterNode.text = tokens[5]
            
            if parsingTwin1:
                if '---' in line:
                    #print '\nExiting parsingTwin1'
                    parsingTwin1 = False
                else:
                    tokens = line.split(':')
                    #print '\n\n Splitting a partingTwin1 Line'
                    #print tokens
                    if len(tokens) > 1:
                        twinningBySymmetryOperatorNode = ET.SubElement(self.twinningNode,'TwinningSymmetryOperator')
                        operatorTokens = tokens[0].split()
                        if len(operatorTokens) > 3:
                            operatorText = operatorTokens[-3]+operatorTokens[-2]+operatorTokens[-1]
                            twinningSymmetryOperatorNode = ET.SubElement(twinningBySymmetryOperatorNode,'SymmetryOperator')
                            twinningSymmetryOperatorNode.text = operatorText
                        otTokens = tokens[1].split('=')
                        if len(otTokens) > 1:
                            operatorText = otTokens[1]
                            if len(operatorTokens) > 1:
                                twinningRmergeNode = ET.SubElement(twinningBySymmetryOperatorNode,'Rmerge')
                                twinningRmergeNode.text = operatorText
            
            if parsingTwin2:
                if '---' in line:
                    parsingTwin2 = False
                else:
                    tokens = line.split(':')
                    if len(tokens) > 3:
                        twinningByTwinOperatorNode = ET.SubElement(self.twinningNode,'TwinOperator')
                        twinningTwinOperatorNode = ET.SubElement(twinningByTwinOperatorNode,'TwinOperator')
                        twinningTwinOperatorNode.text = tokens[1]
                        try:
                            fractionText = tokens[2].split(';')[0].split('=')[1]
                            if len(operatorTokens) > 1:
                                twinningFractionNode = ET.SubElement(twinningByTwinOperatorNode,'Fraction')
                                twinningFractionNode.text = fractionText
                        except:
                            pass
            
            if parsingAlignment:
                if '---' not in line and ':' not in line and '*' not in line:
                    parsingAlignment = False
            
            if lastCycleComplete:
                pass
            
            if parsingAlignment:
                if 'No of aligned' not in line and '---' not in line and 'Alignment results' not in line:
                    tokens = line.split(':')
                    equivalenceNode = ET.SubElement(ncsNode,'equivalence',selection1 = tokens[2],selection2=tokens[3],nEquivalent=tokens[4],score=tokens[5],rms=tokens[6],ave_rmsLoc=tokens[7])
            
            if parsingOutliers:
                outliersLine += 1
                if outliersLine == 1:
                    pass
                elif outliersLine == 2:
                    pass
                elif outliersLine == 3:
                    outliersByCriteria[outliersMode]['Criteria'] = line
                elif outliersLine == 4:
                    pass
                else:
                    tokens = line.split()
                    if len(tokens) == 0:
                        parsingOutliers = False
                    else:
                        outlier = {}
                        outlier['chainId1'] = line[0:1]
                        outlier['resId1'] = line[4:8]
                        outlier['at1'] = line[12:16]
                        outlier['ins1'] = line[17:18]
                        outlier['chainId2'] = line[21:22]
                        outlier['resId2'] = line[25:29]
                        outlier['at2'] = line[33:37]
                        outlier['ins2'] = line[38:39]
                        outlier['mod'] = line[45:51]
                        outlier['ideal'] = line[56:62]
                        if outliersMode is 'Bond':
                            outlier['sigma'] = line[80:86]
                        else:
                            outlier['sigma'] = line[79:84]
                            outlier['operator'] = line[90:102]
                            outlier['type'] = line[111:112]
                        outliersByCriteria[outliersMode]['Outliers'].append(outlier)
            
            if 'ERROR' in line or 'error' in line or 'Error' in line:
                if errorNode is None: errorNode = ET.SubElement(root,'ErrorLines')
                errorLine = ET.SubElement(errorNode,'ErrorLine')
                errorLine.text = line

        if len(list(outliersByCriteria.keys())) is not 0:
            outliersNode = ET.SubElement(root,'OutliersByCriteria')
            for key in list(outliersByCriteria.keys()):
                outlierDict = outliersByCriteria[key]
                outlierCriterionNode = ET.SubElement(outliersNode,key)
                outlierCriterionCriteriaNode = ET.SubElement(outlierCriterionNode,'Criteria')
                outlierCriterionCriteriaNode.text = outlierDict['Criteria']
                outliers = outlierDict['Outliers']
                for outlier in outliers:
                    outlierNode = ET.SubElement(outlierCriterionNode,'Outlier')
                    for key in list(outlier.keys()):
                        outlierNode.set(key,outlier[key])

    def listOfTransgressingSegments(self, outliersByCriteria):
        from sets import Set
        badResidueSet = Set()
        for key in list(outliersByCriteria.keys()):
            outlierDict = outliersByCriteria[key]
            outliers = outlierDict['Outliers']
            for outlier in outliers:
                res1Tuple = outlier.get('chainId1'),outlier.get('resId1')
                res2Tuple = outlier.get('chainId2'),outlier.get('resId2')
                badResidueSet.add(res1Tuple)
                badResidueSet.add(res2Tuple)
        badResidueTuple = [tuple for tuple in badResidueSet]
        #Sort on the basis of a string formed by adding the chain and residue elements of the tuple
        orderedBadResidues = sorted(badResidueTuple, key=lambda residue: residue[0]+residue[1])
        badStretches = []
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
        import os
        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")
        # make refmac command script
        self.appendCommandLine(['XYZIN',self.inputCoordPath])
        self.appendCommandLine(['HKLIN',self.hklin])
        if self.container.inputData.DICT.isSet():
            self.appendCommandLine(['LIBIN',self.container.inputData.DICT.fullPath])
        if self.container.inputData.TLSIN.isSet():
            self.appendCommandLine(['TLSIN',self.container.inputData.TLSIN.fullPath])
        self.appendCommandLine(['XYZOUT',self.container.outputData.XYZOUT.fullPath])
        self.appendCommandLine(['HKLOUT',self.hklout])
                #self.appendCommandLine(['LIBOUT',self.libout])
        self.appendCommandLine(['LIBOUT',self.container.outputData.LIBOUT.fullPath])
        self.appendCommandLine(['XMLOUT',self.container.outputData.XMLOUT.fullPath])
        if self.container.controlParameters.TITLE.isSet():
            self.appendCommandScript("TITLE %s"%(str(self.container.controlParameters.TITLE)))
        if self.container.controlParameters.NCYCLES.isSet():
            self.appendCommandScript("NCYCLES %s"%(str(self.container.controlParameters.NCYCLES)))
        if self.container.controlParameters.WEIGHT.isSet():
            self.appendCommandScript("WEIGHT MATRIX %s"%(str(self.container.controlParameters.WEIGHT)))
        else:
            self.appendCommandScript("WEIGHT AUTO")
        if str(self.container.controlParameters.MAP_SHARP) != 'None':
            if self.container.controlParameters.MAP_SHARP.isSet():
                if self.container.controlParameters.MAP_SHARP:
                    self.appendCommandScript("MAPCALCULATE SHARP")
        if self.container.controlParameters.BSHARP.isSet():
            self.appendCommandScript("MAPCALCULATE SHARP "+str(self.container.controlParameters.BSHARP))
        if self.container.controlParameters.BFACSET.isSet():
            self.appendCommandScript("BFAC SET "+str(self.container.controlParameters.BFACSET))
        if self.container.controlParameters.TWIN_TYPE != "NONE":
            self.appendCommandScript("TWIN")
        if self.container.controlParameters.RESOLUTION:
            self.appendCommandScript("RESOLUTION %s"%(str(self.container.controlParameters.RESOLUTION)))
        if self.container.controlParameters.USE_LOCAL_SYMMETRY:
            self.appendCommandScript("NCSR LOCAL")
        if self.container.controlParameters.MAKE_NEW_LIGAND_EXIT:
            self.appendCommandScript("MAKE NEWLIGAND EXIT")
        else:
            self.appendCommandScript("MAKE NEWLIGAND NOEXIT")
        if self.container.controlParameters.B_REFINEMENT_MODE.isSet():
            self.appendCommandScript("REFI BREF %s"%(str(self.container.controlParameters.B_REFINEMENT_MODE)))
        if self.container.controlParameters.HYDROGENS.isSet():
            self.appendCommandScript("MAKE HYDROGENS %s"%(str(self.container.controlParameters.HYDROGENS)))
        if self.container.controlParameters.NTLSCYCLES.isSet():
            if self.container.controlParameters.NTLSCYCLES > 0:
                self.appendCommandScript("REFI TLSC %s"%(str(self.container.controlParameters.NTLSCYCLES)))
                self.appendCommandLine(['TLSOUT',self.container.outputData.TLSOUT.fullPath])
        if self.container.controlParameters.TLSOUT_ADDU.isSet():
            if self.container.controlParameters.TLSOUT_ADDU:
                self.appendCommandScript("TLSOUT ADDU")
        if self.container.controlParameters.EXTRAREFMACKEYWORDS.isSet():
            for kwLine in str(self.container.controlParameters.EXTRAREFMACKEYWORDS).split('\n'):
                #print 'KwLine','['+str(kwLine)+']'
                self.appendCommandScript(kwLine.rstrip() + '\n')
    
        #if self.container.controlParameters.PHOUT:
        self.appendCommandScript("PHOUT")
        if self.container.inputData.EXTERNALRESTRAINTS.isSet():
            self.appendCommandScript("@%s"%(str(self.container.inputData.EXTERNALRESTRAINTS.fullPath)))
    
        labin = None
        if self.container.controlParameters.TWIN_TYPE.isSet() and self.container.controlParameters.TWIN_TYPE=="ON_I":
            labin = "labin IP=I SIGIP=SIGI"
        else:
            labin = "LABIN FP=F SIGFP=SIGF"

        if self.container.inputData.ABCD.isSet():
            from core import CCP4XtalData
            if  self.container.inputData.ABCD.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
                labin += " HLA=HLA HLB=HLB HLC=HLC HLD=HLD"
            else:
                labin += " PHIB=PHI FOM=FOM"
        if self.container.inputData.FREERFLAG.isSet(): 
            labin += " FREE=FREER"
            print('FreeR flag set')

        self.appendCommandScript(labin)
        self.appendCommandScript('END')
    
        return CPluginScript.SUCCEEDED

    def logToXML(self, refmacInProgressElement):
        nTLSCycles = 0
        nCGMATCycles = 0
        nCycles = 0
        cyclesDone = []
        linesRead = open(self.makeFileName('LOG')).readlines()
        for line in linesRead:
            if line.strip().startswith('CGMAT cycle number ='):
                terms = line.split()
                nCGMATCycles = int(terms[4])
                nCycles = nTLSCycles + nCGMATCycles
                newCycle = {"number":str(nCycles)}
                cyclesDone.append(newCycle)
            elif line.strip().startswith('***TLS refinement cycle***'):
                terms = line.split()
                nTLSCycles = int(terms[3])
                nCycles = nTLSCycles + nCGMATCycles
                newCycle = {"number":str(nCycles)}
                cyclesDone.append(newCycle)
            elif line.startswith('Overall R factor'):
                if len(cyclesDone)>0:
                    currentCycle = cyclesDone[-1]
                    terms = line.split()
                    r_factor = terms[4]
                    currentCycle["r_factor"]=r_factor
            elif line.startswith('Free R factor'):
                if len(cyclesDone)>0:
                    currentCycle = cyclesDone[-1]
                    terms = line.split()
                    r_free = terms[4]
                    currentCycle["r_free"]=r_free
            elif line.startswith('Bond distances: refined atoms'):
                if len(cyclesDone)>0:
                    currentCycle = cyclesDone[-1]
                    terms = line.split()
                    rmsBondsx10 = str(10.*float(terms[5]))
                    currentCycle["rmsBondsx10"]=rmsBondsx10
                    currentCycle["rmsBonds"]= terms[5]
    
        for iCycle in range(0, len(cyclesDone)):
            cycleDone = cyclesDone[iCycle]
            cycleElement=ET.SubElement(refmacInProgressElement, 'Cycle')
            cycleNumberElement = ET.SubElement(cycleElement, "number")
            cycleNumberElement.text = cycleDone["number"]
            
            rfactorElement = ET.SubElement(cycleElement, "r_factor")
            if "r_factor" in cycleDone:
                rfactorElement.text = cycleDone["r_factor"]
            
            rfreeElement = ET.SubElement(cycleElement, "r_free")
            if "r_free" in cycleDone:
                rfreeElement.text = cycleDone["r_free"]
            
            rmsBondsx10Element = ET.SubElement(cycleElement, "rmsBondsx10")
            if "rmsBondsx10" in cycleDone:
                rmsBondsx10Element.text = cycleDone["rmsBondsx10"]
            else:
                rmsBondsx10Element.text = ""

            rmsBondsElement = ET.SubElement(cycleElement, "rmsBonds")
            if "rmsBonds" in cycleDone:
                rmsBondsElement.text = cycleDone["rmsBonds"]
            else:
                rmsBondsElement.text = ""

#=============================================================================================
import unittest
class testRefmac(unittest.TestCase):
    
    def test1(self):
        # Test creation of log file using ../test_data/test1.params.xml input
        from core.CCP4Utils import getCCP4I2Dir
        from core import CCP4Utils
        import os
        workDirectory = CCP4Utils.getTestTmpDir()
        logFile = os.path.join(workDirectory,'refmac_martin_test1.log')
        # Delete any existing log file
        if os.path.exists(logFile): os.remove(logFile)
        self.wrapper = refmac_martin(name='refmac_martin_test1',workDirectory=workDirectory)
        self.wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(),'wrappers','refmac_martin','test_data','test1.params.xml'))
        self.wrapper.setWaitForFinished(1000000)
        pid = self.wrapper.process()
        self.wrapper.setWaitForFinished(-1)
        if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
#self.assertTrue(os.path.exists(logFile),'No log file found')


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testRefmac)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
