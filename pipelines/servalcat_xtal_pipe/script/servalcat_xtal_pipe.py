from __future__ import print_function

"""
    servalcat_xtal_pipe.py: CCP4 GUI Project
    Copyright (C) 2024 University of Southampton, MRC LMB Cambridge

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

from lxml import etree
#from xml.etree import ElementTree as ET
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core import CCP4Utils
import os, sys, shutil, re
import base64

class servalcat_xtal_pipe(CPluginScript):

    TASKMODULE = 'test'
    TASKTITLE = 'Refinement - Servalcat (experimental)'
    TASKNAME = 'servalcat_xtal_pipe'
    TASKVERSION= 0.0
    WHATNEXT = ['servalcat_xtal_pipe','buccaneer_build_refine_mr','coot_rebuild','modelcraft']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 240
    MAXNJOBS = 4
    PERFORMANCECLASS = 'CRefinementPerformance'
    SUBTASKS=['Platonyzer','prosmart','servalcat']
    RUNEXTERNALPROCESS=False
    PURGESEARCHLIST =  [[ 'refmac%*/hklout.mtz', 0, "hklout" ], [ 'refmac%*/hklout.mtz', 7, "hklout" ], [ '*%*/ANOMFPHIOUT.mtz', 1, "ANOMFPHIOUT" ], [ '*%*/DIFANOMFPHIOUT.mtz', 1, "DIFANOMFPHIOUT" ]]


    ERROR_CODES = { 101 : { 'description' : 'Error copying data file from final job to pipeline directory' }
                    }

    def __init__(self, *args, **kws):
        super(servalcat_xtal_pipe, self).__init__(*args, **kws)
        self.pipelinexmlfile = self.makeFileName(format='PROGRAMXML')
        self.refmacMonitors = {}
        self.xmlroot = etree.Element("SERVALCAT") #self.xmlroot = ET.Element("RefmacOptimiseWeight")
        self.xmlroot2 = etree.Element("SERVALCAT") #self.xmlroot2 = ET.Element("RefmacOptimiseWeight")

#    def startProcess(self, processId):
#        if self.container.prosmartProtein.TOGGLE:
#            self.executeProsmartProtein()
#        else:
#            if self.container.controlParameters.WEIGHT_OPT.__str__()=='MANUAL' and self.container.controlParameters.WEIGHT.isSet(): withWeight = float(self.container.controlParameters.WEIGHT)
#            else: withWeight = -1.
#            self.executeFirstServalcat(withWeight)
#        return CPluginScript.SUCCEEDED

    def startProcess(self, processId):
        try:
           self.executeProsmartProtein()
           self.executeProsmartNucleicAcid()
           #self.executePlatonyzer()
           self.executeFirstServalcat()
        except:
           self.reportStatus(CPluginScript.FAILED)
           return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def executeProsmartProtein(self):
       if self.container.prosmartProtein.TOGGLE:
           self.prosmart_protein = self.makePluginObject('prosmart')
           self.prosmart_protein.container.inputData.TARGET_MODEL = self.container.inputData.XYZIN
           self.prosmart_protein.container.inputData.CHAINLIST_1 = self.container.prosmartProtein.CHAINLIST_1
           self.prosmart_protein.container.inputData.REFERENCE_MODELS = self.container.prosmartProtein.REFERENCE_MODELS
           self.prosmart_protein.container.controlParameters.RESTRAIN_ALL_VS_BEST = self.container.prosmartProtein.ALL_BEST
           self.prosmart_protein.container.controlParameters.RESTRAIN_SEQID = self.container.prosmartProtein.SEQID
           self.prosmart_protein.container.controlParameters.RESTRAIN_MAIN_VS_SIDE = self.container.prosmartProtein.SIDE_MAIN
           self.prosmart_protein.container.controlParameters.RESTRAIN_RMIN = self.container.prosmartProtein.RMIN
           self.prosmart_protein.container.controlParameters.RESTRAIN_RMAX = self.container.prosmartProtein.RMAX
           if self.container.prosmartProtein.ADVANCED:
              self.prosmart_protein.container.controlParameters.RESTRAIN_BFAC_FILTER = self.container.prosmartProtein.TOGGLE_BFAC
              self.prosmart_protein.container.controlParameters.RESTRAIN_BFAC_ALPHA = self.container.prosmartProtein.BFAC
              self.prosmart_protein.container.controlParameters.RESTRAIN_ALT = self.container.prosmartProtein.TOGGLE_ALT
              self.prosmart_protein.container.controlParameters.RESTRAIN_OCCUP = self.container.prosmartProtein.OCCUPANCY
              self.prosmart_protein.container.controlParameters.KEYWORDS = self.container.prosmartProtein.KEYWORDS
           self.connectSignal(self.prosmart_protein, 'finished', self.prosmartProteinFinished)
           self.prosmart_protein.waitForFinished = -1
           self.prosmart_protein.process()

    @QtCore.Slot(dict)
    def prosmartProteinFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.reportStatus(status)
            return

    def executeProsmartNucleicAcid(self):
       if self.container.prosmartNucleicAcid.TOGGLE:
           self.prosmart_nucleicacid = self.makePluginObject('prosmart')
           self.prosmart_nucleicacid.container.controlParameters.NUCLEIC_ACID = True
           self.prosmart_nucleicacid.container.inputData.TARGET_MODEL = self.container.inputData.XYZIN
           self.prosmart_nucleicacid.container.inputData.CHAINLIST_1 = self.container.prosmartNucleicAcid.CHAINLIST_1
           self.prosmart_nucleicacid.container.inputData.REFERENCE_MODELS = self.container.prosmartNucleicAcid.REFERENCE_MODELS
           self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_ALL_VS_BEST = self.container.prosmartNucleicAcid.ALL_BEST
           self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_SEQID = self.container.prosmartNucleicAcid.SEQID
           self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_MAIN_VS_SIDE = self.container.prosmartNucleicAcid.SIDE_MAIN
           self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_RMIN = self.container.prosmartNucleicAcid.RMIN
           self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_RMAX = self.container.prosmartNucleicAcid.RMAX
           if self.container.prosmartNucleicAcid.ADVANCED:
              self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_BFAC_FILTER = self.container.prosmartNucleicAcid.TOGGLE_BFAC
              self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_BFAC_ALPHA = self.container.prosmartNucleicAcid.BFAC
              self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_ALT = self.container.prosmartNucleicAcid.TOGGLE_ALT
              self.prosmart_nucleicacid.container.controlParameters.RESTRAIN_OCCUP = self.container.prosmartNucleicAcid.OCCUPANCY
              self.prosmart_nucleicacid.container.controlParameters.KEYWORDS = self.container.prosmartNucleicAcid.KEYWORDS
           self.connectSignal(self.prosmart_nucleicacid,'finished',self.prosmartNucleicAcidFinished)
           self.prosmart_nucleicacid.waitForFinished = -1
           self.prosmart_nucleicacid.process()

    @QtCore.Slot(dict)
    def prosmartNucleicAcidFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.reportStatus(status)
            return

    def executePlatonyzer(self):
       if self.container.platonyzer.TOGGLE:
          self.platonyzer = self.makePluginObject('Platonyzer')
          self.platonyzer.container.inputData.XYZIN = self.container.inputData.XYZIN
          self.platonyzer.container.controlParameters.MODE = self.container.platonyzer.MODE
          self.platonyzer.container.controlParameters.RM_VDW = self.container.platonyzer.RM_VDW
          self.connectSignal(self.platonyzer,'finished',self.platonyzerFinished)
          self.platonyzer.waitForFinished = -1
          self.platonyzer.process()

    @QtCore.Slot(dict)
    def platonyzerFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.reportStatus(status)
            return

    def executeFirstServalcat(self, withWeight=-1):
        #create wrapper
        self.firstServalcat = self.createServalcatJob(withWeight)
        # Run asynchronously ...this is needed so that commands downstream of process launch
        # (i.e. logwatcher) will be calledbefore process completion
        self.firstServalcat.doAsync = self.doAsync
        self.firstServalcat.connectSignal(self.firstServalcat,'finished',self.firstServalcatFinished)
        # Install xml node for an in progress servalcat
        self.xmlLength = 0
        # Start process
        firstServalcatXMLFilename = self.firstServalcat.makeFileName(format='PROGRAMXML')
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("executeFirstServalcat, firstServalcatXMLFilename", firstServalcatXMLFilename)
        self.watchFile(firstServalcatXMLFilename, handler=self.handleXmlChanged, minDeltaSize=34, unwatchWhileHandling=True)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        rv = self.firstServalcat.process()

    @QtCore.Slot(str)
    def handleXmlChanged2(self, xmlFilename):
        self.xmlroot2.clear()
        servalcatEtree = CCP4Utils.openFileToEtree(xmlFilename) #, useLXML=False)
        servalcatXML = servalcatEtree.xpath('//SERVALCAT')
        if len(servalcatXML) == 1:
            servalcatXML[0].tag = "SERVALCAT_WATERS"
            self.xmlroot2.append(servalcatXML[0])
        self.saveXml2()

    @QtCore.Slot(str)
    def handleXmlChanged(self, xmlFilename):
        self.xmlroot.clear()
        servalcatEtree = CCP4Utils.openFileToEtree(xmlFilename) # , useLXML=False)
        servalcatXML = servalcatEtree.xpath("//SERVALCAT")
        if len(servalcatXML) == 1:
            servalcatXML[0].tag = "SERVALCAT_FIRST"
            self.xmlroot.append(servalcatXML[0])
        self.saveXml()

    def saveXml2(self):
        newXml = etree.tostring(self.xmlroot2, pretty_print=True) # newXml = ET.tostring(self.xmlroot2)
        if len(newXml) > self.xmlLength2:
            # Get content of program.xml_first
            firstFileName = self.pipelinexmlfile + '_first'
            with open(firstFileName, 'r') as aFile:
                oldXml = etree.fromstring(aFile.read()) # oldXml = ET.fromstring(aFile.read())
            # Add content of program.xml of the last servalcat subjob
            oldXml.xpath('//SERVALCAT')[0].append(self.xmlroot2.xpath("//SERVALCAT/SERVALCAT_WATERS")[0])
            # Save as program.xml_tmo and then move to program.xml
            tmpFileName = self.pipelinexmlfile + '_tmp'
            with open(tmpFileName, 'w') as aFile:
                CCP4Utils.writeXML(aFile,etree.tostring(oldXml, pretty_print=True) ) # CCP4Utils.writeXML(aFile,ET.tostring(oldXml) )
            shutil.move(tmpFileName, self.pipelinexmlfile)
            self.xmlLength2 = len(newXml)

    def saveXml(self):
        # Save the xml if it has grown
        newXml = etree.tostring(self.xmlroot, pretty_print=True) # newXml = ET.tostring(self.xmlroot)
        if len(newXml) > self.xmlLength:
           tmpFileName = self.pipelinexmlfile + '_tmp'
           with open(tmpFileName, 'w') as aFile:
               CCP4Utils.writeXML(aFile, newXml)
           shutil.move(tmpFileName, self.pipelinexmlfile)
           self.xmlLength = len(newXml)

    def createServalcatJob(self, withWeight=-1, inputCoordinates=None, ncyc=-1):
        result = self.makePluginObject('servalcat_xtal')
        #input data for this servalcat instance is the same as the input data for the program
        result.container.inputData.copyData(self.container.inputData)
        #result.container.inputData = self.container.inputData

        #Copy over most of the control parameters
        for attr in self.container.controlParameters.dataOrder():
            if (attr != "OPTIMISE_WEIGHT"
                  and attr != "REFMAC_CLEANUP"
                  and attr != "VALIDATE_IRIS"
                  and attr != "VALIDATE_BAVERAGE"
                  and attr != "VALIDATE_RAMACHANDRAN"
                  and attr != "VALIDATE_MOLPROBITY"
                  and attr != "RUN_MOLPROBITY"):
                setattr(result.container.controlParameters, attr, getattr(self.container.controlParameters, attr))

        if self.container.prosmartProtein.TOGGLE:
            # result.container.controlParameters.PROSMART_PROTEIN_WEIGHT=self.container.prosmartProtein.WEIGHT
            result.container.controlParameters.PROSMART_PROTEIN_SGMN=self.container.prosmartProtein.SGMN
            result.container.controlParameters.PROSMART_PROTEIN_SGMX=self.container.prosmartProtein.SGMX
            result.container.controlParameters.PROSMART_PROTEIN_ALPHA=self.container.prosmartProtein.ALPHA
            result.container.controlParameters.PROSMART_PROTEIN_DMAX=self.container.prosmartProtein.DMAX
            result.container.inputData.PROSMART_PROTEIN_RESTRAINTS=self.prosmart_protein.container.outputData.RESTRAINTS
        if self.container.prosmartNucleicAcid.TOGGLE:
            # result.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT=self.container.prosmartNucleicAcid.WEIGHT
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMN=self.container.prosmartNucleicAcid.SGMN
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMX=self.container.prosmartNucleicAcid.SGMX
            result.container.controlParameters.PROSMART_NUCLEICACID_ALPHA=self.container.prosmartNucleicAcid.ALPHA
            result.container.controlParameters.PROSMART_NUCLEICACID_DMAX=self.container.prosmartNucleicAcid.DMAX
            result.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS=self.prosmart_nucleicacid.container.outputData.RESTRAINTS
        """if self.container.platonyzer.TOGGLE:
            result.container.inputData.PLATONYZER_RESTRAINTS=self.platonyzer.container.outputData.RESTRAINTS
            result.container.inputData.XYZIN=self.platonyzer.container.outputData.XYZOUT"""

        #specify weight if a meaningful one hsa been offered
        if withWeight>=0.:
            result.container.controlParameters.WEIGHT_OPT='MANUAL'
            result.container.controlParameters.WEIGHT = withWeight

        if inputCoordinates is not None:
            result.container.inputData.XYZIN.set(inputCoordinates)
            if ncyc > 0:
                result.container.controlParameters.NCYCLES.set(ncyc)
        return result

    @QtCore.Slot(dict)
    def firstServalcatFinished(self, statusDict):
        print("AAA1")
        if statusDict['finishStatus'] == CPluginScript.UNSATISFACTORY:
            print("AAA1.UNSATISFACTORY")
            import os
            if os.path.isfile(self.firstServalcat.container.outputData.LIBOUT.__str__()):
                from wrappers.acedrg.script import acedrg
                try:
                    rdkitMol = acedrg.molFromDict(self.firstServalcat.container.outputData.LIBOUT.__str__())
                    from rdkit import Chem
                    molRemovedHs = Chem.RemoveHs(rdkitMol)
                    svgXml = acedrg.svgFromMol(molRemovedHs)
                    self.xmlroot.append(svgXml)
                except:
                    print('Unable to generate svg from DICT')
                shutil.copyfile(self.firstServalcat.container.outputData.LIBOUT.__str__(), self.container.outputData.LIBOUT.__str__())
                self.container.outputData.LIBOUT.annotation = 'Refmac-generated library...use with caution'
            if os.path.isfile(self.firstServalcat.container.outputData.PSOUT.__str__()):
                shutil.copyfile(self.firstServalcat.container.outputData.PSOUT.__str__(), self.container.outputData.PSOUT.__str__())
                self.container.outputData.PSOUT.annotation.set('Pictures of ligand prepared by refmac')
            with open(self.makeFileName('PROGRAMXML'), 'w') as programXML:
                CCP4Utils.writeXML(programXML, etree.tostring(self.xmlroot, pretty_print=True)) # CCP4Utils.writeXML(programXML, ET.tostring(self.xmlroot))
            self.reportStatus(CPluginScript.UNSATISFACTORY)

        elif self.firstServalcat.errorReport.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            print("AAA1.MAXSEVERITY")
            #This gets done in thefirstServalcat.reportStatus() - Liz
            #self.extendErrorReport(self.firstServalcat.errorReport)
            try:
                servalcatEtree = CCP4Utils.openFileToEtree(self.firstServalcat.makeFileName('PROGRAMXML'))
                servalcatXML = servalcatEtree.xpath('//SERVALCAT') # servalcatXML = servalcatEtree.findall('.//SERVALCAT')
                if len(servalcatXML) == 1: self.xmlroot.append(servalcatXML[0])
            except:
                print('Failed attempt to read XML file from first Refmac')
            try:
                newXml = etree.tostring(self.xmlroot,pretty_print=True) # newXml = ET.tostring(self.xmlroot)
                aFile = open(self.pipelinexmlfile, 'w')
                CCP4Utils.writeXML(aFile, newXml)
                aFile.close()
            except:
               print('Failed attempt to write pipeline XML file')
            self.reportStatus(CPluginScript.FAILED)
            return
        print("AAA11")

        self.handleXmlChanged(self.firstServalcat.makeFileName(format='PROGRAMXML'))

        print("AAA10")
        import os
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            # This gets done in the firstServalcat.reportStatus() - Liz
            self.fileSystemWatcher = None
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            print("AAA12")
            # self.addCycleXML(self.firstServalcat) # MM
            aFile=open( self.pipelinexmlfile,'w')
            CCP4Utils.writeXML(aFile, etree.tostring(self.xmlroot, pretty_print=True) ) # CCP4Utils.writeXML(aFile, ET.tostring(self.xmlroot) )
            aFile.close()
            print("AAA13")
            if self.container.controlParameters.OPTIMISE_WEIGHT:
                print("AAA14")
                self.fileSystemWatcher = None
                weightUsed = float(self.xmlroot.xpath('//weight')[-1].text) # weightUsed = float(self.xmlroot.findall('.//weight')[-1].text)
                self.tryVariousRefmacWeightsAround(weightUsed)
            else:
               print("AAA15")
               best_r_free = self.firstServalcat.container.outputData.PERFORMANCEINDICATOR.RFactor
               print("AAA15.1")
               if self.container.controlParameters.ADD_WATERS: # and best_r_free < self.container.controlParameters.REFPRO_RSR_RWORK_LIMIT :
                   # Coot sujob to add waters
                   print("AAA16")
                   self.currentCoordinates = self.firstServalcat.container.outputData.XYZOUT
                   self.cootPlugin = self.makeCootPlugin()
                   self.cootPlugin.doAsync = self.doAsync
                   self.cootPlugin.connectSignal(self.cootPlugin, 'finished', self.cootFinished)
                   print("AAA17")
                   rv = self.cootPlugin.process()
                   if rv == CPluginScript.FAILED:
                        self.reportStatus(rv)
               else:
                     self.fileSystemWatcher = None
                     self.finishUp(self.firstServalcat)
        print('done servalcat_xtal_pipe.firstServalcatFinished')
        return

    def makeCootPlugin(self):
         # FIXME - This is all nonsense - needs to consider best task, etc... *NOT* just firstServalcata?
        cootPlugin = self.makePluginObject('coot_script_lines')
        xyzinList = cootPlugin.container.inputData.XYZIN
        xyzinList.append(xyzinList.makeItem())
        xyzinList[-1].set(self.currentCoordinates)
        fphiinList = cootPlugin.container.inputData.FPHIIN
        fphiinList.append(fphiinList.makeItem())
        fphiinList[-1].set(self.firstServalcat.container.outputData.FPHIOUT)
        #coot_stepped_refine,coot_fit_residues,coot_script_lines
        if self.container.controlParameters.ADD_WATERS:
            cootPlugin.container.controlParameters.SCRIPT = '''coot.set_ligand_water_to_protein_distance_limits(2.2, 3.3)
coot.execute_find_waters_real(MapHandle_1,MolHandle_1,0,1.3)
coot.write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        """
        rso = self.container.controlParameters.REFPRO_COOT_REALSPACE_OPERATION.__str__()
        if rso == "coot_script_lines":
            cootPlugin.container.controlParameters.SCRIPT = self.container.controlParameters.SCRIPT
        elif rso == "coot_fit_residues":
            cootPlugin.container.controlParameters.SCRIPT = '''fill_partial_residues(MolHandle_1)
fit_protein(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif rso == "coot_stepped_refine":
            if self.container.controlParameters.USERAMA.isSet() and self.container.controlParameters.USERAMA:
                cootPlugin.container.controlParameters.SCRIPT = '''fill_partial_residues(MolHandle_1)
stepped_refine_protein_for_rama(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
            else:
                cootPlugin.container.controlParameters.SCRIPT = '''fill_partial_residues(MolHandle_1)
stepped_refine_protein(MolHandle_1)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif rso == "coot_add_waters":
            cootPlugin.container.controlParameters.SCRIPT = '''set_ligand_water_to_protein_distance_limits(2.2, 3.3)
execute_find_waters_real(MapHandle_1,MolHandle_1,0,1.3)
write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        elif rso == "none":
            cootPlugin.container.controlParameters.SCRIPT = '''write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))
'''
        """
        return cootPlugin

    def mapVerdictSuggestionsToi2Params(self,suggestedParameters):

        newSuggestions = {}

        for k,v in suggestedParameters.items():
            if k=="NCYC":
                newSuggestions["NCYCLES"] = v
            elif k=="JELLY_DMAX":
                newSuggestions["JELLY_DIST"] = v
            elif k=="JELLY_SIGMA":
                newSuggestions["JELLY_SIGMA"] = v
            elif k=="WAUTO_VAL":
                newSuggestions["WEIGHT"] = v
            elif k=="TLS_CYCLES":
                newSuggestions["NTLSCYCLES"] = v
            elif k=="BFAC":
                newSuggestions["B_REFINEMENT_MODE"] = v
            elif k=="VDW_VAL":
                #newSuggestions["VDWRESTRAINTS"] = v #FIXME - This option does not currently exist in i2.
                pass
            elif k=="JELLY":
                if v.lower() == "yes":
                    newSuggestions["USE_JELLY"] = "True"
                else:
                    newSuggestions["USE_JELLY"] = "False"
            elif k=="TLS":
                if v.lower() == "auto":
                    newSuggestions["TLSMODE"] = "AUTO"
            elif k=="NCSR":
                if v.lower() == "yes":
                    newSuggestions["USE_NCS"] = "True"
                else:
                    newSuggestions["USE_NCS"] = "False"
            elif k=="RESET_B":
                if v.lower() == "yes":
                    newSuggestions["BFACSETUSE"] = "True" #According to Oleg, we should not get here currently (14/09/2021).
                else:
                    newSuggestions["BFACSETUSE"] = "False"
            elif k=="WAUTO_YES":
                if v.lower() == "yes":
                    newSuggestions["WEIGHT_OPT"] = "AUTO"
                else:
                    newSuggestions["WEIGHT_OPT"] = "MANUAL"
            elif k=="MKHYDR":
                if v.lower() == "all":
                    newSuggestions["HYDR_ALL"] = "ALL"
                    newSuggestions["HYDR_USE"] = "True"

        return newSuggestions

    def checkFinishStatus( self,statusDict,failedErrCode,outputFile = None,noFileErrCode= None):
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        if outputFile is not None and not outputFile.exists():
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile))
            self.reportStatus(CPluginScript.FAILED)
        return

    @QtCore.Slot(dict)
    def cootFinished(self, statusDict={}):
        import functools
        # Check coot status and start servalcat
        if len(self.cootPlugin.container.outputData.XYZOUT) == 0:
            self.appendErrorReport(205,'Coot failed to produce an output file')
            self.reportStatus(CPluginScript.FAILED)
        self.checkFinishStatus(statusDict=statusDict,failedErrCode=204, outputFile= self.cootPlugin.container.outputData.XYZOUT[0],noFileErrCode=205)
        try:
          if self.container.controlParameters.ADD_WATERS:
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = etree.fromstring(aFile.read()) # oldXml = ET.fromstring(aFile.read())
            aFile.close()
            nwaters = "unknown"
            cootLogTxt = os.path.join(os.path.dirname(self.cootPlugin.container.outputData.XYZOUT[0].__str__()), "log.txt")
            with open(cootLogTxt, 'r') as f:
               for l in f:
                   if l.startswith("INFO::") and "found" in l and "water fitting" in l:
                      nwaters = l.strip()
                      numsearch = [ x for x in nwaters.split() if x.isdigit() ]
                      if len(numsearch)>0:
                         nwaters = numsearch[0]
                      break
            postRefmacCoot = etree.Element("CootAddWaters") # postRefmacCoot = ET.Element("CootAddWaters")
            postRefmacCoot.text = "Coot added " + nwaters + " water molecules."
            oldXml.append(postRefmacCoot)
            aFile = open(self.pipelinexmlfile+'_tmpcoot','w')
            CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True)) # CCP4Utils.writeXML(aFile,ET.tostring(oldXml))
            aFile.close()
            shutil.move(self.pipelinexmlfile+'_tmpcoot', self.pipelinexmlfile)
          self.cootPlugin.container.outputData.XYZOUT[0].subType = 1
          self.currentCoordinates = self.cootPlugin.container.outputData.XYZOUT[0]
          self.servalcatPostCootPlugin = self.createServalcatJob(inputCoordinates=self.currentCoordinates, ncyc=int(self.container.controlParameters.NCYCLES_AFTER_ADD_WATERS))
          self.servalcatPostCootPlugin.doAsync = True
          self.servalcatPostCootPlugin.connectSignal(self.servalcatPostCootPlugin,'finished',functools.partial(self.postCootRefmacFinished,self.servalcatPostCootPlugin))
          servalcatPostCootXMLFilename = self.servalcatPostCootPlugin.makeFileName(format='PROGRAMXML')
          self.xmlLength2 = 0
          self.watchFile(servalcatPostCootXMLFilename, handler=self.handleXmlChanged2, minDeltaSize=34, unwatchWhileHandling=True)
          shutil.copyfile(self.pipelinexmlfile, self.pipelinexmlfile+"_first")
          rv = self.servalcatPostCootPlugin.process()
          if rv == CPluginScript.FAILED: self.reportStatus(rv)
        except Exception as e:
          self.appendErrorReport(CPluginScript,39,str(e))

    @QtCore.Slot('CPluginScript', dict)
    def postCootRefmacFinished(self, servalcatJob, statusDict={}):
        self.handleXmlChanged2(servalcatJob.makeFileName(format='PROGRAMXML'))

        import os
        self.fileSystemWatcher = None
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            #This gets done in the firstServalcat.reportStatus() - Liz
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            # self.addCycleXML(servalcatJob,"servalcatPostCoot") # MM
            newXml = etree.tostring(self.xmlroot,pretty_print=True) # newXml = ET.tostring(self.xmlroot)
            servalcatPostCoot = self.xmlroot.xpath("servalcatPostCoot") # servalcatPostCoot = self.xmlroot.findall("servalcatPostCoot")
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = etree.fromstring(aFile.read()) # oldXml = ET.fromstring(aFile.read())
            aFile.close()
            oldXml.extend(servalcatPostCoot)
            aFile = open(self.pipelinexmlfile,'w')
            CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True)) # CCP4Utils.writeXML(aFile,ET.tostring(oldXml))
            aFile.close()

        self.finishUp(servalcatJob)

    def finishUp(self, servalcatJob):
        from core import CCP4ProjectsManager
        print('into servalcat_xtal_pipe.finishUp')
        for attr in self.container.outputData.dataOrder():
            print('servalcat_xtal_pipe.finishUp attr',attr)
            wrappersAttr = getattr(servalcatJob.container.outputData, attr)
            pipelinesAttr = getattr(self.container.outputData, attr)
            if attr in ["PERFORMANCEINDICATOR"]:
                setattr(self.container.outputData, attr, wrappersAttr)
            else:
                if os.path.exists(str(wrappersAttr.fullPath)):
                  try:
                    shutil.copyfile( str(wrappersAttr.fullPath), str(pipelinesAttr.fullPath) )
                  except:
                    self.appendErrorReport(101,str(wrappersAttr.fullPath)+' to '+str(pipelinesAttr.fullPath))
                  #self.container.outputData.copyData(servalcatJob.container.outputData,[attr])
                if attr == "XMLOUT":
                    pass

        print('servalcat_xtal_pipe.finishUp 1')
        from core import CCP4XtalData
        # Apply database annotations
        self.container.outputData.XYZOUT.annotation.set('Model from refinement (PDB format)')
        self.container.outputData.CIFFILE.annotation.set('Model from refinement (mmCIF format)')
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from refinement'
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from refinement'
        self.container.outputData.DIFFPHIOUT.subType = 2
        # self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        # self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        self.container.outputData.TLSOUT.annotation = 'TLS parameters from refinement'
        self.container.outputData.COOTSCRIPTOUT.annotation = 'Coot script written from refinement'
        self.container.outputData.ANOMFPHIOUT.annotation = 'Weighted anomalous difference map from refinement'
        self.container.outputData.DIFANOMFPHIOUT.annotation = 'Weighted differences of anomalous difference map'

        print('servalcat_xtal_pipe.finishUp 2')
        if self.container.outputData.DICT.exists():
           self.container.outputData.DICT.annotation = 'Accumulated ligand geometry dictionary'
        if self.container.outputData.LIBOUT.exists():
          annotation = 'Refmac generated geometry'
          try:
            print('LIBOUT monomerList',self.container.outputData.LIBOUT.fileContent.monomerList)
            if len(self.container.outputData.LIBOUT.fileContent.monomerList)>0:
              annotation = 'Refmac generated geometry for:'
              for item in self.container.outputData.LIBOUT.fileContent.monomerList:
                annotation = annotation + ' ' + str(item.three_letter_code)
              ligxml = etree.SubElement(self.xmlroot,"LIGANDS")
              for item in self.container.outputData.LIBOUT.fileContent.monomerList:
                ligNode = etree.SubElement(ligxml,"ligand")
                ligNode.text = str(item.three_letter_code)
              self.saveXml()
          except:
              print('Error creating LIBOUT annotation')
              self.container.outputData.LIBOUT.annotation = annotation
          try:
              self.mergeDictToProjectLib(fileName=self.container.outputData.LIBOUT.__str__())
          except:
              print('Error merging library to Project Dictionary')
        print('servalcat_xtal_pipe.finishUp 3'); sys.stdout.flush()
#
        cleanUpIntermediate = False
        if hasattr(self.container.controlParameters,"REFMAC_CLEANUP"):
            cleanUpIntermediate = self.container.controlParameters.REFMAC_CLEANUP
            if cleanUpIntermediate:
                print('servalcat_xtal_pipe.finishUp 4'); sys.stdout.flush()
                cleanup = CCP4ProjectsManager.CPurgeProject(self.firstServalcat._dbProjectId)
                print('servalcat_xtal_pipe.finishUp 5'); sys.stdout.flush()
                cleanup.purgeJob(self.firstServalcat.jobId,context="extended_intermediate",reportMode="skip")

                if hasattr(self,"servalcatPostCootPlugin"):
                    print('servalcat_xtal_pipe.finishUp 6'); sys.stdout.flush()
                    cleanup = CCP4ProjectsManager.CPurgeProject(self.servalcatPostCootPlugin._dbProjectId)
                    print('servalcat_xtal_pipe.finishUp 7'); sys.stdout.flush()
                    cleanup.purgeJob(self.servalcatPostCootPlugin.jobId,context="extended_intermediate",reportMode="skip")

        self.reportStatus(CPluginScript.SUCCEEDED)
        return # MM
        if not str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID':
           #Geometry validation not performed in rigid body mode

           runMolprobity = False
           if hasattr(self.container.controlParameters,"RUN_MOLPROBITY"):
              runMolprobity = self.container.controlParameters.RUN_MOLPROBITY
              if runMolprobity:
                 #MN: First order attempt at providing molprobity analysis
                 try:
                     print("Attempting molprobity run after refinement...")
                     from mmtbx.command_line import molprobity
                     coordPath = self.container.outputData.XYZOUT.fullPath.__str__()
                     fileRoot, fileExt = os.path.splitext(coordPath)
                     sanitizedCoordPath = fileRoot + "+asPDB.pdb"

                     #Use mmdb to do some sanitization
                     import ccp4mg
                     import mmdb2
                     mmdb2.InitMatType()
                     m = mmdb2.Manager()

                     #This line ought to just do it, but seems not to...
                     m.SetFlag(mmdb2.MMDBF_IgnoreSegID)
                     m.ReadCoorFile(coordPath)

                     #Remove any SEGIDs (which can confuse molprobity if heterogenous)
                     sel = m.NewSelection()
                     m.SelectAtoms(sel, 0,"*",mmdb2.ANY_RES,"*",mmdb2.ANY_RES,"*","*","*","*","*",mmdb2.SKEY_OR )
                     selindexp = mmdb2.intp()
                     selAtoms = mmdb2.GetAtomSelIndex(m,sel,selindexp)
                     nSelAtoms = selindexp.value()
                     # ... but this certainly does.
                     for i in range(nSelAtoms):
                         at = mmdb2.getPCAtom(selAtoms,i)
                         at.segID = b"    "
                     m.FinishStructEdit()

                     #Looking at the molprobity output, it seems to me that it tres to handle alternates properly...
                     #I am not sure that the following code is appropriate, so am commenting out MN
                     '''
                     #Delete alternate locations using mmdb logic
                     sel = m.NewSelection()
                     m.SelectAtoms(sel, 0,"*",mmdb2.ANY_RES,"*",mmdb2.ANY_RES,"*","*","*","*","! ",mmdb2.SKEY_OR )
                     selindexp = mmdb2.intp()
                     selAtoms = mmdb2.GetAtomSelIndex(m,sel,selindexp)
                     nSelAtoms = selindexp.value()
                     # Edit out the SEGIDs - should have happened at read time ?
                     residuesWithAltes = set()
                     for i in range(nSelAtoms):
                         at = mmdb2.getPCAtom(selAtoms,i)
                         residuesWithAltes.add(at.GetResidue())
                     for residue in residuesWithAltes:
                         residue.DeleteAltLocs()
                     m.FinishStructEdit()
                     '''
                     m.WritePDBASCII(sanitizedCoordPath)

                     fileRoot = os.path.join(self.getWorkDirectory(),"molprobity")
                     molprobity.run(["input.pdb.file_name={}".format(sanitizedCoordPath),
                                     "output.prefix={}".format(fileRoot)])
                     mpCootScriptPath = fileRoot+"_coot.py"
                     mpProbePath = fileRoot+"_probe.txt"
                     if os.path.isfile(mpCootScriptPath):
                         with (open(self.container.outputData.COOTSCRIPTOUT.fullPath.__str__(),"a+")) as cootScript:
                             cootScript.write("\n")
                             with open(mpCootScriptPath,"r") as mpCootScript:
                                 content = mpCootScript.read()
                                 content = content.replace('"molprobity_probe.txt"',
                                                           '"{}"'.format(mpProbePath))
                                 cootScript.write(content)

                     mpOutPath = fileRoot+".out"
                     if os.path.isfile(mpOutPath):
                         etree.SubElement(etree.SubElement(self.xmlroot,"Molprobity"), "Output").text = etree.CDATA(open(mpOutPath).read())
                     self.saveXml()
                     print("...Succeeded molprobity run after refinement :-)")
                 except Exception as err:
                     etree.SubElement(etree.SubElement(self.xmlroot,"Molprobity"), "Output").text = etree.CDATA(str(err))
                     self.saveXml()
                     print("...Failed molprobity run after refinement :-(", err)

           validate_iris = False
           if hasattr(self.container.controlParameters,"VALIDATE_IRIS"):
               validate_iris = self.container.controlParameters.VALIDATE_IRIS

           validate_baverage = False
           if hasattr(self.container.controlParameters,"VALIDATE_BAVERAGE"):
               validate_baverage = self.container.controlParameters.VALIDATE_BAVERAGE

           validate_ramachandran = False
           if hasattr(self.container.controlParameters,"VALIDATE_RAMACHANDRAN"):
              validate_ramachandran = self.container.controlParameters.VALIDATE_RAMACHANDRAN

           validate_molprobity = False
           if hasattr(self.container.controlParameters,"VALIDATE_MOLPROBITY"):
              validate_molprobity = self.container.controlParameters.VALIDATE_MOLPROBITY

           # if validate_baverage or validate_molprobity or validate_ramachandran:
           if True:    
             try:
                print("\n\n\n++++1") 
                self.validate = self.makePluginObject('validate_protein')
                self.validate.container.inputData.XYZIN_1 = self.container.outputData.XYZOUT
                self.validate.container.inputData.F_SIGF_1 = self.container.inputData.HKLIN # I ???
                self.validate.container.inputData.XYZIN_2 = self.container.outputData.XYZOUT
                self.validate.container.inputData.F_SIGF_2 = self.container.inputData.HKLIN
                # self.validate.container.inputData.XYZIN = self.container.inputData.XYZIN
                # self.validate.container.inputData.XYZIN_1.setFullPath(coordPath)
                # self.validate.container.inputData.F_SIGF_1 = self.container.inputData.HKLIN
                #MN...Using "="" to set this is an odd thing and breaks under some circumstances.
                #Specifically, the path for this is now in the prosmart_refmac (i.e. parent job) directory,
                #but it is an output of the validate_protein subjob.  When that subjob seeks to save it to the
                #database, it finds that it is in the wrong directory and croaks.  I have tried to get round this
                #by flagging it as saveToDb=False in the validate_protein .def.xml. But htis needs consideration by jon
                #Agirre and Rob Nicholls
                self.validate.container.outputData.COOTSCRIPTOUT = self.container.outputData.COOTSCRIPTOUT
                # self.validate.container.outputData.COOTSCRIPTOUT.setFullPath("/tmp/coot.tmp")
                print("\n\n\n++++2")

                """self.validate.container.controlParameters.DO_IRIS = validate_iris
                self.validate.container.controlParameters.DO_BFACT = validate_baverage
                self.validate.container.controlParameters.DO_RAMA = validate_ramachandran
                self.validate.container.controlParameters.DO_MOLPROBITY = validate_molprobity"""
                self.validate.container.controlParameters.DO_IRIS = True
                self.validate.container.controlParameters.DO_BFACT = True
                self.validate.container.controlParameters.DO_RAMA = True
                self.validate.container.controlParameters.DO_MOLPROBITY = True

                print(str(validate_baverage), str(validate_ramachandran), str(validate_molprobity))
                self.validate.doAsync = False
                self.validate.waitForFinished = -1
                print("\n\n\n++++3")
                self.validate.process()
                print("\n\n\n++++4")

                validateXMLPath = self.validate.makeFileName('PROGRAMXML')
                validateXML = CCP4Utils.openFileToEtree(validateXMLPath)
                xml_validation = etree.SubElement(self.xmlroot,"Validation")
                if len(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info"))>0:
                   xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info")[0]) 
                """
                if self.validate.container.controlParameters.DO_IRIS:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/Iris"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Iris")[0]) 
                """
                if self.validate.container.controlParameters.DO_BFACT:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/B_factors"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/B_factors")[0])
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/B_averages"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/B_averages")[0])
                if self.validate.container.controlParameters.DO_RAMA:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/Ramachandran"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Ramachandran")[0])
                if self.validate.container.controlParameters.DO_MOLPROBITY:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/Molprobity"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Molprobity")[0])
                   self.saveXml()
                   try:
                       from . import prosmart_refmac_verdict
                       programxml = self.makeFileName('PROGRAMXML')                  #"/Users/stuart/CCP4I2_PROJECTS/5_7_2021/CCP4_JOBS/job_19/program.xml"
                       pdbfile = self.container.outputData.XYZOUT.fullPath.__str__() #"/Users/stuart/CCP4I2_PROJECTS/5_7_2021/CCP4_JOBS/job_19/19_5_7_2021_xyzout_prosmart_refmac.pdb"
                       if hasattr(self,"servalcatPostCootPlugin"):
                           servalcatlog = self.servalcatPostCootPlugin.makeFileName('LOG') #"/Users/stuart/CCP4I2_PROJECTS/5_7_2021/CCP4_JOBS/job_19/job_1/log.txt"
                       else:
                           servalcatlog = self.firstServalcat.makeFileName('LOG')

                       verdict_result = prosmart_refmac_verdict.getJSCOFERefmac5Verdict(programxml=programxml,pdbfile=pdbfile,refmaclog=servalcatlog)
                       verdict_score = verdict_result["score"]
                       verdict_message  = verdict_result["message"]
                       bottomline = verdict_result["bottomLine"]
                       meanRfree = verdict_result["meanRfree"]
                       medianClash = verdict_result["medianClash"]
                       ramaOutliers = verdict_result["ramaOutliers"]
                       suggestedParameters = verdict_result["suggestedParameters"]

                       suggestedParameters = self.mapVerdictSuggestionsToi2Params(suggestedParameters)

                       xml_verdict = etree.SubElement(self.xmlroot,"Verdict")
                       xml_verdict_score = etree.SubElement(xml_verdict,"verdict_score")
                       xml_verdict_score.text = str(verdict_score)
                       xml_verdict_message = etree.SubElement(xml_verdict,"verdict_message")
                       xml_verdict_message.text = etree.CDATA(verdict_message)
                       xml_bottomline = etree.SubElement(xml_verdict,"bottomline")
                       xml_bottomline.text = etree.CDATA(bottomline)
                       xml_meanRfree = etree.SubElement(xml_verdict,"meanRfree")
                       xml_meanRfree.text = str(meanRfree)
                       xml_medianClash = etree.SubElement(xml_verdict,"medianClash")
                       xml_medianClash.text = str(medianClash)
                       xml_ramaOutliers = etree.SubElement(xml_verdict,"ramaOutliers")
                       xml_ramaOutliers.text = str(ramaOutliers)
                       xml_suggestedParameters = etree.SubElement(xml_verdict,"suggestedParameters")
                       for k,v in suggestedParameters.items():
                          xml_suggestedParameters_k = etree.SubElement(xml_suggestedParameters,k)
                          xml_suggestedParameters_k.text = str(v)

                       self.saveXml()
                   except:
                       import traceback
                       print("Some problem with verdict...."); sys.stdout.flush()
                       exc_type, exc_value, exc_tb = sys.exc_info()[:3]
                       sys.stderr.write(str(exc_type) + '\n')
                       sys.stderr.write(str(exc_value) + '\n')
                       traceback.print_tb(exc_tb)
                self.saveXml()
             except Exception as err:
                import traceback
                traceback.print_exc()
                print("...Failed validation run after refinement", err)

        logfiles = []
        if hasattr(self,"firstServalcat"):
            logfiles.append(self.firstServalcat.makeFileName('LOG'))
        if hasattr(self,"servalcatPostCootPlugin"):
            logfiles.append(self.servalcatPostCootPlugin.makeFileName('LOG'))

        # self.createWarningsXML(logfiles) # MM
        self.saveXml()

        print('done servalcat_xtal_pipe.finishUp'); sys.stdout.flush()
        self.reportStatus(CPluginScript.SUCCEEDED)

    def tryVariousRefmacWeightsAround(self, weight):
        import math
        print('Generating jobs with weights around ', weight)
        #make an array to hold the child-jobs
        refmacJobs = []
        for factorExponent in range(-3,4):
            if factorExponent != 0:
                factor = math.pow(2., factorExponent)
                refmacJobs.append(self.createServalcatJob(float(factor)*float(weight)))
                # Set to run asynchronously and set a callback
        self.submitBunchOfJobs(refmacJobs)

    def submitBunchOfJobs(self, jobs):
        self.jobsToSubmit = []
        self.jobsInTrain = {}
        self.jobsCompleted = []
        for job in jobs:
            job.doAsync  = True
            job.connectSignal(job,'finished',self.handleDone)
            self.jobsToSubmit.append(job)
        print('ready to submit from list of length ',len(self.jobsToSubmit))

        for job in self.jobsToSubmit:
            if len(self.jobsInTrain) < prosmart_refmac.MAXNJOBS:
                self.submitJob(job)

    def submitJob(self,job):
        rv = job.process()
        #The mtzdump instance must be saved to keep it in scope and everything else can be got from that.
        self.jobsInTrain[str(job.processId)]=job
        self.jobsToSubmit.remove(job)
        print('submitted job 0')

    @QtCore.Slot(dict)
    def handleDone(self, ret):
        pid =  ret.get('pid',None)
        status,exitStatus,exitCode = self.postProcessCheck(pid)
        if  status == CPluginScript.FAILED:
            self.reportStatus(status)
            return
        import sys
        # callback is passed the jobId (=Non
        # if not in ccp4i2-db context) and processId that
        # can serve at identifier for subProcess
        import sys
        import time
        from copy import deepcopy

        #print 'demo_multi_mtzdump.handleDone',ret

        rtask = self.jobsInTrain[str(pid)]

        # self.addCycleXML(rtask) # MM
        # Save the xml on every cycle
        aFile=open( self.pipelinexmlfile,'w')
        CCP4Utils.writeXML(aFile, etree.tostring(self.xmlroot,pretty_print=True) ) # CCP4Utils.writeXML(aFile, ET.tostring(self.xmlroot) )
        aFile.close()

        #decrement count of running jobs
        self.jobsCompleted.append(rtask)
        del self.jobsInTrain[str(pid)]

        if len(self.jobsToSubmit) != 0:
            job = self.jobsToSubmit[0]
            self.submitJob(job)

        elif len(self.jobsInTrain)==0:
            best_r_free = 9999.
            self.best_rtask = 0
            for rtask in self.jobsCompleted:
                if float(rtask.container.outputData.PERFORMANCEINDICATOR.RFree) < best_r_free:
                    best_r_free = float(rtask.container.outputData.PERFORMANCEINDICATOR.RFree)
                    self.best_rtask = rtask
            self.finishUp(self.best_rtask)
        return

    """def addCycleXML(self, rtask, name="RefmacWeight"):
        xmlcyc = ET.SubElement(self.xmlroot,name)
        cycWeightNode = ET.SubElement(xmlcyc,"weight")
        _rxml = CCP4Utils.openFileToEtree(rtask.makeFileName('PROGRAMXML'))
        rxml=ET.Element('ccp4i2root')
        rxml.append(_rxml)
        try: cycWeightNode.text = rxml.findall(".//Cycle/WeightUsed")[-1].text
        except: pass
        rstats = rxml.findall(".//REFMAC")
        xmlcyc.append(rstats[0])"""

    """def addCycleXML(self, rtask, name="RefmacWeight"):
        # TO DO MM AVOID THIS!
        xmlcyc = etree.SubElement(self.xmlroot,name)
        cycWeightNode = etree.SubElement(xmlcyc,"weight")
        return
        rxml = CCP4Utils.openFileToEtree(rtask.makeFileName('PROGRAMXML'))
        try: cycWeightNode.text = rxml.xpath("//Cycle/WeightUsed")[-1].text
        except: pass
        rstats = rxml.xpath("//REFMAC") # rstats = rxml.findall(".//REFMAC")
        xmlcyc.append (rstats[0])"""

    def handleTimeout(self):
        import sys;sys.stdout.flush()

        for rtask in self.jobsInTrain:
            print('TERMINATING', rtask.processId,sys.stdout.flush())
            try:
                rtask.terminate()
            except:
                pass

        self.appendErrorReport(40,str(self.TIMEOUT_PERIOD))
        self.reportStatus(CPluginScript.FAILED)

def coefficientsToMap(coefficientsPath, mapPath=None, overSample=1.0):
    import clipper
    mtz_file = clipper.CCP4MTZfile()
    hkl_info = clipper.HKL_info()
    mtz_file.open_read (str(coefficientsPath))
    mtz_file.import_hkl_info ( hkl_info )
    sg, cell = hkl_info.spacegroup(), hkl_info.cell()
    fphidata = clipper.HKL_data_F_phi_float(hkl_info)
    mtz_file.import_hkl_data( fphidata, str("/*/*/[F,PHI]") );
    #mtz_file.import_hkl_data( fphidata, str("/*/*/[FC,PHFC]") );
    mtz_file.close_read()
    #Clipper will sample the output map according to Fourier theory and hte nominal resolution
    #for visualisation, it is generally nicer to make things a bit more finely sampled
    fudgedResolution = hkl_info.resolution()
    fudgedResolution.init(hkl_info.resolution().limit()/overSample)
    mygrid=clipper.Grid_sampling ( hkl_info.spacegroup(), hkl_info.cell(), fudgedResolution )
    mymap = clipper.Xmap_float(hkl_info.spacegroup(), hkl_info.cell(), mygrid )
    mymap.fft_from_float(fphidata)

    mapout = clipper.CCP4MAPfile()
    if mapPath is None:
        coefficientsRoot, extension = os.path.splitext(os.path.abspath(coefficientsPath))
        mapPath = coefficientsRoot+".map"

    mapout.open_write( mapPath )
    mapout.export_xmap_float( mymap )
    mapout.close_write()
    return mapPath

# Function called from gui to support exporting MTZ files
def exportJobFile(jobId=None,mode=None,fileInfo={}):
    import os
    from core import CCP4Modules
    from dbapi.CCP4DbApi import CDbApi
    from dbapi.CCP4DbApi import FILE_ROLE_OUT

    theDb = CCP4Modules.PROJECTSMANAGER().db()
    if mode == 'complete_mtz':
        #print 'refmac.exportJobFile',mode
        childJobs = theDb.getChildJobs(jobId=jobId,details=True)
        #print 'exportJobFile childJobs',childJobs
        if childJobs[-1][2] == 'refmac':
          jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=childJobs[-1][1],create=False)
          if os.path.exists(os.path.join(jobDir,'hklout.mtz')):
             return  os.path.join(jobDir,'hklout.mtz')
        elif childJobs[-1][2] == 'validate_protein':
          if childJobs[-2][2] == 'refmac':
             jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=childJobs[-2][1],create=False)
             if os.path.exists(os.path.join(jobDir,'hklout.mtz')):
                return  os.path.join(jobDir,'hklout.mtz')

    elif mode == '2FoFc_as_map' or mode == 'FoFc_as_map':
        files = theDb.getJobFiles(jobId=jobId, role=FILE_ROLE_OUT, searchFileUses=True, fileTypes=[13])
        for fileId in files:
            fileInfo = theDb.getFileInfo(fileId=fileId, mode='jobparamname')
            if (fileInfo == 'FPHIOUT' and mode == '2FoFc_as_map') or (fileInfo == 'DIFFPHIOUT' and mode == 'FoFc_as_map'):
                filePath = theDb.getFullPath(fileId=fileId)
                mapPath = os.path.splitext(os.path.abspath(filePath))[0]+".map"
                if os.path.isfile(mapPath):
                    return mapPath
                return coefficientsToMap(filePath, mapPath=mapPath, overSample=1.5)

    return None

# Function to return list of names of exportable MTZ(s)
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [ [ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ],
            ['2FoFc_as_map', '2FoFc as map', 'application/CCP4-map'],
            ['FoFc_as_map', 'FoFc as map', 'application/CCP4-map'],
            ]

#============================================================================================
import unittest
class testRefmac(unittest.TestCase):

    def test1(self):
        # Test creation of log file using ../test_data/test1.params.xml input
        from core.CCP4Utils import getCCP4I2Dir
        import os
        workDirectory = CCP4Utils.getTestTmpDir()
        logFile = os.path.join(workDirectory,'prosmart_refmac_test1.log')
        # Delete any existing log file
        if os.path.exists(logFile): os.remove(logFile)
        self.wrapper = prosmart_refmac(name='prosmart_refmac_test1',workDirectory=workDirectory)
        self.wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(),'wrappers','prosmart_refmac','test_data','test1.params.xml'))
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
