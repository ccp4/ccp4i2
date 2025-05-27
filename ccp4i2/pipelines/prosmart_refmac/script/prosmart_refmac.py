"""
Copyright (C) 2010 University of York
"""

import functools
import math
import os
import shutil
import sys
import traceback
import xml.etree.ElementTree as ET

from lxml import etree
from mmtbx.command_line import molprobity
from PySide2 import QtCore
from rdkit import Chem
import clipper

from ....core import CCP4ErrorHandling
from ....core import CCP4Utils
from ....core.CCP4MgImports import mmdb2
from ....core.CCP4Modules import PROJECTSMANAGER
from ....core.CCP4PluginScript import CPluginScript


class prosmart_refmac(CPluginScript):

    TASKMODULE = 'refinement'
    TASKTITLE = 'Refine with Refmac & optional restraints from Prosmart & Platonyzer'
    TASKNAME = 'prosmart_refmac'
    TASKVERSION= 0.0
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 240
    MAXNJOBS = 4
    PERFORMANCECLASS = 'CRefinementPerformance'
    SUBTASKS=['Platonyzer','prosmart','refmac']
    RUNEXTERNALPROCESS=False
    PURGESEARCHLIST =  [[ 'refmac%*/hklout.mtz', 0, "hklout" ], [ 'refmac%*/hklout.mtz', 7, "hklout" ], [ '*%*/ANOMFPHIOUT.mtz', 1, "ANOMFPHIOUT" ], [ '*%*/DIFANOMFPHIOUT.mtz', 1, "DIFANOMFPHIOUT" ]]


    ERROR_CODES = { 101 : { 'description' : 'Error copying data file from final job to pipeline directory' }
                    }

    def __init__(self, *args, **kws):
        super(prosmart_refmac,self).__init__(*args, **kws)
        self.pipelinexmlfile = self.makeFileName(format='PROGRAMXML')
        self.refmacMonitors = {}
        self.xmlroot = ET.Element("RefmacOptimiseWeight")
        self.xmlroot2 = ET.Element("RefmacOptimiseWeight")

#    def startProcess(self, processId):
#        if self.container.prosmartProtein.TOGGLE:
#            self.executeProsmartProtein()
#        else:
#            if self.container.controlParameters.WEIGHT_OPT.__str__()=='MANUAL' and self.container.controlParameters.WEIGHT.isSet(): withWeight = float(self.container.controlParameters.WEIGHT)
#            else: withWeight = -1.
#            self.executeFirstRefmac(withWeight)
#        return CPluginScript.SUCCEEDED

    def startProcess(self, processId):
        try:
           self.executeProsmartProtein()
           self.executeProsmartNucleicAcid()
           self.executePlatonyzer()
           self.executeFirstRefmac()
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
           self.connectSignal(self.prosmart_protein,'finished',self.prosmartProteinFinished)
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

    def executeFirstRefmac(self, withWeight=-1):
        #create wrapper
        self.firstRefmac = self.refmacJobWithWeight(withWeight)
        #Run asynchronously ...this is needed so that commands downsrteam ofprocess launch (i.e. logwather) will be called
        #before process completion
        self.firstRefmac.doAsync = self.doAsync
        self.firstRefmac.connectSignal(self.firstRefmac,'finished',self.firstRefmacFinished)
        #Install xml node for an in progress refmac
        self.xmlLength = 0
        #Start process
        firstRefmacXMLFilename = self.firstRefmac.makeFileName(format='PROGRAMXML')
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("executeFirstRefmac, firstRefmacXMLFilename",firstRefmacXMLFilename)
        self.watchFile(firstRefmacXMLFilename, handler=self.handleXmlChanged, minDeltaSize=34, unwatchWhileHandling=True)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        rv = self.firstRefmac.process()

    @QtCore.Slot(str)
    def handleXmlChanged2(self, xmlFilename):
        self.xmlroot2.clear()
        refmacEtree = ET.parse(xmlFilename).getroot()
        refmacXML = refmacEtree.xpath('//REFMAC')
        if len(refmacXML) == 1:
            refmacXML[0].tag="RefmacPostCootInProgress"
            self.xmlroot2.append(refmacXML[0])
        self.saveXml2()

    @QtCore.Slot(str)
    def handleXmlChanged(self, xmlFilename):
        self.xmlroot.clear()
        refmacEtree = ET.parse(xmlFilename).getroot()
        refmacXML = refmacEtree.xpath('//REFMAC')
        if len(refmacXML) == 1:
            refmacXML[0].tag="RefmacInProgress"
            self.xmlroot.append(refmacXML[0])
        self.saveXml()

    def saveXml2(self):
        ET.indent(newXml)
        newXml = ET.tostring(self.xmlroot2)
        if len(newXml) > self.xmlLength2:
           firstFileName = self.pipelinexmlfile+'_first'
           with open(firstFileName,'r') as aFile:
               oldXml = ET.fromstring(aFile.read())
           oldXml.xpath('//RefmacOptimiseWeight')[0].append(self.xmlroot2.xpath("//RefmacOptimiseWeight/RefmacPostCootInProgress")[0])
           tmpFileName = self.pipelinexmlfile+'_tmp'
           with open(tmpFileName,'w') as aFile:
               ET.indent(oldXml)
               CCP4Utils.writeXML(aFile,ET.tostring(oldXml))
           shutil.move(tmpFileName, self.pipelinexmlfile)
           self.xmlLength2 = len(newXml)

    def saveXml(self):
        # Save the xml if it has grown
        ET.indent(self.xmlroot)
        newXml = ET.tostring(self.xmlroot)
        if len(newXml) > self.xmlLength:
           tmpFileName = self.pipelinexmlfile+'_tmp'
           with open(tmpFileName,'w') as aFile:
               CCP4Utils.writeXML(aFile, newXml )
           shutil.move(tmpFileName, self.pipelinexmlfile)
           self.xmlLength = len(newXml)

    def refmacJobWithWeight(self, withWeight=-1, inputCoordinates=None, ncyc=-1):
        result = self.makePluginObject('refmac')
        #input data for this refmac instance is the same as the input data for the program
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

        #Set parameters if there is prior ProSMART or Platonyzer
        if self.container.prosmartProtein.TOGGLE:
            result.container.controlParameters.PROSMART_PROTEIN_WEIGHT=self.container.prosmartProtein.WEIGHT
            result.container.controlParameters.PROSMART_PROTEIN_ALPHA=self.container.prosmartProtein.ALPHA
            result.container.controlParameters.PROSMART_PROTEIN_DMAX=self.container.prosmartProtein.DMAX
            result.container.inputData.PROSMART_PROTEIN_RESTRAINTS=self.prosmart_protein.container.outputData.RESTRAINTS
        if self.container.prosmartNucleicAcid.TOGGLE:
            result.container.controlParameters.PROSMART_NUCLEICACID_WEIGHT=self.container.prosmartNucleicAcid.WEIGHT
            result.container.controlParameters.PROSMART_NUCLEICACID_ALPHA=self.container.prosmartNucleicAcid.ALPHA
            result.container.controlParameters.PROSMART_NUCLEICACID_DMAX=self.container.prosmartNucleicAcid.DMAX
            result.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS=self.prosmart_nucleicacid.container.outputData.RESTRAINTS
        if self.container.platonyzer.TOGGLE:
            result.container.inputData.PLATONYZER_RESTRAINTS=self.platonyzer.container.outputData.RESTRAINTS
            result.container.inputData.XYZIN=self.platonyzer.container.outputData.XYZOUT

        #Specify weight if a meaningful one has been offered
        if withWeight>=0.:
            result.container.controlParameters.WEIGHT_OPT='MANUAL'
            result.container.controlParameters.WEIGHT = withWeight

        if inputCoordinates is not None:
            result.container.inputData.XYZIN.set(inputCoordinates)
            if ncyc > 0:
                result.container.controlParameters.NCYCLES.set(ncyc)
        return result

    @QtCore.Slot(dict)
    def firstRefmacFinished(self, statusDict):
        print("AAA1")
        if statusDict['finishStatus'] == CPluginScript.UNSATISFACTORY:
            print("AAA1.UNSATISFACTORY")
            if os.path.isfile(self.firstRefmac.container.outputData.LIBOUT.__str__()):
                from ....wrappers.acedrg.script import acedrg
                try:
                    rdkitMol = acedrg.molFromDict(self.firstRefmac.container.outputData.LIBOUT.__str__())
                    molRemovedHs = Chem.RemoveHs(rdkitMol)
                    svgXml = acedrg.svgFromMol(molRemovedHs)
                    self.xmlroot.append(svgXml)
                except:
                    print('Unable to generate svg from DICT')
                shutil.copyfile(self.firstRefmac.container.outputData.LIBOUT.__str__(), self.container.outputData.LIBOUT.__str__())
                self.container.outputData.LIBOUT.annotation = 'Refmac-generated library...use with caution'
            if os.path.isfile(self.firstRefmac.container.outputData.PSOUT.__str__()):
                shutil.copyfile(self.firstRefmac.container.outputData.PSOUT.__str__(), self.container.outputData.PSOUT.__str__())
                self.container.outputData.PSOUT.annotation.set('Pictures of ligand prepared by refmac')
            with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
                ET.indent(self.xmlroot)
                CCP4Utils.writeXML(programXML,ET.tostring(self.xmlroot))
            self.reportStatus(CPluginScript.UNSATISFACTORY)

        elif self.firstRefmac.errorReport.maxSeverity() > CCP4ErrorHandling.Severity.WARNING:
            print("AAA1.MAXSEVERITY")
            #This gets done in thefirstRefmac.reportStatus() - Liz
            #self.extendErrorReport(self.firstRefmac.errorReport)
            try:
              refmacEtree = ET.parse(self.firstRefmac.makeFileName('PROGRAMXML')).getroot()
              refmacXML = refmacEtree.xpath('//REFMAC')
              if len(refmacXML) == 1: self.xmlroot.append(refmacXML[0])
            except:
              print('Failed attempt to read XML file from first Refmac')
            try:
              ET.indent(self.xmlroot)
              newXml = ET.tostring(self.xmlroot)
              aFile = open(self.pipelinexmlfile,'w')
              CCP4Utils.writeXML(aFile,newXml)
              aFile.close()
            except:
               print('Failed attempt to write pipeline XML file')
            self.reportStatus(CPluginScript.FAILED)
            return
        print("AAA11")

        self.handleXmlChanged(self.firstRefmac.makeFileName(format='PROGRAMXML'))

        print("AAA10")
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            #This gets done in the firstRefmac.reportStatus() - Liz
            self.fileSystemWatcher = None
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            print("AAA12")
            self.addCycleXML(self.firstRefmac)
            aFile=open( self.pipelinexmlfile,'w')
            ET.indent(self.xmlroot)
            CCP4Utils.writeXML(aFile, ET.tostring(self.xmlroot) )
            aFile.close()
            print("AAA13")
            if self.container.controlParameters.OPTIMISE_WEIGHT:
                print("AAA14")
                self.fileSystemWatcher = None
                weightUsed = float(self.xmlroot.xpath('//weight')[-1].text)
                self.tryVariousRefmacWeightsAround(weightUsed)
            else:
               print("AAA15")
               best_r_free = self.firstRefmac.container.outputData.PERFORMANCEINDICATOR.RFactor
               print("AAA15.1")
               if self.container.controlParameters.ADD_WATERS and best_r_free < self.container.controlParameters.REFPRO_RSR_RWORK_LIMIT :
                   print("AAA16")
                   self.currentCoordinates = self.firstRefmac.container.outputData.CIFFILE
                   self.cootPlugin = self.makeCootPlugin()
                   self.cootPlugin.doAsync = self.doAsync
                   self.cootPlugin.connectSignal(self.cootPlugin,'finished',self.cootFinished)
                   print("AAA17")
                   rv = self.cootPlugin.process()
                   if rv == CPluginScript.FAILED:
                        self.reportStatus(rv)
               else:
                     self.fileSystemWatcher = None
                     self.finishUp(self.firstRefmac)
        print('done prosmart_refmac.firstRefmacFinished')
        return

    def makeCootPlugin(self):
         # FIXME - This is all nonsense - needs to consider best task, etc... *NOT* just firstRefmaca?
        cootPlugin = self.makePluginObject('coot_find_waters')
        cootPlugin.container.inputData.XYZIN.set(self.currentCoordinates)
        cootPlugin.container.inputData.FPHIIN.set(self.firstRefmac.container.outputData.FPHIOUT)
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
        # Check coot status and start refmac
        self.checkFinishStatus(statusDict=statusDict,failedErrCode=204, outputFile= self.cootPlugin.container.outputData.XYZOUT,noFileErrCode=205)
        try:
          if self.container.controlParameters.ADD_WATERS:
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = ET.fromstring(aFile.read())
            aFile.close()
            nwaters = "unknown"
            cootLogTxt = os.path.join(os.path.dirname(self.cootPlugin.container.outputData.XYZOUT.__str__()),"log.txt")
            with open(cootLogTxt, 'r') as f:
               for l in f:
                   if l.startswith("INFO::") and "found" in l and "water fitting" in l:
                      nwaters = l.strip()
                      numsearch = [ x for x in nwaters.split() if x.isdigit() ]
                      if len(numsearch)>0:
                         nwaters = numsearch[0]
                      break
            postRefmacCoot = ET.Element("CootAddWaters")
            postRefmacCoot.text = "Coot added "+nwaters+" waters"
            oldXml.append(postRefmacCoot)
            aFile = open(self.pipelinexmlfile+'_tmpcoot','w')
            ET.indent(oldXml)
            CCP4Utils.writeXML(aFile,ET.tostring(oldXml))
            aFile.close()
            shutil.move(self.pipelinexmlfile+'_tmpcoot', self.pipelinexmlfile)
          self.cootPlugin.container.outputData.XYZOUT.subType = 1
          self.currentCoordinates = self.cootPlugin.container.outputData.XYZOUT
          self.refmacPostCootPlugin = self.refmacJobWithWeight(inputCoordinates=self.currentCoordinates,ncyc=5)
          self.refmacPostCootPlugin.doAsync = True
          self.refmacPostCootPlugin.connectSignal(self.refmacPostCootPlugin,'finished',functools.partial(self.postCootRefmacFinished,self.refmacPostCootPlugin))
          refmacPostCootXMLFilename = self.refmacPostCootPlugin.makeFileName(format='PROGRAMXML')
          self.xmlLength2 = 0
          self.watchFile(refmacPostCootXMLFilename, handler=self.handleXmlChanged2, minDeltaSize=34, unwatchWhileHandling=True)
          shutil.copyfile(self.pipelinexmlfile, self.pipelinexmlfile+"_first")
          rv = self.refmacPostCootPlugin.process()
          if rv == CPluginScript.FAILED: self.reportStatus(rv)
        except Exception as e:
          self.appendErrorReport(CPluginScript,39,str(e))

    @QtCore.Slot('CPluginScript',dict)
    def postCootRefmacFinished(self, refmacJob, statusDict={}):
        self.handleXmlChanged(refmacJob.makeFileName(format='PROGRAMXML'))

        self.fileSystemWatcher = None
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            #This gets done in the firstRefmac.reportStatus() - Liz
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            self.addCycleXML(refmacJob,"refmacPostCoot")
            ET.indent(self.xmlroot)
            newXml = ET.tostring(self.xmlroot)
            refmacPostCoot = self.xmlroot.xpath("refmacPostCoot")
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = ET.fromstring(aFile.read())
            aFile.close()
            oldXml.extend(refmacPostCoot)
            aFile = open(self.pipelinexmlfile,'w')
            ET.indent(oldXml)
            CCP4Utils.writeXML(aFile,ET.tostring(oldXml))
            aFile.close()

        self.finishUp(refmacJob)

    def finishUp(self, refmacJob):
        from ....core import CCP4ProjectsManager
        print('into prosmart_refmac.finishUp')
        for attr in self.container.outputData.dataOrder():
            print('prosmart_refmac.finishUp attr',attr)
            wrappersAttr = getattr(refmacJob.container.outputData, attr)
            pipelinesAttr = getattr(self.container.outputData, attr)
            if attr in [ "PERFORMANCEINDICATOR"]:
                setattr(self.container.outputData, attr, wrappersAttr)
            else:
                if os.path.exists(str(wrappersAttr.fullPath)):
                  try:
                    shutil.copyfile( str(wrappersAttr.fullPath), str(pipelinesAttr.fullPath) )
                  except:
                    self.appendErrorReport(101,str(wrappersAttr.fullPath)+' to '+str(pipelinesAttr.fullPath))
                  #self.container.outputData.copyData(refmacJob.container.outputData,[attr])
                if attr == "XMLOUT":
                    pass

        print('prosmart_refmac.finishUp 1')
        from ....core import CCP4XtalData
        # Apply database annotations
        self.container.outputData.XYZOUT.annotation.set('Model from refinement (PDB format)')
        self.container.outputData.CIFFILE.annotation.set('Model from refinement (mmCIF format)')
        self.container.outputData.CIFFILEDEP.annotation.set('Model from refinement (mmCIF format for deposition)')
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from refinement'
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from refinement'
        self.container.outputData.DIFFPHIOUT.subType = 2
        self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        self.container.outputData.TLSOUT.annotation = 'TLS parameters from refinement'
        self.container.outputData.COOTSCRIPTOUT.annotation = 'Coot script written from refinement'
        self.container.outputData.ANOMFPHIOUT.annotation = 'Weighted anomalous difference map from refinement'
        self.container.outputData.DIFANOMFPHIOUT.annotation = 'Weighted differences of anomalous difference map'

        print('prosmart_refmac.finishUp 2')
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
              ligxml = ET.SubElement(self.xmlroot,"LIGANDS")
              for item in self.container.outputData.LIBOUT.fileContent.monomerList:
                ligNode = ET.SubElement(ligxml,"ligand")
                ligNode.text = str(item.three_letter_code)
              self.saveXml()
          except:
              print('Error creating LIBOUT annotation')
              self.container.outputData.LIBOUT.annotation = annotation
          try:
              self.mergeDictToProjectLib(fileName=self.container.outputData.LIBOUT.__str__())
          except:
              print('Error merging library to Project Dictionary')
        print('prosmart_refmac.finishUp 3'); sys.stdout.flush()
#
        cleanUpIntermediate = False
        if hasattr(self.container.controlParameters,"REFMAC_CLEANUP"):
            cleanUpIntermediate = self.container.controlParameters.REFMAC_CLEANUP
            if cleanUpIntermediate:
                print('prosmart_refmac.finishUp 4'); sys.stdout.flush()
                cleanup = CCP4ProjectsManager.CPurgeProject(self.firstRefmac._dbProjectId)
                print('prosmart_refmac.finishUp 5'); sys.stdout.flush()
                cleanup.purgeJob(self.firstRefmac.jobId,context="extended_intermediate",reportMode="skip")

                if hasattr(self,"refmacPostCootPlugin"):
                    print('prosmart_refmac.finishUp 6'); sys.stdout.flush()
                    cleanup = CCP4ProjectsManager.CPurgeProject(self.refmacPostCootPlugin._dbProjectId)
                    print('prosmart_refmac.finishUp 7'); sys.stdout.flush()
                    cleanup.purgeJob(self.refmacPostCootPlugin.jobId,context="extended_intermediate",reportMode="skip")

        if not str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID':
           #Geometry validation not performed in rigid body mode

           runMolprobity = False
           if hasattr(self.container.controlParameters,"RUN_MOLPROBITY"):
              runMolprobity = self.container.controlParameters.RUN_MOLPROBITY
              if runMolprobity:
                 #MN: First order attempt at providing molprobity analysis
                 try:
                     print("Attempting molprobity run after refinement...")
                     coordPath = self.container.outputData.XYZOUT.fullPath.__str__()
                     fileRoot, fileExt = os.path.splitext(coordPath)
                     sanitizedCoordPath = fileRoot + "+asPDB.pdb"

                     #Use mmdb to do some sanitization
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
                         ET.SubElement(ET.SubElement(self.xmlroot,"Molprobity"), "Output").text = etree.CDATA(open(mpOutPath).read())
                     self.saveXml()
                     print("...Succeeded molprobity run after refinement :-)")
                 except Exception as err:
                     ET.SubElement(ET.SubElement(self.xmlroot,"Molprobity"), "Output").text = etree.CDATA(str(err))
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

           if validate_baverage or validate_molprobity or validate_ramachandran or validate_iris:
             xml_validation = ET.SubElement(self.xmlroot,"Validation")
             xml_validation_status = ET.SubElement(xml_validation,"Success")
             try:
                self.validate = self.makePluginObject('validate_protein')
                self.validate.container.inputData.XYZIN_1.set(self.container.outputData.XYZOUT)
                self.validate.container.inputData.XYZIN_2.set(self.container.inputData.XYZIN)
                if str(self.container.controlParameters.SCATTERING_FACTORS) == "XRAY":
                    self.validate.container.inputData.F_SIGF_1.set(self.container.inputData.F_SIGF)
                    self.validate.container.inputData.F_SIGF_2.set(self.container.inputData.F_SIGF)
                else:
                    self.validate.container.inputData.F_SIGF_1.set(None)
                    self.validate.container.inputData.F_SIGF_2.set(None)
                self.validate.container.inputData.NAME_1 = "Refined"
                self.validate.container.inputData.NAME_2 = "Input"
                #MN...Using "="" to set this is an odd thing and breaks under some circumstances.
                #Specifically, the path for this is now in the prosmart_refmac (i.e. parent job) directory,
                #but it is an output of the validate_protein subjob.  When that subjob seeks to save it to the
                #database, it finds that it is in the wrong directory and croaks.  I have tried to get round this
                #by flagging it as saveToDb=False in the validate_protein .def.xml. But htis needs consideration by jon
                #Agirre and Rob Nicholls
                self.validate.container.outputData.COOTSCRIPTOUT = self.container.outputData.COOTSCRIPTOUT

                self.validate.container.controlParameters.DO_IRIS = validate_iris
                self.validate.container.controlParameters.DO_BFACT = validate_baverage
                self.validate.container.controlParameters.DO_RAMA = validate_ramachandran
                self.validate.container.controlParameters.DO_MOLPROBITY = validate_molprobity

                self.validate.doAsync = False
                self.validate.waitForFinished = -1
                self.validate.process()

                validateXMLPath = self.validate.makeFileName('PROGRAMXML')
                validateXML = ET.parse(validateXMLPath).getroot()
                if len(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info"))>0:
                   xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info")[0]) 
                if self.validate.container.controlParameters.DO_IRIS:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/Iris"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Iris")[0]) 
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
                       if hasattr(self,"refmacPostCootPlugin"):
                           refmaclog = self.refmacPostCootPlugin.makeFileName('LOG') #"/Users/stuart/CCP4I2_PROJECTS/5_7_2021/CCP4_JOBS/job_19/job_1/log.txt"
                       else:
                           refmaclog = self.firstRefmac.makeFileName('LOG')

                       verdict_result = prosmart_refmac_verdict.getJSCOFERefmac5Verdict(programxml=programxml,pdbfile=pdbfile,refmaclog=refmaclog)
                       verdict_score = verdict_result["score"]
                       verdict_message  = verdict_result["message"]
                       bottomline = verdict_result["bottomLine"]
                       meanRfree = verdict_result["meanRfree"]
                       medianClash = verdict_result["medianClash"]
                       ramaOutliers = verdict_result["ramaOutliers"]
                       suggestedParameters = verdict_result["suggestedParameters"]

                       suggestedParameters = self.mapVerdictSuggestionsToi2Params(suggestedParameters)

                       xml_verdict = ET.SubElement(self.xmlroot,"Verdict")
                       xml_verdict_score = ET.SubElement(xml_verdict,"verdict_score")
                       xml_verdict_score.text = str(verdict_score)
                       xml_verdict_message = ET.SubElement(xml_verdict,"verdict_message")
                       xml_verdict_message.text = etree.CDATA(verdict_message)
                       xml_bottomline = ET.SubElement(xml_verdict,"bottomline")
                       xml_bottomline.text = etree.CDATA(bottomline)
                       xml_meanRfree = ET.SubElement(xml_verdict,"meanRfree")
                       xml_meanRfree.text = str(meanRfree)
                       xml_medianClash = ET.SubElement(xml_verdict,"medianClash")
                       xml_medianClash.text = str(medianClash)
                       xml_ramaOutliers = ET.SubElement(xml_verdict,"ramaOutliers")
                       xml_ramaOutliers.text = str(ramaOutliers)
                       xml_suggestedParameters = ET.SubElement(xml_verdict,"suggestedParameters")
                       for k,v in suggestedParameters.items():
                          xml_suggestedParameters_k = ET.SubElement(xml_suggestedParameters,k)
                          xml_suggestedParameters_k.text = str(v)

                       self.saveXml()
                   except:
                       print("Some problem with verdict...."); sys.stdout.flush()
                       exc_type, exc_value, exc_tb = sys.exc_info()[:3]
                       sys.stderr.write(str(exc_type) + '\n')
                       sys.stderr.write(str(exc_value) + '\n')
                       traceback.print_tb(exc_tb)
                xml_validation_status.text = "SUCCESS"
                self.saveXml()
             except Exception as err:
                xml_validation_status.text = "FAILURE"
                self.saveXml()
                traceback.print_exc()
                print("...Failed validation run after refinement", err)

        logfiles = []
        if hasattr(self,"firstRefmac"):
            logfiles.append(self.firstRefmac.makeFileName('LOG'))
        if hasattr(self,"refmacPostCootPlugin"):
            logfiles.append(self.refmacPostCootPlugin.makeFileName('LOG'))

        self.createWarningsXML(logfiles)
        self.saveXml()

        print('done prosmart_refmac.finishUp'); sys.stdout.flush()
        self.reportStatus(CPluginScript.SUCCEEDED)

    def tryVariousRefmacWeightsAround(self, weight):
        print('Generating jobs with weights around ', weight)
        #make an array to hold the child-jobs
        refmacJobs = []
        for factorExponent in range(-3,4):
            if factorExponent != 0:
                factor = math.pow(2., factorExponent)
                refmacJobs.append(self.refmacJobWithWeight(float(factor)*float(weight)))
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
        # callback is passed the jobId (=Non
        # if not in ccp4i2-db context) and processId that
        # can serve at identifier for subProcess

        #print 'demo_multi_mtzdump.handleDone',ret

        rtask = self.jobsInTrain[str(pid)]

        self.addCycleXML(rtask)
        # Save the xml on every cycle
        aFile=open( self.pipelinexmlfile,'w')
        ET.indent(self.xmlroot)
        CCP4Utils.writeXML(aFile, ET.tostring(self.xmlroot) )
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

    def addCycleXML(self, rtask, name="RefmacWeight"):
        xmlcyc = ET.SubElement(self.xmlroot,name)
        cycWeightNode = ET.SubElement(xmlcyc,"weight")
        rxml = ET.parse(rtask.makeFileName('PROGRAMXML')).getroot()
        try: cycWeightNode.text = rxml.xpath("//Cycle/WeightUsed")[-1].text
        except: pass
        rstats = rxml.xpath("//REFMAC")
        xmlcyc.append (rstats[0])

    def handleTimeout(self):
        sys.stdout.flush()

        for rtask in self.jobsInTrain:
            print('TERMINATING', rtask.processId,sys.stdout.flush())
            try:
                rtask.terminate()
            except:
                pass

        self.appendErrorReport(40,str(self.TIMEOUT_PERIOD))
        self.reportStatus(CPluginScript.FAILED)

def coefficientsToMap(coefficientsPath, mapPath=None, overSample=1.0):
    mtz_file = clipper.CCP4MTZfile()
    hkl_info = clipper.HKL_info()
    mtz_file.open_read (str(coefficientsPath))
    mtz_file.import_hkl_info ( hkl_info )
    sg, cell = hkl_info.spacegroup(), hkl_info.cell()
    fphidata = clipper.HKL_data_F_phi_float(hkl_info)
    mtz_file.import_hkl_data( fphidata, str("/*/*/[F,PHI]") );
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
    from ....dbapi.CCP4DbApi import FILE_ROLE_OUT

    theDb = PROJECTSMANAGER().db()
    if mode == 'complete_mtz':
        #print 'refmac.exportJobFile',mode
        childJobs = theDb.getChildJobs(jobId=jobId,details=True)
        #print 'exportJobFile childJobs',childJobs
        if childJobs[-1][2] == 'refmac':
          jobDir = PROJECTSMANAGER().jobDirectory(jobId=childJobs[-1][1],create=False)
          if os.path.exists(os.path.join(jobDir,'hklout.mtz')):
             return  os.path.join(jobDir,'hklout.mtz')
        elif childJobs[-1][2] == 'validate_protein':
          if childJobs[-2][2] == 'refmac':
             jobDir = PROJECTSMANAGER().jobDirectory(jobId=childJobs[-2][1],create=False)
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
