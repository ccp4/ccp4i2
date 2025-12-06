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

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from lxml import etree
from baselayer import QtCore, DJANGO
from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core import CCP4Utils
import os,sys,shutil,re
import traceback
from wrappers.modelASUCheck.script.modelASUCheck import sequenceAlignment


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


    ERROR_CODES = {
        101: {'description': 'Error copying data file from final job to pipeline directory'},
        102: {'description': 'ProSMART protein restraints failed'},
        103: {'description': 'ProSMART nucleic acid restraints failed'},
        104: {'description': 'Platonyzer restraints failed'},
        105: {'description': 'Refmac refinement failed'},
        106: {'description': 'Refmac output file not created'},
        107: {'description': 'Coot find waters failed'},
        108: {'description': 'Post-coot refmac failed'},
        109: {'description': 'Weight optimization failed'},
        110: {'description': 'Validation failed'},
    }

    def __init__(self, *args, **kws):
        super(prosmart_refmac,self).__init__(*args, **kws)
        self.pipelinexmlfile = self.makeFileName(format='PROGRAMXML')
        self.refmacMonitors = {}
        self.xmlroot = etree.Element("RefmacOptimiseWeight")
        self.xmlroot2 = etree.Element("RefmacOptimiseWeight")

    def startProcess(self, processId):
        """
        Main pipeline execution - runs ProSMART, Platonyzer, then Refmac.

        Each phase uses try/except with proper CErrorReport logging including traceback.
        """
        # Phase 1: ProSMART protein restraints (optional)
        try:
            self.executeProsmartProtein()
        except Exception as e:
            tb = traceback.format_exc()
            self.appendErrorReport(102, f'ProSMART protein restraints: {e}\n{tb}')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Phase 2: ProSMART nucleic acid restraints (optional)
        try:
            self.executeProsmartNucleicAcid()
        except Exception as e:
            tb = traceback.format_exc()
            self.appendErrorReport(103, f'ProSMART nucleic acid restraints: {e}\n{tb}')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Phase 3: Platonyzer restraints (optional)
        try:
            self.executePlatonyzer()
        except Exception as e:
            tb = traceback.format_exc()
            self.appendErrorReport(104, f'Platonyzer restraints: {e}\n{tb}')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Phase 4: Refmac refinement
        try:
            self.executeFirstRefmac()
        except Exception as e:
            tb = traceback.format_exc()
            self.appendErrorReport(105, f'Refmac refinement: {e}\n{tb}')
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
        """Execute the main Refmac refinement step."""
        # Create wrapper
        self.firstRefmac = self.refmacJobWithWeight(withWeight)
        # Run asynchronously - needed so logwatcher callbacks work before process completion
        self.firstRefmac.doAsync = self.doAsync
        self.firstRefmac.connectSignal(self.firstRefmac, 'finished', self.firstRefmacFinished)
        # Install xml node for in-progress refmac
        self.xmlLength = 0

        firstRefmacXMLFilename = self.firstRefmac.makeFileName(format='PROGRAMXML')

        if DJANGO():
            # Django mode: Use Signal/Slot for streaming updates (atomic writes, no race conditions)
            self.firstRefmac.connectSignal(self.firstRefmac, 'progressUpdated', self.handleRefmacProgress)
        else:
            # Qt mode: Use file watching as QProcess doesn't support our streaming mechanism
            self.watchFile(firstRefmacXMLFilename, handler=self.handleXmlChanged, minDeltaSize=34, unwatchWhileHandling=True)

        rv = self.firstRefmac.process()

    def handleRefmacProgress(self, progressInfo):
        """Handle progress signal from refmac (streaming mode).

        This is called when refmac emits progressUpdated signal, avoiding
        file-watching race conditions. We read the XML file that refmac
        just wrote atomically.
        """
        try:
            firstRefmacXMLFilename = self.firstRefmac.makeFileName(format='PROGRAMXML')
            self.handleXmlChanged(firstRefmacXMLFilename)
        except Exception as e:
            print(f"Warning: Error handling refmac progress: {e}")

    @QtCore.Slot(str)
    def handleXmlChanged2(self, xmlFilename):
        self.xmlroot2.clear()
        refmacEtree = CCP4Utils.openFileToEtree(xmlFilename)
        refmacXML = refmacEtree.xpath('//REFMAC')
        if len(refmacXML) == 1:
            refmacXML[0].tag="RefmacPostCootInProgress"
            self.xmlroot2.append(refmacXML[0])
        self.saveXml2()

    @QtCore.Slot(str)
    def handleXmlChanged(self, xmlFilename):
        self.xmlroot.clear()
        refmacEtree = CCP4Utils.openFileToEtree(xmlFilename)
        refmacXML = refmacEtree.xpath('//REFMAC')
        if len(refmacXML) == 1:
            refmacXML[0].tag="RefmacInProgress"
            self.xmlroot.append(refmacXML[0])
        self.saveXml()

    def saveXml2(self):
        newXml = etree.tostring(self.xmlroot2,pretty_print=True)
        if len(newXml) > self.xmlLength2:
           firstFileName = self.pipelinexmlfile+'_first'
           with open(firstFileName,'r') as aFile:
               oldXml = etree.fromstring(aFile.read())
           oldXml.xpath('//RefmacOptimiseWeight')[0].append(self.xmlroot2.xpath("//RefmacOptimiseWeight/RefmacPostCootInProgress")[0])
           tmpFileName = self.pipelinexmlfile+'_tmp'
           with open(tmpFileName,'w') as aFile:
               CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True) )
           shutil.move(tmpFileName, self.pipelinexmlfile)
           self.xmlLength2 = len(newXml)

    def saveXml(self):
        # Save the xml if it has grown
        newXml = etree.tostring(self.xmlroot,pretty_print=True)
        if len(newXml) > self.xmlLength:
           tmpFileName = self.pipelinexmlfile+'_tmp'
           with open(tmpFileName,'w') as aFile:
               CCP4Utils.writeXML(aFile, newXml )
           import shutil
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
        """Handle completion of the first Refmac refinement job."""
        # Handle unsatisfactory result (e.g., ligand geometry issues)
        if statusDict['finishStatus'] == CPluginScript.UNSATISFACTORY:
            import os
            if os.path.isfile(self.firstRefmac.container.outputData.LIBOUT.__str__()):
                from wrappers.acedrg.script import acedrg
                try:
                    rdkitMol = acedrg.molFromDict(self.firstRefmac.container.outputData.LIBOUT.__str__())
                    from rdkit import Chem
                    molRemovedHs = Chem.RemoveHs(rdkitMol)
                    svgXml = acedrg.svgFromMol(molRemovedHs)
                    self.xmlroot.append(svgXml)
                except Exception as e:
                    self.appendErrorReport(106, f'Unable to generate SVG from ligand dictionary: {e}')
                import shutil
                shutil.copyfile(self.firstRefmac.container.outputData.LIBOUT.__str__(), self.container.outputData.LIBOUT.__str__())
                self.container.outputData.LIBOUT.annotation = 'Refmac-generated library...use with caution'
            if os.path.isfile(self.firstRefmac.container.outputData.PSOUT.__str__()):
                shutil.copyfile(self.firstRefmac.container.outputData.PSOUT.__str__(), self.container.outputData.PSOUT.__str__())
                self.container.outputData.PSOUT.annotation.set('Pictures of ligand prepared by refmac')
            with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
                CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))
            self.reportStatus(CPluginScript.UNSATISFACTORY)
            return

        # Handle errors with severity > warning
        if self.firstRefmac.errorReport.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            try:
                refmacEtree = CCP4Utils.openFileToEtree(self.firstRefmac.makeFileName('PROGRAMXML'))
                refmacXML = refmacEtree.xpath('//REFMAC')
                if len(refmacXML) == 1:
                    self.xmlroot.append(refmacXML[0])
            except Exception as e:
                self.appendErrorReport(105, f'Failed to read Refmac XML output: {e}')
            try:
                newXml = etree.tostring(self.xmlroot,pretty_print=True)
                with open(self.pipelinexmlfile,'w') as aFile:
                    CCP4Utils.writeXML(aFile,newXml)
            except Exception as e:
                self.appendErrorReport(105, f'Failed to write pipeline XML: {e}')
            self.reportStatus(CPluginScript.FAILED)
            return

        # Update XML with refmac results
        self.handleXmlChanged(self.firstRefmac.makeFileName(format='PROGRAMXML'))

        # Handle explicit failure status
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            self.fileSystemWatcher = None
            self.appendErrorReport(105, 'Refmac reported failure status')
            self.reportStatus(CPluginScript.FAILED)
            return

        # Success path - add cycle XML and continue
        self.addCycleXML(self.firstRefmac)
        with open(self.pipelinexmlfile,'w') as aFile:
            CCP4Utils.writeXML(aFile, etree.tostring(self.xmlroot,pretty_print=True))

        # Check if weight optimization is requested
        if self.container.controlParameters.OPTIMISE_WEIGHT:
            self.fileSystemWatcher = None
            try:
                weightUsed = float(self.xmlroot.xpath('//weight')[-1].text)
                self.tryVariousRefmacWeightsAround(weightUsed)
            except Exception as e:
                self.appendErrorReport(109, f'Weight optimization: {e}')
                self.reportStatus(CPluginScript.FAILED)
        else:
            # Check if water finding should be done
            best_r_free = self.firstRefmac.container.outputData.PERFORMANCEINDICATOR.RFactor
            if self.container.controlParameters.ADD_WATERS and best_r_free < self.container.controlParameters.REFPRO_RSR_RWORK_LIMIT:
                try:
                    self.currentCoordinates = self.firstRefmac.container.outputData.CIFFILE
                    self.cootPlugin = self.makeCootPlugin()
                    self.cootPlugin.doAsync = self.doAsync
                    self.cootPlugin.connectSignal(self.cootPlugin,'finished',self.cootFinished)
                    rv = self.cootPlugin.process()
                    if rv == CPluginScript.FAILED:
                        self.appendErrorReport(107, 'Coot find waters process() returned FAILED')
                        self.reportStatus(rv)
                except Exception as e:
                    self.appendErrorReport(107, f'Coot find waters: {e}')
                    self.reportStatus(CPluginScript.FAILED)
            else:
                self.fileSystemWatcher = None
                self.finishUp(self.firstRefmac)
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
        import functools
        # Check coot status and start refmac
        self.checkFinishStatus(statusDict=statusDict,failedErrCode=204, outputFile= self.cootPlugin.container.outputData.XYZOUT,noFileErrCode=205)
        try:
          if self.container.controlParameters.ADD_WATERS:
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = etree.fromstring(aFile.read())
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
            postRefmacCoot = etree.Element("CootAddWaters")
            postRefmacCoot.text = "Coot added "+nwaters+" waters"
            oldXml.append(postRefmacCoot)
            aFile = open(self.pipelinexmlfile+'_tmpcoot','w')
            CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True))
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

        import os
        self.fileSystemWatcher = None
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            #This gets done in the firstRefmac.reportStatus() - Liz
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            self.addCycleXML(refmacJob,"refmacPostCoot")
            newXml = etree.tostring(self.xmlroot,pretty_print=True)
            refmacPostCoot = self.xmlroot.xpath("refmacPostCoot")
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = etree.fromstring(aFile.read())
            aFile.close()
            oldXml.extend(refmacPostCoot)
            aFile = open(self.pipelinexmlfile,'w')
            CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True))
            aFile.close()

        self.finishUp(refmacJob)

    def finishUp(self, refmacJob):
        """
        Copy output files from refmac job to pipeline output and finalize.

        Handles validation, cleanup, and optional molprobity analysis.
        """
        from core import CCP4ProjectsManager
        import shutil

        # Copy output files from refmac job to pipeline directory
        for attr in self.container.outputData.dataOrder():
            wrappersAttr = getattr(refmacJob.container.outputData, attr)
            pipelinesAttr = getattr(self.container.outputData, attr)
            if attr in ["PERFORMANCEINDICATOR"]:
                setattr(self.container.outputData, attr, wrappersAttr)
            else:
                if os.path.exists(str(wrappersAttr.fullPath)):
                    try:
                        shutil.copyfile(str(wrappersAttr.fullPath), str(pipelinesAttr.fullPath))
                    except Exception as e:
                        self.appendErrorReport(101, f'{wrappersAttr.fullPath} to {pipelinesAttr.fullPath}: {e}')

        from core import CCP4XtalData
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

        # Set ligand dictionary annotations
        if self.container.outputData.DICT.exists():
            self.container.outputData.DICT.annotation = 'Accumulated ligand geometry dictionary'
        if self.container.outputData.LIBOUT.exists():
            annotation = 'Refmac generated geometry'
            try:
                if len(self.container.outputData.LIBOUT.fileContent.monomerList) > 0:
                    annotation = 'Refmac generated geometry for:'
                    for item in self.container.outputData.LIBOUT.fileContent.monomerList:
                        annotation = annotation + ' ' + str(item.three_letter_code)
                    ligxml = etree.SubElement(self.xmlroot, "LIGANDS")
                    for item in self.container.outputData.LIBOUT.fileContent.monomerList:
                        ligNode = etree.SubElement(ligxml, "ligand")
                        ligNode.text = str(item.three_letter_code)
                    self.saveXml()
            except Exception as e:
                # Non-fatal - just use default annotation
                self.container.outputData.LIBOUT.annotation = annotation
            try:
                self.mergeDictToProjectLib(fileName=self.container.outputData.LIBOUT.__str__())
            except Exception as e:
                # Non-fatal - dictionary merge is optional
                pass

        # Optional cleanup of intermediate files
        cleanUpIntermediate = False
        if hasattr(self.container.controlParameters, "REFMAC_CLEANUP"):
            cleanUpIntermediate = self.container.controlParameters.REFMAC_CLEANUP
            if cleanUpIntermediate:
                try:
                    cleanup = CCP4ProjectsManager.CPurgeProject(self.firstRefmac._dbProjectId)
                    cleanup.purgeJob(self.firstRefmac.jobId, context="extended_intermediate", reportMode="skip")

                    if hasattr(self, "refmacPostCootPlugin"):
                        cleanup = CCP4ProjectsManager.CPurgeProject(self.refmacPostCootPlugin._dbProjectId)
                        cleanup.purgeJob(self.refmacPostCootPlugin.jobId, context="extended_intermediate", reportMode="skip")
                except Exception as e:
                    # Non-fatal - cleanup failure shouldn't fail the job
                    pass

        if not str(self.container.controlParameters.REFINEMENT_MODE) == 'RIGID':
           #Geometry validation not performed in rigid body mode

           # Optional molprobity analysis
           runMolprobity = False
           if hasattr(self.container.controlParameters, "RUN_MOLPROBITY"):
               runMolprobity = self.container.controlParameters.RUN_MOLPROBITY
               if runMolprobity:
                   try:
                       from mmtbx.command_line import molprobity
                       coordPath = self.container.outputData.XYZOUT.fullPath.__str__()
                       fileRoot, fileExt = os.path.splitext(coordPath)
                       sanitizedCoordPath = fileRoot + "+asPDB.pdb"

                       # Use mmdb to sanitize the coordinate file
                       import ccp4mg
                       import mmdb2
                       mmdb2.InitMatType()
                       m = mmdb2.Manager()
                       m.SetFlag(mmdb2.MMDBF_IgnoreSegID)
                       m.ReadCoorFile(coordPath)

                       # Remove SEGIDs which can confuse molprobity if heterogeneous
                       sel = m.NewSelection()
                       m.SelectAtoms(sel, 0, "*", mmdb2.ANY_RES, "*", mmdb2.ANY_RES, "*", "*", "*", "*", "*", mmdb2.SKEY_OR)
                       selindexp = mmdb2.intp()
                       selAtoms = mmdb2.GetAtomSelIndex(m, sel, selindexp)
                       nSelAtoms = selindexp.value()
                       for i in range(nSelAtoms):
                           at = mmdb2.getPCAtom(selAtoms, i)
                           at.segID = b"    "
                       m.FinishStructEdit()
                       m.WritePDBASCII(sanitizedCoordPath)

                       fileRoot = os.path.join(self.getWorkDirectory(), "molprobity")
                       molprobity.run(["input.pdb.file_name={}".format(sanitizedCoordPath),
                                       "output.prefix={}".format(fileRoot)])
                       mpCootScriptPath = fileRoot + "_coot.py"
                       mpProbePath = fileRoot + "_probe.txt"
                       if os.path.isfile(mpCootScriptPath):
                           with open(self.container.outputData.COOTSCRIPTOUT.fullPath.__str__(), "a+") as cootScript:
                               cootScript.write("\n")
                               with open(mpCootScriptPath, "r") as mpCootScript:
                                   content = mpCootScript.read()
                                   content = content.replace('"molprobity_probe.txt"', '"{}"'.format(mpProbePath))
                                   cootScript.write(content)

                       mpOutPath = fileRoot + ".out"
                       if os.path.isfile(mpOutPath):
                           etree.SubElement(etree.SubElement(self.xmlroot, "Molprobity"), "Output").text = etree.CDATA(open(mpOutPath).read())
                       self.saveXml()
                   except Exception as err:
                       # Non-fatal - molprobity is optional
                       etree.SubElement(etree.SubElement(self.xmlroot, "Molprobity"), "Output").text = etree.CDATA(str(err))
                       self.saveXml()

           # Stop file watcher before validation to prevent handleXmlChanged from clearing xmlroot
           # (Only in Qt mode - in Django mode we use Signal/Slot, not file watching)
           if not DJANGO() and hasattr(self, 'firstRefmac'):
               firstRefmacXMLFilename = self.firstRefmac.makeFileName(format='PROGRAMXML')
               self.unwatchFile(firstRefmacXMLFilename)

           # Get validation settings - defaults match def.xml (True for all except MOLPROBITY which is test-controlled)
           validate_iris = getattr(self.container.controlParameters, "VALIDATE_IRIS", True)
           validate_baverage = getattr(self.container.controlParameters, "VALIDATE_BAVERAGE", True)
           validate_ramachandran = getattr(self.container.controlParameters, "VALIDATE_RAMACHANDRAN", True)
           validate_molprobity = getattr(self.container.controlParameters, "VALIDATE_MOLPROBITY", True)


           if validate_baverage or validate_molprobity or validate_ramachandran or validate_iris:
             xml_validation = etree.SubElement(self.xmlroot,"Validation")
             xml_validation_status = etree.SubElement(xml_validation,"Success")
             try:
                self.validate = self.makePluginObject('validate_protein')
                self.validate.container.inputData.XYZIN_2.set(self.container.outputData.CIFFILE)
                self.validate.container.inputData.XYZIN_1.set(self.container.inputData.XYZIN)
                if str(self.container.controlParameters.SCATTERING_FACTORS) == "XRAY":
                    self.validate.container.inputData.F_SIGF_2.set(self.container.inputData.F_SIGF)
                    self.validate.container.inputData.F_SIGF_1.set(self.container.inputData.F_SIGF)
                else:
                    self.validate.container.inputData.F_SIGF_2.set(None)
                    self.validate.container.inputData.F_SIGF_1.set(None)
                self.validate.container.inputData.NAME_2 = "Refined"
                self.validate.container.inputData.NAME_1 = "Input"
                #MN...Using "="" to set this is an odd thing and breaks under some circumstances.
                #Specifically, the path for this is now in the prosmart_refmac (i.e. parent job) directory,
                #but it is an output of the validate_protein subjob.  When that subjob seeks to save it to the
                #database, it finds that it is in the wrong directory and croaks.  I have tried to get round this
                #by flagging it as saveToDb=False in the validate_protein .def.xml. But htis needs consideration by jon
                #Agirre and Rob Nicholls
                self.validate.container.outputData.COOTSCRIPTOUT = self.container.outputData.COOTSCRIPTOUT

                self.validate.container.controlParameters.TWO_DATASETS.set(True)
                self.validate.container.controlParameters.DO_IRIS = validate_iris
                self.validate.container.controlParameters.DO_BFACT = validate_baverage
                self.validate.container.controlParameters.DO_RAMA = validate_ramachandran
                self.validate.container.controlParameters.DO_MOLPROBITY = validate_molprobity

                self.validate.doAsync = False
                self.validate.waitForFinished = -1
                self.validate.process()

                validateXMLPath = self.validate.makeFileName('PROGRAMXML')
                validateXML = CCP4Utils.openFileToEtree(validateXMLPath)
                if len(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info"))>0:
                   xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Model_info")[0])
                if self.validate.container.controlParameters.DO_IRIS:
                   if len(validateXML.xpath("//Validate_geometry_CCP4i2/Iris"))>0:
                      xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Iris")[0])
                if self.validate.container.controlParameters.DO_BFACT:
                   bfactors_found = validateXML.xpath("//Validate_geometry_CCP4i2/B_factors")
                   if len(bfactors_found)>0:
                      xml_validation.append(bfactors_found[0])
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
                   except Exception as e:
                       # Non-fatal - verdict analysis is optional
                       pass
                xml_validation_status.text = "SUCCESS"
                self.saveXml()
             except Exception as err:
                xml_validation_status.text = "FAILURE"
                self.saveXml()
                self.appendErrorReport(110, f'Validation failed: {err}')

        logfiles = []
        if hasattr(self,"firstRefmac"):
            logfiles.append(self.firstRefmac.makeFileName('LOG'))
        if hasattr(self,"refmacPostCootPlugin"):
            logfiles.append(self.refmacPostCootPlugin.makeFileName('LOG'))

        asuin = self.container.inputData.ASUIN
        if asuin.isSet():
            self.saveXml()
            try:
                xyzinPath = str(self.container.outputData.XYZOUT)
                self.xmlroot.append(sequenceAlignment(xyzinPath, asuin))
            except Exception as err:
                # Non-fatal - sequence alignment is optional
                pass

        self.createWarningsXML(logfiles)
        self.saveXml()

        self.reportStatus(CPluginScript.SUCCEEDED)

    def tryVariousRefmacWeightsAround(self, weight):
        import math
        # Generate jobs with weights around the initial weight
        # make an array to hold the child-jobs
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

        for job in self.jobsToSubmit:
            if len(self.jobsInTrain) < prosmart_refmac.MAXNJOBS:
                self.submitJob(job)

    def submitJob(self,job):
        rv = job.process()
        # The job instance must be saved to keep it in scope
        self.jobsInTrain[str(job.processId)]=job
        self.jobsToSubmit.remove(job)

    @QtCore.Slot(dict)
    def handleDone(self, ret):
        pid =  ret.get('pid',None)
        status,exitStatus,exitCode = self.postProcessCheck(pid)
        if  status == CPluginScript.FAILED:
            self.reportStatus(status)
            return

        rtask = self.jobsInTrain[str(pid)]

        self.addCycleXML(rtask)
        # Save the xml on every cycle
        aFile=open( self.pipelinexmlfile,'w')
        CCP4Utils.writeXML(aFile, etree.tostring(self.xmlroot,pretty_print=True) )
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
        xmlcyc = etree.SubElement(self.xmlroot,name)
        cycWeightNode = etree.SubElement(xmlcyc,"weight")
        rxml = CCP4Utils.openFileToEtree(rtask.makeFileName('PROGRAMXML'))
        try: cycWeightNode.text = rxml.xpath("//Cycle/WeightUsed")[-1].text
        except: pass
        rstats = rxml.xpath("//REFMAC")
        xmlcyc.append (rstats[0])

    def handleTimeout(self):
        for rtask in self.jobsInTrain:
            try:
                rtask.terminate()
            except Exception:
                pass

        self.appendErrorReport(40,str(self.TIMEOUT_PERIOD))
        self.reportStatus(CPluginScript.FAILED)

# Function called from gui to support exporting MTZ files
def exportJobFile(jobId=None,mode=None,fileInfo={}):
    import os
    from core import CCP4Modules

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

    return None

# Function to return list of names of exportable MTZ(s)
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [[ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ]]

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
