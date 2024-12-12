"""
    servalcat_pipe.py: CCP4 GUI Project
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
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling
from core import CCP4Utils
from . import monitor_differences
from wrappers.servalcat.script.json2xml import json2xml
import os, sys, shutil
import gemmi
import numpy
import json
import statistics
from operator import itemgetter


class servalcat_pipe(CPluginScript):

    TASKMODULE = 'refinement'
    SHORTTASKTITLE = 'Servalcat'
    TASKTITLE = 'Refinement against diffraction data or SPA map & optional restraints from ProSMART & MetalCoord'
    TASKNAME = 'servalcat_pipe'  # Task name - same as class name
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'
    TASKVERSION= 0.1
    WHATNEXT = ['servalcat_pipe','coot_rebuild','modelcraft']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 240
    MAXNJOBS = 4
    PERFORMANCECLASS = 'CServalcatPerformance'
    SUBTASKS=['servalcat','prosmart','metalCoord']
    RUNEXTERNALPROCESS=False
    PURGESEARCHLIST =  [[ 'refmac%*/hklout.mtz', 0, "hklout" ], [ 'refmac%*/hklout.mtz', 7, "hklout" ], [ '*%*/ANOMFPHIOUT.mtz', 1, "ANOMFPHIOUT" ], [ '*%*/DIFANOMFPHIOUT.mtz', 1, "DIFANOMFPHIOUT" ]]


    ERROR_CODES = { 101 : { 'description' : 'Error copying data file from final job to pipeline directory' }
                    }

    def __init__(self, *args, **kws):
        super(servalcat_pipe, self).__init__(*args, **kws)
        self.pipelinexmlfile = self.makeFileName(format='PROGRAMXML')
        self.refmacMonitors = {}
        self.xmlroot = etree.Element("SERVALCAT")
        self.xmlroot2 = etree.Element("SERVALCAT")

    def startProcess(self, processId):
        try:
            self.executeProsmartProtein()
            self.executeProsmartNucleicAcid()
        except Exception as e:
            sys.stderr.write("ERROR while running ProSMART: " + str(e) + "\n")
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        try:
            self.executeMetalCoords()
        except Exception as e:
            sys.stderr.write("ERROR while running MetalCoord: " + str(e) + "\n")
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        try:
            self.executeFirstServalcat()
        except Exception as e:
            sys.stderr.write("ERROR while running Servalcat: " + str(e) + "\n")
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

    def executeMetalCoords(self):
        if self.container.metalCoordPipeline.RUN_METALCOORD and \
                self.container.metalCoordPipeline.GENERATE_OR_USE == "GENERATE":
            if not shutil.which("metalCoord", mode=os.X_OK):
                print("WARNING: MetalCoord will not be executed because it is not installed.")
                return
            self.metalCoordOutputJsonPaths = []
            ligand_codes_selected = self.container.metalCoordPipeline.LIGAND_CODES_SELECTED
            for ligand_code in ligand_codes_selected:
                self.executeMetalCoord(ligand_code)
            print(self.metalCoordOutputJsonPaths)
            if not self.metalCoordOutputJsonPaths:
                print("ERROR: No output from MetalCoord found.")
                return
            # Get all the JSON files together and convert JSON to external restraint keywords
            self.outputRestraintsPrefix = "metal_restraints"
            self.outputRestraintsFilename = self.outputRestraintsPrefix + ".txt"
            self.outputRestraintsMmcifFilename = self.outputRestraintsPrefix + ".mmcif"
            self.outputRestraintsPathPrefix = os.path.join(self.getWorkDirectory(), self.outputRestraintsPrefix)
            self.outputRestraintsPath = os.path.join(self.getWorkDirectory(), self.outputRestraintsFilename)
            self.outputRestraintsMmcifPath = os.path.join(self.getWorkDirectory(), self.outputRestraintsMmcifFilename)
            if self.container.metalCoordPipeline.LINKS == "UPDATE":
                stPath = str(self.container.inputData.XYZIN.fullPath)
                keep_links = False
            elif self.container.metalCoordPipeline.LINKS == "KEEP":
                stPath = str(self.container.inputData.XYZIN.fullPath)
                keep_links = True
            else:  # self.container.metalCoordPipeline.LINKS == "NOTTOUCH":
                stPath = None
                keep_links = True
            print("Converting MetalCoord analyses from JSON files to restraints")
            from wrappers.metalCoord.script import json2restraints
            json2restraints.main(
                jsonPaths=self.metalCoordOutputJsonPaths,
                stPath=stPath,
                outputPrefix=self.outputRestraintsPathPrefix,
                jsonEquivalentsPath=None,
                keep_links=keep_links)
            if os.path.isfile(self.outputRestraintsPath):
                self.container.outputData.METALCOORD_RESTRAINTS.setFullPath(self.outputRestraintsPath)
                self.container.outputData.METALCOORD_RESTRAINTS.annotation = 'Restraints for metal sites'
            if os.path.isfile(self.outputRestraintsMmcifPath) and stPath:
                self.container.outputData.METALCOORD_XYZ.setFullPath(self.outputRestraintsMmcifPath)
                self.container.outputData.METALCOORD_XYZ.annotation = 'Input structure with links from MetalCoord'

    def executeMetalCoord(self, ligand_code):
        self.metalCoordPlugin = self.makePluginObject('metalCoord')
        self.metalCoordPlugin.container.inputData.XYZIN.set(self.container.inputData.XYZIN)
        self.metalCoordPlugin.container.controlParameters.copyData(self.container.metalCoordWrapper.controlParameters)
        self.metalCoordPlugin.container.inputData.LIGAND_CODE.set(ligand_code)
        self.metalCoordPlugin.container.controlParameters.SAVE_PDBMMCIF.set(False)
        if str(self.container.metalCoordPipeline.LINKS) == "KEEP":
            self.metalCoordPlugin.container.controlParameters.KEEP_LINKS.set(True)
        self.connectSignal(self.metalCoordPlugin, 'finished', self.metalCoordFinished)
        self.metalCoordPlugin.waitForFinished = -1
        self.metalCoordPlugin.process()
        self.outputJsonFilename = str(self.metalCoordPlugin.container.inputData.LIGAND_CODE) + ".json"
        self.outputJsonPath = os.path.join(self.metalCoordPlugin.getWorkDirectory(), self.outputJsonFilename)
        if os.path.isfile(self.outputJsonPath):
            self.metalCoordOutputJsonPaths.append(self.outputJsonPath)

    @QtCore.Slot(dict)
    def metalCoordFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.reportStatus(status)

    def executeFirstServalcat(self, withWeight=-1):
        #create wrapper
        self.firstServalcat = self.createServalcatJob(withWeight)
        # Run asynchronously ...this is needed so that commands downstream of process launch
        # (i.e. logwatcher) will be calledbefore process completion
        self.firstServalcat.doAsync = self.doAsync
        self.firstServalcat.connectSignal(self.firstServalcat,'finished',self.firstServalcatFinished)
        # Install xml node for an in progress from core import CCP4ProjectsManager
        self.xmlLength = 0
        # Start process
        firstServalcatXMLFilename = self.firstServalcat.makeFileName(format='PROGRAMXML')
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("executeFirstServalcat, firstServalcatXMLFilename", firstServalcatXMLFilename)
        self.watchFile(firstServalcatXMLFilename, handler=self.handleXmlChanged, minDeltaSize=34, unwatchWhileHandling=True)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.firstServalcat.process()

    @QtCore.Slot(str)
    def handleXmlChanged2(self, xmlFilename2):
        self.xmlroot2.clear()
        servalcatEtree2 = CCP4Utils.openFileToEtree(xmlFilename2) # MM , useLXML=False)
        servalcatXML2 = servalcatEtree2.xpath('//SERVALCAT')
        if len(servalcatXML2) == 1:
            servalcatXML2[0].tag = "SERVALCAT_WATERS"
            self.xmlroot2.append(servalcatXML2[0])
        self.saveXml2()

    @QtCore.Slot(str)
    def handleXmlChanged(self, xmlFilename):
        self.xmlroot.clear()
        servalcatEtree = CCP4Utils.openFileToEtree(xmlFilename) # MM , useLXML=False)
        servalcatXML = servalcatEtree.xpath("//SERVALCAT")
        if len(servalcatXML) == 1:
            servalcatXML[0].tag = "SERVALCAT_FIRST"
            self.xmlroot.append(servalcatXML[0])
        self.saveXml()

    def saveXml2(self):
        newXml = etree.tostring(self.xmlroot2, pretty_print=True)
        if len(newXml) > self.xmlLength2:
            # Get content of program.xml_first
            firstFileName = self.pipelinexmlfile + '_first'
            with open(firstFileName, 'r') as aFile:
                masterXml = etree.fromstring(aFile.read())
            # Add content of program.xml of the last servalcat subjob
            masterXml.xpath('//SERVALCAT')[0].append(self.xmlroot2.xpath("//SERVALCAT/SERVALCAT_WATERS")[0])
            self.xmlLength2 = len(newXml)
            # Save as program.xml_tmp and then move to program.xml
            tmpFileName = self.pipelinexmlfile + '_tmp'
            with open(tmpFileName, 'w') as aFile:
                CCP4Utils.writeXML(aFile,etree.tostring(masterXml, pretty_print=True) )
            shutil.move(tmpFileName, self.pipelinexmlfile)
            self.xmlroot = masterXml

    def saveXml(self):
        # Save the xml if it has grown
        newXml = etree.tostring(self.xmlroot, pretty_print=True)
        if len(newXml) > self.xmlLength:
           # Save as program.xml_tmp and then move to program.xml
           tmpFileName = self.pipelinexmlfile + '_tmp'
           with open(tmpFileName, 'w') as aFile:
               CCP4Utils.writeXML(aFile, newXml)
           shutil.move(tmpFileName, self.pipelinexmlfile)
           self.xmlLength = len(newXml)

    def createServalcatJob(self, withWeight=-1, inputCoordinates=None, ncyc=-1):
        result = self.makePluginObject('servalcat')
        #input data for this servalcat instance is the same as the input data for the program
        result.container.inputData.copyData(self.container.inputData)

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

        if self.container.metalCoordPipeline.RUN_METALCOORD:
            if self.container.metalCoordPipeline.GENERATE_OR_USE == "GENERATE":
                result.container.inputData.METALCOORD_RESTRAINTS=self.container.outputData.METALCOORD_RESTRAINTS
                if self.container.outputData.METALCOORD_XYZ and self.container.metalCoordPipeline.LINKS != "NOTTOUCH":
                    if os.path.isfile(str(self.container.outputData.METALCOORD_XYZ.fullPath)):
                        result.container.inputData.XYZIN.set(self.container.outputData.METALCOORD_XYZ)
                    # else report error?
            else:
                result.container.inputData.METALCOORD_RESTRAINTS=self.container.metalCoordPipeline.METALCOORD_RESTRAINTS
        if self.container.prosmartProtein.TOGGLE:
            result.container.controlParameters.PROSMART_PROTEIN_SGMN=self.container.prosmartProtein.SGMN
            result.container.controlParameters.PROSMART_PROTEIN_SGMX=self.container.prosmartProtein.SGMX
            result.container.controlParameters.PROSMART_PROTEIN_ALPHA=self.container.prosmartProtein.ALPHA
            result.container.controlParameters.PROSMART_PROTEIN_DMAX=self.container.prosmartProtein.DMAX
            result.container.inputData.PROSMART_PROTEIN_RESTRAINTS=self.prosmart_protein.container.outputData.RESTRAINTS
        if self.container.prosmartNucleicAcid.TOGGLE:
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMN=self.container.prosmartNucleicAcid.SGMN
            result.container.controlParameters.PROSMART_NUCLEICACID_SGMX=self.container.prosmartNucleicAcid.SGMX
            result.container.controlParameters.PROSMART_NUCLEICACID_ALPHA=self.container.prosmartNucleicAcid.ALPHA
            result.container.controlParameters.PROSMART_NUCLEICACID_DMAX=self.container.prosmartNucleicAcid.DMAX
            result.container.inputData.PROSMART_NUCLEICACID_RESTRAINTS=self.prosmart_nucleicacid.container.outputData.RESTRAINTS

        #Specify weight if a meaningful one has been offered
        if withWeight>=0.:
            result.container.controlParameters.WEIGHT_OPT='MANUAL'
            result.container.controlParameters.WEIGHT = withWeight

        if inputCoordinates is not None:
            result.container.inputData.XYZIN.set(inputCoordinates)
            if ncyc > 0:
                result.container.controlParameters.NCYCLES.set(ncyc)
        return result


    def multimericValidation(self):
        # Geometry validation
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

        if validate_iris or validate_baverage or validate_molprobity or validate_ramachandran:
            self.validate = self.makePluginObject('validate_protein')
            self.validate.container.inputData.XYZIN_1.set(self.container.outputData.XYZOUT)
            self.validate.container.inputData.F_SIGF_1.set(self.container.inputData.HKLIN)
            self.validate.container.inputData.XYZIN_2.set(self.container.outputData.XYZOUT)
            self.validate.container.inputData.F_SIGF_2.set(self.container.inputData.HKLIN)
            self.validate.container.inputData.NAME_1.set("Refined")
            self.validate.container.inputData.NAME_2.set("Input")

            self.validate.container.controlParameters.DO_IRIS.set(validate_iris)
            self.validate.container.controlParameters.DO_BFACT.set(validate_baverage)
            self.validate.container.controlParameters.DO_RAMA.set(validate_ramachandran)
            self.validate.container.controlParameters.DO_MOLPROBITY.set(validate_molprobity)

            self.validate.doAsync = False
            self.validate.waitForFinished = -1
            self.validate.process()

            validateXMLPath = self.validate.makeFileName('PROGRAMXML')
            validateXML = CCP4Utils.openFileToEtree(validateXMLPath)
            xml_validation = etree.SubElement(self.xmlroot,"Validation")
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
    

    def adp_analysis(self, modelPath, iqrFactor=2.0):
        adp_dict = {}
        adp_per_resi = {}
        adp_dict["All"] = []
        st = gemmi.read_structure(modelPath)
        for model in st:
            for chain in model:
                polymer = chain.get_polymer()
                ptype = polymer.check_polymer_type()
                adp_dict[chain.name] = []
                adp_per_resi[chain.name] = [[],[],[]]  # residue.seqid.num; mean ADP in backbone; mean ADP in side chain
                for residue in chain:
                    adp_this_resi = []
                    adp_this_resi_sidechain = []
                    for atom in residue:
                        if not atom.is_hydrogen() and atom.occ > 0:
                            if atom.aniso.nonzero():
                                adp_atom = gemmi.calculate_b_est(atom)
                            else:
                                adp_atom = atom.b_iso
                            adp_dict["All"].append(adp_atom)
                            adp_dict[chain.name].append(adp_atom)
                            if ptype in [gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD] and \
                                    atom.name not in ["CA", "C", "O", "N", "OXT"]:
                                adp_this_resi_sidechain.append(adp_atom)
                            else:
                                adp_this_resi.append(adp_atom)
                    try:
                        adp_per_resi[chain.name][0].append(residue.seqid.num)  # ignoring insertion code, sorry
                        if adp_this_resi:
                            adp_per_resi[chain.name][1].append(numpy.mean(adp_this_resi))
                        else:
                            adp_per_resi[chain.name][1].append(None)
                        if adp_this_resi_sidechain:
                            adp_per_resi[chain.name][2].append(numpy.mean(adp_this_resi_sidechain))
                        else:
                            adp_per_resi[chain.name][2].append(None)
                    except:
                        pass

        # Find ADP values which are too small or large - outliers
        q1 = numpy.quantile(adp_dict["All"], 0.25)
        q3 = numpy.quantile(adp_dict["All"], 0.75)
        iqr = q3 - q1
        adp_limit_low = q1 - iqrFactor * iqr
        adp_limit_high = q3 + iqrFactor * iqr
        adp_low = []
        adp_high = []
        for model in st:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.element != gemmi.Element('H') and atom.occ > 0:
                            if atom.aniso.nonzero():
                                adp_atom = gemmi.calculate_b_est(atom)
                            else:
                                adp_atom = atom.b_iso
                            if adp_atom < adp_limit_low:
                                adp_low.append({"atom": str(model.get_cra(atom)),
                                                "adp": adp_atom})
                            elif adp_atom > adp_limit_high:
                                adp_high.append({"atom": str(model.get_cra(atom)),
                                                 "adp": adp_atom})
        adp_low = sorted(adp_low, key=itemgetter('adp'))
        adp_high = sorted(adp_high, key=itemgetter('adp'), reverse=True)

        # Write the analysis in XML
        adp_root = etree.Element('ADP_ANALYSIS')
        chains_root = etree.SubElement(adp_root, "chains")
        for ch, values in adp_dict.items():
            chain = etree.SubElement(chains_root, "chain")
            chain_name = etree.SubElement(chain, "name")
            chain_name.text = str(ch)
            chain.set("name", str(ch))
            chain_min = etree.SubElement(chain, "min")
            chain_min.text = "{:.2f}".format(min(values))
            chain_max = etree.SubElement(chain, "max")
            chain_max.text = "{:.2f}".format(max(values))
            chain_med_val = numpy.median(values)
            chain_mad_val = numpy.median(numpy.absolute(values - chain_med_val))
            chain_med = etree.SubElement(chain, "med")
            chain_med.text = "{:.2f}".format(chain_med_val)
            chain_mad = etree.SubElement(chain, "mad")
            chain_mad.text = "{:.2f}".format(chain_mad_val)
            chain_q1 = etree.SubElement(chain, "q1")
            chain_q1.text = "{:.2f}".format(numpy.quantile(values, 0.25))
            chain_q3 = etree.SubElement(chain, "q3")
            chain_q3.text = "{:.2f}".format(numpy.quantile(values, 0.75))
            chain_mean_val = numpy.mean(values)
            chain_std_val = numpy.std(values)
            chain_mean = etree.SubElement(chain, "mean")
            chain_mean.text = "{:.2f}".format(chain_mean_val)
            chain_std = etree.SubElement(chain, "std")
            chain_std.text = "{:.2f}".format(chain_std_val)

            hist, bin_edges = numpy.histogram(values, bins="auto")
            bin_edges = numpy.delete(bin_edges, -1)
            chain_histogram = etree.SubElement(chain, 'histogram')
            for i in range(len(bin_edges)):
                bin_elem = etree.SubElement(chain_histogram, 'bin')
                bin_adp = etree.SubElement(bin_elem, 'adp')
                bin_adp.text = "{:.2f}".format(bin_edges[i])
                bin_count = etree.SubElement(bin_elem, 'count')
                bin_count.text = str(hist[i])

            if ch != "All":
                chain_per_resi = etree.SubElement(chain, "per_resi")
                for i in range(len(adp_per_resi[str(ch)][0])):
                    adp_per_resi_elem = etree.SubElement(chain_per_resi, "data")
                    adp_per_resi_resi = etree.SubElement(adp_per_resi_elem, 'resi')
                    adp_per_resi_resi.text = str(adp_per_resi[str(ch)][0][i])
                    adp_per_resi_adp = etree.SubElement(adp_per_resi_elem, 'adp')
                    if adp_per_resi[str(ch)][1][i]:
                        adp_per_resi_adp.text = "{:.2f}".format(adp_per_resi[str(ch)][1][i])
                    else:
                        adp_per_resi_adp.text = "-"
                    adp_per_resi_adp_sidechain = etree.SubElement(adp_per_resi_elem, 'adp_sidechain')
                    if adp_per_resi[str(ch)][2][i]:
                        adp_per_resi_adp_sidechain.text = "{:.2f}".format(adp_per_resi[str(ch)][2][i])
                    else:
                        adp_per_resi_adp_sidechain.text = "-"

        outliers_root = etree.SubElement(adp_root, "outliers")
        outliers_adp_limit_low = etree.SubElement(outliers_root, "adp_limit_low")
        outliers_adp_limit_low.text = "{:.2f}".format(adp_limit_low)
        outliers_adp_limit_high = etree.SubElement(outliers_root, "adp_limit_high")
        outliers_adp_limit_high.text = "{:.2f}".format(adp_limit_high)
        outliers_adp_iqr_factor = etree.SubElement(outliers_root, "iqr_factor")
        outliers_adp_iqr_factor.text = "{:.2f}".format(iqrFactor)
        outliers_low = etree.SubElement(outliers_root, "low")
        for outlier in adp_low:
            outlier_elem = etree.SubElement(outliers_low, "data")
            outlier_adp = etree.SubElement(outlier_elem, "adp")
            outlier_adp.text = "{:.2f}".format(outlier["adp"])
            outlier_atom = etree.SubElement(outlier_elem, "atom")
            outlier_atom.text = str(outlier["atom"])
        outliers_high = etree.SubElement(outliers_root, "high")
        for outlier in adp_high:
            outlier_elem = etree.SubElement(outliers_high, "data")
            outlier_adp = etree.SubElement(outlier_elem, "adp")
            outlier_adp.text = "{:.2f}".format(outlier["adp"])
            outlier_atom = etree.SubElement(outlier_elem, "atom")
            outlier_atom.text = str(outlier["atom"])

        aFile = open(self.pipelinexmlfile, 'r')
        oldXml = etree.fromstring(aFile.read())
        aFile.close()
        oldXml.append(adp_root)
        aFile = open(self.pipelinexmlfile + '_tmp', 'w')
        CCP4Utils.writeXML(aFile,etree.tostring(oldXml, pretty_print=True))
        aFile.close()
        shutil.move(self.pipelinexmlfile + '_tmp', self.pipelinexmlfile)


    def coord_adp_dev_analysis(self, model1Path, model2Path):
        import io
        import pandas
        try:
            coordDevMinReported = self.container.monitor.MIN_COORDDEV
            ADPAbsDevMinReported = self.container.monitor.MIN_ADPDEV
            csv_string = monitor_differences.main(
                file1=model1Path, file2=model2Path, output=None,
                minCoordDev=float(coordDevMinReported), minADPDev=float(ADPAbsDevMinReported))
            # Load
            csvStringIO = io.StringIO(csv_string)
            df = pandas.read_csv(csvStringIO, sep=",", header=0)
            coordDevMean = df["CoordDev"].mean()
            ADPAbsDevMean = df["ADPDev"].mean()
            # Save csv in program.xml
            xmlText = "\n<COORD_ADP_DEV>"
            xmlText += "\n<STATISTICS>"
            xmlText += "\n<coordDevMean>"
            xmlText += str(round(coordDevMean, 2))
            xmlText += "</coordDevMean>"
            xmlText += "\n<coordDevMinReported>"
            xmlText += str(round(coordDevMinReported, 2))
            xmlText += "</coordDevMinReported>"
            xmlText += "\n<ADPAbsDevMean>"
            xmlText += str(round(ADPAbsDevMean, 2))
            xmlText += "</ADPAbsDevMean>"
            xmlText += "\n<coordADPAbsMinReported>"
            xmlText += str(round(ADPAbsDevMinReported, 2))
            xmlText += "</coordADPAbsMinReported>"
            xmlText += "\n</STATISTICS>"
            xmlText += "\n<CSV><![CDATA[\n" + csv_string + "\n]]></CSV>"
            xmlText += "\n</COORD_ADP_DEV>"
            print(xmlText)
            xmlTree = etree.fromstring(xmlText)
            aFile = open(self.pipelinexmlfile, 'r')
            oldXml = etree.fromstring(aFile.read())
            aFile.close()
            oldXml.append(xmlTree)
            aFile = open(self.pipelinexmlfile + '_tmp', 'w')
            CCP4Utils.writeXML(aFile, etree.tostring(oldXml, pretty_print=True))
            aFile.close()
            shutil.move(self.pipelinexmlfile + '_tmp', self.pipelinexmlfile)
        except Exception as e:
            sys.stderr.write("Monitoring of the changes in coordinates and ADPs was not successful: " + str(e) + "\n")

    @QtCore.Slot(dict)
    def firstServalcatFinished(self, statusDict):
        print("AAA1")
        if statusDict['finishStatus'] == CPluginScript.UNSATISFACTORY:
            print("AAA1.UNSATISFACTORY")
            with open(self.makeFileName('PROGRAMXML'), 'w') as programXML:
                CCP4Utils.writeXML(programXML, etree.tostring(self.xmlroot, pretty_print=True))
            self.reportStatus(CPluginScript.UNSATISFACTORY)

        elif self.firstServalcat.errorReport.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            print("AAA1.MAXSEVERITY")
            #This gets done in thefirstServalcat.reportStatus() - Liz
            try:
                servalcatEtree = CCP4Utils.openFileToEtree(self.firstServalcat.makeFileName('PROGRAMXML'))
                servalcatXML = servalcatEtree.xpath('//SERVALCAT')
                if len(servalcatXML) == 1: self.xmlroot.append(servalcatXML[0])
            except:
                print('Failed attempt to read XML file from first Refmac')
            try:
                newXml = etree.tostring(self.xmlroot,pretty_print=True)
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
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            # This gets done in the firstServalcat.reportStatus() - Liz
            self.fileSystemWatcher = None
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            print("AAA12")
            # self.addCycleXML(self.firstServalcat) # MM
            aFile=open( self.pipelinexmlfile,'w')
            CCP4Utils.writeXML(aFile, etree.tostring(self.xmlroot, pretty_print=True) )
            aFile.close()
            print("AAA13")
            if self.container.controlParameters.OPTIMISE_WEIGHT:
                print("AAA14")
                self.fileSystemWatcher = None
                weightUsed = float(self.xmlroot.xpath('//weight')[-1].text)
                self.tryVariousRefmacWeightsAround(weightUsed)
            else:
               print("AAA15")
               print("AAA15.1")
               if self.container.controlParameters.ADD_WATERS:
                   # Coot sujob to add waters
                   print("AAA16")
                   self.currentCoordinates = self.firstServalcat.container.outputData.CIFFILE
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
        print('done servalcat_pipe.firstServalcatFinished')

    def makeCootPlugin(self):
         # FIXME - This is all nonsense - needs to consider best task, etc... *NOT* just firstRefmaca?
        cootPlugin = self.makePluginObject('coot_find_waters')
        cootPlugin.container.inputData.XYZIN.set(self.currentCoordinates)
        cootPlugin.container.inputData.FPHIIN.set(self.firstServalcat.container.outputData.FPHIOUT)
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

    @QtCore.Slot(dict)
    def cootFinished(self, statusDict={}):
        import functools
        # Check coot status and start servalcat
        self.checkFinishStatus(statusDict=statusDict, failedErrCode=204, outputFile=self.cootPlugin.container.outputData.XYZOUT, noFileErrCode=205)
        try:
          if self.container.controlParameters.ADD_WATERS:
            aFile = open(self.pipelinexmlfile,'r')
            oldXml = etree.fromstring(aFile.read())
            aFile.close()
            nwaters = "unknown"
            cootLogTxt = os.path.join(os.path.dirname(self.cootPlugin.container.outputData.XYZOUT.__str__()), "log.txt")
            with open(cootLogTxt, 'r') as f:
               for l in f:
                   if l.startswith("INFO::") and "found" in l and "water fitting" in l:
                      nwaters = l.strip()
                      numsearch = [ x for x in nwaters.split() if x.isdigit() ]
                      if len(numsearch)>0:
                         nwaters = numsearch[0]
                      break
            postRefmacCoot = etree.Element("CootAddWaters")
            postRefmacCoot.text = "Coot added " + nwaters + " water molecules."
            oldXml.append(postRefmacCoot)
            aFile = open(self.pipelinexmlfile+'_tmpcoot','w')
            CCP4Utils.writeXML(aFile,etree.tostring(oldXml,pretty_print=True))
            aFile.close()
            shutil.move(self.pipelinexmlfile+'_tmpcoot', self.pipelinexmlfile)
          self.cootPlugin.container.outputData.XYZOUT.subType = 1
          self.currentCoordinates = self.cootPlugin.container.outputData.XYZOUT
          self.servalcatPostCootPlugin = self.createServalcatJob(inputCoordinates=self.currentCoordinates, ncyc=int(self.container.controlParameters.NCYCLES_AFTER_ADD_WATERS))
          self.servalcatPostCootPlugin.doAsync = True
          self.servalcatPostCootPlugin.connectSignal(self.servalcatPostCootPlugin,'finished',functools.partial(self.postCootServalcatFinished,self.servalcatPostCootPlugin))
          servalcatPostCootXMLFilename = self.servalcatPostCootPlugin.makeFileName(format='PROGRAMXML')
          self.xmlLength2 = 0
          self.watchFile(servalcatPostCootXMLFilename, handler=self.handleXmlChanged2, minDeltaSize=34, unwatchWhileHandling=True)
          shutil.copyfile(self.pipelinexmlfile, self.pipelinexmlfile + "_first")
          rv = self.servalcatPostCootPlugin.process()
          if rv == CPluginScript.FAILED: self.reportStatus(rv)
        except Exception as e:
          self.appendErrorReport(CPluginScript,39,str(e))

    @QtCore.Slot('CPluginScript', dict)
    def postCootServalcatFinished(self, servalcatJob, statusDict={}):
        self.handleXmlChanged2(servalcatJob.makeFileName(format='PROGRAMXML'))
        self.fileSystemWatcher = None
        if statusDict['finishStatus'] == CPluginScript.FAILED:
            #This gets done in the firstServalcat.reportStatus() - Liz
            self.reportStatus(CPluginScript.FAILED)
            return
        else:
            self.handleXmlChanged2(self.servalcatPostCootPlugin.makeFileName(format='PROGRAMXML'))
        self.finishUp(servalcatJob)

    def finishUp(self, servalcatJob):
        from core import CCP4ProjectsManager
        print('into servalcat_pipe.finishUp')
        for attr in self.container.outputData.dataOrder():
            try:
                wrappersAttr = getattr(servalcatJob.container.outputData, attr)
                pipelinesAttr = getattr(self.container.outputData, attr)
            except:
                print('servalcat_pipe.finishUp attr', attr, 'not copied from wrapper to pipeline')
                continue
            print('servalcat_pipe.finishUp attr', attr)
            if attr in ["PERFORMANCEINDICATOR"]:
                setattr(self.container.outputData, attr, wrappersAttr)
            else:
                if os.path.exists(str(wrappersAttr.fullPath)):
                  try:
                    shutil.copyfile( str(wrappersAttr.fullPath), str(pipelinesAttr.fullPath) )
                  except:
                    self.appendErrorReport(101,str(wrappersAttr.fullPath)+' to '+str(pipelinesAttr.fullPath))
                if attr == "XMLOUT":
                    pass
        print('servalcat_pipe.finishUp 1')
        # Apply database annotations
        self.container.outputData.XYZOUT.annotation.set(servalcatJob.container.outputData.XYZOUT.annotation)
        self.container.outputData.CIFFILE.annotation.set(servalcatJob.container.outputData.CIFFILE.annotation)
        self.container.outputData.FPHIOUT.annotation.set(servalcatJob.container.outputData.FPHIOUT.annotation)
        self.container.outputData.FPHIOUT.subType = servalcatJob.container.outputData.FPHIOUT.subType
        self.container.outputData.DIFFPHIOUT.annotation.set(servalcatJob.container.outputData.DIFFPHIOUT.annotation)
        self.container.outputData.DIFFPHIOUT.subType = servalcatJob.container.outputData.DIFFPHIOUT.subType
        if self.container.outputData.DICT.exists():
            self.container.outputData.DICT.annotation = 'Accumulated monomer dictionary'
        if servalcatJob.container.outputData.COOTSCRIPTOUT.exists():
            if servalcatJob.container.outputData.COOTSCRIPTOUT.annotation.isSet():
                self.container.outputData.COOTSCRIPTOUT.annotation.set(servalcatJob.container.outputData.COOTSCRIPTOUT.annotation)
        if str(servalcatJob.container.controlParameters.DATA_METHOD) == 'xtal':
            if servalcatJob.container.outputData.ANOMFPHIOUT.exists():
                if servalcatJob.container.outputData.ANOMFPHIOUT.annotation.isSet():
                    self.container.outputData.ANOMFPHIOUT.annotation.set(servalcatJob.container.outputData.ANOMFPHIOUT.annotation)
            if servalcatJob.container.outputData.DIFANOMFPHIOUT.exists():
                if servalcatJob.container.outputData.DIFANOMFPHIOUT.annotation.isSet():
                    self.container.outputData.DIFANOMFPHIOUT.annotation.set(servalcatJob.container.outputData.DIFANOMFPHIOUT.annotation)
        elif str(servalcatJob.container.controlParameters.DATA_METHOD) == 'spa':
            self.container.outputData.MAP_FO.annotation.set(servalcatJob.container.outputData.MAP_FO.annotation)
            self.container.outputData.MAP_FO.subType = servalcatJob.container.outputData.MAP_FO.subType
            self.container.outputData.MAP_FOFC.annotation.set(servalcatJob.container.outputData.MAP_FOFC.annotation)
            self.container.outputData.MAP_FOFC.subType = servalcatJob.container.outputData.MAP_FOFC.subType
        print('servalcat_pipe.finishUp 3'); sys.stdout.flush()

        cleanUpIntermediate = False
        if hasattr(self.container.controlParameters,"REFMAC_CLEANUP"):
            cleanUpIntermediate = self.container.controlParameters.REFMAC_CLEANUP
            if cleanUpIntermediate:
                print('servalcat_pipe.finishUp 4'); sys.stdout.flush()
                cleanup = CCP4ProjectsManager.CPurgeProject(self.firstServalcat._dbProjectId)
                print('servalcat_pipe.finishUp 5'); sys.stdout.flush()
                cleanup.purgeJob(self.firstServalcat.jobId,context="extended_intermediate",reportMode="skip")

                if hasattr(self,"servalcatPostCootPlugin"):
                    print('servalcat_pipe.finishUp 6'); sys.stdout.flush()
                    cleanup = CCP4ProjectsManager.CPurgeProject(self.servalcatPostCootPlugin._dbProjectId)
                    print('servalcat_pipe.finishUp 7'); sys.stdout.flush()
                    cleanup.purgeJob(self.servalcatPostCootPlugin.jobId,context="extended_intermediate",reportMode="skip")

        self.multimericValidation()
        if self.container.controlParameters.RUN_ADP_ANALYSIS:
            self.adp_analysis(
                str(self.container.outputData.CIFFILE.fullPath),
                float(self.container.controlParameters.ADP_IQR_FACTOR))
        if self.container.monitor.RUN_COORDADPDEV_ANALYSIS:
            self.coord_adp_dev_analysis(
                str(self.container.inputData.XYZIN.fullPath),
                str(self.container.outputData.CIFFILE.fullPath))
        self.saveXml()

        print('done servalcat_pipe.finishUp'); sys.stdout.flush()
        self.reportStatus(CPluginScript.SUCCEEDED)


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
    import clipper
    mtz_file = clipper.CCP4MTZfile()
    hkl_info = clipper.HKL_info()
    mtz_file.open_read (str(coefficientsPath))
    mtz_file.import_hkl_info ( hkl_info )
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
    from core import CCP4Modules

    theDb = CCP4Modules.PROJECTSMANAGER().db()
    if mode == 'complete_mtz':
        childJobs = theDb.getChildJobs(jobId=jobId,details=True)
        if childJobs[-1][2] == 'servalcat':
          jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=childJobs[-1][1],create=False)
          if os.path.exists(os.path.join(jobDir,'refined.mtz')):
             return  os.path.join(jobDir,'refined.mtz')
        elif childJobs[-1][2] == 'validate_protein':
          if childJobs[-2][2] == 'servalcat':
             jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=childJobs[-2][1],create=False)
             if os.path.exists(os.path.join(jobDir,'refined.mtz')):
                return  os.path.join(jobDir,'refined.mtz')


# Function to return list of names of exportable MTZ(s)
def exportJobFileMenu(jobId=None):
    # Return a list of items to appear on the 'Export' menu - each has three subitems:
    # [ unique identifier - will be mode argument to exportJobFile() , menu item , mime type (see CCP4CustomMimeTypes module) ]
    return [[ 'complete_mtz' ,'MTZ file' , 'application/CCP4-mtz' ]]
