from __future__ import print_function

"""
    buccaneer_build_refine_mr.py: CCP4 GUI Project
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

import inspect
import sys
import time
import os
import shutil
from lxml import etree
from copy import deepcopy
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from core import CCP4ModelData
from core import CCP4ProjectsManager

class buccaneer_build_refine_mr(CPluginScript):

    TASKMODULE = 'model_building'
    TASKTITLE = 'Autobuild protein (Buccaneer pipeline)'
    TASKNAME = 'buccaneer_build_refine_mr'
    TASKVERSION = 202003271051
    INTERRUPTABLE = True
    GUINAME = 'buccaneer_build_refine_mr'
    RESTARTABLE = True
    ASYNCHRONOUS = False
    PERFORMANCECLASS = 'CModelBuildPerformance'
    MAINTAINER = 'paul.bond@york.ac.uk'
    RUNEXTERNALPROCESS = False
    WHATNEXT = ['coot_rebuild', 'prosmart_refmac', 'buccaneer_build_refine_mr']
    PURGESEARCHLIST = [['buccaneer_mr%*/XYZOUT.*', 5, 'XYZOUT'],
                       ['refmac%*/ABCDOUT.mtz', 5, 'ABCDOUT'],
                       ['refmac%*/FPHIOUT.mtz', 5, 'FPHIOUT'],
                       ['refmac%*/DIFFPHIOUT.mtz', 5, 'DIFFPHIOUT'],
                       ['refmac%*/XYZOUT.pdb', 5, 'XYZOUT'],
                       ['refmac%*/COOTSCRIPTOUT.scm', 5, 'COOTSCRIPTOUT'],
                       ['prosmart_refmac%*/ABCDOUT.mtz', 5, 'ABCDOUT'],
                       ['prosmart_refmac%*/FPHIOUT.mtz', 5, 'FPHIOUT'],
                       ['prosmart_refmac%*/DIFFPHIOUT.mtz', 5, 'DIFFPHIOUT'],
                       ['prosmart_refmac%*/XYZOUT.pdb', 5, 'XYZOUT'],
                       ['prosmart_refmac%*/COOTSCRIPTOUT.scm', 5, 'COOTSCRIPTOUT'],
                       ['prosmart_refmac%*/refmac%*/ABCDOUT.mtz', 5, 'ABCDOUT'],
                       ['prosmart_refmac%*/refmac%*/FPHIOUT.mtz', 5, 'FPHIOUT'],
                       ['prosmart_refmac%*/refmac%*/DIFFPHIOUT.mtz', 5, 'DIFFPHIOUT'],
                       ['prosmart_refmac%*/refmac%*/XYZOUT.pdb', 5, 'XYZOUT'],
                       ['prosmart_refmac%*/refmac%*/COOTSCRIPTOUT.scm', 5, 'COOTSCRIPTOUT']]
    ERROR_CODES = {200: {'description': 'Buccaneer task failed'},
                   201: {'description': 'Buccaneer task failed to produce a required output file'},
                   202: {'description': 'Refmac task failed'},
                   203: {'description': 'Refmac task failed to produce a required output file'},
                   204: {'description': 'Coot realspace operation task failed'},
                   205: {'description': 'Coot realspace operation task failed to produce a required output file'}}

    def __init__(self, *args, **kws):
        super(buccaneer_build_refine_mr, self).__init__(*args, **kws)
        if self.container.controlParameters.BUCCANEER_PHSIN_TYPE == 'mr':
            self.container.inputData.ABCD_EXP.setQualifier('allowUndefined', True)
            self.container.inputData.FREERFLAG.setQualifier('allowUndefined', True)

# ===== Pipeline methods =======================================================

    def startProcess(self, processId):
        self.initialiseProperties()
        self.initialiseXML()
        if self.restarted: self.handleRestart()

        if self.shouldCalculateModelPhases():
            self.refmacABCD()

        for self.cycle in range(self.cycle, self.ncycles):
            self.processedPlugins[self.cycle] = []
            self.cleanUpIntermediateFiles()
            self.addXMLCycle()

            print('Starting cycle', self.cycle + 1, 'of', self.ncycles)
            self.buildRefineCycle()
            print('Finished cycle', self.cycle + 1, 'of', self.ncycles)

            self.checkForImprovement()
            if self.auto_stop and self.cycles_without_improvement == 4: break

            if self.testForInterrupt() and self.cycle + 1 < self.ncycles:
                self.alignment()
                self.container.interruptStatus.LASTCYCLE = self.cycle
                self.reportStatus(CPluginScript.INTERRUPTED)

        self.writeFinalXML()
        self.alignment()
        if self.container.controlParameters.BUCCANEER_CLEANUP:
            cleanup = CCP4ProjectsManager.CPurgeProject(self._dbProjectId)
            cleanup.purgeJob(self.jobId, context="extended_intermediate", reportMode="skip")
        self.reportStatus(CPluginScript.SUCCEEDED)

    def buildRefineCycle(self):
        if self.full_prune and self.cycle > 0:
            self.cootPrune()
            self.refmac(cycles=5)
        self.buccaneer()
        if self.chain_prune or self.shouldRunCootRSO():
            self.refmac()
            if self.chain_prune:
                self.cootPrune(chains_only=True)
            if self.shouldRunCootRSO():
                self.cootRSO()
            self.refmac(cycles=5, final=True)
        else:
            self.refmac(final=True)

    def initialiseProperties(self):
        self.cycle = 0
        self.ncycles = int(self.container.controlParameters.ITERATIONS)
        self.auto_stop = self.container.controlParameters.STOP_AUTOMATICALLY
        self.full_prune = self.container.controlParameters.FULL_PRUNE
        self.chain_prune = self.container.controlParameters.CHAIN_PRUNE
        self.best_cycle = 0
        self.rwork = 1
        self.rfree = 1
        self.completeness_by_res = 0
        self.completeness_by_chn = 0
        self.n_fragments = 999
        self.longest_fragment = 1
        self.residues_built = 1
        self.residues_sequenced = 1
        self.min_rwork = self.rwork
        self.min_rfree = self.rfree
        self.min_n_fragments = self.n_fragments
        self.max_longest_fragment = self.longest_fragment
        self.max_residues_built = self.residues_built
        self.max_residues_sequenced = self.residues_sequenced
        self.cycles_without_improvement = 0
        self.essentiallyComplete = False
        self.currentXYZ = self.container.inputData.XYZIN
        self.currentABCD = self.container.inputData.ABCD
        self.processedPlugins = {}
        self.restarted = self.container.interruptStatus.LASTCYCLE.isSet()

    def handleRestart(self):
        self.cycle = int(self.container.interruptStatus.LASTCYCLE) + 1
        self.currentXYZ = os.path.join(self.getSubDirectories()[-1], 'XYZOUT.pdb')
        self.currentABCD = os.path.join(self.getSubDirectories()[-1], 'ABCDOUT.mtz')
        self.currentFPHI = os.path.join(self.getSubDirectories()[-1], 'FPHIOUT.mtz')
        self.currentDIFFPHI = os.path.join(self.getSubDirectories()[-1], 'DIFFPHIOUT.mtz')
        print('Restarting from cycle', self.cycle, 'with coordinates in file:', self.currentXYZ)
        self.parseXmlFromBeforeInterrupt()

    def shouldCalculateModelPhases(self):
        return (self.cycle == 0 and
            self.container.controlParameters.BUCCANEER_PHSIN_TYPE == "mr" and
            not self.container.inputData.ABCD.isSet())

    def shouldRunCootRSO(self):
        return (self.container.controlParameters.COOT_REALSPACE_OPERATION != "none" and
            self.rwork < self.container.controlParameters.BUCCANEER_RSR_RWORK_LIMIT)

    def improved(self):
        required_improvement = 0.02
        improvement = (self.min_rwork - self.rwork) / self.min_rwork
        if improvement > required_improvement: return True
        improvement = (self.min_n_fragments - self.n_fragments) / float(self.min_n_fragments)
        if improvement > required_improvement: return True
        improvement = (self.longest_fragment - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement: return True
        improvement = (self.residues_built - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement: return True
        improvement = (self.residues_sequenced - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement: return True
        return False

    def checkForImprovement(self):
        if self.improved(): self.cycles_without_improvement = 0
        else: self.cycles_without_improvement += 1

        if self.rwork < self.min_rwork: self.min_rwork = self.rwork
        if self.n_fragments < self.min_n_fragments: self.min_n_fragments = self.n_fragments
        if self.longest_fragment > self.max_longest_fragment: self.max_longest_fragment = self.longest_fragment
        if self.residues_built > self.max_residues_built: self.max_residues_built = self.residues_built
        if self.residues_sequenced > self.max_residues_sequenced: self.max_residues_sequenced = self.residues_sequenced
        if self.rfree < self.min_rfree:
            self.min_rfree = self.rfree
            self.setOutput()
            self.best_cycle = self.cycle

    def setOutput(self):
        self.copyPipelineOutputFiles()
        self.setPerformanceData()

# ===== Sub-job methods ========================================================

    def refmacABCD(self):
        if self.container.controlParameters.USE_SHIFTFIELD:
          plugin = self.makeSynchronousPlugin('sheetbend')
          plugin.container.inputData.XYZIN = self.container.inputData.BUCCANEER_MR_MODE_XYZIN
          plugin.container.inputData.F_SIGF = self.container.inputData.F_SIGF
          plugin.container.inputData.FREERFLAG = self.container.inputData.FREERFLAG
          self.processPlugin(plugin, 202, ['XYZOUT'], 203)
          model = plugin.container.outputData.XYZOUT
        else:
          model = self.container.inputData.BUCCANEER_MR_MODE_XYZIN
        plugin = self.makeRefmacPlugin(xyzin=model, cycles=10, localSymmetry=False, usePhases=False, useProsmart=False)
        self.processPlugin(plugin, 202, ['XYZOUT', 'ABCDOUT'], 203)
        self.container.inputData.ABCD = plugin.container.outputData.ABCDOUT
        self.container.inputData.ABCD.annotation.set('Initial MR phases')
        self.container.inputData.BUCCANEER_MR_MODE_XYZIN = plugin.container.outputData.XYZOUT
        self.container.inputData.BUCCANEER_MR_MODE_XYZIN.annotation.set('Refined molecular replacement solution')
        self.currentABCD = plugin.container.outputData.ABCDOUT

    def buccaneer(self):
        plugin = self.makeBuccaneerPlugin()
        self.processPlugin(plugin, 200, ['XYZOUT'], 201)
        xml = self.parseBuccaneerXML(plugin)
        self.appendXMLToCycle(xml)
        self.extractBuccaneerMetrics(xml)
        self.currentXYZ = plugin.container.outputData.XYZOUT

    def coot(self, plugin):
        self.processPlugin(plugin, 204)
        self.checkFileExists(plugin.container.outputData.XYZOUT[0], 205)
        plugin.container.outputData.XYZOUT[0].subType = 1
        self.currentXYZ = plugin.container.outputData.XYZOUT[0]

    def cootRSO(self):
        self.coot(self.makeCootRSOPlugin())

    def cootPrune(self, chains_only=False):
        self.coot(self.makeCootPrunePlugin(chains_only))

    def refmac(self, cycles=None, final=False):
        plugin = self.makeRefmacPlugin(cycles=cycles)
        self.processPlugin(plugin, 202, ['XYZOUT', 'ABCDOUT', 'FPHIOUT', 'DIFFPHIOUT'], 203)
        self.currentXYZ = plugin.container.outputData.XYZOUT
        self.currentABCD = plugin.container.outputData.ABCDOUT
        self.currentFPHI = plugin.container.outputData.FPHIOUT
        self.currentDIFFPHI = plugin.container.outputData.DIFFPHIOUT
        if final:
            xml = self.parseRefmacXML(plugin)
            self.appendXMLToCycle(xml)
            self.extractRefmacMetrics(xml)

# ===== Plugin creation methods ================================================

    def makeBuccaneerPlugin(self):
        plugin = self.makeSynchronousPlugin('buccaneer_mr')
        plugin.container.inputData.copyData(otherContainer=self.container.inputData,
                dataList=('F_SIGF', 'ASUIN', 'XYZIN_SEQ', ('BUCCANEER_MR_MODE_XYZIN', 'MR_MODE_XYZIN')))
        plugin.container.inputData.XYZIN = self.currentXYZ
        plugin.container.inputData.ABCD = self.currentABCD
        if self.container.controlParameters.XYZIN_MODE:
            plugin.container.inputData.XYZIN_MODE = True
            plugin.container.controlParameters.KNOWN_STRUCTURE = self.container.controlParameters.KNOWN_STRUCTURE
        if self.cycle > 0:
            plugin.container.inputData.XYZIN_MODE = True
            plugin.container.inputData.FWT_PHWT_IN = self.currentFPHI
            plugin.container.controlParameters.CYCLES = self.container.controlParameters.BUCCANEER_CYCLES_NEXT
            if self.essentiallyComplete: self.container.controlParameters.BUCCANEER_MR_MODE = 'nothing'
        else:
            plugin.container.controlParameters.CYCLES = self.container.controlParameters.BUCCANEER_CYCLES
        if self.container.controlParameters.BUCCANEER_USE_FREER:
            plugin.container.inputData.copyData(otherContainer=self.container.inputData, dataList=('FREERFLAG',))
        plugin.container.controlParameters.copyData(otherContainer=self.container.controlParameters,
            dataList=(('BUCCANEER_ANISOTROPY_CORRECTION', 'ANISOTROPY_CORRECTION'),
                      ('BUCCANEER_BUILD_SEMET', 'BUILD_SEMET'),
                      ('BUCCANEER_FAST', 'FAST'),
                      ('BUCCANEER_FIX_POSITION', 'FIX_POSITION'),
                      ('BUCCANEER_RESOLUTION', 'RESOLUTION'),
                      ('BUCCANEER_SEQUENCE_RELIABILITY', 'SEQUENCE_RELIABILITY'),
                      ('BUCCANEER_NEW_RESIDUE_NAME', 'NEW_RESIDUE_NAME'),
                      ('BUCCANEER_MODEL_SIGMA', 'MODEL_SIGMA'),
                      ('BUCCANEER_JOBS', 'JOBS'),
                      ('BUCCANEER_VERBOSE', 'VERBOSE'),
                      ('BUCCANEER_PHSIN_TYPE', 'PHSIN_TYPE'),
                      ('BUCCANEER_MR_MODE', 'MR_MODE'),
                      ('BUCCANEER_MR_MODE_SIGMA', 'MR_MODE_SIGMA'),
                      'F_SIGF_REF', 'ABCD_REF', 'XYZIN_REF'))
        return plugin

    def makeRefmacPlugin(self, xyzin=None, cycles=None, localSymmetry=None, usePhases=None, useProsmart=None):
        if xyzin == None: xyzin = self.currentXYZ
        if cycles == None: cycles = self.container.controlParameters.REFMAC_CYCLES
        if localSymmetry == None: localSymmetry = self.container.controlParameters.REFMAC_LOCAL_NCS and self.n_fragments < 14
        if usePhases == None: usePhases = self.refmacShouldUsePhases()
        if useProsmart == None: useProsmart = self.container.controlParameters.USE_PROSMART and self.container.inputData.TARGET.isSet()

        plugin = self.makeSynchronousPlugin('prosmart_refmac') if useProsmart else self.makeSynchronousPlugin('refmac')
        plugin.container.inputData.XYZIN = xyzin
        plugin.container.inputData.ABCD = self.container.inputData.ABCD if usePhases else None
        plugin.container.inputData.copyData(otherContainer=self.container.inputData, dataList=('F_SIGF', 'FREERFLAG', 'DICT'))
        if useProsmart: plugin.container.inputData.REFERENCE_MODEL = self.container.inputData.TARGET
        plugin.container.controlParameters.HYDR_USE = False
        plugin.container.controlParameters.NCYCLES = cycles
        plugin.container.controlParameters.USE_LOCAL_SYMMETRY = localSymmetry
        plugin.container.controlParameters.MAKE_NEW_LIGAND_EXIT = False
        plugin.container.controlParameters.copyData(otherContainer=self.container.controlParameters, dataList=('EXTRAREFMACKEYWORDS',))
        return plugin

    def refmacShouldUsePhases(self):
        if (self.container.controlParameters.BUCCANEER_PHSIN_TYPE == "mr" and
            self.container.controlParameters.REFMAC_MR_USEPHI):
            return True
        if (self.container.controlParameters.BUCCANEER_PHSIN_TYPE == "experimental" and
            self.container.controlParameters.REFMAC_EXP_USEPHI and
            self.rwork > self.container.controlParameters.BUCCANEER_MLHL_RWORK_LIMIT):
            return True
        return False

    def makeCootPlugin(self, script):
        plugin = self.makeSynchronousPlugin('coot_script_lines')
        def appendInput(inputList, item):
            inputList.append(inputList.makeItem())
            inputList[-1].set(item)
        appendInput(plugin.container.inputData.XYZIN, self.currentXYZ)
        appendInput(plugin.container.inputData.FPHIIN, self.currentFPHI)
        appendInput(plugin.container.inputData.DELFPHIIN, self.currentDIFFPHI)
        plugin.container.controlParameters.SCRIPT = script
        return plugin

    def makeCootRSOPlugin(self):
        rso = self.container.controlParameters.COOT_REALSPACE_OPERATION.__str__()
        if rso == "coot_script_lines":
            script = self.container.controlParameters.SCRIPT
        elif rso == "coot_fit_residues":
            script = (
                'fill_partial_residues(MolHandle_1)\n'
                'fit_protein(MolHandle_1)\n'
                'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n')
        elif rso == "coot_stepped_refine":
            if self.container.controlParameters.USERAMA:
                script = (
                    'fill_partial_residues(MolHandle_1)\n'
                    'stepped_refine_protein_for_rama(MolHandle_1)\n'
                    'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n')
            else:
                script = (
                    'fill_partial_residues(MolHandle_1)\n'
                    'stepped_refine_protein(MolHandle_1)\n'
                    'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n')
        elif rso == "coot_add_waters":
            script = (
                'set_ligand_water_to_protein_distance_limits(2.2, 3.3)\n'
                'execute_find_waters_real(MapHandle_1, MolHandle_1, 0, 1.3)\n'
                'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n')
        elif rso == "none":
            script = 'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n'
        return self.makeCootPlugin(script)

    def makeCootPrunePlugin(self, chains_only=False):
        current_path = os.path.abspath(inspect.getfile(inspect.currentframe()))
        script_path = os.path.join(os.path.dirname(current_path), "coot_prune.txt")
        script = open(script_path).read()
        if chains_only:
            script += 'prune(MolHandle_1, MapHandle_1, DifmapHandle_1, residues=False, sidechains=False)\n'
        else:
            script += 'prune(MolHandle_1, MapHandle_1, DifmapHandle_1)\n'
        script += 'write_pdb_file(MolHandle_1, os.path.join(dropDir, "output.pdb"))\n'
        return self.makeCootPlugin(script)

# ===== Generic sub-plugin methods =============================================

    def processPlugin(self, plugin, failedCode, requiredOutputFiles=[], missingFileCode=None):
        status = plugin.process()
        self.checkFinishStatus(status, failedCode)
        for outputFile in requiredOutputFiles:
            path = getattr(plugin.container.outputData, outputFile)
            self.checkFileExists(path, missingFileCode)
        if self.cycle in self.processedPlugins:
            self.processedPlugins[self.cycle].append(plugin)

    def makeSynchronousPlugin(self, name):
        plugin = self.makePluginObject(name)
        plugin.doAsync = False
        return plugin

# ===== XML methods ============================================================

    def initialiseXML(self):
        self.pipelinexmlfile = self.makeFileName('PROGRAMXML')
        if self.restarted and os.path.exists(self.pipelinexmlfile):
            self.xmlroot = etree.fromstring(open(self.pipelinexmlfile).read())
        else:
            self.xmlroot = etree.Element("BuccaneerBuildRefineResult")

    def addXMLCycle(self):
        self.xmlcyc = etree.SubElement(self.xmlroot, "BuildRefineCycle")
        etree.SubElement(self.xmlcyc, "Number").text = str(self.cycle + 1)

    def appendXMLToRoot(self, xml):
        self.xmlroot.append(xml)
        self.writeXMLRoot()

    def appendXMLToCycle(self, xml):
        self.xmlcyc.append(xml)
        self.writeXMLRoot()

    def writeFinalXML(self):
        xml = etree.Element("FinalStatistics")
        etree.SubElement(xml, "BestCycle").text = str(self.best_cycle + 1)
        cycle = self.xmlroot.findall("BuildRefineCycle")[self.best_cycle]
        for node in cycle.find("BuccaneerResult").find("Final"):
            xml.append(deepcopy(node))
        for node in cycle.find("RefmacResult"):
            xml.append(deepcopy(node))
        self.appendXMLToRoot(xml)

    def writeXMLRoot(self):
        with open(self.pipelinexmlfile, 'w') as f:
            CCP4Utils.writeXML(f, etree.tostring(self.xmlroot, pretty_print=True))

    def parseBuccaneerXML(self, plugin):
        return CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))

    def parseRefmacXML(self, plugin):
        raw_xml = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
        stats = raw_xml.xpath("//REFMAC/Overall_stats/stats_vs_cycle")
        new_xml = etree.Element('RefmacResult')
        for node in stats[0].xpath("new_cycle[last()]/r_factor | new_cycle[last()]/r_free | new_cycle[last()]/rmsBOND |  new_cycle[last()]/rmsANGLE"):
            node.text = str(node.text).strip()
            if node.tag == 'rmsBOND':
                node.text = str(100*float(node.text))
                node.tag = 'rmsBONDx100'
            new_xml.append(deepcopy(node))
        return new_xml

    def parseXmlFromBeforeInterrupt(self):
        cycles = self.xmlroot.xpath("//BuccaneerBuildRefineResult/BuildRefineCycle")
        for cycle in cycles:
            bxml = cycle.xpath("//BuildRefineCycle/BuccaneerResult")[-1]
            self.extractBuccaneerMetrics(bxml)
            rxml = cycle.xpath("//BuildRefineCycle/RefmacResult")[-1]
            self.extractRefmacMetrics(rxml)
            self.checkForImprovement()

    def extractBuccaneerMetrics(self, xml):
        self.completeness_by_res = float(xml.xpath('//BuccaneerResult/Final/CompletenessByResiduesBuilt')[-1].text)
        self.completeness_by_chn = float(xml.xpath('//BuccaneerResult/Final/CompletenessByChainsBuilt')[-1].text)
        self.n_fragments = int(xml.xpath('//BuccaneerResult/Final/FragmentsBuilt')[-1].text)
        self.longest_fragment = int(xml.xpath('//BuccaneerResult/Final/ResiduesLongestFragment')[-1].text)
        self.residues_built = int(xml.xpath('//BuccaneerResult/Final/ResiduesBuilt')[-1].text)
        self.residues_sequenced = int(xml.xpath('//BuccaneerResult/Final/ResiduesSequenced')[-1].text)
        print('Number of fragments:', self.n_fragments, " Completeness:", self.completeness_by_res)
        if self.n_fragments == 0: self.reportStatus(CPluginScript.UNSATISFACTORY)

    def extractRefmacMetrics(self, xml):
        self.rwork = float(xml.xpath('//RefmacResult/r_factor')[-1].text)
        self.rfree = float(xml.xpath('//RefmacResult/r_free')[-1].text)
        print('R-work:', self.rwork, " R-free:", self.rfree)
        if self.rwork < 0.30: self.essentiallyComplete = True

# ===== Error handling =========================================================

    def checkFinishStatus(self, status, errorCode):
        if status == CPluginScript.FAILED:
            self.appendErrorReport(errorCode)
            self.reportStatus(status)

    def checkFileExists(self, path, errorCode):
        if not path.exists():
            self.appendErrorReport(errorCode, 'Expected file: ' + str(path))
            self.reportStatus(CPluginScript.FAILED)

# ===== Utility methods ========================================================

    def copyPipelineOutputFiles(self):
        shutil.copyfile(str(self.currentXYZ), str(self.container.outputData.XYZOUT.fullPath))
        shutil.copyfile(str(self.currentABCD), str(self.container.outputData.ABCDOUT.fullPath))
        shutil.copyfile(str(self.currentFPHI), str(self.container.outputData.FPHIOUT.fullPath))
        shutil.copyfile(str(self.currentDIFFPHI), str(self.container.outputData.DIFFPHIOUT.fullPath))
        self.container.outputData.XYZOUT.annotation.set('Model built by Autobuild protein')
        self.container.outputData.XYZOUT.subType.set(1)
        self.container.outputData.FPHIOUT.annotation.set('2mFo-DFc map coefficients')
        self.container.outputData.FPHIOUT.subType.set(1)
        self.container.outputData.DIFFPHIOUT.annotation.set('mFo-DFc map coefficients')
        self.container.outputData.DIFFPHIOUT.subType.set(2)

    def setPerformanceData(self):
        self.container.outputData.PERFORMANCE.RFactor.set(self.rwork)
        self.container.outputData.PERFORMANCE.completeness.set(self.completeness_by_res)

    def alignment(self):
        try:
            print('Trying to create an alignment between output coordinates and input sequence.')
            match = CCP4ModelData.CChainMatch(self.container.outputData.XYZOUT, self.container.inputData.ASUIN)
            self.appendXMLToRoot(match.reportXmlAlignments())
            print('Alignment succeeded:', match)
        except Exception as e:
            print('Alignment failed:', e)

    def cleanUpIntermediateFiles(self):
        if not self.container.controlParameters.BUCCANEER_CLEANUP: return
        if self.cycle - 2 not in self.processedPlugins: return
        for plugin in self.processedPlugins[self.cycle - 2]:
            cleanup = CCP4ProjectsManager.CPurgeProject(plugin._dbProjectId)
            cleanup.purgeJob(plugin.jobId, context="extended_intermediate", reportMode="skip")

# ===== Unit testing ===========================================================

import unittest

class test_buccaneer_build_refine_mr(unittest.TestCase):

    def setUp(self):
        # make all background jobs wait for completion
        from core.CCP4Modules import QTAPPLICATION, PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)

    def tearDown(self):
        from core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)

    def test_1(self):
        from core.CCP4Modules import QTAPPLICATION
        from core.CCP4Utils import getCCP4I2Dir

        # Run the pipeline
        wrapper = buccaneer_build_refine_mr(parent=QTAPPLICATION(), name='buccaneer_build_refine_mr')
        wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(), 'pipelines', 'buccaneer_build_refine_mr', 'test_data', 'test_1.params.xml'))
        # Ensure no output file exists already
        xyzout = wrapper.container.outputData.XYZOUT.fullPath.get()
        if xyzout is not None and os.path.exists(xyzout):
            os.remove(xyzout)
        xmlout = wrapper.makeFileName('PROGRAMXML')
        if xmlout is not None and os.path.exists(xmlout):
            os.remove(xmlout)
        wrapper.process()

        # test if output file created
        self.assertEqual(os.path.exists(xyzout), 1, 'Failed to create copied pdb file '+xyzout)

def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(test_buccaneer_build_refine_mr)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
